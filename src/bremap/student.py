"""Experimental, dataset-trained distilled encoders. Inference requires no torch."""
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from numbers import Real
import warnings

import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.metrics import f1_score

from .preprocessing import normalize_counts
from .expression import values as expression_values, features, effective_contract, check_contract, CANONICAL, zero_rows, get_declaration


def _ids(values, name):
    ids = pd.Index(values)
    if not ids.is_unique or ids.isna().any() or any(not str(x).strip() for x in ids):
        raise ValueError(f'{name} IDs must be unique, nonempty and nonmissing')
    return ids


def _aligned(frame, ids, name):
    if not isinstance(frame, (pd.DataFrame, pd.Series)):
        raise TypeError(f'{name} must be a pandas object indexed by exact cell IDs')
    _ids(frame.index, name + ' cell')
    if set(frame.index) != set(ids):
        raise ValueError(f'{name} cell IDs must match counts exactly; align/filter explicitly')
    return frame.loc[ids]


def _counts(data, gene_case):
    _ids(data.obs_names, 'cell')
    genes = _ids(data.var_names, 'gene')
    if gene_case not in ('exact', 'upper'):
        raise ValueError('gene_case must be exact or upper')
    if gene_case == 'upper': genes = _ids(genes.str.upper(), 'uppercase gene')
    x = expression_values(data).astype(np.float32)
    if not np.isfinite(x.data).all():
        raise ValueError('Expression exceeds float32 numeric range for the student')
    if not data.n_obs or not data.n_vars or zero_rows(x).any():
        raise ValueError('Empty input or zero-expression cells; remove them explicitly')
    return x, genes


def _features(data, x):
    # Preserve historical float32-count conversion before float64 normalization.
    return (normalize_counts(x) if get_declaration(data)['type']=='counts' else features(data)).astype(np.float32)


def _digest(data, matrix=None):
    x = sparse.csr_matrix(data.X if matrix is None else matrix).copy()
    x.sum_duplicates(); x.sort_indices(); x.eliminate_zeros()
    h = hashlib.sha256()
    h.update(json.dumps([data.obs_names.tolist(), data.var_names.tolist()]).encode())
    for part in (np.asarray(x.shape, dtype='<i8'), x.indptr.astype('<i8'),
                 x.indices.astype('<i8'), x.data.astype('<f8')):
        h.update(part.tobytes())
    return h.hexdigest()


@dataclass
class StudentEncoder:
    """Frozen feature selection + SVD + ReLU MLP predicting one teacher's vectors.

    Train with paired declared expression/vectors and donor-disjoint validation. This is an
    experimental surrogate, not a pretrained foundation model. Build a separate
    BReMAP reference with its vectors; never mix teacher and student partitions.
    """
    genes: np.ndarray
    components: np.ndarray
    input_mean: np.ndarray
    input_scale: np.ndarray
    weights: list
    biases: list
    target_mean: np.ndarray
    target_scale: np.ndarray
    metadata: dict

    @classmethod
    def fit(cls, train, teacher, *, validation, validation_teacher,
            teacher_provenance, donor_key='donor_id', n_genes=2048,
            n_components=128, hidden_size=128, epochs=40, batch_size=128,
            learning_rate=.001, loss='mse', relational_weight=.05,
            minimum_gene_fraction=.90, gene_case='exact', seed=0,
            validation_weights=None):
        """Fit on train only; select epoch by validation raw-vector weighted MSE.

        Teacher/validation_teacher are cell-indexed DataFrames produced with the
        SAME teacher configuration. Provenance must include model and checkpoint.
        Optional validation_weights is a positive cell-indexed Series. torch is
        optional and only imported here; training runs on CPU.
        """
        if loss not in ('mse', 'mse_relational'): raise ValueError('loss must be mse or mse_relational')
        for name, val in [('n_genes',n_genes),('n_components',n_components),
                          ('hidden_size',hidden_size),('epochs',epochs),('batch_size',batch_size)]:
            if isinstance(val,bool) or not isinstance(val,int) or val < 1:
                raise ValueError(f'{name} must be a positive integer')
        if not 0 < minimum_gene_fraction <= 1: raise ValueError('minimum_gene_fraction must be in (0,1]')
        if not np.isfinite(learning_rate) or learning_rate <= 0: raise ValueError('learning_rate must be positive')
        if not np.isfinite(relational_weight) or relational_weight < 0: raise ValueError('relational_weight must be nonnegative')
        if not isinstance(teacher_provenance,dict) or any(not teacher_provenance.get(k) for k in ('model','checkpoint')):
            raise ValueError('teacher_provenance requires model and checkpoint')
        for table in (teacher,validation_teacher):
            declared=getattr(table,'attrs',{}).get('bremap_encoder_metadata',{})
            for key in ('model','checkpoint','preprocessing'):
                if key in declared and key in teacher_provenance and declared[key]!=teacher_provenance[key]:
                    raise ValueError(f'Teacher provenance {key} differs from teacher_provenance; use one encoder configuration')
        left=getattr(teacher,'attrs',{}).get('bremap_encoder_metadata',{})
        right=getattr(validation_teacher,'attrs',{}).get('bremap_encoder_metadata',{})
        for key in ('model','checkpoint','preprocessing'):
            if key in left and key in right and left[key]!=right[key]:
                raise ValueError(f'Train/validation teacher provenance {key} differs')
        for a in (train, validation):
            if donor_key not in a.obs or a.obs[donor_key].isna().any() or a.obs[donor_key].astype(str).str.strip().eq('').any():
                raise ValueError(f'Provide nonmissing donor metadata in obs[{donor_key!r}]')
        if set(train.obs[donor_key].astype(str)) & set(validation.obs[donor_key].astype(str)):
            raise ValueError('Train/validation donor overlap; namespace studies and split donors first')
        if set(train.obs_names) & set(validation.obs_names): raise ValueError('Train/validation cell overlap')
        x, genes = _counts(train,gene_case)
        contract = effective_contract(train)
        check_contract(validation, contract)
        # Canonical column-major targets make column statistics independent of
        # DataFrame/CSV memory layout and preserve the historical CSV fit path.
        y = np.asfortranarray(_aligned(teacher, train.obs_names, 'teacher'), dtype=np.float32)
        vy = np.asfortranarray(_aligned(validation_teacher, validation.obs_names, 'validation teacher'),dtype=np.float32)
        if y.ndim != 2 or vy.ndim != 2 or y.shape[1] == 0 or y.shape[1] != vy.shape[1] or not np.isfinite(y).all() or not np.isfinite(vy).all():
            raise ValueError('Teacher vectors must be finite, 2-D and have matching dimensions')
        if not isinstance(teacher,pd.DataFrame) or not teacher.columns.equals(validation_teacher.columns):
            raise ValueError('Teacher embedding columns must match in the same order')
        if float(np.var(y,axis=0).sum()) < 1e-12: raise ValueError('Teacher vectors are constant; cannot distill variation')
        if train.n_obs < 3 or train.n_vars < 3: raise ValueError('At least three training cells/genes required')
        log = _features(train, x)
        variance = np.asarray(log.power(2).mean(0)).ravel()-np.asarray(log.mean(0)).ravel()**2
        ix = np.argsort(-variance,kind='stable')[:min(n_genes,len(genes))]
        svd = TruncatedSVD(n_components=min(n_components,len(ix)-1,train.n_obs-1),random_state=seed)
        z = svd.fit_transform(log[:,ix])
        meta = dict(format_version=1, experimental=True, gene_case=gene_case,
                    minimum_gene_fraction=minimum_gene_fraction, teacher=dict(teacher_provenance),
                    train_sha256=_digest(train), validation_sha256=_digest(validation),
                    teacher_sha256=hashlib.sha256(y.tobytes()).hexdigest(),
                    validation_teacher_sha256=hashlib.sha256(vy.tobytes()).hexdigest(),
                    train_donors=sorted(set(train.obs[donor_key].astype(str))),
                    validation_donors=sorted(set(validation.obs[donor_key].astype(str))),
                    parameters=dict(n_genes=n_genes,n_components=n_components,hidden_size=hidden_size,
                        epochs=epochs,batch_size=batch_size,learning_rate=learning_rate,loss=loss,
                        relational_weight=relational_weight,seed=seed), history=[])
        if contract['normalization'] != CANONICAL:
            meta['expression_contract'] = contract
        obj = cls(genes[ix].to_numpy(dtype=str),svd.components_,z.mean(0),np.maximum(z.std(0),1e-6),
                  [],[],y.mean(0),np.maximum(y.std(0),1e-6),meta)
        z = ((z-obj.input_mean)/obj.input_scale).astype(np.float32)
        vz = obj._project(validation)
        if validation_weights is None:
            counts = validation.obs[donor_key].astype(str).value_counts()
            vw = 1/validation.obs[donor_key].astype(str).map(counts).to_numpy(float)
        else:
            vw = np.asarray(_aligned(validation_weights,validation.obs_names,'validation weights'),float)
        if vw.shape != (validation.n_obs,) or not np.isfinite(vw).all() or (vw <= 0).any():
            raise ValueError('validation_weights must be finite positive cell weights')
        meta['validation_weights_sha256'] = hashlib.sha256(vw.tobytes()).hexdigest()
        try: import torch
        except ImportError as e: raise ImportError('Student training requires PyTorch; install bremap[train]. Inference does not.') from e
        # Local RNGs: do not reset the caller's global NumPy or torch RNG state.
        with torch.random.fork_rng(devices=[]):
            torch.random.default_generator.manual_seed(seed)
            net = torch.nn.Sequential(torch.nn.Linear(z.shape[1],hidden_size),torch.nn.ReLU(),
                                      torch.nn.Linear(hidden_size,y.shape[1]))
        optimizer = torch.optim.Adam(net.parameters(),lr=learning_rate)
        tz, ty = torch.from_numpy(z), torch.from_numpy((y-obj.target_mean)/obj.target_scale)
        tvz = torch.from_numpy(vz.astype(np.float32)); rng = np.random.default_rng(seed)
        best = float('inf')
        for epoch in range(epochs):
            order = rng.permutation(len(z)); net.train(); losses=[]
            for start in range(0,len(z),batch_size):
                batch = order[start:start+batch_size]; pred = net(tz[batch]); target = ty[batch]
                error = (pred-target).square().mean()
                if loss == 'mse_relational' and len(batch)>1:
                    error = error + relational_weight * (torch.pdist(pred)/np.sqrt(y.shape[1])-
                                                         torch.pdist(target)/np.sqrt(y.shape[1])).square().mean()
                optimizer.zero_grad(); error.backward(); optimizer.step(); losses.append(float(error.detach()))
            net.eval()
            with torch.no_grad(): pv = net(tvz).numpy()*obj.target_scale+obj.target_mean
            score = float(np.average(np.mean((pv-vy)**2,axis=1),weights=vw))
            if not np.isfinite(score): raise ValueError('Nonfinite training outcome; inspect target scales/learning rate')
            meta['history'].append({'epoch':epoch+1,'train_loss':float(np.mean(losses)),'validation_mse':score})
            if score < best:
                best=score; meta['selected_epoch']=epoch+1
                layers=[net[0],net[2]]
                obj.weights=[p.weight.detach().numpy().T.copy() for p in layers]
                obj.biases=[p.bias.detach().numpy().copy() for p in layers]
        return obj

    def _project(self, data):
        check_contract(data,self.metadata.get('expression_contract', {'normalization':CANONICAL}))
        x, genes = _counts(data,self.metadata['gene_case'])
        ix=genes.get_indexer(self.genes); present=ix>=0; fraction=float(present.mean())
        if fraction < self.metadata['minimum_gene_fraction']:
            raise ValueError(f'Student gene coverage {fraction:.3f} is below {self.metadata["minimum_gene_fraction"]}; check gene IDs/assay')
        if not present.all(): warnings.warn(f'{sum(~present)} student genes absent; zero-filled, coverage={fraction:.3f}',UserWarning,stacklevel=2)
        log=_features(data,x)
        selected=np.zeros((data.n_obs,len(ix)),dtype=np.float32)
        selected[:,present]=log[:,ix[present]].toarray()
        return (selected@self.components.T-self.input_mean)/self.input_scale

    def transform(self, data):
        """Return cell-indexed vectors. Use representation='custom' in BReMAP."""
        z=self._project(data)
        for i,(w,b) in enumerate(zip(self.weights,self.biases)):
            z=z@w+b
            if i<len(self.weights)-1: z=np.maximum(z,0)
        values=z*self.target_scale+self.target_mean
        if not np.isfinite(values).all(): raise ValueError('Nonfinite student predictions')
        out=pd.DataFrame(values,index=data.obs_names,columns=[f'embedding_{i}' for i in range(values.shape[1])])
        out.attrs['bremap_encoder_metadata']={'model':'custom','checkpoint':self.fingerprint(),
            'preprocessing':self._preprocessing(), 'teacher':self.metadata['teacher']}
        return out

    def _preprocessing(self):
        if 'expression_contract' in self.metadata:
            return 'student provided expression; '+json.dumps(self.metadata['expression_contract'],sort_keys=True)+'; '+self.metadata['gene_case']
        return 'student CP10k/log1p; '+self.metadata['gene_case']

    def _arrays(self):
        return dict(genes=self.genes,components=self.components,input_mean=self.input_mean,
            input_scale=self.input_scale,target_mean=self.target_mean,target_scale=self.target_scale,
            metadata=np.asarray(json.dumps(self.metadata,sort_keys=True,allow_nan=False)),layers=np.asarray(len(self.weights)),
            **{f'w{i}':w for i,w in enumerate(self.weights)},**{f'b{i}':b for i,b in enumerate(self.biases)})

    def attach(self,data,*,key='X_student'):
        """Return counts plus registered custom vectors, ready for BReMAP APIs."""
        from .usability import attach_embeddings
        if key in data.obsm: raise ValueError(f'obsm[{key!r}] already exists; choose a new key')
        vectors=self.transform(data)
        return attach_embeddings(data,vectors.rename_axis('cell_id').reset_index(),key=key,
            model='custom',checkpoint=self.fingerprint(),preprocessing=self._preprocessing())

    def fingerprint(self):
        h=hashlib.sha256()
        for k,v in sorted(self._arrays().items()): h.update(k.encode()); h.update(np.asarray(v).tobytes())
        return h.hexdigest()

    def save(self,path):
        """Write an immutable, pickle-free NPZ artifact."""
        with Path(path).open('xb') as f: np.savez_compressed(f,**self._arrays(),fingerprint=self.fingerprint())

    @classmethod
    def load(cls,path):
        with np.load(path,allow_pickle=False) as a:
            meta=json.loads(str(a['metadata']))
            if meta.get('format_version')!=1: raise ValueError('Unsupported student format_version')
            n=int(a['layers'])
            if n!=2: raise ValueError('Student must have two dense layers')
            obj=cls(a['genes'],a['components'],a['input_mean'],a['input_scale'],
                [a[f'w{i}'] for i in range(n)],[a[f'b{i}'] for i in range(n)],a['target_mean'],a['target_scale'],meta)
            if obj.fingerprint()!=str(a['fingerprint']): raise ValueError('Student integrity fingerprint mismatch')
        return obj


def evaluate_student(teacher,student,*,states,truth=None,agreement_threshold=.85,
                     maximum_coverage_change=.05,maximum_f1_loss=.05,
                     minimum_task_f1=None,minimum_coverage=None):
    """Donor-macro task fidelity; missing truth cannot yield joint acceptance.

    Inputs are cell-indexed mapping tables with donor_id, candidate_state, accepted,
    or CSV paths containing those columns and cell_id. truth is a cell-indexed
    Series or CSV with cell_id, donor_id, label; file donors are checked as well.
    F1 uses all cells and every declared state (absent states contribute zero),
    then averages donors equally. Coverage also uses all cells within each donor.
    joint_pass retains the relative teacher-fidelity checks. case_pass adds the
    optional absolute minimum_task_f1 and minimum_coverage; both are fractions.
    Without truth, joint_pass/case_pass are None. Choose thresholds before testing.
    The thresholds are engineering margins, not statistical noninferiority tests.
    """
    thresholds=dict(agreement_threshold=agreement_threshold,
                    maximum_coverage_change=maximum_coverage_change,maximum_f1_loss=maximum_f1_loss,
                    minimum_task_f1=minimum_task_f1,minimum_coverage=minimum_coverage)
    for key,v in thresholds.items():
        if v is None and key.startswith('minimum_'): continue
        if isinstance(v,(bool,np.bool_)) or not isinstance(v,Real) or not np.isfinite(v) or not 0<=v<=1:
            raise ValueError(f'{key} threshold must be a numeric fraction between 0 and 1')
    def read_table(value, required):
        if not isinstance(value,(str,Path)): return value
        table=pd.read_csv(value,dtype={'cell_id':str,'donor_id':str,'candidate_state':str,'label':str})
        missing={'cell_id',*required}-set(table)
        if missing: raise ValueError(f'{value}: missing required columns {sorted(missing)}')
        return table.set_index('cell_id')
    teacher=read_table(teacher,['donor_id','candidate_state','accepted'])
    student=read_table(student,['donor_id','candidate_state','accepted'])
    if not isinstance(teacher,pd.DataFrame) or not isinstance(student,pd.DataFrame):
        raise TypeError('teacher and student mappings must be cell-indexed DataFrames')
    _ids(teacher.index,'teacher cell'); student=_aligned(student,teacher.index,'student')
    if not states or len(set(states))!=len(states): raise ValueError('states must be a nonempty unique vocabulary')
    for a in (teacher,student):
        if not {'donor_id','candidate_state','accepted'}<=set(a): raise ValueError('Mapping requires donor_id, candidate_state, accepted')
        if a[['donor_id','candidate_state','accepted']].isna().any().any(): raise ValueError('Mapping contains missing values')
        if a.donor_id.astype(str).str.strip().eq('').any(): raise ValueError('Mapping donor_id must be nonempty')
        if a.accepted.dtype != bool: raise ValueError('accepted must contain booleans')
        if not set(a.candidate_state)<=set(states): raise ValueError('candidate_state outside declared states')
    if not teacher.donor_id.equals(student.donor_id): raise ValueError('Teacher/student donor metadata differ')
    if isinstance(truth,(str,Path)):
        labels=_aligned(read_table(truth,['donor_id','label']),teacher.index,'truth')
        if not labels.donor_id.equals(teacher.donor_id): raise ValueError('Truth/mapping donor metadata differ')
        truth=labels.label
    if truth is not None:
        if not isinstance(truth,pd.Series): raise TypeError('truth must be a cell-indexed Series of state labels')
        truth=_aligned(truth,teacher.index,'truth')
        if truth.isna().any() or not set(truth)<=set(states): raise ValueError('truth outside declared states')
    rows=[]
    for donor, t in teacher.groupby('donor_id',sort=True,observed=True):
        s=student.loc[t.index]
        row={'donor':str(donor),'n_cells':len(t),'agreement':float((t.candidate_state==s.candidate_state).mean()),
             'teacher_coverage':float(t.accepted.mean()),'student_coverage':float(s.accepted.mean())}
        if truth is not None:
            for label,p in [('teacher',t),('student',s)]:
                row[label+'_f1']=float(f1_score(truth.loc[p.index],p.candidate_state,labels=states,average='macro',zero_division=0))
        rows.append(row)
    if not rows: raise ValueError('Empty mapping tables')
    d=pd.DataFrame(rows); agreement=float(d.agreement.mean()); delta=float((d.student_coverage-d.teacher_coverage).mean())
    result={'agreement':agreement,'coverage_change':delta,'agreement_threshold':agreement_threshold,
            'agreement_pass':agreement+1e-12>=agreement_threshold,'coverage_pass':abs(delta)<=maximum_coverage_change+1e-12,
            'f1_loss':None,'f1_pass':None,'joint_pass':None,'donors':rows,
            'student_f1':None,'student_coverage':float(d.student_coverage.mean()),
            'task_pass':None,'case_pass':None,'thresholds':thresholds}
    if truth is not None:
        gap=float((d.teacher_f1-d.student_f1).mean()); result.update(f1_loss=gap,f1_pass=gap<=maximum_f1_loss+1e-12)
        result['joint_pass']=bool(result['agreement_pass'] and result['coverage_pass'] and result['f1_pass'])
        result['student_f1']=float(d.student_f1.mean())
        result['task_pass']=bool((minimum_task_f1 is None or result['student_f1']+1e-12>=minimum_task_f1)
            and (minimum_coverage is None or result['student_coverage']+1e-12>=minimum_coverage))
        result['case_pass']=bool(result['joint_pass'] and result['task_pass'])
    return result
