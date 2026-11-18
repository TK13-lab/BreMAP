"""Explicit local teacher inference; optional upstream dependencies, resumable cache."""
import hashlib
import importlib.metadata
import importlib.util
import inspect
import json
import os
from pathlib import Path
import random
import tempfile
import time

import numpy as np
import pandas as pd

from .preprocessing import file_hash
from .student import _counts, _digest


class TeacherEncoder:
    """Generate dataset-specific embeddings from a local pretrained checkpoint.

    scGPT uses CLS with the official binned-count inference helper. Geneformer
    uses mean cell pooling at layer -1, requires Ensembl IDs and explicit V1/V2.
    The constructor uses local files; ``from_pretrained`` can download an official
    checkpoint after confirmation. No automatic GPU-to-CPU fallback.
    A cache belongs to one exact input/configuration; use another directory when
    changing either. Do not run concurrent writers against the same cache.
    """
    def __init__(self,model,checkpoint,*,device='cuda',batch_size=32,chunk_size=256,
                 max_length=None,gene_column=None,gene_case='exact',model_version=None,
                 seed=0,workers=2):
        if model not in ('scgpt','geneformer'): raise ValueError('model must be scgpt or geneformer')
        if device not in ('cpu','cuda'): raise ValueError('device must explicitly be cpu or cuda')
        if model=='geneformer' and model_version not in ('V1','V2'):
            raise ValueError('Geneformer requires model_version=V1 or V2 matching its checkpoint')
        if model=='geneformer' and device!='cuda': raise ValueError('This Geneformer adapter requires device=cuda')
        if gene_case not in ('exact','insensitive'): raise ValueError('gene_case must be exact or insensitive')
        if model=='geneformer' and gene_case!='exact': raise ValueError('Geneformer requires exact Ensembl IDs')
        if max_length is None: max_length=1200 if model=='scgpt' else (2048 if model_version=='V1' else 4096)
        if model=='geneformer' and model_version=='V1' and max_length!=2048:
            raise ValueError('Official Geneformer V1 tokenizer fixes max_length=2048')
        for name,v in [('batch_size',batch_size),('chunk_size',chunk_size),('max_length',max_length),('workers',workers)]:
            if isinstance(v,bool) or not isinstance(v,int) or v<1: raise ValueError(f'{name} must be a positive integer')
        self.model=model;self.checkpoint=Path(checkpoint);self.device=device
        self.options=dict(model=model,device=device,batch_size=batch_size,chunk_size=chunk_size,
            max_length=max_length,gene_column=gene_column,gene_case=gene_case,
            model_version=model_version,seed=seed,workers=workers)
        self._bundle=None
        self._bundle_identity=None

    @classmethod
    def from_pretrained(cls, checkpoint_id, *, cache_dir=None, confirm=None, **kwargs):
        """Select/download an official checkpoint, then construct a teacher.

        See :func:`bremap.list_checkpoints` for IDs. Confirmation and verified
        cache reuse follow :func:`bremap.download_checkpoint`. Family and
        Geneformer V1/V2 are inferred from the selected checkpoint; conflicting
        options fail before downloading. Other options are constructor arguments.
        Upstream inference libraries still need to be installed separately.
        """
        from .checkpoints import _select, download_checkpoint
        if checkpoint_id is None:
            raise ValueError('from_pretrained requires a checkpoint_id; use list_checkpoints()')
        entry = _select(checkpoint_id)
        version = kwargs.pop('model_version', entry['model_version'])
        if version != entry['model_version']:
            raise ValueError('model_version conflicts with the selected checkpoint')
        teacher = cls(entry['model'], '.', model_version=version, **kwargs)
        teacher.checkpoint = download_checkpoint(checkpoint_id, cache_dir=cache_dir, confirm=confirm)
        return teacher

    def _manifest(self):
        names=['args.json','vocab.json','best_model.pt'] if self.model=='scgpt' else ['config.json']
        if self.model=='geneformer':
            weights=list(self.checkpoint.glob('*.safetensors'))
            if not weights: weights=list(self.checkpoint.glob('pytorch_model*.bin'))
            if not weights: raise ValueError('Geneformer checkpoint lacks model weights')
            names += [p.name for p in sorted(weights)]
            cfg=json.loads((self.checkpoint/'config.json').read_text())
            capacity=cfg.get('max_position_embeddings')
            expected=2048 if self.options['model_version']=='V1' else 4096
            if capacity!=expected or self.options['max_length']>expected:
                raise ValueError('Geneformer model_version/max_length conflicts with checkpoint position capacity')
        hashes={}
        for name in names:
            p=self.checkpoint/name
            if not p.is_file(): raise ValueError(f'Missing local {self.model} checkpoint file: {p}')
            hashes[name]=file_hash(p)
        try: version=importlib.metadata.version(self.model)
        except importlib.metadata.PackageNotFoundError: version='not-installed'
        if self.model=='scgpt':
            spec=importlib.util.find_spec('scgpt')
            if spec and spec.origin:
                folder=Path(spec.origin).parent
                for name in ['tasks/cell_emb.py','data_collator.py','preprocess.py','model/model.py','utils/util.py']:
                    if (folder/name).is_file(): hashes['backend/'+name]=file_hash(folder/name)
        if self.model=='geneformer':
            try:
                import geneformer
                from geneformer import tokenizer
            except ImportError as e: raise ImportError('Install the official Geneformer package and its tokenizer assets in this environment') from e
            for key in ('GENE_MEDIAN_FILE','TOKEN_DICTIONARY_FILE','ENSEMBL_MAPPING_FILE'):
                suffix='_30M' if self.options['model_version']=='V1' else ''
                path=Path(getattr(geneformer,key+suffix)); hashes['tokenizer/'+key+suffix]=file_hash(path)
            hashes['tokenizer/source']=file_hash(inspect.getfile(tokenizer))
            folder=Path(inspect.getfile(tokenizer)).parent
            for name in ['emb_extractor.py','perturber_utils.py']:
                hashes['backend/'+name]=file_hash(folder/name)
        return {'checkpoint_files':hashes,'backend_version':version,'adapter_sha256':file_hash(__file__)}

    def transform(self,data,*,cache_dir):
        """Return cell-indexed DataFrame; completed chunks survive interruptions."""
        _counts(data,'exact')
        col=self.options['gene_column']
        genes=data.var_names if col is None else data.var[col] if col in data.var else None
        if genes is None: raise ValueError(f'gene_column {col!r} is missing from data.var')
        if pd.isna(genes).any() or any(not str(x).strip() for x in genes): raise ValueError('Missing gene identifiers')
        if self.model=='geneformer' and not all(str(g).startswith('ENSG') and '.' not in str(g) for g in genes):
            raise ValueError('Geneformer requires unversioned human Ensembl IDs; provide gene_column with an explicit mapping')
        model_manifest=self._manifest()
        identity=hashlib.sha256(json.dumps({'options':self.options,**model_manifest},sort_keys=True).encode()).hexdigest()
        if self._bundle_identity!=identity:
            self._bundle=None; self._bundle_identity=identity
        record={'input_sha256':_digest(data),'gene_identifiers':list(map(str,genes)),
                'options':self.options,**model_manifest}
        stamp=hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()
        folder=Path(cache_dir);folder.mkdir(parents=True,exist_ok=True); manifest=folder/'manifest.json'
        if manifest.exists():
            if json.loads(manifest.read_text())['fingerprint']!=stamp:
                raise ValueError('Teacher cache fingerprint differs: input/checkpoint/options/code changed; use a new cache_dir')
        else:
            with manifest.open('x') as f: json.dump(dict(fingerprint=stamp,**record),f,indent=2)
        arrays=[]; timings=[]; step=self.options['chunk_size']
        for start in range(0,data.n_obs,step):
            path=folder/f'chunk_{start:06d}.npz'; ids=data.obs_names[start:start+step]
            if path.exists():
                try:
                    with np.load(path,allow_pickle=False) as a:
                        values=a['values']; saved_ids=a['cell_ids']
                        if str(a['fingerprint'])!=stamp or not np.array_equal(saved_ids,ids): raise ValueError('identities')
                        if hashlib.sha256(values.tobytes()).hexdigest()!=str(a['values_sha256']): raise ValueError('hash')
                except Exception as e: raise ValueError(f'Invalid teacher cache chunk {path}; remove only this corrupted chunk to recompute') from e
            else:
                begin=time.perf_counter()
                values=np.asarray(self._encode_chunk(data[start:start+step].copy(),self.options['seed']+start),dtype=np.float32)
                self._validate(values,len(ids))
                temp=path.with_suffix('.tmp')
                with temp.open('wb') as f:
                    np.savez_compressed(f,values=values,cell_ids=np.asarray(ids,dtype=str),fingerprint=stamp,
                                        values_sha256=hashlib.sha256(values.tobytes()).hexdigest())
                temp.replace(path);timings.append({'start':start,'seconds':time.perf_counter()-begin})
            self._validate(values,len(ids));arrays.append(values)
        values=np.concatenate(arrays)
        result=pd.DataFrame(values,index=data.obs_names,columns=[f'embedding_{i}' for i in range(values.shape[1])])
        result.attrs['bremap_encoder_metadata']={'model':self.model,'checkpoint':identity,
            'preprocessing':json.dumps(self.options,sort_keys=True),'fingerprint':stamp}
        (folder/'completed.json').write_text(json.dumps({'fingerprint':stamp,'shape':list(values.shape),
            'new_chunks':timings,'values_sha256':hashlib.sha256(values.tobytes()).hexdigest()},indent=2))
        return result

    def attach(self,data,*,cache_dir,key='X_teacher'):
        """Embed and register vectors/provenance for reference or query mapping."""
        from .usability import attach_embeddings
        if key in data.obsm: raise ValueError(f'obsm[{key!r}] already exists; choose a new key')
        vectors=self.transform(data,cache_dir=cache_dir);meta=vectors.attrs['bremap_encoder_metadata']
        return attach_embeddings(data,vectors.rename_axis('cell_id').reset_index(),key=key,
            model=meta['model'],checkpoint=meta['checkpoint'],preprocessing=meta['preprocessing'])

    @staticmethod
    def _validate(values,n):
        if values.ndim!=2 or values.shape[0]!=n or values.shape[1]<1 or not np.isfinite(values).all():
            raise ValueError('Teacher returned nonfinite or mismatched embedding vectors')

    def _encode_chunk(self,data,seed):
        try: import torch
        except ImportError as e: raise ImportError('Teacher inference requires PyTorch and the official scGPT/Geneformer package') from e
        if self.device=='cuda' and not torch.cuda.is_available(): raise RuntimeError('CUDA requested but unavailable; check GPU access/environment')
        old_np=np.random.get_state();old_random=random.getstate()
        affinity=os.sched_getaffinity(0) if hasattr(os,'sched_getaffinity') else None
        try:
            np.random.seed(seed);random.seed(seed)
            if affinity: os.sched_setaffinity(0,sorted(affinity)[:self.options['workers']])
            with torch.random.fork_rng(devices=[torch.cuda.current_device()] if self.device=='cuda' else []):
                torch.random.default_generator.manual_seed(seed)
                if self.device=='cuda': torch.cuda.manual_seed(seed)
                return self._scgpt(data) if self.model=='scgpt' else self._geneformer(data)
        finally:
            np.random.set_state(old_np);random.setstate(old_random)
            if affinity: os.sched_setaffinity(0,affinity)

    def _scgpt(self,data):
        try:
            import torch
            from scgpt.model import TransformerModel
            from scgpt.tokenizer import GeneVocab
            from scgpt.tasks.cell_emb import get_batch_cell_embeddings
            from scgpt.utils import load_pretrained
        except ImportError as e: raise ImportError('Install the official scGPT package in this Python environment') from e
        if not hasattr(os,'sched_getaffinity'): raise RuntimeError('This scGPT adapter uses upstream Linux affinity API; run in Linux/WSL')
        if self._bundle is None:
            vocab=GeneVocab.from_file(self.checkpoint/'vocab.json')
            for token in ('<pad>','<cls>','<eoc>'):
                if token not in vocab: vocab.append_token(token)
            vocab.set_default_index(vocab['<pad>'])
            cfg=json.loads((self.checkpoint/'args.json').read_text())
            model=TransformerModel(ntoken=len(vocab),d_model=cfg['embsize'],nhead=cfg['nheads'],
                d_hid=cfg['d_hid'],nlayers=cfg['nlayers'],nlayers_cls=cfg['n_layers_cls'],n_cls=1,
                vocab=vocab,dropout=cfg['dropout'],pad_token=cfg['pad_token'],pad_value=cfg['pad_value'],
                do_mvc=True,do_dab=False,use_batch_labels=False,domain_spec_batchnorm=False,
                explicit_zero_prob=False,use_fast_transformer=False,pre_norm=False)
            load_pretrained(model,torch.load(self.checkpoint/'best_model.pt',map_location='cpu',weights_only=True),verbose=False)
            model.to(self.device);model.eval();self._bundle=(vocab,cfg,model)
        vocab,cfg,model=self._bundle
        col=self.options['gene_column']; genes=list(map(str,data.var_names if col is None else data.var[col]))
        if self.options['gene_case']=='insensitive':
            lookup={}
            for g in vocab.get_itos():
                key=g.upper()
                if key in lookup: raise ValueError('scGPT vocabulary has ambiguous case-insensitive gene IDs')
                lookup[key]=g
            genes=[lookup.get(g.upper(),g) for g in genes]
        if len(set(genes))!=len(genes): raise ValueError('Duplicate gene IDs after scGPT identifier matching')
        keep=np.array([g in vocab for g in genes]);data=data[:,keep].copy()
        if not keep.any() or (np.asarray(data.X.sum(1)).ravel()<=0).any(): raise ValueError('No expressed vocabulary genes in one or more cells; check gene identifiers')
        ids=np.array(vocab([g for g,k in zip(genes,keep) if k]),dtype=int)
        return get_batch_cell_embeddings(data,cell_embedding_mode='cls',model=model,vocab=vocab,
            max_length=self.options['max_length'],batch_size=self.options['batch_size'],model_configs=cfg,
            gene_ids=ids,use_batch_labels=False)

    def _geneformer(self,data):
        try: from geneformer import TranscriptomeTokenizer, EmbExtractor
        except ImportError as e: raise ImportError('Install the official Geneformer package with tokenizer resources') from e
        col=self.options['gene_column']; data.var['ensembl_id']=list(data.var_names if col is None else data.var[col])
        if data.var.ensembl_id.duplicated().any(): raise ValueError('Duplicate Ensembl IDs; resolve counts before tokenization')
        data.obs=pd.DataFrame({'cell_id':data.obs_names,'n_counts':np.asarray(data.X.sum(1)).ravel(),
                              'filter_pass':1},index=data.obs_names)
        version=self.options['model_version']
        with tempfile.TemporaryDirectory(prefix='bremap-geneformer-') as tmp:
            p=Path(tmp);src=p/'input';src.mkdir();data.write_h5ad(src/'counts.h5ad')
            tokenizer=TranscriptomeTokenizer(custom_attr_name_dict={'cell_id':'cell_id'},nproc=self.options['workers'],
                model_version=version,model_input_size=self.options['max_length'],special_token=version=='V2')
            tokenizer.tokenize_data(src,p,'cells',file_format='h5ad')
            extractor=EmbExtractor(model_type='Pretrained',emb_mode='cell',max_ncells=None,emb_layer=-1,
                emb_label=['cell_id'],forward_batch_size=self.options['batch_size'],nproc=self.options['workers'],model_version=version)
            df=extractor.extract_embs(str(self.checkpoint),str(p/'cells.dataset'),str(p),'vectors').set_index('cell_id')
            if not df.index.is_unique or set(df.index)!=set(data.obs_names): raise ValueError('Geneformer dropped/duplicated cells; check mapped gene coverage')
            return df.loc[data.obs_names].to_numpy(dtype=np.float32)
