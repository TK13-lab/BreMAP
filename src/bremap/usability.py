"""Explicit input preparation and convenient access to workflow outputs."""
from pathlib import Path
from importlib.resources import files
import json
import anndata as ad
import numpy as np
import pandas as pd

from .inputs import prepare_input
from .configuration import inspect_input
from .specification import ReferenceSpec
from .workflow import run_workflow


def read_data(source, *, metadata_columns=None, cell_id='cell_id', counts_layer=None):
    """Read raw-count H5AD/AnnData or a CSV with explicitly declared metadata columns.

    CSV rows are cells; remaining columns are genes. No normalization, donor
    splitting, gene conversion or missing-value imputation is performed.
    """
    if isinstance(source,(str,Path)):
        source=Path(source)
        if not source.is_file():
            raise FileNotFoundError(f"Input file not found: {source}. Supply an existing local count file.")
        if source.suffix.lower() not in ('.h5ad','.csv','.tsv'):
            raise ValueError('Unsupported input file format. Use .h5ad, .csv or .tsv; Excel and normalized matrices are not accepted.')
    if isinstance(source, ad.AnnData):
        data = source.copy()
    elif isinstance(source, (str, Path)) and Path(source).suffix.lower() == '.h5ad':
        data = ad.read_h5ad(source)
    else:
        if not metadata_columns:
            raise ValueError('CSV requires explicit metadata_columns')
        frame = pd.read_csv(source, sep='\t' if isinstance(source,Path) and source.suffix.lower()=='.tsv' else ',', dtype={cell_id: str})
        required = [cell_id, *metadata_columns]
        if len(set(required)) != len(required) or not set(required) <= set(frame):
            raise ValueError('Missing or repeated cell ID/metadata columns')
        if frame[cell_id].isna().any() or frame[cell_id].duplicated().any():
            raise ValueError('Cell IDs must be present and unique')
        frame = frame.set_index(cell_id)
        counts = frame.drop(columns=metadata_columns)
        if not len(counts.columns):
            raise ValueError('No gene columns')
        data = ad.AnnData(counts.to_numpy(dtype=float), obs=frame[list(metadata_columns)].copy(),
                         var=pd.DataFrame(index=counts.columns))
    if counts_layer is not None:
        if counts_layer not in data.layers:
            raise ValueError('Missing counts layer: ' + counts_layer)
        data.X = data.layers[counts_layer].copy()
    if data.X is None:
        raise ValueError('No count matrix in X; select an existing raw layer with counts_layer=')
    values = data.X.data if hasattr(data.X, 'tocsr') else np.asarray(data.X)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.rint(values)).all():
        raise ValueError('Raw finite nonnegative integer counts required')
    if not data.n_obs or not data.n_vars or not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError('Nonempty data with unique cell and gene IDs required')
    return data


def load_example():
    """Return the bundled synthetic counts as AnnData; query truth is not included."""
    with files('bremap').joinpath('resources/demo_counts.csv').open() as stream:
        return read_data(stream, metadata_columns=['role', 'donor_id', 'state'])


def example_annotations(destination=None):
    """Return separate synthetic query labels for evaluation, never model fitting."""
    with files('bremap').joinpath('resources/demo_truth.csv').open() as stream:
        table = pd.read_csv(stream)
    if destination is None:
        return table
    with Path(destination).open('x') as output:
        table.to_csv(output, index=False)
    return Path(destination)


def describe_data(data, *, role_column='role', donor_column='donor_id'):
    """Summarize cell and donor counts per declared role."""
    if role_column not in data.obs or donor_column not in data.obs:
        raise ValueError('Missing role/donor column')
    return data.obs.groupby(role_column, observed=True).agg(
        cells=(donor_column, 'size'), donors=(donor_column, 'nunique'))


def prepare_partitions(data, spec, destination, *, role_column='role',
                       donor_column='donor_id', label_column='state'):
    """Validate explicit donor-disjoint roles and prepare three immutable inputs.

    This does not choose a split. Query labels are removed from prepared inputs.
    """
    spec = spec if isinstance(spec, ReferenceSpec) else ReferenceSpec.load(spec)
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Partition destination already exists')
    for key in (role_column, donor_column, label_column):
        if key not in data.obs:
            raise ValueError('Missing metadata column: '+key)
    roles = ['reference', 'calibration', 'query']
    if data.obs[role_column].isna().any() or set(data.obs[role_column]) != set(roles):
        raise ValueError('Declare exactly reference, calibration and query roles')
    if data.obs[donor_column].isna().any():
        raise ValueError('Missing donor identifiers')
    sets = [set(data.obs.loc[data.obs[role_column] == r, donor_column].astype(str)) for r in roles]
    if any(sets[i] & sets[j] for i in range(3) for j in range(i)):
        raise ValueError('Donors overlap between partitions')
    # Stage all files first; invalid final partitions leave no destination.
    import tempfile
    import shutil
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as tmp:
        for role in roles:
            subset = data[data.obs[role_column] == role].copy()
            if role == 'query':
                subset.obs = subset.obs.drop(columns=list({label_column, spec.data['label_column']} & set(subset.obs)))
            prepare_input(subset, spec, Path(tmp)/f'{role}.h5ad', role=role,
                          donor_column=donor_column, label_column=None if role=='query' else label_column)
        shutil.copytree(tmp, destination)
    return {role: destination/f'{role}.h5ad' for role in roles}


class WorkflowResult:
    """Read tables and visualize an existing completed workflow without refitting."""
    def __init__(self, path):
        self.path = Path(path)
        record = json.loads((self.path/'workflow.json').read_text())
        if record.get('status') != 'complete':
            raise ValueError('Workflow is not complete')

    def table(self, name='mapping'):
        """Read mapping, programs, donor_programs or evaluation as a DataFrame."""
        paths = {'mapping':'mapping/cells.csv', 'programs':'programs/cells.csv',
                 'donor_programs':'programs/donors.csv', 'evaluation':'evaluation/summary.csv'}
        if name not in paths:
            raise ValueError('Unknown table: '+name)
        return pd.read_csv(self.path/paths[name])

    def summary(self):
        """Cell/donor counts and cell-weighted versus donor-macro acceptance."""
        cells = self.table()
        return pd.Series({'cells':len(cells), 'donors':cells.donor_id.nunique(),
                          'accepted':int(cells.accepted.sum()),
                          'acceptance_cell_weighted':cells.accepted.mean(),
                          'acceptance_donor_macro':cells.groupby('donor_id').accepted.mean().mean()})

    def plot(self, kind, **kwargs):
        """Plot programs, affinity, coverage or confusion; return matplotlib Figure."""
        from .plotting import plot_result
        return plot_result(self, kind, **kwargs)


def run_analysis(inputs, spec, output, *, method='affinity', policy='normal_donor_v1', annotations=None):
    """Run the established workflow and return a table/plot accessor.

    Inputs are the mapping returned by prepare_partitions. Annotation paths remain
    separate from counts. Existing lower-level run_workflow still returns a Path.
    """
    recipe = dict(schema='bremap-workflow-v1', name='bremap-analysis', spec=spec,
                  **inputs, output=output, method=method, policy=policy)
    if annotations is not None:
        recipe['annotations'] = annotations
    return WorkflowResult(run_workflow(recipe))


def apply_reference(archive, query, output):
    """Apply a calibrated portable reference to a prepared query, without fitting.

    Query must meet the reference input contract; use inspect_input/prepare_input
    for raw external data. This does not perform biological eligibility selection.
    """
    return WorkflowResult(run_workflow(dict(schema='bremap-workflow-v1',
        name='reference-application', archive=archive, query=query, output=output)))


def compare_measurement(mapping, measurements, *, state, measurement):
    """Join a separate measurement and report Spearman correlation per donor.

    All mapped cells must have matched cell/donor IDs; missing numeric scores are
    counted explicitly. Rejected cells remain included. No p-values or pooled
    cell-level significance is inferred. Return joined data and donor diagnostics.
    """
    from scipy.stats import spearmanr
    labels = pd.read_csv(measurements) if not isinstance(measurements,pd.DataFrame) else measurements.copy()
    score='conditional_affinity::'+state
    keys=['cell_id','donor_id']
    if score not in mapping or not set(keys+[measurement])<=set(labels):
        raise ValueError('Missing declared state/measurement or cell/donor identifiers')
    if measurement in mapping.columns:
        raise ValueError('Measurement name collides with mapping columns')
    if mapping.cell_id.duplicated().any() or labels.cell_id.duplicated().any():
        raise ValueError('Duplicate cell identifiers')
    if labels[keys].isna().any().any() or mapping[keys].isna().any().any():
        raise ValueError('Missing identifiers')
    joined=mapping.merge(labels[keys+[measurement]],on=keys,how='left',validate='one_to_one',indicator=True)
    if not joined['_merge'].eq('both').all():
        raise ValueError('Measurements must cover all mapped cells with matching donors')
    joined=joined.drop(columns='_merge')
    numeric=joined[[score,measurement]].apply(pd.to_numeric,errors='raise')
    rows=[]
    for donor,index in joined.groupby('donor_id').groups.items():
        subset=numeric.loc[index];valid=np.isfinite(subset).all(axis=1);subset=subset.loc[valid]
        reason='too_few_cells' if len(subset)<3 else ('constant_values' if (subset.nunique()<2).any() else '')
        rho=np.nan if reason else float(spearmanr(subset[score],subset[measurement]).statistic)
        rows.append(dict(donor_id=donor,n_total=len(index),n_scored=len(subset),n_missing=int((~valid).sum()),spearman_rho=rho,reason=reason))
    return joined,pd.DataFrame(rows)


def attach_embeddings(data, embeddings, *, key, model, checkpoint, preprocessing, cell_id='cell_id'):
    """Attach precomputed CSV/TSV/DataFrame vectors by exact cell IDs to a copy.

    Model is scgpt, geneformer or custom. These labels record user-supplied
    provenance; BReMAP does not infer encoder identity or run encoder inference.
    The table must contain cell_id plus numeric dimensions. No positional matching,
    truncation, automatic imputation, overwritten keys or partial cell overlap.
    """
    if not isinstance(data,ad.AnnData):
        raise ValueError('data must be AnnData; load it with read_data first')
    if model not in ('scgpt','geneformer','custom'):
        raise ValueError('model must be scgpt, geneformer or custom; SVD uses counts without uploaded embeddings')
    for name,value in [('key',key),('checkpoint',checkpoint),('preprocessing',preprocessing)]:
        if not isinstance(value,str) or not value.strip():
            raise ValueError(f'{name} must be a nonempty string describing the supplied embeddings')
    if key in data.obsm:
        raise ValueError(f'obsm[{key!r}] already exists; choose a new key to avoid overwriting vectors')
    if isinstance(embeddings,pd.DataFrame):
        table=embeddings.copy()
    else:
        path=Path(embeddings)
        if not path.is_file():raise FileNotFoundError(f'Embedding file not found: {path}')
        if path.suffix.lower() not in ('.csv','.tsv'):
            raise ValueError('Embedding file format must be .csv or .tsv with explicit cell_id and numeric dimensions')
        table=pd.read_csv(path,sep='\t' if path.suffix.lower()=='.tsv' else ',',dtype={cell_id:str})
    if cell_id not in table or table[cell_id].isna().any() or table[cell_id].astype(str).duplicated().any():
        raise ValueError(f'Embedding {cell_id!r} column must contain unique, nonmissing cell IDs')
    table[cell_id]=table[cell_id].astype(str)
    if not data.obs_names.is_unique or set(table[cell_id]) != set(data.obs_names):
        missing=set(data.obs_names)-set(table[cell_id]);extra=set(table[cell_id])-set(data.obs_names)
        raise ValueError(f'Embedding cell IDs must exactly match input: {len(missing)} missing, {len(extra)} extra. Re-export vectors for the same cells.')
    table=table.set_index(cell_id).loc[data.obs_names]
    try:values=table.to_numpy(dtype=float)
    except (TypeError,ValueError) as exc:
        raise ValueError('Embedding dimensions must be numeric; remove metadata columns other than cell_id') from exc
    if not values.shape[1] or not np.isfinite(values).all():
        raise ValueError('Embedding dimensions must be nonempty and finite; NaN/Inf cannot be imputed automatically')
    result=data.copy();result.obsm[key]=values
    provenance=dict(model=model,checkpoint=checkpoint,preprocessing=preprocessing)
    result.uns.setdefault('bremap_encoder_metadata',{})[key]=provenance
    result.uns.setdefault('bremap_embedding_provenance',{})[key]=json.dumps(provenance,sort_keys=True)
    return result
