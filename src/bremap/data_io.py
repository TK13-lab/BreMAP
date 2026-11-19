"""Read explicitly selected expression and inspect processed input containers."""
from copy import deepcopy
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

from .expression import KEY, declaration, get_declaration, values


def read_data(source, *, metadata_columns=None, cell_id='cell_id', counts_layer=None,
              matrix=None, expression_type=None, log_base=None, processing_id=None,
              units=None, assay=None):
    """Load counts or explicitly declared processed expression into AnnData.

    ``matrix`` is X, raw.X or layers/<name> for AnnData, and an assay layer for
    RDS. ``counts_layer`` remains an alias for layers/<name>. Library-normalized
    and log1p inputs require a known per-cell scaling transform; integrated or
    residual values require preprocessed mode with processing_id and units.
    No normalization or gene conversion occurs here; source values stay intact.
    """
    if counts_layer is not None and matrix is not None:
        raise ValueError('Choose matrix or counts_layer, not both')
    selection = 'layers/'+counts_layer if counts_layer is not None else ('X' if matrix is None else matrix)
    if not isinstance(selection, str) or not selection.strip():
        raise ValueError('matrix must name X, raw.X or layers/<name>')
    if isinstance(source, (str, Path)):
        source = Path(source)
        if not source.is_file():
            raise FileNotFoundError(f'Input file not found: {source}. Supply an existing local file.')
        if source.suffix.lower() not in ('.h5ad', '.csv', '.tsv', '.rds'):
            raise ValueError('Unsupported input file format. Use .h5ad, .csv, .tsv or supported .rds objects.')
    is_rds = isinstance(source, Path) and source.suffix.lower() == '.rds'
    if assay is not None and not is_rds:
        raise ValueError('assay is for RDS; select AnnData expression with matrix or counts_layer')
    if is_rds:
        if counts_layer is not None:
            raise ValueError('RDS uses assay and matrix, not counts_layer')
        from .rds_io import read_rds
        data = read_rds(source, assay=assay, matrix=matrix)
    else:
        if isinstance(source, ad.AnnData):
            original = source
        elif isinstance(source, Path) and source.suffix.lower() == '.h5ad':
            original = ad.read_h5ad(source)
        else:
            if selection != 'X':
                raise ValueError('CSV/TSV expression has only matrix="X"')
            if not metadata_columns:
                raise ValueError('CSV requires explicit metadata_columns')
            frame = pd.read_csv(source, sep='\t' if isinstance(source, Path) and source.suffix.lower()=='.tsv' else ',', dtype={cell_id: str})
            required = [cell_id, *metadata_columns]
            if len(set(required)) != len(required) or not set(required) <= set(frame):
                raise ValueError('Missing or repeated cell ID/metadata columns')
            if frame[cell_id].isna().any() or frame[cell_id].duplicated().any():
                raise ValueError('Cell IDs must be present and unique')
            frame = frame.set_index(cell_id)
            counts = frame.drop(columns=metadata_columns)
            if not len(counts.columns):
                raise ValueError('No gene columns')
            original = ad.AnnData(counts.to_numpy(dtype=float), obs=frame[list(metadata_columns)].copy(),
                                 var=pd.DataFrame(index=counts.columns))
        if selection == 'raw.X':
            if original.raw is None:
                raise ValueError('Missing raw.X; inspect_data lists available matrices')
            data = ad.AnnData(original.raw.X.copy(), obs=original.obs.copy(),
                var=original.raw.var.copy(), uns=deepcopy(original.uns),
                obsm={k: v.copy() for k,v in original.obsm.items()})
        else:
            data = original.copy()
            if selection != 'X':
                if not selection.startswith('layers/') or selection[7:] not in data.layers:
                    raise ValueError('Missing matrix '+selection+'; use X, raw.X or layers/<name>')
                data.X = data.layers[selection[7:]].copy()
    if expression_type is None and KEY in data.uns and matrix is None and counts_layer is None:
        if log_base is not None or processing_id is not None or units is not None:
            raise ValueError('Existing expression declaration: supply expression_type explicitly to change log_base/processing_id/units')
    log_base = 'e' if log_base is None else log_base
    if expression_type is None:
        # A different selected matrix must not inherit X's declaration.
        expression_type = get_declaration(data)['type'] if matrix is None and counts_layer is None and KEY in data.uns else 'counts'
        if expression_type != 'counts' and KEY in data.uns:
            old = get_declaration(data)
            log_base = old.get('log_base', log_base)
            processing_id = old.get('processing_id', processing_id); units = old.get('units', units)
    info = declaration(expression_type, log_base=log_base, processing_id=processing_id, units=units)
    info['source'] = {'matrix': selection, 'container': str(source) if isinstance(source,Path) else 'AnnData/CSV'}
    if is_rds:
        info['source'].update(data.uns.pop('_bremap_rds_source'))
    data.uns[KEY] = info
    checked = values(data)
    if expression_type == 'counts' and not np.equal(checked.data, np.rint(checked.data)).all():
        raise ValueError('Raw integer counts required; declare the actual expression_type rather than round values')
    if not data.n_obs or not data.n_vars or not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError('Nonempty data with unique cell and gene IDs required')
    from .core import _identifiers
    _identifiers(data.obs_names, 'cell identifiers', data.n_obs)
    _identifiers(data.var_names, 'gene identifiers', data.n_vars)
    return data


def _matrix_report(name, matrix, shape=None, cells=None, genes=None):
    sample = matrix[:min(128, matrix.shape[0])]
    sample = sample.data if sparse.issparse(sample) else np.asarray(sample).ravel()
    finite = bool(np.isfinite(sample).all())
    return {'matrix': name, 'shape': list(matrix.shape if shape is None else shape),
            'sparse': sparse.issparse(matrix) or getattr(matrix, 'format', None) in ('csr','csc'),
            'dtype': str(matrix.dtype),
            'cell_ids_preview': [] if cells is None else list(map(str,cells[:5])),
            'gene_ids_preview': [] if genes is None else list(map(str,genes[:5])),
            'sample_rows': min(128, matrix.shape[0]), 'finite_sample': finite,
            'nonnegative_sample': bool(finite and (sample >= 0).all()),
            'integer_like_sample': bool(finite and np.equal(sample, np.rint(sample)).all())}


def _capabilities(report):
    columns = report['obs_columns']
    report['metadata_candidates'] = {
        'donor': [c for c in columns if any(t in c.lower() for t in ('donor','patient','sample','subject'))],
        'label': [c for c in columns if any(t in c.lower() for t in ('type','state','label','cluster'))],
    }
    report['operations_by_expression'] = {
        kind: {'mapping': True, 'programs': True, 'student': True,
               'teacher': kind != 'preprocessed',
               'requirements': ('Requires matching processing_id/units, a separately fitted reference/student; '
                   'supply compatible teacher targets for student fitting. Arbitrary residuals/integrated values cannot generate FM targets.'
                   if kind == 'preprocessed' else
                   'Requires declared expression semantics, complete validation, compatible genes and metadata; '
                   'teacher adapters additionally require upstream runtimes/checkpoints and model-specific gene IDs.')}
        for kind in ('counts','library_normalized','log1p_library_normalized','preprocessed')}
    report['notes'] += ' Operation availability is conditional, not a validation result. Metadata candidates are names only; select them explicitly.'
    return report


def inspect_data(source):
    """Inventory source matrices/metadata, without inferring normalization.

    H5AD expression is sampled from disk, not loaded in full. Sampling is not
    validation: read_data checks the complete selected matrix before analysis.
    RDS parsing may require the full in-memory object.
    """
    notes = ('Expression semantics require an explicit declaration. Integer-like samples and raw.X '
             'do not prove raw counts; matrix names do not establish preprocessing. '
             'Only the first 128 rows are sampled; complete selected inputs are validated at loading.')
    if isinstance(source, ad.AnnData):
        matrices = [('X', source.X)] if source.X is not None else []
        matrices += [('layers/'+k, v) for k,v in source.layers.items()]
        if source.raw is not None: matrices.append(('raw.X', source.raw.X))
        return _capabilities({'format':'AnnData', 'matrices':[_matrix_report(k,v,cells=source.obs_names,genes=source.raw.var_names if k=='raw.X' else source.var_names) for k,v in matrices],
                'obs_columns': list(source.obs.columns), 'embeddings': list(source.obsm), 'notes': notes})
    source = Path(source)
    if not source.is_file(): raise FileNotFoundError('Input file not found: '+str(source))
    if source.suffix.lower() == '.rds':
        from .rds_io import inspect_rds
        report = inspect_rds(source); report['notes'] = notes+' RDS parsing is not streaming.'
        return _capabilities(report)
    if source.suffix.lower() != '.h5ad':
        raise ValueError('inspect_data supports .h5ad, AnnData or supported .rds objects')
    with h5py.File(source, 'r') as f:
        names = [('X','X')] if 'X' in f else []
        if 'layers' in f: names += [('layers/'+k,'layers/'+k) for k in f['layers']]
        if 'raw/X' in f: names.append(('raw.X','raw/X'))
        matrices = []
        obs = ad.io.read_elem(f['obs'])
        for name, key in names:
            node = f[key]
            matrix = ad.io.sparse_dataset(node) if isinstance(node,h5py.Group) else node
            var = ad.io.read_elem(f['raw/var' if name=='raw.X' else 'var'])
            matrices.append(_matrix_report(name, matrix, cells=obs.index, genes=var.index))
        columns = list(obs.columns)
        return _capabilities({'format':'h5ad', 'matrices':matrices, 'obs_columns':columns,
                'embeddings': list(f['obsm']) if 'obsm' in f else [], 'notes':notes})
