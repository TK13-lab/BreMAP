"""Python-only import of selected in-memory R single-cell matrices.

This is an input adapter, not general Seurat/SCE object conversion. No R
functions are evaluated. External-backed assays need an in-memory export.
"""
import warnings
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def _classes(obj):
    return list(map(str, np.atleast_1d(getattr(obj, 'class', []))))


def _load(path):
    try:
        import rdata
    except ImportError as exc:
        raise ImportError('RDS input requires optional dependencies: pip install "bremap[rds]". Rscript is not required.') from exc
    # These containers keep their slots; only selected numeric matrices and
    # tabular metadata are interpreted below, never commands or R functions.
    classes = ['Seurat','Assay','Assay5','StdAssay','DimReduc','Graph','Neighbor',
        'SeuratCommand','JackStrawData','SingleCellExperiment','RangedSummarizedExperiment',
        'SummarizedExperiment','SimpleAssays','SimpleList','DFrame','GRanges',
        'CompressedGRangesList','PartitioningByEnd','IRanges','Rle','Seqinfo',
        'CompressedList','dgCMatrix','dgRMatrix','dgTMatrix','lgCMatrix',
        'dgeMatrix','package_version','numeric_version','LogMap']
    constructors = dict(rdata.conversion.DEFAULT_CLASS_MAP)
    constructors.update({name: lambda obj, attrs: obj for name in classes})
    try:
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            obj = rdata.read_rds(path, constructor_dict=constructors)
        # Missing constructors for unused auxiliary slots are explicitly disclosed.
        messages = sorted({str(w.message) for w in captured})
        if messages:
            warnings.warn('RDS auxiliary conversion notes: '+'; '.join(messages[:4]), UserWarning, stacklevel=3)
        return obj
    except (NotImplementedError, TypeError, KeyError, RecursionError) as exc:
        raise ValueError('Unsupported serialized RDS structure; export the selected in-memory assay to H5AD or a named sparse matrix. '+str(exc)) from exc


def _names(value):
    if value is None or (isinstance(value,str) and value=='\x01NULL\x01'):
        return None
    if pd.isna(np.asarray(value)).any():
        raise ValueError('RDS names must be present and unique; found a missing identifier')
    return pd.Index(np.asarray(value).astype(str))


def _matrix(obj):
    """Return gene-by-cell values plus names, without densifying sparse inputs."""
    classes = _classes(obj)
    if any(c in classes for c in ['dgCMatrix','lgCMatrix','dgRMatrix','dgTMatrix']):
        shape = tuple(map(int,obj.Dim)); dtype = bool if 'lgCMatrix' in classes else float
        if 'dgTMatrix' in classes:
            x = sparse.coo_matrix((np.asarray(obj.x,dtype=dtype),(obj.i,obj.j)),shape=shape).tocsc()
        else:
            constructor, index = (sparse.csr_matrix,obj.j) if 'dgRMatrix' in classes else (sparse.csc_matrix,obj.i)
            x = constructor((np.asarray(obj.x,dtype=dtype), index, obj.p),shape=shape)
            x.check_format(full_check=True)
        return x, _names(obj.Dimnames[0]), _names(obj.Dimnames[1])
    if 'dgeMatrix' in classes:
        return np.asarray(obj.x).reshape(tuple(obj.Dim),order='F'),_names(obj.Dimnames[0]),_names(obj.Dimnames[1])
    if hasattr(obj,'dims') and hasattr(obj,'coords'):
        if obj.ndim != 2: raise ValueError('R expression must be a two-dimensional matrix')
        names = [_names(obj.coords[d].values) if d in obj.coords else None for d in obj.dims]
        return np.asarray(obj.values), *names
    if isinstance(obj,np.ndarray) and obj.ndim==2:
        return obj, None, None
    raise ValueError('Unsupported R matrix class '+str(classes or type(obj).__name__)+
                     '; export an in-memory numeric matrix. External-backed/delayed assays are not loaded.')


def _frame(obj, fallback_names=None):
    if isinstance(obj,pd.DataFrame):
        result = obj.copy()
    elif hasattr(obj,'listData'):
        names = _names(getattr(obj,'rownames',None))
        names = names if names is not None else fallback_names
        columns = {}
        for key,value in obj.listData.items():
            if 'Rle' in _classes(value):
                value = np.repeat(value.values, value.lengths)
            if isinstance(value,SimpleNamespace):
                raise ValueError('Unsupported R metadata column '+str(key)+'; export ordinary vectors/factors')
            columns[str(key)] = value
        result = pd.DataFrame(columns,index=names)
    else:
        raise ValueError('Unsupported R cell/feature metadata; expected data.frame or DFrame')
    result.index = result.index.astype(str)
    if not result.index.is_unique:
        raise ValueError('R metadata has duplicate identifiers')
    return result


def _choose(name, choices, kind):
    choices = list(map(str,choices))
    if name is None and len(choices)==1: return choices[0]
    if name not in choices:
        raise ValueError(f'Select an exact {kind}; available: {choices}. Split layers are not merged automatically.')
    return name


def _membership(logmap, layer):
    matrix, names, columns = _matrix(logmap)
    if names is None or columns is None or layer not in columns:
        raise ValueError('Malformed Seurat v5 layer membership map')
    selected = np.asarray(matrix[:, columns.get_loc(layer)]).ravel().astype(bool)
    return names[selected]


def _seurat_layers(assay):
    if hasattr(assay,'layers'):
        return dict(assay.layers)
    return {k:getattr(assay,k) for k in ['counts','data','scale.data']
            if hasattr(assay,k) and all(n>0 for n in _matrix(getattr(assay,k))[0].shape)}


def _sce_layers(obj):
    try:
        return dict(obj.assays.data.listData)
    except AttributeError as exc:
        raise ValueError('Unsupported SCE assay storage; export in-memory assays first') from exc


def _sce_genes(obj):
    names = _names(getattr(obj,'NAMES',None))
    if names is None:
        ranges = getattr(obj,'rowRanges',None)
        names = _names(getattr(getattr(ranges,'partitioning',None),'NAMES',None))
        if names is None: names = _names(getattr(ranges,'NAMES',None))
    return names


def _align_names(names, size, label):
    if names is None or len(names)!=size or not names.is_unique or names.isna().any() or any(not n.strip() for n in names):
        raise ValueError(f'RDS {label} names must be present, unique and match the matrix')
    return names


def read_rds(path, *, assay=None, matrix=None):
    obj = _load(path); classes = _classes(obj); reductions = {}; feature_meta = None
    if 'Seurat' in classes:
        assay = _choose(assay,obj.assays,'assay')
        a = obj.assays[assay]; layers = _seurat_layers(a)
        matrix = _choose(matrix,layers,'matrix/layer')
        x,genes,cells = _matrix(layers[matrix])
        if hasattr(a,'layers'):
            mapped_genes = _membership(a.features,matrix); mapped_cells = _membership(a.cells,matrix)
            if genes is not None and not genes.equals(mapped_genes): raise ValueError('Layer gene names disagree with Seurat membership')
            if cells is not None and not cells.equals(mapped_cells): raise ValueError('Layer cell names disagree with Seurat membership')
            genes,cells = mapped_genes,mapped_cells
        obs = _frame(getattr(obj,'meta.data'))
        feature_meta = getattr(a,'meta.features',None)
        if feature_meta is None and hasattr(a,'layers'):
            feature_meta = getattr(a,'meta.data',None)
        reductions = {str(k):getattr(v,'cell.embeddings') for k,v in dict(obj.reductions).items()}
    elif 'SingleCellExperiment' in classes or 'SummarizedExperiment' in classes or 'RangedSummarizedExperiment' in classes:
        layers = _sce_layers(obj)
        if assay is not None and matrix is not None: raise ValueError('SCE assay names are selected with matrix; do not also supply assay')
        matrix = _choose(matrix if matrix is not None else assay,layers,'matrix/assay'); assay = None
        x,genes,cells = _matrix(layers[matrix]); obs = _frame(obj.colData)
        if cells is None: cells = obs.index
        if genes is None: genes = _sce_genes(obj)
        feature_meta = getattr(obj,'elementMetadata',None)
        internal = getattr(getattr(obj,'int_colData',None),'listData',{})
        reductions = dict(getattr(internal.get('reducedDims'),'listData',{}))
    else:
        if assay is not None or matrix not in (None,'X'):
            raise ValueError('Standalone R matrix has no assay/layer; omit assay and matrix')
        x,genes,cells = _matrix(obj); matrix='X'
        obs = pd.DataFrame(index=cells)
    genes = _align_names(genes,x.shape[0],'gene'); cells = _align_names(cells,x.shape[1],'cell')
    if not cells.isin(obs.index).all(): raise ValueError('RDS cell metadata does not cover selected matrix cells')
    var = pd.DataFrame(index=genes)
    if feature_meta is not None:
        fm = _frame(feature_meta,genes)
        if len(fm.columns):
            if not genes.isin(fm.index).all(): raise ValueError('RDS feature metadata does not cover selected genes')
            var = fm.loc[genes].copy()
    data = ad.AnnData(X=x.T.tocsr() if sparse.issparse(x) else x.T.copy(),obs=obs.loc[cells].copy(),var=var)
    for name, reduction in reductions.items():
        v, ids, _ = _matrix(reduction)
        if ids is None:
            raise ValueError('Embedding '+str(name)+' lacks cell IDs; export with row names')
        ids = _align_names(ids,v.shape[0],'embedding cell')
        if not cells.isin(ids).all():
            warnings.warn('Skipped embedding '+str(name)+': it does not cover every selected cell',UserWarning,stacklevel=2)
            continue
        data.obsm[str(name)] = v[ids.get_indexer(cells)].copy()
    data.uns['_bremap_rds_source'] = {'class': ','.join(classes) or 'matrix', 'matrix':matrix,
                                      'assay':assay or 'none'}
    return data


def inspect_rds(path):
    from .data_io import _matrix_report
    obj = _load(path); classes = _classes(obj); matrices=[]
    if 'Seurat' in classes:
        selected = [(str(a),k,v) for a,assay in obj.assays.items() for k,v in _seurat_layers(assay).items()]
        columns = list(_frame(getattr(obj,'meta.data')).columns); embeddings=list(obj.reductions)
    elif any(c in classes for c in ['SingleCellExperiment','SummarizedExperiment','RangedSummarizedExperiment']):
        selected = [(None,k,v) for k,v in _sce_layers(obj).items()]
        columns = list(_frame(obj.colData).columns)
        internal = getattr(getattr(obj,'int_colData',None),'listData',{})
        embeddings = list(getattr(internal.get('reducedDims'),'listData',{}))
    else:
        selected=[(None,'X',obj)]; columns=[]; embeddings=[]
    for assay,name,value in selected:
        x,genes,cells = _matrix(value)
        if 'Seurat' in classes and hasattr(obj.assays[assay],'layers'):
            genes = _membership(obj.assays[assay].features,name)
            cells = _membership(obj.assays[assay].cells,name)
        row = _matrix_report(str(name),x.T,cells=cells,genes=genes)
        row['assay']=assay; matrices.append(row)
    return {'format':'rds','class':classes,'matrices':matrices,'obs_columns':list(map(str,columns)),
            'embeddings':list(map(str,embeddings))}
