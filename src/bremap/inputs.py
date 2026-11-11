"""Prepare declared package inputs from an AnnData object or existing H5AD."""
from pathlib import Path
import tempfile

import anndata as ad
import pandas as pd

from .configuration import inspect_input
from .specification import ReferenceSpec


def prepare_input(source, spec, destination, *, role='query', cell_ids=None,
                  donor_column=None, label_column=None):
    """Copy selected data, declare spec identity, validate, and save an H5AD.

    The caller declares tissue/state eligibility; this function does not infer it.
    X must contain raw integer counts. Supplied embeddings/provenance stay in the
    AnnData. Source objects/files are never modified and existing outputs are refused.
    """
    spec = spec if isinstance(spec, ReferenceSpec) else ReferenceSpec.load(spec)
    destination = Path(destination)
    if destination.exists() or destination.is_symlink():
        raise ValueError('Input-preparation output already exists')
    data = source.copy() if isinstance(source, ad.AnnData) else ad.read_h5ad(source)
    if not data.obs_names.is_unique:
        raise ValueError('Source cell IDs must be unique')
    if cell_ids is not None:
        ids = pd.Index(cell_ids)
        if not len(ids) or not ids.is_unique or not ids.isin(data.obs_names).all():
            raise ValueError('Selected cell IDs must be nonempty, unique and present')
        data = data[ids].copy()
    for supplied, target in [(donor_column, spec.data['donor_column']), (label_column, spec.data['label_column'])]:
        if supplied is not None:
            if supplied not in data.obs:
                raise ValueError('Missing source metadata column: '+supplied)
            if target in data.obs and not data.obs[target].equals(data.obs[supplied]):
                raise ValueError('Conflicting metadata column: '+target)
            data.obs[target] = data.obs[supplied].copy()
    data.uns['bremap'] = {key: spec.data[key] for key in ['tissue', 'spec_id', 'spec_version']}
    # Use broadly supported H5AD strings for metadata imported from Parquet.
    for frame in [data.obs, data.var]:
        frame.index = pd.Index(frame.index.to_numpy(dtype=str), name=frame.index.name, dtype=object)
        for name in frame.columns:
            if isinstance(frame[name].dtype, pd.StringDtype):
                frame[name] = frame[name].astype(object)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as temp:
        staged = Path(temp)/'prepared.h5ad'
        data.write_h5ad(staged, compression='gzip')
        check = inspect_input(staged, spec, role)
        if not check['valid']:
            raise ValueError('; '.join(check['errors']))
        # Exclusive destination creation, including if another caller raced us.
        with staged.open('rb') as reader, destination.open('xb') as writer:
            import shutil
            shutil.copyfileobj(reader, writer)
    return destination
