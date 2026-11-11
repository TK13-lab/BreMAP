"""Portable numeric reference bundles and explicitly selected representation exports."""
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import zipfile

import anndata as ad
import numpy as np
import pandas as pd

from . import __version__
from .api import Artifact, _path, load_reference
from .legacy import _embedding
from .preprocessing import file_hash

_FILES = {'artifact.json', 'model.json', 'arrays.npz', 'spec.json', 'policy.json',
          'checksums.json', 'calibration.json'}


def export_reference(model, destination):
    """Archive a validated v2 artifact. No raw counts or pretrained weights included."""
    source, destination = _path(model), Path(destination)
    meta, _, _, calibration = load_reference(source)
    if not meta.get('spec'):
        raise ValueError('Portable export requires a versioned specification artifact')
    files = sorted(p for p in source.iterdir())
    if any(not p.is_file() or p.is_symlink() or p.name not in _FILES for p in files):
        raise ValueError('Reference contains unexpected files or links')
    if destination.exists():
        raise ValueError('Export destination already exists')
    destination.parent.mkdir(parents=True, exist_ok=True)
    manifest = {'schema': 'bremap-reference-bundle-v1', 'software_version': __version__,
                'tissue': meta['spec']['tissue'], 'spec_id': meta['spec']['spec_id'],
                'representation': meta['representation'], 'calibrated': calibration is not None,
                'files': {p.name: file_hash(p) for p in files},
                'contents': 'Fitted numeric arrays, configuration and calibration; no counts or FM weights',
                'limits': 'Contains reference/calibration cell and donor identifiers; signatures are user-specified, not biologically certified'}
    with tempfile.TemporaryDirectory(dir=destination.parent) as temp:
        archive = Path(temp) / 'bundle.zip'
        with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED) as z:
            z.writestr('manifest.json', json.dumps(manifest, indent=2))
            for path in files:
                z.write(path, 'reference/' + path.name)
        # Exclusive creation prevents accidental replacement of an existing export.
        with archive.open('rb') as src, destination.open('xb') as dst:
            shutil.copyfileobj(src, dst)
    return destination


def import_reference(archive, destination):
    """Validate allowed archive members/checksums before publishing an artifact."""
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Import destination already exists')
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as z:
        names = z.namelist()
        allowed = {'manifest.json'} | {'reference/' + name for name in _FILES}
        if len(set(names)) != len(names) or set(names) - allowed or 'manifest.json' not in names:
            raise ValueError('Unexpected or duplicate reference archive member')
        if any(i.file_size > 1024**3 for i in z.infolist()) or sum(i.file_size for i in z.infolist()) > 2*1024**3:
            raise ValueError('Reference archive exceeds supported unpacked size')
        manifest = json.loads(z.read('manifest.json'))
        files = manifest.get('files', {})
        if (manifest.get('schema') != 'bremap-reference-bundle-v1' or not isinstance(files, dict)
                or not files or set(files) - _FILES
                or set(names) != {'manifest.json'} | {'reference/' + n for n in files}):
            raise ValueError('Invalid reference bundle manifest')
        with tempfile.TemporaryDirectory(dir=destination.parent) as temp:
            staged = Path(temp) / 'reference'; staged.mkdir()
            for name, digest in files.items():
                payload = z.read('reference/' + name)
                if hashlib.sha256(payload).hexdigest() != digest:
                    raise ValueError('Reference bundle checksum mismatch: ' + name)
                (staged / name).write_bytes(payload)
            meta, _, _, cal = load_reference(staged)
            if (manifest['tissue'] != meta['spec']['tissue'] or manifest['spec_id'] != meta['spec']['spec_id']
                    or manifest['representation'] != meta['representation']
                    or manifest['calibrated'] != (cal is not None)):
                raise ValueError('Bundle summary disagrees with its reference artifact')
            if destination.exists():
                raise ValueError('Import destination already exists')
            staged.rename(destination)
    return Artifact(destination, json.loads((destination / 'artifact.json').read_text()))


def export_embeddings(source, destination, *, key, obs_columns=()):
    """Write a count-free H5AD with vectors in X; not a mapping-count input.

    Only cell IDs and explicitly requested observation columns are retained.
    Source gene universe and embedding provenance are retained in uns.
    """
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Embedding export destination already exists')
    a = ad.read_h5ad(source, backed='r')
    try:
        if not a.obs_names.is_unique:
            raise ValueError('Embedding export requires unique cell IDs')
        columns = list(obs_columns)
        if len(set(columns)) != len(columns) or set(columns) - set(a.obs.columns):
            raise ValueError('Metadata columns must be unique and present in the source')
        values, provenance = _embedding(a, key)
        original = np.asarray(a.obsm[key])
        if np.issubdtype(original.dtype, np.number):
            values = original
        result = ad.AnnData(X=values, obs=a.obs.loc[:, columns].copy(),
                           var=pd.DataFrame(index=[f'embedding_{i}' for i in range(values.shape[1])]))
        result.uns['bremap_representation_export'] = {
            'schema': 'bremap-representation-v1', 'software_version': __version__,
            'embedding_key': key, 'encoder_provenance': provenance,
            'source_sha256': file_hash(source), 'source_gene_universe': list(a.var_names),
            'counts_included': False, 'dimensions_are': 'latent coordinates, not named pathways'}
        destination.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=destination.parent) as temp:
            path = Path(temp) / 'vectors.h5ad'; result.write_h5ad(path, compression='gzip')
            with path.open('rb') as src, destination.open('xb') as dst:
                shutil.copyfileobj(src, dst)
    finally:
        a.file.close()
    return destination
