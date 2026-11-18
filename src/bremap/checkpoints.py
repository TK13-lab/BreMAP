"""Official checkpoint selection and consent-based downloads, independent of torch."""
import hashlib
from importlib.resources import files
import json
import os
from pathlib import Path
import zipfile

import pandas as pd

from .preprocessing import file_hash


def _catalog():
    return json.loads(files('bremap').joinpath('resources/checkpoints.json').read_text())


def list_checkpoints(model=None):
    """List supported official checkpoints without network access.

    ``model`` is None, ``"scgpt"`` or ``"geneformer"``. Returns a DataFrame with
    checkpoint IDs, family/version, description, source, revision and download
    size in bytes (missing when the source does not publish a fixed size).
    Listing a checkpoint does not imply it is validated for your biological task.
    """
    if model not in (None, 'scgpt', 'geneformer'):
        raise ValueError('model must be scgpt or geneformer, or None for all')
    return pd.DataFrame([
        {k: e[k] for k in ('checkpoint_id', 'model', 'model_version', 'description',
                           'source', 'revision', 'size_bytes')}
        for e in _catalog() if model is None or e['model'] == model
    ])


def _ask(message):
    try:
        return input(message).strip()
    except (EOFError, OSError) as exc:
        raise RuntimeError('Interactive input unavailable. Choose a checkpoint_id and pass '
                           'confirm=True only after approving its download.') from exc


def _select(checkpoint_id):
    entries = _catalog()
    if checkpoint_id is None:
        for i, entry in enumerate(entries, 1):
            print(f"{i}. {entry['checkpoint_id']}: {entry['description']}")
        choice = _ask('Select checkpoint number or ID (blank cancels): ')
        if not choice:
            raise RuntimeError('Checkpoint download cancelled.')
        if choice.isdigit() and 1 <= int(choice) <= len(entries):
            choice = entries[int(choice)-1]['checkpoint_id']
        checkpoint_id = choice
    for entry in entries:
        if entry['checkpoint_id'] == checkpoint_id:
            return entry
    raise ValueError(f'Unknown checkpoint {checkpoint_id!r}; use list_checkpoints() to see choices.')


def _validate_file(entry, name, path):
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f'Missing or empty checkpoint file: {path}')
    with path.open('rb') as handle:
        prefix = handle.read(256).lstrip().lower()
    if prefix.startswith((b'<!doctype html', b'<html', b'version https://git-lfs')):
        raise ValueError(f'{name} contains HTML or a Git LFS pointer, not model weights/configuration')
    if name == 'best_model.pt':
        try:
            with zipfile.ZipFile(path) as archive:
                if not any(n.endswith('/data.pkl') for n in archive.namelist()) or archive.testzip():
                    raise ValueError('Incomplete PyTorch weights archive')
        except zipfile.BadZipFile as exc:
            raise ValueError('Invalid or truncated PyTorch weights archive') from exc
    if name.endswith('.json'):
        value = json.loads(path.read_text())
        if not isinstance(value, dict) or not value:
            raise ValueError(f'{name} must contain a nonempty JSON object')
        if name == 'args.json':
            required = {'embsize', 'nheads', 'd_hid', 'nlayers', 'n_layers_cls',
                        'dropout', 'pad_token', 'pad_value'}
            if not required <= value.keys():
                raise ValueError(f'args.json missing scGPT configuration keys: {sorted(required-value.keys())}')
        if name == 'vocab.json' and any(type(v) is not int or v < 0 for v in value.values()):
            raise ValueError('vocab.json must map gene/token names to nonnegative integer IDs')
        if name == 'config.json':
            expected = 2048 if entry['model_version'] == 'V1' else 4096
            if value.get('max_position_embeddings') != expected:
                raise ValueError(f'config.json does not match Geneformer {entry["model_version"]}')
    digest = file_hash(path)
    expected_hash = entry.get('sha256', {}).get(name)
    if expected_hash and digest != expected_hash:
        raise ValueError(f'{name} differs from the official pinned checkpoint SHA256')
    return digest


def _download_file(entry, name, destination):
    if entry['model'] == 'scgpt':
        try:
            import gdown
        except ImportError as exc:
            raise ImportError('Install the BReMAP download extra from the source checkout: '
                              'python -m pip install ".[download]". See the installation guide.') from exc
        # Only the three declared files are transferred, even if the folder contains extras.
        listing = gdown.download_folder(url=entry['source'], skip_download=True,
                                        quiet=True, use_cookies=False)
        matches = [f for f in listing or [] if f.path == name]
        if len(matches) != 1:
            raise ValueError(f'Official scGPT folder must contain exactly one {name}')
        result = gdown.download(id=matches[0].id, output=str(destination),
                               quiet=False, resume=True, use_cookies=False)
        if result is None:
            raise OSError('Google Drive did not return a file; check access or download quota')
    else:
        try:
            from huggingface_hub import hf_hub_download
        except ImportError as exc:
            raise ImportError('Install the BReMAP download extra from the source checkout: '
                              'python -m pip install ".[download]". See the installation guide.') from exc
        # local_dir avoids keeping a second full copy in the Hub global cache.
        downloaded = hf_hub_download(repo_id='ctheodoris/Geneformer',
            filename=f"{entry['subfolder']}/{name}", revision=entry['revision'],
            local_dir=str(destination.parent), token=False)
        Path(downloaded).replace(destination)


def _read_manifest(folder, signature, names, complete):
    manifest = folder / '.bremap-checkpoint.json'
    if not manifest.is_file():
        raise ValueError(f'Unmanaged checkpoint directory (no BReMAP manifest): {folder}. '
                         'Use TeacherEncoder with this local path, or choose another cache_dir.')
    record = json.loads(manifest.read_text())
    if record.get('catalog_sha256') != signature:
        raise ValueError(f'Checkpoint catalog changed for {folder}; choose another cache_dir.')
    hashes = record.get('files', {})
    if not isinstance(hashes, dict) or not set(hashes) <= set(names) or (complete and set(hashes) != set(names)):
        raise ValueError(f'Invalid checkpoint manifest in {folder}; inspect or remove that checkpoint directory.')
    for name, digest in hashes.items():
        path = folder/name
        if not path.is_file() or file_hash(path) != digest:
            raise ValueError(f'Cached checkpoint file is corrupt or changed: {path}. '
                             f'Remove only {folder} and retry, or use another cache_dir.')
    return record


def _write_manifest(folder, record):
    temp = folder/'.bremap-checkpoint.json.tmp'
    temp.write_text(json.dumps(record, indent=2))
    temp.replace(folder/'.bremap-checkpoint.json')


def download_checkpoint(checkpoint_id=None, *, cache_dir=None, confirm=None):
    """Download a selected official checkpoint and return its local ``Path``.

    With no ID, show a numbered selection menu. Before a new download, print the
    source, required files, size if known and destination. ``confirm=None`` asks
    [y/N]; ``False`` cancels; ``True`` explicitly authorizes an unattended download
    (with an explicit ID). No network or cache writes occur before consent.

    The default cache is ``$XDG_CACHE_HOME/bremap/checkpoints`` or
    ``~/.cache/bremap/checkpoints``. A completed, hash-verified cache is reused
    offline without asking again. Interrupted transfers retain completed files
    in a hidden staging directory; rerun the same call to retry. A lock prevents
    simultaneous writers. Corrupt or unmanaged directories fail without overwrite.
    Dependencies in ``bremap[download]`` are optional. Model weights and inference
    packages are not bundled; this function does not install inference packages.
    """
    if confirm is not None and type(confirm) is not bool:
        raise ValueError('confirm must be None (ask), True (approved), or False (cancel)')
    if confirm is False:
        raise RuntimeError('Checkpoint download cancelled.')
    if checkpoint_id is None and confirm is True:
        raise ValueError('Provide checkpoint_id with confirm=True; no model is chosen automatically.')
    entry = _select(checkpoint_id)
    root = Path(cache_dir) if cache_dir is not None else Path(
        os.environ.get('XDG_CACHE_HOME', str(Path.home()/'.cache'))) / 'bremap' / 'checkpoints'
    root = root.expanduser().resolve()
    target = root / entry['checkpoint_id']
    signature = hashlib.sha256(json.dumps(entry, sort_keys=True).encode()).hexdigest()
    if target.exists():
        _read_manifest(target, signature, entry['files'], complete=True)
        print(f'Using verified checkpoint: {target}')
        return target
    size = f"{entry['size_bytes']/1_000_000:.1f} MB" if entry['size_bytes'] else 'not published; weights may be large'
    print(f"Checkpoint: {entry['checkpoint_id']}\nSource: {entry['source']}\n"
          f"Files: {', '.join(entry['files'])}\nDownload size: {size}\nSave to: {target}")
    if confirm is None and _ask('Download this checkpoint? [y/N] ').lower() not in ('y', 'yes'):
        raise RuntimeError('Checkpoint download cancelled.')
    root.mkdir(parents=True, exist_ok=True)
    lock = root/f".{entry['checkpoint_id']}.lock"
    try:
        lock.mkdir()
    except FileExistsError as exc:
        raise RuntimeError(f'Download already in progress: {lock}. If no process is running, '
                           'remove this stale lock directory and retry.') from exc
    staging = root/f".{entry['checkpoint_id']}.partial"
    try:
        if target.exists():
            _read_manifest(target, signature, entry['files'], complete=True)
            return target
        if staging.exists():
            record = _read_manifest(staging, signature, entry['files'], complete=False)
        else:
            staging.mkdir()
            record = {'checkpoint_id': entry['checkpoint_id'], 'catalog_sha256': signature,
                      'source': entry['source'], 'revision': entry['revision'], 'files': {}}
            _write_manifest(staging, record)
        for name in entry['files']:
            if name in record['files']:
                continue
            try:
                _download_file(entry, name, staging/name)
                try:
                    digest = _validate_file(entry, name, staging/name)
                except (ValueError, UnicodeError):
                    # This failed payload belongs to this download; never resume it as weights.
                    (staging/name).unlink(missing_ok=True)
                    raise
                record['files'][name] = digest
                _write_manifest(staging, record)
            except ImportError:
                raise
            except Exception as exc:
                raise RuntimeError(f'Download failed for {entry["checkpoint_id"]}/{name}: {exc}. '
                                   f'Check {entry["source"]} and retry the same call; completed files are retained.') from exc
        staging.rename(target)
    finally:
        lock.rmdir()
    print(f'Checkpoint ready: {target}')
    return target
