"""Path-based reference service with immutable, checksummed artifact bundles."""
from dataclasses import dataclass
import json
from pathlib import Path
import shutil
import tempfile
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd

from . import legacy
from .calibration import load_policy, select_settings
from .core import _identifiers
from .preprocessing import file_hash
from .specification import ReferenceSpec


@dataclass(frozen=True)
class Artifact:
    path: Path
    metadata: dict


ReferenceModel = CalibratedReference = MappingResult = Artifact


@dataclass(frozen=True)
class AnalysisResult(Artifact):
    """Saved analysis bundle and the same tables written beneath ``path``."""
    tables: dict


ReliabilityResult = EvaluationResult = AnalysisResult


def _path(value):
    return Path(value.path if isinstance(value, Artifact) else value)


def _json(path):
    return json.loads(Path(path).read_text())


def _seal(path):
    legacy._write_json(path / 'checksums.json', {p.name: file_hash(p) for p in sorted(path.iterdir())
                                               if p.is_file() and p.name != 'checksums.json'})


def _verify(path):
    checks = _json(path / 'checksums.json')
    required = {'artifact.json', 'model.json', 'arrays.npz', 'spec.json', 'policy.json'}
    files = {p.name for p in path.iterdir() if p.is_file() and p.name != 'checksums.json'}
    if not isinstance(checks, dict) or not required <= checks.keys() or set(checks) != files:
        raise ValueError('Artifact checksum manifest is incomplete')
    for name, digest in checks.items():
        if Path(name).name != name or file_hash(path / name) != digest:
            raise ValueError(f'Artifact checksum mismatch: {name}')


def _number(value, name, low=0, high=None, positive=False):
    if isinstance(value, (str, bool)) or not isinstance(value, (int, float)) or not np.isfinite(value):
        raise ValueError(f'Invalid numeric metadata: {name}')
    if value < low or (positive and value <= low) or (high is not None and value > high):
        raise ValueError(f'Invalid numeric metadata: {name}')


def _validate_settings(settings, method):
    _number(settings['beta'], 'beta', high=1)
    _number(settings['temperature'], 'temperature', positive=True)
    _number(settings['threshold'], 'threshold', high=1 if method == 'transport' else None)
    if method == 'transport':
        _number(settings['reject_cost'], 'reject_cost')
        _number(settings['entropy'], 'entropy', positive=True)
        _number(settings['marginal_penalty'], 'marginal_penalty')
        _number(settings['reject_prior'], 'reject_prior', high=1, positive=True)
        if settings['reject_prior'] == 1:
            raise ValueError('Invalid reject_prior')
    elif settings['reject_cost'] is not None:
        raise ValueError('Affinity calibration cannot have a reject cost')
    expected = 'reject_mass' if method == 'transport' else 'reference_distance'
    if settings['threshold_variable'] != expected:
        raise ValueError('Calibration threshold variable does not match method')


def _validate_policy_settings(settings, cost_mode, policy):
    """Check declared fixed/grid constraints; this does not authenticate tuning."""
    legacy_policy = policy['policy_id'] == 'legacy_cli_v1'
    allowed = {
        'beta': [0.] if cost_mode == 'embedding' else ([policy['beta']] if legacy_policy else policy['transfer_beta_grid']),
        'temperature': [policy['temperature']] if legacy_policy else policy['transfer_temperature_grid'],
        'entropy': [policy['entropy']] if legacy_policy else policy['entropy_grid'],
        'marginal_penalty': [policy['marginal_penalty']] if legacy_policy else policy['marginal_penalty_grid'],
        'reject_prior': [policy['reject_prior']],
    }
    for name, values in allowed.items():
        if name in settings and settings[name] not in values:
            raise ValueError(f'Calibration policy/cost mode does not permit {name}={settings[name]}')


def load_reference(model):
    """Return (metadata, numeric arrays, mapper, calibration-or-None).

    Metadata is an independent dict, with effective selected settings applied to
    the returned mapper and transport options. Legacy v1 directories are readable;
    they have no verified tissue identity. No artifact is modified by this reader.
    """
    path = _path(model)
    if not (path / 'artifact.json').exists():
        meta, arrays, mapper = legacy._load_model(path)
        meta['tissue_identity_verified'] = False
        calibration = _json(path / 'calibration.json') if (path / 'calibration.json').exists() else None
        if calibration is not None:
            if calibration.get('schema') != 'bremap-calibration-v1' or calibration.get('model_sha256') != file_hash(path / 'model.json'):
                raise ValueError('Calibration schema or model checksum mismatch')
            _number(calibration['threshold'], 'threshold')
        return meta, arrays, mapper, calibration
    _verify(path)
    artifact = _json(path / 'artifact.json')
    if artifact.get('schema') != 'bremap-artifact-v2' or artifact.get('kind') not in ('reference', 'calibrated'):
        raise ValueError('Unsupported artifact schema or kind')
    spec = ReferenceSpec.load(path / 'spec.json')
    policy = _json(path / 'policy.json')
    if policy != load_policy(artifact.get('policy_id')):
        raise ValueError('Artifact policy does not match the versioned policy')
    if artifact.get('spec_sha256') != file_hash(path / 'spec.json') or artifact.get('policy_sha256') != file_hash(path / 'policy.json'):
        raise ValueError('Artifact spec/policy checksum mismatch')
    meta, arrays, mapper = legacy._load_model(path)
    if (artifact.get('model_sha256') != file_hash(path / 'model.json') or
            artifact.get('arrays_sha256') != meta['arrays_sha256'] or
            artifact.get('reference_sha256') != meta['reference_sha256']):
        raise ValueError('Artifact model/reference checksum linkage mismatch')
    if (meta['classes'] != spec.data['states'] or meta['programs'] != spec.data['programs'] or
            meta['donor_column'] != spec.data['donor_column'] or meta['label_column'] != spec.data['label_column'] or
            meta['representation']['kind'] != spec.data['representation']['kind']):
        raise ValueError('Model metadata does not match the specification')
    _number(meta['beta'], 'beta', high=1)
    _number(meta['temperature'], 'temperature', positive=True)
    _validate_policy_settings({'beta': meta['beta'], 'temperature': meta['temperature'], **meta['transport']},
                              spec.data['cost_mode'], policy)
    if len(set(meta['genes'])) != len(meta['genes']) or meta['genes'] != [g.upper() for g in meta['genes']]:
        raise ValueError('Malformed model gene order')
    for key in ('reference_cells', 'reference_donors'):
        values = artifact.get(key)
        if not isinstance(values, list) or not values or len(set(values)) != len(values):
            raise ValueError(f'Malformed artifact {key}')
    if artifact['reference_donors'] != meta['reference_donors']:
        raise ValueError('Reference donor identities do not match model')
    meta.update(spec=spec.data, cost_mode=spec.data['cost_mode'], policy=policy,
                tissue_identity_verified=True, reference_cells=artifact['reference_cells'], artifact=artifact)
    calibration = None
    if artifact['kind'] == 'calibrated':
        calibration = _json(path / 'calibration.json')
        if (calibration.get('schema') != 'bremap-calibration-v2' or
                calibration.get('model_sha256') != file_hash(path / 'model.json') or
                calibration.get('spec_sha256') != artifact['spec_sha256'] or
                calibration.get('policy_sha256') != artifact['policy_sha256'] or
                calibration.get('source_artifact_sha256') != artifact.get('source_artifact_sha256')):
            raise ValueError('Calibration linkage does not match model/spec/policy')
        selected = calibration['settings']
        _validate_settings(selected, meta['method'])
        _validate_policy_settings(selected, meta['cost_mode'], policy)
        if calibration['threshold'] != selected['threshold'] or calibration['reject_cost'] != selected['reject_cost']:
            raise ValueError('Inconsistent calibration settings')
        legacy._disjoint(calibration['donors'], meta['reference_donors'])
        if set(calibration['cells']) & set(meta['reference_cells']):
            raise ValueError('Cell overlap between reference and calibration')
        mapper.beta, mapper.temperature = selected['beta'], selected['temperature']
        meta.update(beta=mapper.beta, temperature=mapper.temperature)
        if meta['method'] == 'transport':
            meta['transport'] = {k: selected[k] for k in ('entropy', 'marginal_penalty', 'reject_prior')}
    elif (path / 'calibration.json').exists():
        raise ValueError('Reference artifact cannot contain a calibration file')
    return meta, arrays, mapper, calibration


def transform_input(path, metadata, arrays, require_labels=False):
    """Return frozen (X, programs, cell_ids, donors, labels-or-None, zero_mask).

    Never fit query/reference representations. Labels are read only on request.
    Zero rows are returned for callers to reject (reference/calibration) or isolate
    before prediction (queries). Study IDs are deliberately optional here.
    """
    if metadata.get('spec'):
        data = ad.read_h5ad(path)
        ReferenceSpec.from_dict(metadata['spec']).validate_identity(data)
    return legacy._transform(path, metadata, arrays, require_labels)


def fit_reference(reference, spec, model_dir, method='transport', policy=None):
    if policy is None:
        raise ValueError('Spec-based workflows require an explicit calibration policy')
    if method not in ('affinity', 'transport'):
        raise ValueError('Method must be affinity or transport')
    out = _path(model_dir)
    if out.exists():
        raise ValueError('Output already exists; choose a new immutable artifact directory')
    policy_data = load_policy(policy)
    spec = spec if isinstance(spec, ReferenceSpec) else ReferenceSpec.load(spec)
    data = ad.read_h5ad(reference)
    spec.validate_identity(data)
    label = spec.data['label_column']
    if label not in data.obs or set(data.obs[label]) != set(spec.data['states']):
        raise ValueError('Reference must contain exactly the declared states')
    with tempfile.TemporaryDirectory() as tmp:
        programs = Path(tmp) / 'programs.json'
        legacy._write_json(programs, spec.data['programs'])
        legacy.fit_reference(SimpleNamespace(reference=reference, model=out, programs=programs,
            method=method, donor_column=spec.data['donor_column'], label_column=label,
            embedding_key=spec.data['representation'].get('key')))
    metadata = _json(out / 'model.json')
    if spec.data['cost_mode'] == 'embedding':
        metadata['beta'] = 0.
        legacy._write_json(out / 'model.json', metadata)
    legacy._write_json(out / 'spec.json', spec.data)
    legacy._write_json(out / 'policy.json', policy_data)
    artifact = {'schema': 'bremap-artifact-v2', 'kind': 'reference', 'policy_id': policy,
                'spec_sha256': file_hash(out / 'spec.json'), 'policy_sha256': file_hash(out / 'policy.json'),
                'reference_cells': _identifiers(data.obs_names, 'cell IDs', data.n_obs).tolist(),
                'reference_donors': metadata['reference_donors'],
                'model_sha256': file_hash(out / 'model.json'), 'arrays_sha256': file_hash(out / 'arrays.npz'),
                'reference_sha256': file_hash(reference)}
    legacy._write_json(out / 'artifact.json', artifact)
    _seal(out)
    return Artifact(out, artifact)


def calibrate(model, calibration, calibrated_dir):
    source, out = _path(model), _path(calibrated_dir)
    meta, arrays, mapper, existing = load_reference(source)
    if existing is not None:
        raise ValueError('Input is already calibrated; calibrate the base reference instead')
    if out.exists():
        raise ValueError('Output already exists; choose a new immutable artifact directory')
    if not meta.get('spec'):
        shutil.copytree(source, out)
        legacy.calibrate(SimpleNamespace(model=out, calibration=calibration))
        return Artifact(out, _json(out / 'calibration.json'))
    X, programs, ids, donors, labels, zero = transform_input(calibration, meta, arrays, True)
    legacy._disjoint(donors, meta['reference_donors'])
    if set(ids) & set(meta['reference_cells']):
        raise ValueError('Cell overlap between reference and calibration')
    if zero.any():
        raise ValueError('Calibration contains zero-library cells')
    if set(labels) != set(meta['classes']):
        raise ValueError('Calibration must contain exactly the declared known states')
    settings, prediction = select_settings(mapper, X, programs, labels, donors, meta['method'], meta['cost_mode'], meta['policy'])
    _validate_settings(settings, meta['method'])
    _validate_policy_settings(settings, meta['cost_mode'], meta['policy'])
    shutil.copytree(source, out)
    artifact = dict(meta['artifact'], kind='calibrated', source_artifact_sha256=file_hash(source / 'artifact.json'))
    record = {'schema': 'bremap-calibration-v2', 'model_sha256': file_hash(out / 'model.json'),
        'spec_sha256': artifact['spec_sha256'], 'policy_sha256': artifact['policy_sha256'],
        'source_artifact_sha256': artifact['source_artifact_sha256'], 'input_sha256': file_hash(calibration),
        'cells': ids.tolist(), 'donors': sorted(set(donors)), 'settings': settings,
        'threshold': settings['threshold'], 'reject_cost': settings['reject_cost'],
        'known_acceptance_target': .95, 'observed_cell_acceptance': float(np.mean(prediction['divergence'] <= settings['threshold'])),
        'n_cells': len(X), 'versions': legacy._versions(),
        'limits': 'Normal-only selection and thresholds reuse calibration donors; no formal coverage guarantee.'}
    legacy._write_json(out / 'calibration.json', record)
    legacy._write_json(out / 'artifact.json', artifact)
    _seal(out)
    return Artifact(out, artifact)


def map_query(calibrated_model, query, out_dir):
    path, out = _path(calibrated_model), _path(out_dir)
    meta, arrays, mapper, calibration = load_reference(path)
    if calibration is None:
        raise ValueError('Calibrate the model before mapping a query')
    if meta.get('spec'):
        if out.exists():
            raise ValueError('Output already exists; choose a new mapping directory')
        _, _, ids, donors, _, _ = transform_input(query, meta, arrays)
        if set(ids) & (set(meta['reference_cells']) | set(calibration['cells'])):
            raise ValueError('Cell overlap between query and reference/calibration')
    args = SimpleNamespace(model=path, query=query, out=out)
    legacy._map_with_model(args, meta, arrays, mapper, calibration)
    record = _json(out / 'reproducibility.json')
    if meta.get('spec'):
        frame = pd.read_csv(out / 'cells.csv', dtype=str, keep_default_na=False)
        frame = frame.rename(columns={'candidate_analog': 'candidate_state', 'best_analog': 'accepted_analog',
            'analog_strength': 'conditional_affinity_max', 'geometric_divergence': 'reference_distance',
            'abstain_reason': 'abstention_reason', **{f'affinity::{c}': f'conditional_affinity::{c}' for c in meta['classes']}})
        frame['accepted'] = frame.pop('abstain_flag').eq('False')
        frame = frame.drop(columns='divergence')
        frame.to_csv(out / 'cells.csv', index=False)
        donors = pd.read_csv(out / 'donors.csv', dtype=str, keep_default_na=False)
        donors = donors.rename(columns={'coverage': 'acceptance', 'mean_divergence': 'mean_reject_mass' if meta['method'] == 'transport' else 'mean_reference_distance',
            **{f'mean_affinity::{c}': f'mean_conditional_affinity::{c}' for c in meta['classes']}})
        donors.to_csv(out / 'donors.csv', index=False)
        record.update(schema='bremap-mapping-v2', spec=meta['spec'], policy=meta['policy'], cost_mode=meta['cost_mode'],
            settings=calibration['settings'], tissue_identity_verified=True,
            artifact_sha256=file_hash(path / 'artifact.json'),
            transport_group_dependence='Complete donor groups; splitting the same donor changes transport predictions.',
            limits=calibration['limits'],
            output_sha256={name: file_hash(out / name) for name in ('cells.csv', 'donors.csv')})
        for name in ('api.py', 'calibration.py', 'specification.py'):
            record['code_sha256'][name] = file_hash(Path(__file__).with_name(name))
    else:
        record['tissue_identity_verified'] = False
    legacy._write_json(out / 'reproducibility.json', record)
    return Artifact(out, record)


def assess_reliability(calibrated_model, reference, query, out_dir,
                       n_bootstrap=20, seed=42, leave_study_out=False):
    """Assess reference sensitivity with frozen representation and calibration."""
    from .reliability import assess_reliability as assess
    return assess(calibrated_model, reference, query, out_dir, n_bootstrap, seed, leave_study_out)


def evaluate(mapping, annotations, out_dir):
    """Evaluate saved predictions using exact cell_id/donor_id annotation joins."""
    from .evaluation import evaluate as evaluate_frozen
    return evaluate_frozen(mapping, annotations, out_dir)
