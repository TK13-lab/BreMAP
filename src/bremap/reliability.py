"""Conditional reference sensitivity; these diagnostics do not change acceptance.

Representation/SVD, model settings, reject cost and acceptance threshold stay
fixed. Each valid draw refits donor-balanced scalers and state prototypes, as in
the audited study workflow. SDs are population SDs over valid draws; donor tables
average defined per-cell diagnostics, with separate per-draw donor summaries.
"""
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from .api import AnalysisResult, _path, load_reference, transform_input
from .core import ReferenceMapper, _identifiers
from .legacy import _disjoint, _predict, _versions, _write_json
from .preprocessing import file_hash


def _bootstrap_rows(donors, sampled):
    """Distinct draw IDs preserve repeated donor weight in scaling/prototypes."""
    blocks = [np.flatnonzero(donors == donor) for donor in np.unique(donors)]
    selected = [blocks[i] for i in sampled]
    return (np.concatenate(selected),
            np.concatenate([np.repeat(f'draw::{i}', len(block)) for i, block in enumerate(selected)]))


def _summary(ids, donors, zero, baseline, predictions, classes, requested, method):
    valid = ~zero
    undefined = 'not_requested' if requested == 0 else ('' if predictions else 'no_valid_draws')
    frame = pd.DataFrame({'cell_id': ids, 'donor_id': donors, 'qc_zero_library': zero,
        'n_requested': requested, 'n_valid': len(predictions),
        'undefined_reason': np.where(zero, 'zero_library', undefined)})
    metrics = {'candidate_agreement': [], 'accepted_agreement': [], 'acceptance_frequency': [],
               'reference_distance_sd': [], **{f'affinity_sd::{c}': [] for c in classes}}
    if method == 'transport':
        metrics['reject_mass_sd'] = []
    if predictions and valid.any():
        metrics['candidate_agreement'] = np.mean([p['candidate_analog'] == baseline['candidate_analog'] for p in predictions], axis=0)
        metrics['accepted_agreement'] = np.mean([p['accepted'] == baseline['accepted'] for p in predictions], axis=0)
        metrics['acceptance_frequency'] = np.mean([p['accepted'] for p in predictions], axis=0)
        metrics['reference_distance_sd'] = np.std([p['geometric_divergence'] for p in predictions], axis=0)
        sd = np.std([p['affinity'] for p in predictions], axis=0)
        for i, state in enumerate(classes):
            metrics[f'affinity_sd::{state}'] = sd[:, i]
        if method == 'transport':
            metrics['reject_mass_sd'] = np.std([p['reject_mass'] for p in predictions], axis=0)
    for name, values in metrics.items():
        frame[name] = np.nan
        if len(values):
            frame.loc[valid, name] = values
    grouped = frame.groupby('donor_id', sort=True)
    donor = grouped[list(metrics)].mean()
    donor['n_cells'] = grouped.size()
    donor['n_zero_library'] = grouped.qc_zero_library.sum()
    donor['n_requested'], donor['n_valid'] = requested, len(predictions)
    donor['undefined_reason'] = np.where(donor.n_zero_library == donor.n_cells, 'zero_library',
                                         undefined)
    return frame, donor.reset_index()


def assess_reliability(calibrated_model, reference, query, out_dir,
                       n_bootstrap=20, seed=42, leave_study_out=False):
    """Return an AnalysisResult with bootstrap/study tables and a draw ledger.

    The original reference file must match its recorded SHA256. Study IDs come
    from the spec's study_column only when leave_study_out=True. Query labels
    are never read. Zero-library rows retain undefined sensitivity diagnostics.
    """
    if isinstance(n_bootstrap, (bool, np.bool_)) or not isinstance(n_bootstrap, (int, np.integer)) or n_bootstrap < 1:
        raise ValueError('n_bootstrap must be a positive integer')
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError('seed must be a nonnegative integer')
    if not isinstance(leave_study_out, bool):
        raise ValueError('leave_study_out must be boolean')
    out = _path(out_dir)
    if out.exists():
        raise ValueError('Output already exists; choose a new sensitivity directory')
    meta, arrays, mapper, calibration = load_reference(calibrated_model)
    if calibration is None:
        raise ValueError('Calibrate the model before reference sensitivity')
    if file_hash(reference) != meta['reference_sha256']:
        raise ValueError('Original reference checksum mismatch')
    xt, pt, reference_ids, dt, yt, reference_zero = transform_input(reference, meta, arrays, True)
    if reference_zero.any():
        raise ValueError('Reference contains zero-library cells')
    if (not np.array_equal(np.unique(yt), mapper.classes_) or
            sorted(set(dt)) != sorted(meta['reference_donors']) or
            (meta.get('reference_cells') is not None and set(reference_ids) != set(meta['reference_cells']))):
        raise ValueError('Reference cell, donor or label identities do not match model')
    xq, pq, ids, dq, _, zero = transform_input(query, meta, arrays)
    _disjoint(dq, meta['reference_donors'])
    _disjoint(dq, calibration['donors'])
    if set(ids) & (set(reference_ids) | set(calibration.get('cells', []))):
        raise ValueError('Cell overlap between query and reference/calibration')
    studies = None
    if leave_study_out:
        column = meta.get('spec', {}).get('study_column')
        obs = ad.read_h5ad(reference, backed='r')
        try:
            if column is None or column not in obs.obs:
                raise ValueError('Study omission requires the declared study column in reference')
            studies = _identifiers(obs.obs[column], 'study IDs', len(xt))
        finally:
            obs.file.close()
    threshold = calibration['threshold']
    def predict(model):
        if not (~zero).any():
            return None
        result = _predict(meta, model, xq[~zero], pq[~zero], dq[~zero], calibration['reject_cost'])
        result['accepted'] = result['divergence'] <= threshold
        return result
    baseline = predict(mapper)
    rng = np.random.default_rng(seed)
    donor_ids = np.unique(dt)
    draws = [('bootstrap', str(i), rng.choice(len(donor_ids), len(donor_ids), replace=True))
             for i in range(n_bootstrap)]
    if studies is not None:
        draws += [('leave_study_out', str(study), None) for study in np.unique(studies)]
    values = {'bootstrap': [], 'leave_study_out': []}
    records, donor_draws, diagnostics = [], [], []
    for kind, name, sampled in draws:
        if kind == 'bootstrap':
            rows, draw_donors = _bootstrap_rows(dt, sampled)
            unique, counts = np.unique(sampled, return_counts=True)
            provenance = dict(zip(donor_ids[unique].tolist(), counts.tolist()))
        else:
            rows = np.flatnonzero(studies != name)
            draw_donors, provenance = dt[rows], {'omitted_study': name}
        missing = sorted(set(mapper.classes_) - set(yt[rows]))
        valid = not missing
        records.append({'kind': kind, 'draw': name, 'valid': valid,
            'reason': 'missing_states:' + '|'.join(missing) if missing else '',
            'n_cells': len(rows), 'n_source_donors': len(np.unique(dt[rows])),
            'provenance': json.dumps(provenance, sort_keys=True)})
        if not valid:
            continue
        refit = ReferenceMapper(mapper.beta, mapper.temperature).fit(xt[rows], pt[rows], yt[rows], draw_donors)
        pred = predict(refit)
        # None denotes a valid reference draw with no scoreable query cells.
        values[kind].append(pred)
        frame = pd.DataFrame({'donor_id': dq, 'acceptance': False, 'mean_reference_distance': np.nan})
        if meta['method'] == 'transport':
            frame['mean_reject_mass'] = np.nan
        for state in mapper.classes_:
            frame[f'mean_conditional_affinity::{state}'] = np.nan
        if pred is not None:
            frame.loc[~zero, 'acceptance'] = pred['accepted']
            frame.loc[~zero, 'mean_reference_distance'] = pred['geometric_divergence']
            if meta['method'] == 'transport':
                frame.loc[~zero, 'mean_reject_mass'] = pred['reject_mass']
            for i, state in enumerate(mapper.classes_):
                frame.loc[~zero, f'mean_conditional_affinity::{state}'] = pred['affinity'][:, i]
            diagnostics.extend({'kind': kind, 'draw': name, **d} for d in pred.get('group_diagnostics', []))
        grouped = frame.groupby('donor_id', sort=True)
        donor = grouped.mean()
        donor['n_cells'] = grouped.size()
        donor['kind'], donor['draw'] = kind, name
        donor_draws.append(donor.reset_index())
    ledger = pd.DataFrame(records)
    tables = {'draws': ledger}
    counts = {}
    for kind, prefix in (('bootstrap', 'bootstrap'), ('leave_study_out', 'study')):
        requested = int((ledger.kind == kind).sum())
        valid = len(values[kind])
        counts[kind] = {'requested': requested, 'valid': valid, 'skipped': requested - valid}
        tables[prefix + '_cells'], tables[prefix + '_donors'] = _summary(
            ids, dq, zero, baseline, values[kind], mapper.classes_, requested, meta['method'])
    draw_columns = ['donor_id', 'acceptance', 'mean_reference_distance']
    if meta['method'] == 'transport':
        draw_columns += ['mean_reject_mass']
    draw_columns += [f'mean_conditional_affinity::{c}' for c in mapper.classes_]
    draw_columns += ['n_cells', 'kind', 'draw']
    tables['donor_draws'] = pd.concat(donor_draws, ignore_index=True) if donor_draws else pd.DataFrame(columns=draw_columns)
    record = {'schema': 'bremap-reliability-v1', 'method': meta['method'], 'seed': int(seed), 'counts': counts,
        'reference_sha256': file_hash(reference), 'query_sha256': file_hash(query),
        'model_sha256': file_hash(_path(calibrated_model) / 'model.json'),
        'calibration_sha256': file_hash(_path(calibrated_model) / 'calibration.json'),
        'settings': calibration.get('settings', {'threshold': threshold, 'beta': mapper.beta,
                                                'temperature': mapper.temperature, **meta['transport']}),
        'policy': 'Donor block resampling with distinct repeat IDs; representation/settings/thresholds fixed; scalers/prototypes refitted.',
        'denominators': 'Cell sensitivity excludes zero-library rows; donor sensitivity means defined cell values. Donor draw acceptance includes QC abstentions.',
        'limits': 'Small conditional reference sensitivity analysis, not precision-guaranteed intervals, biological truth or a combined reliability score.',
        'transport_group_dependence': 'Complete donor groups; splitting a donor changes transport predictions.',
        'group_diagnostics': diagnostics, 'versions': _versions(),
        'code_sha256': {name: file_hash(Path(__file__).with_name(name)) for name in ('reliability.py', 'core.py', 'transport.py', 'api.py')}}
    out.mkdir(parents=True)
    for name, table in tables.items():
        table.to_csv(out / (name + '.csv'), index=False)
    record['output_sha256'] = {name + '.csv': file_hash(out / (name + '.csv')) for name in tables}
    _write_json(out / 'reproducibility.json', record)
    return AnalysisResult(out, record, tables)
