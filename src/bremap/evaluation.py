"""Evaluate frozen mapping CSVs without fitting, selecting or modifying a model.

Annotations use exactly cell_id, donor_id, label and optional study_id. Labels
outside the recorded state vocabulary are unknown. F1 averages observed known
states (as in the study evaluator); absent-state recall stays undefined. Brier
is the mean sum of squared errors across all reference states. QC rows enter
acceptance and unknown false-assignment denominators, but not score metrics.
Risk/coverage accepts whole equal-score blocks over score-defined cells only.
Donor/study macros give each group equal weight over defined group metrics.
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import f1_score, roc_auc_score

from .api import AnalysisResult, _json, _path
from .core import _identifiers
from .legacy import _versions, _write_json
from .preprocessing import file_hash


def _boolean(frame, column):
    values = frame[column]
    if not values.isin(['True', 'False']).all():
        raise ValueError(f'{column} must contain True or False')
    return values.eq('True').to_numpy()


def _metrics(frame, classes, score_column):
    labels = frame.label.to_numpy()
    known = np.isin(labels, classes)
    scored = ~frame.qc_zero_library.to_numpy()
    accepted = frame.accepted.to_numpy()
    candidate = frame.candidate_state.to_numpy()
    correct = known & (candidate == labels)
    eligible = known & scored
    metrics = {}
    def add(name, value, denominator, reason=''):
        metrics[name] = {'value': float(value), 'denominator': int(denominator), 'undefined_reason': reason}
    add('acceptance', accepted.mean(), len(frame))
    add('known_acceptance', accepted[known].mean() if known.any() else np.nan, known.sum(),
        '' if known.any() else 'no_known_labels')
    add('known_macro_f1', f1_score(labels[eligible], candidate[eligible], labels=np.unique(labels[eligible]),
        average='macro', zero_division=0) if eligible.any() else np.nan, eligible.sum(),
        '' if eligible.any() else 'no_scored_known_cells')
    affinity = frame[[f'conditional_affinity::{c}' for c in classes]].to_numpy(float)
    add('known_brier', np.mean(np.sum((affinity[eligible] - (labels[eligible, None] == np.asarray(classes)[None, :])) ** 2,
        axis=1)) if eligible.any() else np.nan, eligible.sum(), '' if eligible.any() else 'no_scored_known_cells')
    for state in classes:
        mask = eligible & (labels == state)
        add(f'recall::{state}', correct[mask].mean() if mask.any() else np.nan, mask.sum(),
            '' if mask.any() else ('state_absent' if not (known & (labels == state)).any() else 'no_scored_state_cells'))
    unknown = ~known
    add('unknown_false_assignment_rate', accepted[unknown].mean() if unknown.any() else np.nan, unknown.sum(),
        '' if unknown.any() else 'no_unknown_labels')
    identifiable = known[scored].any() and unknown[scored].any()
    add('unknown_auroc', roc_auc_score(unknown[scored], frame.loc[scored, score_column]) if identifiable else np.nan,
        scored.sum(), '' if identifiable else 'requires_known_and_unknown_scored_cells')
    add('selective_risk', (~correct[accepted]).mean() if accepted.any() else np.nan, accepted.sum(),
        '' if accepted.any() else 'no_accepted_cells')
    curve = pd.DataFrame(columns=['score_threshold', 'n_accepted', 'coverage', 'risk'])
    if scored.any():
        # Same blockwise formula as audited mapping_metrics; no arbitrary ordering
        # can split tied scores. QC scores are undefined and excluded explicitly.
        scores, group = np.unique(frame.loc[scored, score_column].to_numpy(float), return_inverse=True)
        totals = np.bincount(group)
        mistakes = np.bincount(group, weights=(~correct[scored]).astype(float))
        cumulative = np.cumsum(totals)
        coverage, risks = cumulative / scored.sum(), np.cumsum(mistakes) / cumulative
        curve = pd.DataFrame({'score_threshold': scores, 'n_accepted': cumulative, 'coverage': coverage, 'risk': risks})
        at80 = np.flatnonzero(coverage >= .8)[0]
        add('risk_coverage_auc', np.sum(risks * totals / scored.sum()), scored.sum())
        add('risk_at_nominal_80pct_coverage', risks[at80], cumulative[at80])
        add('achieved_80pct_coverage', coverage[at80], scored.sum())
    else:
        for metric in ('risk_coverage_auc', 'risk_at_nominal_80pct_coverage', 'achieved_80pct_coverage'):
            add(metric, np.nan, 0, 'no_scored_cells')
    return metrics, curve


def evaluate(mapping, annotations, out_dir):
    """Return EvaluationResult(path, metadata, tables) for a saved v2 mapping.

    ``mapping`` accepts a MappingResult or directory. ``annotations`` is a CSV
    with a one-to-one exact cell/donor join; optional study_id must be complete
    and each donor must belong to one study. Inputs remain byte-unchanged.
    """
    source, out = _path(mapping), _path(out_dir)
    if out.exists():
        raise ValueError('Output already exists; choose a new evaluation directory')
    record = _json(source / 'reproducibility.json')
    if record.get('schema') != 'bremap-mapping-v2':
        raise ValueError('Evaluation requires a frozen bremap-mapping-v2 artifact')
    checksums = record.get('output_sha256', {})
    if 'cells.csv' not in checksums:
        raise ValueError('Mapping cells checksum is missing')
    for name, digest in checksums.items():
        if Path(name).name != name or file_hash(source / name) != digest:
            raise ValueError(f'Mapping checksum mismatch: {name}')
    classes = record['spec']['states']
    if (not isinstance(classes, list) or len(classes) < 2 or len(set(classes)) != len(classes)
            or any(not isinstance(c, str) or not c or c == 'unknown' for c in classes)):
        raise ValueError('Invalid mapping state vocabulary')
    method = record.get('method')
    if method not in ('affinity', 'transport'):
        raise ValueError('Unsupported mapping method')
    score_column = 'reject_mass' if method == 'transport' else 'reference_distance'
    if record['settings']['threshold_variable'] != score_column:
        raise ValueError('Mapping score variable contradicts method')
    cells = pd.read_csv(source / 'cells.csv', dtype=str, keep_default_na=False)
    annotation = pd.read_csv(annotations, dtype=str, keep_default_na=False)
    required = {'cell_id', 'donor_id', 'candidate_state', 'accepted_analog', 'accepted', 'qc_zero_library',
                'abstention_reason', score_column, *[f'conditional_affinity::{c}' for c in classes]}
    if not required <= set(cells) or not {'cell_id', 'donor_id', 'label'} <= set(annotation):
        raise ValueError('Mapping/annotation CSV is missing required cell_id, donor_id, label or prediction columns')
    if len(cells) == 0:
        raise ValueError('Mapping contains no cells')
    # Preserve literal identifiers such as NA and reject ambiguous whitespace.
    for frame, columns in ((cells, ('cell_id', 'donor_id')), (annotation, ('cell_id', 'donor_id', 'label'))):
        for column in columns:
            normalized = _identifiers(frame[column], column, len(frame))
            if not np.array_equal(normalized, frame[column].to_numpy()):
                raise ValueError(f'{column} must not contain surrounding whitespace')
        if frame.cell_id.duplicated().any():
            raise ValueError('Duplicate cell IDs in mapping or annotations')
    if set(cells.cell_id) != set(annotation.cell_id):
        raise ValueError('Cell IDs must match exactly between mapping and annotations')
    annotation = annotation.set_index('cell_id').loc[cells.cell_id].reset_index()
    if not np.array_equal(annotation.donor_id, cells.donor_id):
        raise ValueError('Donor IDs disagree between mapping and annotations')
    cells['label'] = annotation.label.to_numpy()
    if 'study_id' in annotation:
        normalized = _identifiers(annotation.study_id, 'study IDs', len(cells))
        if not np.array_equal(normalized, annotation.study_id):
            raise ValueError('Study IDs must not contain surrounding whitespace')
        cells['study_id'] = annotation.study_id.to_numpy()
        if (cells.groupby('donor_id').study_id.nunique() != 1).any():
            raise ValueError('Each donor must belong to exactly one study')
    cells['accepted'] = _boolean(cells, 'accepted')
    cells['qc_zero_library'] = _boolean(cells, 'qc_zero_library')
    zero = cells.qc_zero_library.to_numpy()
    for name in [score_column, *[f'conditional_affinity::{c}' for c in classes]]:
        # NumPy's float conversion preserves round-trip float64 strings; pandas
        # to_numeric can merge adjacent scores and thereby alter ranking/ties.
        cells[name] = cells[name].mask(cells[name].eq('')).astype(np.float64)
        if not np.isfinite(cells.loc[~zero, name]).all() or not cells.loc[zero, name].isna().all():
            raise ValueError('Prediction scores must be finite on valid cells and undefined on zero-library cells')
    affinity = cells.loc[~zero, [f'conditional_affinity::{c}' for c in classes]].to_numpy(float)
    if np.any(affinity < 0) or np.any(affinity > 1) or not np.allclose(affinity.sum(axis=1), 1, atol=1e-8, rtol=0):
        raise ValueError('Conditional affinities must be probabilities summing to one')
    if (cells.loc[~zero, score_column] < 0).any() or (method == 'transport' and (cells.loc[~zero, score_column] > 1).any()):
        raise ValueError('Invalid mapping score range')
    if (not cells.loc[~zero, 'candidate_state'].isin(classes).all() or cells.loc[zero, 'accepted'].any()
            or not cells.loc[zero, 'abstention_reason'].eq('zero_library').all()
            or not np.array_equal(cells.accepted_analog.to_numpy(), np.where(cells.accepted, cells.candidate_state, 'unknown'))):
        raise ValueError('Mapping candidates, accepted analogs or QC abstentions are inconsistent')
    tables, summaries, curves = {'cells': cells}, [], []
    levels = [('cell_pooled', None), ('donor', 'donor_id')]
    if 'study_id' in cells:
        levels.append(('study', 'study_id'))
    for level, column in levels:
        rows = []
        groups = [('', cells)] if column is None else list(cells.groupby(column, sort=True))
        for group_id, frame in groups:
            metrics, curve = _metrics(frame, classes, score_column)
            for metric, values in metrics.items():
                rows.append({'aggregation': level, 'group_id': group_id, 'metric': metric,
                             'n_cells': len(frame), 'n_scored': int((~frame.qc_zero_library).sum()), **values})
            curve['aggregation'], curve['group_id'] = level, group_id
            curves.append(curve)
        unit = pd.DataFrame(rows)
        if column is None:
            summaries.extend(rows)
        else:
            tables['donors' if level == 'donor' else 'studies'] = unit
            for metric, frame in unit.groupby('metric', sort=False):
                defined = frame.value.notna()
                summaries.append({'aggregation': level + '_macro', 'group_id': '', 'metric': metric,
                    'value': frame.loc[defined, 'value'].mean(), 'denominator': int(defined.sum()),
                    'n_groups_total': len(frame), 'n_groups_defined': int(defined.sum()),
                    'n_cells': len(cells), 'n_scored': int((~zero).sum()),
                    'undefined_reason': '' if defined.any() else 'no_defined_groups'})
    tables['summary'] = pd.DataFrame(summaries)
    tables['risk_coverage'] = pd.concat(curves, ignore_index=True)
    metadata = {'schema': 'bremap-evaluation-v1', 'mapping_sha256': file_hash(source / 'reproducibility.json'),
        'mapping_cells_sha256': file_hash(source / 'cells.csv'), 'annotations_sha256': file_hash(annotations),
        'states': classes, 'unknown_policy': 'Any annotation label outside the frozen reference vocabulary is unknown.',
        'score_variable': score_column, 'n_cells': len(cells), 'n_scored': int((~zero).sum()),
        'denominators': {'acceptance': 'All cells, including zero-library QC abstentions.',
            'known_macro_f1': 'Score-defined known cells; equal mean over observed known states.',
            'known_brier': 'Score-defined known cells; sum over all frozen reference states.',
            'unknown_auroc': 'Score-defined cells; unknown positive; larger divergence means unknown.',
            'unknown_false_assignment_rate': 'All unknown cells, including QC abstentions.',
            'selective_risk': 'Accepted cells; unknown assignments count as errors.',
            'risk_coverage': 'Score-defined cells only; accept complete tied-score blocks, right-endpoint area.',
            'group_macros': 'Equal mean of defined donor/study metrics; denominator is defined groups; group counts explicit.'},
        'limits': 'Frozen-prediction evaluation; no model fitting/tuning, biological truth claim or combined reliability score.',
        'versions': _versions(), 'code_sha256': {'evaluation.py': file_hash(Path(__file__))}}
    out.mkdir(parents=True)
    for name, table in tables.items():
        table.to_csv(out / (name + '.csv'), index=False)
    metadata['output_sha256'] = {name + '.csv': file_hash(out / (name + '.csv')) for name in tables}
    _write_json(out / 'reproducibility.json', metadata)
    return AnalysisResult(out, metadata, tables)
