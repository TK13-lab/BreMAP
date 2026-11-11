"""Frozen calibration policies, sharing the audited study selection helpers."""
from copy import copy
from importlib.resources import files
from itertools import product
import json

import numpy as np
from scipy.special import softmax
from sklearn.metrics import f1_score

from .transport import transport_predict


def donor_quantile(values, donors, quantile):
    values, donors = np.asarray(values, dtype=float), np.asarray(donors)
    if values.ndim != 1 or len(values) != len(donors) or not len(values) or not np.isfinite(values).all():
        raise ValueError('Finite, nonempty matched scores/donors required')
    if not 0 <= quantile <= 1:
        raise ValueError('Quantile must be in [0,1]')
    _, inverse, counts = np.unique(donors, return_inverse=True, return_counts=True)
    weights = 1 / counts[inverse]
    order = np.argsort(values, kind='stable')
    cumulative = np.cumsum(weights[order]) / weights.sum()
    index = min(np.searchsorted(cumulative, quantile - 1e-12), len(values) - 1)
    return float(values[order[index]])

def selection_score(prediction, labels, classes):
    known = np.isin(labels, classes)
    if not known.any():
        raise ValueError('Calibration requires known states')
    f1 = f1_score(labels[known], prediction['candidate_analog'][known], labels=classes,
                  average='macro', zero_division=0)
    target = labels[known, None] == classes[None, :]
    brier = np.mean(np.sum((prediction['affinity'][known] - target) ** 2, axis=1))
    return float(f1), -float(brier)

def select_transport(mapper, xcal, pcal, labels, donors, config):
    base = mapper.predict(xcal, pcal)
    known = np.isin(labels, mapper.classes_)
    reject_cost = donor_quantile(base['divergence'][known], donors[known], config['reject_cost_known_quantile'])
    candidates = []
    for entropy, penalty in product(config['entropy_grid'], config['marginal_penalty_grid']):
        pred = transport_predict(mapper, xcal, pcal, donors, reject_cost,
                                 entropy=entropy, marginal_penalty=penalty,
                                 reject_prior=config['reject_prior'])
        f1, negative_brier = selection_score(pred, labels, mapper.classes_)
        candidates.append((f1, negative_brier, -penalty, -entropy, pred, entropy, penalty))
    selected = max(candidates, key=lambda item: item[:4])
    return selected[4], {'reject_cost': reject_cost, 'entropy': selected[5],
                         'marginal_penalty': selected[6], 'reject_prior': config['reject_prior'],
                         'calibration_macro_f1': selected[0], 'calibration_brier': -selected[1]}

def _temperature(probabilities, temperature, original=.1):
    return softmax(np.log(np.maximum(probabilities, 1e-300)) * original / temperature, axis=1)

def _tune_temperature(cal, y, classes, grid):
    known = np.isin(y, classes)
    target = y[known, None] == classes[None, :]
    return min(grid, key=lambda t: np.mean(np.sum((_temperature(cal['affinity'][known], t) - target) ** 2, axis=1)))


def load_policy(name):
    if name == 'normal_donor_v1':
        return json.loads(files('bremap').joinpath('resources/normal_donor_v1.json').read_text())
    if name == 'legacy_cli_v1':
        return {'policy_id': name, 'version': 1, 'beta': .25, 'temperature': .1,
                'entropy': .1, 'marginal_penalty': .2, 'reject_prior': .1,
                'known_acceptance': .95, 'reject_cost_known_quantile': .95,
                'quantile': 'cell_weighted_numpy_linear', 'acceptance_comparison': 'score <= threshold'}
    raise ValueError(f'Unsupported calibration policy: {name}')


def select_settings(mapper, X, programs, labels, donors, method, cost_mode, policy):
    """Select only calibration settings; do not refit frozen reference parameters."""
    if method not in ('affinity', 'transport') or cost_mode not in ('biology', 'embedding'):
        raise ValueError('Unsupported method or cost mode')
    model = copy(mapper)
    if policy['policy_id'] == 'legacy_cli_v1':
        model.beta = 0. if cost_mode == 'embedding' else policy['beta']
        model.temperature = policy['temperature']
        pred = model.predict(X, programs)
        reject_cost = float(np.quantile(pred['divergence'], .95)) if method == 'transport' else None
        selected = {'beta': model.beta, 'temperature': model.temperature,
                    'reject_cost': reject_cost, 'entropy': policy['entropy'],
                    'marginal_penalty': policy['marginal_penalty'], 'reject_prior': policy['reject_prior']}
        if method == 'transport':
            pred = transport_predict(model, X, programs, donors, **{k: selected[k] for k in
                                     ('reject_cost', 'entropy', 'marginal_penalty', 'reject_prior')})
        threshold = float(np.quantile(pred['divergence'], .95))
    else:
        candidates = []
        for beta in ([0.] if cost_mode == 'embedding' else policy['transfer_beta_grid']):
            candidate = copy(model)
            candidate.beta, candidate.temperature = beta, .1
            pred = candidate.predict(X, programs)
            candidate.temperature = _tune_temperature(pred, labels, candidate.classes_, policy['transfer_temperature_grid'])
            pred = candidate.predict(X, programs)
            f1, negative_brier = selection_score(pred, labels, candidate.classes_)
            candidates.append((f1, negative_brier, -beta, candidate, pred))
        f1, negative_brier, _, model, pred = max(candidates, key=lambda item: item[:3])
        selected = {'beta': model.beta, 'temperature': model.temperature, 'reject_cost': None,
                    'calibration_macro_f1': f1, 'calibration_brier': -negative_brier}
        if method == 'transport':
            pred, transport = select_transport(model, X, programs, labels, donors, policy)
            selected.update(transport)
        threshold = donor_quantile(pred['divergence'], donors, policy['known_acceptance'])
    selected.update(threshold=threshold,
                    threshold_variable='reject_mass' if method == 'transport' else 'reference_distance')
    return selected, pred
