"""Explicit expression semantics, without reconstructing raw counts."""
from numbers import Real
import numpy as np
from scipy import sparse

from .preprocessing import normalize_counts

KEY = 'bremap_expression'
CANONICAL = 'CP10k_log1p_on_reference_gene_universe'
TYPES = ('counts', 'library_normalized', 'log1p_library_normalized', 'preprocessed')


def declaration(kind='counts', *, log_base='e', processing_id=None, units=None):
    if kind not in TYPES:
        raise ValueError('expression_type must be ' + ', '.join(TYPES))
    result = {'schema': 'bremap-expression-v1', 'type': kind}
    if kind == 'log1p_library_normalized':
        base = float(np.e) if log_base == 'e' else log_base
        if isinstance(base, (bool, np.bool_)) or not isinstance(base, Real) or not np.isfinite(base) or base <= 1:
            raise ValueError('log_base must be "e" or a finite number greater than 1 (e.g. 2)')
        result['log_base'] = float(base)
    elif log_base != 'e':
        raise ValueError('log_base only applies to log1p_library_normalized expression')
    if kind == 'preprocessed':
        for name, value in [('processing_id', processing_id), ('units', units)]:
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f'preprocessed expression requires nonempty {name}; describe the actual transform')
            result[name] = value.strip()
    elif processing_id is not None or units is not None:
        raise ValueError('processing_id/units apply only to preprocessed expression')
    return result


def get_declaration(data):
    value = data.uns.get(KEY)
    if value is None:
        return declaration()
    if not isinstance(value, dict) or value.get('schema') != 'bremap-expression-v1':
        raise ValueError('Invalid expression declaration; load with read_data and explicit expression_type')
    expected = declaration(value.get('type'), log_base=value.get('log_base', 'e'),
        processing_id=value.get('processing_id'), units=value.get('units'))
    if set(value) - (set(expected) | {'source'}):
        raise ValueError('Unknown expression declaration fields')
    return expected


def effective_contract(data):
    info = get_declaration(data)
    if info['type'] != 'preprocessed':
        return {'normalization': CANONICAL}
    return {'normalization': 'provided_expression', 'processing_id': info['processing_id'],
            'units': info['units']}


def check_contract(data, expected):
    validate_contract(expected)
    if effective_contract(data) != expected:
        raise ValueError('Expression processing does not match the fitted reference/student: '
                         'use the same declared processing or fit a separate model')


def validate_contract(contract):
    if contract == {'normalization': CANONICAL}:
        return
    if not isinstance(contract, dict) or set(contract) != {'normalization','processing_id','units'} or contract.get('normalization') != 'provided_expression':
        raise ValueError('Unsupported expression processing contract')
    declaration('preprocessed', processing_id=contract['processing_id'], units=contract['units'])


def values(data):
    """Validated counts, relative linear expression, or declared as-supplied values.

    Validate before feature selection. Do not convert normalized values to integers.
    The input object/matrix is never changed.
    """
    info = get_declaration(data)
    if data.X is None:
        raise ValueError('No expression matrix in X; select a layer or raw.X with read_data')
    x = sparse.csr_matrix(data.X, dtype=np.float64, copy=True)
    if not np.isfinite(x.data).all():
        raise ValueError('Expression must contain finite numeric values; NaN/Inf are not imputed')
    if info['type'] != 'preprocessed' and (x.data < 0).any():
        raise ValueError('Counts/library-normalized expression must be nonnegative; residuals need preprocessed mode')
    if info['type'] == 'counts' and not np.allclose(x.data, np.rint(x.data), atol=1e-5, rtol=0):
        raise ValueError('Raw integer counts required; select counts or declare the known expression_type; do not round')
    if info['type'] == 'log1p_library_normalized':
        with np.errstate(over='ignore', invalid='ignore'):
            x.data = np.expm1(x.data * np.log(info['log_base']))
        if not np.isfinite(x.data).all():
            raise ValueError('Inverse log1p overflow: check log_base and expression_type; scaled values are not logcounts')
    x.sum_duplicates(); x.sort_indices(); x.eliminate_zeros()
    return x


def features(data, positions=None):
    """Return program/SVD features on the selected reference gene universe."""
    x = values(data)
    if positions is not None:
        x = x[:, positions]
    # Alignment can reorder CSR column indices. Restore the historical order
    # before row sums; roundoff otherwise changes decisions at a threshold.
    x.sort_indices()
    kind = get_declaration(data)['type']
    if kind == 'preprocessed':
        return x
    if kind == 'counts':
        return normalize_counts(x)  # Keep the historical count path unchanged.
    # Relative libraries may sum to less than one. Unlike integer counts, use
    # every positive total as-is; clipping to 1 would alter their expression.
    totals = np.asarray(x.sum(axis=1)).ravel()
    if not np.isfinite(totals).all():
        raise ValueError('Expression library totals overflow; check input scale')
    factors = np.divide(10000., totals, out=np.zeros_like(totals), where=totals > 0)
    x = (sparse.diags(factors) @ x).tocsr()
    x.data = np.log1p(x.data)
    if not np.isfinite(x.data).all():
        raise ValueError('Nonfinite normalized expression; check the supplied numeric scale')
    return x


def zero_rows(matrix):
    x = sparse.csr_matrix(matrix, copy=True)
    x.eliminate_zeros()
    return np.diff(x.indptr) == 0


def score_units(contract):
    return 'mean ' + (contract['units'] if contract['normalization'] == 'provided_expression' else 'log1p(CP10k)')
