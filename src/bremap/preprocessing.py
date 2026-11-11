"""Shared count normalization, program scoring and reference-only raw SVD.

These functions are extracted unchanged from the historical study workflow.
"""

import hashlib

import numpy as np
from scipy import sparse
from sklearn.decomposition import TruncatedSVD


def file_hash(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def normalize_counts(counts):
    x = sparse.csr_matrix(counts, dtype=np.float64)
    if not np.isfinite(x.data).all() or (x.data < 0).any():
        raise ValueError('Counts must be finite and nonnegative')
    if not np.allclose(x.data, np.rint(x.data), atol=1e-5, rtol=0):
        raise ValueError('Expected raw integer counts, not a normalized matrix')
    totals = np.asarray(x.sum(axis=1)).ravel()
    x = sparse.diags(10000. / np.maximum(totals, 1)) @ x
    x.data = np.log1p(x.data)
    return x.tocsr()


def program_scores(logcounts, genes, config):
    anchors = {g.upper() for v in config['anchors'].values() for g in v}
    heldout = {g.upper() for v in config['heldout_validation'].values() for g in v}
    if anchors & heldout:
        raise ValueError(f'Anchor/validation overlap: {sorted(anchors & heldout)}')
    genes = np.asarray([str(g).upper() for g in genes])
    scores, coverage = [], {}
    for name, symbols in config['anchors'].items():
        mask = np.isin(genes, symbols)
        coverage[name] = sorted(set(genes[mask]))
        if len(coverage[name]) < config['minimum_genes_per_anchor']:
            raise ValueError(f'Insufficient anchor coverage for {name}')
        scores.append(np.asarray(logcounts[:, mask].mean(axis=1)).ravel())
    return np.column_stack(scores), coverage


def fit_raw_svd(train, query, components=32, seed=42):
    k = min(components, train.shape[0] - 1, train.shape[1] - 1)
    if k < 2:
        raise ValueError('SVD requires at least three training cells and genes')
    model = TruncatedSVD(n_components=k, random_state=seed)
    return model.fit_transform(train), model.transform(query), model
