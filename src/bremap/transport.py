"""Donor-local entropic transport with relaxed state marginals and reject mass.

This is a standard generalized Sinkhorn construction, not a novelty claim.
See docs/revision/TRANSPORT_OBJECTIVE.md for the exact objective and limitations.
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy.special import logsumexp, softmax, xlogy

from .core import _finite_number, _identifiers, _matrix


def transport_plan(cost, reject_cost, entropy=.1, marginal_penalty=.2,
                   reject_prior=.1, max_iter=2000, tol=1e-8):
    """Minimize cost + entropy KL(T|a*b) + penalty KL(real columns|b_real).

    Rows have mass 1/n; b=((1-reject_prior)/K,...,reject_prior). The last
    column is rejection and has no marginal KL penalty. Returned reject_mass
    is a per-cell fraction, n*T[:, -1], not an unnormalized joint mass.
    Conditional affinities sum to one over real states even for rejected cells.
    No threshold selection or fitting happens here; nonconvergence raises.
    """
    cost = _matrix(cost, "cost")
    reject_cost = _finite_number(reject_cost, "reject_cost")
    entropy = _finite_number(entropy, "entropy")
    marginal_penalty = _finite_number(marginal_penalty, "marginal_penalty")
    reject_prior = _finite_number(reject_prior, "reject_prior")
    tol = _finite_number(tol, "tol")
    if np.any(cost < 0) or reject_cost < 0:
        raise ValueError("cost and reject_cost must be nonnegative")
    if entropy <= 0 or marginal_penalty < 0 or tol <= 0:
        raise ValueError("entropy and tol must be positive; marginal_penalty nonnegative")
    if not 0 < reject_prior < 1:
        raise ValueError("reject_prior must lie strictly between zero and one")
    if (isinstance(max_iter, (bool, np.bool_))
            or not isinstance(max_iter, (int, np.integer)) or max_iter < 1):
        raise ValueError("max_iter must be a positive integer")

    n, k = cost.shape
    b = np.r_[np.full(k, (1 - reject_prior) / k), reject_prior]
    augmented = np.column_stack([cost, np.full(n, reject_cost)])
    log_a = -np.log(n)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        log_kernel = log_a + np.log(b) - augmented / entropy
    if not np.isfinite(log_kernel).all():
        raise ValueError("Transport costs overflowed relative to entropy")
    # Subtract row constants: exact row constraints make these immaterial.
    log_kernel -= log_kernel.max(axis=1, keepdims=True)
    log_v = np.zeros(k + 1)
    relaxation = 1 / (1 + entropy / marginal_penalty) if marginal_penalty else 0.
    for iteration in range(1, max_iter + 1):
        log_u = log_a - logsumexp(log_kernel + log_v, axis=1)
        updated = np.zeros(k + 1)
        updated[:k] = relaxation * (np.log(b[:k]) -
            logsumexp(log_kernel[:, :k] + log_u[:, None], axis=0))
        error = float(np.max(np.abs(updated - log_v)))
        log_v = updated
        if error <= tol:
            break
    else:
        raise RuntimeError(f"Transport did not converge after {max_iter} iterations "
                           f"(log scaling residual={error:.3g})")

    log_u = log_a - logsumexp(log_kernel + log_v, axis=1)
    plan = np.exp(log_kernel + log_u[:, None] + log_v)
    row_error = float(np.max(np.abs(plan.sum(axis=1) - 1 / n)))
    columns = plan[:, :k].sum(axis=0)
    base = b[None, :] / n
    objective = float((plan * augmented).sum()
        + entropy * (xlogy(plan, plan / base) - plan + base).sum()
        + marginal_penalty * (xlogy(columns, columns / b[:k]) - columns + b[:k]).sum())
    affinity = softmax(log_kernel[:, :k] + log_v[:k], axis=1)
    if (not np.isfinite(objective) or not np.isfinite(plan).all()
            or not np.isfinite(affinity).all() or row_error > tol):
        raise RuntimeError("Transport produced nonfinite output or failed row convergence")
    return {
        "plan": plan,
        "conditional_affinity": affinity,
        "reject_mass": np.clip(n * plan[:, -1], 0, 1),
        "n_iter": iteration,
        "converged": True,
        "row_error": row_error,
        "scaling_error": error,
        "objective": objective,
    }


def _mapper_cost(mapper, X, programs):
    """Reuse the reference-fitted transforms and the existing mapper cost."""
    if not hasattr(mapper, "classes_"):
        raise ValueError("Fit a reference before transport prediction")
    X, programs = _matrix(X, "X"), _matrix(programs, "programs")
    if len(X) != len(programs):
        raise ValueError("X and programs must have the same number of rows")
    if programs.shape[1] != mapper.program_scaler_.n_features_in_:
        raise ValueError("Query program dimensions must match training reference")
    Z = mapper.transform_embedding(X)
    P = mapper.program_scaler_.transform(programs)
    embedding_cost = np.clip(1 - Z @ mapper.embedding_centroids_.T, 0, 2)
    program_cost = cdist(P, mapper.program_prototypes_, metric="sqeuclidean") / P.shape[1]
    cost = (1 - mapper.beta) * embedding_cost + mapper.beta * program_cost
    if not np.isfinite(cost).all():
        raise ValueError("Query distances overflowed; inspect input magnitudes")
    return cost, embedding_cost


def transport_predict(mapper, X, programs, groups, reject_cost, entropy=.1,
                      marginal_penalty=.2, reject_prior=.1):
    """Solve each query donor separately using frozen reference transforms.

    Affinity is conditional on real-state matching; divergence is reject fraction.
    Geometric and embedding distances are separate diagnostics. This function
    does not threshold reject fractions; callers calibrate abstention on held-out
    reference donors. beta=0 gives the same transport without biological costs.
    """
    cost, embedding_cost = _mapper_cost(mapper, X, programs)
    groups = _identifiers(groups, "groups", len(cost))
    affinity = np.empty_like(cost)
    reject_mass = np.empty(len(cost))
    diagnostics = []
    for group in np.unique(groups):
        rows = groups == group
        result = transport_plan(cost[rows], reject_cost, entropy=entropy,
                                marginal_penalty=marginal_penalty, reject_prior=reject_prior)
        affinity[rows] = result["conditional_affinity"]
        reject_mass[rows] = result["reject_mass"]
        diagnostics.append({"group": str(group), "n_cells": int(rows.sum()),
                            **{key: result[key] for key in
                               ("n_iter", "converged", "row_error", "scaling_error", "objective")}})
    candidate = mapper.classes_[affinity.argmax(axis=1)]
    return {
        "affinity": affinity,
        "candidate_analog": candidate,
        "best_analog": candidate.copy(),
        "analog_strength": affinity.max(axis=1),
        "divergence": reject_mass,
        "reject_mass": reject_mass,
        "matched_mass": 1 - reject_mass,
        "geometric_divergence": cost.min(axis=1),
        "embedding_divergence": embedding_cost.min(axis=1),
        "abstain_flag": np.zeros(len(cost), dtype=bool),
        "group_diagnostics": diagnostics,
    }
