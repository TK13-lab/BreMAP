"""Train-only donor-balanced affinity mapping (M1/M2a), not optimal transport.

Affinity is conditional on the available normal reference states. Divergence is
a distance, not malignancy probability. Calibrate rejection on separate donors
outside this module; neither fit nor predict selects weights or thresholds.
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy.special import softmax
from sklearn.preprocessing import StandardScaler, normalize


def _finite_number(value, name):
    if not np.isscalar(value):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not np.isfinite(result):
        raise ValueError(f"{name} must be a finite number")
    return result


def _matrix(values, name):
    try:
        array = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite numeric matrix") from exc
    if array.ndim != 2 or 0 in array.shape or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a nonempty finite 2D matrix")
    return array


def _identifiers(values, name, n):
    array = np.asarray(values, dtype=object)
    if array.ndim != 1 or len(array) != n:
        raise ValueError(f"{name} must be a 1D array with {n} entries")
    result = []
    for value in array:
        text = str(value).strip()
        if value is None or text.lower() in {"", "nan", "none", "<na>", "nat"}:
            raise ValueError(f"{name} contains a missing value")
        result.append(text)
    return np.asarray(result)


class ReferenceMapper:
    """Blend cosine embedding cost and standardized program squared distance.

    Both StandardScalers fit training reference only, with equal total weight per
    donor. Embeddings are then row-normalized; state centroids average donor-state
    means equally, then normalize again. Program prototypes average standardized
    donor-state means equally without L2 normalization.

    Joint cost = (1-beta) * cosine_distance + beta * mean_squared_program_distance.
    Cosine distance is bounded by [0, 2]; program distance is unbounded. Their
    scales deliberately differ: beta is a cost weight, not a fraction of evidence.
    Select beta and temperature on disjoint calibration data. beta=0 is embedding
    only, beta=1 is standardized signature distance only.
    """

    def __init__(self, beta=0.5, temperature=0.1):
        beta = _finite_number(beta, "beta")
        temperature = _finite_number(temperature, "temperature")
        if not 0 <= beta <= 1:
            raise ValueError("beta must be finite and lie in [0, 1]")
        if temperature <= 0:
            raise ValueError("temperature must be finite and positive")
        self.beta = beta
        self.temperature = temperature

    def fit(self, X, programs, labels, donors):
        """Fit training reference; each donor contributes equally to each state."""
        X, programs = _matrix(X, "X"), _matrix(programs, "programs")
        if len(programs) != len(X):
            raise ValueError("X and programs must have the same number of rows")
        labels = _identifiers(labels, "labels", len(X))
        donors = _identifiers(donors, "donors", len(X))
        classes = np.unique(labels)
        if len(classes) < 2:
            raise ValueError("Reference must contain at least two states")
        donor_ids, donor_index, counts = np.unique(donors, return_inverse=True,
                                                   return_counts=True)
        # Equal donor totals also make duplicating a donor's cells invariant.
        weights = 1. / counts[donor_index]
        embedding_scaler = StandardScaler().fit(X, sample_weight=weights)
        program_scaler = StandardScaler().fit(programs, sample_weight=weights)
        Z = normalize(embedding_scaler.transform(X))
        P = program_scaler.transform(programs)
        embedding_centroids, program_prototypes = [], []
        donor_embeddings, donor_programs, state_labels, state_donors = [], [], [], []
        for label in classes:
            state_embedding, state_program = [], []
            in_state = labels == label
            for donor in np.unique(donors[in_state]):
                mask = in_state & (donors == donor)
                state_embedding.append(Z[mask].mean(axis=0))
                state_program.append(P[mask].mean(axis=0))
                state_labels.append(label)
                state_donors.append(donor)
            donor_embeddings.extend(state_embedding)
            donor_programs.extend(state_program)
            embedding_centroids.append(np.mean(state_embedding, axis=0))
            program_prototypes.append(np.mean(state_program, axis=0))
        # Assign only after successful fit, and retain no aliases to input arrays.
        self.classes_ = classes
        self.donors_ = donor_ids
        self.embedding_scaler_ = embedding_scaler
        self.program_scaler_ = program_scaler
        self.embedding_centroids_ = normalize(np.asarray(embedding_centroids))
        self.program_prototypes_ = np.asarray(program_prototypes)
        self.donor_embedding_centroids_ = np.asarray(donor_embeddings)
        self.donor_program_prototypes_ = np.asarray(donor_programs)
        self.donor_state_labels_ = np.asarray(state_labels)
        self.donor_state_donors_ = np.asarray(state_donors)
        return self

    def transform_embedding(self, X):
        """Apply training-fitted embedding scaling and row L2 normalization."""
        if not hasattr(self, "classes_"):
            raise ValueError("Fit a reference before transforming embeddings")
        X = _matrix(X, "X")
        if X.shape[1] != self.embedding_scaler_.n_features_in_:
            raise ValueError("Query embedding dimensions must match training reference")
        return normalize(self.embedding_scaler_.transform(X))

    def predict(self, X, programs, threshold=None):
        """Return conditional affinities and reject if minimum cost > threshold.

        ``None`` disables rejection. ``candidate_analog`` always retains the
        minimum-cost state; ``best_analog`` becomes ``unknown`` on rejection.
        ``embedding_divergence`` is the independent minimum embedding cost, even
        when the state minimizing joint cost differs. No query data are fitted.
        """
        if not hasattr(self, "classes_"):
            raise ValueError("Fit a reference before prediction")
        if threshold is not None:
            threshold = _finite_number(threshold, "threshold")
            if threshold < 0:
                raise ValueError("threshold must be nonnegative, or None")
        X, programs = _matrix(X, "X"), _matrix(programs, "programs")
        if len(X) != len(programs):
            raise ValueError("X and programs must have the same number of rows")
        if (X.shape[1] != self.embedding_scaler_.n_features_in_
                or programs.shape[1] != self.program_scaler_.n_features_in_):
            raise ValueError("Query feature dimensions must match training reference")
        Z = self.transform_embedding(X)
        P = self.program_scaler_.transform(programs)
        embedding_cost = np.clip(1 - Z @ self.embedding_centroids_.T, 0, 2)
        program_cost = cdist(P, self.program_prototypes_, metric="sqeuclidean") / P.shape[1]
        cost = (1 - self.beta) * embedding_cost + self.beta * program_cost
        if not np.isfinite(cost).all():
            raise ValueError("Query distances overflowed; inspect input magnitudes")
        # Subtract the row minimum before division, including for tiny temperatures.
        with np.errstate(over="ignore"):
            affinity = softmax(-(cost - cost.min(axis=1, keepdims=True)) / self.temperature,
                               axis=1)
        best_index = np.argmin(cost, axis=1)
        divergence = cost[np.arange(len(X)), best_index]
        rejected = np.zeros(len(X), dtype=bool) if threshold is None else divergence > threshold
        candidate = self.classes_[best_index]
        best = np.where(rejected, "unknown", candidate)
        return {
            "affinity": affinity,
            "best_analog": best,
            "candidate_analog": candidate,
            "analog_strength": affinity.max(axis=1),
            "divergence": divergence,
            "embedding_divergence": embedding_cost.min(axis=1),
            "abstain_flag": rejected,
        }


def donor_bootstrap_stability(mapper, X_train, programs_train, labels_train,
                              donors_train, X_query, programs_query,
                              n_bootstrap=20, seed=42):
    """Resample reference donor blocks and refit scaling/prototypes per draw.

    Stability is agreement with the original reference's *candidate* analog, not
    confidence in biological origin or rejection. Affinity SD is population SD
    across valid draws. Missing-state draws are skipped and counted explicitly;
    no valid draws returns NaN estimates. Repeated donor draws are separate block
    copies, preserving multiplicity during both scaling and prototype averaging.
    """
    if (isinstance(n_bootstrap, (bool, np.bool_))
            or not isinstance(n_bootstrap, (int, np.integer)) or n_bootstrap < 1):
        raise ValueError("n_bootstrap must be a positive integer")
    baseline = mapper.predict(X_query, programs_query)
    X = _matrix(X_train, "X_train")
    programs = _matrix(programs_train, "programs_train")
    if len(X) != len(programs):
        raise ValueError("Training matrices must have the same number of rows")
    labels = _identifiers(labels_train, "labels_train", len(X))
    donors = _identifiers(donors_train, "donors_train", len(X))
    if not np.array_equal(np.unique(labels), mapper.classes_):
        raise ValueError("Bootstrap training states must match the fitted reference")
    if (X.shape[1] != mapper.embedding_scaler_.n_features_in_
            or programs.shape[1] != mapper.program_scaler_.n_features_in_):
        raise ValueError("Bootstrap training dimensions must match fitted reference")
    donor_ids = np.unique(donors)
    donor_rows = [np.flatnonzero(donors == donor) for donor in donor_ids]
    rng = np.random.default_rng(seed)
    affinities, agreements = [], []
    skipped = 0
    for _ in range(n_bootstrap):
        sampled = rng.choice(len(donor_ids), size=len(donor_ids), replace=True)
        blocks = [donor_rows[i] for i in sampled]
        rows = np.concatenate(blocks)
        if not np.array_equal(np.unique(labels[rows]), mapper.classes_):
            skipped += 1
            continue
        # Draw identifiers deliberately differ when the same donor is repeated.
        draw_donors = np.concatenate([np.repeat(str(i), len(block))
                                      for i, block in enumerate(blocks)])
        boot = ReferenceMapper(beta=mapper.beta, temperature=mapper.temperature)
        boot.fit(X[rows], programs[rows], labels[rows], draw_donors)
        prediction = boot.predict(X_query, programs_query)
        affinities.append(prediction["affinity"])
        agreements.append(prediction["candidate_analog"] == baseline["candidate_analog"])
    n_valid = len(affinities)
    return {
        "mapping_stability": (np.mean(agreements, axis=0) if n_valid
                              else np.full(len(baseline["candidate_analog"]), np.nan)),
        "affinity_sd": (np.std(affinities, axis=0) if n_valid
                        else np.full_like(baseline["affinity"], np.nan)),
        "n_valid_bootstrap": n_valid,
        "n_skipped_missing_states": skipped,
        "n_requested_bootstrap": int(n_bootstrap),
    }
