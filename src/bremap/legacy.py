"""Portable fit/calibrate/map workflow; analogs are resemblance, not lineage."""

import argparse
from importlib.metadata import version
import json
from pathlib import Path
import platform

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.preprocessing import StandardScaler

from . import __version__
from .core import ReferenceMapper, _identifiers, _matrix
from .transport import transport_predict
from .preprocessing import file_hash, fit_raw_svd, normalize_counts, program_scores

SCHEMA = "bremap-model-v1"


def _write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _versions():
    return {"bremap": __version__, "python": platform.python_version(),
            **{name: version(name) for name in
               ("numpy", "scipy", "pandas", "scikit-learn", "anndata")}}


def _new_directory(path):
    path = Path(path)
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError(f"Output already exists and is not empty: {path}")
    path.mkdir(parents=True, exist_ok=True)
    return path


def _program_config(path):
    config = json.loads(Path(path).read_text())
    anchors = config.get("anchors")
    if not isinstance(anchors, dict) or not anchors:
        raise ValueError("Programs JSON requires a nonempty anchors mapping")
    for name, genes in anchors.items():
        if not isinstance(name, str) or not name.strip() or not isinstance(genes, list) or not genes:
            raise ValueError("Each anchor requires a name and a nonempty gene list")
        anchors[name] = sorted(set(_identifiers(genes, "anchor genes", len(genes)).tolist()))
        anchors[name] = [gene.upper() for gene in anchors[name]]
    minimum = config.get("minimum_genes_per_anchor", 1)
    if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 1:
        raise ValueError("minimum_genes_per_anchor must be a positive integer")
    heldout = config.get("heldout_validation", {})
    if not isinstance(heldout, dict) or any(not isinstance(g, list) for g in heldout.values()):
        raise ValueError("heldout_validation must map names to gene lists")
    return {"anchors": {name: anchors[name] for name in sorted(anchors)},
            "heldout_validation": heldout, "minimum_genes_per_anchor": minimum}


def _read_input(path, donor_column, label_column=None, genes=None):
    data = ad.read_h5ad(path)
    if not data.n_obs or not data.n_vars:
        raise ValueError("H5AD input must contain cells and genes")
    cell_ids = _identifiers(data.obs_names, "cell IDs", data.n_obs)
    gene_ids = np.char.upper(_identifiers(data.var_names, "gene IDs", data.n_vars))
    if len(np.unique(cell_ids)) != len(cell_ids) or len(np.unique(gene_ids)) != len(gene_ids):
        raise ValueError("Cell IDs and case-insensitive gene IDs must be unique")
    if donor_column not in data.obs:
        raise ValueError(f"Missing donor column: {donor_column}")
    donors = _identifiers(data.obs[donor_column], "donors", data.n_obs)
    labels = None
    if label_column:
        if label_column not in data.obs:
            raise ValueError(f"Missing label column: {label_column}")
        labels = _identifiers(data.obs[label_column], "labels", data.n_obs)
    counts = sparse.csr_matrix(data.X, dtype=float)
    # Validate the complete matrix, including extra genes, before alignment.
    if (not np.isfinite(counts.data).all() or (counts.data < 0).any()
            or not np.allclose(counts.data, np.rint(counts.data), rtol=0, atol=1e-5)):
        raise ValueError("H5AD X must contain finite nonnegative raw integer counts")
    if genes is not None:
        positions = pd.Index(gene_ids).get_indexer(genes)
        if (positions < 0).any():
            raise ValueError(f"Missing {int((positions < 0).sum())} reference genes")
        counts = counts[:, positions]
        gene_ids = np.asarray(genes)
    counts.sort_indices()
    zero = np.asarray(counts.sum(axis=1)).ravel() == 0
    return data, normalize_counts(counts), gene_ids, cell_ids, donors, labels, zero


def _embedding(data, key):
    if key not in data.obsm:
        raise ValueError(f"Missing frozen embedding: obsm[{key!r}]")
    provenance = data.uns.get("bremap_embedding_provenance", {}).get(key)
    if not isinstance(provenance, str) or not provenance.strip():
        raise ValueError("Frozen embeddings require a nonempty "
                         "uns['bremap_embedding_provenance'][key] identifier")
    return _matrix(data.obsm[key], "embedding"), provenance


def fit_reference(args):
    data, logcounts, genes, _, donors, labels, zero = _read_input(
        args.reference, args.donor_column, args.label_column)
    if zero.any():
        raise ValueError("Reference contains zero-library cells; remove them before fitting")
    if "unknown" in labels:
        raise ValueError("Reference label 'unknown' is reserved for abstention")
    config = _program_config(args.programs)
    programs, coverage = program_scores(logcounts, genes, config)
    arrays = {}
    if args.embedding_key:
        X, provenance = _embedding(data, args.embedding_key)
        representation = {"kind": "frozen_embedding", "key": args.embedding_key,
                          "provenance": provenance, "dimensions": X.shape[1]}
    else:
        X, _, svd = fit_raw_svd(logcounts, logcounts, components=32, seed=42)
        # Null SVD directions otherwise turn roundoff into scaled biology.
        keep = svd.singular_values_ > (svd.singular_values_[0] *
                                      max(logcounts.shape) * np.finfo(float).eps)
        arrays["svd_components"] = svd.components_[keep]
        X = np.asarray(logcounts @ arrays["svd_components"].T)
        representation = {"kind": "raw_svd", "dimensions": X.shape[1],
                          "normalization": "CP10k_log1p_on_reference_gene_universe", "seed": 42}
    mapper = ReferenceMapper(beta=.25, temperature=.1).fit(X, programs, labels, donors)
    for name in ("embedding", "program"):
        scaler = getattr(mapper, name + "_scaler_")
        for attribute in ("mean", "scale", "var"):
            arrays[f"{name}_{attribute}"] = getattr(scaler, attribute + "_")
    arrays["embedding_centroids"] = mapper.embedding_centroids_
    arrays["program_prototypes"] = mapper.program_prototypes_
    out = _new_directory(args.model)
    np.savez_compressed(out / "arrays.npz", **arrays)
    metadata = {"schema": SCHEMA, "software_version": __version__, "method": args.method,
        "genes": genes.tolist(), "classes": mapper.classes_.tolist(),
        "reference_donors": sorted(set(donors)), "label_column": args.label_column,
        "donor_column": args.donor_column, "programs": config, "anchor_coverage": coverage,
        "representation": representation, "beta": mapper.beta, "temperature": mapper.temperature,
        "transport": {"entropy": .1, "marginal_penalty": .2, "reject_prior": .1},
        "arrays_sha256": file_hash(out / "arrays.npz"),
        "reference_sha256": file_hash(args.reference), "programs_sha256": file_hash(args.programs),
        "versions": _versions(), "status": "provisional_research_software"}
    _write_json(out / "model.json", metadata)


def _load_model(path):
    path = Path(path)
    meta = json.loads((path / "model.json").read_text())
    if meta.get("schema") != SCHEMA or meta.get("method") not in ("affinity", "transport"):
        raise ValueError("Unsupported model schema or method")
    if file_hash(path / "arrays.npz") != meta.get("arrays_sha256"):
        raise ValueError("Model array checksum mismatch")
    with np.load(path / "arrays.npz", allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    if any(a.dtype.kind not in "fiu" or not np.isfinite(a).all() for a in arrays.values()):
        raise ValueError("Model arrays must be finite numeric arrays")
    mapper = ReferenceMapper(beta=meta["beta"], temperature=meta["temperature"])
    mapper.classes_ = _identifiers(meta["classes"], "model classes", len(meta["classes"]))
    mapper.donors_ = _identifiers(meta["reference_donors"], "model donors", len(meta["reference_donors"]))
    if len(mapper.classes_) < 2 or len(set(mapper.classes_)) != len(mapper.classes_):
        raise ValueError("Model requires at least two unique classes")
    dimensions = {"embedding": meta["representation"]["dimensions"],
                  "program": len(meta["programs"]["anchors"])}
    for name, size in dimensions.items():
        scaler = StandardScaler()
        for attr in ("mean", "scale", "var"):
            value = arrays[f"{name}_{attr}"]
            if value.shape != (size,) or (attr == "scale" and (value <= 0).any()):
                raise ValueError("Malformed model scaler dimensions or scales")
            setattr(scaler, attr + "_", value)
        scaler.n_features_in_ = size
        setattr(mapper, name + "_scaler_", scaler)
    for name, size in (("embedding_centroids", dimensions["embedding"]),
                       ("program_prototypes", dimensions["program"])):
        if arrays[name].shape != (len(mapper.classes_), size):
            raise ValueError("Malformed model prototype dimensions")
        setattr(mapper, name + "_", arrays[name])
    if meta["representation"]["kind"] == "raw_svd":
        if arrays["svd_components"].shape != (dimensions["embedding"], len(meta["genes"])):
            raise ValueError("Malformed model SVD dimensions")
    elif meta["representation"]["kind"] != "frozen_embedding":
        raise ValueError("Unknown model representation")
    return meta, arrays, mapper


def _transform(path, meta, arrays, require_labels=False):
    data, logcounts, genes, ids, donors, labels, zero = _read_input(
        path, meta["donor_column"], meta["label_column"] if require_labels else None, meta["genes"])
    programs, _ = program_scores(logcounts, genes, meta["programs"])
    representation = meta["representation"]
    if representation["kind"] == "raw_svd":
        X = np.asarray(logcounts @ arrays["svd_components"].T)
    else:
        X, provenance = _embedding(data, representation["key"])
        if provenance != representation["provenance"]:
            raise ValueError("Frozen embedding provenance does not match the reference")
    if X.shape[1] != representation["dimensions"]:
        raise ValueError("Embedding dimensions do not match the reference")
    return X, programs, ids, donors, labels, zero


def _disjoint(donors, previous):
    overlap = set(donors) & set(previous)
    if overlap:
        raise ValueError(f"Donor overlap between partitions: {list(map(str, sorted(overlap)[:5]))}")


def _predict(meta, mapper, X, programs, donors, reject_cost=None):
    if meta["method"] == "transport":
        return transport_predict(mapper, X, programs, donors, reject_cost, **meta["transport"])
    result = mapper.predict(X, programs)
    result["geometric_divergence"] = result["divergence"].copy()
    return result


def calibrate(args):
    directory = Path(args.model)
    if (directory / "artifact.json").exists():
        raise ValueError("Versioned artifacts require separate calibration through bremap calibrate --out")
    if (directory / "calibration.json").exists():
        raise ValueError("Model is already calibrated; fit a new model for different calibration")
    meta, arrays, mapper = _load_model(directory)
    X, programs, _, donors, labels, zero = _transform(args.calibration, meta, arrays, True)
    _disjoint(donors, meta["reference_donors"])
    if zero.any():
        raise ValueError("Calibration contains zero-library cells; remove them before calibration")
    if not np.isin(labels, mapper.classes_).all():
        raise ValueError("This calibration mode requires known reference-state labels only")
    if set(labels) != set(mapper.classes_):
        raise ValueError("Calibration must represent every reference state")
    geometric = mapper.predict(X, programs)["divergence"]
    reject_cost = float(np.quantile(geometric, .95)) if meta["method"] == "transport" else None
    prediction = _predict(meta, mapper, X, programs, donors, reject_cost)
    threshold = float(np.quantile(prediction["divergence"], .95))
    _write_json(directory / "calibration.json", {
        "schema": "bremap-calibration-v1", "model_sha256": file_hash(directory / "model.json"),
        "input_sha256": file_hash(args.calibration), "donors": sorted(set(donors)),
        "threshold": threshold, "reject_cost": reject_cost, "known_acceptance_target": .95,
        "observed_cell_acceptance": float(np.mean(prediction["divergence"] <= threshold)),
        "n_cells": len(X), "policy": "known_only_cell_quantile_no_unknown_detection_guarantee",
        "versions": _versions()})


def map_query(args):
    directory = Path(args.model)
    if not (directory / "calibration.json").is_file():
        raise ValueError("Calibrate the model before mapping a query")
    meta, arrays, mapper = _load_model(directory)
    calibration = json.loads((directory / "calibration.json").read_text())
    if (calibration.get("schema") != "bremap-calibration-v1" or
            calibration.get("model_sha256") != file_hash(directory / "model.json")):
        raise ValueError("Calibration schema or model checksum mismatch")
    threshold = calibration["threshold"]
    if not np.isscalar(threshold) or not np.isfinite(threshold) or threshold < 0:
        raise ValueError("Invalid calibrated threshold")
    _map_with_model(args, meta, arrays, mapper, calibration)


def _map_with_model(args, meta, arrays, mapper, calibration):
    """Write predictions with supplied validated frozen model and calibration."""
    directory = Path(args.model)
    threshold = calibration["threshold"]
    X, programs, ids, donors, _, zero = _transform(args.query, meta, arrays)
    _disjoint(donors, [*meta["reference_donors"], *calibration["donors"]])
    # A known-invalid row must not alter the coupled solution for valid cells.
    # Retain it in the output with undefined scores and an explicit QC reason.
    prediction = {key: np.full(len(ids), np.nan) for key in
                  ("analog_strength", "divergence", "geometric_divergence", "embedding_divergence")}
    prediction["candidate_analog"] = np.full(len(ids), "unknown", dtype=object)
    prediction["affinity"] = np.full((len(ids), len(mapper.classes_)), np.nan)
    prediction["group_diagnostics"] = []
    if meta["method"] == "transport":
        prediction.update({key: np.full(len(ids), np.nan) for key in ("reject_mass", "matched_mass")})
    valid = ~zero
    if valid.any():
        fitted = _predict(meta, mapper, X[valid], programs[valid], donors[valid], calibration["reject_cost"])
        for key in prediction:
            if key == "group_diagnostics":
                prediction[key] = fitted.get(key, [])
            else:
                prediction[key][valid] = fitted[key]
    rejected = (prediction["divergence"] > threshold) | zero
    frame = pd.DataFrame({"cell_id": ids, "donor_id": donors,
        "candidate_analog": prediction["candidate_analog"],
        "best_analog": np.where(rejected, "unknown", prediction["candidate_analog"]),
        "analog_strength": prediction["analog_strength"], "divergence": prediction["divergence"],
        "geometric_divergence": prediction["geometric_divergence"],
        "embedding_divergence": prediction["embedding_divergence"],
        "abstain_flag": rejected, "qc_zero_library": zero,
        "abstain_reason": np.where(zero, "zero_library", np.where(rejected, "divergence", ""))})
    for key in ("reject_mass", "matched_mass"):
        if key in prediction:
            frame[key] = prediction[key]
    for i, state in enumerate(mapper.classes_):
        frame[f"affinity::{state}"] = prediction["affinity"][:, i]
    summaries = []
    for donor, group in frame.groupby("donor_id", sort=True):
        summaries.append({"donor_id": donor, "n_cells": len(group),
            "coverage": float((~group.abstain_flag).mean()), "n_zero_library": int(group.qc_zero_library.sum()),
            "mean_divergence": float(group.divergence.mean()),
            **{f"mean_affinity::{state}": float(group[f"affinity::{state}"].mean())
               for state in mapper.classes_}})
    out = _new_directory(args.out)
    frame.to_csv(out / "cells.csv", index=False)
    pd.DataFrame(summaries).to_csv(out / "donors.csv", index=False)
    _write_json(out / "reproducibility.json", {
        "schema": "bremap-prediction-v1", "method": meta["method"], "versions": _versions(),
        "input_sha256": file_hash(args.query), "model_sha256": file_hash(directory / "model.json"),
        "arrays_sha256": meta["arrays_sha256"],
        "calibration_sha256": file_hash(directory / "calibration.json"),
        "threshold": threshold, "reject_cost": calibration["reject_cost"],
        "representation": meta["representation"], "donor_column": meta["donor_column"],
        "zero_library_policy": "retain_unknown_with_undefined_scores_exclude_from_coupling",
        "n_valid_cells": int(valid.sum()), "n_zero_library": int(zero.sum()),
        "group_diagnostics": prediction.get("group_diagnostics", []),
        "output_sha256": {name: file_hash(out / name) for name in ("cells.csv", "donors.csv")},
        "code_sha256": {name: file_hash(Path(__file__).with_name(name)) for name in
                        ("legacy.py", "core.py", "transport.py", "preprocessing.py")},
        "limits": "Known-only cell-quantile calibration; no patient CI, bootstrap or lineage inference."})


def toy(args):
    root = _new_directory(args.out)
    genes = ["A1", "A2", "B1", "B2", "N1", "N2"]
    counts = np.array([[100, 80, 2, 1, 8, 10], [90, 75, 1, 2, 10, 9],
                       [1, 2, 100, 80, 10, 8], [2, 1, 90, 75, 9, 10]] * 6)
    for split in ("reference", "calibration", "query"):
        values = counts if split != "query" else np.vstack(
            [counts[:8], np.tile([1, 1, 1, 1, 2000, 2000], (4, 1)), np.zeros(6)])
        labels = np.resize(["a", "a", "b", "b"], len(values))
        if split == "query":
            labels = np.array([*labels[:8], *(["outside"] * 4), "zero"], dtype=str)
        obs = pd.DataFrame({"reference_state": labels,
                            "global_donor_key": [f"{split}_d{i % 2}" for i in range(len(values))]},
                           index=[f"{split}_{i}" for i in range(len(values))])
        ad.AnnData(sparse.csr_matrix(values, dtype=float), obs=obs,
                   var=pd.DataFrame(index=genes)).write_h5ad(root / f"{split}.h5ad")
    _write_json(root / "programs.json", {"anchors": {"a": genes[:2], "b": genes[2:4]},
                                        "minimum_genes_per_anchor": 2})
    report = {"deterministic": True, "outside_queries": 4, "methods": {}}
    for method in ("affinity", "transport"):
        model = root / method / "model"
        main(["fit-reference", "--reference", str(root / "reference.h5ad"), "--model", str(model),
              "--programs", str(root / "programs.json"), "--method", method])
        main(["calibrate", "--model", str(model), "--calibration", str(root / "calibration.h5ad")])
        for name in ("mapped", "repeat"):
            main(["map", "--model", str(model), "--query", str(root / "query.h5ad"),
                  "--out", str(root / method / name)])
        same = all(file_hash(root / method / "mapped" / name) ==
                   file_hash(root / method / "repeat" / name) for name in ("cells.csv", "donors.csv"))
        cells = pd.read_csv(root / method / "mapped/cells.csv")
        rejected = int(cells.iloc[8:12].best_analog.eq("unknown").sum())
        if not same or rejected == 0:
            raise RuntimeError("Toy determinism or outside-query rejection check failed")
        report["methods"][method] = {"deterministic": same, "outside_queries_rejected": rejected}
    _write_json(root / "toy.json", report)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    fit = sub.add_parser("fit-reference", help="Fit reference-only preprocessing and prototypes")
    fit.add_argument("--reference", required=True)
    fit.add_argument("--model", required=True)
    fit.add_argument("--programs", required=True)
    fit.add_argument("--label-column", default="reference_state")
    fit.add_argument("--donor-column", default="global_donor_key")
    fit.add_argument("--embedding-key", help="Frozen obsm embedding; otherwise train-only raw SVD")
    fit.add_argument("--method", choices=("affinity", "transport"), default="affinity")
    fit.set_defaults(func=fit_reference)
    cal = sub.add_parser("calibrate", help="Freeze known-state 95th-percentile rejection settings")
    cal.add_argument("--model", required=True)
    cal.add_argument("--calibration", required=True)
    cal.set_defaults(func=calibrate)
    mapping = sub.add_parser("map", help="Map disjoint query donors using a calibrated model")
    mapping.add_argument("--model", required=True)
    mapping.add_argument("--query", required=True)
    mapping.add_argument("--out", required=True)
    mapping.set_defaults(func=map_query)
    example = sub.add_parser("toy", help="Generate and run a deterministic two-method example")
    example.add_argument("--out", required=True)
    example.set_defaults(func=toy)
    args = parser.parse_args(argv)
    try:
        args.func(args)
    except (ValueError, OSError, KeyError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
