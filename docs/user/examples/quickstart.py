"""A complete installed-BReMAP example; all cells and programs are synthetic."""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

import bremap as br
from bremap.inputs import prepare_input
from bremap.workflow import run_workflow


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    spec = br.create_spec(
        spec_id="synthetic-immune-demo", tissue="immune",
        states=["synthetic_quiescent", "synthetic_cytotoxic"],
        programs={"quiescent": ["CCR7", "TCF7", "LEF1"],
                  "cytotoxic": ["NKG7", "PRF1", "GZMB"]},
        reference_source="Invented counts for this software tutorial",
        eligibility="Synthetic cells only; not a validated biological signature",
        minimum_genes_per_anchor=2,
    )
    genes = ["CCR7", "TCF7", "LEF1", "NKG7", "PRF1", "GZMB", "CD3D", "CD3E"]
    states = spec.data["states"]
    rng = np.random.default_rng(42)
    prepared = {}
    for role in ["reference", "calibration", "query"]:
        counts = np.tile([[90, 80, 70, 2, 1, 2, 30, 25],
                          [2, 1, 2, 90, 80, 70, 30, 25]], (9, 1))
        counts += rng.poisson(2, size=counts.shape)
        obs = pd.DataFrame(
            {"patient": [f"{role}-donor-{i // 6}" for i in range(18)]},
            index=[f"{role}-cell-{i}" for i in range(18)],
        )
        if role != "query":
            obs["annotation"] = states * 9
        data = ad.AnnData(sparse.csr_matrix(counts), obs=obs,
                         var=pd.DataFrame(index=genes))
        prepared[role] = prepare_input(
            data, spec, out / "inputs" / f"{role}.h5ad", role=role,
            donor_column="patient",
            label_column="annotation" if role != "query" else None,
        )
        if role == "query":
            # Truth stays separate from the query input; it is used only to evaluate.
            pd.DataFrame({"cell_id": obs.index, "donor_id": obs.patient.to_numpy(),
                          "label": states * 9}).to_csv(out / "annotations.csv", index=False)

    result = run_workflow({
        "schema": "bremap-workflow-v1", "name": "synthetic-immune-demo",
        "spec": spec, **prepared, "method": "affinity",
        "policy": "normal_donor_v1",  # Historical name; no normal-tissue requirement.
        "annotations": out / "annotations.csv", "output": out / "workflow",
    })
    cells = pd.read_csv(result / "mapping/cells.csv")
    programs = pd.read_csv(result / "programs/cells.csv")
    assert len(cells) == len(programs) == 18
    assert cells.cell_id.is_unique and cells.donor_id.nunique() == 3
    assert {"quiescent", "cytotoxic"} <= set(programs.columns)

    # Reuse the exported calibrated reference through individual package APIs.
    restored = br.import_reference(result / "reference.zip", out / "restored")
    mapped = br.map_query(restored, prepared["query"], out / "reused_mapping")
    pd.testing.assert_frame_equal(cells, pd.read_csv(mapped.path / "cells.csv"))
    print(f"BReMAP {br.__version__}: mapped {len(cells)} synthetic query cells / 3 donors")
    print("Reference reuse: identical mapping; program and evaluation tables written")
    print(f"Results: {result}")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("bremap_quickstart"))
    main(parser.parse_args().out)
