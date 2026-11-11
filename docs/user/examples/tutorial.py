# %% load
from pathlib import Path
import anndata as ad
import pandas as pd
import bremap as br
from bremap.inputs import prepare_input
from bremap.workflow import run_workflow

counts = pd.read_csv("demo_counts.csv", index_col="cell_id")
print(counts.groupby("role").agg(cells=("donor_id", "size"), donors=("donor_id", "nunique")))
# %% spec
spec = br.create_spec(
    spec_id="synthetic-immune-demo", tissue="immune",
    states=["synthetic_quiescent", "synthetic_cytotoxic"],
    programs={"quiescent": ["CCR7", "TCF7", "LEF1"],
              "cytotoxic": ["NKG7", "PRF1", "GZMB"]},
    reference_source="Invented counts for this software tutorial",
    eligibility="Synthetic cells only; not a validated biological signature",
    minimum_genes_per_anchor=2,
)
print(spec.data["states"])
# %% prepare
inputs = {}
for role in ["reference", "calibration", "query"]:
    table = counts.loc[counts.role == role]
    genes = table.drop(columns=["role", "donor_id", "state"])
    obs = table[["donor_id", "state"]].copy()
    if role == "query":
        obs = obs.drop(columns="state")
    data = ad.AnnData(genes.to_numpy(), obs=obs, var=pd.DataFrame(index=genes.columns))
    inputs[role] = prepare_input(data, spec, Path("tutorial_inputs") / f"{role}.h5ad", role=role)
    report = br.inspect_input(inputs[role], spec, role=role)
    print(role, "valid:", report["valid"])
# %% run
result = run_workflow({
    "schema": "bremap-workflow-v1", "name": "first-reference-map",
    "spec": spec, **inputs, "method": "affinity", "policy": "normal_donor_v1",
    "annotations": "demo_truth.csv", "output": "tutorial_run",
})
print("Mapping:", (result / "mapping/cells.csv").is_file())
print("Programs:", (result / "programs/cells.csv").is_file())
print("Reusable reference:", (result / "reference.zip").is_file())
# %% mapping
cells = pd.read_csv(result / "mapping/cells.csv")
print(cells[["cell_id", "candidate_state", "accepted"]].head(4).to_string(index=False))
print(f"Accepted: {cells.accepted.sum()}/{len(cells)} cells")
# %% programs
programs = pd.read_csv(result / "programs/cells.csv")
print(programs[["cell_id", "quiescent", "cytotoxic"]].head(4).round(3).to_string(index=False))
# %% evaluate
summary = pd.read_csv(result / "evaluation/summary.csv")
selected = summary.loc[(summary.aggregation == "donor_macro") &
                       summary.metric.isin(["acceptance", "known_macro_f1"])]
print(selected[["metric", "value"]].round(4).to_string(index=False))
# %% reuse
reference = br.import_reference(result / "reference.zip", "tutorial_restored")
reused = br.map_query(reference, inputs["query"], "tutorial_reused")
pd.testing.assert_frame_equal(cells, pd.read_csv(reused.path / "cells.csv"))
print("Identical mapping after reference reuse: True")
# %% end
