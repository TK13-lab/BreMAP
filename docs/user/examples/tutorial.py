# %% load
import bremap as br

data = br.load_example()
print(br.describe_data(data))
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
inputs = br.prepare_partitions(data, spec, "tutorial_inputs")
print(br.inspect_input(inputs["query"], spec, role="query")["valid"])
# %% run
truth = br.example_annotations("tutorial_truth.csv")
result = br.run_analysis(inputs, spec, "tutorial_run", annotations=truth)
print(result.summary().round(4))
# %% mapping
print(result.table("mapping")[["cell_id", "candidate_state", "accepted"]].head(4).to_string(index=False))
result.plot("coverage", save="tutorial_coverage.png")
# %% programs
print(result.table("programs")[["cell_id", "quiescent", "cytotoxic"]].head(4).round(3).to_string(index=False))
result.plot("programs", x="quiescent", y="cytotoxic", save="tutorial_programs.png")
# %% affinity
result.plot("affinity", save="tutorial_affinity.png")
# %% evaluate
result.plot("confusion", annotations=truth, save="tutorial_confusion.png")
print(result.table("evaluation").query(
    "aggregation == 'donor_macro' and metric in ['acceptance', 'known_macro_f1']"
)[["metric", "value"]].round(4).to_string(index=False))
# %% reuse
reference = br.import_reference(result.path / "reference.zip", "tutorial_restored")
reused = br.map_query(reference, inputs["query"], "tutorial_reused")
print("Reusable reference mapped successfully:", reused.path.is_dir())
# %% end
