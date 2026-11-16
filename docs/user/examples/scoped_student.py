# Run in a directory containing the three raw-count H5AD inputs and local model path.
# Package-only code mirrored from scoped-student.rst; no research script imports.
import bremap as br

reference = br.read_data("reference.h5ad")
validation = br.read_data("calibration.h5ad")
query = br.read_data("query.h5ad")

teacher = br.TeacherEncoder(
    "scgpt", "models/scgpt/whole-human", device="cuda",
    gene_case="insensitive", max_length=1200,
    batch_size=32, chunk_size=256, seed=20261116,
)
training_targets = teacher.transform(reference, cache_dir="teacher/reference")
validation_targets = teacher.transform(validation, cache_dir="teacher/calibration")

student = br.StudentEncoder.fit(
    reference, training_targets,
    validation=validation, validation_teacher=validation_targets,
    teacher_provenance=training_targets.attrs["bremap_encoder_metadata"],
    n_genes=4096, n_components=128, hidden_size=128,
    loss="mse", epochs=40, batch_size=128, seed=20261116,
)
student.save("lymphnode_student.npz")

student = br.StudentEncoder.load("lymphnode_student.npz")
reference_s = student.attach(reference, key="X_student")
validation_s = student.attach(validation, key="X_student")
query_s = student.attach(query, key="X_student")

spec = br.create_spec(
    spec_id="my-lymphnode-student", tissue="lymph_node_T_cells",
    states=["Cytotoxic CD8+ T", "Exhausted CD8+ T", "Naive CD4+ T", "Treg"],
    programs={
        "naive_memory": ["CCR7", "TCF7", "LEF1", "IL7R", "LTB", "MAL"],
        "cytotoxic": ["NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "CTSW"],
        "checkpoint": ["PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "TOX"],
        "cycling": ["MKI67", "TOP2A", "UBE2C", "BIRC5"],
    },
    reference_source="Annotated lymph-node reference donors",
    eligibility="Four declared T-cell states; disjoint patient partitions",
    representation="custom", embedding_key="X_student",
    cost_mode="embedding", minimum_genes_per_anchor=2,
)
inputs = {
    "reference": br.prepare_input(reference_s, spec, "inputs/reference.h5ad", role="reference"),
    "calibration": br.prepare_input(validation_s, spec, "inputs/calibration.h5ad", role="calibration"),
    "query": br.prepare_input(query_s, spec, "inputs/query.h5ad", role="query"),
}
result = br.run_analysis(inputs, spec, "student_mapping", method="affinity")
print(result.summary())
print(result.table("mapping").head())
result.plot("coverage", save="student_coverage.png")
result.plot("affinity", group_by="donor_id", save="student_affinity.png")
result.plot("programs", x="naive_memory", y="cytotoxic", save="student_programs.png")
