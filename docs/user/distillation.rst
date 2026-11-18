Experimental teacher and student encoders
=========================================

BReMAP 0.6.0 includes optional teacher adapters and an **experimental user-trained
student API**. This is not a validated general-use student checkpoint. A student learns
from the paired counts and teacher vectors you provide. It is not a compressed
copy of every cell used to pretrain the foundation model.

The intended use is a **student tailored to your data and task**. One successful
within-domain application can justify that use; performance across unrelated
tissues is not required. Choose acceptable task error and coverage before the
final test, and evaluate on donors not used to train or select your recipe.
Only claim cross-cohort transfer when tested. You do not need to retrain the
underlying foundation model.

Choose the route
----------------

* Start with count SVD when you do not need teacher embeddings.
* Use ``TeacherEncoder`` to generate vectors for your own dataset with an explicit
  local scGPT or Geneformer checkpoint. No checkpoint is downloaded automatically.
* Use ``StudentEncoder`` when you have paired counts/teacher vectors and separate
  validation donors. Train once, save, then predict on CPU without torch.
* Evaluate on untouched donors from the intended domain. Use another cohort when
  claiming cross-cohort transfer, rather than for every scoped application.

1. Load counts and generate teacher vectors
-------------------------------------------

Install BReMAP from the source checkout with ``pip install '.[train]'`` for student
training. Install the official scGPT or Geneformer package in the environment
used for that teacher. Their dependencies are optional and are not installed by
the base BReMAP package. This adapter was exercised on Linux/WSL with CUDA;
Geneformer currently requires CUDA, and scGPT uses an upstream Linux affinity API.

**Get the checkpoint once.** Model files are not included in BReMAP. For scGPT,
download a complete checkpoint folder from the official `model zoo
<https://github.com/bowang-lab/scGPT#pretrained-scgpt-model-zoo>`_. It contains
``args.json`` (configuration), ``vocab.json`` (gene vocabulary) and
``best_model.pt`` (trained weights). Use these files from the same checkpoint;
do not create them yourself. Set ``checkpoint`` below to that folder's actual
location. Geneformer uses its own `official checkpoints
<https://huggingface.co/ctheodoris/Geneformer>`_ and tokenizer assets instead.
Count-SVD mapping and loading a saved student require no teacher checkpoint.

.. code-block:: python

   import bremap as br

   train = br.read_data("train_counts.h5ad", counts_layer="counts")
   validation = br.read_data("validation_counts.h5ad", counts_layer="counts")

   teacher = br.TeacherEncoder(
       model="scgpt",
       checkpoint="models/scgpt/whole-human",
       device="cuda", batch_size=32, chunk_size=256,
       max_length=1200, gene_case="exact", seed=42,
   )
   train_vectors = teacher.transform(train, cache_dir="cache/train_teacher")
   validation_vectors = teacher.transform(
       validation, cache_dir="cache/validation_teacher"
   )

Omit ``counts_layer`` when ``X`` already contains raw integer counts. Normalized
expression and zero-library cells fail explicitly. These are local files, not
uploads to an external service. Donor metadata must be present in ``obs``.

The return value is a cell-indexed DataFrame, not a positional array. Completed
chunks are reused after interruption. Changed counts, cell IDs, gene identifiers,
checkpoint files, configuration or adapter code invalidate the cache; choose a
new cache directory deliberately. Do not concurrently write to one cache.
``batch_size`` controls GPU batch memory; ``chunk_size`` controls restart granularity.
Geneformer currently reloads its upstream extractor per chunk, so very small
chunks add startup cost. No unmeasured speedup is promised.

Geneformer uses the same API with explicit identifiers and version:

.. code-block:: python

   teacher = br.TeacherEncoder(
       model="geneformer", checkpoint="models/Geneformer-V2-316M",
       model_version="V2", gene_column="ensembl_id",
       device="cuda", batch_size=4, chunk_size=256, seed=42,
   )

``gene_column`` names a column in ``data.var`` containing unversioned human
Ensembl IDs. This lets the count matrix retain gene-symbol names for mapping.
No implicit symbol-to-Ensembl conversion is performed. The Geneformer adapter
uses mean cell pooling at layer -1; V1 uses length 2048 and V2 defaults to 4096.
The scGPT adapter uses CLS pooling and its official binned-count preprocessing.
These are distinct representations, not interchangeable coordinates.

2. Train your own student
-------------------------

.. code-block:: python

   student = br.StudentEncoder.fit(
       train, train_vectors,
       validation=validation, validation_teacher=validation_vectors,
       donor_key="donor_id",
       teacher_provenance=train_vectors.attrs["bremap_encoder_metadata"],
       n_genes=4096, n_components=128, hidden_size=128,
       loss="mse", epochs=40, seed=42,
       gene_case="exact", minimum_gene_fraction=0.90,
   )
   student.save("my_student.npz")

Both vector tables must come from the **same teacher configuration**. The API
joins them by exact cell IDs and checks dimensions/column order. Training and
validation donors/cells must be disjoint. Namespace donor IDs across studies
only after checking whether the studies contain the same people.

Gene selection, SVD projection and target normalization are fitted on training
cells only. Validation selects the epoch by raw embedding MSE, weighted equally
by donor by default. ``loss="mse_relational"`` adds a pairwise-distance term;
``relational_weight`` controls its weight. This option did not improve the local
ablation and is not an automatic recommendation. The selected model's history,
donors, input hashes, parameters and teacher declaration are in ``student.metadata``.

Saving refuses to overwrite an existing model. The NPZ contains arrays and JSON,
not executable pickle objects; loading verifies its integrity fingerprint.

3. Predict and use the vectors in BReMAP
----------------------------------------

.. code-block:: python

   student = br.StudentEncoder.load("my_student.npz")
   query = br.read_data("query_counts.h5ad")
   vectors = student.transform(query)
   print(vectors.shape)

   # Attach counts, vectors and provenance in one package API call.
   query_with_vectors = student.attach(query, key="X_student")

``attach`` returns a copy and refuses to overwrite an existing embedding key.
It registers a **custom** encoder fingerprint; student vectors are never labeled
as authentic scGPT/Geneformer outputs. A missing-gene fraction above the declared
limit raises an error. Below that limit missing features are explicitly warned
about and zero-filled; adequate coverage alone does not establish transferability.
Counts are normalized by all supplied measured genes before selected features
are projected. Differences in assay coverage can therefore affect predictions.

Use ``student.attach`` on **all** reference, calibration and query partitions,
then create a specification with ``representation="custom"`` and
``embedding_key="X_student"``. Follow :doc:`own-data` to fit and calibrate that
reference. Never put student query vectors into a teacher-fitted reference.
Alternatively, ``teacher.attach(data, cache_dir=..., key="X_teacher")`` registers
genuine teacher vectors through the same input workflow.

4. Check the task, not just vector similarity
---------------------------------------------

For annotation mapping, compare cell-indexed mapping tables from separately
fitted teacher and student references:

.. code-block:: python

   report = br.evaluate_student(
       teacher_mapping, student_mapping,
       states=declared_states, truth=heldout_state_labels,
       agreement_threshold=0.85,
       maximum_coverage_change=0.05,
       maximum_f1_loss=0.05,
       minimum_task_f1=0.75, minimum_coverage=0.80,
   )
   print(report["agreement"], report["student_f1"], report["case_pass"])

Tables require ``donor_id``, ``candidate_state`` and boolean ``accepted`` columns.
Labels are a Series indexed by the same cells. You may instead pass CSV paths
for both mappings (``cell_id, donor_id, candidate_state, accepted``) and truth
(``cell_id, donor_id, label``). The API checks exact IDs and matching donors.
Metrics are averaged over donors. F1 includes every declared state, with zero
for an absent state, and all query cells, including rejected cells.
Without independent task labels, ``joint_pass`` is ``None``; agreement alone
cannot establish task fidelity. ``case_pass`` also remains ``None`` without truth.
``joint_pass`` checks relative fidelity; ``case_pass`` additionally checks any
explicit ``minimum_task_f1`` and ``minimum_coverage``. The defaults for these two
absolute floors are ``None``: there is no universal task threshold. Inspect
``thresholds`` and ``donors`` rather than treating a boolean as scientific proof.
Engineering thresholds are configurable, not statistical noninferiority tests.
For a continuous endpoint, use held-out regression error and correlation instead
of applying an annotation-agreement threshold.

What was actually observed locally
----------------------------------

The 22-fit scGPT experiment compared equal-size Wu-only and Wu+Kim training,
compression and loss. The selected mixed student used 4096 genes, SVD128,
hidden128 and MSE. It passed all three mapping checks in 5/8 Kim case/method
combinations at both 85% and 90% agreement. Kim was development data.

On 1548 CD4 cells from eight donors, excluded from student selection, donor-mean
protein-prediction MAE was **0.876** for the selected student, **0.886** for the
historical student, **0.735** for scGPT and **0.732** for count SVD. The student
remained worse than teacher/SVD. This cohort had been used in earlier BReMAP work;
it was new to student development, not a newly collected dataset.

Real API checks produced scGPT vectors for all 1548 cells and Geneformer vectors
for 16 cells with 1152 dimensions. The latter is a smoke test, **not** a
Geneformer multi-tissue efficacy experiment. A separate 240-cell train / 80-cell
validation smoke verified training and portable inference on Geneformer targets.
These checks justify an experimental user-training API, not a public universal
student checkpoint or foundation-model superiority claim.


A scoped case with a frozen follow-up
--------------------------------------

A new LN T-cell student was trained from 1970 cells and validated on 864 cells
from separate donors. Its first five-donor LN test preserved teacher mappings
(94.4% agreement) but missed the absolute F1 floor (0.7375 versus 0.75).
A single follow-up was declared before inspecting four new biopsy patients:
1435 cells, 92.0% agreement, student F1 0.8170 versus teacher 0.8488 and SVD
0.7970, and 97.96% coverage. The frozen student/reference passed the declared
point-estimate checks. This small, same-source-cohort test does not establish
statistical noninferiority or cross-cohort validity. See :doc:`scoped-student`
for the executable recipe, measured costs and complete interpretation.
