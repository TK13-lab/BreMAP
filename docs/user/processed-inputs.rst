Use a processed AnnData or Seurat object
========================================

You do not need to recover raw counts when your object contains documented
library-normalized expression. BReMAP 0.8 supports this route from loading through
mapping, gene programs, portable references and user-trained students. Other
processed expression has a separate, explicitly declared analysis mode.

1. Inspect the downloaded object
--------------------------------

.. code-block:: python

   import bremap as br

   report = br.inspect_data("processed.h5ad")
   print(report["matrices"])
   print(report["obs_columns"])
   print(report["embeddings"])

The inventory lists matrix names, dimensions, sparse storage, a preview of IDs
and numeric checks. ``metadata_candidates`` suggests column names to inspect;
``operations_by_expression`` explains conditional compatibility. These are not
automatic annotations or validation of the source's processing.

H5AD inspection samples at most 128 rows. ``read_data`` validates the whole selected
matrix. An integer-looking sample or a slot called ``raw.X`` does **not** establish
raw counts. Consult the source paper or processing description.

2. Select the matrix and its actual meaning
-------------------------------------------

For natural-log-transformed library-normalized expression in X:

.. code-block:: python

   data = br.read_data(
       "processed.h5ad", matrix="X",
       expression_type="log1p_library_normalized", log_base="e",
   )
   print(data.shape)
   print(data.uns["bremap_expression"])

``expression_type`` is retained in the object; it survives H5AD save/load. No values
are rounded into counts. For log2(1 + normalized expression), use ``log_base=2``.
For linear library-normalized values use ``expression_type="library_normalized"``.
The declaration means each cell's values are a positive scalar times its counts,
optionally log1p-transformed. It does not cover TPM, imputation, residuals or
batch-corrected expression. Missing genes and numerical precision can affect results.

If counts are present, select them directly:

.. code-block:: python

   data = br.read_data("processed.h5ad", matrix="layers/counts", expression_type="counts")
   # Or matrix="raw.X", only when the source documents raw counts there.

``counts_layer="counts"`` remains a compatible shortcut. Selecting ``raw.X`` keeps
its own gene set, which can differ from X. The original file is unchanged.

For Seurat or SingleCellExperiment RDS, install the optional Python reader from
the BReMAP checkout with ``python -m pip install '.[rds]'``; **R is not required**:

.. code-block:: python

   print(br.inspect_data("processed.rds")["matrices"])
   data = br.read_data(
       "processed.rds", assay="RNA", matrix="data",
       expression_type="log1p_library_normalized",
   )
   # SingleCellExperiment: matrix="logcounts"; omit assay="RNA".

Supported objects are in-memory Seurat Assay/Assay5, SingleCellExperiment and
named numeric sparse/dense gene-by-cell matrices. Cells, genes, metadata and named
reductions are aligned explicitly. Split layers such as ``data.batchA`` require
that exact name; they are not silently merged. A standalone matrix has no donor
annotations: add your verified metadata before partitioning. Imported reductions
do not establish scGPT/Geneformer provenance.

RDS parsing decompresses the whole object and can need substantially more RAM than
the file size. Sparse matrices stay sparse after conversion, but this is not a
streaming atlas reader. Delayed/external-backed classes need an in-memory export
or H5AD from their source workflow. The reader does not execute serialized R code.

3. Bind processing to the reference, then run
---------------------------------------------

Use your actual states, programs and donor-disjoint roles. Assume ``data.obs``
already contains verified ``patient``, ``annotation`` and ``role`` columns:

.. code-block:: python

   spec = br.create_spec(
       spec_id="my-processed-reference", tissue="immune",
       states=["naive", "effector"], programs="programs.gmt",
       reference_source="My accession and selected donors",
       eligibility="My population and inclusion criteria",
       expression=data, representation="svd",
   )
   inputs = br.prepare_partitions(
       data, spec, "prepared", donor_column="patient",
       label_column="annotation", role_column="role",
   )
   result = br.run_analysis(inputs, spec, "analysis", method="affinity")
   print(result.summary())
   print(result.table("donor_programs").head())
   result.plot("programs", x="naive_program", y="effector_program")
   result.plot("coverage")

Use the exact program names in your GMT for x/y. ``role`` must contain reference,
calibration and query with disjoint donors. ``method="transport"`` uses the same
input contract. Both methods export a reusable ``reference.zip``. The package
normalizes count/library-derived values to log1p(CP10k) on the reference gene
universe; it does not claim to recover original UMI totals.

4. If the expression is integrated or scaled
---------------------------------------------

Use a new reference fitted on that explicitly declared scale:

.. code-block:: python

   data = br.read_data(
       "scaled.h5ad", matrix="X", expression_type="preprocessed",
       processing_id="my-fixed-development-scaler-v1", units="gene z score",
   )
   # Continue with create_spec(..., expression=data) and the same APIs above.

Finite negative values are allowed. No CP10k or log transform is applied; program
scores are means in the supplied units, including plot labels. Reference,
calibration and query must declare the same processing identifier and units.
This identifier is your declaration, not proof that independent integrations are
comparable or that preprocessing avoided query leakage. Count-trained references
and students reject this incompatible scale.

5. Optional teacher and student route
-------------------------------------

Declared library-normalized/log1p inputs can use :doc:`distillation`. Pass
``data=data`` to ``TeacherEncoder.from_pretrained`` so incompatible expression is
rejected **before download**. Without that argument, data are checked at transform.
Arbitrary residual/integrated expression cannot generate teacher vectors; use
compatible vectors generated from supported input. Student fitting can use such
paired targets with an explicitly declared processed input contract.

A bounded check used 96 real cells from six HBCA donors with matched counts,
CP10k and natural log1p(CP10k). scGPT whole-human and Geneformer V1-10M were tested
against their actual upstream binning/tokenization, with cached replay and student
save/load checks. This is an interoperability check, not biological validation or
an upstream guarantee for every normalization, checkpoint or gene universe.
