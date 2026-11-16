Choose input files, gene programs and representations
======================================================================

BReMAP receives local paths or Python data objects. It does not upload your data
to a server. These options are included in BReMAP 0.6.0. Reading a file never certifies its biological suitability.

1. Load counts explicitly
------------------------------------------------------------

.. code-block:: python

   import bremap as br

   data = br.read_data("cells.h5ad", counts_layer="counts")
   # Alternative: CSV/TSV rows are cells, other numeric columns are genes.
   data = br.read_data(
       "counts.csv", cell_id="barcode",
       metadata_columns=["patient", "role", "annotation"],
   )

Choose one loader call. H5AD uses X if counts_layer is omitted; a named layer must
exist. CSV/TSV requires explicit metadata_columns so patient numbers cannot be
silently treated as gene counts. Counts must be finite, nonnegative integers;
normalized data must not be repaired by rounding. Unsupported formats, absent
files/layers, duplicate cell IDs and invalid values produce errors.

2. Supply your programs
------------------------------------------------------------

.. list-table:: Supported program inputs
   :header-rows: 1

   * - Input
     - Schema
   * - Python dict / JSON
     - Program name mapped to a list of human gene symbols.
   * - GMT
     - Tab-separated program name, description, then gene symbols.
   * - CSV / TSV
     - Exactly two columns: program and gene; one row per gene.
   * - TXT
     - One symbol per line; program_name is required.

.. code-block:: python

   programs = br.read_gene_programs("programs.gmt")
   # Alternative: a single custom program.
   programs = br.read_gene_programs("my_genes.txt", program_name="my_program")

CSV example:

.. code-block:: text

   program,gene
   naive,CCR7
   naive,TCF7
   cytotoxic,NKG7
   cytotoxic,PRF1

Duplicate genes within a program, empty lists, malformed JSON/GMT and unsupported
identifiers are rejected with the file/program name. Shared genes across different
programs are allowed and reported as a warning. The loader validates syntax; it
is not a current gene-database lookup and does not convert Ensembl IDs.

3. Choose counts or precomputed model features
------------------------------------------------------------

.. list-table:: create_spec representation options
   :header-rows: 1

   * - representation
     - Required inputs
   * - svd
     - Raw counts; omit embedding_key. SVD is fitted on the reference.
   * - scgpt
     - Precomputed scGPT vectors, embedding_key, declared checkpoint/preprocessing.
   * - geneformer
     - Precomputed Geneformer vectors, embedding_key, declared checkpoint/preprocessing.
   * - custom
     - Other precomputed features with the same provenance and alignment requirements.

Choosing a model does **not** download weights or run that encoder. The package
accepts frozen vectors; raw counts remain necessary for gene programs. Encoder
identity is recorded from your declaration, not inferred from numerical vectors.

.. code-block:: python

   data = br.attach_embeddings(
       data, "scgpt_embeddings.csv", key="X_scgpt", model="scgpt",
       checkpoint="my-exact-checkpoint-id",
       preprocessing="my-gene-universe-and-extraction-version",
       cell_id="barcode",
   )

The embedding table contains the cell ID plus numeric dimensions, with exactly
the same cells as the AnnData. Rows are aligned by IDs; missing/extra/duplicate
cells are rejected. NaN/Inf and an existing obsm key are rejected. CSV, TSV and
DataFrame are accepted. The original AnnData is unchanged; the return value is
a copy containing vectors and provenance. Embeddings already in H5AD may be used,
but named model options require the equivalent encoder metadata/provenance.
Use the same checkpoint, preprocessing and dimension ordering across partitions.

4. Bind the options to a reference specification
------------------------------------------------------------

.. code-block:: python

   spec = br.create_spec(
       spec_id="my_reference", tissue="immune",
       states=["naive", "effector"],
       programs="programs.gmt",
       reference_source="My source accession and selection criteria",
       eligibility="My eligible population and donor split",
       representation="scgpt", embedding_key="X_scgpt",
       minimum_genes_per_anchor=2,
   )
   inputs = br.prepare_partitions(
       data, spec, "prepared",
       donor_column="patient", label_column="annotation", role_column="role",
   )
   result = br.run_analysis(inputs, spec, "analysis", method="affinity")

For a TXT program file, pass programs="genes.txt" and program_name="my_program".
For SVD replace the representation choice with representation="svd" and omit
embedding_key. Calibration/reference labels must match the chosen states; query
labels are removed. The examples' states/programs are placeholders to replace,
not automatic annotations. Each program must meet minimum_genes_per_anchor.

Existing calls with only embedding_key retain their legacy generic frozen-vector
behavior. New explicit encoder declarations additionally check model metadata.
When applying an existing archive, its programs, representation and calibrated
settings are frozen: do not override them at map time. Build a new reference to
change those options and evaluate that configuration separately.

5. Understand errors and warnings
------------------------------------------------------------

Examples of errors:

.. code-block:: text

   scgpt requires embedding_key for precomputed embeddings. Use attach_embeddings or an experimental encoder's attach method before creating the specification.
   Embedding cell IDs must exactly match input: 3 missing, 2 extra. Re-export vectors for the same cells.
   Program file my_genes.txt: TXT gene lists require program_name, e.g. program_name="my_program"

Errors stop preparation/fitting; warnings describe issues requiring interpretation
(for example partially missing program genes above the minimum or shared genes).
inspect_input returns valid, errors, warnings and per-program gene coverage.
Reference, calibration and query embeddings must also agree in dimensions and
provenance during fitting/calibration/mapping. Do not relabel provenance to bypass
these checks; re-extract compatible vectors instead.
