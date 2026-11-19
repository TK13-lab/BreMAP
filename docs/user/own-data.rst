Use your own data
=================

See :doc:`input-options` for expression files, model choices, custom program files and
validation errors.

Start with the question
-----------------------

Choose a reference whose states are meaningful for the query. For example, a
researcher may compare immune activation profiles between patient groups, or
describe how query vascular cells resemble states in a declared vascular atlas.
The reference can be task-defined, disease-associated or drawn from a relevant
baseline population. It need not be normal epithelium.

The current mapper requires at least two discrete, declared reference states.
Gene programs need not correspond one-to-one to those states. Overlapping program
genes are allowed but reported; a program is a positive gene list scored by mean
log1p(CP10k) expression, or supplied units in preprocessed mode. It is not a signed regulatory network or a newly inferred pathway.

Prepare the biological partitions
---------------------------------

Assign **reference**, **calibration** and **query** donors before fitting. Do not
split cells randomly across these roles while retaining the same patients.
Reference and calibration must each contain exactly the declared states. Query
labels are optional, and any evaluation labels should remain separate from the
inputs used to select the configuration.

With paired biopsies, the patient is still the biological unit: samples from the
same patient must not cross development and held-out roles. BReMAP validates
identifiers and recorded roles; it does not discover patient identity across
inconsistent study metadata or design your partition for you.

Input contract
--------------

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - AnnData field
     - Required content
   * - ``X``
     - Counts or explicitly declared processed expression, loaded with read_data; sparse or dense.
   * - ``obs_names``
     - Unique, nonempty cell identifiers, disjoint across roles.
   * - ``var_names``
     - Human gene symbols. Inspect missing/duplicate symbols and harmonize deliberately.
   * - ``obs[donor_column]``
     - Consistent biological donor/patient identifiers.
   * - ``obs[label_column]``
     - Declared state labels for reference and calibration; optional for query.
   * - ``obsm[embedding_key]``
     - Optional cell-aligned frozen vectors; see :doc:`representations` for provenance requirements.

Use ``br.read_data("cells.h5ad", counts_layer="counts")`` to select a count layer.
For normalized-only AnnData/Seurat/SCE, follow :doc:`processed-inputs`. Bind the
returned declaration with ``create_spec(..., expression=data)``. Canonical inputs
use CP10k/log1p on the reference gene universe; as-supplied processed expression
retains its declared units.

Define your specification
-------------------------

.. code-block:: python

   import bremap as br

   spec = br.create_spec(
       spec_id="my-immune-reference", tissue="immune",
       states=["naive", "effector"],
       programs={
           "naive_program": ["CCR7", "TCF7", "LEF1"],
           "cytotoxic_program": ["NKG7", "PRF1", "GZMB"],
       },
       reference_source="My annotated cohort and its accession",
       eligibility="My documented donor and cell inclusion criteria",
       donor_column="donor_id", label_column="state",
       minimum_genes_per_anchor=2,
   )

These lists demonstrate syntax, not a validated biological signature. Replace
the states, programs and eligibility with the choices for your task. A GMT or
supported JSON program file can also be passed as ``programs``;
:func:`bremap.read_gene_programs` reads either format. GMT rows contain program
name, description and tab-separated genes.

Prepare three files
-------------------

.. code-block:: python

   from pathlib import Path
   from bremap.inputs import prepare_input

   inputs = Path("prepared")
   reference = prepare_input(
       "reference_counts.h5ad", spec, inputs / "reference.h5ad",
       role="reference", donor_column="patient", label_column="annotation",
   )
   calibration = prepare_input(
       "calibration_counts.h5ad", spec, inputs / "calibration.h5ad",
       role="calibration", donor_column="patient", label_column="annotation",
   )
   query = prepare_input(
       "query_counts.h5ad", spec, inputs / "query.h5ad",
       role="query", donor_column="patient",
   )
   report = br.inspect_input(query, spec, role="query")
   print(report["valid"], report["errors"], report["warnings"])

An AnnData object can replace each source filename. ``cell_ids`` can select an
explicit list of cells. Preparation copies the source, aliases the specified
metadata and records specification identity; it does not modify the source,
infer cell types or choose a donor split. Existing destinations are refused.

Inspect missing program genes before analysis. Passing schema validation means
the data meet the software contract; it does not certify reference suitability.
Proceed to :doc:`reuse` to fit and map, or :doc:`representations` to supply
precomputed embeddings.
