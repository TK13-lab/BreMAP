Representations and mapping methods
===================================

Choose a representation
-----------------------

**Count-derived SVD:** omit ``embedding_key`` when creating the specification.
The representation is fitted on the reference and then applied to calibration
and query inputs using the stored normalization/gene basis. This is a CPU workflow
and a useful starting point before adding an encoder dependency.

**Supplied frozen embeddings:** set ``embedding_key`` to an existing AnnData slot,
for example ``X_scgpt``, ``X_geneformer`` or a user-defined name. BReMAP accepts
numeric vectors; it does not verify biological quality from the slot name.
The package does not download model weights, tokenize cells, extract embeddings
or fine-tune an encoder.

For every input partition, provide:

.. code-block:: python

   # vectors must already be aligned to data.obs_names.
   data.obsm["X_encoder"] = vectors
   data.uns["bremap_embedding_provenance"] = {
       "X_encoder": "encoder/checkpoint/preprocessing/gene-universe identifier"
   }

Create the specification with ``embedding_key="X_encoder"``. Use identical
provenance and dimensions across reference, calibration and query. Keep the
checkpoint, feature universe and extraction recipe in your own provenance record.
Copying the same metadata string onto incompatible vectors does not make them
compatible. Raw counts remain required for RNA-program scores and input checks.

If an endpoint claims genes were held out from the representation, those genes
must actually have been excluded during extraction. The embedding slot cannot
establish that claim retrospectively.

Choose a mapping method
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 25 40 35

   * - Method
     - Interpretation
     - Operational consideration
   * - ``affinity``
     - Reference-state resemblance with a calibrated distance-based acceptance rule.
     - Straightforward starting point; inspect both affinity and acceptance.
   * - ``transport``
     - Coupled assignment with unresolved mass and a calibrated rejection rule.
     - Keep complete donor groups: dividing a donor across runs changes the coupled problem.

Both methods return conditional state affinities. Neither makes these values
validated state proportions, malignancy probabilities or lineage ancestry.
Use a task-appropriate comparator to determine whether the representation helps.

Programs and calibration
------------------------

``cost_mode="biology"`` uses the existing combined geometry/program procedure;
``cost_mode="embedding"`` disables the program contribution to the mapping cost.
Program reporting remains distinct from mapping coordinates. The current
calibration policy selects from predefined settings; changing the number of
programs does not license tuning a configuration on the held-out endpoint.

Pass ``policy="normal_donor_v1"`` explicitly when fitting. This is a historical
identifier for the implemented donor-balanced calibration policy, **not** a rule
that the reference must contain normal cells. Some artifact strings retain older
wording. Calibration donors are separate from reference and query donors, and
the reported acceptance is not a formal coverage guarantee.

For a scientific comparison, retain the tested configurations, declare which
criterion selects one, and evaluate it on samples not used for that selection.
The :doc:`cases` show why foundation-model and mapped-coordinate performance
should be assessed per task rather than assumed.
