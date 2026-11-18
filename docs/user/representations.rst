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
The base workflow accepts precomputed vectors. The optional
:doc:`distillation` API can also extract them using explicitly configured local
scGPT/Geneformer checkpoints, and train a dataset-specific student. It does not
download checkpoints or fine-tune a foundation model.

For every input partition, provide:

.. code-block:: python

   import bremap as br

   data = br.attach_embeddings(
       data, "cell_embeddings.csv", key="X_encoder",
       model="scgpt", checkpoint="my-checkpoint-identifier",
       preprocessing="my-fixed-extraction-recipe",
   )

The CSV must contain unique ``cell_id`` values matching the count input exactly,
followed by numeric embedding columns. Rows are joined by ID, not position.
Create the specification with ``representation="scgpt"`` and
``embedding_key="X_encoder"``. Use identical
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

``cost_mode="biology"`` combines representation and program distances;
``cost_mode="embedding"`` uses representation distances alone. Both report RNA
program scores. ``run_analysis`` handles the calibration policy automatically;
provide calibration donors separate from reference and query donors.
