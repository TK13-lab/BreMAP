Troubleshooting
===============

Start with ``br.inspect_input(path, spec, role=...)`` and read its ``errors``.
For file schemas and model options, see :doc:`input-options` and :doc:`distillation`.

.. list-table:: Common problems
   :header-rows: 1
   :widths: 30 70

   * - Problem
     - Action
   * - Invalid counts
     - Load raw integer counts, using ``counts_layer`` if needed. Do not round normalized expression.
   * - Missing metadata or state labels
     - Check column names and the declared state vocabulary. Reference and calibration require known-state labels.
   * - Donor/cell overlap
     - Split biological donors before preparing inputs; check patient/sample identities.
   * - Missing program or student genes
     - Check identifiers and coverage. Mapping uses human gene symbols; Geneformer requires explicit Ensembl IDs.
   * - Embedding/provenance mismatch
     - Use the same encoder and extraction recipe across partitions. Student vectors need a separate student-derived reference.
   * - Checkpoint or upstream package missing
     - Use ``TeacherEncoder.from_pretrained`` with the download extra, or provide a local checkpoint. Install the official inference package separately. See :doc:`distillation`.
   * - Download interrupted or source unavailable
     - Retry the same checkpoint/cache after checking the source or quota; completed files are retained. No alternative checkpoint is substituted.
   * - Checkpoint cache corrupt or unmanaged
     - Inspect the named directory. Remove only that checkpoint to redownload, or use another ``cache_dir``. For existing local models, use ``TeacherEncoder(model, checkpoint)``.
   * - Cache fingerprint differs
     - Inputs or encoder settings changed. Choose a new cache directory for the changed experiment.
   * - Model not calibrated
     - Pass the artifact returned by ``br.calibrate`` or import a calibrated reference.
   * - Output already exists
     - Choose a fresh destination; fitted artifacts are immutable.
   * - Low acceptance
     - Inspect reference suitability, gene coverage and distance/unresolved mass. Affinity alone does not imply acceptance.
   * - Transport changes after batching
     - Keep each donor's query cells together; transport depends on the group.
   * - Evaluation IDs/donors differ
     - Supply labels for exactly the mapped cells with matching donor identifiers.
   * - ``case_pass`` is ``None``
     - Provide held-out task labels; teacher agreement alone cannot establish task quality.

For unresolved errors, open a `GitHub issue
<https://github.com/TK13-lab/BreMAP/issues>`_ with the package/Python versions,
exact error and a small reproducible example.
