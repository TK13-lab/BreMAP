Python API
==========

Expand a group to inspect signatures and parameters. Start with :doc:`quickstart`;
see :doc:`input-options` for file schemas and :doc:`distillation` for encoders.

.. raw:: html

   <details><summary><strong>Inputs and specifications</strong></summary>

.. autofunction:: bremap.usability.read_data
.. autofunction:: bremap.usability.load_example
.. autofunction:: bremap.usability.example_annotations
.. autofunction:: bremap.usability.describe_data
.. autofunction:: bremap.create_spec
.. autofunction:: bremap.read_gene_programs
.. autoclass:: bremap.ReferenceSpec
   :members: from_dict, load, validate_identity
.. autofunction:: bremap.usability.attach_embeddings
.. autofunction:: bremap.inputs.prepare_input
.. autofunction:: bremap.usability.prepare_partitions
.. autofunction:: bremap.inspect_input

.. raw:: html

   </details>
   <details><summary><strong>Run an analysis and inspect results</strong></summary>

.. autofunction:: bremap.usability.run_analysis
.. autoclass:: bremap.usability.WorkflowResult
   :members: table, summary, plot
.. autofunction:: bremap.usability.apply_reference
.. autofunction:: bremap.usability.compare_measurement

.. raw:: html

   </details>
   <details><summary><strong>Teacher and experimental student encoders</strong></summary>

.. autoclass:: bremap.TeacherEncoder
   :members: transform, attach
.. autoclass:: bremap.StudentEncoder
   :members: fit, transform, attach, save, load, fingerprint
.. autofunction:: bremap.evaluate_student

.. raw:: html

   </details>
   <details><summary><strong>Individual mapping steps and evaluation</strong></summary>

.. autofunction:: bremap.fit_reference

   For this lower-level call, pass ``policy="normal_donor_v1"`` explicitly.
   ``run_analysis`` supplies it automatically.

.. autofunction:: bremap.calibrate
.. autofunction:: bremap.map_query
.. autofunction:: bremap.evaluate
.. autofunction:: bremap.assess_reliability
.. autofunction:: bremap.workflow.run_workflow
.. autoclass:: bremap.Artifact
.. autoclass:: bremap.AnalysisResult

.. raw:: html

   </details>
   <details><summary><strong>Export and reuse</strong></summary>

.. autofunction:: bremap.export_reference
.. autofunction:: bremap.import_reference
.. autofunction:: bremap.export_embeddings

.. raw:: html

   </details>
   <details><summary><strong>Optional plotting</strong></summary>

.. autofunction:: bremap.plotting.plot_result
.. autofunction:: bremap.plot_donor_correlations

.. raw:: html

   </details>
