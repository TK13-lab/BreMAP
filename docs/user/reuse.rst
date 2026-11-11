Compose and reuse a workflow
============================

Individual Python APIs
----------------------

After preparing ``reference``, ``calibration`` and ``query`` as in :doc:`own-data`,
call the public functions in your notebook or analysis program:

.. code-block:: python

   import bremap as br

   model = br.fit_reference(
       reference, spec, "analysis/reference",
       method="affinity", policy="normal_donor_v1",
   )
   calibrated = br.calibrate(model, calibration, "analysis/calibrated")
   mapping = br.map_query(calibrated, query, "analysis/mapping")
   evaluation = br.evaluate(mapping, "query_annotations.csv", "analysis/evaluation")
   print(mapping.path)
   print(evaluation.tables.keys())

Fit/calibrate/map return artifact objects with ``path`` and ``metadata``. Use
their saved CSV files for table access. Evaluation/reliability additionally return
``tables``. Fitting alone is insufficient: calibrate before mapping.

One composed call
-----------------

``run_workflow`` combines those operations, input inspection, program tables,
an archive and a reproducibility record. The dictionary can contain the actual
``ReferenceSpec`` and ``Path`` objects:

.. code-block:: python

   from pathlib import Path
   from bremap.workflow import run_workflow

   result_dir = run_workflow({
       "schema": "bremap-workflow-v1", "name": "my-analysis",
       "spec": spec, "reference": reference, "calibration": calibration,
       "query": query, "method": "affinity", "policy": "normal_donor_v1",
       "output": Path("complete_analysis"),
       # Optional:
       # "annotations": Path("query_annotations.csv"),
       # "reliability": {"n_bootstrap": 20, "seed": 42},
   })

The return value is a **Path**, not an AnnData object or an ``AnalysisResult``.
Relative paths in a Python dictionary resolve against the current working
directory, or an explicit ``base_dir``. A JSON recipe is optional; its relative
paths instead resolve against the recipe file's directory.

Share a calibrated reference
----------------------------

.. code-block:: python

   archive = br.export_reference(calibrated, "analysis/reference.zip")
   restored = br.import_reference(archive, "shared_reference")
   new_mapping = br.map_query(restored, new_prepared_query, "new_mapping")

The archive contains numeric fitted arrays, settings, specification, calibration
and checksums. It also contains reference/calibration identifiers. It contains
neither original count matrices nor scGPT/Geneformer weights. It transfers the
fitted BReMAP reference, not an encoder or a biological certification of the
user's gene programs.

The composed workflow already writes ``reference.zip``. Reuse that archive
without supplying fitting options:

.. code-block:: python

   reused_dir = run_workflow({
       "schema": "bremap-workflow-v1", "name": "new-query",
       "archive": Path("complete_analysis/reference.zip"),
       "query": Path("new_prepared_query.h5ad"),
       "output": Path("new_query_analysis"),
   })

Prepare the new query using the same specification, including any embedding
provenance. The quickstart verifies exact mapping parity after archive reuse.

Export cell representations
---------------------------

.. code-block:: python

   br.export_embeddings(
       "query_counts.h5ad", "shared_vectors.h5ad",
       key="X_encoder", obs_columns=["donor_id"],
   )

This writes vectors into the exported H5AD's ``X`` with cell IDs, selected
metadata and extraction provenance. These dimensions are latent coordinates,
not named biological pathways. The file is count-free and **cannot substitute
for the count input of another mapping run**. The export does not train or extract
new features; it shares vectors already present in the source.

Reruns and output ownership
---------------------------

Individual artifact destinations must be new. Do not edit their internal files:
checksums are part of the reuse contract. ``run_workflow(..., overwrite=True)``
can replace its own recognized, unmodified prior result after the replacement
succeeds. It refuses unrelated or modified output directories. Keep input files
outside the workflow output directory, and keep personal notes outside the
checksummed bundle.
