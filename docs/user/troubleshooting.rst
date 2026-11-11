Troubleshooting
===============

Start by running ``br.inspect_input(path, spec, role=...)``. Read ``errors`` and
``warnings`` separately: a valid schema is not a biological validation. For a
completed workflow, inspect ``input_inspection.json`` and ``workflow.json`` before
changing the analysis.

Counts and metadata
-------------------

``Raw integer counts in X are required`` or invalid count values
   Restore the original count layer into an AnnData copy's ``X``. Scaled/log data
   cannot be repaired by rounding. Recheck sparsity, missing values and nonnegativity.

Missing metadata column or conflicting metadata column
   Supply the correct ``donor_column``/``label_column`` to ``prepare_input``.
   If the source already has a target-named column with different values, resolve
   that ambiguity explicitly instead of overwriting it silently.

Reference/calibration labels do not match the declared states
   Check spelling and cell selection. Both development partitions must contain
   exactly the known-state vocabulary. Do not invent absent labels merely to make
   validation pass; redesign the eligible split or the declared reference task.

Donor overlap or cell overlap between partitions
   Review the patient/sample crosswalk and split donors before preparing files.
   Renaming the same biological patient in another partition does not produce an
   independent test set.

Missing program genes or overlapping gene sets
   Inspect the input gene namespace, coverage and chosen gene universe. BReMAP
   expects human symbols and does not automatically convert Ensembl IDs or mouse
   orthologues. The minimum covered-gene setting is an eligibility parameter,
   not a guarantee of a reliable biological program.

Frozen embeddings
-----------------

Missing embedding or missing provenance
   Supply the named ``obsm`` matrix and a nonempty provenance string at
   ``uns["bremap_embedding_provenance"][key]``. All rows must align to cell IDs.

Embedding dimensions/provenance disagree
   Recheck checkpoint, preprocessing, extraction universe and row alignment.
   Re-extract incompatible inputs with a common procedure; do not merely change
   the metadata label. The package does not run encoder inference for you.

Mapping and output folders
--------------------------

Calibrate the model before mapping
   Use the object returned by ``br.calibrate``, or import a calibrated archive.
   A fitted reference without calibration is not ready for ``map_query``.

Output already exists
   Individual preparation/artifact outputs are immutable. Use a fresh destination.
   For a complete composed run, ``overwrite=True`` works only for an unmodified
   output previously created by ``run_workflow``. Keep your notes outside that folder.

Overwrite requires a recognized, unmodified workflow output
   Check for added or edited files and changed checksums. The refusal protects
   unrelated work. Do not remove checksum checks to force a replacement.

Low acceptance or strong affinity despite rejection
   Inspect reference distance/unresolved mass as well as conditional affinities.
   Check gene overlap, population eligibility and representation compatibility.
   Do not tune acceptance on a held-out outcome simply to increase coverage.

Transport output changes after splitting a donor across batches
   Keep complete donor groups. Coupled transport depends on the query group, so
   cell batching is not guaranteed to preserve the original solution.

Requesting support
------------------

Record the BReMAP/Python versions, method, specification and exact error. A tiny
synthetic example is preferable to sharing a full patient dataset. Include the
relevant inspection errors and input dimensions, with identifiers appropriate
for the sharing context. Report reproducible problems at
`GitHub Issues <https://github.com/TK13-lab/BreMAP/issues>`_.
