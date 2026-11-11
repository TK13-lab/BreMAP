Python API reference
====================

Signatures below are imported from the installed BReMAP version used to build
this documentation. The primary workflow is specification → preparation → fit →
calibrate → map. The convenience orchestrator returns an output path; individual
mapping functions return artifact objects.

Specifications and gene programs
--------------------------------

.. autofunction:: bremap.create_spec

   Returns a ``ReferenceSpec``. ``programs`` accepts a dictionary or a supported
   GMT/JSON file. ``embedding_key=None`` chooses reference-fitted count SVD.
   ``cost_mode`` is ``biology`` or ``embedding``. At least two states are required.

.. autofunction:: bremap.read_gene_programs

   Returns the validated program-name to gene-list mapping.

.. autoclass:: bremap.ReferenceSpec
   :members: from_dict, load, validate_identity
   :undoc-members:

   ``spec.data`` contains the serializable specification. Identity validation
   checks declared metadata; it does not infer biological eligibility.

Preparation and inspection
--------------------------

.. autofunction:: bremap.inputs.prepare_input

   Returns the prepared H5AD ``Path``. ``source`` is an AnnData object or H5AD path;
   ``cell_ids`` optionally selects cells. Existing destinations are refused.

.. autofunction:: bremap.inspect_input

   Returns a dictionary with ``valid``, ``errors``, ``warnings``, shape, program
   coverage and QC information. Reference/calibration roles require state labels.

Fitting, calibration and mapping
--------------------------------

.. autofunction:: bremap.fit_reference

   ``reference`` is a prepared H5AD, ``spec`` a ``ReferenceSpec`` or specification
   path, and ``model_dir`` a new artifact directory. ``policy`` must be explicit;
   use the implemented ``normal_donor_v1`` policy. Returns a reference artifact.

.. autofunction:: bremap.calibrate

   ``model`` is a reference artifact or path. ``calibration`` supplies separate
   donor/cell roles and known labels. Returns a new calibrated artifact;
   the original reference is not changed.

.. autofunction:: bremap.map_query

   Requires a calibrated artifact/path and prepared query H5AD. Returns an
   artifact with saved ``cells.csv`` and ``donors.csv``. See :doc:`outputs`.

.. autoclass:: bremap.Artifact

   Access ``path`` and ``metadata``. ``ReferenceModel``, ``CalibratedReference``
   and ``MappingResult`` are aliases of this container, not distinct estimators.

Evaluation and sensitivity
--------------------------

.. autofunction:: bremap.evaluate

   Accepts a mapping artifact/path and an annotation CSV path. The annotation
   columns are ``cell_id``, ``donor_id`` and ``label``. Returns an ``AnalysisResult``.

.. autofunction:: bremap.assess_reliability

   Assesses donor-bootstrap and optional reference-study omission sensitivity
   with frozen representation/calibration. Returns an ``AnalysisResult``; inspect
   valid draw counts and undefined reasons.

.. autoclass:: bremap.AnalysisResult

   Adds ``tables`` to the artifact's ``path`` and ``metadata``. ``EvaluationResult``
   and ``ReliabilityResult`` are aliases of this container.

Reference and representation sharing
------------------------------------

.. autofunction:: bremap.export_reference

   Returns the archive ``Path``. The export contains fitted arrays and identifiers,
   not counts or encoder weights. Export a calibrated reference for direct reuse.

.. autofunction:: bremap.import_reference

   Validates archive members/checksums and returns a restored artifact. The
   destination must not already exist.

.. autofunction:: bremap.export_embeddings

   Exports existing vectors and selected metadata, not a new encoder model.
   The count-free export is not a valid replacement for mapping count inputs.

Composed workflow
-----------------

.. autofunction:: bremap.workflow.run_workflow

   ``recipe`` is a Python dictionary or JSON filename. Returns a ``Path`` to the
   completed output folder. ``base_dir`` is only for dictionary inputs;
   ``out_dir`` overrides the recipe output. ``overwrite`` is restricted to owned,
   unchanged workflow outputs. See :doc:`reuse` for the two recipe forms.

Scope of this reference
-----------------------

These are the supported user-facing interfaces for the documented workflow.
The package also exposes lower-level numerical/legacy interfaces, which are not
required by these tutorials. This page does not import private research modules
or present them as installed-package capabilities.
