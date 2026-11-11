Read outputs and use them downstream
====================================

Cell-level mapping
------------------

The individual ``map_query`` API writes ``cells.csv`` and ``donors.csv`` beneath
its output path. The composed workflow places them under ``mapping/``.

.. list-table:: Key cell columns
   :header-rows: 1
   :widths: 38 62

   * - Column
     - Meaning
   * - ``cell_id``, ``donor_id``
     - Identifiers for exact metadata joins.
   * - ``candidate_state``
     - The proposed state, including for a cell that is not accepted.
   * - ``accepted``, ``accepted_analog``
     - Whether the acceptance rule was met, and the accepted label when applicable.
   * - ``conditional_affinity::<state>``
     - Relative affinity among the declared known states for scored cells.
   * - ``reference_distance``
     - Geometric mismatch to the reference; its meaning differs from probability.
   * - ``reject_mass``
     - Unresolved mass for transport; not a new cell type.
   * - ``abstention_reason``, ``qc_zero_library``
     - Diagnostic reason for non-acceptance and zero-count-cell flag.

Do not replace missing scores for zero-library cells with biological zeros.
Report the denominator for any summary: all query cells, scored cells or accepted
cells can give different answers. A visually strong affinity can coexist with
rejection because conditional affinity and reference compatibility answer
different questions.

Program profiles
----------------

``run_workflow`` writes one row per cell in ``programs/cells.csv`` and a long-form
donor table in ``programs/donors.csv``. The latter contains ``donor_id``, ``program``,
``n_cells``, ``n_scored``, ``mean`` and ``sd``. Scores use mean log1p(CP10k) expression
on the frozen reference gene universe. The SD is the population SD (``ddof=0``),
not an uncertainty interval for the donor mean.

These program summaries include score-defined cells regardless of mapping
acceptance. An accepted-only summary is a different analysis and should be
identified as such. Missing genes or differences in the gene universe can change
the meaning of a score; inspect ``input_inspection.json`` and reference metadata.

.. code-block:: python

   import pandas as pd

   donor_programs = pd.read_csv(result_dir / "programs/donors.csv")
   features = donor_programs.pivot(index="donor_id", columns="program", values="mean")
   # clinical_metadata must have one row per donor_id.
   analysis = features.merge(clinical_metadata, on="donor_id", validate="one_to_one")

This produces a feature table, not a fitted clinical predictor. For serial
samples, choose and record the sample/patient aggregation appropriate to the
question, retain patient grouping, and use a downstream analysis that accounts
for repeated measurements. Do not treat thousands of cells as thousands of
independent patient outcomes.

Evaluation
----------

``evaluate`` joins saved mappings to a separate annotation table containing
``cell_id``, ``donor_id`` and ``label``. It returns an ``AnalysisResult`` with
``path``, ``metadata`` and ``tables``. Evaluation reports pooled and donor-level
results, eligible counts and reasons for undefined values.

* Brier loss measures squared discrepancy between the conditional affinity
  vector and a one-hot known-state label; **lower is better**. The multiclass
  implementation sums over states without a division by the number of states.
* Macro-F1 summarizes hard candidate-state classification; **higher is better**.
  This does not automatically mean accepted-only performance.
* Acceptance is the fraction passing the operational rule. High acceptance is
  not itself evidence of correct biology.
* Risk/coverage summaries describe changing acceptance; state their denominator
  and inspect unavailable endpoints rather than silently replacing them with zero.

Source RNA annotations can support an annotation-transfer evaluation. They are
not independent protein measurements or functional validation. Different state
vocabularies and patient populations also prevent interpreting raw scores across
tasks as a common difficulty scale.

Reference sensitivity
---------------------

``assess_reliability`` perturbs the reference while keeping the representation
and calibration frozen. Donor bootstrap and optional reference-study omission
summarize sensitivity, including the number of valid perturbations. Agreement
across these runs is not proof of biological accuracy, and these outputs are not
automatically confidence intervals for a downstream clinical effect.

Portable records
----------------

``workflow.json`` records configuration, versions, hashes and stage timing.
Process peak memory includes the process lifetime; isolate runs before comparing
resource measurements. Reference archives contain fitted numerical parameters
and recorded cell/donor identifiers, but no raw count matrices or encoder weights.
Choose identifiers suitable for the intended sharing context.
