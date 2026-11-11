Run and check a prepared case
=========================================

This shared procedure applies to all three case studies. Run commands from the
cloned BReMAP repository, after following :doc:`installation`. The example files
are tutorial clients of the installed package; they do not import research scripts.

1. Check your installation
--------------------------------------

.. code-block:: console

   python -c "import bremap; print(bremap.__version__)"

Expected version for this tutorial release: ``0.4.1``. The original study release
was ``0.4.0``; version 0.4.1 adds documentation and checks without changing the
mapping algorithms. Update your checkout and reinstall to obtain these examples.
A successful installation alone does not provide biological datasets.

2. Assemble the inputs
----------------------------------

Create a directory called ``case_inputs`` containing:

.. code-block:: text

   case_inputs/
     spec.json          # ReferenceSpec saved as JSON
     reference.h5ad     # Raw counts + donor/state metadata
     calibration.h5ad   # Separate donors, same declared state vocabulary
     query.h5ad         # Separate donors; query labels not needed for mapping
     annotations.csv    # Optional: cell_id, donor_id, label, for evaluation only

Use :doc:`own-data` to create the specification and call ``prepare_input`` for
all three roles. Save the specification with:

.. code-block:: python

   import json
   from pathlib import Path
   Path("case_inputs/spec.json").write_text(json.dumps(spec.data, indent=2))

Do not rename source labels to fit a vocabulary without a documented crosswalk.
For supplied embeddings, follow :doc:`representations` before input preparation.

The study accession is a source of data, not a download link for this prepared
folder. The processed paper inputs are not currently bundled in the public repo.
The case pages specify what must be prepared; reported study output is not an
assertion that a new donor split or alternative preprocessing will match it.

3. Validate and run the package APIs
------------------------------------------------

.. code-block:: console

   python docs/user/examples/run_case.py --inputs case_inputs --out case_run --method affinity

Use a fresh output folder. The runner checks each input and donor separation
before calling ``run_workflow``. For transport, change ``--method transport`` and
use another output directory; retain whole query donor groups.

:download:`Download the runner <examples/run_case.py>` or inspect the complete API calls:

.. literalinclude:: examples/run_case.py
   :language: python
   :start-at: def run(
   :end-before: if __name__

Recorded output from the retained Kim immune SVD inputs:

.. literalinclude:: examples/expected/kim_run.txt
   :language: text

4. Check mapping output before interpreting it
----------------------------------------------------------

.. code-block:: console

   python docs/user/examples/check_mapping.py case_run/mapping

:download:`Download the checker <examples/check_mapping.py>`.
It tests unique cell identifiers, donor presence, Boolean acceptance flags,
nonnegative affinities and row sums of one. Zero-library cells may be unscored;
the checker reports them rather than replacing missing values with zeros.

Recorded Kim immune SVD-affinity output:

.. literalinclude:: examples/expected/kim_immune.txt
   :language: text

Cell-weighted and donor-macro acceptance are different summaries. Candidate state
counts include rejected cells and are not estimates of tissue composition.
Passing these checks establishes software consistency, not biological correctness.

5. Inspect programs and annotation performance
----------------------------------------------------------

.. code-block:: python

   from pathlib import Path
   import pandas as pd

   result = Path("case_run")
   programs = pd.read_csv(result / "programs/donors.csv")
   print(programs.head().to_string(index=False))
   print(programs.pivot(index="donor_id", columns="program", values="mean"))
   # Available only when annotations.csv was supplied:
   print(pd.read_csv(result / "evaluation/summary.csv").to_string(index=False))

Program output has columns ``donor_id, program, n_cells, n_scored, mean, sd``.
Read :doc:`outputs` for denominators, Brier loss, F1 and acceptance. Compare methods
on identical query cells; do not select the winning reference using the final test set.

6. Reuse the reference
----------------------------------

``case_run/reference.zip`` stores the calibrated reference. Follow :doc:`reuse`
to apply it to a new compatible query. Input validation, compatible feature
provenance and donor-level assessment are still needed for the new cohort.
