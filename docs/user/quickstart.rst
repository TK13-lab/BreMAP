Quickstart
==========

This example runs without a dataset download. It invents 18 cells in each of three
partitions, with three distinct donors per partition and two synthetic states.
The genes illustrate the API; neither these profiles nor their output are a
biological validation or an approved immune signature.

Run the example
---------------

After :doc:`installation`, download :download:`quickstart.py <examples/quickstart.py>`
and run it from any working folder:

.. code-block:: console

   python quickstart.py --out bremap_quickstart

The chosen output directory must not already exist. You should see:

.. code-block:: text

   BReMAP 0.4.1: mapped 18 synthetic query cells / 3 donors
   Reference reuse: identical mapping; program and evaluation tables written

The script checks that importing the saved reference gives exactly the same
mapping as the initial run. It uses only the installed BReMAP package and its
scientific dependencies, not repository research scripts.

Read the outputs
----------------

Open these files beneath ``bremap_quickstart/workflow``:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - File
     - Purpose
   * - ``mapping/cells.csv``
     - Candidate states, conditional affinities and acceptance for each query cell.
   * - ``mapping/donors.csv``
     - Donor-level mapping summaries.
   * - ``programs/cells.csv`` and ``programs/donors.csv``
     - Expression-program scores and donor summaries.
   * - ``evaluation/summary.csv``
     - Annotation-based evaluation of the saved predictions.
   * - ``reference.zip``
     - Portable fitted and calibrated reference, without source counts.
   * - ``workflow.json`` and ``input_inspection.json``
     - Settings, provenance, checksums, timing and input checks.

See :doc:`outputs` before interpreting results from biological data.

Complete executable source
--------------------------

This is the same file offered for download and used for the installed-package
check, so the displayed example cannot diverge from a separately maintained copy.

.. literalinclude:: examples/quickstart.py
   :language: python
   :linenos:

Next, replace the synthetic matrices using :doc:`own-data`, or read :doc:`reuse`
to perform fitting, calibration and mapping as individual Python calls.
