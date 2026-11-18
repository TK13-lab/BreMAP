Multilineage tutorial: change the reference and programs
====================================================================

1. Choose the source and population
-----------------------------------------------

Wu breast: `GSE176078 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078>`_.
Kim lung: `GSE131907 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907>`_.
Obtain raw counts and the matching cell-label/patient crosswalk. Select immune,
stromal, epithelial or vascular cells before splitting donors. The recorded study
contains eight source/population tasks, not eight independent external cohorts.

2. Define a task-specific specification
---------------------------------------------------

For the Kim immune example:

.. code-block:: python

   import bremap as br

   spec = br.create_spec(
       spec_id="kim_immune", tissue="immune",
       states=["Cytotoxic CD8+ T", "Exhausted CD8+ T", "Naive CD4+ T", "Treg"],
       programs={
           "checkpoint": ["CTLA4", "HAVCR2", "LAG3", "PDCD1", "TIGIT", "TOX"],
           "cycling": ["BIRC5", "MKI67", "TOP2A", "UBE2C"],
           "cytotoxic": ["CTSW", "GNLY", "GZMB", "GZMH", "NKG7", "PRF1"],
           "naive_memory": ["CCR7", "IL7R", "LEF1", "LTB", "MAL", "TCF7"],
       },
       reference_source="Kim source-donor training partition",
       eligibility="Declared source labels; metadata-only donor split",
       minimum_genes_per_anchor=2, cost_mode="biology",
   )

Create separate specifications for other populations; do not apply these immune
genes to stromal cells. Use the source-specific labels and record the crosswalk.
A state vocabulary and a program list need not have the same length.

3. Prepare and run
------------------------------

Use :doc:`own-data` to prepare donor-disjoint ``inputs`` with the specification
above. Compare both methods through the package API:

.. code-block:: python

   affinity_result = br.run_analysis(inputs, spec, "affinity_run", method="affinity")
   transport_result = br.run_analysis(inputs, spec, "transport_run", method="transport")
   print(affinity_result.summary())
   print(transport_result.summary())

Exact paper replication requires the recorded selected-cell and donor manifests,
which are not distributed as public processed inputs.

4. Compare methods on matching cells
------------------------------------------------

.. code-block:: python

   affinity = affinity_result.table("mapping").set_index("cell_id")
   transport = transport_result.table("mapping").set_index("cell_id")
   assert set(affinity.index) == set(transport.index)
   transport = transport.loc[affinity.index]
   assert affinity.donor_id.equals(transport.donor_id)
   affinity_result.plot("coverage")
   transport_result.plot("coverage")

Do not compare Brier losses across different state vocabularies as if task
difficulty were identical. Compare coverage and hard-label F1 as well: the Kim
immune transport result improved Brier against its direct-feature comparator,
but did not improve F1. The direct-feature classifier is a separate downstream
analysis and is not fitted by ``run_workflow``.

5. Read the population profiles
-------------------------------------------

.. image:: _images/multilineage.png
   :alt: Kim immune, stromal, epithelial and vascular program and affinity profiles.
   :width: 100%

This is the retained four-population study figure. Program expression and
conditional affinity have different units and separate color scales. They show
reference-relative organization, not independent functional validation. The
sampling cap per donor/state means plotted counts are not natural cell abundances.
Wu also used supplied scGPT and Geneformer embeddings; Kim used SVD only. To change
representation, follow :doc:`representations` and keep input partitions fixed.
