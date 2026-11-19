BReMAP
======

Configurable reference and gene-program profiling for single-cell RNA data.

BReMAP turns a declared reference population into reusable coordinates for query
cells. Define the states and gene programs relevant to your question, fit and
calibrate a reference, then map new samples using Python APIs. Read the resulting
state affinities, RNA-program profiles and acceptance diagnostics together.

References can describe immune, stromal, epithelial, vascular or other eligible
human cell populations. Neither normal epithelium nor breast tissue is required.
Choose expression SVD, supply compatible embeddings, or generate them with an optional
local scGPT/Geneformer adapter. Train a small student on your own paired data when
you need repeated CPU inference within a validated task. Checkpoint downloads require confirmation; foundation models are not retrained.

Start here
----------

* :doc:`installation`: install the package and check your environment.
* :doc:`quickstart`: run a complete, tiny synthetic example without a dataset download.
* :doc:`own-data`: replace that example with your reference, expression and programs.
* :doc:`distillation`: generate teacher vectors and train your own scoped student.
* :doc:`scoped-student`: follow a measured T-cell application, including its limits.
* :doc:`outputs`: understand what each table can support downstream.

The workflow
------------

.. code-block:: text

   Your question + annotated reference + gene programs
                         |
           Prepare disjoint donor partitions
                         |
              Fit reference -> Calibrate
                         |
                Map query + Score programs
                         |
       Inspect / Evaluate -> Export reference for reuse

These steps are public functions in the installed package. You can compose them
in a notebook or call :func:`bremap.workflow.run_workflow` with a Python dictionary.
No repository research scripts or JSON recipe file are required.

What the output means
---------------------

An affinity describes resemblance to the **states you supplied**. A program score
summarizes expression of the **genes you supplied**. Acceptance reports whether
the query meets the calibrated mapping rule. These are useful reference-relative
descriptions; an application-specific measurement is needed to establish a
functional or clinical conclusion. See :doc:`cases` for both supporting evidence
and unfavorable comparisons.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User guide

   own-data
   processed-inputs
   input-options
   representations
   distillation
   scoped-student
   outputs
   reuse
   cases

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api
   troubleshooting
   release
