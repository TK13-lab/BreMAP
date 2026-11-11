Release, citation and documentation builds
==========================================

Current snapshot
----------------

This guide documents BReMAP 0.4.1. Package capabilities include custom reference
states/programs, count SVD or supplied frozen embeddings, fitting and calibration,
query mapping, annotation evaluation, reference sensitivity and portable resources.
The composed Python API adds input inspection, program summaries and a single
reproducibility record around those operations.

The source repository is `TK13-lab/BreMAP <https://github.com/TK13-lab/BreMAP>`_.
The code is licensed under MIT. The `stable documentation
<https://bremap.readthedocs.io/en/stable/>`_ is hosted on Read the Docs.
A public PyPI release and accepted manuscript are not claimed in this snapshot.

Citation
--------

Record the package version and reference/configuration provenance in your Methods.
An archived software DOI and final manuscript citation are not yet assigned in
this snapshot. Use the citation supplied with the eventual release; this guide
does not invent a publication or DOI. Cite the original data and encoder studies
used in your own analysis as well.

Build the documentation
-----------------------

From the package root in an activated Python environment:

.. code-block:: console

   python -m pip install .
   python -m pip install -r docs/user/requirements.txt
   python -m sphinx -W --keep-going -b html docs/user outputs/documentation/html
   python -m http.server 8000 --bind 127.0.0.1 --directory outputs/documentation/html

Open ``http://127.0.0.1:8000`` in a browser. The static ``index.html`` also opens
directly from disk. Generated HTML stays outside git and outside the manuscript
revision folders. Rebuild the same output directory; do not create numbered copies.

``autodoc`` imports the installed package. It does not load research scripts,
download datasets or execute tutorials while building. Run the downloadable
quickstart separately to check that the documented workflow still works.

Read the Docs deployment
------------------------

The repository includes ``.readthedocs.yaml`` configured to install BReMAP and
the documentation dependencies before building ``docs/user/conf.py``. It targets
Python 3.12 and fails on Sphinx warnings.

For a public release:

1. Use the standalone ``TK13-lab/BreMAP`` repository. Configuration paths assume
   the package is the repository root. Data and model weights are not bundled.
2. Run the strict local build and installed-wheel quickstart checks.
3. Connect that repository to an author-controlled Read the Docs project.
4. Check the hosted build, canonical URL, search and enabled release versions.
5. Only then add the verified documentation URL to the package/manuscript.

The YAML configuration alone does not create a hosted site or reserve a name.
Read the Docs' `Sphinx deployment guide
<https://docs.readthedocs.com/platform/stable/intro/sphinx.html>`_ describes the
hosting setup. The `TCRi documentation <https://tcri.readthedocs.io/en/stable/index.html>`_
provided a user-journey reference for this documentation design.
