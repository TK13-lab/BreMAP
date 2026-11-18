Installation
============

Requirements
------------

BReMAP 0.6.2 requires Python 3.10 or newer. Its runtime dependencies are NumPy,
SciPy, pandas, scikit-learn and AnnData; the package metadata declares compatible
minimum versions. The documentation build uses Python 3.12. A GPU and an R
installation are not required to run BReMAP. R was used separately for manuscript
figures and is not a package runtime dependency.

Install from source
-------------------

Clone the public `BReMAP repository <https://github.com/TK13-lab/BreMAP>`_:

.. code-block:: console

   git clone https://github.com/TK13-lab/BreMAP.git
   cd BreMAP
   python -m venv .venv
   # Linux/macOS:
   source .venv/bin/activate
   # Windows PowerShell, instead:
   # .venv\Scripts\Activate.ps1
   python -m pip install .
   python -c "import bremap; print(bremap.__version__)"

The expected package version for this guide is ``0.6.2``. Use ``pip install -e .``
only when you intend to edit the package source during development.

Install a supplied wheel
------------------------

When the authors provide a wheel, install that file in your activated environment:

.. code-block:: console

   python -m pip install /path/to/bremap-0.6.2-py3-none-any.whl

Use the actual supplied filename. No public PyPI release is claimed; use the
repository installation above rather than assuming ``pip install bremap`` resolves
to this project. The source code is distributed under the MIT License.

Check the public interfaces
---------------------------

.. code-block:: python

   import bremap as br
   from bremap.inputs import prepare_input
   from bremap.workflow import run_workflow

   print(br.__version__)
   print(br.__file__)  # Check that this is the intended installation.

Then follow :doc:`quickstart`. It creates its own small synthetic inputs and
requires no foundation-model weights or external dataset.

Optional encoder routes
------------------------

.. list-table:: Install only what your workflow uses
   :header-rows: 1

   * - Route
     - Requirements
   * - Count SVD / supplied embeddings / saved student inference
     - Base package; CPU; no torch or pretrained weights.
   * - Tutorial plots
     - ``python -m pip install '.[plot]'`` from the checkout.
   * - Train your own student
     - ``python -m pip install '.[train]'``; paired teacher vectors; CPU training.
   * - Generate scGPT or Geneformer vectors
     - Official upstream package and compatible local checkpoint, configured in
       ``TeacherEncoder``. These are not included in ``.[train]``.

For both training and plots use ``python -m pip install '.[train,plot]'``.
The teacher adapters were exercised on Linux/WSL. Geneformer extraction requires
CUDA; scGPT's upstream helper uses a Linux CPU-affinity API. Native Windows teacher
extraction is not claimed. Base mapping and saved-student inference do not use
those upstream dependencies. See :doc:`distillation` for checkpoint names,
identifier options, cache behavior and error remedies.
