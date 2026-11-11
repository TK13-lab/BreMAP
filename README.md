# BReMAP

Configurable reference mapping and gene-program analysis for single-cell RNA-seq.
BReMAP provides Python APIs to prepare inputs, define reference states and gene
programs, fit and calibrate a reference, map query cells, evaluate annotations,
and export reusable reference resources in one reproducible workflow.

References and programs are user-defined. Immune, stromal, epithelial and vascular
populations can be studied when suitable reference data and validation are available.
Use count-derived SVD or compatible precomputed embeddings, including scGPT and
Geneformer representations. Foundation models are optional; BReMAP does not train
or run the encoders and does not assume that their features outperform simpler inputs.

## Install

Python 3.10 or newer is required; Python 3.12 is used for documentation checks.

```bash
git clone https://github.com/TK13-lab/BreMAP.git
cd BreMAP
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install .
```

This installs BReMAP 0.4.0 from source. No PyPI release is currently claimed.
NumPy, SciPy, pandas, scikit-learn and AnnData are installed as dependencies.
Neither R nor a GPU is required for the package.

## Run the tutorial

```bash
python docs/user/examples/quickstart.py --out bremap_quickstart
```

The example creates synthetic inputs, calls only installed BReMAP APIs, maps
18 query cells from three donors, exports program/evaluation tables and verifies
that a reused reference produces identical mappings. Use a new output destination.
The synthetic programs demonstrate software usage, not biological validation.

Start with `bremap.create_spec`, `bremap.inputs.prepare_input` and
`bremap.workflow.run_workflow`; use the individual APIs for finer control.
Read the [quickstart](https://github.com/TK13-lab/BreMAP/blob/main/docs/user/quickstart.rst)
and [own-data tutorial](https://github.com/TK13-lab/BreMAP/blob/main/docs/user/own-data.rst)
before preparing a real dataset. Split donors across development and query data;
keep held-out outcome labels outside model fitting and calibration.

## Documentation

The Sphinx guide covers installation, input schemas, representations, output
interpretation, reference reuse, case studies, troubleshooting and API reference.
To build and view it locally:

```bash
python -m pip install -r docs/user/requirements.txt
python -m sphinx -W --keep-going -b html docs/user outputs/documentation/html
python -m http.server 8000 --bind 127.0.0.1 --directory outputs/documentation/html
```

Open <http://127.0.0.1:8000>. Read the Docs configuration is included; a hosted
URL will be added after deployment is verified. GitHub Actions checks the
installed-package tutorial and strict documentation build on pushes and pull requests.

## Scope and citation

BReMAP summarizes a query relative to a declared reference. Affinity is not proof
of cell identity, malignancy or clinical response. Check reference suitability,
coverage and donor-level evidence for the intended task. Data and pretrained
encoder weights are not bundled; their original licenses and access terms apply.

Record the BReMAP version and configuration in your Methods and cite the original
datasets and encoders you use. A software DOI and final manuscript citation have
not yet been assigned.

## License and support

BReMAP code is released under the [MIT License](https://github.com/TK13-lab/BreMAP/blob/main/LICENSE).
Report bugs with a minimal synthetic example at
[GitHub Issues](https://github.com/TK13-lab/BreMAP/issues).
