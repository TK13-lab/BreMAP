# BReMAP

Configurable reference mapping and gene-program analysis for single-cell RNA-seq.
BReMAP provides Python APIs to prepare inputs, define reference states and gene
programs, fit and calibrate a reference, map query cells, evaluate annotations,
and export reusable reference resources in one reproducible workflow.

References and programs are user-defined. Immune, stromal, epithelial and vascular
populations can be studied when suitable reference data and validation are available.
Use count-derived SVD or compatible precomputed embeddings, including scGPT and
Geneformer representations. Optional TeacherEncoder adapters generate embeddings
from your counts using pretrained scGPT/Geneformer checkpoints. Select an official
checkpoint with `TeacherEncoder.from_pretrained`; BReMAP asks before downloading
and reuses its verified cache. Existing local checkpoints also work. StudentEncoder
trains a compact model from paired counts and teacher vectors for your own task.
Foundation models are optional; no foundation-model retraining or universal student
checkpoint is required or claimed.

## Install

Python 3.10 or newer is required; Python 3.12 is used for documentation checks.

```bash
git clone https://github.com/TK13-lab/BreMAP.git
cd BreMAP
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install ".[plot]"
```

This installs BReMAP 0.7.0 from source. No PyPI release is currently claimed.
NumPy, SciPy, pandas, scikit-learn and AnnData are installed as dependencies.
Neither R nor a GPU is required for count-SVD mapping or student inference.
Student training additionally uses `.[train]`; teacher extraction needs the
upstream encoder environment and checkpoint. Add `.[download]` to enable confirmed
checkpoint downloads. Geneformer extraction requires CUDA.

## Run the tutorial

For an interactive, API-led walkthrough with input loading, four plots and
recorded output, start with the [Quickstart](https://bremap.readthedocs.io/en/stable/quickstart.html).
The following script remains a minimal non-plotting smoke test.

```bash
python docs/user/examples/quickstart.py --out bremap_quickstart
```

The example creates synthetic inputs, calls only installed BReMAP APIs, maps
18 query cells from three donors, exports program/evaluation tables and verifies
that a reused reference produces identical mappings. Use a new output destination.
The synthetic programs demonstrate software usage, not biological validation.

Start with `bremap.read_data`, `bremap.create_spec`, `bremap.prepare_partitions` and
`bremap.run_analysis`; use the individual APIs for finer control.
For encoder targets and your own student, see the
[distillation guide](https://bremap.readthedocs.io/en/stable/distillation.html) and
[scoped T-cell tutorial](https://bremap.readthedocs.io/en/stable/scoped-student.html).
Read the [quickstart](https://github.com/TK13-lab/BreMAP/blob/main/docs/user/quickstart.rst)
and [own-data tutorial](https://github.com/TK13-lab/BreMAP/blob/main/docs/user/own-data.rst)
before preparing a real dataset. Split donors across development and query data;
keep held-out outcome labels outside model fitting and calibration.

## Documentation

Read the **[online documentation](https://bremap.readthedocs.io/en/stable/)**
for the stable release.

The Sphinx guide covers installation, input schemas, representations, output
interpretation, reference reuse, case studies, troubleshooting and API reference.
GitHub Actions checks the installed-package tutorial and documentation build.

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
