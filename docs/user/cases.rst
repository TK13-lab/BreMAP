Application case studies
========================

These cases explain what users can do with reference profiles and how strong
the validation is. They are not three universally validated preset models.
The current package exposes reusable mapping APIs; study-specific regression,
cohort sampling and manuscript rendering are separate analyses.

Step-by-step guides
-------------------

Start with the common runner and checks, then follow your biological case.
Each guide includes commands, input contracts, recorded output and the existing
study visualization. Recorded biological results require the retained study
inputs; they are not outputs of the synthetic quickstart.

.. toctree::
   :maxdepth: 1

   case-workflow
   case-cd4
   case-multilineage
   case-feline

.. list-table:: Choose the example closest to your question
   :header-rows: 1
   :widths: 24 36 40

   * - Case
     - Research question
     - Validation boundary
   * - CD4 CITE-seq
     - Do RNA-derived reference coordinates track a measured protein axis?
     - Protein measurements provide an additional modality; source labels have multimodal information.
   * - Wu/Kim multilineage
     - Can the same API describe different populations using user-defined states and programs?
     - Within-source, donor-disjoint annotation transfer; labels remain RNA-derived.
   * - FELINE
     - How much endpoint information is retained after mapping native representations to reference coordinates?
     - Exploratory, concurrent RNA-score changes; not pretreatment clinical-response prediction.

CD4: compare RNA profiles with measured protein
-----------------------------------------------

**Input:** counts, donor metadata and matched CD45RA/CD45RO protein measurements
from Hao PBMC CITE-seq, `GSE164378 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378>`_.
The evaluated subset contained 1,548 CD4 cells from eight donors.

**Workflow:** prepare donor-held-out reference/calibration/query roles; exclude
``PTPRC`` from RNA inputs; map the query; join saved cell coordinates to the
independently retained protein table by cell ID. Summarize concordance per donor.
For prediction, use a donor-held-out downstream model and compare RNA programs,
RNA plus mapped coordinates, and RNA plus a direct representation.

**Observed result:** mean within-donor correlation with the measured protein axis
was 0.6675, at 88.74% mean acceptance. RNA plus transport/biology coordinates had
MSE 1.0706 versus 1.5484 for RNA programs alone, a 30.86% reduction. RNA plus direct
SVD was better still, at 0.9533. The result supports incremental information beyond
that program baseline, not superiority to SVD or validation in every immune population.

**Adaptation:** keep your protein or functional measurement separate from inputs
used to define the reference. Declare whether program genes should be withheld,
retain donor-level replication and compare against a reasonable simpler readout.
``br.evaluate`` evaluates annotation labels; it does not run this protein-axis
regression automatically.

Wu/Kim: change the population, reference states and programs
------------------------------------------------------------

**Input:** source counts and annotation/patient crosswalks for Wu breast
`GSE176078 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078>`_ and Kim lung
`GSE131907 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907>`_.
The application panel includes immune, stromal, epithelial and vascular tasks
within each source: eight tasks, 14,759 selected cells and 4,828 query cells.

**Workflow:** define one ``ReferenceSpec`` for each population; select relevant
source states and programs; prepare disjoint donors; run affinity and transport
with the same declared inputs. Evaluate saved predictions against the source
annotation and inspect program/affinity profiles. Wu had supplied SVD, scGPT and
Geneformer representations; Kim was evaluated with SVD only. There were 32
mapping runs and 32 comparator fits, not 64 independent cohorts.

**Observed result:** donor-macro acceptance across mapping configurations ranged
from 85.17% to 97.76%. Kim immune SVD transport had Brier loss 0.1605 versus 0.2510
for the direct-SVD-plus-program comparator, but lower F1 (0.8803 versus 0.9052).
Kim vascular affinity likewise improved Brier while slightly reducing F1.
Stromal and epithelial comparisons did not establish a uniform BReMAP advantage.

**Adaptation:** use :doc:`own-data` to substitute your state vocabulary and genes,
then :doc:`reuse` to run the same public API. Keep both successes and losses.
These source annotations and RNA programs are not independent functional assays;
one fixed source-specific donor split is not cross-cohort validation. The study's
donor-state cell cap also prevents interpreting the selected cell frequencies
as natural tissue composition.

FELINE: compare native features with reference coordinates
----------------------------------------------------------

**Input:** serial breast-tumor RNA measurements, patient/visit metadata and
compatible supplied embeddings from FELINE,
`GSE211434 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211434>`_.
The analyzed material contains 30,513 cells, 78 biopsies and 27 patients.

**Workflow:** map the biopsies using a fixed reference configuration. Construct
patient/biopsy-level native-feature summaries and mapped-coordinate summaries.
Keep the downstream endpoint, patient-held-out fitting and comparison procedures
matched. In this study, nine validation genes were excluded from encoder inputs;
that property is specific to the recorded extraction, not guaranteed by BReMAP.

**Observed result:** among 26 day-0/day-14 pairs, direct Geneformer features had
MSE 0.027081 for estrogen-mean change, versus 0.051020 for the RNA baseline and
0.048823 for Geneformer-based mapped coordinates. That gain did not recur against
RNA in the other five tasks. At day 0–180, the ordering reversed: mapped
Geneformer MSE was 0.039268 versus 0.071858 for native Geneformer.

**Adaptation:** use a downstream endpoint to assess whether reference compression
retains the information your question needs. Do not assume that fewer,
interpretable coordinates improve prediction. These analyses use both visits
and same-assay RNA targets; they are not pretreatment treatment-response prediction.

Reproduction boundary
---------------------

The :doc:`quickstart` is fully executable without external data. These biological
case walkthroughs require the specified study inputs, cell/sample joins and
recorded analysis definitions. Public accession links do not imply that processed
inputs, foundation weights, or the paper's complete downstream scripts are
distributed with the installed package. An author-provided reproducibility
release must identify the exact manifests and files before claiming paper-exact
reproduction. Until then, use these pages as application and interpretation guides.
