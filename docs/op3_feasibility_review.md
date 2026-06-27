# OP3 Feasibility Review

## Bottom Line

OP3 should not be treated as a second independent validation dataset for this manuscript. The current OpenProblems NeurIPS 2023 perturbation prediction data used in our benchmark is the same benchmark family and effectively the same biological dataset underlying OP3.

Best use in the paper:

- Cite OP3 as the closest benchmark context.
- Sharpen the difference from OP3.
- Do not present OP3 as an added independent dataset unless using a genuinely distinct OP3-derived raw-data reprocessing question, which would still not solve the "single dataset" critique.

## Evidence

OpenReview/NeurIPS description:

- OP3 is described as a benchmark for predicting small-molecule perturbation effects on cell-type-specific gene expression.
- It is enabled by a new dataset of 146 compounds tested on human blood cells.
- It includes data representations, metrics, and winning methods from the NeurIPS 2023 competition.
- Source: https://openreview.net/forum?id=WTI4RJYSVm

OpenProblems task repository:

- The task is the NeurIPS 2023 Kaggle/OpenProblems perturbation prediction competition.
- It states that the dataset was generated in human PBMCs, using small molecules selected from LINCS CMap, with 24-hour treatment.
- It exposes the same task files used in this project, including `de_train`, `de_test`, and `id_map`.
- Source: https://github.com/openproblems-bio/task_perturbation_prediction

GEO record GSE279945:

- The title is "OP3: single-cell multimodal dataset in PBMCs for perturbation prediction benchmarking."
- The summary says compounds were dosed in PBMCs and perturbation effects were computed across four cell populations: T cells, B cells, NK cells, and myeloid cells.
- The design reports 146 small molecules, three donors, and baseline multiome profiling.
- Supplementary processed files include a 13.2 GB scRNA-seq H5AD and a 1.5 GB multiome H5AD.
- Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945

## Implication For Our Manuscript

The current manuscript's dataset can be described as the OP3/OpenProblems PBMC small-molecule perturbation dataset. The high-score revision should avoid wording that implies OP3 is a separate dataset still to be added.

The correct distinction is:

> OP3 established a standardized perturbation-prediction benchmark over the PBMC small-molecule dataset and evaluated competition methods under benchmark-defined settings. Our study asks a narrower validation question: whether random or weakly constrained splits preserve chemical-neighbor leakage, and whether Bemis-Murcko scaffold-strict and joint cell-plus-scaffold-held-out splits expose a larger generalization failure. This adds explicit scaffold-level leakage quantification and model-ranking instability analysis rather than proposing another general-purpose benchmark leaderboard.

## Practical Options

Option A: Use OP3 as prior-work context only.

- Fastest and cleanest.
- Strengthens the Introduction and Discussion.
- Does not solve the one-dataset limitation.

Option B: Reprocess GSE279945 raw/processed counts.

- Possible but heavy: processed scRNA-seq file is 13.2 GB.
- Could audit donor-level or cell-type-level variants not available in `de_train/de_test`.
- Still remains the same biological study, so it should be framed as an internal robustness analysis rather than a second dataset.

Option C: Choose a true second dataset.

- Best answer to the reviewer risk.
- Candidates remain Sci-Plex/scPerturb or another public chemical perturbation single-cell dataset with SMILES and at least two cellular contexts.
- This is the recommended high-score path.

## Manuscript Action

Add a paragraph in the Introduction:

> The OP3 benchmark provided an important standardized setting for small-molecule perturbation prediction across PBMC cell types. However, benchmark-level performance comparisons do not by themselves quantify how much chemical-neighbor information is retained across train-test splits. We therefore focus on a complementary audit: explicit scaffold-based split construction, nearest-training-drug Tanimoto analysis, and random-to-strict model-ranking instability.

Add a paragraph in the Discussion:

> Our analysis is complementary to OP3 rather than a replacement for it. OP3 evaluates prediction methods in a standardized competition framework; the present study isolates a narrower validation failure mode, showing that random splitting can preserve near-complete drug and scaffold overlap and can select models that do not remain top-ranked under scaffold-strict or joint held-out validation.
