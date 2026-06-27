# 100-seed Upgrade Results

## Completed Analyses

- Repeated OpenProblems split generation and baseline evaluation from 10 to 100 seeds.
- Added leakage overlap diagnostics for same-drug, same-scaffold, same-cell, and same drug-cell-pair overlap.
- Recomputed nearest-training-drug Tanimoto bootstrap intervals with 100 seed-level observations.
- Added model-ranking instability audit using mean DE-gene rowwise Pearson.

## Core Performance Result

For the ridge model using cell context plus drug fingerprint, random validation remained substantially more optimistic than scaffold-aware validation across 100 seeds.

| Split | Mean test records | All-gene Pearson, mean +/- SD | DE-gene Pearson, mean +/- SD |
|---|---:|---:|---:|
| Random | 108.41 | 0.362 +/- 0.016 | 0.568 +/- 0.018 |
| Cell held-out | 136.34 | 0.313 +/- 0.056 | 0.496 +/- 0.089 |
| Scaffold held-out | 106.70 | 0.234 +/- 0.027 | 0.392 +/- 0.037 |
| Cell + scaffold held-out | 26.76 | 0.118 +/- 0.034 | 0.209 +/- 0.059 |

Suggested manuscript wording:

> Across 100 repeated split seeds, the cell-context plus drug-fingerprint ridge model achieved a mean DE-gene Pearson correlation of 0.568 under random splitting, but this decreased to 0.392 under scaffold-held-out splitting and to 0.209 under joint cell-plus-scaffold-held-out splitting. Thus, the apparent gain under random splitting was not preserved when chemical series and cell contexts were both withheld.

## Leakage Overlap Result

Random splitting left almost all test drugs and scaffolds represented in the training set, whereas scaffold-strict and joint splits removed this overlap by construction.

| Split | Same drug in train | Same scaffold in train | Same cell in train | Mean test records |
|---|---:|---:|---:|---:|
| Random | 0.989 | 0.975 | 1.000 | 109.97 |
| Cell held-out | 1.000 | 0.986 | 0.000 | 138.34 |
| Scaffold held-out | 0.000 | 0.000 | 1.000 | 108.30 |
| Cell + scaffold held-out | 0.000 | 0.000 | 0.000 | 27.16 |

Suggested manuscript wording:

> The random split placed a same-drug training example next to 98.9% of test observations and a same-scaffold training example next to 97.5% of test observations. In contrast, the scaffold-held-out and joint held-out splits reduced both same-drug and same-scaffold overlap to zero, demonstrating that the random split was not a chemically independent generalization test.

## Chemical Similarity Result

The nearest training-drug Tanimoto similarity remained close to identity under random splitting and dropped sharply under scaffold-aware splits.

| Split | Mean max training-drug Tanimoto | 95% bootstrap CI |
|---|---:|---:|
| Random | 0.992 | 0.989 to 0.995 |
| Scaffold held-out | 0.278 | 0.275 to 0.281 |
| Cell + scaffold held-out | 0.278 | 0.275 to 0.282 |

Contrasts:

| Contrast | Mean difference | 95% bootstrap CI |
|---|---:|---:|
| Random minus scaffold held-out | 0.714 | 0.710 to 0.718 |
| Random minus cell + scaffold held-out | 0.714 | 0.709 to 0.718 |

Suggested manuscript wording:

> The average maximum Tanimoto similarity between each test drug and its nearest training drug was 0.992 under random splitting, compared with 0.278 under scaffold-held-out and joint held-out splitting. The random-minus-scaffold difference was 0.714, with a 95% seed-level bootstrap interval of 0.710 to 0.718.

## Ranking Instability Result

Random-split model ranking did not transfer to scaffold-strict or joint generalization.

| Split | Best average-ranked baseline | Mean rank | Mean DE-gene Pearson |
|---|---|---:|---:|
| Random | ridge_cell_drug_fp | 1.00 | 0.568 |
| Cell held-out | ridge_cell_drug_fp | 1.24 | 0.496 |
| Scaffold held-out | cell_context_mean | 1.08 | 0.477 |
| Cell + scaffold held-out | cell_context_mean / drug_mean / global_train_mean | 1.24 | 0.288 |

Random-vs-strict Spearman rank correlations:

| Comparison | Mean Spearman rho +/- SD |
|---|---:|
| Random vs cell held-out | 0.531 +/- 0.277 |
| Random vs scaffold held-out | 0.097 +/- 0.436 |
| Random vs cell + scaffold held-out | -0.574 +/- 0.191 |

Suggested manuscript wording:

> Random splitting also distorted model selection. The drug-fingerprint ridge model ranked first in every random split, but its average rank fell to 3.58 under scaffold-held-out splitting and 5.33 under joint held-out splitting. Across seeds, the Spearman rank correlation between random and joint held-out rankings was negative (mean rho = -0.574), indicating that random validation can select models that fail under chemically independent generalization.

## Ridge Hyperparameter Sensitivity

Output:

- `results/ridge_sensitivity_30/all_ridge_sensitivity_summary.csv`
- `results/ridge_sensitivity_30/ridge_sensitivity_by_config.csv`
- `results/ridge_sensitivity_30/ridge_cell_drug_fp_sensitivity.csv`

The sensitivity audit used 30 split seeds and focused on the `ridge_cell_drug_fp` baseline. Alpha was varied at 0.1, 1, 10, 100, and 1000 with 1024-bit Morgan fingerprints. Fingerprint length was also varied at 512, 1024, and 2048 bits with alpha fixed at 10.

| Alpha | Fingerprint bits | Random DE Pearson | Scaffold DE Pearson | Joint DE Pearson | Random - scaffold | Random - joint |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 1024 | 0.555 | 0.384 | 0.200 | 0.171 | 0.354 |
| 1 | 1024 | 0.556 | 0.385 | 0.203 | 0.171 | 0.354 |
| 10 | 512 | 0.564 | 0.383 | 0.214 | 0.181 | 0.351 |
| 10 | 1024 | 0.563 | 0.396 | 0.218 | 0.167 | 0.345 |
| 10 | 2048 | 0.563 | 0.409 | 0.227 | 0.155 | 0.336 |
| 100 | 1024 | 0.557 | 0.424 | 0.272 | 0.133 | 0.286 |
| 1000 | 1024 | 0.456 | 0.411 | 0.304 | 0.046 | 0.152 |

Suggested manuscript wording:

> The random-to-strict performance decay was not specific to the default ridge setting. Across alpha values from 0.1 to 1000 and fingerprint lengths from 512 to 2048 bits, random-split DE-gene Pearson remained higher than joint cell-plus-scaffold-held-out performance. Stronger regularization reduced the magnitude of this gap, indicating that hyperparameters influence absolute scores and should be reported, but did not remove the central validation-regime effect.

## Direct Response to Reviewer Risks

- Small n=10 bootstrap risk: substantially reduced by rerunning the main lightweight analyses with 100 seed-level split replicates.
- Joint held-out uncertainty: still limited by a small mean test set size of 26.76 records, but now characterized across 100 split seeds.
- Leakage mechanism: now explicitly quantified by both categorical overlap and continuous chemical-neighbor similarity.
- Model-selection risk: now explicitly shown by ranking instability, not only performance decay.
- Ridge hyperparameter risk: partially addressed with 30-seed alpha and fingerprint-length sensitivity.

## Remaining High-score Gaps

- A second dataset is still the largest remaining weakness.
- No true foundation-model or perturbation-specific deep baseline has yet been run.
- SVD dimensionality sensitivity is still pending.
- Submission metadata placeholders still need to be completed before journal submission.
