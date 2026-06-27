# Progress Report

## Current status

The project has moved from concept planning to a working feasibility benchmark and an initial Results-section draft.

Current manuscript draft output:

- `manuscript/title_abstract_keywords_zh.md`
- `manuscript/introduction_draft_zh.md`
- `manuscript/results_draft_zh.md`
- `manuscript/methods_draft_zh.md`
- `manuscript/discussion_draft_zh.md`
- `manuscript/figure_table_plan_zh.md`
- `manuscript/manuscript_full_zh.md`
- `manuscript/manuscript_audit_zh.md`
- `manuscript/tables_draft_zh.md`
- `manuscript/figure_legends_zh.md`
- `manuscript/main_figure_inventory.md`

Current main figure output:

- `figures/manuscript_main/figure1_benchmark_design.*`
- `figures/manuscript_main/figure2_baseline_decay.*`
- `figures/manuscript_main/figure3_topk_response_gene_recovery.*`
- `figures/manuscript_main/figure4_chemical_similarity_audit.*`
- `figures/manuscript_main/figure5_svd_ridge_stress_test.*`
- `figures/manuscript_main/source_data/`

Current submission package output:

- `submission_package/manuscript/manuscript_full_zh.docx`
- `submission_package/rendered_docx/manuscript_full_zh.pdf`
- `submission_package/PACKAGE_README.md`
- `submission_package.zip`

## Environment

A project-level Python environment was created at:

- `env/.venv/`

Installed packages include pandas, numpy, scikit-learn, scipy, matplotlib, seaborn, pyarrow, RDKit, anndata, h5py, requests, and AWS CLI.

Because the workspace path contains spaces, AWS CLI should be called through:

```bash
env/.venv/bin/python -m awscli
```

instead of the generated `aws` entry point.

## Data sources checked

### scPerturb / Zenodo

- Zenodo record was queried successfully.
- File manifest was saved to `metadata/scperturb_zenodo_files.csv`.
- Sci-Plex 3 was downloaded to `data/raw/scperturb/SrivatsanTrapnell2020_sciplex3.h5ad`.
- PubChem SMILES mapping was completed for all 188 Sci-Plex 3 drug perturbations.
- RDKit scaffold audit found 165 unique Bemis-Murcko scaffolds.
- Feasibility report: `docs/sciplex3_feasibility_review.md`.

### OpenProblems / NeurIPS 2023 perturbation prediction

The task repository was cloned to:

- `data/raw/openproblems_task_perturbation_prediction/`

Public S3 resources were inspected and downloaded:

- `data/raw/openproblems_neurips2023/de_train.h5ad`
- `data/raw/openproblems_neurips2023/de_test.h5ad`
- `data/raw/openproblems_neurips2023/id_map.csv`
- `data/raw/openproblems_neurips2023/moa_annotations.csv`
- `data/raw/openproblems_neurips2023/dataset_info.yaml`

## OpenProblems feasibility result

Full train+test metadata:

- 553 perturbation rows
- 545 treated rows
- 140 drugs
- 4 cell contexts
- 140 drugs with SMILES
- 136 valid Bemis-Murcko scaffolds
- 136 scaffolds measured in at least 2 cell contexts
- 135 scaffolds measured in at least 3 cell contexts

Decision:

> Proceed with a full Level 4-5 benchmark on this dataset.

## Candidate split summary

Generated file:

- `metadata/openproblems/splits/candidate_splits.csv`

Split counts:

| Split | Train records | Test records | Excluded records |
|---|---:|---:|---:|
| Random | 443 | 110 | 0 |
| Cell held-out | 413 | 140 | 0 |
| Scaffold held-out | 443 | 110 | 0 |
| Cell + scaffold held-out | 331 | 28 | 194 |

## First simple baseline result

Target layer:

- `clipped_sign_log10_pval`

Output:

- `results/baselines_openproblems/simple_baseline_summary.csv`
- `results/baselines_openproblems/simple_baseline_per_row.csv`
- `figures/openproblems/simple_baseline_decay_mean_pearson.png`
- `figures/openproblems/simple_baseline_decay_mean_pearson.pdf`

Main observation:

The simple baseline performance declines under stricter scaffold and joint cell+scaffold held-out splits. In the joint held-out setting, cell-context and drug-specific mean baselines collapse to the same value as the global training mean, because the split removes the relevant matching shortcuts.

This supports the central manuscript logic: strict split design changes the apparent difficulty of perturbation prediction.

## Ridge baseline result

Output:

- `results/ridge_openproblems/ridge_baseline_summary.csv`
- `results/ridge_openproblems/ridge_baseline_per_row.csv`
- `figures/openproblems/combined_baseline_decay_mean_pearson.png`
- `figures/openproblems/combined_baseline_decay_mean_pearson.pdf`
- `figures/openproblems/combined_baseline_decay_plot_data.csv`

Feature sets:

- `ridge_cell`: cell-context one-hot features.
- `ridge_drug_fp`: Morgan drug fingerprints.
- `ridge_cell_drug_fp`: cell-context one-hot plus Morgan fingerprints.

Key preliminary observations:

- Random split: the combined cell+drug fingerprint ridge baseline achieved the highest mean row-wise Pearson among current baselines.
- Cell-held-out split: drug fingerprint ridge performed best among current simple baselines, indicating that drug-level structure can transfer across held-out cell contexts when the same drugs are represented in training.
- Scaffold-held-out split: drug fingerprint performance dropped below cell-context baselines, consistent with scaffold-strict validation reducing chemical-neighbor shortcuts.
- Cell+scaffold-held-out split: all simple models performed poorly, with drug-fingerprint ridge below the global mean baseline in this seed. This is consistent with the intended strictness of the joint held-out task.

These results are preliminary because only one split seed and one ridge alpha have been evaluated.

## Drug-structure similarity audit

Output:

- `results/drug_similarity_audit_ridge/drug_similarity_vs_error_per_row.csv`
- `results/drug_similarity_audit_ridge/drug_similarity_vs_error_summary.csv`
- `figures/openproblems/drug_similarity_audit_ridge_cell_drug_fp.png`
- `figures/openproblems/drug_similarity_audit_ridge_cell_drug_fp.pdf`

Key preliminary observations:

- Random split had very high nearest training-drug similarity, with mean maximum Tanimoto around 0.97.
- Scaffold-held-out split reduced nearest training-drug similarity to around 0.28.
- Joint cell+scaffold-held-out split also had low nearest training-drug similarity, around 0.28.
- In scaffold-strict settings, higher nearest training-drug similarity was positively associated with row-wise Pearson for ridge baselines in this seed.

This directly supports the manuscript's drug-structure leakage audit: random splits permit chemical-neighbor shortcuts that are strongly reduced by scaffold-strict validation.

## Repeated-seed benchmark

Output:

- `results/openproblems_multiseed_10/all_baseline_summary.csv`
- `results/openproblems_multiseed_10/baseline_summary_by_seed_mean.csv`
- `results/openproblems_multiseed_10/all_drug_similarity_summary.csv`
- `results/openproblems_multiseed_10/drug_similarity_summary_by_seed_mean.csv`
- `figures/openproblems_multiseed/multiseed_baseline_decay_mean_pearson.png`
- `figures/openproblems_multiseed/multiseed_baseline_decay_mean_pearson.pdf`
- `results/nearest_drug_multiseed_10/all_nearest_drug_baseline_summary.csv`
- `results/nearest_drug_multiseed_10/nearest_drug_summary_by_seed_mean.csv`
- `results/openproblems_multiseed_10/all_baseline_with_nearest_summary.csv`
- `results/openproblems_multiseed_10/baseline_with_nearest_summary_by_seed_mean.csv`

Ten split seeds were evaluated.

Key repeated-seed results:

| Split | Baseline | Mean row-wise Pearson | SD |
|---|---|---:|---:|
| Random | Ridge: cell + drug FP | 0.365 | 0.024 |
| Cell held-out | Ridge: cell + drug FP | 0.335 | 0.058 |
| Scaffold held-out | Ridge: cell + drug FP | 0.225 | 0.027 |
| Cell + scaffold held-out | Ridge: cell + drug FP | 0.110 | 0.034 |
| Cell + scaffold held-out | Global mean | 0.168 | 0.031 |
| Cell + scaffold held-out | Nearest drug | 0.071 | 0.021 |

Key repeated-seed similarity audit:

| Split | Mean nearest training-drug Tanimoto | SD |
|---|---:|---:|
| Random | 0.993 | 0.012 |
| Cell held-out | 1.000 | 0.000 |
| Scaffold held-out | 0.279 | 0.014 |
| Cell + scaffold held-out | 0.279 | 0.016 |

Interpretation:

The repeated-seed results support the central benchmark claim. Performance is strongest under random validation, where test drugs are nearly identical or highly similar to training drugs. Scaffold-strict validation sharply reduces chemical proximity and lowers drug-fingerprint ridge performance. In the joint cell+scaffold-held-out task, the drug-fingerprint ridge baseline performs below the global mean, indicating that this task removes both cell-context and chemical-neighbor shortcuts.

Nearest-drug baseline:

- Random split nearest-drug any-cell mean Pearson: 0.216.
- Scaffold-held-out nearest-drug any-cell mean Pearson: 0.083.
- Cell+scaffold-held-out nearest-drug any-cell mean Pearson: 0.071.

This indicates that a simple nearest-neighbor retrieval strategy is not sufficient under scaffold-strict and joint held-out validation.

## DE-gene-only metrics

Output:

- `results/openproblems_multiseed_10_de/all_baseline_with_nearest_summary.csv`
- `results/openproblems_multiseed_10_de/baseline_with_nearest_summary_by_seed_mean.csv`
- `results/nearest_drug_multiseed_10_de/all_nearest_drug_baseline_summary.csv`
- `results/nearest_drug_multiseed_10_de/nearest_drug_summary_by_seed_mean.csv`
- `figures/openproblems_multiseed_de/multiseed_de_baseline_decay_mean_pearson.png`
- `figures/openproblems_multiseed_de/multiseed_de_baseline_decay_mean_pearson.pdf`

Metric definition:

- Target layer: `clipped_sign_log10_pval`
- DE mask: `is_de`
- Minimum DE genes per evaluated row: 10

Key repeated-seed DE-gene-only results:

| Split | Baseline | Mean DE-gene row-wise Pearson | SD |
|---|---|---:|---:|
| Random | Ridge: cell + drug FP | 0.570 | 0.027 |
| Cell held-out | Ridge: cell + drug FP | 0.531 | 0.094 |
| Scaffold held-out | Ridge: cell + drug FP | 0.381 | 0.039 |
| Cell + scaffold held-out | Ridge: cell + drug FP | 0.198 | 0.058 |
| Cell + scaffold held-out | Global mean | 0.290 | 0.054 |
| Cell + scaffold held-out | Nearest drug | 0.131 | 0.029 |

Interpretation:

The DE-gene-only metric gives the same qualitative conclusion as all-gene row-wise Pearson, but with a stronger biological focus. Performance for the cell+drug fingerprint ridge baseline drops from 0.570 under random split to 0.198 under joint cell+scaffold held-out validation. In the strictest split, the global mean baseline remains higher than the drug-structure-informed ridge and nearest-drug baselines, supporting the claim that chemical-neighbor information does not solve unseen drug plus unseen cell-context prediction.

## Top-k DE gene recovery

Output:

- `results/topk_de_overlap_multiseed_10/all_topk_de_overlap_summary.csv`
- `results/topk_de_overlap_multiseed_10/topk_de_overlap_summary_by_seed_mean.csv`
- `figures/openproblems_topk/top100_de_overlap_direction.png`
- `figures/openproblems_topk/top100_de_overlap_direction.pdf`
- `figures/openproblems_topk/top200_de_overlap_direction.png`
- `figures/openproblems_topk/top200_de_overlap_direction.pdf`

Metric definition:

- Top-k overlap compares predicted and observed sets of strongest-response genes.
- Direction agreement evaluates whether overlapping top-k genes have matching response signs.
- Evaluated k values: 50, 100, and 200.

Key repeated-seed top-100 results:

| Split | Baseline | Mean top-100 overlap | Mean direction agreement |
|---|---|---:|---:|
| Random | Ridge: cell + drug FP | 0.138 | 0.975 |
| Scaffold held-out | Ridge: cell + drug FP | 0.090 | 0.864 |
| Cell + scaffold held-out | Ridge: cell + drug FP | 0.050 | 0.762 |
| Cell + scaffold held-out | Global mean | 0.065 | 0.787 |
| Cell + scaffold held-out | Nearest drug | 0.030 | 0.670 |

Key repeated-seed top-200 results:

| Split | Baseline | Mean top-200 overlap | Mean direction agreement |
|---|---|---:|---:|
| Random | Ridge: cell + drug FP | 0.165 | 0.957 |
| Scaffold held-out | Ridge: cell + drug FP | 0.118 | 0.845 |
| Cell + scaffold held-out | Ridge: cell + drug FP | 0.073 | 0.718 |
| Cell + scaffold held-out | Global mean | 0.089 | 0.787 |
| Cell + scaffold held-out | Nearest drug | 0.052 | 0.663 |

Interpretation:

Top-k DE gene recovery confirms the main benchmark pattern from a response-gene perspective. The cell+drug fingerprint ridge baseline recovers more top response genes in the random split than in scaffold-strict and joint held-out splits. In the strictest joint split, global training mean is higher than the ridge and nearest-drug baselines for top-k overlap, indicating that response-gene recovery also deteriorates when both cell-context reuse and chemical-neighbor shortcuts are removed.

## Low-rank SVD-ridge baseline

Output:

- `results/svd_ridge_multiseed_10/all_svd_ridge_baseline_summary.csv`
- `results/svd_ridge_multiseed_10/svd_ridge_summary_by_seed_mean.csv`
- `figures/openproblems_svd_ridge/svd_ridge_comparison_all_genes.png`
- `figures/openproblems_svd_ridge/svd_ridge_comparison_all_genes.pdf`
- `figures/openproblems_svd_ridge/svd_ridge_comparison_de_genes.png`
- `figures/openproblems_svd_ridge/svd_ridge_comparison_de_genes.pdf`

Method:

- Training response matrices were centered and decomposed by SVD within each training split.
- The first 50 response components were retained.
- Ridge regression mapped cell-context one-hot features and Morgan drug fingerprints to the SVD response coordinates.
- Predictions were reconstructed back to gene space before evaluation.

Key repeated-seed all-gene results:

| Split | Baseline | Mean row-wise Pearson | SD |
|---|---|---:|---:|
| Random | SVD50 + ridge: cell + drug FP | 0.404 | 0.025 |
| Cell held-out | SVD50 + ridge: cell + drug FP | 0.366 | 0.059 |
| Scaffold held-out | SVD50 + ridge: cell + drug FP | 0.237 | 0.029 |
| Cell + scaffold held-out | SVD50 + ridge: cell + drug FP | 0.119 | 0.038 |
| Cell + scaffold held-out | SVD50 + ridge: cell only | 0.168 | 0.030 |

Key repeated-seed DE-gene-only results:

| Split | Baseline | Mean DE-gene row-wise Pearson | SD |
|---|---|---:|---:|
| Random | SVD50 + ridge: cell + drug FP | 0.610 | 0.025 |
| Cell held-out | SVD50 + ridge: cell + drug FP | 0.565 | 0.092 |
| Scaffold held-out | SVD50 + ridge: cell + drug FP | 0.393 | 0.040 |
| Cell + scaffold held-out | SVD50 + ridge: cell + drug FP | 0.210 | 0.062 |
| Cell + scaffold held-out | SVD50 + ridge: cell only | 0.290 | 0.054 |

Interpretation:

The low-rank SVD-ridge baseline improves performance in easier validation settings, especially under random validation and cell-held-out validation. However, the same model still loses much of its advantage under scaffold-held-out and joint cell+scaffold-held-out validation. In the strictest joint split, adding drug fingerprints to the SVD-ridge model performs worse than the cell-only SVD-ridge variant and remains below the global mean pattern observed in the earlier baseline set. This strengthens the benchmark claim: the observed failure is not only a weakness of direct per-gene ridge regression, but a broader failure of chemical-neighbor transfer under strict unseen-drug and unseen-cell-context validation.

## Bootstrap confidence intervals for similarity audit

Output:

- `results/drug_similarity_bootstrap_ci/similarity_metric_seed_bootstrap_ci.csv`
- `results/drug_similarity_bootstrap_ci/similarity_metric_seed_bootstrap_contrasts.csv`
- `figures/openproblems_similarity_ci/nearest_train_drug_tanimoto_bootstrap_ci.png`
- `figures/openproblems_similarity_ci/nearest_train_drug_tanimoto_bootstrap_ci.pdf`

Method:

- Bootstrap confidence intervals were computed over split seeds.
- The primary audited baseline was `ridge_cell_drug_fp`.
- The primary leakage metric was mean maximum Tanimoto similarity from each test drug to the nearest training drug.

Key seed-level bootstrap estimates:

| Split | Mean nearest training-drug Tanimoto | 95% CI |
|---|---:|---:|
| Random | 0.993 | 0.986 to 1.000 |
| Cell held-out | 1.000 | 1.000 to 1.000 |
| Scaffold held-out | 0.279 | 0.271 to 0.288 |
| Cell + scaffold held-out | 0.279 | 0.269 to 0.288 |

Key paired contrasts:

| Contrast | Mean difference | 95% CI |
|---|---:|---:|
| Random minus scaffold held-out | 0.714 | 0.702 to 0.726 |
| Random minus cell+scaffold held-out | 0.714 | 0.701 to 0.728 |

Interpretation:

The bootstrap confidence intervals show that the reduction in nearest training-drug similarity is stable across repeated split seeds. Random validation and cell-held-out validation keep test drugs chemically identical or near-identical to training drugs, whereas scaffold-held-out and joint held-out validation sharply reduce this shortcut. Correlations between similarity and prediction performance are useful as secondary evidence, but are less stable in random validation because nearest-drug similarity is nearly constant there.

## 100-seed statistical hardening

Output:

- `results/openproblems_multiseed_100_de/all_baseline_summary.csv`
- `results/openproblems_multiseed_100_de/baseline_summary_by_seed_mean.csv`
- `results/leakage_overlap_100/leakage_overlap_summary.csv`
- `results/drug_similarity_bootstrap_ci_100/similarity_metric_seed_bootstrap_ci.csv`
- `results/drug_similarity_bootstrap_ci_100/similarity_metric_seed_bootstrap_contrasts.csv`
- `results/ranking_instability_100/ranking_instability_summary.md`
- `results/ridge_sensitivity_30/ridge_cell_drug_fp_sensitivity.csv`
- `docs/high_score_upgrade_results_100seed.md`

The main lightweight OpenProblems analyses were rerun over 100 split seeds to reduce the dependence on a small n=10 seed-level bootstrap.

Key 100-seed DE-gene-only results for `ridge_cell_drug_fp`:

| Split | Mean test records | Mean DE-gene row-wise Pearson | SD |
|---|---:|---:|---:|
| Random | 108.41 | 0.568 | 0.018 |
| Cell held-out | 136.34 | 0.496 | 0.089 |
| Scaffold held-out | 106.70 | 0.392 | 0.037 |
| Cell + scaffold held-out | 26.76 | 0.209 | 0.059 |

Leakage overlap audit:

| Split | Same drug in train | Same scaffold in train | Same cell in train |
|---|---:|---:|---:|
| Random | 0.989 | 0.975 | 1.000 |
| Cell held-out | 1.000 | 0.986 | 0.000 |
| Scaffold held-out | 0.000 | 0.000 | 1.000 |
| Cell + scaffold held-out | 0.000 | 0.000 | 0.000 |

Updated nearest training-drug Tanimoto estimates:

| Split | Mean max Tanimoto | 95% CI |
|---|---:|---:|
| Random | 0.992 | 0.989 to 0.995 |
| Scaffold held-out | 0.278 | 0.275 to 0.281 |
| Cell + scaffold held-out | 0.278 | 0.275 to 0.282 |

Model-ranking instability:

- Under random splitting, `ridge_cell_drug_fp` ranked first in every seed.
- Under scaffold-held-out splitting, the best average-ranked model became `cell_context_mean`, and `ridge_cell_drug_fp` fell to average rank 3.58.
- Under joint cell+scaffold-held-out splitting, `ridge_cell_drug_fp` fell to average rank 5.33.
- Mean Spearman rank correlation was 0.097 for random vs scaffold-held-out and -0.574 for random vs joint held-out.

Interpretation:

The 100-seed analyses strengthen the original manuscript claim. Random validation is not only associated with higher prediction scores; it also retains near-complete drug and scaffold overlap with training data and can select a model ranking that reverses under joint cell+scaffold-held-out validation.

## Ridge hyperparameter sensitivity

Output:

- `scripts/run_ridge_sensitivity.py`
- `results/ridge_sensitivity_30/all_ridge_sensitivity_summary.csv`
- `results/ridge_sensitivity_30/ridge_sensitivity_by_config.csv`
- `results/ridge_sensitivity_30/ridge_cell_drug_fp_sensitivity.csv`

The `ridge_cell_drug_fp` baseline was evaluated over 30 split seeds across alpha values 0.1, 1, 10, 100, and 1000 with 1024-bit Morgan fingerprints. Fingerprint length was also varied at 512, 1024, and 2048 bits with alpha fixed at 10.

Key DE-gene results:

| Alpha | Bits | Random | Scaffold held-out | Cell + scaffold held-out | Random - joint |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 1024 | 0.555 | 0.384 | 0.200 | 0.354 |
| 1 | 1024 | 0.556 | 0.385 | 0.203 | 0.354 |
| 10 | 512 | 0.564 | 0.383 | 0.214 | 0.351 |
| 10 | 1024 | 0.563 | 0.396 | 0.218 | 0.345 |
| 10 | 2048 | 0.563 | 0.409 | 0.227 | 0.336 |
| 100 | 1024 | 0.557 | 0.424 | 0.272 | 0.286 |
| 1000 | 1024 | 0.456 | 0.411 | 0.304 | 0.152 |

Interpretation:

The random-to-strict performance gap is robust across the tested ridge hyperparameters. Stronger regularization reduces the magnitude of the joint held-out gap, but does not remove it. This should be reported as a sensitivity result rather than hidden, because it shows that hyperparameters affect absolute performance while the validation-regime effect remains.

## Sci-Plex 3 second-dataset feasibility

Output:

- `data/raw/scperturb/SrivatsanTrapnell2020_sciplex3.h5ad`
- `metadata/sciplex3_pubchem_smiles.csv`
- `metadata/sciplex3_pubchem_smiles_with_scaffolds.csv`
- `metadata/sciplex3_candidate_records.csv`
- `docs/sciplex3_feasibility_review.md`

Downloaded Sci-Plex 3 data:

- 799,317 cells x 110,983 genes/features
- 188 drug perturbations excluding control
- 3 major cell lines: MCF7, A549, and K562
- 4 dose levels: 10, 100, 1000, and 10000 nM
- 2 time points: 24 h and 72 h
- 17,578 control cells

SMILES and scaffold result:

- PubChem matched 188/188 drug perturbations.
- RDKit parsed 188/188 molecules.
- Bemis-Murcko scaffold count: 165 unique scaffolds.

Candidate record feasibility:

| Minimum cells per treated group | Records | Drugs | Cell lines | Scaffolds |
|---:|---:|---:|---:|---:|
| 10 | 2,444 | 188 | 3 | 165 |
| 20 | 2,440 | 188 | 3 | 165 |
| 50 | 2,423 | 188 | 3 | 165 |
| 100 | 2,377 | 188 | 3 | 165 |

Interpretation:

Sci-Plex 3 is a strong true second-dataset candidate. It is independent of the OpenProblems/OP3 PBMC benchmark and has enough drugs, cell lines, doses, cells per group, and valid scaffolds to support a second leakage-aware benchmark.

## Sci-Plex 3 smoke-test benchmark

Output:

- `scripts/build_sciplex3_pseudobulk.py`
- `scripts/evaluate_sciplex3_response_baselines.py`
- `data/processed/sciplex3/sciplex3_24h_top2000_response.h5ad`
- `metadata/sciplex3_24h_top2000_response_metadata.csv`
- `metadata/sciplex3_24h_top2000_split_metadata.csv`
- `results/sciplex3_24h_top2000_seed_001/splits/candidate_splits.csv`
- `results/sciplex3_24h_top2000_seed_001/baselines/baseline_summary.csv`
- `results/sciplex3_24h_top2000_seed_001/drug_similarity/drug_similarity_vs_error_summary.csv`

Processed response matrix:

- 2,200 records x 2,000 genes
- 188 drugs
- 3 cell lines
- 164 scaffolds represented after 24 h and minimum-cell filtering
- Response: treated mean log1p(CP10K) minus matched cell-line control mean log1p(CP10K)

Seed 1 split sizes:

| Split | Train records | Test records | Excluded records |
|---|---:|---:|---:|
| Random | 1,767 | 433 | 0 |
| Cell-line held-out | 1,463 | 737 | 0 |
| Scaffold held-out | 1,766 | 434 | 0 |
| Cell-line + scaffold held-out | 1,173 | 144 | 883 |

Key smoke-test performance for `ridge_cell_dose_drug_fp`:

| Split | Mean row-wise Pearson |
|---|---:|
| Random | 0.278 |
| Cell-line held-out | 0.149 |
| Scaffold held-out | 0.167 |
| Cell-line + scaffold held-out | 0.061 |

Chemical similarity audit:

| Split | Mean nearest training-drug Tanimoto |
|---|---:|
| Random | 1.000 |
| Cell-line held-out | 1.000 |
| Scaffold held-out | 0.283 |
| Cell-line + scaffold held-out | 0.284 |

Interpretation:

The Sci-Plex 3 smoke test independently reproduces the main leakage pattern. Random validation keeps identical drugs in training and gives the strongest drug-fingerprint ridge performance. Scaffold-strict and joint validation sharply reduce nearest training-drug similarity and reduce performance. This directly addresses the earlier single-dataset critique, although repeated seeds and a numerically stable ridge implementation are still needed before treating the Sci-Plex results as final manuscript evidence.

## Next steps

1. Replace the Sci-Plex ridge solver with a numerically stable implementation.
2. Run repeated Sci-Plex 3 split seeds.
3. Add SVD dimensionality sensitivity for OpenProblems.
4. Fold the 100-seed OpenProblems and Sci-Plex smoke-test findings into the manuscript Results, Discussion, and supplementary tables.
