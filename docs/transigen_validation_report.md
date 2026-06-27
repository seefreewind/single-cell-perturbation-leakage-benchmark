# TranSiGen-Style Sci-Plex 3 Validation Report

Date: 2026-06-24

## Scope

This report validates the completed `transigen_adapted_sciplex3` analysis before any further manuscript expansion. The completed model should be described as a TranSiGen-style pseudobulk neural adaptation, not as a full reproduction of the original TranSiGen paper. The official third-party repository was cloned for provenance and inspection, but the benchmark run used a local auditable adapter because the official environment is old, Linux/CUDA-oriented, and not directly compatible with the present macOS CPU runtime.

## Third-Party Repository Provenance

- Official repository: `https://github.com/myzheng-SIMM/TranSiGen.git`
- Local path: `external/TranSiGen/`
- Commit: `8ec2218e2fe4fbb5f3a2c14271a8640a62c1af3a`
- License: MIT
- Source modifications: none

## Runtime

- Environment: `envs/transigen_venv/`
- Torch: 2.8.0
- CUDA available: false
- NumPy: 2.0.2
- pandas: 2.3.3
- scikit-learn: 1.6.1
- SciPy: 1.13.1
- Execution mode: CPU

## Input Data

- Dataset: Sci-Plex 3, 24 h subset, top 2,000 genes.
- Input directory: `data/processed/transigen/sciplex3/`
- Response records: 2,200 drug-cell-line-dose pseudobulk records.
- Basal/control matrix: 2,200 x 2,000.
- Treated matrix: 2,200 x 2,000.
- Response target: treated mean log1p(CP10K) minus matched control mean log1p(CP10K).
- Drug representation: Morgan fingerprint, 1,024 dimensions.
- Dose representation: log10(dose + 1).
- Cell-context representation: one-hot cell-line encoding.
- Split manifest: `data/processed/transigen/sciplex3/split_assignments_long.csv`.

## Model Adapter

The adapter is implemented in `src/models/deep/transigen_adapter.py`. It uses a two-hidden-layer MLP response head with ReLU activations and dropout to predict the 2,000-gene response vector from benchmark-compatible features.

The full current input vector has 3,028 dimensions:

- 2,000 basal/control expression features;
- 1,024 Morgan fingerprint features;
- 1 dose feature;
- 3 cell-line one-hot features.

The model standardizes training features and targets within each split. Test labels are not used during feature or target scaling. Each training split is further divided into an internal 90/10 train-validation split for early stopping.

## Training Protocol

- Smoke run: seed 1, random split, 2 epochs.
- Pilot run: seeds 1-3, random/scaffold/joint splits, 5 epochs.
- Main run: seeds 1-30, random/cell-line/scaffold/joint splits, 10-12 epoch budget depending on run script defaults.
- Ablation run: seeds 1-10, all four split types, 10 epochs.
- Optimizer: AdamW.
- Learning rate: 1e-3.
- Batch size: 128 or 256 depending on run script.
- Early stopping: internal validation only, patience `max(2, min(5, epochs // 3 + 1))`.

## Main 30-Seed Performance

| Split | Seeds | Mean Pearson | SD | Median Pearson | IQR | Mean RMSE | Mean top-200 overlap |
|---|---:|---:|---:|---:|---:|---:|---:|
| Random | 30 | 0.2329 | 0.0128 | 0.2346 | 0.0215 | 0.0804 | 0.2465 |
| Cell-line held-out | 30 | 0.0270 | 0.0210 | 0.0279 | 0.0289 | 0.6630 | 0.1377 |
| Scaffold held-out | 30 | 0.2118 | 0.0198 | 0.2099 | 0.0266 | 0.0835 | 0.2395 |
| Joint cell-line+scaffold held-out | 30 | 0.0169 | 0.0170 | 0.0181 | 0.0213 | 0.5808 | 0.1329 |

## Paired Random-to-Strict Contrasts

| Contrast | Seeds | Mean Pearson difference | Bootstrap 95% CI |
|---|---:|---:|---:|
| Random minus scaffold held-out | 30 | 0.0211 | 0.0131 to 0.0289 |
| Random minus joint held-out | 30 | 0.2159 | 0.2083 to 0.2239 |
| Scaffold held-out minus joint held-out | 30 | 0.1948 | 0.1859 to 0.2043 |

These paired contrasts support the main benchmark point: the TranSiGen-style adaptation retains moderate performance when chemical scaffold is held out but collapses when both cell-line context and scaffold are held out.

## Comparison With Sci-Plex 3 Baselines

Paired comparisons were calculated against existing Sci-Plex 3 baseline tables using matched seeds and split types.

| Baseline | Split | Mean TranSiGen-style Pearson | Mean baseline Pearson | Difference | Bootstrap 95% CI |
|---|---|---:|---:|---:|---:|
| `ridge_cell_dose_drug_fp` | Random | 0.2329 | 0.2916 | -0.0587 | -0.0617 to -0.0557 |
| `ridge_cell_dose_drug_fp` | Scaffold held-out | 0.2118 | 0.1885 | 0.0233 | 0.0194 to 0.0272 |
| `ridge_cell_dose_drug_fp` | Joint held-out | 0.0169 | 0.0597 | -0.0427 | -0.0537 to -0.0319 |
| `ridge_drug_fp` | Random | 0.2329 | 0.2431 | -0.0103 | -0.0135 to -0.0069 |
| `ridge_drug_fp` | Scaffold held-out | 0.2118 | 0.1264 | 0.0854 | 0.0809 to 0.0901 |
| `ridge_drug_fp` | Joint held-out | 0.0169 | 0.0582 | -0.0413 | -0.0518 to -0.0308 |

The adaptation therefore should not be presented as uniformly stronger than ridge. Its clearest favorable comparison is under scaffold-held-out validation, while ridge remains stronger under random and joint held-out settings.

## Ten-Seed Input Ablation

| Model | Input dimension | Random | Cell-line | Scaffold | Joint |
|---|---:|---:|---:|---:|---:|
| Basal only | 2,000 | 0.2229 | 0.0030 | 0.2297 | -0.0019 |
| Drug only | 1,024 | 0.2323 | 0.1038 | 0.1656 | 0.0769 |
| Basal+drug | 3,024 | 0.2353 | 0.0267 | 0.2149 | 0.0160 |
| Basal+drug+dose | 3,025 | 0.2351 | 0.0226 | 0.2129 | 0.0134 |
| Basal+drug+dose+cell | 3,028 | 0.2347 | 0.0307 | 0.2152 | 0.0130 |
| Full current | 3,028 | 0.2347 | 0.0307 | 0.2152 | 0.0130 |

The ablation indicates that the full input set does not rescue joint held-out generalization. The drug-only model has the highest joint held-out mean among the ablations, but it is weaker on scaffold-held-out validation and is not the main architecture used for the completed TranSiGen-style analysis.

## Output Files

- Main metrics: `results/deep_model_panel/sciplex3/transigen_metrics_long.csv`
- Main prediction manifest: `results/deep_model_panel/sciplex3/transigen_predictions_manifest.csv`
- Leakage audit: `results/deep_model_panel/sciplex3/transigen_leakage_audit.csv`
- Robustness summary: `results/deep_model_panel/sciplex3/transigen_summary_stats.csv`
- Paired strict-split contrasts: `results/deep_model_panel/sciplex3/transigen_paired_contrasts.csv`
- Baseline paired differences: `results/deep_model_panel/sciplex3/transigen_vs_baseline_paired_differences.csv`
- Ablation metrics: `results/deep_model_panel/sciplex3/transigen_ablation_metrics_long.csv`
- Ablation contrasts: `results/deep_model_panel/sciplex3/transigen_ablation_contrasts.csv`
- Ablation input dimensions: `results/deep_model_panel/sciplex3/transigen_ablation_input_dims.csv`
- Figure source data: `results/deep_model_panel/source_data/figure6_transigen_ablation_sciplex3.csv`
- Additional source data: `results/deep_model_panel/source_data/supp_transigen_seed_distribution.csv`, `results/deep_model_panel/source_data/supp_transigen_random_to_strict_ci.csv`, `results/deep_model_panel/source_data/supp_transigen_vs_ridge_paired.csv`

## Figures

- `figures/figure6_transigen_ablation_sciplex3.png`
- `figures/figure6_transigen_ablation_sciplex3.pdf`
- `figures/supp_transigen_seed_distribution.png`
- `figures/supp_transigen_seed_distribution.pdf`
- `figures/supp_transigen_random_to_strict_ci.png`
- `figures/supp_transigen_random_to_strict_ci.pdf`
- `figures/supp_transigen_vs_ridge_paired.png`
- `figures/supp_transigen_vs_ridge_paired.pdf`

## Validation Tests

The following checks were run successfully with `envs/transigen_venv/bin/python`:

- `tests/test_transigen_input_shapes.py`
- `tests/test_transigen_no_test_leakage.py`
- `tests/test_transigen_prediction_alignment.py`
- `tests/test_transigen_split_integrity.py`
- `tests/test_transigen_input_ablation_shapes.py`

## Manuscript Wording Recommendation

Use: "TranSiGen-style pseudobulk neural adaptation" or "TranSiGen-inspired benchmark adapter."

Avoid: "full TranSiGen reproduction", "official TranSiGen performance", or any claim that the result reproduces the original TranSiGen training recipe.

The manuscript can use this result to support a resource-oriented claim: a modern chemical-response neural architecture can be plugged into the leakage-aware Sci-Plex 3 split manifests, and its apparent random-split performance does not transfer to joint cell-line and scaffold extrapolation.
