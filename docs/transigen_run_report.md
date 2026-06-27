# TranSiGen Sci-Plex 3 Run Report

Date: 2026-06-24

## Installation and Repository

- TranSiGen repository: `https://github.com/myzheng-SIMM/TranSiGen.git`
- Local path: `external/TranSiGen/`
- Commit: `8ec2218e2fe4fbb5f3a2c14271a8640a62c1af3a`
- License: MIT
- Third-party source modifications: none

## Runtime

- Execution environment: `envs/transigen_venv/`
- Python: 3.9 virtual environment
- Torch: 2.8.0
- CUDA/GPU: no NVIDIA GPU detected; run used CPU.
- Install check: passed with `scripts/check_transigen_install.py`.

## Adaptation Used

The completed model is named `transigen_adapted_sciplex3`. It is a Sci-Plex 3 pseudobulk adaptation that uses:

- matched basal/control expression profile;
- Morgan drug fingerprint;
- log-transformed dose;
- cell-line one-hot encoding;
- MLP response head predicting 2,000-gene treated-control response vectors.

This is not a full reproduction of the original TranSiGen training setting. It is a benchmark-compatible adaptation designed to test whether a TranSiGen-style basal-expression plus chemical-representation model can be evaluated under the fixed leakage-aware split manifests.

## Input Matrix

- Input package: `data/processed/transigen/sciplex3/`
- Response records: 2,200
- Genes: 2,000
- Basal/control matrix: `(2200, 2000)`
- Treated matrix: `(2200, 2000)`
- Response target matrix: `(2200, 2000)`
- Split assignment file: `data/processed/transigen/sciplex3/split_assignments_long.csv`

## Training Runs

- Smoke test: seed 1, random split, completed.
- Pilot: seeds 1-3, random/scaffold/joint splits, completed.
- Main: seeds 1-30, random/cell-line/scaffold/joint splits, completed.
- Failed main combinations: 0.

## Core Performance

| Split | Seeds | Mean all-gene Pearson | SD | Mean RMSE | SD |
|---|---:|---:|---:|---:|---:|
| Random | 30 | 0.233 | 0.013 | 0.080 | 0.002 |
| Cell-line held-out | 30 | 0.027 | 0.021 | 0.663 | 0.602 |
| Scaffold held-out | 30 | 0.212 | 0.020 | 0.084 | 0.006 |
| Joint cell-line+scaffold held-out | 30 | 0.017 | 0.017 | 0.581 | 0.436 |

## Random-to-Strict Drops

- Random minus scaffold held-out Pearson difference: 0.021.
- Random minus joint cell-line+scaffold held-out Pearson difference: 0.216.

## Comparison With Ridge Baseline

The existing Sci-Plex 3 `ridge_cell_dose_drug_fp` baseline had random, scaffold held-out, and joint held-out mean Pearson values of approximately 0.292, 0.189, and 0.060. The TranSiGen adaptation was lower in random and joint held-out settings, but slightly higher than this ridge baseline in scaffold held-out validation. Both models showed a strong drop under joint held-out validation.

## Output Files

- `results/deep_model_panel/sciplex3/transigen_metrics_long.csv`
- `results/deep_model_panel/sciplex3/transigen_predictions_manifest.csv`
- `results/deep_model_panel/sciplex3/transigen_leakage_audit.csv`
- `results/deep_model_panel/sciplex3/transigen_random_to_strict_contrasts.csv`
- `figures/transigen_pilot_sciplex3_performance.png`
- `figures/figure6_modern_model_panel.png`
- `figures/figure7_rank_transfer_with_transigen.png`

## QC Status

The following checks passed:

- input matrix shape;
- prediction-target shape and record alignment;
- no train/test/excluded overlap;
- scaffold and joint split integrity;
- no-test-tuning configuration;
- deep-model-panel output completeness.

## Recommendation

TranSiGen should enter the manuscript as a completed Sci-Plex 3-only stress test, with careful wording. The result supports the main benchmark claim that modern chemical-response models can be integrated into the same leakage-aware split/audit framework, and that strict joint held-out validation remains substantially harder.

## Next Model Recommendation

The next model should be chemCPA/CPA if the goal is to evaluate a canonical single-cell perturbation model on raw Sci-Plex 3 data. PRnet is also relevant, but chemCPA/CPA is more directly tied to single-cell perturbation benchmarking and dose-aware covariates.

## Human Inputs Still Needed

- Confirm whether the pseudobulk adaptation should be described as `TranSiGen-style` or `TranSiGen-adapted` in the final manuscript.
- Decide whether to run chemCPA/CPA next.
- Confirm final public GitHub URL and Zenodo DOI once the resource package is archived.
