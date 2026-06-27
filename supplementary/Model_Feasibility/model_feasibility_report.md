# Model Feasibility Report for the v8 Deep-Model Panel

Date: 2026-06-24

This report records whether each planned modern perturbation model or frozen foundation-model representation was runnable in the current local benchmark environment. It is part of the benchmark resource and should be updated whenever a new adapter is installed, configured, or run.

## Environment Check

Local imports were tested in the project virtual environment. `anndata`, `RDKit`, and `scikit-learn` were available. The following model or training dependencies were not available locally: `torch`, `prnet`, `transigen`, `chemCPA`, `cpa`, `scgpt`, `scfoundation`, and `geneformer`.

## Feasibility Table

| Model | Dependency status | Intended input | Uses drug structure | Uses dose | Raw single-cell required | Pseudobulk/DGE adaptation | OpenProblems status | Sci-Plex 3 status | Result files | Included in main performance results | Main limitation |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PRnet | Not installed | Raw or model-specific perturbation data with chemical information | Yes | Pending implementation | Likely yes for original workflow | Not completed | Not run | Not run | Status rows in `results/deep_model_panel/all_model_metrics_long.csv` | No | External implementation and required input contract are not configured locally. |
| TranSiGen | Repository cloned; CPU adaptation completed | Basal expression plus molecular representation | Yes | Yes | Not necessarily | Completed as Sci-Plex 3 pseudobulk response adaptation | Not run | Completed | `results/deep_model_panel/sciplex3/transigen_metrics_long.csv` | Yes, Sci-Plex 3 only | This is a pseudobulk `transigen_adapted_sciplex3` stress test, not a full reproduction of the original TranSiGen training setting. |
| chemCPA/CPA | Not installed | Raw single-cell counts, perturbation labels, dose, covariates | Yes for chemCPA | Yes | Yes | Matrix-level adaptation would not be equivalent to the original model | Not run | Not run | Status rows in `results/deep_model_panel/all_model_metrics_long.csv` | No | Neither chemCPA nor CPA is installed; original workflows require raw single-cell training data and model-specific preprocessing. |
| scGPT frozen ridge | Not installed | Precomputed scGPT embeddings plus response head inputs | Optional via concatenated drug fingerprint | Yes when available | Raw single-cell or precomputed embeddings | Feasible after embedding extraction | Not run | Not run | Status rows in `results/deep_model_panel/all_model_metrics_long.csv` | No | scGPT package and pretrained embedding workflow are not configured. |
| scGPT frozen MLP | Not installed | Precomputed scGPT embeddings plus MLP response head | Optional via concatenated drug fingerprint | Yes when available | Raw single-cell or precomputed embeddings | Feasible after embedding extraction | Not run | Not run | Status rows in `results/deep_model_panel/all_model_metrics_long.csv` | No | Same as scGPT frozen ridge; `torch` is also absent. |
| scFoundation frozen ridge | Not installed | Precomputed scFoundation embeddings plus response head inputs | Optional via concatenated drug fingerprint | Yes when available | Raw single-cell or precomputed embeddings | Feasible after embedding extraction | Not run | Not run | Status rows in `results/deep_model_panel/all_model_metrics_long.csv` | No | scFoundation package and pretrained embedding workflow are not configured. |
| scFoundation frozen MLP | Not installed | Precomputed scFoundation embeddings plus MLP response head | Optional via concatenated drug fingerprint | Yes when available | Raw single-cell or precomputed embeddings | Feasible after embedding extraction | Not run | Not run | Status rows in `results/deep_model_panel/all_model_metrics_long.csv` | No | Same as scFoundation frozen ridge; `torch` is also absent. |

## Completed Model Families

The panel currently includes completed existing results for mean baselines, ridge baselines, nearest-drug retrieval baselines, SVD-ridge stress tests, OpenProblems MoA-held-out baselines, Sci-Plex 3 pseudobulk baselines, and a completed `transigen_adapted_sciplex3` 30-seed Sci-Plex 3 stress test.

## Required Next Steps Before Claiming Deep-Model Benchmark Results

1. Install and pin a `torch` runtime compatible with the target models.
2. Add external model repositories through a documented script without committing large weights.
3. Define model-specific input builders that reuse the existing split manifests without changing test assignments.
4. Run one-seed smoke tests, then 10-seed pilots, then 30-seed main panels.
5. Update this report and the long-format result CSVs with completed metrics or explicit failure rows.
