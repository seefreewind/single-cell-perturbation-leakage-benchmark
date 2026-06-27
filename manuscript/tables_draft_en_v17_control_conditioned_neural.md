# Tables for v17 control-conditioned neural manuscript

## Table 1. Dataset and benchmark coverage

| Dataset | Records | Drugs | Cell contexts | Scaffold annotation | Main role |
|---|---:|---:|---:|---|---|
| OpenProblems NeurIPS 2023 | 545 treated records | 140 | 4 | 136 valid Bemis-Murcko scaffolds | Primary perturbation transcriptomics benchmark |
| Sci-Plex 3 24 h pseudobulk | 2,200 drug-cell-line-dose records | 188 mapped drugs | 3 cell lines | 165 valid Bemis-Murcko scaffolds | External chemical-transcriptomics validation |

## Table 2. Leakage-aware split definitions

| Split | Held-out unit | Training-test shortcut removed | Interpretation |
|---|---|---|---|
| Random | Record | None by design | Local interpolation benchmark |
| Cell held-out | Cell context or cell line | Same cell context | Cellular-context extrapolation |
| Scaffold held-out | Bemis-Murcko scaffold | Same drug and same scaffold | Chemical-scaffold extrapolation |
| Joint cell+scaffold held-out | Cell context plus scaffold | Same drug, scaffold, and cell context | Joint cellular and chemical extrapolation |
| MoA-held-out | Mechanism annotation | Same MoA | Mechanism-level secondary stress test |

## Table 3. Completed model performance summary

| Dataset | Split | Mean test records | Model family | Seeds | Main result |
|---|---|---:|---|---:|---|
| OpenProblems | Random | 108.4 | ridge cell+drug FP | 100 | all-gene Pearson 0.362; DE-gene Pearson 0.568 |
| OpenProblems | Scaffold held-out | 106.7 | ridge cell+drug FP | 100 | all-gene Pearson 0.234; DE-gene Pearson 0.392 |
| OpenProblems | Joint held-out | 26.8 | ridge cell+drug FP | 100 | all-gene Pearson 0.118; DE-gene Pearson 0.209 |
| OpenProblems | MoA-held-out | 110.8 | ridge cell+drug FP | 100 | all-gene Pearson 0.212; DE-gene Pearson 0.363 |
| Sci-Plex 3 | Random | 437.5 | ridge cell+dose+drug FP | 30 | all-gene Pearson 0.292 |
| Sci-Plex 3 | Scaffold held-out | 427.1 | ridge cell+dose+drug FP | 30 | all-gene Pearson 0.189 |
| Sci-Plex 3 | Joint held-out | 141.9 | ridge cell+dose+drug FP | 30 | all-gene Pearson 0.060 |
| Sci-Plex 3 | Random | 437.5 | Control-conditioned response regressor | 30 | all-gene Pearson 0.233; RMSE 0.080 |
| Sci-Plex 3 | Scaffold held-out | 427.1 | Control-conditioned response regressor | 30 | all-gene Pearson 0.212; RMSE 0.084 |
| Sci-Plex 3 | Joint held-out | 141.9 | Control-conditioned response regressor | 30 | all-gene Pearson 0.017; RMSE 0.581 |
| Sci-Plex 3 | Random | 437.5 | PRnet-interface pseudobulk adapter | 30 | all-gene Pearson 0.017; bootstrap 95% CI 0.013 to 0.022 |
| Sci-Plex 3 | Scaffold held-out | 427.1 | PRnet-interface pseudobulk adapter | 30 | all-gene Pearson 0.015; bootstrap 95% CI 0.012 to 0.018 |
| Sci-Plex 3 | Joint held-out | 141.9 | PRnet-interface pseudobulk adapter | 30 | all-gene Pearson 0.001; bootstrap 95% CI -0.013 to 0.014 |

## Table 4. Modern model feasibility and adaptation status

| Model | Status | Sci-Plex 3 use | Input contract | Adaptation level | Main limitation |
|---|---|---|---|---|---|
| Control-conditioned response regressor | Completed | Repeated-seed benchmark | Basal, drug FP, dose, and cell-line features | Local pseudobulk neural architecture | Not a full original TranSiGen reproduction; internal file ID remains transigen_adapted_sciplex3 |
| PRnet-interface pseudobulk adapter | Completed | Repeated-seed benchmark | Basal plus dose-scaled drug FP | Pseudobulk adapter using official PGM module | Not a full original PRnet reproduction |
| chemCPA/CPA | Not run | Not used | Raw cells with perturbation covariates | Adapter/status row | Requires raw single-cell workflow |
| scGPT frozen ridge/MLP | Not run | Not used | Frozen embeddings plus response head | Adapter/status row | Embeddings not generated |
| scFoundation frozen ridge/MLP | Not run | Not used | Frozen embeddings plus response head | Adapter/status row | Embeddings not generated |

## Table 5. Random-to-strict model ranking transfer

| Dataset | Contrast | Model count | Aggregate or mean Spearman | Random-best model | Random-best rank in joint |
|---|---|---:|---:|---|---:|
| OpenProblems | random vs scaffold held-out | 11 | 0.424 | SVD50 ridge cell+drug FP | 4 |
| OpenProblems | random vs joint held-out | 11 | 0.060 | SVD50 ridge cell+drug FP | 6 |
| Sci-Plex 3 | random vs scaffold held-out | 5 | 0.530 | seed-level variable | NA |
| Sci-Plex 3 | random vs joint held-out | 5 | 0.143 | seed-level variable | 2.43 |

## Table 6. Minimum reporting checklist for leakage-aware perturbation transcriptomics benchmarks

| Reporting item | Required output |
|---|---|
| Split type | Random, cell held-out, scaffold held-out, joint held-out, and mechanism-aware split where feasible |
| Leakage overlap | Same drug, scaffold, MoA, cell context, and dose overlap |
| Chemical similarity | Nearest training-drug Tanimoto |
| Mechanism-neighbor audit | Same-MoA or target-class overlap where annotations are available |
| Model adaptation level | Original workflow, benchmark-compatible adaptation, pseudobulk adaptation, frozen embedding head, failed, or not run |
| Performance | Pearson, RMSE, DE-gene metrics, top-k overlap, direction agreement |
| Model ranking | Random-to-strict rank transfer |
| Failure reporting | Explicit status rows for not-run or failed models |
| External validation | At least one independent dataset when available |

## Table 7. Control-conditioned response regressor specification

| Component | Specification |
|---|---|
| Dataset and splits | Sci-Plex 3 24 h pseudobulk top-2,000-gene benchmark; fixed manifest seeds 1-30; random, cell-line held-out, scaffold held-out, and joint cell-line-plus-scaffold held-out splits |
| Target | Treated-minus-control response vector with 2,000 genes |
| Inputs | Matched 24 h basal/control expression vector (2,000), Morgan fingerprint (1,024), log10(dose + 1), and cell-line one-hot vector (3) |
| Morgan fingerprint | Existing local 1,024-bit Morgan fingerprint matrix from `compound_fingerprints.npz`; radius not recorded in the current output metadata |
| Input dimension | 3,028 for the full model |
| Network | Fully connected MLP: Linear(input, 128), ReLU, Dropout(0.1), Linear(128, 128), ReLU, Dropout(0.1), Linear(128, 2,000) |
| Normalization | Feature and target standardization fitted on fitting records only; validation and test records transformed with training-fitted scalers |
| Validation | 10% of training records selected internally by split seed plus 1729; no test records used for validation or early stopping |
| Loss and optimizer | Mean-squared error loss; AdamW optimizer; learning rate 1e-3; weight decay 1e-5 |
| Training schedule | Maximum 12 epochs for main30; batch size 128; early-stopping patience 5 |
| Seed policy | Model initialization seed equals fixed split seed; seeds 1-30 only |
| Adaptation note | Local auditable architecture inspired by basal-plus-compound perturbation modeling; not a line-by-line reproduction of a published TranSiGen workflow |
