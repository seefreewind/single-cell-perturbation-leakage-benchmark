# PRnet Smoke-to-Main Audit

Date: 2026-06-24

## Purpose

This audit records the changes needed to upgrade `prnet_adapted_sciplex3` from the earlier two-epoch smoke test to a repeated-seed Sci-Plex 3 benchmark. The main correction is that all main10 runs must use real split-manifest seeds 1-10, with the model initialization seed equal to the split seed.

## Repository and Provenance

- Official repository: `https://github.com/Perturbation-Response-Prediction/PRnet`
- Recorded upstream commit: `f19174bde3ed2633f54c7831799cc38c4ffc7a0d`
- Local path: `external/PRnet/`
- Local checkout type: GitHub source archive, because direct git clone failed during pack transfer.
- License: Apache-2.0.
- Third-party source modifications: none.

## Official PRnet Component Used

The local adapter imports the official PRnet `PGM` module from:

- `external/PRnet/models/PRnet.py`

The adapter does not edit official PRnet code. It wraps the official PGM module in a Sci-Plex 3 pseudobulk benchmark interface.

## Current Adapter Input Contract

The main10 benchmark uses `prnet_adapted_sciplex3` with:

- basal/control expression: 2,000 genes;
- treated expression: 2,000 genes, used as training target;
- response target for evaluation: treated minus basal/control expression;
- dose-scaled drug fingerprint: 1,024-dimensional Morgan fingerprint multiplied by `log10(dose + 1)`;
- cell-line features: not used by the current PRnet adapter.

The adapter predicts treated expression and then converts it to a response vector by subtracting the matched basal/control vector. The evaluated output is therefore a 2,000-gene response vector aligned with the existing Sci-Plex 3 pipeline.

## Metric Compatibility

The smoke test already computed the same metric family used by the Sci-Plex 3 benchmark:

- all-gene row-wise Pearson;
- all-gene RMSE;
- top-50, top-100, and top-200 overlap;
- top-k direction agreement.

The main10 runner keeps the same metric functions from `src.metrics.vector_metrics`.

## Smoke-Test Seed Issue

The earlier smoke test used:

- `requested_seed=0` for model initialization;
- `split_seed=1` because the fixed split manifest contains seeds 1-30 and does not contain seed 0.

This is acceptable for a smoke test but not for a repeated-seed benchmark. The main10 benchmark corrects this by using real split-manifest seeds 1-10 only. For every main10 combination:

- `seed` equals the split manifest seed;
- the model initialization seed equals `seed`;
- no `split_seed` substitution is used;
- output rows report only real manifest seed values.

## No-Leakage Rules for Main10

The main10 runner must preserve the fixed split assignments:

- no split assignment is modified for PRnet;
- validation records are sampled only from the training records;
- test records are used only for final evaluation;
- excluded records are not used in train, validation, or test;
- no test records are used for normalization, early stopping, or hyperparameter selection.

The current PRnet adapter performs no feature normalization. This avoids a normalization leakage channel. The validation set is used only for early stopping within each training split.

## Main10 Decision

Proceed to `prnet_adapted_sciplex3` main10 with:

- dataset: `sciplex3_24h_top2000`;
- seeds: 1-10;
- splits: `random`, `scaffold_heldout`, and `cell_scaffold_heldout`;
- epochs: 20;
- early stopping: true;
- patience: 5;
- batch size: 128;
- device: CPU unless CUDA is detected by PyTorch.
