# BMC Master Result Provenance

This directory is the single numerical source for the BMC Bioinformatics revision. It was generated from immutable fixed-manifest, completed seed-level outputs only; no model training occurs in this build step.

## Added controls

- OpenProblems: 100 fixed seeds for test-size-matched, train-and-test-size-matched, and cell-context plus DE-gene-count-quartile-matched random controls of the cell-plus-drug-fingerprint ridge.
- Sci-Plex 3: 30 fixed seeds for test-size-matched and train-and-test-size-matched random controls of the control-conditioned nonlinear benchmark. Training uses only sampled random training records; validation is a deterministic partition of that training set.
- Joint manifests supply target counts only for these controls. Their test, training, and excluded record IDs are never reused by a random-control split.

## Primary source files

- `results/tcbb_master/`: pre-existing fixed-manifest baseline and audit source tables.
- `results/openproblems_train_size_matched_100/`: new OpenProblems controls and reconstructed audit fields.
- `results/deep_model_panel/sciplex3/control_conditioned_nonlinear_train_size_matched_30/`: new Sci-Plex nonlinear controls, prediction manifests, audit fields, and training-loss traces.

All contrasts are paired by the same manifest seed and use seed-level bootstrap resampling (10,000 resamples).
