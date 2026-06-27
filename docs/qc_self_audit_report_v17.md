# v17 QC and self-audit report

## Scope

This report covers the v17 manuscript revision that renames the completed Sci-Plex 3 neural arm as a control-conditioned response regressor, updates figure labels, and adds a Methods-ready architecture specification table. No new performance value was fabricated or manually tuned.

## Commands run

| Check | Command | Result |
|---|---|---|
| Deep-model panel QC | `python3 scripts/qc_deep_model_panel.py` | Passed. The script reported metric toy tests, split integrity checks, output completeness, and no-test-tuning config check as passed. |
| Control-conditioned input shape | `python3 tests/test_transigen_input_shapes.py` | Passed. |
| Control-conditioned no-test-leakage | `python3 tests/test_transigen_no_test_leakage.py` | Passed. |
| Control-conditioned prediction alignment | `python3 tests/test_transigen_prediction_alignment.py` | Passed. |
| Control-conditioned split integrity | `python3 tests/test_transigen_split_integrity.py` | Passed. |
| Control-conditioned ablation shapes | `python3 tests/test_transigen_input_ablation_shapes.py` | Passed. |
| PRnet-interface prediction alignment | `python3 tests/test_prnet_main30_prediction_alignment.py` | Passed. |
| PRnet-interface split integrity | `python3 tests/test_prnet_main30_split_integrity.py` | Passed. |
| Figure regeneration | `envs/transigen_venv/bin/python scripts/summarize_prnet_main30.py` | Passed after switching from system Python, which lacked SciPy, to the project environment. |

## Split and leakage assertions

| Model output | Split | Same drug max | Same scaffold max | Same cell-line max | Status |
|---|---|---:|---:|---:|---|
| Control-conditioned response regressor | scaffold held-out | 0.0 | 0.0 | 1.0 | Pass |
| Control-conditioned response regressor | joint cell-line+scaffold held-out | 0.0 | 0.0 | 0.0 | Pass |
| PRnet-interface adapter | scaffold held-out | 0.0 | 0.0 | 1.0 | Pass |
| PRnet-interface adapter | joint cell-line+scaffold held-out | 0.0 | 0.0 | 0.0 | Pass |

## Prediction-vector checks

| Manifest | Unique NPZ files checked | Required vector length | Failed files | Status |
|---|---:|---:|---:|---|
| `results/deep_model_panel/sciplex3/transigen_main30_predictions_manifest.csv` | 120 | 2,000 genes | 0 | Pass |
| `results/deep_model_panel/sciplex3/prnet_main30_predictions_manifest.csv` | 120 | 2,000 genes | 0 | Pass |

All checked NPZ files had matching `prediction` and `target` shapes, and prediction vectors had 2,000 genes.

## No-test-tuning audit

- The control-conditioned runner builds `fit_idx` and `val_idx` only from records assigned to `train` in the fixed split manifest.
- Feature and target standardization are fitted only on the fitting records in `TranSiGenAdapter.fit`.
- Validation records are transformed with training-fitted scalers and used only for early stopping.
- Test records are passed only to `predict` after fitting and are used only for final evaluation.
- Excluded records are not included in `train_idx` or `test_idx`.

## Items not fully machine-verifiable from saved artifacts

- Historical per-epoch train/validation loss curves were not saved during the 30-seed neural runs.
- The current fingerprint NPZ records a 1,024-bit Morgan matrix but does not store the radius as metadata.
- Sci-Plex 3 DE-gene masks are not present in the current pseudobulk artifact, so DE-gene-only neural metrics were not computed.

## Manuscript integrity notes

- The manuscript now describes the completed neural arm by architecture, not as a third-party model reproduction.
- The internal file identifier `transigen_adapted_sciplex3` is retained only for provenance.
- The 10-seed fingerprint/input ablation is described as diagnostic and is not presented as a 30-seed main result.
