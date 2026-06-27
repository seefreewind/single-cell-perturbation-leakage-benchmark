# PRnet Sci-Plex 3 Smoke Report

Date: 2026-06-24

## Scope

This report documents a minimal PRnet-style Sci-Plex 3 smoke test. The goal was only to verify that the model code can read the benchmark inputs, fit, predict, emit 2,000-gene response vectors, and compute metrics. The result should not be treated as a stable PRnet benchmark estimate.

## Model Name and Adaptation Level

- Model name: `prnet_adapted_sciplex3`
- Adaptation level: PRnet-style pseudobulk adapter.
- Official source used: PRnet `PGM` module from `external/PRnet/models/PRnet.py`.
- Official source modifications: none.
- Full PRnet reproduction: no.

## Input and Split

- Dataset: Sci-Plex 3, 24 h, top 2,000 genes.
- Requested seed: 0.
- Effective split seed: 1, because the existing split manifest contains seeds 1-30 and does not include seed 0.
- Split type: random.
- Training records: 1,767.
- Internal validation records: 177.
- Test records: 433.
- Test records were not used for validation or early stopping.

## Training Configuration

- Epochs: 2.
- Batch size: 128.
- Learning rate: 1e-3.
- Device: CPU.
- Input control expression dimension: 2,000.
- Input drug-dose dimension: 1,024.
- Output response vector length: 2,000 genes.

## Smoke-Test Output

| Metric | Value |
|---|---:|
| Status | completed-smoke |
| Mean all-gene row-wise Pearson | 0.0014 |
| Mean all-gene RMSE | 0.5203 |
| Mean top-50 overlap | 0.1578 |
| Mean top-100 overlap | 0.1868 |
| Mean top-200 overlap | 0.2311 |
| Mean top-200 direction agreement | 0.5065 |

The near-zero Pearson is not interpreted as PRnet performance because this was a two-epoch smoke test on a pseudobulk adapter. It only confirms that the pipeline can complete fit, predict, and metric computation.

## Output Files

- Metrics: `results/deep_model_panel/sciplex3/prnet_smoke_metrics.csv`
- Prediction manifest: `results/deep_model_panel/sciplex3/prnet_smoke_predictions_manifest.csv`
- Prediction array: `results/deep_model_panel/sciplex3/predictions/prnet_smoke_seed0_random.npz`
- Log: `logs/prnet_smoke_seed0.log`
- Registry metrics with retained placeholder rows: `results/deep_model_panel/sciplex3/prnet_metrics.csv`
- Registry manifest with retained placeholder rows: `results/deep_model_panel/sciplex3/prnet_predictions_manifest.csv`

## Shape Checks

- Prediction array shape: 433 x 2,000.
- Target array shape: 433 x 2,000.
- Record IDs: 433.
- Manifest rows: 433.

## Recommendation

PRnet can be kept in the resource-readiness section as a completed smoke test, not as a main model-comparison result. A 10-seed or 30-seed PRnet main run should be considered only after deciding whether to continue with this pseudobulk adapter or to build a fuller AnnData bridge that more closely follows PRnet's official single-cell input contract.
