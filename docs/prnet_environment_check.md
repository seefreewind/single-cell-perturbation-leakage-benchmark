# PRnet Environment Check

Date: 2026-06-24

## Repository Access

- Official repository: `https://github.com/Perturbation-Response-Prediction/PRnet`
- Access status: accessible.
- Initial `git clone` and shallow clone attempts failed because the network connection ended during pack transfer.
- Fallback archive download succeeded from GitHub codeload.
- Local path: `external/PRnet/`
- Recorded upstream commit: `f19174bde3ed2633f54c7831799cc38c4ffc7a0d`
- Local checkout type: downloaded source archive, not a live git clone.
- Source modifications: none.

## License

- License file: `external/PRnet/LICENSE`
- License: Apache-2.0.

## Requirements

Official `requirements.txt` lists:

- numpy
- pandas
- tqdm
- scikit-learn
- torch
- scanpy
- rdkit
- anndata
- matplotlib
- seaborn

## Local Runtime

- Environment: `envs/transigen_venv/`
- Torch: 2.8.0
- CUDA available: false
- scanpy: 1.10.3
- anndata: 0.10.9
- RDKit: 2025.09.2
- scikit-learn: 1.6.1
- NumPy: 2.0.2
- pandas: 2.3.3
- Execution mode: CPU.

`scanpy` and its dependencies were installed into `envs/transigen_venv/` because the official PRnet demo imports scanpy directly.

## GPU Requirement

The official README recommends a Docker image based on PyTorch CUDA, and the official test scripts set `CUDA_VISIBLE_DEVICES=0`. However, the trainer chooses `cuda` only if `torch.cuda.is_available()` is true. Therefore, CPU execution is possible for small smoke tests, but full PRnet training or repeated benchmark runs would likely benefit from a GPU.

## Official Demo Status

The unmodified official `test_demo.py` failed on this CPU-only machine because the bundled checkpoint was serialized from CUDA and the script calls `torch.load(model_path)` without `map_location`.

Observed failure:

```text
RuntimeError: Attempting to deserialize object on a CUDA device but torch.cuda.is_available() is False.
```

A non-invasive CPU wrapper that monkeypatched `torch.load(..., map_location='cpu')` succeeded without editing official PRnet source.

Official demo CPU-wrapper output:

- x_true shape: 661 x 978
- y_true shape: 661 x 978
- y_pre shape: 661 x 978
- log: `logs/prnet_official_demo_cpu_wrapper.log`
- output directory: `external/PRnet/results/demo_cpu/`

## Sci-Plex 3 Input Compatibility

The official PRnet data path expects raw or profile-level AnnData with single-cell-like rows, control rows, `SMILES`, `dose`, `paired_control_index`, and split columns. The current benchmark input for TranSiGen is a pseudobulk response matrix with matched basal/control and treated means:

- `basal`: 2,200 x 2,000
- `treated`: 2,200 x 2,000
- `response`: 2,200 x 2,000
- Morgan fingerprint: 1,024 dimensions
- split manifest: seeds 1-30, four split types

This means the current input is not a drop-in replacement for official PRnet single-cell training. A benchmark-compatible adapter is required.

## Adapter Decision

The local smoke test therefore uses `prnet_adapted_sciplex3`:

- It imports the official PRnet `PGM` module from `external/PRnet/models/PRnet.py`.
- It does not modify official PRnet source code.
- It trains on the local pseudobulk Sci-Plex 3 fixed split.
- It predicts treated expression from basal/control expression and dose-scaled 1,024-bit drug fingerprint.
- It evaluates response as predicted treated expression minus basal/control expression.

This is a PRnet-style pseudobulk adapter and should not be described as a full PRnet reproduction.

## Smoke-Test Readiness

Status: ready and completed.

Output files:

- adapter: `src/models/deep/prnet_adapter.py`
- config: `configs/prnet_sciplex3_smoke.yaml`
- runner: `scripts/run_prnet_sciplex3_smoke.py`
- metrics: `results/deep_model_panel/sciplex3/prnet_smoke_metrics.csv`
- manifest: `results/deep_model_panel/sciplex3/prnet_smoke_predictions_manifest.csv`
- log: `logs/prnet_smoke_seed0.log`

One caveat: the user-requested smoke seed was 0, but the existing split manifest contains seeds 1-30 only. The runner used `requested_seed=0` for model initialization and `split_seed=1` for the fixed split assignment, without modifying the split manifest.
