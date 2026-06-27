# TranSiGen Environment Check

Date: 2026-06-24

## Hardware and CUDA

- `nvidia-smi`: not available on this machine.
- NVIDIA GPU: not detected.
- CUDA runtime: not detected through `nvidia-smi`.
- Planned execution mode: CPU smoke/pilot/main unless a GPU environment is later provided.

## Python and Project Environment

- System Python: 3.9.6.
- Current project virtual environment Python: 3.12.13.
- Current project environment has `anndata`, `RDKit`, and `scikit-learn`.
- Current project environment did not have `torch` or `pyarrow` at the time of the check.
- A separate TranSiGen environment was created under `envs/transigen_venv/` to avoid polluting `.venv`.
- TranSiGen adaptation runtime after setup: torch 2.8.0, numpy 2.0.2, pandas 2.3.3, scikit-learn 1.6.1, anndata 0.10.9, RDKit 2025.09.2.

## TranSiGen Repository

- Repository URL: `https://github.com/myzheng-SIMM/TranSiGen.git`
- Local path: `external/TranSiGen/`
- Commit checked out: `8ec2218e2fe4fbb5f3a2c14271a8640a62c1af3a`
- License: MIT.
- Official requirements: Linux/conda-oriented environment with Python 3.6.13 and torch 1.5.1+CUDA, plus older numpy/scikit-learn/RDKit versions.
- Source modifications: none. Third-party source is not edited.

## Sci-Plex 3 Input Availability

- Processed response matrix exists: `data/processed/sciplex3/sciplex3_24h_top2000_response.h5ad`
- Processed response matrix shape: 2,200 records by 2,000 genes.
- Available layers:
  - `treated_mean_log1p_cp10k`
  - `control_mean_log1p_cp10k`
- Response target definition: treated mean log1p(CP10K) minus matched cell-line control mean log1p(CP10K).
- Raw Sci-Plex 3 file exists: `data/raw/scperturb/SrivatsanTrapnell2020_sciplex3.h5ad`
- Split manifest exists: `benchmark_resource/split_manifests/split_manifest_sciplex3.csv`
- Metadata exists: `metadata/sciplex3_24h_top2000_split_metadata.csv`
- Drug SMILES and scaffold metadata exists: `metadata/sciplex3_pubchem_smiles_with_scaffolds.csv`

## Smoke-Test Readiness

The project entered and completed CPU smoke, pilot, 10-seed, and 30-seed runs for the pseudobulk adaptation named `transigen_adapted_sciplex3`.

## Blockers and Caveats

- No NVIDIA GPU is available.
- The official TranSiGen environment is old and Linux/CUDA-specific. Running the official code exactly is not appropriate on this macOS/CPU environment.
- The current run will therefore use a minimal TranSiGen-style pseudobulk adaptation: basal/control expression plus compound representation, dose, and cell line are used to predict the treated-control response vector.
- This adaptation is suitable for integration testing and Sci-Plex 3 stress testing, but it should not be described as a full reproduction of the original TranSiGen paper.
