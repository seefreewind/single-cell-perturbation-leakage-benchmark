#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_DIR="${ROOT}/envs/transigen_venv"

if [ ! -d "${ENV_DIR}" ]; then
  python3 -m venv "${ENV_DIR}"
fi

"${ENV_DIR}/bin/python" -m pip install --upgrade pip setuptools wheel
"${ENV_DIR}/bin/python" -m pip install numpy pandas scipy scikit-learn matplotlib anndata rdkit torch

echo "TranSiGen adaptation environment ready: ${ENV_DIR}"
