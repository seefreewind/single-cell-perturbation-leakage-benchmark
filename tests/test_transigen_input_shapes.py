from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/processed/transigen/sciplex3"


def test_transigen_input_shapes():
    mats = np.load(INPUT / "matrices.npz", allow_pickle=True)
    basal = mats["basal"]
    treated = mats["treated"]
    response = mats["response"]
    record_id = mats["record_id"]
    gene = mats["gene"]
    assert basal.shape == treated.shape == response.shape
    assert basal.shape == (2200, 2000)
    assert len(record_id) == basal.shape[0]
    assert len(gene) == basal.shape[1]


if __name__ == "__main__":
    test_transigen_input_shapes()
