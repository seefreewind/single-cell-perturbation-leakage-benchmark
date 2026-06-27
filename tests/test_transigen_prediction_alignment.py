from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "results/deep_model_panel/sciplex3/transigen_predictions_manifest.csv"


def test_prediction_alignment():
    manifest = pd.read_csv(MANIFEST)
    assert not manifest.empty
    for path, df in manifest.groupby("pred_response_vector_path"):
        arr = np.load(path, allow_pickle=True)
        ids = arr["record_id"].astype(str)
        pred = arr["prediction"]
        target = arr["target"]
        assert pred.shape == target.shape
        assert pred.shape[1] == 2000
        assert ids.tolist() == df["record_id"].astype(str).tolist()


if __name__ == "__main__":
    test_prediction_alignment()
