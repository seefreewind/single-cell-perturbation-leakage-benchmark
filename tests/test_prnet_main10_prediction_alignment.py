from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCI = ROOT / "results/deep_model_panel/sciplex3"


def test_prnet_main10_prediction_record_alignment_and_shape():
    metrics = pd.read_csv(SCI / "prnet_main10_metrics_long.csv")
    manifest = pd.read_csv(SCI / "prnet_main10_predictions_manifest.csv")
    completed = metrics[metrics["status"].eq("completed")]
    assert not completed.empty

    for row in completed.itertuples():
        sub = manifest[
            manifest["seed"].eq(row.seed)
            & manifest["split_type"].eq(row.split_type)
            & manifest["status"].eq("completed")
        ]
        assert len(sub) == row.n_test
        path = Path(sub["pred_response_vector_path"].iloc[0])
        assert path.exists()
        arr = np.load(path, allow_pickle=True)
        assert arr["prediction"].shape == arr["target"].shape
        assert arr["prediction"].shape[1] == 2000
        assert len(arr["record_id"]) == arr["prediction"].shape[0]
        assert list(arr["record_id"].astype(str)) == sub["record_id"].astype(str).tolist()


if __name__ == "__main__":
    test_prnet_main10_prediction_record_alignment_and_shape()
