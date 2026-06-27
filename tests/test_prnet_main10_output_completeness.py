from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCI = ROOT / "results/deep_model_panel/sciplex3"

EXPECTED_SEEDS = set(range(1, 11))
EXPECTED_SPLITS = {"random", "scaffold_heldout", "cell_scaffold_heldout"}


def test_prnet_main10_has_status_row_for_every_seed_split():
    metrics = pd.read_csv(SCI / "prnet_main10_metrics_long.csv")
    observed = set(zip(metrics["seed"], metrics["split_type"]))
    expected = {(seed, split) for seed in EXPECTED_SEEDS for split in EXPECTED_SPLITS}
    assert observed == expected
    assert metrics["status"].isin({"completed", "failure"}).all()


def test_prnet_main10_required_outputs_exist():
    required = [
        SCI / "prnet_main10_metrics_long.csv",
        SCI / "prnet_main10_predictions_manifest.csv",
        SCI / "prnet_main10_leakage_audit.csv",
        SCI / "prnet_main10_random_to_strict_contrasts.csv",
    ]
    for path in required:
        assert path.exists()

    metrics = pd.read_csv(SCI / "prnet_main10_metrics_long.csv")
    for split, sub in metrics.groupby("split_type"):
        assert split in EXPECTED_SPLITS
        assert set(sub["seed"]) == EXPECTED_SEEDS


if __name__ == "__main__":
    test_prnet_main10_has_status_row_for_every_seed_split()
    test_prnet_main10_required_outputs_exist()
