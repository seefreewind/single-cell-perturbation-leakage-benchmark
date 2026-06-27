from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCI = ROOT / "results/deep_model_panel/sciplex3"


EXPECTED_DIMS = {
    "transigen_basal_only": 2000,
    "transigen_drug_only": 1024,
    "transigen_basal_drug": 3024,
    "transigen_basal_drug_dose": 3025,
    "transigen_basal_drug_dose_cell": 3028,
    "transigen_full_current": 3028,
}
EXPECTED_SPLITS = {"random", "cell_heldout", "scaffold_heldout", "cell_scaffold_heldout"}


def test_transigen_ablation_input_dimensions():
    dims = pd.read_csv(SCI / "transigen_ablation_input_dims.csv")
    observed = dict(zip(dims["model_name"], dims["input_dim"]))
    assert observed == EXPECTED_DIMS


def test_transigen_ablation_seed_split_completeness():
    metrics = pd.read_csv(SCI / "transigen_ablation_metrics_long.csv")
    completed = metrics[metrics["status"].eq("completed")].copy()
    assert set(completed["model_name"]) == set(EXPECTED_DIMS)
    assert set(completed["split_type"]) == EXPECTED_SPLITS
    assert completed["seed"].nunique() == 10
    assert len(completed) == len(EXPECTED_DIMS) * len(EXPECTED_SPLITS) * 10

    counts = completed.groupby(["model_name", "split_type"])["seed"].nunique()
    assert counts.min() == counts.max() == 10
    assert completed["input_dim"].notna().all()


if __name__ == "__main__":
    test_transigen_ablation_input_dimensions()
    test_transigen_ablation_seed_split_completeness()
