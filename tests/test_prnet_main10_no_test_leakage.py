from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCI = ROOT / "results/deep_model_panel/sciplex3"
SPLITS = ROOT / "data/processed/transigen/sciplex3/split_assignments_long.csv"


def test_prnet_main10_reports_train_only_validation_and_no_test_tuning():
    metrics = pd.read_csv(SCI / "prnet_main10_metrics_long.csv")
    assert metrics["validation_source"].eq("train_only").all()
    assert metrics["test_usage"].eq("final_evaluation_only").all()
    assert metrics["feature_normalization"].eq("none").all()
    completed = metrics[metrics["status"].eq("completed")]
    assert (completed["n_train"] == completed["n_fit"] + completed["n_val"]).all()


def test_prnet_main10_split_assignment_groups_are_disjoint():
    metrics = pd.read_csv(SCI / "prnet_main10_metrics_long.csv")
    splits = pd.read_csv(SPLITS)
    for row in metrics.itertuples():
        df = splits[(splits.seed == row.seed) & (splits.split_name == row.split_type)]
        train = set(df.loc[df.assignment.eq("train"), "record_id"])
        test = set(df.loc[df.assignment.eq("test"), "record_id"])
        excluded = set(df.loc[df.assignment.eq("excluded"), "record_id"])
        assert not (train & test)
        assert not (train & excluded)
        assert not (test & excluded)


if __name__ == "__main__":
    test_prnet_main10_reports_train_only_validation_and_no_test_tuning()
    test_prnet_main10_split_assignment_groups_are_disjoint()
