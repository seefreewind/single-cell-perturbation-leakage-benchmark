from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPLITS = ROOT / "data/processed/transigen/sciplex3/split_assignments_long.csv"


def test_no_test_leakage_and_excluded_handling():
    splits = pd.read_csv(SPLITS)
    for (seed, split_name), df in splits.groupby(["seed", "split_name"]):
        train = set(df.loc[df.assignment == "train", "record_id"])
        test = set(df.loc[df.assignment == "test", "record_id"])
        excluded = set(df.loc[df.assignment == "excluded", "record_id"])
        assert train.isdisjoint(test), (seed, split_name, "train/test overlap")
        assert train.isdisjoint(excluded), (seed, split_name, "train/excluded overlap")
        assert test.isdisjoint(excluded), (seed, split_name, "test/excluded overlap")


if __name__ == "__main__":
    test_no_test_leakage_and_excluded_handling()
