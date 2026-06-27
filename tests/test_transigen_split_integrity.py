from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SPLITS = ROOT / "data/processed/transigen/sciplex3/split_assignments_long.csv"


def test_transigen_split_integrity():
    splits = pd.read_csv(SPLITS)
    for (seed, split_name), df in splits.groupby(["seed", "split_name"]):
        train = df[df.assignment == "train"]
        test = df[df.assignment == "test"]
        if split_name == "scaffold_heldout":
            assert test["drug_name"].isin(set(train["drug_name"])).mean() == 0
            assert test["scaffold"].isin(set(train["scaffold"])).mean() == 0
        if split_name == "cell_scaffold_heldout":
            assert test["drug_name"].isin(set(train["drug_name"])).mean() == 0
            assert test["scaffold"].isin(set(train["scaffold"])).mean() == 0
            assert test["cell_line"].isin(set(train["cell_line"])).mean() == 0
        assert pd.api.types.is_numeric_dtype(df["dose_value"])


if __name__ == "__main__":
    test_transigen_split_integrity()
