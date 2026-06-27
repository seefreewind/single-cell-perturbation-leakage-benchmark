from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCI = ROOT / "results/deep_model_panel/sciplex3"


def test_prnet_main10_scaffold_and_joint_leakage_zero():
    leakage = pd.read_csv(SCI / "prnet_main10_leakage_audit.csv")
    completed = leakage[leakage["status"].eq("completed")].copy()

    scaffold = completed[completed["split_type"].eq("scaffold_heldout")]
    assert not scaffold.empty
    assert scaffold["same_drug_overlap"].fillna(1).eq(0).all()
    assert scaffold["same_scaffold_overlap"].fillna(1).eq(0).all()

    joint = completed[completed["split_type"].eq("cell_scaffold_heldout")]
    assert not joint.empty
    assert joint["same_drug_overlap"].fillna(1).eq(0).all()
    assert joint["same_scaffold_overlap"].fillna(1).eq(0).all()
    assert joint["same_cell_line_overlap"].fillna(1).eq(0).all()


if __name__ == "__main__":
    test_prnet_main10_scaffold_and_joint_leakage_zero()
