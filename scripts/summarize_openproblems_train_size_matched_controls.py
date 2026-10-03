#!/usr/bin/env python3
"""Summarize OP3 matched-random controls with seed-paired inference and audits."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
CONTROL = ROOT / "results/openproblems_train_size_matched_100"
MASTER = ROOT / "results/tcbb_master"
OUT = CONTROL
N_BOOT = 10_000

RDLogger.DisableLog("rdApp.*")


def bootstrap(values: pd.Series, seed: int) -> tuple[float, float, float]:
    x = values.to_numpy(dtype=float)
    rng = np.random.default_rng(seed)
    samples = rng.choice(x, size=(N_BOOT, len(x)), replace=True).mean(axis=1)
    return float(x.mean()), float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def fingerprint(smiles: object):
    if not isinstance(smiles, str) or not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return None if mol is None else AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def comp_draw(random_test: pd.DataFrame, joint_test: pd.DataFrame, seed: int) -> np.ndarray:
    values = pd.concat([random_test["n_de_genes"], joint_test["n_de_genes"]], ignore_index=True)
    try:
        quartile = pd.qcut(values, q=4, labels=False, duplicates="drop")
        random = random_test.assign(de_quartile=quartile.iloc[: len(random_test)].to_numpy())
        joint = joint_test.assign(de_quartile=quartile.iloc[len(random_test) :].to_numpy())
    except ValueError:
        random = random_test.assign(de_quartile=0)
        joint = joint_test.assign(de_quartile=0)
    groups = random.groupby(["cell_context", "de_quartile"], dropna=False).groups
    rng = np.random.default_rng(91000 + seed)
    selected = []
    for row in joint.itertuples(index=False):
        candidates = np.asarray(groups.get((row.cell_context, row.de_quartile), []))
        if not len(candidates):
            candidates = random.index[random.cell_context.eq(row.cell_context)].to_numpy()
        if not len(candidates):
            candidates = random.index.to_numpy()
        selected.append(int(rng.choice(candidates)))
    return np.asarray(selected)


def audit(train: pd.DataFrame, test: pd.DataFrame, fmap: dict[str, object]) -> dict:
    train_drugs = set(train.drug_name.dropna())
    train_scaffold = set(train.scaffold.dropna())
    train_cells = set(train.cell_context.dropna())
    fps = [fmap.get(drug) for drug in train_drugs if fmap.get(drug) is not None]
    similarities = []
    for drug in test.drug_name.dropna().unique():
        current = fmap.get(drug)
        if current is not None and fps:
            similarities.append(max(DataStructs.BulkTanimotoSimilarity(current, fps)))
    return {
        "n_train": int(len(train)),
        "n_test": int(len(test)),
        "same_drug_overlap": float(test.drug_name.isin(train_drugs).mean()),
        "same_scaffold_overlap": float(test.scaffold.isin(train_scaffold).mean()),
        "same_cell_overlap": float(test.cell_context.isin(train_cells).mean()),
        "mean_nearest_train_tanimoto": float(np.mean(similarities)) if similarities else np.nan,
        "cell_composition": json.dumps(test.cell_context.value_counts(normalize=True).sort_index().round(6).to_dict()),
        "scaffold_composition": json.dumps(test.scaffold.fillna("MISSING").value_counts(normalize=True).sort_index().round(6).to_dict()),
    }


def main() -> None:
    metrics = pd.read_csv(CONTROL / "train_size_matched_random_metrics_by_seed.csv")
    performance = pd.read_csv(MASTER / "master_performance_seed_level.csv")
    joint = performance[
        performance.dataset.eq("OpenProblems")
        & performance.model.eq("ridge_cell_drug_fp")
        & performance.metric.eq("rowwise_pearson")
        & performance.split.eq("split_cell_scaffold_heldout")
    ][["seed", "value", "n_train", "n_test"]].rename(columns={"value": "joint_value", "n_train": "joint_n_train", "n_test": "joint_n_test"})
    contrasts = []
    for control, group in metrics.groupby("control"):
        merged = group.merge(joint, on="seed", validate="one_to_one")
        merged["difference"] = merged["value"] - merged["joint_value"]
        mean, low, high = bootstrap(merged["difference"], 20260722 + len(contrasts))
        contrasts.append({
            "control": control,
            "model": "ridge_cell_drug_fp",
            "comparison": f"{control} minus joint held-out",
            "n_paired_seeds": len(merged),
            "mean_difference": mean,
            "bootstrap_95ci_low": low,
            "bootstrap_95ci_high": high,
            "inference_unit": "seed",
        })
        merged.to_csv(OUT / f"{control}_paired_seed_level.csv", index=False)
    pd.DataFrame(contrasts).to_csv(OUT / "train_size_matched_random_contrasts.csv", index=False)

    # Reconstruct audit fields from the immutable seed manifests using the same
    # deterministic RNG sequences as the evaluator.
    audit_rows = []
    for seed in range(1, 101):
        path = ROOT / f"results/openproblems_multiseed_100_de/seed_{seed:03d}/splits/candidate_splits.csv"
        meta = pd.read_csv(path)
        data = meta[meta.condition.eq("treated")].copy()
        data["n_de_genes"] = 0  # Only used by the deterministic composition draw; actual quartiles are reconstructed below.
        # Load exact DE counts from the completed composition draws is not required
        # for size controls. For composition, use the stored selected IDs.
        random_train = data[data.split_random.eq("train")]
        random_test = data[data.split_random.eq("test")]
        joint_train = data[data.split_cell_scaffold_heldout.eq("train")]
        joint_test = data[data.split_cell_scaffold_heldout.eq("test")]
        drug_smiles = data[["drug_name", "smiles"]].drop_duplicates("drug_name")
        fmap = dict(zip(drug_smiles.drug_name, drug_smiles.smiles.map(fingerprint)))
        rng = np.random.default_rng(81000 + seed)
        test_size = rng.choice(random_test.index.to_numpy(), size=len(joint_test), replace=False)
        train_size = rng.choice(random_train.index.to_numpy(), size=len(joint_train), replace=False)
        selected_comp = pd.read_csv(CONTROL / "train_size_matched_random_draws.csv")
        selected_comp = selected_comp[(selected_comp.seed == seed) & selected_comp.control.eq("composition_train_test_size_matched_random")].sample_id.astype(str)
        controls = {
            "joint_heldout": (joint_train, joint_test),
            "test_size_matched_random": (random_train, random_test.loc[test_size]),
            "train_test_size_matched_random": (random_train.loc[train_size], random_test.loc[test_size]),
            "composition_train_test_size_matched_random": (
                random_train.loc[train_size],
                random_test.set_index(random_test["sample_id"].astype(str)).loc[selected_comp].reset_index(drop=True),
            ),
        }
        for control, (train, test) in controls.items():
            result = audit(train, test, fmap)
            result.update({"seed": seed, "control": control})
            audit_rows.append(result)
    audit_df = pd.DataFrame(audit_rows)
    audit_df.to_csv(OUT / "train_size_matched_random_leakage_audit.csv", index=False)
    summary = audit_df.groupby("control", as_index=False).agg(
        mean_n_train=("n_train", "mean"),
        mean_n_test=("n_test", "mean"),
        same_drug_overlap=("same_drug_overlap", "mean"),
        same_scaffold_overlap=("same_scaffold_overlap", "mean"),
        same_cell_overlap=("same_cell_overlap", "mean"),
        mean_nearest_train_tanimoto=("mean_nearest_train_tanimoto", "mean"),
    )
    summary.to_csv(OUT / "train_size_matched_random_audit_summary.csv", index=False)


if __name__ == "__main__":
    main()
