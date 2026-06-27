"""Run input ablations for the TranSiGen-style Sci-Plex 3 adaptation."""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_transigen_sciplex3 import evaluate, leakage_for_split, split_train_val  # noqa: E402
from src.models.deep.transigen_adapter import TranSiGenAdapter  # noqa: E402


INPUT = ROOT / "data/processed/transigen/sciplex3"
OUTROOT = ROOT / "results/deep_model_panel/sciplex3"
LOG = ROOT / "logs/transigen_ablation_input_dims.csv"

FEATURE_SETS = {
    "transigen_basal_only": ["basal"],
    "transigen_drug_only": ["drug"],
    "transigen_basal_drug": ["basal", "drug"],
    "transigen_basal_drug_dose": ["basal", "drug", "dose"],
    "transigen_basal_drug_dose_cell": ["basal", "drug", "dose", "cell"],
    "transigen_full_current": ["basal", "drug", "dose", "cell"],
}

SPLITS = ["random", "cell_heldout", "scaffold_heldout", "cell_scaffold_heldout"]
SPLIT_TO_CANONICAL = {
    "random": "split_random",
    "cell_heldout": "split_cell_heldout",
    "scaffold_heldout": "split_scaffold_heldout",
    "cell_scaffold_heldout": "split_cell_scaffold_heldout",
}


def load_base() -> dict[str, object]:
    mats = np.load(INPUT / "matrices.npz", allow_pickle=True)
    meta = pd.read_csv(INPUT / "record_metadata.csv")
    splits = pd.read_csv(INPUT / "split_assignments_long.csv")
    fp = np.load(INPUT / "compound_fingerprints.npz", allow_pickle=True)
    fp_df = pd.DataFrame(fp["morgan1024"])
    fp_df["drug_name"] = fp["drug_name"]
    meta = meta.merge(fp_df, on="drug_name", how="left")
    cell = pd.get_dummies(meta["cell_line"], prefix="cell", dtype=float).to_numpy(dtype=np.float32)
    dose = np.log10(meta[["dose_value"]].to_numpy(dtype=float) + 1.0).astype(np.float32)
    return {
        "basal": mats["basal"].astype(np.float32),
        "drug": meta[list(range(1024))].to_numpy(dtype=np.float32),
        "dose": dose,
        "cell": cell,
        "y": mats["response"].astype(np.float32),
        "record_id": mats["record_id"].astype(str),
        "splits": splits,
    }


def make_x(data: dict[str, object], parts: list[str]) -> np.ndarray:
    return np.hstack([data[p] for p in parts]).astype(np.float32)


def run_one(data: dict[str, object], x: np.ndarray, model_name: str, seed: int, split_name: str, epochs: int) -> tuple[dict, pd.DataFrame]:
    y = data["y"]
    record_ids = np.asarray(data["record_id"])
    splits = data["splits"]
    split_df = splits[(splits.seed == seed) & (splits.split_name == split_name)].copy()
    split_df["record_order"] = pd.Categorical(split_df["record_id"], categories=record_ids, ordered=True)
    split_df = split_df.sort_values("record_order")
    train_ids = set(split_df.loc[split_df.assignment == "train", "record_id"])
    test_ids = set(split_df.loc[split_df.assignment == "test", "record_id"])
    train_idx = np.array([i for i, rid in enumerate(record_ids) if rid in train_ids], dtype=int)
    test_idx = np.array([i for i, rid in enumerate(record_ids) if rid in test_ids], dtype=int)
    fit_idx, val_idx = split_train_val(train_idx, seed)
    start = time.time()
    model = TranSiGenAdapter(input_dim=x.shape[1], output_dim=y.shape[1], hidden_dim=128, dropout=0.1, seed=seed)
    model.fit(
        {"x": x[fit_idx], "y": y[fit_idx]},
        {"x": x[val_idx], "y": y[val_idx]},
        {"epochs": epochs, "batch_size": 256, "learning_rate": 1e-3, "patience": max(2, min(5, epochs // 3 + 1))},
    )
    pred = model.predict({"x": x[test_idx]})
    runtime = time.time() - start
    metric = evaluate(y[test_idx], pred)
    metric.update(
        {
            "dataset": "sciplex3_24h_top2000",
            "seed": seed,
            "split_type": split_name,
            "split": SPLIT_TO_CANONICAL[split_name],
            "model_name": model_name,
            "input_dim": x.shape[1],
            "n_train": len(train_idx),
            "n_val": len(val_idx),
            "n_test": len(test_idx),
            "status": "completed",
            "runtime_seconds": runtime,
            "notes": "TranSiGen-style input ablation; same architecture and split protocol as full adaptation.",
        }
    )
    leakage = leakage_for_split(split_df)
    metric.update(leakage)
    return metric, pd.DataFrame(
        {
            "record_id": record_ids[test_idx],
            "dataset": "sciplex3_24h_top2000",
            "split_type": split_name,
            "seed": seed,
            "model_name": model_name,
            "status": "completed",
        }
    )


def main() -> None:
    OUTROOT.mkdir(parents=True, exist_ok=True)
    data = load_base()
    rows = []
    manifest_rows = []
    dim_rows = []
    failures = []
    for model_name, parts in FEATURE_SETS.items():
        x = make_x(data, parts)
        dim_rows.append({"model_name": model_name, "feature_parts": "+".join(parts), "input_dim": x.shape[1]})
        for seed in range(1, 11):
            for split in SPLITS:
                try:
                    metric, manifest = run_one(data, x, model_name, seed, split, epochs=10)
                    rows.append(metric)
                    manifest_rows.append(manifest)
                    print(f"completed {model_name} seed={seed} split={split} pearson={metric['mean_all_gene_pearson']:.4f}")
                except Exception as exc:
                    failures.append({"model_name": model_name, "seed": seed, "split_type": split, "status": "failure", "notes": f"{type(exc).__name__}: {exc}"})
                    print(f"failure {model_name} seed={seed} split={split}: {exc}", file=sys.stderr)
    metrics = pd.DataFrame(rows + failures)
    metrics.to_csv(OUTROOT / "transigen_ablation_metrics_long.csv", index=False)
    (pd.concat(manifest_rows, ignore_index=True) if manifest_rows else pd.DataFrame()).to_csv(
        OUTROOT / "transigen_ablation_predictions_manifest.csv", index=False
    )
    dims = pd.DataFrame(dim_rows)
    dims.to_csv(LOG, index=False)
    dims.to_csv(OUTROOT / "transigen_ablation_input_dims.csv", index=False)
    if failures:
        pd.DataFrame(failures).to_csv(ROOT / "logs/transigen_ablation_failures.csv", index=False)


if __name__ == "__main__":
    main()
