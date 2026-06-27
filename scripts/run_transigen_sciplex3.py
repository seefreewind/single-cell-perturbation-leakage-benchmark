"""Run TranSiGen-style Sci-Plex 3 pseudobulk adaptation on fixed split manifests."""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.metrics.vector_metrics import rmse, rowwise_pearson, topk_direction_agreement, topk_overlap  # noqa: E402
from src.models.deep.transigen_adapter import TranSiGenAdapter  # noqa: E402


INPUT = ROOT / "data/processed/transigen/sciplex3"
OUTROOT = ROOT / "results/deep_model_panel/sciplex3"

SPLIT_MAP = {
    "random": "split_random",
    "cell_heldout": "split_cell_heldout",
    "scaffold_heldout": "split_scaffold_heldout",
    "cell_scaffold_heldout": "split_cell_scaffold_heldout",
}


def load_inputs() -> dict[str, object]:
    mats = np.load(INPUT / "matrices.npz", allow_pickle=True)
    compound = pd.read_csv(INPUT / "compound_metadata.csv")
    fp = np.load(INPUT / "compound_fingerprints.npz", allow_pickle=True)
    fp_df = pd.DataFrame(fp["morgan1024"])
    fp_df["drug_name"] = fp["drug_name"]
    meta = pd.read_csv(INPUT / "record_metadata.csv")
    splits = pd.read_csv(INPUT / "split_assignments_long.csv")
    meta = meta.merge(fp_df, on="drug_name", how="left")
    cell_dummies = pd.get_dummies(meta["cell_line"], prefix="cell", dtype=float)
    dose = np.log10(meta[["dose_value"]].to_numpy(dtype=float) + 1.0)
    basal = mats["basal"].astype(np.float32)
    response = mats["response"].astype(np.float32)
    fp_cols = list(range(1024))
    x = np.hstack([basal, meta[fp_cols].to_numpy(dtype=np.float32), dose.astype(np.float32), cell_dummies.to_numpy(dtype=np.float32)])
    return {
        "x": x.astype(np.float32),
        "y": response,
        "record_id": mats["record_id"].astype(str),
        "gene": mats["gene"].astype(str),
        "meta": meta,
        "splits": splits,
    }


def split_train_val(train_idx: np.ndarray, seed: int, val_fraction: float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed + 1729)
    order = rng.permutation(train_idx)
    n_val = max(1, int(round(len(order) * val_fraction)))
    return order[n_val:], order[:n_val]


def evaluate(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    out = {
        "mean_all_gene_pearson": float(np.nanmean(rowwise_pearson(y_true, y_pred))),
        "mean_all_gene_rmse": float(np.nanmean(rmse(y_true, y_pred))),
    }
    for k in [50, 100, 200]:
        out[f"mean_top{k}_overlap"] = float(np.nanmean(topk_overlap(y_true, y_pred, k)))
        out[f"mean_top{k}_direction_agreement"] = float(np.nanmean(topk_direction_agreement(y_true, y_pred, k)))
    return out


def leakage_for_split(split_df: pd.DataFrame) -> dict[str, float]:
    train = split_df[split_df.assignment == "train"]
    test = split_df[split_df.assignment == "test"]
    train_drugs = set(train.drug_name)
    train_scaffolds = set(train.scaffold)
    train_cells = set(train.cell_line)
    train_drug_dose = set(zip(train.drug_name, train.dose_value))
    train_drug_cell_dose = set(zip(train.drug_name, train.cell_line, train.dose_value))
    train_dose = set(train.dose_value)
    return {
        "same_drug_overlap": float(test.drug_name.isin(train_drugs).mean()) if len(test) else np.nan,
        "same_scaffold_overlap": float(test.scaffold.isin(train_scaffolds).mean()) if len(test) else np.nan,
        "same_cell_line_overlap": float(test.cell_line.isin(train_cells).mean()) if len(test) else np.nan,
        "same_drug_dose_overlap": float(pd.Series(list(zip(test.drug_name, test.dose_value))).isin(train_drug_dose).mean()) if len(test) else np.nan,
        "same_drug_cell_line_dose_overlap": float(pd.Series(list(zip(test.drug_name, test.cell_line, test.dose_value))).isin(train_drug_cell_dose).mean()) if len(test) else np.nan,
        "same_numeric_dose_overlap": float(test.dose_value.isin(train_dose).mean()) if len(test) else np.nan,
    }


def run_one(data: dict[str, object], seed: int, split_name: str, epochs: int, hidden_dim: int, batch_size: int) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    x = data["x"]
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
    if len(train_idx) == 0 or len(test_idx) == 0:
        raise ValueError(f"Empty train/test for seed={seed} split={split_name}")
    fit_idx, val_idx = split_train_val(train_idx, seed)
    start = time.time()
    model = TranSiGenAdapter(input_dim=x.shape[1], output_dim=y.shape[1], hidden_dim=hidden_dim, dropout=0.1, seed=seed)
    model.fit(
        {"x": x[fit_idx], "y": y[fit_idx]},
        {"x": x[val_idx], "y": y[val_idx]},
        {"epochs": epochs, "batch_size": batch_size, "learning_rate": 1e-3, "patience": max(2, min(5, epochs // 3 + 1))},
    )
    pred = model.predict({"x": x[test_idx]})
    runtime = time.time() - start
    metric = evaluate(y[test_idx], pred)
    metric.update(
        {
            "dataset": "sciplex3_24h_top2000",
            "seed": seed,
            "split_name": split_name,
            "split": SPLIT_MAP[split_name],
            "split_type": split_name,
            "model_name": "transigen_adapted_sciplex3",
            "n_train": len(train_idx),
            "n_val": len(val_idx),
            "n_test": len(test_idx),
            "status": "completed",
            "runtime_seconds": runtime,
            "notes": "TranSiGen-style pseudobulk adaptation; not full original TranSiGen reproduction.",
        }
    )
    pred_path = OUTROOT / "predictions" / f"transigen_seed{seed:03d}_{split_name}.npz"
    pred_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(pred_path, record_id=record_ids[test_idx], prediction=pred.astype(np.float32), target=y[test_idx].astype(np.float32))
    manifest = pd.DataFrame(
        {
            "record_id": record_ids[test_idx],
            "dataset": "sciplex3_24h_top2000",
            "split_type": split_name,
            "seed": seed,
            "model_name": "transigen_adapted_sciplex3",
            "true_response_vector_path": str(pred_path),
            "true_response_hash": "",
            "pred_response_vector_path": str(pred_path),
            "pred_response_array": "",
            "status": "completed",
            "runtime_seconds": runtime,
            "notes": "Rows are aligned to record_id order stored in the NPZ file.",
        }
    )
    leakage = leakage_for_split(split_df)
    leakage.update({"dataset": "sciplex3_24h_top2000", "seed": seed, "split_type": split_name, "model_name": "transigen_adapted_sciplex3", "status": "completed"})
    return metric, manifest, pd.DataFrame([leakage])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "pilot", "main10", "main30"], required=True)
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--hidden-dim", type=int, default=128)
    parser.add_argument("--batch-size", type=int, default=128)
    args = parser.parse_args()
    plan = {
        "smoke": ([1], ["random"], args.epochs or 2),
        "pilot": ([1, 2, 3], ["random", "scaffold_heldout", "cell_scaffold_heldout"], args.epochs or 5),
        "main10": (list(range(1, 11)), ["random", "cell_heldout", "scaffold_heldout", "cell_scaffold_heldout"], args.epochs or 12),
        "main30": (list(range(1, 31)), ["random", "cell_heldout", "scaffold_heldout", "cell_scaffold_heldout"], args.epochs or 12),
    }
    seeds, splits, epochs = plan[args.mode]
    data = load_inputs()
    metrics, manifests, leakages, failures = [], [], [], []
    OUTROOT.mkdir(parents=True, exist_ok=True)
    for seed in seeds:
        for split in splits:
            try:
                m, pred_manifest, leakage = run_one(data, seed, split, epochs, args.hidden_dim, args.batch_size)
                metrics.append(m)
                manifests.append(pred_manifest)
                leakages.append(leakage)
                print(f"completed seed={seed} split={split} pearson={m['mean_all_gene_pearson']:.4f}")
            except Exception as exc:
                failures.append({"seed": seed, "split_type": split, "status": "failure", "notes": f"{type(exc).__name__}: {exc}"})
                print(f"failure seed={seed} split={split}: {exc}", file=sys.stderr)
    metrics_df = pd.DataFrame(metrics + failures)
    manifest_df = pd.concat(manifests, ignore_index=True) if manifests else pd.DataFrame()
    leakage_df = pd.concat(leakages, ignore_index=True) if leakages else pd.DataFrame()
    metrics_df.to_csv(OUTROOT / f"transigen_{args.mode}_metrics_long.csv", index=False)
    manifest_df.to_csv(OUTROOT / f"transigen_{args.mode}_predictions_manifest.csv", index=False)
    leakage_df.to_csv(OUTROOT / f"transigen_{args.mode}_leakage_audit.csv", index=False)
    pd.DataFrame(failures).to_csv(ROOT / "logs/transigen_failures.csv", index=False)
    if args.mode in {"main10", "main30"}:
        metrics_df.to_csv(OUTROOT / "transigen_metrics_long.csv", index=False)
        manifest_df.to_csv(OUTROOT / "transigen_predictions_manifest.csv", index=False)
        leakage_df.to_csv(OUTROOT / "transigen_leakage_audit.csv", index=False)


if __name__ == "__main__":
    main()
