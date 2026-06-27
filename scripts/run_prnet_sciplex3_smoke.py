"""Run a minimal PRnet-style Sci-Plex 3 pseudobulk smoke test."""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.metrics.vector_metrics import rmse, rowwise_pearson, topk_direction_agreement, topk_overlap  # noqa: E402
from src.models.deep.prnet_adapter import PRnetAdapter  # noqa: E402


INPUT = ROOT / "data/processed/transigen/sciplex3"
OUT = ROOT / "results/deep_model_panel/sciplex3"
LOG = ROOT / "logs/prnet_smoke_seed0.log"


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


def load_inputs() -> dict[str, object]:
    mats = np.load(INPUT / "matrices.npz", allow_pickle=True)
    meta = pd.read_csv(INPUT / "record_metadata.csv")
    splits = pd.read_csv(INPUT / "split_assignments_long.csv")
    fp = np.load(INPUT / "compound_fingerprints.npz", allow_pickle=True)
    fp_df = pd.DataFrame(fp["morgan1024"])
    fp_df["drug_name"] = fp["drug_name"]
    meta = meta.merge(fp_df, on="drug_name", how="left")
    dose = np.log10(meta[["dose_value"]].to_numpy(dtype=float) + 1.0).astype(np.float32)
    fp_cols = list(range(1024))
    drug_dose = meta[fp_cols].to_numpy(dtype=np.float32) * dose
    return {
        "control": mats["basal"].astype(np.float32),
        "treated": mats["treated"].astype(np.float32),
        "response": mats["response"].astype(np.float32),
        "record_id": mats["record_id"].astype(str),
        "splits": splits,
        "drug_dose": drug_dose.astype(np.float32),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    requested_seed = 0
    split_name = "random"
    epochs = 2
    batch_size = 128
    status = "failure"
    notes = []
    metric: dict[str, object] = {
        "dataset": "sciplex3_24h_top2000",
        "model_name": "prnet_adapted_sciplex3",
        "seed": requested_seed,
        "requested_seed": requested_seed,
        "split_type": split_name,
        "status": status,
    }
    manifest = pd.DataFrame()
    start = time.time()
    try:
        data = load_inputs()
        splits = data["splits"]
        available = sorted(int(s) for s in splits["seed"].unique())
        split_seed = requested_seed if requested_seed in available else available[0]
        if split_seed != requested_seed:
            notes.append(f"split manifest has seeds {available[0]}-{available[-1]}; used split_seed={split_seed} while keeping requested_seed=0 for model initialization")
        record_ids = np.asarray(data["record_id"])
        split_df = splits[(splits.seed == split_seed) & (splits.split_name == split_name)].copy()
        split_df["record_order"] = pd.Categorical(split_df["record_id"], categories=record_ids, ordered=True)
        split_df = split_df.sort_values("record_order")
        train_ids = set(split_df.loc[split_df.assignment == "train", "record_id"])
        test_ids = set(split_df.loc[split_df.assignment == "test", "record_id"])
        train_idx = np.array([i for i, rid in enumerate(record_ids) if rid in train_ids], dtype=int)
        test_idx = np.array([i for i, rid in enumerate(record_ids) if rid in test_ids], dtype=int)
        fit_idx, val_idx = split_train_val(train_idx, requested_seed)
        model = PRnetAdapter(x_dimension=data["control"].shape[1], drug_dimension=data["drug_dose"].shape[1], seed=requested_seed)
        model.fit(
            {
                "control": data["control"][fit_idx],
                "treated": data["treated"][fit_idx],
                "drug_dose": data["drug_dose"][fit_idx],
            },
            {
                "control": data["control"][val_idx],
                "treated": data["treated"][val_idx],
                "drug_dose": data["drug_dose"][val_idx],
            },
            {"epochs": epochs, "batch_size": batch_size, "learning_rate": 1e-3, "patience": 2},
        )
        pred_treated = model.predict_treated(
            {
                "control": data["control"][test_idx],
                "treated": data["treated"][test_idx],
                "drug_dose": data["drug_dose"][test_idx],
            },
            batch_size=batch_size,
        )
        pred_response = pred_treated - data["control"][test_idx]
        y_true = data["response"][test_idx]
        pred_path = OUT / "predictions/prnet_smoke_seed0_random.npz"
        pred_path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(pred_path, record_id=record_ids[test_idx], prediction=pred_response.astype(np.float32), target=y_true.astype(np.float32))
        metric.update(evaluate(y_true, pred_response))
        metric.update(
            {
                "status": "completed-smoke",
                "split_seed": split_seed,
                "n_train": len(train_idx),
                "n_val": len(val_idx),
                "n_test": len(test_idx),
                "epochs": epochs,
                "input_control_dim": data["control"].shape[1],
                "input_drug_dose_dim": data["drug_dose"].shape[1],
                "output_vector_length": pred_response.shape[1],
                "runtime_seconds": time.time() - start,
                "notes": "; ".join(notes) if notes else "PRnet-style pseudobulk smoke test completed.",
            }
        )
        manifest = pd.DataFrame(
            {
                "record_id": record_ids[test_idx],
                "dataset": "sciplex3_24h_top2000",
                "split_type": split_name,
                "seed": requested_seed,
                "split_seed": split_seed,
                "model_name": "prnet_adapted_sciplex3",
                "pred_response_vector_path": str(pred_path),
                "status": "completed-smoke",
                "notes": "Rows are aligned to record_id order stored in the NPZ file.",
            }
        )
    except Exception as exc:
        metric.update({"status": "failure", "runtime_seconds": time.time() - start, "notes": f"{type(exc).__name__}: {exc}"})
        notes.append(f"failure: {type(exc).__name__}: {exc}")

    pd.DataFrame([metric]).to_csv(OUT / "prnet_smoke_metrics.csv", index=False)
    manifest.to_csv(OUT / "prnet_smoke_predictions_manifest.csv", index=False)
    LOG.write_text("\n".join([f"{k}: {v}" for k, v in metric.items()]) + "\n", encoding="utf-8")
    print(metric)


if __name__ == "__main__":
    main()
