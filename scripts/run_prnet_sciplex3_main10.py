"""Run PRnet-style Sci-Plex 3 pseudobulk main10 benchmark."""

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
PRED = OUT / "predictions"
LOGDIR = ROOT / "logs/prnet_main10"

DATASET = "sciplex3_24h_top2000"
MODEL = "prnet_adapted_sciplex3"
SEEDS = list(range(1, 11))
SPLITS = ["random", "scaffold_heldout", "cell_scaffold_heldout"]
SPLIT_TO_CANONICAL = {
    "random": "split_random",
    "scaffold_heldout": "split_scaffold_heldout",
    "cell_scaffold_heldout": "split_cell_scaffold_heldout",
}
EPOCHS = 20
BATCH_SIZE = 128
PATIENCE = 5
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-8


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
        "gene": mats["gene"].astype(str),
        "splits": splits,
        "drug_dose": drug_dose.astype(np.float32),
    }


def split_indices(data: dict[str, object], seed: int, split_name: str) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    record_ids = np.asarray(data["record_id"])
    splits = data["splits"]
    split_df = splits[(splits.seed == seed) & (splits.split_name == split_name)].copy()
    if split_df.empty:
        raise ValueError(f"No split rows for seed={seed} split={split_name}")
    split_df["record_order"] = pd.Categorical(split_df["record_id"], categories=record_ids, ordered=True)
    split_df = split_df.sort_values("record_order")
    train_ids = set(split_df.loc[split_df.assignment == "train", "record_id"])
    test_ids = set(split_df.loc[split_df.assignment == "test", "record_id"])
    excluded_ids = set(split_df.loc[split_df.assignment == "excluded", "record_id"])
    if train_ids & test_ids or train_ids & excluded_ids or test_ids & excluded_ids:
        raise ValueError(f"Overlapping assignment groups for seed={seed} split={split_name}")
    train_idx = np.array([i for i, rid in enumerate(record_ids) if rid in train_ids], dtype=int)
    test_idx = np.array([i for i, rid in enumerate(record_ids) if rid in test_ids], dtype=int)
    excluded_idx = np.array([i for i, rid in enumerate(record_ids) if rid in excluded_ids], dtype=int)
    fit_idx, val_idx = split_train_val(train_idx, seed)
    return split_df, train_idx, fit_idx, val_idx, test_idx


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


def run_one(data: dict[str, object], seed: int, split_name: str) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame, str]:
    start = time.time()
    record_ids = np.asarray(data["record_id"])
    split_df, train_idx, fit_idx, val_idx, test_idx = split_indices(data, seed, split_name)
    model = PRnetAdapter(
        x_dimension=data["control"].shape[1],
        drug_dimension=data["drug_dose"].shape[1],
        seed=seed,
    )
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
        {
            "epochs": EPOCHS,
            "batch_size": BATCH_SIZE,
            "learning_rate": LEARNING_RATE,
            "weight_decay": WEIGHT_DECAY,
            "patience": PATIENCE,
        },
    )
    pred_treated = model.predict_treated(
        {
            "control": data["control"][test_idx],
            "treated": data["treated"][test_idx],
            "drug_dose": data["drug_dose"][test_idx],
        },
        batch_size=BATCH_SIZE,
    )
    pred_response = pred_treated - data["control"][test_idx]
    y_true = data["response"][test_idx]
    pred_path = PRED / f"prnet_main10_seed{seed:03d}_{split_name}.npz"
    np.savez_compressed(pred_path, record_id=record_ids[test_idx], prediction=pred_response.astype(np.float32), target=y_true.astype(np.float32))
    metric = evaluate(y_true, pred_response)
    metric.update(
        {
            "dataset": DATASET,
            "model_name": MODEL,
            "seed": seed,
            "split_name": split_name,
            "split": SPLIT_TO_CANONICAL[split_name],
            "split_type": split_name,
            "n_train": int(len(train_idx)),
            "n_fit": int(len(fit_idx)),
            "n_val": int(len(val_idx)),
            "n_test": int(len(test_idx)),
            "epochs": EPOCHS,
            "batch_size": BATCH_SIZE,
            "learning_rate": LEARNING_RATE,
            "patience": PATIENCE,
            "validation_source": "train_only",
            "test_usage": "final_evaluation_only",
            "feature_normalization": "none",
            "input_control_dim": int(data["control"].shape[1]),
            "input_drug_dose_dim": int(data["drug_dose"].shape[1]),
            "output_vector_length": int(pred_response.shape[1]),
            "status": "completed",
            "runtime_seconds": float(time.time() - start),
            "notes": "PRnet-style pseudobulk adaptation using official PGM module; not a full PRnet reproduction.",
        }
    )
    manifest = pd.DataFrame(
        {
            "record_id": record_ids[test_idx],
            "dataset": DATASET,
            "split_type": split_name,
            "seed": seed,
            "model_name": MODEL,
            "pred_response_vector_path": str(pred_path),
            "status": "completed",
            "notes": "Rows are aligned to record_id order stored in the NPZ file.",
        }
    )
    leakage = leakage_for_split(split_df)
    leakage.update(
        {
            "dataset": DATASET,
            "seed": seed,
            "split_type": split_name,
            "model_name": MODEL,
            "status": "completed",
        }
    )
    log = (
        f"status=completed\nseed={seed}\nsplit={split_name}\n"
        f"n_train={len(train_idx)}\nn_val={len(val_idx)}\nn_test={len(test_idx)}\n"
        f"mean_all_gene_pearson={metric['mean_all_gene_pearson']}\n"
        f"runtime_seconds={metric['runtime_seconds']}\n"
    )
    return metric, manifest, pd.DataFrame([leakage]), log


def failure_row(seed: int, split_name: str, exc: Exception) -> tuple[dict[str, object], pd.DataFrame, pd.DataFrame, str]:
    msg = f"{type(exc).__name__}: {exc}"
    metric = {
        "dataset": DATASET,
        "model_name": MODEL,
        "seed": seed,
        "split_name": split_name,
        "split": SPLIT_TO_CANONICAL.get(split_name, split_name),
        "split_type": split_name,
        "status": "failure",
        "notes": msg,
        "validation_source": "train_only",
        "test_usage": "final_evaluation_only",
        "feature_normalization": "none",
    }
    leakage = pd.DataFrame(
        [
            {
                "dataset": DATASET,
                "seed": seed,
                "split_type": split_name,
                "model_name": MODEL,
                "status": "failure",
                "notes": msg,
            }
        ]
    )
    return metric, pd.DataFrame(), leakage, f"status=failure\nseed={seed}\nsplit={split_name}\nnotes={msg}\n"


def contrast_table(metrics: pd.DataFrame) -> pd.DataFrame:
    ok = metrics[metrics["status"].eq("completed")].copy()
    rows = []
    wide = ok.pivot(index="seed", columns="split_type", values="mean_all_gene_pearson")
    for a, b in [("random", "scaffold_heldout"), ("random", "cell_scaffold_heldout"), ("scaffold_heldout", "cell_scaffold_heldout")]:
        if a in wide and b in wide:
            diff = wide[a] - wide[b]
            rows.append(
                {
                    "dataset": DATASET,
                    "model_name": MODEL,
                    "contrast": f"{a}_minus_{b}",
                    "n_seeds": int(diff.notna().sum()),
                    "mean_all_gene_pearson_difference": float(diff.mean()),
                    "sd_all_gene_pearson_difference": float(diff.std(ddof=1)),
                    "status": "completed" if diff.notna().sum() else "failure",
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PRED.mkdir(parents=True, exist_ok=True)
    LOGDIR.mkdir(parents=True, exist_ok=True)
    data = load_inputs()
    metrics: list[dict[str, object]] = []
    manifests: list[pd.DataFrame] = []
    leakages: list[pd.DataFrame] = []
    for seed in SEEDS:
        for split_name in SPLITS:
            try:
                metric, manifest, leakage, log = run_one(data, seed, split_name)
                print(f"completed seed={seed} split={split_name} pearson={metric['mean_all_gene_pearson']:.4f}")
            except Exception as exc:
                metric, manifest, leakage, log = failure_row(seed, split_name, exc)
                print(f"failure seed={seed} split={split_name}: {exc}", file=sys.stderr)
            metrics.append(metric)
            if not manifest.empty:
                manifests.append(manifest)
            leakages.append(leakage)
            (LOGDIR / f"seed{seed:03d}_{split_name}.log").write_text(log, encoding="utf-8")
            pd.DataFrame(metrics).to_csv(OUT / "prnet_main10_metrics_long.csv", index=False)
            (pd.concat(manifests, ignore_index=True) if manifests else pd.DataFrame()).to_csv(
                OUT / "prnet_main10_predictions_manifest.csv", index=False
            )
            pd.concat(leakages, ignore_index=True).to_csv(OUT / "prnet_main10_leakage_audit.csv", index=False)
    metrics_df = pd.DataFrame(metrics)
    contrast_table(metrics_df).to_csv(OUT / "prnet_main10_random_to_strict_contrasts.csv", index=False)


if __name__ == "__main__":
    main()
