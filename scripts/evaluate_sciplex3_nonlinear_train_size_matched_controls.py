"""Evaluate train-size-matched random controls for the nonlinear Sci-Plex probe.

The control draws use only the fixed random manifest for each seed.  Joint
manifests provide target train/test counts only; their record identities are
never introduced into a random-control training, validation, or test set.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_transigen_sciplex3 import evaluate, leakage_for_split, load_inputs, split_train_val  # noqa: E402
from src.models.deep.transigen_adapter import TranSiGenAdapter  # noqa: E402


OUT = ROOT / "results/deep_model_panel/sciplex3/control_conditioned_nonlinear_train_size_matched_30"
MODEL = "control_conditioned_nonlinear_benchmark"
JOINT = "cell_scaffold_heldout"


def indices(record_ids: np.ndarray, ids: set[str]) -> np.ndarray:
    return np.array([i for i, value in enumerate(record_ids) if value in ids], dtype=int)


def run_one(data: dict[str, object], seed: int, control: str, epochs: int, hidden_dim: int, batch_size: int) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    splits = data["splits"]
    record_ids = np.asarray(data["record_id"])
    x = np.asarray(data["x"])
    y = np.asarray(data["y"])
    seed_rows = splits[splits.seed.eq(seed)]
    random_rows = seed_rows[seed_rows.split_name.eq("random")]
    joint_rows = seed_rows[seed_rows.split_name.eq(JOINT)]
    random_train = random_rows[random_rows.assignment.eq("train")]
    random_test = random_rows[random_rows.assignment.eq("test")]
    joint_train = joint_rows[joint_rows.assignment.eq("train")]
    joint_test = joint_rows[joint_rows.assignment.eq("test")]
    if min(len(random_train), len(random_test), len(joint_train), len(joint_test)) < 1:
        raise ValueError("empty fixed split partition")

    rng = np.random.default_rng(71000 + seed)
    sampled_test = random_test.iloc[rng.choice(len(random_test), size=len(joint_test), replace=False)].copy()
    if control == "test_size_matched_random":
        sampled_train = random_train.copy()
    elif control == "train_test_size_matched_random":
        sampled_train = random_train.iloc[rng.choice(len(random_train), size=len(joint_train), replace=False)].copy()
    else:
        raise ValueError(control)

    train_idx = indices(record_ids, set(sampled_train.record_id.astype(str)))
    test_idx = indices(record_ids, set(sampled_test.record_id.astype(str)))
    fit_idx, val_idx = split_train_val(train_idx, seed)
    start = time.time()
    model = TranSiGenAdapter(input_dim=x.shape[1], output_dim=y.shape[1], hidden_dim=hidden_dim, dropout=0.1, seed=seed)
    history = model.fit(
        {"x": x[fit_idx], "y": y[fit_idx]},
        {"x": x[val_idx], "y": y[val_idx]},
        {"epochs": epochs, "batch_size": batch_size, "learning_rate": 1e-3, "patience": max(2, min(5, epochs // 3 + 1))},
    )
    pred = model.predict({"x": x[test_idx]})
    elapsed = time.time() - start
    metric = evaluate(y[test_idx], pred)
    metric.update({
        "dataset": "sciplex3_24h_top2000",
        "seed": seed,
        "split_type": control,
        "model_name": MODEL,
        "n_train": len(train_idx),
        "n_fit": len(fit_idx),
        "n_val": len(val_idx),
        "n_test": len(test_idx),
        "status": "completed",
        "runtime_seconds": elapsed,
        "notes": "Benchmark-compatible control-conditioned nonlinear adaptation; not an original TranSiGen reproduction.",
    })
    output = OUT / "predictions" / f"{control}_seed{seed:03d}.npz"
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, record_id=record_ids[test_idx], prediction=pred.astype(np.float32), target=y[test_idx].astype(np.float32))
    prediction_manifest = pd.DataFrame({
        "record_id": record_ids[test_idx], "seed": seed, "split_type": control,
        "model_name": MODEL, "prediction_path": str(output), "status": "completed",
        "vector_length": y.shape[1],
    })
    split_for_audit = pd.concat([
        sampled_train.assign(assignment="train"), sampled_test.assign(assignment="test"),
    ], ignore_index=True)
    leakage = leakage_for_split(split_for_audit)
    leakage.update({"dataset": "sciplex3_24h_top2000", "seed": seed, "split_type": control, "model_name": MODEL, "status": "completed"})
    loss = pd.DataFrame(history)
    loss.insert(0, "seed", seed)
    loss.insert(1, "split_type", control)
    loss.insert(2, "model_name", MODEL)
    return metric, prediction_manifest, pd.DataFrame([leakage]), loss


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seeds", type=int, nargs="+", default=list(range(1, 31)))
    parser.add_argument("--epochs", type=int, default=12)
    parser.add_argument("--hidden-dim", type=int, default=128)
    parser.add_argument("--batch-size", type=int, default=128)
    args = parser.parse_args()
    data = load_inputs()
    metrics: list[dict] = []
    manifests: list[pd.DataFrame] = []
    audits: list[pd.DataFrame] = []
    losses: list[pd.DataFrame] = []
    failures: list[dict] = []
    OUT.mkdir(parents=True, exist_ok=True)
    for seed in args.seeds:
        for control in ("test_size_matched_random", "train_test_size_matched_random"):
            try:
                metric, manifest, audit, loss = run_one(data, seed, control, args.epochs, args.hidden_dim, args.batch_size)
                metrics.append(metric)
                manifests.append(manifest)
                audits.append(audit)
                losses.append(loss)
                print(f"completed seed={seed} control={control} pearson={metric['mean_all_gene_pearson']:.4f}", flush=True)
            except Exception as exc:
                failures.append({"seed": seed, "split_type": control, "model_name": MODEL, "status": "failure", "notes": f"{type(exc).__name__}: {exc}"})
                print(f"failure seed={seed} control={control}: {exc}", file=sys.stderr, flush=True)
    pd.DataFrame(metrics + failures).to_csv(OUT / "metrics_long.csv", index=False)
    pd.concat(manifests, ignore_index=True).to_csv(OUT / "predictions_manifest.csv", index=False) if manifests else pd.DataFrame().to_csv(OUT / "predictions_manifest.csv", index=False)
    pd.concat(audits, ignore_index=True).to_csv(OUT / "leakage_audit.csv", index=False) if audits else pd.DataFrame().to_csv(OUT / "leakage_audit.csv", index=False)
    pd.concat(losses, ignore_index=True).to_csv(OUT / "training_loss_by_epoch.csv", index=False) if losses else pd.DataFrame().to_csv(OUT / "training_loss_by_epoch.csv", index=False)
    pd.DataFrame(failures).to_csv(OUT / "failures.csv", index=False)


if __name__ == "__main__":
    main()
