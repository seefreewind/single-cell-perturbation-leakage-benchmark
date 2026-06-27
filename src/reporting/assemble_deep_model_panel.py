"""Assemble existing and planned model outputs into deep-model-panel tables.

This script is intentionally conservative. It imports completed baseline
results from existing CSVs and emits explicit status rows for modern/foundation
models that are not yet runnable in the local environment.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/deep_model_panel"

OPENPROBLEMS_AGG = ROOT / "results/openproblems_multiseed_100_de/baseline_summary_by_seed_mean.csv"
OPENPROBLEMS_SEED = ROOT / "results/openproblems_multiseed_100_de/all_baseline_summary.csv"
OPENPROBLEMS_SVD = ROOT / "results/svd_ridge_multiseed_10/svd_ridge_summary_by_seed_mean.csv"
OPENPROBLEMS_NEAREST = ROOT / "results/nearest_drug_multiseed_10_de/nearest_drug_summary_by_seed_mean.csv"
OPENPROBLEMS_MOA = ROOT / "results/openproblems_moa_heldout_100_fast/baseline_summary_aggregate.csv"
SCIPLEX_AGG = ROOT / "results/sciplex3_24h_top2000_multiseed_30/baseline_summary_aggregate.csv"
SCIPLEX_SEED = ROOT / "results/sciplex3_24h_top2000_multiseed_30/baseline_summary_by_seed.csv"
SCIPLEX_TRANSIGEN = ROOT / "results/deep_model_panel/sciplex3/transigen_metrics_long.csv"

OPENPROBLEMS_LEAKAGE = ROOT / "results/leakage_overlap_100/all_leakage_overlap_by_seed.csv"
SCIPLEX_LEAKAGE = ROOT / "results/sciplex3_24h_top2000_multiseed_30/leakage_overlap_by_seed.csv"
MOA_LEAKAGE = ROOT / "results/openproblems_moa_heldout_100_fast/leakage_summary_by_seed.csv"
SCIPLEX_TRANSIGEN_LEAKAGE = ROOT / "results/deep_model_panel/sciplex3/transigen_leakage_audit.csv"

PLANNED_MODELS = [
    "prnet",
    "chemcpa_or_cpa",
    "scgpt_frozen_ridge",
    "scgpt_frozen_mlp",
    "scfoundation_frozen_ridge",
    "scfoundation_frozen_mlp",
]

SPLIT_LABELS = {
    "split_random": "random",
    "split_cell_heldout": "cell_heldout",
    "split_scaffold_heldout": "scaffold_heldout",
    "split_cell_scaffold_heldout": "joint_cell_scaffold_heldout",
    "split_moa_heldout": "moa_heldout",
}

MODEL_RENAME = {
    "global_train_mean": "global_mean",
    "svd50_ridge": "svd50_ridge_cell_drug_fp",
    "nearest_drug": "nearest_drug",
}


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def normalize_model(name: str) -> str:
    return MODEL_RENAME.get(str(name), str(name))


def normalize_aggregate(df: pd.DataFrame, dataset: str, status: str = "completed") -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    out = df.copy()
    out["dataset"] = dataset
    out["split_type"] = out["split"].map(SPLIT_LABELS).fillna(out["split"])
    out["model_name"] = out["baseline"].map(normalize_model)
    out["status"] = status
    out["notes"] = "Imported from existing benchmark result CSV."
    rename = {
        "n_seeds": "seeds",
        "mean_n_test": "mean_test_records",
        "mean_n_train": "mean_train_records",
        "mean_rowwise_pearson": "mean_all_gene_pearson",
        "sd_rowwise_pearson": "sd_all_gene_pearson",
        "mean_de_rowwise_pearson": "mean_de_gene_pearson",
        "sd_de_rowwise_pearson": "sd_de_gene_pearson",
        "mean_rmse": "mean_all_gene_rmse",
        "sd_rmse": "sd_all_gene_rmse",
        "mean_de_rmse": "mean_de_gene_rmse",
        "sd_de_rmse": "sd_de_gene_rmse",
    }
    out = out.rename(columns=rename)
    for col in [
        "seeds",
        "mean_train_records",
        "mean_test_records",
        "mean_all_gene_pearson",
        "sd_all_gene_pearson",
        "mean_de_gene_pearson",
        "sd_de_gene_pearson",
        "mean_all_gene_rmse",
        "sd_all_gene_rmse",
        "mean_de_gene_rmse",
        "sd_de_gene_rmse",
    ]:
        if col not in out.columns:
            out[col] = np.nan
    return out[
        [
            "dataset",
            "split_type",
            "model_name",
            "seeds",
            "mean_train_records",
            "mean_test_records",
            "mean_all_gene_pearson",
            "sd_all_gene_pearson",
            "mean_de_gene_pearson",
            "sd_de_gene_pearson",
            "mean_all_gene_rmse",
            "sd_all_gene_rmse",
            "mean_de_gene_rmse",
            "sd_de_gene_rmse",
            "status",
            "notes",
        ]
    ]


def planned_status_rows() -> pd.DataFrame:
    rows = []
    for dataset, splits, seeds in [
        ("openproblems_neurips2023", list(SPLIT_LABELS.values()), 0),
        (
            "sciplex3_24h_top2000",
            ["random", "cell_heldout", "scaffold_heldout", "joint_cell_scaffold_heldout"],
            0,
        ),
    ]:
        for split in splits:
            for model in PLANNED_MODELS:
                rows.append(
                    {
                        "dataset": dataset,
                        "split_type": split,
                        "model_name": model,
                        "seeds": seeds,
                        "mean_train_records": np.nan,
                        "mean_test_records": np.nan,
                        "mean_all_gene_pearson": np.nan,
                        "sd_all_gene_pearson": np.nan,
                        "mean_de_gene_pearson": np.nan,
                        "sd_de_gene_pearson": np.nan,
                        "mean_all_gene_rmse": np.nan,
                        "sd_all_gene_rmse": np.nan,
                        "mean_de_gene_rmse": np.nan,
                        "sd_de_gene_rmse": np.nan,
                        "status": "not_run_dependency_or_input_pending",
                        "notes": (
                            "External implementation, torch runtime, or required raw single-cell/foundation "
                            "embedding input is not configured locally; see docs/model_feasibility_report.md."
                        ),
                    }
                )
    return pd.DataFrame(rows)


def aggregate_transigen() -> pd.DataFrame:
    df = read_csv(SCIPLEX_TRANSIGEN)
    if df.empty:
        return pd.DataFrame()
    ok = df[df["status"].eq("completed")].copy()
    if ok.empty:
        return normalize_aggregate(pd.DataFrame(), "sciplex3_24h_top2000")
    agg = (
        ok.groupby(["split", "split_type", "model_name"], as_index=False)
        .agg(
            seeds=("seed", "nunique"),
            mean_train_records=("n_train", "mean"),
            mean_test_records=("n_test", "mean"),
            mean_all_gene_pearson=("mean_all_gene_pearson", "mean"),
            sd_all_gene_pearson=("mean_all_gene_pearson", "std"),
            mean_all_gene_rmse=("mean_all_gene_rmse", "mean"),
            sd_all_gene_rmse=("mean_all_gene_rmse", "std"),
        )
    )
    agg["split_type"] = agg["split_type"].replace({"cell_scaffold_heldout": "joint_cell_scaffold_heldout"})
    agg["dataset"] = "sciplex3_24h_top2000"
    agg["mean_de_gene_pearson"] = np.nan
    agg["sd_de_gene_pearson"] = np.nan
    agg["mean_de_gene_rmse"] = np.nan
    agg["sd_de_gene_rmse"] = np.nan
    agg["status"] = "completed"
    agg["notes"] = "TranSiGen-style Sci-Plex 3 pseudobulk adaptation, 30 seeds."
    return agg[
        [
            "dataset",
            "split_type",
            "model_name",
            "seeds",
            "mean_train_records",
            "mean_test_records",
            "mean_all_gene_pearson",
            "sd_all_gene_pearson",
            "mean_de_gene_pearson",
            "sd_de_gene_pearson",
            "mean_all_gene_rmse",
            "sd_all_gene_rmse",
            "mean_de_gene_rmse",
            "sd_de_gene_rmse",
            "status",
            "notes",
        ]
    ]


def assemble_metrics() -> pd.DataFrame:
    frames = [
        normalize_aggregate(read_csv(OPENPROBLEMS_AGG), "openproblems_neurips2023"),
        normalize_aggregate(read_csv(OPENPROBLEMS_SVD), "openproblems_neurips2023"),
        normalize_aggregate(read_csv(OPENPROBLEMS_NEAREST), "openproblems_neurips2023"),
        normalize_aggregate(read_csv(OPENPROBLEMS_MOA), "openproblems_neurips2023"),
        normalize_aggregate(read_csv(SCIPLEX_AGG), "sciplex3_24h_top2000"),
        aggregate_transigen(),
        planned_status_rows(),
    ]
    metrics = pd.concat([f for f in frames if not f.empty], ignore_index=True)
    metrics = metrics.drop_duplicates(["dataset", "split_type", "model_name", "status"], keep="first")
    return metrics.sort_values(["dataset", "split_type", "status", "model_name"]).reset_index(drop=True)


def assemble_leakage() -> pd.DataFrame:
    frames = []
    if OPENPROBLEMS_LEAKAGE.exists():
        df = pd.read_csv(OPENPROBLEMS_LEAKAGE)
        df["dataset"] = "openproblems_neurips2023"
        frames.append(df)
    if SCIPLEX_LEAKAGE.exists():
        df = pd.read_csv(SCIPLEX_LEAKAGE)
        df["dataset"] = "sciplex3_24h_top2000"
        frames.append(df)
    if SCIPLEX_TRANSIGEN_LEAKAGE.exists():
        df = pd.read_csv(SCIPLEX_TRANSIGEN_LEAKAGE)
        frames.append(df)
    if MOA_LEAKAGE.exists():
        df = pd.read_csv(MOA_LEAKAGE)
        df["dataset"] = "openproblems_neurips2023"
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    leakage = pd.concat(frames, ignore_index=True, sort=False)
    leakage["split_type"] = leakage["split"].map(SPLIT_LABELS).fillna(leakage["split"])
    return leakage


def rank_stability(metrics: pd.DataFrame) -> pd.DataFrame:
    completed = metrics[metrics["status"] == "completed"].copy()
    completed["ranking_metric"] = completed["mean_de_gene_pearson"]
    completed.loc[completed["ranking_metric"].isna(), "ranking_metric"] = completed.loc[
        completed["ranking_metric"].isna(), "mean_all_gene_pearson"
    ]
    completed = completed[completed["ranking_metric"].notna()].copy()
    rows = []
    for dataset, df in completed.groupby("dataset"):
        pivot = df.pivot_table(index="model_name", columns="split_type", values="ranking_metric", aggfunc="mean")
        for contrast in [("random", "scaffold_heldout"), ("random", "joint_cell_scaffold_heldout")]:
            a, b = contrast
            if a in pivot.columns and b in pivot.columns:
                pair = pivot[[a, b]].dropna()
                rho = pair[a].rank(ascending=False).corr(pair[b].rank(ascending=False), method="spearman")
                random_best = pair[a].idxmax() if not pair.empty else ""
                random_best_rank_in_joint = (
                    float(pair[b].rank(ascending=False).loc[random_best]) if random_best and b in pair.columns else np.nan
                )
                rows.append(
                    {
                        "dataset": dataset,
                        "contrast": f"{a}_vs_{b}",
                        "model_count": len(pair),
                        "mean_spearman": rho,
                        "sd_spearman": np.nan,
                        "percent_negative": float(rho < 0) if pd.notna(rho) else np.nan,
                        "random_best_model": random_best,
                        "random_best_rank_in_joint": random_best_rank_in_joint,
                        "status": "aggregate_existing_results",
                    }
                )
    return pd.DataFrame(rows)


def random_to_strict(metrics: pd.DataFrame) -> pd.DataFrame:
    completed = metrics[(metrics["status"] == "completed") & metrics["mean_all_gene_pearson"].notna()].copy()
    rows = []
    for (dataset, model), df in completed.groupby(["dataset", "model_name"]):
        vals = df.set_index("split_type")
        for strict in ["scaffold_heldout", "joint_cell_scaffold_heldout", "moa_heldout"]:
            if "random" in vals.index and strict in vals.index:
                rows.append(
                    {
                        "dataset": dataset,
                        "model_name": model,
                        "contrast": f"random_minus_{strict}",
                        "mean_all_gene_pearson_difference": vals.loc["random", "mean_all_gene_pearson"]
                        - vals.loc[strict, "mean_all_gene_pearson"],
                        "mean_de_gene_pearson_difference": vals.loc["random", "mean_de_gene_pearson"]
                        - vals.loc[strict, "mean_de_gene_pearson"]
                        if pd.notna(vals.loc[strict, "mean_de_gene_pearson"])
                        else np.nan,
                        "status": "aggregate_existing_results",
                    }
                )
    return pd.DataFrame(rows)


def bootstrap_ci() -> pd.DataFrame:
    frames = []
    paths = [
        ROOT / "results/manuscript_v4_stats/openproblems_ridge_performance_delta_ci.csv",
        ROOT / "results/manuscript_v4_stats/sciplex3_ridge_performance_delta_ci.csv",
    ]
    for path, dataset in zip(paths, ["openproblems_neurips2023", "sciplex3_24h_top2000"]):
        if path.exists():
            df = pd.read_csv(path)
            df["dataset"] = dataset
            df["model_name"] = "ridge_cell_drug_fp"
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    metrics = assemble_metrics()
    leakage = assemble_leakage()
    ranks = rank_stability(metrics)
    contrasts = random_to_strict(metrics)
    ci = bootstrap_ci()
    metrics.to_csv(OUT / "all_model_metrics_long.csv", index=False)
    leakage.to_csv(OUT / "all_model_leakage_audit_long.csv", index=False)
    ranks.to_csv(OUT / "model_rank_stability.csv", index=False)
    contrasts.to_csv(OUT / "random_to_strict_contrasts.csv", index=False)
    ci.to_csv(OUT / "bootstrap_ci.csv", index=False)
    (OUT / "logs").mkdir(exist_ok=True)
    (OUT / "logs/assemble_deep_model_panel.log").write_text(
        "Assembled existing completed baselines and explicit not-run status rows for modern/foundation models.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
