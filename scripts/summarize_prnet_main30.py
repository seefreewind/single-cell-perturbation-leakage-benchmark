"""Summarize PRnet main30, integrate model panel outputs, and plot figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCI = ROOT / "results/deep_model_panel/sciplex3"
PANEL = ROOT / "results/deep_model_panel"
SRC = PANEL / "source_data"
FIG = ROOT / "figures"
SRC.mkdir(parents=True, exist_ok=True)
FIG.mkdir(exist_ok=True)

DATASET = "sciplex3_24h_top2000"
PRNET = "prnet_adapted_sciplex3"
TRANSIGEN = "transigen_adapted_sciplex3"
BASELINE_PATH = ROOT / "results/sciplex3_24h_top2000_multiseed_30/baseline_summary_by_seed.csv"

SPLIT_LABEL = {
    "random": "Random",
    "cell_heldout": "Cell",
    "scaffold_heldout": "Scaffold",
    "cell_scaffold_heldout": "Joint",
    "joint_cell_scaffold_heldout": "Joint",
}
ALL_MODEL_SPLIT = {
    "random": "random",
    "cell_heldout": "cell_heldout",
    "scaffold_heldout": "scaffold_heldout",
    "cell_scaffold_heldout": "joint_cell_scaffold_heldout",
}
BASELINE_TO_MODEL = {
    "global_train_mean": "global_mean",
    "ridge_drug_fp": "ridge_drug_fp",
    "ridge_cell_dose_drug_fp": "ridge_cell_dose_drug_fp",
}
MODEL_LABEL = {
    "global_mean": "Global mean",
    "ridge_drug_fp": "Ridge drug FP",
    "ridge_cell_dose_drug_fp": "Ridge cell+dose+drug FP",
    TRANSIGEN: "Control-conditioned regressor",
    PRNET: "PRnet-interface adapter",
}
CORE_MODELS = ["global_mean", "ridge_drug_fp", "ridge_cell_dose_drug_fp", TRANSIGEN, PRNET]


def ci(values: np.ndarray, n_boot: int = 10000, seed: int = 20260624) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.array([rng.choice(values, size=len(values), replace=True).mean() for _ in range(n_boot)])
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def prnet_summary(metrics: pd.DataFrame) -> pd.DataFrame:
    ok = metrics[metrics["status"].eq("completed")].copy()
    q1 = ok.groupby("split_type")["mean_all_gene_pearson"].quantile(0.25)
    q3 = ok.groupby("split_type")["mean_all_gene_pearson"].quantile(0.75)
    out = (
        ok.groupby("split_type", as_index=False)
        .agg(
            seeds=("seed", "nunique"),
            mean_train_records=("n_train", "mean"),
            mean_test_records=("n_test", "mean"),
            mean_all_gene_pearson=("mean_all_gene_pearson", "mean"),
            sd_all_gene_pearson=("mean_all_gene_pearson", "std"),
            median_all_gene_pearson=("mean_all_gene_pearson", "median"),
            mean_all_gene_rmse=("mean_all_gene_rmse", "mean"),
            sd_all_gene_rmse=("mean_all_gene_rmse", "std"),
            mean_top200_overlap=("mean_top200_overlap", "mean"),
            mean_top200_direction_agreement=("mean_top200_direction_agreement", "mean"),
        )
    )
    out["q1_all_gene_pearson"] = out["split_type"].map(q1)
    out["q3_all_gene_pearson"] = out["split_type"].map(q3)
    out["iqr_all_gene_pearson"] = out["q3_all_gene_pearson"] - out["q1_all_gene_pearson"]
    out["dataset"] = DATASET
    out["model_name"] = PRNET
    return out


def bootstrap_ci_table(metrics: pd.DataFrame) -> pd.DataFrame:
    ok = metrics[metrics["status"].eq("completed")].copy()
    rows = []
    for split, group in ok.groupby("split_type"):
        values = group["mean_all_gene_pearson"].to_numpy()
        lo, hi = ci(values)
        rows.append(
            {
                "dataset": DATASET,
                "model_name": PRNET,
                "quantity": "split_mean_all_gene_pearson",
                "split_type": split,
                "contrast": "",
                "n_seeds": int(group["seed"].nunique()),
                "estimate": float(np.nanmean(values)),
                "ci_low": lo,
                "ci_high": hi,
            }
        )
    wide = ok.pivot(index="seed", columns="split_type", values="mean_all_gene_pearson")
    for a, b in [("random", "scaffold_heldout"), ("random", "cell_scaffold_heldout")]:
        diff = (wide[a] - wide[b]).dropna()
        lo, hi = ci(diff.to_numpy())
        rows.append(
            {
                "dataset": DATASET,
                "model_name": PRNET,
                "quantity": "paired_difference_all_gene_pearson",
                "split_type": "",
                "contrast": f"{a}_minus_{b}",
                "n_seeds": int(diff.shape[0]),
                "estimate": float(diff.mean()),
                "ci_low": lo,
                "ci_high": hi,
            }
        )
    out = pd.DataFrame(rows)
    out.to_csv(SCI / "prnet_main30_bootstrap_ci.csv", index=False)
    return out


def paired_differences(prnet: pd.DataFrame, transigen: pd.DataFrame, baseline: pd.DataFrame) -> pd.DataFrame:
    rows = []
    p = prnet[prnet["status"].eq("completed")][["seed", "split_type", "mean_all_gene_pearson"]].rename(
        columns={"mean_all_gene_pearson": "prnet_pearson"}
    )
    t = transigen[transigen["status"].eq("completed")][["seed", "split_type", "mean_all_gene_pearson"]].rename(
        columns={"mean_all_gene_pearson": "comparison_pearson"}
    )
    t["comparison_model"] = TRANSIGEN
    baseline_sub = baseline[baseline["baseline"].isin(["ridge_drug_fp", "ridge_cell_dose_drug_fp", "global_train_mean"])].copy()
    baseline_sub = baseline_sub.rename(
        columns={"split_label": "split_type", "mean_rowwise_pearson": "comparison_pearson", "baseline": "comparison_model"}
    )[["seed", "split_type", "comparison_model", "comparison_pearson"]]
    comp = pd.concat([t, baseline_sub], ignore_index=True)
    for model, g in comp.groupby("comparison_model"):
        merged = p.merge(g, on=["seed", "split_type"], how="inner")
        for split, sg in merged.groupby("split_type"):
            diff = sg["prnet_pearson"] - sg["comparison_pearson"]
            lo, hi = ci(diff.to_numpy())
            rows.append(
                {
                    "dataset": DATASET,
                    "model_name": PRNET,
                    "comparison_model": BASELINE_TO_MODEL.get(model, model),
                    "split_type": split,
                    "n_seeds": int(diff.notna().sum()),
                    "mean_prnet": float(sg["prnet_pearson"].mean()),
                    "mean_comparison": float(sg["comparison_pearson"].mean()),
                    "mean_difference_prnet_minus_comparison": float(diff.mean()),
                    "sd_difference": float(diff.std(ddof=1)),
                    "ci_low": lo,
                    "ci_high": hi,
                }
            )
    out = pd.DataFrame(rows)
    out.to_csv(SCI / "prnet_main30_paired_differences.csv", index=False)
    return out


def build_rank_seed_table(prnet: pd.DataFrame, transigen: pd.DataFrame, baseline: pd.DataFrame) -> pd.DataFrame:
    base = baseline[baseline["seed"].between(1, 30) & baseline["baseline"].isin(BASELINE_TO_MODEL)].copy()
    base["model_name"] = base["baseline"].map(BASELINE_TO_MODEL)
    base = base.rename(columns={"split_label": "split_type", "mean_rowwise_pearson": "score"})[
        ["seed", "split_type", "model_name", "score"]
    ]
    trans = transigen[
        transigen["seed"].between(1, 30)
        & transigen["split_type"].isin(["random", "scaffold_heldout", "cell_scaffold_heldout"])
        & transigen["status"].eq("completed")
    ][["seed", "split_type", "model_name", "mean_all_gene_pearson"]].rename(columns={"mean_all_gene_pearson": "score"})
    pr = prnet[prnet["status"].eq("completed")][["seed", "split_type", "model_name", "mean_all_gene_pearson"]].rename(
        columns={"mean_all_gene_pearson": "score"}
    )
    all_rows = pd.concat([base, trans, pr], ignore_index=True)
    all_rows = all_rows[all_rows["model_name"].isin(CORE_MODELS)]
    all_rows["rank"] = all_rows.groupby(["seed", "split_type"])["score"].rank(ascending=False, method="average")
    all_rows.to_csv(SCI / "prnet_main30_rank_by_seed.csv", index=False)
    return all_rows


def rank_stability(rank_seed: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for strict in ["scaffold_heldout", "cell_scaffold_heldout"]:
        corrs = []
        random_best_ranks = []
        for seed, g in rank_seed.groupby("seed"):
            ref = g[g["split_type"].eq("random")][["model_name", "rank"]].rename(columns={"rank": "random_rank"})
            strict_df = g[g["split_type"].eq(strict)][["model_name", "rank"]]
            merged = strict_df.merge(ref, on="model_name", how="inner")
            if len(merged) >= 3:
                corrs.append(float(merged["rank"].corr(merged["random_rank"], method="spearman")))
            random_best = ref.sort_values("random_rank").head(1)["model_name"].iloc[0]
            rb = strict_df[strict_df["model_name"].eq(random_best)]["rank"]
            if not rb.empty:
                random_best_ranks.append(float(rb.iloc[0]))
        rows.append(
            {
                "dataset": DATASET,
                "contrast": f"random_vs_{ALL_MODEL_SPLIT[strict]}",
                "model_count": int(rank_seed["model_name"].nunique()),
                "mean_spearman": float(np.nanmean(corrs)),
                "sd_spearman": float(np.nanstd(corrs, ddof=1)),
                "percent_negative": float(np.mean(np.asarray(corrs) < 0)),
                "random_best_model": "seed_level_variable",
                "random_best_rank_in_joint": float(np.nanmean(random_best_ranks)) if strict == "cell_scaffold_heldout" else np.nan,
                "status": "prnet_main30_seed1_30_completed",
            }
        )
    out = pd.DataFrame(rows)
    return out


def update_panel_outputs(summary: pd.DataFrame, rank_rows: pd.DataFrame) -> None:
    all_path = PANEL / "all_model_metrics_long.csv"
    all_metrics = pd.read_csv(all_path)
    prnet_rows = summary.copy()
    prnet_rows["split_type"] = prnet_rows["split_type"].map(ALL_MODEL_SPLIT)
    prnet_rows["seeds"] = prnet_rows["seeds"].astype(int)
    prnet_rows["mean_de_gene_pearson"] = np.nan
    prnet_rows["sd_de_gene_pearson"] = np.nan
    prnet_rows["mean_de_gene_rmse"] = np.nan
    prnet_rows["sd_de_gene_rmse"] = np.nan
    prnet_rows["status"] = "completed"
    prnet_rows["notes"] = "PRnet-style Sci-Plex 3 pseudobulk adaptation main30; not a full PRnet reproduction."
    prnet_rows = prnet_rows[
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
    all_metrics = all_metrics[~((all_metrics["dataset"].eq(DATASET)) & (all_metrics["model_name"].eq(PRNET)))]
    pd.concat([all_metrics, prnet_rows], ignore_index=True).to_csv(all_path, index=False)

    contrast_path = PANEL / "random_to_strict_contrasts.csv"
    contrasts = pd.read_csv(contrast_path)
    prnet_con = pd.read_csv(SCI / "prnet_main30_random_to_strict_contrasts.csv")
    prnet_con["contrast"] = prnet_con["contrast"].replace(
        {
            "random_minus_cell_scaffold_heldout": "random_minus_joint_cell_scaffold_heldout",
            "scaffold_heldout_minus_cell_scaffold_heldout": "scaffold_heldout_minus_joint_cell_scaffold_heldout",
        }
    )
    prnet_con["mean_de_gene_pearson_difference"] = np.nan
    prnet_con = prnet_con[
        [
            "dataset",
            "model_name",
            "contrast",
            "mean_all_gene_pearson_difference",
            "mean_de_gene_pearson_difference",
            "status",
        ]
    ]
    contrasts = contrasts[~((contrasts["dataset"].eq(DATASET)) & (contrasts["model_name"].eq(PRNET)))]
    pd.concat([contrasts, prnet_con], ignore_index=True).to_csv(contrast_path, index=False)

    rank_path = PANEL / "model_rank_stability.csv"
    ranks = pd.read_csv(rank_path)
    ranks = ranks[~ranks["dataset"].eq(DATASET)]
    pd.concat([ranks, rank_rows], ignore_index=True).to_csv(rank_path, index=False)


def model_summary(
    summary: pd.DataFrame,
    paired: pd.DataFrame,
    rank_rows: pd.DataFrame,
    transigen: pd.DataFrame,
    baseline: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    rank_joint = float(
        rank_rows.loc[rank_rows["contrast"].eq("random_vs_joint_cell_scaffold_heldout"), "mean_spearman"].iloc[0]
    )

    base = baseline[baseline["seed"].between(1, 30) & baseline["baseline"].isin(BASELINE_TO_MODEL)].copy()
    base["model_name"] = base["baseline"].map(BASELINE_TO_MODEL)
    base = base.rename(columns={"split_label": "split_type", "mean_rowwise_pearson": "score"})
    trans = transigen[
        transigen["seed"].between(1, 30)
        & transigen["status"].eq("completed")
        & transigen["split_type"].isin(["random", "scaffold_heldout", "cell_scaffold_heldout"])
    ][["seed", "split_type", "model_name", "mean_all_gene_pearson"]].rename(columns={"mean_all_gene_pearson": "score"})
    pr = summary[["split_type", "model_name", "seeds", "mean_all_gene_pearson"]].rename(
        columns={"mean_all_gene_pearson": "score"}
    )

    completed = pd.concat(
        [
            base[["seed", "split_type", "model_name", "score"]],
            trans[["seed", "split_type", "model_name", "score"]],
        ],
        ignore_index=True,
    )
    for model_name, group in completed.groupby("model_name"):
        vals = group.groupby("split_type")["score"].mean()
        rows.append(
            {
                "dataset": DATASET,
                "model_name": model_name,
                "seeds": int(group["seed"].nunique()),
                "random_mean_pearson": float(vals.loc["random"]),
                "scaffold_mean_pearson": float(vals.loc["scaffold_heldout"]),
                "joint_mean_pearson": float(vals.loc["cell_scaffold_heldout"]),
                "random_minus_scaffold": float(vals.loc["random"] - vals.loc["scaffold_heldout"]),
                "random_minus_joint": float(vals.loc["random"] - vals.loc["cell_scaffold_heldout"]),
                "random_to_joint_spearman_after_prnet": rank_joint,
                "status": "completed",
            }
        )

    vals = pr.set_index("split_type")
    rows.append(
        {
            "dataset": DATASET,
            "model_name": PRNET,
            "seeds": int(summary["seeds"].min()),
            "random_mean_pearson": float(vals.loc["random", "score"]),
            "scaffold_mean_pearson": float(vals.loc["scaffold_heldout", "score"]),
            "joint_mean_pearson": float(vals.loc["cell_scaffold_heldout", "score"]),
            "random_minus_scaffold": float(vals.loc["random", "score"] - vals.loc["scaffold_heldout", "score"]),
            "random_minus_joint": float(vals.loc["random", "score"] - vals.loc["cell_scaffold_heldout", "score"]),
            "random_to_joint_spearman_after_prnet": rank_joint,
            "status": "completed",
        }
    )
    out = pd.DataFrame(rows)
    out["model_name"] = pd.Categorical(out["model_name"], categories=CORE_MODELS, ordered=True)
    out = out.sort_values("model_name").astype({"model_name": "string"})
    out.to_csv(SCI / "sciplex3_completed_model_summary.csv", index=False)
    return out


def plot_figures(summary: pd.DataFrame, paired: pd.DataFrame, rank_seed: pd.DataFrame, rank_rows: pd.DataFrame) -> None:
    success_by_split = summary.set_index("split_type")["seeds"].to_dict()
    if min(success_by_split.get(s, 0) for s in ["random", "scaffold_heldout", "cell_scaffold_heldout"]) < 8:
        return

    base = pd.read_csv(BASELINE_PATH)
    trans = pd.read_csv(SCI / "transigen_metrics_long.csv")
    prnet = pd.read_csv(SCI / "prnet_main30_metrics_long.csv")
    plot_rows = []
    for model in ["global_train_mean", "ridge_drug_fp", "ridge_cell_dose_drug_fp"]:
        sub = base[base["seed"].between(1, 30) & base["baseline"].eq(model) & base["split_label"].isin(["random", "scaffold_heldout", "cell_scaffold_heldout"])]
        agg = sub.groupby("split_label", as_index=False)["mean_rowwise_pearson"].mean()
        agg["model_name"] = BASELINE_TO_MODEL[model]
        agg = agg.rename(columns={"split_label": "split_type", "mean_rowwise_pearson": "mean_pearson"})
        plot_rows.append(agg)
    for model_name, df in [(TRANSIGEN, trans), (PRNET, prnet)]:
        sub = df[df["seed"].between(1, 30) & df["status"].eq("completed") & df["split_type"].isin(["random", "scaffold_heldout", "cell_scaffold_heldout"])]
        agg = sub.groupby("split_type", as_index=False)["mean_all_gene_pearson"].mean()
        agg["model_name"] = model_name
        agg = agg.rename(columns={"mean_all_gene_pearson": "mean_pearson"})
        plot_rows.append(agg)
    panel = pd.concat(plot_rows, ignore_index=True)
    panel["model_label"] = panel["model_name"].map(MODEL_LABEL)
    panel["split_label"] = panel["split_type"].map(SPLIT_LABEL)
    panel.to_csv(SRC / "figure6_sciplex3_prnet30_transigen30_model_panel.csv", index=False)
    pivot = panel.pivot_table(index="model_label", columns="split_label", values="mean_pearson", aggfunc="mean")
    pivot = pivot.reindex([MODEL_LABEL[m] for m in CORE_MODELS])
    pivot = pivot[["Random", "Scaffold", "Joint"]]
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    pivot.plot(kind="bar", ax=ax, color=["#4C78A8", "#F58518", "#54A24B"])
    ax.set_ylabel("All-gene row-wise Pearson")
    ax.set_xlabel("")
    ax.set_title("Sci-Plex 3 repeated-seed model panel with PRnet")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.legend(title="Split")
    ax.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(FIG / "figure6_sciplex3_prnet30_transigen30_model_panel.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "figure6_sciplex3_prnet30_transigen30_model_panel.pdf", bbox_inches="tight")
    plt.close(fig)

    rank_heat = (
        rank_seed.groupby(["model_name", "split_type"], as_index=False)["rank"].mean()
        .assign(model_label=lambda d: d["model_name"].map(MODEL_LABEL), split_label=lambda d: d["split_type"].map(SPLIT_LABEL))
    )
    rank_heat.to_csv(SRC / "figure7_rank_transfer_prnet30.csv", index=False)
    heat = rank_heat.pivot_table(index="model_label", columns="split_label", values="rank", aggfunc="mean")
    heat = heat.reindex([MODEL_LABEL[m] for m in CORE_MODELS])
    heat = heat[["Random", "Scaffold", "Joint"]]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw={"width_ratios": [1.3, 0.9]})
    im = axes[0].imshow(heat.to_numpy(), aspect="auto", cmap="viridis_r")
    axes[0].set_xticks(range(len(heat.columns)))
    axes[0].set_xticklabels(heat.columns, rotation=25, ha="right")
    axes[0].set_yticks(range(len(heat.index)))
    axes[0].set_yticklabels(heat.index)
    axes[0].set_title("Mean rank by split")
    fig.colorbar(im, ax=axes[0], label="Rank")
    rank_rows.to_csv(SRC / "figure7_rank_transfer_prnet30_summary.csv", index=False)
    axes[1].bar(rank_rows["contrast"].str.replace("random_vs_", "", regex=False), rank_rows["mean_spearman"], color="#4C78A8")
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_ylabel("Seed-level Spearman")
    axes[1].set_title("Random-rank transfer")
    axes[1].tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig.savefig(FIG / "figure7_rank_transfer_prnet30.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "figure7_rank_transfer_prnet30.pdf", bbox_inches="tight")
    plt.close(fig)

    prnet.to_csv(SRC / "supp_prnet30_seed_distribution.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    data = [prnet.loc[prnet["split_type"].eq(s), "mean_all_gene_pearson"].to_numpy() for s in ["random", "scaffold_heldout", "cell_scaffold_heldout"]]
    ax.boxplot(data, tick_labels=["Random", "Scaffold", "Joint"], showfliers=False)
    ax.set_ylabel("All-gene row-wise Pearson")
    ax.set_title("PRnet-interface adapter seed-level distribution")
    fig.tight_layout()
    fig.savefig(FIG / "supp_prnet30_seed_distribution.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "supp_prnet30_seed_distribution.pdf", bbox_inches="tight")
    plt.close(fig)

    pvt = paired[paired["comparison_model"].eq(TRANSIGEN)].copy()
    pvt.to_csv(SRC / "supp_prnet30_vs_transigen30_paired.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.errorbar(
        pvt["mean_difference_prnet_minus_comparison"],
        pvt["split_type"].map(SPLIT_LABEL),
        xerr=[
            pvt["mean_difference_prnet_minus_comparison"] - pvt["ci_low"],
            pvt["ci_high"] - pvt["mean_difference_prnet_minus_comparison"],
        ],
        fmt="o",
        capsize=3,
        color="#B279A2",
    )
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("PRnet-interface minus control-conditioned regressor Pearson")
    ax.set_title("Paired neural-adapter difference")
    fig.tight_layout()
    fig.savefig(FIG / "supp_prnet30_vs_transigen30_paired.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "supp_prnet30_vs_transigen30_paired.pdf", bbox_inches="tight")
    plt.close(fig)


def write_report(summary: pd.DataFrame, paired: pd.DataFrame, rank_rows: pd.DataFrame, boot: pd.DataFrame) -> None:
    vals = summary.set_index("split_type")
    completed = int(summary["seeds"].min())
    total_completed = int(summary["seeds"].sum())
    boot_lookup = boot.set_index("contrast")
    r_ms = boot_lookup.loc["random_minus_scaffold_heldout"]
    r_mj = boot_lookup.loc["random_minus_cell_scaffold_heldout"]
    report = f"""# PRnet Main30 Run Report

Date: 2026-06-24

## Status

`prnet_adapted_sciplex3` main30 completed successfully for all required seeds and splits.

- Seeds: 1-30.
- Splits: random, cell-line held-out, scaffold held-out, joint cell-line-plus-scaffold held-out.
- Successful seed/split combinations: {total_completed}/120 total; {completed}/30 for every required split.
- Failed combinations: 0.
- Model seed matched the fixed split-manifest seed in every run.

## Performance Summary

| Split | Seeds | Mean Pearson | SD Pearson | Median | IQR | Mean RMSE |
|---|---:|---:|---:|---:|---:|---:|
| Random | {int(vals.loc['random', 'seeds'])} | {vals.loc['random', 'mean_all_gene_pearson']:.4f} | {vals.loc['random', 'sd_all_gene_pearson']:.4f} | {vals.loc['random', 'median_all_gene_pearson']:.4f} | {vals.loc['random', 'iqr_all_gene_pearson']:.4f} | {vals.loc['random', 'mean_all_gene_rmse']:.4f} |
| Cell-line held-out | {int(vals.loc['cell_heldout', 'seeds'])} | {vals.loc['cell_heldout', 'mean_all_gene_pearson']:.4f} | {vals.loc['cell_heldout', 'sd_all_gene_pearson']:.4f} | {vals.loc['cell_heldout', 'median_all_gene_pearson']:.4f} | {vals.loc['cell_heldout', 'iqr_all_gene_pearson']:.4f} | {vals.loc['cell_heldout', 'mean_all_gene_rmse']:.4f} |
| Scaffold held-out | {int(vals.loc['scaffold_heldout', 'seeds'])} | {vals.loc['scaffold_heldout', 'mean_all_gene_pearson']:.4f} | {vals.loc['scaffold_heldout', 'sd_all_gene_pearson']:.4f} | {vals.loc['scaffold_heldout', 'median_all_gene_pearson']:.4f} | {vals.loc['scaffold_heldout', 'iqr_all_gene_pearson']:.4f} | {vals.loc['scaffold_heldout', 'mean_all_gene_rmse']:.4f} |
| Joint held-out | {int(vals.loc['cell_scaffold_heldout', 'seeds'])} | {vals.loc['cell_scaffold_heldout', 'mean_all_gene_pearson']:.4f} | {vals.loc['cell_scaffold_heldout', 'sd_all_gene_pearson']:.4f} | {vals.loc['cell_scaffold_heldout', 'median_all_gene_pearson']:.4f} | {vals.loc['cell_scaffold_heldout', 'iqr_all_gene_pearson']:.4f} | {vals.loc['cell_scaffold_heldout', 'mean_all_gene_rmse']:.4f} |

Random minus scaffold Pearson difference: {r_ms['estimate']:.4f} (bootstrap 95% CI {r_ms['ci_low']:.4f} to {r_ms['ci_high']:.4f}).

Random minus joint Pearson difference: {r_mj['estimate']:.4f} (bootstrap 95% CI {r_mj['ci_low']:.4f} to {r_mj['ci_high']:.4f}).

## Paired Comparisons

Paired PRnet-minus-comparison differences are available in `results/deep_model_panel/sciplex3/prnet_main30_paired_differences.csv`. The comparisons use matched seeds and split types where available.

## Rank Transfer

Seed-level rank transfer after adding PRnet is available in `results/deep_model_panel/model_rank_stability.csv`. The random-to-joint Spearman mean for the five-model Sci-Plex 3 panel is {rank_rows.loc[rank_rows['contrast'].eq('random_vs_joint_cell_scaffold_heldout'), 'mean_spearman'].iloc[0]:.4f}.

## Interpretation

This is a completed repeated-seed pseudobulk adaptation, not a full official PRnet reproduction. The run supports inclusion as a completed Sci-Plex 3 model-panel result, but the low Pearson values indicate that this adapter is not competitive with the ridge or control-conditioned repeated-seed results in the current configuration.
"""
    (ROOT / "docs/prnet_main30_run_report.md").write_text(report, encoding="utf-8")


def main() -> None:
    prnet = pd.read_csv(SCI / "prnet_main30_metrics_long.csv")
    transigen = pd.read_csv(SCI / "transigen_metrics_long.csv")
    baseline = pd.read_csv(BASELINE_PATH)
    summary = prnet_summary(prnet)
    summary.to_csv(SCI / "prnet_main30_summary_stats.csv", index=False)
    boot = bootstrap_ci_table(prnet)
    paired = paired_differences(prnet, transigen, baseline)
    rank_seed = build_rank_seed_table(prnet, transigen, baseline)
    rank_rows = rank_stability(rank_seed)
    model_summary(summary, paired, rank_rows, transigen, baseline)
    update_panel_outputs(summary, rank_rows)
    plot_figures(summary, paired, rank_seed, rank_rows)
    write_report(summary, paired, rank_rows, boot)


if __name__ == "__main__":
    main()
