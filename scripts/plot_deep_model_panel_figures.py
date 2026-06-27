"""Create v8 deep-model-panel figures from completed and feasibility outputs."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/deep_model_panel/figures"
SRC = ROOT / "results/deep_model_panel/source_data"
OUT.mkdir(parents=True, exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)


def clean_label(text: str) -> str:
    return (
        str(text)
        .replace("openproblems_neurips2023", "OpenProblems")
        .replace("sciplex3_24h_top2000", "Sci-Plex 3")
        .replace("ridge_cell_dose_drug_fp", "ridge cell+dose+FP")
        .replace("ridge_cell_drug_fp", "ridge cell+drug FP")
        .replace("svd50_ridge_cell_drug_fp", "SVD50 ridge cell+drug FP")
        .replace("nearest_drug_same_cell", "nearest drug same cell")
        .replace("nearest_drug_any_cell", "nearest drug any cell")
        .replace("cell_context_mean", "cell-context mean")
        .replace("global_mean", "global mean")
        .replace("split_", "")
        .replace("_", " ")
        .replace("cell scaffold", "cell+scaffold")
        .replace("cell drug fp", "cell+drug FP")
        .replace("svd50", "SVD50")
    )


def save(fig: plt.Figure, name: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def figure6(metrics: pd.DataFrame) -> None:
    completed = metrics[(metrics.status == "completed") & metrics.mean_all_gene_pearson.notna()].copy()
    completed = completed[completed.split_type.isin(["random", "scaffold_heldout", "joint_cell_scaffold_heldout"])]
    completed["plot_metric"] = completed["mean_de_gene_pearson"].fillna(completed["mean_all_gene_pearson"])
    completed.to_csv(SRC / "figure6_completed_model_performance.csv", index=False)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={"width_ratios": [1.2, 1.2, 1.05]})
    for ax, dataset in zip(axes[:2], ["openproblems_neurips2023", "sciplex3_24h_top2000"]):
        df = completed[completed.dataset == dataset].copy()
        top_models = (
            df[df.split_type == "random"]
            .sort_values("plot_metric", ascending=False)
            .model_name.head(7)
            .tolist()
        )
        df = df[df.model_name.isin(top_models)]
        pivot = df.pivot_table(index="model_name", columns="split_type", values="plot_metric", aggfunc="mean")
        pivot = pivot.reindex(top_models)
        x = np.arange(len(pivot.index))
        width = 0.25
        colors = {"random": "#4C78A8", "scaffold_heldout": "#F58518", "joint_cell_scaffold_heldout": "#54A24B"}
        for i, split in enumerate(["random", "scaffold_heldout", "joint_cell_scaffold_heldout"]):
            if split in pivot:
                ax.bar(x + (i - 1) * width, pivot[split], width=width, label=clean_label(split), color=colors[split])
        ax.set_title(dataset.replace("_", " "))
        ax.set_ylabel("Pearson (DE-gene when available)")
        ax.set_xticks(x)
        ax.set_xticklabels([clean_label(m) for m in pivot.index], rotation=45, ha="right", fontsize=8)
        ax.axhline(0, color="black", linewidth=0.6)
    drops = completed.pivot_table(index=["dataset", "model_name"], columns="split_type", values="plot_metric", aggfunc="mean")
    drops["random_minus_joint"] = drops.get("random") - drops.get("joint_cell_scaffold_heldout")
    drop_df = drops.reset_index().dropna(subset=["random_minus_joint"]).sort_values("random_minus_joint", ascending=False)
    drop_df.to_csv(SRC / "figure6_random_minus_joint_drop.csv", index=False)
    plot_drop = drop_df.head(10).copy()
    plot_drop["label"] = plot_drop.apply(
        lambda r: f"{clean_label(r['dataset'])}: {clean_label(r['model_name'])}", axis=1
    )
    axes[2].barh(
        plot_drop["label"],
        plot_drop.random_minus_joint,
        color="#B279A2",
    )
    axes[2].invert_yaxis()
    axes[2].set_xlabel("Random minus joint Pearson")
    axes[2].set_title("Largest drops")
    axes[2].tick_params(axis="y", labelsize=8)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.05))
    save(fig, "figure6_completed_model_stress_test")


def figure7(metrics: pd.DataFrame, ranks: pd.DataFrame) -> None:
    completed = metrics[(metrics.status == "completed") & metrics.mean_all_gene_pearson.notna()].copy()
    completed["plot_metric"] = completed["mean_de_gene_pearson"].fillna(completed["mean_all_gene_pearson"])
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, dataset in zip(axes[:2], ["openproblems_neurips2023", "sciplex3_24h_top2000"]):
        df = completed[completed.dataset == dataset].copy()
        df["rank"] = df.groupby("split_type")["plot_metric"].rank(ascending=False, method="average")
        pivot = df.pivot_table(index="model_name", columns="split_type", values="rank", aggfunc="mean")
        keep = [c for c in ["random", "cell_heldout", "scaffold_heldout", "joint_cell_scaffold_heldout", "moa_heldout"] if c in pivot.columns]
        pivot = pivot[keep].sort_values("random")
        pivot.to_csv(SRC / f"figure7_{dataset}_rank_heatmap.csv")
        im = ax.imshow(pivot.to_numpy(), aspect="auto", cmap="viridis_r")
        ax.set_title(dataset.replace("_", " "))
        ax.set_xticks(np.arange(len(keep)))
        ax.set_xticklabels([clean_label(c) for c in keep], rotation=45, ha="right", fontsize=8)
        ax.set_yticks(np.arange(len(pivot.index)))
        ax.set_yticklabels([clean_label(m) for m in pivot.index], fontsize=7)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Rank")
    ranks.to_csv(SRC / "figure7_rank_transfer_summary.csv", index=False)
    axes[2].bar(
        [f"{r.dataset.split('_')[0]}\n{r.contrast.replace('random_vs_', '')}" for r in ranks.itertuples()],
        ranks["mean_spearman"],
        color=["#4C78A8" if x >= 0 else "#E45756" for x in ranks["mean_spearman"]],
    )
    axes[2].axhline(0, color="black", linewidth=0.8)
    axes[2].set_ylabel("Rank Spearman")
    axes[2].set_title("Random-to-strict rank transfer")
    axes[2].tick_params(axis="x", labelrotation=45)
    save(fig, "figure7_ranking_transfer_completed_models")


def figure8(metrics: pd.DataFrame) -> None:
    planned = metrics[metrics.status.str.contains("not_run", na=False)].copy()
    completed = metrics[metrics.status == "completed"].copy()
    counts = pd.DataFrame(
        {
            "category": ["completed existing baselines", "planned modern/foundation status rows"],
            "count": [completed.model_name.nunique(), planned.model_name.nunique()],
        }
    )
    counts.to_csv(SRC / "figure8_model_feasibility_counts.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.bar(counts["category"], counts["count"], color=["#4C78A8", "#E45756"])
    ax.set_ylabel("Number of model names")
    ax.set_title("Model panel status")
    ax.tick_params(axis="x", labelrotation=20)
    ax.text(
        0.5,
        0.88,
        "Foundation and deep-model rows are retained as explicit not-run statuses;\nno performance values are imputed.",
        transform=ax.transAxes,
        ha="center",
        va="top",
        fontsize=9,
    )
    save(fig, "figure8_foundation_embedding_feasibility")


def main() -> None:
    metrics = pd.read_csv(ROOT / "results/deep_model_panel/all_model_metrics_long.csv")
    ranks = pd.read_csv(ROOT / "results/deep_model_panel/model_rank_stability.csv")
    figure6(metrics)
    figure7(metrics, ranks)
    figure8(metrics)


if __name__ == "__main__":
    main()
