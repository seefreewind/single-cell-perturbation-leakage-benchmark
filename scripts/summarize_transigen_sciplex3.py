"""Summarize TranSiGen Sci-Plex 3 runs and create pilot/main figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCI = ROOT / "results/deep_model_panel/sciplex3"
FIG = ROOT / "figures"


SPLIT_ORDER = ["random", "cell_heldout", "scaffold_heldout", "cell_scaffold_heldout"]
SPLIT_LABELS = {
    "random": "Random",
    "cell_heldout": "Cell-line held-out",
    "scaffold_heldout": "Scaffold held-out",
    "cell_scaffold_heldout": "Joint held-out",
    "joint_cell_scaffold_heldout": "Joint held-out",
}

MODEL_LABELS = {
    "global_mean": "Global mean",
    "cell_context_mean": "Cell-line mean",
    "cell_dose_mean": "Cell+dose mean",
    "drug_mean": "Drug mean",
    "ridge_cell_dose": "Ridge cell+dose",
    "ridge_drug_fp": "Ridge drug FP",
    "ridge_cell_dose_drug_fp": "Ridge cell+dose+drug FP",
    "transigen_adapted_sciplex3": "TranSiGen adapted",
}


def summarize(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    ok = df[df["status"].eq("completed")].copy()
    return (
        ok.groupby(["split_type"], as_index=False)
        .agg(
            seeds=("seed", "nunique"),
            mean_all_gene_pearson=("mean_all_gene_pearson", "mean"),
            sd_all_gene_pearson=("mean_all_gene_pearson", "std"),
            mean_all_gene_rmse=("mean_all_gene_rmse", "mean"),
            sd_all_gene_rmse=("mean_all_gene_rmse", "std"),
            mean_top50_overlap=("mean_top50_overlap", "mean"),
            mean_top100_overlap=("mean_top100_overlap", "mean"),
            mean_top200_overlap=("mean_top200_overlap", "mean"),
            mean_top50_direction_agreement=("mean_top50_direction_agreement", "mean"),
            mean_top100_direction_agreement=("mean_top100_direction_agreement", "mean"),
            mean_top200_direction_agreement=("mean_top200_direction_agreement", "mean"),
        )
    )


def contrasts(summary: pd.DataFrame) -> pd.DataFrame:
    vals = summary.set_index("split_type")
    rows = []
    for strict in ["scaffold_heldout", "cell_scaffold_heldout"]:
        if "random" in vals.index and strict in vals.index:
            rows.append(
                {
                    "dataset": "sciplex3_24h_top2000",
                    "model_name": "transigen_adapted_sciplex3",
                    "contrast": f"random_minus_{strict}",
                    "mean_all_gene_pearson_difference": vals.loc["random", "mean_all_gene_pearson"]
                    - vals.loc[strict, "mean_all_gene_pearson"],
                    "mean_all_gene_rmse_difference": vals.loc["random", "mean_all_gene_rmse"]
                    - vals.loc[strict, "mean_all_gene_rmse"],
                    "status": "completed",
                }
            )
    return pd.DataFrame(rows)


def plot_transigen(summary: pd.DataFrame, out_png: Path, out_pdf: Path, title: str) -> None:
    plot = summary.copy()
    plot["order"] = plot["split_type"].map({s: i for i, s in enumerate(SPLIT_ORDER)})
    plot = plot.sort_values("order")
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.bar(
        [SPLIT_LABELS.get(x, x) for x in plot["split_type"]],
        plot["mean_all_gene_pearson"],
        yerr=plot["sd_all_gene_pearson"].fillna(0),
        color="#4C78A8",
        capsize=3,
    )
    ax.set_ylabel("All-gene row-wise Pearson")
    ax.set_title(title)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.tick_params(axis="x", labelrotation=25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def formal_figures() -> None:
    metrics = pd.read_csv(ROOT / "results/deep_model_panel/all_model_metrics_long.csv")
    sci = metrics[(metrics.dataset == "sciplex3_24h_top2000") & (metrics.status == "completed")].copy()
    keep_models = [
        "global_mean",
        "nearest_drug_any_cell",
        "ridge_cell_dose_drug_fp",
        "ridge_drug_fp",
        "transigen_adapted_sciplex3",
    ]
    sci = sci[sci.model_name.isin(keep_models)]
    sci = sci[sci.split_type.isin(["random", "scaffold_heldout", "joint_cell_scaffold_heldout"])]
    sci["plot_split"] = sci["split_type"].map(SPLIT_LABELS)
    sci["plot_model"] = sci["model_name"].map(MODEL_LABELS).fillna(sci["model_name"])
    pivot = sci.pivot_table(index="plot_model", columns="plot_split", values="mean_all_gene_pearson", aggfunc="mean")
    pivot.to_csv(ROOT / "results/deep_model_panel/source_data/figure6_modern_model_panel.csv")
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    pivot = pivot[[c for c in ["Random", "Scaffold held-out", "Joint held-out"] if c in pivot.columns]]
    pivot.plot(kind="bar", ax=ax, color=["#4C78A8", "#F58518", "#54A24B"])
    ax.set_ylabel("All-gene row-wise Pearson")
    ax.set_xlabel("")
    ax.set_title("Sci-Plex 3 model stress test with TranSiGen adaptation")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.tick_params(axis="x", labelrotation=35)
    fig.tight_layout()
    fig.savefig(FIG / "figure6_modern_model_panel.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "figure6_modern_model_panel.pdf", bbox_inches="tight")
    plt.close(fig)

    ranks = sci.copy()
    ranks["rank"] = ranks.groupby("split_type")["mean_all_gene_pearson"].rank(ascending=False, method="average")
    heat = ranks.pivot_table(index="plot_model", columns="split_type", values="rank", aggfunc="mean")
    heat = heat[[c for c in ["random", "scaffold_heldout", "joint_cell_scaffold_heldout"] if c in heat.columns]]
    heat.to_csv(ROOT / "results/deep_model_panel/source_data/figure7_rank_transfer_with_transigen.csv")
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    im = ax.imshow(heat.to_numpy(), aspect="auto", cmap="viridis_r")
    ax.set_xticks(range(len(heat.columns)))
    ax.set_xticklabels([SPLIT_LABELS.get(c, c) for c in heat.columns], rotation=25, ha="right")
    ax.set_yticks(range(len(heat.index)))
    ax.set_yticklabels(heat.index)
    ax.set_title("Rank transfer after adding TranSiGen")
    ax.set_xlabel("")
    fig.colorbar(im, ax=ax, label="Rank")
    fig.tight_layout()
    fig.savefig(FIG / "figure7_rank_transfer_with_transigen.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "figure7_rank_transfer_with_transigen.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    FIG.mkdir(exist_ok=True)
    (ROOT / "results/deep_model_panel/source_data").mkdir(parents=True, exist_ok=True)
    pilot_summary = summarize(SCI / "transigen_pilot_metrics_long.csv")
    pilot_summary.to_csv(SCI / "transigen_pilot_summary.csv", index=False)
    plot_transigen(
        pilot_summary,
        FIG / "transigen_pilot_sciplex3_performance.png",
        FIG / "transigen_pilot_sciplex3_performance.pdf",
        "TranSiGen adaptation Sci-Plex 3 pilot",
    )
    main_summary = summarize(SCI / "transigen_metrics_long.csv")
    main_summary.to_csv(SCI / "transigen_summary.csv", index=False)
    contrast_df = contrasts(main_summary)
    contrast_df.to_csv(SCI / "transigen_random_to_strict_contrasts.csv", index=False)
    formal_figures()


if __name__ == "__main__":
    main()
