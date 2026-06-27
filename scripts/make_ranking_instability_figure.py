#!/usr/bin/env python3
"""Create Fig. 6 ranking-instability manuscript figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
    }
)


ROOT = Path(".")
OUT = ROOT / "figures/manuscript_main"
SRC = OUT / "source_data"
OUT.mkdir(parents=True, exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)

SPLIT_LABELS = {
    "split_random": "Random",
    "split_cell_heldout": "Cell\nheld-out",
    "split_scaffold_heldout": "Scaffold\nheld-out",
    "split_cell_scaffold_heldout": "Cell + scaffold\nheld-out",
}

BASELINE_LABELS = {
    "global_train_mean": "Global mean",
    "cell_context_mean": "Cell mean",
    "cell_dose_mean": "Cell+dose mean",
    "drug_mean": "Drug mean",
    "ridge_cell": "Ridge cell",
    "ridge_cell_dose": "Ridge cell+dose",
    "ridge_drug_fp": "Ridge drug FP",
    "ridge_cell_drug_fp": "Ridge cell+drug FP",
    "ridge_cell_dose_drug_fp": "Ridge cell+dose+drug FP",
}


def load_mean_ranks(path: Path, dataset: str, baselines: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df[df["baseline"].isin(baselines)].copy()
    df["dataset"] = dataset
    df["split_label"] = df["split"].map(SPLIT_LABELS)
    df["baseline_label"] = df["baseline"].map(BASELINE_LABELS)
    return df


def rank_matrix(df: pd.DataFrame, baselines: list[str]) -> tuple[np.ndarray, list[str], list[str]]:
    splits = ["split_random", "split_cell_heldout", "split_scaffold_heldout", "split_cell_scaffold_heldout"]
    mat = np.full((len(baselines), len(splits)), np.nan)
    for i, baseline in enumerate(baselines):
        for j, split in enumerate(splits):
            vals = df.loc[df["baseline"].eq(baseline) & df["split"].eq(split), "mean_rank"]
            if not vals.empty:
                mat[i, j] = vals.iloc[0]
    return mat, [BASELINE_LABELS[b] for b in baselines], [SPLIT_LABELS[s] for s in splits]


def draw_rank_heatmap(ax, df: pd.DataFrame, baselines: list[str], title: str, vmax: float) -> None:
    mat, ylabels, xlabels = rank_matrix(df, baselines)
    im = ax.imshow(mat, cmap="YlGnBu_r", vmin=1, vmax=vmax, aspect="auto")
    ax.set_title(title, loc="left", fontweight="bold", pad=6)
    ax.set_xticks(np.arange(len(xlabels)), xlabels)
    ax.set_yticks(np.arange(len(ylabels)), ylabels)
    ax.tick_params(axis="both", length=0)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if np.isfinite(mat[i, j]):
                color = "white" if mat[i, j] <= 2.0 else "black"
                ax.text(j, i, f"{mat[i, j]:.1f}", ha="center", va="center", fontsize=6.5, color=color)
    for spine in ax.spines.values():
        spine.set_visible(False)
    return im


def correlation_summary(path: Path, dataset: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    if {"mean_spearman", "sd_spearman"}.issubset(df.columns):
        out = df.copy()
    else:
        out = (
            df.groupby("comparison_split")
            .agg(
                mean_spearman=("spearman_rank_correlation", "mean"),
                sd_spearman=("spearman_rank_correlation", "std"),
            )
            .reset_index()
        )
    out["dataset"] = dataset
    out["comparison_label"] = out["comparison_split"].map(SPLIT_LABELS)
    return out


def slope_data(open_df: pd.DataFrame, sci_df: pd.DataFrame) -> pd.DataFrame:
    specs = [
        ("OpenProblems", open_df, "ridge_cell_drug_fp"),
        ("Sci-Plex 3", sci_df, "ridge_cell_dose_drug_fp"),
    ]
    rows = []
    for dataset, df, baseline in specs:
        for split in ["split_random", "split_scaffold_heldout", "split_cell_scaffold_heldout"]:
            row = df[(df["baseline"].eq(baseline)) & (df["split"].eq(split))].iloc[0]
            rows.append(
                {
                    "dataset": dataset,
                    "baseline": baseline,
                    "split": split,
                    "split_label": SPLIT_LABELS[split],
                    "mean_rank": row["mean_rank"],
                    "mean_score": row["mean_score"],
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    open_baselines = [
        "ridge_cell_drug_fp",
        "ridge_drug_fp",
        "drug_mean",
        "cell_context_mean",
        "global_train_mean",
    ]
    sci_baselines = [
        "ridge_cell_dose_drug_fp",
        "cell_dose_mean",
        "ridge_drug_fp",
        "drug_mean",
        "cell_context_mean",
        "global_train_mean",
        "ridge_cell_dose",
    ]
    open_df = load_mean_ranks(ROOT / "results/ranking_instability_100/mean_rank_summary.csv", "OpenProblems", open_baselines)
    sci_df = load_mean_ranks(ROOT / "results/sciplex3_24h_top2000_multiseed_30/mean_rank_summary.csv", "Sci-Plex 3", sci_baselines)
    corr = pd.concat(
        [
            correlation_summary(ROOT / "results/ranking_instability_100/split_rank_correlations.csv", "OpenProblems"),
            correlation_summary(ROOT / "results/sciplex3_24h_top2000_multiseed_30/rank_correlation_summary.csv", "Sci-Plex 3"),
        ],
        ignore_index=True,
    )
    slopes = slope_data(open_df, sci_df)

    open_df.to_csv(SRC / "figure6a_openproblems_mean_ranks.csv", index=False)
    sci_df.to_csv(SRC / "figure6b_sciplex3_mean_ranks.csv", index=False)
    corr.to_csv(SRC / "figure6c_rank_correlations.csv", index=False)
    slopes.to_csv(SRC / "figure6d_key_model_rank_shift.csv", index=False)

    fig = plt.figure(figsize=(7.2, 6.4), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], width_ratios=[1.0, 1.0])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    im = draw_rank_heatmap(ax_a, open_df, open_baselines, "a  OpenProblems model ranks", vmax=6)
    draw_rank_heatmap(ax_b, sci_df, sci_baselines, "b  Sci-Plex 3 model ranks", vmax=7)
    cbar = fig.colorbar(im, ax=[ax_a, ax_b], location="top", shrink=0.55, pad=0.02, aspect=30)
    cbar.set_label("Mean rank (1 = best)")
    cbar.ax.invert_xaxis()

    order = ["split_cell_heldout", "split_scaffold_heldout", "split_cell_scaffold_heldout"]
    x = np.arange(len(order))
    width = 0.34
    colors = {"OpenProblems": "#6f8fb7", "Sci-Plex 3": "#c78f5a"}
    for offset, dataset in [(-width / 2, "OpenProblems"), (width / 2, "Sci-Plex 3")]:
        sub = corr[corr["dataset"].eq(dataset)].set_index("comparison_split").loc[order]
        ax_c.bar(x + offset, sub["mean_spearman"], width=width, color=colors[dataset], label=dataset)
        ax_c.errorbar(
            x + offset,
            sub["mean_spearman"],
            yerr=sub["sd_spearman"],
            fmt="none",
            ecolor="#333333",
            elinewidth=0.8,
            capsize=2,
        )
    ax_c.axhline(0, color="#333333", linewidth=0.8)
    ax_c.set_title("c  Random-rank transfer weakens under strict splits", loc="left", fontweight="bold", pad=6)
    ax_c.set_ylabel("Spearman rho vs random")
    ax_c.set_xticks(x, [SPLIT_LABELS[s] for s in order])
    ax_c.legend(loc="upper right")
    ax_c.set_ylim(-0.9, 0.9)

    x_map = {"split_random": 0, "split_scaffold_heldout": 1, "split_cell_scaffold_heldout": 2}
    for dataset, sub in slopes.groupby("dataset"):
        sub = sub.sort_values("split", key=lambda s: s.map(x_map))
        xs = [x_map[s] for s in sub["split"]]
        ax_d.plot(xs, sub["mean_rank"], marker="o", linewidth=1.6, color=colors[dataset], label=dataset)
        for _, row in sub.iterrows():
            ax_d.text(x_map[row["split"]], row["mean_rank"] + 0.18, f"{row['mean_rank']:.1f}", ha="center", fontsize=6.5)
    ax_d.invert_yaxis()
    ax_d.set_title("d  Random-best drug-aware ridge loses rank", loc="left", fontweight="bold", pad=6)
    ax_d.set_ylabel("Mean rank")
    ax_d.set_xticks([0, 1, 2], [SPLIT_LABELS[s] for s in ["split_random", "split_scaffold_heldout", "split_cell_scaffold_heldout"]])
    ax_d.set_ylim(7.2, 0.5)
    ax_d.legend(loc="lower left")

    for ax in [ax_c, ax_d]:
        ax.grid(axis="y", color="#e5e5e5", linewidth=0.7)

    fig.suptitle("Random validation can select models that do not transfer to strict generalization", x=0.01, ha="left", fontsize=10, fontweight="bold")
    base = OUT / "figure6_ranking_instability"
    fig.savefig(base.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    print(base)


if __name__ == "__main__":
    main()
