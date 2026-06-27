"""Summarize and plot TranSiGen input ablations."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCI = ROOT / "results/deep_model_panel/sciplex3"
FIG = ROOT / "figures"
SRC = ROOT / "results/deep_model_panel/source_data"
FIG.mkdir(exist_ok=True)
SRC.mkdir(parents=True, exist_ok=True)


SPLIT_LABELS = {
    "random": "Random",
    "cell_heldout": "Cell-line",
    "scaffold_heldout": "Scaffold",
    "cell_scaffold_heldout": "Joint",
}
MODEL_LABELS = {
    "transigen_basal_only": "Basal only",
    "transigen_drug_only": "Drug only",
    "transigen_basal_drug": "Basal+drug",
    "transigen_basal_drug_dose": "Basal+drug+dose",
    "transigen_basal_drug_dose_cell": "Basal+drug+dose+cell",
    "transigen_full_current": "Full current",
}


def ci(values: np.ndarray, n_boot: int = 10000, seed: int = 1729) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    boot = np.array([rng.choice(values, size=len(values), replace=True).mean() for _ in range(n_boot)])
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def main() -> None:
    df = pd.read_csv(SCI / "transigen_ablation_metrics_long.csv")
    ok = df[df["status"].eq("completed")].copy()
    summary = (
        ok.groupby(["model_name", "split_type"], as_index=False)
        .agg(
            seeds=("seed", "nunique"),
            mean_pearson=("mean_all_gene_pearson", "mean"),
            sd_pearson=("mean_all_gene_pearson", "std"),
            median_pearson=("mean_all_gene_pearson", "median"),
            mean_rmse=("mean_all_gene_rmse", "mean"),
            sd_rmse=("mean_all_gene_rmse", "std"),
            input_dim=("input_dim", "first"),
        )
    )
    summary["model_label"] = summary["model_name"].map(MODEL_LABELS)
    summary["split_label"] = summary["split_type"].map(SPLIT_LABELS)
    summary.to_csv(SRC / "figure6_transigen_ablation_sciplex3.csv", index=False)

    rows = []
    for model, g in ok.groupby("model_name"):
        wide = g.pivot(index="seed", columns="split_type", values="mean_all_gene_pearson")
        for a, b in [
            ("random", "scaffold_heldout"),
            ("random", "cell_scaffold_heldout"),
            ("scaffold_heldout", "cell_scaffold_heldout"),
        ]:
            if a in wide and b in wide:
                diff = wide[a] - wide[b]
                lo, hi = ci(diff.to_numpy())
                rows.append(
                    {
                        "model_name": model,
                        "contrast": f"{a}_minus_{b}",
                        "n_seeds": int(diff.notna().sum()),
                        "mean_difference": float(diff.mean()),
                        "sd_difference": float(diff.std(ddof=1)),
                        "ci_low": lo,
                        "ci_high": hi,
                    }
                )
    contrasts = pd.DataFrame(rows)
    contrasts.to_csv(SCI / "transigen_ablation_contrasts.csv", index=False)

    rank = summary.copy()
    rank["rank"] = rank.groupby("split_type")["mean_pearson"].rank(ascending=False, method="average")
    rank_table = rank.pivot_table(index="model_name", columns="split_type", values="rank", aggfunc="mean")
    rank_table.to_csv(SCI / "transigen_ablation_rank_table.csv")

    plot = summary[summary["split_type"].isin(["random", "scaffold_heldout", "cell_scaffold_heldout"])].copy()
    pivot = plot.pivot_table(index="model_label", columns="split_label", values="mean_pearson", aggfunc="mean")
    pivot = pivot[[c for c in ["Random", "Scaffold", "Joint"] if c in pivot.columns]]
    fig, ax = plt.subplots(figsize=(9, 4.8))
    pivot.plot(kind="bar", ax=ax, color=["#4C78A8", "#F58518", "#54A24B"])
    ax.set_ylabel("All-gene row-wise Pearson")
    ax.set_xlabel("")
    ax.set_title("TranSiGen-style Sci-Plex 3 input ablation")
    ax.legend(title="Split")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(FIG / "figure6_transigen_ablation_sciplex3.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "figure6_transigen_ablation_sciplex3.pdf", bbox_inches="tight")
    plt.close(fig)

    ok.to_csv(SRC / "supp_transigen_seed_distribution.csv", index=False)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    order = ["random", "cell_heldout", "scaffold_heldout", "cell_scaffold_heldout"]
    data = [ok.loc[ok["split_type"].eq(s), "mean_all_gene_pearson"].to_numpy() for s in order]
    ax.boxplot(data, tick_labels=[SPLIT_LABELS[s] for s in order], showfliers=False)
    ax.set_ylabel("All-gene row-wise Pearson")
    ax.set_title("TranSiGen ablation seed-level distributions")
    fig.tight_layout()
    fig.savefig(FIG / "supp_transigen_seed_distribution.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "supp_transigen_seed_distribution.pdf", bbox_inches="tight")
    plt.close(fig)

    contrasts.to_csv(SRC / "supp_transigen_random_to_strict_ci.csv", index=False)
    focus = contrasts[contrasts["contrast"].eq("random_minus_cell_scaffold_heldout")].copy()
    focus["model_label"] = focus["model_name"].map(MODEL_LABELS)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.errorbar(
        focus["mean_difference"],
        focus["model_label"],
        xerr=[focus["mean_difference"] - focus["ci_low"], focus["ci_high"] - focus["mean_difference"]],
        fmt="o",
        color="#4C78A8",
        ecolor="#9ecae9",
        capsize=3,
    )
    ax.set_xlabel("Random minus joint Pearson")
    ax.set_title("Random-to-joint drop by input ablation")
    fig.tight_layout()
    fig.savefig(FIG / "supp_transigen_random_to_strict_ci.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "supp_transigen_random_to_strict_ci.pdf", bbox_inches="tight")
    plt.close(fig)

    pair = pd.read_csv(SCI / "transigen_vs_baseline_paired_differences.csv")
    pair.to_csv(SRC / "supp_transigen_vs_ridge_paired.csv", index=False)
    pair_focus = pair[pair["baseline"].isin(["ridge_cell_dose_drug_fp", "ridge_drug_fp", "global_train_mean"])].copy()
    pair_focus["label"] = pair_focus["baseline"] + " / " + pair_focus["split_key"]
    fig, ax = plt.subplots(figsize=(8, 5.2))
    ax.errorbar(
        pair_focus["mean_difference_transigen_minus_baseline"],
        pair_focus["label"],
        xerr=[
            pair_focus["mean_difference_transigen_minus_baseline"] - pair_focus["ci_low"],
            pair_focus["ci_high"] - pair_focus["mean_difference_transigen_minus_baseline"],
        ],
        fmt="o",
        color="#B279A2",
        ecolor="#d6b1ca",
        capsize=3,
    )
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("TranSiGen minus baseline Pearson")
    ax.set_title("Paired comparison against Sci-Plex 3 baselines")
    fig.tight_layout()
    fig.savefig(FIG / "supp_transigen_vs_ridge_paired.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG / "supp_transigen_vs_ridge_paired.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
