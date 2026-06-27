#!/usr/bin/env python3
"""Inspect AnnData .h5ad observation and variable metadata without running full analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def summarize_h5ad(path: Path, outdir: Path) -> None:
    import anndata as ad

    outdir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(path, backed="r")

    obs = adata.obs.copy()
    var = adata.var.copy()
    obs.to_csv(outdir / f"{path.stem}_obs_head.csv", index=True)
    var.head(200).to_csv(outdir / f"{path.stem}_var_head.csv", index=True)

    summary_rows = []
    for col in obs.columns:
        series = obs[col]
        n_unique = series.nunique(dropna=True)
        examples = "; ".join(map(str, series.dropna().astype(str).unique()[:8]))
        summary_rows.append(
            {
                "column": col,
                "dtype": str(series.dtype),
                "n_unique": int(n_unique),
                "examples": examples,
            }
        )

    summary = pd.DataFrame(summary_rows).sort_values("n_unique", ascending=False)
    summary.to_csv(outdir / f"{path.stem}_obs_column_summary.csv", index=False)

    shape = pd.DataFrame(
        [
            {
                "file": path.name,
                "n_obs": adata.n_obs,
                "n_vars": adata.n_vars,
                "n_obs_columns": len(obs.columns),
                "n_var_columns": len(var.columns),
            }
        ]
    )
    shape.to_csv(outdir / f"{path.stem}_shape.csv", index=False)

    print(shape.to_string(index=False))
    print(summary.head(40).to_string(index=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("--outdir", type=Path, default=Path("results/h5ad_metadata"))
    args = parser.parse_args()
    summarize_h5ad(args.h5ad, args.outdir)


if __name__ == "__main__":
    main()

