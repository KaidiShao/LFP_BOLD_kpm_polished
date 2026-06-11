#!/usr/bin/env python3
"""Plot SOP figures for P5 theta/RG density correlation.

Main figure:
    dataset x method-k heatmap of canonical theta/RG session-demeaned
    density correlation.

Optional supplement figures:
    canonical-pair distribution by dataset and all-pair distribution.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WORKSPACE = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = (
    WORKSPACE
    / "results"
    / "standardized_csplit_k03_k16_all_current_20260607"
    / "p5_theta_rg_density_correlation"
)
DEFAULT_FIGURE_DIR = DEFAULT_INPUT_DIR / "sop_figures"

METHOD_ORDER = ["svd", "nmf", "mds", "umap"]
DEFAULT_KS = list(range(3, 17))
DEFAULT_LOW_ABS_THRESHOLD = 0.2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=None,
        help="Default: <input-dir>/sop_figures.",
    )
    parser.add_argument(
        "--metric",
        choices=["session_demeaned_density_corr", "density_corr"],
        default="session_demeaned_density_corr",
        help="Correlation metric to plot. Session-demeaned is the SOP default.",
    )
    parser.add_argument(
        "--summary-threshold",
        type=float,
        default=DEFAULT_LOW_ABS_THRESHOLD,
        help="Threshold for the per-dataset fraction of low absolute correlation.",
    )
    parser.add_argument(
        "--vlim",
        type=float,
        default=0.6,
        help="Symmetric color limit. Values outside are clipped visually only.",
    )
    parser.add_argument("--include-cell-text", action="store_true")
    parser.add_argument(
        "--include-supplements",
        action="store_true",
        help="Also write optional distribution plots. Default SOP writes only the main heatmap.",
    )
    parser.add_argument("--png-dpi", type=int, default=220)
    return parser.parse_args()


def natural_dataset_order(values: pd.Series) -> list[str]:
    preferred = [
        "e10gb1",
        "e10fV1",
        "e10gh1",
        "e10gw1",
        "f12m01",
        "f12m02",
        "f12m03",
        "f12m05",
        "k13m17",
        "k13m18",
        "k13m19",
        "k13m20",
        "k13m21",
        "k13m23",
    ]
    seen = set(str(v) for v in values.dropna().unique())
    ordered = [ds for ds in preferred if ds in seen]
    ordered.extend(sorted(seen.difference(ordered)))
    return ordered


def method_k_order(df: pd.DataFrame) -> list[str]:
    observed = set(str(v) for v in df["method_k"].dropna().unique())
    ordered = [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in DEFAULT_KS]
    ordered = [mk for mk in ordered if mk in observed]
    extras = sorted(observed.difference(ordered))
    return ordered + extras


def read_inputs(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    canonical_path = input_dir / "theta_rg_density_canonical_pair_correlations.csv"
    all_pairs_path = input_dir / "theta_rg_density_all_pair_correlations.csv"
    if not canonical_path.exists():
        raise FileNotFoundError(canonical_path)
    canonical = pd.read_csv(canonical_path)
    all_pairs = pd.read_csv(all_pairs_path) if all_pairs_path.exists() else None
    return canonical, all_pairs


def prepare_pivot(canonical: pd.DataFrame, metric: str) -> pd.DataFrame:
    required = {"dataset", "method_k", metric}
    missing = required.difference(canonical.columns)
    if missing:
        raise ValueError(f"Missing required columns in canonical CSV: {sorted(missing)}")

    ds_order = natural_dataset_order(canonical["dataset"])
    mk_order = method_k_order(canonical)
    pivot = canonical.pivot_table(
        index="dataset",
        columns="method_k",
        values=metric,
        aggfunc="median",
        observed=False,
    )
    pivot = pivot.reindex(index=ds_order, columns=mk_order)
    return pivot


def draw_method_separators(ax: plt.Axes, columns: list[str]) -> None:
    last_method = None
    for idx, col in enumerate(columns):
        method = str(col).split("_k", 1)[0]
        if last_method is not None and method != last_method:
            ax.axvline(idx - 0.5, color="#222222", lw=1.0)
        last_method = method


def plot_main_heatmap(
    canonical: pd.DataFrame,
    figure_dir: Path,
    metric: str,
    low_abs_threshold: float,
    vlim: float,
    include_cell_text: bool,
    dpi: int,
) -> Path:
    pivot = prepare_pivot(canonical, metric)
    data = pivot.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(data)

    cmap = plt.get_cmap("coolwarm").copy()
    cmap.set_bad("#e6e6e6")

    n_rows, n_cols = data.shape
    fig_w = max(16, 0.32 * n_cols + 3.5)
    fig_h = max(5.5, 0.42 * n_rows + 2.2)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=-vlim, vmax=vlim)
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels(pivot.index, fontsize=9)
    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels(pivot.columns, rotation=90, fontsize=7)
    ax.set_xlabel("P5 dimred method-k")
    ax.set_ylabel("dataset")
    ax.set_title(
        "P5 canonical theta/RG density correlation | "
        f"{metric.replace('_', ' ')} | near 0 = independent subprocesses"
    )

    draw_method_separators(ax, list(pivot.columns))
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.45)
    ax.tick_params(which="minor", bottom=False, left=False)

    if include_cell_text:
        for row_idx in range(n_rows):
            for col_idx in range(n_cols):
                val = data[row_idx, col_idx]
                if math.isfinite(val):
                    ax.text(
                        col_idx,
                        row_idx,
                        f"{val:+.2f}",
                        ha="center",
                        va="center",
                        fontsize=5.5,
                        color="black",
                    )

    # Add per-dataset low-correlation fraction as a compact side annotation.
    for row_idx, ds in enumerate(pivot.index):
        vals = pd.to_numeric(pivot.loc[ds], errors="coerce").dropna()
        if vals.empty:
            label = "n=0"
        else:
            frac = float((vals.abs() < low_abs_threshold).mean())
            label = f"{frac:.0%} |r|<{low_abs_threshold:g}"
        ax.text(
            n_cols + 0.35,
            row_idx,
            label,
            va="center",
            ha="left",
            fontsize=8,
            color="#333333",
        )
    ax.set_xlim(-0.5, n_cols + 2.8)

    cbar = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.015)
    cbar.set_label("canonical theta/RG density correlation")

    figure_dir.mkdir(parents=True, exist_ok=True)
    out = figure_dir / f"01_main_dataset_methodk_canonical_theta_rg_density_corr_heatmap__{metric}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_canonical_distribution(
    canonical: pd.DataFrame,
    figure_dir: Path,
    metric: str,
    low_abs_threshold: float,
    dpi: int,
) -> Path:
    ds_order = natural_dataset_order(canonical["dataset"])
    values = [
        pd.to_numeric(canonical.loc[canonical["dataset"].eq(ds), metric], errors="coerce").dropna().to_numpy()
        for ds in ds_order
    ]

    fig, ax = plt.subplots(figsize=(12, max(5, 0.42 * len(ds_order) + 1.2)))
    parts = ax.violinplot(values, vert=False, showmeans=False, showextrema=False, showmedians=True)
    for body in parts["bodies"]:
        body.set_facecolor("#80b1d3")
        body.set_edgecolor("#4c78a8")
        body.set_alpha(0.55)
    if "cmedians" in parts:
        parts["cmedians"].set_color("#1f1f1f")
        parts["cmedians"].set_linewidth(1.2)

    for idx, vals in enumerate(values, start=1):
        if vals.size:
            jitter = np.linspace(-0.08, 0.08, vals.size)
            ax.scatter(vals, np.full(vals.size, idx) + jitter, s=12, alpha=0.45, color="#4c78a8")

    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(low_abs_threshold, color="#777777", lw=0.8, ls="--")
    ax.axvline(-low_abs_threshold, color="#777777", lw=0.8, ls="--")
    ax.set_xlim(-1, 1)
    ax.set_yticks(np.arange(1, len(ds_order) + 1))
    ax.set_yticklabels(ds_order)
    ax.set_xlabel("canonical theta/RG density correlation")
    ax.set_ylabel("dataset")
    ax.set_title("Supplement: canonical theta/RG density correlation distribution by dataset")
    ax.grid(axis="x", alpha=0.25)

    out = figure_dir / f"02_supplement_canonical_theta_rg_density_corr_distribution_by_dataset__{metric}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def plot_all_pair_distribution(
    all_pairs: pd.DataFrame | None,
    figure_dir: Path,
    metric: str,
    dpi: int,
) -> Path | None:
    if all_pairs is None or all_pairs.empty or metric not in all_pairs.columns:
        return None
    vals = pd.to_numeric(all_pairs[metric], errors="coerce").dropna()
    if vals.empty:
        return None

    fig, ax = plt.subplots(figsize=(10, 4.5))
    bins = np.linspace(-1, 1, 45)
    ax.hist(vals, bins=bins, color="#4c78a8", alpha=0.75)
    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(vals.median(), color="#d62728", lw=1.2, label=f"median={vals.median():+.3f}")
    ax.set_xlabel("all theta/RG pair density correlation")
    ax.set_ylabel("component pairs")
    ax.set_title("Supplement: all theta/RG component-pair density correlation")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    out = figure_dir / f"03_supplement_all_pair_theta_rg_density_corr_distribution__{metric}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    return out


def write_readme(paths: list[Path], args: argparse.Namespace) -> None:
    lines = [
        "# P5 theta/RG density correlation SOP figures",
        "",
        "Main figure:",
        "",
        "- `01_main_dataset_methodk_canonical_theta_rg_density_corr_heatmap__*.png`",
        "",
        "Interpretation:",
        "",
        "```text",
        "Each cell is one dataset x method-k canonical theta/RG density correlation.",
        "The SOP metric is session-demeaned Pearson correlation.",
        "Near 0 means theta-like and RG-no-theta densities are relatively independent.",
        "Large positive or negative magnitude means the two labeled subprocess densities",
        "may not be well separated in time for that dataset/method-k.",
        "```",
        "",
        "Parameters:",
        "",
        "```text",
        f"input_dir = {args.input_dir}",
        f"metric = {args.metric}",
        f"summary_threshold = {args.summary_threshold}",
        f"vlim = {args.vlim}",
        "```",
        "",
        "Generated files:",
        "",
    ]
    if args.include_supplements:
        lines[lines.index("Parameters:"):lines.index("Parameters:")] = [
            "Supplement figures:",
            "",
            "- `02_supplement_canonical_theta_rg_density_corr_distribution_by_dataset__*.png`",
            "- `03_supplement_all_pair_theta_rg_density_corr_distribution__*.png`",
            "",
        ]
    for path in paths:
        lines.append(f"- `{path.name}`")
    (args.figure_dir / "README_p5_theta_rg_density_correlation_sop_figures.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def main() -> None:
    args = parse_args()
    args.input_dir = args.input_dir.resolve()
    if args.figure_dir is None:
        args.figure_dir = args.input_dir / "sop_figures"
    args.figure_dir = args.figure_dir.resolve()
    canonical, all_pairs = read_inputs(args.input_dir)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    paths: list[Path] = []
    paths.append(
        plot_main_heatmap(
            canonical,
            args.figure_dir,
            args.metric,
            args.summary_threshold,
            args.vlim,
            args.include_cell_text,
            args.png_dpi,
        )
    )
    if args.include_supplements:
        paths.append(
            plot_canonical_distribution(
                canonical,
                args.figure_dir,
                args.metric,
                args.summary_threshold,
                args.png_dpi,
            )
        )
        all_pair_path = plot_all_pair_distribution(all_pairs, args.figure_dir, args.metric, args.png_dpi)
        if all_pair_path is not None:
            paths.append(all_pair_path)
    write_readme(paths, args)

    print("Wrote:")
    for path in paths:
        print(path)


if __name__ == "__main__":
    main()
