"""Focused summary panel for P8 raw_csplit_q070 efun/deconv ROI checks."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from summarize_pipeline8_cross_session_consistency import DEFAULT_PROCESSED_ROOT


DEFAULT_RESULTS_DIR = Path("results") / "pipeline_roi_profile_consistency_current"
DEFAULT_SOURCE_DIR = DEFAULT_RESULTS_DIR / "p8_raw_density_efun_vs_deconv_roi_similarity_20260603"
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p8_raw_csplit_q070_focus_20260603"
)
OBS_ORDER = ("global_svd100", "gsvd100_ds", "roi_mean")
DATASET_ORDER = ("e10gb1", "e10fv1", "e10gh1", "e10gw1", "f12m01")
COMPARISON_ORDER = (
    "same_dataset_efun_vs_deconv",
    "cross_dataset_efun_within",
    "cross_dataset_deconv_efun_within",
    "cross_dataset_efun_vs_deconv",
)
COMPARISON_LABELS = (
    "same ds\nE vs D",
    "cross ds\nE within",
    "cross ds\nD within",
    "cross ds\nE vs D",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--raw-density", default="raw_csplit_q070")
    return parser.parse_args()


def as_float(value: object) -> float:
    try:
        out = float(value) if value not in (None, "") else math.nan
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def safe_name(text: object) -> str:
    return "".join(ch if ch.isalnum() or ch in ("-", "_") else "_" for ch in str(text)).strip("_") or "blank"


def plot_focus_panel(source_dir: Path, figure_dir: Path, raw_density: str) -> Path:
    roi_summary = [
        row
        for row in read_csv(source_dir / "p8_raw_density_efun_vs_deconv_roi_summary.csv")
        if row.get("raw_density") == raw_density
    ]
    roi_pairwise = [
        row
        for row in read_csv(source_dir / "p8_raw_density_efun_vs_deconv_roi_pairwise.csv")
        if row.get("raw_density") == raw_density
    ]
    index_summary = [
        row
        for row in read_csv(source_dir / "p8_raw_density_efun_vs_deconv_mode_index_overlap_summary.csv")
        if row.get("raw_density") == raw_density
    ]

    summary_lookup = {
        (row.get("observable", ""), row.get("comparison", "")): as_float(row.get("mean"))
        for row in roi_summary
    }
    pair_lookup = {
        (row.get("observable", ""), row.get("dataset_a", "").lower()): as_float(row.get("roi_corr"))
        for row in roi_pairwise
        if row.get("comparison") == "same_dataset_efun_vs_deconv"
    }
    index_lookup = {
        (row.get("observable", ""), row.get("metric", "")): as_float(row.get("mean"))
        for row in index_summary
    }

    matrix = np.full((len(OBS_ORDER), len(COMPARISON_ORDER)), np.nan, dtype=float)
    for i, obs in enumerate(OBS_ORDER):
        for j, comp in enumerate(COMPARISON_ORDER):
            matrix[i, j] = summary_lookup.get((obs, comp), math.nan)

    dataset_matrix = np.full((len(OBS_ORDER), len(DATASET_ORDER)), np.nan, dtype=float)
    for i, obs in enumerate(OBS_ORDER):
        for j, dataset in enumerate(DATASET_ORDER):
            dataset_matrix[i, j] = pair_lookup.get((obs, dataset), math.nan)

    overlap_matrix = np.full((len(OBS_ORDER), 2), np.nan, dtype=float)
    for i, obs in enumerate(OBS_ORDER):
        overlap_matrix[i, 0] = index_lookup.get((obs, "unique_jaccard"), math.nan)
        overlap_matrix[i, 1] = index_lookup.get((obs, "weighted_cosine"), math.nan)

    fig = plt.figure(figsize=(14.8, 8.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1.0], height_ratios=[1.0, 1.0], hspace=0.35, wspace=0.3)
    ax_summary = fig.add_subplot(gs[:, 0])
    ax_dataset = fig.add_subplot(gs[0, 1])
    ax_overlap = fig.add_subplot(gs[1, 1])

    im = ax_summary.imshow(matrix, cmap="coolwarm", vmin=0.0, vmax=1.0, aspect="auto")
    ax_summary.set_xticks(range(len(COMPARISON_ORDER)))
    ax_summary.set_xticklabels(COMPARISON_LABELS)
    ax_summary.set_yticks(range(len(OBS_ORDER)))
    ax_summary.set_yticklabels(OBS_ORDER)
    ax_summary.set_title("ROI profile corr means")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if math.isfinite(float(matrix[i, j])):
                ax_summary.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=11)
    cbar = fig.colorbar(im, ax=ax_summary, shrink=0.72)
    cbar.set_label("Pearson corr")

    im2 = ax_dataset.imshow(dataset_matrix, cmap="coolwarm", vmin=0.0, vmax=1.0, aspect="auto")
    ax_dataset.set_xticks(range(len(DATASET_ORDER)))
    ax_dataset.set_xticklabels(DATASET_ORDER, rotation=35, ha="right")
    ax_dataset.set_yticks(range(len(OBS_ORDER)))
    ax_dataset.set_yticklabels(OBS_ORDER)
    ax_dataset.set_title("Same-dataset E vs D ROI corr")
    for i in range(dataset_matrix.shape[0]):
        for j in range(dataset_matrix.shape[1]):
            if math.isfinite(float(dataset_matrix[i, j])):
                ax_dataset.text(j, i, f"{dataset_matrix[i, j]:.2f}", ha="center", va="center", fontsize=10)
            else:
                ax_dataset.text(j, i, "NA", ha="center", va="center", fontsize=9, color="#606060")
    fig.colorbar(im2, ax=ax_dataset, shrink=0.8)

    im3 = ax_overlap.imshow(overlap_matrix, cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
    ax_overlap.set_xticks([0, 1])
    ax_overlap.set_xticklabels(["mode-index\nJaccard", "weighted-index\ncosine"])
    ax_overlap.set_yticks(range(len(OBS_ORDER)))
    ax_overlap.set_yticklabels(OBS_ORDER)
    ax_overlap.set_title("Selected BOLD mode-index overlap")
    for i in range(overlap_matrix.shape[0]):
        for j in range(overlap_matrix.shape[1]):
            if math.isfinite(float(overlap_matrix[i, j])):
                color = "white" if overlap_matrix[i, j] < 0.55 else "black"
                ax_overlap.text(j, i, f"{overlap_matrix[i, j]:.2f}", ha="center", va="center", fontsize=10, color=color)
    fig.colorbar(im3, ax=ax_overlap, shrink=0.8)

    fig.suptitle(
        f"P8 {raw_density}: efun_real vs deconv_real top-xcorr ROI footprint",
        fontsize=16,
        weight="bold",
    )
    fig.text(
        0.5,
        0.025,
        "Interpretation: low same-dataset E-vs-D corr plus low mode overlap suggests different ROI/mode footprints; high values suggest the same BOLD spatial subspace.",
        ha="center",
        fontsize=10,
        color="#404040",
    )
    figure_dir.mkdir(parents=True, exist_ok=True)
    out_file = figure_dir / f"p8_{safe_name(raw_density)}_efun_vs_deconv_roi_focus_panel.png"
    fig.tight_layout(rect=[0, 0.05, 1, 0.94])
    fig.savefig(out_file, dpi=180)
    plt.close(fig)
    return out_file


def main() -> None:
    args = parse_args()
    out_file = plot_focus_panel(args.source_dir, args.figure_dir, args.raw_density)
    print(f"Figure: {out_file}")


if __name__ == "__main__":
    main()
