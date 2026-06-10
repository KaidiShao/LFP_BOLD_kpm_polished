#!/usr/bin/env python3
"""Pairwise signed ROI and selected-unit tests for roi_mean SOP targets.

This is a control analysis for the focused efun-theta vs deconv-RG test.
It compares all four target combinations:

    efun_theta, efun_rg, deconv_theta, deconv_rg

For each topN it summarizes pairwise ROI profile similarity, top ROI-set
overlap, and selected BOLD unit overlap. P8 selected-mode overlap is
interpretable because it refers to sorted P7 BOLD modes. P10 overlap is kept as
QC only, because efun and deconv_efun P9 components live in separate bases.
"""

from __future__ import annotations

import argparse
import math
from itertools import combinations
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import analyze_roi_mean_signed_target_modes as signed_base


TARGETS = {
    "efun_theta": ("efun", "theta_selective"),
    "efun_rg": ("efun", "ripple_gamma_no_theta"),
    "deconv_theta": ("deconv_efun", "theta_selective"),
    "deconv_rg": ("deconv_efun", "ripple_gamma_no_theta"),
}
TARGET_ORDER = tuple(TARGETS)
ROI_VALUE_MODES = signed_base.ROI_VALUE_MODES
TOP_NS = (3, 5, 10, 20)
TOP_ROIS = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hits-csv",
        type=Path,
        default=Path("results")
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "p8_p10_strict_band_coupling"
        / "p8_p10_strict_band_hits_long.csv",
    )
    parser.add_argument(
        "--p7-roi-csv",
        type=Path,
        default=Path("results")
        / "pipeline_roi_profile_consistency_k03_k16_20260607"
        / "p7_roi_mean_profiles_direct_from_bold_post.csv",
    )
    parser.add_argument("--processed-root", type=Path, default=Path("/mnt/e/DataPons_processed"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results")
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "roi_mean_signed_target_pairwise_tests",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("/mnt/e/DataPons_processed")
        / "summary_figures"
        / "pipeline11_current_analysis_summary"
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "roi_mean_signed_target_pairwise_tests",
    )
    parser.add_argument("--datasets", nargs="+", default=list(signed_base.DATASETS))
    parser.add_argument("--top-ns", nargs="+", type=int, default=list(TOP_NS))
    parser.add_argument("--top-rois", type=int, default=TOP_ROIS)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    signed_base.TARGETS = TARGETS
    datasets = {d.lower() for d in args.datasets}
    roi_order, dataset_roi_order, bold_posts = signed_base.read_p7_metadata(args.p7_roi_csv, datasets)

    p8_cache: dict[str, tuple[dict[int, np.ndarray], dict[int, np.ndarray], str]] = {}
    p10_cache: dict[tuple[str, str, str, str], tuple[dict[int, np.ndarray], str]] = {}

    contrast_rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, object]] = []
    profile_rows: list[dict[str, object]] = []

    for top_n in args.top_ns:
        hits = signed_base.read_hits(args.hits_csv, datasets, top_n)
        unit_sets = signed_base.collect_selected_units(hits)
        overlap_rows.extend(compute_pairwise_unit_overlap(unit_sets, top_n))

        profiles = signed_base.aggregate_target_profiles(
            hits,
            top_n,
            roi_order,
            dataset_roi_order,
            bold_posts,
            args.processed_root,
            p8_cache,
            p10_cache,
        )
        contrast_rows.extend(compute_pairwise_roi_contrasts(profiles, top_n, args.top_rois))

        for (pipeline, dataset, target, mode), payload in profiles.items():
            profile_rows.append(
                {
                    "top_n": top_n,
                    "pipeline": pipeline,
                    "dataset": dataset,
                    "target": target,
                    "roi_value_mode": mode,
                    "weight_sum": payload["weight_sum"],
                    "n_hit_rows": payload["n_hit_rows"],
                }
            )

    signed_base.write_csv(args.output_dir / "pairwise_roi_contrasts.csv", contrast_rows)
    signed_base.write_csv(args.output_dir / "pairwise_selected_unit_overlap.csv", overlap_rows)
    signed_base.write_csv(args.output_dir / "target_profile_availability.csv", profile_rows)

    contrast = pd.DataFrame(contrast_rows)
    overlap = pd.DataFrame(overlap_rows)
    plot_pairwise_heatmaps(contrast, overlap, args.figure_dir)
    write_readme(args.output_dir, args.figure_dir, contrast, overlap)
    print(f"Wrote {args.output_dir}")
    print(f"Figures {args.figure_dir}")


def compute_pairwise_roi_contrasts(
    profiles: dict[tuple[str, str, str, str], dict[str, object]],
    top_n: int,
    top_rois: int,
) -> list[dict[str, object]]:
    rows = []
    contexts = sorted({(pipeline, dataset, mode) for pipeline, dataset, _, mode in profiles})
    for pipeline, dataset, mode in contexts:
        for target_a, target_b in combinations(TARGET_ORDER, 2):
            a_key = (pipeline, dataset, target_a, mode)
            b_key = (pipeline, dataset, target_b, mode)
            if a_key not in profiles or b_key not in profiles:
                continue
            a = np.asarray(profiles[a_key]["profile"], dtype=float)
            b = np.asarray(profiles[b_key]["profile"], dtype=float)
            set_a = signed_base.top_roi_set(a, mode, top_rois)
            set_b = signed_base.top_roi_set(b, mode, top_rois)
            rows.append(
                {
                    "top_n": top_n,
                    "pipeline": pipeline,
                    "dataset": dataset,
                    "roi_value_mode": mode,
                    "target_a": target_a,
                    "target_b": target_b,
                    "profile_corr": signed_base.pearson(a, b),
                    "top_roi_jaccard": signed_base.jaccard(set_a, set_b),
                    "n_top_roi_overlap": len(set_a & set_b),
                    "top_rois": top_rois,
                    "weight_sum_a": profiles[a_key]["weight_sum"],
                    "weight_sum_b": profiles[b_key]["weight_sum"],
                    "n_hit_rows_a": profiles[a_key]["n_hit_rows"],
                    "n_hit_rows_b": profiles[b_key]["n_hit_rows"],
                }
            )
    return rows


def compute_pairwise_unit_overlap(
    unit_sets: dict[tuple[str, str, str], set[tuple]],
    top_n: int,
) -> list[dict[str, object]]:
    rows = []
    contexts = sorted({(pipeline, dataset) for pipeline, dataset, _ in unit_sets})
    for pipeline, dataset in contexts:
        for target_a, target_b in combinations(TARGET_ORDER, 2):
            a = unit_sets.get((pipeline, dataset, target_a), set())
            b = unit_sets.get((pipeline, dataset, target_b), set())
            if not a or not b:
                continue
            inter = a & b
            union = a | b
            adjacent = count_adjacent_units(pipeline, a, b)
            rows.append(
                {
                    "top_n": top_n,
                    "pipeline": pipeline,
                    "dataset": dataset,
                    "target_a": target_a,
                    "target_b": target_b,
                    "n_units_a": len(a),
                    "n_units_b": len(b),
                    "n_exact_overlap": len(inter),
                    "exact_jaccard": len(inter) / len(union) if union else math.nan,
                    "adjacent_or_exact_fraction_of_a": adjacent / len(a) if a else math.nan,
                }
            )
    return rows


def count_adjacent_units(pipeline: str, a: set[tuple], b: set[tuple]) -> int:
    if pipeline == "P8":
        by_run: dict[str, list[int]] = {}
        for _, run, mode in b:
            by_run.setdefault(str(run), []).append(int(mode))
        return sum(
            any(abs(int(mode) - other) <= 1 for other in by_run.get(str(run), []))
            for _, run, mode in a
        )

    by_context: dict[tuple[str, str, str], list[int]] = {}
    for _, run, feature, method_k, comp in b:
        by_context.setdefault((str(run), str(feature), str(method_k)), []).append(int(comp))
    return sum(
        any(abs(int(comp) - other) <= 1 for other in by_context.get((str(run), str(feature), str(method_k)), []))
        for _, run, feature, method_k, comp in a
    )


def plot_pairwise_heatmaps(contrast: pd.DataFrame, overlap: pd.DataFrame, figure_dir: Path) -> None:
    for top_n in sorted(contrast["top_n"].dropna().unique()):
        for pipeline in ["P8", "P10"]:
            for mode in ROI_VALUE_MODES:
                sub = contrast[
                    contrast["top_n"].eq(top_n)
                    & contrast["pipeline"].eq(pipeline)
                    & contrast["roi_value_mode"].eq(mode)
                ]
                if sub.empty:
                    continue
                mat_corr = pairwise_matrix(sub, "profile_corr")
                mat_jacc = pairwise_matrix(sub, "top_roi_jaccard")
                save_matrix(
                    mat_corr,
                    figure_dir / f"01_pairwise_roi_profile_corr__top{int(top_n)}__{pipeline}__{mode}.png",
                    f"{pipeline} ROI profile corr | {mode} | top{int(top_n)}",
                    vmin=-1,
                    vmax=1,
                    cmap="coolwarm",
                )
                save_matrix(
                    mat_jacc,
                    figure_dir / f"02_pairwise_top_roi_jaccard__top{int(top_n)}__{pipeline}__{mode}.png",
                    f"{pipeline} top ROI Jaccard | {mode} | top{int(top_n)}",
                    vmin=0,
                    vmax=1,
                    cmap="viridis",
                )

        sub_o = overlap[overlap["top_n"].eq(top_n) & overlap["pipeline"].eq("P8")]
        if not sub_o.empty:
            save_matrix(
                pairwise_matrix(sub_o, "exact_jaccard"),
                figure_dir / f"03_pairwise_p8_selected_mode_exact_jaccard__top{int(top_n)}.png",
                f"P8 selected BOLD mode exact Jaccard | top{int(top_n)}",
                vmin=0,
                vmax=1,
                cmap="viridis",
            )
            save_matrix(
                pairwise_matrix(sub_o, "adjacent_or_exact_fraction_of_a"),
                figure_dir / f"04_pairwise_p8_selected_mode_adjacent_fraction__top{int(top_n)}.png",
                f"P8 selected BOLD mode same/adjacent fraction | top{int(top_n)}",
                vmin=0,
                vmax=1,
                cmap="viridis",
            )


def pairwise_matrix(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    mat = pd.DataFrame(np.eye(len(TARGET_ORDER)), index=TARGET_ORDER, columns=TARGET_ORDER, dtype=float)
    for target_a, target_b in combinations(TARGET_ORDER, 2):
        vals = df[df["target_a"].eq(target_a) & df["target_b"].eq(target_b)][value_col]
        vals = pd.to_numeric(vals, errors="coerce").dropna()
        value = float(vals.median()) if not vals.empty else np.nan
        mat.loc[target_a, target_b] = value
        mat.loc[target_b, target_a] = value
    return mat


def save_matrix(
    mat: pd.DataFrame,
    path: Path,
    title: str,
    vmin: float,
    vmax: float,
    cmap: str,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(mat.values, vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_yticks(range(len(mat.index)))
    ax.set_xticklabels(mat.columns, rotation=35, ha="right")
    ax.set_yticklabels(mat.index)
    ax.set_title(title)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            value = mat.iat[i, j]
            text = "NA" if not np.isfinite(value) else f"{value:.2f}"
            ax.text(j, i, text, ha="center", va="center", color="black", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_readme(output_dir: Path, figure_dir: Path, contrast: pd.DataFrame, overlap: pd.DataFrame) -> None:
    lines = [
        "# ROI mean signed target pairwise tests",
        "",
        "Targets: efun_theta, efun_rg, deconv_theta, deconv_rg.",
        "",
        "P8 selected-mode overlap refers to sorted P7 BOLD modes and is interpretable.",
        "P10 selected-unit overlap is QC only because efun and deconv_efun are separate P9 bases.",
        "",
        f"Figures: `{figure_dir}`",
        "",
    ]
    if not contrast.empty:
        lines.append("## Median pairwise ROI profile correlation")
        for top_n in sorted(contrast["top_n"].unique()):
            sub = contrast[contrast["top_n"].eq(top_n)]
            lines.append(f"top{int(top_n)}")
            med = (
                sub.groupby(["pipeline", "roi_value_mode", "target_a", "target_b"], dropna=False)["profile_corr"]
                .median()
                .reset_index()
                .sort_values(["pipeline", "roi_value_mode", "profile_corr"])
            )
            for _, row in med.iterrows():
                lines.append(
                    f"- {row['pipeline']} {row['roi_value_mode']} {row['target_a']} vs {row['target_b']}: "
                    f"{float(row['profile_corr']):.3f}"
                )
            lines.append("")
    output_dir.joinpath("README_roi_mean_signed_target_pairwise_tests.md").write_text(
        "\n".join(lines),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
