"""Plot current P8/P10 cross-session ROI profile consistency.

Inputs are ROI-vector CSVs exported by
``export_p8_p10_roi_profile_consistency_sources.m``.  Each consistency score is
the mean off-diagonal Pearson correlation of the same aggregate ROI profile
across datasets.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Iterable, Sequence

from summarize_pipeline8_cross_session_consistency import (
    COMPONENT_COUNTS,
    DEFAULT_PROCESSED_ROOT,
    METHOD_ORDER,
    RUN_TAG_LABELS,
    plot_heatmap,
    slug,
    write_csv,
)


DEFAULT_RESULTS_DIR = Path("results") / "pipeline_roi_profile_consistency_current"
DEFAULT_P8_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p8"
    / "roi_profile_consistency"
)
DEFAULT_P10_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p10"
    / "roi_profile_consistency"
)
DEFAULT_RESULTS_FIGURE_DIR = DEFAULT_RESULTS_DIR / "figures"
DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "f12m01")
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--p8-profile-csv", type=Path, default=DEFAULT_RESULTS_DIR / "p8_roi_profiles_long.csv")
    parser.add_argument("--p10-profile-csv", type=Path, default=DEFAULT_RESULTS_DIR / "p10_roi_profiles_long.csv")
    parser.add_argument("--p8-figure-dir", type=Path, default=DEFAULT_P8_FIGURE_DIR)
    parser.add_argument("--p10-figure-dir", type=Path, default=DEFAULT_P10_FIGURE_DIR)
    parser.add_argument("--results-figure-dir", type=Path, default=DEFAULT_RESULTS_FIGURE_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--value-column", default="roi_value_weighted")
    parser.add_argument("--skip-figures", action="store_true")
    return parser.parse_args()


def as_float(value: str | None) -> float:
    try:
        out = float(value) if value not in (None, "") else math.nan
    except ValueError:
        return math.nan
    return out if math.isfinite(out) else math.nan


def finite(values: Iterable[float]) -> list[float]:
    return [v for v in values if math.isfinite(v)]


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return math.nan
    xvals = [p[0] for p in pairs]
    yvals = [p[1] for p in pairs]
    xmean = mean(xvals)
    ymean = mean(yvals)
    num = sum((x - xmean) * (y - ymean) for x, y in pairs)
    xden = math.sqrt(sum((x - xmean) ** 2 for x in xvals))
    yden = math.sqrt(sum((y - ymean) ** 2 for y in yvals))
    if xden == 0 or yden == 0:
        return math.nan
    return num / (xden * yden)


def method_k_order() -> list[str]:
    return [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]


def read_profile_groups(
    path: Path,
    pipeline: str,
    datasets: Sequence[str],
    run_tags: Sequence[str],
    value_column: str,
) -> tuple[dict[tuple[str, ...], dict[str, dict[str, float]]], dict[tuple[str, ...], dict[str, str]]]:
    dataset_set = {d.lower() for d in datasets}
    run_tag_set = {r.lower() for r in run_tags}
    groups: dict[tuple[str, ...], dict[str, dict[str, float]]] = defaultdict(lambda: defaultdict(dict))
    meta: dict[tuple[str, ...], dict[str, str]] = {}
    if not path.is_file():
        return groups, meta

    with path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("pipeline", "").upper() != pipeline.upper():
                continue
            dataset = row.get("dataset", "")
            run_tag = row.get("run_tag", "")
            if dataset.lower() not in dataset_set or run_tag.lower() not in run_tag_set:
                continue
            value = as_float(row.get(value_column))
            if not math.isfinite(value):
                continue
            roi_label = row.get("roi_label", "")
            if not roi_label:
                continue
            if pipeline.upper() == "P8":
                key = (
                    row.get("run_tag", ""),
                    row.get("observable", ""),
                    row.get("feature_family", ""),
                    row.get("bold_feature", ""),
                    row.get("density_condition", ""),
                    row.get("density_method_k", ""),
                    row.get("density_name", ""),
                )
            else:
                key = (
                    row.get("run_tag", ""),
                    row.get("observable", ""),
                    row.get("p9_feature", ""),
                    row.get("p9_method_k", ""),
                    row.get("feature_family", ""),
                    row.get("bold_feature", ""),
                    row.get("density_condition", ""),
                    row.get("density_method_k", ""),
                    row.get("density_name", ""),
                )
            groups[key][dataset][roi_label] = value
            if key not in meta:
                meta[key] = {
                    "pipeline": pipeline.upper(),
                    "run_tag": row.get("run_tag", ""),
                    "observable": row.get("observable", ""),
                    "p9_feature": row.get("p9_feature", ""),
                    "p9_method": row.get("p9_method", ""),
                    "p9_k": row.get("p9_k", ""),
                    "p9_method_k": row.get("p9_method_k", ""),
                    "feature_family": row.get("feature_family", ""),
                    "bold_feature": row.get("bold_feature", ""),
                    "density_condition": row.get("density_condition", ""),
                    "density_method": row.get("density_method", ""),
                    "density_k": row.get("density_k", ""),
                    "density_method_k": row.get("density_method_k", ""),
                    "density_name": row.get("density_name", ""),
                    "density_display": row.get("density_display", ""),
                    "roi_value_mode": row.get("roi_value_mode", ""),
                    "n_selected": row.get("n_selected", ""),
                }
    return groups, meta


def summarize_group_consistency(
    groups: dict[tuple[str, ...], dict[str, dict[str, float]]],
    meta: dict[tuple[str, ...], dict[str, str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key, dataset_vectors in groups.items():
        datasets = sorted(dataset_vectors)
        if len(datasets) < 2:
            continue
        common_rois = sorted(set.intersection(*(set(dataset_vectors[d]) for d in datasets)))
        if len(common_rois) < 2:
            continue
        corrs = []
        pair_text = []
        for i, d1 in enumerate(datasets):
            v1 = [dataset_vectors[d1][roi] for roi in common_rois]
            for d2 in datasets[i + 1 :]:
                v2 = [dataset_vectors[d2][roi] for roi in common_rois]
                r = pearson(v1, v2)
                if math.isfinite(r):
                    corrs.append(r)
                    pair_text.append(f"{d1}-{d2}:{r:.4f}")
        if not corrs:
            continue
        info = dict(meta.get(key, {}))
        info.update(
            {
                "n_datasets": len(datasets),
                "datasets": ";".join(datasets),
                "n_rois_common": len(common_rois),
                "mean_pairwise_roi_corr": f"{mean(corrs):.10g}",
                "min_pairwise_roi_corr": f"{min(corrs):.10g}",
                "max_pairwise_roi_corr": f"{max(corrs):.10g}",
                "pairwise_roi_corrs": ";".join(pair_text),
            }
        )
        rows.append(info)
    rows.sort(key=lambda r: (r.get("pipeline", ""), r.get("run_tag", ""), r.get("p9_method_k", ""), r.get("density_method_k", ""), r.get("feature_family", "")))
    return rows


def average_matrix(rows: Sequence[dict[str, object]], row_field: str, col_field: str, row_labels: Sequence[str], col_labels: Sequence[str]) -> dict[tuple[str, str], float]:
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        rlabel = str(row.get(row_field, ""))
        clabel = str(row.get(col_field, ""))
        value = as_float(str(row.get("mean_pairwise_roi_corr", "")))
        if math.isfinite(value):
            grouped[(rlabel, clabel)].append(value)
    return {key: mean(vals) for key, vals in grouped.items() if finite(vals)}


def with_observable_labels(matrix: dict[tuple[str, str], float]) -> dict[tuple[str, str], float]:
    return {
        (RUN_TAG_LABELS.get(row, row), col): value
        for (row, col), value in matrix.items()
    }


def plot_p8(rows: Sequence[dict[str, object]], figure_dirs: Sequence[Path], run_tags: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    row_labels = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]
    cols = method_k_order()
    for feature_family in ("efun", "deconv_efun"):
        for condition in ("abs", "csplit"):
            subset = [
                row
                for row in rows
                if row.get("feature_family") == feature_family
                and row.get("density_condition") == condition
            ]
            matrix = with_observable_labels(average_matrix(subset, "run_tag", "density_method_k", run_tags, cols))
            if not matrix:
                continue
            for figure_dir in figure_dirs:
                path = figure_dir / f"p8_roi_profile_consistency__{condition}__{feature_family}.png"
                title = f"P8 ROI profile consistency | {condition} | {feature_family} | peak-weighted top5"
                plot_heatmap(matrix, row_labels, cols, title, path, value_format=".2f")
                paths.append(path)
    return paths


def plot_p10(rows: Sequence[dict[str, object]], figure_dirs: Sequence[Path], run_tags: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    row_labels = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]
    cols = method_k_order()
    for condition in ("abs", "csplit"):
        subset = [row for row in rows if row.get("density_condition") == condition]
        matrix = with_observable_labels(average_matrix(subset, "run_tag", "p9_method_k", run_tags, cols))
        if matrix:
            for figure_dir in figure_dirs:
                path = figure_dir / f"p10_roi_profile_consistency_by_p9_method_k__{condition}.png"
                title = f"P10 ROI profile consistency by P9 method-k | density {condition} | peak-weighted top5"
                plot_heatmap(matrix, row_labels, cols, title, path, value_format=".2f")
                paths.append(path)

        matrix = with_observable_labels(average_matrix(subset, "run_tag", "density_method_k", run_tags, cols))
        if matrix:
            for figure_dir in figure_dirs:
                path = figure_dir / f"p10_roi_profile_consistency_by_lfp_density_method_k__{condition}.png"
                title = f"P10 ROI profile consistency by LFP density method-k | {condition} | peak-weighted top5"
                plot_heatmap(matrix, row_labels, cols, title, path, value_format=".2f")
                paths.append(path)

        for run_tag in run_tags:
            rt_subset = [row for row in subset if row.get("run_tag") == run_tag]
            matrix = average_matrix(rt_subset, "p9_method_k", "density_method_k", cols, cols)
            if not matrix:
                continue
            for figure_dir in figure_dirs:
                path = figure_dir / f"p10_roi_profile_consistency_p9_by_lfp_method_k__{condition}__{slug(run_tag)}.png"
                title = f"P10 ROI profile consistency | P9 x LFP density | {condition} | {RUN_TAG_LABELS.get(run_tag, run_tag)}"
                plot_heatmap(matrix, cols, cols, title, path, value_format=".2f")
                paths.append(path)
    return paths


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)

    p8_groups, p8_meta = read_profile_groups(args.p8_profile_csv, "P8", args.datasets, args.run_tags, args.value_column)
    p10_groups, p10_meta = read_profile_groups(args.p10_profile_csv, "P10", args.datasets, args.run_tags, args.value_column)
    p8_rows = summarize_group_consistency(p8_groups, p8_meta)
    p10_rows = summarize_group_consistency(p10_groups, p10_meta)

    p8_csv = args.results_dir / "p8_roi_profile_consistency_summary.csv"
    p10_csv = args.results_dir / "p10_roi_profile_consistency_summary.csv"
    write_csv(p8_csv, p8_rows)
    write_csv(p10_csv, p10_rows)
    write_csv(
        args.results_dir / "p8_roi_profile_consistency_ranked.csv",
        sorted(p8_rows, key=lambda r: as_float(str(r.get("mean_pairwise_roi_corr", ""))), reverse=True),
    )
    write_csv(
        args.results_dir / "p10_roi_profile_consistency_ranked.csv",
        sorted(p10_rows, key=lambda r: as_float(str(r.get("mean_pairwise_roi_corr", ""))), reverse=True),
    )

    figure_paths: list[Path] = []
    if not args.skip_figures:
        p8_dirs = [args.p8_figure_dir, args.results_figure_dir / "p8_roi_profile_consistency"]
        p10_dirs = [args.p10_figure_dir, args.results_figure_dir / "p10_roi_profile_consistency"]
        figure_paths.extend(plot_p8(p8_rows, p8_dirs, args.run_tags))
        figure_paths.extend(plot_p10(p10_rows, p10_dirs, args.run_tags))

    print(f"P8 profile groups: {len(p8_groups)}")
    print(f"P8 consistency rows: {len(p8_rows)}")
    print(f"P10 profile groups: {len(p10_groups)}")
    print(f"P10 consistency rows: {len(p10_rows)}")
    print(f"Value column: {args.value_column}")
    print(f"Results dir: {args.results_dir}")
    print(f"P8 figure dir: {args.p8_figure_dir}")
    print(f"P10 figure dir: {args.p10_figure_dir}")
    print(f"Figures: {len(figure_paths)}")


if __name__ == "__main__":
    main()
