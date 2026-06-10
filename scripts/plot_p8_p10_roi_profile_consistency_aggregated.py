"""Aggregate P8/P10 ROI profiles into cross-dataset consistency views.

This script consumes the long ROI-vector CSVs exported by
``export_p8_p10_roi_profile_consistency_sources.m``.  The older plotting
script matches exact density/component names across datasets; that is too
strict for P10 because component identities are dataset-local.  Here we
aggregate exact hit groups into broader analysis keys such as:

    pipeline x observable x efun/deconv x density method-k

and compare the resulting ROI profile shape across datasets.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np

from summarize_pipeline8_cross_session_consistency import (
    COMPONENT_COUNTS,
    METHOD_ORDER,
    RUN_TAG_LABELS,
    plot_heatmap,
    write_csv,
)


DEFAULT_RESULTS_DIR = Path("results") / "pipeline_roi_profile_consistency_current"
DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01", "k13m17", "k13m23")
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")
STABLE_METHOD_K = (
    "mds_k04",
    "mds_k05",
    "mds_k06",
    "mds_k07",
    "mds_k08",
    "nmf_k04",
    "umap_k04",
    "umap_k05",
    "umap_k06",
    "umap_k08",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--p8-profile-csv", type=Path, default=DEFAULT_RESULTS_DIR / "p8_roi_profiles_long.csv")
    parser.add_argument("--p10-profile-csv", type=Path, default=DEFAULT_RESULTS_DIR / "p10_roi_profiles_long.csv")
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_RESULTS_DIR / "figures" / "aggregated_roi_profile_consistency")
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--conditions", nargs="+", default=["csplit"])
    parser.add_argument("--value-column", default="roi_value_weighted")
    parser.add_argument("--weight-column", default="mean_peak_abs_corr")
    parser.add_argument("--method-k-scope", choices=["all", "p5-stable"], default="p5-stable")
    parser.add_argument("--top-example-groups", type=int, default=24)
    return parser.parse_args()


def method_k_order() -> list[str]:
    return [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]


def as_float(value: object) -> float:
    try:
        out = float(value) if value not in (None, "") else math.nan
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return math.nan
    x = np.asarray([p[0] for p in pairs], dtype=float)
    y = np.asarray([p[1] for p in pairs], dtype=float)
    x = x - np.mean(x)
    y = y - np.mean(y)
    denom = float(np.sqrt(np.sum(x * x) * np.sum(y * y)))
    if denom <= 0:
        return math.nan
    return float(np.sum(x * y) / denom)


def safe_slug(text: str) -> str:
    out = []
    for ch in str(text):
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_") or "blank"


def row_method_k_ok(row: dict[str, str], scope: str) -> bool:
    if scope == "all":
        return True
    density_mk = row.get("density_method_k", "")
    p9_mk = row.get("p9_method_k", "")
    if density_mk and density_mk not in STABLE_METHOD_K:
        return False
    if p9_mk and p9_mk not in STABLE_METHOD_K:
        return False
    return True


def aggregate_profiles(
    csv_path: Path,
    pipeline: str,
    datasets: set[str],
    run_tags: set[str],
    conditions: set[str],
    value_column: str,
    weight_column: str,
    method_k_scope: str,
) -> tuple[
    dict[tuple[str, ...], dict[str, dict[str, float]]],
    dict[tuple[str, ...], dict[str, str]],
]:
    # broad_key -> dataset -> roi -> (weighted sum, weight sum)
    sums: dict[tuple[str, ...], dict[str, dict[str, list[float]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: [0.0, 0.0]))
    )
    meta: dict[tuple[str, ...], dict[str, str]] = {}

    if not csv_path.is_file():
        return {}, {}

    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("pipeline", "").upper() != pipeline.upper():
                continue
            dataset = row.get("dataset", "")
            run_tag = row.get("run_tag", "")
            condition = row.get("density_condition", "")
            if dataset.lower() not in datasets:
                continue
            if run_tag.lower() not in run_tags:
                continue
            if condition.lower() not in conditions:
                continue
            if row.get("density_source_kind", "") != "dimred":
                continue
            if not row_method_k_ok(row, method_k_scope):
                continue
            roi = row.get("roi_label", "")
            value = as_float(row.get(value_column))
            if not roi or not math.isfinite(value):
                continue
            weight = as_float(row.get(weight_column))
            if not math.isfinite(weight) or weight <= 0:
                weight = 1.0

            if pipeline.upper() == "P8":
                broad_key = (
                    "P8",
                    run_tag,
                    row.get("observable", ""),
                    row.get("feature_family", ""),
                    condition,
                    row.get("density_method_k", ""),
                )
                info = {
                    "pipeline": "P8",
                    "run_tag": run_tag,
                    "observable": row.get("observable", ""),
                    "p9_feature": "",
                    "p9_method_k": "",
                    "feature_family": row.get("feature_family", ""),
                    "density_condition": condition,
                    "density_method_k": row.get("density_method_k", ""),
                    "method_k_scope": method_k_scope,
                }
            else:
                broad_key = (
                    "P10",
                    run_tag,
                    row.get("observable", ""),
                    row.get("p9_feature", ""),
                    row.get("p9_method_k", ""),
                    row.get("feature_family", ""),
                    condition,
                    row.get("density_method_k", ""),
                )
                info = {
                    "pipeline": "P10",
                    "run_tag": run_tag,
                    "observable": row.get("observable", ""),
                    "p9_feature": row.get("p9_feature", ""),
                    "p9_method_k": row.get("p9_method_k", ""),
                    "feature_family": row.get("feature_family", ""),
                    "density_condition": condition,
                    "density_method_k": row.get("density_method_k", ""),
                    "method_k_scope": method_k_scope,
                }

            slot = sums[broad_key][dataset][roi]
            slot[0] += value * weight
            slot[1] += weight
            meta.setdefault(broad_key, info)

    vectors: dict[tuple[str, ...], dict[str, dict[str, float]]] = defaultdict(dict)
    for broad_key, dataset_map in sums.items():
        for dataset, roi_map in dataset_map.items():
            vectors[broad_key][dataset] = {
                roi: sw[0] / sw[1]
                for roi, sw in roi_map.items()
                if sw[1] > 0 and math.isfinite(sw[0])
            }
    return dict(vectors), meta


def summarize_consistency(
    vectors: dict[tuple[str, ...], dict[str, dict[str, float]]],
    meta: dict[tuple[str, ...], dict[str, str]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key, dataset_vectors in vectors.items():
        dataset_names = sorted(dataset_vectors)
        if len(dataset_names) < 2:
            continue
        common_rois = sorted(set.intersection(*(set(dataset_vectors[d]) for d in dataset_names)))
        if len(common_rois) < 2:
            continue
        corrs: list[float] = []
        pair_text: list[str] = []
        for i, d1 in enumerate(dataset_names):
            v1 = [dataset_vectors[d1][roi] for roi in common_rois]
            for d2 in dataset_names[i + 1 :]:
                v2 = [dataset_vectors[d2][roi] for roi in common_rois]
                r = pearson(v1, v2)
                if math.isfinite(r):
                    corrs.append(r)
                    pair_text.append(f"{d1}-{d2}:{r:.4f}")
        if not corrs:
            continue
        info = dict(meta[key])
        info.update(
            {
                "n_datasets": len(dataset_names),
                "datasets": ";".join(dataset_names),
                "n_rois_common": len(common_rois),
                "mean_pairwise_roi_corr": f"{mean(corrs):.10g}",
                "min_pairwise_roi_corr": f"{min(corrs):.10g}",
                "max_pairwise_roi_corr": f"{max(corrs):.10g}",
                "pairwise_roi_corrs": ";".join(pair_text),
            }
        )
        rows.append(info)
    rows.sort(
        key=lambda r: (
            r.get("pipeline", ""),
            r.get("observable", ""),
            r.get("feature_family", ""),
            r.get("p9_method_k", ""),
            r.get("density_method_k", ""),
        )
    )
    return rows


def average_matrix(
    rows: Sequence[dict[str, object]],
    row_fn,
    col_fn,
    row_labels: Sequence[str],
    col_labels: Sequence[str],
) -> dict[tuple[str, str], float]:
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        r = row_fn(row)
        c = col_fn(row)
        value = as_float(row.get("mean_pairwise_roi_corr"))
        if r in row_labels and c in col_labels and math.isfinite(value):
            grouped[(r, c)].append(value)
    return {key: mean(vals) for key, vals in grouped.items() if vals}


def plot_summary_heatmaps(rows: list[dict[str, object]], figure_dir: Path) -> list[Path]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    obs_labels = list(dict.fromkeys(RUN_TAG_LABELS.values()))
    feature_cols = ["efun", "deconv_efun"]
    method_cols = method_k_order()
    obs_feature_rows = [f"{obs} | {feature}" for obs in obs_labels for feature in feature_cols]

    for pipeline in ("P8", "P10"):
        subset = [r for r in rows if r.get("pipeline") == pipeline]
        matrix = average_matrix(
            subset,
            lambda r: RUN_TAG_LABELS.get(str(r.get("run_tag", "")), str(r.get("observable", ""))),
            lambda r: str(r.get("feature_family", "")),
            obs_labels,
            feature_cols,
        )
        if matrix:
            path = figure_dir / f"{pipeline.lower()}_roi_consistency_by_observable_feature.png"
            plot_heatmap(matrix, obs_labels, feature_cols, f"{pipeline} aggregated ROI consistency | observable x efun type", path)
            paths.append(path)

        matrix = average_matrix(
            subset,
            lambda r: f"{RUN_TAG_LABELS.get(str(r.get('run_tag', '')), str(r.get('observable', '')))} | {r.get('feature_family', '')}",
            lambda r: str(r.get("density_method_k", "")),
            obs_feature_rows,
            method_cols,
        )
        if matrix:
            path = figure_dir / f"{pipeline.lower()}_roi_consistency_by_lfp_density_method_k.png"
            plot_heatmap(matrix, obs_feature_rows, method_cols, f"{pipeline} aggregated ROI consistency | LFP density method-k", path)
            paths.append(path)

        if pipeline == "P10":
            matrix = average_matrix(
                subset,
                lambda r: f"{RUN_TAG_LABELS.get(str(r.get('run_tag', '')), str(r.get('observable', '')))} | {r.get('feature_family', '')}",
                lambda r: str(r.get("p9_method_k", "")),
                obs_feature_rows,
                method_cols,
            )
            if matrix:
                path = figure_dir / "p10_roi_consistency_by_bold_p9_method_k.png"
                plot_heatmap(matrix, obs_feature_rows, method_cols, "P10 aggregated ROI consistency | BOLD P9 method-k", path)
                paths.append(path)
    return paths


def row_normalize(values: list[float]) -> list[float]:
    arr = np.asarray(values, dtype=float)
    mask = np.isfinite(arr)
    if mask.sum() < 2:
        return [math.nan for _ in values]
    mu = float(np.nanmean(arr[mask]))
    sd = float(np.nanstd(arr[mask]))
    if sd <= 0:
        sd = 1.0
    out = (arr - mu) / sd
    return [float(x) if math.isfinite(float(x)) else math.nan for x in out]


def plot_example_profiles(
    rows: list[dict[str, object]],
    vectors: dict[tuple[str, ...], dict[str, dict[str, float]]],
    meta: dict[tuple[str, ...], dict[str, str]],
    figure_dir: Path,
    limit: int,
) -> list[Path]:
    example_dir = figure_dir / "roi_profile_examples"
    example_dir.mkdir(parents=True, exist_ok=True)
    by_signature = {
        tuple(sorted(info.items())): key
        for key, info in meta.items()
    }
    ranked = sorted(
        rows,
        key=lambda r: (int(r.get("n_datasets", 0)), as_float(r.get("mean_pairwise_roi_corr"))),
        reverse=True,
    )
    paths: list[Path] = []
    for row in ranked[:limit]:
        signature = tuple(sorted((k, str(row.get(k, ""))) for k in meta[next(iter(meta))].keys())) if meta else ()
        key = None
        for candidate_key, info in meta.items():
            if all(str(row.get(k, "")) == str(v) for k, v in info.items()):
                key = candidate_key
                break
        if key is None:
            continue
        dataset_vectors = vectors.get(key, {})
        datasets = sorted(dataset_vectors)
        if len(datasets) < 2:
            continue
        rois = sorted(set.intersection(*(set(dataset_vectors[d]) for d in datasets)))
        if len(rois) < 2:
            continue
        matrix = np.asarray([row_normalize([dataset_vectors[d].get(roi, math.nan) for roi in rois]) for d in datasets])
        fig_w = max(9, 0.38 * len(rois))
        fig_h = max(3.2, 0.45 * len(datasets) + 1.4)
        fig, ax = plt.subplots(figsize=(fig_w, fig_h))
        im = ax.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=-2.5, vmax=2.5)
        ax.set_yticks(range(len(datasets)))
        ax.set_yticklabels(datasets)
        ax.set_xticks(range(len(rois)))
        ax.set_xticklabels(rois, rotation=70, ha="right", fontsize=7)
        obs = RUN_TAG_LABELS.get(str(row.get("run_tag", "")), str(row.get("observable", "")))
        title_parts = [
            str(row.get("pipeline", "")),
            obs,
            str(row.get("feature_family", "")),
        ]
        if row.get("p9_method_k"):
            title_parts.append(f"P9 {row.get('p9_method_k')}")
        title_parts.append(f"LFP {row.get('density_method_k')}")
        title_parts.append(f"r={as_float(row.get('mean_pairwise_roi_corr')):.2f}")
        ax.set_title(" | ".join(title_parts), fontsize=11)
        ax.set_xlabel("ROI")
        fig.colorbar(im, ax=ax, shrink=0.82, label="within-dataset ROI profile z")
        fig.tight_layout()
        filename = "__".join(
            safe_slug(x)
            for x in [
                row.get("pipeline", ""),
                obs,
                row.get("feature_family", ""),
                row.get("p9_method_k", ""),
                row.get("density_method_k", ""),
            ]
            if x
        )
        path = example_dir / f"{len(paths)+1:02d}_{filename}.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(path)
    return paths


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    datasets = {d.lower() for d in args.datasets}
    run_tags = {r.lower() for r in args.run_tags}
    conditions = {c.lower() for c in args.conditions}

    p8_vectors, p8_meta = aggregate_profiles(
        args.p8_profile_csv,
        "P8",
        datasets,
        run_tags,
        conditions,
        args.value_column,
        args.weight_column,
        args.method_k_scope,
    )
    p10_vectors, p10_meta = aggregate_profiles(
        args.p10_profile_csv,
        "P10",
        datasets,
        run_tags,
        conditions,
        args.value_column,
        args.weight_column,
        args.method_k_scope,
    )
    p8_rows = summarize_consistency(p8_vectors, p8_meta)
    p10_rows = summarize_consistency(p10_vectors, p10_meta)
    all_rows = p8_rows + p10_rows

    suffix = f"_{args.method_k_scope.replace('-', '_')}"
    write_csv(args.results_dir / f"p8_roi_profile_consistency_aggregated{suffix}.csv", p8_rows)
    write_csv(args.results_dir / f"p10_roi_profile_consistency_aggregated{suffix}.csv", p10_rows)
    write_csv(
        args.results_dir / f"p8_p10_roi_profile_consistency_aggregated_ranked{suffix}.csv",
        sorted(all_rows, key=lambda r: as_float(r.get("mean_pairwise_roi_corr")), reverse=True),
    )

    figure_dir = args.figure_dir / args.method_k_scope.replace("-", "_")
    heatmap_paths = plot_summary_heatmaps(all_rows, figure_dir)
    example_paths = plot_example_profiles(
        all_rows,
        {**p8_vectors, **p10_vectors},
        {**p8_meta, **p10_meta},
        figure_dir,
        args.top_example_groups,
    )
    print(f"P8 aggregated groups: {len(p8_vectors)}")
    print(f"P8 consistency rows: {len(p8_rows)}")
    print(f"P10 aggregated groups: {len(p10_vectors)}")
    print(f"P10 consistency rows: {len(p10_rows)}")
    print(f"Method-k scope: {args.method_k_scope}")
    print(f"Results dir: {args.results_dir}")
    print(f"Figure dir: {figure_dir}")
    print(f"Heatmaps: {len(heatmap_paths)}")
    print(f"Example profile figures: {len(example_paths)}")


if __name__ == "__main__":
    main()
