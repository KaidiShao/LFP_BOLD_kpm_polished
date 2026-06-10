"""Summarize roi_mean ROI profiles for P8 RG-no-theta dimred hits.

This is a deliberately narrow readout:

    BLP: standardized complex-split dimred density
    BOLD: roi_mean
    xcorr: P8 dimred hits with strict_label == ripple_gamma_no_theta

For each dataset x BOLD feature family x BLP method-k, the script selects
topN hits by peak_abs_corr, maps the selected BOLD mode indices to P7 roi_mean
ROI vectors, and reports which ROI labels dominate.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np


METHOD_ORDER = ("svd", "nmf", "mds", "umap")
K_ORDER = (3, 4, 5, 6, 7, 8)
DEFAULT_TOP_NS = (1, 5, 10)
DEFAULT_FEATURES = ("deconv_efun", "efun")
CORTICAL_ROIS = {
    "ParPrec",
    "ParIntra",
    "ParLat",
    "TmpVis",
    "TmpSTS",
    "TmpPol",
    "TmpAu",
    "V5",
    "V4",
    "V2V3",
    "pV1",
    "pfV1",
    "fV1",
    "V1",
    "dlPFC",
    "medPFC",
    "orbPFC",
    "RetroSp",
    "PCC",
    "dACC",
    "ACC",
    "aIns",
    "Ins",
    "S2",
    "S1",
    "Motor",
    "Premotor",
    "PirFo",
    "Ent",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hits-csv",
        type=Path,
        nargs="+",
        default=[
            Path("results")
            / "pipeline8_10_strict_band_coupling_k13m18_m21_20260604"
            / "p8_p10_strict_band_hits_long.csv"
        ],
    )
    parser.add_argument(
        "--roi-profile-csv",
        type=Path,
        nargs="+",
        default=[
            Path("results")
            / "pipeline_roi_profile_consistency_k13m18_m21_20260604"
            / "p7_intrinsic_bold_efun_roi_profiles_long.csv"
        ],
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results") / "k13_roi_mean_rg_no_theta_roi_summary_20260604",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("E:/DataPons_processed/summary_figures/pipeline11_current_analysis_summary/k13_roi_mean_rg_no_theta_roi_summary_20260604"),
    )
    parser.add_argument("--datasets", nargs="+", default=["k13m18", "k13m19", "k13m20", "k13m21"])
    parser.add_argument("--pipeline", default="P8")
    parser.add_argument("--bold-observable", default="roi_mean")
    parser.add_argument("--features", nargs="+", default=list(DEFAULT_FEATURES))
    parser.add_argument("--strict-label", default="ripple_gamma_no_theta")
    parser.add_argument("--density-class", default="dimred_efun_density")
    parser.add_argument("--top-ns", type=int, nargs="+", default=list(DEFAULT_TOP_NS))
    parser.add_argument("--main-top-n", type=int, default=5)
    parser.add_argument("--top-rois", type=int, default=20)
    parser.add_argument("--min-datasets-for-top-roi", type=int, default=1)
    parser.add_argument(
        "--roi-order",
        choices=["p7_plot", "p7_native", "median"],
        default="p7_plot",
        help=(
            "ROI order for all-ROI plots and ROI-profile correlations. "
            "p7_plot matches the default flipped P7 ROI bar-summary display; "
            "p7_native preserves the source roiTs/P7 CSV order; median is only for debugging."
        ),
    )
    parser.add_argument("--dataset-panel-cols", type=int, default=12)
    parser.add_argument(
        "--dataset-panel-scale",
        choices=["per_dataset", "shared"],
        default="per_dataset",
        help=(
            "Color scaling for all-ROI-by-dataset panels. per_dataset is the "
            "default for pattern inspection; shared keeps absolute values "
            "visually comparable across datasets."
        ),
    )
    return parser.parse_args()


def as_float(value: object) -> float:
    try:
        out = float(value) if value not in (None, "") else math.nan
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def as_int(value: object) -> int | None:
    value_f = as_float(value)
    if not math.isfinite(value_f):
        return None
    return int(value_f)


def safe_name(text: object) -> str:
    out = []
    for ch in str(text):
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_") or "blank"


def write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys: list[str] = []
        for row in rows:
            for key in row:
                if key not in keys:
                    keys.append(key)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def method_k_order() -> list[str]:
    return [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in K_ORDER]


def is_standardized_csplit(row: dict[str, str]) -> bool:
    pieces = [
        row.get("density_condition", ""),
        row.get("density_name", ""),
        row.get("density_suffix", ""),
        row.get("source_csv", ""),
    ]
    text = " ".join(str(piece).lower() for piece in pieces)
    return "csplit" in text and ("standardize" in text or "standardized" in text)


def read_p7_roi_profiles(
    csv_paths: Sequence[Path],
    *,
    datasets: set[str],
    observable: str,
) -> dict[tuple[str, int], dict[str, float]]:
    out: dict[tuple[str, int], dict[str, float]] = defaultdict(dict)
    for csv_path in csv_paths:
        with csv_path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                dataset = row.get("dataset", "").lower()
                if dataset not in datasets:
                    continue
                if row.get("observable", "") != observable:
                    continue
                pos = as_int(row.get("sorted_position"))
                roi = row.get("roi_label", "")
                value = as_float(row.get("roi_value"))
                if pos is None or not roi or not math.isfinite(value):
                    continue
                out[(dataset, pos)][roi] = value
    return dict(out)


def read_p7_roi_order(
    csv_paths: Sequence[Path],
    *,
    datasets: Sequence[str],
    observable: str,
    order_mode: str,
) -> list[str]:
    """Read the native ROI order from the P7 profile CSV, preserving anatomy order.

    The CSV has one ROI row per BOLD mode. sorted_position is the BOLD mode rank,
    not an ROI order field, so the ROI order is the row order within the first
    available mode for each dataset. Multiple datasets are merged with stable
    append semantics because not every dataset exposes exactly the same ROI set.
    """
    if order_mode == "median":
        return []
    dataset_order = [dataset.lower() for dataset in datasets]
    ranks_by_roi: dict[str, list[int]] = defaultdict(list)
    for csv_path in csv_paths:
        by_dataset_mode: dict[tuple[str, str], list[str]] = defaultdict(list)
        with csv_path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                dataset = row.get("dataset", "").lower()
                if dataset not in dataset_order:
                    continue
                if row.get("observable", "") != observable:
                    continue
                mode = row.get("sorted_position", "")
                roi = row.get("roi_label", "")
                if not mode or not roi:
                    continue
                by_dataset_mode[(dataset, mode)].append(roi)
        for dataset in dataset_order:
            mode_keys = [
                key for key in by_dataset_mode
                if key[0] == dataset and len(by_dataset_mode[key]) > 1
            ]
            mode_keys = sorted(mode_keys, key=lambda key: as_int(key[1]) or 10**9)
            if not mode_keys:
                continue
            for rank, roi in enumerate(by_dataset_mode[mode_keys[0]], start=1):
                ranks_by_roi[roi].append(rank)
    out = sorted(
        ranks_by_roi,
        key=lambda roi: (float(np.median(ranks_by_roi[roi])), roi),
    )
    if order_mode == "p7_plot":
        out = list(reversed(out))
    return out


def read_filtered_hits(args: argparse.Namespace) -> list[dict[str, object]]:
    datasets = {dataset.lower() for dataset in args.datasets}
    features = set(args.features)
    rows: list[dict[str, object]] = []
    for csv_path in args.hits_csv:
        with csv_path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                dataset = row.get("dataset", "").lower()
                if dataset not in datasets:
                    continue
                if row.get("pipeline", "") != args.pipeline:
                    continue
                if row.get("bold_observable", "") != args.bold_observable:
                    continue
                if row.get("bold_feature_family", "") not in features:
                    continue
                if row.get("density_class", "") != args.density_class:
                    continue
                if row.get("strict_label", "") != args.strict_label:
                    continue
                if not is_standardized_csplit(row):
                    continue
                if row.get("finite_peak", "True").lower() not in {"true", "1", "yes"}:
                    continue
                method_k = row.get("density_method_k", "")
                mode = as_int(row.get("bold_mode_index"))
                peak = as_float(row.get("peak_abs_corr"))
                if not method_k or mode is None or not math.isfinite(peak):
                    continue
                rows.append(
                    {
                        "pipeline": row.get("pipeline", ""),
                        "dataset": dataset,
                        "bold_observable": row.get("bold_observable", ""),
                        "bold_feature_family": row.get("bold_feature_family", ""),
                        "bold_feature": row.get("bold_feature", ""),
                        "density_method_k": method_k,
                        "density_method": row.get("density_method", ""),
                        "density_k": row.get("density_k", ""),
                        "density_index": row.get("density_index", ""),
                        "density_label": row.get("density_label", ""),
                        "strict_label": row.get("strict_label", ""),
                        "theta_effect_z": row.get("theta_effect_z", ""),
                        "gamma_effect_z": row.get("gamma_effect_z", ""),
                        "ripple_effect_z": row.get("ripple_effect_z", ""),
                        "ripple_gamma_effect_z": row.get("ripple_gamma_effect_z", ""),
                        "component_timescale_sec": row.get("component_timescale_sec", ""),
                        "bold_mode_index": mode,
                        "peak_abs_corr": peak,
                        "peak_corr": row.get("peak_corr", ""),
                        "peak_lag_sec": row.get("peak_lag_sec", ""),
                        "source_csv": row.get("source_csv", ""),
                    }
                )
    return rows


def top_hits_by_context(rows: Sequence[dict[str, object]], top_n: int) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = (
            str(row["dataset"]),
            str(row["pipeline"]),
            str(row["bold_feature_family"]),
            str(row["density_method_k"]),
        )
        groups[key].append(row)

    selected: list[dict[str, object]] = []
    for key_rows in groups.values():
        key_rows = sorted(key_rows, key=lambda row: float(row["peak_abs_corr"]), reverse=True)
        for rank, row in enumerate(key_rows[:top_n], start=1):
            out = dict(row)
            out["top_n"] = top_n
            out["rank_within_context"] = rank
            selected.append(out)
    return selected


def build_dataset_roi_profiles(
    selected_hits: Sequence[dict[str, object]],
    roi_profiles: dict[tuple[str, int], dict[str, float]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    accum: dict[tuple[str, str, str, str, int], dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(lambda: [0.0, 0.0])
    )
    used_hits: list[dict[str, object]] = []
    for row in selected_hits:
        dataset = str(row["dataset"])
        mode = int(row["bold_mode_index"])
        peak = float(row["peak_abs_corr"])
        roi_vec = roi_profiles.get((dataset, mode), {})
        if not roi_vec:
            continue
        key = (
            dataset,
            str(row["pipeline"]),
            str(row["bold_feature_family"]),
            str(row["density_method_k"]),
            int(row["top_n"]),
        )
        for roi, value in roi_vec.items():
            accum[key][roi][0] += value * peak
            accum[key][roi][1] += peak
        used = dict(row)
        used["n_rois_for_bold_mode"] = len(roi_vec)
        used_hits.append(used)

    profile_rows: list[dict[str, object]] = []
    for key, roi_accum in accum.items():
        dataset, pipeline, feature, method_k, top_n = key
        for roi, (weighted_sum, weight_sum) in roi_accum.items():
            if weight_sum <= 0:
                continue
            profile_rows.append(
                {
                    "dataset": dataset,
                    "pipeline": pipeline,
                    "bold_observable": "roi_mean",
                    "bold_feature_family": feature,
                    "density_method_k": method_k,
                    "top_n": top_n,
                    "roi_label": roi,
                    "roi_value_weighted": weighted_sum / weight_sum,
                    "weight_sum_peak_abs_corr": weight_sum,
                }
            )
    return profile_rows, used_hits


def aggregate_roi_profiles(profile_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, int, str], list[dict[str, object]]] = defaultdict(list)
    for row in profile_rows:
        key = (
            str(row["pipeline"]),
            str(row["bold_feature_family"]),
            str(row["density_method_k"]),
            int(row["top_n"]),
            str(row["roi_label"]),
        )
        grouped[key].append(row)

    out_rows: list[dict[str, object]] = []
    for key, rows in grouped.items():
        pipeline, feature, method_k, top_n, roi = key
        vals = np.asarray([float(row["roi_value_weighted"]) for row in rows], dtype=float)
        datasets = sorted({str(row["dataset"]) for row in rows})
        out_rows.append(
            {
                "pipeline": pipeline,
                "bold_observable": "roi_mean",
                "bold_feature_family": feature,
                "density_method_k": method_k,
                "top_n": top_n,
                "roi_label": roi,
                "n_datasets": len(datasets),
                "datasets": ";".join(datasets),
                "mean_roi_value": float(np.mean(vals)),
                "median_roi_value": float(np.median(vals)),
                "std_roi_value": float(np.std(vals)) if len(vals) > 1 else 0.0,
                "min_roi_value": float(np.min(vals)),
                "max_roi_value": float(np.max(vals)),
            }
        )
    return out_rows


def top_roi_rows(agg_rows: Sequence[dict[str, object]], top_rois: int, min_datasets: int) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, int], list[dict[str, object]]] = defaultdict(list)
    for row in agg_rows:
        if int(row["n_datasets"]) < min_datasets:
            continue
        key = (
            str(row["pipeline"]),
            str(row["bold_feature_family"]),
            str(row["density_method_k"]),
            int(row["top_n"]),
        )
        grouped[key].append(row)

    out: list[dict[str, object]] = []
    for key, rows in grouped.items():
        rows = sorted(rows, key=lambda row: float(row["mean_roi_value"]), reverse=True)
        for rank, row in enumerate(rows[:top_rois], start=1):
            item = dict(row)
            item["roi_rank_within_method"] = rank
            out.append(item)
    return out


def overall_roi_rows(agg_rows: Sequence[dict[str, object]], *, top_n: int, feature: str, top_rois: int) -> list[dict[str, object]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in agg_rows:
        if int(row["top_n"]) != top_n:
            continue
        if row["bold_feature_family"] != feature:
            continue
        grouped[str(row["roi_label"])].append(row)
    out: list[dict[str, object]] = []
    for roi, rows in grouped.items():
        vals = np.asarray([float(row["mean_roi_value"]) for row in rows], dtype=float)
        methods = sorted({str(row["density_method_k"]) for row in rows})
        datasets = sorted(set(";".join(str(row["datasets"]) for row in rows).split(";")) - {""})
        out.append(
            {
                "bold_feature_family": feature,
                "top_n": top_n,
                "roi_label": roi,
                "n_methods": len(methods),
                "methods": ";".join(methods),
                "n_datasets_any_method": len(datasets),
                "mean_across_methods": float(np.mean(vals)),
                "median_across_methods": float(np.median(vals)),
                "max_across_methods": float(np.max(vals)),
            }
        )
    out = sorted(out, key=lambda row: float(row["mean_across_methods"]), reverse=True)
    for i, row in enumerate(out[:top_rois], start=1):
        row["overall_roi_rank"] = i
    return out[:top_rois]


def plot_heatmap(
    agg_rows: Sequence[dict[str, object]],
    figure_dir: Path,
    *,
    top_n: int,
    feature: str,
    top_rois: int,
) -> Path | None:
    rows = [
        row
        for row in agg_rows
        if int(row["top_n"]) == top_n and row["bold_feature_family"] == feature
    ]
    if not rows:
        return None
    overall = overall_roi_rows(rows, top_n=top_n, feature=feature, top_rois=top_rois)
    roi_order = [row["roi_label"] for row in overall]
    method_order = [mk for mk in method_k_order() if any(row["density_method_k"] == mk for row in rows)]
    if not roi_order or not method_order:
        return None

    lookup = {
        (row["roi_label"], row["density_method_k"]): float(row["mean_roi_value"])
        for row in rows
    }
    matrix = np.full((len(roi_order), len(method_order)), np.nan, dtype=float)
    for i, roi in enumerate(roi_order):
        for j, method_k in enumerate(method_order):
            matrix[i, j] = lookup.get((roi, method_k), math.nan)

    fig_h = max(6.0, 0.32 * len(roi_order) + 2.0)
    fig_w = max(18.0, 0.55 * len(method_order) + 5.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(matrix, cmap="viridis", aspect="auto")
    ax.set_title(
        f"P8 roi_mean RG-no-theta ROI profile by BLP method-k | {feature} | top{top_n}",
        fontsize=13,
        weight="bold",
    )
    ax.set_xticks(range(len(method_order)))
    ax.set_xticklabels(method_order, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(roi_order)))
    ax.set_yticklabels(roi_order, fontsize=8)
    ax.set_xlabel("BLP dimred density method-k")
    ax.set_ylabel("ROI label")
    for i in range(len(roi_order)):
        for j in range(len(method_order)):
            val = matrix[i, j]
            if math.isfinite(val):
                ax.text(j, i, f"{val:.2g}", ha="center", va="center", fontsize=6, color="white" if val > np.nanmax(matrix) * 0.55 else "black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cbar.set_label("weighted mean roi_value")
    fig.tight_layout()
    figure_dir.mkdir(parents=True, exist_ok=True)
    out_file = figure_dir / f"01_p8_roi_mean_rg_no_theta_roi_by_method__{safe_name(feature)}__top{top_n}.png"
    fig.savefig(out_file, dpi=180)
    plt.close(fig)
    return out_file


def plot_top_roi_bars(
    overall_rows: Sequence[dict[str, object]],
    figure_dir: Path,
    *,
    top_n: int,
    feature: str,
) -> Path | None:
    rows = [row for row in overall_rows if row["bold_feature_family"] == feature and int(row["top_n"]) == top_n]
    if not rows:
        return None
    labels = [str(row["roi_label"]) for row in rows][::-1]
    vals = [float(row["mean_across_methods"]) for row in rows][::-1]
    fig_h = max(5.0, 0.32 * len(labels) + 1.4)
    fig, ax = plt.subplots(figsize=(14.0, fig_h))
    ax.barh(labels, vals, color="#5b8cc0")
    ax.set_title(f"P8 roi_mean RG-no-theta overall top ROIs | {feature} | top{top_n}", weight="bold")
    ax.set_xlabel("mean weighted roi_value across BLP method-k")
    ax.grid(axis="x", color="#dddddd", linewidth=0.7)
    fig.tight_layout()
    figure_dir.mkdir(parents=True, exist_ok=True)
    out_file = figure_dir / f"02_p8_roi_mean_rg_no_theta_overall_top_rois__{safe_name(feature)}__top{top_n}.png"
    fig.savefig(out_file, dpi=180)
    plt.close(fig)
    return out_file


def finite_corr(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if int(np.sum(mask)) < 3:
        return math.nan
    aa = a[mask]
    bb = b[mask]
    if float(np.std(aa)) == 0.0 or float(np.std(bb)) == 0.0:
        return math.nan
    return float(np.corrcoef(aa, bb)[0, 1])


def roi_order_for_display(
    profile_rows: Sequence[dict[str, object]],
    *,
    top_n: int,
    feature: str,
    canonical_roi_order: Sequence[str] | None = None,
) -> list[str]:
    """Return all ROI labels, preserving the requested anatomical order."""
    present = {
        str(row["roi_label"])
        for row in profile_rows
        if int(row["top_n"]) == top_n and row["bold_feature_family"] == feature
    }
    if canonical_roi_order:
        ordered = [roi for roi in canonical_roi_order if roi in present]
        extras = sorted(present - set(ordered))
        return ordered + extras

    values: dict[str, list[float]] = defaultdict(list)
    for row in profile_rows:
        if int(row["top_n"]) != top_n:
            continue
        if row["bold_feature_family"] != feature:
            continue
        values[str(row["roi_label"])].append(float(row["roi_value_weighted"]))
    return sorted(
        values,
        key=lambda roi: (float(np.nanmedian(values[roi])), roi),
        reverse=True,
    )


def plot_all_roi_dataset_panels(
    profile_rows: Sequence[dict[str, object]],
    figure_dir: Path,
    *,
    top_n: int,
    feature: str,
    datasets: Sequence[str],
    canonical_roi_order: Sequence[str],
    max_cols: int = 4,
    scale_mode: str = "per_dataset",
) -> Path | None:
    rows = [
        row
        for row in profile_rows
        if int(row["top_n"]) == top_n and row["bold_feature_family"] == feature
    ]
    if not rows:
        return None

    roi_order = roi_order_for_display(
        profile_rows,
        top_n=top_n,
        feature=feature,
        canonical_roi_order=canonical_roi_order,
    )
    method_order = [mk for mk in method_k_order() if any(row["density_method_k"] == mk for row in rows)]
    dataset_order = [dataset.lower() for dataset in datasets if any(row["dataset"] == dataset.lower() for row in rows)]
    if not roi_order or not method_order or not dataset_order:
        return None

    lookup = {
        (row["dataset"], row["roi_label"], row["density_method_k"]): float(row["roi_value_weighted"])
        for row in rows
    }
    matrices: list[np.ndarray] = []
    for dataset in dataset_order:
        matrix = np.full((len(roi_order), len(method_order)), np.nan, dtype=float)
        for i, roi in enumerate(roi_order):
            for j, method_k in enumerate(method_order):
                matrix[i, j] = lookup.get((dataset, roi, method_k), math.nan)
        matrices.append(matrix)

    finite_by_matrix = [matrix[np.isfinite(matrix)] for matrix in matrices if np.any(np.isfinite(matrix))]
    if not finite_by_matrix:
        return None
    all_values = np.concatenate(finite_by_matrix)

    def color_limits(values: np.ndarray) -> tuple[float, float]:
        vmax_i = float(np.nanpercentile(values, 98))
        if not math.isfinite(vmax_i) or vmax_i <= 0:
            vmax_i = float(np.nanmax(values))
        vmin_i = float(np.nanpercentile(values, 2))
        if not math.isfinite(vmin_i):
            vmin_i = 0.0
        if vmin_i > 0:
            vmin_i = 0.0
        if vmax_i <= vmin_i:
            vmax_i = vmin_i + 1.0
        return vmin_i, vmax_i

    shared_vmin, shared_vmax = color_limits(all_values)

    n_cols = max(1, min(int(max_cols), len(dataset_order)))
    n_rows = int(math.ceil(len(dataset_order) / n_cols))
    fig_h = max(7.5 * n_rows, 0.18 * len(roi_order) * n_rows + 2.4)
    fig_w = max(16.0, 3.7 * n_cols + 4.0)
    fig, axes_arr = plt.subplots(n_rows, n_cols, figsize=(fig_w, fig_h), sharey=True)
    axes_flat = np.asarray(axes_arr).reshape(-1)

    cortical_idx = [idx for idx, roi in enumerate(roi_order) if roi in CORTICAL_ROIS]
    separator_y = (max(cortical_idx) + 0.5) if cortical_idx and max(cortical_idx) < len(roi_order) - 1 else None
    for ax, dataset, matrix in zip(axes_flat, dataset_order, matrices):
        finite_values = matrix[np.isfinite(matrix)]
        if scale_mode == "shared":
            panel_vmin, panel_vmax = shared_vmin, shared_vmax
        else:
            panel_vmin, panel_vmax = color_limits(finite_values)
        im = ax.imshow(matrix, cmap="viridis", aspect="auto", vmin=panel_vmin, vmax=panel_vmax)
        if scale_mode == "shared":
            title = dataset
        else:
            title = f"{dataset}\nscale 0-{panel_vmax:.2g}"
        ax.set_title(title, fontsize=9 if scale_mode == "per_dataset" else 11, weight="bold")
        ax.set_xticks(range(len(method_order)))
        ax.set_xticklabels(method_order, rotation=90, fontsize=6)
        ax.set_yticks(range(len(roi_order)))
        ax.set_yticklabels(roi_order, fontsize=7)
        ax.set_xlabel("method-k")
        if separator_y is not None:
            ax.axhline(separator_y, color="white", linewidth=1.2)
        for j, method_k in enumerate(method_order):
            if method_k.endswith("_k03") and j > 0:
                ax.axvline(j - 0.5, color="white", linewidth=1.0)
    for ax in axes_flat[len(dataset_order):]:
        ax.axis("off")
    for idx, ax in enumerate(axes_flat[:len(dataset_order)]):
        if idx % n_cols == 0:
            ax.set_ylabel("ROI label, all ROIs")
        else:
            ax.set_ylabel("")
    fig.suptitle(
        f"P8 roi_mean RG-no-theta all ROI profiles by dataset | {feature} | top{top_n}",
        fontsize=14,
        weight="bold",
    )
    fig.text(
        0.01,
        0.01,
        (
            "No across-dataset averaging. ROI order preserves the P7 ROI bar-summary anatomical order. "
            "Colors use per-dataset scales, so color intensity is not comparable across datasets."
            if scale_mode == "per_dataset"
            else "No across-dataset averaging. ROI order preserves the P7 ROI bar-summary anatomical order. Colors use one shared scale across datasets."
        ),
        fontsize=8,
        color="#555555",
    )
    if scale_mode == "shared":
        cbar = fig.colorbar(im, ax=axes_flat[:len(dataset_order)].tolist(), fraction=0.015, pad=0.01)
        cbar.set_label("dataset-level weighted roi_value")
    figure_dir.mkdir(parents=True, exist_ok=True)
    out_file = figure_dir / f"03_p8_roi_mean_rg_no_theta_all_roi_by_dataset__{safe_name(feature)}__top{top_n}.png"
    fig.savefig(out_file, dpi=180)
    plt.close(fig)
    return out_file


def dataset_pair_roi_correlations(
    profile_rows: Sequence[dict[str, object]],
    *,
    top_n: int,
    feature: str,
    datasets: Sequence[str],
    canonical_roi_order: Sequence[str],
) -> list[dict[str, object]]:
    rows = [
        row
        for row in profile_rows
        if int(row["top_n"]) == top_n and row["bold_feature_family"] == feature
    ]
    roi_order = roi_order_for_display(
        profile_rows,
        top_n=top_n,
        feature=feature,
        canonical_roi_order=canonical_roi_order,
    )
    method_order = [mk for mk in method_k_order() if any(row["density_method_k"] == mk for row in rows)]
    dataset_order = [dataset.lower() for dataset in datasets if any(row["dataset"] == dataset.lower() for row in rows)]
    lookup = {
        (row["dataset"], row["density_method_k"], row["roi_label"]): float(row["roi_value_weighted"])
        for row in rows
    }

    out: list[dict[str, object]] = []
    for method_k in method_order:
        vectors: dict[str, np.ndarray] = {}
        for dataset in dataset_order:
            vectors[dataset] = np.asarray(
                [lookup.get((dataset, method_k, roi), math.nan) for roi in roi_order],
                dtype=float,
            )
        for i, dataset_a in enumerate(dataset_order):
            for dataset_b in dataset_order[i + 1 :]:
                corr = finite_corr(vectors[dataset_a], vectors[dataset_b])
                out.append(
                    {
                        "pipeline": "P8",
                        "bold_observable": "roi_mean",
                        "bold_feature_family": feature,
                        "top_n": top_n,
                        "density_method_k": method_k,
                        "dataset_a": dataset_a,
                        "dataset_b": dataset_b,
                        "dataset_pair": f"{dataset_a}__{dataset_b}",
                        "n_rois": int(np.sum(np.isfinite(vectors[dataset_a]) & np.isfinite(vectors[dataset_b]))),
                        "roi_profile_pearson": corr,
                    }
                )
    return out


def plot_dataset_pair_consistency(
    corr_rows: Sequence[dict[str, object]],
    figure_dir: Path,
    *,
    top_n: int,
    feature: str,
) -> Path | None:
    rows = [
        row
        for row in corr_rows
        if int(row["top_n"]) == top_n and row["bold_feature_family"] == feature
    ]
    if not rows:
        return None
    pair_order = sorted({str(row["dataset_pair"]) for row in rows})
    method_order = [mk for mk in method_k_order() if any(row["density_method_k"] == mk for row in rows)]
    if not pair_order or not method_order:
        return None

    lookup = {
        (row["dataset_pair"], row["density_method_k"]): float(row["roi_profile_pearson"])
        for row in rows
    }
    matrix = np.full((len(pair_order), len(method_order)), np.nan, dtype=float)
    for i, pair in enumerate(pair_order):
        for j, method_k in enumerate(method_order):
            matrix[i, j] = lookup.get((pair, method_k), math.nan)

    # Ultrawide-friendly orientation: many dataset pairs belong on x, not y.
    # This keeps the figure readable on a wide monitor instead of producing a
    # tall, skinny matrix.
    matrix = matrix.T
    x_order = pair_order
    y_order = method_order

    fig_w = max(22.0, 0.34 * len(x_order) + 6.0)
    fig_h = max(7.0, 0.33 * len(y_order) + 2.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(matrix, cmap="coolwarm", aspect="auto", vmin=-1.0, vmax=1.0)
    ax.set_title(
        f"P8 roi_mean RG-no-theta all-ROI dataset-pair consistency | {feature} | top{top_n}",
        fontsize=13,
        weight="bold",
    )
    ax.set_xticks(range(len(x_order)))
    ax.set_xticklabels(x_order, rotation=60, ha="right", fontsize=7)
    ax.set_yticks(range(len(y_order)))
    ax.set_yticklabels(y_order, fontsize=8)
    ax.set_xlabel("dataset pair")
    ax.set_ylabel("BLP dimred density method-k")
    for i in range(len(y_order)):
        for j in range(len(x_order)):
            val = matrix[i, j]
            if math.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color="black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Pearson corr across all ROI values")
    fig.tight_layout()
    figure_dir.mkdir(parents=True, exist_ok=True)
    out_file = figure_dir / f"04_p8_roi_mean_rg_no_theta_dataset_pair_consistency__{safe_name(feature)}__top{top_n}.png"
    fig.savefig(out_file, dpi=180)
    plt.close(fig)
    return out_file


def main() -> None:
    args = parse_args()
    datasets = {dataset.lower() for dataset in args.datasets}

    roi_profiles = read_p7_roi_profiles(args.roi_profile_csv, datasets=datasets, observable=args.bold_observable)
    canonical_roi_order = read_p7_roi_order(
        args.roi_profile_csv,
        datasets=args.datasets,
        observable=args.bold_observable,
        order_mode=args.roi_order,
    )
    hits = read_filtered_hits(args)

    selected_all: list[dict[str, object]] = []
    for top_n in args.top_ns:
        selected_all.extend(top_hits_by_context(hits, top_n))

    profile_rows, used_hits = build_dataset_roi_profiles(selected_all, roi_profiles)
    agg_rows = aggregate_roi_profiles(profile_rows)
    top_rows = top_roi_rows(agg_rows, args.top_rois, args.min_datasets_for_top_roi)
    corr_rows: list[dict[str, object]] = []
    for top_n in args.top_ns:
        for feature in args.features:
            corr_rows.extend(
                dataset_pair_roi_correlations(
                    profile_rows,
                    top_n=top_n,
                    feature=feature,
                    datasets=args.datasets,
                    canonical_roi_order=canonical_roi_order,
                )
            )

    overall_rows: list[dict[str, object]] = []
    for top_n in args.top_ns:
        for feature in args.features:
            overall_rows.extend(overall_roi_rows(agg_rows, top_n=top_n, feature=feature, top_rois=args.top_rois))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "p8_roi_mean_rg_no_theta_filtered_hits_all.csv", hits)
    write_csv(args.output_dir / "p8_roi_mean_rg_no_theta_selected_hits.csv", used_hits)
    write_csv(args.output_dir / "p8_roi_mean_rg_no_theta_dataset_roi_profiles_long.csv", profile_rows)
    write_csv(args.output_dir / "p8_roi_mean_rg_no_theta_roi_profiles_by_method.csv", agg_rows)
    write_csv(args.output_dir / "p8_roi_mean_rg_no_theta_top_rois_by_method.csv", top_rows)
    write_csv(args.output_dir / "p8_roi_mean_rg_no_theta_overall_top_rois.csv", overall_rows)
    write_csv(args.output_dir / "p8_roi_mean_rg_no_theta_dataset_pair_roi_correlations.csv", corr_rows)

    figure_paths: list[Path] = []
    for feature in args.features:
        heatmap = plot_heatmap(agg_rows, args.figure_dir, top_n=args.main_top_n, feature=feature, top_rois=args.top_rois)
        if heatmap is not None:
            figure_paths.append(heatmap)
        bars = plot_top_roi_bars(overall_rows, args.figure_dir, top_n=args.main_top_n, feature=feature)
        if bars is not None:
            figure_paths.append(bars)
        all_roi = plot_all_roi_dataset_panels(
            profile_rows,
            args.figure_dir,
            top_n=args.main_top_n,
            feature=feature,
            datasets=args.datasets,
            canonical_roi_order=canonical_roi_order,
            max_cols=args.dataset_panel_cols,
            scale_mode=args.dataset_panel_scale,
        )
        if all_roi is not None:
            figure_paths.append(all_roi)
        consistency = plot_dataset_pair_consistency(corr_rows, args.figure_dir, top_n=args.main_top_n, feature=feature)
        if consistency is not None:
            figure_paths.append(consistency)

    readme = args.output_dir / "README_p8_roi_mean_rg_no_theta_roi_summary.md"
    readme.write_text(
        "\n".join(
            [
                "# P8 roi_mean RG-no-theta ROI summary",
                "",
                "Filters:",
                "",
                "```text",
                f"datasets = {', '.join(args.datasets)}",
                f"pipeline = {args.pipeline}",
                f"BOLD observable = {args.bold_observable}",
                f"BOLD feature families = {', '.join(args.features)}",
                f"density class = {args.density_class}",
                "BLP condition = standardized complex-split",
                f"strict label = {args.strict_label}",
                f"topN = {', '.join(str(x) for x in args.top_ns)}",
                "```",
                "",
                f"Filtered hit rows before topN: `{len(hits)}`",
                f"Selected hit rows with ROI vectors: `{len(used_hits)}`",
                f"Dataset-level ROI profile rows: `{len(profile_rows)}`",
                f"Method-level ROI rows: `{len(agg_rows)}`",
                f"Dataset-pair ROI correlation rows: `{len(corr_rows)}`",
                "",
                "Main figures:",
                "",
                *[f"- `{path}`" for path in figure_paths],
                "",
                "Interpretation note:",
                "",
                "This is a P8 direct BOLD-mode ROI analysis. P10 ROI needs a separate P9-component backprojection step before it can be interpreted as brain-region ROI.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    print(f"P7 ROI mode profiles loaded: {len(roi_profiles)}")
    print(f"Filtered hit rows before topN: {len(hits)}")
    print(f"Selected hit rows with ROI vectors: {len(used_hits)}")
    print(f"Output dir: {args.output_dir}")
    print(f"Figure dir: {args.figure_dir}")
    for path in figure_paths:
        print(f"Figure: {path}")


if __name__ == "__main__":
    main()
