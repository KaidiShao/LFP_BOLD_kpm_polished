"""Analyze whether P8 subprocess ROI profiles are spatially separated.

The by-subprocess ROI consistency plots test whether each P5 subprocess has a
repeatable BOLD ROI profile across datasets.  This script asks the stricter
question:

    Is corr(theta ROI, ripple-gamma ROI) lower than the within-subprocess
    cross-dataset ROI correlation?

If the cross-subprocess correlations are similar to the within-subprocess
correlations, ROI profiles do not support spatial separation between the
subprocesses, even if each subprocess is individually consistent.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean, median
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np

from summarize_pipeline8_cross_session_consistency import (
    DEFAULT_PROCESSED_ROOT,
    RUN_TAG_LABELS,
    write_csv,
)


DEFAULT_RESULTS_DIR = Path("results") / "pipeline_roi_profile_consistency_current"
DEFAULT_PROFILE_CSV = DEFAULT_RESULTS_DIR / "p8_roi_profiles_by_subprocess_long.csv"
DEFAULT_OUTPUT_DIR = DEFAULT_RESULTS_DIR / "p8_roi_subprocess_separation_20260603"
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p8_roi_subprocess_separation_20260603"
)
DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01", "k13m17", "k13m23")
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")
LABEL_ORDER = ("theta", "ripple_gamma", "mixed", "inactive")
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
METHOD_ORDER = ("svd", "nmf", "mds", "umap")
K_ORDER = ("k03", "k04", "k05", "k06", "k07", "k08")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-csv", type=Path, default=DEFAULT_PROFILE_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--conditions", nargs="+", default=["csplit"])
    parser.add_argument("--label-groups", nargs="+", default=list(LABEL_ORDER))
    parser.add_argument("--method-k-scope", choices=["all", "p5-stable"], default="p5-stable")
    parser.add_argument("--min-datasets", type=int, default=5)
    parser.add_argument("--min-rois", type=int, default=5)
    parser.add_argument("--value-column", default="roi_value_weighted")
    parser.add_argument("--weight-column", default="mean_peak_abs_corr")
    return parser.parse_args()


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


def safe_slug(text: object) -> str:
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
    return row.get("density_method_k", "") in STABLE_METHOD_K


def read_aggregated_vectors(
    csv_path: Path,
    datasets: set[str],
    run_tags: set[str],
    conditions: set[str],
    label_groups: set[str],
    method_k_scope: str,
    value_column: str,
    weight_column: str,
) -> tuple[
    dict[tuple[str, ...], dict[str, dict[str, dict[str, float]]]],
    dict[tuple[str, ...], dict[str, str]],
]:
    # key -> label -> dataset -> roi -> [weighted sum, weight sum]
    sums: dict[tuple[str, ...], dict[str, dict[str, dict[str, list[float]]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: [0.0, 0.0])))
    )
    meta: dict[tuple[str, ...], dict[str, str]] = {}
    if not csv_path.is_file():
        return {}, {}

    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("pipeline", "").upper() != "P8":
                continue
            dataset = row.get("dataset", "")
            run_tag = row.get("run_tag", "")
            condition = row.get("density_condition", "")
            label = row.get("strict_label_group", "")
            if dataset.lower() not in datasets:
                continue
            if run_tag.lower() not in run_tags:
                continue
            if condition.lower() not in conditions:
                continue
            if label not in label_groups:
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

            key = (
                run_tag,
                row.get("observable", ""),
                row.get("feature_family", ""),
                condition,
                row.get("density_method_k", ""),
            )
            slot = sums[key][label][dataset][roi]
            slot[0] += value * weight
            slot[1] += weight
            meta.setdefault(
                key,
                {
                    "pipeline": "P8",
                    "run_tag": run_tag,
                    "observable": row.get("observable", ""),
                    "feature_family": row.get("feature_family", ""),
                    "density_condition": condition,
                    "density_method_k": row.get("density_method_k", ""),
                    "method_k_scope": method_k_scope,
                },
            )

    vectors: dict[tuple[str, ...], dict[str, dict[str, dict[str, float]]]] = defaultdict(lambda: defaultdict(dict))
    for key, label_map in sums.items():
        for label, dataset_map in label_map.items():
            for dataset, roi_map in dataset_map.items():
                vectors[key][label][dataset] = {
                    roi: sw[0] / sw[1]
                    for roi, sw in roi_map.items()
                    if sw[1] > 0 and math.isfinite(sw[0])
                }
    return dict(vectors), meta


def corr_for_vectors(v1: dict[str, float], v2: dict[str, float], min_rois: int) -> tuple[float, int]:
    rois = sorted(set(v1) & set(v2))
    if len(rois) < min_rois:
        return math.nan, len(rois)
    return pearson([v1[r] for r in rois], [v2[r] for r in rois]), len(rois)


def collect_correlation_records(
    vectors: dict[tuple[str, ...], dict[str, dict[str, dict[str, float]]]],
    meta: dict[tuple[str, ...], dict[str, str]],
    labels: Sequence[str],
    min_datasets: int,
    min_rois: int,
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for key, label_map in vectors.items():
        info = meta[key]
        for label in labels:
            dataset_map = label_map.get(label, {})
            datasets = sorted(dataset_map)
            if len(datasets) < min_datasets:
                continue
            for d1, d2 in combinations(datasets, 2):
                r, n_rois = corr_for_vectors(dataset_map[d1], dataset_map[d2], min_rois)
                if not math.isfinite(r):
                    continue
                records.append(
                    {
                        **info,
                        "comparison": "within_subprocess_cross_dataset",
                        "label_a": label,
                        "label_b": label,
                        "dataset_a": d1,
                        "dataset_b": d2,
                        "same_dataset": "false",
                        "n_rois": n_rois,
                        "roi_corr": r,
                    }
                )

        for i, label_a in enumerate(labels):
            map_a = label_map.get(label_a, {})
            for label_b in labels[i + 1 :]:
                map_b = label_map.get(label_b, {})
                common_datasets = sorted(set(map_a) & set(map_b))
                if len(common_datasets) >= min_datasets:
                    for dataset in common_datasets:
                        r, n_rois = corr_for_vectors(map_a[dataset], map_b[dataset], min_rois)
                        if not math.isfinite(r):
                            continue
                        records.append(
                            {
                                **info,
                                "comparison": "cross_subprocess_same_dataset",
                                "label_a": label_a,
                                "label_b": label_b,
                                "dataset_a": dataset,
                                "dataset_b": dataset,
                                "same_dataset": "true",
                                "n_rois": n_rois,
                                "roi_corr": r,
                            }
                        )
                if len(map_a) >= min_datasets and len(map_b) >= min_datasets:
                    for d1, v1 in sorted(map_a.items()):
                        for d2, v2 in sorted(map_b.items()):
                            if d1 == d2:
                                continue
                            r, n_rois = corr_for_vectors(v1, v2, min_rois)
                            if not math.isfinite(r):
                                continue
                            records.append(
                                {
                                    **info,
                                    "comparison": "cross_subprocess_cross_dataset",
                                    "label_a": label_a,
                                    "label_b": label_b,
                                    "dataset_a": d1,
                                    "dataset_b": d2,
                                    "same_dataset": "false",
                                    "n_rois": n_rois,
                                    "roi_corr": r,
                                }
                            )
    return records


def summarize_records(records: Sequence[dict[str, object]], group_fields: Sequence[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[object, ...], list[float]] = defaultdict(list)
    example: dict[tuple[object, ...], dict[str, object]] = {}
    for row in records:
        key = tuple(row.get(field, "") for field in group_fields)
        val = as_float(row.get("roi_corr"))
        if not math.isfinite(val):
            continue
        grouped[key].append(val)
        example.setdefault(key, {field: row.get(field, "") for field in group_fields})

    out: list[dict[str, object]] = []
    for key, vals in grouped.items():
        vals = [v for v in vals if math.isfinite(v)]
        if not vals:
            continue
        item = dict(example[key])
        item.update(
            {
                "n_corr": len(vals),
                "mean_roi_corr": f"{mean(vals):.10g}",
                "median_roi_corr": f"{median(vals):.10g}",
                "min_roi_corr": f"{min(vals):.10g}",
                "max_roi_corr": f"{max(vals):.10g}",
            }
        )
        out.append(item)
    out.sort(key=lambda r: tuple(str(r.get(f, "")) for f in group_fields))
    return out


def mean_lookup(rows: Sequence[dict[str, object]], comparison: str) -> dict[tuple[str, str], float]:
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        if row.get("comparison") != comparison:
            continue
        a = str(row.get("label_a", ""))
        b = str(row.get("label_b", ""))
        value = as_float(row.get("roi_corr"))
        if math.isfinite(value):
            grouped[(a, b)].append(value)
            grouped[(b, a)].append(value)
    return {k: mean(v) for k, v in grouped.items() if v}


def filter_records(records: Sequence[dict[str, object]], field: str, value: str) -> list[dict[str, object]]:
    return [row for row in records if str(row.get(field, "")) == value]


def safe_name(text: object) -> str:
    return (
        str(text)
        .replace("\\", "_")
        .replace("/", "_")
        .replace(" ", "_")
        .replace("|", "_")
        .replace(":", "_")
    )


def display_observable(row: dict[str, object]) -> str:
    return RUN_TAG_LABELS.get(str(row.get("run_tag", "")), str(row.get("observable", "")))


def split_method_k(method_k: object) -> tuple[str, str]:
    text = str(method_k)
    if "_" not in text:
        return text, ""
    method, kval = text.split("_", 1)
    return method, kval


def label_sort_key(label: str, preferred: Sequence[str]) -> tuple[int, str]:
    try:
        return (preferred.index(label), label)
    except ValueError:
        return (len(preferred), label)


def plot_feature_confusion_matrices(
    records: Sequence[dict[str, object]],
    label_groups: Sequence[str],
    feature_family: str,
    figure_dir: Path,
) -> None:
    feature_records = filter_records(records, "feature_family", feature_family)
    if not feature_records:
        return

    cross_dataset_values = mean_lookup(
        [
            r
            for r in feature_records
            if r.get("comparison") in {"within_subprocess_cross_dataset", "cross_subprocess_cross_dataset"}
        ],
        "within_subprocess_cross_dataset",
    )
    cross_values = mean_lookup(feature_records, "cross_subprocess_cross_dataset")
    for key, value in cross_values.items():
        if key[0] != key[1]:
            cross_dataset_values[key] = value
    plot_matrix(
        cross_dataset_values,
        label_groups,
        f"P8 ROI subprocess confusion matrix | {feature_family}",
        figure_dir / f"01_p8_roi_subprocess_confusion_cross_dataset__{feature_family}.png",
        vmin=0.0,
        vmax=1.0,
    )

    same_values = mean_lookup(feature_records, "cross_subprocess_same_dataset")
    plot_matrix(
        same_values,
        label_groups,
        f"P8 ROI cross-subprocess similarity | same dataset | {feature_family}",
        figure_dir / f"02_p8_roi_cross_subprocess_same_dataset__{feature_family}.png",
        vmin=0.0,
        vmax=1.0,
    )


def plot_confusion_for_subset(
    records: Sequence[dict[str, object]],
    label_groups: Sequence[str],
    title_context: str,
    out_prefix: str,
    figure_dir: Path,
) -> None:
    if not records:
        return

    cross_dataset_values = mean_lookup(
        [
            r
            for r in records
            if r.get("comparison") in {"within_subprocess_cross_dataset", "cross_subprocess_cross_dataset"}
        ],
        "within_subprocess_cross_dataset",
    )
    cross_values = mean_lookup(records, "cross_subprocess_cross_dataset")
    for key, value in cross_values.items():
        if key[0] != key[1]:
            cross_dataset_values[key] = value
    plot_matrix(
        cross_dataset_values,
        label_groups,
        f"P8 ROI subprocess confusion | {title_context}",
        figure_dir / f"01_{out_prefix}_confusion_cross_dataset.png",
        vmin=0.0,
        vmax=1.0,
    )

    same_values = mean_lookup(records, "cross_subprocess_same_dataset")
    plot_matrix(
        same_values,
        label_groups,
        f"P8 ROI cross-subprocess same-dataset | {title_context}",
        figure_dir / f"02_{out_prefix}_same_dataset.png",
        vmin=0.0,
        vmax=1.0,
    )


def plot_observable_feature_confusion_matrices(
    records: Sequence[dict[str, object]],
    label_groups: Sequence[str],
    figure_dir: Path,
) -> None:
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in records:
        key = (
            str(row.get("run_tag", "")),
            display_observable(row),
            str(row.get("feature_family", "")),
        )
        groups[key].append(row)

    out_dir = figure_dir / "by_observable_feature"
    for (_run_tag, obs, feature), subset in sorted(groups.items(), key=lambda item: (item[0][1], item[0][2])):
        title_context = f"{obs} | {feature}"
        out_prefix = f"p8_roi_subprocess__{safe_name(obs)}__{safe_name(feature)}"
        plot_confusion_for_subset(subset, label_groups, title_context, out_prefix, out_dir)


def plot_matrix(
    values: dict[tuple[str, str], float],
    labels: Sequence[str],
    title: str,
    out_file: Path,
    vmin: float = 0.0,
    vmax: float = 1.0,
) -> None:
    matrix = np.full((len(labels), len(labels)), np.nan, dtype=float)
    for i, row_label in enumerate(labels):
        for j, col_label in enumerate(labels):
            matrix[i, j] = values.get((row_label, col_label), math.nan)

    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    im = ax.imshow(matrix, cmap="coolwarm", vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title(title)
    for i in range(len(labels)):
        for j in range(len(labels)):
            value = matrix[i, j]
            if math.isfinite(float(value)):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.82, label="ROI profile correlation")
    fig.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def plot_method_k_separation_heatmaps(rows: Sequence[dict[str, object]], figure_dir: Path) -> None:
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = (
            str(row.get("run_tag", "")),
            RUN_TAG_LABELS.get(str(row.get("run_tag", "")), str(row.get("observable", ""))),
            str(row.get("feature_family", "")),
        )
        groups[key].append(row)

    metrics = [
        ("within_mean", "within theta/RG"),
        ("theta_rg_cross_dataset_mean", "theta vs RG cross-dataset"),
        ("separation_cross_dataset", "within - cross"),
    ]
    out_dir = figure_dir / "by_observable_feature_method_k"
    out_dir.mkdir(parents=True, exist_ok=True)

    for (_run_tag, obs, feature), subset in sorted(groups.items(), key=lambda item: (item[0][1], item[0][2])):
        method_labels = sorted(
            {split_method_k(row.get("density_method_k", ""))[0] for row in subset},
            key=lambda label: label_sort_key(label, METHOD_ORDER),
        )
        k_labels = sorted(
            {split_method_k(row.get("density_method_k", ""))[1] for row in subset},
            key=lambda label: label_sort_key(label, K_ORDER),
        )
        if not method_labels or not k_labels:
            continue

        fig, axes = plt.subplots(1, len(metrics), figsize=(5.2 * len(metrics), 4.2), constrained_layout=True)
        if len(metrics) == 1:
            axes = [axes]
        for ax, (metric, metric_title) in zip(axes, metrics):
            matrix = np.full((len(method_labels), len(k_labels)), np.nan, dtype=float)
            for row in subset:
                method, kval = split_method_k(row.get("density_method_k", ""))
                if method not in method_labels or kval not in k_labels:
                    continue
                value = as_float(row.get(metric))
                if math.isfinite(value):
                    matrix[method_labels.index(method), k_labels.index(kval)] = value

            if metric == "separation_cross_dataset":
                im = ax.imshow(matrix, cmap="coolwarm", vmin=-0.08, vmax=0.08)
                cbar_label = "ROI separation"
            else:
                im = ax.imshow(matrix, cmap="coolwarm", vmin=0.0, vmax=1.0)
                cbar_label = "ROI correlation"
            ax.set_title(metric_title)
            ax.set_xticks(range(len(k_labels)))
            ax.set_xticklabels(k_labels)
            ax.set_yticks(range(len(method_labels)))
            ax.set_yticklabels(method_labels)
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    value = matrix[i, j]
                    if math.isfinite(float(value)):
                        ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=9)
            fig.colorbar(im, ax=ax, shrink=0.78, label=cbar_label)

        fig.suptitle(f"P8 theta/RG ROI separation by method-k | {obs} | {feature}", fontsize=13)
        out_file = out_dir / f"p8_theta_rg_roi_separation_method_k__{safe_name(obs)}__{safe_name(feature)}.png"
        fig.savefig(out_file, dpi=180)
        plt.close(fig)


def plot_theta_rg_bars(rows: Sequence[dict[str, object]], out_file: Path) -> None:
    # One row per observable/feature, comparing within-subprocess vs theta-RG.
    grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in rows:
        obs = RUN_TAG_LABELS.get(str(row.get("run_tag", "")), str(row.get("observable", "")))
        feature = str(row.get("feature_family", ""))
        comparison = str(row.get("comparison", ""))
        labels = {str(row.get("label_a", "")), str(row.get("label_b", ""))}
        value = as_float(row.get("roi_corr"))
        if not math.isfinite(value):
            continue
        if comparison == "within_subprocess_cross_dataset" and row.get("label_a") in {"theta", "ripple_gamma"}:
            metric = f"within_{row.get('label_a')}"
        elif comparison == "cross_subprocess_same_dataset" and labels == {"theta", "ripple_gamma"}:
            metric = "theta_rg_same_dataset"
        elif comparison == "cross_subprocess_cross_dataset" and labels == {"theta", "ripple_gamma"}:
            metric = "theta_rg_cross_dataset"
        else:
            continue
        grouped[(obs, feature, metric)].append(value)

    row_labels = sorted({(obs, feature) for obs, feature, _metric in grouped})
    metrics = ["within_theta", "within_ripple_gamma", "theta_rg_same_dataset", "theta_rg_cross_dataset"]
    matrix = np.full((len(row_labels), len(metrics)), np.nan, dtype=float)
    for i, (obs, feature) in enumerate(row_labels):
        for j, metric in enumerate(metrics):
            vals = grouped.get((obs, feature, metric), [])
            if vals:
                matrix[i, j] = mean(vals)

    fig_h = max(4.5, 0.5 * len(row_labels) + 1.5)
    fig, ax = plt.subplots(figsize=(10, fig_h))
    im = ax.imshow(matrix, aspect="auto", cmap="coolwarm", vmin=0.0, vmax=1.0)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels([f"{obs} | {feature}" for obs, feature in row_labels])
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(metrics, rotation=25, ha="right")
    ax.set_title("Theta/RG ROI separation check | mean correlations")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if math.isfinite(float(value)):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.82, label="ROI profile correlation")
    fig.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def summarize_theta_rg_separation(records: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in records:
        labels = {str(row.get("label_a", "")), str(row.get("label_b", ""))}
        key = (
            str(row.get("run_tag", "")),
            str(row.get("observable", "")),
            str(row.get("feature_family", "")),
            str(row.get("density_condition", "")),
            str(row.get("density_method_k", "")),
        )
        value = as_float(row.get("roi_corr"))
        if not math.isfinite(value):
            continue
        if row.get("comparison") == "within_subprocess_cross_dataset" and row.get("label_a") in {"theta", "ripple_gamma"}:
            grouped[key][f"within_{row.get('label_a')}"].append(value)
        elif row.get("comparison") == "cross_subprocess_same_dataset" and labels == {"theta", "ripple_gamma"}:
            grouped[key]["theta_rg_same_dataset"].append(value)
        elif row.get("comparison") == "cross_subprocess_cross_dataset" and labels == {"theta", "ripple_gamma"}:
            grouped[key]["theta_rg_cross_dataset"].append(value)

    rows: list[dict[str, object]] = []
    for key, values in grouped.items():
        wt = values.get("within_theta", [])
        wrg = values.get("within_ripple_gamma", [])
        same = values.get("theta_rg_same_dataset", [])
        cross = values.get("theta_rg_cross_dataset", [])
        if not wt or not wrg or not same:
            continue
        within_mean = mean([mean(wt), mean(wrg)])
        same_mean = mean(same)
        cross_mean = mean(cross) if cross else math.nan
        rows.append(
            {
                "run_tag": key[0],
                "observable": key[1],
                "feature_family": key[2],
                "density_condition": key[3],
                "density_method_k": key[4],
                "within_theta_mean": f"{mean(wt):.10g}",
                "within_ripple_gamma_mean": f"{mean(wrg):.10g}",
                "within_mean": f"{within_mean:.10g}",
                "theta_rg_same_dataset_mean": f"{same_mean:.10g}",
                "theta_rg_cross_dataset_mean": f"{cross_mean:.10g}" if math.isfinite(cross_mean) else "",
                "separation_same_dataset": f"{within_mean - same_mean:.10g}",
                "separation_cross_dataset": f"{within_mean - cross_mean:.10g}" if math.isfinite(cross_mean) else "",
                "n_within_theta": len(wt),
                "n_within_ripple_gamma": len(wrg),
                "n_theta_rg_same_dataset": len(same),
                "n_theta_rg_cross_dataset": len(cross),
            }
        )
    rows.sort(key=lambda r: as_float(r.get("separation_same_dataset")), reverse=True)
    return rows


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    vectors, meta = read_aggregated_vectors(
        args.profile_csv,
        {d.lower() for d in args.datasets},
        {r.lower() for r in args.run_tags},
        {c.lower() for c in args.conditions},
        set(args.label_groups),
        args.method_k_scope,
        args.value_column,
        args.weight_column,
    )
    records = collect_correlation_records(
        vectors,
        meta,
        args.label_groups,
        args.min_datasets,
        args.min_rois,
    )
    write_csv(args.output_dir / "p8_roi_subprocess_pairwise_correlations.csv", records)
    write_csv(
        args.output_dir / "p8_roi_subprocess_pairwise_summary.csv",
        summarize_records(records, ["comparison", "label_a", "label_b"]),
    )
    write_csv(
        args.output_dir / "p8_roi_subprocess_pairwise_by_observable.csv",
        summarize_records(records, ["comparison", "run_tag", "observable", "feature_family", "label_a", "label_b"]),
    )
    write_csv(
        args.output_dir / "p8_roi_subprocess_pairwise_by_observable_method_k.csv",
        summarize_records(
            records,
            ["comparison", "run_tag", "observable", "feature_family", "density_method_k", "label_a", "label_b"],
        ),
    )
    write_csv(
        args.output_dir / "p8_roi_subprocess_pairwise_summary_by_feature.csv",
        summarize_records(records, ["comparison", "feature_family", "label_a", "label_b"]),
    )
    theta_rg_rows = summarize_theta_rg_separation(records)
    write_csv(args.output_dir / "p8_theta_rg_roi_separation_by_method_k.csv", theta_rg_rows)

    feature_families = sorted({str(row.get("feature_family", "")) for row in records if row.get("feature_family")})
    for feature_family in feature_families:
        plot_feature_confusion_matrices(records, args.label_groups, feature_family, args.figure_dir)
    plot_observable_feature_confusion_matrices(records, args.label_groups, args.figure_dir)
    plot_method_k_separation_heatmaps(theta_rg_rows, args.figure_dir)
    plot_theta_rg_bars(records, args.figure_dir / "03_p8_theta_rg_roi_separation_by_observable.png")

    print(f"Keys: {len(vectors)}")
    print(f"Pairwise correlation rows: {len(records)}")
    print(f"Theta/RG separation rows: {len(theta_rg_rows)}")
    print(f"Output dir: {args.output_dir}")
    print(f"Figure dir: {args.figure_dir}")


if __name__ == "__main__":
    main()
