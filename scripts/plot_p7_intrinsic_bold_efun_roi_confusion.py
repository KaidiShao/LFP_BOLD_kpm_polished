"""Plot all-mode P7 intrinsic BOLD efun ROI confusion matrices.

This is independent of P8/P10 xcorr and BLP density. It compares the original
P7 BOLD efun ROI profiles, split by BOLD observable, across all retained sorted
BOLD efun modes exported from the current P7 BOLD_POST runs.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean, median
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from summarize_pipeline8_cross_session_consistency import DEFAULT_PROCESSED_ROOT, write_csv


DEFAULT_RESULTS_DIR = Path("results") / "pipeline_roi_profile_consistency_current"
DEFAULT_PROFILE_CSV = DEFAULT_RESULTS_DIR / "p7_intrinsic_bold_efun_roi_profiles_long.csv"
DEFAULT_OUTPUT_DIR = DEFAULT_RESULTS_DIR / "p7_intrinsic_bold_efun_roi_confusion_20260603"
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p7_intrinsic_bold_efun_roi_confusion_20260603"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-csv", type=Path, default=DEFAULT_PROFILE_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--observables", nargs="*", default=[])
    parser.add_argument("--min-datasets", type=int, default=3)
    parser.add_argument("--min-rois", type=int, default=5)
    parser.add_argument("--block-size", type=int, default=10)
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
    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def safe_name(text: object) -> str:
    return (
        str(text)
        .replace("\\", "_")
        .replace("/", "_")
        .replace(" ", "_")
        .replace("|", "_")
        .replace(":", "_")
    )


def mode_sort_key(label: str) -> tuple[int, str]:
    digits = "".join(ch for ch in label if ch.isdigit())
    if digits:
        return (int(digits), label)
    return (10_000, label)


def mode_position(label: str) -> int:
    digits = "".join(ch for ch in label if ch.isdigit())
    return int(digits) if digits else 10_000


def block_label(label: str, block_size: int) -> str:
    pos = mode_position(label)
    if pos >= 10_000 or block_size <= 0:
        return label
    start = ((pos - 1) // block_size) * block_size + 1
    end = start + block_size - 1
    return f"{start:03d}-{end:03d}"


def read_profiles(path: Path, observables: set[str]) -> dict[tuple[str, str, str], dict[str, object]]:
    accum: dict[tuple[str, str, str], dict[str, object]] = {}
    roi_values: dict[tuple[str, str, str], dict[str, float]] = defaultdict(dict)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            obs = row.get("observable", "")
            if observables and obs not in observables:
                continue
            dataset = row.get("dataset", "")
            mode = row.get("mode_label", "")
            roi = row.get("roi_label", "")
            value = as_float(row.get("roi_value"))
            if not dataset or not obs or not mode or not roi or not math.isfinite(value):
                continue
            key = (obs, dataset, mode)
            roi_values[key][roi] = value
            accum.setdefault(
                key,
                {
                    "observable": obs,
                    "dataset": dataset,
                    "mode_label": mode,
                    "run_name": row.get("run_name", ""),
                    "raw_index": row.get("raw_index", ""),
                    "sorted_position": row.get("sorted_position", ""),
                },
            )

    for key, values in roi_values.items():
        accum[key]["roi_values"] = values
    return accum


def correlate_profiles(a: dict[str, float], b: dict[str, float], min_rois: int) -> tuple[float, int]:
    common = sorted(set(a) & set(b))
    if len(common) < min_rois:
        return math.nan, len(common)
    r = pearson([a[x] for x in common], [b[x] for x in common])
    return r, len(common)


def collect_records(profiles: dict[tuple[str, str, str], dict[str, object]], min_rois: int) -> list[dict[str, object]]:
    by_obs: dict[str, list[dict[str, object]]] = defaultdict(list)
    for item in profiles.values():
        by_obs[str(item["observable"])].append(item)

    records: list[dict[str, object]] = []
    for obs, items in by_obs.items():
        for a, b in combinations(items, 2):
            mode_a = str(a["mode_label"])
            mode_b = str(b["mode_label"])
            dataset_a = str(a["dataset"])
            dataset_b = str(b["dataset"])
            value, n_rois = correlate_profiles(
                a["roi_values"],  # type: ignore[arg-type]
                b["roi_values"],  # type: ignore[arg-type]
                min_rois,
            )
            if not math.isfinite(value):
                continue
            if dataset_a == dataset_b and mode_a == mode_b:
                continue
            if dataset_a == dataset_b:
                comparison = "cross_mode_same_dataset"
            elif mode_a == mode_b:
                comparison = "within_mode_cross_dataset"
            else:
                comparison = "cross_mode_cross_dataset"
            records.append(
                {
                    "observable": obs,
                    "comparison": comparison,
                    "dataset_a": dataset_a,
                    "dataset_b": dataset_b,
                    "mode_a": mode_a,
                    "mode_b": mode_b,
                    "mode_a_pos": mode_position(mode_a),
                    "mode_b_pos": mode_position(mode_b),
                    "n_rois": n_rois,
                    "roi_corr": value,
                }
            )
    return records


def summarize_records(records: Sequence[dict[str, object]], group_fields: Sequence[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[object, ...], list[float]] = defaultdict(list)
    example: dict[tuple[object, ...], dict[str, object]] = {}
    for row in records:
        key = tuple(row.get(field, "") for field in group_fields)
        value = as_float(row.get("roi_corr"))
        if math.isfinite(value):
            grouped[key].append(value)
            example.setdefault(key, {field: row.get(field, "") for field in group_fields})

    out: list[dict[str, object]] = []
    for key, vals in grouped.items():
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
    out.sort(key=lambda r: tuple(str(r.get(field, "")) for field in group_fields))
    return out


def mean_lookup(records: Sequence[dict[str, object]], comparison: str) -> dict[tuple[str, str], float]:
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in records:
        if row.get("comparison") != comparison:
            continue
        a = str(row.get("mode_a", ""))
        b = str(row.get("mode_b", ""))
        value = as_float(row.get("roi_corr"))
        if math.isfinite(value):
            grouped[(a, b)].append(value)
            grouped[(b, a)].append(value)
    return {key: mean(vals) for key, vals in grouped.items() if vals}


def plot_matrix(
    values: dict[tuple[str, str], float],
    labels: Sequence[str],
    title: str,
    out_file: Path,
    *,
    annotate_max_labels: int = 15,
) -> None:
    matrix = np.full((len(labels), len(labels)), np.nan, dtype=float)
    for i, row_label in enumerate(labels):
        for j, col_label in enumerate(labels):
            matrix[i, j] = values.get((row_label, col_label), math.nan)

    n = len(labels)
    fig_w = max(18.0, min(30.0, 0.28 * n + 7.0))
    fig_h = max(8.0, min(18.0, 0.16 * n + 5.0))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(matrix, cmap="coolwarm", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=65 if n > 20 else 30, ha="right")
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title(title)
    if len(labels) <= annotate_max_labels:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix[i, j]
                if math.isfinite(float(value)):
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=10)
    if n > 20:
        tick_step = max(1, int(math.ceil(n / 20)))
        keep = np.arange(n) % tick_step == 0
        ax.set_xticks(np.where(keep)[0])
        ax.set_xticklabels([labels[i] for i in np.where(keep)[0]], rotation=65, ha="right")
        ax.set_yticks(np.where(keep)[0])
        ax.set_yticklabels([labels[i] for i in np.where(keep)[0]])
    fig.colorbar(im, ax=ax, shrink=0.82, label="ROI profile correlation")
    fig.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def grouped_mean_lookup(
    records: Sequence[dict[str, object]],
    *,
    dataset_relation: str,
    block_size: int,
) -> dict[tuple[str, str], float]:
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in records:
        same_dataset = str(row.get("dataset_a", "")) == str(row.get("dataset_b", ""))
        if dataset_relation == "cross_dataset" and same_dataset:
            continue
        if dataset_relation == "same_dataset" and not same_dataset:
            continue
        a = block_label(str(row.get("mode_a", "")), block_size)
        b = block_label(str(row.get("mode_b", "")), block_size)
        value = as_float(row.get("roi_corr"))
        if math.isfinite(value):
            grouped[(a, b)].append(value)
            grouped[(b, a)].append(value)
    return {key: mean(vals) for key, vals in grouped.items() if vals}


def plot_observable_matrices(
    records: Sequence[dict[str, object]],
    figure_dir: Path,
    min_datasets: int,
    block_size: int,
) -> None:
    by_obs: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in records:
        by_obs[str(row.get("observable", ""))].append(row)

    for obs, obs_records in sorted(by_obs.items()):
        datasets = sorted(
            {
                str(row.get("dataset_a", ""))
                for row in obs_records
                if row.get("dataset_a")
            }
            | {
                str(row.get("dataset_b", ""))
                for row in obs_records
                if row.get("dataset_b")
            }
        )
        if len(datasets) < min_datasets:
            continue
        labels = sorted(
            {
                str(row.get("mode_a", ""))
                for row in obs_records
            }
            | {
                str(row.get("mode_b", ""))
                for row in obs_records
            },
            key=mode_sort_key,
        )
        cross_dataset_values = mean_lookup(obs_records, "within_mode_cross_dataset")
        cross_mode_values = mean_lookup(obs_records, "cross_mode_cross_dataset")
        for key, value in cross_mode_values.items():
            if key[0] != key[1]:
                cross_dataset_values[key] = value
        plot_matrix(
            cross_dataset_values,
            labels,
            f"P7 BOLD ROI confusion | {obs} | full | cross-dataset",
            figure_dir / f"p7_intrinsic_bold_efun_roi_confusion_cross_dataset__{safe_name(obs)}__full.png",
        )

        same_dataset_values = mean_lookup(obs_records, "cross_mode_same_dataset")
        plot_matrix(
            same_dataset_values,
            labels,
            f"P7 BOLD ROI similarity | {obs} | full | same dataset",
            figure_dir / f"p7_intrinsic_bold_efun_roi_confusion_same_dataset__{safe_name(obs)}__full.png",
        )

        block_labels = sorted(
            {block_label(label, block_size) for label in labels},
            key=mode_sort_key,
        )
        plot_matrix(
            grouped_mean_lookup(obs_records, dataset_relation="cross_dataset", block_size=block_size),
            block_labels,
            f"P7 BOLD ROI confusion | {obs} | blocks{block_size} | cross-dataset",
            figure_dir / f"p7_intrinsic_bold_efun_roi_confusion_cross_dataset__{safe_name(obs)}__blocks{block_size}.png",
        )
        plot_matrix(
            grouped_mean_lookup(obs_records, dataset_relation="same_dataset", block_size=block_size),
            block_labels,
            f"P7 BOLD ROI similarity | {obs} | blocks{block_size} | same dataset",
            figure_dir / f"p7_intrinsic_bold_efun_roi_confusion_same_dataset__{safe_name(obs)}__blocks{block_size}.png",
        )


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    profiles = read_profiles(args.profile_csv, set(args.observables))
    records = collect_records(profiles, args.min_rois)
    write_csv(args.output_dir / "p7_intrinsic_bold_efun_roi_pairwise_correlations.csv", records)
    write_csv(
        args.output_dir / "p7_intrinsic_bold_efun_roi_pairwise_summary.csv",
        summarize_records(records, ["observable", "comparison", "mode_a", "mode_b"]),
    )
    write_csv(
        args.output_dir / "p7_intrinsic_bold_efun_roi_comparison_summary.csv",
        summarize_records(records, ["observable", "comparison"]),
    )
    plot_observable_matrices(records, args.figure_dir, args.min_datasets, args.block_size)
    print(f"Profiles: {len(profiles)}")
    print(f"Pairwise records: {len(records)}")
    print(f"Output dir: {args.output_dir}")
    print(f"Figure dir: {args.figure_dir}")


if __name__ == "__main__":
    main()
