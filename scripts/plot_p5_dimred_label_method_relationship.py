#!/usr/bin/env python3
"""Plot P5 dimred method/k versus component process-label relationships."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

from summarize_pipeline8_cross_session_consistency import (
    COMPONENT_COUNTS,
    METHOD_ORDER,
    plot_heatmap,
)


DEFAULT_LABEL_TABLE = (
    Path("results")
    / "pipeline5_dimred_component_process_labels_current_rmsenv_adaptive_all_components"
    / "dimred_efun_process_labels.csv"
)
DEFAULT_RESULTS_DIR = Path("results") / "p5_dimred_label_method_relationship_current_rmsenv_adaptive"
DEFAULT_FIGURE_DIR = (
    Path(r"E:\DataPons_processed")
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p5"
    / "dimred_label_method_relationship_v1"
)

LABEL_ORDER = (
    "theta_selective_similar",
    "theta_selective_unequal",
    "ripple_selective_similar",
    "ripple_selective_unequal",
    "theta_ripple_joint",
    "mixed_theta_ripple",
    "mixed_theta_gamma",
    "mixed_ripple_gamma",
    "gamma_selective",
    "pan_event",
    "partial_or_inactive",
    "no_event_activity",
    "label_missing",
    "nonselective",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label-table", type=Path, default=DEFAULT_LABEL_TABLE)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def truthy(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def condition_short(row: dict[str, str]) -> str:
    existing = str(row.get("condition_short", "")).strip().lower()
    if existing:
        return existing
    condition = str(row.get("condition", "")).lower()
    if "complex" in condition or "csplit" in condition:
        return "csplit"
    if "abs" in condition:
        return "abs"
    return condition


def method_k(row: dict[str, str]) -> str:
    method = str(row.get("method", "")).lower()
    k = as_int(row.get("k") or row.get("component_count"))
    return f"{method}_k{k:02d}"


def method_k_rows() -> list[str]:
    return [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]


def normalize_label(label: object) -> str:
    text = str(label or "").strip()
    if text == "unlabeled":
        return "no_event_activity"
    return text or "label_missing"


def labels_present(rows: Iterable[dict[str, str]]) -> list[str]:
    found = {normalize_label(row.get("primary_process_label")) for row in rows}
    ordered = [label for label in LABEL_ORDER if label in found]
    extras = sorted(found.difference(ordered))
    return ordered + extras


def datasets_present(rows: Iterable[dict[str, str]]) -> list[str]:
    return sorted({str(row.get("dataset", "")) for row in rows if str(row.get("dataset", ""))})


def conditions_present(rows: Iterable[dict[str, str]]) -> list[str]:
    preferred = ["abs", "csplit"]
    found = {condition_short(row) for row in rows}
    return [c for c in preferred if c in found] + sorted(found.difference(preferred))


def coverage_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, object]]:
    counts: Counter[tuple[str, str, str, int]] = Counter()
    for row in rows:
        dataset = str(row.get("dataset", ""))
        condition = condition_short(row)
        method = str(row.get("method", "")).lower()
        k = as_int(row.get("k") or row.get("component_count"))
        if dataset and condition and method and k:
            counts[(dataset, condition, method, k)] += 1

    out: list[dict[str, object]] = []
    for dataset in datasets_present(rows):
        for condition in conditions_present(rows):
            for method in METHOD_ORDER:
                for k in COMPONENT_COUNTS:
                    n = counts[(dataset, condition, method, k)]
                    expected = k
                    fraction = n / expected if expected else math.nan
                    out.append(
                        {
                            "dataset": dataset,
                            "condition_short": condition,
                            "method": method,
                            "k": k,
                            "method_k": f"{method}_k{k:02d}",
                            "n_label_rows": n,
                            "expected_components": expected,
                            "coverage_fraction": min(fraction, 1.0),
                            "missing_fraction": max(0.0, 1.0 - fraction),
                        }
                    )
    return out


def label_composition_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, object]]:
    counts: Counter[tuple[str, str]] = Counter()
    totals: Counter[str] = Counter()
    for row in rows:
        condition = condition_short(row)
        mk = method_k(row)
        label = normalize_label(row.get("primary_process_label"))
        counts[(condition, mk, label)] += 1
        totals[(condition, mk)] += 1

    out: list[dict[str, object]] = []
    all_labels = labels_present(rows)
    for condition in conditions_present(rows):
        for mk in method_k_rows():
            total = totals[(condition, mk)]
            for label in all_labels:
                count = counts[(condition, mk, label)]
                out.append(
                    {
                        "condition_short": condition,
                        "method_k": mk,
                        "label": label,
                        "count": count,
                        "n_components": total,
                        "fraction": count / total if total else math.nan,
                    }
                )
    return out


def selective_yield_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, object]]:
    metric_names = (
        "theta_selective",
        "ripple_selective",
        "theta_ripple_joint",
        "any_theta_or_ripple_selective",
        "partial_or_inactive",
        "no_event_activity",
    )
    counts: Counter[tuple[str, str, int, str, str]] = Counter()
    totals: Counter[tuple[str, str, int, str]] = Counter()
    datasets_by_grid: dict[tuple[str, str, int, str], set[str]] = defaultdict(set)
    for row in rows:
        dataset = str(row.get("dataset", ""))
        condition = condition_short(row)
        method = str(row.get("method", "")).lower()
        k = as_int(row.get("k") or row.get("component_count"))
        label = normalize_label(row.get("primary_process_label"))
        theta = truthy(row.get("theta_selective"))
        ripple = truthy(row.get("ripple_selective"))
        metric_values = {
            "theta_selective": theta,
            "ripple_selective": ripple,
            "theta_ripple_joint": theta and ripple,
            "any_theta_or_ripple_selective": theta or ripple,
            "partial_or_inactive": label == "partial_or_inactive",
            "no_event_activity": label == "no_event_activity",
        }
        grid_key = (condition, method, k, dataset)
        totals[grid_key] += 1
        datasets_by_grid[(condition, method, k, "datasets")].add(dataset)
        for metric, value in metric_values.items():
            if value:
                counts[(condition, method, k, dataset, metric)] += 1

    out: list[dict[str, object]] = []
    for condition in conditions_present(rows):
        for method in METHOD_ORDER:
            for k in COMPONENT_COUNTS:
                datasets = sorted(datasets_by_grid[(condition, method, k, "datasets")])
                for metric in metric_names:
                    per_dataset_counts = [counts[(condition, method, k, dataset, metric)] for dataset in datasets]
                    per_dataset_totals = [totals[(condition, method, k, dataset)] for dataset in datasets]
                    n_components = sum(per_dataset_totals)
                    count = sum(per_dataset_counts)
                    out.append(
                        {
                            "condition_short": condition,
                            "method": method,
                            "k": k,
                            "method_k": f"{method}_k{k:02d}",
                            "metric": metric,
                            "count": count,
                            "n_components": n_components,
                            "fraction": count / n_components if n_components else math.nan,
                            "mean_count_per_dataset": (
                                sum(per_dataset_counts) / len(per_dataset_counts)
                                if per_dataset_counts
                                else math.nan
                            ),
                            "n_datasets": len(datasets),
                        }
                    )
    return out


def label_dataset_prevalence_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, object]]:
    label_dataset_sets: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    available_datasets: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        dataset = str(row.get("dataset", ""))
        condition = condition_short(row)
        mk = method_k(row)
        label = normalize_label(row.get("primary_process_label"))
        available_datasets[(condition, mk)].add(dataset)
        label_dataset_sets[(condition, mk, label)].add(dataset)

    out: list[dict[str, object]] = []
    all_labels = labels_present(rows)
    for condition in conditions_present(rows):
        for mk in method_k_rows():
            denom = len(available_datasets[(condition, mk)])
            for label in all_labels:
                n = len(label_dataset_sets[(condition, mk, label)])
                out.append(
                    {
                        "condition_short": condition,
                        "method_k": mk,
                        "label": label,
                        "n_datasets_with_label": n,
                        "n_available_datasets": denom,
                        "dataset_prevalence_fraction": n / denom if denom else math.nan,
                    }
                )
    return out


def plot_coverage(rows: Sequence[dict[str, object]], figure_dir: Path) -> None:
    for condition in sorted({str(r["condition_short"]) for r in rows}):
        matrix = {
            (str(row["dataset"]), str(row["method_k"])): float(row["coverage_fraction"])
            for row in rows
            if row["condition_short"] == condition
        }
        row_labels = sorted({str(row["dataset"]) for row in rows})
        plot_heatmap(
            matrix,
            row_labels,
            method_k_rows(),
            f"P5 all-component label coverage | {condition}",
            figure_dir / "01_coverage_qc" / f"p5_label_coverage__{condition}.png",
            value_format=".2f",
        )


def plot_label_composition(rows: Sequence[dict[str, object]], labels: Sequence[str], figure_dir: Path) -> None:
    for condition in sorted({str(r["condition_short"]) for r in rows}):
        matrix = {
            (str(row["method_k"]), str(row["label"])): float(row["fraction"])
            for row in rows
            if row["condition_short"] == condition
        }
        plot_heatmap(
            matrix,
            method_k_rows(),
            labels,
            f"P5 label composition by dimred method-k | {condition}",
            figure_dir / "02_label_composition" / f"p5_label_composition__{condition}.png",
            value_format=".2f",
        )


def plot_selective_yield(rows: Sequence[dict[str, object]], figure_dir: Path) -> None:
    metrics = sorted({str(row["metric"]) for row in rows})
    for condition in sorted({str(r["condition_short"]) for r in rows}):
        for metric in metrics:
            matrix_fraction: dict[tuple[str, str], float] = {}
            matrix_count: dict[tuple[str, str], float] = {}
            for row in rows:
                if row["condition_short"] != condition or row["metric"] != metric:
                    continue
                matrix_fraction[(str(row["method"]), f"k{int(row['k']):02d}")] = float(row["fraction"])
                matrix_count[(str(row["method"]), f"k{int(row['k']):02d}")] = float(row["mean_count_per_dataset"])
            plot_heatmap(
                matrix_fraction,
                list(METHOD_ORDER),
                [f"k{k:02d}" for k in COMPONENT_COUNTS],
                f"P5 selective yield fraction | {metric} | {condition}",
                figure_dir / "03_selective_yield" / f"p5_selective_yield_fraction__{metric}__{condition}.png",
                value_format=".2f",
            )
            plot_heatmap(
                matrix_count,
                list(METHOD_ORDER),
                [f"k{k:02d}" for k in COMPONENT_COUNTS],
                f"P5 selective yield mean count/dataset | {metric} | {condition}",
                figure_dir / "03_selective_yield" / f"p5_selective_yield_mean_count__{metric}__{condition}.png",
                value_format=".2f",
            )


def plot_label_prevalence(rows: Sequence[dict[str, object]], labels: Sequence[str], figure_dir: Path) -> None:
    for condition in sorted({str(r["condition_short"]) for r in rows}):
        matrix = {
            (str(row["method_k"]), str(row["label"])): float(row["dataset_prevalence_fraction"])
            for row in rows
            if row["condition_short"] == condition
        }
        plot_heatmap(
            matrix,
            method_k_rows(),
            labels,
            f"P5 label dataset prevalence by method-k | {condition}",
            figure_dir / "04_label_dataset_prevalence" / f"p5_label_dataset_prevalence__{condition}.png",
            value_format=".2f",
        )


def main() -> int:
    args = parse_args()
    rows = read_csv(args.label_table)
    if not rows:
        raise SystemExit(f"No label rows found: {args.label_table}")
    labels = labels_present(rows)

    cov = coverage_rows(rows)
    comp = label_composition_rows(rows)
    yld = selective_yield_rows(rows)
    prev = label_dataset_prevalence_rows(rows)

    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.results_dir / "p5_dimred_label_coverage.csv", cov)
    write_csv(args.results_dir / "p5_dimred_label_composition.csv", comp)
    write_csv(args.results_dir / "p5_dimred_selective_yield.csv", yld)
    write_csv(args.results_dir / "p5_dimred_label_dataset_prevalence.csv", prev)

    plot_coverage(cov, args.figure_dir)
    plot_label_composition(comp, labels, args.figure_dir)
    plot_selective_yield(yld, args.figure_dir)
    plot_label_prevalence(prev, labels, args.figure_dir)

    summary = [
        "# P5 dimred label-method relationship",
        "",
        f"- Label table: `{args.label_table}`",
        f"- Label rows: `{len(rows)}`",
        f"- Datasets: `{', '.join(datasets_present(rows))}`",
        f"- Conditions: `{', '.join(conditions_present(rows))}`",
        f"- Labels: `{', '.join(labels)}`",
        f"- Figure root: `{args.figure_dir}`",
        "",
        "These figures use all current P5 dimred component-label rows, not only P8/P10 top hits.",
    ]
    (args.results_dir / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
