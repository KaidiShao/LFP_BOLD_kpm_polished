#!/usr/bin/env python3
"""Compare e10gb1 P5 raw-vs-standardized dimred ripple selectivity."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Sequence

from summarize_pipeline8_cross_session_consistency import (
    COMPONENT_COUNTS,
    METHOD_ORDER,
    plot_heatmap,
)


DEFAULT_RAW_LABEL_TABLE = (
    Path("results")
    / "pipeline5_dimred_component_process_labels_current_rmsenv_adaptive_all_components"
    / "dimred_efun_process_labels.csv"
)
DEFAULT_STD_LABEL_TABLE = (
    Path("results")
    / "pipeline5_dimred_component_process_labels_e10gb1_standardize_rmsenv_adaptive_all_components"
    / "dimred_efun_process_labels.csv"
)
DEFAULT_RESULTS_DIR = Path("results") / "e10gb1_p5_standardize_ripple_probe"
DEFAULT_FIGURE_DIR = (
    Path(r"E:\DataPons_processed")
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p5"
    / "e10gb1_standardize_ripple_probe"
)

METRICS = (
    "ripple_active",
    "ripple_gamma_joint",
    "ripple_active_allow_gamma_theta_inactive",
    "ripple_active_allow_gamma_theta_not_selective",
    "ripple_selective",
    "ripple_selective_similar",
    "ripple_selective_unequal",
    "ripple_active_not_selective",
    "theta_active",
    "theta_selective",
    "gamma_active",
    "gamma_selective",
    "any_theta_or_ripple_selective",
    "partial_or_inactive",
    "no_event_activity",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-label-table", type=Path, default=DEFAULT_RAW_LABEL_TABLE)
    parser.add_argument("--std-label-table", type=Path, default=DEFAULT_STD_LABEL_TABLE)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--dataset", default="e10gb1")
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
    text = str(row.get("condition_short") or row.get("condition") or "").lower()
    if "complex" in text or "csplit" in text:
        return "csplit"
    if "abs" in text:
        return "abs"
    return text


def normalize_label(value: object) -> str:
    label = str(value or "").strip()
    if label == "unlabeled":
        return "no_event_activity"
    return label or "label_missing"


def row_metrics(row: dict[str, str]) -> dict[str, bool]:
    label = normalize_label(row.get("primary_process_label"))
    ripple_active = truthy(row.get("ripple_active"))
    ripple = truthy(row.get("ripple_selective"))
    theta_active = truthy(row.get("theta_active"))
    theta = truthy(row.get("theta_selective"))
    gamma_active = truthy(row.get("gamma_active"))
    gamma = truthy(row.get("gamma_selective"))
    ripple_similar = ripple and truthy(row.get("ripple_similar"))
    return {
        "ripple_active": ripple_active,
        "ripple_gamma_joint": ripple_active and gamma_active,
        "ripple_active_allow_gamma_theta_inactive": ripple_active and not theta_active,
        "ripple_active_allow_gamma_theta_not_selective": ripple_active and not theta,
        "ripple_selective": ripple,
        "ripple_selective_similar": ripple_similar,
        "ripple_selective_unequal": ripple and not ripple_similar,
        "ripple_active_not_selective": ripple_active and not ripple,
        "theta_active": theta_active,
        "theta_selective": theta,
        "gamma_active": gamma_active,
        "gamma_selective": gamma,
        "any_theta_or_ripple_selective": theta or ripple,
        "partial_or_inactive": label == "partial_or_inactive",
        "no_event_activity": label == "no_event_activity",
    }


def summarize_branch(rows: Sequence[dict[str, str]], branch: str, dataset: str) -> list[dict[str, object]]:
    counts: Counter[tuple[str, str, int, str]] = Counter()
    totals: Counter[tuple[str, str, int]] = Counter()
    for row in rows:
        if str(row.get("dataset", "")).lower() != dataset.lower():
            continue
        condition = condition_short(row)
        method = str(row.get("method", "")).lower()
        k = as_int(row.get("k") or row.get("component_count"))
        if condition not in {"abs", "csplit"} or method not in METHOD_ORDER or k not in COMPONENT_COUNTS:
            continue
        totals[(condition, method, k)] += 1
        metrics = row_metrics(row)
        for metric, value in metrics.items():
            if value:
                counts[(condition, method, k, metric)] += 1

    out: list[dict[str, object]] = []
    for condition in ("abs", "csplit"):
        for method in METHOD_ORDER:
            for k in COMPONENT_COUNTS:
                total = totals[(condition, method, k)]
                for metric in METRICS:
                    count = counts[(condition, method, k, metric)]
                    out.append(
                        {
                            "dataset": dataset,
                            "branch": branch,
                            "condition_short": condition,
                            "method": method,
                            "k": k,
                            "method_k": f"{method}_k{k:02d}",
                            "metric": metric,
                            "count": count,
                            "n_components": total,
                            "fraction": count / total if total else math.nan,
                        }
                    )
    return out


def comparison_rows(summary_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    by_key: dict[tuple[str, str, int, str, str], dict[str, object]] = {}
    for row in summary_rows:
        key = (
            str(row["condition_short"]),
            str(row["method"]),
            int(row["k"]),
            str(row["metric"]),
            str(row["branch"]),
        )
        by_key[key] = row

    out: list[dict[str, object]] = []
    for condition in ("abs", "csplit"):
        for method in METHOD_ORDER:
            for k in COMPONENT_COUNTS:
                for metric in METRICS:
                    raw = by_key.get((condition, method, k, metric, "raw"), {})
                    std = by_key.get((condition, method, k, metric, "standardize"), {})
                    raw_fraction = float(raw.get("fraction", math.nan))
                    std_fraction = float(std.get("fraction", math.nan))
                    raw_count = int(raw.get("count", 0) or 0)
                    std_count = int(std.get("count", 0) or 0)
                    out.append(
                        {
                            "condition_short": condition,
                            "method": method,
                            "k": k,
                            "method_k": f"{method}_k{k:02d}",
                            "metric": metric,
                            "raw_count": raw_count,
                            "standardize_count": std_count,
                            "raw_fraction": raw_fraction,
                            "standardize_fraction": std_fraction,
                            "delta_fraction": (
                                std_fraction - raw_fraction
                                if math.isfinite(raw_fraction) and math.isfinite(std_fraction)
                                else math.nan
                            ),
                        }
                    )
    return out


def plot_metric(summary_rows: Sequence[dict[str, object]], compare_rows: Sequence[dict[str, object]], metric: str, figure_dir: Path) -> None:
    for condition in ("abs", "csplit"):
        for branch in ("raw", "standardize"):
            matrix = {
                (str(row["method"]), f"k{int(row['k']):02d}"): float(row["fraction"])
                for row in summary_rows
                if row["condition_short"] == condition and row["branch"] == branch and row["metric"] == metric
            }
            plot_heatmap(
                matrix,
                list(METHOD_ORDER),
                [f"k{k:02d}" for k in COMPONENT_COUNTS],
                f"e10gb1 P5 {metric} fraction | {branch} | {condition}",
                figure_dir / f"{metric}_fraction__{branch}__{condition}.png",
                value_format=".2f",
            )

        delta_matrix = {
            (str(row["method"]), f"k{int(row['k']):02d}"): float(row["delta_fraction"])
            for row in compare_rows
            if row["condition_short"] == condition and row["metric"] == metric
        }
        plot_heatmap(
            delta_matrix,
            list(METHOD_ORDER),
            [f"k{k:02d}" for k in COMPONENT_COUNTS],
            f"e10gb1 P5 {metric}: standardize - raw | {condition}",
            figure_dir / f"{metric}_delta_standardize_minus_raw__{condition}.png",
            value_format=".2f",
        )


def top_delta_lines(
    compare: Sequence[dict[str, object]],
    metric: str,
    n: int = 8,
    reverse: bool = True,
) -> list[str]:
    rows = [r for r in compare if r["metric"] == metric and math.isfinite(float(r["delta_fraction"]))]
    rows = sorted(rows, key=lambda r: float(r["delta_fraction"]), reverse=reverse)
    lines = []
    for row in rows[:n]:
        lines.append(
            f"- {row['condition_short']} {row['method_k']}: "
            f"raw={float(row['raw_fraction']):.2f}, "
            f"std={float(row['standardize_fraction']):.2f}, "
            f"delta={float(row['delta_fraction']):+.2f} "
            f"(count {row['raw_count']} -> {row['standardize_count']})"
        )
    return lines


def main() -> int:
    args = parse_args()
    raw_rows = read_csv(args.raw_label_table)
    std_rows = read_csv(args.std_label_table)
    if not raw_rows:
        raise SystemExit(f"Missing raw label table: {args.raw_label_table}")
    if not std_rows:
        raise SystemExit(f"Missing standardized label table: {args.std_label_table}")

    summary = summarize_branch(raw_rows, "raw", args.dataset) + summarize_branch(std_rows, "standardize", args.dataset)
    compare = comparison_rows(summary)
    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.results_dir / "e10gb1_raw_vs_standardize_p5_label_summary.csv", summary)
    write_csv(args.results_dir / "e10gb1_raw_vs_standardize_p5_label_delta.csv", compare)

    for metric in (
        "ripple_active",
        "ripple_gamma_joint",
        "ripple_active_allow_gamma_theta_inactive",
        "ripple_active_allow_gamma_theta_not_selective",
        "ripple_selective",
        "ripple_selective_similar",
        "ripple_active_not_selective",
        "theta_active",
        "theta_selective",
        "gamma_active",
        "any_theta_or_ripple_selective",
        "no_event_activity",
    ):
        plot_metric(summary, compare, metric, args.figure_dir)

    md = [
        "# e10gb1 P5 standardized ripple-selectivity probe",
        "",
        f"- Raw label table: `{args.raw_label_table}`",
        f"- Standardized label table: `{args.std_label_table}`",
        f"- Figure root: `{args.figure_dir}`",
        "",
        "## Largest ripple_selective gains",
        *top_delta_lines(compare, "ripple_selective"),
        "",
        "## Largest ripple_selective losses",
        *top_delta_lines(compare, "ripple_selective", reverse=False),
        "",
        "## Largest ripple_active gains",
        *top_delta_lines(compare, "ripple_active"),
        "",
        "## Largest ripple_active_allow_gamma_theta_inactive gains",
        *top_delta_lines(compare, "ripple_active_allow_gamma_theta_inactive"),
        "",
        "## Largest ripple_gamma_joint gains",
        *top_delta_lines(compare, "ripple_gamma_joint"),
        "",
        "## Largest ripple_active_not_selective gains",
        *top_delta_lines(compare, "ripple_active_not_selective"),
        "",
        "## Largest ripple_selective_similar gains",
        *top_delta_lines(compare, "ripple_selective_similar"),
        "",
        "Positive delta means standardized branch has a larger fraction than raw for the same condition/method/k.",
    ]
    (args.results_dir / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
