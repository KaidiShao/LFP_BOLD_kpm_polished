#!/usr/bin/env python3
"""Compare e10gb1 raw-vs-standardized P5 with relaxed ripple/gamma criteria."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Sequence

from summarize_pipeline8_cross_session_consistency import COMPONENT_COUNTS, METHOD_ORDER, plot_heatmap


DEFAULT_RAW_ACTIVITY = (
    Path("results")
    / "peak_event_family_component_activity_current_rmsenv_adaptive"
    / "component_family_activity.csv"
)
DEFAULT_STD_ACTIVITY = (
    Path("results")
    / "peak_event_family_component_activity_e10gb1_standardize_rmsenv_adaptive"
    / "component_family_activity.csv"
)
DEFAULT_RESULTS_DIR = Path("results") / "e10gb1_p5_relaxed_ripple_gamma_probe"
DEFAULT_FIGURE_DIR = (
    Path(r"E:\DataPons_processed")
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p5"
    / "e10gb1_relaxed_ripple_gamma_probe"
)

METRICS = (
    "strict_ripple_selective",
    "relaxed_ripple_allow_gamma_no_pure_theta",
    "ripple_gamma_joint",
    "ripple_gamma_joint_no_pure_theta",
    "ripple_active",
    "pure_theta_active",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-activity", type=Path, default=DEFAULT_RAW_ACTIVITY)
    parser.add_argument("--std-activity", type=Path, default=DEFAULT_STD_ACTIVITY)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--dataset", default="e10gb1")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
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


def split_items(value: object) -> set[str]:
    return {item.strip().lower() for item in str(value or "").split(";") if item.strip()}


def condition_short(condition: str) -> str:
    text = condition.lower()
    if "complex" in text or "csplit" in text:
        return "csplit"
    if "abs" in text:
        return "abs"
    return text


def component_metrics(family_rows: Sequence[dict[str, str]]) -> dict[str, bool]:
    by_family = {str(row.get("event_family", "")).lower(): row for row in family_rows}
    ripple = by_family.get("ripple", {})
    gamma = by_family.get("gamma", {})
    theta = by_family.get("theta", {})

    ripple_active = truthy(ripple.get("all_family_events_active"))
    strict_ripple = ripple_active and truthy(ripple.get("selective_against_forbidden_events"))
    gamma_active = truthy(gamma.get("all_family_events_active"))
    theta_active_events = split_items(theta.get("active_events"))
    pure_theta_active = "theta" in theta_active_events

    return {
        "strict_ripple_selective": strict_ripple,
        "relaxed_ripple_allow_gamma_no_pure_theta": ripple_active and not pure_theta_active,
        "ripple_gamma_joint": ripple_active and gamma_active,
        "ripple_gamma_joint_no_pure_theta": ripple_active and gamma_active and not pure_theta_active,
        "ripple_active": ripple_active,
        "pure_theta_active": pure_theta_active,
    }


def summarize_branch(rows: Sequence[dict[str, str]], branch: str, dataset: str) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, int, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if str(row.get("dataset", "")).lower() != dataset.lower():
            continue
        condition = condition_short(str(row.get("condition", "")))
        method = str(row.get("method", "")).lower()
        k = as_int(row.get("component_count"))
        comp = as_int(row.get("component_idx"))
        if condition and method and k and comp:
            groups[(condition, method, k, comp)].append(row)

    counts: Counter[tuple[str, str, int, str]] = Counter()
    totals: Counter[tuple[str, str, int]] = Counter()
    for (condition, method, k, _comp), group_rows in groups.items():
        totals[(condition, method, k)] += 1
        metrics = component_metrics(group_rows)
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


def compare_rows(summary: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    by_key = {
        (
            row["condition_short"],
            row["method"],
            int(row["k"]),
            row["metric"],
            row["branch"],
        ): row
        for row in summary
    }
    out: list[dict[str, object]] = []
    for condition in ("abs", "csplit"):
        for method in METHOD_ORDER:
            for k in COMPONENT_COUNTS:
                for metric in METRICS:
                    raw = by_key.get((condition, method, k, metric, "raw"), {})
                    std = by_key.get((condition, method, k, metric, "standardize"), {})
                    raw_fraction = float(raw.get("fraction", math.nan))
                    std_fraction = float(std.get("fraction", math.nan))
                    out.append(
                        {
                            "condition_short": condition,
                            "method": method,
                            "k": k,
                            "method_k": f"{method}_k{k:02d}",
                            "metric": metric,
                            "raw_count": raw.get("count", 0),
                            "standardize_count": std.get("count", 0),
                            "raw_n_components": raw.get("n_components", 0),
                            "standardize_n_components": std.get("n_components", 0),
                            "raw_fraction": raw_fraction,
                            "standardize_fraction": std_fraction,
                            "delta_fraction": std_fraction - raw_fraction
                            if math.isfinite(raw_fraction) and math.isfinite(std_fraction)
                            else math.nan,
                        }
                    )
    return out


def plot_metric(summary: Sequence[dict[str, object]], compare: Sequence[dict[str, object]], metric: str, fig_dir: Path) -> None:
    for branch in ("raw", "standardize"):
        for condition in ("abs", "csplit"):
            matrix = {
                (str(row["method"]), f"k{int(row['k']):02d}"): float(row["fraction"])
                for row in summary
                if row["branch"] == branch
                and row["condition_short"] == condition
                and row["metric"] == metric
            }
            plot_heatmap(
                matrix,
                list(METHOD_ORDER),
                [f"k{k:02d}" for k in COMPONENT_COUNTS],
                f"e10gb1 P5 {metric} | {branch} | {condition}",
                fig_dir / f"{metric}_fraction__{branch}__{condition}.png",
                value_format=".2f",
            )
    for condition in ("abs", "csplit"):
        matrix = {
            (str(row["method"]), f"k{int(row['k']):02d}"): float(row["delta_fraction"])
            for row in compare
            if row["condition_short"] == condition and row["metric"] == metric
        }
        plot_heatmap(
            matrix,
            list(METHOD_ORDER),
            [f"k{k:02d}" for k in COMPONENT_COUNTS],
            f"e10gb1 P5 {metric}: standardize - raw | {condition}",
            fig_dir / f"{metric}_delta_standardize_minus_raw__{condition}.png",
            value_format=".2f",
        )


def top_lines(compare: Sequence[dict[str, object]], metric: str, n: int = 10) -> list[str]:
    rows = [r for r in compare if r["metric"] == metric and math.isfinite(float(r["delta_fraction"]))]
    rows = sorted(rows, key=lambda r: float(r["delta_fraction"]), reverse=True)
    return [
        (
            f"- {row['condition_short']} {row['method_k']}: "
            f"raw={float(row['raw_fraction']):.2f}, "
            f"std={float(row['standardize_fraction']):.2f}, "
            f"delta={float(row['delta_fraction']):+.2f} "
            f"(count {row['raw_count']} -> {row['standardize_count']})"
        )
        for row in rows[:n]
    ]


def main() -> int:
    args = parse_args()
    raw_rows = read_csv(args.raw_activity)
    std_rows = read_csv(args.std_activity)
    summary = summarize_branch(raw_rows, "raw", args.dataset) + summarize_branch(std_rows, "standardize", args.dataset)
    compare = compare_rows(summary)
    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.results_dir / "e10gb1_relaxed_ripple_gamma_summary.csv", summary)
    write_csv(args.results_dir / "e10gb1_relaxed_ripple_gamma_delta.csv", compare)
    for metric in METRICS:
        plot_metric(summary, compare, metric, args.figure_dir)

    md = [
        "# e10gb1 P5 relaxed ripple/gamma probe",
        "",
        f"- Raw activity table: `{args.raw_activity}`",
        f"- Standardized activity table: `{args.std_activity}`",
        f"- Figure root: `{args.figure_dir}`",
        "",
        "Definitions:",
        "- `relaxed_ripple_allow_gamma_no_pure_theta`: ripple and sharp-wave-ripple active; pure theta event inactive; gamma/theta-gamma allowed.",
        "- `ripple_gamma_joint`: ripple family and gamma family both active.",
        "- `ripple_gamma_joint_no_pure_theta`: ripple/gamma joint with pure theta inactive.",
        "",
        "## Largest relaxed ripple gains",
        *top_lines(compare, "relaxed_ripple_allow_gamma_no_pure_theta"),
        "",
        "## Largest ripple-gamma joint gains",
        *top_lines(compare, "ripple_gamma_joint"),
        "",
        "## Largest ripple-gamma joint no-pure-theta gains",
        *top_lines(compare, "ripple_gamma_joint_no_pure_theta"),
    ]
    (args.results_dir / "summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
