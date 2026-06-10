#!/usr/bin/env python3
"""Summarize P5 component activity across event-name families.

For each Pipeline 5 peak-state stats CSV, this script asks:
1. Is any single component active in every event whose label contains
   theta/gamma/ripple?
2. If so, are the active amplitudes similar across those events?

Default event families:
  theta  -> labels containing "theta"  (for example theta, theta-gamma)
  gamma  -> labels containing "gamma"  (for example gamma, theta-gamma)
  ripple -> labels containing "ripple" (for example ripple, sharp-wave-ripple)
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import mean, pstdev
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_OUTPUT_DIR = Path("results") / "peak_event_family_component_activity_current_maxabs"
DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01")
DEFAULT_CONDITIONS = ("abs_projected_vlambda", "complex_split_projected_vlambda")
DEFAULT_METHODS = ("svd", "nmf", "mds", "umap")
DEFAULT_COMPONENT_COUNTS = tuple(range(3, 9))
DEFAULT_EVENT_KEYWORDS = ("theta", "gamma", "ripple")
DEFAULT_EXPLICIT_FAMILY_EVENTS = {
    "theta": ("theta", "theta-gamma", "sharp-wave-ripple"),
    "gamma": ("gamma", "theta-gamma", "sharp-wave-ripple"),
    "ripple": ("ripple", "sharp-wave-ripple"),
}


@dataclass(frozen=True)
class RunSpec:
    dataset: str
    condition: str
    method: str
    component_count: int
    method_tag: str
    stats_file: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--peak-stage", default="pipeline5_eigenfunction_peaks_by_state_maxabs")
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--conditions", nargs="+", default=list(DEFAULT_CONDITIONS))
    parser.add_argument("--methods", nargs="+", default=list(DEFAULT_METHODS))
    parser.add_argument("--component-counts", nargs="+", type=int, default=list(DEFAULT_COMPONENT_COUNTS))
    parser.add_argument("--event-keywords", nargs="+", default=list(DEFAULT_EVENT_KEYWORDS))
    parser.add_argument(
        "--family-mode",
        choices=("explicit", "substring"),
        default="explicit",
        help=(
            "explicit uses curated family membership where sharp-wave-ripple is allowed "
            "for theta/gamma as well as ripple; substring uses the older name-contains rule."
        ),
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-effect", type=float, default=0.0)
    parser.add_argument("--min-cohen-d", type=float, default=0.2)
    parser.add_argument(
        "--amplitude-column",
        default="mean_event_minus_baseline",
        choices=("mean_event_minus_baseline", "event_mean_peak", "event_median_peak"),
        help="Column used when judging whether amplitudes are similar.",
    )
    parser.add_argument("--equal-cv-threshold", type=float, default=0.20)
    parser.add_argument("--equal-relative-range-threshold", type=float, default=0.25)
    return parser.parse_args()


def safe_float(value: object, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_int(value: object, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: List[str] = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fieldnames.append(key)
                seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def method_tag(method: str, component_count: int) -> str:
    return f"{method.lower()}_k{component_count:02d}"


def discover_runs(
    processed_root: Path,
    peak_stage: str,
    datasets: Sequence[str],
    conditions: Sequence[str],
    methods: Sequence[str],
    component_counts: Sequence[int],
) -> Tuple[List[RunSpec], List[Dict[str, object]]]:
    runs: List[RunSpec] = []
    missing: List[Dict[str, object]] = []
    for dataset in datasets:
        for condition in conditions:
            for method in methods:
                for k in component_counts:
                    tag = method_tag(method, k)
                    root = processed_root / dataset / peak_stage / condition / tag
                    stats_file = root / f"{dataset}_{condition}_{tag}_peaks_stats.csv"
                    if stats_file.is_file():
                        runs.append(RunSpec(dataset, condition, method.lower(), k, tag, stats_file))
                    else:
                        missing.append(
                            {
                                "dataset": dataset,
                                "condition": condition,
                                "method": method.lower(),
                                "component_count": k,
                                "method_tag": tag,
                                "missing_stats_file": str(stats_file),
                            }
                        )
    return runs, missing


def event_matches(label: str, keyword: str) -> bool:
    return keyword.lower() in label.lower()


def family_labels_for(
    state_labels: Sequence[str],
    family: str,
    family_mode: str,
) -> List[str]:
    family_key = family.lower()
    if family_mode == "substring":
        return [label for label in state_labels if event_matches(label, family_key)]

    explicit = DEFAULT_EXPLICIT_FAMILY_EVENTS.get(family_key)
    if explicit is None:
        return [label for label in state_labels if event_matches(label, family_key)]
    explicit_set = {label.lower() for label in explicit}
    return [label for label in state_labels if label.lower() in explicit_set]


def is_active(row: Dict[str, str], alpha: float, min_effect: float, min_cohen_d: float) -> bool:
    q_value = safe_float(row.get("q_vs_baseline_paired_ttest_two_sided"))
    effect = safe_float(row.get("mean_event_minus_baseline"))
    cohen_d = safe_float(row.get("cohen_d_paired_vs_baseline"))
    return (
        math.isfinite(q_value)
        and q_value <= alpha
        and math.isfinite(effect)
        and effect > min_effect
        and math.isfinite(cohen_d)
        and cohen_d >= min_cohen_d
    )


def summarize_values(values: Sequence[float]) -> Dict[str, float]:
    clean = [value for value in values if math.isfinite(value)]
    if not clean:
        return {
            "amplitude_mean": math.nan,
            "amplitude_min": math.nan,
            "amplitude_max": math.nan,
            "amplitude_range": math.nan,
            "amplitude_cv": math.nan,
            "amplitude_relative_range": math.nan,
        }
    avg = mean(clean)
    amp_min = min(clean)
    amp_max = max(clean)
    amp_range = amp_max - amp_min
    denom = max(abs(avg), 1e-12)
    cv = pstdev(clean) / denom if len(clean) > 1 else 0.0
    rel_range = amp_range / denom
    return {
        "amplitude_mean": avg,
        "amplitude_min": amp_min,
        "amplitude_max": amp_max,
        "amplitude_range": amp_range,
        "amplitude_cv": cv,
        "amplitude_relative_range": rel_range,
    }


def joined(items: Iterable[object]) -> str:
    return ";".join(str(item) for item in items)


def component_rows_for_run(
    run: RunSpec,
    stats_rows: Sequence[Dict[str, str]],
    event_keywords: Sequence[str],
    family_mode: str,
    alpha: float,
    min_effect: float,
    min_cohen_d: float,
    amplitude_column: str,
    equal_cv_threshold: float,
    equal_relative_range_threshold: float,
) -> List[Dict[str, object]]:
    by_component: Dict[int, List[Dict[str, str]]] = defaultdict(list)
    for row in stats_rows:
        by_component[safe_int(row.get("component_idx"))].append(row)

    state_labels = sorted(
        {str(row.get("state_label", "")) for row in stats_rows},
        key=lambda label: min(
            safe_int(row.get("state_code"))
            for row in stats_rows
            if str(row.get("state_label", "")) == label
        ),
    )

    out: List[Dict[str, object]] = []
    for keyword in event_keywords:
        family_labels = family_labels_for(state_labels, keyword, family_mode)
        if not family_labels:
            continue
        family_set = set(family_labels)
        forbidden_labels = [label for label in state_labels if label not in family_set]
        forbidden_set = set(forbidden_labels)
        for comp_idx, rows in sorted(by_component.items()):
            family_rows = [row for row in rows if str(row.get("state_label", "")) in family_set]
            if not family_rows:
                continue
            forbidden_rows = [row for row in rows if str(row.get("state_label", "")) in forbidden_set]

            active_rows = [
                row
                for row in family_rows
                if is_active(row, alpha=alpha, min_effect=min_effect, min_cohen_d=min_cohen_d)
            ]
            active_labels = [str(row.get("state_label", "")) for row in active_rows]
            forbidden_active_rows = [
                row
                for row in forbidden_rows
                if is_active(row, alpha=alpha, min_effect=min_effect, min_cohen_d=min_cohen_d)
            ]
            forbidden_active_labels = [str(row.get("state_label", "")) for row in forbidden_active_rows]
            missing_active = [label for label in family_labels if label not in set(active_labels)]
            amplitudes = [safe_float(row.get(amplitude_column)) for row in family_rows]
            value_summary = summarize_values(amplitudes)
            all_active = len(missing_active) == 0
            selective = len(forbidden_active_labels) == 0
            similar = (
                all_active
                and value_summary["amplitude_cv"] <= equal_cv_threshold
                and value_summary["amplitude_relative_range"] <= equal_relative_range_threshold
            )
            selective_similar = similar and selective
            if not all_active:
                profile = "partial_or_inactive"
            elif not selective:
                profile = "nonselective"
            elif selective_similar:
                profile = "selective_similar"
            else:
                profile = "selective_unequal"
            sorted_by_amp = sorted(
                family_rows,
                key=lambda row: safe_float(row.get(amplitude_column), -math.inf),
            )
            weak_event = str(sorted_by_amp[0].get("state_label", "")) if sorted_by_amp else ""
            strong_event = str(sorted_by_amp[-1].get("state_label", "")) if sorted_by_amp else ""
            out.append(
                {
                    "dataset": run.dataset,
                    "condition": run.condition,
                    "method": run.method,
                    "component_count": run.component_count,
                    "method_tag": run.method_tag,
                    "event_family": keyword.lower(),
                    "family_mode": family_mode,
                    "component_idx": comp_idx,
                    "family_event_count": len(family_labels),
                    "family_events": joined(family_labels),
                    "forbidden_event_count": len(forbidden_labels),
                    "forbidden_events": joined(forbidden_labels),
                    "active_event_count": len(set(active_labels)),
                    "active_events": joined(label for label in family_labels if label in set(active_labels)),
                    "inactive_or_weak_events": joined(missing_active),
                    "forbidden_active_event_count": len(set(forbidden_active_labels)),
                    "forbidden_active_events": joined(label for label in forbidden_labels if label in set(forbidden_active_labels)),
                    "all_family_events_active": int(all_active),
                    "selective_against_forbidden_events": int(selective),
                    "all_active_and_selective": int(all_active and selective),
                    "amplitude_similar_by_tolerance": int(similar),
                    "all_active_selective_and_similar": int(selective_similar),
                    "amplitude_profile": profile,
                    "amplitude_column": amplitude_column,
                    "weakest_event_by_amplitude": weak_event,
                    "strongest_event_by_amplitude": strong_event,
                    **value_summary,
                    "min_q_vs_baseline": min(
                        safe_float(row.get("q_vs_baseline_paired_ttest_two_sided"))
                        for row in family_rows
                    ),
                    "max_q_vs_baseline": max(
                        safe_float(row.get("q_vs_baseline_paired_ttest_two_sided"))
                        for row in family_rows
                    ),
                    "min_cohen_d": min(safe_float(row.get("cohen_d_paired_vs_baseline")) for row in family_rows),
                    "max_cohen_d": max(safe_float(row.get("cohen_d_paired_vs_baseline")) for row in family_rows),
                    "stats_file": str(run.stats_file),
                }
            )
    return out


def best_rows(component_rows: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str, str, int, str], List[Dict[str, object]]] = defaultdict(list)
    for row in component_rows:
        key = (
            str(row["dataset"]),
            str(row["condition"]),
            str(row["method"]),
            safe_int(row["component_count"]),
            str(row["event_family"]),
        )
        grouped[key].append(row)

    best: List[Dict[str, object]] = []
    for _key, rows in grouped.items():
        rows_sorted = sorted(
            rows,
            key=lambda row: (
                -safe_int(row["all_family_events_active"]),
                -safe_int(row["selective_against_forbidden_events"]),
                -safe_int(row["all_active_and_selective"]),
                -safe_int(row["active_event_count"]),
                -safe_int(row["all_active_selective_and_similar"]),
                safe_float(row["amplitude_cv"], math.inf),
                safe_float(row["amplitude_relative_range"], math.inf),
                -safe_float(row["amplitude_mean"], -math.inf),
            ),
        )
        best.append(dict(rows_sorted[0]))
    return best


def aggregate_rows(best: Sequence[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    for row in best:
        grouped[(str(row["event_family"]), str(row["method"]))].append(row)

    out: List[Dict[str, object]] = []
    for (family, method), rows in sorted(grouped.items()):
        n = len(rows)
        all_active = sum(safe_int(row["all_family_events_active"]) for row in rows)
        selective = sum(safe_int(row["selective_against_forbidden_events"]) for row in rows)
        all_active_selective = sum(safe_int(row["all_active_and_selective"]) for row in rows)
        similar = sum(safe_int(row["amplitude_similar_by_tolerance"]) for row in rows)
        selective_similar = sum(safe_int(row["all_active_selective_and_similar"]) for row in rows)
        forbidden_any = sum(1 for row in rows if safe_int(row["forbidden_active_event_count"]) > 0)
        out.append(
            {
                "event_family": family,
                "method": method,
                "run_family_count": n,
                "best_component_all_active_count": all_active,
                "best_component_all_active_fraction": all_active / n if n else math.nan,
                "best_component_selective_count": selective,
                "best_component_selective_fraction": selective / n if n else math.nan,
                "best_component_all_active_selective_count": all_active_selective,
                "best_component_all_active_selective_fraction": all_active_selective / n if n else math.nan,
                "best_component_forbidden_active_count": forbidden_any,
                "best_component_forbidden_active_fraction": forbidden_any / n if n else math.nan,
                "best_component_similar_amplitude_count": similar,
                "best_component_similar_amplitude_fraction": similar / n if n else math.nan,
                "best_component_selective_similar_amplitude_count": selective_similar,
                "best_component_selective_similar_amplitude_fraction": selective_similar / n if n else math.nan,
                "median_best_cv": median(safe_float(row["amplitude_cv"]) for row in rows),
                "median_best_relative_range": median(
                    safe_float(row["amplitude_relative_range"]) for row in rows
                ),
            }
        )
    return out


def median(values: Iterable[float]) -> float:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return math.nan
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return (clean[mid - 1] + clean[mid]) / 2.0


def write_markdown(
    path: Path,
    args: argparse.Namespace,
    run_count: int,
    missing_count: int,
    best: Sequence[Dict[str, object]],
    aggregate: Sequence[Dict[str, object]],
) -> None:
    lines: List[str] = []
    lines.append("# Peak Event-Family Component Activity")
    lines.append("")
    lines.append(f"- Generated at: `{datetime.now().isoformat(timespec='seconds')}`")
    lines.append(f"- Datasets: `{', '.join(args.datasets)}`")
    lines.append(f"- Conditions: `{', '.join(args.conditions)}`")
    lines.append(f"- Methods: `{', '.join(args.methods)}`")
    lines.append(f"- Component counts: `{', '.join(str(k) for k in args.component_counts)}`")
    lines.append(f"- Event families: `{', '.join(args.event_keywords)}`")
    lines.append(f"- Family mode: `{args.family_mode}`")
    if args.family_mode == "explicit":
        explicit_text = "; ".join(
            f"{family}: {', '.join(DEFAULT_EXPLICIT_FAMILY_EVENTS[family])}"
            for family in args.event_keywords
            if family in DEFAULT_EXPLICIT_FAMILY_EVENTS
        )
        lines.append(f"- Explicit family membership: `{explicit_text}`")
    lines.append(f"- Runs found: `{run_count}`; missing stats files: `{missing_count}`")
    lines.append(
        "- Active definition: "
        f"`q <= {args.alpha}`, `mean_event_minus_baseline > {args.min_effect}`, "
        f"`cohen_d >= {args.min_cohen_d}`"
    )
    lines.append(
        "- Similar-amplitude definition: "
        f"`CV <= {args.equal_cv_threshold}` and "
        f"`relative range <= {args.equal_relative_range_threshold}` using "
        f"`{args.amplitude_column}`."
    )
    lines.append("")
    lines.append("## Method Summary")
    lines.append("")
    lines.append(
        "| Event family | Method | Runs | All-active best | All-active fraction | "
        "All-active selective | Selective fraction | Selective similar | Similar fraction | "
        "Forbidden-active fraction | Median CV | Median relative range |"
    )
    lines.append("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for row in aggregate:
        lines.append(
            "| {event_family} | {method} | {run_family_count} | "
            "{best_component_all_active_count} | {best_component_all_active_fraction:.3f} | "
            "{best_component_all_active_selective_count} | "
            "{best_component_all_active_selective_fraction:.3f} | "
            "{best_component_selective_similar_amplitude_count} | "
            "{best_component_selective_similar_amplitude_fraction:.3f} | "
            "{best_component_forbidden_active_fraction:.3f} | "
            "{median_best_cv:.3f} | {median_best_relative_range:.3f} |".format(**row)
        )
    lines.append("")
    lines.append("## Outputs")
    lines.append("")
    lines.append("- `component_family_activity.csv`: every component candidate.")
    lines.append("- `best_component_by_run_family.csv`: best component per dataset/condition/method/k/family.")
    lines.append("- `method_family_summary.csv`: aggregate fractions by event family and method.")
    lines.append("- `missing_stats_files.csv`: missing inputs, if any.")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    runs, missing = discover_runs(
        args.processed_root.resolve(),
        args.peak_stage,
        args.datasets,
        args.conditions,
        args.methods,
        args.component_counts,
    )
    component_rows: List[Dict[str, object]] = []
    for run in runs:
        stats_rows = read_csv(run.stats_file)
        component_rows.extend(
            component_rows_for_run(
                run,
                stats_rows,
                args.event_keywords,
                args.family_mode,
                args.alpha,
                args.min_effect,
                args.min_cohen_d,
                args.amplitude_column,
                args.equal_cv_threshold,
                args.equal_relative_range_threshold,
            )
        )

    best = best_rows(component_rows)
    aggregate = aggregate_rows(best)

    write_csv(output_dir / "component_family_activity.csv", component_rows)
    write_csv(output_dir / "best_component_by_run_family.csv", best)
    write_csv(output_dir / "method_family_summary.csv", aggregate)
    write_csv(output_dir / "missing_stats_files.csv", missing)
    write_markdown(output_dir / "summary.md", args, len(runs), len(missing), best, aggregate)

    print(f"Runs found          : {len(runs)}")
    print(f"Missing stats files : {len(missing)}")
    print(f"Component rows      : {len(component_rows)}")
    print(f"Output dir          : {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
