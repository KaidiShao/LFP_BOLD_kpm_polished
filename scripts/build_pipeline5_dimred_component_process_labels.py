#!/usr/bin/env python3
"""Build P5 dimred component process labels from event-family activity rows.

This is the P11 bridge between P5 selectivity and P8/P10 interpretation.  The
current default input is the existing maxabs selectivity table, so the output is
marked as transitional.  Once activity-magnitude peak statistics exist, point
``--input`` at that table and change ``--activity-transform`` accordingly.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean
from typing import Iterable, Sequence


DEFAULT_INPUT = (
    Path("results")
    / "peak_event_family_component_activity_current_maxabs"
    / "component_family_activity.csv"
)
DEFAULT_OUTPUT_DIR = Path("results") / "pipeline5_dimred_component_process_labels_current"
DEFAULT_DATASETS = (
    "e10aw1",
    "e10bv1",
    "e10fV1",
    "e10gb1",
    "e10gh1",
    "e10gw1",
    "f12m01",
    "k13m17",
    "k13m23",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument(
        "--expected-conditions",
        nargs="+",
        default=None,
        help="Optional condition grid used only for missing-grid reporting. Defaults to conditions observed in the input.",
    )
    parser.add_argument("--activity-transform", default="maxabs_transitional")
    parser.add_argument("--activity-window-policy", default="not_available_maxabs_transitional")
    parser.add_argument("--threshold-ratio", default="070")
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
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


def truthy(value: object) -> bool:
    text = str(value or "").strip().lower()
    return text in {"1", "true", "yes", "y"}


def as_float(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def condition_short(condition: str) -> str:
    text = condition.lower()
    if "complex" in text or "csplit" in text:
        return "csplit"
    if "abs" in text:
        return "abs"
    return text


def density_name(condition: str, method: str, k: int, threshold_ratio: str) -> str:
    return f"dim_{condition_short(condition)}_{method.lower()}{k}_q{threshold_ratio}"


def joined(items: Iterable[object]) -> str:
    return ";".join(str(item) for item in items if str(item) != "")


def family_state(row: dict[str, str]) -> dict[str, object]:
    active = truthy(row.get("all_family_events_active"))
    selective = truthy(row.get("selective_against_forbidden_events")) and active
    similar = truthy(row.get("amplitude_similar_by_tolerance")) and selective
    return {
        "active": active,
        "selective": selective,
        "similar": similar,
        "amplitude_cv": as_float(row.get("amplitude_cv")),
        "relative_range": as_float(row.get("amplitude_relative_range")),
        "amplitude_mean": as_float(row.get("amplitude_mean")),
        "forbidden_active_events": row.get("forbidden_active_events", ""),
        "active_events": row.get("active_events", ""),
        "family_events": row.get("family_events", ""),
        "source_peak_stats_file": row.get("stats_file", ""),
    }


def choose_label(states: dict[str, dict[str, object]]) -> tuple[str, str, str]:
    theta = states.get("theta", {})
    ripple = states.get("ripple", {})
    gamma = states.get("gamma", {})
    theta_sel = bool(theta.get("selective"))
    ripple_sel = bool(ripple.get("selective"))
    gamma_sel = bool(gamma.get("selective"))
    theta_sim = bool(theta.get("similar"))
    ripple_sim = bool(ripple.get("similar"))
    gamma_active = bool(gamma.get("active"))
    active_families = [fam for fam, state in states.items() if state.get("active")]
    selective_families = [fam for fam, state in states.items() if state.get("selective")]

    if theta_sel and ripple_sel:
        if theta_sim and ripple_sim:
            return "theta_ripple_joint", "high", "theta and ripple both selective with similar amplitudes"
        return "mixed_theta_ripple", "medium", "theta and ripple both selective but amplitude similarity is incomplete"
    if theta_sel:
        if gamma_active:
            return "mixed_theta_gamma", "medium", "theta selective with gamma-family activity"
        return (
            "theta_selective_similar" if theta_sim else "theta_selective_unequal",
            "high" if theta_sim else "medium",
            "theta selective",
        )
    if ripple_sel:
        if gamma_active:
            return "mixed_ripple_gamma", "medium", "ripple selective with gamma-family activity"
        return (
            "ripple_selective_similar" if ripple_sim else "ripple_selective_unequal",
            "high" if ripple_sim else "medium",
            "ripple selective",
        )
    if gamma_sel:
        return "gamma_selective", "low", "gamma is optional/secondary for current mainline"
    if len(active_families) >= 3:
        return "pan_event", "low", "active across all tracked families but not selective"
    if active_families:
        return "partial_or_inactive", "low", "some event-family activity but no clean theta/ripple selectivity"
    if selective_families:
        return "nonselective", "low", "selectivity fields present but family activity incomplete"
    return "unlabeled", "low", "no event-family activity detected by current criteria"


def summarize_component(rows: Sequence[dict[str, str]], activity_transform: str, window_policy: str, threshold_ratio: str) -> dict[str, object]:
    first = rows[0]
    states = {str(row.get("event_family", "")).lower(): family_state(row) for row in rows}
    label, confidence, note = choose_label(states)
    dataset = str(first.get("dataset", ""))
    condition = str(first.get("condition", ""))
    method = str(first.get("method", "")).lower()
    k = as_int(first.get("component_count"))
    component_idx = as_int(first.get("component_idx"))
    active_set = [fam for fam in ("theta", "ripple", "gamma") if states.get(fam, {}).get("active")]
    selective_set = [fam for fam in ("theta", "ripple", "gamma") if states.get(fam, {}).get("selective")]
    cvs = [float(states[fam]["amplitude_cv"]) for fam in states if math.isfinite(float(states[fam].get("amplitude_cv", math.nan)))]
    rel_ranges = [
        float(states[fam]["relative_range"])
        for fam in states
        if math.isfinite(float(states[fam].get("relative_range", math.nan)))
    ]
    source_files = sorted({str(state.get("source_peak_stats_file", "")) for state in states.values() if state.get("source_peak_stats_file")})
    forbidden_events = sorted({str(state.get("forbidden_active_events", "")) for state in states.values() if state.get("forbidden_active_events")})
    density = density_name(condition, method, k, threshold_ratio)
    return {
        "dataset": dataset,
        "condition": condition,
        "condition_short": condition_short(condition),
        "method": method,
        "component_count": k,
        "k": k,
        "method_tag": str(first.get("method_tag", f"{method}_k{k:02d}")),
        "component_idx": component_idx,
        "density_name_for_p8_p10": density,
        "density_name": density,
        "density_index": component_idx,
        "theta_active": int(states.get("theta", {}).get("active", False)),
        "theta_selective": int(states.get("theta", {}).get("selective", False)),
        "theta_similar": int(states.get("theta", {}).get("similar", False)),
        "ripple_active": int(states.get("ripple", {}).get("active", False)),
        "ripple_selective": int(states.get("ripple", {}).get("selective", False)),
        "ripple_similar": int(states.get("ripple", {}).get("similar", False)),
        "gamma_active": int(states.get("gamma", {}).get("active", False)),
        "gamma_selective": int(states.get("gamma", {}).get("selective", False)),
        "gamma_similar": int(states.get("gamma", {}).get("similar", False)),
        "active_family_set": joined(active_set),
        "selective_family_set": joined(selective_set),
        "primary_process_label": label,
        "label_confidence": confidence,
        "label_note": note,
        "amplitude_cv_by_family": joined(f"{fam}:{states[fam]['amplitude_cv']:.6g}" for fam in sorted(states) if math.isfinite(float(states[fam].get("amplitude_cv", math.nan)))),
        "amplitude_relative_range_by_family": joined(f"{fam}:{states[fam]['relative_range']:.6g}" for fam in sorted(states) if math.isfinite(float(states[fam].get("relative_range", math.nan)))),
        "amplitude_cv_mean": f"{mean(cvs):.10g}" if cvs else "",
        "amplitude_relative_range_mean": f"{mean(rel_ranges):.10g}" if rel_ranges else "",
        "forbidden_active_events_by_family": joined(forbidden_events),
        "source_peak_stats_file": joined(source_files),
        "lfp_activity_transform": activity_transform,
        "lfp_activity_window_policy": window_policy,
        "label_source_status": "transitional_maxabs" if "maxabs" in activity_transform else "activity_magnitude",
    }


def build_labels(
    rows: Sequence[dict[str, str]],
    datasets: Sequence[str],
    activity_transform: str,
    window_policy: str,
    threshold_ratio: str,
    expected_conditions: Sequence[str] | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    dataset_set = {d.lower() for d in datasets}
    grouped: dict[tuple[str, str, str, int, int], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        dataset = str(row.get("dataset", ""))
        if dataset.lower() not in dataset_set:
            continue
        key = (
            dataset,
            str(row.get("condition", "")),
            str(row.get("method", "")).lower(),
            as_int(row.get("component_count")),
            as_int(row.get("component_idx")),
        )
        grouped[key].append(row)

    labels: list[dict[str, object]] = []
    for _key, group_rows in sorted(grouped.items()):
        labels.append(summarize_component(group_rows, activity_transform, window_policy, threshold_ratio))

    present = {(str(row["dataset"]).lower(), str(row["condition"]), str(row["method"]), int(row["component_count"])) for row in labels}
    if expected_conditions is None:
        expected_conditions = sorted(
            {
                str(row.get("condition", ""))
                for row in rows
                if str(row.get("condition", "")).strip()
            }
        )
    expected_methods = ("svd", "nmf", "mds", "umap")
    expected_ks = tuple(range(3, 9))
    missing: list[dict[str, object]] = []
    for dataset in datasets:
        for condition in expected_conditions:
            for method in expected_methods:
                for k in expected_ks:
                    if (dataset.lower(), condition, method, k) not in present:
                        missing.append(
                            {
                                "dataset": dataset,
                                "condition": condition,
                                "method": method,
                                "component_count": k,
                                "method_tag": f"{method}_k{k:02d}",
                                "missing_reason": "no component-family rows for this dataset/condition/method/k",
                            }
                        )
    return labels, missing


def main() -> int:
    args = parse_args()
    rows = read_csv(args.input)
    labels, missing = build_labels(
        rows,
        args.datasets,
        args.activity_transform,
        args.activity_window_policy,
        args.threshold_ratio,
        args.expected_conditions,
    )
    write_csv(args.output_dir / "dimred_efun_process_labels.csv", labels)
    write_csv(args.output_dir / "missing_dimred_efun_process_labels.csv", missing)
    summary = [
        "# P5 Dimred Component Process Labels",
        "",
        f"- Input: `{args.input}`",
        f"- Output rows: `{len(labels)}`",
        f"- Missing run grids: `{len(missing)}`",
        f"- Activity transform: `{args.activity_transform}`",
        f"- Window policy: `{args.activity_window_policy}`",
        "",
        "Current maxabs-derived labels are transitional until the activity-magnitude P5 branch is regenerated.",
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"Input rows : {len(rows)}")
    print(f"Label rows : {len(labels)}")
    print(f"Missing    : {len(missing)}")
    print(f"Output dir : {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
