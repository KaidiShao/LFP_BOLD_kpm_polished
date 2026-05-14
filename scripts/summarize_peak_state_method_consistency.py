#!/usr/bin/env python3
"""Summarize cross-method consistency for pipeline5 event-to-component matches.

This script scans `*_peaks_stats.csv` files under one pipeline5 root, then:
1. Builds an event x component score table for each run.
2. Solves a one-to-one maximum-score assignment with unmatched events allowed.
3. Writes assignment and consistency summaries to a results directory.

Default score:
    cohen_d_paired_vs_baseline

Default validity filter:
    q_vs_baseline_paired_ttest_two_sided <= alpha
    cohen_d_paired_vs_baseline > 0
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


DEFAULT_ROOT = Path(r"E:\DataPons_processed\e10gb1\pipeline5_eigenfunction_peaks_by_state")
DEFAULT_OUTPUT_DIR = Path("results") / "peak_state_method_consistency_e10gb1"
DEFAULT_METHODS = ("svd", "nmf", "mds", "umap", "logsvd")
DEFAULT_ALPHA = 0.05
DEFAULT_SCORE_COLUMN = "cohen_d_paired_vs_baseline"
DEFAULT_Q_COLUMN = "q_vs_baseline_paired_ttest_two_sided"


@dataclass(frozen=True)
class PairStat:
    state_code: int
    state_label: str
    component_idx: int
    score_value: float
    q_value: float
    mean_event_minus_baseline: float
    significant_vs_baseline: int
    n_paired_windows: int


@dataclass(frozen=True)
class RunSpec:
    variant: str
    method: str
    run_id: str
    stats_file: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="Root folder containing variant/method/*_peaks_stats.csv files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for summary outputs.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_ALPHA,
        help="FDR threshold for valid event-component pairs.",
    )
    parser.add_argument(
        "--methods",
        nargs="*",
        default=list(DEFAULT_METHODS),
        help="Methods to include. Default: %(default)s",
    )
    parser.add_argument(
        "--include-shared",
        action="store_true",
        help="Include the 'shared' subfolder when present.",
    )
    parser.add_argument(
        "--score-column",
        default=DEFAULT_SCORE_COLUMN,
        help="Column used as the assignment objective.",
    )
    parser.add_argument(
        "--q-column",
        default=DEFAULT_Q_COLUMN,
        help="Adjusted-p column used in the significance filter.",
    )
    parser.add_argument(
        "--require-positive-mean-effect",
        action="store_true",
        help="Also require mean_event_minus_baseline > 0.",
    )
    parser.add_argument(
        "--require-significant-flag",
        action="store_true",
        help="Also require significant_vs_baseline == 1.",
    )
    return parser.parse_args()


def safe_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def safe_float(value: str, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def discover_runs(
    root: Path,
    methods: Sequence[str],
    include_shared: bool,
) -> List[RunSpec]:
    allowed_methods = set(methods)
    if include_shared:
        allowed_methods.add("shared")

    runs: List[RunSpec] = []
    for variant_dir in sorted(root.iterdir()):
        if not variant_dir.is_dir():
            continue
        variant = variant_dir.name
        for method_dir in sorted(variant_dir.iterdir()):
            if not method_dir.is_dir():
                continue
            method = method_dir.name
            if method not in allowed_methods:
                continue
            stats_files = sorted(method_dir.glob("*_peaks_stats.csv"))
            if not stats_files:
                continue
            stats_file = stats_files[0]
            run_id = f"{variant}/{method}"
            runs.append(RunSpec(variant, method, run_id, stats_file))
    return runs


def load_pair_stats(
    path: Path,
    score_column: str,
    q_column: str,
) -> List[PairStat]:
    rows: List[PairStat] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                PairStat(
                    state_code=safe_int(row.get("state_code", "")),
                    state_label=str(row.get("state_label", "")).strip(),
                    component_idx=safe_int(row.get("component_idx", "")),
                    score_value=safe_float(row.get(score_column, "")),
                    q_value=safe_float(row.get(q_column, "")),
                    mean_event_minus_baseline=safe_float(
                        row.get("mean_event_minus_baseline", "")
                    ),
                    significant_vs_baseline=safe_int(
                        row.get("significant_vs_baseline", "")
                    ),
                    n_paired_windows=safe_int(row.get("n_paired_windows", "")),
                )
            )
    rows.sort(key=lambda item: (item.state_code, item.component_idx))
    return rows


def is_valid_pair(
    stat: PairStat,
    alpha: float,
    require_positive_mean_effect: bool,
    require_significant_flag: bool,
) -> bool:
    if not math.isfinite(stat.q_value) or stat.q_value > alpha:
        return False
    if not math.isfinite(stat.score_value) or stat.score_value <= 0:
        return False
    if require_positive_mean_effect:
        if not math.isfinite(stat.mean_event_minus_baseline):
            return False
        if stat.mean_event_minus_baseline <= 0:
            return False
    if require_significant_flag and stat.significant_vs_baseline != 1:
        return False
    return True


def summarize_run(
    run: RunSpec,
    alpha: float,
    score_column: str,
    q_column: str,
    require_positive_mean_effect: bool,
    require_significant_flag: bool,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    pair_stats = load_pair_stats(run.stats_file, score_column, q_column)
    if not pair_stats:
        return [], {}

    events: List[Tuple[int, str]] = sorted(
        {(stat.state_code, stat.state_label) for stat in pair_stats},
        key=lambda item: item[0],
    )
    components: List[int] = sorted({stat.component_idx for stat in pair_stats})
    stat_lookup: Dict[Tuple[int, int], PairStat] = {
        (stat.state_code, stat.component_idx): stat for stat in pair_stats
    }

    valid_options: Dict[int, List[Tuple[int, PairStat]]] = {code: [] for code, _ in events}
    best_by_event: Dict[int, Optional[PairStat]] = {}
    for event_code, _event_label in events:
        best_stat: Optional[PairStat] = None
        for component_idx in components:
            stat = stat_lookup.get((event_code, component_idx))
            if stat is None or not is_valid_pair(
                stat,
                alpha,
                require_positive_mean_effect=require_positive_mean_effect,
                require_significant_flag=require_significant_flag,
            ):
                continue
            valid_options[event_code].append((component_idx, stat))
            if best_stat is None:
                best_stat = stat
            elif stat.score_value > best_stat.score_value:
                best_stat = stat
            elif stat.score_value == best_stat.score_value and stat.q_value < best_stat.q_value:
                best_stat = stat
        best_by_event[event_code] = best_stat

    @lru_cache(maxsize=None)
    def solve(event_idx: int, used_components: Tuple[int, ...]) -> Tuple[float, int, float, Tuple[Optional[int], ...]]:
        if event_idx >= len(events):
            return 0.0, 0, 0.0, tuple()

        used = set(used_components)
        event_code, _event_label = events[event_idx]

        best = solve(event_idx + 1, tuple(sorted(used)))
        best_total, best_matched, best_qsum, best_assignments = best
        best_candidate = (best_total, best_matched, best_qsum, (None,) + best_assignments)

        for component_idx, stat in valid_options[event_code]:
            if component_idx in used:
                continue
            next_used = tuple(sorted((*used_components, component_idx)))
            next_total, next_matched, next_qsum, next_assignments = solve(event_idx + 1, next_used)
            candidate = (
                next_total + stat.score_value,
                next_matched + 1,
                next_qsum + stat.q_value,
                (component_idx,) + next_assignments,
            )
            if is_better_solution(candidate, best_candidate):
                best_candidate = candidate

        return best_candidate

    total_score, matched_count, q_sum, assignments = solve(0, tuple())
    assignment_rows: List[Dict[str, object]] = []
    unmatched_labels: List[str] = []
    matched_labels: List[str] = []

    for (event_code, event_label), selected_component in zip(events, assignments):
        selected_stat = (
            stat_lookup.get((event_code, selected_component))
            if selected_component is not None
            else None
        )
        best_stat = best_by_event[event_code]
        matched = int(selected_component is not None and selected_stat is not None)
        if matched:
            matched_labels.append(event_label)
        else:
            unmatched_labels.append(event_label)

        assignment_rows.append(
            {
                "variant": run.variant,
                "method": run.method,
                "run_id": run.run_id,
                "stats_file": str(run.stats_file),
                "state_code": event_code,
                "state_label": event_label,
                "selected_component_idx": "" if selected_component is None else selected_component,
                "matched": matched,
                "assignment_score": "" if selected_stat is None else selected_stat.score_value,
                "assignment_q_value": "" if selected_stat is None else selected_stat.q_value,
                "assignment_mean_event_minus_baseline": ""
                if selected_stat is None
                else selected_stat.mean_event_minus_baseline,
                "assignment_n_paired_windows": ""
                if selected_stat is None
                else selected_stat.n_paired_windows,
                "best_component_idx": "" if best_stat is None else best_stat.component_idx,
                "best_score": "" if best_stat is None else best_stat.score_value,
                "best_q_value": "" if best_stat is None else best_stat.q_value,
                "assignment_is_event_best": int(
                    best_stat is not None
                    and selected_component is not None
                    and best_stat.component_idx == selected_component
                ),
                "run_total_score": total_score,
                "run_matched_event_count": matched_count,
                "run_event_count": len(events),
                "run_component_count": len(components),
            }
        )

    run_summary = {
        "variant": run.variant,
        "method": run.method,
        "run_id": run.run_id,
        "stats_file": str(run.stats_file),
        "event_count": len(events),
        "component_count": len(components),
        "matched_event_count": matched_count,
        "unmatched_event_count": len(events) - matched_count,
        "matched_state_labels": "; ".join(matched_labels),
        "unmatched_state_labels": "; ".join(unmatched_labels),
        "total_assignment_score": total_score,
        "assignment_q_sum": q_sum,
        "alpha": alpha,
        "score_column": score_column,
        "q_column": q_column,
    }
    return assignment_rows, run_summary


def is_better_solution(
    candidate: Tuple[float, int, float, Tuple[Optional[int], ...]],
    incumbent: Tuple[float, int, float, Tuple[Optional[int], ...]],
) -> bool:
    cand_total, cand_matched, cand_qsum, cand_assignments = candidate
    inc_total, inc_matched, inc_qsum, inc_assignments = incumbent

    tol = 1e-12
    if cand_total > inc_total + tol:
        return True
    if inc_total > cand_total + tol:
        return False

    if cand_matched > inc_matched:
        return True
    if inc_matched > cand_matched:
        return False

    if cand_qsum + tol < inc_qsum:
        return True
    if inc_qsum + tol < cand_qsum:
        return False

    return stringify_assignments(cand_assignments) < stringify_assignments(inc_assignments)


def stringify_assignments(assignments: Sequence[Optional[int]]) -> str:
    parts = []
    for item in assignments:
        parts.append("NA" if item is None else f"{item:04d}")
    return "|".join(parts)


def mean_or_blank(values: Iterable[float]) -> object:
    values = [value for value in values if math.isfinite(value)]
    if not values:
        return ""
    return sum(values) / len(values)


def median_or_blank(values: Iterable[float]) -> object:
    values = [value for value in values if math.isfinite(value)]
    if not values:
        return ""
    return statistics.median(values)


def join_sorted(values: Iterable[str]) -> str:
    clean = sorted({value for value in values if value})
    return "; ".join(clean)


def build_event_summary(
    assignments: List[Dict[str, object]],
    group_keys: Sequence[str],
) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[object, ...], List[Dict[str, object]]] = {}
    for row in assignments:
        key = tuple(row[group_key] for group_key in group_keys)
        grouped.setdefault(key, []).append(row)

    summary_rows: List[Dict[str, object]] = []
    for key in sorted(grouped):
        rows = grouped[key]
        matched_rows = [row for row in rows if int(row["matched"]) == 1]
        scores = [safe_float(str(row["assignment_score"])) for row in matched_rows]
        base: Dict[str, object] = {group_key: key[idx] for idx, group_key in enumerate(group_keys)}
        base.update(
            {
                "run_count": len(rows),
                "matched_run_count": len(matched_rows),
                "unmatched_run_count": len(rows) - len(matched_rows),
                "match_rate": len(matched_rows) / len(rows) if rows else "",
                "mean_assignment_score_matched": mean_or_blank(scores),
                "median_assignment_score_matched": median_or_blank(scores),
                "matched_methods": join_sorted(str(row["method"]) for row in matched_rows),
                "matched_variants": join_sorted(str(row["variant"]) for row in matched_rows),
                "selected_components": join_sorted(
                    str(row["selected_component_idx"]) for row in matched_rows
                ),
                "source_runs": join_sorted(str(row["run_id"]) for row in rows),
            }
        )
        summary_rows.append(base)
    return summary_rows


def build_method_state_summary(assignments: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str, int], List[Dict[str, object]]] = {}
    for row in assignments:
        key = (str(row["method"]), str(row["state_label"]), int(row["state_code"]))
        grouped.setdefault(key, []).append(row)

    summary_rows: List[Dict[str, object]] = []
    for key in sorted(grouped):
        rows = grouped[key]
        matched_rows = [row for row in rows if int(row["matched"]) == 1]
        scores = [safe_float(str(row["assignment_score"])) for row in matched_rows]
        summary_rows.append(
            {
                "method": key[0],
                "state_label": key[1],
                "state_code": key[2],
                "variant_count": len(rows),
                "matched_variant_count": len(matched_rows),
                "match_rate": len(matched_rows) / len(rows) if rows else "",
                "mean_assignment_score_matched": mean_or_blank(scores),
                "median_assignment_score_matched": median_or_blank(scores),
                "matched_variants": join_sorted(str(row["variant"]) for row in matched_rows),
                "unmatched_variants": join_sorted(
                    str(row["variant"]) for row in rows if int(row["matched"]) == 0
                ),
            }
        )
    return summary_rows


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_summary_markdown(
    path: Path,
    metadata: Dict[str, object],
    run_summaries: List[Dict[str, object]],
    overall_event_summary: List[Dict[str, object]],
    variant_event_summary: List[Dict[str, object]],
    method_state_summary: List[Dict[str, object]],
) -> None:
    filter_terms = [f"`{metadata['q_column']} <= {metadata['alpha']}`", f"`{metadata['score_column']} > 0`"]
    if metadata.get("require_positive_mean_effect"):
        filter_terms.append("`mean_event_minus_baseline > 0`")
    if metadata.get("require_significant_flag"):
        filter_terms.append("`significant_vs_baseline == 1`")

    lines: List[str] = []
    lines.append("# Peak-State Method Consistency Summary")
    lines.append("")
    lines.append(f"- Generated at: `{metadata['generated_at']}`")
    lines.append(f"- Input root: `{metadata['input_root']}`")
    lines.append(f"- Included methods: `{', '.join(metadata['included_methods'])}`")
    lines.append(f"- Alpha: `{metadata['alpha']}`")
    lines.append(f"- Score column: `{metadata['score_column']}`")
    lines.append(f"- Valid pair filter: {', '.join(filter_terms)}")
    lines.append("")
    lines.append("## Run Summary")
    lines.append("")
    lines.append("| Run | Matched Events | Unmatched Events | Total Score | Unmatched States |")
    lines.append("| --- | ---: | ---: | ---: | --- |")
    for row in run_summaries:
        lines.append(
            "| {run_id} | {matched_event_count} | {unmatched_event_count} | {total_assignment_score:.6f} | {unmatched_state_labels} |".format(
                **row
            )
        )
    lines.append("")
    lines.append("## Overall Event Stability")
    lines.append("")
    lines.append("| Event | Matched Runs | Match Rate | Median Score | Matched Variants |")
    lines.append("| --- | ---: | ---: | ---: | --- |")
    for row in overall_event_summary:
        median_score = row["median_assignment_score_matched"]
        median_text = "" if median_score == "" else f"{median_score:.6f}"
        match_rate = row["match_rate"]
        match_rate_text = "" if match_rate == "" else f"{match_rate:.3f}"
        lines.append(
            f"| {row['state_label']} | {row['matched_run_count']}/{row['run_count']} | "
            f"{match_rate_text} | {median_text} | {row['matched_variants']} |"
        )
    lines.append("")
    lines.append("## Variant-by-Variant Event Stability")
    lines.append("")
    lines.append("| Variant | Event | Matched Methods | Match Rate | Median Score |")
    lines.append("| --- | --- | ---: | ---: | ---: |")
    for row in variant_event_summary:
        median_score = row["median_assignment_score_matched"]
        median_text = "" if median_score == "" else f"{median_score:.6f}"
        match_rate = row["match_rate"]
        match_rate_text = "" if match_rate == "" else f"{match_rate:.3f}"
        lines.append(
            f"| {row['variant']} | {row['state_label']} | {row['matched_run_count']}/{row['run_count']} | "
            f"{match_rate_text} | {median_text} |"
        )
    lines.append("")
    lines.append("## Method-Level Event Stability")
    lines.append("")
    lines.append("| Method | Event | Matched Variants | Match Rate | Median Score |")
    lines.append("| --- | --- | ---: | ---: | ---: |")
    for row in method_state_summary:
        median_score = row["median_assignment_score_matched"]
        median_text = "" if median_score == "" else f"{median_score:.6f}"
        match_rate = row["match_rate"]
        match_rate_text = "" if match_rate == "" else f"{match_rate:.3f}"
        lines.append(
            f"| {row['method']} | {row['state_label']} | {row['matched_variant_count']}/{row['variant_count']} | "
            f"{match_rate_text} | {median_text} |"
        )
    lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    runs = discover_runs(root, args.methods, args.include_shared)
    if not runs:
        raise SystemExit(f"No stats files found under {root}")

    all_assignments: List[Dict[str, object]] = []
    run_summaries: List[Dict[str, object]] = []
    for run in runs:
        assignment_rows, run_summary = summarize_run(
            run=run,
            alpha=args.alpha,
            score_column=args.score_column,
            q_column=args.q_column,
            require_positive_mean_effect=args.require_positive_mean_effect,
            require_significant_flag=args.require_significant_flag,
        )
        all_assignments.extend(assignment_rows)
        if run_summary:
            run_summaries.append(run_summary)

    overall_event_summary = build_event_summary(
        all_assignments,
        group_keys=("state_code", "state_label"),
    )
    variant_event_summary = build_event_summary(
        all_assignments,
        group_keys=("variant", "state_code", "state_label"),
    )
    method_state_summary = build_method_state_summary(all_assignments)

    metadata = {
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "input_root": str(root),
        "output_dir": str(output_dir),
        "alpha": args.alpha,
        "score_column": args.score_column,
        "q_column": args.q_column,
        "require_positive_mean_effect": args.require_positive_mean_effect,
        "require_significant_flag": args.require_significant_flag,
        "included_methods": list(args.methods) + (["shared"] if args.include_shared else []),
        "run_count": len(run_summaries),
    }

    write_csv(output_dir / "assignment_summary.csv", all_assignments)
    write_csv(output_dir / "run_summary.csv", run_summaries)
    write_csv(output_dir / "overall_event_consistency.csv", overall_event_summary)
    write_csv(output_dir / "variant_event_consistency.csv", variant_event_summary)
    write_csv(output_dir / "method_event_consistency.csv", method_state_summary)
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    write_summary_markdown(
        output_dir / "summary.md",
        metadata=metadata,
        run_summaries=run_summaries,
        overall_event_summary=overall_event_summary,
        variant_event_summary=variant_event_summary,
        method_state_summary=method_state_summary,
    )

    print(f"Wrote consistency summaries to: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
