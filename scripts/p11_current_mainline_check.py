#!/usr/bin/env python3
"""One-command P11 mainline check for P1-P10 current analysis.

This wrapper is intentionally audit/summary only. It refreshes repo-local P11
CSV products and the P11 audit plan, but it does not launch MATLAB backfills,
delete files, or move legacy outputs unless explicit flags are passed through
to the underlying audit.
"""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATASETS_CURRENT5 = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01")
DEFAULT_DATASETS_ALL9 = (
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


@dataclass(frozen=True)
class StepResult:
    name: str
    returncode: int
    command: str


P1_P10_PIPELINES = {f"P{i}" for i in range(1, 11)} | {f"PIPELINE{i}" for i in range(1, 11)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-scope", choices=("current5", "all9"), default="all9")
    parser.add_argument("--datasets", nargs="+", default=None)
    parser.add_argument(
        "--audit-output-dir",
        type=Path,
        default=None,
        help="Output dir for p11_current_analysis_audit.py. Defaults to a scope-specific results folder.",
    )
    parser.add_argument("--skip-scorecard-figures", action="store_true")
    parser.add_argument("--copy-flat-figures", action="store_true")
    parser.add_argument("--archive-existing-summary", action="store_true")
    parser.add_argument("--quarantine-legacy", action="store_true")
    parser.add_argument("--overwrite-flat", action="store_true")
    parser.add_argument("--max-flat-copies", type=int, default=None)
    return parser.parse_args()


def default_datasets(scope: str) -> list[str]:
    return list(DEFAULT_DATASETS_ALL9 if scope == "all9" else DEFAULT_DATASETS_CURRENT5)


def run_step(name: str, cmd: Sequence[str]) -> StepResult:
    print("")
    print(f"== {name} ==")
    print(" ".join(cmd))
    completed = subprocess.run(cmd, cwd=REPO_ROOT)
    return StepResult(name=name, returncode=completed.returncode, command=" ".join(cmd))


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


def p1_p10_rows(rows: Sequence[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if row.get("pipeline") in P1_P10_PIPELINES]


def write_public_outputs(audit_output_dir: Path, output_dir: Path) -> dict[str, Path]:
    status_rows = read_csv(audit_output_dir / "p11_stage_status_latest.csv")
    fill_rows = read_csv(audit_output_dir / "p11_fill_plan_latest.csv")
    flat_rows = read_csv(audit_output_dir / "p11_flat_figure_sources_latest.csv")
    legacy_rows = read_csv(audit_output_dir / "p11_legacy_candidates_latest.csv")

    current_status = p1_p10_rows(status_rows)
    recompute_plan = p1_p10_rows(fill_rows)
    flat_manifest = p1_p10_rows(flat_rows)
    legacy_candidates = p1_p10_rows(legacy_rows)

    outputs = {
        "current_status": output_dir / "P1-P10_current_status.csv",
        "recompute_plan": output_dir / "P1-P10_recompute_plan.csv",
        "flat_figure_manifest": output_dir / "P1-P10_flat_figure_manifest.csv",
        "legacy_candidates": output_dir / "P1-P10_legacy_candidates.csv",
    }
    write_csv(outputs["current_status"], current_status)
    write_csv(outputs["recompute_plan"], recompute_plan)
    write_csv(outputs["flat_figure_manifest"], flat_manifest)
    write_csv(outputs["legacy_candidates"], legacy_candidates)
    return outputs


def summarize_status_counts(status_path: Path) -> list[str]:
    rows = read_csv(status_path)
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        counts[row.get("pipeline", "")][row.get("status", "")] += 1
    lines: list[str] = []
    for pipeline in sorted(counts):
        if pipeline not in P1_P10_PIPELINES:
            continue
        summary = ", ".join(f"{status}={count}" for status, count in sorted(counts[pipeline].items()))
        lines.append(f"- `{pipeline}`: {summary}")
    return lines


def summarize_recompute_counts(plan_path: Path) -> list[str]:
    rows = read_csv(plan_path)
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        pipeline = row.get("pipeline", "")
        if pipeline not in P1_P10_PIPELINES:
            continue
        counts[pipeline][row.get("action", "")] += 1
    lines: list[str] = []
    for pipeline in sorted(counts):
        actions = "; ".join(f"{action}={count}" for action, count in counts[pipeline].most_common())
        lines.append(f"- `{pipeline}`: {actions}")
    return lines


def write_summary(
    path: Path,
    datasets: Sequence[str],
    results: Sequence[StepResult],
    public_outputs: dict[str, Path],
) -> None:
    lines = [
        "# P11 Current Mainline Check",
        "",
        f"- Timestamp: `{datetime.now().isoformat(timespec='seconds')}`",
        f"- Datasets: `{', '.join(datasets)}`",
        "",
        "## Steps",
        "",
    ]
    for result in results:
        status = "ok" if result.returncode == 0 else f"failed:{result.returncode}"
        lines.append(f"- `{result.name}`: `{status}`")
        lines.append(f"  - `{result.command}`")
    lines.extend(["", "## Public P1-P10 Outputs", ""])
    for name, output_path in public_outputs.items():
        lines.append(f"- `{name}`: `{output_path}`")
    lines.extend(["", "## P1-P10 Status Counts", ""])
    lines.extend(summarize_status_counts(public_outputs["current_status"]) or ["- No P1-P10 status rows found."])
    lines.extend(["", "## P1-P10 Recompute Actions", ""])
    lines.extend(summarize_recompute_counts(public_outputs["recompute_plan"]) or ["- No P1-P10 recompute actions found."])
    lines.extend(
        [
            "",
            "## Main Outputs",
            "",
            "- `results/pipeline11_p1_p4_readiness_current/`",
            "- `results/pipeline5_dimred_component_process_labels_current/`",
            "- `results/pipeline11_parameter_selection_scorecard_current/`",
            "- `results/pipeline11_current_analysis_audit_all9/` or the selected audit output dir",
            "",
            "This wrapper does not run missing MATLAB/P4/P5-P10 backfills. Use the generated fill plan before launching heavy jobs.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    datasets = args.datasets or default_datasets(args.dataset_scope)
    audit_output_dir = args.audit_output_dir
    if audit_output_dir is None:
        audit_output_dir = Path("results") / f"pipeline11_current_analysis_audit_{args.dataset_scope}"

    py = sys.executable
    results: list[StepResult] = []

    results.append(
        run_step(
            "p1_p4_readiness",
            [py, "scripts/audit_p1_p4_readiness.py", "--datasets", *datasets],
        )
    )
    results.append(
        run_step(
            "p5_dimred_component_labels",
            [py, "scripts/build_pipeline5_dimred_component_process_labels.py", "--datasets", *datasets],
        )
    )
    scorecard_cmd = [
        py,
        "scripts/summarize_p11_parameter_selection_scorecard.py",
        "--datasets",
        *datasets,
    ]
    if args.skip_scorecard_figures:
        scorecard_cmd.append("--skip-figures")
    results.append(run_step("p11_parameter_selection_scorecard", scorecard_cmd))

    audit_cmd = [
        py,
        "scripts/p11_current_analysis_audit.py",
        "--dataset-scope",
        args.dataset_scope,
        "--datasets",
        *datasets,
        "--output-dir",
        str(audit_output_dir),
    ]
    if args.copy_flat_figures:
        audit_cmd.append("--copy-flat-figures")
    if args.archive_existing_summary:
        audit_cmd.append("--archive-existing-summary")
    if args.quarantine_legacy:
        audit_cmd.append("--quarantine-legacy")
    if args.overwrite_flat:
        audit_cmd.append("--overwrite-flat")
    if args.max_flat_copies is not None:
        audit_cmd.extend(["--max-flat-copies", str(args.max_flat_copies)])
    results.append(run_step("p11_current_analysis_audit", audit_cmd))

    public_output_dir = Path("results") / "pipeline11_current_mainline_check"
    public_outputs = write_public_outputs(audit_output_dir, public_output_dir)
    summary_path = public_output_dir / "summary.md"
    write_summary(summary_path, datasets, results, public_outputs)
    print("")
    print(f"Summary: {summary_path}")
    return 0 if all(result.returncode == 0 for result in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
