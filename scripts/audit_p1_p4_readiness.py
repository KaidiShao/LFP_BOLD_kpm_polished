#!/usr/bin/env python3
"""Audit P1-P4 readiness for the current multi-dataset mainline.

This complements the downstream P11/P5-P10 audit.  It checks the front half of
the pipeline at product level instead of only checking whether a stage folder
exists.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Sequence


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_AUTODL_BLP_ROOT = Path(r"E:\autodl_results_new")
DEFAULT_AUTODL_BOLD_ROOT = Path(r"E:\autodl_results_local\bold_wsl")
DEFAULT_OUTPUT_DIR = Path("results") / "pipeline11_p1_p4_readiness_current"
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
CURRENT_BOLD_OBSERVABLES = ("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean")
CURRENT_BLP_CONDITIONS = ("abs", "complex_split")
CURRENT_P4_PARAM_MODES = ("raw", "standardize")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--autodl-blp-root", type=Path, default=DEFAULT_AUTODL_BLP_ROOT)
    parser.add_argument("--autodl-bold-root", type=Path, default=DEFAULT_AUTODL_BOLD_ROOT)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS_ALL9))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


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


def count_files(path: Path, pattern: str = "*") -> int:
    if not path.exists():
        return 0
    return sum(1 for p in path.glob(pattern) if p.is_file())


def has_file(path: Path, pattern: str) -> tuple[bool, str]:
    if not path.is_dir():
        return False, ""
    matches = sorted(p for p in path.glob(pattern) if p.is_file() and ".corrupt_" not in p.name.lower())
    return bool(matches), str(matches[-1]) if matches else ""


def has_recursive_file(path: Path, pattern: str = "*") -> tuple[bool, str, int]:
    if not path.is_dir():
        return False, "", 0
    matches = sorted(p for p in path.rglob(pattern) if p.is_file() and ".corrupt_" not in p.name.lower())
    return bool(matches), str(matches[-1]) if matches else "", len(matches)


def row(dataset: str, pipeline: str, item: str, status: str, path: str = "", detail: str = "", count: int = 0) -> dict[str, object]:
    return {
        "dataset": dataset,
        "pipeline": pipeline,
        "item": item,
        "status": status,
        "path": path,
        "count": count,
        "detail": detail,
    }


def p4_run_param_mode(path_or_name: Path | str) -> str:
    name = path_or_name.name if isinstance(path_or_name, Path) else str(path_or_name)
    lower = name.lower()
    if (
        re.search(r"(^|[_-])std(t|complex|[_-]|$)", lower)
        or "standardize" in lower
        or "standardise" in lower
        or "time_standard" in lower
    ):
        return "standardize"
    return "raw"


def p4_matches_mode(path: Path, param_mode: str) -> bool:
    return p4_run_param_mode(path) == param_mode


def audit_p1(dataset: str, ds_root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    spec_root = ds_root / "pipeline1_spectrograms"
    dict_root = ds_root / "pipeline1_reskoopnet_dictionary"
    for condition, spec_token, dict_token in (
        ("abs", "spectrograms_abs.mat", "abs_single"),
        ("complex_split", "spectrograms_complex.mat", "complex_split_single"),
    ):
        ok, path = has_file(spec_root, f"*{spec_token}")
        rows.append(row(dataset, "P1", f"spectrogram:{condition}", "ok" if ok else "missing", path, "required P1 spectrogram source"))
        ok_mat, mat_path = has_file(dict_root, f"*{dict_token}.mat")
        ok_csv, csv_path = has_file(dict_root, f"*{dict_token}_obs_info.csv")
        status = "ok" if ok_mat and ok_csv else "missing"
        rows.append(
            row(
                dataset,
                "P1",
                f"dictionary:{condition}",
                status,
                mat_path or csv_path,
                f"mat={int(ok_mat)}; obs_info={int(ok_csv)}",
                int(ok_mat) + int(ok_csv),
            )
        )
    return rows


def audit_p2(dataset: str, ds_root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    required = (
        ("event_detection", "pipeline2_event_detection"),
        ("event_density", "pipeline2_event_density"),
        ("consensus_states", "pipeline2_consensus_states"),
        ("consensus_state_summary", "pipeline2_consensus_state_summary"),
        ("consensus_state_diversity_windows", "pipeline2_consensus_state_diversity_windows"),
    )
    for item, stage in required:
        ok, path, count = has_recursive_file(ds_root / stage)
        rows.append(row(dataset, "P2", item, "ok" if ok else "missing", path or str(ds_root / stage), "required for current P2/P8/P10 event-density context", count))
    # Event-diversity windows are explicitly minor/nonblocking.
    ok, path, count = has_recursive_file(ds_root / "pipeline2_event_diversity_windows")
    rows.append(row(dataset, "P2", "event_diversity_windows_minor", "minor_ok" if ok else "minor_missing", path or str(ds_root / "pipeline2_event_diversity_windows"), "minor/nonblocking by current policy", count))
    return rows


def audit_p3(dataset: str, ds_root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    obs_root = ds_root / "pipeline3_bold_observables"
    qc_root = ds_root / "pipeline3_figures_bold_pre_reskoopnet_qc"
    for observable in CURRENT_BOLD_OBSERVABLES:
        ok, path = has_file(obs_root, f"*bold_observables_{observable}.mat")
        rows.append(row(dataset, "P3", f"observable_mat:{observable}", "ok" if ok else "missing", path or str(obs_root), "current four-observable BOLD source"))
        ok_qc, qc_path, qc_count = has_recursive_file(qc_root, f"*{observable}*.png")
        rows.append(row(dataset, "P3", f"observable_qc:{observable}", "ok" if ok_qc else "missing", qc_path or str(qc_root), "current P3 QC figure", qc_count))
    return rows


def audit_p4(dataset: str, blp_root: Path, bold_root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    blp_outputs = blp_root / dataset / "mlp" / "outputs"
    for condition in CURRENT_BLP_CONDITIONS:
        all_matches = sorted(
            p
            for p in blp_outputs.glob(f"*projected_vlambda_{condition}")
            if p.is_dir()
        ) if blp_outputs.is_dir() else []
        matches = all_matches
        rows.append(
            row(
                dataset,
                "P4",
                f"blp_p4_output:{condition}",
                "ok" if matches else "missing",
                str(matches[-1]) if matches else str(blp_outputs),
                "current BLP P4 output root expected by P5/P6",
                len(matches),
            )
        )
        for param_mode in CURRENT_P4_PARAM_MODES:
            mode_matches = [p for p in all_matches if p4_matches_mode(p, param_mode)]
            rows.append(
                row(
                    dataset,
                    "P4",
                    f"blp_p4_output:{condition}:{param_mode}",
                    "ok" if mode_matches else "missing",
                    str(mode_matches[-1]) if mode_matches else str(blp_outputs),
                    "P4 BLP output split by raw vs standardize parameter mode",
                    len(mode_matches),
                )
            )
    bold_outputs = bold_root / dataset / "mlp" / "outputs"
    for observable in CURRENT_BOLD_OBSERVABLES:
        all_matches = sorted(
            p
            for p in bold_outputs.glob(f"*projected_vlambda_{observable}")
            if p.is_dir()
        ) if bold_outputs.is_dir() else []
        matches = all_matches
        rows.append(
            row(
                dataset,
                "P4",
                f"bold_p4_output:{observable}",
                "ok" if matches else "missing",
                str(matches[-1]) if matches else str(bold_outputs),
                "current BOLD P4 output root expected by P7",
                len(matches),
            )
        )
        for param_mode in CURRENT_P4_PARAM_MODES:
            mode_matches = [p for p in all_matches if p4_matches_mode(p, param_mode)]
            rows.append(
                row(
                    dataset,
                    "P4",
                    f"bold_p4_output:{observable}:{param_mode}",
                    "ok" if mode_matches else "missing",
                    str(mode_matches[-1]) if mode_matches else str(bold_outputs),
                    "P4 BOLD output split by raw vs standardize parameter mode",
                    len(mode_matches),
                )
            )
    return rows


def summarize(rows: Sequence[dict[str, object]], datasets: Sequence[str]) -> list[dict[str, object]]:
    by_dataset: dict[str, list[dict[str, object]]] = defaultdict(list)
    for r in rows:
        by_dataset[str(r["dataset"])].append(r)
    out: list[dict[str, object]] = []
    for dataset in datasets:
        ds_rows = by_dataset.get(dataset, [])
        missing = [r for r in ds_rows if r["status"] == "missing"]
        p1_missing = [r["item"] for r in missing if r["pipeline"] == "P1"]
        p2_missing = [r["item"] for r in missing if r["pipeline"] == "P2"]
        p3_missing = [r["item"] for r in missing if r["pipeline"] == "P3"]
        p4_missing = [r["item"] for r in missing if r["pipeline"] == "P4"]
        actions: list[str] = []
        if p1_missing:
            actions.append("run_or_rebuild_p1_spectrogram_dictionary")
        if p2_missing:
            actions.append("run_p2_event_density_consensus")
        if p3_missing:
            actions.append("run_p3_bold_observables_or_qc")
        if any("blp_p4" in item for item in p4_missing):
            actions.append("run_p4_blp_training_export")
        if any("bold_p4" in item for item in p4_missing):
            actions.append("run_p4_bold_training_export")
        out.append(
            {
                "dataset": dataset,
                "p1_missing": ";".join(p1_missing),
                "p2_missing": ";".join(p2_missing),
                "p3_missing": ";".join(p3_missing),
                "p4_missing": ";".join(p4_missing),
                "needs_frontend_recompute": int(bool(actions)),
                "recommended_actions": ";".join(actions),
                "n_missing_required_items": len(missing),
            }
        )
    return out


def main() -> int:
    args = parse_args()
    rows: list[dict[str, object]] = []
    for dataset in args.datasets:
        ds_root = args.processed_root / dataset
        rows.extend(audit_p1(dataset, ds_root))
        rows.extend(audit_p2(dataset, ds_root))
        rows.extend(audit_p3(dataset, ds_root))
        rows.extend(audit_p4(dataset, args.autodl_blp_root, args.autodl_bold_root))
    summary = summarize(rows, args.datasets)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "p1_p4_readiness_long.csv", rows)
    write_csv(args.output_dir / "p1_p4_recompute_requirements_by_dataset.csv", summary)
    lines = [
        "# P1-P4 Readiness Audit",
        "",
        f"- Datasets: `{', '.join(args.datasets)}`",
        f"- Required item rows: `{len(rows)}`",
        f"- Datasets with required front-end recompute: `{sum(int(r['needs_frontend_recompute']) for r in summary)}` / `{len(summary)}`",
        "",
        "## Missing Required Items",
        "",
    ]
    for r in summary:
        if int(r["needs_frontend_recompute"]):
            lines.append(f"- `{r['dataset']}`: {r['recommended_actions']} ({r['n_missing_required_items']} missing)")
    (args.output_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Rows: {len(rows)}")
    print(f"Datasets needing P1-P4 recompute: {sum(int(r['needs_frontend_recompute']) for r in summary)} / {len(summary)}")
    print(f"Output dir: {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
