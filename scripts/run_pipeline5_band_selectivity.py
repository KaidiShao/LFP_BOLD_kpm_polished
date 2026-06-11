#!/usr/bin/env python3
"""Run the formal P5 band-selectivity analysis for one dataset.

This wrapper turns the SOP into a single runnable entry point:

1. Build P2-band event response labels for P5 dimred components.
2. Generate the required P5 band-selectivity summary figures.
3. Generate component-pair diagnostic plates for method-k settings that have
   both strict theta-like and RG-like candidates.
4. Optionally generate paired top-window validation sheets.
5. Optionally generate P2/P5 top-30 window QC sheets.

Run from WSL/conda Python, for example:

python scripts/run_pipeline5_band_selectivity.py --dataset e10gb1
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path


DEFAULT_CONDITION = "complex_split_projected_vlambda_standardize"
DEFAULT_LABEL_SUFFIX = "p5_p2_band_event_response_v2_selective_envelope_20260528"
AVAILABLE_TRANSFORMS = ["abs", "adaptive_envelope"]
DEFAULT_TRANSFORMS = ["adaptive_envelope"]
DEFAULT_METHODS = ["svd", "nmf", "mds", "umap"]
DEFAULT_KS = list(range(3, 17))


@dataclass
class StepResult:
    name: str
    status: str
    seconds: float
    command: list[str]
    returncode: int | None = None
    note: str = ""


def main() -> None:
    args = parse_args()
    workspace = Path(args.workspace)
    processed_root = Path(args.processed_root)
    script_dir = Path(__file__).resolve().parent
    dataset_root = processed_root / args.dataset
    label_root = workspace / "results" / f"{args.dataset}_{DEFAULT_LABEL_SUFFIX}"
    manifest_file = label_root / "pipeline5_band_selectivity_run_manifest.json"

    label_root.mkdir(parents=True, exist_ok=True)

    run_manifest: dict[str, object] = {
        "dataset": args.dataset,
        "condition": args.condition,
        "component_counts": args.component_counts,
        "workspace": str(workspace),
        "processed_root": str(processed_root),
        "started_unix": time.time(),
        "steps": [],
        "input_check": check_inputs(dataset_root, args.dataset, args.condition, args.component_counts),
    }

    if args.dry_run:
        print(json.dumps(run_manifest["input_check"], indent=2))
        write_manifest(manifest_file, run_manifest)
        print(f"Dry-run manifest: {manifest_file}")
        return

    if not args.skip_labels:
        cmd = [
            sys.executable,
            str(script_dir / "analyze_e10gb1_p5_dimred_band_event_response_v2.py"),
            "--dataset",
            args.dataset,
            "--condition",
            args.condition,
            "--workspace",
            str(workspace),
            "--processed-root",
            str(processed_root),
            "--component-counts",
            *[str(k) for k in args.component_counts],
        ]
        if args.summary_only:
            cmd.append("--summary-only")
        append_result(run_manifest, run_step("p5_band_labels_and_summary_figures", cmd, required=True))

    if not args.skip_diagnostic_plates:
        for transform in args.transforms:
            for ribbon in args.ribbons:
                cmd = [
                    sys.executable,
                    str(script_dir / "plot_p5_event_diagnostic_plate.py"),
                    "--dataset",
                    args.dataset,
                    "--condition",
                    args.condition,
                    "--workspace",
                    str(workspace),
                    "--processed-root",
                    str(processed_root),
                    "--transform",
                    transform,
                    "--ribbon",
                    ribbon,
                    "--all-method-k",
                ]
                append_result(
                    run_manifest,
                    run_step(
                        f"p5_event_diagnostic_plates_{transform}_{ribbon}",
                        cmd,
                        required=False,
                    ),
                )

    if args.include_paired_top_windows and not args.skip_paired_top_windows:
        cmd = [
            sys.executable,
            str(script_dir / "plot_p5_paired_subprocess_traces.py"),
            "--dataset",
            args.dataset,
            "--condition",
            args.condition,
            "--workspace",
            str(workspace),
            "--processed-root",
            str(processed_root),
            "--transforms",
            *args.transforms,
            "--all-method-k",
        ]
        append_result(run_manifest, run_step("p5_paired_subprocess_top_windows", cmd, required=False))

    if not args.skip_p2_p5_top30:
        cmd = [
            sys.executable,
            str(script_dir / "plot_p2_p5_top30_windows.py"),
            "--dataset",
            args.dataset,
            "--condition",
            args.condition,
            "--workspace",
            str(workspace),
            "--processed-root",
            str(processed_root),
            "--transforms",
            *args.transforms,
        ]
        append_result(run_manifest, run_step("p2_p5_top30_window_qc", cmd, required=False))

    run_manifest["finished_unix"] = time.time()
    write_manifest(manifest_file, run_manifest)
    print(f"P5 band-selectivity pipeline manifest: {manifest_file}")

    failed_required = [
        step for step in run_manifest["steps"]
        if isinstance(step, dict) and step.get("status") == "failed_required"
    ]
    if failed_required:
        raise SystemExit(1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run P5 band-selectivity pipeline for one dataset.")
    parser.add_argument("--dataset", required=True, help="Dataset folder name under processed root.")
    parser.add_argument("--condition", default=DEFAULT_CONDITION)
    parser.add_argument("--workspace", default="/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
    parser.add_argument("--processed-root", default="/mnt/e/DataPons_processed")
    parser.add_argument("--transforms", nargs="+", default=DEFAULT_TRANSFORMS, choices=AVAILABLE_TRANSFORMS)
    parser.add_argument("--component-counts", nargs="+", type=int, default=DEFAULT_KS)
    parser.add_argument("--ribbons", nargs="+", default=["sem"], choices=["sem", "std", "iqr", "none"])
    parser.add_argument("--summary-only", action="store_true", help="Only generate minimal label summary figures.")
    parser.add_argument("--skip-labels", action="store_true")
    parser.add_argument("--skip-diagnostic-plates", action="store_true")
    parser.add_argument(
        "--include-paired-top-windows",
        action="store_true",
        help=(
            "Generate the heavy paired subprocess top-window validation sheets. "
            "Default is off; use this for selected method-k or suspicious datasets."
        ),
    )
    parser.add_argument("--skip-paired-top-windows", action="store_true")
    parser.add_argument("--skip-p2-p5-top30", action="store_true")
    parser.add_argument("--dry-run", action="store_true", help="Only report input availability.")
    return parser.parse_args()


def check_inputs(
    dataset_root: Path,
    dataset: str,
    condition: str,
    component_counts: list[int] | None = None,
) -> dict[str, object]:
    if component_counts is None:
        component_counts = DEFAULT_KS
    p2_event_file = dataset_root / "pipeline2_event_detection" / f"{dataset}_bandpass_events_3bands.mat"
    p2_top_csv = (
        dataset_root
        / "pipeline2_consensus_state_diversity_windows"
        / f"{dataset}_consensus_state_diversity_windows_6000samp_globalwin_top.csv"
    )
    p5_root = dataset_root / "pipeline5_eigenfunction_reduction" / condition
    method_k: dict[str, dict[str, object]] = {}
    for method in DEFAULT_METHODS:
        for k in component_counts:
            key = f"{method}_k{k:02d}"
            mat_dir = p5_root / key / "mat"
            mats = sorted(mat_dir.glob("*.mat")) if mat_dir.exists() else []
            method_k[key] = {
                "mat_dir": str(mat_dir),
                "exists": bool(mats),
                "n_mat": len(mats),
            }
    return {
        "p2_event_file": str(p2_event_file),
        "p2_event_file_exists": p2_event_file.exists(),
        "p2_top_window_csv": str(p2_top_csv),
        "p2_top_window_csv_exists": p2_top_csv.exists(),
        "p5_condition_root": str(p5_root),
        "p5_condition_root_exists": p5_root.exists(),
        "method_k": method_k,
    }


def run_step(name: str, cmd: list[str], required: bool) -> StepResult:
    start = time.time()
    print(f"\n=== {name} ===")
    print(" ".join(cmd))
    try:
        proc = subprocess.run(cmd, check=False)
    except Exception as exc:  # noqa: BLE001 - manifest should capture unexpected launch failures.
        seconds = time.time() - start
        status = "failed_required" if required else "failed_optional"
        note = f"{type(exc).__name__}: {exc}"
        print(f"{name}: {status}: {note}")
        return StepResult(name=name, status=status, seconds=seconds, command=cmd, returncode=None, note=note)

    seconds = time.time() - start
    if proc.returncode == 0:
        status = "completed"
    else:
        status = "failed_required" if required else "failed_optional"
    print(f"{name}: {status} ({seconds:.1f}s, returncode={proc.returncode})")
    return StepResult(name=name, status=status, seconds=seconds, command=cmd, returncode=proc.returncode)


def append_result(manifest: dict[str, object], result: StepResult) -> None:
    steps = manifest.setdefault("steps", [])
    if not isinstance(steps, list):
        raise TypeError("Manifest steps field is not a list.")
    steps.append(asdict(result))


def write_manifest(path: Path, manifest: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
