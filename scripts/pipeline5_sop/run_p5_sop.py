#!/usr/bin/env python3
"""One-command runner for the current P5 SOP figures.

Default scope:
    - standardized complex-split P5
    - k03:k16
    - adaptive_envelope only
    - strict P2-band labels

For each dataset, this runner generates:
    1. P5 strict gate CSV and summary figures.
    2. Diagnostic plates for method-k settings that have both theta-like and
       RG-like strict candidates.

Optionally, after all datasets, it can refresh the cross-dataset theta/RG
density-correlation CSVs and SOP heatmap.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
WORKSPACE = SCRIPT_DIR.parents[1]
DEFAULT_WORKSPACE = str(WORKSPACE)
DEFAULT_PROCESSED_ROOT = "/mnt/e/DataPons_processed" if str(WORKSPACE).startswith("/mnt/") else "E:/DataPons_processed"

DEFAULT_CONDITION = "complex_split_projected_vlambda_standardize"
DEFAULT_TRANSFORM = "adaptive_envelope"
DEFAULT_COMPONENT_COUNTS = list(range(3, 17))
DEFAULT_RIBBON = "sem"


@dataclass
class StepResult:
    name: str
    status: str
    seconds: float
    command: list[str]
    returncode: int | None = None
    note: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run current P5 SOP figures for one or more datasets.")
    parser.add_argument("--datasets", nargs="+", required=True, help="Dataset names, e.g. e10gb1 f12m05.")
    parser.add_argument("--condition", default=DEFAULT_CONDITION)
    parser.add_argument("--workspace", default=DEFAULT_WORKSPACE)
    parser.add_argument("--processed-root", default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--component-counts", nargs="+", type=int, default=DEFAULT_COMPONENT_COUNTS)
    parser.add_argument("--transform", default=DEFAULT_TRANSFORM, choices=["adaptive_envelope", "abs"])
    parser.add_argument("--ribbon", default=DEFAULT_RIBBON, choices=["sem", "std", "iqr", "none"])
    parser.add_argument("--skip-labels", action="store_true")
    parser.add_argument("--skip-diagnostic-plates", action="store_true")
    parser.add_argument("--refresh-density-correlation", action="store_true")
    parser.add_argument("--skip-density-correlation-plot", action="store_true")
    parser.add_argument(
        "--manifest",
        default="",
        help="Optional manifest JSON path. Default writes under results/p5_sop_runs.",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    workspace = Path(args.workspace)
    run_id = time.strftime("%Y%m%d_%H%M%S")
    manifest_path = Path(args.manifest) if args.manifest else (
        workspace / "results" / "p5_sop_runs" / f"p5_sop_run_{run_id}.json"
    )
    manifest: dict[str, object] = {
        "started_unix": time.time(),
        "datasets": args.datasets,
        "condition": args.condition,
        "component_counts": args.component_counts,
        "transform": args.transform,
        "ribbon": args.ribbon,
        "workspace": str(workspace),
        "processed_root": args.processed_root,
        "steps": [],
    }

    for dataset in args.datasets:
        dataset_steps = build_dataset_steps(args, dataset)
        for name, cmd, required in dataset_steps:
            if args.dry_run:
                result = StepResult(name=name, status="dry_run", seconds=0.0, command=cmd)
                print(f"[dry-run] {name}: {' '.join(cmd)}")
            else:
                result = run_step(name, cmd, required=required)
            append_result(manifest, result)
            write_manifest(manifest_path, manifest)

    if args.refresh_density_correlation:
        density_steps = build_density_correlation_steps(args)
        for name, cmd, required in density_steps:
            if args.dry_run:
                result = StepResult(name=name, status="dry_run", seconds=0.0, command=cmd)
                print(f"[dry-run] {name}: {' '.join(cmd)}")
            else:
                result = run_step(name, cmd, required=required)
            append_result(manifest, result)
            write_manifest(manifest_path, manifest)

    manifest["finished_unix"] = time.time()
    write_manifest(manifest_path, manifest)
    print(f"P5 SOP manifest: {manifest_path}")

    failed_required = [
        step for step in manifest["steps"]
        if isinstance(step, dict) and step.get("status") == "failed_required"
    ]
    if failed_required:
        raise SystemExit(1)


def build_dataset_steps(args: argparse.Namespace, dataset: str) -> list[tuple[str, list[str], bool]]:
    steps: list[tuple[str, list[str], bool]] = []
    component_counts = [str(k) for k in args.component_counts]

    if not args.skip_labels:
        steps.append(
            (
                f"{dataset}:p5_gate_labels_and_figures",
                [
                    sys.executable,
                    str(SCRIPT_DIR / "analyze_band_event_response.py"),
                    "--dataset",
                    dataset,
                    "--condition",
                    args.condition,
                    "--workspace",
                    args.workspace,
                    "--processed-root",
                    args.processed_root,
                    "--component-counts",
                    *component_counts,
                    "--transforms",
                    args.transform,
                ],
                True,
            )
        )

    if not args.skip_diagnostic_plates:
        steps.append(
            (
                f"{dataset}:p5_pair_pass_diagnostic_plates",
                [
                    sys.executable,
                    str(SCRIPT_DIR / "plot_event_diagnostic_plate.py"),
                    "--dataset",
                    dataset,
                    "--condition",
                    args.condition,
                    "--workspace",
                    args.workspace,
                    "--processed-root",
                    args.processed_root,
                    "--transform",
                    args.transform,
                    "--ribbon",
                    args.ribbon,
                    "--all-method-k",
                ],
                False,
            )
        )
    return steps


def build_density_correlation_steps(args: argparse.Namespace) -> list[tuple[str, list[str], bool]]:
    steps = [
        (
            "cross_dataset:theta_rg_density_correlation_tables",
            [
                sys.executable,
                str(SCRIPT_DIR / "analyze_theta_rg_density_correlation.py"),
                "--workspace",
                args.workspace,
                "--processed-root",
                args.processed_root,
            ],
            True,
        )
    ]
    if not args.skip_density_correlation_plot:
        input_dir = (
            Path(args.workspace)
            / "results"
            / "standardized_csplit_k03_k16_all_current_20260607"
            / "p5_theta_rg_density_correlation"
        )
        steps.append(
            (
                "cross_dataset:theta_rg_density_correlation_sop_figures",
                [
                    sys.executable,
                    str(SCRIPT_DIR / "plot_theta_rg_density_correlation.py"),
                    "--input-dir",
                    str(input_dir),
                ],
                True,
            )
        )
    return steps


def run_step(name: str, cmd: list[str], required: bool) -> StepResult:
    start = time.time()
    print(f"\n=== {name} ===", flush=True)
    print(" ".join(cmd), flush=True)
    try:
        proc = subprocess.run(cmd, check=False)
    except Exception as exc:  # noqa: BLE001 - manifest should capture launch failures.
        seconds = time.time() - start
        status = "failed_required" if required else "failed_optional"
        note = f"{type(exc).__name__}: {exc}"
        print(f"{name}: {status}: {note}", flush=True)
        return StepResult(name=name, status=status, seconds=seconds, command=cmd, returncode=None, note=note)

    seconds = time.time() - start
    status = "completed" if proc.returncode == 0 else ("failed_required" if required else "failed_optional")
    print(f"{name}: {status} ({seconds:.1f}s, returncode={proc.returncode})", flush=True)
    return StepResult(name=name, status=status, seconds=seconds, command=cmd, returncode=proc.returncode)


def append_result(manifest: dict[str, object], result: StepResult) -> None:
    steps = manifest.setdefault("steps", [])
    if not isinstance(steps, list):
        raise TypeError("manifest['steps'] is not a list")
    steps.append(asdict(result))


def write_manifest(path: Path, manifest: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
