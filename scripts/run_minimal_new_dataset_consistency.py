#!/usr/bin/env python3
"""Run the current minimal cross-dataset consistency readout.

This runner does not recompute P1-P10.  It assumes the current mainline outputs
already exist and generates the small set of analyses used by
docs/current_minimal_cross_dataset_consistency_sop_2026-06-04_zh.md.

Run from WSL/conda Python, for example:

    /home/kdshao/anaconda3/bin/python scripts/run_minimal_new_dataset_consistency.py \
        --datasets e10gb1 e10gh1 k13m18 \
        --tag standardized_csplit_minimal_20260604
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path


DEFAULT_WORKSPACE = Path("/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
DEFAULT_PROCESSED_ROOT = Path("/mnt/e/DataPons_processed")
DEFAULT_SUMMARY_ROOT = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
)
DEFAULT_ROI_PROFILE_CSVS = [
    Path("results/pipeline_roi_profile_consistency_current/p7_intrinsic_bold_efun_roi_profiles_long.csv"),
    Path("results/pipeline_roi_profile_consistency_current/p7_roi_mean_profiles_direct_from_bold_post_k13m23.csv"),
    Path("results/pipeline_roi_profile_consistency_k13m18_m21_20260604/p7_intrinsic_bold_efun_roi_profiles_long.csv"),
]
DEFAULT_P8_SUBPROCESS_ROI_PROFILE_CSV = Path(
    "results/pipeline_roi_profile_consistency_current/p8_roi_profiles_by_subprocess_long.csv"
)
DEFAULT_RUN_TAGS = ["pv_roi"]


@dataclass
class StepResult:
    name: str
    status: str
    seconds: float
    command: list[str]
    returncode: int | None = None
    note: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="+", required=True)
    parser.add_argument("--tag", default="minimal_new_dataset_consistency")
    parser.add_argument("--workspace", type=Path, default=DEFAULT_WORKSPACE)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--summary-root", type=Path, default=DEFAULT_SUMMARY_ROOT)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument(
        "--run-tags",
        nargs="+",
        default=list(DEFAULT_RUN_TAGS),
        help="BOLD P8/P10 run tags to include. The minimal SOP defaults to roi_mean only: pv_roi.",
    )
    parser.add_argument("--roi-profile-csv", type=Path, nargs="*", default=list(DEFAULT_ROI_PROFILE_CSVS))
    parser.add_argument("--p8-subprocess-roi-profile-csv", type=Path, default=DEFAULT_P8_SUBPROCESS_ROI_PROFILE_CSV)
    parser.add_argument(
        "--roi-confusion-min-datasets",
        type=int,
        default=0,
        help="Minimum datasets required for subprocess ROI confusion. 0 uses min(4, n requested datasets).",
    )
    parser.add_argument("--skip-p5", action="store_true")
    parser.add_argument("--skip-coupling", action="store_true")
    parser.add_argument("--skip-traces", action="store_true")
    parser.add_argument("--skip-roi", action="store_true")
    parser.add_argument("--component-counts", nargs="+", type=int, default=list(range(3, 9)))
    parser.add_argument(
        "--full-p5-qc",
        action="store_true",
        help="Also generate P2/P5 top30 QC sheets; default keeps P5 gate lighter.",
    )
    parser.add_argument(
        "--include-p5-paired-top-windows",
        action="store_true",
        help=(
            "Also generate heavy paired subprocess top-window sheets. "
            "Default is off; use only for selected method-k or suspicious datasets."
        ),
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    workspace = args.workspace
    script_dir = workspace / "scripts"
    result_root = workspace / "results" / args.tag
    figure_root = args.summary_root / args.tag
    manifest_path = result_root / "minimal_consistency_run_manifest.json"

    result_root.mkdir(parents=True, exist_ok=True)
    figure_root.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, object] = {
        "tag": args.tag,
        "datasets": args.datasets,
        "run_tags": args.run_tags,
        "workspace": str(workspace),
        "processed_root": str(args.processed_root),
        "result_root": str(result_root),
        "figure_root": str(figure_root),
        "started_unix": time.time(),
        "input_audit": audit_inputs(args),
        "steps": [],
    }

    commands = build_commands(args, script_dir, result_root, figure_root)
    if args.dry_run:
        manifest["planned_commands"] = commands
        write_manifest(manifest_path, manifest)
        print(json.dumps(manifest, indent=2))
        return

    for name, cmd, required in commands:
        result = run_step(name, cmd, required=required, cwd=workspace)
        manifest["steps"].append(asdict(result))
        write_manifest(manifest_path, manifest)
        if required and result.status != "ok":
            manifest["finished_unix"] = time.time()
            write_manifest(manifest_path, manifest)
            raise SystemExit(result.returncode or 1)

    manifest["finished_unix"] = time.time()
    write_manifest(manifest_path, manifest)
    print(f"Minimal consistency manifest: {manifest_path}")


def audit_inputs(args: argparse.Namespace) -> dict[str, object]:
    audit: dict[str, object] = {}
    for dataset in args.datasets:
        dataset_root = args.processed_root / dataset
        audit[dataset] = {
            "p2_events": str(dataset_root / "pipeline2_event_detection" / f"{dataset}_bandpass_events_3bands.mat"),
            "p2_events_exists": (dataset_root / "pipeline2_event_detection" / f"{dataset}_bandpass_events_3bands.mat").is_file(),
            "p5_std_csplit_root": str(dataset_root / "pipeline5_eigenfunction_reduction" / "complex_split_projected_vlambda_standardize"),
            "p5_std_csplit_root_exists": (dataset_root / "pipeline5_eigenfunction_reduction" / "complex_split_projected_vlambda_standardize").is_dir(),
            "p8_root": str(dataset_root / "pipeline8_xcorr"),
            "p8_root_exists": (dataset_root / "pipeline8_xcorr").is_dir(),
            "p10_root": str(dataset_root / "pipeline10_dimred_xcorr"),
            "p10_root_exists": (dataset_root / "pipeline10_dimred_xcorr").is_dir(),
            "p7_root": str(dataset_root / "pipeline7_bold_reskoopnet_postprocessing"),
            "p7_root_exists": (dataset_root / "pipeline7_bold_reskoopnet_postprocessing").is_dir(),
        }
    audit["roi_profile_csv"] = [
        {"path": str(path), "exists": path.is_file()}
        for path in args.roi_profile_csv
    ]
    audit["p8_subprocess_roi_profile_csv"] = {
        "path": str(args.p8_subprocess_roi_profile_csv),
        "exists": args.p8_subprocess_roi_profile_csv.is_file(),
    }
    return audit


def build_commands(
    args: argparse.Namespace,
    script_dir: Path,
    result_root: Path,
    figure_root: Path,
) -> list[tuple[str, list[str], bool]]:
    commands: list[tuple[str, list[str], bool]] = []

    if not args.skip_p5:
        for dataset in args.datasets:
            cmd = [
                args.python,
                str(script_dir / "run_pipeline5_band_selectivity.py"),
                "--dataset",
                dataset,
                "--workspace",
                str(args.workspace),
                        "--processed-root",
                        str(args.processed_root),
                        "--component-counts",
                        *[str(k) for k in args.component_counts],
                    ]
            if not args.full_p5_qc:
                cmd.append("--skip-p2-p5-top30")
            if args.include_p5_paired_top_windows:
                cmd.append("--include-paired-top-windows")
            commands.append((f"p5_band_selectivity_{dataset}", cmd, True))

    coupling_results = result_root / "p8_p10_strict_band_coupling"
    coupling_figures = figure_root / "p8_p10_strict_band_coupling"
    coupling_hits = coupling_results / "p8_p10_strict_band_hits_long.csv"
    if not args.skip_coupling:
        commands.append(
            (
                "p8_p10_strict_band_coupling",
                [
                    args.python,
                    str(script_dir / "analyze_p8_p10_strict_band_coupling.py"),
                    "--processed-root",
                    str(args.processed_root),
                    "--workspace-root",
                    str(args.workspace),
                    "--datasets",
                    *args.datasets,
                    "--run-tags",
                    *args.run_tags,
                    "--component-counts",
                    *[str(k) for k in args.component_counts],
                    "--results-dir",
                    str(coupling_results),
                    "--figure-dir",
                    str(coupling_figures),
                ],
                True,
            )
        )

    if not args.skip_traces:
        for dataset in args.datasets:
            for pipeline in ("P8", "P10"):
                commands.append(
                    (
                        f"{pipeline.lower()}_density_trace_{dataset}",
                        [
                            args.python,
                            str(script_dir / "plot_bold_density_trace_by_dataset_observable.py"),
                            "--hits-csv",
                            str(coupling_hits),
                            "--dataset",
                            dataset,
                            "--pipeline",
                            pipeline,
                            "--density-classes",
                            "raw_efun_density",
                            "dimred_efun_density",
                            "--results-dir",
                            str(result_root / f"{pipeline.lower()}_{dataset}_density_traces"),
                        ],
                        False,
                    )
                )

    existing_roi_csvs = [path for path in args.roi_profile_csv if path.is_file()]
    if not args.skip_roi and existing_roi_csvs:
        commands.append(
            (
                "p8_roi_mean_rg_no_theta_roi_summary",
                [
                    args.python,
                    str(script_dir / "analyze_p8_roi_mean_rg_no_theta_roi_summary.py"),
                    "--hits-csv",
                    str(coupling_hits),
                    "--roi-profile-csv",
                    *[str(path) for path in existing_roi_csvs],
                    "--datasets",
                    *args.datasets,
                    "--output-dir",
                    str(result_root / "p8_roi_mean_rg_no_theta_roi_summary"),
                    "--figure-dir",
                    str(figure_root / "p8_roi_mean_rg_no_theta_roi_summary"),
                    "--roi-order",
                    "p7_plot",
                    "--dataset-panel-cols",
                    "12",
                    "--dataset-panel-scale",
                    "per_dataset",
                ],
                False,
            )
        )

        commands.append(
            (
                "p8_raw_csplit_roi_efun_vs_deconv",
                [
                    args.python,
                    str(script_dir / "plot_p8_raw_density_efun_vs_deconv_roi_similarity.py"),
                    "--processed-root",
                    str(args.processed_root),
                    "--p7-profile-csv",
                    *[str(path) for path in existing_roi_csvs],
                    "--datasets",
                    *args.datasets,
                    "--run-tags",
                    *args.run_tags,
                    "--raw-densities",
                    "raw_csplit_q070",
                    "--output-dir",
                    str(result_root / "p8_raw_csplit_roi_efun_vs_deconv"),
                    "--figure-dir",
                    str(figure_root / "p8_raw_csplit_roi_efun_vs_deconv"),
                ],
                False,
            )
        )

    if not args.skip_roi and args.p8_subprocess_roi_profile_csv.is_file():
        roi_confusion_min_datasets = (
            args.roi_confusion_min_datasets
            if args.roi_confusion_min_datasets > 0
            else min(4, max(1, len(args.datasets)))
        )
        commands.append(
            (
                "p8_roi_subprocess_confusion",
                [
                    args.python,
                    str(script_dir / "analyze_p8_roi_subprocess_separation.py"),
                    "--profile-csv",
                    str(args.p8_subprocess_roi_profile_csv),
                    "--datasets",
                    *args.datasets,
                    "--run-tags",
                    *args.run_tags,
                    "--conditions",
                    "csplit",
                    "--label-groups",
                    "theta",
                    "ripple_gamma",
                    "mixed",
                    "inactive",
                    "--method-k-scope",
                    "p5-stable",
                    "--min-datasets",
                    str(roi_confusion_min_datasets),
                    "--output-dir",
                    str(result_root / "p8_roi_subprocess_confusion"),
                    "--figure-dir",
                    str(figure_root / "p8_roi_subprocess_confusion"),
                ],
                False,
            )
        )

        commands.append(
            (
                "p8_dimred_rg_roi_efun_vs_deconv",
                [
                    args.python,
                    str(script_dir / "plot_p8_efun_vs_deconv_top_xcorr_roi_difference.py"),
                    "--profile-csv",
                    str(args.p8_subprocess_roi_profile_csv),
                    "--datasets",
                    *args.datasets,
                    "--run-tags",
                    *args.run_tags,
                    "--observables",
                    "roi_mean",
                    "--conditions",
                    "csplit",
                    "--label-scopes",
                    "ripple_gamma",
                    "--method-k-scope",
                    "p5-stable",
                    "--output-dir",
                    str(result_root / "p8_dimred_rg_roi_efun_vs_deconv"),
                    "--figure-dir",
                    str(figure_root / "p8_dimred_rg_roi_efun_vs_deconv"),
                ],
                False,
            )
        )

    return commands


def run_step(name: str, cmd: list[str], *, required: bool, cwd: Path) -> StepResult:
    print(f"\n=== {name} ===")
    print(" ".join(cmd))
    start = time.time()
    try:
        proc = subprocess.run(cmd, cwd=cwd, check=False)
        seconds = time.time() - start
        if proc.returncode == 0:
            return StepResult(name=name, status="ok", seconds=seconds, command=cmd, returncode=0)
        status = "failed_required" if required else "failed_optional"
        return StepResult(name=name, status=status, seconds=seconds, command=cmd, returncode=proc.returncode)
    except Exception as exc:  # noqa: BLE001 - manifest should capture launch failures.
        seconds = time.time() - start
        status = "failed_required" if required else "failed_optional"
        return StepResult(name=name, status=status, seconds=seconds, command=cmd, note=str(exc))


def write_manifest(path: Path, manifest: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
