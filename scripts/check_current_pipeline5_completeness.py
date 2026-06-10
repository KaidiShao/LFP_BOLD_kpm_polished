#!/usr/bin/env python3
"""Check current Pipeline 5 completeness for the mainline dataset/method/k grid.

The checker treats Pipeline 5 reduction MAT files as source outputs and
figures/statistics as derived inspection outputs. Scatter trajectory figures
are considered current only if the summary PNG is newer than the scatter
renderer source file, unless --no-require-fresh-scatter is passed.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, List, Optional, Sequence


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_RESULTS_ROOT = Path("results")
DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "f12m01", "e10gw1")
DEFAULT_CONDITIONS = ("abs_projected_vlambda", "complex_split_projected_vlambda")
DEFAULT_METHODS = ("svd", "nmf", "mds", "umap")
DEFAULT_COMPONENT_COUNTS = tuple(range(3, 9))


@dataclass(frozen=True)
class CheckRow:
    dataset: str
    condition: str
    method: str
    component_count: str
    method_tag: str
    category: str
    item: str
    status: str
    path: str
    detail: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--conditions", nargs="+", default=list(DEFAULT_CONDITIONS))
    parser.add_argument("--methods", nargs="+", default=list(DEFAULT_METHODS))
    parser.add_argument("--component-counts", nargs="+", type=int, default=list(DEFAULT_COMPONENT_COUNTS))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Default: <repo>/results/current_pipeline5_completeness",
    )
    parser.add_argument(
        "--no-require-fresh-scatter",
        action="store_true",
        help="Count any existing scatter summary PNG as OK, even if it predates the current scatter renderer.",
    )
    parser.add_argument(
        "--skip-derived",
        action="store_true",
        help="Check reduction/scatter/default-density outputs only; skip peak/consistency derived layers.",
    )
    return parser.parse_args()


def method_tag(method: str, component_count: int) -> str:
    return f"{method}_k{component_count:02d}"


def latest_file(root: Path, pattern: str) -> Optional[Path]:
    if not root.is_dir():
        return None
    files = [path for path in root.glob(pattern) if path.is_file()]
    if not files:
        return None
    return max(files, key=lambda path: path.stat().st_mtime)


def latest_recursive(root: Path, pattern: str) -> Optional[Path]:
    if not root.is_dir():
        return None
    files = [path for path in root.rglob(pattern) if path.is_file()]
    if not files:
        return None
    return max(files, key=lambda path: path.stat().st_mtime)


def row(
    dataset: str,
    condition: str,
    method: str,
    component_count: object,
    tag: str,
    category: str,
    item: str,
    status: str,
    path: Optional[Path],
    detail: str = "",
) -> CheckRow:
    return CheckRow(
        dataset=dataset,
        condition=condition,
        method=method,
        component_count="" if component_count is None else str(component_count),
        method_tag=tag,
        category=category,
        item=item,
        status=status,
        path="" if path is None else str(path),
        detail=detail,
    )


def freshness_status(path: Optional[Path], fresh_after: Optional[float]) -> tuple[str, str]:
    if path is None:
        return "missing", ""
    if fresh_after is None:
        return "ok", ""
    mtime = path.stat().st_mtime
    if mtime >= fresh_after:
        return "ok", "fresh"
    return "stale", "older than current scatter renderer"


def check_method_outputs(
    processed_root: Path,
    dataset: str,
    condition: str,
    method: str,
    k: int,
    state_renderer_mtime: Optional[float],
    consensus_renderer_mtime: Optional[float],
) -> List[CheckRow]:
    rows: List[CheckRow] = []
    tag = method_tag(method, k)
    dataset_root = processed_root / dataset
    method_root = dataset_root / "pipeline5_eigenfunction_reduction" / condition / tag
    mat_dir = method_root / "mat"
    fig_dir = method_root / "fig"
    mat_file = latest_file(mat_dir, "*_efun__time__*.mat") or latest_file(mat_dir, "*_efun__spectrum__*.mat")
    rows.append(
        row(dataset, condition, method, k, tag, "source", "reduction_mat", "ok" if mat_file else "missing", mat_file)
    )

    summary_root = dataset_root / "pipeline5_summary_figures" / condition / "eigenfunction_reduction"
    state_png = latest_file(summary_root / "state_space", f"{tag}__*.png")
    status, detail = freshness_status(state_png, state_renderer_mtime)
    rows.append(row(dataset, condition, method, k, tag, "scatter", "state_space_scatter_png", status, state_png, detail))

    consensus_png = latest_file(summary_root / "consensus_state_space", f"{tag}__*.png")
    status, detail = freshness_status(consensus_png, consensus_renderer_mtime)
    rows.append(
        row(
            dataset,
            condition,
            method,
            k,
            tag,
            "scatter",
            "consensus_state_space_scatter_png",
            status,
            consensus_png,
            detail,
        )
    )

    if method in {"mds", "umap"}:
        spec_png = latest_file(summary_root / "spectrum_diagnostics", f"{tag}__*.png")
        rows.append(
            row(
                dataset,
                condition,
                method,
                k,
                tag,
                "default_figure",
                "spectrum_diagnostics_png",
                "ok" if spec_png else "missing",
                spec_png,
            )
        )

    dimred_density_root = dataset_root / "pipeline5_dimred_thresholded_density" / condition / tag
    dimred_mat = latest_recursive(dimred_density_root / "mat", f"*ratio_070*{condition}_{tag}*.mat")
    dimred_fig = latest_recursive(dimred_density_root / "fig", f"*ratio_070*{condition}_{tag}*.png")
    rows.append(
        row(dataset, condition, method, k, tag, "density", "dimred_thresholded_density_mat", "ok" if dimred_mat else "missing", dimred_mat)
    )
    rows.append(
        row(dataset, condition, method, k, tag, "density", "dimred_thresholded_density_png", "ok" if dimred_fig else "missing", dimred_fig)
    )

    return rows


def check_condition_outputs(processed_root: Path, dataset: str, condition: str) -> List[CheckRow]:
    rows: List[CheckRow] = []
    dataset_root = processed_root / dataset
    raw_density_root = dataset_root / "pipeline5_raw_thresholded_density" / condition
    raw_mat = latest_recursive(raw_density_root / "mat", f"*ratio_070*{condition}*.mat")
    raw_fig = latest_recursive(raw_density_root / "fig", f"*ratio_070*{condition}*.png")
    rows.append(row(dataset, condition, "", None, "", "density", "raw_thresholded_density_mat", "ok" if raw_mat else "missing", raw_mat))
    rows.append(row(dataset, condition, "", None, "", "density", "raw_thresholded_density_png", "ok" if raw_fig else "missing", raw_fig))
    return rows


def check_peak_outputs(processed_root: Path, dataset: str, condition: str, method: str, k: int) -> List[CheckRow]:
    rows: List[CheckRow] = []
    tag = method_tag(method, k)
    peak_root = processed_root / dataset / "pipeline5_eigenfunction_peaks_by_state_maxabs" / condition / tag
    save_tag = f"{dataset}_{condition}_{tag}_peaks"
    expected = {
        "peak_stats_csv": peak_root / f"{save_tag}_stats.csv",
        "peak_dist_png": peak_root / f"{save_tag}_dist.png",
        "peak_mean_png": peak_root / f"{save_tag}_mean.png",
        "peak_effect_png": peak_root / f"{save_tag}_effect.png",
    }
    for item, path in expected.items():
        rows.append(row(dataset, condition, method, k, tag, "peak_statistics", item, "ok" if path.is_file() else "missing", path))

    summary_png = (
        processed_root
        / "summary_figures"
        / "pipeline5_peak_statistics_maxabs"
        / f"{dataset}__{condition}__{tag}__peak_stats_summary.png"
    )
    rows.append(
        row(
            dataset,
            condition,
            method,
            k,
            tag,
            "peak_statistics",
            "peak_summary_png",
            "ok" if summary_png.is_file() else "missing",
            summary_png,
        )
    )
    return rows


def check_consistency_outputs(results_root: Path, dataset: str, k: int) -> List[CheckRow]:
    rows: List[CheckRow] = []
    kk = f"{k:02d}"
    per_dataset_root = results_root / f"peak_state_method_consistency_current_maxabs_k{kk}_{dataset}"
    for item in ("event_match_rate.png", "run_matched_events.png", "method_event_match_rate.png"):
        path = per_dataset_root / item
        rows.append(row(dataset, "", "", k, "", "method_consistency", item, "ok" if path.is_file() else "missing", path))
    return rows


def check_cross_dataset_outputs(results_root: Path, k: int) -> List[CheckRow]:
    rows: List[CheckRow] = []
    kk = f"{k:02d}"
    root = results_root / f"peak_state_method_consistency_current_all_datasets_maxabs_k{kk}"
    for item in (
        "dataset_event_match_rate.png",
        "dataset_event_median_score.png",
        "dataset_method_avg_matched_events.png",
        "method_event_avg_match_rate.png",
    ):
        path = root / item
        rows.append(row("ALL", "", "", k, "", "cross_dataset_consistency", item, "ok" if path.is_file() else "missing", path))
    return rows


def write_rows(path: Path, rows: Sequence[CheckRow]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(CheckRow.__dataclass_fields__.keys()))
        writer.writeheader()
        for item in rows:
            writer.writerow(item.__dict__)


def summarize(rows: Iterable[CheckRow]) -> str:
    by_category: dict[str, Counter[str]] = defaultdict(Counter)
    for item in rows:
        by_category[item.category][item.status] += 1
    lines = []
    for category in sorted(by_category):
        counts = by_category[category]
        total = sum(counts.values())
        parts = ", ".join(f"{key}={counts[key]}" for key in sorted(counts))
        lines.append(f"{category}: total={total}, {parts}")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    processed_root = args.processed_root.resolve()
    results_root = args.results_root.resolve()
    output_dir = args.output_dir or (results_root / "current_pipeline5_completeness")
    output_dir = output_dir.resolve()

    script_root = Path(__file__).resolve().parent.parent
    state_renderer = script_root / "functions" / "plottings" / "plot_eigenfunction_state_space_trajectory.m"
    consensus_renderer = script_root / "functions" / "plottings" / "plot_eigenfunction_state_space_consensus_trajectory.m"
    state_mtime = None if args.no_require_fresh_scatter else state_renderer.stat().st_mtime
    consensus_mtime = None if args.no_require_fresh_scatter else consensus_renderer.stat().st_mtime

    rows: List[CheckRow] = []
    for dataset in args.datasets:
        for condition in args.conditions:
            rows.extend(check_condition_outputs(processed_root, dataset, condition))
            for k in args.component_counts:
                for method in args.methods:
                    rows.extend(
                        check_method_outputs(
                            processed_root,
                            dataset,
                            condition,
                            method,
                            k,
                            state_mtime,
                            consensus_mtime,
                        )
                    )
                    if not args.skip_derived:
                        rows.extend(check_peak_outputs(processed_root, dataset, condition, method, k))
        if not args.skip_derived:
            for k in args.component_counts:
                rows.extend(check_consistency_outputs(results_root, dataset, k))

    if not args.skip_derived:
        for k in args.component_counts:
            rows.extend(check_cross_dataset_outputs(results_root, k))

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report = output_dir / f"pipeline5_completeness_{stamp}.csv"
    latest = output_dir / "pipeline5_completeness_latest.csv"
    missing = output_dir / "pipeline5_completeness_missing_latest.csv"
    scatter_missing = output_dir / "pipeline5_scatter_missing_latest.csv"
    write_rows(report, rows)
    write_rows(latest, rows)
    write_rows(missing, [item for item in rows if item.status != "ok"])
    write_rows(scatter_missing, [item for item in rows if item.category == "scatter" and item.status != "ok"])

    print(f"Processed root: {processed_root}")
    print(f"Results root  : {results_root}")
    print(f"Report        : {report}")
    print(f"Latest report : {latest}")
    print(f"Missing report: {missing}")
    print(f"Scatter miss  : {scatter_missing}")
    print("")
    print(summarize(rows))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
