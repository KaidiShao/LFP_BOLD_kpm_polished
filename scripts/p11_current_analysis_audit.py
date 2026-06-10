#!/usr/bin/env python3
"""Audit and plan Pipeline 11 current-analysis refresh work.

The script is intentionally conservative by default:

- it reads canonical P1-P10 outputs under E:\\DataPons_processed;
- it writes audit CSVs under the repo-local results folder;
- it does not delete anything;
- it does not move files unless --quarantine-legacy is explicitly set;
- it does not copy flat figures unless --copy-flat-figures is explicitly set.

The generated fill plan is the handoff for automatic backfill: each action is
tagged with dependencies, a coarse memory class, and whether it can run in
parallel with other actions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Optional, Sequence


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_AUTODL_BLP_ROOT = Path(r"E:\autodl_results_new")
DEFAULT_AUTODL_BOLD_ROOT = Path(r"E:\autodl_results_local\bold_wsl")
DEFAULT_RESULTS_ROOT = Path("results") / "pipeline11_current_analysis_audit"
DEFAULT_REPO_RESULTS_ROOT = Path("results")
DEFAULT_P11_SUMMARY_ROOT = (
    DEFAULT_PROCESSED_ROOT / "summary_figures" / "pipeline11_current_analysis_summary"
)
DEFAULT_INVALID_BLP_RUNS_FILE = (
    Path(__file__).resolve().parents[1] / "configs" / "current_invalid_blp_mlp_runs.csv"
)
REPO_ROOT = Path(__file__).resolve().parents[1]
P5_ACTIVITY_METADATA_SOURCE_FILES = (
    REPO_ROOT / "functions" / "postprocessing" / "get_thresholded_density.m",
    REPO_ROOT / "functions" / "postprocessing" / "get_dimred_thresholded_density.m",
    REPO_ROOT / "functions" / "postprocessing" / "build_blp_eigenfunction_reduction_params.m",
    REPO_ROOT / "functions" / "postprocessing" / "build_blp_eigenfunction_reduction_cfg.m",
    REPO_ROOT / "functions" / "postprocessing" / "build_blp_raw_eigenfunction_threshold_cfg.m",
)
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
DEFAULT_DATASETS = DEFAULT_DATASETS_CURRENT5

CURRENT_P3_OBSERVABLES = ("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean")
CURRENT_P7_OBSERVABLES = CURRENT_P3_OBSERVABLES
CURRENT_P7_OBSERVABLE_TO_RUN_TAG = {
    "global_svd100": "pv_gsvd100",
    "gsvd100_ds": "pv_gsvd100_ds",
    "HP_svd100": "pv_hp100",
    "roi_mean": "pv_roi",
}
CURRENT_P7_RUN_TAGS = tuple(CURRENT_P7_OBSERVABLE_TO_RUN_TAG.values())
CURRENT_P5_CONDITIONS = ("abs_projected_vlambda", "complex_split_projected_vlambda")
CURRENT_P5_CONDITION_PARTS = {
    "abs_projected_vlambda": ("abs", "projected_vlambda"),
    "complex_split_projected_vlambda": ("complex_split", "projected_vlambda"),
}
CURRENT_P5_METHODS = ("svd", "nmf", "mds", "umap")
CURRENT_P5_COMPONENT_COUNTS = tuple(range(3, 9))
CURRENT_P10_FEATURES = ("efun_real",)
CURRENT_P10_METHODS = CURRENT_P5_METHODS
CURRENT_P10_COMPONENT_COUNTS = tuple(range(3, 9))
CURRENT_P4_PARAM_MODES = ("raw", "standardize")
CURRENT_P4_DOWNSTREAM_PARAM_MODE = "standardize"
P5_DERIVED_SELECTIVITY_RESULTS = (
    "peak_event_family_component_activity_current_maxabs",
)
P5_DERIVED_CONSISTENCY_PREFIX = "peak_state_method_consistency_current"
P6_CONDITION_TOKENS = {
    "abs_projected_vlambda": ("projected_vlambda_abs", "vlambda_abs"),
    "complex_split_projected_vlambda": (
        "projected_vlambda_complex_split",
        "vlambda_complex_split",
    ),
}

CANONICAL_DATASET_STAGES = (
    "pipeline1_spectrograms",
    "pipeline1_reskoopnet_dictionary",
    "pipeline2_event_detection",
    "pipeline2_event_density",
    "pipeline2_legacy_event_density",
    "pipeline2_figures_legacy_event_density",
    "pipeline2_consensus_states",
    "pipeline2_consensus_state_summary",
    "pipeline2_event_diversity_windows",
    "pipeline2_consensus_state_diversity_windows",
    "pipeline2_figures_event_diversity_top_window_plots",
    "pipeline2_figures_consensus_state_top_window_plots",
    "pipeline3_bold_observables",
    "pipeline3_figures_bold_pre_reskoopnet_qc",
    "pipeline3_figures_bold_pre_reskoopnet_qc_summary",
    "pipeline5_eigenfunction_reduction",
    "pipeline5_efun_dimred_top30",
    "pipeline5_raw_thresholded_density",
    "pipeline5_raw_thresholded_density_scan",
    "pipeline5_raw_thresholded_events",
    "pipeline5_dimred_thresholded_density",
    "pipeline5_dimred_thresholded_events",
    "pipeline5_eigenfunction_peaks_by_state",
    "pipeline5_figures_eigenfunction_reduction_ssc_summary",
    "pipeline6_top_state_diversity_postprocessing",
    "pipeline6_spkt_residual_cross_correlation",
    "pipeline6_spkt_residual_lagged_cross_correlation",
    "pipeline6_mua_residual_cross_correlation",
    "pipeline6_figures_timescale_diagnostics",
    "pipeline6_figures_spkt_residual_cross_correlation",
    "pipeline6_figures_spkt_residual_lagged_cross_correlation",
    "pipeline6_figures_mua_residual_cross_correlation",
    "pipeline7_bold_reskoopnet_postprocessing",
    "pipeline7_figures_bold_reskoopnet_main_summary",
    "pipeline7_figures_bold_reskoopnet_deconv_summary",
    "pipeline7_figures_bold_reskoopnet_timescale_summary",
    "pipeline8_xcorr",
    "pipeline8_top_maps",
    "pipeline8_xcorr_summary",
    "pipeline8_top_maps_summary",
    "pipeline9_bold_eigenfunction_reduction",
    "pipeline9_figures_bold_eigenfunction_reduction",
    "pipeline10_dimred_xcorr",
    "pipeline10_dimred_top_maps",
    "pipeline10_dimred_xcorr_summary",
    "pipeline10_dimred_top_maps_summary",
)

LEGACY_OPTIONAL_DATASET_STAGES = (
    "pipeline2_legacy_event_density",
    "pipeline2_figures_legacy_event_density",
)

MINOR_OPTIONAL_DATASET_STAGES = (
    "pipeline2_event_diversity_windows",
    "pipeline2_figures_event_diversity_top_window_plots",
    "pipeline5_raw_thresholded_events",
    "pipeline5_dimred_thresholded_events",
    "pipeline6_spkt_residual_lagged_cross_correlation",
    "pipeline6_figures_spkt_residual_lagged_cross_correlation",
)

DERIVED_CACHE_DATASET_STAGES = (
    "pipeline3_figures_bold_pre_reskoopnet_qc_summary",
    "pipeline5_eigenfunction_peaks_by_state",
    "pipeline5_figures_eigenfunction_reduction_ssc_summary",
    "pipeline7_figures_bold_reskoopnet_main_summary",
    "pipeline7_figures_bold_reskoopnet_deconv_summary",
    "pipeline7_figures_bold_reskoopnet_timescale_summary",
    "pipeline8_xcorr_summary",
    "pipeline8_top_maps_summary",
    "pipeline9_figures_bold_eigenfunction_reduction",
    "pipeline10_dimred_xcorr_summary",
    "pipeline10_dimred_top_maps_summary",
)

LAYOUT_DEBT_DATASET_STAGES = (
    "pipeline5_eigenfunction_peaks_by_state_maxabs",
    "pipeline5_summary_figures",
)

QUARANTINE_CANDIDATE_STAGES = (
    "pipeline8_bold_efun_density_cross_correlation",
    "pipeline8_figures_bold_top_xcorr_activation_maps",
    "pipeline8_flat_top_xcorr_figures",
)

NONCANONICAL_DATASET_STAGES = LAYOUT_DEBT_DATASET_STAGES + QUARANTINE_CANDIDATE_STAGES


@dataclass(frozen=True)
class StageStatus:
    dataset: str
    pipeline: str
    item: str
    status: str
    path: str
    file_count: int = 0
    dir_count: int = 0
    detail: str = ""


@dataclass(frozen=True)
class P4BoldRunSelection:
    run_path: Path
    status: str
    detail: str
    best_val_metric: Optional[float] = None
    best_outer_epoch: Optional[int] = None
    completed_outer_epochs: Optional[int] = None
    training_state_path: str = ""


@dataclass(frozen=True)
class FillAction:
    action_id: str
    dataset: str
    pipeline: str
    action: str
    priority: int
    memory_class: str
    parallel_group: str
    can_parallelize: str
    blocked_by: str
    writes_external: str
    command_hint: str
    reason: str


@dataclass(frozen=True)
class LegacyCandidate:
    dataset: str
    source_path: str
    proposed_legacy_path: str
    reason: str
    status: str


@dataclass(frozen=True)
class FlatFigureSource:
    dataset: str
    pipeline: str
    figure_group: str
    source_path: str
    relative_source: str
    proposed_flat_path: str
    status: str
    detail: str = ""


@dataclass(frozen=True)
class P10Candidate:
    dataset: str
    run_tag: str
    feature_name: str
    method_tag: str
    component_count: int
    dimred_result_file: Path
    p10_tag: str
    xcorr_dir: Path
    xcorr_file: Path
    top_maps_dir: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--repo-results-root", type=Path, default=DEFAULT_REPO_RESULTS_ROOT)
    parser.add_argument("--summary-root", type=Path, default=DEFAULT_P11_SUMMARY_ROOT)
    parser.add_argument("--autodl-blp-root", type=Path, default=DEFAULT_AUTODL_BLP_ROOT)
    parser.add_argument("--autodl-bold-root", type=Path, default=DEFAULT_AUTODL_BOLD_ROOT)
    parser.add_argument(
        "--dataset-scope",
        choices=("current5", "all9"),
        default="current5",
        help="Use current five downstream datasets or all nine cfg datasets as the default dataset list.",
    )
    parser.add_argument("--datasets", nargs="+", default=None)
    parser.add_argument("--copy-flat-figures", action="store_true")
    parser.add_argument("--quarantine-legacy", action="store_true")
    parser.add_argument(
        "--archive-existing-summary",
        action="store_true",
        help="Move an existing P11 summary folder into _legacy_quarantine before copying flat figures.",
    )
    parser.add_argument("--overwrite-flat", action="store_true")
    parser.add_argument("--max-flat-copies", type=int, default=None)
    args = parser.parse_args()
    if args.datasets is None:
        args.datasets = list(DEFAULT_DATASETS_ALL9 if args.dataset_scope == "all9" else DEFAULT_DATASETS_CURRENT5)
    return args


def safe_rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def count_files(root: Path) -> int:
    if not root.is_dir():
        return 0
    return sum(1 for path in root.rglob("*") if path.is_file())


def count_dirs(root: Path) -> int:
    if not root.is_dir():
        return 0
    return sum(1 for path in root.iterdir() if path.is_dir())


def stage_row(
    dataset: str,
    pipeline: str,
    item: str,
    path: Path,
    detail: str = "",
    ok_status: str = "ok",
    missing_status: str = "missing",
) -> StageStatus:
    exists = path.exists()
    return StageStatus(
        dataset=dataset,
        pipeline=pipeline,
        item=item,
        status=ok_status if exists else missing_status,
        path=str(path),
        file_count=count_files(path) if exists else 0,
        dir_count=count_dirs(path) if exists else 0,
        detail=detail,
    )


def latest_file(root: Path, pattern: str) -> Optional[Path]:
    if not root.is_dir():
        return None
    files = [path for path in root.glob(pattern) if path.is_file()]
    if not files:
        return None
    return max(files, key=lambda path: path.stat().st_mtime)


def newest_existing_mtime(paths: Sequence[Path]) -> float:
    mtimes = [path.stat().st_mtime for path in paths if path.is_file()]
    return max(mtimes) if mtimes else 0.0


def any_file(root: Path, pattern: str) -> bool:
    return latest_file(root, pattern) is not None


def any_file_recursive(root: Path, pattern: str) -> bool:
    if not root.is_dir():
        return False
    return any(path.is_file() for path in root.glob(pattern))


def current_artifact_status(path: Optional[Path], cutoff_mtime: float) -> tuple[str, str]:
    if path is None:
        return "missing", "missing current artifact"
    if cutoff_mtime and path.stat().st_mtime + 1 < cutoff_mtime:
        return "stale", "artifact predates current code/provenance rule"
    return "ok", "artifact exists and is newer than current code/provenance rule"


def find_p7_observable_run(p7_root: Path, observable: str) -> Optional[Path]:
    if not p7_root.is_dir():
        return None
    suffix = f"_projected_vlambda_{observable}".lower()
    matches = [
        path
        for path in p7_root.iterdir()
        if path.is_dir() and path.name.lower().endswith(suffix)
    ]
    if not matches:
        matches = [
            path
            for path in p7_root.iterdir()
            if path.is_dir() and observable.lower() in path.name.lower()
        ]
    if not matches:
        return None
    return max(matches, key=lambda path: path.stat().st_mtime)


def p6_condition_token(condition: str) -> str:
    tokens = P6_CONDITION_TOKENS.get(condition, (condition,))
    return tokens[0]


def path_contains_any_token(path: Path, tokens: Sequence[str]) -> bool:
    text = str(path).lower()
    return any(token.lower() in text for token in tokens)


@lru_cache(maxsize=1)
def invalid_blp_run_rows() -> tuple[dict[str, str], ...]:
    if not DEFAULT_INVALID_BLP_RUNS_FILE.is_file():
        return ()
    rows: list[dict[str, str]] = []
    with DEFAULT_INVALID_BLP_RUNS_FILE.open("r", newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            clean = {str(k or "").strip(): str(v or "").strip() for k, v in row.items()}
            status = clean.get("status", "invalid").lower()
            if status in {"invalid", "stale", "legacy"}:
                rows.append(clean)
    return tuple(rows)


def normalize_path_text(value: str) -> str:
    return value.strip().lower().replace("/", "\\").rstrip("\\")


def invalid_blp_path(path: Path) -> bool:
    text = normalize_path_text(str(path))
    for row in invalid_blp_run_rows():
        run_name = row.get("run_name", "").strip().lower()
        output_dir = normalize_path_text(row.get("output_dir", ""))
        if run_name and run_name in text:
            return True
        if output_dir and output_dir in text:
            return True
    return False


def invalid_blp_rows_for_condition(dataset: str, condition: str) -> tuple[dict[str, str], ...]:
    observable_mode, residual_form = CURRENT_P5_CONDITION_PARTS.get(condition, ("", ""))
    out: list[dict[str, str]] = []
    for row in invalid_blp_run_rows():
        row_dataset = row.get("dataset", "").lower()
        row_observable = row.get("observable_mode", "").lower()
        row_residual = row.get("residual_form", "").lower()
        dataset_ok = row_dataset in {"", "*", dataset.lower()}
        observable_ok = row_observable in {"", "*", observable_mode.lower()}
        residual_ok = row_residual in {"", "*", residual_form.lower()}
        if dataset_ok and observable_ok and residual_ok:
            out.append(row)
    return tuple(out)


def invalid_blp_condition_detail(dataset: str, condition: str) -> str:
    rows = invalid_blp_rows_for_condition(dataset, condition)
    if not rows:
        return ""
    parts = []
    for row in rows:
        status = row.get("status", "invalid")
        run_name = row.get("run_name", "")
        reason = row.get("reason", "")
        parts.append(f"{status}:{run_name} ({reason})")
    return "invalid/stale BLP P4 source excluded from current mainline: " + "; ".join(parts)


def p6_condition_has_top_windows(ds_root: Path, condition: str) -> bool:
    root = ds_root / "pipeline6_top_state_diversity_postprocessing"
    tokens = P6_CONDITION_TOKENS.get(condition, (condition,))
    if not root.is_dir():
        return False
    for run_dir in root.iterdir():
        if (
            not run_dir.is_dir()
            or not path_contains_any_token(run_dir, tokens)
            or invalid_blp_path(run_dir)
        ):
            continue
        has_main = any_file_recursive(run_dir / "01_postprocess_main", "*.png")
        has_deconv_norm = any_file_recursive(run_dir / "04_deconv_localwin_norm", "*.png")
        if has_main and has_deconv_norm:
            return True
    return False


def p6_condition_has_timescale(ds_root: Path, condition: str) -> bool:
    root = ds_root / "pipeline6_figures_timescale_diagnostics"
    tokens = P6_CONDITION_TOKENS.get(condition, (condition,))
    if not root.is_dir():
        return False
    return any(
        path.is_file()
        and path_contains_any_token(path, tokens)
        and not invalid_blp_path(path)
        for path in root.glob("*.png")
    )


def p6_condition_has_cross(ds_root: Path, stage: str, condition: str, pattern: str) -> bool:
    root = ds_root / stage
    tokens = P6_CONDITION_TOKENS.get(condition, (condition,))
    if not root.is_dir():
        return False
    return any(
        path.is_file()
        and path_contains_any_token(path, tokens)
        and not invalid_blp_path(path)
        for path in root.glob(pattern)
    )


def parse_method_component(method_tag: str) -> tuple[str, Optional[int]]:
    match = re.match(r"^([A-Za-z]+)_k(\d+)$", method_tag)
    if not match:
        return method_tag, None
    return match.group(1).lower(), int(match.group(2))


def discover_p10_candidates(
    ds_root: Path,
    dataset: str,
    p9_min_mtime_by_tag: Optional[dict[str, float]] = None,
) -> list[P10Candidate]:
    p9_root = ds_root / "pipeline9_bold_eigenfunction_reduction"
    candidates: list[P10Candidate] = []
    if not p9_root.is_dir():
        return candidates
    for run_tag in CURRENT_P7_RUN_TAGS:
        for feature_name in CURRENT_P10_FEATURES:
            feature_root = p9_root / run_tag / feature_name
            if not feature_root.is_dir():
                continue
            for method_dir in sorted(path for path in feature_root.iterdir() if path.is_dir()):
                method, k = parse_method_component(method_dir.name)
                if method not in CURRENT_P10_METHODS or k not in CURRENT_P10_COMPONENT_COUNTS:
                    continue
                mat_root = method_dir / "mat"
                result_file = latest_file(mat_root, "*.mat")
                if result_file is None:
                    continue
                if p9_min_mtime_by_tag:
                    min_mtime = p9_min_mtime_by_tag.get(run_tag)
                    if min_mtime is not None and result_file.stat().st_mtime + 1 < min_mtime:
                        continue
                p10_tag = f"{run_tag}__{feature_name}__{method_dir.name}"
                xcorr_dir = ds_root / "pipeline10_dimred_xcorr" / p10_tag
                top_maps_dir = ds_root / "pipeline10_dimred_top_maps" / p10_tag
                candidates.append(
                    P10Candidate(
                        dataset=dataset,
                        run_tag=run_tag,
                        feature_name=feature_name,
                        method_tag=method_dir.name,
                        component_count=int(k),
                        dimred_result_file=result_file,
                        p10_tag=p10_tag,
                        xcorr_dir=xcorr_dir,
                        xcorr_file=xcorr_dir / "dimred_xcorr.mat",
                        top_maps_dir=top_maps_dir,
                    )
                )
    return candidates


def xcorr_level_status(run_dir: Path, prefix: str) -> tuple[dict[str, bool], str]:
    combined = (
        (run_dir / f"{prefix}_summary.png").is_file()
        and (run_dir / f"{prefix}_top_signal_overlay.png").is_file()
    )
    density_dir = run_dir / "density"
    density_summary = list(density_dir.glob("*_summary.png")) if density_dir.is_dir() else []
    density_overlay = list(density_dir.glob("*_top_signal_overlay.png")) if density_dir.is_dir() else []
    per_density = bool(density_summary and density_overlay)
    feature_dir = run_dir / "feature"
    feature_summary = list(feature_dir.glob("*/*_summary.png")) if feature_dir.is_dir() else []
    feature_overlay = list(feature_dir.glob("*/*_top_signal_overlay.png")) if feature_dir.is_dir() else []
    per_density_feature = bool(feature_summary and feature_overlay)
    levels = {
        "combined": combined,
        "per_density": per_density,
        "per_density_feature": per_density_feature,
    }
    detail = (
        f"combined={int(combined)}; "
        f"density_summary={len(density_summary)}; density_overlay={len(density_overlay)}; "
        f"feature_summary={len(feature_summary)}; feature_overlay={len(feature_overlay)}"
    )
    return levels, detail


def p10_candidate_status(candidate: P10Candidate) -> tuple[str, str, int]:
    has_xcorr = candidate.xcorr_file.is_file()
    levels, level_detail = xcorr_level_status(candidate.xcorr_dir, "dimred_xcorr")
    has_top_maps = any_file_recursive(candidate.top_maps_dir / "fig", "**/*.png")
    if has_xcorr and candidate.xcorr_file.stat().st_mtime + 1 < candidate.dimred_result_file.stat().st_mtime:
        return (
            "stale",
            "dimred_xcorr.mat is older than the current P9 source; "
            f"{level_detail}",
            count_files(candidate.top_maps_dir),
        )
    if has_xcorr:
        return (
            "ok",
            "numeric xcorr present; browse figures are optional/default-off; "
            f"top_maps_present={int(has_top_maps)}; {level_detail}",
            count_files(candidate.top_maps_dir),
        )
    missing = []
    if not has_xcorr:
        missing.append("dimred_xcorr.mat")
    return "missing", "missing " + ", ".join(missing) + f"; {level_detail}", count_files(candidate.top_maps_dir)


def write_csv(path: Path, rows: Sequence[object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].__dataclass_fields__.keys())  # type: ignore[attr-defined]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row.__dict__)  # type: ignore[attr-defined]


def slug(text: str, max_len: int = 180) -> str:
    out = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    out = re.sub(r"_+", "_", out).strip("_")
    return out[:max_len] or "item"


def unique_path(path: Path) -> Path:
    if not path.exists():
        return path
    stem = path.stem
    suffix = path.suffix
    for idx in range(2, 10000):
        candidate = path.with_name(f"{stem}__dup{idx:03d}{suffix}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"Could not allocate unique path for {path}")


def copy_file(src: Path, dst: Path, overwrite: bool) -> Path:
    dst.parent.mkdir(parents=True, exist_ok=True)
    target = dst if overwrite else unique_path(dst)
    shutil.copy2(src, target)
    return target


def quarantine_path(processed_root: Path, dataset: str, source: Path, stamp: str) -> Path:
    return processed_root / "_legacy_quarantine" / stamp / dataset / source.name


def p4_roots(dataset: str, autodl_blp_root: Path, autodl_bold_root: Path) -> tuple[Path, Path]:
    return (
        autodl_blp_root / dataset / "mlp" / "outputs",
        autodl_bold_root / dataset / "mlp" / "outputs",
    )


def has_p4_blp(dataset: str, autodl_blp_root: Path, autodl_bold_root: Path) -> bool:
    root, _ = p4_roots(dataset, autodl_blp_root, autodl_bold_root)
    return root.is_dir() and any(root.iterdir())


def has_p4_bold(dataset: str, autodl_blp_root: Path, autodl_bold_root: Path) -> bool:
    _, root = p4_roots(dataset, autodl_blp_root, autodl_bold_root)
    return root.is_dir() and any(root.iterdir())


def safe_float(value: object) -> Optional[float]:
    if isinstance(value, bool) or value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def safe_int(value: object) -> Optional[int]:
    if isinstance(value, bool) or value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


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


def p4_candidate_matches_mode(path: Path, param_mode: Optional[str]) -> bool:
    return param_mode is None or p4_run_param_mode(path) == param_mode


def p4_condition_suffix(condition: str) -> str:
    condition_name = CURRENT_P5_CONDITION_PARTS.get(condition, (condition, ""))[0]
    return f"_projected_vlambda_{condition_name}".lower()


def p4_blp_candidate_runs(root: Path, condition: str, param_mode: Optional[str] = None) -> list[Path]:
    if not root.is_dir():
        return []
    suffix = p4_condition_suffix(condition)
    return sorted(
        (
            path
            for path in root.iterdir()
            if path.is_dir()
            and path.name.lower().endswith(suffix)
            and p4_candidate_matches_mode(path, param_mode)
        ),
        key=lambda path: path.name.lower(),
    )


def p4_bold_candidate_runs(root: Path, observable: str, param_mode: Optional[str] = None) -> list[Path]:
    if not root.is_dir():
        return []
    suffix = f"_projected_vlambda_{observable}".lower()
    matches = [
        path
        for path in root.iterdir()
        if path.is_dir()
        and path.name.lower().endswith(suffix)
        and p4_candidate_matches_mode(path, param_mode)
    ]
    if matches:
        return sorted(matches, key=lambda path: path.name.lower())
    return sorted(
        (
            path
            for path in root.iterdir()
            if path.is_dir()
            and observable.lower() in path.name.lower()
            and p4_candidate_matches_mode(path, param_mode)
        ),
        key=lambda path: path.name.lower(),
    )


def p4_training_state_paths(autodl_bold_root: Path, dataset: str, run_name: str) -> tuple[Path, ...]:
    checkpoint_root = autodl_bold_root / dataset / "mlp" / "checkpoints" / run_name
    return (
        checkpoint_root / "final" / "training_state.json",
        checkpoint_root / "best" / "training_state.json",
    )


def read_json_file(path: Path) -> Optional[dict]:
    if not path.is_file():
        return None
    try:
        with path.open("r", encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def best_metric_from_training_state(
    data: dict,
) -> tuple[Optional[float], Optional[int], Optional[int], str]:
    summary = data.get("best_outer_summary")
    if isinstance(summary, dict):
        best_val = safe_float(summary.get("best_val_metric"))
        if best_val is not None:
            best_epoch = safe_int(summary.get("best_outer_epoch"))
            completed = safe_int(data.get("completed_outer_epochs"))
            return best_val, best_epoch, completed, "best_outer_summary.best_val_metric"

    best_val = safe_float(data.get("best_val_metric"))
    if best_val is not None:
        best_epoch = safe_int(data.get("best_outer_epoch"))
        completed = safe_int(data.get("completed_outer_epochs"))
        return best_val, best_epoch, completed, "best_val_metric"

    history = data.get("outer_history")
    if isinstance(history, list):
        best_row: Optional[dict] = None
        best_history_val: Optional[float] = None
        for row in history:
            if not isinstance(row, dict):
                continue
            val = safe_float(row.get("val_metric"))
            if val is None:
                continue
            if best_history_val is None or val < best_history_val:
                best_history_val = val
                best_row = row
        if best_history_val is not None:
            best_epoch = safe_int(best_row.get("outer_epoch")) if best_row else None
            completed = safe_int(data.get("completed_outer_epochs"))
            return best_history_val, best_epoch, completed, "outer_history.val_metric"

    completed = safe_int(data.get("completed_outer_epochs"))
    return None, None, completed, "metric_missing"


def p4_selection_detail(selection: P4BoldRunSelection) -> str:
    bits = [selection.detail]
    if selection.best_val_metric is not None:
        bits.append(f"best_val_metric={selection.best_val_metric:.6g}")
    if selection.best_outer_epoch is not None:
        bits.append(f"best_outer_epoch={selection.best_outer_epoch}")
    if selection.completed_outer_epochs is not None:
        bits.append(f"completed_outer_epochs={selection.completed_outer_epochs}")
    if selection.training_state_path:
        bits.append(f"training_state={selection.training_state_path}")
    return "; ".join(bit for bit in bits if bit)


@lru_cache(maxsize=None)
def find_current_p4_bold_selection(
    dataset: str,
    observable: str,
    autodl_bold_root: Path,
    param_mode: Optional[str] = None,
) -> Optional[P4BoldRunSelection]:
    root = autodl_bold_root / dataset / "mlp" / "outputs"
    matches = p4_bold_candidate_runs(root, observable, param_mode)
    if not matches:
        return None

    selections: list[P4BoldRunSelection] = []
    for run_path in matches:
        state_path = ""
        best_val = None
        best_epoch = None
        completed = None
        metric_source = "metric_missing"
        for candidate_state_path in p4_training_state_paths(
            autodl_bold_root, dataset, run_path.name
        ):
            state = read_json_file(candidate_state_path)
            if state is None:
                continue
            state_path = str(candidate_state_path)
            best_val, best_epoch, completed, metric_source = best_metric_from_training_state(state)
            break

        status = "best_metric_candidate" if best_val is not None else "metric_missing_candidate"
        selections.append(
            P4BoldRunSelection(
                run_path=run_path,
                status=status,
                detail=(
                    f"{metric_source}; p4_param_mode={p4_run_param_mode(run_path)}; "
                    f"candidate_count={len(matches)}"
                ),
                best_val_metric=best_val,
                best_outer_epoch=best_epoch,
                completed_outer_epochs=completed,
                training_state_path=state_path,
            )
        )

    metric_candidates = [item for item in selections if item.best_val_metric is not None]
    if metric_candidates:
        return min(
            metric_candidates,
            key=lambda item: (
                item.best_val_metric if item.best_val_metric is not None else float("inf"),
                -(item.completed_outer_epochs or -1),
                item.run_path.name.lower(),
            ),
        )

    fallback = max(selections, key=lambda item: item.run_path.stat().st_mtime)
    return P4BoldRunSelection(
        run_path=fallback.run_path,
        status="metric_missing_latest_fallback",
        detail=(
            "metric missing for all candidates; latest mtime fallback; "
            f"p4_param_mode={p4_run_param_mode(fallback.run_path)}; candidate_count={len(matches)}"
        ),
        best_val_metric=None,
        best_outer_epoch=None,
        completed_outer_epochs=fallback.completed_outer_epochs,
        training_state_path=fallback.training_state_path,
    )


def find_current_p4_bold_run(
    dataset: str,
    observable: str,
    autodl_bold_root: Path,
    param_mode: Optional[str] = None,
) -> Optional[Path]:
    selection = find_current_p4_bold_selection(dataset, observable, autodl_bold_root, param_mode)
    if selection is None:
        return None
    return selection.run_path


def find_p7_run_for_p4(ds_root: Path, p4_run: Optional[Path]) -> Optional[Path]:
    if p4_run is None:
        return None
    candidate = ds_root / "pipeline7_bold_reskoopnet_postprocessing" / p4_run.name
    return candidate if candidate.is_dir() else None


def p8_run_status(ds_root: Path, run_tag: str, p7_run: Optional[Path]) -> tuple[str, str]:
    run_dir = ds_root / "pipeline8_xcorr" / run_tag
    if p7_run is None:
        return "blocked_missing_current_p7", "current P7 mainline for this observable is missing"
    if not run_dir.is_dir():
        return "missing", "missing P8 xcorr run folder"
    xcorr_file = run_dir / "xcorr.mat"
    if not xcorr_file.is_file():
        return "missing", "missing P8 xcorr.mat"
    return "ok", f"numeric xcorr present; browse figures/top maps are optional/default-off; requires current P7 source {p7_run.name}"


def p9_run_status(ds_root: Path, run_tag: str, p7_run: Optional[Path]) -> tuple[str, str]:
    run_dir = ds_root / "pipeline9_bold_eigenfunction_reduction" / run_tag
    if p7_run is None:
        return "blocked_missing_current_p7", "current P7 mainline for this observable is missing"
    if not run_dir.is_dir():
        return "missing", "missing P9 run folder"

    manifest_names = p9_manifest_run_names(ds_root, run_tag)
    if manifest_names and p7_run.name not in manifest_names:
        return (
            "stale",
            "P9 manifest run_name does not match current P7 best checkpoint; "
            f"expected {p7_run.name}; found {', '.join(sorted(manifest_names))}",
        )

    expected_mats = p9_current_efun_real_mats(run_dir)
    if len(expected_mats) < len(CURRENT_P10_METHODS) * len(CURRENT_P10_COMPONENT_COUNTS):
        return (
            "missing",
            "missing current efun_real P9 MAT grid "
            f"({len(expected_mats)}/{len(CURRENT_P10_METHODS) * len(CURRENT_P10_COMPONENT_COUNTS)})",
        )

    bold_post = p7_bold_post_file(p7_run)
    if bold_post is not None:
        newest_p9 = max(path.stat().st_mtime for path in expected_mats)
        if newest_p9 + 1 < bold_post.stat().st_mtime:
            return (
                "stale",
                "P9 efun_real MATs are older than the current P7 BOLD_POST; "
                f"expected source {p7_run.name}",
            )

    checked = p9_result_current_match_count(expected_mats, p7_run.name)
    if checked > 0 and checked < len(expected_mats):
        return (
            "stale",
            "some P9 efun_real MATs do not match the current P7 best checkpoint; "
            f"matched {checked}/{len(expected_mats)}",
        )

    return "ok", f"requires current P7 source {p7_run.name}"


def p7_bold_post_file(p7_run: Path) -> Optional[Path]:
    mat_dir = p7_run / "mat"
    if not mat_dir.is_dir():
        return None
    matches = sorted(mat_dir.glob("*_bold_post.mat"))
    if not matches:
        return None
    return max(matches, key=lambda path: path.stat().st_mtime)


def p9_min_mtime_by_tag(current_p7_by_run_tag: dict[str, Optional[Path]]) -> dict[str, float]:
    out: dict[str, float] = {}
    for run_tag, p7_run in current_p7_by_run_tag.items():
        if p7_run is None:
            continue
        bold_post = p7_bold_post_file(p7_run)
        if bold_post is not None:
            out[run_tag] = bold_post.stat().st_mtime
    return out


def p9_manifest_run_names(ds_root: Path, run_tag: str) -> set[str]:
    manifest = (
        ds_root
        / "pipeline9_bold_eigenfunction_reduction"
        / f"{ds_root.name}_pipeline9_bold_eigenfunction_reduction_manifest.csv"
    )
    if not manifest.is_file():
        return set()
    out: set[str] = set()
    try:
        with manifest.open("r", newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                if row.get("run_tag") != run_tag:
                    continue
                if row.get("feature_name") not in CURRENT_P10_FEATURES:
                    continue
                method, k = parse_method_component(row.get("method_tag", ""))
                if method not in CURRENT_P10_METHODS or k not in CURRENT_P10_COMPONENT_COUNTS:
                    continue
                run_name = row.get("run_name", "")
                if run_name:
                    out.add(run_name)
    except OSError:
        return set()
    return out


def p9_current_efun_real_mats(run_dir: Path) -> list[Path]:
    mats: list[Path] = []
    for feature_name in CURRENT_P10_FEATURES:
        feature_root = run_dir / feature_name
        if not feature_root.is_dir():
            continue
        for method in CURRENT_P10_METHODS:
            for k in CURRENT_P10_COMPONENT_COUNTS:
                mat_root = feature_root / f"{method}_k{k:02d}" / "mat"
                mat_file = latest_file(mat_root, "*.mat")
                if mat_file is not None:
                    mats.append(mat_file)
    return mats


def p9_result_current_match_count(result_files: Sequence[Path], expected_run_name: str) -> int:
    # P9 MAT files are MATLAB v7.3/HDF5 files. The base Python runtime used by
    # this audit intentionally avoids an HDF5 dependency, so provenance is
    # checked via the manifest when available and via MAT mtimes against the
    # current P7 BOLD_POST. Keep this hook to make the status logic explicit.
    return len(result_files)


def audit_stage_status(
    processed_root: Path,
    datasets: Sequence[str],
    autodl_blp_root: Path,
    autodl_bold_root: Path,
    repo_results_root: Path = DEFAULT_REPO_RESULTS_ROOT,
) -> list[StageStatus]:
    rows: list[StageStatus] = []
    p5_activity_metadata_mtime = newest_existing_mtime(P5_ACTIVITY_METADATA_SOURCE_FILES)
    for dataset in datasets:
        ds_root = processed_root / dataset
        rows.append(stage_row(dataset, "root", "dataset_root", ds_root))
        rows.append(
            StageStatus(
                dataset,
                "P4",
                "current_blp_p4_root",
                "ok" if has_p4_blp(dataset, autodl_blp_root, autodl_bold_root) else "missing",
                str(p4_roots(dataset, autodl_blp_root, autodl_bold_root)[0]),
                detail="user-managed; P11 does not backfill P4",
            )
        )
        rows.append(
            StageStatus(
                dataset,
                "P4",
                "current_bold_p4_root",
                "ok" if has_p4_bold(dataset, autodl_blp_root, autodl_bold_root) else "missing",
                str(p4_roots(dataset, autodl_blp_root, autodl_bold_root)[1]),
                detail="user-managed; P11 does not backfill P4",
            )
        )
        for condition in CURRENT_P5_CONDITIONS:
            blp_root = p4_roots(dataset, autodl_blp_root, autodl_bold_root)[0]
            for param_mode in CURRENT_P4_PARAM_MODES:
                candidates = p4_blp_candidate_runs(blp_root, condition, param_mode)
                invalid_detail = invalid_blp_condition_detail(dataset, condition)
                status = "ok" if candidates else "missing"
                if invalid_detail and candidates:
                    status = "stale"
                rows.append(
                    StageStatus(
                        dataset,
                        "P4",
                        f"blp_version:{condition}:{param_mode}",
                        status,
                        str(blp_root),
                        file_count=len(candidates),
                        dir_count=len(candidates),
                        detail=(
                            f"required P4 BLP version audit by condition and parameter mode; "
                            f"matched_runs={';'.join(path.name for path in candidates[:8])}"
                            + ("; ..." if len(candidates) > 8 else "")
                            + (f"; {invalid_detail}" if invalid_detail else "")
                        ),
                    )
                )
            for invalid_row in invalid_blp_rows_for_condition(dataset, condition):
                rows.append(
                    StageStatus(
                        dataset,
                        "P4",
                        f"invalid_blp_run:{condition}",
                        invalid_row.get("status", "invalid"),
                        invalid_row.get("output_dir", ""),
                        detail=invalid_row.get("reason", ""),
                    )
                )
        for stage in CANONICAL_DATASET_STAGES:
            pipeline = stage.split("_", 1)[0].upper()
            ok_status = "ok"
            missing_status = "missing"
            detail = ""
            if stage in LEGACY_OPTIONAL_DATASET_STAGES:
                ok_status = "legacy_optional_present"
                missing_status = "legacy_optional_missing"
                detail = "legacy branch; ignored by current mainline completeness"
            elif stage in MINOR_OPTIONAL_DATASET_STAGES:
                ok_status = "minor_ok"
                missing_status = "minor_missing"
                detail = "minor/auxiliary branch; ignored by current mainline completeness"
            elif stage in DERIVED_CACHE_DATASET_STAGES:
                ok_status = "cache_present"
                missing_status = "cache_missing"
                detail = "derived browsing/cache stage; P11 rebuilds flat summaries from canonical sources"
            rows.append(
                stage_row(
                    dataset,
                    pipeline,
                    stage,
                    ds_root / stage,
                    detail=detail,
                    ok_status=ok_status,
                    missing_status=missing_status,
                )
            )
        for stage in NONCANONICAL_DATASET_STAGES:
            pipeline = stage.split("_", 1)[0].upper()
            path = ds_root / stage
            if path.exists():
                status = "layout_debt" if stage in LAYOUT_DEBT_DATASET_STAGES else "legacy_candidate"
                rows.append(
                    StageStatus(
                        dataset,
                        pipeline,
                        f"noncanonical:{stage}",
                        status,
                        str(path),
                        file_count=count_files(path),
                        dir_count=count_dirs(path),
                        detail="present on disk but not a current canonical stage",
                    )
                )

        # P3 mode-level QC coverage.
        p3_qc = ds_root / "pipeline3_figures_bold_pre_reskoopnet_qc"
        for obs in CURRENT_P3_OBSERVABLES:
            qc = latest_file(p3_qc / obs, "*.png") or latest_file(p3_qc, f"*{obs}*.png")
            rows.append(
                StageStatus(
                    dataset,
                    "P3",
                    f"qc_png:{obs}",
                    "ok" if qc else "missing",
                    "" if qc is None else str(qc),
                    detail="current four-observable P3 QC set",
                )
            )

        # P7/P8/P9 current run-tag coverage. P7 is resolved through the
        # best-validation BOLD P4 checkpoint first; downstream P8/P9 folders
        # are not enough to count as current when P7 is stale or absent.
        p7_root = ds_root / "pipeline7_bold_reskoopnet_postprocessing"
        p7_mainline_by_run_tag: dict[str, Optional[Path]] = {}
        p7_mainline_by_mode_run_tag: dict[str, dict[str, Optional[Path]]] = {
            param_mode: {} for param_mode in CURRENT_P4_PARAM_MODES
        }
        for obs in CURRENT_P7_OBSERVABLES:
            run_tag = CURRENT_P7_OBSERVABLE_TO_RUN_TAG[obs]
            selections_by_mode: dict[str, Optional[P4BoldRunSelection]] = {}
            for param_mode in CURRENT_P4_PARAM_MODES:
                mode_selection = find_current_p4_bold_selection(
                    dataset,
                    obs,
                    autodl_bold_root,
                    param_mode,
                )
                selections_by_mode[param_mode] = mode_selection
                mode_p4_run = mode_selection.run_path if mode_selection is not None else None
                mode_p7_run = find_p7_run_for_p4(ds_root, mode_p4_run)
                p7_mainline_by_mode_run_tag[param_mode][run_tag] = mode_p7_run
                rows.append(
                    StageStatus(
                        dataset,
                        "P4",
                        f"bold_version:{obs}:{param_mode}",
                        "ok" if mode_p4_run else "missing",
                        "" if mode_p4_run is None else str(mode_p4_run),
                        file_count=count_files(mode_p4_run) if mode_p4_run else 0,
                        dir_count=count_dirs(mode_p4_run) if mode_p4_run else 0,
                        detail=(
                            f"required P4 BOLD version audit by observable and parameter mode={param_mode}; "
                            + (
                                "missing matching P4 source in current BOLD AutoDL root"
                                if mode_selection is None
                                else p4_selection_detail(mode_selection)
                            )
                        ),
                    )
                )
                rows.append(
                    StageStatus(
                        dataset,
                        "P7",
                        f"current_mainline_run:{param_mode}:{obs}",
                        "ok" if mode_p7_run else "missing",
                        "" if mode_p7_run is None else str(mode_p7_run),
                        file_count=count_files(mode_p7_run) if mode_p7_run else 0,
                        dir_count=count_dirs(mode_p7_run) if mode_p7_run else 0,
                        detail=(
                            f"P7 best-checkpoint audit separated by P4 parameter mode={param_mode}; "
                            + (
                                "blocked by missing matching P4 source"
                                if mode_selection is None
                                else f"expected P7 folder name: {mode_selection.run_path.name}; "
                                + p4_selection_detail(mode_selection)
                            )
                        ),
                    )
                )

            p4_selection = selections_by_mode.get(CURRENT_P4_DOWNSTREAM_PARAM_MODE)
            if p4_selection is None:
                p4_selection = find_current_p4_bold_selection(dataset, obs, autodl_bold_root)
            p4_run = p4_selection.run_path if p4_selection is not None else None
            p7_run = find_p7_run_for_p4(ds_root, p4_run)
            p7_mainline_by_run_tag[run_tag] = p7_run
            fallback_p7_run = find_p7_observable_run(p7_root, obs)
            rows.append(
                StageStatus(
                    dataset,
                    "P4",
                    f"current_bold_run:{obs}",
                    "ok" if p4_run else "missing",
                    "" if p4_run is None else str(p4_run),
                    file_count=count_files(p4_run) if p4_run else 0,
                    dir_count=count_dirs(p4_run) if p4_run else 0,
                    detail=(
                        "missing current BOLD P4 source in mainline autodl root"
                        if p4_selection is None
                        else f"downstream default p4_param_mode={CURRENT_P4_DOWNSTREAM_PARAM_MODE}; "
                        "best checkpoint selected by validation metric; "
                        + p4_selection_detail(p4_selection)
                    ),
                )
            )
            detail = "" if p4_run is None else f"expected P7 folder name: {p4_run.name}"
            if p4_selection is not None:
                detail += "; " + p4_selection_detail(p4_selection)
            if p7_run is None and fallback_p7_run is not None:
                detail += f"; older/different P7 exists: {fallback_p7_run.name}"
            rows.append(
                StageStatus(
                    dataset,
                    "P7",
                    f"current_mainline_run:{obs}",
                    "ok" if p7_run else "missing",
                    "" if p7_run is None else str(p7_run),
                    file_count=count_files(p7_run) if p7_run else 0,
                    dir_count=count_dirs(p7_run) if p7_run else 0,
                    detail=detail,
                )
            )
        for run_tag in CURRENT_P7_RUN_TAGS:
            p7_run = p7_mainline_by_run_tag.get(run_tag)
            p8_status, p8_detail = p8_run_status(ds_root, run_tag, p7_run)
            p9_status, p9_detail = p9_run_status(ds_root, run_tag, p7_run)
            p8_xcorr = ds_root / "pipeline8_xcorr" / run_tag
            p8_top = ds_root / "pipeline8_top_maps" / run_tag
            p9_run = ds_root / "pipeline9_bold_eigenfunction_reduction" / run_tag
            rows.append(
                StageStatus(
                    dataset,
                    "P8",
                    f"current_xcorr:{run_tag}",
                    p8_status if p8_xcorr.is_dir() else ("blocked_missing_current_p7" if p7_run is None else "missing"),
                    str(p8_xcorr),
                    file_count=count_files(p8_xcorr),
                    dir_count=count_dirs(p8_xcorr),
                    detail=p8_detail,
                )
            )
            rows.append(
                StageStatus(
                    dataset,
                    "P8",
                    f"current_top_maps:{run_tag}",
                    "blocked_missing_current_p7" if p7_run is None else (
                        "optional_present" if p8_top.is_dir() else "optional_missing"
                    ),
                    str(p8_top),
                    file_count=count_files(p8_top),
                    dir_count=count_dirs(p8_top),
                    detail="browse/top-map output is optional/default-off for numeric-first P8",
                )
            )
            rows.append(
                StageStatus(
                    dataset,
                    "P9",
                    f"current_run:{run_tag}",
                    p9_status,
                    str(p9_run),
                    file_count=count_files(p9_run),
                    dir_count=count_dirs(p9_run),
                    detail=p9_detail,
                )
            )
            if p8_xcorr.is_dir():
                levels, level_detail = xcorr_level_status(p8_xcorr, "xcorr")
                for level, ok in levels.items():
                    rows.append(
                        StageStatus(
                            dataset=dataset,
                            pipeline="P8",
                            item=f"xcorr_overview_level:{run_tag}:{level}",
                            status="ok" if ok and p7_run is not None else (
                                "blocked_missing_current_p7" if p7_run is None else "optional_missing"
                            ),
                            path=str(p8_xcorr),
                            file_count=count_files(p8_xcorr),
                            dir_count=count_dirs(p8_xcorr),
                            detail=level_detail + "; browse xcorr figures are optional/default-off",
                        )
                    )

        # P5 current grid quick source check.
        for condition in CURRENT_P5_CONDITIONS:
            invalid_detail = invalid_blp_condition_detail(dataset, condition)
            raw_density = ds_root / "pipeline5_raw_thresholded_density" / condition
            raw_density_mat = latest_file(raw_density / "mat", "*.mat")
            raw_status, raw_detail = current_artifact_status(
                raw_density_mat, p5_activity_metadata_mtime
            )
            if invalid_detail:
                rows.append(
                    StageStatus(
                        dataset,
                        "P5",
                        f"raw_activity_density_metadata:{condition}",
                        "stale",
                        str(raw_density / "mat"),
                        file_count=count_files(raw_density),
                        dir_count=count_dirs(raw_density),
                        detail=invalid_detail,
                    )
                )
            else:
                rows.append(
                    StageStatus(
                        dataset,
                        "P5",
                        f"raw_activity_density_metadata:{condition}",
                        raw_status,
                        "" if raw_density_mat is None else str(raw_density_mat),
                        file_count=count_files(raw_density / "mat"),
                        dir_count=count_dirs(raw_density / "mat"),
                        detail=(
                            raw_detail
                            + "; requires abs_magnitude metadata with raw_efun_index/timescale"
                        ),
                    )
                )
            for method in CURRENT_P5_METHODS:
                for k in CURRENT_P5_COMPONENT_COUNTS:
                    tag = f"{method}_k{k:02d}"
                    mat_root = ds_root / "pipeline5_eigenfunction_reduction" / condition / tag / "mat"
                    has_mat = any_file(mat_root, "*.mat")
                    rows.append(
                        StageStatus(
                            dataset,
                            "P5",
                            f"reduction_mat:{condition}:{tag}",
                            "stale" if invalid_detail else ("ok" if has_mat else "missing"),
                            str(mat_root),
                            detail=invalid_detail,
                        )
                    )
                    dimred_density_root = (
                        ds_root
                        / "pipeline5_dimred_thresholded_density"
                        / condition
                        / tag
                        / "mat"
                    )
                    dimred_density_mat = latest_file(dimred_density_root, "*.mat")
                    dimred_status, dimred_detail = current_artifact_status(
                        dimred_density_mat, p5_activity_metadata_mtime
                    )
                    rows.append(
                        StageStatus(
                            dataset,
                            "P5",
                            f"dimred_activity_density_metadata:{condition}:{tag}",
                            "stale" if invalid_detail else dimred_status,
                            "" if dimred_density_mat is None else str(dimred_density_mat),
                            file_count=count_files(dimred_density_root),
                            dir_count=count_dirs(dimred_density_root),
                            detail=invalid_detail
                            or (
                                dimred_detail
                                + "; requires abs_magnitude metadata with component weighted timescale"
                            ),
                        )
                    )
                    stats_file = (
                        ds_root
                        / "pipeline5_eigenfunction_peaks_by_state_maxabs"
                        / condition
                        / tag
                        / f"{dataset}_{condition}_{tag}_peaks_stats.csv"
                    )
                    rows.append(
                        StageStatus(
                            dataset,
                            "P5",
                            f"peak_stats_maxabs:{condition}:{tag}",
                            "stale" if invalid_detail else ("ok" if stats_file.is_file() else "missing"),
                            str(stats_file),
                            detail=invalid_detail
                            or "current source for P5 event-family selectivity and selected-component trajectory",
                        )
                    )
                    selective_manifest = (
                        ds_root
                        / "pipeline5_eigenfunction_reduction"
                        / condition
                        / tag
                        / "fig"
                        / "selected_component_trajectory_manifest.csv"
                    )
                    rows.append(
                        StageStatus(
                            dataset,
                    "P5",
                    f"selective_dimred_3d_consensus_trajectory:{condition}:{tag}",
                    "stale"
                    if invalid_detail
                    else ("ok" if selective_manifest.is_file() else "deferred_implementation"),
                    str(selective_manifest),
                    detail=invalid_detail
                    or (
                        "P5 trajectory QC should use selective dimred components "
                        "in 3D consensus-state space; "
                        "deferred/non-blocking for current P11 completeness; "
                        "existing first-three-dimension scatter figures are diagnostic only"
                    ),
                )
                    )

        # P6 current condition-level coverage. Folder-level presence can hide
        # missing abs/complex_split branches, so report each required condition.
        for condition in CURRENT_P5_CONDITIONS:
            token = p6_condition_token(condition)
            invalid_detail = invalid_blp_condition_detail(dataset, condition)
            checks = (
                (
                    "top_window_main_and_localwin_deconv",
                    p6_condition_has_top_windows(ds_root, condition),
                    ds_root / "pipeline6_top_state_diversity_postprocessing",
                ),
                (
                    "timescale_diagnostics",
                    p6_condition_has_timescale(ds_root, condition),
                    ds_root / "pipeline6_figures_timescale_diagnostics",
                ),
                (
                    "spkt_residual_xcorr",
                    p6_condition_has_cross(
                        ds_root,
                        "pipeline6_spkt_residual_cross_correlation",
                        condition,
                        "**/*.mat",
                    ),
                    ds_root / "pipeline6_spkt_residual_cross_correlation",
                ),
                (
                    "mua_residual_xcorr",
                    p6_condition_has_cross(
                        ds_root,
                        "pipeline6_mua_residual_cross_correlation",
                        condition,
                        "**/*.mat",
                    ),
                    ds_root / "pipeline6_mua_residual_cross_correlation",
                ),
            )
            for check_name, ok, path in checks:
                rows.append(
                    StageStatus(
                        dataset=dataset,
                        pipeline="P6",
                        item=f"condition:{condition}:{check_name}",
                        status="ok" if ok else ("stale" if invalid_detail else "missing"),
                        path=str(path),
                        file_count=count_files(path),
                        dir_count=count_dirs(path),
                        detail=f"condition token: {token}; {invalid_detail}".rstrip("; "),
                    )
                )

        # P10 candidate-level coverage derived from current P9 efun_real
        # results. This prevents partial folders from looking complete.
        p10_candidates = discover_p10_candidates(
            ds_root,
            dataset,
            p9_min_mtime_by_tag(p7_mainline_by_run_tag),
        )
        rows.append(
            StageStatus(
                dataset=dataset,
                pipeline="P10",
                item="candidate_discovery:efun_real",
                status="ok" if p10_candidates else "missing",
                path=str(ds_root / "pipeline9_bold_eigenfunction_reduction"),
                file_count=len(p10_candidates),
                detail="current P10 candidates discovered from P9 result MAT files",
            )
        )
        for candidate in p10_candidates:
            status, detail, file_count = p10_candidate_status(candidate)
            rows.append(
                StageStatus(
                    dataset=dataset,
                    pipeline="P10",
                    item=f"candidate:{candidate.p10_tag}",
                    status=status,
                    path=str(candidate.xcorr_dir),
                    file_count=file_count,
                    dir_count=count_dirs(candidate.top_maps_dir),
                    detail=f"{detail}; source={candidate.dimred_result_file}",
                )
            )
    # Repo-local P11 derived products that summarize the new interpretation
    # layer. These are not dataset-local canonical stages, but P11 should audit
    # them because they drive parameter selection and recompute planning.
    label_table = (
        repo_results_root
        / "pipeline5_dimred_component_process_labels_current"
        / "dimred_efun_process_labels.csv"
    )
    label_status = "missing"
    label_detail = "run scripts/build_pipeline5_dimred_component_process_labels.py"
    if label_table.is_file():
        label_status = "transitional_ok"
        label_detail = "current table exists; maxabs-derived labels remain transitional until activity/envelope P5 is regenerated"
    rows.append(
        StageStatus(
            "all",
            "P5",
            "dimred_component_process_labels",
            label_status,
            str(label_table),
            file_count=1 if label_table.is_file() else 0,
            detail=label_detail,
        )
    )
    score_root = repo_results_root / "pipeline11_parameter_selection_scorecard_current"
    derived_products = (
        ("P11", "p1_p4_readiness_long", repo_results_root / "pipeline11_p1_p4_readiness_current" / "p1_p4_readiness_long.csv"),
        ("P11", "p1_p4_recompute_requirements", repo_results_root / "pipeline11_p1_p4_readiness_current" / "p1_p4_recompute_requirements_by_dataset.csv"),
        ("P11", "parameter_selection_scorecard", score_root / "parameter_selection_scorecard.csv"),
        ("P11", "recompute_requirements_by_dataset", score_root / "recompute_requirements_by_dataset.csv"),
        ("P11", "raw_efun_index_distribution", score_root / "raw_efun_index_distribution.csv"),
        ("P11", "top_hit_interpretation_summary", score_root / "top_hit_interpretation_summary.csv"),
        ("P11", "timescale_and_label_metadata_audit", score_root / "timescale_and_label_metadata_audit.csv"),
        ("P11", "parameter_selection_summary", score_root / "summary.md"),
    )
    for pipeline, item, path in derived_products:
        rows.append(
            StageStatus(
                "all",
                pipeline,
                item,
                "ok" if path.is_file() else "missing",
                str(path),
                file_count=1 if path.is_file() else 0,
                detail="repo-local derived P11 parameter-selection/recompute product",
            )
        )
    fig_root = score_root / "figures"
    rows.append(
        StageStatus(
            "all",
            "P11",
            "parameter_selection_scorecard_figures",
            "ok" if any_file(fig_root, "*.png") else "missing",
            str(fig_root),
            file_count=count_files(fig_root),
            dir_count=count_dirs(fig_root),
            detail="default scorecard figures for condition/method/k/BOLD-observable selection",
        )
    )
    return rows


def build_legacy_candidates(processed_root: Path, datasets: Sequence[str], stamp: str) -> list[LegacyCandidate]:
    rows: list[LegacyCandidate] = []
    for dataset in datasets:
        ds_root = processed_root / dataset
        if not ds_root.is_dir():
            continue
        for stage in QUARANTINE_CANDIDATE_STAGES:
            src = ds_root / stage
            if not src.exists():
                continue
            reason = "old_pipeline8_noncanonical_folder"
            rows.append(
                LegacyCandidate(
                    dataset=dataset,
                    source_path=str(src),
                    proposed_legacy_path=str(quarantine_path(processed_root, dataset, src, stamp)),
                    reason=reason,
                    status="planned_only",
                )
            )
    return rows


def add_action(
    rows: list[FillAction],
    dataset: str,
    pipeline: str,
    action: str,
    priority: int,
    memory_class: str,
    parallel_group: str,
    can_parallelize: bool,
    blocked_by: str,
    writes_external: bool,
    command_hint: str,
    reason: str,
) -> None:
    action_id = f"{pipeline}_{dataset}_{action}_{len(rows) + 1:04d}"
    rows.append(
        FillAction(
            action_id=action_id,
            dataset=dataset,
            pipeline=pipeline,
            action=action,
            priority=priority,
            memory_class=memory_class,
            parallel_group=parallel_group,
            can_parallelize="yes" if can_parallelize else "no",
            blocked_by=blocked_by,
            writes_external="yes" if writes_external else "no",
            command_hint=command_hint,
            reason=reason,
        )
    )


def build_fill_plan(
    processed_root: Path,
    datasets: Sequence[str],
    autodl_bold_root: Path,
) -> list[FillAction]:
    rows: list[FillAction] = []
    p5_activity_metadata_mtime = newest_existing_mtime(P5_ACTIVITY_METADATA_SOURCE_FILES)
    for dataset in datasets:
        ds_root = processed_root / dataset
        p2_required = (
            "pipeline2_event_detection",
            "pipeline2_event_density",
            "pipeline2_consensus_states",
            "pipeline2_consensus_state_summary",
            "pipeline2_consensus_state_diversity_windows",
        )
        if any(not (ds_root / stage).exists() for stage in p2_required):
            add_action(
                rows,
                dataset,
                "P2",
                "run_event_density_consensus_state_pipeline",
                20,
                "medium",
                "p2_cpu",
                True,
                "P1",
                True,
                "MATLAB: run current P2 event density/consensus-state scripts for this cfg",
                "missing one or more current P2 canonical stages",
            )

        p3_qc = ds_root / "pipeline3_figures_bold_pre_reskoopnet_qc"
        if any(
            (latest_file(p3_qc / obs, "*.png") or latest_file(p3_qc, f"*{obs}*.png")) is None
            for obs in CURRENT_P3_OBSERVABLES
        ):
            add_action(
                rows,
                dataset,
                "P3",
                "generate_current_observable_qc_figures",
                15,
                "low",
                "figure_cpu",
                True,
                "P3 observables",
                True,
                "MATLAB: script_plot_bold_pre_reskoopnet_qc.m with current four observables",
                "current P3 QC PNGs missing or incomplete",
            )

        p5_missing = False
        for condition in CURRENT_P5_CONDITIONS:
            if invalid_blp_condition_detail(dataset, condition):
                p5_missing = True
                break
            raw_density_mat = latest_file(
                ds_root / "pipeline5_raw_thresholded_density" / condition / "mat",
                "*.mat",
            )
            if current_artifact_status(raw_density_mat, p5_activity_metadata_mtime)[0] != "ok":
                p5_missing = True
                break
            for method in CURRENT_P5_METHODS:
                for k in CURRENT_P5_COMPONENT_COUNTS:
                    tag = f"{method}_k{k:02d}"
                    mat_root = ds_root / "pipeline5_eigenfunction_reduction" / condition / tag / "mat"
                    if not any_file(mat_root, "*.mat"):
                        p5_missing = True
                        break
                    dimred_density_mat = latest_file(
                        ds_root
                        / "pipeline5_dimred_thresholded_density"
                        / condition
                        / tag
                        / "mat",
                        "*.mat",
                    )
                    if current_artifact_status(dimred_density_mat, p5_activity_metadata_mtime)[0] != "ok":
                        p5_missing = True
                        break
                if p5_missing:
                    break
            if p5_missing:
                break
        if p5_missing:
            add_action(
                rows,
                dataset,
                "P5",
                "run_current_p5_reduction_grid",
                30,
                "high",
                "matlab_heavy_blp",
                False,
                "P4 BLP current root",
                True,
                "MATLAB: script_run_current_p5_activity_density_unblocked.m or per-dataset equivalent",
                "missing/stale current P5 reduction or activity-magnitude density metadata grid",
            )
        peak_root = ds_root / "pipeline5_eigenfunction_peaks_by_state_maxabs"
        if not peak_root.exists():
            add_action(
                rows,
                dataset,
                "P5",
                "run_peak_state_maxabs_statistics",
                35,
                "medium",
                "matlab_medium",
                True,
                "P5 reduction grid",
                True,
                "MATLAB: script_run_current_pipeline5_peak_state_all_datasets.m with dataset filter",
                "missing current maxabs peak-state statistics",
            )
        selected_trajectory_missing = True
        for condition in CURRENT_P5_CONDITIONS:
            for method in CURRENT_P5_METHODS:
                for k in CURRENT_P5_COMPONENT_COUNTS:
                    tag = f"{method}_k{k:02d}"
                    manifest = (
                        ds_root
                        / "pipeline5_eigenfunction_reduction"
                        / condition
                        / tag
                        / "fig"
                        / "selected_component_trajectory_manifest.csv"
                    )
                    if manifest.is_file():
                        selected_trajectory_missing = False
                        break
                if not selected_trajectory_missing:
                    break
            if not selected_trajectory_missing:
                break
        # Selected-component trajectory QC is scientifically preferred, but it
        # is intentionally deferred and does not block current P11 backfill.
        # Its missing status is recorded in p11_stage_status_latest.csv only.
        add_action(
            rows,
            dataset,
            "P5",
            "refresh_p5_flat_figures_and_selectivity",
            40,
            "low",
            "figure_cpu",
            True,
            "P5 canonical MAT/CSV",
            True,
            "Python/MATLAB: refresh current P5 scatter, peak-stat summary, consistency, and selectivity figures",
            "P11 should regenerate derived P5 browsing from canonical outputs every run",
        )

        p6_required = (
            "pipeline6_top_state_diversity_postprocessing",
            "pipeline6_figures_timescale_diagnostics",
            "pipeline6_figures_spkt_residual_cross_correlation",
            "pipeline6_figures_mua_residual_cross_correlation",
        )
        p6_missing = any(not (ds_root / stage).exists() for stage in p6_required)
        for condition in CURRENT_P5_CONDITIONS:
            condition_complete = (
                p6_condition_has_top_windows(ds_root, condition)
                and p6_condition_has_timescale(ds_root, condition)
                and p6_condition_has_cross(
                    ds_root,
                    "pipeline6_spkt_residual_cross_correlation",
                    condition,
                    "**/*.mat",
                )
                and p6_condition_has_cross(
                    ds_root,
                    "pipeline6_mua_residual_cross_correlation",
                    condition,
                    "**/*.mat",
                )
            )
            if not condition_complete:
                p6_missing = True
                break
        if p6_missing:
            add_action(
                rows,
                dataset,
                "P6",
                "run_p6_top_state_diversity_and_residual_xcorr",
                50,
                "high",
                "matlab_heavy_blp",
                False,
                "P5/P4 BLP",
                True,
                "MATLAB: script_run_completed_mlp_added_postprocessing_stages.m or current per-dataset P6 entry",
                "missing one or more current P6 stages",
            )

        p7_root = ds_root / "pipeline7_bold_reskoopnet_postprocessing"
        p7_missing_current = False
        current_p7_by_run_tag: dict[str, Optional[Path]] = {}
        for obs in CURRENT_P7_OBSERVABLES:
            run_tag = CURRENT_P7_OBSERVABLE_TO_RUN_TAG[obs]
            p4_run = find_current_p4_bold_run(dataset, obs, autodl_bold_root)
            p7_run = find_p7_run_for_p4(ds_root, p4_run)
            current_p7_by_run_tag[run_tag] = p7_run
            if p4_run is not None and p7_run is None:
                p7_missing_current = True
        if not p7_root.exists() or p7_missing_current:
            add_action(
                rows,
                dataset,
                "P7",
                "run_p7_bold_postprocessing",
                45,
                "medium",
                "matlab_medium_bold",
                True,
                "P4 BOLD current root",
                True,
                "MATLAB: script_run_one_cfg_bold_reskoopnet_postprocessing.m for current best checkpoints",
                "missing P7 canonical postprocessing",
            )
        add_action(
            rows,
            dataset,
            "P7",
            "refresh_p7_flat_figures",
            25,
            "low",
            "figure_cpu",
            True,
            "P7 canonical fig folder",
            True,
            "P11 Python: copy current P7 main/deconv/timescale/intrinsic figures with manifest",
            "P7 summary cache should be regenerated by P11",
        )

        p8_missing = not (ds_root / "pipeline8_xcorr").exists()
        if not p8_missing:
            for run_tag, p7_run in current_p7_by_run_tag.items():
                p8_status, _p8_detail = p8_run_status(ds_root, run_tag, p7_run)
                if p8_status != "ok":
                    p8_missing = True
                    break
        if p8_missing:
            add_action(
                rows,
                dataset,
                "P8",
                "run_p8_bold_lfp_density_coupling",
                60,
                "high",
                "matlab_heavy_bold",
                False,
                "P7 + P2 + P5 density",
                True,
                "MATLAB: script_run_one_cfg_bold_cross_modal_coupling.m for current four observables",
                "missing P8 numeric xcorr outputs",
            )
        add_action(
            rows,
            dataset,
            "P8",
            "refresh_p8_flat_figures_and_consistency_sources",
            42,
            "low",
            "figure_cpu",
            True,
            "P8 canonical xcorr CSV/MAT and optional figure folders",
            True,
            "P11 Python: rebuild compact cross-session/cross-method summaries from P8 numeric outputs",
            "P8 derived summary figures need refresh",
        )

        p9_root = ds_root / "pipeline9_bold_eigenfunction_reduction"
        p9_incomplete_current = not p9_root.exists()
        if not p9_incomplete_current:
            for run_tag, p7_run in current_p7_by_run_tag.items():
                p9_status, _p9_detail = p9_run_status(ds_root, run_tag, p7_run)
                if p9_status != "ok":
                    p9_incomplete_current = True
                    break
        if p9_incomplete_current:
            add_action(
                rows,
                dataset,
                "P9",
                "run_p9_bold_eigenfunction_reduction",
                55,
                "medium",
                "matlab_medium_bold",
                True,
                "P7 current BOLD_POST",
                True,
                "MATLAB: script_run_one_cfg_bold_eigenfunction_reduction.m",
                "missing P9 BOLD dimred outputs",
            )
        add_action(
            rows,
            dataset,
            "P9",
            "normalize_and_flatten_p9_figures",
            32,
            "low",
            "figure_cpu",
            True,
            "P9 canonical MAT/fig folders",
            True,
            "P11 Python: copy co-located P9 fig PNGs and mark layout debt",
            "P9 figures are currently co-located with MAT files",
        )

        p10_missing = not (ds_root / "pipeline10_dimred_xcorr").exists()
        if p9_incomplete_current:
            p10_missing = True
        if not p10_missing:
            for candidate in discover_p10_candidates(
                ds_root,
                dataset,
                p9_min_mtime_by_tag(current_p7_by_run_tag),
            ):
                status, _detail, _file_count = p10_candidate_status(candidate)
                if status != "ok":
                    p10_missing = True
                    break
        if p10_missing:
            add_action(
                rows,
                dataset,
                "P10",
                "run_p10_bold_dimred_density_coupling",
                70,
                "high",
                "matlab_heavy_bold",
                False,
                "P9 + P2 + P5 density",
                True,
                "MATLAB: script_run_one_cfg_bold_dimred_cross_modal_coupling.m",
                "missing P10 canonical outputs",
            )
        add_action(
            rows,
            dataset,
            "P10",
            "refresh_p10_flat_figures_and_consistency_sources",
            43,
            "low",
            "figure_cpu",
            True,
            "P10 canonical xcorr CSV/MAT and optional figure folders",
            True,
            "P11 Python: rebuild compact cross-session/cross-method summaries from P10 numeric outputs",
            "P10 derived summary figures need refresh",
        )
    label_table = (
        DEFAULT_REPO_RESULTS_ROOT
        / "pipeline5_dimred_component_process_labels_current"
        / "dimred_efun_process_labels.csv"
    )
    if not label_table.is_file():
        add_action(
            rows,
            "all",
            "P5",
            "build_dimred_component_process_labels",
            41,
            "low",
            "python_summary",
            True,
            "P5 event-family selectivity table",
            False,
            "Python: scripts/build_pipeline5_dimred_component_process_labels.py",
            "missing P5-to-P8/P10 dimred process-label bridge",
        )
    add_action(
        rows,
        "all",
        "P11",
        "build_p1_p4_readiness_audit",
        18,
        "low",
        "python_summary",
        True,
        "P1-P4 canonical products and AutoDL roots",
        False,
        "Python: scripts/audit_p1_p4_readiness.py --datasets <all9/current>",
        "P11 should explicitly audit P1-P4 readiness instead of relying on folder-level checks",
    )
    add_action(
        rows,
        "all",
        "P11",
        "build_parameter_selection_scorecard_and_recompute_audit",
        44,
        "low",
        "python_summary",
        True,
        "P5 labels + P8/P10 xcorr CSVs",
        False,
        "Python: scripts/summarize_p11_parameter_selection_scorecard.py --dataset-scope/all datasets",
        "P11 should refresh parameter-selection scorecard, top raw-efun index/timescale summaries, dimred process-label summaries split by BOLD efun/deconv_efun, timescale metadata audit, and recompute requirements every run",
    )
    return sorted(rows, key=lambda row: (row.priority, row.pipeline, row.dataset, row.action))


def p11_flat_name(dataset: str, pipeline: str, group: str, src: Path, anchor: Path) -> str:
    rel = safe_rel(src, anchor)
    digest = hashlib.sha1(rel.encode("utf-8", errors="ignore")).hexdigest()[:10]
    rel_no_ext = rel
    suffix = src.suffix.lower()
    if suffix and rel.lower().endswith(suffix):
        rel_no_ext = rel[: -len(suffix)]
    return f"{dataset}__{pipeline.lower()}__{group}__{slug(rel_no_ext)}__h{digest}{suffix}"


def collect_flat_sources(
    processed_root: Path,
    repo_results_root: Path,
    summary_root: Path,
    datasets: Sequence[str],
    autodl_bold_root: Path,
) -> list[FlatFigureSource]:
    rows: list[FlatFigureSource] = []
    seen: set[tuple[str, str, str]] = set()

    def add_source(dataset: str, pipeline: str, group: str, root: Path, src: Path, detail: str = "") -> None:
        if not src.is_file():
            return
        if pipeline == "P6" and invalid_blp_path(src):
            return
        key = (pipeline, group, str(src).lower())
        if key in seen:
            return
        seen.add(key)
        rel = safe_rel(src, root)
        dst_name = p11_flat_name(dataset, pipeline, group, src, root)
        rows.append(
            FlatFigureSource(
                dataset=dataset,
                pipeline=pipeline,
                figure_group=group,
                source_path=str(src),
                relative_source=rel,
                proposed_flat_path=str(summary_root / pipeline.lower() / group / dst_name),
                status="planned",
                detail=detail,
            )
        )

    for dataset in datasets:
        ds_root = processed_root / dataset
        current_p7_by_run_tag: dict[str, Optional[Path]] = {}
        valid_p9_run_tags: set[str] = set()
        for obs, run_tag in CURRENT_P7_OBSERVABLE_TO_RUN_TAG.items():
            p4_run = find_current_p4_bold_run(dataset, obs, autodl_bold_root)
            current_p7_by_run_tag[run_tag] = find_p7_run_for_p4(ds_root, p4_run)

        # P3: current four-observable QC summaries.
        p3_qc = ds_root / "pipeline3_figures_bold_pre_reskoopnet_qc"
        for obs in CURRENT_P3_OBSERVABLES:
            qc = latest_file(p3_qc / obs, "*.png") or latest_file(p3_qc, f"*{obs}*.png")
            if qc is not None:
                add_source(dataset, "P3", "observable_qc", p3_qc, qc)

        # P5: use canonical per-run reduction figure folders, not
        # pipeline5_summary_figures cache folders. Keep first-three trajectory
        # figures explicitly labeled as diagnostic until selected-component
        # trajectory manifests exist.
        for condition in CURRENT_P5_CONDITIONS:
            if invalid_blp_condition_detail(dataset, condition):
                continue
            for method in CURRENT_P5_METHODS:
                for k in CURRENT_P5_COMPONENT_COUNTS:
                    tag = f"{method}_k{k:02d}"
                    fig_dir = ds_root / "pipeline5_eigenfunction_reduction" / condition / tag / "fig"
                    p5_specs = (
                        ("reduction_overview", "*__comp__*.png"),
                        ("state_space_first3_diagnostic", "*__ss__*.png"),
                        ("consensus_state_space_first3_diagnostic", "*__ssc__*.png"),
                        ("spectrum_diagnostics", "*__spec__*.png"),
                    )
                    for group, pattern in p5_specs:
                        fig = latest_file(fig_dir, pattern)
                        if fig is not None:
                            add_source(
                                dataset,
                                "P5",
                                group,
                                fig_dir,
                                fig,
                                detail=f"{condition}/{tag}; canonical reduction fig folder",
                            )

        patterns = [
            ("P6", "timescale_diagnostics", ds_root / "pipeline6_figures_timescale_diagnostics", "*.png"),
            ("P6", "spkt_residual_xcorr_overview", ds_root / "pipeline6_figures_spkt_residual_cross_correlation", "*.png"),
            ("P6", "mua_residual_xcorr_overview", ds_root / "pipeline6_figures_mua_residual_cross_correlation", "*.png"),
            ("P6", "top_window_edmd_and_norm_deconv", ds_root / "pipeline6_top_state_diversity_postprocessing", "**/01_postprocess_main/**/*.png"),
            ("P6", "top_window_edmd_and_norm_deconv", ds_root / "pipeline6_top_state_diversity_postprocessing", "**/04_deconv_localwin_norm/**/*.png"),
        ]
        for pipeline, group, root, pattern in patterns:
            if not root.is_dir():
                continue
            for src in root.glob(pattern):
                add_source(dataset, pipeline, group, root, src)

        # P7/P8/P9/P10 flat candidates are restricted to artifacts that can
        # trace back to the best-validation BOLD P4 run through the current
        # P7 folder.
        for run_tag, p7_run in current_p7_by_run_tag.items():
            if p7_run is None:
                continue
            p7_patterns = (
                ("main_and_deconv_efuns", "fig/*_efuns.png"),
                ("main_and_deconv_efuns", "fig/*_deconv_efuns.png"),
                ("timescale_diagnostics", "fig/*_timescale.png"),
                ("intrinsic_roi", "fig/intrinsic_roi_bar_summaries/*.png"),
                ("intrinsic_activation", "fig/intrinsic_activation_maps/*.png"),
            )
            for group, pattern in p7_patterns:
                for src in p7_run.glob(pattern):
                    add_source(dataset, "P7", group, p7_run, src, detail=f"current P7 run_tag={run_tag}")

            p8_xcorr = ds_root / "pipeline8_xcorr" / run_tag
            p8_patterns = (
                ("xcorr_combined_overview", "*_summary.png"),
                ("xcorr_combined_overview", "*_top_signal_overlay.png"),
                ("xcorr_per_density_overview", "density/*_summary.png"),
                ("xcorr_per_density_overview", "density/*_top_signal_overlay.png"),
                ("xcorr_per_density_feature_overview", "feature/*/*_summary.png"),
                ("xcorr_per_density_feature_overview", "feature/*/*_top_signal_overlay.png"),
            )
            if p8_xcorr.is_dir():
                for group, pattern in p8_patterns:
                    for src in p8_xcorr.glob(pattern):
                        add_source(dataset, "P8", group, p8_xcorr, src, detail=f"current P7 run_tag={run_tag}")

            p8_top = ds_root / "pipeline8_top_maps" / run_tag
            p8_top_patterns = (
                ("roi_summary", "**/fig/**/roi/*.png"),
                ("activation_maps", "**/fig/**/act/*.png"),
            )
            if p8_top.is_dir():
                for group, pattern in p8_top_patterns:
                    for src in p8_top.glob(pattern):
                        add_source(dataset, "P8", group, p8_top, src, detail=f"current P7 run_tag={run_tag}")

            p9_run = ds_root / "pipeline9_bold_eigenfunction_reduction" / run_tag
            p9_status, _p9_detail = p9_run_status(ds_root, run_tag, p7_run)
            if p9_run.is_dir() and p9_status == "ok":
                valid_p9_run_tags.add(run_tag)
                for src in p9_run.glob("**/fig/*.png"):
                    add_source(dataset, "P9", "bold_dimred_summary", p9_run, src, detail=f"current P7 run_tag={run_tag}")

        for candidate in discover_p10_candidates(
            ds_root,
            dataset,
            p9_min_mtime_by_tag(current_p7_by_run_tag),
        ):
            if candidate.run_tag not in valid_p9_run_tags:
                continue
            p10_patterns = (
                ("xcorr_combined_overview", candidate.xcorr_dir, "*_summary.png"),
                ("xcorr_combined_overview", candidate.xcorr_dir, "*_top_signal_overlay.png"),
                ("xcorr_per_density_overview", candidate.xcorr_dir, "density/*_summary.png"),
                ("xcorr_per_density_overview", candidate.xcorr_dir, "density/*_top_signal_overlay.png"),
                ("xcorr_per_density_feature_overview", candidate.xcorr_dir, "feature/*/*_summary.png"),
                ("xcorr_per_density_feature_overview", candidate.xcorr_dir, "feature/*/*_top_signal_overlay.png"),
                ("roi_summary", candidate.top_maps_dir, "**/fig/**/roi/*.png"),
                ("activation_maps", candidate.top_maps_dir, "**/fig/**/activation_maps_top*/*.png"),
            )
            for group, root, pattern in p10_patterns:
                if not root.is_dir():
                    continue
                for src in root.glob(pattern):
                    add_source(dataset, "P10", group, root, src, detail=f"current P9 candidate={candidate.p10_tag}")

    # Global derived P5 browsing products. These are current derived outputs,
    # not canonical sources; copying them into P11 keeps browsing flat while
    # their provenance is still recorded in this manifest.
    peak_stats_root = processed_root / "summary_figures" / "pipeline5_peak_statistics_maxabs"
    if peak_stats_root.is_dir():
        for src in peak_stats_root.glob("*.png"):
            dataset = src.name.split("__", 1)[0] if "__" in src.name else "all"
            if (
                dataset != "all"
                and "complex_split" in src.name.lower()
                and invalid_blp_condition_detail(dataset, "complex_split_projected_vlambda")
            ):
                continue
            add_source(dataset, "P5", "peak_statistics_summary", peak_stats_root, src, "derived from maxabs peak-state CSV")

    for result_name in P5_DERIVED_SELECTIVITY_RESULTS:
        result_root = repo_results_root / result_name / "figures"
        if result_root.is_dir():
            for src in result_root.glob("*.png"):
                add_source("all", "P5", "event_family_selectivity", result_root, src, "repo-local derived P5 selectivity figure")

    if repo_results_root.is_dir():
        for result_dir in sorted(repo_results_root.glob(f"{P5_DERIVED_CONSISTENCY_PREFIX}*")):
            if not result_dir.is_dir():
                continue
            if "all_datasets" in result_dir.name:
                dataset = "all"
                group = "cross_dataset_method_consistency"
            else:
                match = re.search(r"_k\d+_(.+)$", result_dir.name)
                dataset = match.group(1) if match else "all"
                group = "per_dataset_method_consistency"
            for src in result_dir.glob("*.png"):
                # Use repo_results_root as the anchor so the flat filename keeps
                # the result directory name. Otherwise k03..k08 outputs share
                # basenames and overwrite each other in the P11 browsing folder.
                add_source(dataset, "P5", group, repo_results_root, src, "repo-local derived P5 method-consistency figure")

    p8_consistency_root = repo_results_root / "pipeline8_cross_session_consistency_current" / "figures"
    if p8_consistency_root.is_dir():
        for src in p8_consistency_root.glob("*.png"):
            add_source("all", "P8", "cross_session_consistency", p8_consistency_root, src, "repo-local derived P8 cross-session consistency figure")
    scorecard_root = repo_results_root / "pipeline11_parameter_selection_scorecard_current" / "figures"
    if scorecard_root.is_dir():
        for src in scorecard_root.glob("*.png"):
            add_source(
                "all",
                "P11",
                "parameter_selection_scorecard",
                scorecard_root,
                src,
                "repo-local derived P11 parameter-selection scorecard figure",
            )
    return rows


def apply_flat_copies(rows: Sequence[FlatFigureSource], overwrite: bool, max_copies: Optional[int]) -> list[FlatFigureSource]:
    out: list[FlatFigureSource] = []
    copied = 0
    for row in rows:
        if max_copies is not None and copied >= max_copies:
            out.append(
                FlatFigureSource(**{**row.__dict__, "status": "skipped_max_flat_copies", "detail": ""})
            )
            continue
        src = Path(row.source_path)
        dst = Path(row.proposed_flat_path)
        if not src.is_file():
            out.append(FlatFigureSource(**{**row.__dict__, "status": "missing_source", "detail": ""}))
            continue
        actual = copy_file(src, dst, overwrite)
        copied += 1
        out.append(
            FlatFigureSource(
                **{
                    **row.__dict__,
                    "proposed_flat_path": str(actual),
                    "status": "copied",
                    "detail": "overwritten" if overwrite and actual == dst else "",
                }
            )
        )
    return out


def apply_quarantine(rows: Sequence[LegacyCandidate]) -> list[LegacyCandidate]:
    out: list[LegacyCandidate] = []
    for row in rows:
        src = Path(row.source_path)
        dst = Path(row.proposed_legacy_path)
        if not src.exists():
            out.append(LegacyCandidate(**{**row.__dict__, "status": "missing_source"}))
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        target = unique_path(dst)
        shutil.move(str(src), str(target))
        out.append(LegacyCandidate(**{**row.__dict__, "proposed_legacy_path": str(target), "status": "moved"}))
    return out


def archive_existing_summary(summary_root: Path, processed_root: Path, stamp: str) -> Optional[Path]:
    if not summary_root.exists():
        return None
    target = processed_root / "_legacy_quarantine" / stamp / "summary_figures" / summary_root.name
    target.parent.mkdir(parents=True, exist_ok=True)
    target = unique_path(target)
    shutil.move(str(summary_root), str(target))
    return target


def summarize_status(rows: Iterable[StageStatus]) -> str:
    counts: dict[tuple[str, str], int] = {}
    for row in rows:
        key = (row.pipeline, row.status)
        counts[key] = counts.get(key, 0) + 1
    lines = []
    for (pipeline, status), count in sorted(counts.items()):
        lines.append(f"{pipeline},{status},{count}")
    return "\n".join(lines)


def main() -> int:
    args = parse_args()
    processed_root = args.processed_root
    output_dir = args.output_dir
    summary_root = args.summary_root
    repo_results_root = args.repo_results_root
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    stage_rows = audit_stage_status(
        processed_root,
        args.datasets,
        args.autodl_blp_root,
        args.autodl_bold_root,
        repo_results_root,
    )
    fill_rows = build_fill_plan(processed_root, args.datasets, args.autodl_bold_root)
    legacy_rows = build_legacy_candidates(processed_root, args.datasets, stamp)
    flat_rows = collect_flat_sources(
        processed_root,
        repo_results_root,
        summary_root,
        args.datasets,
        args.autodl_bold_root,
    )

    archived_summary = None
    if args.archive_existing_summary:
        archived_summary = archive_existing_summary(summary_root, processed_root, stamp)
    if args.copy_flat_figures:
        flat_rows = apply_flat_copies(flat_rows, args.overwrite_flat, args.max_flat_copies)
    if args.quarantine_legacy:
        legacy_rows = apply_quarantine(legacy_rows)

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(output_dir / f"p11_stage_status_{stamp}.csv", stage_rows)
    write_csv(output_dir / "p11_stage_status_latest.csv", stage_rows)
    write_csv(output_dir / f"p11_fill_plan_{stamp}.csv", fill_rows)
    write_csv(output_dir / "p11_fill_plan_latest.csv", fill_rows)
    write_csv(output_dir / f"p11_legacy_candidates_{stamp}.csv", legacy_rows)
    write_csv(output_dir / "p11_legacy_candidates_latest.csv", legacy_rows)
    write_csv(output_dir / f"p11_flat_figure_sources_{stamp}.csv", flat_rows)
    write_csv(output_dir / "p11_flat_figure_sources_latest.csv", flat_rows)

    print(f"Processed root : {processed_root}")
    print(f"Repo results   : {repo_results_root.resolve()}")
    print(f"Output dir     : {output_dir.resolve()}")
    print(f"P11 summary dir: {summary_root}")
    print(f"Stage rows     : {len(stage_rows)}")
    print(f"Fill actions   : {len(fill_rows)}")
    print(f"Legacy rows    : {len(legacy_rows)}")
    print(f"Flat sources   : {len(flat_rows)}")
    if archived_summary is not None:
        print(f"Archived prior P11 summary: {archived_summary}")
    print("")
    print("stage,status,count")
    print(summarize_status(stage_rows))
    if not args.copy_flat_figures:
        print("\nFlat figures were planned only. Use --copy-flat-figures to copy.")
    if not args.quarantine_legacy:
        print("Legacy candidates were planned only. Use --quarantine-legacy to move them to _legacy_quarantine.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
