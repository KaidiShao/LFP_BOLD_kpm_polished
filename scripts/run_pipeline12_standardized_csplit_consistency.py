#!/usr/bin/env python3
"""Pipeline12 standardized complex-split consistency check.

This script consumes existing P5/P8/P10 outputs for one dataset and summarizes
whether standardized complex-split BLP density sources reproduce the E10gb1
slow-variable / intrinsic-trigger pattern.  It is intentionally strict: if the
standardized branch is missing, it reports missing requirements instead of
falling back to legacy/nonstandard results.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import math
import re
from collections import Counter, defaultdict
from itertools import count
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_RESULT_ROOT = REPO_ROOT / "results" / "pipeline12_standardized_csplit_consistency"
DEFAULT_LABEL_TABLE = (
    REPO_ROOT
    / "results"
    / "pipeline5_dimred_component_process_labels_standardized_csplit_rmsenv_adaptive_7datasets"
    / "dimred_efun_process_labels.csv"
)
DEFAULT_RAW_METADATA = (
    REPO_ROOT
    / "results"
    / "p8_p10_parameter_selection_current_rmsenv_adaptive_v2"
    / "p5_raw_density_mode_metadata.csv"
)

METHODS = ("svd", "nmf", "mds", "umap")
KS = tuple(range(3, 9))
RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")
RUN_TAG_LABELS = {
    "pv_gsvd100": "global_svd100",
    "pv_gsvd100_ds": "gsvd100_ds",
    "pv_hp100": "HP_svd100",
    "pv_roi": "roi_mean",
    "pv_roi_mean": "roi_mean",
}
FEATURE_FAMILIES = ("efun", "deconv_efun")
DENSITY_CLASSES = ("event_density", "raw_efun_density", "dimred_efun_density")
LABEL_ORDER = (
    "theta_strict",
    "theta_selective_similar",
    "theta_selective_unequal",
    "theta_relaxed",
    "ripple_gamma_no_pure_theta",
    "ripple_selective_similar",
    "ripple_selective_unequal",
    "mixed_ripple_gamma",
    "mixed_theta_ripple",
    "partial_or_inactive",
    "other",
    "unlabeled",
    "label_missing",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--figure-dir", type=Path, default=None)
    parser.add_argument("--p4-summary-mat", type=Path, default=None)
    parser.add_argument(
        "--p4-full-output-dir",
        type=Path,
        default=None,
        help="Optional directory containing exported P4 *_outputs_*.mat chunks when the summary MAT does not store efuns.",
    )
    parser.add_argument(
        "--condition-tag",
        default="complex_split_projected_vlambda_standardize_rmsenv_adaptive",
    )
    parser.add_argument("--p8-save-tag", default="xcorr_csplit_standardize_rmsenv_adaptive")
    parser.add_argument("--p10-save-tag", default="dimred_xcorr_csplit_standardize_rmsenv_adaptive")
    parser.add_argument("--label-table", type=Path, default=DEFAULT_LABEL_TABLE)
    parser.add_argument("--raw-metadata-csv", type=Path, default=DEFAULT_RAW_METADATA)
    parser.add_argument("--run-tags", nargs="+", default=list(RUN_TAGS))
    parser.add_argument("--top-n-values", nargs="+", type=int, default=[1, 3, 5, 10, 20])
    parser.add_argument("--max-rows-per-group", type=int, default=50)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: object) -> float:
    try:
        out = float(str(value))
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def as_int(value: object, default: int = -1) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def safe_median(values: Iterable[float]) -> float | str:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return ""
    return float(median(clean))


def feature_family(text: str) -> str:
    return "deconv_efun" if "deconv" in str(text).lower() else "efun"


def observable_mode(run_tag: str) -> str:
    return RUN_TAG_LABELS.get(run_tag, run_tag)


def normalize_density_name(name: str) -> str:
    text = str(name or "").lower()
    text = re.sub(r"_standardize_rmsenv_adaptive$", "", text)
    text = re.sub(r"_rmsenv_adaptive$", "", text)
    return text


def parse_density_name(name: str) -> dict[str, object]:
    text = str(name or "").lower()
    base = normalize_density_name(text)
    if base.startswith("blp_evt"):
        return {
            "density_class": "event_density",
            "condition": "event",
            "method": "",
            "k": "",
        }
    if base.startswith("raw_csplit"):
        return {
            "density_class": "raw_efun_density",
            "condition": "csplit",
            "method": "",
            "k": "",
        }
    match = re.match(r"dim_csplit_([a-z]+)(\d+)_q\d+", base)
    if match:
        return {
            "density_class": "dimred_efun_density",
            "condition": "csplit",
            "method": match.group(1),
            "k": int(match.group(2)),
        }
    return {
        "density_class": "other",
        "condition": "",
        "method": "",
        "k": "",
    }


def norm_path(value: object) -> str:
    return str(value or "").replace("/", "\\").lower()


def label_from_custom_row(row: dict[str, str]) -> str:
    def truthy(key: str) -> bool:
        return str(row.get(key, "")).strip().lower() in {"1", "true", "yes", "y"}

    theta_strict = truthy("theta_strict")
    theta_relaxed = truthy("theta_relaxed")
    ripple_gamma = truthy("ripple_gamma_no_pure_theta")
    ripple_no_theta = truthy("ripple_no_pure_theta")
    if theta_strict and ripple_gamma:
        return "mixed_theta_ripple"
    if theta_strict:
        return "theta_strict"
    if theta_relaxed and ripple_gamma:
        return "mixed_theta_ripple"
    if theta_relaxed:
        return "theta_relaxed"
    if ripple_gamma:
        return "ripple_gamma_no_pure_theta"
    if ripple_no_theta:
        return "ripple_selective_unequal"
    return "other"


def load_labels(path: Path) -> tuple[dict[tuple[str, str, int], str], dict[tuple[str, str, int, int], str]]:
    by_density: dict[tuple[str, str, int], str] = {}
    by_method: dict[tuple[str, str, int, int], str] = {}
    for row in read_csv(path):
        if "primary_process_label" in row:
            dataset = str(row.get("dataset", "")).lower()
            density = normalize_density_name(row.get("density_name_for_p8_p10") or row.get("density_name") or "")
            idx = as_int(row.get("density_index") or row.get("component_idx"))
            label = str(row.get("primary_process_label") or "label_missing")
            if dataset and density and idx >= 0:
                by_density[(dataset, density, idx)] = label
            method = str(row.get("method", "")).lower()
            k = as_int(row.get("k") or row.get("component_count"))
            condition = "csplit" if "csplit" in density else str(row.get("condition_short", "")).lower()
            if dataset and method and k >= 0 and idx >= 0:
                by_method[(dataset, f"{condition}:{method}", k, idx)] = label
        elif "theta_strict" in row:
            branch = str(row.get("branch", "")).lower()
            condition = str(row.get("condition", "")).lower()
            if branch and branch != "standardize":
                continue
            if condition and condition != "csplit":
                continue
            method = str(row.get("method", "")).lower()
            k = as_int(row.get("k"))
            idx = as_int(row.get("component_idx"))
            label = label_from_custom_row(row)
            if method and k >= 0 and idx >= 0:
                by_method[("*", f"csplit:{method}", k, idx)] = label
    return by_density, by_method


def lookup_label(
    dataset: str,
    density_name: str,
    density_index: int,
    method: str,
    k: int | str,
    by_density: dict[tuple[str, str, int], str],
    by_method: dict[tuple[str, str, int, int], str],
) -> str:
    density = normalize_density_name(density_name)
    dataset_lc = dataset.lower()
    label = by_density.get((dataset_lc, density, density_index))
    if label:
        return label
    k_int = as_int(k)
    for key in (
        (dataset_lc, f"csplit:{method}", k_int, density_index),
        (dataset_lc, "csplit", k_int, density_index),
        ("*", f"csplit:{method}", k_int, density_index),
        ("*", "csplit", k_int, density_index),
    ):
        label = by_method.get(key)
        if label:
            return label
    return "label_missing"


def load_raw_metadata(path: Path) -> dict[tuple[str, int], dict[str, str]]:
    out: dict[tuple[str, int], dict[str, str]] = {}
    for row in read_csv(path):
        density_file = norm_path(row.get("density_file"))
        idx = as_int(row.get("density_index"))
        if density_file and idx >= 0:
            out[(density_file, idx)] = row
    return out


def raw_meta_for(raw_metadata: dict[tuple[str, int], dict[str, str]], row: dict[str, object]) -> dict[str, str]:
    idx = as_int(row.get("density_index"))
    return raw_metadata.get((norm_path(row.get("density_file")), idx), {})


class TopKStore:
    def __init__(self, k: int) -> None:
        self.k = k
        self._counter = count()
        self._data: dict[tuple[object, ...], list[tuple[float, int, dict[str, object]]]] = defaultdict(list)

    def add(self, key: tuple[object, ...], row: dict[str, object], score: float) -> None:
        if not math.isfinite(score):
            return
        heap = self._data[key]
        item = (score, next(self._counter), row)
        if len(heap) < self.k:
            heapq.heappush(heap, item)
        elif score > heap[0][0]:
            heapq.heapreplace(heap, item)

    def rows(self) -> list[dict[str, object]]:
        out: list[dict[str, object]] = []
        for key, heap in self._data.items():
            ranked = sorted(heap, key=lambda item: item[0], reverse=True)
            for rank, (_, _, row) in enumerate(ranked, start=1):
                item = dict(row)
                item["rank_within_group"] = rank
                item["top_group_key"] = "|".join(str(x) for x in key)
                out.append(item)
        return out

    def rows_for_key(self, key: tuple[object, ...]) -> list[dict[str, object]]:
        ranked = sorted(self._data.get(key, []), key=lambda item: item[0], reverse=True)
        return [dict(item[2], rank_within_group=i + 1) for i, item in enumerate(ranked)]


def p10_dir_parts(name: str) -> tuple[str, str, str, str]:
    parts = name.split("__")
    if len(parts) < 3:
        return name, "", "", ""
    run_tag = parts[0]
    bold_feature = parts[1]
    method_k = parts[2]
    method = method_k.split("_k")[0].lower() if "_k" in method_k else method_k.lower()
    return run_tag, bold_feature, method_k, method


def enrich_peak_row(
    *,
    pipeline: str,
    dataset: str,
    run_tag: str,
    source_file: Path,
    source_row_rank: int,
    row: dict[str, str],
    label_maps: tuple[dict[tuple[str, str, int], str], dict[tuple[str, str, int, int], str]],
    raw_metadata: dict[tuple[str, int], dict[str, str]],
    p10_bold_feature: str = "",
    p10_method_k: str = "",
) -> dict[str, object]:
    density_name = str(row.get("density_name", ""))
    meta = parse_density_name(density_name)
    bold_feature = str(row.get("bold_feature") or p10_bold_feature)
    fam = feature_family(bold_feature or p10_bold_feature)
    density_index = as_int(row.get("density_index"))
    score = as_float(row.get("peak_abs_corr"))
    method = str(meta["method"])
    k = meta["k"]
    label = ""
    if meta["density_class"] == "dimred_efun_density":
        label = lookup_label(dataset, density_name, density_index, method, k, label_maps[0], label_maps[1])
    out: dict[str, object] = {
        "pipeline": pipeline,
        "dataset": dataset,
        "run_tag": run_tag,
        "observable_mode": observable_mode(run_tag),
        "bold_feature": bold_feature,
        "feature_family": fam,
        "p10_method_k": p10_method_k,
        "source_csv": str(source_file),
        "source_row_rank": source_row_rank,
        "density_name": density_name,
        "density_class": meta["density_class"],
        "density_condition": meta["condition"],
        "blp_method": method,
        "blp_k": k,
        "density_file": row.get("density_file", ""),
        "density_index": density_index,
        "density_label": row.get("density_label", ""),
        "bold_mode_index": row.get("bold_mode_index", ""),
        "bold_component_index": row.get("bold_component_index", ""),
        "peak_abs_corr": score,
        "peak_corr": as_float(row.get("peak_corr")),
        "peak_lag_sec": as_float(row.get("peak_lag_sec")),
        "two_subprocess_label": label,
    }
    if meta["density_class"] == "raw_efun_density":
        raw_meta = raw_meta_for(raw_metadata, out)
        for key in (
            "raw_efun_index",
            "source_mode_index",
            "timescale_sec_preferred",
            "frequency_hz_preferred",
            "timescale_sec_discrete_log",
            "frequency_hz_discrete_angle",
            "timescale_source_preferred",
            "envelope_window_sec",
        ):
            out[key] = raw_meta.get(key, "")
    return out


def collect_rows(args: argparse.Namespace, label_maps, raw_metadata):
    max_k = max(args.top_n_values + [args.max_rows_per_group])
    top_all = TopKStore(max_k)
    top_class = TopKStore(max_k)
    top_dimred = TopKStore(max_k)
    top_raw = TopKStore(max_k)
    top_method_k = TopKStore(1)
    file_rows: list[dict[str, object]] = []

    def handle_file(pipeline: str, run_tag: str, source_file: Path, p10_bold_feature: str = "", p10_method_k: str = "") -> None:
        n_rows = 0
        with source_file.open("r", newline="", encoding="utf-8-sig") as handle:
            for source_row_rank, row in enumerate(csv.DictReader(handle), start=1):
                peak = enrich_peak_row(
                    pipeline=pipeline,
                    dataset=args.dataset,
                    run_tag=run_tag,
                    source_file=source_file,
                    source_row_rank=source_row_rank,
                    row=row,
                    label_maps=label_maps,
                    raw_metadata=raw_metadata,
                    p10_bold_feature=p10_bold_feature,
                    p10_method_k=p10_method_k,
                )
                if peak["density_class"] not in DENSITY_CLASSES:
                    continue
                n_rows += 1
                context = (
                    pipeline,
                    run_tag,
                    peak["observable_mode"],
                    peak["bold_feature"],
                    peak["feature_family"],
                )
                score = as_float(peak["peak_abs_corr"])
                top_all.add(context, peak, score)
                top_class.add(context + (peak["density_class"],), peak, score)
                if peak["density_class"] == "dimred_efun_density":
                    top_dimred.add((pipeline, peak["observable_mode"], peak["feature_family"]), peak, score)
                    top_method_k.add(
                        (
                            pipeline,
                            peak["observable_mode"],
                            peak["feature_family"],
                            peak["blp_method"],
                            peak["blp_k"],
                        ),
                        peak,
                        score,
                    )
                elif peak["density_class"] == "raw_efun_density":
                    top_raw.add((pipeline, peak["observable_mode"], peak["feature_family"]), peak, score)
        file_rows.append(
            {
                "pipeline": pipeline,
                "run_tag": run_tag,
                "source_csv": str(source_file),
                "rows_read": n_rows,
            }
        )

    p8_root = args.processed_root / args.dataset / "pipeline8_xcorr"
    for run_tag in args.run_tags:
        density_dir = p8_root / run_tag / "density"
        for source_file in sorted(density_dir.glob(f"{args.p8_save_tag}_peaks__*.csv")):
            if keep_density_file(source_file.name):
                handle_file("P8", run_tag, source_file)

    p10_root = args.processed_root / args.dataset / "pipeline10_dimred_xcorr"
    for run_dir in sorted(p10_root.glob("*")):
        if not run_dir.is_dir():
            continue
        run_tag, bold_feature, method_k, _ = p10_dir_parts(run_dir.name)
        if run_tag not in set(args.run_tags):
            continue
        density_dir = run_dir / "density"
        for source_file in sorted(density_dir.glob(f"{args.p10_save_tag}_peaks__*.csv")):
            if keep_density_file(source_file.name):
                handle_file("P10", run_tag, source_file, p10_bold_feature=bold_feature, p10_method_k=method_k)

    return {
        "file_rows": file_rows,
        "top_rows": top_all.rows(),
        "top_class_rows": top_class.rows(),
        "top_dimred_rows": top_dimred.rows(),
        "top_raw_rows": top_raw.rows(),
        "method_k_rows": top_method_k.rows(),
    }


def keep_density_file(name: str) -> bool:
    text = name.lower()
    return "__blp_evt" in text or "__raw_csplit" in text or "__dim_csplit" in text


def audit_inputs(
    args: argparse.Namespace,
    label_count: int,
    raw_meta_count: int,
    file_rows: list[dict[str, object]],
    top_rows: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    missing: list[dict[str, object]] = []

    def add(item: str, status: str, detail: object = "", required: bool = True) -> None:
        rows.append({"item": item, "status": status, "detail": detail, "required": int(required)})
        if required and status not in {"ok", "present"}:
            missing.append({"item": item, "status": status, "detail": detail})

    contains_efuns = False
    p4_chunk_dirs: list[Path] = []
    if args.p4_summary_mat:
        if args.p4_summary_mat.is_file():
            contains_efuns, fields, err = inspect_p4_summary(args.p4_summary_mat)
            p4_chunk_dirs.append(args.p4_summary_mat.parent)
            add("p4_summary_mat", "ok", args.p4_summary_mat)
            add("p4_summary_contains_efuns", "ok" if contains_efuns else "missing", f"fields={';'.join(fields[:20])}; error={err}", required=False)
        else:
            add("p4_summary_mat", "missing", args.p4_summary_mat)
    else:
        add("p4_summary_mat", "not_checked", "no --p4-summary-mat provided", required=False)

    if args.p4_full_output_dir:
        if args.p4_full_output_dir.is_dir():
            p4_chunk_dirs.append(args.p4_full_output_dir)
            add("p4_full_output_dir", "ok", args.p4_full_output_dir, required=False)
        else:
            add("p4_full_output_dir", "missing", args.p4_full_output_dir, required=False)

    seen_dirs: set[str] = set()
    chunk_details: list[str] = []
    chunk_count = 0
    for chunk_dir in p4_chunk_dirs:
        key = str(chunk_dir).lower()
        if key in seen_dirs:
            continue
        seen_dirs.add(key)
        count_here = len(list(chunk_dir.glob("*_outputs_*.mat"))) if chunk_dir.is_dir() else 0
        chunk_count += count_here
        chunk_details.append(f"{count_here} in {chunk_dir}")
    add("p4_full_output_chunks", "ok" if chunk_count > 0 else "missing", "; ".join(chunk_details) or "no chunk dirs checked", required=False)
    add("p4_full_efun_available", "ok" if contains_efuns or chunk_count > 0 else "missing", "summary efuns or *_outputs_*.mat chunks")

    raw_dir = args.processed_root / args.dataset / "pipeline5_raw_thresholded_density" / args.condition_tag / "mat"
    raw_count = len(list(raw_dir.glob("*.mat"))) if raw_dir.is_dir() else 0
    add("p5_raw_density_standardized_csplit", "ok" if raw_count > 0 else "missing", f"{raw_count} mat files in {raw_dir}")

    dim_root = args.processed_root / args.dataset / "pipeline5_dimred_thresholded_density" / args.condition_tag
    dim_expected = len(METHODS) * len(KS)
    dim_present = 0
    for method in METHODS:
        for k in KS:
            mat_dir = dim_root / f"{method}_k{k:02d}" / "mat"
            if mat_dir.is_dir() and any(mat_dir.glob("*.mat")):
                dim_present += 1
    add("p5_dimred_density_standardized_csplit", "ok" if dim_present == dim_expected else "missing", f"{dim_present}/{dim_expected} method-k density dirs in {dim_root}")

    add("p5_dimred_component_labels", "ok" if label_count > 0 else "missing", f"{label_count} labels loaded from {args.label_table}", required=False)
    add("raw_efun_timescale_metadata", "ok" if raw_meta_count > 0 else "missing", f"{raw_meta_count} metadata rows from {args.raw_metadata_csv}", required=False)

    p8_files = [r for r in file_rows if r["pipeline"] == "P8"]
    p10_files = [r for r in file_rows if r["pipeline"] == "P10"]
    add("p8_standardized_csplit_xcorr_peak_files", "ok" if p8_files else "missing", f"{len(p8_files)} peak files matched save tag {args.p8_save_tag}")
    add("p10_standardized_csplit_xcorr_peak_files", "ok" if p10_files else "missing", f"{len(p10_files)} peak files matched save tag {args.p10_save_tag}")
    expected_families = set(FEATURE_FAMILIES)
    for pipeline in ("P8", "P10"):
        families = {
            str(row.get("feature_family", ""))
            for row in top_rows
            if str(row.get("pipeline", "")) == pipeline
        }
        status = "ok" if expected_families.issubset(families) else "missing"
        add(
            f"{pipeline.lower()}_efun_and_deconv_feature_families",
            status,
            ",".join(sorted(families)) or "none",
        )
        for family in FEATURE_FAMILIES:
            present_run_tags = {
                str(row.get("run_tag", ""))
                for row in top_rows
                if str(row.get("pipeline", "")) == pipeline
                and str(row.get("feature_family", "")) == family
            }
            missing_run_tags = [tag for tag in args.run_tags if tag not in present_run_tags]
            status = "ok" if not missing_run_tags else "missing"
            detail = (
                f"present={','.join(sorted(present_run_tags)) or 'none'}; "
                f"missing={','.join(missing_run_tags) or 'none'}"
            )
            add(f"{pipeline.lower()}_{family}_observable_coverage", status, detail)
    return rows, missing


def inspect_p4_summary(path: Path) -> tuple[bool, list[str], str]:
    try:
        import scipy.io

        mat = scipy.io.loadmat(path, simplify_cells=True)
        edmd = mat.get("EDMD_outputs", {})
        if isinstance(edmd, dict):
            fields = list(edmd.keys())
            return "efuns" in edmd, fields, ""
        fields = []
        if hasattr(edmd, "dtype") and getattr(edmd.dtype, "names", None):
            fields = list(edmd.dtype.names or [])
        return "efuns" in fields, fields, ""
    except Exception as exc:  # pragma: no cover - audit path
        return False, [], f"{type(exc).__name__}: {exc}"


def density_class_best_scores(top_class_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in top_class_rows:
        if as_int(row.get("rank_within_group")) != 1:
            continue
        out.append(row)
    return out


def topn_membership(top_rows: list[dict[str, object]], top_ns: Sequence[int]) -> list[dict[str, object]]:
    by_context: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in top_rows:
        by_context[str(row.get("top_group_key", ""))].append(row)
    out: list[dict[str, object]] = []
    for _, rows in by_context.items():
        rows_sorted = sorted(rows, key=lambda r: as_float(r.get("peak_abs_corr")), reverse=True)
        first = rows_sorted[0]
        for top_n in top_ns:
            subset = rows_sorted[:top_n]
            counts = Counter(str(r.get("density_class")) for r in subset)
            for density_class in DENSITY_CLASSES:
                out.append(
                    {
                        "pipeline": first.get("pipeline"),
                        "dataset": first.get("dataset"),
                        "run_tag": first.get("run_tag"),
                        "observable_mode": first.get("observable_mode"),
                        "bold_feature": first.get("bold_feature"),
                        "feature_family": first.get("feature_family"),
                        "top_n": top_n,
                        "density_class": density_class,
                        "count": counts.get(density_class, 0),
                        "fraction": counts.get(density_class, 0) / float(max(1, len(subset))),
                    }
                )
    return out


def save_placeholder(path: Path, title: str, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9, 4.8))
    ax.axis("off")
    ax.text(0.02, 0.82, title, fontsize=16, weight="bold", transform=ax.transAxes)
    ax.text(0.02, 0.55, message, fontsize=11, transform=ax.transAxes, va="top", wrap=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_availability(rows: list[dict[str, object]], out: Path) -> None:
    labels = [str(r["item"]) for r in rows]
    vals = [1 if str(r["status"]) == "ok" else 0 for r in rows]
    colors = ["#2f7d32" if v else "#b23b32" for v in vals]
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, max(4, 0.42 * len(labels))))
    ax.barh(labels, vals, color=colors)
    ax.set_xlim(0, 1)
    ax.set_xlabel("available")
    ax.set_title("Pipeline12 standardized csplit availability")
    for i, row in enumerate(rows):
        ax.text(0.03, i, str(row["status"]), va="center", color="white" if vals[i] else "black")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_best_score_heatmaps(rows: list[dict[str, object]], fig_dir: Path) -> None:
    if not rows:
        save_placeholder(fig_dir / "density_class_best_score__missing.png", "Density class best score", "No standardized csplit P8/P10 rows available.")
        return
    for pipeline in ("P8", "P10"):
        for fam in FEATURE_FAMILIES:
            subset = [r for r in rows if r.get("pipeline") == pipeline and r.get("feature_family") == fam]
            out = fig_dir / f"density_class_best_score__{pipeline.lower()}__{fam}.png"
            if not subset:
                save_placeholder(out, f"{pipeline} {fam}", "No rows available.")
                continue
            observables = sorted({str(r.get("observable_mode")) for r in subset})
            matrix: list[list[float]] = []
            for obs in observables:
                row_vals = []
                for cls in DENSITY_CLASSES:
                    vals = [as_float(r.get("peak_abs_corr")) for r in subset if r.get("observable_mode") == obs and r.get("density_class") == cls]
                    med = safe_median(vals)
                    row_vals.append(float(med) if med != "" else math.nan)
                matrix.append(row_vals)
            plot_heatmap(matrix, observables, list(DENSITY_CLASSES), out, f"{pipeline} {fam}: best |xcorr| by density class", "median best |r|")


def plot_heatmap(matrix: list[list[float]], ylabels: list[str], xlabels: list[str], out: Path, title: str, cbar_label: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    if not matrix:
        save_placeholder(out, title, "No data.")
        return
    arr = matrix
    fig, ax = plt.subplots(figsize=(max(7, 1.5 * len(xlabels)), max(3.6, 0.55 * len(ylabels) + 1.8)))
    im = ax.imshow(arr, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(xlabels)), xlabels, rotation=25, ha="right")
    ax.set_yticks(range(len(ylabels)), ylabels)
    ax.set_title(title)
    for i, row in enumerate(arr):
        for j, value in enumerate(row):
            if math.isfinite(value):
                ax.text(j, i, f"{value:.3f}", ha="center", va="center", color="white" if value > 0.45 else "black", fontsize=8)
            else:
                ax.text(j, i, "NA", ha="center", va="center", color="black", fontsize=8)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label(cbar_label)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_topn_membership(rows: list[dict[str, object]], out: Path) -> None:
    if not rows:
        save_placeholder(out, "Density class topN membership", "No standardized csplit rows available.")
        return
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharey=True)
    top_ns = sorted({as_int(r.get("top_n")) for r in rows})
    for ax, (pipeline, fam) in zip(axes.ravel(), [(p, f) for p in ("P8", "P10") for f in FEATURE_FAMILIES]):
        subset = [r for r in rows if r.get("pipeline") == pipeline and r.get("feature_family") == fam]
        bottoms = [0.0] * len(top_ns)
        for cls in DENSITY_CLASSES:
            vals = []
            for top_n in top_ns:
                these = [as_float(r.get("fraction")) for r in subset if as_int(r.get("top_n")) == top_n and r.get("density_class") == cls]
                med = safe_median(these)
                vals.append(float(med) if med != "" else 0.0)
            ax.bar([str(n) for n in top_ns], vals, bottom=bottoms, label=cls)
            bottoms = [b + v for b, v in zip(bottoms, vals)]
        ax.set_title(f"{pipeline} {fam}")
        ax.set_xlabel("topN")
        ax.set_ylim(0, 1.02)
        ax.grid(axis="y", alpha=0.2)
    axes[0, 0].set_ylabel("median fraction")
    axes[1, 0].set_ylabel("median fraction")
    axes[0, 1].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig.suptitle("Density class membership inside global topN")
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def rank_topn_subset(rows: list[dict[str, object]], top_n: int) -> list[dict[str, object]]:
    out = []
    for row in rows:
        rank = as_int(row.get("rank_within_group"))
        if 1 <= rank <= top_n:
            out.append(row)
    return out


def plot_dimred_label_fraction(rows: list[dict[str, object]], out: Path, top_n: int) -> None:
    if not rows:
        save_placeholder(out, f"Dimred label fraction top{top_n}", "No dimred standardized csplit rows available.")
        return
    contexts = [(p, f, obs) for p in ("P8", "P10") for f in FEATURE_FAMILIES for obs in sorted({str(r.get("observable_mode")) for r in rows})]
    contexts = [ctx for ctx in contexts if any(r.get("pipeline") == ctx[0] and r.get("feature_family") == ctx[1] and r.get("observable_mode") == ctx[2] for r in rows)]
    if not contexts:
        save_placeholder(out, f"Dimred label fraction top{top_n}", "No contexts available.")
        return
    labels = [lab for lab in LABEL_ORDER if any((r.get("two_subprocess_label") or "label_missing") == lab for r in rows)]
    if not labels:
        labels = ["label_missing"]
    fig, ax = plt.subplots(figsize=(12, max(4, 0.45 * len(contexts))))
    y = list(range(len(contexts)))
    left = [0.0] * len(contexts)
    for label in labels:
        vals = []
        for ctx in contexts:
            subset = [
                r
                for r in rows
                if r.get("pipeline") == ctx[0] and r.get("feature_family") == ctx[1] and r.get("observable_mode") == ctx[2]
            ]
            subset = sorted(subset, key=lambda r: as_float(r.get("peak_abs_corr")), reverse=True)[:top_n]
            denom = max(1, len(subset))
            vals.append(sum(1 for r in subset if (r.get("two_subprocess_label") or "label_missing") == label) / denom)
        ax.barh(y, vals, left=left, label=label)
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks(y, [f"{p} | {f} | {obs}" for p, f, obs in contexts])
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction")
    ax.set_title(f"Dimred BLP component label fraction in top{top_n}")
    ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_dimred_label_fraction_by_topn(rows: list[dict[str, object]], fig_dir: Path, top_ns: Sequence[int]) -> None:
    if not rows:
        save_placeholder(fig_dir / "dimred_label_fraction_by_topN__missing.png", "Dimred label fraction by topN", "No dimred standardized csplit rows available.")
        return
    labels = [lab for lab in LABEL_ORDER if any((r.get("two_subprocess_label") or "label_missing") == lab for r in rows)]
    if not labels:
        labels = ["label_missing"]
    for pipeline in ("P8", "P10"):
        for fam in FEATURE_FAMILIES:
            subset = [r for r in rows if r.get("pipeline") == pipeline and r.get("feature_family") == fam]
            out = fig_dir / f"dimred_label_fraction_by_topN__{pipeline.lower()}__{fam}.png"
            if not subset:
                save_placeholder(out, f"{pipeline} {fam} dimred labels by topN", "No rows available.")
                continue
            matrix: list[list[float]] = []
            for label in labels:
                row_vals = []
                for top_n in top_ns:
                    selected = rank_topn_subset(subset, top_n)
                    denom = max(1, len(selected))
                    row_vals.append(sum(1 for r in selected if (r.get("two_subprocess_label") or "label_missing") == label) / denom)
                matrix.append(row_vals)
            plot_fraction_heatmap(
                matrix,
                labels,
                [f"top{n}" for n in top_ns],
                out,
                f"{pipeline} {fam}: dimred label fraction by topN",
            )


def plot_fraction_heatmap(matrix: list[list[float]], ylabels: list[str], xlabels: list[str], out: Path, title: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    if not matrix:
        save_placeholder(out, title, "No data.")
        return
    fig, ax = plt.subplots(figsize=(max(7, 1.25 * len(xlabels)), max(3.6, 0.45 * len(ylabels) + 1.8)))
    im = ax.imshow(matrix, aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(len(xlabels)), xlabels)
    ax.set_yticks(range(len(ylabels)), ylabels)
    ax.set_title(title)
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white" if value > 0.45 else "black", fontsize=8)
    cbar = fig.colorbar(im, ax=ax, shrink=0.82)
    cbar.set_label("fraction")
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_method_k(rows: list[dict[str, object]], fig_dir: Path) -> None:
    if not rows:
        save_placeholder(fig_dir / "method_k_best_score__missing.png", "Method-k best score", "No dimred standardized csplit rows available.")
        return
    for pipeline in ("P8", "P10"):
        for fam in FEATURE_FAMILIES:
            subset = [r for r in rows if r.get("pipeline") == pipeline and r.get("feature_family") == fam]
            out = fig_dir / f"method_k_best_score__{pipeline.lower()}__{fam}.png"
            matrix: list[list[float]] = []
            for method in METHODS:
                vals_row = []
                for k in KS:
                    vals = [
                        as_float(r.get("peak_abs_corr"))
                        for r in subset
                        if r.get("blp_method") == method and as_int(r.get("blp_k")) == k
                    ]
                    med = safe_median(vals)
                    vals_row.append(float(med) if med != "" else math.nan)
                matrix.append(vals_row)
            plot_heatmap(matrix, list(METHODS), [f"k{k:02d}" for k in KS], out, f"{pipeline} {fam}: dimred method-k best |xcorr|", "median best |r|")


def plot_raw_timescale_ecdf(rows: list[dict[str, object]], out: Path, top_ns: Sequence[int]) -> None:
    with_tau = [r for r in rows if as_float(r.get("timescale_sec_preferred")) > 0]
    if not with_tau:
        save_placeholder(out, "Raw efun timescale ECDF", "No raw efun timescale metadata was available for standardized csplit top hits.")
        return
    fig, axes = plt.subplots(2, len(top_ns), figsize=(4.7 * len(top_ns), 7), sharey=True)
    for i, pipeline in enumerate(("P8", "P10")):
        for j, top_n in enumerate(top_ns):
            ax = axes[i, j]
            for fam in FEATURE_FAMILIES:
                subset = [r for r in rows if r.get("pipeline") == pipeline and r.get("feature_family") == fam]
                subset = rank_topn_subset(subset, top_n)
                vals = sorted(math.log10(as_float(r.get("timescale_sec_preferred"))) for r in subset if as_float(r.get("timescale_sec_preferred")) > 0)
                if vals:
                    y = [(idx + 1) / len(vals) for idx in range(len(vals))]
                    ax.step(vals, y, where="post", label=f"{fam} n={len(vals)}")
            ax.set_title(f"{pipeline} top{top_n}")
            ax.set_xlabel("log10(timescale sec)")
            ax.grid(alpha=0.2)
            if j == 0:
                ax.set_ylabel("ECDF")
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Raw BLP efun preferred timescale: efun vs deconv_efun")
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_raw_timescale_summary_by_topn(rows: list[dict[str, object]], out: Path, top_ns: Sequence[int]) -> None:
    with_tau = [r for r in rows if as_float(r.get("timescale_sec_preferred")) > 0]
    if not with_tau:
        save_placeholder(out, "Raw efun timescale by topN", "No raw efun timescale metadata was available for standardized csplit top hits.")
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    markers = {"efun": "o", "deconv_efun": "s"}
    for pipeline in ("P8", "P10"):
        for fam in FEATURE_FAMILIES:
            medians = []
            slow_fracs = []
            for top_n in top_ns:
                subset = [
                    r for r in rows
                    if r.get("pipeline") == pipeline and r.get("feature_family") == fam
                ]
                subset = rank_topn_subset(subset, top_n)
                vals = [as_float(r.get("timescale_sec_preferred")) for r in subset if as_float(r.get("timescale_sec_preferred")) > 0]
                if vals:
                    log_vals = sorted(math.log10(v) for v in vals)
                    medians.append(log_vals[len(log_vals) // 2])
                    slow_fracs.append(sum(1 for v in vals if v >= 1.0) / len(vals))
                else:
                    medians.append(math.nan)
                    slow_fracs.append(math.nan)
            label = f"{pipeline} {fam}"
            axes[0].plot(top_ns, medians, marker=markers[fam], label=label)
            axes[1].plot(top_ns, slow_fracs, marker=markers[fam], label=label)
    axes[0].set_title("Median raw efun timescale in topN")
    axes[0].set_ylabel("median log10(timescale sec)")
    axes[1].set_title("Slow-mode fraction in topN")
    axes[1].set_ylabel("fraction with timescale >= 1 sec")
    for ax in axes:
        ax.set_xlabel("topN per observable context")
        ax.set_xticks(list(top_ns))
        ax.grid(alpha=0.25)
    axes[1].set_ylim(-0.02, 1.02)
    axes[1].legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def write_summary(
    args: argparse.Namespace,
    availability: list[dict[str, object]],
    missing: list[dict[str, object]],
    top_rows: list[dict[str, object]],
    best_rows: list[dict[str, object]],
    output_dir: Path,
) -> None:
    lines = [
        f"# Pipeline12 standardized csplit consistency: {args.dataset}",
        "",
        f"- dataset: `{args.dataset}`",
        f"- condition tag: `{args.condition_tag}`",
        f"- P8 save tag: `{args.p8_save_tag}`",
        f"- P10 save tag: `{args.p10_save_tag}`",
        f"- top rows retained: `{len(top_rows)}`",
        f"- density-class best rows: `{len(best_rows)}`",
        "",
        "## Availability",
        "",
    ]
    for row in availability:
        lines.append(f"- {row['item']}: `{row['status']}` - {row['detail']}")
    if missing:
        lines.extend(["", "## Blocking / missing requirements", ""])
        for row in missing:
            lines.append(f"- {row['item']}: `{row['status']}` - {row['detail']}")
    else:
        lines.extend(["", "## Blocking / missing requirements", "", "- none"])
    if best_rows:
        lines.extend(["", "## Density class best-score snapshot", ""])
        grouped: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
        for row in best_rows:
            grouped[(str(row.get("pipeline")), str(row.get("feature_family")))].update([str(row.get("density_class"))])
        for key in sorted(grouped):
            counts = grouped[key]
            lines.append(f"- {key[0]} {key[1]}: " + ", ".join(f"{c}={counts[c]}" for c in DENSITY_CLASSES))
    else:
        lines.extend(["", "## Density class best-score snapshot", "", "- no standardized csplit xcorr rows available"])
    (output_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir or (DEFAULT_RESULT_ROOT / args.dataset)
    fig_dir = args.figure_dir or (output_dir / "figures")
    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    label_maps = load_labels(args.label_table)
    label_count = len(label_maps[0]) + len(label_maps[1])
    raw_metadata = load_raw_metadata(args.raw_metadata_csv)
    collected = collect_rows(args, label_maps, raw_metadata)
    availability, missing = audit_inputs(args, label_count, len(raw_metadata), collected["file_rows"], collected["top_rows"])

    top_rows = collected["top_rows"]
    best_rows = density_class_best_scores(collected["top_class_rows"])
    membership_rows = topn_membership(top_rows, args.top_n_values)

    write_csv(output_dir / "p12_file_rows_read.csv", collected["file_rows"])
    write_csv(output_dir / "availability.csv", availability)
    write_csv(output_dir / "missing_requirements.csv", missing)
    write_csv(output_dir / "p12_top_rows.csv", top_rows)
    write_csv(output_dir / "p12_density_class_best_score.csv", best_rows)
    write_csv(output_dir / "p12_density_class_topN_membership.csv", membership_rows)
    write_csv(output_dir / "p12_dimred_top_rows_with_labels.csv", collected["top_dimred_rows"])
    write_csv(output_dir / "p12_raw_top_rows_with_timescale.csv", collected["top_raw_rows"])
    write_csv(output_dir / "p12_method_k_best_score.csv", collected["method_k_rows"])

    plot_availability(availability, fig_dir / "availability.png")
    plot_best_score_heatmaps(best_rows, fig_dir)
    plot_topn_membership(membership_rows, fig_dir / "density_class_topN_membership.png")
    for top_n in args.top_n_values:
        plot_dimred_label_fraction(collected["top_dimred_rows"], fig_dir / f"dimred_label_fraction_top{top_n:02d}.png", top_n)
    plot_dimred_label_fraction_by_topn(collected["top_dimred_rows"], fig_dir, args.top_n_values)
    plot_method_k(collected["method_k_rows"], fig_dir)
    plot_raw_timescale_ecdf(collected["top_raw_rows"], fig_dir / "raw_timescale_ecdf_by_topN.png", args.top_n_values)
    plot_raw_timescale_summary_by_topn(collected["top_raw_rows"], fig_dir / "raw_timescale_summary_by_topN.png", args.top_n_values)
    write_summary(args, availability, missing, top_rows, best_rows, output_dir)

    print(f"Pipeline12 output: {output_dir}")
    print(f"Pipeline12 figures: {fig_dir}")
    print(f"Top rows: {len(top_rows)}")
    print(f"Missing requirements: {len(missing)}")


if __name__ == "__main__":
    main()
