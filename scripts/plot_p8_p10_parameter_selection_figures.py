#!/usr/bin/env python3
"""Plot P8/P10 parameter-selection figures from current xcorr peak tables.

This script is figure-first.  It reads full ``*_peaks*.csv`` files, rebuilds
topN views for multiple N values, and writes a compact set of PNG figures for
available/current data only.  Missing datasets/conditions are shown in
availability figures rather than being silently filled.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

try:
    from PIL import Image, ImageDraw, ImageFont
except Exception:  # pragma: no cover
    Image = None
    ImageDraw = None
    ImageFont = None

from summarize_pipeline8_cross_session_consistency import (
    COMPONENT_COUNTS,
    METHOD_ORDER,
    RUN_TAG_LABELS,
    as_float,
    density_group_name,
    density_metadata_fields,
    feature_family,
    plot_heatmap,
    slug,
)


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_RESULTS_DIR = Path("results") / "p8_p10_parameter_selection_current_rmsenv_adaptive_v2"
DEFAULT_RAW_METADATA_CSV = (
    Path("results")
    / "p8_p10_parameter_selection_current_rmsenv_adaptive_v2"
    / "p5_raw_density_mode_metadata.csv"
)
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p8_p10_parameter_selection_v2"
)
DEFAULT_LABEL_TABLE = (
    Path("results")
    / "pipeline5_dimred_component_process_labels_current_rmsenv_adaptive"
    / "dimred_efun_process_labels.csv"
)
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
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")
DEFAULT_TOP_N = (1, 3, 5, 10, 20, 50)
DEFAULT_METHOD_TOP_K = (1, 3, 5, 10, 20)
LABEL_ORDER = (
    "theta_selective_similar",
    "theta_selective_unequal",
    "ripple_selective_similar",
    "ripple_selective_unequal",
    "theta_ripple_joint",
    "mixed_theta_ripple",
    "mixed_theta_gamma",
    "mixed_ripple_gamma",
    "gamma_selective",
    "pan_event",
    "partial_or_inactive",
    "no_event_activity",
    "label_missing",
    "nonselective",
)
DENSITY_CLASS_ORDER = ("event_density", "raw_efun_density", "dimred_efun_density")
DENSITY_CLASS_LABEL = {
    "event_density": "event",
    "raw_efun_density": "raw efun",
    "dimred_efun_density": "dimred efun",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--label-table", type=Path, default=DEFAULT_LABEL_TABLE)
    parser.add_argument("--raw-metadata-csv", type=Path, default=DEFAULT_RAW_METADATA_CSV)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--top-n-values", nargs="+", type=int, default=list(DEFAULT_TOP_N))
    parser.add_argument("--method-top-k-values", nargs="+", type=int, default=list(DEFAULT_METHOD_TOP_K))
    parser.add_argument("--p8-save-tag", default="xcorr_rmsenv_adaptive")
    parser.add_argument("--p10-save-tag", default="dimred_xcorr_rmsenv_adaptive")
    parser.add_argument("--heatmap-top-n", type=int, default=10)
    parser.add_argument("--max-p10-dirs-per-dataset", type=int, default=0)
    return parser.parse_args()


def read_csv_limited(path: Path, limit: int | None = None) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    out: list[dict[str, str]] = []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            out.append(row)
            if limit is not None and len(out) >= limit:
                break
    return out


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


def as_int(value: object, default: int = -1) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def finite(values: Iterable[float]) -> list[float]:
    return [v for v in values if math.isfinite(v)]


def parse_p10_dir(name: str) -> tuple[str, str, str, int, str]:
    match = re.match(r"^(pv_.+?)__(.+?)__([a-z]+)_k(\d+)$", name)
    if not match:
        return "", "", "", -1, ""
    run_tag, p9_feature, method, k = match.groups()
    k_int = int(k)
    return run_tag, p9_feature, method, k_int, f"{method}_k{k_int:02d}"


def normalized_density_name(density_name: str) -> str:
    meta = density_metadata_fields(density_name)
    kind = meta.get("density_source_kind", "")
    if kind == "event":
        return "blp_evt"
    if kind == "raw":
        return f"raw_{meta['density_condition']}_q{meta['density_threshold']}"
    if kind == "dimred":
        return (
            f"dim_{meta['density_condition']}_"
            f"{meta['density_method']}{int(meta['density_k'])}_"
            f"q{meta['density_threshold']}"
        )
    return re.sub(r"_rmsenv_adaptive$", "", density_name)


def load_label_map(path: Path) -> dict[tuple[str, str, int], dict[str, str]]:
    rows = read_csv_limited(path, None)
    labels: dict[tuple[str, str, int], dict[str, str]] = {}
    for row in rows:
        dataset = str(row.get("dataset", "")).lower()
        density = str(row.get("density_name_for_p8_p10") or row.get("density_name") or "")
        idx = as_int(row.get("density_index") or row.get("component_idx"))
        labels[(dataset, density, idx)] = row
    return labels


def normalize_path_for_key(path: object) -> str:
    return str(path or "").replace("/", "\\").lower()


def load_raw_metadata(path: Path) -> dict[tuple[str, int], dict[str, str]]:
    rows = read_csv_limited(path, None)
    out: dict[tuple[str, int], dict[str, str]] = {}
    for row in rows:
        density_file = normalize_path_for_key(row.get("density_file"))
        idx = as_int(row.get("density_index"))
        if density_file and idx >= 0:
            out[(density_file, idx)] = row
    return out


def label_for(labels: dict[tuple[str, str, int], dict[str, str]], row: dict[str, object]) -> str:
    if row.get("density_group") != "dimred_efun_density":
        return ""
    key = (
        str(row.get("dataset", "")).lower(),
        normalized_density_name(str(row.get("density_name", ""))),
        as_int(row.get("density_index")),
    )
    item = labels.get(key, {})
    if not item:
        return "label_missing"
    label = str(item.get("primary_process_label") or "").strip()
    if label == "unlabeled":
        return "no_event_activity"
    return label or "label_missing"


def raw_metadata_for(
    raw_metadata: dict[tuple[str, int], dict[str, str]],
    row: dict[str, str],
) -> dict[str, str]:
    key = (normalize_path_for_key(row.get("density_file")), as_int(row.get("density_index")))
    return raw_metadata.get(key, {})


def condition_scopes(row: dict[str, object]) -> list[str]:
    condition = str(row.get("density_condition", ""))
    if row.get("density_group") == "event_density":
        return ["all", "abs", "csplit"]
    if condition in {"abs", "csplit"}:
        return ["all", condition]
    return ["all"]


def enrich_row(
    *,
    pipeline: str,
    dataset: str,
    run_tag: str,
    source_csv: Path,
    rank_in_source: int,
    row: dict[str, str],
    labels: dict[tuple[str, str, int], dict[str, str]],
    raw_metadata: dict[tuple[str, int], dict[str, str]],
    p9_feature: str = "",
    p9_method: str = "",
    p9_k: int = -1,
    p9_method_k: str = "",
) -> dict[str, object]:
    density = str(row.get("density_name", "")).strip()
    meta = density_metadata_fields(density)
    group = density_group_name(density)
    raw_meta = raw_metadata_for(raw_metadata, row) if group == "raw_efun_density" else {}
    out: dict[str, object] = {
        "pipeline": pipeline,
        "dataset": dataset,
        "run_tag": run_tag,
        "bold_observable": RUN_TAG_LABELS.get(run_tag, run_tag),
        "p9_feature": p9_feature,
        "p9_method": p9_method,
        "p9_k": p9_k if p9_k >= 0 else "",
        "p9_method_k": p9_method_k,
        "source_csv": str(source_csv),
        "rank_in_source_csv": rank_in_source,
        "density_name": density,
        "density_group": group,
        **meta,
        "density_index": as_int(row.get("density_index")),
        "density_label": row.get("density_label", ""),
        "bold_feature": row.get("bold_feature", ""),
        "bold_feature_family": feature_family(str(row.get("bold_feature", ""))),
        "bold_mode_index": row.get("bold_mode_index", ""),
        "bold_component_index": row.get("bold_component_index", ""),
        "peak_abs_corr": as_float(row.get("peak_abs_corr")),
        "peak_corr": as_float(row.get("peak_corr")),
        "peak_lag_sec": as_float(row.get("peak_lag_sec")),
        "raw_efun_index": as_int(raw_meta.get("raw_efun_index") or row.get("density_index")),
        "raw_timescale_sec": as_float(
            raw_meta.get("timescale_sec_preferred") or raw_meta.get("timescale_sec")
        ),
        "raw_frequency_hz": as_float(
            raw_meta.get("frequency_hz_preferred") or raw_meta.get("frequency_hz")
        ),
        "raw_timescale_source": raw_meta.get("timescale_source_preferred", ""),
        "raw_frequency_source": raw_meta.get("frequency_source_preferred", ""),
        "raw_lambda_continuous_real": as_float(raw_meta.get("lambda_continuous_real")),
        "raw_lambda_continuous_imag": as_float(raw_meta.get("lambda_continuous_imag")),
        "raw_envelope_window_sec": as_float(raw_meta.get("envelope_window_sec")),
        "raw_metadata_status": "present" if raw_meta else ("not_raw" if group != "raw_efun_density" else "missing"),
    }
    out["selectivity_label"] = label_for(labels, out)
    return out


def collect_p8(
    args: argparse.Namespace,
    labels: dict[tuple[str, str, int], dict[str, str]],
    raw_metadata: dict[tuple[str, int], dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    hits: list[dict[str, object]] = []
    availability: list[dict[str, object]] = []
    max_n = max(args.top_n_values)
    for dataset in args.datasets:
        for run_tag in args.run_tags:
            root = args.processed_root / dataset / "pipeline8_xcorr" / run_tag
            feature_files = sorted((root / "feature").glob(f"*/*{args.p8_save_tag}_peaks__*.csv"))
            density_files = sorted((root / "density").glob(f"{args.p8_save_tag}_peaks__*.csv"))
            files = feature_files if feature_files else density_files
            status = "available" if files else ("missing" if not root.is_dir() else "missing_peaks")
            availability.append(
                {
                    "pipeline": "P8",
                    "dataset": dataset,
                    "run_tag": run_tag,
                    "bold_observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                    "status": status,
                    "n_peak_csv": len(files),
                }
            )
            for csv_path in files:
                for rank, row in enumerate(read_csv_limited(csv_path, max_n), start=1):
                    hits.append(
                        enrich_row(
                            pipeline="P8",
                            dataset=dataset,
                            run_tag=run_tag,
                            source_csv=csv_path,
                            rank_in_source=rank,
                            row=row,
                            labels=labels,
                            raw_metadata=raw_metadata,
                        )
                    )
    return hits, availability


def collect_p10(
    args: argparse.Namespace,
    labels: dict[tuple[str, str, int], dict[str, str]],
    raw_metadata: dict[tuple[str, int], dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    hits: list[dict[str, object]] = []
    availability: list[dict[str, object]] = []
    max_n = max(args.top_n_values)
    for dataset in args.datasets:
        root = args.processed_root / dataset / "pipeline10_dimred_xcorr"
        p10_dirs = [p for p in sorted(root.iterdir()) if p.is_dir()] if root.is_dir() else []
        if args.max_p10_dirs_per_dataset > 0:
            p10_dirs = p10_dirs[: args.max_p10_dirs_per_dataset]
        seen_run_tags: Counter[str] = Counter()
        for p10_dir in p10_dirs:
            run_tag, p9_feature, p9_method, p9_k, p9_method_k = parse_p10_dir(p10_dir.name)
            if run_tag not in args.run_tags:
                continue
            feature_files = sorted((p10_dir / "feature").glob(f"*/*{args.p10_save_tag}_peaks__*.csv"))
            density_files = sorted((p10_dir / "density").glob(f"{args.p10_save_tag}_peaks__*.csv"))
            files = feature_files if feature_files else density_files
            if files:
                seen_run_tags[run_tag] += 1
            for csv_path in files:
                for rank, row in enumerate(read_csv_limited(csv_path, max_n), start=1):
                    hits.append(
                        enrich_row(
                            pipeline="P10",
                            dataset=dataset,
                            run_tag=run_tag,
                            source_csv=csv_path,
                            rank_in_source=rank,
                            row=row,
                            labels=labels,
                            raw_metadata=raw_metadata,
                            p9_feature=p9_feature,
                            p9_method=p9_method,
                            p9_k=p9_k,
                            p9_method_k=p9_method_k,
                        )
                    )
        for run_tag in args.run_tags:
            status = "available" if seen_run_tags[run_tag] else ("missing" if not root.is_dir() else "missing_peaks")
            availability.append(
                {
                    "pipeline": "P10",
                    "dataset": dataset,
                    "run_tag": run_tag,
                    "bold_observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                    "status": status,
                    "n_p10_dirs_with_peaks": seen_run_tags[run_tag],
                }
            )
    return hits, availability


def context_key(row: dict[str, object], condition_scope: str) -> tuple[object, ...]:
    return (
        row["pipeline"],
        row["dataset"],
        row["run_tag"],
        row["bold_feature_family"],
        condition_scope,
        row.get("p9_method_k", ""),
    )


def build_topn_contexts(hits: Sequence[dict[str, object]], top_n_values: Sequence[int]) -> dict[int, dict[tuple[object, ...], list[dict[str, object]]]]:
    max_n = max(top_n_values)
    grouped: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in hits:
        for scope in condition_scopes(row):
            grouped[context_key(row, scope)].append(row)
    out: dict[int, dict[tuple[object, ...], list[dict[str, object]]]] = {n: {} for n in top_n_values}
    for key, rows in grouped.items():
        rows_sorted = sorted(rows, key=lambda r: (r.get("peak_abs_corr") if math.isfinite(float(r.get("peak_abs_corr", math.nan))) else -1), reverse=True)
        rows_sorted = rows_sorted[:max_n]
        ranked_rows = []
        for rank, row in enumerate(rows_sorted, start=1):
            row_i = dict(row)
            row_i["rank_in_context"] = rank
            ranked_rows.append(row_i)
        for n in top_n_values:
            out[n][key] = ranked_rows[: min(n, len(ranked_rows))]
    return out


def density_membership_summary(contexts_by_n: dict[int, dict[tuple[object, ...], list[dict[str, object]]]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for top_n, contexts in sorted(contexts_by_n.items()):
        grouped: dict[tuple[str, str, str, str], Counter[str]] = defaultdict(Counter)
        totals: Counter[tuple[str, str, str, str]] = Counter()
        for key, hits in contexts.items():
            pipeline, _dataset, _run_tag, family, condition_scope, _p9 = key
            group_key = (str(pipeline), str(family), str(condition_scope), "all_p9")
            for row in hits:
                group = str(row.get("density_group"))
                if group in DENSITY_CLASS_ORDER:
                    grouped[group_key][group] += 1
                    totals[group_key] += 1
        for group_key, counts in grouped.items():
            pipeline, family, condition_scope, p9 = group_key
            total = max(totals[group_key], 1)
            for density_group in DENSITY_CLASS_ORDER:
                rows.append(
                    {
                        "pipeline": pipeline,
                        "bold_feature_family": family,
                        "condition_scope": condition_scope,
                        "p9_method_k": p9,
                        "top_n": top_n,
                        "density_group": density_group,
                        "count": counts[density_group],
                        "fraction": counts[density_group] / total,
                        "n_context_hits": total,
                    }
                )
    return rows


def winner_summary(contexts: dict[tuple[object, ...], list[dict[str, object]]], top_n: int) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], Counter[str]] = defaultdict(Counter)
    dataset_sets: dict[tuple[str, str, str, str], set[str]] = defaultdict(set)
    for key, rows in contexts.items():
        if not rows:
            continue
        pipeline, dataset, _run_tag, family, condition_scope, _p9 = key
        winner = str(rows[0].get("density_group"))
        if winner not in DENSITY_CLASS_ORDER:
            continue
        group_key = (str(pipeline), str(family), str(condition_scope), "all_p9")
        grouped[group_key][winner] += 1
        dataset_sets[group_key].add(str(dataset))
    out: list[dict[str, object]] = []
    for group_key, counts in sorted(grouped.items()):
        pipeline, family, condition_scope, p9 = group_key
        total = sum(counts.values())
        for density_group in DENSITY_CLASS_ORDER:
            out.append(
                {
                    "pipeline": pipeline,
                    "bold_feature_family": family,
                    "condition_scope": condition_scope,
                    "p9_method_k": p9,
                    "top_n": top_n,
                    "density_group": density_group,
                    "winner_count": counts[density_group],
                    "winner_fraction": counts[density_group] / total if total else math.nan,
                    "n_available_contexts": total,
                    "n_available_datasets": len(dataset_sets[group_key]),
                    "datasets": ";".join(sorted(dataset_sets[group_key])),
                }
            )
    return out


def method_k_robustness(contexts: dict[tuple[object, ...], list[dict[str, object]]], pipeline_filter: str) -> list[dict[str, object]]:
    context_totals: Counter[tuple[str, str]] = Counter()
    context_has_mk: dict[tuple[str, str, str], set[tuple[object, ...]]] = defaultdict(set)
    for key, rows in contexts.items():
        pipeline, _dataset, _run_tag, family, condition_scope, _p9 = key
        if pipeline != pipeline_filter or condition_scope not in {"abs", "csplit"}:
            continue
        group_key = (str(family), str(condition_scope))
        context_totals[group_key] += 1
        for row in rows:
            if row.get("density_group") != "dimred_efun_density":
                continue
            method_k = str(row.get("density_method_k", ""))
            if method_k:
                context_has_mk[(str(family), str(condition_scope), method_k)].add(key)
    out: list[dict[str, object]] = []
    for (family, condition, method_k), keys in sorted(context_has_mk.items()):
        total = context_totals[(family, condition)]
        out.append(
            {
                "pipeline": pipeline_filter,
                "bold_feature_family": family,
                "condition_scope": condition,
                "method_k": method_k,
                "context_count": len(keys),
                "n_available_contexts": total,
                "fraction": len(keys) / total if total else math.nan,
            }
        )
    return out


def blp_method_k_v2_stats(
    contexts_by_n: dict[int, dict[tuple[object, ...], list[dict[str, object]]]],
    pipeline_filter: str,
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for top_n, contexts in sorted(contexts_by_n.items()):
        denominators: Counter[tuple[str, str]] = Counter()
        membership: dict[tuple[str, str, str], set[tuple[object, ...]]] = defaultdict(set)
        winners: Counter[tuple[str, str, str]] = Counter()
        best_ranks: dict[tuple[str, str, str], list[int]] = defaultdict(list)
        for key, rows in contexts.items():
            pipeline, _dataset, _run_tag, family, condition_scope, _p9 = key
            if pipeline != pipeline_filter or condition_scope not in {"abs", "csplit"}:
                continue
            dimred_rows = [r for r in rows if r.get("density_group") == "dimred_efun_density"]
            if not dimred_rows:
                continue
            group_key = (str(family), str(condition_scope))
            denominators[group_key] += 1
            seen: dict[str, int] = {}
            for row in dimred_rows:
                method_k = str(row.get("density_method_k", ""))
                if not method_k:
                    continue
                rank = as_int(row.get("rank_in_context"), 999999)
                if method_k not in seen or rank < seen[method_k]:
                    seen[method_k] = rank
            if not seen:
                continue
            winner_mk = min(seen.items(), key=lambda kv: kv[1])[0]
            winners[(str(family), str(condition_scope), winner_mk)] += 1
            for method_k, rank in seen.items():
                membership[(str(family), str(condition_scope), method_k)].add(key)
                best_ranks[(str(family), str(condition_scope), method_k)].append(rank)
        for (family, condition, method_k), keys in sorted(membership.items()):
            denom = denominators[(family, condition)]
            ranks = best_ranks[(family, condition, method_k)]
            out.append(
                {
                    "pipeline": pipeline_filter,
                    "bold_feature_family": family,
                    "condition_scope": condition,
                    "top_n": top_n,
                    "method_k": method_k,
                    "membership_context_count": len(keys),
                    "winner_context_count": winners[(family, condition, method_k)],
                    "n_available_contexts": denom,
                    "membership_fraction": len(keys) / denom if denom else math.nan,
                    "winner_fraction": winners[(family, condition, method_k)] / denom if denom else math.nan,
                    "median_best_rank": median(ranks) if ranks else math.nan,
                }
            )
    return out


def p10_bold_method_k_robustness(contexts: dict[tuple[object, ...], list[dict[str, object]]]) -> list[dict[str, object]]:
    context_totals: Counter[tuple[str, str]] = Counter()
    has_p9: dict[tuple[str, str, str], set[tuple[object, ...]]] = defaultdict(set)
    for key, rows in contexts.items():
        pipeline, _dataset, _run_tag, family, condition_scope, p9_method_k = key
        if pipeline != "P10" or condition_scope not in {"all", "abs", "csplit"}:
            continue
        group_key = (str(family), str(condition_scope))
        context_totals[group_key] += 1
        if p9_method_k:
            has_p9[(str(family), str(condition_scope), str(p9_method_k))].add(key)
    out: list[dict[str, object]] = []
    for (family, condition, p9_method_k), keys in sorted(has_p9.items()):
        total = context_totals[(family, condition)]
        out.append(
            {
                "pipeline": "P10",
                "bold_feature_family": family,
                "condition_scope": condition,
                "p9_method_k": p9_method_k,
                "context_count": len(keys),
                "n_available_contexts": total,
                "fraction": len(keys) / total if total else math.nan,
            }
        )
    return out


def p10_bold_method_k_v2_stats(
    contexts_by_n: dict[int, dict[tuple[object, ...], list[dict[str, object]]]],
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for top_n, contexts in sorted(contexts_by_n.items()):
        parent_scores: dict[tuple[object, ...], dict[str, float]] = defaultdict(dict)
        for key, rows in contexts.items():
            pipeline, dataset, run_tag, family, condition_scope, p9_method_k = key
            if pipeline != "P10" or not p9_method_k:
                continue
            score_vals = finite(float(r.get("peak_abs_corr", math.nan)) for r in rows)
            if not score_vals:
                continue
            parent_key = (dataset, run_tag, family, condition_scope)
            parent_scores[parent_key][str(p9_method_k)] = max(score_vals)

        denominators: Counter[tuple[str, str]] = Counter()
        present: dict[tuple[str, str, str], set[tuple[object, ...]]] = defaultdict(set)
        winners: Counter[tuple[str, str, str]] = Counter()
        ranks_by_method: dict[tuple[str, str, str], list[int]] = defaultdict(list)
        scores_by_method: dict[tuple[str, str, str], list[float]] = defaultdict(list)
        for parent_key, score_by_method in parent_scores.items():
            _dataset, _run_tag, family, condition_scope = parent_key
            if condition_scope not in {"all", "abs", "csplit"}:
                continue
            ranked = sorted(score_by_method.items(), key=lambda kv: kv[1], reverse=True)
            if not ranked:
                continue
            denom_key = (str(family), str(condition_scope))
            denominators[denom_key] += 1
            for rank, (method_k, _score) in enumerate(ranked, start=1):
                stat_key = (str(family), str(condition_scope), method_k)
                present[stat_key].add(parent_key)
                ranks_by_method[stat_key].append(rank)
                scores_by_method[stat_key].append(_score)
                if rank == 1:
                    winners[stat_key] += 1
        for (family, condition, method_k), ranks in sorted(ranks_by_method.items()):
            denom = denominators[(family, condition)]
            scores = scores_by_method[(family, condition, method_k)]
            out.append(
                {
                    "pipeline": "P10",
                    "bold_feature_family": family,
                    "condition_scope": condition,
                    "top_n": top_n,
                    "p9_method_k": method_k,
                    "membership_context_count": len(present[(family, condition, method_k)]),
                    "winner_context_count": winners[(family, condition, method_k)],
                    "n_available_contexts": denom,
                    "membership_fraction": len(present[(family, condition, method_k)]) / denom if denom else math.nan,
                    "winner_fraction": winners[(family, condition, method_k)] / denom if denom else math.nan,
                    "median_best_rank": median(ranks) if ranks else math.nan,
                    "median_score": median(scores) if scores else math.nan,
                }
            )
    return out


def p10_bold_method_k_rank_stats(
    contexts_by_n: dict[int, dict[tuple[object, ...], list[dict[str, object]]]],
    density_top_n: int,
    method_top_k_values: Sequence[int],
) -> list[dict[str, object]]:
    """Rank P9/BOLD method-k candidates using P10 xcorr strength.

    ``density_top_n`` controls how many density hits are used inside each P10
    run to define that run's score. ``method_top_k_values`` controls the rank
    threshold across P9/BOLD method-k candidates.
    """
    contexts = contexts_by_n.get(density_top_n) or contexts_by_n[max(contexts_by_n)]
    parent_scores: dict[tuple[object, ...], dict[str, float]] = defaultdict(dict)
    for key, rows in contexts.items():
        pipeline, dataset, run_tag, family, condition_scope, p9_method_k = key
        if pipeline != "P10" or not p9_method_k or not is_current_method_k(p9_method_k):
            continue
        score_vals = finite(as_float(r.get("peak_abs_corr")) for r in rows)
        if not score_vals:
            continue
        parent_key = (dataset, run_tag, family, condition_scope)
        parent_scores[parent_key][str(p9_method_k)] = max(score_vals)

    denominators: Counter[tuple[str, str]] = Counter()
    method_counts: dict[tuple[str, str], list[int]] = defaultdict(list)
    winners: Counter[tuple[str, str, str]] = Counter()
    ranks_by_method: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    scores_by_method: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    topk_hits: dict[tuple[str, str, str, int], set[tuple[object, ...]]] = defaultdict(set)
    present: dict[tuple[str, str, str], set[tuple[object, ...]]] = defaultdict(set)

    for parent_key, score_by_method in parent_scores.items():
        _dataset, _run_tag, family, condition_scope = parent_key
        if condition_scope not in {"all", "abs", "csplit"}:
            continue
        ranked = sorted(score_by_method.items(), key=lambda kv: kv[1], reverse=True)
        if not ranked:
            continue
        denom_key = (str(family), str(condition_scope))
        denominators[denom_key] += 1
        method_counts[denom_key].append(len(ranked))
        for rank, (method_k, score) in enumerate(ranked, start=1):
            stat_key = (str(family), str(condition_scope), method_k)
            present[stat_key].add(parent_key)
            ranks_by_method[stat_key].append(rank)
            scores_by_method[stat_key].append(score)
            if rank == 1:
                winners[stat_key] += 1
            for method_top_k in method_top_k_values:
                if rank <= method_top_k:
                    topk_hits[(str(family), str(condition_scope), method_k, int(method_top_k))].add(parent_key)

    out: list[dict[str, object]] = []
    for (family, condition, method_k), ranks in sorted(ranks_by_method.items()):
        denom = denominators[(family, condition)]
        scores = scores_by_method[(family, condition, method_k)]
        available_counts = method_counts[(family, condition)]
        for method_top_k in method_top_k_values:
            topk_count = len(topk_hits[(family, condition, method_k, int(method_top_k))])
            out.append(
                {
                    "pipeline": "P10",
                    "bold_feature_family": family,
                    "condition_scope": condition,
                    "density_top_n": density_top_n,
                    "method_top_k": int(method_top_k),
                    "p9_method_k": method_k,
                    "method_topk_context_count": topk_count,
                    "present_context_count": len(present[(family, condition, method_k)]),
                    "winner_context_count": winners[(family, condition, method_k)],
                    "n_available_contexts": denom,
                    "method_topk_fraction": topk_count / denom if denom else math.nan,
                    "presence_fraction": len(present[(family, condition, method_k)]) / denom if denom else math.nan,
                    "winner_fraction": winners[(family, condition, method_k)] / denom if denom else math.nan,
                    "median_method_rank": median(ranks) if ranks else math.nan,
                    "median_score": median(scores) if scores else math.nan,
                    "median_available_method_count": median(available_counts) if available_counts else math.nan,
                }
            )
    return out


def label_composition(contexts: dict[tuple[object, ...], list[dict[str, object]]], pipeline_filter: str) -> list[dict[str, object]]:
    counts: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    for key, rows in contexts.items():
        pipeline, _dataset, _run_tag, family, condition_scope, _p9 = key
        if pipeline != pipeline_filter or condition_scope not in {"abs", "csplit"}:
            continue
        for row in rows:
            if row.get("density_group") != "dimred_efun_density":
                continue
            method_k = str(row.get("density_method_k", ""))
            label = str(row.get("selectivity_label") or "unlabeled")
            counts[(str(family), str(condition_scope), method_k)][label] += 1
    out: list[dict[str, object]] = []
    for (family, condition, method_k), counter in sorted(counts.items()):
        total = sum(counter.values())
        for label in LABEL_ORDER:
            out.append(
                {
                    "pipeline": pipeline_filter,
                    "bold_feature_family": family,
                    "condition_scope": condition,
                    "method_k": method_k,
                    "label": label,
                    "count": counter[label],
                    "fraction": counter[label] / total if total else math.nan,
                    "n_dimred_hits": total,
                }
            )
    return out


def raw_index_counts(contexts: dict[tuple[object, ...], list[dict[str, object]]], pipeline_filter: str) -> Counter[int]:
    counts: Counter[int] = Counter()
    for key, rows in contexts.items():
        pipeline = key[0]
        if pipeline != pipeline_filter:
            continue
        for row in rows:
            if row.get("density_group") == "raw_efun_density":
                idx = as_int(row.get("density_index"))
                if idx >= 0:
                    counts[idx] += 1
    return counts


def raw_values_by_family(
    contexts: dict[tuple[object, ...], list[dict[str, object]]],
    pipeline_filter: str,
    value_field: str,
) -> dict[str, list[float]]:
    values: dict[str, list[float]] = {"efun": [], "deconv_efun": []}
    seen: set[tuple[object, ...]] = set()
    for key, rows in contexts.items():
        pipeline = key[0]
        if pipeline != pipeline_filter:
            continue
        for row in rows:
            if row.get("density_group") != "raw_efun_density":
                continue
            family = str(row.get("bold_feature_family", ""))
            if family not in values:
                continue
            value = as_float(row.get(value_field))
            if not math.isfinite(value):
                continue
            unique_key = (
                pipeline_filter,
                row.get("dataset"),
                row.get("run_tag"),
                family,
                row.get("density_name"),
                row.get("density_index"),
                row.get("bold_feature"),
                row.get("bold_mode_index"),
                row.get("peak_lag_sec"),
            )
            if unique_key in seen:
                continue
            seen.add(unique_key)
            values[family].append(value)
    return values


def percentile(sorted_values: Sequence[float], p: float) -> float:
    if not sorted_values:
        return math.nan
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = (len(sorted_values) - 1) * p
    lo = int(math.floor(pos))
    hi = int(math.ceil(pos))
    if lo == hi:
        return sorted_values[lo]
    frac = pos - lo
    return sorted_values[lo] * (1 - frac) + sorted_values[hi] * frac


def raw_timescale_topn_summary(
    contexts_by_n: dict[int, dict[tuple[object, ...], list[dict[str, object]]]],
    pipeline_filter: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    by_topn_family: dict[tuple[int, str], dict[str, float]] = {}
    for top_n, contexts in sorted(contexts_by_n.items()):
        values = raw_values_by_family(contexts, pipeline_filter, "raw_timescale_sec")
        for family in ("efun", "deconv_efun"):
            clean = sorted(v for v in values.get(family, []) if math.isfinite(v) and v > 0)
            med = median(clean) if clean else math.nan
            q25 = percentile(clean, 0.25)
            q75 = percentile(clean, 0.75)
            mean = sum(clean) / len(clean) if clean else math.nan
            med_log = math.log10(med) if math.isfinite(med) and med > 0 else math.nan
            row = {
                "pipeline": pipeline_filter,
                "top_n": top_n,
                "bold_feature_family": family,
                "n_raw_hits_with_timescale": len(clean),
                "median_timescale_sec": med,
                "q25_timescale_sec": q25,
                "q75_timescale_sec": q75,
                "mean_timescale_sec": mean,
                "median_log10_timescale_sec": med_log,
            }
            rows.append(row)
            by_topn_family[(top_n, family)] = row
        efun_log = by_topn_family[(top_n, "efun")]["median_log10_timescale_sec"]
        deconv_log = by_topn_family[(top_n, "deconv_efun")]["median_log10_timescale_sec"]
        rows.append(
            {
                "pipeline": pipeline_filter,
                "top_n": top_n,
                "bold_feature_family": "deconv_minus_efun",
                "n_raw_hits_with_timescale": "",
                "median_timescale_sec": "",
                "q25_timescale_sec": "",
                "q75_timescale_sec": "",
                "mean_timescale_sec": "",
                "median_log10_timescale_sec": "",
                "delta_median_log10_timescale_sec": (
                    deconv_log - efun_log
                    if math.isfinite(float(deconv_log)) and math.isfinite(float(efun_log))
                    else math.nan
                ),
            }
        )
    return rows


def parse_method_k(value: object) -> tuple[str, str]:
    text = str(value or "")
    match = re.match(r"^([a-z]+)_k(\d+)$", text)
    if not match:
        return "", ""
    method, k = match.groups()
    return method, f"k{int(k):02d}"


def method_k_rows() -> list[str]:
    return [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]


def is_current_method_k(value: object) -> bool:
    method, k_label = parse_method_k(value)
    if method not in METHOD_ORDER or not k_label:
        return False
    return as_int(k_label.lstrip("k")) in set(COMPONENT_COUNTS)


def load_font(size: int, bold: bool = False):
    if ImageFont is None:
        return None
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except Exception:
            pass
    return ImageFont.load_default() if ImageFont is not None else None


def draw_text(draw, xy, text: str, font, fill=(30, 30, 30), anchor=None):
    if anchor:
        draw.text(xy, text, fill=fill, font=font, anchor=anchor)
    else:
        draw.text(xy, text, fill=fill, font=font)


def write_text_figure(path: Path, title: str, lines: Sequence[str]) -> None:
    if Image is None:
        return
    width = 1280
    height = max(360, 110 + 28 * len(lines))
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw_text(draw, (24, 24), title, load_font(24, True), fill=(20, 30, 40))
    y = 78
    for line in lines:
        draw_text(draw, (30, y), line, load_font(16, False), fill=(70, 70, 70))
        y += 28
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def plot_topn_lines(rows: Sequence[dict[str, object]], pipeline: str, family: str, condition: str, path: Path) -> None:
    if Image is None:
        return
    subset = [
        r
        for r in rows
        if r["pipeline"] == pipeline
        and r["bold_feature_family"] == family
        and r["condition_scope"] == condition
    ]
    if not subset:
        write_text_figure(path, f"{pipeline} topN density-class membership", ["No available data for this condition."])
        return
    by_group: dict[str, dict[int, float]] = defaultdict(dict)
    for row in subset:
        by_group[str(row["density_group"])][int(row["top_n"])] = float(row["fraction"])
    top_ns = sorted({int(r["top_n"]) for r in subset})
    width, height = 980, 620
    left, top, right, bottom = 82, 80, 32, 72
    plot_w = width - left - right
    plot_h = height - top - bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title = f"{pipeline} density-class topN membership | {family} | {condition}"
    draw_text(draw, (24, 22), title, load_font(22, True), fill=(20, 30, 40))
    draw.line((left, top + plot_h, left + plot_w, top + plot_h), fill=(80, 80, 80), width=2)
    draw.line((left, top, left, top + plot_h), fill=(80, 80, 80), width=2)
    colors = {
        "event_density": (80, 80, 80),
        "raw_efun_density": (52, 126, 184),
        "dimred_efun_density": (77, 175, 74),
    }
    min_x, max_x = min(top_ns), max(top_ns)
    for yv in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = top + plot_h - yv * plot_h
        draw.line((left, y, left + plot_w, y), fill=(230, 230, 230))
        draw_text(draw, (left - 12, y), f"{yv:.2f}", load_font(12), anchor="rm", fill=(80, 80, 80))
    for i, xval in enumerate(top_ns):
        x = left + ((xval - min_x) / max(max_x - min_x, 1)) * plot_w
        draw_text(draw, (x, top + plot_h + 20), str(xval), load_font(12), anchor="mm", fill=(60, 60, 60))
    legend_x = left + plot_w - 260
    legend_y = top + 10
    for j, group in enumerate(DENSITY_CLASS_ORDER):
        pts = []
        for xval in top_ns:
            val = by_group.get(group, {}).get(xval, math.nan)
            if math.isfinite(val):
                x = left + ((xval - min_x) / max(max_x - min_x, 1)) * plot_w
                y = top + plot_h - val * plot_h
                pts.append((x, y))
        color = colors[group]
        if len(pts) >= 2:
            draw.line(pts, fill=color, width=4)
        for x, y in pts:
            draw.ellipse((x - 5, y - 5, x + 5, y + 5), fill=color)
        draw.rectangle((legend_x, legend_y + 22 * j, legend_x + 14, legend_y + 22 * j + 14), fill=color)
        draw_text(draw, (legend_x + 20, legend_y + 22 * j - 1), DENSITY_CLASS_LABEL[group], load_font(13), fill=(40, 40, 40))
    draw_text(draw, (left + plot_w / 2, height - 30), "topN", load_font(14, True), anchor="mm")
    draw_text(draw, (18, top + plot_h / 2), "fraction", load_font(14, True), anchor="lm")
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def plot_bar_fraction(rows: Sequence[dict[str, object]], pipeline: str, family: str, condition: str, path: Path) -> None:
    if Image is None:
        return
    subset = [
        r
        for r in rows
        if r["pipeline"] == pipeline
        and r["bold_feature_family"] == family
        and r["condition_scope"] == condition
    ]
    if not subset:
        write_text_figure(path, f"{pipeline} density-class winner fraction", ["No available data for this condition."])
        return
    vals = {str(r["density_group"]): float(r["winner_fraction"]) for r in subset}
    counts = {str(r["density_group"]): str(r.get("winner_count", "")) for r in subset}
    width, height = 760, 480
    left, top, bottom = 110, 86, 70
    plot_h = height - top - bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw_text(draw, (24, 22), f"{pipeline} winner fraction | {family} | {condition}", load_font(22, True), fill=(20, 30, 40))
    colors = [(80, 80, 80), (52, 126, 184), (77, 175, 74)]
    bar_w = 120
    gap = 55
    for i, group in enumerate(DENSITY_CLASS_ORDER):
        val = vals.get(group, 0.0)
        x0 = left + i * (bar_w + gap)
        y1 = top + plot_h
        y0 = y1 - val * plot_h
        draw.rectangle((x0, y0, x0 + bar_w, y1), fill=colors[i], outline=(60, 60, 60))
        draw_text(draw, (x0 + bar_w / 2, y0 - 18), f"{val:.2f}", load_font(14, True), anchor="mm")
        draw_text(draw, (x0 + bar_w / 2, y1 + 22), DENSITY_CLASS_LABEL[group], load_font(13), anchor="mm")
        draw_text(draw, (x0 + bar_w / 2, y1 + 42), f"n={counts.get(group, '0')}", load_font(12), anchor="mm", fill=(90, 90, 90))
    draw.line((left - 30, top + plot_h, width - 50, top + plot_h), fill=(80, 80, 80), width=2)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def plot_raw_hist(counts: Counter[int], title: str, path: Path) -> None:
    if Image is None:
        return
    if not counts:
        write_text_figure(path, title, ["No raw efun density hits found in current available data."])
        return
    width, height = 1220, 520
    left, top, right, bottom = 72, 72, 32, 66
    plot_w = width - left - right
    plot_h = height - top - bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw_text(draw, (24, 22), title, load_font(22, True), fill=(20, 30, 40))
    max_idx = max(counts)
    max_count = max(counts.values())
    bar_w = max(1, plot_w / max(max_idx, 1))
    draw.line((left, top + plot_h, left + plot_w, top + plot_h), fill=(80, 80, 80), width=2)
    draw.line((left, top, left, top + plot_h), fill=(80, 80, 80), width=2)
    for idx, count in counts.items():
        x0 = left + (idx - 1) * bar_w
        x1 = x0 + max(1, bar_w * 0.9)
        y1 = top + plot_h
        y0 = y1 - (count / max_count) * plot_h
        draw.rectangle((x0, y0, x1, y1), fill=(52, 126, 184))
    draw_text(draw, (left + plot_w / 2, height - 28), "raw efun index", load_font(14, True), anchor="mm")
    draw_text(draw, (24, top + 10), "count", load_font(13), fill=(70, 70, 70))
    draw_text(draw, (left, height - 50), f"unique indices={len(counts)}; max count={max_count}", load_font(12), fill=(90, 90, 90))
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def plot_raw_index_hist_by_family(values: dict[str, list[float]], title: str, path: Path) -> None:
    if Image is None:
        return
    counts_by_family = {
        family: Counter(int(v) for v in vals if math.isfinite(v) and v >= 0)
        for family, vals in values.items()
    }
    if not any(counts_by_family.values()):
        write_text_figure(path, title, ["No raw efun index values found."])
        return
    width, height = 1220, 560
    left, top, right, bottom = 72, 78, 32, 74
    plot_w = width - left - right
    plot_h = height - top - bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw_text(draw, (24, 22), title, load_font(22, True), fill=(20, 30, 40))
    max_idx = max(max(c.keys()) for c in counts_by_family.values() if c)
    max_count = max(max(c.values()) for c in counts_by_family.values() if c)
    bar_w = max(1, plot_w / max(max_idx, 1))
    colors = {"efun": (52, 126, 184), "deconv_efun": (228, 108, 10)}
    offsets = {"efun": -0.18, "deconv_efun": 0.18}
    draw.line((left, top + plot_h, left + plot_w, top + plot_h), fill=(80, 80, 80), width=2)
    draw.line((left, top, left, top + plot_h), fill=(80, 80, 80), width=2)
    for family, counts in counts_by_family.items():
        for idx, count in counts.items():
            x_mid = left + (idx - 0.5) * bar_w + offsets[family] * bar_w
            x0 = x_mid - max(1, bar_w * 0.18)
            x1 = x_mid + max(1, bar_w * 0.18)
            y1 = top + plot_h
            y0 = y1 - (count / max_count) * plot_h
            draw.rectangle((x0, y0, x1, y1), fill=colors[family])
    legend_x = width - 260
    legend_y = top + 8
    for i, family in enumerate(("efun", "deconv_efun")):
        draw.rectangle((legend_x, legend_y + 24 * i, legend_x + 15, legend_y + 24 * i + 15), fill=colors[family])
        draw_text(draw, (legend_x + 22, legend_y + 24 * i - 1), f"{family} n={len(values.get(family, []))}", load_font(13))
    draw_text(draw, (left + plot_w / 2, height - 30), "raw efun index", load_font(14, True), anchor="mm")
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def plot_timescale_ecdf(values: dict[str, list[float]], title: str, path: Path) -> None:
    if Image is None:
        return
    clean = {
        family: sorted(math.log10(v) for v in vals if math.isfinite(v) and v > 0)
        for family, vals in values.items()
    }
    if not any(clean.values()):
        write_text_figure(
            path,
            title,
            [
                "No positive finite timescale_sec values were available after metadata join.",
                "Check p5_raw_density_mode_metadata.csv and raw density_file path matching.",
            ],
        )
        return
    all_vals = [v for vals in clean.values() for v in vals]
    xmin, xmax = min(all_vals), max(all_vals)
    if xmin == xmax:
        xmin -= 0.5
        xmax += 0.5
    width, height = 900, 600
    left, top, right, bottom = 88, 78, 40, 76
    plot_w = width - left - right
    plot_h = height - top - bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw_text(draw, (24, 22), title, load_font(22, True), fill=(20, 30, 40))
    draw.line((left, top + plot_h, left + plot_w, top + plot_h), fill=(80, 80, 80), width=2)
    draw.line((left, top, left, top + plot_h), fill=(80, 80, 80), width=2)
    for yv in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = top + plot_h - yv * plot_h
        draw.line((left, y, left + plot_w, y), fill=(230, 230, 230))
        draw_text(draw, (left - 12, y), f"{yv:.2f}", load_font(12), anchor="rm", fill=(80, 80, 80))
    colors = {"efun": (52, 126, 184), "deconv_efun": (228, 108, 10)}
    for family in ("efun", "deconv_efun"):
        vals = clean.get(family, [])
        if not vals:
            continue
        pts = []
        n = len(vals)
        for i, val in enumerate(vals, start=1):
            x = left + ((val - xmin) / (xmax - xmin)) * plot_w
            y = top + plot_h - (i / n) * plot_h
            pts.append((x, y))
        if len(pts) >= 2:
            draw.line(pts, fill=colors[family], width=4)
        for x, y in pts[:: max(1, len(pts) // 40)]:
            draw.ellipse((x - 3, y - 3, x + 3, y + 3), fill=colors[family])
    legend_x = width - 250
    legend_y = top + 12
    for i, family in enumerate(("efun", "deconv_efun")):
        draw.rectangle((legend_x, legend_y + 24 * i, legend_x + 15, legend_y + 24 * i + 15), fill=colors[family])
        draw_text(draw, (legend_x + 22, legend_y + 24 * i - 1), f"{family} n={len(clean.get(family, []))}", load_font(13))
    draw_text(draw, (left + plot_w / 2, height - 32), "log10(timescale_sec)", load_font(14, True), anchor="mm")
    draw_text(draw, (24, top + plot_h / 2), "ECDF", load_font(14, True), anchor="lm")
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def plot_timescale_topn_median_lines(
    rows: Sequence[dict[str, object]],
    pipeline: str,
    title: str,
    path: Path,
) -> None:
    if Image is None:
        return
    subset = [
        r for r in rows
        if r.get("pipeline") == pipeline and r.get("bold_feature_family") in {"efun", "deconv_efun"}
    ]
    top_ns = sorted({as_int(r.get("top_n")) for r in subset if as_int(r.get("top_n")) >= 0})
    values: dict[str, dict[int, float]] = {"efun": {}, "deconv_efun": {}}
    for row in subset:
        family = str(row.get("bold_feature_family"))
        top_n = as_int(row.get("top_n"))
        value = as_float(row.get("median_log10_timescale_sec"))
        if family in values and top_n >= 0 and math.isfinite(value):
            values[family][top_n] = value
    all_vals = [v for fam_vals in values.values() for v in fam_vals.values() if math.isfinite(v)]
    if not all_vals:
        write_text_figure(path, title, ["No finite raw efun timescale medians found."])
        return
    ymin, ymax = min(all_vals), max(all_vals)
    if ymin == ymax:
        ymin -= 0.5
        ymax += 0.5
    pad = max(0.05, (ymax - ymin) * 0.12)
    ymin -= pad
    ymax += pad

    width, height = 940, 600
    left, top, right, bottom = 92, 80, 42, 78
    plot_w = width - left - right
    plot_h = height - top - bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    draw_text(draw, (24, 22), title, load_font(22, True), fill=(20, 30, 40))
    draw.line((left, top + plot_h, left + plot_w, top + plot_h), fill=(80, 80, 80), width=2)
    draw.line((left, top, left, top + plot_h), fill=(80, 80, 80), width=2)
    min_x, max_x = min(top_ns), max(top_ns)
    for yfrac in (0.0, 0.25, 0.5, 0.75, 1.0):
        value = ymin + yfrac * (ymax - ymin)
        y = top + plot_h - yfrac * plot_h
        draw.line((left, y, left + plot_w, y), fill=(230, 230, 230))
        draw_text(draw, (left - 12, y), f"{value:.2f}", load_font(12), anchor="rm", fill=(70, 70, 70))
    for xval in top_ns:
        x = left + ((xval - min_x) / max(max_x - min_x, 1)) * plot_w
        draw_text(draw, (x, top + plot_h + 20), f"top{xval}", load_font(12), anchor="mm", fill=(70, 70, 70))
    colors = {"efun": (52, 126, 184), "deconv_efun": (228, 108, 10)}
    for family in ("efun", "deconv_efun"):
        pts = []
        for xval in top_ns:
            value = values[family].get(xval, math.nan)
            if not math.isfinite(value):
                continue
            x = left + ((xval - min_x) / max(max_x - min_x, 1)) * plot_w
            y = top + plot_h - ((value - ymin) / (ymax - ymin)) * plot_h
            pts.append((x, y))
        if len(pts) >= 2:
            draw.line(pts, fill=colors[family], width=4)
        for x, y in pts:
            draw.ellipse((x - 5, y - 5, x + 5, y + 5), fill=colors[family])
    legend_x = width - 265
    legend_y = top + 10
    for i, family in enumerate(("efun", "deconv_efun")):
        draw.rectangle((legend_x, legend_y + 24 * i, legend_x + 15, legend_y + 24 * i + 15), fill=colors[family])
        draw_text(draw, (legend_x + 22, legend_y + 24 * i - 1), family, load_font(13))
    draw_text(draw, (left + plot_w / 2, height - 34), "density-hit topN", load_font(14, True), anchor="mm")
    draw_text(draw, (22, top + plot_h / 2), "median log10(timescale_sec)", load_font(14, True), anchor="lm")
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def plot_method_k_topn_membership(
    rows: Sequence[dict[str, object]],
    *,
    method_field: str,
    family: str,
    condition: str,
    title: str,
    path: Path,
) -> None:
    matrix: dict[tuple[str, str], float] = {}
    top_ns = sorted({as_int(r.get("top_n")) for r in rows if as_int(r.get("top_n")) >= 0})
    for row in rows:
        if row.get("bold_feature_family") != family or row.get("condition_scope") != condition:
            continue
        method_k = str(row.get(method_field, ""))
        top_n = as_int(row.get("top_n"))
        if not method_k or top_n < 0:
            continue
        matrix[(method_k, f"top{top_n}")] = as_float(row.get("membership_fraction"))
    plot_heatmap(
        matrix,
        method_k_rows(),
        [f"top{n}" for n in top_ns],
        title,
        path,
        value_format=".2f",
    )


def plot_method_k_method_topk_membership(
    rows: Sequence[dict[str, object]],
    *,
    method_field: str,
    family: str,
    condition: str,
    title: str,
    path: Path,
) -> None:
    matrix: dict[tuple[str, str], float] = {}
    method_top_ks = sorted({as_int(r.get("method_top_k")) for r in rows if as_int(r.get("method_top_k")) >= 0})
    for row in rows:
        if row.get("bold_feature_family") != family or row.get("condition_scope") != condition:
            continue
        method_k = str(row.get(method_field, ""))
        method_top_k = as_int(row.get("method_top_k"))
        if not method_k or method_top_k < 0:
            continue
        matrix[(method_k, f"method_top{method_top_k}")] = as_float(row.get("method_topk_fraction"))
    plot_heatmap(
        matrix,
        method_k_rows(),
        [f"method_top{k}" for k in method_top_ks],
        title,
        path,
        value_format=".2f",
    )


def plot_method_k_grid_metric(
    rows: Sequence[dict[str, object]],
    *,
    method_field: str,
    metric_field: str,
    top_n: int,
    family: str,
    condition: str,
    title: str,
    path: Path,
    value_format: str = ".2f",
) -> None:
    matrix: dict[tuple[str, str], float] = {}
    for row in rows:
        if row.get("bold_feature_family") != family or row.get("condition_scope") != condition:
            continue
        row_top_marker = row.get("top_n") if "top_n" in row else row.get("method_top_k")
        if as_int(row_top_marker) != top_n:
            continue
        method, k_label = parse_method_k(row.get(method_field))
        if not method or not k_label:
            continue
        matrix[(method, k_label)] = as_float(row.get(metric_field))
    plot_heatmap(
        matrix,
        list(METHOD_ORDER),
        [f"k{k:02d}" for k in COMPONENT_COUNTS],
        title,
        path,
        value_format=value_format,
    )


def availability_plot(availability: Sequence[dict[str, object]], datasets: Sequence[str], path: Path) -> None:
    rows = []
    matrix: dict[tuple[str, str], float] = {}
    status_value = {"missing": 0.0, "missing_peaks": 0.25, "available": 1.0}
    for pipeline in ("P8", "P10"):
        for run_tag in DEFAULT_RUN_TAGS:
            row_label = f"{pipeline} {RUN_TAG_LABELS.get(run_tag, run_tag)}"
            rows.append(row_label)
            for item in availability:
                if item["pipeline"] == pipeline and item["run_tag"] == run_tag:
                    matrix[(row_label, str(item["dataset"]))] = status_value.get(str(item["status"]), 0.0)
    plot_heatmap(matrix, rows, list(datasets), "P8/P10 condition availability | 0 missing, 1 available", path, value_format=".0f")


def write_recompute_plan(
    path: Path,
    availability: Sequence[dict[str, object]],
    hits: Sequence[dict[str, object]],
) -> None:
    rows: list[dict[str, object]] = []
    for item in availability:
        if item["status"] != "available":
            rows.append(
                {
                    "dataset": item["dataset"],
                    "pipeline": item["pipeline"],
                    "bold_observable": item["bold_observable"],
                    "status": item["status"],
                    "recommended_action": "run_current_xcorr_or_export_peaks",
                }
            )
    raw_rows = [r for r in hits if r.get("density_group") == "raw_efun_density"]
    raw_present = [r for r in raw_rows if r.get("raw_metadata_status") == "present"]
    raw_missing = [r for r in raw_rows if r.get("raw_metadata_status") != "present"]
    if raw_rows and not raw_present:
        rows.append(
            {
                "dataset": "all_available",
                "pipeline": "P8/P10",
                "bold_observable": "all",
                "status": "missing_joined_raw_timescale_metadata",
                "recommended_action": "export_or_join_P5_raw_density_mode_metadata_before_timescale_ECDF",
            }
        )
    elif raw_missing:
        rows.append(
            {
                "dataset": "all_available",
                "pipeline": "P8/P10",
                "bold_observable": "all",
                "status": "partial_joined_raw_timescale_metadata",
                "recommended_action": "inspect raw_metadata_status counts and missing density_file path matches",
            }
        )
    write_csv(path, rows)


def main() -> int:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    labels = load_label_map(args.label_table)
    raw_metadata = load_raw_metadata(args.raw_metadata_csv)
    p8_hits, p8_availability = collect_p8(args, labels, raw_metadata)
    p10_hits, p10_availability = collect_p10(args, labels, raw_metadata)
    hits = p8_hits + p10_hits
    availability = p8_availability + p10_availability
    contexts_by_n = build_topn_contexts(hits, args.top_n_values)
    heatmap_contexts = contexts_by_n.get(args.heatmap_top_n) or contexts_by_n[max(args.top_n_values)]

    membership = density_membership_summary(contexts_by_n)
    winners = winner_summary(heatmap_contexts, args.heatmap_top_n)
    p8_mk_v2 = blp_method_k_v2_stats(contexts_by_n, "P8")
    p10_blp_mk_v2 = blp_method_k_v2_stats(contexts_by_n, "P10")
    p10_bold_mk_v2 = p10_bold_method_k_rank_stats(
        contexts_by_n,
        args.heatmap_top_n,
        args.method_top_k_values,
    )
    p8_labels = label_composition(heatmap_contexts, "P8")
    p10_labels = label_composition(heatmap_contexts, "P10")
    raw_timescale_summary = (
        raw_timescale_topn_summary(contexts_by_n, "P8")
        + raw_timescale_topn_summary(contexts_by_n, "P10")
    )

    write_csv(args.results_dir / "p8_p10_parameter_selection_hits_long.csv", hits)
    write_csv(args.results_dir / "availability.csv", availability)
    write_csv(args.results_dir / "density_class_topN_membership.csv", membership)
    write_csv(args.results_dir / "density_class_winner_fraction.csv", winners)
    write_csv(args.results_dir / "p8_blp_method_k_v2.csv", p8_mk_v2)
    write_csv(args.results_dir / "p10_blp_method_k_v2.csv", p10_blp_mk_v2)
    write_csv(args.results_dir / "p10_bold_method_k_v2.csv", p10_bold_mk_v2)
    write_csv(args.results_dir / "p8_dimred_selectivity_label_composition.csv", p8_labels)
    write_csv(args.results_dir / "p10_dimred_selectivity_label_composition.csv", p10_labels)
    write_csv(args.results_dir / "raw_efun_timescale_topN_summary.csv", raw_timescale_summary)
    write_recompute_plan(args.results_dir / "recompute_or_metadata_needed.csv", availability, hits)

    availability_plot(availability, args.datasets, args.figure_dir / "00_availability" / "condition_availability_heatmap.png")

    for pipeline, dirname in (("P8", "01_p8_density_class_competition"), ("P10", "02_p10_density_class_competition")):
        for family in ("efun", "deconv_efun"):
            for condition in ("all", "abs", "csplit"):
                plot_topn_lines(
                    membership,
                    pipeline,
                    family,
                    condition,
                    args.figure_dir / dirname / f"{pipeline.lower()}_topN_density_class__{family}__{condition}.png",
                )
                plot_bar_fraction(
                    winners,
                    pipeline,
                    family,
                    condition,
                    args.figure_dir / dirname / f"{pipeline.lower()}_winner_density_class__{family}__{condition}.png",
                )

    for pipeline, rows, dirname, title_prefix in (
        ("P8", p8_mk_v2, "03_p8_blp_method_k_v2", "P8 BLP method-k v2"),
        ("P10", p10_blp_mk_v2, "04_p10_blp_method_k_v2", "P10 BLP method-k v2"),
    ):
        for family in ("efun", "deconv_efun"):
            for condition in ("abs", "csplit"):
                plot_method_k_topn_membership(
                    rows,
                    method_field="method_k",
                    family=family,
                    condition=condition,
                    title=f"{title_prefix}: topN membership | {family} | {condition}",
                    path=args.figure_dir / dirname / f"{pipeline.lower()}_blp_method_k_topN_membership__{family}__{condition}.png",
                )
                plot_method_k_grid_metric(
                    rows,
                    method_field="method_k",
                    metric_field="winner_fraction",
                    top_n=args.heatmap_top_n,
                    family=family,
                    condition=condition,
                    title=f"{title_prefix}: top1 winner fraction after top{args.heatmap_top_n} xcorr | {family} | {condition}",
                    path=args.figure_dir / dirname / f"{pipeline.lower()}_blp_method_k_winner_fraction__{family}__{condition}.png",
                )
                plot_method_k_grid_metric(
                    rows,
                    method_field="method_k",
                    metric_field="median_best_rank",
                    top_n=args.heatmap_top_n,
                    family=family,
                    condition=condition,
                    title=f"{title_prefix}: median best xcorr rank | top{args.heatmap_top_n} | {family} | {condition}",
                    path=args.figure_dir / dirname / f"{pipeline.lower()}_blp_method_k_median_rank__{family}__{condition}.png",
                    value_format=".1f",
                )

    for family in ("efun", "deconv_efun"):
        for condition in ("all", "abs", "csplit"):
            plot_method_k_method_topk_membership(
                p10_bold_mk_v2,
                method_field="p9_method_k",
                family=family,
                condition=condition,
                title=f"P10 BOLD-side P9 method-k v2.1: method-rank topK fraction | density top{args.heatmap_top_n} | {family} | {condition}",
                path=args.figure_dir / "05_p10_bold_method_k_v2_1" / f"p10_bold_method_k_method_topK_fraction__{family}__{condition}.png",
            )
            plot_method_k_grid_metric(
                p10_bold_mk_v2,
                method_field="p9_method_k",
                metric_field="winner_fraction",
                top_n=1,
                family=family,
                condition=condition,
                title=f"P10 BOLD-side P9 method-k v2.1: winner fraction across P9 methods | density top{args.heatmap_top_n} | {family} | {condition}",
                path=args.figure_dir / "05_p10_bold_method_k_v2_1" / f"p10_bold_method_k_winner_fraction__{family}__{condition}.png",
            )
            plot_method_k_grid_metric(
                p10_bold_mk_v2,
                method_field="p9_method_k",
                metric_field="median_method_rank",
                top_n=1,
                family=family,
                condition=condition,
                title=f"P10 BOLD-side P9 method-k v2.1: median rank among P9 methods | density top{args.heatmap_top_n} | {family} | {condition}",
                path=args.figure_dir / "05_p10_bold_method_k_v2_1" / f"p10_bold_method_k_median_method_rank__{family}__{condition}.png",
                value_format=".1f",
            )
            plot_method_k_grid_metric(
                p10_bold_mk_v2,
                method_field="p9_method_k",
                metric_field="median_score",
                top_n=1,
                family=family,
                condition=condition,
                title=f"P10 BOLD-side P9 method-k v2.1: median max |r| | density top{args.heatmap_top_n} | {family} | {condition}",
                path=args.figure_dir / "05_p10_bold_method_k_v2_1" / f"p10_bold_method_k_median_score__{family}__{condition}.png",
                value_format=".3f",
            )

    for pipeline, rows, dirname in (
        ("P8", p8_labels, "06_p8_dimred_selectivity_labels_v2"),
        ("P10", p10_labels, "07_p10_dimred_selectivity_labels_v2"),
    ):
        for family in ("efun", "deconv_efun"):
            for condition in ("abs", "csplit"):
                matrix = {
                    (str(r["method_k"]), str(r["label"])): float(r["fraction"])
                    for r in rows
                    if r["bold_feature_family"] == family and r["condition_scope"] == condition
                }
                plot_heatmap(
                    matrix,
                    method_k_rows(),
                    list(LABEL_ORDER),
                    f"{pipeline} dimred selectivity labels | top{args.heatmap_top_n} | {family} | {condition}",
                    args.figure_dir / dirname / f"{pipeline.lower()}_dimred_label_composition__{family}__{condition}.png",
                    value_format=".2f",
                )

    for pipeline, dirname in (("P8", "08_p8_raw_efun_timescales_v2"), ("P10", "09_p10_raw_efun_timescales_v2")):
        plot_raw_index_hist_by_family(
            raw_values_by_family(heatmap_contexts, pipeline, "raw_efun_index"),
            f"{pipeline} raw efun index distribution: efun vs deconv | top{args.heatmap_top_n}",
            args.figure_dir / dirname / f"{pipeline.lower()}_raw_efun_index_hist_by_family_top{args.heatmap_top_n}.png",
        )
        plot_timescale_ecdf(
            raw_values_by_family(heatmap_contexts, pipeline, "raw_timescale_sec"),
            f"{pipeline} raw efun timescale ECDF: efun vs deconv | top{args.heatmap_top_n}",
            args.figure_dir / dirname / f"{pipeline.lower()}_raw_efun_timescale_ecdf_by_family_top{args.heatmap_top_n}.png",
        )
        plot_timescale_topn_median_lines(
            raw_timescale_summary,
            pipeline,
            f"{pipeline} raw efun median timescale by density-hit topN",
            args.figure_dir / dirname / f"{pipeline.lower()}_raw_efun_timescale_median_by_topN.png",
        )
        for top_n in sorted(contexts_by_n):
            plot_timescale_ecdf(
                raw_values_by_family(contexts_by_n[top_n], pipeline, "raw_timescale_sec"),
                f"{pipeline} raw efun timescale ECDF: efun vs deconv | top{top_n}",
                args.figure_dir / dirname / "topN_sensitivity" / f"{pipeline.lower()}_raw_efun_timescale_ecdf_by_family_top{top_n}.png",
            )

    write_text_figure(
        args.figure_dir / "10_p8_roi_activation_consistency" / "p8_roi_activation_candidate_selection_needed.png",
        "P8 ROI / activation consistency",
        [
            "This second-stage figure should be generated after selecting candidate groups from density/label/timescale plots.",
            "Current run did not redraw all activation maps by default.",
        ],
    )
    write_text_figure(
        args.figure_dir / "11_p10_roi_activation_consistency" / "p10_roi_activation_candidate_selection_needed.png",
        "P10 ROI / activation consistency",
        [
            "This second-stage figure should be generated after selecting candidate groups from density/label/timescale plots.",
            "P10 must keep P9/BOLD method-k in the candidate definition.",
        ],
    )
    write_text_figure(
        args.figure_dir / "12_p8_vs_p10_comparison" / "p8_vs_p10_comparison_after_individual_figures.png",
        "P8 vs P10 comparison",
        [
            "P8 and P10 individual figure sets were generated first.",
            "Agreement plots should only use paired available conditions.",
            "Use after deciding which condition groups are scientifically relevant.",
        ],
    )

    summary = [
        "# P8/P10 parameter-selection figures",
        "",
        f"- P8 hit rows from peaks CSV: `{len(p8_hits)}`",
        f"- P10 hit rows from peaks CSV: `{len(p10_hits)}`",
        f"- Raw P5 mode metadata rows loaded: `{len(raw_metadata)}`",
        f"- Raw density hit rows with joined finite timescale: `{sum(1 for r in hits if r.get('density_group') == 'raw_efun_density' and r.get('raw_metadata_status') == 'present' and math.isfinite(as_float(r.get('raw_timescale_sec'))))}`",
        f"- TopN values: `{', '.join(str(n) for n in args.top_n_values)}`",
        f"- Heatmap topN: `{args.heatmap_top_n}`",
        f"- Figure root: `{args.figure_dir}`",
        "",
        "TopN values greater than 5 were rebuilt from full peaks CSV files, not from the old top5 CSV files.",
        "V2 method-k figures separate topN membership, top1 winner fraction, and median best-rank.",
        "Raw efun timescale figures now compare BOLD efun vs BOLD deconv_efun hit distributions.",
    ]
    (args.results_dir / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
