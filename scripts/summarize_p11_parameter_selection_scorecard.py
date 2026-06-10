#!/usr/bin/env python3
"""Build P11 parameter-selection scorecards and recompute requirements.

The script is deliberately numeric-first.  It reads current P8/P10 xcorr top
tables and optional P5 component labels, then writes:

- top-hit tables with interpreted density metadata;
- raw efun index distribution summaries;
- raw/dimred timescale metadata audit rows;
- a parameter-selection scorecard for condition/method/k/BOLD observable;
- lightweight PNG figures for the scorecard.

When newly requested metadata is absent from existing outputs, the script marks
the row as ``metadata_missing`` and records a recompute requirement instead of
failing.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median, pstdev
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
    density_display_name,
    density_group_name,
    density_metadata_fields,
    feature_family,
    finite,
    plot_heatmap,
    slug,
)


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
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
DEFAULT_CURRENT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")
DEFAULT_RESULTS_DIR = Path("results") / "pipeline11_parameter_selection_scorecard_current"
DEFAULT_LABEL_TABLE = (
    Path("results")
    / "pipeline5_dimred_component_process_labels_current"
    / "dimred_efun_process_labels.csv"
)
DEFAULT_TOP_N = 5
RAW_TIMESCALE_FIELDS = (
    "raw_efun_index",
    "source_mode_index",
    "lambda_discrete",
    "lambda_continuous_real",
    "lambda_continuous_imag",
    "timescale_sec",
    "frequency_hz",
    "lfp_activity_transform",
    "lfp_activity_window_policy",
)
DIMRED_TIMESCALE_FIELDS = (
    "component_timescale_median_sec",
    "component_timescale_mean_sec",
    "component_timescale_weighted_median_sec",
    "component_timescale_weighted_mean_sec",
    "component_timescale_iqr_sec",
    "component_timescale_source",
)
CURRENT_P5_CONDITIONS = ("abs_projected_vlambda", "complex_split_projected_vlambda")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS_ALL9))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_CURRENT_RUN_TAGS))
    parser.add_argument("--label-table", type=Path, default=DEFAULT_LABEL_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--top-n", type=int, default=DEFAULT_TOP_N)
    parser.add_argument("--skip-figures", action="store_true")
    return parser.parse_args()


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


def as_int(value: object, default: int = -1) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def fmt(value: float) -> str:
    return f"{value:.10g}" if math.isfinite(value) else ""


def lag_sign(value: float, threshold_sec: float = 1.0) -> str:
    if not math.isfinite(value):
        return "missing"
    if abs(value) <= threshold_sec:
        return "near_zero"
    return "positive" if value > 0 else "negative"


def majority_fraction(values: Sequence[str]) -> tuple[str, float]:
    clean = [v for v in values if v and v != "missing"]
    if not clean:
        return "missing", math.nan
    counts = Counter(clean)
    category, count = counts.most_common(1)[0]
    return category, count / len(clean)


def parse_p10_tag(name: str) -> tuple[str, str, str, int, str]:
    match = re.match(r"^(pv_.+?)__(.+?)__([a-z]+)_k(\d+)$", name)
    if not match:
        return "", "", "", -1, ""
    run_tag, feature, method, k = match.groups()
    k_int = int(k)
    return run_tag, feature, method, k_int, f"{method}_k{k_int:02d}"


def label_key(dataset: str, density_name: str, density_index: object) -> tuple[str, str, int]:
    return (dataset.lower(), density_name, as_int(density_index))


def load_label_map(path: Path) -> tuple[dict[tuple[str, str, int], dict[str, str]], str, dict[str, str]]:
    rows = read_csv(path)
    out: dict[tuple[str, str, int], dict[str, str]] = {}
    dataset_status: dict[str, str] = {}
    for row in rows:
        density = row.get("density_name_for_p8_p10") or row.get("density_name") or ""
        idx = row.get("density_index") or row.get("component_idx") or ""
        out[label_key(str(row.get("dataset", "")), density, idx)] = row
        dataset = str(row.get("dataset", "")).lower()
        status = str(row.get("label_source_status", "") or "present")
        if dataset:
            if status == "transitional_maxabs":
                dataset_status[dataset] = "transitional_maxabs"
            elif dataset not in dataset_status:
                dataset_status[dataset] = status
    return out, "present" if rows else "missing", dataset_status


def extract_raw_index(row: dict[str, str]) -> int:
    explicit = row.get("raw_efun_index")
    if explicit not in (None, ""):
        return as_int(explicit)
    return as_int(row.get("density_index"))


def metadata_presence(row: dict[str, str], fields: Sequence[str]) -> str:
    present = [field for field in fields if str(row.get(field, "")).strip()]
    return "present" if len(present) == len(fields) else ("partial" if present else "missing")


def enrich_hit(
    *,
    analysis_source: str,
    dataset: str,
    run_tag: str,
    source_level: str,
    source_csv: Path,
    row: dict[str, str],
    rank: int,
    top_n: int,
    label_map: dict[tuple[str, str, int], dict[str, str]],
    p9_feature: str = "",
    p9_method: str = "",
    p9_k: int = -1,
    p9_method_k: str = "",
) -> dict[str, object]:
    density = str(row.get("density_name") or "").strip()
    meta = density_metadata_fields(density)
    density_group = density_group_name(density)
    density_index = row.get("density_index", "")
    raw_index = extract_raw_index(row) if density_group == "raw_efun_density" else -1
    component_idx = as_int(density_index) if density_group == "dimred_efun_density" else -1
    peak_abs = as_float(row.get("peak_abs_corr"))
    lag = as_float(row.get("peak_lag_sec"))
    label = label_map.get(label_key(dataset, density, density_index), {})
    raw_metadata_status = metadata_presence(row, RAW_TIMESCALE_FIELDS) if density_group == "raw_efun_density" else "not_raw"
    dimred_metadata_status = metadata_presence(row, DIMRED_TIMESCALE_FIELDS) if density_group == "dimred_efun_density" else "not_dimred"
    out = {
        "analysis_source": analysis_source,
        "dataset": dataset,
        "run_tag": run_tag,
        "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
        "p9_feature": p9_feature,
        "p9_method": p9_method,
        "p9_k": p9_k if p9_k >= 0 else "",
        "p9_method_k": p9_method_k,
        "source_level": source_level,
        "top_rank": rank,
        "top_n": top_n,
        "density_name": density,
        "density_display": density_display_name(density),
        "density_group": density_group,
        **meta,
        "condition": meta.get("density_condition", ""),
        "method": meta.get("density_method", ""),
        "k": meta.get("density_k", ""),
        "method_k": meta.get("density_method_k", ""),
        "density_index": density_index,
        "raw_efun_index": raw_index if raw_index >= 0 else "",
        "component_idx": component_idx if component_idx >= 0 else "",
        "density_label": row.get("density_label", ""),
        "bold_feature": row.get("bold_feature", ""),
        "feature_family": feature_family(str(row.get("bold_feature", ""))),
        "bold_mode_index": row.get("bold_mode_index", ""),
        "bold_component_index": row.get("bold_component_index", ""),
        "peak_abs_corr": fmt(peak_abs),
        "peak_corr": row.get("peak_corr", ""),
        "peak_lag_sec": fmt(lag),
        "lag_sign": lag_sign(lag),
        "raw_timescale_metadata_status": raw_metadata_status,
        "dimred_timescale_metadata_status": dimred_metadata_status,
        "raw_timescale_sec": row.get("timescale_sec", ""),
        "raw_frequency_hz": row.get("frequency_hz", ""),
        "dimred_component_timescale_median_sec": row.get("component_timescale_median_sec", ""),
        "dimred_component_timescale_weighted_median_sec": row.get("component_timescale_weighted_median_sec", ""),
        "lfp_activity_transform": row.get("lfp_activity_transform", ""),
        "lfp_activity_window_policy": row.get("lfp_activity_window_policy", ""),
        "lfp_process_label": label.get("primary_process_label", "unlabeled" if density_group == "dimred_efun_density" else ""),
        "lfp_selective_family_set": label.get("selective_family_set", ""),
        "lfp_active_family_set": label.get("active_family_set", ""),
        "lfp_label_confidence": label.get("label_confidence", ""),
        "lfp_label_source_status": label.get("label_source_status", ""),
        "lfp_label_source_peak_stats_file": label.get("source_peak_stats_file", ""),
        "source_csv": str(source_csv),
    }
    return out


def collect_p8_hits(processed_root: Path, datasets: Sequence[str], run_tags: Sequence[str], top_n: int, label_map: dict[tuple[str, str, int], dict[str, str]]) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    for dataset in datasets:
        for run_tag in run_tags:
            root = processed_root / dataset / "pipeline8_xcorr" / run_tag
            for csv_path in sorted((root / "feature").glob("*/*_top__*.csv")):
                rows = read_csv(csv_path)[:top_n]
                for rank, row in enumerate(rows, start=1):
                    hits.append(
                        enrich_hit(
                            analysis_source="P8",
                            dataset=dataset,
                            run_tag=run_tag,
                            source_level="per_density_feature",
                            source_csv=csv_path,
                            row=row,
                            rank=rank,
                            top_n=top_n,
                            label_map=label_map,
                        )
                    )
            if not (root / "feature").is_dir():
                for csv_path in sorted((root / "density").glob("*_top__*.csv")):
                    rows = read_csv(csv_path)[:top_n]
                    for rank, row in enumerate(rows, start=1):
                        hits.append(
                            enrich_hit(
                                analysis_source="P8",
                                dataset=dataset,
                                run_tag=run_tag,
                                source_level="per_density",
                                source_csv=csv_path,
                                row=row,
                                rank=rank,
                                top_n=top_n,
                                label_map=label_map,
                            )
                        )
    return hits


def collect_p10_hits(processed_root: Path, datasets: Sequence[str], run_tags: Sequence[str], top_n: int, label_map: dict[tuple[str, str, int], dict[str, str]]) -> list[dict[str, object]]:
    hits: list[dict[str, object]] = []
    for dataset in datasets:
        root = processed_root / dataset / "pipeline10_dimred_xcorr"
        if not root.is_dir():
            continue
        for p10_dir in sorted(path for path in root.iterdir() if path.is_dir()):
            run_tag, p9_feature, p9_method, p9_k, p9_method_k = parse_p10_tag(p10_dir.name)
            if run_tag not in run_tags:
                continue
            for csv_path in sorted((p10_dir / "feature").glob("*/*_top__*.csv")):
                rows = read_csv(csv_path)[:top_n]
                for rank, row in enumerate(rows, start=1):
                    hits.append(
                        enrich_hit(
                            analysis_source="P10",
                            dataset=dataset,
                            run_tag=run_tag,
                            source_level="per_density_feature",
                            source_csv=csv_path,
                            row=row,
                            rank=rank,
                            top_n=top_n,
                            label_map=label_map,
                            p9_feature=p9_feature,
                            p9_method=p9_method,
                            p9_k=p9_k,
                            p9_method_k=p9_method_k,
                        )
                    )
            if not (p10_dir / "feature").is_dir():
                for csv_path in sorted((p10_dir / "density").glob("*_top__*.csv")):
                    rows = read_csv(csv_path)[:top_n]
                    for rank, row in enumerate(rows, start=1):
                        hits.append(
                            enrich_hit(
                                analysis_source="P10",
                                dataset=dataset,
                                run_tag=run_tag,
                                source_level="per_density",
                                source_csv=csv_path,
                                row=row,
                                rank=rank,
                                top_n=top_n,
                                label_map=label_map,
                                p9_feature=p9_feature,
                                p9_method=p9_method,
                                p9_k=p9_k,
                                p9_method_k=p9_method_k,
                            )
                        )
    return hits


def weighted_median(pairs: Sequence[tuple[float, float]]) -> float:
    clean = sorted((x, max(w, 0.0)) for x, w in pairs if math.isfinite(x) and math.isfinite(w) and w >= 0)
    if not clean:
        return math.nan
    total = sum(w for _x, w in clean)
    if total <= 0:
        return median(x for x, _w in clean)
    target = total / 2
    running = 0.0
    for x, w in clean:
        running += w
        if running >= target:
            return x
    return clean[-1][0]


def summarize_raw_index_distribution(hits: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in hits:
        if row.get("density_group") != "raw_efun_density":
            continue
        key = (
            str(row.get("analysis_source")),
            str(row.get("dataset")),
            str(row.get("observable")),
            str(row.get("density_name")),
            str(row.get("feature_family")),
            str(row.get("p9_method_k")),
        )
        grouped[key].append(row)
    out: list[dict[str, object]] = []
    for key, rows in sorted(grouped.items()):
        source, dataset, obs, density, family, p9_method_k = key
        pairs: list[tuple[float, float]] = []
        counts: Counter[int] = Counter()
        q_counts: Counter[str] = Counter()
        for row in rows:
            idx = as_int(row.get("raw_efun_index"))
            if idx < 0:
                continue
            weight = as_float(row.get("peak_abs_corr"))
            pairs.append((float(idx), weight))
            counts[idx] += 1
            if idx <= 50:
                q_counts["idx_001_050"] += 1
            elif idx <= 100:
                q_counts["idx_051_100"] += 1
            elif idx <= 200:
                q_counts["idx_101_200"] += 1
            else:
                q_counts["idx_201_plus"] += 1
        values = [idx for idx, _w in pairs]
        weights = [w for _idx, w in pairs]
        out.append(
            {
                "analysis_source": source,
                "dataset": dataset,
                "observable": obs,
                "p9_method_k": p9_method_k,
                "density_name": density,
                "feature_family": family,
                "n_hits": len(rows),
                "n_unique_raw_efun_index": len(counts),
                "raw_efun_index_set": ";".join(str(i) for i in sorted(counts)),
                "raw_efun_index_histogram": ";".join(f"{idx}:{count}" for idx, count in sorted(counts.items())),
                "raw_efun_index_bin_histogram": ";".join(f"{name}:{q_counts[name]}" for name in ("idx_001_050", "idx_051_100", "idx_101_200", "idx_201_plus")),
                "weighted_median_raw_efun_index": fmt(weighted_median(pairs)),
                "mean_raw_efun_index": fmt(mean(values)) if values else "",
                "weighted_mean_raw_efun_index": fmt(sum(x * w for x, w in pairs) / sum(weights)) if weights and sum(weights) > 0 else "",
                "raw_timescale_metadata_status": "present" if all(str(row.get("raw_timescale_metadata_status")) == "present" for row in rows) else "missing_or_partial",
            }
        )
    return out


def summarize_top_hit_interpretation(hits: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str, str, str, str, str, str, int, str], list[dict[str, object]]] = defaultdict(list)
    for row in hits:
        density_group = str(row.get("density_group", ""))
        if density_group not in {"raw_efun_density", "dimred_efun_density"}:
            continue
        k = as_int(row.get("k"), -1)
        key = (
            str(row.get("analysis_source", "")),
            str(row.get("dataset", "")),
            str(row.get("observable", "")),
            str(row.get("feature_family", "")),
            str(row.get("p9_feature", "")),
            str(row.get("p9_method_k", "")),
            density_group,
            str(row.get("condition", "")),
            str(row.get("method", "")),
            k,
            str(row.get("density_name", "")),
        )
        grouped[key].append(row)

    out: list[dict[str, object]] = []
    for key, rows in sorted(grouped.items()):
        (
            source,
            dataset,
            observable,
            feat_family,
            p9_feature,
            p9_method_k,
            density_group,
            condition,
            method,
            k,
            density_name,
        ) = key
        values = finite(as_float(row.get("peak_abs_corr")) for row in rows)
        lag_values = finite(as_float(row.get("peak_lag_sec")) for row in rows)
        raw_counts: Counter[int] = Counter()
        comp_counts: Counter[int] = Counter()
        label_counts: Counter[str] = Counter()
        raw_timescale_pairs: list[tuple[float, float]] = []
        dim_timescale_pairs: list[tuple[float, float]] = []
        for row in rows:
            weight = as_float(row.get("peak_abs_corr"))
            if density_group == "raw_efun_density":
                raw_idx = as_int(row.get("raw_efun_index"), -1)
                if raw_idx >= 0:
                    raw_counts[raw_idx] += 1
                ts = as_float(row.get("raw_timescale_sec"))
                if math.isfinite(ts) and math.isfinite(weight):
                    raw_timescale_pairs.append((ts, weight))
            elif density_group == "dimred_efun_density":
                comp_idx = as_int(row.get("component_idx"), -1)
                if comp_idx >= 0:
                    comp_counts[comp_idx] += 1
                label_counts[str(row.get("lfp_process_label") or "unlabeled")] += 1
                ts = as_float(
                    row.get("dimred_component_timescale_weighted_median_sec")
                    or row.get("dimred_component_timescale_median_sec")
                )
                if math.isfinite(ts) and math.isfinite(weight):
                    dim_timescale_pairs.append((ts, weight))
        out.append(
            {
                "analysis_source": source,
                "dataset": dataset,
                "observable": observable,
                "bold_feature_family": feat_family,
                "p9_feature": p9_feature,
                "p9_method_k": p9_method_k,
                "density_group": density_group,
                "condition": condition,
                "method": method,
                "k": k if k >= 0 else "",
                "density_name": density_name,
                "n_top_hits": len(rows),
                "mean_peak_abs_corr": fmt(mean(values)) if values else "",
                "max_peak_abs_corr": fmt(max(values)) if values else "",
                "median_lag_sec": fmt(median(lag_values)) if lag_values else "",
                "top_raw_efun_index_histogram": ";".join(f"{idx}:{count}" for idx, count in sorted(raw_counts.items())),
                "top_raw_efun_index_set": ";".join(str(idx) for idx in sorted(raw_counts)),
                "weighted_median_raw_timescale_sec": fmt(weighted_median(raw_timescale_pairs)),
                "top_dimred_component_histogram": ";".join(f"{idx}:{count}" for idx, count in sorted(comp_counts.items())),
                "top_dimred_component_set": ";".join(str(idx) for idx in sorted(comp_counts)),
                "top_dimred_process_label_histogram": ";".join(f"{label}:{count}" for label, count in label_counts.most_common()),
                "weighted_median_dimred_timescale_sec": fmt(weighted_median(dim_timescale_pairs)),
                "metadata_status": (
                    "raw_timescale:"
                    + ("present" if raw_timescale_pairs else "missing_or_partial")
                    if density_group == "raw_efun_density"
                    else "dimred_label:"
                    + ("present" if label_counts else "missing")
                    + ";dimred_timescale:"
                    + ("present" if dim_timescale_pairs else "missing_or_partial")
                ),
            }
        )
    return out


def label_score(label: str) -> float:
    if label in {"theta_selective_similar", "ripple_selective_similar", "theta_ripple_joint"}:
        return 1.0
    if label in {"theta_selective_unequal", "ripple_selective_unequal", "mixed_theta_ripple"}:
        return 0.8
    if label in {"mixed_theta_gamma", "mixed_ripple_gamma"}:
        return 0.5
    if label in {"pan_event", "gamma_selective"}:
        return 0.25
    return 0.0


def complexity_score(method: str, k: int, condition: str) -> float:
    method_penalty = {"svd": 0.05, "nmf": 0.10, "mds": 0.20, "umap": 0.25}.get(method, 0.25)
    condition_penalty = 0.0 if condition == "abs" else 0.10
    k_penalty = max(k - 3, 0) * 0.04
    return max(0.0, 1.0 - method_penalty - condition_penalty - k_penalty)


def build_scorecard(hits: Sequence[dict[str, object]], max_dataset_count: int) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str, str, int, str], list[dict[str, object]]] = defaultdict(list)
    for row in hits:
        if row.get("density_group") != "dimred_efun_density":
            continue
        condition = str(row.get("condition", ""))
        method = str(row.get("method", ""))
        k = as_int(row.get("k"))
        if condition not in {"abs", "csplit"} or method not in METHOD_ORDER or k not in COMPONENT_COUNTS:
            continue
        key = (
            str(row.get("analysis_source")),
            condition,
            method,
            f"k{k:02d}",
            str(row.get("observable")),
            k,
            str(row.get("feature_family")),
        )
        grouped[key].append(row)

    raw_rows: list[dict[str, object]] = []
    for key, rows in sorted(grouped.items()):
        source, condition, method, k_tag, observable, k, family = key
        datasets = sorted({str(r.get("dataset")) for r in rows})
        values = finite(as_float(r.get("peak_abs_corr")) for r in rows)
        lags = [lag_sign(as_float(r.get("peak_lag_sec"))) for r in rows]
        _lag_cat, lag_agree = majority_fraction(lags)
        labels = [str(r.get("lfp_process_label") or "unlabeled") for r in rows]
        label_values = [label_score(label) for label in labels]
        label_counts = Counter(labels)
        coverage = len(datasets) / max(max_dataset_count, 1)
        strength = mean(values) if values else math.nan
        label_avg = mean(label_values) if label_values else 0.0
        lag_score = lag_agree if math.isfinite(lag_agree) else 0.0
        comp_score = complexity_score(method, int(k), condition)
        has_time = any(str(r.get("dimred_timescale_metadata_status")) == "present" for r in rows)
        timescale_score = math.nan if not has_time else 1.0
        final = (
            0.25 * coverage
            + 0.30 * (strength if math.isfinite(strength) else 0.0)
            + 0.20 * label_avg
            + 0.15 * lag_score
            + 0.10 * comp_score
        )
        tier = "needs_metadata_recompute" if not has_time else ("candidate" if coverage >= 0.5 else "low_coverage")
        raw_rows.append(
            {
                "analysis_source": source,
                "condition": condition,
                "method": method,
                "k": int(k),
                "method_k": f"{method}_k{int(k):02d}",
                "bold_observable": observable,
                "feature_family": family,
                "n_datasets": len(datasets),
                "datasets": ";".join(datasets),
                "coverage_score": fmt(coverage),
                "mean_topN_abs_corr": fmt(strength),
                "median_topN_abs_corr": fmt(median(values)) if values else "",
                "std_topN_abs_corr": fmt(pstdev(values)) if len(values) > 1 else "0",
                "label_agreement": ";".join(f"{label}:{count}" for label, count in label_counts.most_common()),
                "theta_ripple_label_score": fmt(label_avg),
                "preference_agreement": "",
                "lag_sign_agreement": fmt(lag_score),
                "roi_profile_consistency": "",
                "raw_weighted_median_timescale_sec": "",
                "dimred_weighted_median_timescale_sec": "",
                "raw_vs_dimred_timescale_delta_sec": "",
                "timescale_agreement": "" if not has_time else "1",
                "complexity_score": fmt(comp_score),
                "final_score": fmt(final),
                "recommendation_tier": tier,
                "notes": "score excludes timescale agreement until P5 density metadata is regenerated" if not has_time else "",
            }
        )
    return sorted(raw_rows, key=lambda r: -as_float(r.get("final_score")))


def dataset_recompute_requirements(
    processed_root: Path,
    datasets: Sequence[str],
    hits: Sequence[dict[str, object]],
    labels_status: str,
    label_dataset_status: dict[str, str],
) -> list[dict[str, object]]:
    by_dataset = defaultdict(list)
    for row in hits:
        by_dataset[str(row.get("dataset"))].append(row)
    out: list[dict[str, object]] = []
    for dataset in datasets:
        ds_root = processed_root / dataset
        rows = by_dataset.get(dataset, [])
        raw_rows = [r for r in rows if r.get("density_group") == "raw_efun_density"]
        dim_rows = [r for r in rows if r.get("density_group") == "dimred_efun_density"]
        p8_root = ds_root / "pipeline8_xcorr"
        p10_root = ds_root / "pipeline10_dimred_xcorr"
        p5_raw = ds_root / "pipeline5_raw_thresholded_density"
        p5_dim = ds_root / "pipeline5_dimred_thresholded_density"
        dataset_label_status = label_dataset_status.get(dataset.lower(), "missing")
        label_needed = (
            labels_status != "present"
            or dataset_label_status in {"missing", "transitional_maxabs"}
            or any(str(r.get("lfp_label_source_status")) == "transitional_maxabs" for r in dim_rows)
        )
        raw_meta_missing = bool(raw_rows) and any(str(r.get("raw_timescale_metadata_status")) != "present" for r in raw_rows)
        dim_meta_missing = bool(dim_rows) and any(str(r.get("dimred_timescale_metadata_status")) != "present" for r in dim_rows)
        actions: list[str] = []
        reasons: list[str] = []
        if not ds_root.is_dir():
            actions.append("create_processed_dataset_outputs")
            reasons.append("processed dataset folder missing")
        if not p5_raw.is_dir() or not p5_dim.is_dir():
            actions.append("run_p5_density_grid")
            reasons.append("P5 raw/dimred density folders missing")
        if label_needed:
            actions.append("rerun_p5_activity_envelope_peak_stats_and_labels")
            reasons.append("P5 labels are missing or maxabs-transitional; need abs/envelope activity labels")
        if raw_meta_missing:
            actions.append("rerun_p5_raw_density_with_mode_timescale_metadata")
            reasons.append("raw density xcorr rows lack explicit raw_efun_index/timescale/activity metadata")
        if dim_meta_missing:
            actions.append("rerun_p5_dimred_density_with_component_timescale_metadata")
            reasons.append("dimred density rows lack component timescale metadata")
        if not p8_root.is_dir():
            actions.append("run_p8_current_xcorr")
            reasons.append("P8 xcorr folder missing")
        elif raw_meta_missing or dim_meta_missing or label_needed:
            actions.append("rerun_p8_after_p5_density_metadata_refresh")
            reasons.append("P8 top rows need refreshed P5 metadata/labels")
        if not p10_root.is_dir():
            actions.append("run_p10_current_xcorr")
            reasons.append("P10 xcorr folder missing")
        elif raw_meta_missing or dim_meta_missing or label_needed:
            actions.append("rerun_p10_after_p5_density_metadata_refresh")
            reasons.append("P10 top rows need refreshed P5 metadata/labels")
        out.append(
            {
                "dataset": dataset,
                "processed_folder_present": int(ds_root.is_dir()),
                "p5_raw_density_present": int(p5_raw.is_dir()),
                "p5_dimred_density_present": int(p5_dim.is_dir()),
                "p8_present": int(p8_root.is_dir()),
                "p10_present": int(p10_root.is_dir()),
                "observed_top_hit_rows": len(rows),
                "observed_raw_top_hit_rows": len(raw_rows),
                "observed_dimred_top_hit_rows": len(dim_rows),
                "p5_process_label_status": dataset_label_status,
                "needs_recompute": int(bool(actions)),
                "recommended_actions": ";".join(dict.fromkeys(actions)),
                "reasons": " | ".join(dict.fromkeys(reasons)),
            }
        )
    return out


def load_font(size: int, bold: bool = False):
    if ImageFont is None:
        return None
    candidates = [
        "arialbd.ttf" if bold else "arial.ttf",
        "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except Exception:
            pass
    return ImageFont.load_default()


def write_text_figure(path: Path, title: str, lines: Sequence[str]) -> None:
    if Image is None:
        return
    width = 1200
    height = max(360, 110 + 30 * len(lines))
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(24, True)
    body_font = load_font(16, False)
    draw.text((24, 22), title, fill=(20, 30, 40), font=title_font)
    y = 76
    for line in lines:
        draw.text((30, y), line, fill=(60, 60, 60), font=body_font)
        y += 30
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def plot_raw_index_histogram(summary_rows: Sequence[dict[str, object]], path: Path) -> None:
    if Image is None:
        return
    counts: Counter[int] = Counter()
    for row in summary_rows:
        hist = str(row.get("raw_efun_index_histogram", ""))
        for item in hist.split(";"):
            if ":" not in item:
                continue
            key, value = item.split(":", 1)
            idx = as_int(key)
            count = as_int(value, 0)
            if idx >= 0:
                counts[idx] += count
    if not counts:
        write_text_figure(path, "Top xcorr raw efun index distribution", ["No raw efun top-hit rows were found."])
        return
    width, height = 1200, 520
    margin_l, margin_r, margin_t, margin_b = 70, 30, 65, 70
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(22, True)
    label_font = load_font(13, False)
    draw.text((22, 20), "Top xcorr raw efun index distribution", fill=(20, 30, 40), font=title_font)
    max_idx = max(counts)
    max_count = max(counts.values())
    plot_w = width - margin_l - margin_r
    plot_h = height - margin_t - margin_b
    bar_w = max(1, plot_w / max(max_idx, 1))
    draw.line((margin_l, margin_t + plot_h, margin_l + plot_w, margin_t + plot_h), fill=(80, 80, 80))
    draw.line((margin_l, margin_t, margin_l, margin_t + plot_h), fill=(80, 80, 80))
    for idx, count in counts.items():
        x0 = margin_l + (idx - 1) * bar_w
        x1 = x0 + max(1, bar_w * 0.9)
        y1 = margin_t + plot_h
        y0 = y1 - (count / max_count) * plot_h
        draw.rectangle((x0, y0, x1, y1), fill=(52, 126, 184))
    draw.text((margin_l, height - 45), "raw efun index", fill=(50, 50, 50), font=label_font)
    draw.text((12, margin_t + 5), "count", fill=(50, 50, 50), font=label_font)
    draw.text((margin_l, height - 25), f"n indices={len(counts)}; max index={max_idx}; max count={max_count}", fill=(90, 90, 90), font=label_font)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def scorecard_matrix_rows(scorecard: Sequence[dict[str, object]], limit: int = 30) -> tuple[dict[tuple[str, str], float], list[str], list[str]]:
    cols = [
        "coverage_score",
        "mean_topN_abs_corr",
        "theta_ripple_label_score",
        "lag_sign_agreement",
        "complexity_score",
        "final_score",
    ]
    rows = []
    matrix: dict[tuple[str, str], float] = {}
    for idx, row in enumerate(scorecard[:limit], start=1):
        label = (
            f"{idx:02d} {row['analysis_source']} {row['condition']} "
            f"{row['method']} k{int(row['k']):02d} {row['bold_observable']} {row['feature_family']}"
        )
        rows.append(label)
        for col in cols:
            matrix[(label, col)] = as_float(row.get(col))
    return matrix, rows, cols


def plot_scorecard_figures(scorecard: Sequence[dict[str, object]], raw_summary: Sequence[dict[str, object]], figure_dir: Path) -> list[Path]:
    paths: list[Path] = []
    figure_dir.mkdir(parents=True, exist_ok=True)
    if not scorecard:
        path = figure_dir / "parameter_scorecard_metric_matrix.png"
        write_text_figure(path, "Parameter scorecard", ["No dimred P8/P10 hits were found."])
        return [path]

    matrix, rows, cols = scorecard_matrix_rows(scorecard)
    path = figure_dir / "parameter_scorecard_metric_matrix.png"
    plot_heatmap(matrix, rows, cols, "P11 parameter scorecard metric matrix", path, value_format=".2f")
    paths.append(path)

    path = figure_dir / "candidate_decision_panel.png"
    plot_heatmap(matrix, rows[:20], cols, "P11 candidate decision panel | top 20", path, value_format=".2f")
    paths.append(path)

    cond_group: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in scorecard:
        for metric in ("final_score", "mean_topN_abs_corr", "theta_ripple_label_score", "lag_sign_agreement"):
            value = as_float(row.get(metric))
            if math.isfinite(value):
                cond_group[(str(row["condition"]), metric)].append(value)
    cond_matrix = {key: mean(vals) for key, vals in cond_group.items() if vals}
    path = figure_dir / "condition_comparison_abs_vs_csplit.png"
    plot_heatmap(cond_matrix, ["abs", "csplit"], ["final_score", "mean_topN_abs_corr", "theta_ripple_label_score", "lag_sign_agreement"], "Condition comparison: abs vs csplit", path, value_format=".2f")
    paths.append(path)

    for condition in ("abs", "csplit"):
        subset = [r for r in scorecard if r.get("condition") == condition]
        matrix = {
            (str(r["bold_observable"]), f"{r['method']}_k{int(r['k']):02d}"): as_float(r.get("final_score"))
            for r in subset
        }
        rows_obs = sorted({str(r["bold_observable"]) for r in subset}) or list(RUN_TAG_LABELS.values())
        cols_mk = [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]
        path = figure_dir / f"method_k_selection_heatmap__{condition}.png"
        plot_heatmap(matrix, rows_obs, cols_mk, f"Method x k selection score | {condition}", path, value_format=".2f")
        paths.append(path)

    obs_group: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in scorecard:
        for metric in ("final_score", "mean_topN_abs_corr", "theta_ripple_label_score", "lag_sign_agreement", "complexity_score"):
            value = as_float(row.get(metric))
            if math.isfinite(value):
                obs_group[(str(row["bold_observable"]), metric)].append(value)
    obs_matrix = {key: mean(vals) for key, vals in obs_group.items() if vals}
    path = figure_dir / "bold_observable_metric_profile.png"
    plot_heatmap(obs_matrix, list(RUN_TAG_LABELS.values()), ["final_score", "mean_topN_abs_corr", "theta_ripple_label_score", "lag_sign_agreement", "complexity_score"], "BOLD observable metric profile", path, value_format=".2f")
    paths.append(path)

    k_group: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in scorecard:
        label = f"{row['condition']} {row['method']}"
        value = as_float(row.get("final_score"))
        if math.isfinite(value):
            k_group[(label, f"k{int(row['k']):02d}")].append(value)
    k_matrix = {key: mean(vals) for key, vals in k_group.items() if vals}
    path = figure_dir / "k_sensitivity_curves.png"
    plot_heatmap(k_matrix, sorted({key[0] for key in k_matrix}), [f"k{k:02d}" for k in COMPONENT_COUNTS], "k sensitivity: mean final score", path, value_format=".2f")
    paths.append(path)

    label_group: dict[tuple[str, str], int] = defaultdict(int)
    totals: Counter[str] = Counter()
    for row in scorecard:
        label = f"{row['condition']} {row['method']}_k{int(row['k']):02d}"
        for item in str(row.get("label_agreement", "")).split(";"):
            if ":" not in item:
                continue
            cat, count = item.rsplit(":", 1)
            n = as_int(count, 0)
            label_group[(label, cat)] += n
            totals[label] += n
    label_cols = ["theta_selective_similar", "theta_selective_unequal", "ripple_selective_similar", "ripple_selective_unequal", "theta_ripple_joint", "mixed_theta_ripple", "pan_event", "unlabeled"]
    label_matrix = {
        key: (count / totals[key[0]] if totals[key[0]] else math.nan)
        for key, count in label_group.items()
    }
    path = figure_dir / "process_label_composition.png"
    plot_heatmap(label_matrix, sorted(totals), label_cols, "Process-label composition among scorecard hits", path, value_format=".2f")
    paths.append(path)

    path = figure_dir / "raw_vs_dimred_timescale_agreement.png"
    write_text_figure(
        path,
        "Raw-vs-dimred timescale agreement",
        [
            "Current outputs do not yet contain raw/dimred timescale metadata.",
            "This figure becomes numeric after P5 raw/dimred density is regenerated with timescale provenance.",
        ],
    )
    paths.append(path)

    path = figure_dir / "top_xcorr_raw_efun_index_distribution.png"
    plot_raw_index_histogram(raw_summary, path)
    paths.append(path)
    return paths


def write_summary(path: Path, scorecard: Sequence[dict[str, object]], recompute: Sequence[dict[str, object]], hits: Sequence[dict[str, object]]) -> None:
    lines = [
        "# P11 Parameter Selection Scorecard",
        "",
        f"- Top-hit rows scanned: `{len(hits)}`",
        f"- Scorecard rows: `{len(scorecard)}`",
        f"- Datasets needing at least one recompute action: `{sum(1 for r in recompute if int(r.get('needs_recompute', 0)))}` / `{len(recompute)}`",
        "",
        "## Top Candidates",
        "",
    ]
    for row in scorecard[:20]:
        lines.append(
            "- "
            f"{row['analysis_source']} | {row['condition']} | {row['method']}_k{int(row['k']):02d} | "
            f"{row['bold_observable']} | {row['feature_family']} | "
            f"score={row['final_score']} | n={row['n_datasets']} | label={row['label_agreement']}"
        )
    lines.extend(["", "## Recompute Requirements", ""])
    for row in recompute:
        if int(row.get("needs_recompute", 0)):
            lines.append(f"- `{row['dataset']}`: {row['recommended_actions']} ({row['reasons']})")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    label_map, label_status, label_dataset_status = load_label_map(args.label_table)
    p8_hits = collect_p8_hits(args.processed_root, args.datasets, args.run_tags, args.top_n, label_map)
    p10_hits = collect_p10_hits(args.processed_root, args.datasets, args.run_tags, args.top_n, label_map)
    hits = p8_hits + p10_hits
    raw_summary = summarize_raw_index_distribution(hits)
    interpretation_summary = summarize_top_hit_interpretation(hits)
    scorecard = build_scorecard(hits, max_dataset_count=len(args.datasets))
    recompute = dataset_recompute_requirements(args.processed_root, args.datasets, hits, label_status, label_dataset_status)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "p11_parameter_selection_top_hits.csv", hits)
    write_csv(args.output_dir / "raw_efun_index_distribution.csv", raw_summary)
    write_csv(args.output_dir / "top_hit_interpretation_summary.csv", interpretation_summary)
    write_csv(args.output_dir / "parameter_selection_scorecard.csv", scorecard)
    write_csv(args.output_dir / "recompute_requirements_by_dataset.csv", recompute)

    metadata_rows = []
    for row in hits:
        if row.get("density_group") in {"raw_efun_density", "dimred_efun_density"}:
            metadata_rows.append(
                {
                    "analysis_source": row.get("analysis_source"),
                    "dataset": row.get("dataset"),
                    "observable": row.get("observable"),
                    "bold_feature": row.get("bold_feature"),
                    "feature_family": row.get("feature_family"),
                    "source_level": row.get("source_level"),
                    "top_rank": row.get("top_rank"),
                    "top_n": row.get("top_n"),
                    "p9_feature": row.get("p9_feature"),
                    "p9_method_k": row.get("p9_method_k"),
                    "density_name": row.get("density_name"),
                    "density_group": row.get("density_group"),
                    "condition": row.get("condition"),
                    "method": row.get("method"),
                    "k": row.get("k"),
                    "density_index": row.get("density_index"),
                    "raw_efun_index": row.get("raw_efun_index"),
                    "component_idx": row.get("component_idx"),
                    "peak_abs_corr": row.get("peak_abs_corr"),
                    "peak_lag_sec": row.get("peak_lag_sec"),
                    "raw_timescale_metadata_status": row.get("raw_timescale_metadata_status"),
                    "dimred_timescale_metadata_status": row.get("dimred_timescale_metadata_status"),
                    "raw_timescale_sec": row.get("raw_timescale_sec"),
                    "raw_frequency_hz": row.get("raw_frequency_hz"),
                    "dimred_component_timescale_median_sec": row.get("dimred_component_timescale_median_sec"),
                    "dimred_component_timescale_weighted_median_sec": row.get("dimred_component_timescale_weighted_median_sec"),
                    "lfp_process_label": row.get("lfp_process_label"),
                    "lfp_selective_family_set": row.get("lfp_selective_family_set"),
                    "lfp_active_family_set": row.get("lfp_active_family_set"),
                    "lfp_label_source_status": row.get("lfp_label_source_status"),
                    "source_csv": row.get("source_csv"),
                }
            )
    write_csv(args.output_dir / "timescale_and_label_metadata_audit.csv", metadata_rows)

    figures: list[Path] = []
    if not args.skip_figures:
        figures = plot_scorecard_figures(scorecard, raw_summary, args.output_dir / "figures")
    write_summary(args.output_dir / "summary.md", scorecard, recompute, hits)

    print(f"Datasets      : {', '.join(args.datasets)}")
    print(f"P8 hits       : {len(p8_hits)}")
    print(f"P10 hits      : {len(p10_hits)}")
    print(f"Scorecard rows: {len(scorecard)}")
    print(f"Raw summaries : {len(raw_summary)}")
    print(f"Interpretation: {len(interpretation_summary)}")
    print(f"Recompute rows: {len(recompute)}")
    print(f"Figures       : {len(figures)}")
    print(f"Output dir    : {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
