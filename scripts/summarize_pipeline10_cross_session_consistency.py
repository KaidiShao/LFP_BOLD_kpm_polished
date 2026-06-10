"""Summarize current P10 cross-session density/component consistency."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, pstdev
from typing import Iterable, Sequence

from summarize_pipeline8_cross_session_consistency import (
    COMPONENT_COUNTS,
    DEFAULT_PROCESSED_ROOT,
    METHOD_ORDER,
    RUN_TAG_LABELS,
    as_float,
    density_display_name,
    density_group_name,
    density_metadata,
    density_metadata_fields,
    feature_family,
    finite,
    plot_heatmap,
    slug,
    write_csv,
)


DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "f12m01")
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")
DEFAULT_FEATURES = ("efun_real",)
DEFAULT_RESULTS_DIR = Path("results") / "pipeline10_cross_session_consistency_current"
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p10"
    / "cross_session_consistency"
)
DEFAULT_RESULTS_FIGURE_DIR = DEFAULT_RESULTS_DIR / "figures"


@dataclass(frozen=True)
class P10ScoreRow:
    dataset: str
    run_tag: str
    observable: str
    p9_feature: str
    p9_method: str
    p9_k: int
    p9_method_k: str
    density_name: str
    bold_feature: str
    source_level: str
    n_top_rows: int
    mean_peak_abs_corr: float
    max_peak_abs_corr: float
    mean_peak_corr: float
    mean_peak_lag_sec: float
    source_csv: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--features", nargs="+", default=list(DEFAULT_FEATURES))
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--xcorr-save-tag", default="dimred_xcorr")
    parser.add_argument("--results-figure-dir", type=Path, default=DEFAULT_RESULTS_FIGURE_DIR)
    parser.add_argument("--skip-figures", action="store_true")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def parse_p10_tag(name: str) -> tuple[str, str, str, int, str]:
    m = re.match(r"^(pv_[^_]+(?:_[^_]+)*)__([^_]+(?:_[^_]+)*)__([a-z]+)_k(\d+)$", name)
    if not m:
        return "", "", "", -1, ""
    run_tag, feature, method, k = m.groups()
    k_int = int(k)
    return run_tag, feature, method, k_int, f"{method}_k{k_int:02d}"


def summarize_rows(
    dataset: str,
    run_tag: str,
    p9_feature: str,
    p9_method: str,
    p9_k: int,
    p9_method_k: str,
    source_level: str,
    csv_path: Path,
    rows: Sequence[dict[str, str]],
) -> list[P10ScoreRow]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        density = (row.get("density_name") or "unknown_density").strip()
        feature = (row.get("bold_feature") or "all_components").strip()
        grouped[(density, feature)].append(row)

    out: list[P10ScoreRow] = []
    for (density, feature), group_rows in sorted(grouped.items()):
        abs_vals = finite(as_float(r.get("peak_abs_corr")) for r in group_rows)
        corr_vals = finite(as_float(r.get("peak_corr")) for r in group_rows)
        lag_vals = finite(as_float(r.get("peak_lag_sec")) for r in group_rows)
        if not abs_vals:
            continue
        out.append(
            P10ScoreRow(
                dataset=dataset,
                run_tag=run_tag,
                observable=RUN_TAG_LABELS.get(run_tag, run_tag),
                p9_feature=p9_feature,
                p9_method=p9_method,
                p9_k=p9_k,
                p9_method_k=p9_method_k,
                density_name=density,
                bold_feature=feature,
                source_level=source_level,
                n_top_rows=len(group_rows),
                mean_peak_abs_corr=mean(abs_vals),
                max_peak_abs_corr=max(abs_vals),
                mean_peak_corr=mean(corr_vals) if corr_vals else math.nan,
                mean_peak_lag_sec=mean(lag_vals) if lag_vals else math.nan,
                source_csv=str(csv_path),
            )
        )
    return out


def make_hit_row(
    dataset: str,
    run_tag: str,
    p9_feature: str,
    p9_method: str,
    p9_k: int,
    p9_method_k: str,
    source_level: str,
    csv_path: Path,
    row: dict[str, str],
    rank: int,
    top_n: int,
) -> dict[str, object]:
    density = (row.get("density_name") or "").strip()
    bold_feature = (row.get("bold_feature") or "").strip()
    peak_abs = as_float(row.get("peak_abs_corr"))
    peak = as_float(row.get("peak_corr"))
    lag = as_float(row.get("peak_lag_sec"))
    out = {
        "dataset": dataset,
        "run_tag": run_tag,
        "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
        "p9_feature": p9_feature,
        "p9_method": p9_method,
        "p9_k": p9_k,
        "p9_method_k": p9_method_k,
        "source_level": source_level,
        "top_rank": rank,
        "top_n": top_n,
        "density_name": density,
        "density_display": density_display_name(density),
        "density_group": density_group_name(density),
        **density_metadata_fields(density),
        "density_index": row.get("density_index", ""),
        "density_label": row.get("density_label", ""),
        "bold_feature": bold_feature,
        "feature_family": feature_family(bold_feature),
        "bold_component_index": row.get("bold_component_index") or row.get("bold_mode_index", ""),
        "peak_abs_corr": f"{peak_abs:.10g}" if math.isfinite(peak_abs) else "",
        "peak_corr": f"{peak:.10g}" if math.isfinite(peak) else "",
        "peak_lag_sec": f"{lag:.10g}" if math.isfinite(lag) else "",
        "source_csv": str(csv_path),
    }
    return out


def collect_scores(
    processed_root: Path,
    datasets: Sequence[str],
    run_tags: Sequence[str],
    features: Sequence[str],
    xcorr_save_tag: str,
) -> tuple[list[P10ScoreRow], list[dict[str, object]]]:
    scores: list[P10ScoreRow] = []
    hits: list[dict[str, object]] = []
    method_tags = [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]
    for dataset in datasets:
        root = processed_root / dataset / "pipeline10_dimred_xcorr"
        for run_tag in run_tags:
            for p9_feature in features:
                for p9_method_k in method_tags:
                    p10_tag = f"{run_tag}__{p9_feature}__{p9_method_k}"
                    p10_dir = root / p10_tag
                    parsed_run, parsed_feature, p9_method, p9_k, parsed_method_k = parse_p10_tag(p10_tag)
                    if parsed_run != run_tag or parsed_feature != p9_feature:
                        continue
                    csv_path = p10_dir / f"{xcorr_save_tag}_top.csv"
                    rows = read_csv_rows(csv_path)
                    for i, row in enumerate(rows, start=1):
                        hits.append(make_hit_row(dataset, run_tag, p9_feature, p9_method, p9_k, parsed_method_k, "combined", csv_path, row, i, len(rows)))
                    scores.extend(summarize_rows(dataset, run_tag, p9_feature, p9_method, p9_k, parsed_method_k, "combined", csv_path, rows))

                    for csv_path in sorted((p10_dir / "density").glob(f"{xcorr_save_tag}_top__*.csv")):
                        rows = read_csv_rows(csv_path)
                        for i, row in enumerate(rows, start=1):
                            hits.append(make_hit_row(dataset, run_tag, p9_feature, p9_method, p9_k, parsed_method_k, "per_density", csv_path, row, i, len(rows)))
                        scores.extend(summarize_rows(dataset, run_tag, p9_feature, p9_method, p9_k, parsed_method_k, "per_density", csv_path, rows))

                    for csv_path in sorted((p10_dir / "feature").glob(f"*/{xcorr_save_tag}_top__*.csv")):
                        rows = read_csv_rows(csv_path)
                        for i, row in enumerate(rows, start=1):
                            hits.append(make_hit_row(dataset, run_tag, p9_feature, p9_method, p9_k, parsed_method_k, "per_density_feature", csv_path, row, i, len(rows)))
                        scores.extend(summarize_rows(dataset, run_tag, p9_feature, p9_method, p9_k, parsed_method_k, "per_density_feature", csv_path, rows))
    return scores, hits


def row_to_dict(row: P10ScoreRow) -> dict[str, object]:
    return {
        "dataset": row.dataset,
        "run_tag": row.run_tag,
        "observable": row.observable,
        "p9_feature": row.p9_feature,
        "p9_method": row.p9_method,
        "p9_k": row.p9_k,
        "p9_method_k": row.p9_method_k,
        "density_name": row.density_name,
        "density_display": density_display_name(row.density_name),
        "density_group": density_group_name(row.density_name),
        **density_metadata_fields(row.density_name),
        "bold_feature": row.bold_feature,
        "feature_family": feature_family(row.bold_feature),
        "source_level": row.source_level,
        "n_top_rows": row.n_top_rows,
        "mean_peak_abs_corr": f"{row.mean_peak_abs_corr:.10g}",
        "max_peak_abs_corr": f"{row.max_peak_abs_corr:.10g}",
        "mean_peak_corr": f"{row.mean_peak_corr:.10g}" if math.isfinite(row.mean_peak_corr) else "",
        "mean_peak_lag_sec": f"{row.mean_peak_lag_sec:.10g}" if math.isfinite(row.mean_peak_lag_sec) else "",
        "source_csv": row.source_csv,
    }


def summarize_cross_session(scores: Sequence[P10ScoreRow]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str, str], list[P10ScoreRow]] = defaultdict(list)
    for row in scores:
        if row.source_level != "per_density_feature":
            continue
        grouped[(row.run_tag, row.p9_method_k, row.density_name, row.bold_feature, row.p9_feature)].append(row)

    out: list[dict[str, object]] = []
    for (run_tag, p9_method_k, density, bold_feature, p9_feature), rows in sorted(grouped.items()):
        vals = [r.mean_peak_abs_corr for r in rows if math.isfinite(r.mean_peak_abs_corr)]
        datasets = sorted({r.dataset for r in rows})
        if not vals:
            continue
        out.append(
            {
                "run_tag": run_tag,
                "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                "p9_feature": p9_feature,
                "p9_method_k": p9_method_k,
                "density_name": density,
                "density_display": density_display_name(density),
                **density_metadata_fields(density),
                "bold_feature": bold_feature,
                "feature_family": feature_family(bold_feature),
                "n_datasets": len(datasets),
                "datasets": ";".join(datasets),
                "mean_peak_abs_corr_across_datasets": f"{mean(vals):.10g}",
                "std_peak_abs_corr_across_datasets": f"{pstdev(vals):.10g}" if len(vals) > 1 else "0",
                "min_peak_abs_corr": f"{min(vals):.10g}",
                "max_peak_abs_corr": f"{max(vals):.10g}",
                "range_peak_abs_corr": f"{(max(vals) - min(vals)):.10g}",
            }
        )
    return out


def best_hit_rows(hits: Sequence[dict[str, object]], dimred_only: bool = False) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in hits:
        if row.get("source_level") != "per_density_feature":
            continue
        if dimred_only and row.get("density_source_kind") != "dimred":
            continue
        key = (str(row.get("dataset")), str(row.get("run_tag")), str(row.get("p9_method_k")), str(row.get("feature_family")))
        grouped[key].append(dict(row))
    out: list[dict[str, object]] = []
    for _, rows in sorted(grouped.items()):
        valid = [r for r in rows if math.isfinite(as_float(str(r.get("peak_abs_corr", ""))))]
        if valid:
            out.append(max(valid, key=lambda r: as_float(str(r.get("peak_abs_corr", "")))))
    return out


def plot_score_overviews(scores: Sequence[P10ScoreRow], figure_dir: Path, run_tags: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    figure_dir.mkdir(parents=True, exist_ok=True)
    p9_cols = [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]
    row_labels = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]

    for condition in ("all", "abs", "csplit"):
        grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in scores:
            if row.source_level != "per_density_feature":
                continue
            meta = density_metadata(row.density_name)
            if condition != "all" and meta["density_condition"] != condition:
                continue
            grouped[(row.observable, row.p9_method_k)].append(row.mean_peak_abs_corr)
        matrix = {key: mean(finite(vals)) for key, vals in grouped.items() if finite(vals)}
        if matrix:
            path = figure_dir / f"p10_p9_method_k_score_overview__{condition}.png"
            plot_heatmap(matrix, row_labels, p9_cols, f"P10 P9 method x k mean top5 score | density {condition}", path, value_format=".3f")
            paths.append(path)

    lfp_cols = [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]
    for condition in ("abs", "csplit"):
        grouped = defaultdict(list)
        for row in scores:
            if row.source_level != "per_density_feature":
                continue
            meta = density_metadata(row.density_name)
            if meta["density_source_kind"] != "dimred" or meta["density_condition"] != condition:
                continue
            grouped[(row.observable, meta["density_method_k"])].append(row.mean_peak_abs_corr)
        matrix = {key: mean(finite(vals)) for key, vals in grouped.items() if finite(vals)}
        if matrix:
            path = figure_dir / f"p10_lfp_density_method_k_score_overview__{condition}.png"
            plot_heatmap(matrix, row_labels, lfp_cols, f"P10 LFP density method x k mean top5 score | {condition}", path, value_format=".3f")
            paths.append(path)

    for run_tag in run_tags:
        for condition in ("abs", "csplit"):
            grouped = defaultdict(list)
            for row in scores:
                if row.run_tag != run_tag or row.source_level != "per_density_feature":
                    continue
                meta = density_metadata(row.density_name)
                if meta["density_source_kind"] != "dimred" or meta["density_condition"] != condition:
                    continue
                grouped[(row.p9_method_k, meta["density_method_k"])].append(row.mean_peak_abs_corr)
            matrix = {key: mean(finite(vals)) for key, vals in grouped.items() if finite(vals)}
            if matrix:
                path = figure_dir / f"p10_p9_by_lfp_method_k_score__{condition}__{slug(run_tag)}.png"
                title = f"P10 P9 method-k x LFP density method-k | {condition} | {RUN_TAG_LABELS.get(run_tag, run_tag)}"
                plot_heatmap(matrix, p9_cols, lfp_cols, title, path, value_format=".3f")
                paths.append(path)
    return paths


def main() -> None:
    args = parse_args()
    scores, hits = collect_scores(
        args.processed_root,
        args.datasets,
        args.run_tags,
        args.features,
        args.xcorr_save_tag,
    )
    args.results_dir.mkdir(parents=True, exist_ok=True)
    score_dicts = [row_to_dict(row) for row in scores]
    write_csv(args.results_dir / "p10_xcorr_topn_scores_long.csv", score_dicts)
    write_csv(args.results_dir / "p10_xcorr_density_feature_scores.csv", [row_to_dict(r) for r in scores if r.source_level == "per_density_feature"])
    write_csv(args.results_dir / "p10_top_xcorr_hits_readable.csv", hits)
    write_csv(args.results_dir / "p10_cross_session_density_feature_summary.csv", summarize_cross_session(scores))
    write_csv(args.results_dir / "p10_best_hit_by_dataset_observable_p9_method_family.csv", best_hit_rows(hits, dimred_only=False))
    write_csv(args.results_dir / "p10_best_dimred_density_hit_by_dataset_observable_p9_method_family.csv", best_hit_rows(hits, dimred_only=True))

    figure_paths: list[Path] = []
    if not args.skip_figures:
        figure_paths.extend(plot_score_overviews(scores, args.figure_dir, args.run_tags))
        if args.results_figure_dir.resolve() != args.figure_dir.resolve():
            figure_paths.extend(plot_score_overviews(scores, args.results_figure_dir, args.run_tags))

    print(f"Datasets: {', '.join(args.datasets)}")
    print(f"Run tags: {', '.join(args.run_tags)}")
    print(f"Features: {', '.join(args.features)}")
    print(f"Score rows: {len(scores)}")
    print(f"Readable hit rows: {len(hits)}")
    print("P10 topN: 5")
    print(f"Results dir: {args.results_dir}")
    print(f"Figure dir: {args.figure_dir}")
    print(f"Results figure dir: {args.results_figure_dir}")
    print(f"Figures: {len(figure_paths)}")


if __name__ == "__main__":
    main()
