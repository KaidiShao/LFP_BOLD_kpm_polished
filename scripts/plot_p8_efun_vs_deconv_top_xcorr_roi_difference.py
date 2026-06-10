"""Compare P8 efun vs deconv_efun top-xcorr ROI profiles.

This analysis ignores theta-vs-RG spatial separation and asks a cleaner
question: do BOLD efun top hits and BOLD deconv_efun top hits select different
BOLD ROI profiles / BOLD sorted-mode indices?
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from statistics import mean, median
from typing import Sequence

import matplotlib.pyplot as plt
import numpy as np

from summarize_pipeline8_cross_session_consistency import DEFAULT_PROCESSED_ROOT, write_csv


DEFAULT_RESULTS_DIR = Path("results") / "pipeline_roi_profile_consistency_current"
DEFAULT_PROFILE_CSV = DEFAULT_RESULTS_DIR / "p8_roi_profiles_by_subprocess_long.csv"
DEFAULT_OUTPUT_DIR = DEFAULT_RESULTS_DIR / "p8_efun_vs_deconv_top_xcorr_roi_difference_20260603"
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p8_efun_vs_deconv_top_xcorr_roi_difference_20260603"
)

DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01")
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_roi")
DEFAULT_OBSERVABLES = ("global_svd100", "gsvd100_ds", "roi_mean")
FEATURES = ("efun", "deconv_efun")
LABEL_SCOPES = ("ripple_gamma", "all_subprocess")
STABLE_METHOD_K = (
    "mds_k04",
    "mds_k05",
    "mds_k06",
    "mds_k07",
    "mds_k08",
    "nmf_k04",
    "umap_k04",
    "umap_k05",
    "umap_k06",
    "umap_k08",
)
LABEL_SOURCE_ORDER = ("theta", "ripple_gamma", "mixed", "inactive")
FEATURE_COLORS = {"efun": "#2b83ba", "deconv_efun": "#d95f02"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-csv", type=Path, default=DEFAULT_PROFILE_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--observables", nargs="+", default=list(DEFAULT_OBSERVABLES))
    parser.add_argument("--conditions", nargs="+", default=["csplit"])
    parser.add_argument("--label-scopes", nargs="+", default=list(LABEL_SCOPES))
    parser.add_argument("--method-k-scope", choices=["all", "p5-stable"], default="p5-stable")
    parser.add_argument("--min-rois", type=int, default=5)
    parser.add_argument("--value-column", default="roi_value_weighted")
    parser.add_argument("--weight-column", default="subprocess_peak_abs_corr_sum")
    return parser.parse_args()


def as_float(value: object) -> float:
    try:
        out = float(value) if value not in (None, "") else math.nan
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def as_int(value: object) -> int | None:
    number = as_float(value)
    if not math.isfinite(number):
        return None
    return int(number)


def safe_name(text: object) -> str:
    out = []
    for ch in str(text):
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_") or "blank"


def pearson(xs: Sequence[float], ys: Sequence[float]) -> float:
    pairs = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if len(pairs) < 2:
        return math.nan
    x = np.asarray([p[0] for p in pairs], dtype=float)
    y = np.asarray([p[1] for p in pairs], dtype=float)
    x = x - np.mean(x)
    y = y - np.mean(y)
    denom = float(np.sqrt(np.sum(x * x) * np.sum(y * y)))
    if denom <= 0:
        return math.nan
    return float(np.sum(x * y) / denom)


def corr_for_vectors(a: dict[str, float], b: dict[str, float], min_rois: int) -> tuple[float, int]:
    common = sorted(set(a) & set(b))
    if len(common) < min_rois:
        return math.nan, len(common)
    return pearson([a[x] for x in common], [b[x] for x in common]), len(common)


def parse_indices(text: object) -> list[int]:
    out: list[int] = []
    for part in str(text or "").replace(",", ";").split(";"):
        idx = as_int(part.strip())
        if idx is not None and idx >= 1:
            out.append(idx)
    return out


def row_method_k_ok(row: dict[str, str], scope: str) -> bool:
    if scope == "all":
        return True
    return row.get("density_method_k", "") in STABLE_METHOD_K


def row_label_scopes(label: str, allowed_scopes: set[str]) -> list[str]:
    scopes = []
    if "all_subprocess" in allowed_scopes and label in LABEL_SOURCE_ORDER:
        scopes.append("all_subprocess")
    if label in allowed_scopes:
        scopes.append(label)
    return scopes


def add_weighted_value(
    store: dict[tuple[str, ...], dict[str, list[float]]],
    key: tuple[str, ...],
    roi: str,
    value: float,
    weight: float,
) -> None:
    slot = store[key][roi]
    slot[0] += value * weight
    slot[1] += weight


def read_profiles_and_indices(
    csv_path: Path,
    *,
    datasets: set[str],
    run_tags: set[str],
    observables: set[str],
    conditions: set[str],
    label_scopes: set[str],
    method_k_scope: str,
    value_column: str,
    weight_column: str,
) -> tuple[
    dict[tuple[str, ...], dict[str, float]],
    dict[tuple[str, ...], dict[int, float]],
    list[dict[str, object]],
]:
    roi_sums: dict[tuple[str, ...], dict[str, list[float]]] = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0]))
    index_weights: dict[tuple[str, ...], dict[int, float]] = defaultdict(lambda: defaultdict(float))
    source_rows: list[dict[str, object]] = []
    seen_sources: set[tuple[str, ...]] = set()

    with csv_path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row.get("pipeline", "").upper() != "P8":
                continue
            dataset = row.get("dataset", "").lower()
            run_tag = row.get("run_tag", "").lower()
            obs = row.get("observable", "")
            feature = row.get("feature_family", "")
            condition = row.get("density_condition", "").lower()
            label = row.get("strict_label_group", "")
            density_mk = row.get("density_method_k", "")
            if dataset not in datasets:
                continue
            if run_tag not in run_tags:
                continue
            if obs not in observables:
                continue
            if feature not in FEATURES:
                continue
            if condition not in conditions:
                continue
            if row.get("density_source_kind", "") != "dimred":
                continue
            if not row_method_k_ok(row, method_k_scope):
                continue
            scopes = row_label_scopes(label, label_scopes)
            if not scopes:
                continue

            roi = row.get("roi_label", "")
            value = as_float(row.get(value_column))
            if not roi or not math.isfinite(value):
                continue
            weight = as_float(row.get(weight_column))
            if not math.isfinite(weight) or weight <= 0:
                weight = as_float(row.get("mean_peak_abs_corr"))
            if not math.isfinite(weight) or weight <= 0:
                weight = 1.0

            for scope in scopes:
                all_key = (dataset, obs, scope, feature, "all_p5_stable_method_k")
                mk_key = (dataset, obs, scope, feature, density_mk)
                add_weighted_value(roi_sums, all_key, roi, value, weight)
                add_weighted_value(roi_sums, mk_key, roi, value, weight)

            source_key = (
                dataset,
                run_tag,
                obs,
                feature,
                row.get("bold_feature", ""),
                row.get("density_name", ""),
                density_mk,
                row.get("density_threshold", ""),
                label,
                row.get("strict_label", ""),
                row.get("subprocess_rank", ""),
                row.get("subprocess_density_indices", ""),
                row.get("selected_indices", ""),
            )
            if source_key in seen_sources:
                continue
            seen_sources.add(source_key)

            indices = parse_indices(row.get("selected_indices", ""))
            if indices:
                source_weight = as_float(row.get(weight_column))
                if not math.isfinite(source_weight) or source_weight <= 0:
                    source_weight = as_float(row.get("mean_peak_abs_corr"))
                if not math.isfinite(source_weight) or source_weight <= 0:
                    source_weight = 1.0
                for scope in scopes:
                    for mk in ("all_p5_stable_method_k", density_mk):
                        hist_key = (dataset, obs, scope, feature, mk)
                        for idx in indices:
                            index_weights[hist_key][idx] += source_weight / max(1, len(indices))
                source_rows.append(
                    {
                        "dataset": dataset,
                        "observable": obs,
                        "feature_family": feature,
                        "strict_label_group": label,
                        "density_method_k": density_mk,
                        "n_indices": len(indices),
                        "unique_n_indices": len(set(indices)),
                        "selected_indices": row.get("selected_indices", ""),
                    }
                )

    vectors: dict[tuple[str, ...], dict[str, float]] = {}
    for key, roi_map in roi_sums.items():
        vectors[key] = {
            roi: sw[0] / sw[1]
            for roi, sw in roi_map.items()
            if sw[1] > 0 and math.isfinite(sw[0])
        }
    return vectors, index_weights, source_rows


def weighted_cosine(a: dict[int, float], b: dict[int, float]) -> float:
    keys = set(a) | set(b)
    if not keys:
        return math.nan
    av = np.asarray([a.get(k, 0.0) for k in keys], dtype=float)
    bv = np.asarray([b.get(k, 0.0) for k in keys], dtype=float)
    denom = float(np.sqrt(np.sum(av * av) * np.sum(bv * bv)))
    if denom <= 0:
        return math.nan
    return float(np.sum(av * bv) / denom)


def jaccard(a: Iterable[int], b: Iterable[int]) -> float:
    aa = set(a)
    bb = set(b)
    union = aa | bb
    if not union:
        return math.nan
    return len(aa & bb) / len(union)


def compare_features(
    vectors: dict[tuple[str, ...], dict[str, float]],
    index_weights: dict[tuple[str, ...], dict[int, float]],
    *,
    datasets: Sequence[str],
    observables: Sequence[str],
    label_scopes: Sequence[str],
    min_rois: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    roi_rows: list[dict[str, object]] = []
    index_rows: list[dict[str, object]] = []
    method_ks = sorted({key[4] for key in vectors if len(key) == 5})

    for obs in observables:
        for scope in label_scopes:
            for mk in method_ks:
                for dataset in datasets:
                    ds = dataset.lower()
                    e_key = (ds, obs, scope, "efun", mk)
                    d_key = (ds, obs, scope, "deconv_efun", mk)
                    if e_key in vectors and d_key in vectors:
                        r, n_rois = corr_for_vectors(vectors[e_key], vectors[d_key], min_rois)
                        if math.isfinite(r):
                            roi_rows.append(
                                {
                                    "comparison": "same_dataset_efun_vs_deconv",
                                    "dataset_a": ds,
                                    "dataset_b": ds,
                                    "observable": obs,
                                    "label_scope": scope,
                                    "density_method_k": mk,
                                    "feature_a": "efun",
                                    "feature_b": "deconv_efun",
                                    "n_rois": n_rois,
                                    "roi_corr": r,
                                }
                            )

                    eh = index_weights.get(e_key, {})
                    dh = index_weights.get(d_key, {})
                    if eh or dh:
                        index_rows.append(
                            {
                                "comparison": "same_dataset_efun_vs_deconv",
                                "dataset": ds,
                                "observable": obs,
                                "label_scope": scope,
                                "density_method_k": mk,
                                "unique_jaccard": jaccard(eh.keys(), dh.keys()),
                                "weighted_cosine": weighted_cosine(eh, dh),
                                "n_efun_indices": len(eh),
                                "n_deconv_indices": len(dh),
                                "n_shared_indices": len(set(eh) & set(dh)),
                            }
                        )

                # Cross-dataset reference distributions.
                for feature in FEATURES:
                    for d1, d2 in combinations([d.lower() for d in datasets], 2):
                        key_a = (d1, obs, scope, feature, mk)
                        key_b = (d2, obs, scope, feature, mk)
                        if key_a not in vectors or key_b not in vectors:
                            continue
                        r, n_rois = corr_for_vectors(vectors[key_a], vectors[key_b], min_rois)
                        if math.isfinite(r):
                            roi_rows.append(
                                {
                                    "comparison": f"cross_dataset_{feature}_within",
                                    "dataset_a": d1,
                                    "dataset_b": d2,
                                    "observable": obs,
                                    "label_scope": scope,
                                    "density_method_k": mk,
                                    "feature_a": feature,
                                    "feature_b": feature,
                                    "n_rois": n_rois,
                                    "roi_corr": r,
                                }
                            )
                for d1, d2 in combinations([d.lower() for d in datasets], 2):
                    for feature_a, feature_b in (("efun", "deconv_efun"), ("deconv_efun", "efun")):
                        key_a = (d1, obs, scope, feature_a, mk)
                        key_b = (d2, obs, scope, feature_b, mk)
                        if key_a not in vectors or key_b not in vectors:
                            continue
                        r, n_rois = corr_for_vectors(vectors[key_a], vectors[key_b], min_rois)
                        if math.isfinite(r):
                            roi_rows.append(
                                {
                                    "comparison": "cross_dataset_efun_vs_deconv",
                                    "dataset_a": d1,
                                    "dataset_b": d2,
                                    "observable": obs,
                                    "label_scope": scope,
                                    "density_method_k": mk,
                                    "feature_a": feature_a,
                                    "feature_b": feature_b,
                                    "n_rois": n_rois,
                                    "roi_corr": r,
                                }
                            )
    return roi_rows, index_rows


def summarize_rows(rows: Sequence[dict[str, object]], value_field: str, group_fields: Sequence[str]) -> list[dict[str, object]]:
    grouped: dict[tuple[object, ...], list[float]] = defaultdict(list)
    for row in rows:
        value = as_float(row.get(value_field))
        if not math.isfinite(value):
            continue
        key = tuple(row.get(field, "") for field in group_fields)
        grouped[key].append(value)
    out: list[dict[str, object]] = []
    for key, vals in sorted(grouped.items()):
        item = {field: value for field, value in zip(group_fields, key)}
        arr = np.asarray(vals, dtype=float)
        item.update(
            {
                "n": len(vals),
                "mean": f"{float(np.mean(arr)):.10g}",
                "median": f"{float(np.median(arr)):.10g}",
                "q25": f"{float(np.percentile(arr, 25)):.10g}",
                "q75": f"{float(np.percentile(arr, 75)):.10g}",
                "min": f"{float(np.min(arr)):.10g}",
                "max": f"{float(np.max(arr)):.10g}",
            }
        )
        out.append(item)
    return out


def plot_roi_summary_heatmap(summary_rows: Sequence[dict[str, object]], out_file: Path) -> None:
    wanted = (
        "same_dataset_efun_vs_deconv",
        "cross_dataset_efun_within",
        "cross_dataset_deconv_efun_within",
        "cross_dataset_efun_vs_deconv",
    )
    rows = sorted({(str(r["observable"]), str(r["label_scope"])) for r in summary_rows})
    if not rows:
        return
    lookup = {
        (str(r["observable"]), str(r["label_scope"]), str(r["comparison"])): r
        for r in summary_rows
    }
    matrix = np.full((len(rows), len(wanted)), np.nan, dtype=float)
    annot = [["" for _ in wanted] for _ in rows]
    for i, (obs, scope) in enumerate(rows):
        for j, comp in enumerate(wanted):
            row = lookup.get((obs, scope, comp))
            if not row:
                continue
            value = as_float(row.get("mean"))
            if math.isfinite(value):
                matrix[i, j] = value
                annot[i][j] = f"{value:.2f}"

    fig, ax = plt.subplots(figsize=(15.5, max(4.2, 0.62 * len(rows) + 2.0)))
    im = ax.imshow(matrix, cmap="coolwarm", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(wanted)))
    ax.set_xticklabels(
        ["same dataset\nefun vs deconv", "cross dataset\nefun", "cross dataset\ndeconv", "cross dataset\nefun vs deconv"],
        rotation=0,
    )
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{obs} | {scope}" for obs, scope in rows])
    ax.set_title("P8 top-xcorr ROI profile similarity: efun vs deconv_efun")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if math.isfinite(float(matrix[i, j])):
                ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.82, label="ROI profile Pearson corr")
    fig.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def plot_roi_boxplots(rows: Sequence[dict[str, object]], out_file: Path) -> None:
    panel_keys = sorted({(str(r["observable"]), str(r["label_scope"])) for r in rows})
    comparisons = (
        "same_dataset_efun_vs_deconv",
        "cross_dataset_efun_within",
        "cross_dataset_deconv_efun_within",
        "cross_dataset_efun_vs_deconv",
    )
    if not panel_keys:
        return
    n_cols = 2
    n_rows = int(math.ceil(len(panel_keys) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18.0, max(4.0, 3.6 * n_rows)), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, (obs, scope) in zip(axes.ravel(), panel_keys):
        ax.axis("on")
        grouped: dict[str, list[float]] = defaultdict(list)
        for row in rows:
            if row["observable"] == obs and row["label_scope"] == scope:
                grouped[str(row["comparison"])].append(float(row["roi_corr"]))
        data = [grouped.get(comp, []) for comp in comparisons]
        bp = ax.boxplot(data, patch_artist=True, showfliers=False, medianprops={"color": "black", "linewidth": 1.1})
        colors = ("#d95f02", "#2b83ba", "#d95f02", "#999999")
        for patch, color in zip(bp["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.72)
            patch.set_edgecolor("#404040")
        for pos, vals in enumerate(data, start=1):
            clean = [v for v in vals if math.isfinite(v)]
            if clean:
                ax.text(pos, 0.98, f"n={len(clean)}", ha="center", va="top", fontsize=7, rotation=90)
        ax.axhline(0.5, color="#303030", linestyle="--", linewidth=0.9, alpha=0.7)
        ax.set_xticks(range(1, len(comparisons) + 1))
        ax.set_xticklabels(["same ds\nE-D", "cross ds\nE-E", "cross ds\nD-D", "cross ds\nE-D"])
        ax.set_ylim(-0.55, 1.02)
        ax.set_title(f"{obs} | {scope}")
        ax.set_ylabel("ROI profile corr")
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("P8 efun vs deconv_efun top-xcorr ROI profile distributions", fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def plot_index_overlap_heatmap(summary_rows: Sequence[dict[str, object]], out_file: Path) -> None:
    rows = sorted({(str(r["observable"]), str(r["label_scope"])) for r in summary_rows})
    cols = ("unique_jaccard", "weighted_cosine")
    if not rows:
        return
    lookup = {
        (str(r["observable"]), str(r["label_scope"]), str(r["metric"])): r
        for r in summary_rows
    }
    matrix = np.full((len(rows), len(cols)), np.nan, dtype=float)
    for i, (obs, scope) in enumerate(rows):
        for j, metric in enumerate(cols):
            row = lookup.get((obs, scope, metric))
            if row:
                matrix[i, j] = as_float(row.get("mean"))
    fig, ax = plt.subplots(figsize=(12.5, max(4.0, 0.58 * len(rows) + 1.8)))
    im = ax.imshow(matrix, cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(["mode-index\nJaccard", "weighted-index\ncosine"])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{obs} | {scope}" for obs, scope in rows])
    ax.set_title("Efun vs deconv selected BOLD mode-index overlap")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if math.isfinite(float(matrix[i, j])):
                ax.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=9, color="white" if matrix[i, j] < 0.55 else "black")
    fig.colorbar(im, ax=ax, shrink=0.82, label="overlap")
    fig.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def plot_index_histograms(
    index_weights: dict[tuple[str, ...], dict[int, float]],
    observables: Sequence[str],
    label_scopes: Sequence[str],
    figure_dir: Path,
) -> None:
    for obs in observables:
        for scope in label_scopes:
            e_total: dict[int, float] = defaultdict(float)
            d_total: dict[int, float] = defaultdict(float)
            for key, hist in index_weights.items():
                dataset, obs_i, scope_i, feature, mk = key
                if obs_i != obs or scope_i != scope or mk != "all_p5_stable_method_k":
                    continue
                target = e_total if feature == "efun" else d_total
                for idx, weight in hist.items():
                    target[idx] += weight
            if not e_total and not d_total:
                continue
            positions = sorted(set(e_total) | set(d_total))
            e_vals = np.asarray([e_total.get(pos, 0.0) for pos in positions], dtype=float)
            d_vals = np.asarray([d_total.get(pos, 0.0) for pos in positions], dtype=float)
            fig, ax = plt.subplots(figsize=(18.0, 3.8))
            x = np.arange(len(positions))
            ax.bar(x, e_vals, color=FEATURE_COLORS["efun"], alpha=0.72, label="efun")
            ax.bar(x, -d_vals, color=FEATURE_COLORS["deconv_efun"], alpha=0.72, label="deconv_efun")
            ax.axhline(0, color="#303030", linewidth=0.8)
            tick_step = max(1, int(math.ceil(len(positions) / 24)))
            ticks = np.arange(0, len(positions), tick_step)
            ax.set_xticks(ticks)
            ax.set_xticklabels([f"{positions[i]:03d}" for i in ticks], rotation=45, ha="right")
            ax.set_ylabel("selected weight\npositive efun / negative deconv")
            ax.set_title(f"P8 selected BOLD mode positions | {obs} | {scope}")
            ax.legend(frameon=False)
            ax.grid(axis="y", alpha=0.22)
            fig.tight_layout()
            figure_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(figure_dir / f"04_selected_mode_histogram__{safe_name(obs)}__{safe_name(scope)}.png", dpi=180)
            plt.close(fig)


def make_index_metric_summary(index_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for metric in ("unique_jaccard", "weighted_cosine"):
        metric_summary = summarize_rows(
            index_rows,
            metric,
            ("observable", "label_scope"),
        )
        for row in metric_summary:
            row["metric"] = metric
            rows.append(row)
    rows.sort(key=lambda r: (str(r.get("observable", "")), str(r.get("label_scope", "")), str(r.get("metric", ""))))
    return rows


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    datasets = {d.lower() for d in args.datasets}
    run_tags = {r.lower() for r in args.run_tags}
    observables = set(args.observables)
    conditions = {c.lower() for c in args.conditions}
    label_scopes = set(args.label_scopes)

    vectors, index_weights, source_rows = read_profiles_and_indices(
        args.profile_csv,
        datasets=datasets,
        run_tags=run_tags,
        observables=observables,
        conditions=conditions,
        label_scopes=label_scopes,
        method_k_scope=args.method_k_scope,
        value_column=args.value_column,
        weight_column=args.weight_column,
    )
    roi_rows, index_rows = compare_features(
        vectors,
        index_weights,
        datasets=args.datasets,
        observables=args.observables,
        label_scopes=args.label_scopes,
        min_rois=args.min_rois,
    )
    roi_summary = summarize_rows(
        roi_rows,
        "roi_corr",
        ("observable", "label_scope", "comparison"),
    )
    index_summary = make_index_metric_summary(index_rows)

    write_csv(args.output_dir / "p8_efun_vs_deconv_source_groups.csv", source_rows)
    write_csv(args.output_dir / "p8_efun_vs_deconv_roi_pairwise.csv", roi_rows)
    write_csv(args.output_dir / "p8_efun_vs_deconv_roi_summary.csv", roi_summary)
    write_csv(args.output_dir / "p8_efun_vs_deconv_mode_index_overlap.csv", index_rows)
    write_csv(args.output_dir / "p8_efun_vs_deconv_mode_index_overlap_summary.csv", index_summary)

    plot_roi_summary_heatmap(
        roi_summary,
        args.figure_dir / "01_p8_efun_vs_deconv_roi_similarity_summary.png",
    )
    plot_roi_boxplots(
        roi_rows,
        args.figure_dir / "02_p8_efun_vs_deconv_roi_similarity_boxplots.png",
    )
    plot_index_overlap_heatmap(
        index_summary,
        args.figure_dir / "03_p8_efun_vs_deconv_selected_mode_overlap.png",
    )
    plot_index_histograms(
        index_weights,
        args.observables,
        args.label_scopes,
        args.figure_dir,
    )

    print(f"Vectors: {len(vectors)}")
    print(f"Source groups: {len(source_rows)}")
    print(f"ROI comparison rows: {len(roi_rows)}")
    print(f"Mode-index overlap rows: {len(index_rows)}")
    print(f"Output dir: {args.output_dir}")
    print(f"Figure dir: {args.figure_dir}")


if __name__ == "__main__":
    main()
