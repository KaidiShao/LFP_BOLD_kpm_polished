"""Compare P8 raw-density top-xcorr ROI profiles for efun vs deconv_efun.

The existing ROI profile exports were built mainly for dimred/P5 subprocess
density sources.  This script reconstructs raw-density ROI profiles directly:

1. read P8 per-feature raw-density xcorr top rows;
2. use `bold_mode_index` and `peak_abs_corr` as selected BOLD-mode weights;
3. fetch each selected BOLD mode ROI vector from the P7 all-mode ROI export;
4. compare efun_real-selected ROI profiles with deconv_real-selected profiles.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np

from summarize_pipeline8_cross_session_consistency import DEFAULT_PROCESSED_ROOT, write_csv


DEFAULT_RESULTS_DIR = Path("results") / "pipeline_roi_profile_consistency_current"
DEFAULT_P7_PROFILE = DEFAULT_RESULTS_DIR / "p7_intrinsic_bold_efun_roi_profiles_long.csv"
DEFAULT_OUTPUT_DIR = DEFAULT_RESULTS_DIR / "p8_raw_density_efun_vs_deconv_roi_similarity_20260603"
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p8_raw_density_efun_vs_deconv_roi_similarity_20260603"
)

DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01")
RUN_TAG_TO_OBSERVABLE = {
    "pv_gsvd100": "global_svd100",
    "pv_gsvd100_ds": "gsvd100_ds",
    "pv_roi": "roi_mean",
    "pv_roi_mean": "roi_mean",
    "pv_hp100": "HP_svd100",
    "pv_elehp": "HP_svd100",
}
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_roi")
DEFAULT_RAW_DENSITIES = ("raw_abs_q070", "raw_csplit_q070")
FEATURE_FAMILIES = ("efun", "deconv_efun")
FEATURE_FILES = {
    "efun": "efun_real",
    "deconv_efun": "deconv_real",
}
FEATURE_COLORS = {"efun": "#2b83ba", "deconv_efun": "#d95f02"}
COMPARISON_ORDER = (
    "same_dataset_efun_vs_deconv",
    "cross_dataset_efun_within",
    "cross_dataset_deconv_efun_within",
    "cross_dataset_efun_vs_deconv",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument(
        "--p7-profile-csv",
        type=Path,
        nargs="+",
        default=[DEFAULT_P7_PROFILE],
        help="One or more P7 all-mode ROI profile CSVs. Multiple files are merged by dataset/observable/mode.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--raw-densities", nargs="+", default=list(DEFAULT_RAW_DENSITIES))
    parser.add_argument("--min-rois", type=int, default=5)
    return parser.parse_args()


def as_float(value: object) -> float:
    try:
        out = float(value) if value not in (None, "") else math.nan
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def as_int(value: object) -> int | None:
    value_f = as_float(value)
    if not math.isfinite(value_f):
        return None
    return int(value_f)


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


def read_p7_roi_profiles(
    csv_paths: Sequence[Path],
    *,
    datasets: set[str],
    observables: set[str],
) -> dict[tuple[str, str, int], dict[str, float]]:
    out: dict[tuple[str, str, int], dict[str, float]] = defaultdict(dict)
    for csv_path in csv_paths:
        if not csv_path.is_file():
            continue
        with csv_path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                dataset = row.get("dataset", "").lower()
                obs = row.get("observable", "")
                pos = as_int(row.get("sorted_position"))
                roi = row.get("roi_label", "")
                value = as_float(row.get("roi_value"))
                if dataset not in datasets:
                    continue
                if obs not in observables:
                    continue
                if pos is None or not roi or not math.isfinite(value):
                    continue
                out[(dataset, obs, pos)][roi] = value
    return dict(out)


def raw_top_file_candidates(
    processed_root: Path,
    dataset: str,
    run_tag: str,
    feature_family: str,
    density_name: str,
) -> list[Path]:
    feature_name = FEATURE_FILES[feature_family]
    base = (
        processed_root
        / dataset
        / "pipeline8_xcorr"
        / run_tag
        / "feature"
        / feature_family
    )
    density_names = [density_name]
    if density_name == "raw_csplit_q070":
        density_names.append("raw_csplit_q070_standardize_rmsenv_adaptive")
    if density_name == "raw_abs_q070":
        density_names.append("raw_abs_q070_standardize_rmsenv_adaptive")

    candidates: list[Path] = []
    for name in dict.fromkeys(density_names):
        candidates.extend(
            [
                base / f"xcorr_top__{name}__{feature_name}.csv",
                base / f"xcorr_csplit_standardize_rmsenv_adaptive_top__{name}__{feature_name}.csv",
                base / f"xcorr_abs_standardize_rmsenv_adaptive_top__{name}__{feature_name}.csv",
            ]
        )
    return candidates


def read_raw_top_rows(
    processed_root: Path,
    *,
    datasets: Sequence[str],
    run_tags: Sequence[str],
    raw_densities: Sequence[str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    missing: list[dict[str, object]] = []
    for dataset in datasets:
        ds = dataset.lower()
        for run_tag in run_tags:
            obs = RUN_TAG_TO_OBSERVABLE.get(run_tag, run_tag)
            for density in raw_densities:
                for feature_family in FEATURE_FAMILIES:
                    candidates = raw_top_file_candidates(processed_root, ds, run_tag, feature_family, density)
                    path = next((candidate for candidate in candidates if candidate.is_file()), candidates[0])
                    if not path.is_file():
                        missing.append(
                            {
                                "dataset": ds,
                                "run_tag": run_tag,
                                "observable": obs,
                                "raw_density": density,
                                "feature_family": feature_family,
                                "path": str(path),
                                "candidate_paths": ";".join(str(candidate) for candidate in candidates),
                            }
                        )
                        continue
                    with path.open(newline="", encoding="utf-8-sig") as handle:
                        for rank, row in enumerate(csv.DictReader(handle), start=1):
                            mode = as_int(row.get("bold_mode_index"))
                            peak = as_float(row.get("peak_abs_corr"))
                            density_index = as_int(row.get("density_index"))
                            if mode is None or not math.isfinite(peak):
                                continue
                            rows.append(
                                {
                                    "dataset": ds,
                                    "run_tag": run_tag,
                                    "observable": obs,
                                    "raw_density": density,
                                    "feature_family": feature_family,
                                    "bold_feature": row.get("bold_feature", FEATURE_FILES[feature_family]),
                                    "bold_mode_index": mode,
                                    "density_index": density_index if density_index is not None else "",
                                    "density_label": row.get("density_label", ""),
                                    "peak_abs_corr": peak,
                                    "peak_corr": row.get("peak_corr", ""),
                                    "peak_lag_sec": row.get("peak_lag_sec", ""),
                                    "rank_in_top_file": rank,
                                    "source_csv": str(path),
                                }
                            )
    return rows, missing


def build_weighted_roi_vectors(
    top_rows: Sequence[dict[str, object]],
    p7_profiles: dict[tuple[str, str, int], dict[str, float]],
) -> tuple[
    dict[tuple[str, str, str, str], dict[str, float]],
    dict[tuple[str, str, str, str], dict[int, float]],
    list[dict[str, object]],
]:
    roi_accum: dict[tuple[str, str, str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0]))
    index_weights: dict[tuple[str, str, str, str], dict[int, float]] = defaultdict(lambda: defaultdict(float))
    used_rows: list[dict[str, object]] = []
    for row in top_rows:
        dataset = str(row["dataset"])
        obs = str(row["observable"])
        density = str(row["raw_density"])
        feature = str(row["feature_family"])
        mode = int(row["bold_mode_index"])
        peak = float(row["peak_abs_corr"])
        roi_vec = p7_profiles.get((dataset, obs, mode), {})
        if not roi_vec:
            continue
        key = (dataset, obs, density, feature)
        index_weights[key][mode] += peak
        for roi, value in roi_vec.items():
            slot = roi_accum[key][roi]
            slot[0] += value * peak
            slot[1] += peak
        item = dict(row)
        item["n_rois_from_p7"] = len(roi_vec)
        used_rows.append(item)

    vectors: dict[tuple[str, str, str, str], dict[str, float]] = {}
    for key, roi_map in roi_accum.items():
        vectors[key] = {
            roi: sw[0] / sw[1]
            for roi, sw in roi_map.items()
            if sw[1] > 0 and math.isfinite(sw[0])
        }
    return vectors, index_weights, used_rows


def compare_vectors(
    vectors: dict[tuple[str, str, str, str], dict[str, float]],
    index_weights: dict[tuple[str, str, str, str], dict[int, float]],
    *,
    datasets: Sequence[str],
    observables: Sequence[str],
    raw_densities: Sequence[str],
    min_rois: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    roi_rows: list[dict[str, object]] = []
    index_rows: list[dict[str, object]] = []
    ds_list = [d.lower() for d in datasets]

    for obs in observables:
        for density in raw_densities:
            for ds in ds_list:
                e_key = (ds, obs, density, "efun")
                d_key = (ds, obs, density, "deconv_efun")
                if e_key in vectors and d_key in vectors:
                    r, n_rois = corr_for_vectors(vectors[e_key], vectors[d_key], min_rois)
                    if math.isfinite(r):
                        roi_rows.append(
                            {
                                "comparison": "same_dataset_efun_vs_deconv",
                                "dataset_a": ds,
                                "dataset_b": ds,
                                "observable": obs,
                                "raw_density": density,
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
                            "raw_density": density,
                            "unique_jaccard": jaccard(eh.keys(), dh.keys()),
                            "weighted_cosine": weighted_cosine(eh, dh),
                            "n_efun_indices": len(eh),
                            "n_deconv_indices": len(dh),
                            "n_shared_indices": len(set(eh) & set(dh)),
                        }
                    )

            for feature in FEATURE_FAMILIES:
                for d1, d2 in combinations(ds_list, 2):
                    key_a = (d1, obs, density, feature)
                    key_b = (d2, obs, density, feature)
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
                                "raw_density": density,
                                "feature_a": feature,
                                "feature_b": feature,
                                "n_rois": n_rois,
                                "roi_corr": r,
                            }
                        )

            for d1, d2 in combinations(ds_list, 2):
                for feature_a, feature_b in (("efun", "deconv_efun"), ("deconv_efun", "efun")):
                    key_a = (d1, obs, density, feature_a)
                    key_b = (d2, obs, density, feature_b)
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
                                "raw_density": density,
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
        arr = np.asarray(vals, dtype=float)
        item = {field: value for field, value in zip(group_fields, key)}
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
    rows = sorted({(str(r["observable"]), str(r["raw_density"])) for r in summary_rows})
    if not rows:
        return
    lookup = {
        (str(r["observable"]), str(r["raw_density"]), str(r["comparison"])): r
        for r in summary_rows
    }
    matrix = np.full((len(rows), len(COMPARISON_ORDER)), np.nan, dtype=float)
    annot = [["" for _ in COMPARISON_ORDER] for _ in rows]
    for i, (obs, density) in enumerate(rows):
        for j, comp in enumerate(COMPARISON_ORDER):
            row = lookup.get((obs, density, comp))
            if not row:
                continue
            value = as_float(row.get("mean"))
            if math.isfinite(value):
                matrix[i, j] = value
                annot[i][j] = f"{value:.2f}"

    fig, ax = plt.subplots(figsize=(15.5, max(4.3, 0.6 * len(rows) + 2.0)))
    im = ax.imshow(matrix, cmap="coolwarm", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(COMPARISON_ORDER)))
    ax.set_xticklabels(
        ["same dataset\nefun vs deconv", "cross dataset\nefun", "cross dataset\ndeconv", "cross dataset\nefun vs deconv"]
    )
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{obs} | {density}" for obs, density in rows])
    ax.set_title("P8 raw-density top-xcorr ROI profile similarity: efun_real vs deconv_real")
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
    panels = sorted({(str(r["observable"]), str(r["raw_density"])) for r in rows})
    if not panels:
        return
    n_cols = 2
    n_rows = int(math.ceil(len(panels) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18.0, max(4.0, 3.5 * n_rows)), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")
    for ax, (obs, density) in zip(axes.ravel(), panels):
        ax.axis("on")
        grouped: dict[str, list[float]] = defaultdict(list)
        for row in rows:
            if row["observable"] == obs and row["raw_density"] == density:
                grouped[str(row["comparison"])].append(float(row["roi_corr"]))
        data = [grouped.get(comp, []) for comp in COMPARISON_ORDER]
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
        ax.set_xticks(range(1, len(COMPARISON_ORDER) + 1))
        ax.set_xticklabels(["same ds\nE-D", "cross ds\nE-E", "cross ds\nD-D", "cross ds\nE-D"])
        ax.set_ylim(-0.55, 1.02)
        ax.set_title(f"{obs} | {density}")
        ax.set_ylabel("ROI profile corr")
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("P8 raw-density top-xcorr ROI profile distributions", fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def make_index_summary(index_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for metric in ("unique_jaccard", "weighted_cosine"):
        metric_rows = summarize_rows(index_rows, metric, ("observable", "raw_density"))
        for row in metric_rows:
            row["metric"] = metric
            rows.append(row)
    rows.sort(key=lambda r: (str(r.get("observable", "")), str(r.get("raw_density", "")), str(r.get("metric", ""))))
    return rows


def plot_index_overlap_heatmap(summary_rows: Sequence[dict[str, object]], out_file: Path) -> None:
    rows = sorted({(str(r["observable"]), str(r["raw_density"])) for r in summary_rows})
    metrics = ("unique_jaccard", "weighted_cosine")
    if not rows:
        return
    lookup = {
        (str(r["observable"]), str(r["raw_density"]), str(r["metric"])): r
        for r in summary_rows
    }
    matrix = np.full((len(rows), len(metrics)), np.nan, dtype=float)
    for i, (obs, density) in enumerate(rows):
        for j, metric in enumerate(metrics):
            row = lookup.get((obs, density, metric))
            if row:
                matrix[i, j] = as_float(row.get("mean"))
    fig, ax = plt.subplots(figsize=(12.5, max(4.0, 0.58 * len(rows) + 1.8)))
    im = ax.imshow(matrix, cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(range(len(metrics)))
    ax.set_xticklabels(["mode-index\nJaccard", "weighted-index\ncosine"])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{obs} | {density}" for obs, density in rows])
    ax.set_title("Raw-density efun vs deconv selected BOLD mode-index overlap")
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
    index_weights: dict[tuple[str, str, str, str], dict[int, float]],
    observables: Sequence[str],
    raw_densities: Sequence[str],
    figure_dir: Path,
) -> None:
    for obs in observables:
        for density in raw_densities:
            e_total: dict[int, float] = defaultdict(float)
            d_total: dict[int, float] = defaultdict(float)
            for (dataset, obs_i, density_i, feature), hist in index_weights.items():
                if obs_i != obs or density_i != density:
                    continue
                target = e_total if feature == "efun" else d_total
                for idx, weight in hist.items():
                    target[idx] += weight
            if not e_total and not d_total:
                continue
            positions = sorted(set(e_total) | set(d_total))
            x = np.arange(len(positions))
            e_vals = np.asarray([e_total.get(pos, 0.0) for pos in positions], dtype=float)
            d_vals = np.asarray([d_total.get(pos, 0.0) for pos in positions], dtype=float)
            fig, ax = plt.subplots(figsize=(18.0, 3.8))
            ax.bar(x, e_vals, color=FEATURE_COLORS["efun"], alpha=0.72, label="efun_real")
            ax.bar(x, -d_vals, color=FEATURE_COLORS["deconv_efun"], alpha=0.72, label="deconv_real")
            ax.axhline(0, color="#303030", linewidth=0.8)
            tick_step = max(1, int(math.ceil(len(positions) / 24)))
            ticks = np.arange(0, len(positions), tick_step)
            ax.set_xticks(ticks)
            ax.set_xticklabels([f"{positions[i]:03d}" for i in ticks], rotation=45, ha="right")
            ax.set_ylabel("selected weight\npositive efun / negative deconv")
            ax.set_title(f"P8 raw-density selected BOLD mode positions | {obs} | {density}")
            ax.legend(frameon=False)
            ax.grid(axis="y", alpha=0.22)
            fig.tight_layout()
            figure_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(figure_dir / f"04_raw_density_selected_mode_histogram__{safe_name(obs)}__{safe_name(density)}.png", dpi=180)
            plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    datasets = {d.lower() for d in args.datasets}
    run_tags = list(args.run_tags)
    observables = {RUN_TAG_TO_OBSERVABLE.get(tag, tag) for tag in run_tags}
    p7_profiles = read_p7_roi_profiles(args.p7_profile_csv, datasets=datasets, observables=observables)
    top_rows, missing = read_raw_top_rows(
        args.processed_root,
        datasets=args.datasets,
        run_tags=run_tags,
        raw_densities=args.raw_densities,
    )
    vectors, index_weights, used_rows = build_weighted_roi_vectors(top_rows, p7_profiles)
    roi_rows, index_rows = compare_vectors(
        vectors,
        index_weights,
        datasets=args.datasets,
        observables=sorted(observables),
        raw_densities=args.raw_densities,
        min_rois=args.min_rois,
    )
    roi_summary = summarize_rows(roi_rows, "roi_corr", ("observable", "raw_density", "comparison"))
    index_summary = make_index_summary(index_rows)

    write_csv(args.output_dir / "p8_raw_density_top_rows_used_for_roi_profiles.csv", used_rows)
    write_csv(args.output_dir / "p8_raw_density_top_rows_missing.csv", missing)
    write_csv(args.output_dir / "p8_raw_density_efun_vs_deconv_roi_pairwise.csv", roi_rows)
    write_csv(args.output_dir / "p8_raw_density_efun_vs_deconv_roi_summary.csv", roi_summary)
    write_csv(args.output_dir / "p8_raw_density_efun_vs_deconv_mode_index_overlap.csv", index_rows)
    write_csv(args.output_dir / "p8_raw_density_efun_vs_deconv_mode_index_overlap_summary.csv", index_summary)

    plot_roi_summary_heatmap(roi_summary, args.figure_dir / "01_p8_raw_density_efun_vs_deconv_roi_similarity_summary.png")
    plot_roi_boxplots(roi_rows, args.figure_dir / "02_p8_raw_density_efun_vs_deconv_roi_similarity_boxplots.png")
    plot_index_overlap_heatmap(index_summary, args.figure_dir / "03_p8_raw_density_efun_vs_deconv_selected_mode_overlap.png")
    plot_index_histograms(index_weights, sorted(observables), args.raw_densities, args.figure_dir)

    print(f"P7 mode ROI vectors: {len(p7_profiles)}")
    print(f"Raw top rows read: {len(top_rows)}")
    print(f"Raw top rows used: {len(used_rows)}")
    print(f"Missing raw top files: {len(missing)}")
    print(f"ROI comparison rows: {len(roi_rows)}")
    print(f"Mode-index overlap rows: {len(index_rows)}")
    print(f"Output dir: {args.output_dir}")
    print(f"Figure dir: {args.figure_dir}")


if __name__ == "__main__":
    main()
