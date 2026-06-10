"""Compare BLP-selected BOLD ROI profiles against the all-BOLD-mode baseline.

The baseline is the P7 intrinsic all-mode BOLD efun ROI profile distribution:
all cross-dataset ROI correlations between original BOLD efun modes, split by
BOLD observable.  The selected distributions come from P8/P10 subprocess ROI
profiles grouped by P5 strict band-selectivity labels.
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

from summarize_pipeline8_cross_session_consistency import DEFAULT_PROCESSED_ROOT, RUN_TAG_LABELS, write_csv


DEFAULT_RESULTS_DIR = Path("results") / "pipeline_roi_profile_consistency_current"
DEFAULT_P7_PAIRWISE = (
    DEFAULT_RESULTS_DIR
    / "p7_intrinsic_bold_efun_roi_confusion_20260603"
    / "p7_intrinsic_bold_efun_roi_pairwise_correlations.csv"
)
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "bold_roi_selected_vs_allmode_baseline_20260603"
)

DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01", "k13m17", "k13m23")
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")
LABEL_ORDER = ("theta", "ripple_gamma", "mixed", "inactive")
LABEL_COLORS = {
    "allmode": "#8c8c8c",
    "theta": "#2ca25f",
    "ripple_gamma": "#756bb1",
    "mixed": "#b76e1f",
    "inactive": "#d9d9d9",
}
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
METHOD_ORDER = ("svd", "nmf", "mds", "umap")
K_ORDER = ("k03", "k04", "k05", "k06", "k07", "k08")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p7-pairwise-csv", type=Path, default=DEFAULT_P7_PAIRWISE)
    parser.add_argument("--p8-profile-csv", type=Path, default=DEFAULT_RESULTS_DIR / "p8_roi_profiles_by_subprocess_long.csv")
    parser.add_argument(
        "--p10-profile-csv",
        type=Path,
        default=DEFAULT_RESULTS_DIR / "p10_roi_profiles_by_subprocess_long.partial_20260603.csv",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_RESULTS_DIR / "bold_roi_selected_vs_allmode_baseline_20260603")
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--conditions", nargs="+", default=["csplit"])
    parser.add_argument("--label-groups", nargs="+", default=list(LABEL_ORDER))
    parser.add_argument("--method-k-scope", choices=["all", "p5-stable"], default="p5-stable")
    parser.add_argument("--min-datasets", type=int, default=5)
    parser.add_argument("--min-rois", type=int, default=5)
    parser.add_argument("--value-column", default="roi_value_weighted")
    parser.add_argument("--weight-column", default="mean_peak_abs_corr")
    return parser.parse_args()


def as_float(value: object) -> float:
    try:
        out = float(value) if value not in (None, "") else math.nan
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


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


def safe_name(text: object) -> str:
    out = []
    for ch in str(text):
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_") or "blank"


def split_method_k(method_k: object) -> tuple[str, str]:
    text = str(method_k)
    if "_" not in text:
        return text, ""
    return tuple(text.split("_", 1))  # type: ignore[return-value]


def label_sort_key(label: str, preferred: Sequence[str]) -> tuple[int, str]:
    try:
        return (preferred.index(label), label)
    except ValueError:
        return (len(preferred), label)


def row_method_k_ok(row: dict[str, str], scope: str) -> bool:
    if scope == "all":
        return True
    density_mk = row.get("density_method_k", "")
    p9_mk = row.get("p9_method_k", "")
    if density_mk and density_mk not in STABLE_METHOD_K:
        return False
    if p9_mk and p9_mk not in STABLE_METHOD_K:
        return False
    return True


def read_baseline_distribution(path: Path, datasets: set[str], min_rois: int) -> dict[str, list[float]]:
    values: dict[str, list[float]] = defaultdict(list)
    if not path.is_file():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row.get("comparison") not in {"within_mode_cross_dataset", "cross_mode_cross_dataset"}:
                continue
            if row.get("dataset_a", "").lower() not in datasets or row.get("dataset_b", "").lower() not in datasets:
                continue
            if int(as_float(row.get("n_rois")) or 0) < min_rois:
                continue
            value = as_float(row.get("roi_corr"))
            if math.isfinite(value):
                values[row.get("observable", "")].append(value)
    return dict(values)


def read_selected_vectors(
    path: Path,
    pipeline: str,
    datasets: set[str],
    run_tags: set[str],
    conditions: set[str],
    labels: set[str],
    method_k_scope: str,
    value_column: str,
    weight_column: str,
) -> tuple[dict[tuple[str, ...], dict[str, dict[str, float]]], dict[tuple[str, ...], dict[str, str]]]:
    sums: dict[tuple[str, ...], dict[str, dict[str, list[float]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(lambda: [0.0, 0.0]))
    )
    meta: dict[tuple[str, ...], dict[str, str]] = {}
    if not path.is_file():
        return {}, {}

    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row.get("pipeline", "").upper() != pipeline.upper():
                continue
            dataset = row.get("dataset", "")
            run_tag = row.get("run_tag", "")
            condition = row.get("density_condition", "")
            label = row.get("strict_label_group", "")
            if dataset.lower() not in datasets:
                continue
            if run_tag.lower() not in run_tags:
                continue
            if condition.lower() not in conditions:
                continue
            if label not in labels:
                continue
            if row.get("density_source_kind", "") != "dimred":
                continue
            if not row_method_k_ok(row, method_k_scope):
                continue

            roi = row.get("roi_label", "")
            value = as_float(row.get(value_column))
            if not roi or not math.isfinite(value):
                continue
            weight = as_float(row.get(weight_column))
            if not math.isfinite(weight) or weight <= 0:
                weight = 1.0

            if pipeline.upper() == "P8":
                key = (
                    "P8",
                    row.get("run_tag", ""),
                    row.get("observable", ""),
                    row.get("feature_family", ""),
                    row.get("density_method_k", ""),
                    label,
                )
                info = {
                    "pipeline": "P8",
                    "run_tag": row.get("run_tag", ""),
                    "observable": row.get("observable", ""),
                    "feature_family": row.get("feature_family", ""),
                    "p9_feature": "",
                    "p9_method_k": "",
                    "density_method_k": row.get("density_method_k", ""),
                    "strict_label_group": label,
                    "method_k_scope": method_k_scope,
                }
            else:
                key = (
                    "P10",
                    row.get("run_tag", ""),
                    row.get("observable", ""),
                    row.get("p9_feature", ""),
                    row.get("p9_method_k", ""),
                    row.get("feature_family", ""),
                    row.get("density_method_k", ""),
                    label,
                )
                info = {
                    "pipeline": "P10",
                    "run_tag": row.get("run_tag", ""),
                    "observable": row.get("observable", ""),
                    "feature_family": row.get("feature_family", ""),
                    "p9_feature": row.get("p9_feature", ""),
                    "p9_method_k": row.get("p9_method_k", ""),
                    "density_method_k": row.get("density_method_k", ""),
                    "strict_label_group": label,
                    "method_k_scope": method_k_scope,
                }

            slot = sums[key][dataset][roi]
            slot[0] += value * weight
            slot[1] += weight
            meta.setdefault(key, info)

    vectors: dict[tuple[str, ...], dict[str, dict[str, float]]] = {}
    for key, dataset_map in sums.items():
        vectors[key] = {}
        for dataset, roi_map in dataset_map.items():
            vectors[key][dataset] = {roi: sw[0] / sw[1] for roi, sw in roi_map.items() if sw[1] > 0}
    return vectors, meta


def corr_for_vectors(a: dict[str, float], b: dict[str, float], min_rois: int) -> tuple[float, int]:
    common = sorted(set(a) & set(b))
    if len(common) < min_rois:
        return math.nan, len(common)
    return pearson([a[r] for r in common], [b[r] for r in common]), len(common)


def selected_pairwise_records(
    vectors: dict[tuple[str, ...], dict[str, dict[str, float]]],
    meta: dict[tuple[str, ...], dict[str, str]],
    min_datasets: int,
    min_rois: int,
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for key, dataset_map in vectors.items():
        datasets = sorted(dataset_map)
        if len(datasets) < min_datasets:
            continue
        info = meta[key]
        for d1, d2 in combinations(datasets, 2):
            value, n_rois = corr_for_vectors(dataset_map[d1], dataset_map[d2], min_rois)
            if not math.isfinite(value):
                continue
            records.append(
                {
                    **info,
                    "dataset_a": d1,
                    "dataset_b": d2,
                    "n_datasets": len(datasets),
                    "n_rois": n_rois,
                    "roi_corr": value,
                }
            )
    return records


def summarize_values(values: Sequence[float]) -> dict[str, object]:
    clean = [v for v in values if math.isfinite(v)]
    if not clean:
        return {
            "n": 0,
            "mean": "",
            "median": "",
            "q25": "",
            "q75": "",
            "min": "",
            "max": "",
        }
    arr = np.asarray(clean, dtype=float)
    return {
        "n": len(clean),
        "mean": f"{float(np.mean(arr)):.10g}",
        "median": f"{float(np.median(arr)):.10g}",
        "q25": f"{float(np.percentile(arr, 25)):.10g}",
        "q75": f"{float(np.percentile(arr, 75)):.10g}",
        "min": f"{float(np.min(arr)):.10g}",
        "max": f"{float(np.max(arr)):.10g}",
    }


def make_summary_rows(
    baseline: dict[str, list[float]],
    selected: Sequence[dict[str, object]],
    labels: Sequence[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for obs, vals in sorted(baseline.items()):
        item = {
            "pipeline": "P7",
            "observable": obs,
            "feature_family": "allmode_baseline",
            "strict_label_group": "allmode",
            "density_method_k": "all_modes",
            "baseline_mean": "",
            "delta_vs_baseline_mean": "",
        }
        item.update(summarize_values(vals))
        rows.append(item)

    groups: dict[tuple[str, str, str, str], list[float]] = defaultdict(list)
    for row in selected:
        key = (
            str(row.get("pipeline", "")),
            str(row.get("observable", "")),
            str(row.get("feature_family", "")),
            str(row.get("strict_label_group", "")),
        )
        value = as_float(row.get("roi_corr"))
        if math.isfinite(value):
            groups[key].append(value)
    for key, vals in sorted(groups.items()):
        pipeline, obs, feature, label = key
        base_vals = baseline.get(obs, [])
        base_mean = mean(base_vals) if base_vals else math.nan
        item = {
            "pipeline": pipeline,
            "observable": obs,
            "feature_family": feature,
            "strict_label_group": label,
            "density_method_k": "all_p5_stable_method_k",
            "baseline_mean": f"{base_mean:.10g}" if math.isfinite(base_mean) else "",
        }
        item.update(summarize_values(vals))
        selected_mean = as_float(item["mean"])
        item["delta_vs_baseline_mean"] = (
            f"{selected_mean - base_mean:.10g}"
            if math.isfinite(selected_mean) and math.isfinite(base_mean)
            else ""
        )
        rows.append(item)
    return rows


def plot_distribution_panels(
    pipeline: str,
    baseline: dict[str, list[float]],
    selected: Sequence[dict[str, object]],
    labels: Sequence[str],
    figure_dir: Path,
) -> None:
    groups: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in selected:
        if row.get("pipeline") != pipeline:
            continue
        obs = str(row.get("observable", ""))
        feature = str(row.get("feature_family", ""))
        label = str(row.get("strict_label_group", ""))
        value = as_float(row.get("roi_corr"))
        if math.isfinite(value):
            groups[(obs, feature)][label].append(value)

    panel_keys = [(obs, feature) for obs, feature in sorted(groups) if obs in baseline]
    if not panel_keys:
        return
    n_cols = 2
    n_rows = int(math.ceil(len(panel_keys) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12.5, max(3.8, 3.3 * n_rows)), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")

    for ax, (obs, feature) in zip(axes.ravel(), panel_keys):
        ax.axis("on")
        names = ["allmode"] + list(labels)
        data = [baseline.get(obs, [])] + [groups[(obs, feature)].get(label, []) for label in labels]
        positions = np.arange(1, len(names) + 1)
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=0.58,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.2},
        )
        for patch, name in zip(bp["boxes"], names):
            patch.set_facecolor(LABEL_COLORS.get(name, "#cccccc"))
            patch.set_alpha(0.78)
            patch.set_edgecolor("#404040")
        for pos, vals, name in zip(positions, data, names):
            clean = [v for v in vals if math.isfinite(v)]
            if not clean:
                continue
            jitter_count = min(220, len(clean))
            idx = np.linspace(0, len(clean) - 1, jitter_count).astype(int)
            sampled = np.asarray(clean)[idx]
            jitter = np.linspace(-0.15, 0.15, jitter_count)
            ax.scatter(
                pos + jitter,
                sampled,
                s=6,
                color=LABEL_COLORS.get(name, "#999999"),
                alpha=0.18,
                linewidths=0,
            )
        ax.set_xticks(positions)
        ax.set_xticklabels(["all-mode"] + [label.replace("_", "\n") for label in labels], rotation=0)
        ax.set_ylim(-0.55, 1.02)
        ax.set_ylabel("cross-dataset ROI corr")
        ax.set_title(f"{obs} | {feature}")
        ax.grid(axis="y", alpha=0.25)
        base_mean = mean(baseline.get(obs, [])) if baseline.get(obs) else math.nan
        if math.isfinite(base_mean):
            ax.axhline(base_mean, color="#404040", linestyle="--", linewidth=1.0, alpha=0.8)
            ax.text(0.02, base_mean + 0.015, f"baseline mean {base_mean:.2f}", transform=ax.get_yaxis_transform(), fontsize=8)

    fig.suptitle(f"{pipeline} selected BOLD ROI profiles vs all-mode BOLD baseline", fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"01_{pipeline.lower()}_selected_vs_allmode_baseline_boxplots.png", dpi=180)
    plt.close(fig)


def plot_delta_heatmap(
    pipeline: str,
    baseline: dict[str, list[float]],
    selected: Sequence[dict[str, object]],
    labels: Sequence[str],
    figure_dir: Path,
) -> None:
    groups: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in selected:
        if row.get("pipeline") != pipeline:
            continue
        obs = str(row.get("observable", ""))
        feature = str(row.get("feature_family", ""))
        label = str(row.get("strict_label_group", ""))
        value = as_float(row.get("roi_corr"))
        if obs in baseline and math.isfinite(value):
            groups[(obs, feature)][label].append(value)
    row_keys = sorted(groups)
    if not row_keys:
        return
    matrix = np.full((len(row_keys), len(labels)), np.nan, dtype=float)
    annot = [["" for _ in labels] for _ in row_keys]
    for i, (obs, feature) in enumerate(row_keys):
        base_mean = mean(baseline.get(obs, []))
        for j, label in enumerate(labels):
            vals = groups[(obs, feature)].get(label, [])
            if not vals:
                continue
            val_mean = mean(vals)
            matrix[i, j] = val_mean - base_mean
            annot[i][j] = f"{val_mean:.2f}\n+{matrix[i,j]:.2f}" if matrix[i, j] >= 0 else f"{val_mean:.2f}\n{matrix[i,j]:.2f}"

    fig_h = max(4.5, 0.55 * len(row_keys) + 2.2)
    fig, ax = plt.subplots(figsize=(8.8, fig_h))
    im = ax.imshow(matrix, cmap="coolwarm", vmin=-0.15, vmax=0.35)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([label.replace("_", "\n") for label in labels])
    ax.set_yticks(range(len(row_keys)))
    ax.set_yticklabels([f"{obs} | {feature}" for obs, feature in row_keys])
    ax.set_title(f"{pipeline}: selected mean ROI corr minus all-mode baseline")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if math.isfinite(float(matrix[i, j])):
                ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.82, label="selected mean - all-mode baseline mean")
    fig.tight_layout()
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"02_{pipeline.lower()}_selected_minus_allmode_baseline_heatmap.png", dpi=180)
    plt.close(fig)


def plot_method_k_delta_heatmaps(
    pipeline: str,
    baseline: dict[str, list[float]],
    selected: Sequence[dict[str, object]],
    labels: Sequence[str],
    figure_dir: Path,
) -> None:
    groups: dict[tuple[str, str, str], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in selected:
        if row.get("pipeline") != pipeline:
            continue
        obs = str(row.get("observable", ""))
        feature = str(row.get("feature_family", ""))
        method_k = str(row.get("density_method_k", ""))
        label = str(row.get("strict_label_group", ""))
        value = as_float(row.get("roi_corr"))
        if obs in baseline and math.isfinite(value):
            groups[(obs, feature, method_k)][label].append(value)

    out_dir = figure_dir / "by_method_k"
    out_dir.mkdir(parents=True, exist_ok=True)
    panel_keys = sorted({(obs, feature) for obs, feature, _mk in groups})
    for obs, feature in panel_keys:
        subset_mks = sorted(
            {mk for (obs_i, feature_i, mk) in groups if obs_i == obs and feature_i == feature},
            key=lambda mk: (label_sort_key(split_method_k(mk)[0], METHOD_ORDER), label_sort_key(split_method_k(mk)[1], K_ORDER)),
        )
        if not subset_mks:
            continue
        base_mean = mean(baseline[obs])
        matrix = np.full((len(subset_mks), len(labels)), np.nan, dtype=float)
        for i, mk in enumerate(subset_mks):
            for j, label in enumerate(labels):
                vals = groups[(obs, feature, mk)].get(label, [])
                if vals:
                    matrix[i, j] = mean(vals) - base_mean

        fig_h = max(5.0, 0.4 * len(subset_mks) + 1.8)
        fig, ax = plt.subplots(figsize=(7.8, fig_h))
        im = ax.imshow(matrix, cmap="coolwarm", vmin=-0.15, vmax=0.35)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels([label.replace("_", "\n") for label in labels])
        ax.set_yticks(range(len(subset_mks)))
        ax.set_yticklabels(subset_mks)
        ax.set_title(f"{pipeline}: selected - all-mode baseline | {obs} | {feature}")
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix[i, j]
                if math.isfinite(float(value)):
                    ax.text(j, i, f"{value:+.2f}", ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.78, label="delta ROI corr")
        fig.tight_layout()
        fig.savefig(out_dir / f"03_{pipeline.lower()}_method_k_delta__{safe_name(obs)}__{safe_name(feature)}.png", dpi=180)
        plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    datasets = {d.lower() for d in args.datasets}
    run_tags = {r.lower() for r in args.run_tags}
    conditions = {c.lower() for c in args.conditions}
    labels = list(args.label_groups)
    label_set = set(labels)

    baseline = read_baseline_distribution(args.p7_pairwise_csv, datasets, args.min_rois)
    all_selected_records: list[dict[str, object]] = []
    coverage_rows: list[dict[str, object]] = []
    for pipeline, path in (("P8", args.p8_profile_csv), ("P10", args.p10_profile_csv)):
        vectors, meta = read_selected_vectors(
            path,
            pipeline,
            datasets,
            run_tags,
            conditions,
            label_set,
            args.method_k_scope,
            args.value_column,
            args.weight_column,
        )
        for key, dataset_map in vectors.items():
            info = meta[key]
            coverage_rows.append(
                {
                    **info,
                    "n_datasets_available": len(dataset_map),
                    "datasets_available": ";".join(sorted(dataset_map)),
                    "passes_min_datasets": int(len(dataset_map) >= args.min_datasets),
                }
            )
        records = selected_pairwise_records(vectors, meta, args.min_datasets, args.min_rois)
        all_selected_records.extend(records)
        plot_distribution_panels(pipeline, baseline, records, labels, args.figure_dir)
        plot_delta_heatmap(pipeline, baseline, records, labels, args.figure_dir)
        plot_method_k_delta_heatmaps(pipeline, baseline, records, labels, args.figure_dir)

    write_csv(args.output_dir / "selected_vs_allmode_pairwise_correlations.csv", all_selected_records)
    write_csv(args.output_dir / "selected_vs_allmode_coverage.csv", coverage_rows)
    write_csv(args.output_dir / "selected_vs_allmode_summary.csv", make_summary_rows(baseline, all_selected_records, labels))

    print(f"Baseline observables: {sorted(baseline)}")
    print(f"Selected pairwise records: {len(all_selected_records)}")
    print(f"Coverage rows: {len(coverage_rows)}")
    print(f"Output dir: {args.output_dir}")
    print(f"Figure dir: {args.figure_dir}")


if __name__ == "__main__":
    main()
