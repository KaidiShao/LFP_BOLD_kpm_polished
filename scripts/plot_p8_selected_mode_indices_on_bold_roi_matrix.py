"""Overlay P8 subprocess-selected BOLD mode indices on P7 ROI matrices.

This tests whether P8/P5 subprocess coupling is acting as a selector over
high-consistency regions of the intrinsic P7 BOLD efun ROI profile matrix.

The P7 matrix is independent of BLP density.  P8 contributes only the selected
BOLD sorted-mode indices for each dataset / observable / BOLD feature family /
strict subprocess label.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np

from summarize_pipeline8_cross_session_consistency import DEFAULT_PROCESSED_ROOT, write_csv


DEFAULT_RESULTS_DIR = Path("results") / "pipeline_roi_profile_consistency_current"
DEFAULT_P7_PAIRWISE = (
    DEFAULT_RESULTS_DIR
    / "p7_intrinsic_bold_efun_roi_confusion_20260603"
    / "p7_intrinsic_bold_efun_roi_pairwise_correlations.csv"
)
DEFAULT_P8_PROFILE = DEFAULT_RESULTS_DIR / "p8_roi_profiles_by_subprocess_long.csv"
DEFAULT_OUTPUT_DIR = DEFAULT_RESULTS_DIR / "p8_selected_mode_index_on_p7_roi_matrix_20260603"
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p8_selected_mode_index_on_p7_roi_matrix_20260603"
)

DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01")
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_roi")
DEFAULT_OBSERVABLES = ("global_svd100", "gsvd100_ds", "roi_mean")
DEFAULT_FEATURES = ("efun", "deconv_efun")
LABEL_ORDER = ("theta", "ripple_gamma")
LABEL_COLORS = {
    "theta": "#2ca25f",
    "ripple_gamma": "#756bb1",
    "theta_within": "#2ca25f",
    "rg_within": "#756bb1",
    "theta_vs_rg": "#d95f02",
    "allmode": "#8c8c8c",
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
RESTRICTED_GROUPS = ("allmode", "theta_within", "rg_within", "theta_vs_rg")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p7-pairwise-csv", type=Path, default=DEFAULT_P7_PAIRWISE)
    parser.add_argument("--p8-profile-csv", type=Path, default=DEFAULT_P8_PROFILE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--observables", nargs="+", default=list(DEFAULT_OBSERVABLES))
    parser.add_argument("--features", nargs="+", default=list(DEFAULT_FEATURES))
    parser.add_argument("--conditions", nargs="+", default=["csplit"])
    parser.add_argument("--label-groups", nargs="+", default=list(LABEL_ORDER))
    parser.add_argument("--method-k-scope", choices=["all", "p5-stable"], default="p5-stable")
    parser.add_argument("--min-rois", type=int, default=5)
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
    density_mk = row.get("density_method_k", "")
    p9_mk = row.get("p9_method_k", "")
    if density_mk and density_mk not in STABLE_METHOD_K:
        return False
    if p9_mk and p9_mk not in STABLE_METHOD_K:
        return False
    return True


def summarize_values(values: Sequence[float]) -> dict[str, object]:
    clean = [float(v) for v in values if math.isfinite(float(v))]
    if not clean:
        return {"n": 0, "mean": "", "median": "", "q25": "", "q75": "", "min": "", "max": ""}
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


def read_selected_indices(
    path: Path,
    *,
    datasets: set[str],
    run_tags: set[str],
    observables: set[str],
    features: set[str],
    conditions: set[str],
    labels: set[str],
    method_k_scope: str,
) -> tuple[
    dict[tuple[str, str, str, str], set[int]],
    dict[tuple[str, str, str], dict[int, float]],
    list[dict[str, object]],
]:
    """Read P8 profile rows and return dataset-specific selected mode sets.

    The input CSV has one row per ROI, so rows are deduplicated at the
    dataset/source/subprocess level before index counts are accumulated.
    """

    selected_sets: dict[tuple[str, str, str, str], set[int]] = defaultdict(set)
    index_weights: dict[tuple[str, str, str], dict[int, float]] = defaultdict(lambda: defaultdict(float))
    source_rows: list[dict[str, object]] = []
    seen_sources: set[tuple[str, ...]] = set()

    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row.get("pipeline", "").upper() != "P8":
                continue
            dataset = row.get("dataset", "").lower()
            run_tag = row.get("run_tag", "").lower()
            obs = row.get("observable", "")
            feature = row.get("feature_family", "")
            condition = row.get("density_condition", "").lower()
            label = row.get("strict_label_group", "")
            if dataset not in datasets:
                continue
            if run_tag not in run_tags:
                continue
            if obs not in observables:
                continue
            if feature not in features:
                continue
            if condition not in conditions:
                continue
            if label not in labels:
                continue
            if row.get("density_source_kind", "") != "dimred":
                continue
            if not row_method_k_ok(row, method_k_scope):
                continue

            indices = parse_indices(row.get("selected_indices", ""))
            if not indices:
                continue

            source_key = (
                dataset,
                run_tag,
                obs,
                feature,
                row.get("bold_feature", ""),
                row.get("density_name", ""),
                row.get("density_method_k", ""),
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

            peak_weight = as_float(row.get("subprocess_peak_abs_corr_sum"))
            if not math.isfinite(peak_weight) or peak_weight <= 0:
                peak_weight = as_float(row.get("mean_peak_abs_corr"))
            if not math.isfinite(peak_weight) or peak_weight <= 0:
                peak_weight = 1.0

            selected_sets[(obs, feature, label, dataset)].update(indices)
            denom = max(1, len(indices))
            for idx in indices:
                index_weights[(obs, feature, label)][idx] += peak_weight / denom

            source_rows.append(
                {
                    "dataset": dataset,
                    "run_tag": row.get("run_tag", ""),
                    "observable": obs,
                    "feature_family": feature,
                    "bold_feature": row.get("bold_feature", ""),
                    "density_name": row.get("density_name", ""),
                    "density_method_k": row.get("density_method_k", ""),
                    "density_threshold": row.get("density_threshold", ""),
                    "strict_label_group": label,
                    "strict_label": row.get("strict_label", ""),
                    "n_indices": len(indices),
                    "unique_n_indices": len(set(indices)),
                    "selected_indices": row.get("selected_indices", ""),
                    "subprocess_peak_abs_corr_sum": row.get("subprocess_peak_abs_corr_sum", ""),
                    "mean_peak_abs_corr": row.get("mean_peak_abs_corr", ""),
                }
            )

    return selected_sets, index_weights, source_rows


def read_p7_pairwise(
    path: Path,
    *,
    datasets: set[str],
    observables: set[str],
    min_rois: int,
) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row.get("comparison") not in {"within_mode_cross_dataset", "cross_mode_cross_dataset"}:
                continue
            dataset_a = row.get("dataset_a", "").lower()
            dataset_b = row.get("dataset_b", "").lower()
            obs = row.get("observable", "")
            if dataset_a not in datasets or dataset_b not in datasets:
                continue
            if obs not in observables:
                continue
            n_rois = as_int(row.get("n_rois")) or 0
            if n_rois < min_rois:
                continue
            pos_a = as_int(row.get("mode_a_pos"))
            pos_b = as_int(row.get("mode_b_pos"))
            value = as_float(row.get("roi_corr"))
            if pos_a is None or pos_b is None or not math.isfinite(value):
                continue
            records.append(
                {
                    "observable": obs,
                    "comparison": row.get("comparison", ""),
                    "dataset_a": dataset_a,
                    "dataset_b": dataset_b,
                    "mode_a": row.get("mode_a", ""),
                    "mode_b": row.get("mode_b", ""),
                    "mode_a_pos": pos_a,
                    "mode_b_pos": pos_b,
                    "n_rois": n_rois,
                    "roi_corr": value,
                }
            )
    return records


def mean_matrix(records: Sequence[dict[str, object]], obs: str) -> tuple[np.ndarray, list[int]]:
    positions = sorted(
        {
            int(row["mode_a_pos"])
            for row in records
            if row.get("observable") == obs
        }
        | {
            int(row["mode_b_pos"])
            for row in records
            if row.get("observable") == obs
        }
    )
    pos_to_i = {pos: i for i, pos in enumerate(positions)}
    buckets: dict[tuple[int, int], list[float]] = defaultdict(list)
    for row in records:
        if row.get("observable") != obs:
            continue
        a = int(row["mode_a_pos"])
        b = int(row["mode_b_pos"])
        value = float(row["roi_corr"])
        buckets[(a, b)].append(value)
        buckets[(b, a)].append(value)
    matrix = np.full((len(positions), len(positions)), np.nan, dtype=float)
    for (a, b), vals in buckets.items():
        matrix[pos_to_i[a], pos_to_i[b]] = mean(vals)
    return matrix, positions


def collect_restricted_pairwise(
    p7_records: Sequence[dict[str, object]],
    selected_sets: dict[tuple[str, str, str, str], set[int]],
    *,
    observables: Sequence[str],
    features: Sequence[str],
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for row in p7_records:
        obs = str(row["observable"])
        if obs not in observables:
            continue
        da = str(row["dataset_a"])
        db = str(row["dataset_b"])
        pa = int(row["mode_a_pos"])
        pb = int(row["mode_b_pos"])
        value = float(row["roi_corr"])
        for feature in features:
            theta_a = selected_sets.get((obs, feature, "theta", da), set())
            theta_b = selected_sets.get((obs, feature, "theta", db), set())
            rg_a = selected_sets.get((obs, feature, "ripple_gamma", da), set())
            rg_b = selected_sets.get((obs, feature, "ripple_gamma", db), set())
            memberships = {
                "allmode": True,
                "theta_within": pa in theta_a and pb in theta_b,
                "rg_within": pa in rg_a and pb in rg_b,
                "theta_vs_rg": (pa in theta_a and pb in rg_b) or (pa in rg_a and pb in theta_b),
            }
            for group, keep in memberships.items():
                if not keep:
                    continue
                out.append(
                    {
                        "observable": obs,
                        "feature_family": feature,
                        "restriction_group": group,
                        "dataset_a": da,
                        "dataset_b": db,
                        "mode_a": row["mode_a"],
                        "mode_b": row["mode_b"],
                        "mode_a_pos": pa,
                        "mode_b_pos": pb,
                        "roi_corr": value,
                        "n_rois": row["n_rois"],
                    }
                )
    return out


def make_index_frequency_rows(
    index_weights: dict[tuple[str, str, str], dict[int, float]],
    selected_sets: dict[tuple[str, str, str, str], set[int]],
    datasets: Sequence[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (obs, feature, label), weights in sorted(index_weights.items()):
        for idx, weight in sorted(weights.items()):
            n_datasets = sum(1 for ds in datasets if idx in selected_sets.get((obs, feature, label, ds.lower()), set()))
            rows.append(
                {
                    "observable": obs,
                    "feature_family": feature,
                    "strict_label_group": label,
                    "mode_pos": idx,
                    "weighted_count": f"{weight:.10g}",
                    "n_datasets_selected": n_datasets,
                    "datasets_selected": ";".join(
                        ds for ds in datasets if idx in selected_sets.get((obs, feature, label, ds.lower()), set())
                    ),
                }
            )
    return rows


def make_summary_rows(records: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    for row in records:
        key = (
            str(row["observable"]),
            str(row["feature_family"]),
            str(row["restriction_group"]),
        )
        grouped[key].append(float(row["roi_corr"]))

    rows: list[dict[str, object]] = []
    for key, values in sorted(grouped.items()):
        obs, feature, group = key
        item: dict[str, object] = {
            "observable": obs,
            "feature_family": feature,
            "restriction_group": group,
        }
        item.update(summarize_values(values))
        rows.append(item)

    # Add explicit same/different block contrasts.
    by_key = {(r["observable"], r["feature_family"], r["restriction_group"]): r for r in rows}
    for obs in sorted({str(r["observable"]) for r in rows}):
        for feature in sorted({str(r["feature_family"]) for r in rows if r["observable"] == obs}):
            theta = as_float(by_key.get((obs, feature, "theta_within"), {}).get("mean"))
            rg = as_float(by_key.get((obs, feature, "rg_within"), {}).get("mean"))
            cross = as_float(by_key.get((obs, feature, "theta_vs_rg"), {}).get("mean"))
            if not (math.isfinite(theta) and math.isfinite(rg) and math.isfinite(cross)):
                continue
            rows.append(
                {
                    "observable": obs,
                    "feature_family": feature,
                    "restriction_group": "within_minus_theta_vs_rg",
                    "n": "",
                    "mean": f"{((theta + rg) / 2.0 - cross):.10g}",
                    "median": "",
                    "q25": "",
                    "q75": "",
                    "min": "",
                    "max": "",
                }
            )
    return rows


def plot_overlay_matrix(
    matrix: np.ndarray,
    positions: Sequence[int],
    weights_by_label: dict[str, dict[int, float]],
    title: str,
    out_file: Path,
) -> None:
    n = len(positions)
    if n == 0:
        return
    theta = np.asarray([weights_by_label.get("theta", {}).get(pos, 0.0) for pos in positions], dtype=float)
    rg = np.asarray([weights_by_label.get("ripple_gamma", {}).get(pos, 0.0) for pos in positions], dtype=float)
    max_hist = max(float(np.max(theta)) if theta.size else 0.0, float(np.max(rg)) if rg.size else 0.0, 1.0)

    fig = plt.figure(figsize=(18.0, 10.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 8.0], width_ratios=[8.0, 1.2], hspace=0.04, wspace=0.04)
    ax_hist = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0])
    ax_side = fig.add_subplot(gs[1, 1], sharey=ax)
    ax_legend = fig.add_subplot(gs[0, 1])
    ax_legend.axis("off")

    im = ax.imshow(matrix, cmap="coolwarm", vmin=0.0, vmax=1.0, origin="upper", aspect="auto")
    ax.set_title(title, fontsize=12, weight="bold")
    tick_step = max(1, int(math.ceil(n / 16)))
    ticks = np.arange(0, n, tick_step)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{positions[i]:03d}" for i in ticks], rotation=45, ha="right", fontsize=8)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{positions[i]:03d}" for i in ticks], fontsize=8)
    ax.set_xlabel("P7 sorted BOLD mode position")
    ax.set_ylabel("P7 sorted BOLD mode position")

    x = np.arange(n)
    ax_hist.bar(x, theta, color=LABEL_COLORS["theta"], alpha=0.72, width=0.86, label="theta selected")
    ax_hist.bar(x, -rg, color=LABEL_COLORS["ripple_gamma"], alpha=0.72, width=0.86, label="RG selected")
    ax_hist.axhline(0, color="#303030", linewidth=0.8)
    ax_hist.set_xlim(-0.5, n - 0.5)
    ax_hist.set_ylim(-max_hist * 1.08, max_hist * 1.08)
    ax_hist.set_xticks([])
    ax_hist.set_ylabel("selected\nweight", fontsize=9)
    ax_hist.legend(loc="upper right", fontsize=8, frameon=False)
    ax_hist.grid(axis="y", alpha=0.2)

    ax_side.barh(x, theta, color=LABEL_COLORS["theta"], alpha=0.72, height=0.86)
    ax_side.barh(x, -rg, color=LABEL_COLORS["ripple_gamma"], alpha=0.72, height=0.86)
    ax_side.axvline(0, color="#303030", linewidth=0.8)
    ax_side.set_xlim(-max_hist * 1.08, max_hist * 1.08)
    ax_side.set_xticks([])
    ax_side.set_yticks([])
    ax_side.set_xlabel("weight", fontsize=8)

    cbar = fig.colorbar(im, ax=ax_legend, fraction=0.75)
    cbar.set_label("mean cross-dataset ROI corr")
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def plot_restricted_boxplots(records: Sequence[dict[str, object]], out_file: Path) -> None:
    keys = sorted({(str(r["observable"]), str(r["feature_family"])) for r in records})
    if not keys:
        return
    n_cols = 2
    n_rows = int(math.ceil(len(keys) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18.0, max(4.2, 3.6 * n_rows)), squeeze=False)
    for ax in axes.ravel():
        ax.axis("off")

    for ax, (obs, feature) in zip(axes.ravel(), keys):
        ax.axis("on")
        groups: dict[str, list[float]] = defaultdict(list)
        for row in records:
            if row["observable"] == obs and row["feature_family"] == feature:
                groups[str(row["restriction_group"])].append(float(row["roi_corr"]))
        data = [groups.get(group, []) for group in RESTRICTED_GROUPS]
        positions = np.arange(1, len(RESTRICTED_GROUPS) + 1)
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=0.56,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.1},
        )
        for patch, group in zip(bp["boxes"], RESTRICTED_GROUPS):
            patch.set_facecolor(LABEL_COLORS.get(group, "#cccccc"))
            patch.set_alpha(0.78)
            patch.set_edgecolor("#404040")
        for pos, vals, group in zip(positions, data, RESTRICTED_GROUPS):
            clean = [float(v) for v in vals if math.isfinite(float(v))]
            if not clean:
                continue
            idx = np.linspace(0, len(clean) - 1, min(250, len(clean))).astype(int)
            jitter = np.linspace(-0.15, 0.15, len(idx))
            ax.scatter(
                pos + jitter,
                np.asarray(clean)[idx],
                s=6,
                color=LABEL_COLORS.get(group, "#777777"),
                alpha=0.16,
                linewidths=0,
            )
            ax.text(pos, 0.98, f"n={len(clean)}", ha="center", va="top", fontsize=7, rotation=90)
        ax.axhline(0.5, color="#303030", linestyle="--", linewidth=0.9, alpha=0.7)
        ax.set_xticks(positions)
        ax.set_xticklabels(["all\nmode", "theta\nwithin", "RG\nwithin", "theta\nvs RG"])
        ax.set_ylim(-0.55, 1.02)
        ax.set_ylabel("cross-dataset ROI corr")
        ax.set_title(f"{obs} | {feature}")
        ax.grid(axis="y", alpha=0.25)

    fig.suptitle("P8 selected mode-index restricted pairs on P7 BOLD ROI matrix", fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def plot_summary_heatmap(summary_rows: Sequence[dict[str, object]], out_file: Path) -> None:
    rows = sorted({(str(r["observable"]), str(r["feature_family"])) for r in summary_rows})
    cols = list(RESTRICTED_GROUPS) + ["within_minus_theta_vs_rg"]
    if not rows:
        return
    matrix = np.full((len(rows), len(cols)), np.nan, dtype=float)
    annot = [["" for _ in cols] for _ in rows]
    lookup = {
        (str(r["observable"]), str(r["feature_family"]), str(r["restriction_group"])): r
        for r in summary_rows
    }
    for i, (obs, feature) in enumerate(rows):
        for j, col in enumerate(cols):
            row = lookup.get((obs, feature, col))
            if not row:
                continue
            value = as_float(row.get("mean"))
            if math.isfinite(value):
                matrix[i, j] = value
                annot[i][j] = f"{value:.2f}"

    fig, ax = plt.subplots(figsize=(14.5, max(4.5, 0.55 * len(rows) + 2.2)))
    im = ax.imshow(matrix, cmap="coolwarm", vmin=0.0, vmax=0.8, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(["allmode", "theta\nwithin", "RG\nwithin", "theta\nvs RG", "within\nminus cross"])
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{obs} | {feature}" for obs, feature in rows])
    ax.set_title("Mean P7 ROI corr after restricting to P8-selected BOLD mode indices")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if math.isfinite(float(matrix[i, j])):
                ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.8, label="mean ROI corr / contrast")
    fig.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    datasets = {d.lower() for d in args.datasets}
    run_tags = {r.lower() for r in args.run_tags}
    observables = set(args.observables)
    features = set(args.features)
    conditions = {c.lower() for c in args.conditions}
    labels = set(args.label_groups)

    selected_sets, index_weights, source_rows = read_selected_indices(
        args.p8_profile_csv,
        datasets=datasets,
        run_tags=run_tags,
        observables=observables,
        features=features,
        conditions=conditions,
        labels=labels,
        method_k_scope=args.method_k_scope,
    )
    p7_records = read_p7_pairwise(
        args.p7_pairwise_csv,
        datasets=datasets,
        observables=observables,
        min_rois=args.min_rois,
    )
    restricted = collect_restricted_pairwise(
        p7_records,
        selected_sets,
        observables=args.observables,
        features=args.features,
    )
    summary = make_summary_rows(restricted)
    index_rows = make_index_frequency_rows(index_weights, selected_sets, args.datasets)

    write_csv(args.output_dir / "p8_selected_bold_mode_index_source_groups.csv", source_rows)
    write_csv(args.output_dir / "p8_selected_bold_mode_index_frequency.csv", index_rows)
    write_csv(args.output_dir / "p8_selected_mode_restricted_pairwise_roi_corr.csv", restricted)
    write_csv(args.output_dir / "p8_selected_mode_restricted_summary.csv", summary)

    by_obs_matrix = {obs: mean_matrix(p7_records, obs) for obs in args.observables}
    for obs, (matrix, positions) in by_obs_matrix.items():
        if not positions:
            continue
        for feature in args.features:
            weights_by_label = {
                label: index_weights.get((obs, feature, label), {})
                for label in args.label_groups
            }
            if not any(weights_by_label.values()):
                continue
            plot_overlay_matrix(
                matrix,
                positions,
                weights_by_label,
                f"P7 BOLD ROI matrix with P8-selected mode indices | {obs} | {feature}",
                args.figure_dir / f"01_p8_selected_index_overlay_on_p7_matrix__{safe_name(obs)}__{safe_name(feature)}.png",
            )

    plot_restricted_boxplots(
        restricted,
        args.figure_dir / "02_p8_selected_index_restricted_roi_corr_boxplots.png",
    )
    plot_summary_heatmap(
        summary,
        args.figure_dir / "03_p8_selected_index_restricted_roi_corr_summary.png",
    )

    print(f"P8 source groups: {len(source_rows)}")
    print(f"P7 pairwise records: {len(p7_records)}")
    print(f"Restricted pairwise records: {len(restricted)}")
    print(f"Output dir: {args.output_dir}")
    print(f"Figure dir: {args.figure_dir}")


if __name__ == "__main__":
    main()
