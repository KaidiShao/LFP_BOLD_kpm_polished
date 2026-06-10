"""Targeted ROI footprint plots for the current roi_mean SOP narrative.

This script compares the ROI footprints selected by:

1. BOLD efun similar theta-selective BLP dimred density
2. BOLD efun similar RG-no-theta BLP dimred density
3. BOLD deconv_efun similar RG-no-theta BLP dimred density

The goal is to test whether the proposed slow-state footprint and fast
perturbation footprint are spatially distinct, rather than merely selecting
the same high-consistency BOLD ROI subspace.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATASETS = (
    "e10gb1",
    "e10fV1",
    "e10gh1",
    "e10gw1",
    "f12m01",
    "f12m02",
    "f12m03",
    "f12m05",
    "k13m17",
    "k13m18",
    "k13m19",
    "k13m20",
    "k13m21",
    "k13m23",
)
TARGETS = (
    ("efun_theta", "efun", "theta_selective", "efun x theta"),
    ("efun_rg", "efun", "ripple_gamma_no_theta", "efun x RG-no-theta"),
    ("deconv_rg", "deconv_efun", "ripple_gamma_no_theta", "deconv x RG-no-theta"),
)
CONTRASTS = (
    ("efun_rg_minus_efun_theta", "efun_rg", "efun_theta", "efun RG - efun theta"),
    ("deconv_rg_minus_efun_theta", "deconv_rg", "efun_theta", "deconv RG - efun theta"),
    ("deconv_rg_minus_efun_rg", "deconv_rg", "efun_rg", "deconv RG - efun RG"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--roi-test-dir",
        type=Path,
        default=Path("results")
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "roi_mean_rg_no_theta_formal_roi_tests",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("/mnt/e/DataPons_processed")
        / "summary_figures"
        / "pipeline11_current_analysis_summary"
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "roi_mean_target_footprint_tests",
    )
    parser.add_argument("--top-ns", type=int, nargs="+", default=[5, 10, 20])
    parser.add_argument("--main-top-n", type=int, default=10)
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS))
    return parser.parse_args()


def safe_name(text: object) -> str:
    out = []
    for ch in str(text):
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_") or "blank"


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
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def row_normalize(mat: np.ndarray) -> np.ndarray:
    out = np.asarray(mat, dtype=float).copy()
    for i in range(out.shape[0]):
        row = out[i]
        finite = np.isfinite(row)
        if not finite.any():
            continue
        lo = float(np.nanmin(row))
        hi = float(np.nanmax(row))
        if hi > lo:
            out[i, finite] = (row[finite] - lo) / (hi - lo)
        else:
            out[i, finite] = 0.0
    return out


def row_zscore(mat: np.ndarray) -> np.ndarray:
    out = np.asarray(mat, dtype=float).copy()
    for i in range(out.shape[0]):
        row = out[i]
        finite = np.isfinite(row)
        if finite.sum() < 3:
            continue
        mu = float(np.nanmean(row))
        sd = float(np.nanstd(row))
        if sd > 0:
            out[i, finite] = (row[finite] - mu) / sd
    return out


def load_profile_vectors(path: Path, datasets: Sequence[str]) -> tuple[list[str], dict[tuple[str, str, str], dict[str, np.ndarray]]]:
    df = pd.read_csv(path)
    df["dataset"] = df["dataset"].str.lower()
    df["roi_index"] = pd.to_numeric(df["roi_index"], errors="coerce")
    df["roi_value_weighted"] = pd.to_numeric(df["roi_value_weighted"], errors="coerce")
    roi_order = (
        df.sort_values("roi_index")[["roi_label", "roi_index"]]
        .drop_duplicates("roi_label", keep="first")
        .sort_values("roi_index")["roi_label"]
        .astype(str)
        .tolist()
    )
    roi_to_pos = {roi: i for i, roi in enumerate(roi_order)}
    out: dict[tuple[str, str, str], dict[str, np.ndarray]] = {}
    for (pipeline, feature, label, dataset), g in df.groupby(
        ["pipeline", "feature_family", "strict_label", "dataset"], sort=False
    ):
        if dataset not in {d.lower() for d in datasets}:
            continue
        vec = np.full(len(roi_order), np.nan, dtype=float)
        for _, row in g.iterrows():
            pos = roi_to_pos.get(str(row["roi_label"]))
            if pos is not None:
                vec[pos] = row["roi_value_weighted"]
        out.setdefault((pipeline, feature, label), {})[dataset] = vec
    return roi_order, out


def target_matrix(vectors: dict[tuple[str, str, str], dict[str, np.ndarray]], pipeline: str, target_id: str, datasets: Sequence[str]) -> np.ndarray:
    _, feature, label, _ = next(t for t in TARGETS if t[0] == target_id)
    vmap = vectors.get((pipeline, feature, label), {})
    n_roi = len(next(iter(vmap.values()))) if vmap else 0
    rows = []
    for ds in datasets:
        rows.append(vmap.get(ds.lower(), np.full(n_roi, np.nan)))
    return np.vstack(rows) if rows and n_roi else np.zeros((len(datasets), 0))


def plot_target_heatmaps(
    roi_order: Sequence[str],
    vectors: dict[tuple[str, str, str], dict[str, np.ndarray]],
    datasets: Sequence[str],
    figure_dir: Path,
    pipeline: str,
    top_n: int,
) -> Path:
    mats = [(target_id, label, row_normalize(target_matrix(vectors, pipeline, target_id, datasets))) for target_id, _, _, label in TARGETS]
    fig_w = max(22, len(roi_order) * 0.55)
    fig_h = max(8, len(datasets) * 0.65)
    fig, axes = plt.subplots(len(TARGETS), 1, figsize=(fig_w, fig_h), sharex=True, constrained_layout=True)
    if len(TARGETS) == 1:
        axes = [axes]
    for ax, (target_id, label, mat) in zip(axes, mats):
        im = ax.imshow(mat, aspect="auto", cmap="viridis", vmin=0, vmax=1)
        ax.set_title(f"{pipeline} {label} ROI footprint | top{top_n} | row-normalized", weight="bold")
        ax.set_yticks(range(len(datasets)))
        ax.set_yticklabels(datasets)
        ax.set_ylabel("dataset")
        fig.colorbar(im, ax=ax, shrink=0.55, label="within-row normalized ROI value")
    axes[-1].set_xticks(range(len(roi_order)))
    axes[-1].set_xticklabels(roi_order, rotation=75, ha="right", fontsize=7)
    axes[-1].set_xlabel("ROI, original P7 roi_mean order")
    figure_dir.mkdir(parents=True, exist_ok=True)
    path = figure_dir / f"01_target_roi_footprints_by_dataset__{pipeline.lower()}__top{top_n}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_contrast_heatmaps(
    roi_order: Sequence[str],
    vectors: dict[tuple[str, str, str], dict[str, np.ndarray]],
    datasets: Sequence[str],
    figure_dir: Path,
    pipeline: str,
    top_n: int,
) -> Path:
    z_mats = {
        target_id: row_zscore(target_matrix(vectors, pipeline, target_id, datasets))
        for target_id, _, _, _ in TARGETS
    }
    contrasts = []
    for contrast_id, a, b, label in CONTRASTS:
        contrasts.append((contrast_id, label, z_mats[a] - z_mats[b]))
    fig_w = max(22, len(roi_order) * 0.55)
    fig_h = max(8, len(datasets) * 0.65)
    fig, axes = plt.subplots(len(CONTRASTS), 1, figsize=(fig_w, fig_h), sharex=True, constrained_layout=True)
    if len(CONTRASTS) == 1:
        axes = [axes]
    for ax, (_, label, mat) in zip(axes, contrasts):
        vmax = np.nanpercentile(np.abs(mat), 97) if np.isfinite(mat).any() else 1
        vmax = max(float(vmax), 0.5)
        im = ax.imshow(mat, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax)
        ax.set_title(f"{pipeline} ROI contrast: {label} | top{top_n} | per-target row z-scored", weight="bold")
        ax.set_yticks(range(len(datasets)))
        ax.set_yticklabels(datasets)
        ax.set_ylabel("dataset")
        fig.colorbar(im, ax=ax, shrink=0.55, label="relative ROI contrast")
    axes[-1].set_xticks(range(len(roi_order)))
    axes[-1].set_xticklabels(roi_order, rotation=75, ha="right", fontsize=7)
    axes[-1].set_xlabel("ROI, original P7 roi_mean order")
    path = figure_dir / f"02_target_roi_footprint_contrasts__{pipeline.lower()}__top{top_n}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def extract_similarity_rows(roi_test_dir: Path, top_n: int) -> list[dict[str, object]]:
    pair_path = roi_test_dir / f"roi_pairwise_label_feature_confusion__top{top_n}.csv"
    base_path = roi_test_dir / f"roi_selected_vs_allmode_baseline_summary__top{top_n}.csv"
    rows: list[dict[str, object]] = []
    pair = pd.read_csv(pair_path)
    base = pd.read_csv(base_path)
    target_defs = [
        ("cross_dataset", "efun_theta_within", "efun", "theta_selective", "efun", "theta_selective"),
        ("cross_dataset", "efun_rg_within", "efun", "ripple_gamma_no_theta", "efun", "ripple_gamma_no_theta"),
        ("cross_dataset", "deconv_rg_within", "deconv_efun", "ripple_gamma_no_theta", "deconv_efun", "ripple_gamma_no_theta"),
        ("same_dataset", "efun_theta_vs_efun_rg", "efun", "theta_selective", "efun", "ripple_gamma_no_theta"),
        ("same_dataset", "efun_theta_vs_deconv_rg", "efun", "theta_selective", "deconv_efun", "ripple_gamma_no_theta"),
        ("same_dataset", "efun_rg_vs_deconv_rg", "efun", "ripple_gamma_no_theta", "deconv_efun", "ripple_gamma_no_theta"),
    ]
    for pipeline in ("P8", "P10"):
        for comparison, metric, fa, la, fb, lb in target_defs:
            sub = pair[
                (pair["pipeline"] == pipeline)
                & (pair["comparison"] == comparison)
                & (pair["feature_a"] == fa)
                & (pair["label_a"] == la)
                & (pair["feature_b"] == fb)
                & (pair["label_b"] == lb)
            ]
            if sub.empty:
                continue
            r = sub.iloc[0].to_dict()
            rows.append(
                {
                    "top_n": top_n,
                    "pipeline": pipeline,
                    "metric": metric,
                    "comparison": comparison,
                    "n_corr": r.get("n_corr"),
                    "mean_corr": r.get("mean_corr"),
                    "median_corr": r.get("median_corr"),
                }
            )
    for _, r in base.iterrows():
        rows.append(
            {
                "top_n": top_n,
                "pipeline": "baseline",
                "metric": r["distribution"],
                "comparison": "cross_dataset",
                "n_corr": r["n_corr"],
                "mean_corr": r["mean_corr"],
                "median_corr": r["median_corr"],
            }
        )
    return rows


def plot_similarity_summary(summary_rows: Sequence[dict[str, object]], figure_dir: Path) -> Path:
    df = pd.DataFrame(summary_rows)
    metrics = [
        "baseline_all_mode_pairs_cross_dataset",
        "baseline_same_sorted_position_cross_dataset",
        "efun_theta_within",
        "efun_rg_within",
        "deconv_rg_within",
        "efun_theta_vs_efun_rg",
        "efun_theta_vs_deconv_rg",
        "efun_rg_vs_deconv_rg",
    ]
    labels = {
        "baseline_all_mode_pairs_cross_dataset": "baseline all\nmode pairs",
        "baseline_same_sorted_position_cross_dataset": "baseline same\nsorted pos",
        "efun_theta_within": "efun theta\ncross-dataset",
        "efun_rg_within": "efun RG\ncross-dataset",
        "deconv_rg_within": "deconv RG\ncross-dataset",
        "efun_theta_vs_efun_rg": "same dataset\nefun theta vs efun RG",
        "efun_theta_vs_deconv_rg": "same dataset\nefun theta vs deconv RG",
        "efun_rg_vs_deconv_rg": "same dataset\nefun RG vs deconv RG",
    }
    fig, axes = plt.subplots(1, 3, figsize=(20, 5.5), sharey=True, constrained_layout=True)
    for ax, top_n in zip(axes, sorted(df["top_n"].unique())):
        sub = df[df["top_n"] == top_n]
        vals = []
        for metric in metrics:
            m = sub[sub["metric"] == metric]
            if metric.startswith("baseline"):
                val = m["median_corr"].astype(float).mean() if not m.empty else np.nan
            else:
                # Plot P8 and P10 separately by putting two bars side by side.
                val = np.nan
            vals.append(val)
        x = np.arange(len(metrics))
        base_vals = [vals[i] if metrics[i].startswith("baseline") else np.nan for i in range(len(metrics))]
        ax.bar(x[:2], base_vals[:2], width=0.55, color="#a6a6a6", label="baseline")
        for offset, pipeline, color in [(-0.18, "P8", "#4c78a8"), (0.18, "P10", "#f58518")]:
            y = []
            xpos = []
            for i, metric in enumerate(metrics[2:], start=2):
                m = sub[(sub["metric"] == metric) & (sub["pipeline"] == pipeline)]
                if m.empty:
                    continue
                y.append(float(m.iloc[0]["median_corr"]))
                xpos.append(i + offset)
            ax.bar(xpos, y, width=0.34, color=color, label=pipeline)
        ax.axhline(0, color="black", lw=0.8)
        ax.set_title(f"top{top_n}", weight="bold")
        ax.set_ylim(-0.25, 1.05)
        ax.set_xticks(x)
        ax.set_xticklabels([labels[m] for m in metrics], rotation=45, ha="right", fontsize=8)
        ax.grid(axis="y", alpha=0.2)
        if top_n == sorted(df["top_n"].unique())[0]:
            ax.set_ylabel("median ROI profile correlation")
        ax.legend(fontsize=8, loc="upper left")
    path = figure_dir / "03_target_roi_similarity_summary_by_topN.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main() -> None:
    args = parse_args()
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    output_dir = args.roi_test_dir / "target_footprint_tests"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_summary_rows: list[dict[str, object]] = []
    figure_paths: list[Path] = []
    for top_n in args.top_ns:
        profile_path = args.roi_test_dir / f"roi_mean_rg_no_theta_selected_roi_profiles_long__top{top_n}.csv"
        roi_order, vectors = load_profile_vectors(profile_path, args.datasets)
        all_summary_rows.extend(extract_similarity_rows(args.roi_test_dir, top_n))
        if top_n == args.main_top_n:
            for pipeline in ("P8", "P10"):
                figure_paths.append(plot_target_heatmaps(roi_order, vectors, args.datasets, args.figure_dir, pipeline, top_n))
                figure_paths.append(plot_contrast_heatmaps(roi_order, vectors, args.datasets, args.figure_dir, pipeline, top_n))
    write_csv(output_dir / "target_roi_similarity_summary_by_topN.csv", all_summary_rows)
    figure_paths.append(plot_similarity_summary(all_summary_rows, args.figure_dir))

    with (output_dir / "README_target_footprint_tests.md").open("w", encoding="utf-8") as handle:
        handle.write("# Target ROI footprint tests\n\n")
        handle.write("Targets:\n")
        for target_id, feature, label, title in TARGETS:
            handle.write(f"- `{target_id}` = `{feature}` x `{label}` ({title})\n")
        handle.write("\nMain figures:\n")
        for path in figure_paths:
            handle.write(f"- `{path}`\n")
        handle.write("\nMain CSV:\n")
        handle.write(f"- `{output_dir / 'target_roi_similarity_summary_by_topN.csv'}`\n")
    print(f"Wrote target footprint tests: {output_dir}")
    print(f"Figures: {args.figure_dir}")


if __name__ == "__main__":
    main()
