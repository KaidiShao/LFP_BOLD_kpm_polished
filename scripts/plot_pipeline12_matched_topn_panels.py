#!/usr/bin/env python3
"""Matched Pipeline12 topN panels for side-by-side dataset comparison."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_P12_ROOT = REPO_ROOT / "results" / "pipeline12_standardized_csplit_consistency"
PIPELINE_FAMILIES = (
    ("P8", "efun"),
    ("P8", "deconv_efun"),
    ("P10", "efun"),
    ("P10", "deconv_efun"),
)
DENSITY_CLASSES = ("event_density", "raw_efun_density", "dimred_efun_density")
DENSITY_COLORS = {
    "event_density": "#9aa0a6",
    "raw_efun_density": "#2f80ed",
    "dimred_efun_density": "#d95f02",
}
PROCESS_LABELS = ("theta_selective", "ripple_gamma_no_pure_theta", "other")
PROCESS_LABEL_COLORS = {
    "theta_selective": "#31a354",
    "ripple_gamma_no_pure_theta": "#3182bd",
    "other": "#8f8f8f",
}
DATASET_COLORS = {
    "e10gb1": "#1b9e77",
    "e10gh1": "#d95f02",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="+", default=["e10gb1", "e10gh1"])
    parser.add_argument("--p12-root", type=Path, default=DEFAULT_P12_ROOT)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--top-n-values", nargs="+", type=int, default=[1, 3, 5, 10, 20])
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value: object) -> float:
    try:
        out = float(str(value))
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def as_int(value: object) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return -1


def mean(values: Iterable[float]) -> float:
    clean = [v for v in values if math.isfinite(v)]
    return sum(clean) / len(clean) if clean else math.nan


def median(values: Iterable[float]) -> float:
    clean = sorted(v for v in values if math.isfinite(v))
    return clean[len(clean) // 2] if clean else math.nan


def rank_topn(rows: Sequence[dict[str, str]], top_n: int) -> list[dict[str, str]]:
    return [row for row in rows if 1 <= as_int(row.get("rank_within_group")) <= top_n]


def load_tables(root: Path, datasets: Sequence[str]) -> dict[str, dict[str, list[dict[str, str]]]]:
    tables = {}
    for dataset in datasets:
        ds_root = root / dataset
        tables[dataset] = {
            "membership": read_csv(ds_root / "p12_density_class_topN_membership.csv"),
            "dimred": read_csv(ds_root / "p12_dimred_top_rows_with_labels.csv"),
            "raw": read_csv(ds_root / "p12_raw_top_rows_with_timescale.csv"),
            "best": read_csv(ds_root / "p12_density_class_best_score.csv"),
        }
    return tables


def style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.18)


def density_fraction(rows: Sequence[dict[str, str]], pipeline: str, family: str, top_n: int, density_class: str) -> float:
    vals = [
        as_float(row.get("fraction"))
        for row in rows
        if row.get("pipeline") == pipeline
        and row.get("feature_family") == family
        and as_int(row.get("top_n")) == top_n
        and row.get("density_class") == density_class
    ]
    return mean(vals)


def plot_density_class_matched(tables, datasets: Sequence[str], top_ns: Sequence[int], out: Path) -> None:
    fig, axes = plt.subplots(len(datasets), len(PIPELINE_FAMILIES), figsize=(16, 3.5 * len(datasets)), sharey=True)
    if len(datasets) == 1:
        axes = [axes]
    for row_axes, dataset in zip(axes, datasets):
        rows = tables[dataset]["membership"]
        for ax, (pipeline, family) in zip(row_axes, PIPELINE_FAMILIES):
            bottom = [0.0] * len(top_ns)
            for density_class in DENSITY_CLASSES:
                vals = [density_fraction(rows, pipeline, family, top_n, density_class) for top_n in top_ns]
                vals = [0.0 if not math.isfinite(v) else v for v in vals]
                ax.bar([str(n) for n in top_ns], vals, bottom=bottom, color=DENSITY_COLORS[density_class], label=density_class)
                bottom = [a + b for a, b in zip(bottom, vals)]
            ax.set_title(f"{dataset} | {pipeline} {family}")
            ax.set_ylim(0, 1.02)
            ax.set_xlabel("topN")
            style_axes(ax)
        row_axes[0].set_ylabel("mean fraction")
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, loc="upper center", ncol=len(DENSITY_CLASSES), bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Matched density-class fraction by topN", y=1.06, fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def collapse_process_label(label: str) -> str:
    text = str(label or "").lower()
    if text in {
        "theta_strict",
        "theta_selective_similar",
        "theta_selective_unequal",
        "theta_relaxed",
    }:
        return "theta_selective"
    if text in {
        "ripple_gamma_no_pure_theta",
        "ripple_selective_similar",
        "ripple_selective_unequal",
        "gamma_selective",
    }:
        return "ripple_gamma_no_pure_theta"
    return "other"


def plot_dimred_label_matched(tables, datasets: Sequence[str], top_ns: Sequence[int], out: Path) -> None:
    fig, axes = plt.subplots(len(datasets), len(PIPELINE_FAMILIES), figsize=(17, 3.9 * len(datasets)), sharey=True)
    if len(datasets) == 1:
        axes = [axes]
    for row_axes, dataset in zip(axes, datasets):
        rows_all = tables[dataset]["dimred"]
        for ax, (pipeline, family) in zip(row_axes, PIPELINE_FAMILIES):
            rows = [row for row in rows_all if row.get("pipeline") == pipeline and row.get("feature_family") == family]
            bottom = [0.0] * len(top_ns)
            for label in PROCESS_LABELS:
                vals = []
                for top_n in top_ns:
                    selected = rank_topn(rows, top_n)
                    denom = max(1, len(selected))
                    vals.append(sum(1 for row in selected if collapse_process_label(row.get("two_subprocess_label")) == label) / denom)
                if any(v > 0 for v in vals):
                    ax.bar([str(n) for n in top_ns], vals, bottom=bottom, color=PROCESS_LABEL_COLORS[label], label=label)
                    bottom = [a + b for a, b in zip(bottom, vals)]
            ax.set_title(f"{dataset} | {pipeline} {family}")
            ax.set_ylim(0, 1.02)
            ax.set_xlabel("topN")
            style_axes(ax)
        row_axes[0].set_ylabel("fraction")
    handles, labels_out = axes[0][0].get_legend_handles_labels()
    # Collect legend entries from all axes because some labels appear only in one dataset/context.
    seen_labels = set(labels_out)
    for ax in [ax for row in axes for ax in row]:
        h, l = ax.get_legend_handles_labels()
        for handle, label in zip(h, l):
            if label not in seen_labels:
                handles.append(handle)
                labels_out.append(label)
                seen_labels.add(label)
    fig.legend(handles, labels_out, frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
    fig.suptitle("Matched dimred component process-label fraction by topN", y=1.04, fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def raw_timescale_values(rows: Sequence[dict[str, str]], pipeline: str, family: str, top_n: int) -> list[float]:
    selected = [
        row for row in rows
        if row.get("pipeline") == pipeline and row.get("feature_family") == family
    ]
    selected = rank_topn(selected, top_n)
    vals = [as_float(row.get("timescale_sec_preferred")) for row in selected]
    return [v for v in vals if v > 0 and v < 1e9]


def plot_raw_timescale_lines(tables, datasets: Sequence[str], top_ns: Sequence[int], out_median: Path, out_slow: Path) -> None:
    fig_m, axes_m = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    fig_s, axes_s = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    all_medians = []
    for ax_m, ax_s, (pipeline, family) in zip(axes_m.ravel(), axes_s.ravel(), PIPELINE_FAMILIES):
        for dataset in datasets:
            rows = tables[dataset]["raw"]
            medians = []
            slow_fracs = []
            for top_n in top_ns:
                vals = raw_timescale_values(rows, pipeline, family, top_n)
                log_vals = [math.log10(v) for v in vals]
                medians.append(median(log_vals))
                slow_fracs.append(sum(1 for v in vals if v >= 1.0) / len(vals) if vals else math.nan)
            all_medians.extend(v for v in medians if math.isfinite(v))
            color = DATASET_COLORS.get(dataset, None)
            ax_m.plot(top_ns, medians, marker="o", linewidth=2, label=dataset, color=color)
            ax_s.plot(top_ns, slow_fracs, marker="o", linewidth=2, label=dataset, color=color)
        ax_m.set_title(f"{pipeline} {family}")
        ax_s.set_title(f"{pipeline} {family}")
        ax_m.set_ylabel("median log10(sec)")
        ax_s.set_ylabel("fraction >= 1 sec")
        ax_s.set_ylim(-0.02, 1.02)
        style_axes(ax_m)
        style_axes(ax_s)
    if all_medians:
        ymin = min(all_medians) - 0.2
        ymax = max(all_medians) + 0.2
        for ax in axes_m.ravel():
            ax.set_ylim(ymin, ymax)
    for ax in axes_m[-1, :]:
        ax.set_xlabel("topN")
    for ax in axes_s[-1, :]:
        ax.set_xlabel("topN")
    axes_m[0, 1].legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    axes_s[0, 1].legend(frameon=False, fontsize=9, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig_m.suptitle("Matched raw BLP efun median timescale by topN", fontsize=15, weight="bold")
    fig_s.suptitle("Matched raw BLP efun slow-mode fraction by topN", fontsize=15, weight="bold")
    fig_m.tight_layout()
    fig_s.tight_layout()
    fig_m.savefig(out_median, dpi=180, bbox_inches="tight")
    fig_s.savefig(out_slow, dpi=180, bbox_inches="tight")
    plt.close(fig_m)
    plt.close(fig_s)


def plot_best_score_matched(tables, datasets: Sequence[str], out: Path) -> None:
    fig, axes = plt.subplots(len(datasets), len(PIPELINE_FAMILIES), figsize=(16, 3.4 * len(datasets)), sharey=True)
    if len(datasets) == 1:
        axes = [axes]
    all_values = []
    for dataset in datasets:
        for row in tables[dataset]["best"]:
            value = as_float(row.get("peak_abs_corr"))
            if math.isfinite(value):
                all_values.append(value)
    ymax = max(all_values) * 1.12 if all_values else 1.0
    for row_axes, dataset in zip(axes, datasets):
        rows_all = tables[dataset]["best"]
        for ax, (pipeline, family) in zip(row_axes, PIPELINE_FAMILIES):
            vals = []
            for density_class in DENSITY_CLASSES:
                these = [
                    as_float(row.get("peak_abs_corr"))
                    for row in rows_all
                    if row.get("pipeline") == pipeline
                    and row.get("feature_family") == family
                    and row.get("density_class") == density_class
                ]
                vals.append(median(these))
            ax.bar(DENSITY_CLASSES, vals, color=[DENSITY_COLORS[c] for c in DENSITY_CLASSES])
            ax.set_title(f"{dataset} | {pipeline} {family}")
            ax.set_ylim(0, ymax)
            ax.tick_params(axis="x", labelrotation=28)
            style_axes(ax)
        row_axes[0].set_ylabel("median best |r|")
    fig.suptitle("Matched best score by density class", y=1.04, fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_summary(out_dir: Path, datasets: Sequence[str], top_ns: Sequence[int]) -> None:
    lines = [
        "# Matched Pipeline12 topN panels",
        "",
        "- datasets: `" + "`, `".join(datasets) + "`",
        "- topN values: `" + "`, `".join(str(n) for n in top_ns) + "`",
        "- Each figure uses the same row/column layout: rows are datasets, columns are `P8 efun`, `P8 deconv_efun`, `P10 efun`, `P10 deconv_efun` where applicable.",
        "- Raw-timescale plots use the same four contexts as separate subplots, with datasets overlaid as lines.",
        "",
        "## Figures",
        "",
    ]
    for path in sorted(out_dir.glob("*.png")):
        lines.append(f"- `{path.name}`")
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    out_dir = args.output_dir or (args.p12_root / ("matched_topn_panels_" + "_".join(args.datasets)))
    out_dir.mkdir(parents=True, exist_ok=True)
    tables = load_tables(args.p12_root, args.datasets)
    plot_density_class_matched(tables, args.datasets, args.top_n_values, out_dir / "01_density_class_fraction_by_topN_matched.png")
    plot_dimred_label_matched(tables, args.datasets, args.top_n_values, out_dir / "02_dimred_label_fraction_by_topN_matched.png")
    plot_dimred_label_matched(tables, args.datasets, args.top_n_values, out_dir / "02_dimred_three_class_fraction_by_topN_matched.png")
    plot_raw_timescale_lines(
        tables,
        args.datasets,
        args.top_n_values,
        out_dir / "03_raw_timescale_median_by_topN_matched.png",
        out_dir / "04_raw_timescale_slow_fraction_by_topN_matched.png",
    )
    plot_best_score_matched(tables, args.datasets, out_dir / "05_best_score_by_density_class_matched.png")
    write_summary(out_dir, args.datasets, args.top_n_values)
    print(f"Matched Pipeline12 figures: {out_dir}")


if __name__ == "__main__":
    main()
