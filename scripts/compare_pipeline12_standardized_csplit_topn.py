#!/usr/bin/env python3
"""Compare Pipeline12 standardized-csplit results across datasets by topN."""

from __future__ import annotations

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_P12_ROOT = REPO_ROOT / "results" / "pipeline12_standardized_csplit_consistency"
DENSITY_CLASSES = ("event_density", "raw_efun_density", "dimred_efun_density")
FEATURE_FAMILIES = ("efun", "deconv_efun")
PIPELINES = ("P8", "P10")
PROCESS_LABELS = ("theta_selective", "ripple_gamma_no_pure_theta", "other")
DENSITY_COLORS = {
    "event_density": "#9aa0a6",
    "raw_efun_density": "#2f80ed",
    "dimred_efun_density": "#d95f02",
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


def as_int(value: object, default: int = -1) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def mean(values: Iterable[float]) -> float:
    clean = [v for v in values if math.isfinite(v)]
    if not clean:
        return math.nan
    return sum(clean) / len(clean)


def median(values: Iterable[float]) -> float:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return math.nan
    return clean[len(clean) // 2]


def rank_topn(rows: Sequence[dict[str, str]], top_n: int) -> list[dict[str, str]]:
    return [row for row in rows if 1 <= as_int(row.get("rank_within_group")) <= top_n]


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


def label_color(label: str) -> str:
    palette = {
        "theta_selective": "#31a354",
        "ripple_gamma_no_pure_theta": "#3182bd",
        "other": "#8f8f8f",
    }
    return palette.get(label, "#bbbbbb")


def load_dataset_tables(root: Path, datasets: Sequence[str]) -> dict[str, dict[str, list[dict[str, str]]]]:
    tables: dict[str, dict[str, list[dict[str, str]]]] = {}
    for dataset in datasets:
        ds_root = root / dataset
        tables[dataset] = {
            "membership": read_csv(ds_root / "p12_density_class_topN_membership.csv"),
            "dimred": read_csv(ds_root / "p12_dimred_top_rows_with_labels.csv"),
            "raw": read_csv(ds_root / "p12_raw_top_rows_with_timescale.csv"),
        }
    return tables


def plot_density_membership(tables, datasets: Sequence[str], top_ns: Sequence[int], out_dir: Path) -> None:
    for pipeline in PIPELINES:
        for family in FEATURE_FAMILIES:
            fig, axes = plt.subplots(len(datasets), 1, figsize=(9.5, max(3.2, 2.6 * len(datasets))), sharex=True, sharey=True)
            if len(datasets) == 1:
                axes = [axes]
            for ax, dataset in zip(axes, datasets):
                rows = [
                    row for row in tables[dataset]["membership"]
                    if row.get("pipeline") == pipeline and row.get("feature_family") == family
                ]
                bottom = [0.0] * len(top_ns)
                for density_class in DENSITY_CLASSES:
                    vals = []
                    for top_n in top_ns:
                        these = [
                            as_float(row.get("fraction"))
                            for row in rows
                            if as_int(row.get("top_n")) == top_n and row.get("density_class") == density_class
                        ]
                        vals.append(mean(these) if these else 0.0)
                    ax.bar([str(n) for n in top_ns], vals, bottom=bottom, label=density_class, color=DENSITY_COLORS[density_class])
                    bottom = [a + b for a, b in zip(bottom, vals)]
                ax.set_title(dataset)
                ax.set_ylim(0, 1.02)
                ax.set_ylabel("mean fraction")
                ax.grid(axis="y", alpha=0.22)
            axes[-1].set_xlabel("topN per observable context")
            axes[0].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
            fig.suptitle(f"{pipeline} {family}: density class fraction by topN")
            fig.tight_layout()
            fig.savefig(out_dir / f"density_class_fraction_by_topN__{pipeline.lower()}__{family}.png", dpi=180)
            plt.close(fig)


def plot_dimred_labels(tables, datasets: Sequence[str], top_ns: Sequence[int], out_dir: Path) -> None:
    for pipeline in PIPELINES:
        for family in FEATURE_FAMILIES:
            fig, axes = plt.subplots(len(datasets), 1, figsize=(11, max(3.2, 2.8 * len(datasets))), sharex=True, sharey=True)
            if len(datasets) == 1:
                axes = [axes]
            for ax, dataset in zip(axes, datasets):
                rows = [
                    row for row in tables[dataset]["dimred"]
                    if row.get("pipeline") == pipeline and row.get("feature_family") == family
                ]
                bottom = [0.0] * len(top_ns)
                for label in PROCESS_LABELS:
                    vals = []
                    for top_n in top_ns:
                        selected = rank_topn(rows, top_n)
                        denom = max(1, len(selected))
                        vals.append(sum(1 for row in selected if collapse_process_label(row.get("two_subprocess_label")) == label) / denom)
                    if any(v > 0 for v in vals):
                        ax.bar([str(n) for n in top_ns], vals, bottom=bottom, label=label, color=label_color(label))
                        bottom = [a + b for a, b in zip(bottom, vals)]
                ax.set_title(dataset)
                ax.set_ylim(0, 1.02)
                ax.set_ylabel("fraction")
                ax.grid(axis="y", alpha=0.22)
            axes[-1].set_xlabel("topN per observable context")
            axes[0].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0))
            fig.suptitle(f"{pipeline} {family}: dimred process-label fraction by topN")
            fig.tight_layout()
            fig.savefig(out_dir / f"dimred_label_fraction_by_topN__{pipeline.lower()}__{family}.png", dpi=180)
            plt.close(fig)


def plot_raw_timescale(tables, datasets: Sequence[str], top_ns: Sequence[int], out_dir: Path) -> None:
    fig_median, axes_median = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    fig_slow, axes_slow = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    contexts = [(pipeline, family) for pipeline in PIPELINES for family in FEATURE_FAMILIES]
    for ax_m, ax_s, (pipeline, family) in zip(axes_median.ravel(), axes_slow.ravel(), contexts):
        for dataset in datasets:
            rows = [
                row for row in tables[dataset]["raw"]
                if row.get("pipeline") == pipeline and row.get("feature_family") == family
            ]
            medians = []
            slow_fracs = []
            for top_n in top_ns:
                selected = rank_topn(rows, top_n)
                vals = [as_float(row.get("timescale_sec_preferred")) for row in selected]
                vals = [value for value in vals if value > 0 and value < 1e9]
                log_vals = [math.log10(value) for value in vals]
                medians.append(median(log_vals))
                slow_fracs.append(sum(1 for value in vals if value >= 1.0) / len(vals) if vals else math.nan)
            ax_m.plot(top_ns, medians, marker="o", label=dataset)
            ax_s.plot(top_ns, slow_fracs, marker="o", label=dataset)
        ax_m.set_title(f"{pipeline} {family}")
        ax_s.set_title(f"{pipeline} {family}")
        ax_m.grid(alpha=0.25)
        ax_s.grid(alpha=0.25)
        ax_m.set_ylabel("median log10(sec)")
        ax_s.set_ylabel("fraction >= 1 sec")
        ax_s.set_ylim(-0.02, 1.02)
    for ax in axes_median[-1, :]:
        ax.set_xlabel("topN per observable context")
    for ax in axes_slow[-1, :]:
        ax.set_xlabel("topN per observable context")
    axes_median[0, 1].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    axes_slow[0, 1].legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1.0))
    fig_median.suptitle("Raw BLP efun topN median timescale")
    fig_slow.suptitle("Raw BLP efun topN slow-mode fraction")
    fig_median.tight_layout()
    fig_slow.tight_layout()
    fig_median.savefig(out_dir / "raw_timescale_median_by_topN.png", dpi=180)
    fig_slow.savefig(out_dir / "raw_timescale_slow_fraction_by_topN.png", dpi=180)
    plt.close(fig_median)
    plt.close(fig_slow)


def write_summary(tables, datasets: Sequence[str], top_ns: Sequence[int], out_dir: Path) -> None:
    lines = ["# Pipeline12 standardized csplit topN comparison", ""]
    lines.append("- datasets: `" + "`, `".join(datasets) + "`")
    lines.append("- topN values: `" + "`, `".join(str(n) for n in top_ns) + "`")
    lines.extend(["", "## Outputs", ""])
    for name in sorted(path.name for path in out_dir.glob("*.png")):
        lines.append(f"- `{name}`")
    lines.extend(["", "## Notes", ""])
    lines.append("- Fractions are computed within each dataset and topN, using the current Pipeline12 CSV outputs.")
    lines.append("- Raw-timescale plots use rank-within-observable-context topN, not a global topN pooled across contexts.")
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    out_dir = args.output_dir or (args.p12_root / ("topn_comparison_" + "_".join(args.datasets)))
    out_dir.mkdir(parents=True, exist_ok=True)
    tables = load_dataset_tables(args.p12_root, args.datasets)
    plot_density_membership(tables, args.datasets, args.top_n_values, out_dir)
    plot_dimred_labels(tables, args.datasets, args.top_n_values, out_dir)
    plot_raw_timescale(tables, args.datasets, args.top_n_values, out_dir)
    write_summary(tables, args.datasets, args.top_n_values, out_dir)
    print(f"Pipeline12 topN comparison figures: {out_dir}")


if __name__ == "__main__":
    main()
