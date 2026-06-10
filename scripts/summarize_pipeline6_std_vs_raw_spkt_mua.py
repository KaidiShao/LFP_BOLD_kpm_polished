#!/usr/bin/env python3
"""Compare raw vs standardized BLP Pipeline 6 SPKT/MUA residual xcorr outputs."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean, median


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_OUTPUT_DIR = Path("results") / "pipeline6_std_vs_raw_spkt_mua_20260525"
DATASETS = ("e10gb1", "e10gh1")
MODALITY_STAGES = {
    "spkt": "pipeline6_spkt_residual_cross_correlation",
    "mua": "pipeline6_mua_residual_cross_correlation",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS))
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def as_int(value: object) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def observable_from_name(name: str) -> str:
    low = name.lower()
    if "complex_split" in low:
        return "complex_split"
    if re.search(r"(^|_)abs(_|$)", low) or "projected_vlambda_abs" in low:
        return "abs"
    return "unknown"


def variant_from_name(name: str) -> str:
    low = name.lower()
    if "stdcomplexpair" in low or "_stdt_" in low:
        return "standardized"
    if "mainline_torchlike" in low or "torchlike_longterm" in low:
        return "raw"
    return "other"


def infer_source_from_csv(dataset: str, csv_path: Path) -> tuple[str, str]:
    observable = observable_from_name(csv_path.name)
    variant = variant_from_name(csv_path.name)
    low = csv_path.name.lower()

    # Long standardized complex_split run names can exceed Windows path limits
    # in the P6 saver, which intentionally falls back to a compact
    # <dataset>_<modality>... tag. In this probe those compact tags are the
    # standardized complex_split outputs.
    compact_prefix = f"{dataset.lower()}_"
    if observable == "unknown" and low.startswith(compact_prefix):
        observable = "complex_split"
        variant = "standardized"
    return observable, variant


def feature_name(modality: str, row: dict[str, str]) -> str:
    if modality == "spkt":
        return row.get("residual_feature", "")
    return row.get("pairing_label", "")


def discover_session_top(processed_root: Path, datasets: list[str]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for dataset in datasets:
        for modality, stage in MODALITY_STAGES.items():
            stage_root = processed_root / dataset / stage
            for csv_path in sorted(stage_root.glob("*_session_xcorr.csv")):
                observable, variant = infer_source_from_csv(dataset, csv_path)
                if observable == "unknown" or variant == "other":
                    continue

                grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
                for row in read_csv(csv_path):
                    session_id = row.get("session_id", "")
                    if not session_id or session_id.lower() == "nan":
                        continue
                    abs_corr = as_float(row.get("abs_corr"))
                    if not math.isfinite(abs_corr):
                        continue
                    grouped[(session_id, feature_name(modality, row))].append(row)

                for (session_id, feature), rows in grouped.items():
                    top = max(rows, key=lambda r: as_float(r.get("abs_corr")))
                    out.append(
                        {
                            "dataset": dataset,
                            "modality": modality,
                            "observable": observable,
                            "variant": variant,
                            "feature": feature,
                            "session_id": session_id,
                            "source_file": str(csv_path),
                            "channel_index": as_int(top.get("channel_index")),
                            "channel_site": top.get("channel_site", ""),
                            "channel_label": top.get("channel_label", ""),
                            "mode_rank": as_int(top.get("mode_rank")),
                            "corr": as_float(top.get("corr")),
                            "abs_corr": as_float(top.get("abs_corr")),
                        }
                    )
    return out


def discover_pooled_top(processed_root: Path, datasets: list[str]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for dataset in datasets:
        for modality, stage in MODALITY_STAGES.items():
            stage_root = processed_root / dataset / stage
            for csv_path in sorted(stage_root.glob("*_pooled_top_xcorr.csv")):
                observable, variant = infer_source_from_csv(dataset, csv_path)
                if observable == "unknown" or variant == "other":
                    continue
                rows = [
                    row for row in read_csv(csv_path)
                    if math.isfinite(as_float(row.get("abs_corr")))
                ]
                if not rows:
                    continue
                top = max(rows, key=lambda r: as_float(r.get("abs_corr")))
                out.append(
                    {
                        "dataset": dataset,
                        "modality": modality,
                        "observable": observable,
                        "variant": variant,
                        "feature": feature_name(modality, top),
                        "source_file": str(csv_path),
                        "channel_index": as_int(top.get("channel_index")),
                        "channel_site": top.get("channel_site", ""),
                        "channel_label": top.get("channel_label", ""),
                        "mode_rank": as_int(top.get("mode_rank")),
                        "corr": as_float(top.get("corr")),
                        "abs_corr": as_float(top.get("abs_corr")),
                    }
                )
    return out


def attach_pooled_top_summary(
    summary_rows: list[dict[str, object]],
    pooled_rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    pooled_lookup: dict[tuple[str, str, str, str], dict[str, object]] = {}
    for row in pooled_rows:
        key = (
            str(row["dataset"]),
            str(row["modality"]),
            str(row["observable"]),
            str(row["variant"]),
        )
        current = pooled_lookup.get(key)
        if current is None or float(row["abs_corr"]) > float(current["abs_corr"]):
            pooled_lookup[key] = row

    for row in summary_rows:
        key = (
            str(row["dataset"]),
            str(row["modality"]),
            str(row["observable"]),
            str(row["variant"]),
        )
        pooled = pooled_lookup.get(key)
        if pooled is None:
            row["pooled_top_abs_corr"] = math.nan
            row["pooled_top_channel"] = ""
            row["pooled_top_site"] = ""
            row["pooled_top_feature"] = ""
            row["pooled_top_mode_rank"] = ""
            row["pooled_top_source_file"] = ""
            continue
        row["pooled_top_abs_corr"] = float(pooled["abs_corr"])
        row["pooled_top_channel"] = pooled["channel_label"]
        row["pooled_top_site"] = pooled["channel_site"]
        row["pooled_top_feature"] = pooled["feature"]
        row["pooled_top_mode_rank"] = pooled["mode_rank"]
        row["pooled_top_source_file"] = pooled["source_file"]
    return summary_rows


def summarize(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[
            (
                str(row["dataset"]),
                str(row["modality"]),
                str(row["observable"]),
                str(row["variant"]),
            )
        ].append(row)

    out: list[dict[str, object]] = []
    for (dataset, modality, observable, variant), items in sorted(groups.items()):
        vals = [float(r["abs_corr"]) for r in items if math.isfinite(float(r["abs_corr"]))]
        chan_counts = Counter(str(r["channel_label"]) for r in items)
        site_counts = Counter(str(r["channel_site"]) for r in items)
        feature_counts = Counter(str(r["feature"]) for r in items)
        out.append(
            {
                "dataset": dataset,
                "modality": modality,
                "observable": observable,
                "variant": variant,
                "n_session_feature_tops": len(items),
                "mean_session_top_abs_corr": mean(vals) if vals else math.nan,
                "median_session_top_abs_corr": median(vals) if vals else math.nan,
                "max_session_top_abs_corr": max(vals) if vals else math.nan,
                "top_channels": ";".join(f"{k}:{v}" for k, v in chan_counts.most_common(5)),
                "top_sites": ";".join(f"{k}:{v}" for k, v in site_counts.most_common(5)),
                "feature_counts": ";".join(f"{k}:{v}" for k, v in feature_counts.most_common()),
            }
        )
    return out


def compare_std_raw(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    lookup = {
        (
            str(r["dataset"]),
            str(r["modality"]),
            str(r["observable"]),
            str(r["variant"]),
        ): r
        for r in summary_rows
    }
    out: list[dict[str, object]] = []
    keys = sorted({(k[0], k[1], k[2]) for k in lookup})
    for dataset, modality, observable in keys:
        raw = lookup.get((dataset, modality, observable, "raw"))
        std = lookup.get((dataset, modality, observable, "standardized"))
        if raw is None or std is None:
            out.append(
                {
                    "dataset": dataset,
                    "modality": modality,
                    "observable": observable,
                    "status": "missing_raw_or_standardized",
                    "raw_present": raw is not None,
                    "standardized_present": std is not None,
                }
            )
            continue
        raw_mean = float(raw["mean_session_top_abs_corr"])
        std_mean = float(std["mean_session_top_abs_corr"])
        raw_median = float(raw["median_session_top_abs_corr"])
        std_median = float(std["median_session_top_abs_corr"])
        raw_max = float(raw["max_session_top_abs_corr"])
        std_max = float(std["max_session_top_abs_corr"])
        raw_pooled = float(raw["pooled_top_abs_corr"])
        std_pooled = float(std["pooled_top_abs_corr"])
        out.append(
            {
                "dataset": dataset,
                "modality": modality,
                "observable": observable,
                "status": "ok",
                "raw_mean_session_top_abs_corr": raw_mean,
                "std_mean_session_top_abs_corr": std_mean,
                "std_minus_raw_mean": std_mean - raw_mean,
                "std_over_raw_mean": std_mean / raw_mean if raw_mean else math.nan,
                "raw_median_session_top_abs_corr": raw_median,
                "std_median_session_top_abs_corr": std_median,
                "std_minus_raw_median": std_median - raw_median,
                "std_over_raw_median": std_median / raw_median if raw_median else math.nan,
                "raw_max_session_top_abs_corr": raw_max,
                "std_max_session_top_abs_corr": std_max,
                "std_minus_raw_max": std_max - raw_max,
                "std_over_raw_max": std_max / raw_max if raw_max else math.nan,
                "raw_pooled_top_abs_corr": raw_pooled,
                "std_pooled_top_abs_corr": std_pooled,
                "std_minus_raw_pooled_top": std_pooled - raw_pooled,
                "std_over_raw_pooled_top": std_pooled / raw_pooled if raw_pooled else math.nan,
                "std_top_channels": std["top_channels"],
                "raw_top_channels": raw["top_channels"],
                "std_top_sites": std["top_sites"],
                "raw_top_sites": raw["top_sites"],
                "raw_pooled_top_channel": raw["pooled_top_channel"],
                "std_pooled_top_channel": std["pooled_top_channel"],
                "raw_pooled_top_feature": raw["pooled_top_feature"],
                "std_pooled_top_feature": std["pooled_top_feature"],
            }
        )
    return out


METRICS = [
    ("mean_session_top_abs_corr", "Mean session-top |corr|"),
    ("median_session_top_abs_corr", "Median session-top |corr|"),
    ("max_session_top_abs_corr", "Max session-top |corr|"),
    ("pooled_top_abs_corr", "Pooled top |corr|"),
]


def plot_summary(summary_rows: list[dict[str, object]], out_file: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        plot_summary_svg(summary_rows, out_file.with_suffix(".svg"))
        return

    groups = defaultdict(dict)
    for row in summary_rows:
        key = (row["dataset"], row["modality"])
        label = f"{row['observable']}\n{row['variant'].replace('standardized', 'std')}"
        groups[key][label] = float(row["mean_session_top_abs_corr"])

    labels = ["abs\nraw", "abs\nstd", "complex_split\nraw", "complex_split\nstd"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharey=False)
    axes = axes.ravel()
    for ax, key in zip(axes, sorted(groups)):
        vals = [groups[key].get(label, math.nan) for label in labels]
        colors = ["#7a7a7a", "#2c7fb8", "#9e9e9e", "#41ab5d"]
        ax.bar(range(len(labels)), vals, color=colors)
        ax.set_title(f"{key[0]} {key[1]}")
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=0, fontsize=8)
        ax.set_ylabel("mean session-top |corr|")
        ax.grid(axis="y", alpha=0.25)
    for ax in axes[len(groups):]:
        ax.axis("off")
    fig.suptitle("P6 SPKT/MUA residual xcorr: raw vs standardized")
    fig.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def plot_metric_grid(summary_rows: list[dict[str, object]], out_file: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ModuleNotFoundError:
        return

    groups = defaultdict(dict)
    for row in summary_rows:
        key = (row["dataset"], row["modality"])
        label = f"{row['observable']}\n{row['variant'].replace('standardized', 'std')}"
        groups[key][label] = row

    labels = ["abs\nraw", "abs\nstd", "complex_split\nraw", "complex_split\nstd"]
    colors = ["#7a7a7a", "#2c7fb8", "#9e9e9e", "#41ab5d"]
    keys = sorted(groups)
    fig, axes = plt.subplots(len(METRICS), len(keys), figsize=(4.1 * len(keys), 2.55 * len(METRICS)))
    if len(METRICS) == 1:
        axes = [axes]
    for i_metric, (metric, metric_title) in enumerate(METRICS):
        vals_all = []
        for key in keys:
            for label in labels:
                row = groups[key].get(label)
                if row is not None:
                    val = float(row.get(metric, math.nan))
                    if math.isfinite(val):
                        vals_all.append(val)
        y_max = max(vals_all or [1.0]) * 1.12
        for i_key, key in enumerate(keys):
            ax = axes[i_metric][i_key] if len(keys) > 1 else axes[i_metric]
            vals = []
            for label in labels:
                row = groups[key].get(label)
                vals.append(float(row.get(metric, math.nan)) if row is not None else math.nan)
            ax.bar(range(len(labels)), vals, color=colors)
            ax.set_ylim(0, y_max)
            ax.set_title(f"{key[0]} {key[1]}" if i_metric == 0 else "")
            ax.set_ylabel(metric_title)
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, fontsize=8)
            ax.grid(axis="y", alpha=0.25)
            for j, val in enumerate(vals):
                if math.isfinite(val):
                    ax.text(j, val + y_max * 0.015, f"{val:.3g}", ha="center", va="bottom", fontsize=7)
    fig.suptitle("P6 SPKT/MUA residual xcorr metrics: raw vs standardized")
    fig.tight_layout()
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def plot_summary_svg(summary_rows: list[dict[str, object]], out_file: Path) -> None:
    groups = defaultdict(dict)
    for row in summary_rows:
        key = (row["dataset"], row["modality"])
        label = f"{row['observable']} {row['variant'].replace('standardized', 'std')}"
        groups[key][label] = float(row["mean_session_top_abs_corr"])

    labels = ["abs raw", "abs std", "complex_split raw", "complex_split std"]
    keys = sorted(groups)
    width = 980
    panel_h = 180
    height = 70 + panel_h * max(1, len(keys))
    max_val = max(
        [v for vals in groups.values() for v in vals.values() if math.isfinite(v)] or [1.0]
    )
    colors = ["#7a7a7a", "#2c7fb8", "#9e9e9e", "#41ab5d"]

    def esc(text: object) -> str:
        return (
            str(text)
            .replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
        )

    body = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="24" y="32" font-family="Arial" font-size="20" font-weight="700">P6 SPKT/MUA residual xcorr: raw vs standardized</text>',
    ]
    for i, key in enumerate(keys):
        y0 = 60 + i * panel_h
        x0 = 190
        plot_w = 710
        plot_h = 105
        body.append(f'<text x="24" y="{y0 + 18}" font-family="Arial" font-size="15" font-weight="700">{esc(key[0])} {esc(key[1])}</text>')
        body.append(f'<line x1="{x0}" y1="{y0 + plot_h}" x2="{x0 + plot_w}" y2="{y0 + plot_h}" stroke="#777"/>')
        body.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y0 + plot_h}" stroke="#777"/>')
        for j, label in enumerate(labels):
            val = groups[key].get(label, math.nan)
            bar_h = 0 if not math.isfinite(val) else (val / max_val) * (plot_h - 10)
            bar_w = 90
            x = x0 + 35 + j * 160
            y = y0 + plot_h - bar_h
            body.append(f'<rect x="{x}" y="{y:.2f}" width="{bar_w}" height="{bar_h:.2f}" fill="{colors[j]}"/>')
            shown = "NA" if not math.isfinite(val) else f"{val:.3g}"
            body.append(f'<text x="{x + bar_w / 2}" y="{y - 5:.2f}" text-anchor="middle" font-family="Arial" font-size="11">{shown}</text>')
            body.append(f'<text x="{x + bar_w / 2}" y="{y0 + plot_h + 34}" text-anchor="middle" font-family="Arial" font-size="11">{esc(label)}</text>')
    body.append("</svg>")
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text("\n".join(body) + "\n", encoding="utf-8")


def write_markdown(path: Path, summary_rows: list[dict[str, object]], delta_rows: list[dict[str, object]]) -> None:
    lines = [
        "# P6 Standardized vs Raw SPKT/MUA Summary",
        "",
        "Metric: for each session and residual feature, take the channel with maximum `abs_corr`; then summarize those session-feature top correlations.",
        "",
        "## Mean Session-Top |corr|",
        "",
        "| dataset | modality | observable | raw | standardized | std/raw | std-raw |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in delta_rows:
        if row.get("status") != "ok":
            lines.append(
                f"| {row['dataset']} | {row['modality']} | {row['observable']} | missing | missing |  |  |"
            )
            continue
        lines.append(
            "| {dataset} | {modality} | {observable} | {raw:.4g} | {std:.4g} | {ratio:.3g} | {delta:.4g} |".format(
                dataset=row["dataset"],
                modality=row["modality"],
                observable=row["observable"],
                raw=float(row["raw_mean_session_top_abs_corr"]),
                std=float(row["std_mean_session_top_abs_corr"]),
                ratio=float(row["std_over_raw_mean"]),
                delta=float(row["std_minus_raw_mean"]),
            )
        )
    lines.extend(["", "## Median / Max / Pooled Top |corr|", ""])
    lines.append("| dataset | modality | observable | metric | raw | standardized | std/raw | std-raw |")
    lines.append("|---|---|---|---|---:|---:|---:|---:|")
    metric_specs = [
        ("median", "raw_median_session_top_abs_corr", "std_median_session_top_abs_corr", "std_over_raw_median", "std_minus_raw_median"),
        ("max session-top", "raw_max_session_top_abs_corr", "std_max_session_top_abs_corr", "std_over_raw_max", "std_minus_raw_max"),
        ("pooled top", "raw_pooled_top_abs_corr", "std_pooled_top_abs_corr", "std_over_raw_pooled_top", "std_minus_raw_pooled_top"),
    ]
    for row in delta_rows:
        if row.get("status") != "ok":
            continue
        for label, raw_key, std_key, ratio_key, delta_key in metric_specs:
            lines.append(
                "| {dataset} | {modality} | {observable} | {label} | {raw:.4g} | {std:.4g} | {ratio:.3g} | {delta:.4g} |".format(
                    dataset=row["dataset"],
                    modality=row["modality"],
                    observable=row["observable"],
                    label=label,
                    raw=float(row[raw_key]),
                    std=float(row[std_key]),
                    ratio=float(row[ratio_key]),
                    delta=float(row[delta_key]),
                )
            )
    lines.extend(["", "## Recurrent Top Channels", ""])
    lines.append("| dataset | modality | observable | raw top channels | standardized top channels |")
    lines.append("|---|---|---|---|---|")
    for row in delta_rows:
        if row.get("status") != "ok":
            continue
        lines.append(
            f"| {row['dataset']} | {row['modality']} | {row['observable']} | "
            f"{row['raw_top_channels']} | {row['std_top_channels']} |"
        )
    lines.extend(["", "## Pooled Top Channel", ""])
    lines.append("| dataset | modality | observable | raw pooled top | standardized pooled top |")
    lines.append("|---|---|---|---|---|")
    for row in delta_rows:
        if row.get("status") != "ok":
            continue
        raw_desc = f"{row['raw_pooled_top_channel']} / {row['raw_pooled_top_feature']}"
        std_desc = f"{row['std_pooled_top_channel']} / {row['std_pooled_top_feature']}"
        lines.append(
            f"| {row['dataset']} | {row['modality']} | {row['observable']} | {raw_desc} | {std_desc} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    rows = discover_session_top(args.processed_root, list(args.datasets))
    pooled_rows = discover_pooled_top(args.processed_root, list(args.datasets))
    summary_rows = summarize(rows)
    summary_rows = attach_pooled_top_summary(summary_rows, pooled_rows)
    delta_rows = compare_std_raw(summary_rows)

    write_csv(output_dir / "session_feature_top_channels.csv", rows)
    write_csv(output_dir / "pooled_top_channels.csv", pooled_rows)
    write_csv(output_dir / "summary_by_variant.csv", summary_rows)
    write_csv(output_dir / "std_vs_raw_delta.csv", delta_rows)
    plot_summary(summary_rows, output_dir / "figures" / "mean_session_top_abs_corr.png")
    plot_metric_grid(summary_rows, output_dir / "figures" / "metrics_raw_vs_standardized.png")
    write_markdown(output_dir / "README.md", summary_rows, delta_rows)

    print(f"session-feature rows: {len(rows)}")
    print(f"summary rows: {len(summary_rows)}")
    print(f"wrote: {output_dir}")


if __name__ == "__main__":
    main()
