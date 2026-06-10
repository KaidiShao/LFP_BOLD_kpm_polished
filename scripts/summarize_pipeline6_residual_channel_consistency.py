#!/usr/bin/env python3
"""Summarize Pipeline 6 spike/MUA residual xcorr channel consistency.

This script reads the lightweight CSVs saved next to the large P6 residual
cross-correlation MAT files. It does not load the MAT files.

The main question is spatial: which spike/MUA channels repeatedly carry the
largest residual cross-correlation across sessions, conditions, and modalities?
Spatial information currently comes from cfg channel order and site labels:
channel_index, channel_label, and channel_site.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import mean, median
from typing import Iterable, Sequence


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_OUTPUT_DIR = Path("results") / "pipeline6_residual_channel_consistency_current"
DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "f12m01", "e10gw1")
MODALITY_STAGES = {
    "spkt": "pipeline6_spkt_residual_cross_correlation",
    "mua": "pipeline6_mua_residual_cross_correlation",
}
CONDITIONS = ("abs", "complex_split")


@dataclass(frozen=True)
class SessionTop:
    dataset: str
    modality: str
    condition: str
    source_file: str
    session_id: str
    session_note: str
    feature: str
    channel_index: int
    channel_site: str
    channel_label: str
    mode_rank: int
    corr: float
    abs_corr: float


@dataclass(frozen=True)
class PooledTop:
    dataset: str
    modality: str
    condition: str
    source_file: str
    feature: str
    channel_index: int
    channel_site: str
    channel_label: str
    mode_rank: int
    corr: float
    abs_corr: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--top-n-overlap", type=int, default=3)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def as_float(value: object, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def condition_from_name(name: str) -> str:
    low = name.lower()
    if "complex_split" in low:
        return "complex_split"
    if "projected_vlambda_abs" in low or re.search(r"(^|_)abs(_|$)", low):
        return "abs"
    return "unknown"


def feature_name(modality: str, row: dict[str, str]) -> str:
    if modality == "spkt":
        return row.get("residual_feature", "")
    return row.get("pairing_label", "")


def discover_rows(processed_root: Path, datasets: Sequence[str]) -> tuple[list[SessionTop], list[PooledTop], list[dict[str, object]]]:
    session_top: list[SessionTop] = []
    pooled_top: list[PooledTop] = []
    missing: list[dict[str, object]] = []

    for dataset in datasets:
        ds_root = processed_root / dataset
        for modality, stage in MODALITY_STAGES.items():
            stage_root = ds_root / stage
            if not stage_root.is_dir():
                missing.append(
                    {
                        "dataset": dataset,
                        "modality": modality,
                        "missing": "stage_dir",
                        "path": str(stage_root),
                    }
                )
                continue

            session_csvs = sorted(stage_root.glob("*_session_xcorr.csv"))
            pooled_csvs = sorted(stage_root.glob("*_pooled_top_xcorr.csv"))
            if not session_csvs:
                missing.append(
                    {
                        "dataset": dataset,
                        "modality": modality,
                        "missing": "session_xcorr_csv",
                        "path": str(stage_root),
                    }
                )
            if not pooled_csvs:
                missing.append(
                    {
                        "dataset": dataset,
                        "modality": modality,
                        "missing": "pooled_top_xcorr_csv",
                        "path": str(stage_root),
                    }
                )

            for csv_path in session_csvs:
                condition = condition_from_name(csv_path.name)
                groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
                for row in read_csv(csv_path):
                    session_id = row.get("session_id", "")
                    if not session_id or session_id.lower() == "nan":
                        continue
                    abs_corr = as_float(row.get("abs_corr"))
                    if not math.isfinite(abs_corr):
                        continue
                    groups[(session_id, feature_name(modality, row))].append(row)

                for (session_id, feat), rows in groups.items():
                    row = max(rows, key=lambda r: as_float(r.get("abs_corr")))
                    session_top.append(
                        SessionTop(
                            dataset=dataset,
                            modality=modality,
                            condition=condition,
                            source_file=str(csv_path),
                            session_id=session_id,
                            session_note=row.get("session_note", ""),
                            feature=feat,
                            channel_index=as_int(row.get("channel_index")),
                            channel_site=row.get("channel_site", ""),
                            channel_label=row.get("channel_label", ""),
                            mode_rank=as_int(row.get("mode_rank")),
                            corr=as_float(row.get("corr")),
                            abs_corr=as_float(row.get("abs_corr")),
                        )
                    )

            for csv_path in pooled_csvs:
                condition = condition_from_name(csv_path.name)
                for row in read_csv(csv_path):
                    pooled_top.append(
                        PooledTop(
                            dataset=dataset,
                            modality=modality,
                            condition=condition,
                            source_file=str(csv_path),
                            feature=feature_name(modality, row),
                            channel_index=as_int(row.get("channel_index")),
                            channel_site=row.get("channel_site", ""),
                            channel_label=row.get("channel_label", ""),
                            mode_rank=as_int(row.get("mode_rank")),
                            corr=as_float(row.get("corr")),
                            abs_corr=as_float(row.get("abs_corr")),
                        )
                    )

    return session_top, pooled_top, missing


def dataclass_rows(items: Iterable[object]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for item in items:
        out.append(dict(item.__dict__))
    return out


def summarize_channel_consistency(session_top: Sequence[SessionTop]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str, str, int, str, str], list[SessionTop]] = defaultdict(list)
    totals: Counter[tuple[str, str, str, str]] = Counter()
    for row in session_top:
        key_total = (row.dataset, row.modality, row.condition, row.feature)
        totals[key_total] += 1
        key = (
            row.dataset,
            row.modality,
            row.condition,
            row.feature,
            row.channel_index,
            row.channel_site,
            row.channel_label,
        )
        groups[key].append(row)

    out: list[dict[str, object]] = []
    for key, rows in groups.items():
        dataset, modality, condition, feature, channel_index, channel_site, channel_label = key
        total = totals[(dataset, modality, condition, feature)]
        abs_vals = [row.abs_corr for row in rows if math.isfinite(row.abs_corr)]
        mode_counts = Counter(row.mode_rank for row in rows)
        out.append(
            {
                "dataset": dataset,
                "modality": modality,
                "condition": condition,
                "feature": feature,
                "channel_index": channel_index,
                "channel_site": channel_site,
                "channel_label": channel_label,
                "top_count": len(rows),
                "session_feature_count": total,
                "top_fraction": len(rows) / total if total else math.nan,
                "mean_top_abs_corr": mean(abs_vals) if abs_vals else math.nan,
                "median_top_abs_corr": median(abs_vals) if abs_vals else math.nan,
                "dominant_mode_rank": mode_counts.most_common(1)[0][0] if mode_counts else "",
                "mode_rank_counts": ";".join(f"m{k}:{v}" for k, v in sorted(mode_counts.items())),
            }
        )
    return sorted(
        out,
        key=lambda r: (
            str(r["dataset"]),
            str(r["modality"]),
            str(r["condition"]),
            str(r["feature"]),
            -float(r["top_fraction"]),
            int(r["channel_index"]),
        ),
    )


def summarize_site_consistency(session_top: Sequence[SessionTop]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str, str, str], list[SessionTop]] = defaultdict(list)
    totals: Counter[tuple[str, str, str, str]] = Counter()
    for row in session_top:
        totals[(row.dataset, row.modality, row.condition, row.feature)] += 1
        groups[(row.dataset, row.modality, row.condition, row.feature, row.channel_site)].append(row)

    out: list[dict[str, object]] = []
    for key, rows in groups.items():
        dataset, modality, condition, feature, channel_site = key
        total = totals[(dataset, modality, condition, feature)]
        abs_vals = [row.abs_corr for row in rows if math.isfinite(row.abs_corr)]
        channel_counts = Counter(row.channel_label for row in rows)
        out.append(
            {
                "dataset": dataset,
                "modality": modality,
                "condition": condition,
                "feature": feature,
                "channel_site": channel_site,
                "top_count": len(rows),
                "session_feature_count": total,
                "top_fraction": len(rows) / total if total else math.nan,
                "mean_top_abs_corr": mean(abs_vals) if abs_vals else math.nan,
                "median_top_abs_corr": median(abs_vals) if abs_vals else math.nan,
                "dominant_channel": channel_counts.most_common(1)[0][0] if channel_counts else "",
                "channel_counts": ";".join(f"{k}:{v}" for k, v in channel_counts.most_common()),
            }
        )
    return sorted(
        out,
        key=lambda r: (
            str(r["dataset"]),
            str(r["modality"]),
            str(r["condition"]),
            str(r["feature"]),
            -float(r["top_fraction"]),
            str(r["channel_site"]),
        ),
    )


def summarize_overall_channel_distribution(session_top: Sequence[SessionTop]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str, int, str, str], list[SessionTop]] = defaultdict(list)
    totals: Counter[tuple[str, str, str]] = Counter()
    for row in session_top:
        totals[(row.dataset, row.modality, row.condition)] += 1
        groups[
            (
                row.dataset,
                row.modality,
                row.condition,
                row.channel_index,
                row.channel_site,
                row.channel_label,
            )
        ].append(row)

    out: list[dict[str, object]] = []
    for key, rows in groups.items():
        dataset, modality, condition, channel_index, channel_site, channel_label = key
        total = totals[(dataset, modality, condition)]
        abs_vals = [row.abs_corr for row in rows if math.isfinite(row.abs_corr)]
        feature_counts = Counter(row.feature for row in rows)
        mode_counts = Counter(row.mode_rank for row in rows)
        out.append(
            {
                "dataset": dataset,
                "modality": modality,
                "condition": condition,
                "channel_index": channel_index,
                "channel_site": channel_site,
                "channel_label": channel_label,
                "top_count": len(rows),
                "session_feature_count": total,
                "top_fraction": len(rows) / total if total else math.nan,
                "mean_top_abs_corr": mean(abs_vals) if abs_vals else math.nan,
                "median_top_abs_corr": median(abs_vals) if abs_vals else math.nan,
                "feature_counts": ";".join(f"{k}:{v}" for k, v in feature_counts.most_common()),
                "mode_rank_counts": ";".join(f"m{k}:{v}" for k, v in sorted(mode_counts.items())),
            }
        )
    return sorted(
        out,
        key=lambda r: (
            str(r["dataset"]),
            str(r["modality"]),
            str(r["condition"]),
            -float(r["top_fraction"]),
            int(r["channel_index"]),
        ),
    )


def summarize_overall_site_distribution(session_top: Sequence[SessionTop]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str, str], list[SessionTop]] = defaultdict(list)
    totals: Counter[tuple[str, str, str]] = Counter()
    for row in session_top:
        totals[(row.dataset, row.modality, row.condition)] += 1
        groups[(row.dataset, row.modality, row.condition, row.channel_site)].append(row)

    out: list[dict[str, object]] = []
    for key, rows in groups.items():
        dataset, modality, condition, channel_site = key
        total = totals[(dataset, modality, condition)]
        abs_vals = [row.abs_corr for row in rows if math.isfinite(row.abs_corr)]
        channel_counts = Counter(row.channel_label for row in rows)
        out.append(
            {
                "dataset": dataset,
                "modality": modality,
                "condition": condition,
                "channel_site": channel_site,
                "top_count": len(rows),
                "session_feature_count": total,
                "top_fraction": len(rows) / total if total else math.nan,
                "mean_top_abs_corr": mean(abs_vals) if abs_vals else math.nan,
                "median_top_abs_corr": median(abs_vals) if abs_vals else math.nan,
                "channel_counts": ";".join(f"{k}:{v}" for k, v in channel_counts.most_common()),
            }
        )
    return sorted(
        out,
        key=lambda r: (
            str(r["dataset"]),
            str(r["modality"]),
            str(r["condition"]),
            -float(r["top_fraction"]),
            str(r["channel_site"]),
        ),
    )


def top_channel_sets(session_top: Sequence[SessionTop], top_n: int) -> dict[tuple[str, str, str], list[str]]:
    groups: dict[tuple[str, str, str], Counter[str]] = defaultdict(Counter)
    for row in session_top:
        groups[(row.dataset, row.modality, row.condition)][row.channel_label] += 1
    return {
        key: [channel for channel, _count in counter.most_common(top_n)]
        for key, counter in groups.items()
    }


def summarize_cross_condition_overlap(session_top: Sequence[SessionTop], top_n: int) -> list[dict[str, object]]:
    sets = top_channel_sets(session_top, top_n)
    rows: list[dict[str, object]] = []
    keys = {(dataset, modality) for dataset, modality, _condition in sets}
    for dataset, modality in sorted(keys):
        abs_set = set(sets.get((dataset, modality, "abs"), []))
        complex_set = set(sets.get((dataset, modality, "complex_split"), []))
        union = abs_set | complex_set
        inter = abs_set & complex_set
        rows.append(
            {
                "dataset": dataset,
                "modality": modality,
                "top_n": top_n,
                "abs_top_channels": ";".join(sets.get((dataset, modality, "abs"), [])),
                "complex_split_top_channels": ";".join(sets.get((dataset, modality, "complex_split"), [])),
                "overlap_channels": ";".join(sorted(inter)),
                "jaccard": len(inter) / len(union) if union else math.nan,
                "same_first_channel": (
                    bool(sets.get((dataset, modality, "abs")))
                    and bool(sets.get((dataset, modality, "complex_split")))
                    and sets[(dataset, modality, "abs")][0] == sets[(dataset, modality, "complex_split")][0]
                ),
            }
        )
    return rows


def summarize_cross_modality_overlap(session_top: Sequence[SessionTop], top_n: int) -> list[dict[str, object]]:
    sets = top_channel_sets(session_top, top_n)
    rows: list[dict[str, object]] = []
    keys = {(dataset, condition) for dataset, _modality, condition in sets}
    for dataset, condition in sorted(keys):
        spkt_set = set(sets.get((dataset, "spkt", condition), []))
        mua_set = set(sets.get((dataset, "mua", condition), []))
        union = spkt_set | mua_set
        inter = spkt_set & mua_set
        rows.append(
            {
                "dataset": dataset,
                "condition": condition,
                "top_n": top_n,
                "spkt_top_channels": ";".join(sets.get((dataset, "spkt", condition), [])),
                "mua_top_channels": ";".join(sets.get((dataset, "mua", condition), [])),
                "overlap_channels": ";".join(sorted(inter)),
                "jaccard": len(inter) / len(union) if union else math.nan,
                "same_first_channel": (
                    bool(sets.get((dataset, "spkt", condition)))
                    and bool(sets.get((dataset, "mua", condition)))
                    and sets[(dataset, "spkt", condition)][0] == sets[(dataset, "mua", condition)][0]
                ),
            }
        )
    return rows


def esc(value: object) -> str:
    return html.escape(str(value), quote=True)


def color_scale(value: float, vmin: float = 0.0, vmax: float = 1.0) -> str:
    if not math.isfinite(value):
        return "#f2f4f7"
    t = max(0.0, min(1.0, (value - vmin) / (vmax - vmin if vmax > vmin else 1.0)))
    low = (240, 245, 249)
    high = (39, 107, 150)
    rgb = tuple(round(low[i] + t * (high[i] - low[i])) for i in range(3))
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def write_svg(path: Path, width: int, height: int, body: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
                '<rect width="100%" height="100%" fill="#ffffff"/>',
                *body,
                "</svg>",
                "",
            ]
        ),
        encoding="utf-8",
    )


def text(x: float, y: float, value: object, size: int = 12, fill: str = "#172033", anchor: str = "start", weight: str = "400") -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" fill="{fill}" text-anchor="{anchor}" font-weight="{weight}">{esc(value)}</text>'
    )


def rect(x: float, y: float, w: float, h: float, fill: str, stroke: str = "#ffffff") -> str:
    return f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" fill="{fill}" stroke="{stroke}" stroke-width="1"/>'


def plot_site_heatmap(site_rows: Sequence[dict[str, object]], output_path: Path) -> None:
    sites = sorted({str(row["channel_site"]) for row in site_rows})
    combos = sorted({(str(row["dataset"]), str(row["modality"]), str(row["condition"])) for row in site_rows})
    lookup: dict[tuple[str, str, str, str], float] = {}
    for row in site_rows:
        key = (
            str(row["dataset"]),
            str(row["modality"]),
            str(row["condition"]),
            str(row["channel_site"]),
        )
        lookup[key] = lookup.get(key, 0.0) + float(row["top_count"]) / float(row["session_feature_count"])

    cell_w = 64
    cell_h = 24
    left = 260
    top = 70
    width = left + cell_w * len(sites) + 40
    height = top + cell_h * len(combos) + 60
    body: list[str] = [text(20, 30, "P6 residual xcorr: top-channel site fraction", 18, weight="700")]
    for j, site in enumerate(sites):
        body.append(text(left + j * cell_w + cell_w / 2, 55, site, 11, anchor="middle", weight="700"))
    for i, combo in enumerate(combos):
        y = top + i * cell_h
        label = " | ".join(combo)
        body.append(text(20, y + 16, label, 11))
        for j, site in enumerate(sites):
            val = lookup.get((*combo, site), 0.0)
            body.append(rect(left + j * cell_w, y, cell_w, cell_h, color_scale(val)))
            if val >= 0.08:
                body.append(text(left + j * cell_w + cell_w / 2, y + 16, f"{val:.2f}", 10, "#ffffff" if val > 0.45 else "#172033", "middle"))
    write_svg(output_path, width, height, body)


def plot_channel_bars(channel_rows: Sequence[dict[str, object]], output_path: Path, max_channels: int = 4) -> None:
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in channel_rows:
        groups[(str(row["dataset"]), str(row["modality"]), str(row["condition"]))].append(row)
    combos = sorted(groups)

    left = 235
    top = 58
    row_h = 38
    bar_w = 260
    width = 760
    height = top + row_h * len(combos) + 50
    body: list[str] = [text(20, 30, "P6 residual xcorr: recurrent top channels", 18, weight="700")]
    for i, combo in enumerate(combos):
        y = top + i * row_h
        body.append(text(20, y + 22, " | ".join(combo), 11))
        rows = sorted(groups[combo], key=lambda r: (-float(r["top_fraction"]), int(r["channel_index"])))[:max_channels]
        x = left
        for row in rows:
            frac = float(row["top_fraction"])
            w = max(8, frac * bar_w)
            body.append(rect(x, y + 5, w, 18, color_scale(frac)))
            label = f"{row['channel_label']} {frac:.2f}"
            body.append(text(x + w + 6, y + 19, label, 10))
            x += 125
    write_svg(output_path, width, height, body)


def plot_overlap(rows: Sequence[dict[str, object]], output_path: Path, title: str, label_fields: tuple[str, ...]) -> None:
    left = 250
    top = 58
    row_h = 30
    bar_w = 300
    width = 640
    height = top + row_h * len(rows) + 45
    body: list[str] = [text(20, 30, title, 18, weight="700")]
    for i, row in enumerate(rows):
        y = top + i * row_h
        label = " | ".join(str(row[field]) for field in label_fields)
        val = float(row["jaccard"]) if row["jaccard"] != "" else math.nan
        body.append(text(20, y + 18, label, 11))
        body.append(rect(left, y + 4, bar_w, 18, "#eef2f6", "#d2d9e2"))
        body.append(rect(left, y + 4, bar_w * (val if math.isfinite(val) else 0), 18, color_scale(val)))
        body.append(text(left + bar_w + 10, y + 18, f"{val:.2f}" if math.isfinite(val) else "NA", 11))
    write_svg(output_path, width, height, body)


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    session_top, pooled_top, missing = discover_rows(args.processed_root, args.datasets)
    channel_summary = summarize_channel_consistency(session_top)
    site_summary = summarize_site_consistency(session_top)
    overall_channel = summarize_overall_channel_distribution(session_top)
    overall_site = summarize_overall_site_distribution(session_top)
    cross_condition = summarize_cross_condition_overlap(session_top, args.top_n_overlap)
    cross_modality = summarize_cross_modality_overlap(session_top, args.top_n_overlap)

    write_csv(output_dir / "session_top_channels.csv", dataclass_rows(session_top))
    write_csv(output_dir / "pooled_top_channels.csv", dataclass_rows(pooled_top))
    write_csv(output_dir / "channel_session_consistency.csv", channel_summary)
    write_csv(output_dir / "site_session_consistency.csv", site_summary)
    write_csv(output_dir / "overall_channel_distribution.csv", overall_channel)
    write_csv(output_dir / "overall_site_distribution.csv", overall_site)
    write_csv(output_dir / "cross_condition_channel_overlap.csv", cross_condition)
    write_csv(output_dir / "cross_modality_channel_overlap.csv", cross_modality)
    write_csv(output_dir / "missing_inputs.csv", missing)

    figure_dir = output_dir / "figures"
    plot_site_heatmap(overall_site, figure_dir / "site_top_fraction_heatmap.svg")
    plot_channel_bars(overall_channel, figure_dir / "channel_top_fraction_bars.svg")
    plot_overlap(
        cross_condition,
        figure_dir / "cross_condition_top_channel_overlap.svg",
        "P6 top-channel overlap: abs vs complex_split",
        ("dataset", "modality"),
    )
    plot_overlap(
        cross_modality,
        figure_dir / "cross_modality_top_channel_overlap.svg",
        "P6 top-channel overlap: spike vs MUA",
        ("dataset", "condition"),
    )

    summary = [
        "# P6 Residual Channel Consistency",
        "",
        f"- Generated at: `{datetime.now().isoformat(timespec='seconds')}`",
        f"- Datasets: `{', '.join(args.datasets)}`",
        f"- Session-top rows: `{len(session_top)}`",
        f"- Pooled-top rows: `{len(pooled_top)}`",
        f"- Missing input rows: `{len(missing)}`",
        "",
        "## Key Outputs",
        "",
        "- `session_top_channels.csv`: one top channel per dataset/modality/condition/session/feature.",
        "- `channel_session_consistency.csv`: recurrent top-channel fraction by channel.",
        "- `site_session_consistency.csv`: recurrent top-channel fraction by site label.",
        "- `overall_channel_distribution.csv`: channel distribution across all session-feature winners.",
        "- `overall_site_distribution.csv`: site distribution across all session-feature winners.",
        "- `cross_condition_channel_overlap.csv`: top-N channel overlap for abs vs complex_split.",
        "- `cross_modality_channel_overlap.csv`: top-N channel overlap for spike vs MUA.",
        "- `figures/site_top_fraction_heatmap.svg`: site-level spatial summary.",
        "- `figures/channel_top_fraction_bars.svg`: top channel labels by dataset/modality/condition.",
        "",
        "## Notes",
        "",
        "Spatial location is represented by `channel_index`, `channel_label`, and `channel_site` from cfg files.",
        "Physical depth/XYZ coordinates are not present in the current P6 CSVs.",
    ]
    (output_dir / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    print(f"Output dir          : {output_dir}")
    print(f"Session-top rows    : {len(session_top)}")
    print(f"Pooled-top rows     : {len(pooled_top)}")
    print(f"Missing input rows  : {len(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
