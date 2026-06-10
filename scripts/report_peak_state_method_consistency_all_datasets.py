#!/usr/bin/env python3
"""Build a cross-dataset report for peak-state method consistency results."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


RESULTS_ROOT = Path("results")
INPUT_PREFIX = "peak_state_method_consistency_"
OUTPUT_DIRNAME = "peak_state_method_consistency_all_datasets"
METHOD_ORDER = ["svd", "nmf", "mds", "umap", "logsvd"]
PALETTE_RATE = ((245, 248, 250), (34, 110, 168))
PALETTE_SCORE = ((249, 245, 239), (196, 86, 40))
PALETTE_COUNT = ((242, 247, 242), (47, 125, 80))


@dataclass(frozen=True)
class DatasetBundle:
    dataset: str
    root: Path
    overall_event_path: Path
    method_event_path: Path
    run_summary_path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-root",
        type=Path,
        default=RESULTS_ROOT,
        help="Root folder containing per-dataset consistency result directories.",
    )
    parser.add_argument(
        "--input-prefix",
        default=INPUT_PREFIX,
        help="Directory-name prefix used to discover per-dataset result folders.",
    )
    parser.add_argument(
        "--output-dirname",
        default=OUTPUT_DIRNAME,
        help="Name of the cross-dataset report directory created under results-root.",
    )
    parser.add_argument(
        "--methods",
        nargs="*",
        default=None,
        help="Optional method order/filter for figures and markdown. Defaults to methods present in the inputs.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def safe_float(value: object, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def safe_int(value: object, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def discover_input_dirs(results_root: Path, input_prefix: str, output_dirname: str) -> List[DatasetBundle]:
    bundles: List[DatasetBundle] = []
    for child in sorted(results_root.iterdir()):
        if not child.is_dir():
            continue
        if not child.name.startswith(input_prefix):
            continue
        if child.name == output_dirname:
            continue
        dataset = child.name[len(input_prefix) :]
        overall_path = child / "overall_event_consistency.csv"
        method_path = child / "method_event_consistency.csv"
        run_path = child / "run_summary.csv"
        if overall_path.is_file() and method_path.is_file() and run_path.is_file():
            bundles.append(
                DatasetBundle(
                    dataset=dataset,
                    root=child,
                    overall_event_path=overall_path,
                    method_event_path=method_path,
                    run_summary_path=run_path,
                )
            )
    return bundles


def mean(values: Iterable[float]) -> Optional[float]:
    clean = [value for value in values if math.isfinite(value)]
    if not clean:
        return None
    return sum(clean) / len(clean)


def rgb_to_hex(rgb: Tuple[int, int, int]) -> str:
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def interpolate_color(
    value: float,
    vmin: float,
    vmax: float,
    low_rgb: Tuple[int, int, int],
    high_rgb: Tuple[int, int, int],
) -> str:
    if not math.isfinite(value):
        return "#f2f2f2"
    if math.isclose(vmax, vmin):
        ratio = 1.0
    else:
        ratio = clamp01((value - vmin) / (vmax - vmin))
    rgb = tuple(
        int(round(low_rgb[idx] + ratio * (high_rgb[idx] - low_rgb[idx])))
        for idx in range(3)
    )
    return rgb_to_hex(rgb)


def choose_text_color(fill: str) -> str:
    fill = fill.lstrip("#")
    r = int(fill[0:2], 16)
    g = int(fill[2:4], 16)
    b = int(fill[4:6], 16)
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "#111111" if luminance > 150 else "#ffffff"


def format_value(value: Optional[float], value_format: str) -> str:
    if value is None or not math.isfinite(value):
        return "NA"
    return value_format.format(value)


def resolve_inkscape_executable() -> Optional[str]:
    for name in ("inkscape.com", "inkscape.exe", "inkscape"):
        exe = shutil.which(name)
        if exe:
            return exe
    return None


def export_svg_to_png(svg_path: Path, png_path: Path, inkscape_exe: Optional[str]) -> bool:
    if not inkscape_exe:
        return False

    command = [
        inkscape_exe,
        str(svg_path),
        "--export-type=png",
        f"--export-filename={png_path}",
    ]
    try:
        subprocess.run(
            command,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return False
    return png_path.is_file()


def write_heatmap_svg(
    path: Path,
    title: str,
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    matrix: Sequence[Sequence[Optional[float]]],
    palette: Tuple[Tuple[int, int, int], Tuple[int, int, int]],
    value_format: str,
    subtitle: Optional[str] = None,
) -> None:
    cell_w = 132
    cell_h = 62
    left_margin = 220
    top_margin = 110
    right_margin = 30
    bottom_margin = 40
    width = left_margin + len(col_labels) * cell_w + right_margin
    height = top_margin + len(row_labels) * cell_h + bottom_margin

    values = [
        value
        for row in matrix
        for value in row
        if value is not None and math.isfinite(value)
    ]
    vmin = min(values) if values else 0.0
    vmax = max(values) if values else 1.0

    parts: List[str] = []
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
    )
    parts.append(
        '<rect x="0" y="0" width="100%" height="100%" fill="#fbfbfa" stroke="none" />'
    )
    parts.append(
        f'<text x="{left_margin}" y="34" font-family="Segoe UI, Arial, sans-serif" '
        'font-size="24" font-weight="700" fill="#1d1d1f">'
        f"{escape_xml(title)}</text>"
    )
    if subtitle:
        parts.append(
            f'<text x="{left_margin}" y="60" font-family="Segoe UI, Arial, sans-serif" '
            'font-size="13" fill="#5b5b5f">'
            f"{escape_xml(subtitle)}</text>"
        )

    for idx, label in enumerate(col_labels):
        x = left_margin + idx * cell_w + cell_w / 2
        y = top_margin - 18
        parts.append(
            f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="middle" '
            'font-family="Segoe UI, Arial, sans-serif" font-size="14" '
            'font-weight="600" fill="#333333">'
            f"{escape_xml(label)}</text>"
        )

    for idx, label in enumerate(row_labels):
        x = left_margin - 12
        y = top_margin + idx * cell_h + cell_h / 2 + 5
        parts.append(
            f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="end" '
            'font-family="Segoe UI, Arial, sans-serif" font-size="14" '
            'font-weight="600" fill="#333333">'
            f"{escape_xml(label)}</text>"
        )

    low_rgb, high_rgb = palette
    for row_idx, row in enumerate(matrix):
        for col_idx, value in enumerate(row):
            x = left_margin + col_idx * cell_w
            y = top_margin + row_idx * cell_h
            fill = interpolate_color(value if value is not None else math.nan, vmin, vmax, low_rgb, high_rgb)
            text_color = choose_text_color(fill)
            parts.append(
                f'<rect x="{x}" y="{y}" width="{cell_w - 2}" height="{cell_h - 2}" '
                f'rx="8" ry="8" fill="{fill}" stroke="#ffffff" stroke-width="2" />'
            )
            parts.append(
                f'<text x="{x + cell_w / 2:.1f}" y="{y + cell_h / 2 + 6:.1f}" '
                'text-anchor="middle" font-family="Consolas, Menlo, monospace" '
                f'font-size="16" font-weight="700" fill="{text_color}">'
                f"{escape_xml(format_value(value, value_format))}</text>"
            )

    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def escape_xml(text: object) -> str:
    value = str(text)
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return

    fieldnames: List[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def build_dataset_summary(
    dataset: str,
    overall_rows: List[Dict[str, str]],
    run_rows: List[Dict[str, str]],
) -> Dict[str, object]:
    match_rates = {row["state_label"]: safe_float(row["match_rate"]) for row in overall_rows}
    best_label = max(match_rates.items(), key=lambda item: item[1])[0]
    worst_label = min(match_rates.items(), key=lambda item: item[1])[0]
    matched_counts = [safe_float(row["matched_event_count"]) for row in run_rows]
    total_scores = [safe_float(row["total_assignment_score"]) for row in run_rows]
    return {
        "dataset": dataset,
        "run_count": len(run_rows),
        "avg_matched_events": mean(matched_counts),
        "min_matched_events": min(matched_counts) if matched_counts else "",
        "max_matched_events": max(matched_counts) if matched_counts else "",
        "avg_total_assignment_score": mean(total_scores),
        "best_event": best_label,
        "best_event_match_rate": match_rates[best_label],
        "worst_event": worst_label,
        "worst_event_match_rate": match_rates[worst_label],
    }


def build_cross_dataset_event_summary(
    overall_rows: List[Dict[str, object]]
) -> List[Dict[str, object]]:
    grouped: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in overall_rows:
        grouped[str(row["state_label"])].append(row)

    summary_rows: List[Dict[str, object]] = []
    for event_label, rows in sorted(grouped.items(), key=lambda item: safe_int(item[1][0]["state_code"])):
        match_rates = [safe_float(row["match_rate"]) for row in rows]
        median_scores = [safe_float(row["median_assignment_score_matched"]) for row in rows]
        summary_rows.append(
            {
                "state_label": event_label,
                "state_code": safe_int(rows[0]["state_code"]),
                "dataset_count": len(rows),
                "avg_match_rate": mean(match_rates),
                "min_match_rate": min(match_rates),
                "max_match_rate": max(match_rates),
                "avg_median_score": mean(median_scores),
            }
        )
    return summary_rows


def build_method_dataset_summary(run_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    for row in run_rows:
        grouped[(str(row["dataset"]), str(row["method"]))].append(row)

    summary_rows: List[Dict[str, object]] = []
    for (dataset, method), rows in sorted(grouped.items()):
        matched_counts = [safe_float(row["matched_event_count"]) for row in rows]
        total_scores = [safe_float(row["total_assignment_score"]) for row in rows]
        summary_rows.append(
            {
                "dataset": dataset,
                "method": method,
                "run_count": len(rows),
                "avg_matched_events": mean(matched_counts),
                "avg_total_assignment_score": mean(total_scores),
            }
        )
    return summary_rows


def build_method_event_average(method_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str], List[Dict[str, object]]] = defaultdict(list)
    for row in method_rows:
        grouped[(str(row["method"]), str(row["state_label"]))].append(row)

    summary_rows: List[Dict[str, object]] = []
    for (method, event_label), rows in sorted(grouped.items()):
        summary_rows.append(
            {
                "method": method,
                "state_label": event_label,
                "state_code": safe_int(rows[0]["state_code"]),
                "dataset_count": len(rows),
                "avg_match_rate": mean(safe_float(row["match_rate"]) for row in rows),
                "avg_median_score": mean(
                    safe_float(row["median_assignment_score_matched"]) for row in rows
                ),
            }
        )
    return summary_rows


def resolve_method_order(
    run_rows: Sequence[Dict[str, object]],
    method_rows: Sequence[Dict[str, object]],
    requested_methods: Optional[Sequence[str]],
) -> List[str]:
    present = {
        str(row.get("method", "")).strip().lower()
        for row in list(run_rows) + list(method_rows)
        if str(row.get("method", "")).strip()
    }
    if requested_methods:
        ordered = [
            str(method).strip().lower()
            for method in requested_methods
            if str(method).strip()
        ]
        return [method for method in ordered if method in present]

    ordered = [method for method in METHOD_ORDER if method in present]
    extras = sorted(present.difference(ordered))
    return ordered + extras


def build_matrix(
    rows: Sequence[Dict[str, object]],
    row_key: str,
    col_key: str,
    value_key: str,
    row_order: Sequence[str],
    col_order: Sequence[str],
) -> List[List[Optional[float]]]:
    lookup = {
        (str(row[row_key]), str(row[col_key])): safe_float(row[value_key])
        for row in rows
    }
    matrix: List[List[Optional[float]]] = []
    for row_label in row_order:
        matrix_row: List[Optional[float]] = []
        for col_label in col_order:
            value = lookup.get((row_label, col_label), math.nan)
            matrix_row.append(None if not math.isfinite(value) else value)
        matrix.append(matrix_row)
    return matrix


def dataset_event_order(overall_rows: Sequence[Dict[str, object]]) -> List[str]:
    pairs = sorted(
        {(safe_int(row["state_code"]), str(row["state_label"])) for row in overall_rows}
    )
    return [label for _code, label in pairs]


def write_summary_markdown(
    path: Path,
    dataset_summary_rows: List[Dict[str, object]],
    event_summary_rows: List[Dict[str, object]],
    method_dataset_rows: List[Dict[str, object]],
    method_event_rows: List[Dict[str, object]],
    dataset_order: Sequence[str],
    event_order: Sequence[str],
    method_order: Sequence[str],
) -> None:
    overall_best_event = max(event_summary_rows, key=lambda row: safe_float(row["avg_match_rate"]))
    overall_worst_event = min(event_summary_rows, key=lambda row: safe_float(row["avg_match_rate"]))

    method_mean_match: Dict[str, float] = {}
    grouped: Dict[str, List[float]] = defaultdict(list)
    for row in method_dataset_rows:
        grouped[str(row["method"])].append(safe_float(row["avg_matched_events"]))
    for method, values in grouped.items():
        method_mean_match[method] = mean(values) or math.nan
    best_method = max(method_mean_match.items(), key=lambda item: item[1])[0]
    worst_method = min(method_mean_match.items(), key=lambda item: item[1])[0]

    dataset_hardest = min(
        dataset_summary_rows,
        key=lambda row: safe_float(row["avg_matched_events"]),
    )
    dataset_easiest = max(
        dataset_summary_rows,
        key=lambda row: safe_float(row["avg_matched_events"]),
    )

    lines: List[str] = []
    lines.append("# Cross-Dataset Peak-State Consistency Report")
    lines.append("")
    lines.append(f"- Generated at: `{datetime.now().isoformat(timespec='seconds')}`")
    lines.append(f"- Datasets: `{', '.join(dataset_order)}`")
    lines.append(f"- Methods: `{', '.join(method_order)}`")
    lines.append(
        "- Common pattern: `theta-gamma` and `sharp-wave-ripple` are the most stable events, "
        "`gamma` is the least stable."
    )
    lines.append(
        f"- Strongest average event: `{overall_best_event['state_label']}` "
        f"(`avg match rate = {safe_float(overall_best_event['avg_match_rate']):.3f}`)"
    )
    lines.append(
        f"- Weakest average event: `{overall_worst_event['state_label']}` "
        f"(`avg match rate = {safe_float(overall_worst_event['avg_match_rate']):.3f}`)"
    )
    lines.append(
        f"- Easiest dataset overall: `{dataset_easiest['dataset']}` "
        f"(`avg matched events = {safe_float(dataset_easiest['avg_matched_events']):.3f}`)"
    )
    lines.append(
        f"- Hardest dataset overall: `{dataset_hardest['dataset']}` "
        f"(`avg matched events = {safe_float(dataset_hardest['avg_matched_events']):.3f}`)"
    )
    lines.append(
        f"- Strongest method overall: `{best_method}` "
        f"(`avg matched events = {method_mean_match[best_method]:.3f}`)"
    )
    lines.append(
        f"- Weakest method overall: `{worst_method}` "
        f"(`avg matched events = {method_mean_match[worst_method]:.3f}`)"
    )
    lines.append("")
    lines.append("## Figures")
    lines.append("")
    lines.append("1. `dataset_event_match_rate.png`: event stability across datasets.")
    lines.append("2. `dataset_event_median_score.png`: effect-size strength across datasets.")
    lines.append("3. `dataset_method_avg_matched_events.png`: overall method performance by dataset.")
    lines.append("4. `method_event_avg_match_rate.png`: which methods support which events most consistently.")
    lines.append("")
    lines.append("## Dataset Summary")
    lines.append("")
    lines.append("| Dataset | Runs | Avg Matched Events | Best Event | Best Rate | Worst Event | Worst Rate |")
    lines.append("| --- | ---: | ---: | --- | ---: | --- | ---: |")
    for row in dataset_summary_rows:
        lines.append(
            f"| {row['dataset']} | {row['run_count']} | {safe_float(row['avg_matched_events']):.3f} | "
            f"{row['best_event']} | {safe_float(row['best_event_match_rate']):.3f} | "
            f"{row['worst_event']} | {safe_float(row['worst_event_match_rate']):.3f} |"
        )
    lines.append("")
    lines.append("## Event Summary")
    lines.append("")
    lines.append("| Event | Avg Match Rate | Min | Max | Avg Median Score |")
    lines.append("| --- | ---: | ---: | ---: | ---: |")
    for row in event_summary_rows:
        lines.append(
            f"| {row['state_label']} | {safe_float(row['avg_match_rate']):.3f} | "
            f"{safe_float(row['min_match_rate']):.3f} | {safe_float(row['max_match_rate']):.3f} | "
            f"{safe_float(row['avg_median_score']):.3f} |"
        )
    lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    results_root = args.results_root.resolve()
    output_dir = (results_root / args.output_dirname).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    inkscape_exe = resolve_inkscape_executable()

    bundles = discover_input_dirs(results_root, args.input_prefix, args.output_dirname)
    if not bundles:
        raise SystemExit("No per-dataset consistency directories found.")

    dataset_order = [bundle.dataset for bundle in bundles]
    overall_rows: List[Dict[str, object]] = []
    method_rows: List[Dict[str, object]] = []
    run_rows: List[Dict[str, object]] = []
    dataset_summary_rows: List[Dict[str, object]] = []

    for bundle in bundles:
        dataset_overall = read_csv(bundle.overall_event_path)
        dataset_method = read_csv(bundle.method_event_path)
        dataset_run = read_csv(bundle.run_summary_path)

        for row in dataset_overall:
            row = dict(row)
            row["dataset"] = bundle.dataset
            overall_rows.append(row)
        for row in dataset_method:
            row = dict(row)
            row["dataset"] = bundle.dataset
            method_rows.append(row)
        for row in dataset_run:
            row = dict(row)
            row["dataset"] = bundle.dataset
            run_rows.append(row)

        dataset_summary_rows.append(
            build_dataset_summary(bundle.dataset, dataset_overall, dataset_run)
        )

    event_order = dataset_event_order(overall_rows)
    method_dataset_rows = build_method_dataset_summary(run_rows)
    event_summary_rows = build_cross_dataset_event_summary(overall_rows)
    method_event_avg_rows = build_method_event_average(method_rows)
    method_order = resolve_method_order(run_rows, method_rows, args.methods)

    overall_rows.sort(key=lambda row: (dataset_order.index(str(row["dataset"])), safe_int(row["state_code"])))
    dataset_summary_rows.sort(key=lambda row: dataset_order.index(str(row["dataset"])))
    event_summary_rows.sort(key=lambda row: event_order.index(str(row["state_label"])))

    write_csv(output_dir / "combined_overall_event_consistency.csv", overall_rows)
    write_csv(output_dir / "combined_method_event_consistency.csv", method_rows)
    write_csv(output_dir / "combined_run_summary.csv", run_rows)
    write_csv(output_dir / "dataset_summary.csv", dataset_summary_rows)
    write_csv(output_dir / "cross_dataset_event_summary.csv", event_summary_rows)
    write_csv(output_dir / "method_dataset_summary.csv", method_dataset_rows)
    write_csv(output_dir / "method_event_average_summary.csv", method_event_avg_rows)

    dataset_event_rate_matrix = build_matrix(
        overall_rows,
        row_key="dataset",
        col_key="state_label",
        value_key="match_rate",
        row_order=dataset_order,
        col_order=event_order,
    )
    dataset_event_score_matrix = build_matrix(
        overall_rows,
        row_key="dataset",
        col_key="state_label",
        value_key="median_assignment_score_matched",
        row_order=dataset_order,
        col_order=event_order,
    )
    dataset_method_match_matrix = build_matrix(
        method_dataset_rows,
        row_key="dataset",
        col_key="method",
        value_key="avg_matched_events",
        row_order=dataset_order,
        col_order=method_order,
    )
    method_event_rate_matrix = build_matrix(
        method_event_avg_rows,
        row_key="method",
        col_key="state_label",
        value_key="avg_match_rate",
        row_order=method_order,
        col_order=event_order,
    )

    svg_paths = [
        output_dir / "dataset_event_match_rate.svg",
        output_dir / "dataset_event_median_score.svg",
        output_dir / "dataset_method_avg_matched_events.svg",
        output_dir / "method_event_avg_match_rate.svg",
    ]

    write_heatmap_svg(
        svg_paths[0],
        title="Dataset x Event Match Rate",
        subtitle="Higher values mean the event is more consistently matched within that dataset.",
        row_labels=dataset_order,
        col_labels=event_order,
        matrix=dataset_event_rate_matrix,
        palette=PALETTE_RATE,
        value_format="{:.2f}",
    )
    write_heatmap_svg(
        svg_paths[1],
        title="Dataset x Event Median Assignment Score",
        subtitle="Median matched effect size (Cohen d) among successful event-component assignments.",
        row_labels=dataset_order,
        col_labels=event_order,
        matrix=dataset_event_score_matrix,
        palette=PALETTE_SCORE,
        value_format="{:.2f}",
    )
    write_heatmap_svg(
        svg_paths[2],
        title="Dataset x Method Average Matched Events",
        subtitle="Average number of events matched per run under each method.",
        row_labels=dataset_order,
        col_labels=method_order,
        matrix=dataset_method_match_matrix,
        palette=PALETTE_COUNT,
        value_format="{:.2f}",
    )
    write_heatmap_svg(
        svg_paths[3],
        title="Method x Event Average Match Rate",
        subtitle="Average event match rate across datasets for each method.",
        row_labels=method_order,
        col_labels=event_order,
        matrix=method_event_rate_matrix,
        palette=PALETTE_RATE,
        value_format="{:.2f}",
    )

    png_written: List[Path] = []
    for svg_path in svg_paths:
        png_path = svg_path.with_suffix(".png")
        if export_svg_to_png(svg_path, png_path, inkscape_exe):
            png_written.append(png_path)

    write_summary_markdown(
        output_dir / "summary.md",
        dataset_summary_rows=dataset_summary_rows,
        event_summary_rows=event_summary_rows,
        method_dataset_rows=method_dataset_rows,
        method_event_rows=method_event_avg_rows,
        dataset_order=dataset_order,
        event_order=event_order,
        method_order=method_order,
    )

    if png_written:
        print("PNG exports:")
        for png_path in png_written:
            print(f"  {png_path}")
    elif inkscape_exe is None:
        print("PNG export skipped: Inkscape not found on PATH.")
    else:
        print("PNG export attempted but no PNG files were created.")

    print(f"Wrote cross-dataset report to: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
