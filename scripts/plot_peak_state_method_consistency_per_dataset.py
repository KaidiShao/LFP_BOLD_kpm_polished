#!/usr/bin/env python3
"""Draw per-dataset figures for peak-state method consistency summaries."""

from __future__ import annotations

import argparse
import csv
import html
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


DEFAULT_RESULTS_ROOT = Path("results")
DEFAULT_INPUT_GLOB = "peak_state_method_consistency_current_maxabs_k*_*"
EVENT_ORDER = ["theta", "gamma", "ripple", "theta-gamma", "sharp-wave-ripple"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", type=Path, default=DEFAULT_RESULTS_ROOT)
    parser.add_argument("--input-glob", default=DEFAULT_INPUT_GLOB)
    parser.add_argument(
        "--skip-png",
        action="store_true",
        help="Write SVG only. By default PNG is exported with Inkscape when available.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: str, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def sort_events(rows: Iterable[Dict[str, str]]) -> List[Dict[str, str]]:
    order = {label: idx for idx, label in enumerate(EVENT_ORDER)}
    return sorted(
        rows,
        key=lambda row: (order.get(row.get("state_label", ""), 999), as_int(row.get("state_code", "0"))),
    )


def color_lerp(v: float, low: Tuple[int, int, int], high: Tuple[int, int, int]) -> str:
    v = max(0.0, min(1.0, v))
    rgb = tuple(round(a + (b - a) * v) for a, b in zip(low, high))
    return f"rgb({rgb[0]},{rgb[1]},{rgb[2]})"


def svg_header(width: int, height: int) -> List[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text { font-family: Arial, Helvetica, sans-serif; fill: #1f2933; }",
        ".title { font-size: 20px; font-weight: 700; }",
        ".axis { font-size: 12px; fill: #52616b; }",
        ".label { font-size: 13px; }",
        ".value { font-size: 12px; font-weight: 700; }",
        ".grid { stroke: #d9e2ec; stroke-width: 1; }",
        ".axisline { stroke: #9fb3c8; stroke-width: 1; }",
        "</style>",
        '<rect x="0" y="0" width="100%" height="100%" fill="#fbfbfa" stroke="none"/>',
    ]


def write_svg(path: Path, lines: Sequence[str]) -> None:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def draw_event_match_rate(rows: List[Dict[str, str]], output_path: Path, title: str) -> None:
    rows = sort_events(rows)
    width = 900
    left = 185
    right = 70
    top = 64
    row_h = 44
    height = top + row_h * len(rows) + 48
    bar_w = width - left - right
    lines = svg_header(width, height)
    lines.append(f'<text x="24" y="34" class="title">{html.escape(title)}</text>')

    for tick in (0.0, 0.25, 0.5, 0.75, 1.0):
        x = left + bar_w * tick
        lines.append(f'<line x1="{x:.1f}" y1="{top - 18}" x2="{x:.1f}" y2="{height - 42}" class="grid"/>')
        lines.append(f'<text x="{x:.1f}" y="{height - 18}" class="axis" text-anchor="middle">{tick:.2f}</text>')
    lines.append(f'<line x1="{left}" y1="{height - 42}" x2="{left + bar_w}" y2="{height - 42}" class="axisline"/>')

    for idx, row in enumerate(rows):
        y = top + idx * row_h
        label = row.get("state_label", "")
        rate = as_float(row.get("match_rate", "0"))
        score = as_float(row.get("median_assignment_score_matched", "0"))
        fill = color_lerp(rate, (229, 240, 247), (32, 112, 162))
        lines.append(f'<text x="{left - 12}" y="{y + 23}" class="label" text-anchor="end">{html.escape(label)}</text>')
        lines.append(f'<rect x="{left}" y="{y + 6}" width="{bar_w * rate:.1f}" height="24" rx="3" fill="{fill}"/>')
        lines.append(f'<text x="{left + bar_w * rate + 8:.1f}" y="{y + 23}" class="value">{rate:.2f}</text>')
        lines.append(f'<text x="{width - 24}" y="{y + 23}" class="axis" text-anchor="end">median d={score:.2f}</text>')

    lines.append("</svg>")
    write_svg(output_path, lines)


def draw_run_matched_events(rows: List[Dict[str, str]], output_path: Path, title: str) -> None:
    width = 980
    left = 330
    right = 84
    top = 64
    row_h = 38
    height = top + row_h * len(rows) + 48
    bar_w = width - left - right
    max_events = max([as_int(row.get("event_count", "5"), 5) for row in rows] + [5])
    lines = svg_header(width, height)
    lines.append(f'<text x="24" y="34" class="title">{html.escape(title)}</text>')

    for tick in range(0, max_events + 1):
        x = left + bar_w * tick / max_events
        lines.append(f'<line x1="{x:.1f}" y1="{top - 18}" x2="{x:.1f}" y2="{height - 42}" class="grid"/>')
        lines.append(f'<text x="{x:.1f}" y="{height - 18}" class="axis" text-anchor="middle">{tick}</text>')
    lines.append(f'<line x1="{left}" y1="{height - 42}" x2="{left + bar_w}" y2="{height - 42}" class="axisline"/>')

    rows = sorted(rows, key=lambda row: row.get("run_id", ""))
    for idx, row in enumerate(rows):
        y = top + idx * row_h
        run_id = row.get("run_id", "")
        matched = as_float(row.get("matched_event_count", "0"))
        total = as_float(row.get("event_count", str(max_events)), max_events)
        fill = color_lerp(matched / max(total, 1.0), (232, 245, 233), (46, 125, 50))
        lines.append(f'<text x="{left - 12}" y="{y + 22}" class="label" text-anchor="end">{html.escape(run_id)}</text>')
        lines.append(f'<rect x="{left}" y="{y + 7}" width="{bar_w * matched / max_events:.1f}" height="22" rx="3" fill="{fill}"/>')
        lines.append(f'<text x="{left + bar_w * matched / max_events + 8:.1f}" y="{y + 22}" class="value">{matched:.0f}/{total:.0f}</text>')

    lines.append("</svg>")
    write_svg(output_path, lines)


def draw_method_event_heatmap(rows: List[Dict[str, str]], output_path: Path, title: str) -> None:
    methods = sorted({row.get("method", "") for row in rows if row.get("method", "")})
    event_rows = sort_events(
        [{"state_label": label, "state_code": str(idx + 1)} for idx, label in enumerate(EVENT_ORDER)]
    )
    events = [row["state_label"] for row in event_rows]
    values = {
        (row.get("method", ""), row.get("state_label", "")): as_float(row.get("match_rate", "0"))
        for row in rows
    }

    width = 900
    left = 150
    top = 92
    cell_w = 126
    cell_h = 48
    height = top + cell_h * len(methods) + 56
    lines = svg_header(width, height)
    lines.append(f'<text x="24" y="34" class="title">{html.escape(title)}</text>')
    lines.append('<text x="24" y="58" class="axis">cell value = method-specific event match rate</text>')

    for j, event in enumerate(events):
        x = left + j * cell_w + cell_w / 2
        lines.append(f'<text x="{x:.1f}" y="{top - 16}" class="axis" text-anchor="middle">{html.escape(event)}</text>')

    for i, method in enumerate(methods):
        y = top + i * cell_h
        lines.append(f'<text x="{left - 12}" y="{y + 30}" class="label" text-anchor="end">{html.escape(method)}</text>')
        for j, event in enumerate(events):
            x = left + j * cell_w
            value = values.get((method, event), 0.0)
            fill = color_lerp(value, (248, 244, 238), (188, 83, 39))
            lines.append(f'<rect x="{x}" y="{y}" width="{cell_w - 3}" height="{cell_h - 3}" rx="3" fill="{fill}"/>')
            lines.append(f'<text x="{x + cell_w / 2:.1f}" y="{y + 30}" class="value" text-anchor="middle">{value:.2f}</text>')

    lines.append("</svg>")
    write_svg(output_path, lines)


def resolve_inkscape() -> Optional[str]:
    for name in ("inkscape.com", "inkscape.exe", "inkscape"):
        exe = shutil.which(name)
        if exe:
            return exe
    return None


def export_png(svg_path: Path, inkscape_exe: Optional[str]) -> bool:
    if not inkscape_exe:
        return False
    png_path = svg_path.with_suffix(".png")
    command = [
        inkscape_exe,
        str(svg_path),
        "--export-type=png",
        f"--export-filename={png_path}",
    ]
    try:
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except (OSError, subprocess.CalledProcessError):
        return False
    return png_path.is_file()


def plot_one_dir(result_dir: Path, inkscape_exe: Optional[str], skip_png: bool) -> int:
    overall_path = result_dir / "overall_event_consistency.csv"
    run_path = result_dir / "run_summary.csv"
    method_event_path = result_dir / "method_event_consistency.csv"
    if not (overall_path.is_file() and run_path.is_file() and method_event_path.is_file()):
        return 0

    title_prefix = result_dir.name.replace("peak_state_method_consistency_current_maxabs_", "")
    svg_paths = [
        result_dir / "event_match_rate.svg",
        result_dir / "run_matched_events.svg",
        result_dir / "method_event_match_rate.svg",
    ]
    draw_event_match_rate(read_csv(overall_path), svg_paths[0], f"{title_prefix}: event match rate")
    draw_run_matched_events(read_csv(run_path), svg_paths[1], f"{title_prefix}: matched events by run")
    draw_method_event_heatmap(read_csv(method_event_path), svg_paths[2], f"{title_prefix}: method x event match rate")

    png_count = 0
    if not skip_png:
        for svg_path in svg_paths:
            png_count += int(export_png(svg_path, inkscape_exe))
    return png_count


def main() -> int:
    args = parse_args()
    results_root = args.results_root.resolve()
    inkscape_exe = None if args.skip_png else resolve_inkscape()

    dirs = sorted(path for path in results_root.glob(args.input_glob) if path.is_dir())
    plotted = 0
    png_count = 0
    skipped = 0
    for result_dir in dirs:
        before = png_count
        made_png = plot_one_dir(result_dir, inkscape_exe, args.skip_png)
        if made_png or (result_dir / "event_match_rate.svg").is_file():
            plotted += 1
            png_count += made_png
        else:
            skipped += 1

    print(f"Plotted per-dataset dirs: {plotted}")
    print(f"Skipped dirs without complete CSV inputs: {skipped}")
    if args.skip_png:
        print("PNG export skipped by --skip-png.")
    elif inkscape_exe:
        print(f"PNG files written: {png_count}")
    else:
        print("PNG export skipped: Inkscape not found on PATH.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
