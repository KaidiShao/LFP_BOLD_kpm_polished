#!/usr/bin/env python3
"""Create flat summary PNGs from pipeline5 peak-statistics CSV files."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from PIL import Image, ImageDraw, ImageFont


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "f12m01", "e10gw1")
DEFAULT_PEAK_STAGE = "pipeline5_eigenfunction_peaks_by_state_maxabs"
REQUIRED_COLUMNS = (
    "state_code",
    "state_label",
    "component_idx",
    "event_mean_peak",
    "mean_event_minus_baseline",
    "cohen_d_paired_vs_baseline",
    "q_vs_baseline_paired_ttest_two_sided",
)


@dataclass(frozen=True)
class StatsFile:
    dataset: str
    variant: str
    method_tag: str
    path: Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--peak-stage", default=DEFAULT_PEAK_STAGE)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Default: <processed-root>/summary_figures/pipeline5_peak_statistics_maxabs",
    )
    parser.add_argument("--max-figures", type=int, default=None)
    return parser.parse_args()


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        Path(r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf"),
        Path(r"C:\Windows\Fonts\segoeuib.ttf" if bold else r"C:\Windows\Fonts\segoeui.ttf"),
    ]
    for path in candidates:
        if path.is_file():
            return ImageFont.truetype(str(path), size=size)
    return ImageFont.load_default()


FONT_TITLE = load_font(30, bold=True)
FONT_SUBTITLE = load_font(21, bold=True)
FONT_LABEL = load_font(17)
FONT_SMALL = load_font(14)
FONT_CELL = load_font(15, bold=True)


def find_stats_files(processed_root: Path, datasets: Sequence[str], peak_stage: str) -> List[StatsFile]:
    files: List[StatsFile] = []
    for dataset in datasets:
        root = processed_root / dataset / peak_stage
        if not root.is_dir():
            continue
        for stats_path in root.glob("*/*/*_peaks_stats.csv"):
            method_tag = stats_path.parent.name
            variant = stats_path.parent.parent.name
            files.append(StatsFile(dataset, variant, method_tag, stats_path))
    return sorted(files, key=lambda x: (x.dataset, x.variant, x.method_tag, str(x.path)))


def read_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No rows in {path}")
    missing = [col for col in REQUIRED_COLUMNS if col not in rows[0]]
    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")
    return rows


def to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def to_int(value: str) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return 0


def text_size(draw: ImageDraw.ImageDraw, text: str, font: ImageFont.ImageFont) -> Tuple[int, int]:
    bbox = draw.textbbox((0, 0), text, font=font)
    return bbox[2] - bbox[0], bbox[3] - bbox[1]


def lerp_color(a: Tuple[int, int, int], b: Tuple[int, int, int], t: float) -> Tuple[int, int, int]:
    t = max(0.0, min(1.0, t))
    return tuple(round(x + (y - x) * t) for x, y in zip(a, b))


def sequential_color(value: float, vmin: float, vmax: float) -> Tuple[int, int, int]:
    if not math.isfinite(value):
        return (238, 242, 246)
    if vmax <= vmin:
        t = 0.5
    else:
        t = (value - vmin) / (vmax - vmin)
    stops = [
        (0.00, (245, 248, 250)),
        (0.35, (186, 218, 232)),
        (0.70, (73, 144, 185)),
        (1.00, (24, 78, 119)),
    ]
    return interpolate_stops(stops, t)


def diverging_color(value: float, limit: float) -> Tuple[int, int, int]:
    if not math.isfinite(value):
        return (238, 242, 246)
    if limit <= 0:
        t = 0.5
    else:
        t = (value + limit) / (2 * limit)
    stops = [
        (0.00, (44, 92, 168)),
        (0.50, (247, 247, 247)),
        (1.00, (190, 65, 52)),
    ]
    return interpolate_stops(stops, t)


def q_color(value: float, vmax: float) -> Tuple[int, int, int]:
    if not math.isfinite(value):
        return (238, 242, 246)
    t = 0.0 if vmax <= 0 else value / vmax
    stops = [
        (0.00, (255, 255, 245)),
        (0.35, (254, 217, 118)),
        (0.70, (241, 105, 19)),
        (1.00, (127, 0, 0)),
    ]
    return interpolate_stops(stops, t)


def interpolate_stops(stops: Sequence[Tuple[float, Tuple[int, int, int]]], t: float) -> Tuple[int, int, int]:
    t = max(0.0, min(1.0, t))
    for idx in range(len(stops) - 1):
        left_t, left_c = stops[idx]
        right_t, right_c = stops[idx + 1]
        if left_t <= t <= right_t:
            local = 0.0 if right_t == left_t else (t - left_t) / (right_t - left_t)
            return lerp_color(left_c, right_c, local)
    return stops[-1][1]


def luminance(color: Tuple[int, int, int]) -> float:
    r, g, b = [x / 255.0 for x in color]
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def build_matrix(
    rows: List[Dict[str, str]],
    column: str,
    state_codes: Sequence[int],
    components: Sequence[int],
) -> List[List[float]]:
    matrix = [[math.nan for _ in components] for _ in state_codes]
    state_pos = {code: idx for idx, code in enumerate(state_codes)}
    comp_pos = {comp: idx for idx, comp in enumerate(components)}
    for row in rows:
        r = state_pos.get(to_int(row["state_code"]))
        c = comp_pos.get(to_int(row["component_idx"]))
        if r is not None and c is not None:
            matrix[r][c] = to_float(row[column])
    return matrix


def finite_values(matrix: Sequence[Sequence[float]]) -> List[float]:
    return [x for row in matrix for x in row if math.isfinite(x)]


def draw_colorbar(
    draw: ImageDraw.ImageDraw,
    x: int,
    y: int,
    width: int,
    height: int,
    mode: str,
    vmin: float,
    vmax: float,
) -> None:
    for dx in range(width):
        t = dx / max(1, width - 1)
        if mode == "diverging":
            color = diverging_color(vmin + (vmax - vmin) * t, max(abs(vmin), abs(vmax)))
        elif mode == "q":
            color = q_color(vmin + (vmax - vmin) * t, vmax)
        else:
            color = sequential_color(vmin + (vmax - vmin) * t, vmin, vmax)
        draw.line((x + dx, y, x + dx, y + height), fill=color)
    draw.rectangle((x, y, x + width, y + height), outline=(142, 154, 164), width=1)
    draw.text((x, y + height + 4), f"{vmin:.2g}", fill=(74, 85, 104), font=FONT_SMALL)
    max_label = f"{vmax:.2g}"
    tw, _ = text_size(draw, max_label, FONT_SMALL)
    draw.text((x + width - tw, y + height + 4), max_label, fill=(74, 85, 104), font=FONT_SMALL)


def draw_heatmap(
    draw: ImageDraw.ImageDraw,
    box: Tuple[int, int, int, int],
    title: str,
    matrix: Sequence[Sequence[float]],
    state_labels: Sequence[str],
    components: Sequence[int],
    mode: str,
    q_matrix: Optional[Sequence[Sequence[float]]] = None,
) -> None:
    x0, y0, x1, y1 = box
    draw.text((x0, y0), title, fill=(31, 41, 55), font=FONT_SUBTITLE)
    plot_x = x0 + 145
    plot_y = y0 + 50
    plot_w = x1 - plot_x - 32
    plot_h = y1 - plot_y - 58
    n_rows = len(state_labels)
    n_cols = len(components)
    cell_w = max(1, plot_w // max(1, n_cols))
    cell_h = max(1, plot_h // max(1, n_rows))
    heat_w = cell_w * n_cols
    heat_h = cell_h * n_rows
    values = finite_values(matrix)
    if values:
        raw_min, raw_max = min(values), max(values)
    else:
        raw_min, raw_max = 0.0, 1.0
    if mode == "diverging":
        limit = max(abs(raw_min), abs(raw_max), 1e-9)
        vmin, vmax = -limit, limit
    elif mode == "q":
        vmin, vmax = 0.0, max(1.0, min(16.0, raw_max))
    else:
        vmin, vmax = raw_min, raw_max if raw_max > raw_min else raw_min + 1.0

    for c, comp in enumerate(components):
        label = f"c{comp}"
        tw, _ = text_size(draw, label, FONT_SMALL)
        draw.text((plot_x + c * cell_w + (cell_w - tw) / 2, plot_y - 24), label, fill=(74, 85, 104), font=FONT_SMALL)
    for r, label in enumerate(state_labels):
        tw, th = text_size(draw, label, FONT_SMALL)
        draw.text((plot_x - tw - 10, plot_y + r * cell_h + (cell_h - th) / 2), label, fill=(74, 85, 104), font=FONT_SMALL)

    for r, row in enumerate(matrix):
        for c, value in enumerate(row):
            if mode == "diverging":
                color = diverging_color(value, max(abs(vmin), abs(vmax)))
            elif mode == "q":
                color = q_color(value, vmax)
            else:
                color = sequential_color(value, vmin, vmax)
            x = plot_x + c * cell_w
            y = plot_y + r * cell_h
            draw.rectangle((x, y, x + cell_w - 2, y + cell_h - 2), fill=color, outline=(255, 255, 255))
            if math.isfinite(value):
                label = f"{value:.1f}" if mode == "q" else f"{value:.2f}"
                text_color = (255, 255, 255) if luminance(color) < 0.48 else (17, 24, 39)
                tw, th = text_size(draw, label, FONT_CELL)
                draw.text((x + (cell_w - tw) / 2, y + (cell_h - th) / 2), label, fill=text_color, font=FONT_CELL)
                if q_matrix is not None:
                    q_val = q_matrix[r][c]
                    if math.isfinite(q_val) and q_val <= 0.05:
                        draw.text((x + cell_w - 15, y + 3), "*", fill=(17, 24, 39), font=FONT_LABEL)

    draw.rectangle((plot_x, plot_y, plot_x + heat_w, plot_y + heat_h), outline=(142, 154, 164), width=1)
    draw_colorbar(draw, plot_x, y1 - 35, min(260, heat_w), 12, mode, vmin, vmax)


def render_one(stats: StatsFile, output_dir: Path) -> Path:
    rows = read_rows(stats.path)
    state_codes = sorted({to_int(row["state_code"]) for row in rows})
    components = sorted({to_int(row["component_idx"]) for row in rows})
    labels_by_code: Dict[int, str] = {}
    for row in rows:
        labels_by_code.setdefault(to_int(row["state_code"]), row["state_label"])
    state_labels = [labels_by_code.get(code, f"state {code}") for code in state_codes]

    mean_peak = build_matrix(rows, "event_mean_peak", state_codes, components)
    delta = build_matrix(rows, "mean_event_minus_baseline", state_codes, components)
    cohen = build_matrix(rows, "cohen_d_paired_vs_baseline", state_codes, components)
    q_value = build_matrix(rows, "q_vs_baseline_paired_ttest_two_sided", state_codes, components)
    neg_log_q = [
        [-math.log10(max(x, 1e-16)) if math.isfinite(x) else math.nan for x in row]
        for row in q_value
    ]

    width, height = 1700, 1040
    image = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(image)
    title = f"{stats.dataset} / {stats.variant} / {stats.method_tag}"
    draw.text((34, 28), title, fill=(17, 24, 39), font=FONT_TITLE)
    draw.text((34, 68), str(stats.path), fill=(82, 96, 109), font=FONT_SMALL)
    draw.text((34, 94), "* marks q <= 0.05 in paired event-vs-baseline test", fill=(82, 96, 109), font=FONT_SMALL)

    panels = [
        (34, 130, 835, 555, "Event mean peak", mean_peak, "sequential", None),
        (880, 130, 1666, 555, "Mean event minus baseline", delta, "diverging", q_value),
        (34, 600, 835, 1010, "Cohen d paired vs baseline", cohen, "diverging", q_value),
        (880, 600, 1666, 1010, "-log10 q paired t-test", neg_log_q, "q", None),
    ]
    for x0, y0, x1, y1, panel_title, matrix, mode, q_matrix in panels:
        draw_heatmap(draw, (x0, y0, x1, y1), panel_title, matrix, state_labels, components, mode, q_matrix)

    output_name = f"{stats.dataset}__{stats.variant}__{stats.method_tag}__peak_stats_summary.png"
    output_path = output_dir / output_name
    image.save(output_path)
    return output_path


def main() -> int:
    args = parse_args()
    processed_root = args.processed_root.resolve()
    output_dir = args.output_dir or (
        processed_root / "summary_figures" / "pipeline5_peak_statistics_maxabs"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    stats_files = find_stats_files(processed_root, args.datasets, args.peak_stage)
    if args.max_figures is not None:
        stats_files = stats_files[: args.max_figures]

    print(f"Peak statistics summary dir: {output_dir}")
    print(f"Discovered {len(stats_files)} peak-stats CSV files.")

    index_rows: List[Dict[str, str]] = []
    ok = 0
    failed = 0
    for idx, stats in enumerate(stats_files, start=1):
        print(f"[{idx:03d}/{len(stats_files):03d}] {stats.dataset} / {stats.variant} / {stats.method_tag}")
        status = "OK"
        message = ""
        figure_path = ""
        try:
            figure_path = str(render_one(stats, output_dir))
            ok += 1
        except Exception as exc:  # noqa: BLE001 - keep batch running and record the failure.
            status = "FAILED"
            message = str(exc)
            failed += 1
            print(f"  FAILED: {message}")
        index_rows.append(
            {
                "dataset": stats.dataset,
                "variant": stats.variant,
                "method_tag": stats.method_tag,
                "source_csv": str(stats.path),
                "figure_png": figure_path,
                "status": status,
                "message": message,
            }
        )

    index_path = output_dir / "index.csv"
    with index_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["dataset", "variant", "method_tag", "source_csv", "figure_png", "status", "message"],
        )
        writer.writeheader()
        writer.writerows(index_rows)

    print(f"Done. OK: {ok} | failed: {failed}")
    print(f"Wrote index: {index_path}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
