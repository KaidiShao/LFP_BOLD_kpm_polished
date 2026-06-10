#!/usr/bin/env python3
"""Plot P5 event-family component activity summary figures.

The script intentionally uses only the Python standard library to avoid adding
plotting dependencies. It writes SVG files directly and exports PNGs with
Inkscape when available.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


DEFAULT_INPUT_DIR = Path("results") / "peak_event_family_component_activity_current_maxabs"
METHOD_ORDER = ("svd", "nmf", "mds", "umap")
FAMILY_ORDER = ("theta", "gamma", "ripple")
METHOD_COLORS = {
    "svd": "#2b6cb0",
    "nmf": "#dd6b20",
    "mds": "#2f855a",
    "umap": "#805ad5",
}
STATUS_COLORS = {
    "inactive": "#d9dee6",
    "nonselective": "#f2a65a",
    "all_active": "#4f8fcf",
    "similar": "#35a46b",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--no-png", action="store_true")
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


def esc(text: object) -> str:
    return html.escape(str(text), quote=True)


def rgb_to_hex(rgb: Tuple[int, int, int]) -> str:
    return "#{:02x}{:02x}{:02x}".format(*rgb)


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def interpolate_color(value: float, low: Tuple[int, int, int], high: Tuple[int, int, int]) -> str:
    t = clamp01(value)
    rgb = tuple(round(low[i] + (high[i] - low[i]) * t) for i in range(3))
    return rgb_to_hex(rgb)  # type: ignore[arg-type]


def svg_doc(width: int, height: int, body: Sequence[str]) -> str:
    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}">',
            '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
            *body,
            "</svg>",
            "",
        ]
    )


def text(
    x: float,
    y: float,
    value: object,
    size: int = 12,
    fill: str = "#172033",
    anchor: str = "start",
    weight: str = "400",
    rotate: Optional[float] = None,
) -> str:
    transform = ""
    if rotate is not None:
        transform = f' transform="rotate({rotate} {x} {y})"'
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, Helvetica, sans-serif" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" '
        f'fill="{fill}"{transform}>{esc(value)}</text>'
    )


def rect(
    x: float,
    y: float,
    w: float,
    h: float,
    fill: str,
    stroke: str = "#ffffff",
    sw: float = 1,
    rx: float = 0,
) -> str:
    return (
        f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
        f'rx="{rx:.1f}" fill="{fill}" stroke="{stroke}" stroke-width="{sw:.1f}"/>'
    )


def line(
    x1: float,
    y1: float,
    x2: float,
    y2: float,
    stroke: str = "#a0a8b8",
    sw: float = 1,
    dash: Optional[str] = None,
) -> str:
    dash_attr = "" if dash is None else f' stroke-dasharray="{dash}"'
    return (
        f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
        f'stroke="{stroke}" stroke-width="{sw:.1f}"{dash_attr}/>'
    )


def circle(cx: float, cy: float, r: float, fill: str, stroke: str = "#ffffff", opacity: float = 1.0) -> str:
    return (
        f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.1f}" fill="{fill}" '
        f'stroke="{stroke}" stroke-width="0.8" opacity="{opacity:.3f}"/>'
    )


def write_svg(path: Path, width: int, height: int, body: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(svg_doc(width, height, body), encoding="utf-8")


def find_inkscape() -> Optional[str]:
    exe = shutil.which("inkscape") or shutil.which("inkscape.com")
    if exe:
        return exe
    candidates = [
        Path(r"C:\Program Files\Inkscape\bin\inkscape.com"),
        Path(r"C:\Program Files\Inkscape\bin\inkscape.exe"),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return str(candidate)
    return None


def export_png(svg_path: Path, png_path: Path, inkscape: Optional[str]) -> bool:
    if inkscape is None:
        return False
    cmd = [inkscape, str(svg_path), "--export-type=png", f"--export-filename={png_path}"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception:
        return False
    return png_path.is_file()


def plot_heatmap(
    summary_rows: Sequence[Dict[str, str]],
    value_key: str,
    title: str,
    out_path: Path,
    subtitle: str = "Rows are event-name families; columns are P5 methods.",
) -> Path:
    lookup = {
        (row["event_family"], row["method"]): safe_float(row[value_key])
        for row in summary_rows
    }
    width, height = 760, 390
    left, top = 145, 112
    cell_w, cell_h = 120, 58
    body: List[str] = [
        text(32, 38, title, 22, weight="700"),
        text(32, 62, subtitle, 12, "#536076"),
    ]
    for j, method in enumerate(METHOD_ORDER):
        x = left + j * cell_w + cell_w / 2
        body.append(text(x, top - 18, method.upper(), 13, "#172033", "middle", "700"))
    for i, family in enumerate(FAMILY_ORDER):
        y = top + i * cell_h + cell_h / 2
        body.append(text(left - 18, y + 5, family, 13, "#172033", "end", "700"))
        for j, method in enumerate(METHOD_ORDER):
            x0 = left + j * cell_w
            y0 = top + i * cell_h
            value = lookup.get((family, method), math.nan)
            color = "#eef1f6" if not math.isfinite(value) else interpolate_color(
                value, (238, 244, 249), (35, 112, 174)
            )
            body.append(rect(x0, y0, cell_w - 2, cell_h - 2, color, "#ffffff", 1, 4))
            label = "NA" if not math.isfinite(value) else f"{value:.2f}"
            body.append(text(x0 + cell_w / 2, y0 + cell_h / 2 + 5, label, 17, "#111827", "middle", "700"))
    body.extend(draw_colorbar(610, top, 26, cell_h * len(FAMILY_ORDER) - 2, "0", "1"))
    write_svg(out_path, width, height, body)
    return out_path


def aggregate_rows_by_method_family(rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    grouped: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[(row["event_family"], row["method"])].append(row)

    out: List[Dict[str, str]] = []
    for (family, method), group in sorted(grouped.items()):
        n = len(group)
        if n == 0:
            continue
        all_active = sum(safe_int(row.get("all_family_events_active", "0")) for row in group)
        selective = sum(safe_int(row.get("selective_against_forbidden_events", "0")) for row in group)
        all_active_selective = sum(safe_int(row.get("all_active_and_selective", "0")) for row in group)
        similar = sum(safe_int(row.get("amplitude_similar_by_tolerance", "0")) for row in group)
        selective_similar = sum(safe_int(row.get("all_active_selective_and_similar", "0")) for row in group)
        forbidden_any = sum(1 for row in group if safe_int(row.get("forbidden_active_event_count", "0")) > 0)
        out.append(
            {
                "event_family": family,
                "method": method,
                "run_family_count": str(n),
                "best_component_all_active_count": str(all_active),
                "best_component_all_active_fraction": f"{all_active / n:.12g}",
                "best_component_selective_count": str(selective),
                "best_component_selective_fraction": f"{selective / n:.12g}",
                "best_component_all_active_selective_count": str(all_active_selective),
                "best_component_all_active_selective_fraction": f"{all_active_selective / n:.12g}",
                "best_component_forbidden_active_count": str(forbidden_any),
                "best_component_forbidden_active_fraction": f"{forbidden_any / n:.12g}",
                "best_component_similar_amplitude_count": str(similar),
                "best_component_similar_amplitude_fraction": f"{similar / n:.12g}",
                "best_component_selective_similar_amplitude_count": str(selective_similar),
                "best_component_selective_similar_amplitude_fraction": f"{selective_similar / n:.12g}",
            }
        )
    return out


def plot_dataset_selective_heatmaps(best_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    svg_paths: List[Path] = []
    for dataset in sorted({row["dataset"] for row in best_rows}):
        dataset_rows = [row for row in best_rows if row["dataset"] == dataset]
        summary_rows = aggregate_rows_by_method_family(dataset_rows)
        per_cell_denominators = {
            len(
                {
                    (row["condition"], row["component_count"])
                    for row in dataset_rows
                    if row["method"] == method and row["event_family"] == family
                }
            )
            for method in METHOD_ORDER
            for family in FAMILY_ORDER
        }
        denominator_text = (
            str(per_cell_denominators.pop())
            if len(per_cell_denominators) == 1
            else "variable"
        )
        svg_paths.append(
            plot_heatmap(
                summary_rows,
                "best_component_all_active_selective_fraction",
                f"{dataset}: All-Active Components That Stay Selective",
                output_dir / f"family_method_selective_heatmap_{dataset}.svg",
                f"Fraction across this dataset only; each cell denominator: {denominator_text} runs.",
            )
        )
    return svg_paths


def condition_short_name(condition: str) -> str:
    if condition.startswith("abs_"):
        return "abs"
    if condition.startswith("complex_split_"):
        return "complex"
    return condition.replace("_projected_vlambda", "")


def plot_cross_session_heatmap(
    best_rows: Sequence[Dict[str, str]],
    condition: str,
    family: str,
    value_key: str,
    title: str,
    out_path: Path,
) -> Path:
    rows = [
        row
        for row in best_rows
        if row["condition"] == condition and row["event_family"] == family
    ]
    width, height = 760, 520
    left, top = 145, 118
    cell_w, cell_h = 120, 52
    body: List[str] = [
        text(32, 38, title, 22, weight="700"),
        text(
            32,
            62,
            "Rows are component counts; each cell is datasets passing the criterion.",
            12,
            "#536076",
        ),
    ]
    for j, method in enumerate(METHOD_ORDER):
        x = left + j * cell_w + cell_w / 2
        body.append(text(x, top - 18, method.upper(), 13, "#172033", "middle", "700"))

    for i, k in enumerate(range(3, 9)):
        y = top + i * cell_h + cell_h / 2
        body.append(text(left - 18, y + 5, f"k{k:02d}", 13, "#172033", "end", "700"))
        for j, method in enumerate(METHOD_ORDER):
            cell_rows = [
                row
                for row in rows
                if row["method"] == method and safe_int(row["component_count"]) == k
            ]
            denom = len({row["dataset"] for row in cell_rows})
            count = sum(safe_int(row.get(value_key, "0")) for row in cell_rows)
            fraction = count / denom if denom else math.nan
            x0 = left + j * cell_w
            y0 = top + i * cell_h
            color = "#eef1f6" if not math.isfinite(fraction) else interpolate_color(
                fraction, (238, 244, 249), (35, 112, 174)
            )
            body.append(rect(x0, y0, cell_w - 2, cell_h - 2, color, "#ffffff", 1, 4))
            label = "NA" if denom == 0 else f"{count}/{denom}"
            body.append(text(x0 + cell_w / 2, y0 + cell_h / 2 + 5, label, 17, "#111827", "middle", "700"))
    body.extend(draw_colorbar(610, top, 26, cell_h * 6 - 2, "0", "1"))
    write_svg(out_path, width, height, body)
    return out_path


def plot_cross_session_joint_family_heatmap(
    best_rows: Sequence[Dict[str, str]],
    condition: str,
    families: Sequence[str],
    value_key: str,
    title: str,
    out_path: Path,
) -> Path:
    rows = [row for row in best_rows if row["condition"] == condition and row["event_family"] in families]
    family_set = set(families)
    width, height = 760, 520
    left, top = 145, 118
    cell_w, cell_h = 120, 52
    family_label = "+".join(families)
    body: List[str] = [
        text(32, 38, title, 22, weight="700"),
        text(
            32,
            62,
            f"Rows are component counts; each cell is datasets where all families pass: {family_label}.",
            12,
            "#536076",
        ),
    ]
    for j, method in enumerate(METHOD_ORDER):
        x = left + j * cell_w + cell_w / 2
        body.append(text(x, top - 18, method.upper(), 13, "#172033", "middle", "700"))

    for i, k in enumerate(range(3, 9)):
        y = top + i * cell_h + cell_h / 2
        body.append(text(left - 18, y + 5, f"k{k:02d}", 13, "#172033", "end", "700"))
        for j, method in enumerate(METHOD_ORDER):
            cell_rows = [
                row
                for row in rows
                if row["method"] == method and safe_int(row["component_count"]) == k
            ]
            by_dataset: Dict[str, Dict[str, Dict[str, str]]] = defaultdict(dict)
            for row in cell_rows:
                by_dataset[row["dataset"]][row["event_family"]] = row
            complete_dataset_rows = [
                dataset_rows
                for dataset_rows in by_dataset.values()
                if family_set.issubset(dataset_rows.keys())
            ]
            denom = len(complete_dataset_rows)
            count = sum(
                1
                for dataset_rows in complete_dataset_rows
                if all(safe_int(dataset_rows[family].get(value_key, "0")) for family in families)
            )
            fraction = count / denom if denom else math.nan
            x0 = left + j * cell_w
            y0 = top + i * cell_h
            color = "#eef1f6" if not math.isfinite(fraction) else interpolate_color(
                fraction, (238, 244, 249), (35, 112, 174)
            )
            body.append(rect(x0, y0, cell_w - 2, cell_h - 2, color, "#ffffff", 1, 4))
            label = "NA" if denom == 0 else f"{count}/{denom}"
            body.append(text(x0 + cell_w / 2, y0 + cell_h / 2 + 5, label, 17, "#111827", "middle", "700"))
    body.extend(draw_colorbar(610, top, 26, cell_h * 6 - 2, "0", "1"))
    write_svg(out_path, width, height, body)
    return out_path


def plot_cross_session_breakdown_heatmaps(best_rows: Sequence[Dict[str, str]], output_dir: Path) -> List[Path]:
    svg_paths: List[Path] = []
    conditions = sorted({row["condition"] for row in best_rows})
    metrics = [
        (
            "all_active_and_selective",
            "cross_session_selective_by_k",
            "Cross-Session Selectivity",
        ),
        (
            "all_active_selective_and_similar",
            "cross_session_selective_similar_by_k",
            "Cross-Session Selective + Similar Amplitude",
        ),
    ]
    for condition in conditions:
        condition_short = condition_short_name(condition)
        for family in FAMILY_ORDER:
            for value_key, prefix, title_prefix in metrics:
                svg_paths.append(
                    plot_cross_session_heatmap(
                        best_rows,
                        condition,
                        family,
                        value_key,
                        f"{title_prefix}: {family} | {condition_short}",
                        output_dir / f"{prefix}_{condition_short}_{family}.svg",
                    )
                )
        for value_key, prefix, title_prefix in metrics:
            svg_paths.append(
                plot_cross_session_joint_family_heatmap(
                    best_rows,
                    condition,
                    ("theta", "ripple"),
                    value_key,
                    f"{title_prefix}: theta+ripple | {condition_short}",
                    output_dir / f"{prefix}_{condition_short}_theta_ripple_joint.svg",
                )
            )
    return svg_paths


def row_status(row: Optional[Dict[str, str]]) -> Tuple[str, str]:
    if row is None:
        return "inactive", ""
    component_idx = row.get("component_idx", "")
    if safe_int(row.get("all_active_and_selective", "0")):
        if safe_int(row.get("all_active_selective_and_similar", "0")):
            return "similar", component_idx
        return "all_active", component_idx
    if safe_int(row.get("all_family_events_active", "0")):
        return "nonselective", component_idx
    return "inactive", component_idx


def plot_dataset_condition_event_selectivity(
    best_rows: Sequence[Dict[str, str]],
    dataset: str,
    condition: str,
    out_path: Path,
) -> Path:
    rows = [
        row
        for row in best_rows
        if row["dataset"] == dataset and row["condition"] == condition
    ]
    lookup = {
        (row["event_family"], row["method"], safe_int(row["component_count"])): row
        for row in rows
    }

    width, height = 1120, 520
    top = 122
    panel_w, panel_h = 286, 312
    left0 = 78
    gap = 58
    cell_w, cell_h = 58, 42
    condition_short = condition_short_name(condition)
    body: List[str] = [
        text(32, 38, f"{dataset} | {condition_short}: Event-Family Selectivity By k", 22, weight="700"),
        text(
            32,
            62,
            "Gray: not all active; orange: all active but nonselective; blue: selective; green: selective + similar amplitude.",
            12,
            "#536076",
        ),
    ]

    for p, family in enumerate(FAMILY_ORDER):
        x_panel = left0 + p * (panel_w + gap)
        body.append(text(x_panel + panel_w / 2, top - 26, family, 17, "#172033", "middle", "700"))
        body.append(rect(x_panel, top, panel_w, panel_h, "#f8fafc", "#cad2df", 1, 4))
        grid_left = x_panel + 48
        grid_top = top + 38
        for j, method in enumerate(METHOD_ORDER):
            body.append(text(grid_left + j * cell_w + cell_w / 2, grid_top - 12, method.upper(), 10, "#253044", "middle", "700"))
        for i, k in enumerate(range(3, 9)):
            y = grid_top + i * cell_h
            body.append(text(grid_left - 12, y + cell_h / 2 + 4, f"k{k:02d}", 10, "#253044", "end", "700"))
            for j, method in enumerate(METHOD_ORDER):
                x = grid_left + j * cell_w
                row = lookup.get((family, method, k))
                status, comp = row_status(row)
                body.append(rect(x, y, cell_w - 4, cell_h - 4, STATUS_COLORS[status], "#ffffff", 0.8, 3))
                if comp:
                    body.append(text(x + cell_w / 2 - 2, y + cell_h / 2 + 4, f"c{comp}", 10, "#111827", "middle", "700"))

    legend_y = height - 46
    legend_items = [
        ("inactive", "not all active"),
        ("nonselective", "all active, nonselective"),
        ("all_active", "selective"),
        ("similar", "selective + similar"),
    ]
    x = 32
    for status, label in legend_items:
        body.append(rect(x, legend_y, 18, 12, STATUS_COLORS[status], "#ffffff", 0.5, 2))
        body.append(text(x + 24, legend_y + 10, label, 11, "#253044"))
        x += 190
    write_svg(out_path, width, height, body)
    return out_path


def plot_dataset_condition_event_selectivity_heatmaps(
    best_rows: Sequence[Dict[str, str]],
    output_dir: Path,
) -> List[Path]:
    svg_paths: List[Path] = []
    datasets = sorted({row["dataset"] for row in best_rows})
    conditions = sorted({row["condition"] for row in best_rows})
    for dataset in datasets:
        for condition in conditions:
            condition_short = condition_short_name(condition)
            svg_paths.append(
                plot_dataset_condition_event_selectivity(
                    best_rows,
                    dataset,
                    condition,
                    output_dir / f"dataset_condition_event_selectivity_by_k_{dataset}_{condition_short}.svg",
                )
            )
    return svg_paths


def draw_colorbar(x: float, y: float, w: float, h: float, low_label: str, high_label: str) -> List[str]:
    body: List[str] = []
    steps = 40
    for i in range(steps):
        t0 = i / steps
        yy = y + h - (i + 1) * h / steps
        body.append(
            rect(
                x,
                yy,
                w,
                h / steps + 0.5,
                interpolate_color(t0, (238, 244, 249), (35, 112, 174)),
                "none",
                0,
            )
        )
    body.append(rect(x, y, w, h, "none", "#7b8798", 0.8))
    body.append(text(x + w + 8, y + h, low_label, 10, "#536076"))
    body.append(text(x + w + 8, y + 8, high_label, 10, "#536076"))
    return body


def plot_status_tilemap(best_rows: Sequence[Dict[str, str]], out_path: Path) -> Path:
    grouped: Dict[Tuple[str, str, str, int], Dict[str, Dict[str, str]]] = defaultdict(dict)
    for row in best_rows:
        key = (row["dataset"], row["condition"], row["method"], safe_int(row["component_count"]))
        grouped[key][row["event_family"]] = row

    dataset_order = sorted({row["dataset"] for row in best_rows})
    condition_order = sorted({row["condition"] for row in best_rows})
    row_keys = [
        (dataset, condition, method, k)
        for dataset in dataset_order
        for condition in condition_order
        for method in METHOD_ORDER
        for k in range(3, 9)
        if (dataset, condition, method, k) in grouped
    ]

    width = 1120
    row_h = 10
    top = 92
    left = 350
    tile_w = 78
    height = top + row_h * len(row_keys) + 82
    body: List[str] = [
        text(32, 36, "Run-Level Event-Family Component Status", 22, weight="700"),
        text(32, 60, "Gray: not all active; orange: all active but nonselective; blue: selective; green: selective + similar amplitude.", 12, "#536076"),
    ]
    for j, family in enumerate(FAMILY_ORDER):
        body.append(text(left + j * tile_w + tile_w / 2, top - 18, family, 12, "#172033", "middle", "700"))

    previous_group = None
    for i, key in enumerate(row_keys):
        dataset, condition, method, k = key
        y = top + i * row_h
        group = (dataset, condition)
        if group != previous_group:
            body.append(line(24, y - 2, width - 42, y - 2, "#c7ced8", 1))
            previous_group = group
        label = f"{dataset} | {condition.replace('_projected_vlambda', '')} | {method.upper()} k{k:02d}"
        body.append(text(32, y + 8, label, 7, "#253044"))
        for j, family in enumerate(FAMILY_ORDER):
            row = grouped[key].get(family)
            if row is None:
                status = "inactive"
                comp = ""
            elif not safe_int(row.get("all_active_and_selective", "0")):
                status = "nonselective"
                comp = row["component_idx"]
            elif safe_int(row.get("all_active_selective_and_similar", "0")):
                status = "similar"
                comp = row["component_idx"]
            elif safe_int(row.get("all_active_and_selective", "0")):
                status = "all_active"
                comp = row["component_idx"]
            else:
                status = "inactive"
                comp = row["component_idx"]
            x = left + j * tile_w
            body.append(rect(x, y, tile_w - 4, row_h - 1, STATUS_COLORS[status], "#ffffff", 0.5, 1))
            if comp:
                body.append(text(x + tile_w / 2 - 2, y + 7.3, comp, 6, "#111827", "middle", "700"))

    legend_y = height - 44
    legend_items = [
        ("inactive", "not all active"),
        ("nonselective", "all active, nonselective"),
        ("all_active", "all active + selective"),
        ("similar", "selective + similar"),
    ]
    x = 32
    for status, label in legend_items:
        body.append(rect(x, legend_y, 18, 12, STATUS_COLORS[status], "#ffffff", 0.5, 2))
        body.append(text(x + 24, legend_y + 10, label, 11, "#253044"))
        x += 190
    write_svg(out_path, width, height, body)
    return out_path


def load_stats_lookup(stats_file: str) -> Dict[Tuple[str, int], Dict[str, str]]:
    rows = read_csv(Path(stats_file))
    return {
        (row["state_label"], safe_int(row["component_idx"])): row
        for row in rows
    }


def equality_points(best_rows: Sequence[Dict[str, str]]) -> List[Dict[str, object]]:
    cache: Dict[str, Dict[Tuple[str, int], Dict[str, str]]] = {}
    points: List[Dict[str, object]] = []
    for row in best_rows:
        if not safe_int(row.get("all_active_and_selective", "0")):
            continue
        events = [item for item in row["family_events"].split(";") if item]
        if len(events) != 2:
            continue
        stats_file = row["stats_file"]
        if stats_file not in cache:
            cache[stats_file] = load_stats_lookup(stats_file)
        comp_idx = safe_int(row["component_idx"])
        amp_col = row["amplitude_column"]
        left_row = cache[stats_file].get((events[0], comp_idx))
        right_row = cache[stats_file].get((events[1], comp_idx))
        if not left_row or not right_row:
            continue
        x = safe_float(left_row.get(amp_col))
        y = safe_float(right_row.get(amp_col))
        if not (math.isfinite(x) and math.isfinite(y)):
            continue
        point = dict(row)
        point["x_event"] = events[0]
        point["y_event"] = events[1]
        point["x"] = x
        point["y"] = y
        points.append(point)
    return points


def nice_max(values: Sequence[float]) -> float:
    max_value = max([value for value in values if math.isfinite(value)] + [1.0])
    if max_value <= 0:
        return 1.0
    exponent = math.floor(math.log10(max_value))
    base = 10 ** exponent
    for mult in (1, 2, 5, 10):
        candidate = mult * base
        if candidate >= max_value:
            return candidate
    return 10 * base


def plot_equality_scatter(best_rows: Sequence[Dict[str, str]], out_path: Path) -> Path:
    pts = equality_points(best_rows)
    by_family: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for pt in pts:
        by_family[str(pt["event_family"])].append(pt)

    width, height = 1180, 480
    panel_w, panel_h = 300, 260
    left0, top = 80, 92
    gap = 62
    body: List[str] = [
        text(32, 36, "Amplitude Equality For Best All-Active Components", 22, weight="700"),
        text(32, 60, "Each point is one run-family best component. Points near y=x have similar event-family amplitudes.", 12, "#536076"),
    ]
    for i, family in enumerate(FAMILY_ORDER):
        panel_pts = by_family.get(family, [])
        x0 = left0 + i * (panel_w + gap)
        y0 = top
        max_axis = nice_max([float(pt["x"]) for pt in panel_pts] + [float(pt["y"]) for pt in panel_pts])
        body.append(text(x0 + panel_w / 2, y0 - 22, family, 16, "#172033", "middle", "700"))
        body.append(rect(x0, y0, panel_w, panel_h, "#f8fafc", "#cad2df", 1))
        for t in (0, 0.25, 0.5, 0.75, 1.0):
            gx = x0 + t * panel_w
            gy = y0 + panel_h - t * panel_h
            body.append(line(gx, y0, gx, y0 + panel_h, "#e3e8ef", 0.8))
            body.append(line(x0, gy, x0 + panel_w, gy, "#e3e8ef", 0.8))
        body.append(line(x0, y0 + panel_h, x0 + panel_w, y0, "#667085", 1.2, "5 4"))
        for pt in panel_pts:
            xx = x0 + float(pt["x"]) / max_axis * panel_w
            yy = y0 + panel_h - float(pt["y"]) / max_axis * panel_h
            method = str(pt["method"])
            similar = safe_int(pt["amplitude_similar_by_tolerance"])
            body.append(circle(xx, yy, 4.1 if similar else 3.2, METHOD_COLORS.get(method, "#666666"), opacity=0.78 if similar else 0.42))
        body.append(text(x0, y0 + panel_h + 22, "0", 10, "#536076", "middle"))
        body.append(text(x0 + panel_w, y0 + panel_h + 22, f"{max_axis:g}", 10, "#536076", "middle"))
        body.append(text(x0 + panel_w / 2, y0 + panel_h + 42, "first family event amplitude", 11, "#536076", "middle"))
        body.append(text(x0 - 42, y0 + panel_h / 2, "second event amplitude", 11, "#536076", "middle", rotate=-90))

    legend_y = 438
    x = 80
    for method in METHOD_ORDER:
        body.append(circle(x, legend_y, 5, METHOD_COLORS[method]))
        body.append(text(x + 12, legend_y + 4, method.upper(), 11, "#253044"))
        x += 82
    write_svg(out_path, width, height, body)
    return out_path


def quantiles(values: Sequence[float]) -> Optional[Tuple[float, float, float, float, float]]:
    clean = sorted(value for value in values if math.isfinite(value))
    if not clean:
        return None
    return (
        clean[0],
        percentile(clean, 0.25),
        median(clean),
        percentile(clean, 0.75),
        clean[-1],
    )


def percentile(clean_sorted: Sequence[float], p: float) -> float:
    if len(clean_sorted) == 1:
        return clean_sorted[0]
    pos = p * (len(clean_sorted) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return clean_sorted[lo]
    frac = pos - lo
    return clean_sorted[lo] * (1 - frac) + clean_sorted[hi] * frac


def plot_cv_box(best_rows: Sequence[Dict[str, str]], out_path: Path) -> Path:
    values: Dict[Tuple[str, str], List[float]] = defaultdict(list)
    for row in best_rows:
        if safe_int(row.get("all_active_and_selective", "0")):
            values[(row["event_family"], row["method"])].append(safe_float(row["amplitude_cv"]))

    width, height = 980, 470
    left, top = 90, 92
    panel_w, panel_h = 240, 260
    gap = 52
    ymax = 0.8
    body: List[str] = [
        text(32, 36, "Amplitude CV By Method And Event Family", 22, weight="700"),
        text(32, 60, "Lower CV means more equal amplitudes across events in the same family. All-active best components only.", 12, "#536076"),
    ]
    for i, family in enumerate(FAMILY_ORDER):
        x0 = left + i * (panel_w + gap)
        y0 = top
        body.append(text(x0 + panel_w / 2, y0 - 22, family, 16, "#172033", "middle", "700"))
        body.append(rect(x0, y0, panel_w, panel_h, "#f8fafc", "#cad2df", 1))
        for tick in (0, 0.2, 0.4, 0.6, 0.8):
            yy = y0 + panel_h - tick / ymax * panel_h
            body.append(line(x0, yy, x0 + panel_w, yy, "#e3e8ef", 0.8))
            body.append(text(x0 - 8, yy + 4, f"{tick:.1f}", 9, "#536076", "end"))
        for j, method in enumerate(METHOD_ORDER):
            q = quantiles(values.get((family, method), []))
            cx = x0 + 34 + j * 52
            body.append(text(cx, y0 + panel_h + 22, method.upper(), 9, "#253044", "middle"))
            if q is None:
                continue
            vmin, q1, med, q3, vmax = q
            y_min = y0 + panel_h - min(vmin, ymax) / ymax * panel_h
            y_q1 = y0 + panel_h - min(q1, ymax) / ymax * panel_h
            y_med = y0 + panel_h - min(med, ymax) / ymax * panel_h
            y_q3 = y0 + panel_h - min(q3, ymax) / ymax * panel_h
            y_max = y0 + panel_h - min(vmax, ymax) / ymax * panel_h
            color = METHOD_COLORS[method]
            body.append(line(cx, y_max, cx, y_min, color, 1.3))
            body.append(line(cx - 8, y_max, cx + 8, y_max, color, 1.3))
            body.append(line(cx - 8, y_min, cx + 8, y_min, color, 1.3))
            body.append(rect(cx - 14, y_q3, 28, max(3, y_q1 - y_q3), "#ffffff", color, 1.4, 2))
            body.append(line(cx - 14, y_med, cx + 14, y_med, color, 2))
    write_svg(out_path, width, height, body)
    return out_path


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = (args.output_dir or (input_dir / "figures")).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = read_csv(input_dir / "method_family_summary.csv")
    best_rows = read_csv(input_dir / "best_component_by_run_family.csv")

    svg_paths = [
        plot_heatmap(
            summary_rows,
            "best_component_all_active_fraction",
            "Single Component Active Across Each Event Family",
            output_dir / "family_method_all_active_heatmap.svg",
        ),
        plot_heatmap(
            summary_rows,
            "best_component_all_active_selective_fraction",
            "All-Active Components That Stay Selective",
            output_dir / "family_method_selective_heatmap.svg",
        ),
        plot_heatmap(
            summary_rows,
            "best_component_selective_similar_amplitude_fraction",
            "Selective All-Active Components With Similar Amplitudes",
            output_dir / "family_method_similar_amplitude_heatmap.svg",
        ),
        plot_status_tilemap(best_rows, output_dir / "run_family_status_tilemap.svg"),
        plot_equality_scatter(best_rows, output_dir / "amplitude_equality_scatter_by_family.svg"),
        plot_cv_box(best_rows, output_dir / "amplitude_cv_box_by_method_family.svg"),
    ]
    svg_paths.extend(plot_dataset_selective_heatmaps(best_rows, output_dir))
    svg_paths.extend(plot_cross_session_breakdown_heatmaps(best_rows, output_dir))
    svg_paths.extend(plot_dataset_condition_event_selectivity_heatmaps(best_rows, output_dir))

    png_paths: List[Path] = []
    inkscape = None if args.no_png else find_inkscape()
    if not args.no_png:
        for svg_path in svg_paths:
            png_path = svg_path.with_suffix(".png")
            if export_png(svg_path, png_path, inkscape):
                png_paths.append(png_path)

    print(f"SVG files written: {len(svg_paths)}")
    print(f"PNG files written: {len(png_paths)}")
    print(f"Output dir       : {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
