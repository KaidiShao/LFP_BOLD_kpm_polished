"""Summarize current P8 cross-session density/feature consistency.

This script reads canonical P8 xcorr CSV outputs:

    <processed-root>/<dataset>/pipeline8_xcorr/<run_tag>/

It produces CSV summaries plus simple PNG heatmaps.  The default scope is the
current four-session set used while E10gW1 is still excluded from downstream
backfill.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median, pstdev
from typing import Iterable, Sequence

try:
    from PIL import Image, ImageDraw, ImageFont
except Exception:  # pragma: no cover - SVG fallback keeps the script usable.
    Image = None
    ImageDraw = None
    ImageFont = None


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "f12m01")
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")
RUN_TAG_LABELS = {
    "pv_gsvd100": "global_svd100",
    "pv_gsvd100_ds": "gsvd100_ds",
    "pv_hp100": "HP_svd100",
    "pv_roi": "roi_mean",
}
DEFAULT_RESULTS_DIR = Path("results") / "pipeline8_cross_session_consistency_current"
DEFAULT_RESULTS_FIGURE_DIR = DEFAULT_RESULTS_DIR / "figures"
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p8"
    / "cross_session_consistency"
)
METHOD_ORDER = ("svd", "nmf", "mds", "umap")
COMPONENT_COUNTS = tuple(range(3, 9))


@dataclass(frozen=True)
class ScoreRow:
    dataset: str
    run_tag: str
    observable: str
    density_name: str
    bold_feature: str
    source_level: str
    n_top_rows: int
    mean_peak_abs_corr: float
    max_peak_abs_corr: float
    mean_peak_corr: float
    mean_peak_lag_sec: float
    source_csv: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--xcorr-save-tag", default="xcorr")
    parser.add_argument(
        "--results-figure-dir",
        type=Path,
        default=DEFAULT_RESULTS_FIGURE_DIR,
        help="Repo-local derived figure source for P11 flat-copy manifests.",
    )
    parser.add_argument("--skip-figures", action="store_true")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def as_float(value: str | None) -> float:
    if value is None:
        return math.nan
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def finite(values: Iterable[float]) -> list[float]:
    return [v for v in values if math.isfinite(v)]


def density_display_name(density_name: str) -> str:
    meta = density_metadata(density_name)
    if meta["density_source_kind"] == "event":
        return "event_density"
    if meta["density_source_kind"] == "raw":
        return f"raw_{meta['density_condition']}_q{meta['density_threshold']}"
    if meta["density_source_kind"] == "dimred":
        return (
            f"dim_{meta['density_condition']}_"
            f"{meta['density_method']}_k{int(meta['density_k']):02d}_"
            f"q{meta['density_threshold']}"
        )
    return density_name


def density_metadata(density_name: str) -> dict[str, object]:
    if density_name == "blp_evt":
        return {
            "density_source_kind": "event",
            "density_condition": "event",
            "density_method": "",
            "density_k": "",
            "density_method_k": "",
            "density_threshold": "",
        }
    raw_match = re.match(r"^raw_(abs|csplit)_q(\d+)(?:_(.+))?$", density_name)
    if raw_match:
        condition, threshold, activity_suffix = raw_match.groups()
        return {
            "density_source_kind": "raw",
            "density_condition": condition,
            "density_method": "raw",
            "density_k": "",
            "density_method_k": "raw",
            "density_threshold": threshold,
            "density_activity_suffix": activity_suffix or "",
        }
    dim_match = re.match(r"^dim_(abs|csplit)_([a-z]+)(\d+)_q(\d+)(?:_(.+))?$", density_name)
    if dim_match:
        condition, method, k, threshold, activity_suffix = dim_match.groups()
        k_int = int(k)
        return {
            "density_source_kind": "dimred",
            "density_condition": condition,
            "density_method": method,
            "density_k": k_int,
            "density_method_k": f"{method}_k{k_int:02d}",
            "density_threshold": threshold,
            "density_activity_suffix": activity_suffix or "",
        }
    return {
        "density_source_kind": "other",
        "density_condition": "",
        "density_method": "",
        "density_k": "",
        "density_method_k": "",
        "density_threshold": "",
        "density_activity_suffix": "",
    }


def density_group_name(density_name: str) -> str:
    if density_name.startswith("dim_"):
        return "dimred_efun_density"
    if density_name.startswith("raw_"):
        return "raw_efun_density"
    if density_name == "blp_evt":
        return "event_density"
    return "other_density"


def is_dimred_density(density_name: str) -> bool:
    return density_group_name(density_name) == "dimred_efun_density"


def density_metadata_fields(density_name: str) -> dict[str, object]:
    meta = density_metadata(density_name)
    return {
        "density_source_kind": meta["density_source_kind"],
        "density_condition": meta["density_condition"],
        "density_method": meta["density_method"],
        "density_k": meta["density_k"],
        "density_method_k": meta["density_method_k"],
        "density_threshold": meta["density_threshold"],
        "density_activity_suffix": meta.get("density_activity_suffix", ""),
    }


def summarize_top_rows(
    dataset: str,
    run_tag: str,
    level: str,
    csv_path: Path,
    rows: Sequence[dict[str, str]],
) -> list[ScoreRow]:
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        density = (row.get("density_name") or "unknown_density").strip()
        feature = (row.get("bold_feature") or "all_features").strip()
        grouped[(density, feature)].append(row)

    out: list[ScoreRow] = []
    for (density, feature), group_rows in sorted(grouped.items()):
        abs_vals = finite(as_float(r.get("peak_abs_corr")) for r in group_rows)
        corr_vals = finite(as_float(r.get("peak_corr")) for r in group_rows)
        lag_vals = finite(as_float(r.get("peak_lag_sec")) for r in group_rows)
        if not abs_vals:
            continue
        out.append(
            ScoreRow(
                dataset=dataset,
                run_tag=run_tag,
                observable=RUN_TAG_LABELS.get(run_tag, run_tag),
                density_name=density,
                bold_feature=feature,
                source_level=level,
                n_top_rows=len(group_rows),
                mean_peak_abs_corr=mean(abs_vals),
                max_peak_abs_corr=max(abs_vals),
                mean_peak_corr=mean(corr_vals) if corr_vals else math.nan,
                mean_peak_lag_sec=mean(lag_vals) if lag_vals else math.nan,
                source_csv=str(csv_path),
            )
        )
    return out


def make_hit_row(
    dataset: str,
    run_tag: str,
    source_level: str,
    csv_path: Path,
    row: dict[str, str],
    rank: int,
    top_n: int,
) -> dict[str, object]:
    density_name = (row.get("density_name") or "").strip()
    bold_feature = (row.get("bold_feature") or "").strip()
    density_index = row.get("density_index", "")
    density_label = row.get("density_label", "")
    bold_mode_index = row.get("bold_mode_index", "")
    peak_abs_corr = as_float(row.get("peak_abs_corr"))
    peak_corr = as_float(row.get("peak_corr"))
    peak_lag_sec = as_float(row.get("peak_lag_sec"))
    return {
        "dataset": dataset,
        "run_tag": run_tag,
        "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
        "source_level": source_level,
        "top_rank": rank,
        "top_n": top_n,
        "density_name": density_name,
        "density_display": density_display_name(density_name),
        "density_group": density_group_name(density_name),
        "is_dimred_efun_density": int(is_dimred_density(density_name)),
        **density_metadata_fields(density_name),
        "density_index": density_index,
        "density_label": density_label,
        "bold_feature": bold_feature,
        "feature_family": feature_family(bold_feature),
        "bold_mode_index": bold_mode_index,
        "peak_abs_corr": f"{peak_abs_corr:.10g}" if math.isfinite(peak_abs_corr) else "",
        "peak_corr": f"{peak_corr:.10g}" if math.isfinite(peak_corr) else "",
        "peak_lag_sec": f"{peak_lag_sec:.10g}" if math.isfinite(peak_lag_sec) else "",
        "source_csv": str(csv_path),
    }


def extend_hit_rows(
    hit_rows: list[dict[str, object]],
    dataset: str,
    run_tag: str,
    source_level: str,
    csv_path: Path,
    rows: Sequence[dict[str, str]],
) -> None:
    for i, row in enumerate(rows, start=1):
        hit_rows.append(make_hit_row(dataset, run_tag, source_level, csv_path, row, i, len(rows)))


def collect_scores(
    processed_root: Path,
    datasets: Sequence[str],
    run_tags: Sequence[str],
    xcorr_save_tag: str,
) -> tuple[list[ScoreRow], list[dict[str, str]], list[dict[str, object]]]:
    score_rows: list[ScoreRow] = []
    raw_top_rows: list[dict[str, str]] = []
    hit_rows: list[dict[str, object]] = []
    for dataset in datasets:
        for run_tag in run_tags:
            run_dir = processed_root / dataset / "pipeline8_xcorr" / run_tag
            combined_csv = run_dir / f"{xcorr_save_tag}_top.csv"
            combined_rows = read_csv_rows(combined_csv)
            extend_hit_rows(hit_rows, dataset, run_tag, "combined", combined_csv, combined_rows)
            for row in combined_rows:
                row_out = dict(row)
                density_name = (row.get("density_name") or "").strip()
                bold_feature = (row.get("bold_feature") or "").strip()
                row_out.update(
                    dataset=dataset,
                    run_tag=run_tag,
                    observable=RUN_TAG_LABELS.get(run_tag, run_tag),
                    source_level="combined",
                    density_display=density_display_name(density_name),
                    density_group=density_group_name(density_name),
                    is_dimred_efun_density=int(is_dimred_density(density_name)),
                    **density_metadata_fields(density_name),
                    feature_family=feature_family(bold_feature),
                    source_csv=str(combined_csv),
                )
                raw_top_rows.append(row_out)
            score_rows.extend(summarize_top_rows(dataset, run_tag, "combined", combined_csv, combined_rows))

            density_dir = run_dir / "density"
            for csv_path in sorted(density_dir.glob(f"{xcorr_save_tag}_top__*.csv")):
                rows = read_csv_rows(csv_path)
                extend_hit_rows(hit_rows, dataset, run_tag, "per_density", csv_path, rows)
                score_rows.extend(summarize_top_rows(dataset, run_tag, "per_density", csv_path, rows))

            feature_dir = run_dir / "feature"
            for csv_path in sorted(feature_dir.glob(f"*/{xcorr_save_tag}_top__*.csv")):
                rows = read_csv_rows(csv_path)
                extend_hit_rows(hit_rows, dataset, run_tag, "per_density_feature", csv_path, rows)
                score_rows.extend(summarize_top_rows(dataset, run_tag, "per_density_feature", csv_path, rows))
    return score_rows, raw_top_rows, hit_rows


def row_to_dict(row: ScoreRow) -> dict[str, object]:
    return {
        "dataset": row.dataset,
        "run_tag": row.run_tag,
        "observable": row.observable,
        "density_name": row.density_name,
        "density_display": density_display_name(row.density_name),
        "density_group": density_group_name(row.density_name),
        "is_dimred_efun_density": int(is_dimred_density(row.density_name)),
        **density_metadata_fields(row.density_name),
        "bold_feature": row.bold_feature,
        "feature_family": feature_family(row.bold_feature),
        "source_level": row.source_level,
        "n_top_rows": row.n_top_rows,
        "mean_peak_abs_corr": f"{row.mean_peak_abs_corr:.10g}",
        "max_peak_abs_corr": f"{row.max_peak_abs_corr:.10g}",
        "mean_peak_corr": f"{row.mean_peak_corr:.10g}" if math.isfinite(row.mean_peak_corr) else "",
        "mean_peak_lag_sec": f"{row.mean_peak_lag_sec:.10g}" if math.isfinite(row.mean_peak_lag_sec) else "",
        "source_csv": row.source_csv,
    }


def write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys: list[str] = []
        seen = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    seen.add(key)
                    keys.append(key)
        fieldnames = keys
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def summarize_cross_session(score_rows: Sequence[ScoreRow], level: str) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[ScoreRow]] = defaultdict(list)
    for row in score_rows:
        if row.source_level != level:
            continue
        grouped[(row.run_tag, row.density_name, row.bold_feature)].append(row)

    out: list[dict[str, object]] = []
    for (run_tag, density, feature), rows in sorted(grouped.items()):
        values = [r.mean_peak_abs_corr for r in rows if math.isfinite(r.mean_peak_abs_corr)]
        datasets = sorted({r.dataset for r in rows})
        if not values:
            continue
        out.append(
            {
                "run_tag": run_tag,
                "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                "density_name": density,
                "density_display": density_display_name(density),
                **density_metadata_fields(density),
                "bold_feature": feature,
                "source_level": level,
                "n_datasets": len(datasets),
                "datasets": ";".join(datasets),
                "mean_peak_abs_corr_across_datasets": f"{mean(values):.10g}",
                "std_peak_abs_corr_across_datasets": f"{pstdev(values):.10g}" if len(values) > 1 else "0",
                "min_peak_abs_corr": f"{min(values):.10g}",
                "max_peak_abs_corr": f"{max(values):.10g}",
                "range_peak_abs_corr": f"{(max(values) - min(values)):.10g}",
            }
        )
    return out


def dominant_rows(score_rows: Sequence[ScoreRow], level: str) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[ScoreRow]] = defaultdict(list)
    for row in score_rows:
        if row.source_level == level:
            grouped[(row.dataset, row.run_tag)].append(row)

    out: list[dict[str, object]] = []
    for (dataset, run_tag), rows in sorted(grouped.items()):
        if not rows:
            continue
        best = max(rows, key=lambda r: r.mean_peak_abs_corr)
        out.append(
            {
                "dataset": dataset,
                "run_tag": run_tag,
                "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                "source_level": level,
                "dominant_density_name": best.density_name,
                "dominant_density_display": density_display_name(best.density_name),
                **{
                    f"dominant_{key}": value
                    for key, value in density_metadata_fields(best.density_name).items()
                },
                "dominant_bold_feature": best.bold_feature,
                "dominant_density_feature": f"{best.density_name} | {best.bold_feature}",
                "mean_peak_abs_corr": f"{best.mean_peak_abs_corr:.10g}",
                "max_peak_abs_corr": f"{best.max_peak_abs_corr:.10g}",
            }
        )
    return out


def categorical_agreement(dominant: Sequence[dict[str, object]], category_field: str) -> list[dict[str, object]]:
    grouped: dict[str, Counter[str]] = defaultdict(Counter)
    for row in dominant:
        grouped[str(row["run_tag"])][str(row[category_field])] += 1
    out: list[dict[str, object]] = []
    for run_tag, counts in sorted(grouped.items()):
        total = sum(counts.values())
        for category, count in counts.most_common():
            out.append(
                {
                    "run_tag": run_tag,
                    "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                    "category_field": category_field,
                    "category": category,
                    "n_datasets": count,
                    "fraction_datasets": f"{(count / total if total else math.nan):.10g}",
                }
            )
    return out


def hit_value(row: dict[str, object], key: str = "peak_abs_corr") -> float:
    return as_float(str(row.get(key, "")))


def best_hit_rows(
    hit_rows: Sequence[dict[str, object]],
    *,
    dimred_only: bool = False,
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in hit_rows:
        if row.get("source_level") != "per_density_feature":
            continue
        family = str(row.get("feature_family", ""))
        if family not in {"efun", "deconv_efun"}:
            continue
        if dimred_only and str(row.get("density_group", "")) != "dimred_efun_density":
            continue
        key = (str(row.get("dataset")), str(row.get("run_tag")), family)
        grouped[key].append(dict(row))

    out: list[dict[str, object]] = []
    for key, rows in sorted(grouped.items()):
        valid = [r for r in rows if math.isfinite(hit_value(r))]
        if not valid:
            continue
        best = max(valid, key=hit_value)
        out.append(best)
    return out


def short_density_for_cell(row: dict[str, object]) -> str:
    density = str(row.get("density_display") or row.get("density_name") or "")
    density = density.replace("_density", "")
    density = density.replace("_q070", "")
    idx = str(row.get("density_index") or "")
    if idx:
        return f"{density} d{idx}"
    return density


def short_feature_for_cell(row: dict[str, object]) -> str:
    feature = str(row.get("bold_feature") or "")
    mode = str(row.get("bold_mode_index") or "")
    if mode:
        return f"{feature} m{mode}"
    return feature


def plot_best_hit_table(
    rows: Sequence[dict[str, object]],
    figure_dir: Path,
    datasets: Sequence[str],
    run_tags: Sequence[str],
    *,
    family: str,
    dimred_only: bool,
) -> Path | None:
    rows = [r for r in rows if str(r.get("feature_family")) == family]
    if not rows:
        return None

    by_key = {(str(r.get("dataset")), str(r.get("run_tag"))): r for r in rows}
    row_labels = list(datasets)
    col_labels = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]
    cell_w = 260
    cell_h = 92
    left = 105
    top = 110
    width = left + 30 + cell_w * len(run_tags)
    height = top + 50 + cell_h * len(row_labels)
    if Image is None:
        return None
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(22, True)
    label_font = load_font(14, False)
    small_font = load_font(12, False)
    value_font = load_font(13, True)

    title = f"P8 top5 best hit | {family}"
    if dimred_only:
        title += " | dimred density only"
    draw_text(draw, (20, 22), title, title_font, fill=(15, 25, 35))
    draw_text(draw, (20, 52), "Each cell: strongest row among per-density-feature top5 tables. d=density component, m=BOLD mode.", small_font, fill=(80, 80, 80))

    for j, col in enumerate(col_labels):
        x = left + j * cell_w + cell_w / 2
        draw_text(draw, (x, top - 14), col, label_font, anchor="ms")
    for i, dataset in enumerate(row_labels):
        y = top + i * cell_h + cell_h / 2
        draw_text(draw, (left - 12, y), dataset, label_font, anchor="rm")
        for j, run_tag in enumerate(run_tags):
            x0 = left + j * cell_w
            y0 = top + i * cell_h
            row = by_key.get((dataset, run_tag))
            value = hit_value(row) if row else math.nan
            draw.rectangle([x0, y0, x0 + cell_w, y0 + cell_h], fill=heat_color(value, 0.0, 0.65 if family == "efun" else 0.35), outline=(205, 205, 205))
            if not row:
                draw_text(draw, (x0 + cell_w / 2, y0 + cell_h / 2), "missing", value_font, anchor="mm")
                continue
            lines = [
                short_density_for_cell(row),
                short_feature_for_cell(row),
                f"|r|={hit_value(row):.3f} lag={hit_value(row, 'peak_lag_sec'):.1f}s",
            ]
            for k, text in enumerate(lines):
                draw_text(draw, (x0 + 8, y0 + 18 + 22 * k), text, value_font if k == 0 else small_font, fill=(10, 10, 10))

    prefix = "p8_dimred_density_best_hit_table" if dimred_only else "p8_best_hit_table"
    path = figure_dir / f"{prefix}__{family}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)
    return path


def plot_readable_hit_tables(
    hit_rows: Sequence[dict[str, object]],
    figure_dir: Path,
    datasets: Sequence[str],
    run_tags: Sequence[str],
) -> list[Path]:
    paths: list[Path] = []
    for dimred_only in (False, True):
        rows = best_hit_rows(hit_rows, dimred_only=dimred_only)
        for family in ("efun", "deconv_efun"):
            path = plot_best_hit_table(
                rows,
                figure_dir,
                datasets,
                run_tags,
                family=family,
                dimred_only=dimred_only,
            )
            if path is not None:
                paths.append(path)
    return paths


def slug(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", text.strip())
    text = text.strip("_")
    return text or "value"


def feature_family(feature_name: str) -> str:
    name = feature_name.lower()
    if "deconv" in name:
        return "deconv_efun"
    if "efun" in name:
        return "efun"
    return "other"


def display_label(text: str, max_len: int = 28) -> str:
    if len(text) <= max_len:
        return text
    return text[: max_len - 1] + "..."


def plot_label(text: str) -> str:
    parts = text.split(" | ")
    if parts:
        parts[0] = density_display_name(parts[0])
    return " | ".join(parts)


def heat_color(value: float, vmin: float, vmax: float) -> tuple[int, int, int]:
    if not math.isfinite(value):
        return (235, 235, 235)
    if vmax <= vmin:
        t = 0.5
    else:
        t = max(0.0, min(1.0, (value - vmin) / (vmax - vmin)))
    # White -> blue -> red, readable for sparse scientific heatmaps.
    if t < 0.5:
        a = t / 0.5
        r = int(245 * (1 - a) + 80 * a)
        g = int(247 * (1 - a) + 150 * a)
        b = int(250 * (1 - a) + 210 * a)
    else:
        a = (t - 0.5) / 0.5
        r = int(80 * (1 - a) + 190 * a)
        g = int(150 * (1 - a) + 55 * a)
        b = int(210 * (1 - a) + 45 * a)
    return (r, g, b)


def load_font(size: int, bold: bool = False):
    if ImageFont is None:
        return None
    candidates = [
        r"C:\Windows\Fonts\arialbd.ttf" if bold else r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\calibrib.ttf" if bold else r"C:\Windows\Fonts\calibri.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except Exception:
            pass
    return ImageFont.load_default()


def draw_text(draw, xy, text: str, font, fill=(20, 20, 20), anchor=None):
    if anchor is None:
        draw.text(xy, text, font=font, fill=fill)
    else:
        draw.text(xy, text, font=font, fill=fill, anchor=anchor)


def plot_heatmap(
    matrix: dict[tuple[str, str], float],
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    title: str,
    path: Path,
    *,
    value_format: str = ".2f",
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if Image is None:
        write_svg_heatmap(matrix, row_labels, col_labels, title, path.with_suffix(".svg"), value_format=value_format)
        return

    cell_w = 150
    cell_h = 42
    max_row_len = max([len(r) for r in row_labels] + [1])
    left = max(210, min(520, 28 + max_row_len * 8))
    top = 110
    right = 40
    bottom = 55
    max_col_len = max([len(c) for c in col_labels] + [1])
    if max_col_len > 18:
        top = 160
    width = left + right + cell_w * max(1, len(col_labels))
    height = top + bottom + cell_h * max(1, len(row_labels))
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(22, True)
    label_font = load_font(14, False)
    small_font = load_font(12, False)
    cell_font = load_font(13, True)

    values = [v for v in matrix.values() if math.isfinite(v)]
    vmin = min(values) if values else 0.0
    vmax = max(values) if values else 1.0

    draw_text(draw, (20, 20), title, title_font, fill=(15, 25, 35))
    draw_text(draw, (20, 52), f"color scale: {vmin:.3g} to {vmax:.3g}", small_font, fill=(80, 80, 80))

    for j, col in enumerate(col_labels):
        x = left + j * cell_w + cell_w / 2
        draw_text(draw, (x, top - 10), display_label(col, 34), label_font, fill=(20, 20, 20), anchor="ms")

    for i, row in enumerate(row_labels):
        y = top + i * cell_h + cell_h / 2
        draw_text(draw, (left - 10, y), row, label_font, fill=(20, 20, 20), anchor="rm")
        for j, col in enumerate(col_labels):
            x0 = left + j * cell_w
            y0 = top + i * cell_h
            value = matrix.get((row, col), math.nan)
            draw.rectangle([x0, y0, x0 + cell_w, y0 + cell_h], fill=heat_color(value, vmin, vmax), outline=(210, 210, 210))
            text = "" if not math.isfinite(value) else format(value, value_format)
            if text:
                draw_text(draw, (x0 + cell_w / 2, y0 + cell_h / 2), text, cell_font, fill=(10, 10, 10), anchor="mm")

    image.save(path)


def write_svg_heatmap(
    matrix: dict[tuple[str, str], float],
    row_labels: Sequence[str],
    col_labels: Sequence[str],
    title: str,
    path: Path,
    *,
    value_format: str = ".2f",
) -> None:
    cell_w = 150
    cell_h = 42
    max_row_len = max([len(r) for r in row_labels] + [1])
    left = max(210, min(520, 28 + max_row_len * 8))
    top = 125
    width = left + 40 + cell_w * max(1, len(col_labels))
    height = top + 55 + cell_h * max(1, len(row_labels))
    values = [v for v in matrix.values() if math.isfinite(v)]
    vmin = min(values) if values else 0.0
    vmax = max(values) if values else 1.0

    def color(value: float) -> str:
        r, g, b = heat_color(value, vmin, vmax)
        return f"rgb({r},{g},{b})"

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="20" y="35" font-family="Arial" font-size="22" font-weight="700">{escape_xml(title)}</text>',
        f'<text x="20" y="62" font-family="Arial" font-size="12" fill="#666">color scale: {vmin:.3g} to {vmax:.3g}</text>',
    ]
    for j, col in enumerate(col_labels):
        x = left + j * cell_w + cell_w / 2
        lines.append(f'<text x="{x}" y="{top - 12}" font-family="Arial" font-size="13" text-anchor="middle">{escape_xml(display_label(col, 34))}</text>')
    for i, row in enumerate(row_labels):
        y = top + i * cell_h + cell_h / 2
        lines.append(f'<text x="{left - 10}" y="{y + 4}" font-family="Arial" font-size="13" text-anchor="end">{escape_xml(row)}</text>')
        for j, col in enumerate(col_labels):
            x0 = left + j * cell_w
            y0 = top + i * cell_h
            value = matrix.get((row, col), math.nan)
            lines.append(f'<rect x="{x0}" y="{y0}" width="{cell_w}" height="{cell_h}" fill="{color(value)}" stroke="#d0d0d0"/>')
            if math.isfinite(value):
                lines.append(f'<text x="{x0 + cell_w/2}" y="{y0 + cell_h/2 + 5}" font-family="Arial" font-size="13" font-weight="700" text-anchor="middle">{format(value, value_format)}</text>')
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")


def escape_xml(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
        .replace("'", "&apos;")
    )


def make_heatmap_matrix(rows: Sequence[ScoreRow], row_field: str, col_field: str) -> tuple[dict[tuple[str, str], float], list[str], list[str]]:
    row_labels: list[str] = []
    col_labels: list[str] = []
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in rows:
        row_label = getattr(row, row_field)
        col_label = plot_label(getattr(row, col_field))
        if row_label not in row_labels:
            row_labels.append(row_label)
        if col_label not in col_labels:
            col_labels.append(col_label)
        values[(row_label, col_label)].append(row.mean_peak_abs_corr)
    matrix = {key: mean(finite(vals)) for key, vals in values.items() if finite(vals)}
    return matrix, row_labels, col_labels


def dimred_method_k_summary(score_rows: Sequence[ScoreRow]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str, str, int], list[ScoreRow]] = defaultdict(list)
    for row in score_rows:
        if row.source_level != "per_density_feature":
            continue
        meta = density_metadata(row.density_name)
        if meta["density_source_kind"] != "dimred":
            continue
        family = feature_family(row.bold_feature)
        family_values = ("all", family) if family in {"efun", "deconv_efun"} else ("all",)
        for family_value in family_values:
            grouped[
                (
                    row.run_tag,
                    row.observable,
                    str(meta["density_condition"]),
                    str(meta["density_method"]),
                    family_value,
                    int(meta["density_k"]),
                )
            ].append(row)

    out: list[dict[str, object]] = []
    for (run_tag, observable, condition, method, family, k), rows in sorted(grouped.items()):
        values = [r.mean_peak_abs_corr for r in rows if math.isfinite(r.mean_peak_abs_corr)]
        datasets = sorted({r.dataset for r in rows})
        if not values:
            continue
        out.append(
            {
                "run_tag": run_tag,
                "observable": observable,
                "density_condition": condition,
                "density_method": method,
                "density_k": k,
                "density_method_k": f"{method}_k{k:02d}",
                "feature_family": family,
                "n_rows": len(rows),
                "n_datasets": len(datasets),
                "datasets": ";".join(datasets),
                "mean_peak_abs_corr": f"{mean(values):.10g}",
                "std_peak_abs_corr": f"{pstdev(values):.10g}" if len(values) > 1 else "0",
                "min_peak_abs_corr": f"{min(values):.10g}",
                "max_peak_abs_corr": f"{max(values):.10g}",
                "range_peak_abs_corr": f"{(max(values) - min(values)):.10g}",
            }
        )
    return out


def plot_dimred_method_k_overviews(
    score_rows: Sequence[ScoreRow],
    figure_dir: Path,
    run_tags: Sequence[str],
) -> list[Path]:
    paths: list[Path] = []
    method_k_cols = [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]
    row_labels = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]

    for condition in ("abs", "csplit"):
        for family in ("all", "efun", "deconv_efun"):
            grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
            for row in score_rows:
                if row.source_level != "per_density_feature":
                    continue
                meta = density_metadata(row.density_name)
                if meta["density_source_kind"] != "dimred":
                    continue
                if meta["density_condition"] != condition:
                    continue
                row_family = feature_family(row.bold_feature)
                if family != "all" and row_family != family:
                    continue
                method = str(meta["density_method"])
                k = int(meta["density_k"])
                if method not in METHOD_ORDER or k not in COMPONENT_COUNTS:
                    continue
                observable = RUN_TAG_LABELS.get(row.run_tag, row.run_tag)
                grouped[(observable, f"{method}_k{k:02d}")].append(row.mean_peak_abs_corr)

            matrix = {key: mean(finite(vals)) for key, vals in grouped.items() if finite(vals)}
            if not matrix:
                continue
            path = figure_dir / f"p8_dimred_method_k_score_overview__{condition}__{family}.png"
            title = f"P8 dimred {condition} method x k mean top5 score | {family}"
            plot_heatmap(matrix, row_labels, method_k_cols, title, path, value_format=".3f")
            paths.append(path)
    return paths


TOP_N = 5
P8_FAMILIES = ("all", "efun", "deconv_efun")
P8_PREFERENCE_SCOPES = ("all", "dimred_only", "abs_only", "csplit_only")
LAG_ZERO_THRESHOLD_SEC = 1.0


def density_source_order() -> list[str]:
    ordered = ["blp_evt", "raw_abs_q070"]
    for method in METHOD_ORDER:
        for k in COMPONENT_COUNTS:
            ordered.append(f"dim_abs_{method}{k}_q070")
    ordered.append("raw_csplit_q070")
    for method in METHOD_ORDER:
        for k in COMPONENT_COUNTS:
            ordered.append(f"dim_csplit_{method}{k}_q070")
    return ordered


def density_source_labels() -> list[str]:
    return [density_display_name(name) for name in density_source_order()]


def lag_sign(value: float, threshold_sec: float = LAG_ZERO_THRESHOLD_SEC) -> str:
    if not math.isfinite(value):
        return "nan"
    if value > threshold_sec:
        return "positive"
    if value < -threshold_sec:
        return "negative"
    return "near_zero"


def category_entropy(counts: Counter[str]) -> float:
    total = sum(counts.values())
    if total <= 0:
        return math.nan
    out = 0.0
    for count in counts.values():
        if count <= 0:
            continue
        p = count / total
        out -= p * math.log2(p)
    return out


def most_common_value(values: Sequence[str]) -> tuple[str, int, float]:
    cleaned = [v for v in values if v != ""]
    if not cleaned:
        return "", 0, math.nan
    counts = Counter(cleaned)
    value, count = counts.most_common(1)[0]
    return value, count, count / len(cleaned)


def family_score_rows_from_hits(hit_rows: Sequence[dict[str, object]], top_n: int = TOP_N) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in hit_rows:
        if row.get("source_level") != "per_density_feature":
            continue
        row_family = str(row.get("feature_family", ""))
        if row_family not in {"efun", "deconv_efun"}:
            continue
        for family in ("all", row_family):
            key = (
                str(row.get("dataset", "")),
                str(row.get("run_tag", "")),
                family,
                str(row.get("density_name", "")),
            )
            grouped[key].append(dict(row))

    out: list[dict[str, object]] = []
    for (dataset, run_tag, family, density), rows in sorted(grouped.items()):
        valid = [r for r in rows if math.isfinite(hit_value(r))]
        valid.sort(key=hit_value, reverse=True)
        selected = valid[:top_n]
        if not selected:
            continue
        abs_vals = finite(hit_value(r) for r in selected)
        corr_vals = finite(hit_value(r, "peak_corr") for r in selected)
        lag_vals = finite(hit_value(r, "peak_lag_sec") for r in selected)
        features = [str(r.get("bold_feature", "")) for r in selected]
        density_indices = [str(r.get("density_index", "")) for r in selected if str(r.get("density_index", ""))]
        bold_indices = [str(r.get("bold_mode_index", "")) for r in selected if str(r.get("bold_mode_index", ""))]
        majority_feature, majority_feature_n, majority_feature_fraction = most_common_value(features)
        meta = density_metadata(density)
        out.append(
            {
                "dataset": dataset,
                "run_tag": run_tag,
                "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                "feature_family": family,
                "density_name": density,
                "density_display": density_display_name(density),
                "density_group": density_group_name(density),
                **density_metadata_fields(density),
                "n_top_rows": len(selected),
                "top_n": top_n,
                "mean_peak_abs_corr": mean(abs_vals) if abs_vals else math.nan,
                "std_peak_abs_corr": pstdev(abs_vals) if len(abs_vals) > 1 else 0.0,
                "min_peak_abs_corr": min(abs_vals) if abs_vals else math.nan,
                "max_peak_abs_corr": max(abs_vals) if abs_vals else math.nan,
                "range_peak_abs_corr": (max(abs_vals) - min(abs_vals)) if abs_vals else math.nan,
                "mean_peak_corr": mean(corr_vals) if corr_vals else math.nan,
                "mean_peak_lag_sec": mean(lag_vals) if lag_vals else math.nan,
                "median_peak_lag_sec": median(lag_vals) if lag_vals else math.nan,
                "lag_sign": lag_sign(mean(lag_vals) if lag_vals else math.nan),
                "majority_bold_feature": majority_feature,
                "majority_bold_feature_n": majority_feature_n,
                "majority_bold_feature_fraction": majority_feature_fraction,
                "selected_bold_features": ";".join(features),
                "selected_density_indices": ";".join(density_indices),
                "selected_bold_mode_indices": ";".join(bold_indices),
                "source_kind_order": {"event": 0, "raw": 1, "dimred": 2}.get(str(meta["density_source_kind"]), 9),
            }
        )
    return out


def family_strength_summary(family_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in family_rows:
        grouped[(str(row["run_tag"]), str(row["feature_family"]), str(row["density_name"]))].append(row)

    out: list[dict[str, object]] = []
    for (run_tag, family, density), rows in sorted(grouped.items()):
        vals = finite(as_float(str(r.get("mean_peak_abs_corr", ""))) for r in rows)
        lags = finite(as_float(str(r.get("mean_peak_lag_sec", ""))) for r in rows)
        datasets = sorted({str(r["dataset"]) for r in rows})
        if not vals:
            continue
        signs = [lag_sign(as_float(str(r.get("mean_peak_lag_sec", "")))) for r in rows]
        sign_counts = Counter(s for s in signs if s != "nan")
        majority_sign, majority_sign_n, majority_sign_fraction = most_common_value(signs)
        positive = sign_counts.get("positive", 0)
        negative = sign_counts.get("negative", 0)
        zero = sign_counts.get("near_zero", 0)
        feature_values = [str(r.get("majority_bold_feature", "")) for r in rows]
        majority_feature, majority_feature_n, majority_feature_fraction = most_common_value(feature_values)
        meta_fields = density_metadata_fields(density)
        out.append(
            {
                "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                "run_tag": run_tag,
                "feature_family": family,
                "density_name": density,
                "density_display": density_display_name(density),
                "majority_density": density_display_name(density),
                "density_group": density_group_name(density),
                **meta_fields,
                "condition": meta_fields["density_condition"],
                "method": meta_fields["density_method"],
                "k": meta_fields["density_k"],
                "method_k": meta_fields["density_method_k"],
                "n_datasets": len(datasets),
                "datasets": ";".join(datasets),
                "mean_score": f"{mean(vals):.10g}",
                "std_score": f"{pstdev(vals):.10g}" if len(vals) > 1 else "0",
                "min_score": f"{min(vals):.10g}",
                "max_score": f"{max(vals):.10g}",
                "range_score": f"{(max(vals) - min(vals)):.10g}",
                "lag_median_sec": f"{median(lags):.10g}" if lags else "",
                "lag_mean_sec": f"{mean(lags):.10g}" if lags else "",
                "lag_min_sec": f"{min(lags):.10g}" if lags else "",
                "lag_max_sec": f"{max(lags):.10g}" if lags else "",
                "lag_range_sec": f"{(max(lags) - min(lags)):.10g}" if lags else "",
                "lag_sign_majority": majority_sign,
                "lag_sign_agreement": f"{majority_sign_fraction:.10g}" if math.isfinite(majority_sign_fraction) else "",
                "lag_n_positive": positive,
                "lag_n_negative": negative,
                "lag_n_near_zero": zero,
                "lag_has_opposite_sign": int(positive > 0 and negative > 0),
                "majority_bold_feature": majority_feature,
                "majority_bold_feature_n": majority_feature_n,
                "majority_bold_feature_fraction": f"{majority_feature_fraction:.10g}" if math.isfinite(majority_feature_fraction) else "",
            }
        )
    return out


def main_ranking_rows(summary_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    rows = [dict(r) for r in summary_rows]

    def sort_key(row: dict[str, object]) -> tuple[object, ...]:
        n_datasets = int(as_float(str(row.get("n_datasets", "0"))) or 0)
        mean_score = as_float(str(row.get("mean_score", "")))
        std_score = as_float(str(row.get("std_score", "")))
        range_score = as_float(str(row.get("range_score", "")))
        lag_agreement = as_float(str(row.get("lag_sign_agreement", "")))
        return (
            -n_datasets,
            -(mean_score if math.isfinite(mean_score) else -1),
            std_score if math.isfinite(std_score) else 999,
            range_score if math.isfinite(range_score) else 999,
            -(lag_agreement if math.isfinite(lag_agreement) else -1),
            str(row.get("observable", "")),
            str(row.get("feature_family", "")),
            str(row.get("density_display", "")),
        )

    return sorted(rows, key=sort_key)


def preference_scope_match(row: dict[str, object], scope: str) -> bool:
    condition = str(row.get("density_condition", ""))
    kind = str(row.get("density_source_kind", ""))
    if scope == "all":
        return True
    if scope == "dimred_only":
        return kind == "dimred"
    if scope == "abs_only":
        return condition == "abs"
    if scope == "csplit_only":
        return condition == "csplit"
    return False


def preference_winner_rows(family_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in family_rows:
        for scope in P8_PREFERENCE_SCOPES:
            if preference_scope_match(row, scope):
                grouped[(str(row["dataset"]), str(row["run_tag"]), str(row["feature_family"]), scope)].append(row)

    out: list[dict[str, object]] = []
    for (dataset, run_tag, family, scope), rows in sorted(grouped.items()):
        valid = [r for r in rows if math.isfinite(as_float(str(r.get("mean_peak_abs_corr", ""))))]
        if not valid:
            continue
        best = max(valid, key=lambda r: as_float(str(r.get("mean_peak_abs_corr", ""))))
        out.append(
            {
                "dataset": dataset,
                "run_tag": run_tag,
                "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                "feature_family": family,
                "preference_scope": scope,
                "winner_density_name": best["density_name"],
                "winner_density_display": best["density_display"],
                "winner_density_source_kind": best["density_source_kind"],
                "winner_density_condition": best["density_condition"],
                "winner_density_method": best["density_method"],
                "winner_density_k": best["density_k"],
                "winner_density_method_k": best["density_method_k"],
                "winner_bold_feature": best["majority_bold_feature"],
                "winner_score": f"{as_float(str(best.get('mean_peak_abs_corr', ''))):.10g}",
                "winner_lag_sec": f"{as_float(str(best.get('mean_peak_lag_sec', ''))):.10g}",
                "winner_lag_sign": best["lag_sign"],
                "winner_selected_bold_features": best["selected_bold_features"],
                "winner_selected_density_indices": best["selected_density_indices"],
                "winner_selected_bold_mode_indices": best["selected_bold_mode_indices"],
            }
        )
    return out


def preference_agreement_rows(winner_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    category_fields = (
        "winner_density_source_kind",
        "winner_density_condition",
        "winner_density_method",
        "winner_density_k",
        "winner_density_method_k",
        "winner_density_display",
        "winner_bold_feature",
        "winner_lag_sign",
    )
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in winner_rows:
        grouped[(str(row["run_tag"]), str(row["feature_family"]), str(row["preference_scope"]))].append(row)

    out: list[dict[str, object]] = []
    for (run_tag, family, scope), rows in sorted(grouped.items()):
        datasets = sorted({str(r["dataset"]) for r in rows})
        for field in category_fields:
            values = [str(r.get(field, "")) if str(r.get(field, "")) else "none" for r in rows]
            counts = Counter(values)
            majority, majority_n = counts.most_common(1)[0]
            total = sum(counts.values())
            out.append(
                {
                    "run_tag": run_tag,
                    "observable": RUN_TAG_LABELS.get(run_tag, run_tag),
                    "feature_family": family,
                    "preference_scope": scope,
                    "category_field": field,
                    "majority_category": majority,
                    "majority_count": majority_n,
                    "majority_fraction": f"{(majority_n / total if total else math.nan):.10g}",
                    "entropy": f"{category_entropy(counts):.10g}",
                    "n_datasets": len(datasets),
                    "datasets": ";".join(datasets),
                    "distribution": ";".join(f"{cat}:{count}" for cat, count in counts.most_common()),
                }
            )
    return out


def plot_strength_heatmaps(family_rows: Sequence[dict[str, object]], figure_dir: Path, datasets: Sequence[str], run_tags: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    out_dir = figure_dir / "strength"
    source_order = density_source_order()
    source_labels = density_source_labels()
    source_label_by_name = dict(zip(source_order, source_labels))
    for run_tag in run_tags:
        obs = RUN_TAG_LABELS.get(run_tag, run_tag)
        for family in P8_FAMILIES:
            rows = [r for r in family_rows if r["run_tag"] == run_tag and r["feature_family"] == family]
            grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
            for row in rows:
                label = source_label_by_name.get(str(row["density_name"]), density_display_name(str(row["density_name"])))
                grouped[(str(row["dataset"]), label)].append(as_float(str(row.get("mean_peak_abs_corr", ""))))
            matrix = {key: mean(finite(vals)) for key, vals in grouped.items() if finite(vals)}
            path = out_dir / f"p8_strength_dataset_by_density_source__{slug(run_tag)}__{family}.png"
            title = f"P8 strength consistency | {obs} | {family} | dataset x 51 density sources | top{TOP_N}"
            plot_heatmap(matrix, list(datasets), source_labels, title, path, value_format=".3f")
            paths.append(path)
    return paths


def dimred_method_k_summary_from_family_rows(family_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    summary = family_strength_summary(family_rows)
    out: list[dict[str, object]] = []
    for row in summary:
        if row.get("density_source_kind") != "dimred":
            continue
        out.append(dict(row))
    return out


def plot_method_k_heatmaps(summary_rows: Sequence[dict[str, object]], figure_dir: Path, run_tags: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    out_dir = figure_dir / "method_k"
    method_k_cols = [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]
    row_labels = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]
    for condition in ("abs", "csplit"):
        for family in ("efun", "deconv_efun"):
            subset = [
                r for r in summary_rows
                if r.get("density_condition") == condition and r.get("feature_family") == family
            ]
            matrix = {
                (str(r["observable"]), str(r["density_method_k"])): as_float(str(r.get("mean_score", "")))
                for r in subset
                if math.isfinite(as_float(str(r.get("mean_score", ""))))
            }
            path = out_dir / f"p8_dimred_method_k_score_overview__{condition}__{family}.png"
            title = f"P8 dimred method x k mean top{TOP_N} score | {condition} | {family}"
            plot_heatmap(matrix, row_labels, method_k_cols, title, path, value_format=".3f")
            paths.append(path)
    return paths


def plot_lag_heatmaps(summary_rows: Sequence[dict[str, object]], figure_dir: Path, run_tags: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    out_dir = figure_dir / "lag"
    method_k_cols = [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in COMPONENT_COUNTS]
    row_labels = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]
    for condition in ("abs", "csplit"):
        for family in ("efun", "deconv_efun"):
            subset = [
                r for r in summary_rows
                if r.get("density_condition") == condition and r.get("feature_family") == family
            ]
            lag_matrix = {
                (str(r["observable"]), str(r["density_method_k"])): as_float(str(r.get("lag_median_sec", "")))
                for r in subset
                if math.isfinite(as_float(str(r.get("lag_median_sec", ""))))
            }
            sign_matrix = {
                (str(r["observable"]), str(r["density_method_k"])): as_float(str(r.get("lag_sign_agreement", "")))
                for r in subset
                if math.isfinite(as_float(str(r.get("lag_sign_agreement", ""))))
            }
            lag_path = out_dir / f"p8_dimred_method_k_median_lag_sec__{condition}__{family}.png"
            sign_path = out_dir / f"p8_dimred_method_k_lag_sign_agreement__{condition}__{family}.png"
            plot_heatmap(lag_matrix, row_labels, method_k_cols, f"P8 median lag sec | {condition} | {family}", lag_path, value_format=".1f")
            plot_heatmap(sign_matrix, row_labels, method_k_cols, f"P8 lag sign agreement | {condition} | {family}", sign_path, value_format=".2f")
            paths.extend([lag_path, sign_path])
    return paths


def plot_preference_agreement_heatmaps(winner_rows: Sequence[dict[str, object]], figure_dir: Path, run_tags: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    out_dir = figure_dir / "preference"
    specs = {
        "winner_density_method": ["event", "raw", "svd", "nmf", "mds", "umap", "none"],
        "winner_density_k": ["none"] + [str(k) for k in COMPONENT_COUNTS],
        "winner_density_condition": ["event", "abs", "csplit", "none"],
        "winner_density_source_kind": ["event", "raw", "dimred", "other", "none"],
    }
    row_order = [
        f"{RUN_TAG_LABELS.get(rt, rt)} | {family} | {scope}"
        for rt in run_tags
        for family in P8_FAMILIES
        for scope in P8_PREFERENCE_SCOPES
    ]
    for field, cols in specs.items():
        grouped: dict[tuple[str, str], int] = defaultdict(int)
        totals: Counter[str] = Counter()
        for row in winner_rows:
            rlabel = f"{row['observable']} | {row['feature_family']} | {row['preference_scope']}"
            category = str(row.get(field, "")) or "none"
            grouped[(rlabel, category)] += 1
            totals[rlabel] += 1
        matrix = {
            key: (count / totals[key[0]] if totals[key[0]] else math.nan)
            for key, count in grouped.items()
        }
        path = out_dir / f"p8_preference_distribution__{field}.png"
        title = f"P8 preference consistency distribution | {field} | fraction of sessions"
        plot_heatmap(matrix, row_order, cols, title, path, value_format=".2f")
        paths.append(path)
    return paths


def plot_preference_winner_table(
    winner_rows: Sequence[dict[str, object]],
    figure_dir: Path,
    datasets: Sequence[str],
    run_tags: Sequence[str],
    *,
    family: str,
    scope: str,
) -> Path | None:
    if Image is None:
        return None
    rows = [
        r for r in winner_rows
        if r.get("feature_family") == family and r.get("preference_scope") == scope
    ]
    if not rows:
        return None
    by_key = {(str(r["dataset"]), str(r["run_tag"])): r for r in rows}
    col_labels = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]
    cell_w = 260
    cell_h = 92
    left = 105
    top = 112
    width = left + 30 + cell_w * len(run_tags)
    height = top + 55 + cell_h * len(datasets)
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = load_font(22, True)
    label_font = load_font(14, False)
    small_font = load_font(12, False)
    value_font = load_font(13, True)
    title = f"P8 preference winner | {family} | {scope} | top{TOP_N}"
    draw_text(draw, (20, 22), title, title_font, fill=(15, 25, 35))
    draw_text(draw, (20, 52), "Each cell: winner density, winner BOLD feature, mean |r|, mean lag.", small_font, fill=(80, 80, 80))
    for j, col in enumerate(col_labels):
        x = left + j * cell_w + cell_w / 2
        draw_text(draw, (x, top - 14), col, label_font, anchor="ms")
    values = [hit_value({"peak_abs_corr": r.get("winner_score")}) for r in rows]
    vmax = max(finite(values)) if finite(values) else 1.0
    for i, dataset in enumerate(datasets):
        y = top + i * cell_h + cell_h / 2
        draw_text(draw, (left - 12, y), dataset, label_font, anchor="rm")
        for j, run_tag in enumerate(run_tags):
            x0 = left + j * cell_w
            y0 = top + i * cell_h
            row = by_key.get((dataset, run_tag))
            value = as_float(str(row.get("winner_score", ""))) if row else math.nan
            draw.rectangle([x0, y0, x0 + cell_w, y0 + cell_h], fill=heat_color(value, 0.0, vmax), outline=(205, 205, 205))
            if not row:
                draw_text(draw, (x0 + cell_w / 2, y0 + cell_h / 2), "missing", value_font, anchor="mm")
                continue
            lines = [
                str(row.get("winner_density_display", "")),
                str(row.get("winner_bold_feature", "")),
                f"|r|={value:.3f} lag={as_float(str(row.get('winner_lag_sec', ''))):.1f}s",
            ]
            for k, text in enumerate(lines):
                draw_text(draw, (x0 + 8, y0 + 18 + 22 * k), display_label(text, 32), value_font if k == 0 else small_font, fill=(10, 10, 10))
    path = figure_dir / "preference" / f"p8_preference_winner_table__{scope}__{family}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)
    return path


def plot_preference_winner_tables(winner_rows: Sequence[dict[str, object]], figure_dir: Path, datasets: Sequence[str], run_tags: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    for scope in P8_PREFERENCE_SCOPES:
        for family in P8_FAMILIES:
            path = plot_preference_winner_table(
                winner_rows,
                figure_dir,
                datasets,
                run_tags,
                family=family,
                scope=scope,
            )
            if path is not None:
                paths.append(path)
    return paths


def plot_main_p8_consistency(
    family_rows: Sequence[dict[str, object]],
    summary_rows: Sequence[dict[str, object]],
    winner_rows: Sequence[dict[str, object]],
    figure_dir: Path,
    datasets: Sequence[str],
    run_tags: Sequence[str],
) -> list[Path]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    paths.extend(plot_strength_heatmaps(family_rows, figure_dir, datasets, run_tags))
    paths.extend(plot_method_k_heatmaps(summary_rows, figure_dir, run_tags))
    paths.extend(plot_preference_agreement_heatmaps(winner_rows, figure_dir, run_tags))
    paths.extend(plot_preference_winner_tables(winner_rows, figure_dir, datasets, run_tags))
    paths.extend(plot_lag_heatmaps(summary_rows, figure_dir, run_tags))
    return paths


def plot_all(score_rows: Sequence[ScoreRow], figure_dir: Path, datasets: Sequence[str], run_tags: Sequence[str]) -> list[Path]:
    paths: list[Path] = []
    figure_dir.mkdir(parents=True, exist_ok=True)

    for run_tag in run_tags:
        obs = RUN_TAG_LABELS.get(run_tag, run_tag)
        density_rows = [
            r for r in score_rows
            if r.run_tag == run_tag and r.source_level == "per_density"
        ]
        if density_rows:
            matrix, _, cols = make_heatmap_matrix(density_rows, "dataset", "density_name")
            path = figure_dir / f"p8_density_score_heatmap__{slug(run_tag)}.png"
            plot_heatmap(matrix, list(datasets), cols, f"P8 cross-session density score | {obs}", path)
            paths.append(path)

        feature_rows = [
            ScoreRow(
                dataset=r.dataset,
                run_tag=r.run_tag,
                observable=r.observable,
                density_name=f"{r.density_name} | {r.bold_feature}",
                bold_feature=r.bold_feature,
                source_level=r.source_level,
                n_top_rows=r.n_top_rows,
                mean_peak_abs_corr=r.mean_peak_abs_corr,
                max_peak_abs_corr=r.max_peak_abs_corr,
                mean_peak_corr=r.mean_peak_corr,
                mean_peak_lag_sec=r.mean_peak_lag_sec,
                source_csv=r.source_csv,
            )
            for r in score_rows
            if r.run_tag == run_tag and r.source_level == "per_density_feature"
        ]
        if feature_rows:
            matrix, _, cols = make_heatmap_matrix(feature_rows, "dataset", "density_name")
            path = figure_dir / f"p8_density_feature_score_heatmap__{slug(run_tag)}.png"
            plot_heatmap(matrix, list(datasets), cols, f"P8 density x feature score | {obs}", path)
            paths.append(path)

        for family in ("efun", "deconv_efun"):
            family_source_rows = [
                r for r in score_rows
                if r.run_tag == run_tag
                and r.source_level == "per_density_feature"
                and feature_family(r.bold_feature) == family
            ]
            if family_source_rows:
                matrix, _, cols = make_heatmap_matrix(family_source_rows, "dataset", "density_name")
                path = figure_dir / f"p8_{family}_density_score_heatmap__{slug(run_tag)}.png"
                plot_heatmap(matrix, list(datasets), cols, f"P8 {family} density score | {obs}", path)
                paths.append(path)

                family_feature_rows = [
                    ScoreRow(
                        dataset=r.dataset,
                        run_tag=r.run_tag,
                        observable=r.observable,
                        density_name=f"{r.density_name} | {r.bold_feature}",
                        bold_feature=r.bold_feature,
                        source_level=r.source_level,
                        n_top_rows=r.n_top_rows,
                        mean_peak_abs_corr=r.mean_peak_abs_corr,
                        max_peak_abs_corr=r.max_peak_abs_corr,
                        mean_peak_corr=r.mean_peak_corr,
                        mean_peak_lag_sec=r.mean_peak_lag_sec,
                        source_csv=r.source_csv,
                    )
                    for r in family_source_rows
                ]
                matrix, _, cols = make_heatmap_matrix(family_feature_rows, "dataset", "density_name")
                path = figure_dir / f"p8_{family}_density_feature_score_heatmap__{slug(run_tag)}.png"
                plot_heatmap(matrix, list(datasets), cols, f"P8 {family} density x feature score | {obs}", path)
                paths.append(path)

    overview_rows = [
        ScoreRow(
            dataset=r.observable,
            run_tag=r.run_tag,
            observable=r.observable,
            density_name=r.density_name,
            bold_feature=r.bold_feature,
            source_level=r.source_level,
            n_top_rows=r.n_top_rows,
            mean_peak_abs_corr=r.mean_peak_abs_corr,
            max_peak_abs_corr=r.max_peak_abs_corr,
            mean_peak_corr=r.mean_peak_corr,
            mean_peak_lag_sec=r.mean_peak_lag_sec,
            source_csv=r.source_csv,
        )
        for r in score_rows
        if r.source_level == "per_density"
    ]
    if overview_rows:
        grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in overview_rows:
            grouped[(row.observable, row.density_name)].append(row.mean_peak_abs_corr)
        matrix = {key: mean(finite(vals)) for key, vals in grouped.items() if finite(vals)}
        cols = sorted({r.density_name for r in overview_rows})
        rows = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]
        path = figure_dir / "p8_cross_session_density_mean_score_overview.png"
        plot_heatmap(matrix, rows, cols, "P8 cross-session mean density score", path)
        paths.append(path)

    feature_overview_rows = [
        (
            r.observable,
            f"{r.density_name} | {r.bold_feature}",
            r.mean_peak_abs_corr,
        )
        for r in score_rows
        if r.source_level == "per_density_feature"
    ]
    if feature_overview_rows:
        grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
        for obs, group, value in feature_overview_rows:
            grouped[(obs, group)].append(value)
        matrix = {key: mean(finite(vals)) for key, vals in grouped.items() if finite(vals)}
        cols = sorted({group for _, group, _ in feature_overview_rows})
        rows = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]
        path = figure_dir / "p8_cross_session_density_feature_mean_score_overview.png"
        plot_heatmap(matrix, rows, cols, "P8 cross-session mean density x feature score", path)
        paths.append(path)

    for family in ("efun", "deconv_efun"):
        family_rows = [
            r for r in score_rows
            if r.source_level == "per_density_feature"
            and feature_family(r.bold_feature) == family
        ]
        if not family_rows:
            continue

        grouped_density: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in family_rows:
            grouped_density[(row.observable, row.density_name)].append(row.mean_peak_abs_corr)
        matrix = {key: mean(finite(vals)) for key, vals in grouped_density.items() if finite(vals)}
        cols = sorted({r.density_name for r in family_rows})
        rows = [RUN_TAG_LABELS.get(rt, rt) for rt in run_tags]
        path = figure_dir / f"p8_{family}_cross_session_density_mean_score_overview.png"
        plot_heatmap(matrix, rows, cols, f"P8 {family} cross-session mean density score", path)
        paths.append(path)

        grouped_feature: dict[tuple[str, str], list[float]] = defaultdict(list)
        for row in family_rows:
            grouped_feature[(row.observable, f"{row.density_name} | {row.bold_feature}")].append(row.mean_peak_abs_corr)
        matrix = {key: mean(finite(vals)) for key, vals in grouped_feature.items() if finite(vals)}
        cols = sorted({f"{r.density_name} | {r.bold_feature}" for r in family_rows})
        path = figure_dir / f"p8_{family}_cross_session_density_feature_mean_score_overview.png"
        plot_heatmap(matrix, rows, cols, f"P8 {family} cross-session mean density x feature score", path)
        paths.append(path)

    paths.extend(plot_dimred_method_k_overviews(score_rows, figure_dir, run_tags))

    return paths


def main() -> None:
    args = parse_args()
    score_rows, raw_top_rows, hit_rows = collect_scores(
        args.processed_root,
        args.datasets,
        args.run_tags,
        args.xcorr_save_tag,
    )
    args.results_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.results_dir / "p8_xcorr_combined_top_rows.csv", raw_top_rows)
    write_csv(args.results_dir / "p8_top_xcorr_hits_readable.csv", hit_rows)

    score_dicts = [row_to_dict(r) for r in score_rows]
    write_csv(args.results_dir / "p8_xcorr_topn_scores_long.csv", score_dicts)
    write_csv(
        args.results_dir / "p8_xcorr_density_scores.csv",
        [row_to_dict(r) for r in score_rows if r.source_level == "per_density"],
    )
    write_csv(
        args.results_dir / "p8_xcorr_density_feature_scores.csv",
        [row_to_dict(r) for r in score_rows if r.source_level == "per_density_feature"],
    )

    density_summary = summarize_cross_session(score_rows, "per_density")
    feature_summary = summarize_cross_session(score_rows, "per_density_feature")
    method_k_summary = dimred_method_k_summary(score_rows)
    family_rows = family_score_rows_from_hits(hit_rows, TOP_N)
    strength_summary = family_strength_summary(family_rows)
    main_ranking = main_ranking_rows(strength_summary)
    preference_winners = preference_winner_rows(family_rows)
    preference_agreement = preference_agreement_rows(preference_winners)
    main_method_k_summary = dimred_method_k_summary_from_family_rows(family_rows)
    write_csv(args.results_dir / "p8_cross_session_density_summary.csv", density_summary)
    write_csv(args.results_dir / "p8_cross_session_density_feature_summary.csv", feature_summary)
    write_csv(args.results_dir / "p8_dimred_method_k_summary.csv", method_k_summary)
    write_csv(args.results_dir / "p8_strength_family_density_scores.csv", family_rows)
    write_csv(args.results_dir / "p8_strength_cross_session_summary.csv", strength_summary)
    write_csv(args.results_dir / "p8_main_consistency_ranking.csv", main_ranking)
    write_csv(args.results_dir / "p8_method_k_consistency_summary.csv", main_method_k_summary)
    write_csv(args.results_dir / "p8_preference_winners_by_dataset.csv", preference_winners)
    write_csv(args.results_dir / "p8_preference_agreement_summary.csv", preference_agreement)

    dominant = dominant_rows(score_rows, "per_density_feature")
    agreement = (
        categorical_agreement(dominant, "dominant_density_name")
        + categorical_agreement(dominant, "dominant_bold_feature")
        + categorical_agreement(dominant, "dominant_density_feature")
    )
    write_csv(args.results_dir / "p8_dominant_density_feature_by_dataset.csv", dominant)
    write_csv(args.results_dir / "p8_categorical_agreement_summary.csv", agreement)
    write_csv(
        args.results_dir / "p8_best_hit_by_dataset_observable_family.csv",
        best_hit_rows(hit_rows, dimred_only=False),
    )
    write_csv(
        args.results_dir / "p8_best_dimred_density_hit_by_dataset_observable_family.csv",
        best_hit_rows(hit_rows, dimred_only=True),
    )

    figure_paths: list[Path] = []
    if not args.skip_figures:
        figure_paths = plot_main_p8_consistency(
            family_rows,
            strength_summary,
            preference_winners,
            args.figure_dir,
            args.datasets,
            args.run_tags,
        )
        if args.results_figure_dir.resolve() != args.figure_dir.resolve():
            figure_paths.extend(
                plot_main_p8_consistency(
                    family_rows,
                    strength_summary,
                    preference_winners,
                    args.results_figure_dir,
                    args.datasets,
                    args.run_tags,
                )
            )

    print(f"Datasets: {', '.join(args.datasets)}")
    print(f"Run tags: {', '.join(args.run_tags)}")
    print(f"Score rows: {len(score_rows)}")
    print(f"Combined top rows: {len(raw_top_rows)}")
    print(f"Readable hit rows: {len(hit_rows)}")
    print("P8 topN: 5")
    print(f"Results dir: {args.results_dir}")
    print(f"Figure dir: {args.figure_dir}")
    print(f"Results figure dir: {args.results_figure_dir}")
    print(f"Figures: {len(figure_paths)}")


if __name__ == "__main__":
    main()
