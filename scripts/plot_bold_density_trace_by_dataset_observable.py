"""Plot BOLD efun/deconv traces against BLP density by dataset and observable.

This is a preview/diagnostic layer built on top of
analyze_bold_efun_top_raw_density_hits.py.  It reuses the trace readers,
normalization, xcorr-sign alignment, mask shading, and fixed zoom-window
helpers, but changes the selection unit to:

    dataset x BOLD observable x BOLD feature family x density class

The goal is to avoid hiding BOLD-observable differences behind one winner per
dataset.
"""

from __future__ import annotations

import argparse
import csv
import math
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import MatplotlibDeprecationWarning

warnings.filterwarnings("ignore", category=MatplotlibDeprecationWarning)

from analyze_bold_efun_top_raw_density_hits import (
    FAMILY_ORDER,
    OBS_ORDER,
    TraceCache,
    annotate_trace_mask_and_borders,
    apply_xcorr_valid_mask,
    as_float,
    as_int,
    color_for_family,
    combined_ylim,
    normalize_trace,
    read_xcorr_bold_post_file,
    select_fixed_summary_windows,
    time_axis,
    win_to_wsl,
    write_csv,
    xcorr_alignment_factor,
)


DEFAULT_HITS = (
    Path("results")
    / "pipeline8_10_strict_band_coupling_current"
    / "p8_p10_strict_band_hits_long.csv"
)
DEFAULT_RESULTS_DIR = Path("results") / "e10gb1_bold_density_trace_by_observable_preview_20260529"
DEFAULT_STABLE_METHOD_K = (
    "mds_k04",
    "mds_k05",
    "mds_k06",
    "mds_k07",
    "mds_k08",
    "nmf_k04",
    "umap_k04",
    "umap_k05",
    "umap_k06",
    "umap_k08",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits-csv", type=Path, default=DEFAULT_HITS)
    parser.add_argument("--dataset", default="e10gb1")
    parser.add_argument("--pipeline", default="P8", choices=["P8", "P10"])
    parser.add_argument(
        "--density-classes",
        nargs="+",
        default=["raw_efun_density", "dimred_efun_density"],
        choices=["raw_efun_density", "dimred_efun_density", "event_density"],
    )
    parser.add_argument("--bold-feature-filter", choices=["real", "abs", "all"], default="real")
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--window-sec", type=float, default=120.0)
    parser.add_argument(
        "--top-hits-per-group",
        type=int,
        default=1,
        help=(
            "Number of top xcorr pairs to keep for each "
            "BOLD observable x efun/deconv family x density class."
        ),
    )
    parser.add_argument(
        "--selection-mode",
        choices=["pair", "unique_density", "unique_bold_mode"],
        default="pair",
        help=(
            "pair keeps top xcorr pairs; unique_density keeps each BLP "
            "density/component once with its best BOLD mode; unique_bold_mode "
            "keeps each BOLD mode once with its best BLP density."
        ),
    )
    parser.add_argument(
        "--method-k-scope",
        choices=["all_method_k", "p5_stable_method_k"],
        default="all_method_k",
        help="For dimred density traces, optionally restrict to the current P5-stable method-k set.",
    )
    parser.add_argument("--stable-method-k", nargs="+", default=list(DEFAULT_STABLE_METHOD_K))
    parser.add_argument("--whole-trace-time-unit", choices=["sec", "min", "hour"], default="min")
    parser.add_argument("--align-xcorr-sign", action="store_true", default=True)
    parser.add_argument("--fig-width", type=float, default=32.0)
    parser.add_argument("--full-row-height", type=float, default=2.8)
    parser.add_argument("--zoom-row-height", type=float, default=2.4)
    return parser.parse_args()


def read_csv_rows(path: Path | str) -> list[dict[str, str]]:
    with win_to_wsl(path).open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def feature_allowed(row: dict[str, str], feature_filter: str) -> bool:
    if feature_filter == "all":
        return True
    return str(row.get("bold_feature", "")).lower().endswith(f"_{feature_filter}")


def row_sort_key(row: dict[str, str]) -> tuple[int, str]:
    try:
        obs_idx = OBS_ORDER.index(str(row.get("bold_observable")))
    except ValueError:
        obs_idx = len(OBS_ORDER)
    try:
        fam_idx = FAMILY_ORDER.index(str(row.get("bold_feature_family")))
    except ValueError:
        fam_idx = len(FAMILY_ORDER)
    return (obs_idx * 10 + fam_idx, str(row.get("bold_observable")))


@lru_cache(maxsize=512)
def source_rows(source_csv: str) -> list[dict[str, str]]:
    return read_csv_rows(source_csv)


def resolve_density_file(row: dict[str, object]) -> tuple[str, str]:
    density_name = str(row.get("density_name", ""))
    density_index = as_int(row.get("density_index"))
    bold_feature = str(row.get("bold_feature", ""))
    bold_mode_index = as_int(row.get("bold_mode_index"))
    source_csv = str(row.get("source_csv", ""))
    for src in source_rows(source_csv):
        if src.get("density_name") != density_name:
            continue
        if as_int(src.get("density_index")) != density_index:
            continue
        if bold_feature and src.get("bold_feature") != bold_feature:
            continue
        if bold_mode_index > 0 and as_int(src.get("bold_mode_index")) != bold_mode_index:
            continue
        return src.get("density_file", ""), src.get("density_field", "")
    for src in source_rows(source_csv):
        if src.get("density_name") == density_name and as_int(src.get("density_index")) == density_index:
            return src.get("density_file", ""), src.get("density_field", "")
    return "", ""


def density_label_for_class(density_class: str) -> str:
    if density_class == "raw_efun_density":
        return "raw BLP efun density"
    if density_class == "dimred_efun_density":
        return "dimred BLP efun density"
    return "event density"


def add_density_trace_legend(fig: plt.Figure, density_class: str) -> None:
    handles: list[object] = [
        Line2D([0], [0], color="#7c3aed", lw=1.4, label=density_label_for_class(density_class)),
        Line2D([0], [0], color=color_for_family("efun"), lw=1.4, label="BOLD efun_real"),
        Line2D([0], [0], color=color_for_family("deconv_efun"), lw=1.4, label="BOLD deconv_efun_real"),
        Patch(facecolor="#9ca3af", alpha=0.18, label="masked samples"),
        Line2D([0], [0], color="#6b7280", lw=0.8, ls="--", label="session border"),
    ]
    fig.legend(handles=handles, loc="upper right", fontsize=8)


def select_top_rows(
    hits_csv: Path,
    dataset: str,
    pipeline: str,
    density_class: str,
    bold_feature_filter: str,
    top_hits_per_group: int,
    selection_mode: str,
    method_k_scope: str = "all_method_k",
    stable_method_k: Sequence[str] = DEFAULT_STABLE_METHOD_K,
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = {}
    dataset_l = dataset.lower()
    stable_method_k_set = set(stable_method_k)
    with hits_csv.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            if row.get("pipeline") != pipeline:
                continue
            if row.get("dataset", "").lower() != dataset_l:
                continue
            if row.get("density_class") != density_class:
                continue
            if (
                density_class == "dimred_efun_density"
                and method_k_scope == "p5_stable_method_k"
                and row.get("density_method_k") not in stable_method_k_set
            ):
                continue
            if row.get("bold_feature_family") not in FAMILY_ORDER:
                continue
            if not feature_allowed(row, bold_feature_filter):
                continue
            peak = as_float(row.get("peak_abs_corr"))
            if not math.isfinite(peak):
                continue
            key = (str(row.get("bold_observable")), str(row.get("bold_feature_family")))
            grouped.setdefault(key, []).append(dict(row))

    selected: list[dict[str, object]] = []
    for key, rows in grouped.items():
        rows.sort(key=lambda r: as_float(r.get("peak_abs_corr")), reverse=True)
        unique_rows: list[dict[str, object]] = []
        seen_pairs: set[tuple[object, ...]] = set()
        for row in rows:
            if selection_mode == "unique_density":
                pair_key = (
                    str(row.get("density_name", "")),
                    as_int(row.get("density_index")),
                    str(row.get("density_method_k", "")),
                    str(row.get("strict_label", "")),
                )
            elif selection_mode == "unique_bold_mode":
                pair_key = (
                    str(row.get("bold_feature", "")),
                    as_int(row.get("bold_mode_index")),
                )
            else:
                pair_key = (
                    str(row.get("density_name", "")),
                    as_int(row.get("density_index")),
                    str(row.get("density_method_k", "")),
                    as_int(row.get("bold_mode_index")),
                    str(row.get("bold_feature", "")),
                    str(row.get("strict_label", "")),
                )
            if pair_key in seen_pairs:
                continue
            seen_pairs.add(pair_key)
            unique_rows.append(row)
            if len(unique_rows) >= max(1, top_hits_per_group):
                break
        for rank, row in enumerate(unique_rows, start=1):
            row = dict(row)
            row["selection_group"] = "|".join(key)
            row["selection_rank"] = rank
            density_file, density_field = resolve_density_file(row)
            row["resolved_density_file"] = density_file
            row["resolved_density_field"] = density_field
            selected.append(row)
    selected.sort(key=row_sort_key)
    return selected


def trace_data_for_row(
    row: dict[str, object],
    cache: TraceCache,
    align_xcorr_sign: bool,
) -> dict[str, object]:
    xcorr_mat, bold_post = read_xcorr_bold_post_file(str(row.get("source_csv", "")))
    if bold_post is None:
        raise FileNotFoundError("BOLD_POST path not resolved")
    density_file = str(row.get("resolved_density_file", ""))
    if not density_file:
        raise FileNotFoundError("density file not resolved from source CSV")
    density_index = as_int(row.get("density_index"))
    bold_mode_index = as_int(row.get("bold_mode_index"))
    t_den, den = cache.density_trace(density_file, density_index)
    t_bold, bold = cache.bold_trace(bold_post, str(row.get("bold_feature", "")), bold_mode_index)
    n = min(den.size, bold.size, t_den.size, t_bold.size)
    t = t_den[:n]
    den_z = normalize_trace(den[:n])
    bold_z = normalize_trace(bold[:n])
    align_factor = xcorr_alignment_factor(row, align_xcorr_sign)
    bold_z = align_factor * bold_z
    mask = cache.border_mask(xcorr_mat)
    borders = cache.session_border_times(bold_post)
    valid = np.isfinite(den_z) & np.isfinite(bold_z)
    valid = apply_xcorr_valid_mask(valid, mask, n)
    if int(np.sum(valid)) < 8:
        valid = np.isfinite(den_z) & np.isfinite(bold_z)
    return {
        "xcorr_mat": xcorr_mat,
        "bold_post": bold_post,
        "t": t,
        "den_z": den_z,
        "bold_z": bold_z,
        "mask": mask,
        "borders": borders,
        "valid": valid,
        "align_factor": align_factor,
    }


def compact_density_title(row: dict[str, object]) -> str:
    cls = str(row.get("density_class", ""))
    if cls == "dimred_efun_density":
        label = str(row.get("strict_label") or "unlabeled")
        return (
            f"{row.get('density_method_k')} C{as_int(row.get('density_index'))} "
            f"{label}"
        )
    if cls == "raw_efun_density":
        return f"raw {as_int(row.get('density_index'))}"
    return f"{row.get('density_name')}:{as_int(row.get('density_index'))}"


def panel_title(prefix: str, row: dict[str, object], align_factor: float) -> str:
    peak_abs = as_float(row.get("peak_abs_corr"))
    peak = as_float(row.get("peak_corr"))
    lag = as_float(row.get("peak_lag_sec"))
    return (
        f"{prefix} | {row.get('bold_observable')} | {row.get('bold_feature_family')} | "
        f"BOLD m{as_int(row.get('bold_mode_index'))} | {compact_density_title(row)} | "
        f"|r|={peak_abs:.2f} r={peak:.2f} "
        f"align={align_factor:+.0f} lag={lag:.1f}s"
    )


def plot_full_trace(
    rows: Sequence[dict[str, object]],
    out_file: Path,
    dataset: str,
    pipeline: str,
    density_class: str,
    time_unit: str,
    align_xcorr_sign: bool,
    fig_width: float,
    row_height: float,
) -> list[dict[str, object]]:
    max_rank = max((as_int(row.get("selection_rank")) for row in rows), default=1)
    if max_rank > 1:
        return plot_full_trace_top_hits(
            rows,
            out_file,
            dataset,
            pipeline,
            density_class,
            time_unit,
            align_xcorr_sign,
            fig_width,
            row_height,
            max_rank,
        )

    observables = [obs for obs in OBS_ORDER if any(row.get("bold_observable") == obs for row in rows)]
    n_rows = max(1, len(observables))
    n_cols = len(FAMILY_ORDER)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(fig_width, max(row_height * n_rows, 4.0)),
        sharex=False,
        sharey=False,
    )
    axes = np.asarray(axes)
    if axes.ndim == 1:
        axes = axes.reshape(n_rows, n_cols)
    for ax in axes.ravel():
        ax.axis("off")

    row_lookup = {
        (str(row.get("bold_observable")), str(row.get("bold_feature_family"))): row
        for row in rows
    }
    cache = TraceCache()
    records: list[dict[str, object]] = []
    for i_obs, obs in enumerate(observables):
        for j_fam, family in enumerate(FAMILY_ORDER):
            ax = axes[i_obs, j_fam]
            ax.axis("on")
            row = row_lookup.get((obs, family))
            if row is None:
                ax.text(0.5, 0.5, "missing", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{obs} | {family}", fontsize=8)
                continue
            try:
                data = trace_data_for_row(row, cache, align_xcorr_sign)
                t = data["t"]
                den_z = data["den_z"]
                bold_z = data["bold_z"]
                x, unit_label = time_axis(t, time_unit)
                annotate_trace_mask_and_borders(
                    ax,
                    t,
                    data["mask"],
                    time_unit,
                    origin=0.0,
                    session_borders=data["borders"],
                )
                ax.plot(x, den_z, color="#7c3aed", lw=0.7, alpha=0.75)
                ax.plot(x, bold_z, color=color_for_family(family), lw=1.1, alpha=0.92)
                ax.grid(True, color="#e5e7eb", lw=0.5)
                if i_obs == n_rows - 1:
                    ax.set_xlabel(f"time in concatenated session ({unit_label})")
                if j_fam == 0:
                    ax.set_ylabel(f"{obs}\nrobust z", fontsize=8)
                ax.set_title(panel_title("full", row, float(data["align_factor"])), fontsize=8)
                records.append(
                    {
                        **row,
                        "trace_kind": "full_trace",
                        "bold_trace_alignment_factor": data["align_factor"],
                        "xcorr_mat": str(data["xcorr_mat"]),
                        "bold_post_file": str(data["bold_post"]),
                    }
                )
            except Exception as exc:
                ax.text(0.5, 0.5, f"trace error:\n{exc}", ha="center", va="center", transform=ax.transAxes, fontsize=8)
                records.append({**row, "trace_kind": "full_trace", "trace_error": str(exc)})

    add_density_trace_legend(fig, density_class)
    fig.suptitle(
        f"{dataset} {pipeline} {density_class} vs BOLD observable full traces | signed aligned",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_file, dpi=220)
    plt.close(fig)
    return records


def plot_full_trace_top_hits(
    rows: Sequence[dict[str, object]],
    out_file: Path,
    dataset: str,
    pipeline: str,
    density_class: str,
    time_unit: str,
    align_xcorr_sign: bool,
    fig_width: float,
    row_height: float,
    max_rank: int,
) -> list[dict[str, object]]:
    groups: list[tuple[str, str]] = []
    for obs in OBS_ORDER:
        for family in FAMILY_ORDER:
            if any(row.get("bold_observable") == obs and row.get("bold_feature_family") == family for row in rows):
                groups.append((obs, family))
    for row in rows:
        key = (str(row.get("bold_observable")), str(row.get("bold_feature_family")))
        if key not in groups:
            groups.append(key)

    n_rows = max(1, len(groups))
    n_cols = max(1, max_rank)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(max(fig_width, 6.0 * n_cols), max(row_height * n_rows, 4.0)),
        sharex=False,
        sharey=False,
    )
    axes = np.asarray(axes)
    if axes.ndim == 1:
        axes = axes.reshape(n_rows, n_cols)
    for ax in axes.ravel():
        ax.axis("off")

    row_lookup = {
        (
            str(row.get("bold_observable")),
            str(row.get("bold_feature_family")),
            as_int(row.get("selection_rank")),
        ): row
        for row in rows
    }
    cache = TraceCache()
    records: list[dict[str, object]] = []
    for i_group, (obs, family) in enumerate(groups):
        for j_rank in range(n_cols):
            rank = j_rank + 1
            ax = axes[i_group, j_rank]
            ax.axis("on")
            row = row_lookup.get((obs, family, rank))
            if row is None:
                ax.text(0.5, 0.5, "missing", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"top{rank}", fontsize=8)
                continue
            try:
                data = trace_data_for_row(row, cache, align_xcorr_sign)
                t = data["t"]
                den_z = data["den_z"]
                bold_z = data["bold_z"]
                x, unit_label = time_axis(t, time_unit)
                annotate_trace_mask_and_borders(
                    ax,
                    t,
                    data["mask"],
                    time_unit,
                    origin=0.0,
                    session_borders=data["borders"],
                )
                ax.plot(x, den_z, color="#7c3aed", lw=0.65, alpha=0.75)
                ax.plot(x, bold_z, color=color_for_family(family), lw=1.0, alpha=0.92)
                ax.grid(True, color="#e5e7eb", lw=0.5)
                if i_group == n_rows - 1:
                    ax.set_xlabel(f"time ({unit_label})")
                if j_rank == 0:
                    ax.set_ylabel(f"{obs}\n{family}\nrobust z", fontsize=8)
                ax.set_title(panel_title(f"top{rank}", row, float(data["align_factor"])), fontsize=7)
                records.append(
                    {
                        **row,
                        "trace_kind": "full_trace",
                        "bold_trace_alignment_factor": data["align_factor"],
                        "xcorr_mat": str(data["xcorr_mat"]),
                        "bold_post_file": str(data["bold_post"]),
                    }
                )
            except Exception as exc:
                ax.text(0.5, 0.5, f"trace error:\n{exc}", ha="center", va="center", transform=ax.transAxes, fontsize=8)
                records.append({**row, "trace_kind": "full_trace", "trace_error": str(exc)})

    add_density_trace_legend(fig, density_class)
    fig.suptitle(
        f"{dataset} {pipeline} {density_class} vs BOLD observable full traces | top{max_rank} per group | signed aligned",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_file, dpi=220)
    plt.close(fig)
    return records


def plot_fixed_zoom(
    rows: Sequence[dict[str, object]],
    out_file: Path,
    dataset: str,
    pipeline: str,
    density_class: str,
    window_sec: float,
    align_xcorr_sign: bool,
    fig_width: float,
    row_height: float,
) -> list[dict[str, object]]:
    cache = TraceCache()
    n_rows = max(1, len(rows))
    n_cols = 3
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(fig_width, max(row_height * n_rows, 4.0)),
        sharex=True,
        sharey=False,
    )
    axes = np.asarray(axes)
    if axes.ndim == 1:
        axes = axes.reshape(n_rows, n_cols)
    records: list[dict[str, object]] = []

    for i_row, row in enumerate(rows):
        try:
            data = trace_data_for_row(row, cache, align_xcorr_sign)
            t = data["t"]
            den_z = data["den_z"]
            bold_z = data["bold_z"]
            row_ylim = combined_ylim([den_z, bold_z])
            windows = select_fixed_summary_windows(t, den_z, bold_z, data["valid"], window_sec)
            for j_col, win in enumerate(windows[:n_cols]):
                ax = axes[i_row, j_col]
                lo = int(win["start"])
                hi = int(win["stop"])
                tx, unit_label = time_axis(t[lo:hi], "min", origin=float(t[lo]))
                segment_mask = data["mask"][lo:hi] if data["mask"].size >= hi else np.asarray([], dtype=bool)
                annotate_trace_mask_and_borders(
                    ax,
                    t[lo:hi],
                    segment_mask,
                    "min",
                    origin=float(t[lo]),
                    session_borders=data["borders"],
                )
                ax.plot(tx, den_z[lo:hi], color="#7c3aed", lw=0.9, alpha=0.85)
                ax.plot(tx, bold_z[lo:hi], color=color_for_family(str(row.get("bold_feature_family"))), lw=1.0, alpha=0.9)
                if row_ylim is not None:
                    ax.set_ylim(*row_ylim)
                ax.set_xlim(0, window_sec / 60.0)
                ax.grid(True, color="#e5e7eb", lw=0.5)
                if j_col == 0:
                    ax.set_ylabel(f"{row.get('bold_observable')}\n{row.get('bold_feature_family')}\nrobust z", fontsize=8)
                if i_row == n_rows - 1:
                    ax.set_xlabel(f"time within selected window ({unit_label})")
                label = str(win.get("label", "window"))
                start_min = float(win["start_sec"]) / 60.0
                ax.set_title(panel_title(label, row, float(data["align_factor"])) + f" | start={start_min:.1f}m", fontsize=8)
                records.append(
                    {
                        **row,
                        "trace_kind": "fixed_zoom",
                        "zoom_window_label": label,
                        "zoom_start_sec": win["start_sec"],
                        "zoom_stop_sec": win["stop_sec"],
                        "bold_trace_alignment_factor": data["align_factor"],
                    }
                )
        except Exception as exc:
            for j_col in range(n_cols):
                ax = axes[i_row, j_col]
                ax.text(0.5, 0.5, f"trace error:\n{exc}", ha="center", va="center", transform=ax.transAxes, fontsize=8)
            records.append({**row, "trace_kind": "fixed_zoom", "trace_error": str(exc)})

    add_density_trace_legend(fig, density_class)
    fig.suptitle(
        f"{dataset} {pipeline} {density_class} vs BOLD observable fixed zoom | signed aligned | y=full_trace",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_file, dpi=220)
    plt.close(fig)
    return records


def safe_slug(text: str) -> str:
    return "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in text)


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    all_records: list[dict[str, object]] = []
    for density_class in args.density_classes:
        rows = select_top_rows(
            args.hits_csv,
            args.dataset,
            args.pipeline,
            density_class,
            args.bold_feature_filter,
            args.top_hits_per_group,
            args.selection_mode,
            args.method_k_scope,
            args.stable_method_k,
        )
        top_slug = f"top{max(1, args.top_hits_per_group)}"
        select_slug = safe_slug(args.selection_mode)
        scope_slug = safe_slug(args.method_k_scope)
        write_csv(args.results_dir / f"{args.pipeline.lower()}_{args.dataset}_{density_class}_{scope_slug}_{select_slug}_{top_slug}_selected_rows.csv", rows)
        if not rows:
            continue
        slug = safe_slug(density_class)
        full_records = plot_full_trace(
            rows,
            args.results_dir / f"{args.pipeline.lower()}_{args.dataset}_{slug}_{scope_slug}_{select_slug}_{top_slug}_full_trace_by_observable.png",
            args.dataset,
            args.pipeline,
            density_class,
            args.whole_trace_time_unit,
            args.align_xcorr_sign,
            args.fig_width,
            args.full_row_height,
        )
        zoom_records = plot_fixed_zoom(
            rows,
            args.results_dir / f"{args.pipeline.lower()}_{args.dataset}_{slug}_{scope_slug}_{select_slug}_{top_slug}_fixed_zoom_by_observable.png",
            args.dataset,
            args.pipeline,
            density_class,
            args.window_sec,
            args.align_xcorr_sign,
            args.fig_width,
            args.zoom_row_height,
        )
        all_records.extend(full_records)
        all_records.extend(zoom_records)
    write_csv(
        args.results_dir
        / f"{args.pipeline.lower()}_{args.dataset}_{safe_slug(args.method_k_scope)}_{safe_slug(args.selection_mode)}_top{max(1, args.top_hits_per_group)}_trace_records.csv",
        all_records,
    )
    print(f"Wrote {args.results_dir}")
    print(f"Trace records: {len(all_records)}")


if __name__ == "__main__":
    main()
