"""Probe BOLD efun top raw-density hits.

This diagnostic asks whether the strongest raw BLP density couplings to BOLD
efun/deconv-efun features are concentrated in slow or fast raw Koopman modes.

Inputs are existing numeric outputs:
- P8/P10 strict-band coupling hit table
- P5 raw-density mode metadata
- P5 raw-density MAT files, for trace examples
- P8 xcorr MAT + P7 BOLD_POST MAT, for BOLD efun/deconv traces

It does not recompute xcorr.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


DEFAULT_HITS = Path("results") / "pipeline8_10_strict_band_coupling_current" / "p8_p10_strict_band_hits_long.csv"
DEFAULT_RAW_META = (
    Path("results")
    / "p8_p10_parameter_selection_current_rmsenv_adaptive_v2"
    / "p5_raw_density_mode_metadata.csv"
)
DEFAULT_RESULTS_DIR = Path("results") / "bold_efun_raw_density_slow_state_probe_current"
DEFAULT_TOP_N = (1, 3, 5, 10, 20, 50)
DEFAULT_DENSITY = "raw_csplit_q070_standardize_rmsenv_adaptive"
DATASET_ORDER = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01", "k13m17", "k13m23")
FAMILY_ORDER = ("efun", "deconv_efun")
PIPELINE_ORDER = ("P8", "P10")
OBS_ORDER = ("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits-csv", type=Path, default=DEFAULT_HITS)
    parser.add_argument("--raw-metadata-csv", type=Path, default=DEFAULT_RAW_META)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--figure-dir", type=Path, default=None)
    parser.add_argument("--density-name", default=DEFAULT_DENSITY)
    parser.add_argument(
        "--bold-feature-filter",
        choices=["all", "real", "abs"],
        default="all",
        help="Restrict BOLD-side features before ranking. Use real for efun_real/deconv_real vs density.",
    )
    parser.add_argument("--top-n-values", nargs="+", type=int, default=list(DEFAULT_TOP_N))
    parser.add_argument("--trace-top-n", type=int, default=10)
    parser.add_argument("--trace-pipeline", default="P8", choices=list(PIPELINE_ORDER))
    parser.add_argument("--trace-max-datasets", type=int, default=7)
    parser.add_argument(
        "--exclude-datasets",
        nargs="*",
        default=[],
        help="Datasets to exclude from raw-timescale summaries and trace examples, e.g. k13m17.",
    )
    parser.add_argument("--trace-half-window-bins", type=int, default=90)
    parser.add_argument("--whole-trace-fig-width", type=float, default=30.0)
    parser.add_argument("--whole-trace-row-height", type=float, default=3.0)
    parser.add_argument("--whole-trace-time-unit", choices=["sec", "min", "hour"], default="min")
    parser.add_argument("--zoom-window-sec", type=float, default=300.0)
    parser.add_argument("--zoom-top-k", type=int, default=3)
    parser.add_argument("--zoom-family", choices=["efun", "deconv_efun", "both"], default="deconv_efun")
    parser.add_argument(
        "--zoom-window-mode",
        choices=["fixed_summary", "local_corr"],
        default="fixed_summary",
        help=(
            "fixed_summary shows raw-density-high, BOLD-feature-high, and representative middle windows. "
            "local_corr is the older exploratory mode that selects windows by local correlation."
        ),
    )
    parser.add_argument(
        "--zoom-y-scale",
        choices=["full_trace", "window_auto"],
        default="full_trace",
        help="Use full_trace to keep each zoom row on the same y-scale as the whole selected trace.",
    )
    parser.add_argument(
        "--align-xcorr-sign",
        action="store_true",
        help="Multiply BOLD traces by sign(peak_corr), keeping raw density fixed, so negative xcorr hits are polarity-aligned.",
    )
    parser.add_argument("--skip-traces", action="store_true")
    parser.add_argument(
        "--legacy-local-peak-traces",
        action="store_true",
        help="Also write the older local-peak-aligned trace sheet. Not recommended for whole-trace coupling interpretation.",
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
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
        out = float(str(value))
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def as_int(value: object, default: int = -1) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def finite(values: Iterable[float]) -> list[float]:
    return [v for v in values if math.isfinite(v)]


def win_to_wsl(path: str | Path) -> Path:
    text = str(path)
    if os.name != "nt" and len(text) > 2 and text[1] == ":":
        drive = text[0].lower()
        rest = text[2:].lstrip("\\/").replace("\\", "/")
        return Path(f"/mnt/{drive}/{rest}")
    return Path(text)


def wsl_to_win(path: str | Path) -> str:
    text = str(path)
    m = re.match(r"^/mnt/([a-zA-Z])/(.*)$", text)
    if m:
        rest = m.group(2).replace("/", "\\")
        return f"{m.group(1).upper()}:\\{rest}"
    return text


def matlab_char(ds: h5py.Dataset) -> str:
    arr = np.asarray(ds)
    if arr.dtype == np.uint16:
        return "".join(chr(int(x)) for x in arr.ravel(order="F") if int(x) != 0)
    if arr.dtype.kind in {"S", "U"}:
        return "".join(str(x) for x in arr.ravel(order="F"))
    return str(arr)


def read_h5_numeric(ds: h5py.Dataset) -> np.ndarray:
    arr = np.asarray(ds)
    if arr.dtype.fields and {"real", "imag"}.issubset(arr.dtype.fields):
        return arr["real"] + 1j * arr["imag"]
    return arr


def normalize_trace(x: np.ndarray) -> np.ndarray:
    y = np.asarray(x, dtype=float).copy()
    ok = np.isfinite(y)
    if not np.any(ok):
        return y
    med = float(np.nanmedian(y[ok]))
    p05, p95 = np.nanpercentile(y[ok], [5, 95])
    scale = float(p95 - p05)
    if not math.isfinite(scale) or scale <= 0:
        scale = float(np.nanstd(y[ok]))
    if not math.isfinite(scale) or scale <= 0:
        scale = 1.0
    return (y - med) / scale


def time_axis(t: np.ndarray, unit: str, origin: float = 0.0) -> tuple[np.ndarray, str]:
    scale = {"sec": 1.0, "min": 60.0, "hour": 3600.0}[unit]
    label = {"sec": "sec", "min": "min", "hour": "hours"}[unit]
    return (np.asarray(t, dtype=float) - float(origin)) / scale, label


def xcorr_alignment_factor(row: dict[str, object], align_xcorr_sign: bool) -> float:
    if not align_xcorr_sign:
        return 1.0
    signed_peak = as_float(row.get("peak_corr"))
    if math.isfinite(signed_peak) and signed_peak < 0:
        return -1.0
    return 1.0


def combined_ylim(arrays: Sequence[np.ndarray], margin_fraction: float = 0.04) -> tuple[float, float] | None:
    vals: list[np.ndarray] = []
    for arr in arrays:
        x = np.asarray(arr, dtype=float).ravel()
        x = x[np.isfinite(x)]
        if x.size:
            vals.append(x)
    if not vals:
        return None
    all_vals = np.concatenate(vals)
    lo = float(np.nanmin(all_vals))
    hi = float(np.nanmax(all_vals))
    if not math.isfinite(lo) or not math.isfinite(hi):
        return None
    if lo == hi:
        pad = max(1.0, abs(lo) * 0.1)
    else:
        pad = max(1e-6, (hi - lo) * margin_fraction)
    return lo - pad, hi + pad


def order_index(value: str, order: Sequence[str]) -> int:
    try:
        return order.index(value)
    except ValueError:
        return len(order)


def context_key(row: dict[str, object]) -> tuple[object, ...]:
    return (
        row.get("pipeline", ""),
        row.get("dataset", ""),
        row.get("run_tag", ""),
        row.get("bold_feature_family", ""),
        row.get("p9_feature", ""),
        row.get("p9_method_k", ""),
        row.get("source_level", ""),
    )


def feature_allowed(row: dict[str, object], feature_filter: str) -> bool:
    if feature_filter == "all":
        return True
    feature = str(row.get("bold_feature", "")).lower()
    return feature.endswith(f"_{feature_filter}")


def load_metadata(path: Path, density_name: str) -> dict[tuple[str, str, int], dict[str, object]]:
    out: dict[tuple[str, str, int], dict[str, object]] = {}
    for row in read_csv(path):
        if row.get("density_name_for_p8_p10") != density_name:
            continue
        dataset = str(row.get("dataset", "")).lower()
        density_index = as_int(row.get("density_index"))
        if not dataset or density_index < 1:
            continue
        item: dict[str, object] = dict(row)
        for key in (
            "density_index",
            "source_mode_index",
            "raw_efun_index",
            "timescale_sec_preferred",
            "frequency_hz_preferred",
            "envelope_window_sec",
            "timescale_sec_discrete_log",
            "frequency_hz_discrete_angle",
            "timescale_sec",
            "frequency_hz",
        ):
            item[key] = as_float(row.get(key))
        out[(dataset, density_name, density_index)] = item
    return out


def enrich_hit(row: dict[str, str], metadata: dict[tuple[str, str, int], dict[str, object]]) -> dict[str, object]:
    out: dict[str, object] = dict(row)
    dataset = str(row.get("dataset", "")).lower()
    density_name = str(row.get("density_name", ""))
    density_index = as_int(row.get("density_index"))
    meta = metadata.get((dataset, density_name, density_index), {})
    for key in (
        "density_index",
        "rank_in_source_csv",
        "bold_mode_index",
        "peak_abs_corr",
        "peak_corr",
        "peak_lag_sec",
        "p9_k",
    ):
        if key in out:
            out[key] = as_float(out[key])
    out["raw_efun_index"] = meta.get("raw_efun_index", math.nan)
    out["source_mode_index"] = meta.get("source_mode_index", math.nan)
    out["raw_timescale_sec"] = meta.get("timescale_sec_preferred", math.nan)
    out["raw_frequency_hz"] = meta.get("frequency_hz_preferred", math.nan)
    out["raw_envelope_window_sec"] = meta.get("envelope_window_sec", math.nan)
    out["raw_density_file"] = meta.get("density_file", "")
    out["timescale_bin"] = timescale_bin(as_float(out.get("raw_timescale_sec")))
    return out


def timescale_bin(tau: float) -> str:
    if not math.isfinite(tau) or tau <= 0:
        return "NA"
    if tau < 0.1:
        return "<0.1s"
    if tau < 1:
        return "0.1-1s"
    if tau < 10:
        return "1-10s"
    if tau < 100:
        return "10-100s"
    return ">=100s"


def build_ranked_selections(
    rows: Sequence[dict[str, object]],
    top_n_values: Sequence[int],
    density_name: str,
    bold_feature_filter: str,
) -> dict[str, dict[int, list[dict[str, object]]]]:
    grouped_all: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    grouped_raw: dict[tuple[object, ...], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        if str(row.get("bold_feature_family")) not in FAMILY_ORDER:
            continue
        if not feature_allowed(row, bold_feature_filter):
            continue
        if str(row.get("density_condition")) not in {"event", "csplit"}:
            continue
        if not math.isfinite(as_float(row.get("peak_abs_corr"))):
            continue
        key = context_key(row)
        grouped_all[key].append(row)
        if str(row.get("density_name")) == density_name:
            grouped_raw[key].append(row)

    selections: dict[str, dict[int, list[dict[str, object]]]] = {
        "global_topN_raw_subset": {int(n): [] for n in top_n_values},
        "raw_topN": {int(n): [] for n in top_n_values},
    }
    max_n = max(int(n) for n in top_n_values)
    for key, group_rows in grouped_all.items():
        ranked = sorted(group_rows, key=lambda r: as_float(r.get("peak_abs_corr")), reverse=True)[:max_n]
        for rank, row in enumerate(ranked, start=1):
            row["rank_in_context"] = rank
        for n in top_n_values:
            top = ranked[: int(n)]
            selections["global_topN_raw_subset"][int(n)].extend(
                dict(row, selection_scope="global_topN_raw_subset", top_n=int(n))
                for row in top
                if str(row.get("density_name")) == density_name
            )
    for key, group_rows in grouped_raw.items():
        ranked = sorted(group_rows, key=lambda r: as_float(r.get("peak_abs_corr")), reverse=True)[:max_n]
        for rank, row in enumerate(ranked, start=1):
            row["rank_in_raw_context"] = rank
        for n in top_n_values:
            selections["raw_topN"][int(n)].extend(
                dict(row, selection_scope="raw_topN", top_n=int(n))
                for row in ranked[: int(n)]
            )
    return selections


def summarize_timescales(selections: dict[str, dict[int, list[dict[str, object]]]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for scope, by_n in selections.items():
        for top_n, rows in by_n.items():
            grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
            for row in rows:
                grouped[(str(row.get("pipeline")), str(row.get("bold_feature_family")), str(row.get("bold_observable")))].append(row)
            for (pipeline, family, observable), group_rows in sorted(grouped.items()):
                taus = finite(as_float(row.get("raw_timescale_sec")) for row in group_rows)
                peaks = finite(as_float(row.get("peak_abs_corr")) for row in group_rows)
                out.append(
                    {
                        "selection_scope": scope,
                        "top_n": top_n,
                        "pipeline": pipeline,
                        "bold_feature_family": family,
                        "bold_observable": observable,
                        "n_hits": len(group_rows),
                        "n_timescale": len(taus),
                        "median_tau_sec": float(np.median(taus)) if taus else math.nan,
                        "mean_tau_sec": float(np.mean(taus)) if taus else math.nan,
                        "p25_tau_sec": float(np.percentile(taus, 25)) if taus else math.nan,
                        "p75_tau_sec": float(np.percentile(taus, 75)) if taus else math.nan,
                        "slow_ge_1s_fraction": float(np.mean(np.asarray(taus) >= 1.0)) if taus else math.nan,
                        "slow_ge_10s_fraction": float(np.mean(np.asarray(taus) >= 10.0)) if taus else math.nan,
                        "median_peak_abs_corr": float(np.median(peaks)) if peaks else math.nan,
                    }
                )
    return out


def write_selection_csv(path: Path, selections: dict[str, dict[int, list[dict[str, object]]]]) -> None:
    rows: list[dict[str, object]] = []
    for scope, by_n in selections.items():
        for _top_n, items in by_n.items():
            rows.extend(items)
    write_csv(path, rows)


def color_for_family(family: str) -> str:
    return {"efun": "#2b6cb0", "deconv_efun": "#c2410c"}.get(family, "#555555")


def add_trace_legend(fig: plt.Figure, families: Sequence[str]) -> None:
    family_set = set(families)
    handles: list[object] = [
        Line2D([0], [0], color="#7c3aed", lw=1.4, label="raw BLP efun density"),
    ]
    if "efun" in family_set:
        handles.append(Line2D([0], [0], color=color_for_family("efun"), lw=1.4, label="BOLD efun_real"))
    if "deconv_efun" in family_set:
        handles.append(Line2D([0], [0], color=color_for_family("deconv_efun"), lw=1.4, label="BOLD deconv_efun_real"))
    handles.extend(
        [
            Patch(facecolor="#9ca3af", alpha=0.18, label="masked samples"),
            Line2D([0], [0], color="#6b7280", lw=0.8, ls="--", label="session border"),
        ]
    )
    fig.legend(handles=handles, loc="upper right", fontsize=8)


def plot_ecdf(selections: dict[str, dict[int, list[dict[str, object]]]], figure_dir: Path) -> None:
    scope = "raw_topN"
    top_values = [1, 3, 5, 10, 20, 50]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True, sharey=True)
    axes = axes.ravel()
    panels = [("P8", "efun"), ("P8", "deconv_efun"), ("P10", "efun"), ("P10", "deconv_efun")]
    for ax, (pipeline, family) in zip(axes, panels):
        for top_n in top_values:
            rows = [
                row
                for row in selections.get(scope, {}).get(top_n, [])
                if row.get("pipeline") == pipeline and row.get("bold_feature_family") == family
            ]
            taus = np.asarray(finite(as_float(row.get("raw_timescale_sec")) for row in rows), dtype=float)
            if taus.size == 0:
                continue
            x = np.sort(np.log10(np.clip(taus, 1e-3, None)))
            y = np.arange(1, x.size + 1) / x.size
            ax.step(x, y, where="post", label=f"top{top_n}", lw=1.8)
        ax.set_title(f"{pipeline} | {family}")
        ax.grid(True, color="#e5e7eb")
        ax.axvline(0, color="#999999", lw=1, ls="--")
        ax.axvline(1, color="#999999", lw=1, ls=":")
    for ax in axes[2:]:
        ax.set_xlabel("log10(raw efun timescale sec)")
    for ax in axes[::2]:
        ax.set_ylabel("ECDF")
    axes[0].legend(loc="lower right", fontsize=8)
    fig.suptitle("Raw BLP efun density top-hit timescale ECDF | raw-topN within each context")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(figure_dir / "01_raw_topN_timescale_ecdf_by_pipeline_family.png", dpi=180)
    plt.close(fig)


def plot_median_by_topn(summary_rows: Sequence[dict[str, object]], figure_dir: Path) -> None:
    rows = [row for row in summary_rows if row.get("selection_scope") == "raw_topN"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, pipeline in zip(axes, PIPELINE_ORDER):
        for family in FAMILY_ORDER:
            xs: list[int] = []
            ys: list[float] = []
            for top_n in DEFAULT_TOP_N:
                vals = [
                    as_float(row.get("median_tau_sec"))
                    for row in rows
                    if row.get("pipeline") == pipeline
                    and row.get("bold_feature_family") == family
                    and int(row.get("top_n", -1)) == int(top_n)
                ]
                vals = finite(vals)
                if not vals:
                    continue
                xs.append(int(top_n))
                ys.append(float(np.median(vals)))
            if xs:
                ax.plot(xs, ys, marker="o", lw=2, color=color_for_family(family), label=family)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(list(DEFAULT_TOP_N))
        ax.set_xticklabels([str(x) for x in DEFAULT_TOP_N])
        ax.grid(True, which="both", color="#e5e7eb")
        ax.axhline(1, color="#999999", lw=1, ls="--")
        ax.axhline(10, color="#999999", lw=1, ls=":")
        ax.set_title(pipeline)
        ax.set_xlabel("raw topN")
    axes[0].set_ylabel("median raw efun timescale (sec)")
    axes[0].legend(loc="best")
    fig.suptitle("Median raw-mode timescale by topN, pooled over available contexts")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(figure_dir / "02_raw_topN_median_timescale_by_topN.png", dpi=180)
    plt.close(fig)


def plot_slow_fraction(summary_rows: Sequence[dict[str, object]], figure_dir: Path) -> None:
    rows = [row for row in summary_rows if row.get("selection_scope") == "raw_topN"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for ax, pipeline in zip(axes, PIPELINE_ORDER):
        for family in FAMILY_ORDER:
            xs: list[int] = []
            ys: list[float] = []
            for top_n in DEFAULT_TOP_N:
                vals = [
                    as_float(row.get("slow_ge_10s_fraction"))
                    for row in rows
                    if row.get("pipeline") == pipeline
                    and row.get("bold_feature_family") == family
                    and int(row.get("top_n", -1)) == int(top_n)
                ]
                vals = finite(vals)
                if not vals:
                    continue
                xs.append(int(top_n))
                ys.append(float(np.mean(vals)))
            if xs:
                ax.plot(xs, ys, marker="o", lw=2, color=color_for_family(family), label=family)
        ax.set_xscale("log")
        ax.set_xticks(list(DEFAULT_TOP_N))
        ax.set_xticklabels([str(x) for x in DEFAULT_TOP_N])
        ax.set_ylim(0, 1)
        ax.grid(True, color="#e5e7eb")
        ax.set_title(pipeline)
        ax.set_xlabel("raw topN")
    axes[0].set_ylabel("fraction with tau >= 10 sec")
    axes[0].legend(loc="best")
    fig.suptitle("Slow raw-mode fraction among top raw-density hits")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(figure_dir / "03_raw_topN_slow_ge10s_fraction_by_topN.png", dpi=180)
    plt.close(fig)


def plot_index_hist(selections: dict[str, dict[int, list[dict[str, object]]]], figure_dir: Path, top_n: int = 10) -> None:
    rows = selections["raw_topN"].get(top_n, [])
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=False, sharey=False)
    axes = axes.ravel()
    panels = [("P8", "efun"), ("P8", "deconv_efun"), ("P10", "efun"), ("P10", "deconv_efun")]
    for ax, (pipeline, family) in zip(axes, panels):
        vals = [
            as_int(row.get("raw_efun_index"))
            for row in rows
            if row.get("pipeline") == pipeline and row.get("bold_feature_family") == family
        ]
        vals = [v for v in vals if v > 0]
        if vals:
            counts = Counter(vals)
            top = counts.most_common(30)
            ax.bar([str(k) for k, _ in top], [v for _, v in top], color=color_for_family(family), alpha=0.85)
            ax.tick_params(axis="x", labelrotation=90)
        ax.set_title(f"{pipeline} | {family} | top raw indices")
        ax.set_ylabel("hit count")
        ax.grid(True, axis="y", color="#e5e7eb")
    fig.suptitle(f"Raw efun index distribution among raw-top{top_n} hits")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(figure_dir / f"04_raw_efun_index_hist_raw_top{top_n}.png", dpi=180)
    plt.close(fig)


def read_xcorr_bold_post_file(source_csv: str) -> tuple[Path, Path | None]:
    csv_path = win_to_wsl(source_csv)
    parts = list(csv_path.parts)
    # .../<run>/feature/<family>/<tag>_peaks__*.csv
    if "feature" in parts:
        root = Path(*parts[: parts.index("feature")])
    elif "density" in parts:
        root = Path(*parts[: parts.index("density")])
    else:
        root = csv_path.parent
    m = re.match(r"(.+?)_peaks__", csv_path.name)
    if not m:
        return root, None
    xcorr_mat = root / f"{m.group(1)}.mat"
    if not xcorr_mat.is_file():
        return xcorr_mat, None
    try:
        with h5py.File(xcorr_mat, "r") as handle:
            bold_post_file = matlab_char(handle["out/bold_post_file"])
    except Exception:
        return xcorr_mat, None
    return xcorr_mat, win_to_wsl(bold_post_file)


class TraceCache:
    def __init__(self) -> None:
        self._density: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        self._bold: dict[tuple[str, str, str], tuple[np.ndarray, np.ndarray]] = {}
        self._border: dict[str, np.ndarray] = {}
        self._session_borders: dict[str, np.ndarray] = {}

    def density_trace(self, density_file: str, density_index: int) -> tuple[np.ndarray, np.ndarray]:
        path = str(win_to_wsl(density_file))
        if path not in self._density:
            with h5py.File(path, "r") as handle:
                if "D/density_time_by_mode" in handle:
                    density = np.asarray(handle["D/density_time_by_mode"], dtype=float)
                else:
                    density = np.asarray(handle["density"], dtype=float)
                if "D/t_centers" in handle:
                    t = np.asarray(handle["D/t_centers"], dtype=float).ravel()
                else:
                    t = np.asarray(handle["t_centers"], dtype=float).ravel()
            self._density[path] = (density, t)
        density, t = self._density[path]
        idx = int(density_index) - 1
        if idx < 0 or idx >= density.shape[0]:
            raise IndexError(f"density_index {density_index} out of bounds for {path}")
        return t, density[idx, :].astype(float)

    def border_mask(self, xcorr_mat: Path) -> np.ndarray:
        key = str(xcorr_mat)
        if key not in self._border:
            try:
                with h5py.File(xcorr_mat, "r") as handle:
                    mask = np.asarray(handle["out/border_mask"], dtype=bool).ravel()
            except Exception:
                mask = np.zeros(0, dtype=bool)
            self._border[key] = mask
        return self._border[key]

    def bold_trace(self, bold_post_file: Path, feature: str, mode_index: int) -> tuple[np.ndarray, np.ndarray]:
        feature_key = "deconv_abs" if feature == "deconv_abs" else "deconv_real" if feature == "deconv_real" else "efun_abs" if feature == "efun_abs" else "efun_real"
        key = (str(bold_post_file), feature_key, str(mode_index))
        if key in self._bold:
            return self._bold[key]
        with h5py.File(bold_post_file, "r") as handle:
            t = np.asarray(handle["BOLD_POST/time_vec"], dtype=float).ravel()
            if feature_key == "efun_abs":
                mat = read_h5_numeric(handle["BOLD_POST/EDMD_outputs/norm_efuns/abs"])
            elif feature_key == "efun_real":
                mat = read_h5_numeric(handle["BOLD_POST/EDMD_outputs/norm_efuns/real"])
            elif feature_key == "deconv_abs":
                mat = read_h5_numeric(handle["BOLD_POST/deconv/norm_efuns/abs_all"])
            else:
                mat = read_h5_numeric(handle["BOLD_POST/deconv/norm_efuns/real_all"])
        if np.iscomplexobj(mat):
            mat = np.abs(mat) if feature_key.endswith("abs") else np.real(mat)
        mat = np.asarray(mat, dtype=float)
        idx = int(mode_index) - 1
        if idx < 0 or idx >= mat.shape[0]:
            raise IndexError(f"bold mode {mode_index} out of bounds for {bold_post_file}")
        trace = mat[idx, :].astype(float)
        self._bold[key] = (t, trace)
        return t, trace

    def session_border_times(self, bold_post_file: Path) -> np.ndarray:
        key = str(bold_post_file)
        if key in self._session_borders:
            return self._session_borders[key]
        borders = np.asarray([], dtype=float)
        try:
            with h5py.File(bold_post_file, "r") as handle:
                if "BOLD_POST/session_border_t" in handle:
                    borders = np.asarray(handle["BOLD_POST/session_border_t"], dtype=float).ravel()
                elif "BOLD_POST/session/start_idx" in handle and "BOLD_POST/time_vec" in handle:
                    starts = np.asarray(handle["BOLD_POST/session/start_idx"], dtype=float).ravel().astype(int)
                    t = np.asarray(handle["BOLD_POST/time_vec"], dtype=float).ravel()
                    valid = starts[(starts > 1) & (starts <= t.size)]
                    borders = t[valid - 1]
        except Exception:
            borders = np.asarray([], dtype=float)
        self._session_borders[key] = borders[np.isfinite(borders)]
        return self._session_borders[key]


def choose_trace_rows(
    selections: dict[str, dict[int, list[dict[str, object]]]],
    top_n: int,
    pipeline: str,
    dataset_order: Sequence[str] = DATASET_ORDER,
) -> list[dict[str, object]]:
    rows = [
        row
        for row in selections["raw_topN"].get(top_n, [])
        if row.get("pipeline") == pipeline
        and row.get("source_level") == "per_feature"
        and row.get("bold_feature_family") in FAMILY_ORDER
    ]
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row.get("dataset")), str(row.get("bold_feature_family")))].append(row)
    selected: list[dict[str, object]] = []
    for dataset in dataset_order:
        for family in FAMILY_ORDER:
            group = grouped.get((dataset, family), [])
            if not group:
                continue
            group = sorted(group, key=lambda r: as_float(r.get("peak_abs_corr")), reverse=True)
            selected.append(dict(group[0], trace_selection="top_peak"))
    return selected


def plot_trace_examples(
    selections: dict[str, dict[int, list[dict[str, object]]]],
    figure_dir: Path,
    top_n: int,
    pipeline: str,
    max_datasets: int,
    half_window_bins: int,
    align_xcorr_sign: bool,
    dataset_order: Sequence[str] = DATASET_ORDER,
) -> list[dict[str, object]]:
    rows = choose_trace_rows(selections, top_n, pipeline, dataset_order)
    rows = [
        row
        for row in rows
        if order_index(str(row.get("dataset")), dataset_order) < max_datasets
    ]
    cache = TraceCache()
    trace_records: list[dict[str, object]] = []
    if not rows:
        return trace_records

    visible_datasets = list(dataset_order[:max_datasets])
    n_rows = len(visible_datasets)
    fig, axes = plt.subplots(n_rows, 2, figsize=(15, max(3.0 * n_rows, 6)), sharex=False, sharey=False)
    if n_rows == 1:
        axes = np.asarray([axes])
    for ax in axes.ravel():
        ax.axis("off")

    row_lookup = {(str(row.get("dataset")), str(row.get("bold_feature_family"))): row for row in rows}
    for i_ds, dataset in enumerate(visible_datasets):
        for j_family, family in enumerate(FAMILY_ORDER):
            ax = axes[i_ds, j_family]
            ax.axis("on")
            row = row_lookup.get((dataset, family))
            if row is None:
                ax.text(0.5, 0.5, "missing", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{dataset} | {family}")
                continue
            try:
                xcorr_mat, bold_post = read_xcorr_bold_post_file(str(row.get("source_csv", "")))
                if bold_post is None:
                    raise FileNotFoundError("BOLD_POST path not resolved")
                density_index = as_int(row.get("density_index"))
                bold_mode_index = as_int(row.get("bold_mode_index"))
                t_den, den = cache.density_trace(str(row.get("raw_density_file", "")), density_index)
                t_bold, bold = cache.bold_trace(bold_post, str(row.get("bold_feature", "")), bold_mode_index)
                n = min(den.size, bold.size, t_den.size, t_bold.size)
                den = den[:n]
                bold = bold[:n]
                t = t_den[:n]
                mask = cache.border_mask(xcorr_mat)
                valid = np.isfinite(den)
                valid = apply_xcorr_valid_mask(valid, mask, n)
                if not np.any(valid):
                    valid = np.isfinite(den)
                center = int(np.nanargmax(np.where(valid, den, np.nan)))
                lo = max(0, center - half_window_bins)
                hi = min(n, center + half_window_bins + 1)
                tx = t[lo:hi] - t[center]
                den_z = normalize_trace(den[lo:hi])
                bold_z = normalize_trace(bold[lo:hi])
                align_factor = xcorr_alignment_factor(row, align_xcorr_sign)
                bold_z = align_factor * bold_z
                segment_mask = mask[lo:hi] if mask.size >= hi else np.asarray([], dtype=bool)
                annotate_trace_mask_and_borders(
                    ax,
                    t[lo:hi],
                    segment_mask,
                    "sec",
                    origin=float(t[center]),
                    session_borders=cache.session_border_times(bold_post),
                )
                ax.plot(tx, den_z, color="#6b46c1", lw=1.5, label="raw density")
                ax.plot(tx, bold_z, color=color_for_family(family), lw=1.5, label=str(row.get("bold_feature")))
                lag = as_float(row.get("peak_lag_sec"))
                if math.isfinite(lag):
                    ax.axvline(lag, color="#555555", lw=1, ls="--", alpha=0.8)
                ax.axvline(0, color="#aaaaaa", lw=0.8)
                tau = as_float(row.get("raw_timescale_sec"))
                peak = as_float(row.get("peak_abs_corr"))
                signed_peak = as_float(row.get("peak_corr"))
                raw_idx = as_int(row.get("raw_efun_index"))
                mode_idx = as_int(row.get("bold_mode_index"))
                ax.set_title(
                    f"{dataset} | {family} | {row.get('bold_observable')} | raw {raw_idx} tau={tau:.2g}s | BOLD mode {mode_idx} | |r|={peak:.3f} r={signed_peak:.3f} align={align_factor:+.0f}",
                    fontsize=8,
                )
                ax.grid(True, color="#e5e7eb")
                if i_ds == n_rows - 1:
                    ax.set_xlabel("seconds from raw-density local max")
                if j_family == 0:
                    ax.set_ylabel("robust z")
                trace_records.append(
                    dict(
                        row,
                        xcorr_mat=str(xcorr_mat),
                        bold_post_file=str(bold_post),
                        trace_center_time_sec=float(t[center]),
                        trace_window_start_sec=float(t[lo]),
                        trace_window_end_sec=float(t[hi - 1]),
                        bold_trace_alignment_factor=float(align_factor),
                    )
                )
            except Exception as exc:
                ax.text(0.5, 0.5, f"trace error:\n{exc}", ha="center", va="center", transform=ax.transAxes, fontsize=8)
                ax.set_title(f"{dataset} | {family}")
                trace_records.append(dict(row, trace_error=str(exc)))
    add_trace_legend(fig, FAMILY_ORDER)
    fig.suptitle(f"{pipeline} BOLD feature vs top raw BLP density trace examples | raw-top{top_n}")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(figure_dir / f"05_{pipeline.lower()}_bold_vs_raw_density_trace_examples_raw_top{top_n}.png", dpi=180)
    plt.close(fig)
    return trace_records


def smooth_trace(x: np.ndarray, window: int = 31) -> np.ndarray:
    if window <= 1 or x.size < 3:
        return x
    window = int(max(3, window))
    if window % 2 == 0:
        window += 1
    if x.size < window:
        window = x.size if x.size % 2 == 1 else x.size - 1
    if window < 3:
        return x
    kernel = np.ones(window, dtype=float) / window
    valid = np.isfinite(x).astype(float)
    y = np.where(np.isfinite(x), x, 0.0)
    num = np.convolve(y, kernel, mode="same")
    den = np.convolve(valid, kernel, mode="same")
    out = num / np.maximum(den, 1e-12)
    out[den <= 0] = np.nan
    return out


def apply_xcorr_valid_mask(valid: np.ndarray, mask: np.ndarray, n: int) -> np.ndarray:
    """Apply P8/P10 xcorr masks without assuming whether 1 means keep or reject."""
    if mask.size < n:
        return valid
    m = mask[:n].astype(bool)
    keep_true = valid & m
    keep_false = valid & ~m
    if int(np.sum(keep_true)) >= int(np.sum(keep_false)):
        return keep_true
    return keep_false


def contiguous_true_runs(flag: np.ndarray) -> Iterable[tuple[int, int]]:
    flag = np.asarray(flag, dtype=bool)
    if flag.size == 0 or not np.any(flag):
        return []
    edges = np.flatnonzero(np.diff(np.r_[False, flag, False]))
    return [(int(start), int(stop)) for start, stop in zip(edges[0::2], edges[1::2])]


def annotate_trace_mask_and_borders(
    ax: plt.Axes,
    t_abs: np.ndarray,
    mask: np.ndarray,
    time_unit: str,
    origin: float = 0.0,
    session_borders: np.ndarray | None = None,
) -> None:
    """Shade masked samples and draw session borders for trace comparison plots.

    t_abs is always in absolute concatenated-session seconds. The visible x-axis
    is converted here with the same time_unit/origin used for plotting, so this
    helper can be reused for both whole-trace and zoom-window panels.
    """
    t_abs = np.asarray(t_abs, dtype=float)
    n = int(t_abs.size)
    if n == 0:
        return
    x, _ = time_axis(t_abs, time_unit, origin=origin)
    finite_t = np.isfinite(t_abs) & np.isfinite(x)
    valid = finite_t.copy()
    valid = apply_xcorr_valid_mask(valid, np.asarray(mask), n)
    masked = finite_t & ~valid
    dx = 0.0
    finite_x = x[np.isfinite(x)]
    if finite_x.size > 1:
        diffs = np.diff(finite_x)
        diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
        if diffs.size:
            dx = float(np.nanmedian(diffs))
    pad = 0.5 * dx
    for start, stop in contiguous_true_runs(masked):
        left = float(x[start] - pad)
        right = float(x[stop - 1] + pad)
        ax.axvspan(left, right, color="#9ca3af", alpha=0.18, lw=0, zorder=-10)
    if session_borders is None:
        return
    borders = np.asarray(session_borders, dtype=float)
    borders = borders[np.isfinite(borders)]
    if borders.size == 0:
        return
    t_min = float(np.nanmin(t_abs[finite_t])) if np.any(finite_t) else float(np.nanmin(t_abs))
    t_max = float(np.nanmax(t_abs[finite_t])) if np.any(finite_t) else float(np.nanmax(t_abs))
    bx, _ = time_axis(borders[(borders >= t_min) & (borders <= t_max)], time_unit, origin=origin)
    for b in bx:
        ax.axvline(float(b), color="#6b7280", lw=0.8, ls="--", alpha=0.65, zorder=-5)


def select_zoom_windows(
    t: np.ndarray,
    den_z: np.ndarray,
    bold_z: np.ndarray,
    valid_mask: np.ndarray,
    window_sec: float,
    top_k: int,
) -> list[dict[str, float]]:
    if t.size < 4:
        return []
    dt = float(np.nanmedian(np.diff(t)))
    if not math.isfinite(dt) or dt <= 0:
        dt = 1.0
    window_bins = max(8, int(round(window_sec / dt)))
    window_bins = min(window_bins, t.size)
    step = max(1, window_bins // 4)
    candidates: list[dict[str, float]] = []
    for start in range(0, max(1, t.size - window_bins + 1), step):
        stop = min(t.size, start + window_bins)
        ok = (
            valid_mask[start:stop]
            & np.isfinite(den_z[start:stop])
            & np.isfinite(bold_z[start:stop])
        )
        if int(np.sum(ok)) < max(6, window_bins // 4):
            continue
        x = den_z[start:stop][ok]
        y = bold_z[start:stop][ok]
        if np.nanstd(x) <= 1e-9 or np.nanstd(y) <= 1e-9:
            continue
        corr = float(np.corrcoef(x, y)[0, 1])
        if not math.isfinite(corr):
            continue
        candidates.append(
            {
                "start": float(start),
                "stop": float(stop),
                "score": abs(corr),
                "corr": corr,
                "start_sec": float(t[start]),
                "stop_sec": float(t[stop - 1]),
            }
        )
    candidates.sort(key=lambda item: item["score"], reverse=True)
    selected: list[dict[str, float]] = []
    for cand in candidates:
        start = int(cand["start"])
        stop = int(cand["stop"])
        overlaps = False
        for prev in selected:
            pstart = int(prev["start"])
            pstop = int(prev["stop"])
            overlap = max(0, min(stop, pstop) - max(start, pstart))
            if overlap > 0.25 * window_bins:
                overlaps = True
                break
        if overlaps:
            continue
        selected.append(cand)
        if len(selected) >= top_k:
            break
    if selected:
        return selected
    stop = min(t.size, window_bins)
    return [
        {
            "start": 0.0,
            "stop": float(stop),
            "score": math.nan,
            "corr": math.nan,
            "start_sec": float(t[0]),
            "stop_sec": float(t[stop - 1]),
        }
    ]


def fixed_window_bounds(center: int, n: int, window_bins: int) -> tuple[int, int]:
    if n <= 0:
        return 0, 0
    window_bins = max(1, min(int(window_bins), int(n)))
    start = int(center) - window_bins // 2
    start = max(0, min(start, n - window_bins))
    stop = start + window_bins
    return start, stop


def finite_argmax_index(x: np.ndarray, valid: np.ndarray, fallback: int) -> int:
    score = np.where(valid & np.isfinite(x), x, -np.inf)
    if not np.any(np.isfinite(score)):
        return int(fallback)
    return int(np.nanargmax(score))


def select_fixed_summary_windows(
    t: np.ndarray,
    den_z: np.ndarray,
    bold_z: np.ndarray,
    valid_mask: np.ndarray,
    window_sec: float,
) -> list[dict[str, float | str]]:
    if t.size < 1:
        return []
    dt = float(np.nanmedian(np.diff(t))) if t.size > 1 else 1.0
    if not math.isfinite(dt) or dt <= 0:
        dt = 1.0
    window_bins = max(8, int(round(window_sec / dt)))
    window_bins = min(window_bins, t.size)
    valid_idx = np.flatnonzero(valid_mask & np.isfinite(den_z) & np.isfinite(bold_z))
    middle = int(valid_idx[len(valid_idx) // 2]) if valid_idx.size else int(t.size // 2)
    picks = [
        ("raw_density_high", finite_argmax_index(den_z, valid_mask, middle)),
        ("bold_feature_high", finite_argmax_index(bold_z, valid_mask, middle)),
        ("middle_representative", middle),
    ]
    windows: list[dict[str, float | str]] = []
    for label, center in picks:
        start, stop = fixed_window_bounds(center, t.size, window_bins)
        windows.append(
            {
                "label": label,
                "start": float(start),
                "stop": float(stop),
                "score": math.nan,
                "corr": math.nan,
                "start_sec": float(t[start]),
                "stop_sec": float(t[stop - 1]),
                "center_sec": float(t[center]),
            }
        )
    return windows


def plot_zoom_windows(
    selections: dict[str, dict[int, list[dict[str, object]]]],
    figure_dir: Path,
    top_n: int,
    pipeline: str,
    max_datasets: int,
    family_filter: str,
    window_sec: float,
    top_k: int,
    align_xcorr_sign: bool,
    zoom_window_mode: str,
    zoom_y_scale: str,
    dataset_order: Sequence[str] = DATASET_ORDER,
) -> list[dict[str, object]]:
    rows = choose_trace_rows(selections, top_n, pipeline, dataset_order)
    rows = [
        row
        for row in rows
        if order_index(str(row.get("dataset")), dataset_order) < max_datasets
        and (family_filter == "both" or str(row.get("bold_feature_family")) == family_filter)
    ]
    cache = TraceCache()
    zoom_records: list[dict[str, object]] = []
    if not rows:
        return zoom_records

    n_rows = len(rows)
    n_cols = 3 if zoom_window_mode == "fixed_summary" else max(1, int(top_k))
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(7.2 * n_cols, max(2.3 * n_rows, 5)),
        sharex=True,
        sharey=False,
    )
    axes = np.asarray(axes)
    if axes.ndim == 1:
        axes = axes.reshape(n_rows, n_cols)

    for i_row, row in enumerate(rows):
        try:
            xcorr_mat, bold_post = read_xcorr_bold_post_file(str(row.get("source_csv", "")))
            if bold_post is None:
                raise FileNotFoundError("BOLD_POST path not resolved")
            density_index = as_int(row.get("density_index"))
            bold_mode_index = as_int(row.get("bold_mode_index"))
            t_den, den = cache.density_trace(str(row.get("raw_density_file", "")), density_index)
            t_bold, bold = cache.bold_trace(bold_post, str(row.get("bold_feature", "")), bold_mode_index)
            n = min(den.size, bold.size, t_den.size, t_bold.size)
            t = t_den[:n]
            den_z = normalize_trace(den[:n])
            bold_z = normalize_trace(bold[:n])
            align_factor = xcorr_alignment_factor(row, align_xcorr_sign)
            bold_z = align_factor * bold_z
            row_ylim = combined_ylim([den_z, bold_z]) if zoom_y_scale == "full_trace" else None
            mask = cache.border_mask(xcorr_mat)
            borders = cache.session_border_times(bold_post)
            valid = np.isfinite(den_z) & np.isfinite(bold_z)
            valid = apply_xcorr_valid_mask(valid, mask, n)
            if int(np.sum(valid)) < max(6, min(n, int(round(window_sec))) // 4):
                valid = np.isfinite(den_z) & np.isfinite(bold_z)
            if zoom_window_mode == "fixed_summary":
                windows = select_fixed_summary_windows(t, den_z, bold_z, valid, window_sec)
            else:
                windows = select_zoom_windows(t, den_z, bold_z, valid, window_sec, top_k)
            for j_col in range(n_cols):
                ax = axes[i_row, j_col]
                if j_col >= len(windows):
                    ax.axis("off")
                    continue
                win = windows[j_col]
                lo = int(win["start"])
                hi = int(win["stop"])
                tx, unit_label = time_axis(t[lo:hi], "min", origin=float(t[lo]))
                segment_mask = mask[lo:hi] if mask.size >= hi else np.asarray([], dtype=bool)
                annotate_trace_mask_and_borders(
                    ax,
                    t[lo:hi],
                    segment_mask,
                    "min",
                    origin=float(t[lo]),
                    session_borders=borders,
                )
                ax.plot(tx, den_z[lo:hi], color="#7c3aed", lw=0.9, alpha=0.85, label="raw density")
                ax.plot(
                    tx,
                    bold_z[lo:hi],
                    color=color_for_family(str(row.get("bold_feature_family"))),
                    lw=1.0,
                    alpha=0.9,
                    label=str(row.get("bold_feature")),
                )
                if row_ylim is not None:
                    ax.set_ylim(*row_ylim)
                ax.grid(True, color="#e5e7eb", lw=0.5)
                ax.set_xlim(0, window_sec / 60.0)
                if j_col == 0:
                    ax.set_ylabel(f"{row.get('dataset')}\nrobust z", fontsize=8)
                if i_row == n_rows - 1:
                    ax.set_xlabel(f"time within selected whole-trace window ({unit_label})")
                raw_idx = as_int(row.get("raw_efun_index"))
                tau = as_float(row.get("raw_timescale_sec"))
                signed_peak = as_float(row.get("peak_corr"))
                peak_abs = as_float(row.get("peak_abs_corr"))
                lag = as_float(row.get("peak_lag_sec"))
                start_min = float(win["start_sec"]) / 60.0
                window_label = str(win.get("label", f"local_corr_top{j_col + 1}"))
                ax.set_title(
                    f"{window_label} | {row.get('bold_feature_family')} | {row.get('bold_observable')} | "
                    f"raw {raw_idx} tau={tau:.2g}s | |r|={peak_abs:.2f} r={signed_peak:.2f} "
                    f"align={align_factor:+.0f} lag={lag:.1f}s | start={start_min:.1f}m",
                    fontsize=8,
                )
                zoom_records.append(
                    dict(
                        row,
                        zoom_rank=j_col + 1,
                        zoom_window_mode=zoom_window_mode,
                        zoom_y_scale=zoom_y_scale,
                        zoom_window_label=window_label,
                        zoom_start_sec=win["start_sec"],
                        zoom_stop_sec=win["stop_sec"],
                        zoom_center_sec=win.get("center_sec", ""),
                        zoom_local_corr=win["corr"],
                        zoom_score=win["score"],
                        bold_trace_alignment_factor=float(align_factor),
                    )
                )
        except Exception as exc:
            for j_col in range(n_cols):
                ax = axes[i_row, j_col]
                ax.text(0.5, 0.5, f"zoom error:\n{exc}", ha="center", va="center", transform=ax.transAxes, fontsize=8)
                ax.set_title(f"{row.get('dataset')} | {row.get('bold_feature_family')}", fontsize=8)
            zoom_records.append(dict(row, zoom_error=str(exc)))
    plotted_families = sorted({str(row.get("bold_feature_family")) for row in rows})
    add_trace_legend(fig, plotted_families)
    family_slug = family_filter.replace("_", "")
    fig.suptitle(
        f"{pipeline} {family_filter} real/raw-density zoom windows | raw-top{top_n} | {window_sec / 60:.1f} min {zoom_window_mode} windows | y={zoom_y_scale}"
        + (" | BOLD trace xcorr-sign aligned" if align_xcorr_sign else ""),
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    out_file = figure_dir / f"06_{pipeline.lower()}_{family_slug}_real_vs_raw_density_zoom_windows_{zoom_window_mode}_yscale_{zoom_y_scale}_raw_top{top_n}_{int(window_sec)}s.png"
    fig.savefig(out_file, dpi=220)
    plt.close(fig)
    return zoom_records


def plot_whole_trace_examples(
    selections: dict[str, dict[int, list[dict[str, object]]]],
    figure_dir: Path,
    top_n: int,
    pipeline: str,
    max_datasets: int,
    fig_width: float,
    row_height: float,
    time_unit: str,
    align_xcorr_sign: bool,
    dataset_order: Sequence[str] = DATASET_ORDER,
) -> list[dict[str, object]]:
    rows = choose_trace_rows(selections, top_n, pipeline, dataset_order)
    rows = [
        row
        for row in rows
        if order_index(str(row.get("dataset")), dataset_order) < max_datasets
    ]
    cache = TraceCache()
    trace_records: list[dict[str, object]] = []
    if not rows:
        return trace_records

    visible_datasets = list(dataset_order[:max_datasets])
    n_rows = len(visible_datasets)
    fig, axes = plt.subplots(
        n_rows,
        2,
        figsize=(fig_width, max(row_height * n_rows, 6)),
        sharex=False,
        sharey=False,
    )
    if n_rows == 1:
        axes = np.asarray([axes])
    for ax in axes.ravel():
        ax.axis("off")

    row_lookup = {(str(row.get("dataset")), str(row.get("bold_feature_family"))): row for row in rows}
    for i_ds, dataset in enumerate(visible_datasets):
        for j_family, family in enumerate(FAMILY_ORDER):
            ax = axes[i_ds, j_family]
            ax.axis("on")
            row = row_lookup.get((dataset, family))
            if row is None:
                ax.text(0.5, 0.5, "missing", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{dataset} | {family}", fontsize=8)
                continue
            try:
                xcorr_mat, bold_post = read_xcorr_bold_post_file(str(row.get("source_csv", "")))
                if bold_post is None:
                    raise FileNotFoundError("BOLD_POST path not resolved")
                density_index = as_int(row.get("density_index"))
                bold_mode_index = as_int(row.get("bold_mode_index"))
                t_den, den = cache.density_trace(str(row.get("raw_density_file", "")), density_index)
                t_bold, bold = cache.bold_trace(bold_post, str(row.get("bold_feature", "")), bold_mode_index)
                n = min(den.size, bold.size, t_den.size, t_bold.size)
                den = den[:n]
                bold = bold[:n]
                t = t_den[:n]
                x, unit_label = time_axis(t, time_unit)
                x_label = f"time in concatenated session ({unit_label})"
                den_z = normalize_trace(den)
                bold_z = normalize_trace(bold)
                align_factor = xcorr_alignment_factor(row, align_xcorr_sign)
                bold_z = align_factor * bold_z
                mask = cache.border_mask(xcorr_mat)
                borders = cache.session_border_times(bold_post)
                annotate_trace_mask_and_borders(
                    ax,
                    t,
                    mask,
                    time_unit,
                    origin=0.0,
                    session_borders=borders,
                )
                ax.plot(x, den_z, color="#7c3aed", lw=0.7, alpha=0.75, label="raw density")
                ax.plot(x, bold_z, color=color_for_family(family), lw=1.2, alpha=0.92, label=str(row.get("bold_feature")))
                tau = as_float(row.get("raw_timescale_sec"))
                peak = as_float(row.get("peak_abs_corr"))
                signed_peak = as_float(row.get("peak_corr"))
                raw_idx = as_int(row.get("raw_efun_index"))
                mode_idx = as_int(row.get("bold_mode_index"))
                lag = as_float(row.get("peak_lag_sec"))
                ax.set_title(
                    f"{dataset} | {family} | {row.get('bold_observable')} | raw {raw_idx} tau={tau:.2g}s | BOLD mode {mode_idx} | |r|={peak:.3f} r={signed_peak:.3f} align={align_factor:+.0f} lag={lag:.1f}s",
                    fontsize=8,
                )
                ax.grid(True, color="#e5e7eb", lw=0.5)
                if i_ds == n_rows - 1:
                    ax.set_xlabel(x_label)
                if j_family == 0:
                    ax.set_ylabel("robust z")
                trace_records.append(
                    dict(
                        row,
                        xcorr_mat=str(xcorr_mat),
                        bold_post_file=str(bold_post),
                        trace_time_start_sec=float(t[0]),
                        trace_time_end_sec=float(t[n - 1]),
                        trace_n_samples=int(n),
                        trace_kind="whole_trace",
                        bold_trace_alignment_factor=float(align_factor),
                    )
                )
            except Exception as exc:
                ax.text(0.5, 0.5, f"trace error:\n{exc}", ha="center", va="center", transform=ax.transAxes, fontsize=8)
                ax.set_title(f"{dataset} | {family}", fontsize=8)
                trace_records.append(dict(row, trace_error=str(exc), trace_kind="whole_trace"))
    add_trace_legend(fig, FAMILY_ORDER)
    fig.suptitle(
        f"{pipeline} BOLD feature vs top raw BLP density whole traces | raw-top{top_n} | x-axis is original concatenated time"
        + (" | BOLD trace xcorr-sign aligned" if align_xcorr_sign else ""),
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(figure_dir / f"05_{pipeline.lower()}_bold_vs_raw_density_whole_trace_examples_raw_top{top_n}.png", dpi=220)
    fig.savefig(figure_dir / f"05_{pipeline.lower()}_bold_vs_raw_density_whole_trace_examples_raw_top{top_n}__wide.png", dpi=220)
    plt.close(fig)
    return trace_records


def audit_current_map_availability(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    audit_rows: list[dict[str, object]] = []
    seen: set[tuple[str, str, str, str]] = set()
    for row in rows:
        dataset = str(row.get("dataset", ""))
        pipeline = str(row.get("pipeline", ""))
        run_tag = str(row.get("run_tag", ""))
        family = str(row.get("bold_feature_family", ""))
        key = (pipeline, dataset, run_tag, family)
        if key in seen:
            continue
        seen.add(key)
        if pipeline == "P8":
            root = win_to_wsl(Path("E:/DataPons_processed") / dataset / "pipeline8_top_maps" / run_tag)
        else:
            root = win_to_wsl(Path("E:/DataPons_processed") / dataset / "pipeline10_dimred_top_maps")
        current_matches: list[Path] = []
        if root.is_dir():
            tokens = ("standardize", "rmsenv", DEFAULT_DENSITY)
            for path in root.rglob("*.png"):
                name = path.name.lower()
                if any(str(token).lower() in name for token in tokens):
                    current_matches.append(path)
                    if len(current_matches) >= 5:
                        break
        audit_rows.append(
            {
                "pipeline": pipeline,
                "dataset": dataset,
                "run_tag": run_tag,
                "bold_feature_family": family,
                "map_root": wsl_to_win(root),
                "current_standardized_map_found": bool(current_matches),
                "n_example_matches": len(current_matches),
                "example_map_1": wsl_to_win(current_matches[0]) if current_matches else "",
                "note": "Current standardized csplit map figures are expected to be missing if P8/P10 was run numeric-first with maps disabled.",
            }
        )
    return audit_rows


def main() -> None:
    args = parse_args()
    results_dir: Path = args.results_dir
    figure_dir: Path = args.figure_dir or (results_dir / "figures")
    results_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    excluded_datasets = {str(name).lower() for name in args.exclude_datasets}
    dataset_order = tuple(ds for ds in DATASET_ORDER if ds.lower() not in excluded_datasets)
    metadata = load_metadata(args.raw_metadata_csv, args.density_name)
    all_rows = [
        enrich_hit(row, metadata)
        for row in read_csv(args.hits_csv)
        if row.get("dataset", "").lower() not in excluded_datasets
        if row.get("density_name") == args.density_name
        or row.get("density_condition") in {"event", "csplit"}
    ]
    selections = build_ranked_selections(
        all_rows,
        args.top_n_values,
        args.density_name,
        args.bold_feature_filter,
    )
    write_selection_csv(results_dir / "bold_top_raw_density_hits_with_timescale_long.csv", selections)
    summary_rows = summarize_timescales(selections)
    write_csv(results_dir / "bold_top_raw_density_timescale_summary_by_topN.csv", summary_rows)

    plot_ecdf(selections, figure_dir)
    plot_median_by_topn(summary_rows, figure_dir)
    plot_slow_fraction(summary_rows, figure_dir)
    plot_index_hist(selections, figure_dir, top_n=args.trace_top_n)

    trace_records: list[dict[str, object]] = []
    if not args.skip_traces:
        trace_records = plot_whole_trace_examples(
            selections,
            figure_dir,
            args.trace_top_n,
            args.trace_pipeline,
            args.trace_max_datasets,
            args.whole_trace_fig_width,
            args.whole_trace_row_height,
            args.whole_trace_time_unit,
            args.align_xcorr_sign,
            dataset_order,
        )
        write_csv(results_dir / "selected_whole_trace_examples.csv", trace_records)
        zoom_records = plot_zoom_windows(
            selections,
            figure_dir,
            args.trace_top_n,
            args.trace_pipeline,
            args.trace_max_datasets,
            args.zoom_family,
            args.zoom_window_sec,
            args.zoom_top_k,
            args.align_xcorr_sign,
            args.zoom_window_mode,
            args.zoom_y_scale,
            dataset_order,
        )
        write_csv(results_dir / "selected_zoom_window_examples.csv", zoom_records)
        if args.legacy_local_peak_traces:
            local_trace_records = plot_trace_examples(
                selections,
                figure_dir,
                args.trace_top_n,
                args.trace_pipeline,
                args.trace_max_datasets,
                args.trace_half_window_bins,
                args.align_xcorr_sign,
                dataset_order,
            )
            write_csv(results_dir / "selected_local_peak_trace_examples_legacy.csv", local_trace_records)

    map_audit_rows = audit_current_map_availability(selections["raw_topN"].get(args.trace_top_n, []))
    write_csv(results_dir / "current_roi_activation_map_availability.csv", map_audit_rows)

    manifest_rows = [
        {
            "result": "exclude_datasets",
            "path": " ".join(args.exclude_datasets),
        },
        {
            "result": "bold_feature_filter",
            "path": args.bold_feature_filter,
        },
        {
            "result": "align_xcorr_sign",
            "path": str(bool(args.align_xcorr_sign)),
        },
        {
            "result": "zoom_window_mode",
            "path": args.zoom_window_mode,
        },
        {
            "result": "zoom_y_scale",
            "path": args.zoom_y_scale,
        },
        {
            "result": "hits_with_timescale",
            "path": str((results_dir / "bold_top_raw_density_hits_with_timescale_long.csv").resolve()),
        },
        {
            "result": "timescale_summary",
            "path": str((results_dir / "bold_top_raw_density_timescale_summary_by_topN.csv").resolve()),
        },
        {"result": "figure_dir", "path": str(figure_dir.resolve())},
        {"result": "whole_trace_examples", "path": str((results_dir / "selected_whole_trace_examples.csv").resolve())},
        {"result": "zoom_window_examples", "path": str((results_dir / "selected_zoom_window_examples.csv").resolve())},
        {"result": "map_availability", "path": str((results_dir / "current_roi_activation_map_availability.csv").resolve())},
    ]
    write_csv(results_dir / "manifest.csv", manifest_rows)
    print(f"Wrote {results_dir}")
    print(f"Wrote figures {figure_dir}")
    print(f"Trace examples: {len(trace_records)}")


if __name__ == "__main__":
    main()
