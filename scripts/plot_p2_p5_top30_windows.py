#!/usr/bin/env python3
"""Plot P2 and P5 top-30 windows for one standardized P5 condition.

This is a lightweight diagnostic bundle for the "new dataset" SOP:

1. P2 consensus-state-diversity top windows:
   show the P2 theta/gamma/ripple event masks inside the saved top-30 windows.

2. P5 dimred component top windows:
   choose representative dimred components from the P2-band response label CSV,
   then show each component's highest-activity 6000-sample global windows with
   the same P2 event masks overlaid.

The P5 windows here are component-activity top windows. They are intentionally
not an overlap/rank statistic against P2 top windows.
"""

from __future__ import annotations

import argparse
import csv
import math
import warnings
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BAND_ORDER = ["theta", "gamma", "ripple"]
BAND_COLORS = {
    "theta": "#2ca25f",
    "gamma": "#f0b429",
    "ripple": "#3182bd",
}
DEFAULT_CONDITION = "complex_split_projected_vlambda_standardize"
DEFAULT_LABEL_SUFFIX = "p5_p2_band_event_response_v2_selective_envelope_20260528"


@dataclass
class BandEvents:
    label: str
    mask: np.ndarray
    centers: np.ndarray
    passband_hz: tuple[float, float]
    n_raw_channel_events: int
    coverage_fraction: float


@dataclass
class WindowRank:
    rank: int
    global_window_idx: int
    start0: int
    end0: int
    score: float
    session_id: int
    crosses_session_boundary: bool


def main() -> None:
    args = parse_args()

    workspace = Path(args.workspace)
    processed_root = Path(args.processed_root)
    dataset_root = processed_root / args.dataset
    p2_out = dataset_root / "pipeline2_figures_consensus_state_top30_event_raster"
    p5_out = (
        dataset_root
        / "pipeline5_component_activity_top30_windows"
        / args.condition
    )
    p2_out.mkdir(parents=True, exist_ok=True)
    p5_out.mkdir(parents=True, exist_ok=True)

    bands, dt, session_start, session_end, n_time = load_p2_events(dataset_root, args.dataset)
    top_table = load_p2_top_window_table(dataset_root, args.dataset)
    top_table = top_table.head(args.top_n_windows).copy()

    p2_manifest = plot_p2_top_windows(top_table, bands, dt, p2_out, args.dataset)
    label_csv = workspace / "results" / f"{args.dataset}_{DEFAULT_LABEL_SUFFIX}" / "component_band_event_response_v2.csv"
    if not label_csv.exists():
        raise FileNotFoundError(f"P5/P2 label CSV not found: {label_csv}")

    selected = select_p5_components(
        pd.read_csv(label_csv),
        n_per_group=args.n_components_per_group,
        transforms=args.transforms,
    )
    selected_csv = p5_out / "selected_components_for_top30_windows.csv"
    selected.to_csv(selected_csv, index=False)

    p5_manifest_rows: list[dict[str, object]] = []
    for row in selected.to_dict("records"):
        activity = load_component_activity(row, session_start, session_end)
        ranks = rank_component_windows(
            activity,
            n_time=n_time,
            session_start=session_start,
            session_end=session_end,
            window_len=args.window_length_samples,
            top_n=args.top_n_windows,
            include_boundary_windows=args.include_boundary_windows,
        )
        fig_file = plot_component_top_windows(
            row,
            activity,
            ranks,
            bands,
            dt,
            p5_out,
            args.dataset,
        )
        summary = summarize_ranked_window_events(ranks, bands)
        p5_manifest_rows.append(
            {
                "dataset": args.dataset,
                "condition": args.condition,
                "selection_group": row["selection_group"],
                "activity_transform": row["activity_transform"],
                "method": row["method"],
                "k": int(row["k"]),
                "component": int(row["component"]),
                "strict_label": row["strict_label"],
                "theta_effect_z": float(row["theta_effect_z"]),
                "gamma_effect_z": float(row["gamma_effect_z"]),
                "ripple_effect_z": float(row["ripple_effect_z"]),
                "n_ranked_windows": len(ranks),
                "theta_top30_mean_fraction": summary["theta"],
                "gamma_top30_mean_fraction": summary["gamma"],
                "ripple_top30_mean_fraction": summary["ripple"],
                "figure_file": str(fig_file),
                "source_file": row["source_file"],
            }
        )

    p5_manifest = p5_out / "component_activity_top30_window_manifest.csv"
    with p5_manifest.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(p5_manifest_rows[0].keys()))
        writer.writeheader()
        writer.writerows(p5_manifest_rows)

    plot_p5_top_window_event_fraction_summary(pd.DataFrame(p5_manifest_rows), p5_out, args.dataset)

    print("P2 top30 event-raster manifest:")
    print(p2_manifest)
    print("P5 selected components:")
    print(selected_csv)
    print("P5 component top30 manifest:")
    print(p5_manifest)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot P2 and P5 top30 diagnostic windows.")
    parser.add_argument("--dataset", default="k13m17")
    parser.add_argument("--condition", default=DEFAULT_CONDITION)
    parser.add_argument("--workspace", default="/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
    parser.add_argument("--processed-root", default="/mnt/e/DataPons_processed")
    parser.add_argument("--window-length-samples", type=int, default=6000)
    parser.add_argument("--top-n-windows", type=int, default=30)
    parser.add_argument("--n-components-per-group", type=int, default=4)
    parser.add_argument(
        "--transforms",
        nargs="+",
        default=["abs", "adaptive_envelope"],
        choices=["abs", "adaptive_envelope"],
    )
    parser.add_argument("--include-boundary-windows", action="store_true")
    return parser.parse_args()


def load_p2_events(dataset_root: Path, dataset: str) -> tuple[dict[str, BandEvents], float, np.ndarray, np.ndarray, int]:
    event_file = dataset_root / "pipeline2_event_detection" / f"{dataset}_bandpass_events_3bands.mat"
    if not event_file.exists():
        raise FileNotFoundError(f"P2 event file not found: {event_file}")

    with h5py.File(event_file, "r") as f:
        detect = f["R/DetectResults"]
        n_time = int(np.asarray(f["R/session_end_idx"]).squeeze()[-1])
        dt = float(np.asarray(f["R/dx"]).squeeze())
        session_start = np.asarray(f["R/session_start_idx"]).squeeze().astype(int) - 1
        session_end = np.asarray(f["R/session_end_idx"]).squeeze().astype(int)
        bands: dict[str, BandEvents] = {}

        for band_col in range(detect.shape[1]):
            mask = np.zeros(n_time, dtype=bool)
            centers_all: list[np.ndarray] = []
            label = ""
            passband = (math.nan, math.nan)
            n_raw = 0

            for ch in range(detect.shape[0]):
                group = f[detect[ch, band_col]]
                label = decode_matlab_char(group["band_label"])
                hz = np.asarray(group["bandpass_hz"]).squeeze().astype(float)
                passband = (float(hz[0]), float(hz[1]))
                wins = np.asarray(group["event_win"], dtype=float).T
                peaks = np.asarray(group["loc_peak"], dtype=float).squeeze()
                n_raw += wins.shape[0]

                if wins.size:
                    starts = np.maximum(0, wins[:, 0].astype(int) - 1)
                    ends = np.minimum(n_time, wins[:, 1].astype(int))
                    for start, end in zip(starts, ends):
                        if end > start:
                            mask[start:end] = True
                if peaks.size:
                    centers_all.append(peaks.astype(int) - 1)

            centers = np.concatenate(centers_all) if centers_all else np.zeros(0, dtype=int)
            bands[label] = BandEvents(
                label=label,
                mask=mask,
                centers=np.unique(centers.astype(int)),
                passband_hz=passband,
                n_raw_channel_events=int(n_raw),
                coverage_fraction=float(mask.mean()),
            )

    return bands, dt, session_start, session_end, n_time


def decode_matlab_char(dataset) -> str:
    arr = np.asarray(dataset).squeeze()
    return "".join(chr(int(x)) for x in arr if int(x) > 0)


def load_p2_top_window_table(dataset_root: Path, dataset: str) -> pd.DataFrame:
    top_csv = (
        dataset_root
        / "pipeline2_consensus_state_diversity_windows"
        / f"{dataset}_consensus_state_diversity_windows_6000samp_globalwin_top.csv"
    )
    if not top_csv.exists():
        raise FileNotFoundError(f"P2 top-window CSV not found: {top_csv}")
    return pd.read_csv(top_csv)


def plot_p2_top_windows(
    top_table: pd.DataFrame,
    bands: dict[str, BandEvents],
    dt: float,
    out_dir: Path,
    dataset: str,
) -> Path:
    overview = out_dir / f"{dataset}_p2_consensus_state_diversity_top30_event_raster.png"
    plot_window_grid(
        top_table,
        bands,
        dt,
        overview,
        title=f"{dataset}: P2 consensus-state diversity top30 windows",
        score_column="normalized_state_entropy",
    )

    detail_dir = out_dir / "individual_windows"
    detail_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    for _, row in top_table.iterrows():
        rank = int(row["state_diversity_rank"])
        fig_file = detail_dir / f"p2_top_window_rank{rank:02d}_global{int(row['global_window_idx']):05d}.png"
        plot_single_window(
            int(row["global_start_idx"]) - 1,
            int(row["global_end_idx"]),
            bands,
            dt,
            fig_file,
            title=(
                f"{dataset} P2 rank {rank:02d} | global window {int(row['global_window_idx'])} | "
                f"dominant={row.get('dominant_state', '')} | entropy={float(row.get('normalized_state_entropy', np.nan)):.3f}"
            ),
        )
        manifest_rows.append(
            {
                "state_diversity_rank": rank,
                "global_window_idx": int(row["global_window_idx"]),
                "global_start_idx": int(row["global_start_idx"]),
                "global_end_idx": int(row["global_end_idx"]),
                "dominant_state": row.get("dominant_state", ""),
                "normalized_state_entropy": float(row.get("normalized_state_entropy", np.nan)),
                "figure_file": str(fig_file),
            }
        )

    manifest = out_dir / "p2_top30_event_raster_manifest.csv"
    pd.DataFrame(manifest_rows).to_csv(manifest, index=False)
    return manifest


def plot_window_grid(
    table: pd.DataFrame,
    bands: dict[str, BandEvents],
    dt: float,
    out_file: Path,
    title: str,
    score_column: str,
) -> None:
    n = len(table)
    n_cols = 5
    n_rows = int(math.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 2.2 * n_rows), sharex=False, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax_idx, ax in enumerate(axes):
        if ax_idx >= n:
            ax.axis("off")
            continue
        row = table.iloc[ax_idx]
        start0 = int(row["global_start_idx"]) - 1
        end0 = int(row["global_end_idx"])
        draw_event_raster(ax, start0, end0, bands, dt)
        score = float(row.get(score_column, np.nan))
        ax.set_title(
            f"#{ax_idx + 1} gw{int(row['global_window_idx'])} | {score_column}={score:.2f}",
            fontsize=8,
        )
    fig.suptitle(title, fontsize=16, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out_file, dpi=180)
    plt.close(fig)


def plot_single_window(
    start0: int,
    end0: int,
    bands: dict[str, BandEvents],
    dt: float,
    out_file: Path,
    title: str,
    activity: np.ndarray | None = None,
    activity_label: str = "",
) -> None:
    fig, ax = plt.subplots(figsize=(12, 3.2 if activity is None else 4.2))
    if activity is not None:
        x = np.arange(end0 - start0) * dt
        y = robust_scale_trace(activity[start0:end0])
        ax.plot(x, 2.2 + y, color="#252525", linewidth=0.8)
        ax.text(0.0, 3.3, activity_label, ha="left", va="center", fontsize=8, color="#252525")
        raster_offset = 0.0
        ax.set_ylim(get_window_ylim(y, raster_top=2.8))
    else:
        raster_offset = 0.0
        ax.set_ylim(-0.4, 2.8)
    draw_event_raster(ax, start0, end0, bands, dt, y_offset=raster_offset)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("time within window (s)")
    fig.tight_layout()
    fig.savefig(out_file, dpi=170)
    plt.close(fig)


def draw_event_raster(
    ax,
    start0: int,
    end0: int,
    bands: dict[str, BandEvents],
    dt: float,
    y_offset: float = 0.0,
) -> None:
    duration = (end0 - start0) * dt
    for y, band_name in enumerate(BAND_ORDER):
        band = bands[band_name]
        mask = band.mask[start0:end0]
        y0 = y_offset + y
        for seg_start, seg_end in mask_segments(mask):
            ax.axvspan(seg_start * dt, seg_end * dt, ymin=(y0 + 0.08) / 3.8, ymax=(y0 + 0.82) / 3.8, color=BAND_COLORS[band_name], alpha=0.85, lw=0)
        ax.text(-0.03 * duration, y0 + 0.45, band_name, ha="right", va="center", fontsize=7, color=BAND_COLORS[band_name])
    ax.set_xlim(0, duration)
    ax.set_yticks([])
    ax.grid(axis="x", color="0.9", linewidth=0.5)


def mask_segments(mask: np.ndarray) -> list[tuple[int, int]]:
    if mask.size == 0 or not mask.any():
        return []
    x = mask.astype(np.int8)
    edges = np.diff(np.r_[0, x, 0])
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)
    return list(zip(starts, ends))


def select_p5_components(df: pd.DataFrame, n_per_group: int, transforms: list[str]) -> pd.DataFrame:
    df = df[df["activity_transform"].isin(transforms)].copy()
    df["rg_effect_z"] = df[["gamma_effect_z", "ripple_effect_z"]].max(axis=1)
    groups = []
    for transform in transforms:
        dft = df[df["activity_transform"] == transform].copy()
        groups.append(
            pick_group(
                dft.sort_values("theta_effect_z", ascending=False),
                "theta_proxy_high_theta_effect",
                n_per_group,
            )
        )
        rg_df = dft[dft["strict_label"].isin(["ripple_gamma_no_theta", "gamma_selective", "ripple_selective"])]
        if rg_df.empty:
            rg_df = dft
        groups.append(
            pick_group(
                rg_df.sort_values("rg_effect_z", ascending=False),
                "ripple_gamma_candidate_high_rg_effect",
                n_per_group,
            )
        )
    selected = pd.concat(groups, ignore_index=True)
    selected = selected.drop_duplicates(
        subset=["selection_group", "activity_transform", "method", "k", "component"],
        keep="first",
    )
    return selected


def pick_group(df: pd.DataFrame, name: str, n: int) -> pd.DataFrame:
    out = df.head(n).copy()
    out.insert(0, "selection_group", name)
    return out


def load_component_activity(row: dict[str, object], session_start: np.ndarray, session_end: np.ndarray) -> np.ndarray:
    source_file = str(row["source_file"])
    comp_idx = int(row["component"]) - 1
    with h5py.File(source_file, "r") as f:
        activity = np.abs(
            np.asarray(
                f["result/core/temporal_components_time_by_comp"][comp_idx, :],
                dtype=np.float32,
            )
        ).astype(np.float64)
    if row["activity_transform"] == "adaptive_envelope":
        activity = rms_envelope_sessionwise(
            activity,
            session_start,
            session_end,
            int(row["envelope_window_samples"]),
        )
    return activity


def rms_envelope_sessionwise(
    activity: np.ndarray,
    session_start: np.ndarray,
    session_end: np.ndarray,
    window_samples: int,
) -> np.ndarray:
    window_samples = max(1, int(window_samples))
    out = np.empty_like(activity, dtype=np.float64)
    for start, end in zip(session_start, session_end):
        start = int(start)
        end = int(end)
        if end <= start:
            continue
        out[start:end] = moving_rms(activity[start:end], window_samples)
    return out


def moving_rms(x: np.ndarray, window_samples: int) -> np.ndarray:
    if window_samples <= 1:
        return np.abs(x).astype(np.float64)
    left = window_samples // 2
    right = window_samples - 1 - left
    x2 = np.square(x.astype(np.float64))
    padded = np.pad(x2, (left, right), mode="edge")
    csum = np.concatenate(([0.0], np.cumsum(padded)))
    mean = (csum[window_samples:] - csum[:-window_samples]) / float(window_samples)
    return np.sqrt(np.maximum(mean, 0.0))


def rank_component_windows(
    activity: np.ndarray,
    n_time: int,
    session_start: np.ndarray,
    session_end: np.ndarray,
    window_len: int,
    top_n: int,
    include_boundary_windows: bool,
) -> list[WindowRank]:
    n = min(int(n_time), int(activity.size))
    n_windows = n // window_len
    trimmed = activity[: n_windows * window_len].reshape(n_windows, window_len)
    scores = trimmed.mean(axis=1)
    ranks = []
    order = np.argsort(scores)[::-1]
    for idx in order:
        start0 = int(idx) * window_len
        end0 = start0 + window_len
        session_id, crosses = resolve_session(start0, end0, session_start, session_end)
        if crosses and not include_boundary_windows:
            continue
        ranks.append(
            WindowRank(
                rank=len(ranks) + 1,
                global_window_idx=int(idx) + 1,
                start0=start0,
                end0=end0,
                score=float(scores[idx]),
                session_id=session_id,
                crosses_session_boundary=crosses,
            )
        )
        if len(ranks) >= top_n:
            break
    return ranks


def resolve_session(start0: int, end0: int, session_start: np.ndarray, session_end: np.ndarray) -> tuple[int, bool]:
    for i, (s0, e0) in enumerate(zip(session_start, session_end), start=1):
        if int(s0) <= start0 and end0 <= int(e0):
            return i, False
    return -1, True


def plot_component_top_windows(
    row: dict[str, object],
    activity: np.ndarray,
    ranks: list[WindowRank],
    bands: dict[str, BandEvents],
    dt: float,
    out_root: Path,
    dataset: str,
) -> Path:
    group = safe_name(str(row["selection_group"]))
    transform = safe_name(str(row["activity_transform"]))
    method = safe_name(str(row["method"]))
    k = int(row["k"])
    comp = int(row["component"])
    out_dir = out_root / transform / group
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{dataset}_{transform}_{group}_{method}_k{k:02d}_C{comp:02d}_top30_windows.png"

    n = len(ranks)
    n_cols = 5
    n_rows = int(math.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(18, 2.25 * n_rows), sharex=False, sharey=False)
    axes = np.atleast_1d(axes).ravel()
    for ax_idx, ax in enumerate(axes):
        if ax_idx >= n:
            ax.axis("off")
            continue
        r = ranks[ax_idx]
        x = np.arange(r.end0 - r.start0) * dt
        y = robust_scale_trace(activity[r.start0 : r.end0])
        ax.plot(x, 2.25 + y, color="#252525", linewidth=0.75)
        draw_event_raster(ax, r.start0, r.end0, bands, dt)
        ax.set_ylim(get_window_ylim(y, raster_top=2.8))
        ax.set_title(
            f"#{r.rank} gw{r.global_window_idx} s{r.session_id} score={r.score:.3g}",
            fontsize=8,
        )
    title = (
        f"{dataset} P5 component-activity top30 | {row['activity_transform']} | "
        f"{row['selection_group']} | {row['method']}_k{k:02d} C{comp}\n"
        f"strict={row['strict_label']} | effect_z T/G/R="
        f"{float(row['theta_effect_z']):.2f}/{float(row['gamma_effect_z']):.2f}/{float(row['ripple_effect_z']):.2f}"
    )
    fig.suptitle(title, fontsize=14, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_file, dpi=175)
    plt.close(fig)
    return out_file


def robust_scale_trace(x: np.ndarray) -> np.ndarray:
    lo, hi = np.nanpercentile(x, [2, 98])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(x, dtype=np.float64)
    y = (x - lo) / (hi - lo)
    return y


def get_window_ylim(y: np.ndarray, raster_top: float) -> tuple[float, float]:
    finite = y[np.isfinite(y)]
    if finite.size == 0:
        return (-0.4, raster_top + 0.9)
    ymin = min(-0.4, 2.25 + float(np.nanmin(finite)))
    ymax = max(raster_top + 0.9, 2.25 + float(np.nanmax(finite)))
    pad = 0.08 * max(1e-9, ymax - ymin)
    return (ymin - pad, ymax + pad)


def summarize_ranked_window_events(ranks: list[WindowRank], bands: dict[str, BandEvents]) -> dict[str, float]:
    summary = {}
    for band_name in BAND_ORDER:
        vals = []
        mask = bands[band_name].mask
        for r in ranks:
            vals.append(float(mask[r.start0 : r.end0].mean()))
        summary[band_name] = float(np.mean(vals)) if vals else math.nan
    return summary


def plot_p5_top_window_event_fraction_summary(manifest: pd.DataFrame, out_root: Path, dataset: str) -> None:
    if manifest.empty:
        return
    labels = [
        f"{r.activity_transform}\n{r.selection_group.replace('_', ' ')}\n{r.method}_k{int(r.k):02d}_C{int(r.component):02d}"
        for r in manifest.itertuples()
    ]
    x = np.arange(len(labels))
    width = 0.24
    fig, ax = plt.subplots(figsize=(max(12, 0.85 * len(labels)), 5.5))
    for j, band in enumerate(BAND_ORDER):
        values = manifest[f"{band}_top30_mean_fraction"].to_numpy(dtype=float)
        ax.bar(x + (j - 1) * width, values, width=width, color=BAND_COLORS[band], label=band)
    ax.set_xticks(x, labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("mean fraction of samples inside P5 top30 windows")
    ax.set_title(f"{dataset}: P2 event coverage inside selected P5 component top30 windows")
    ax.grid(axis="y", color="0.9")
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(out_root / f"{dataset}_p5_selected_component_top30_event_fraction_summary.png", dpi=180)
    plt.close(fig)


def safe_name(text: str) -> str:
    return "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in text)


if __name__ == "__main__":
    main()
