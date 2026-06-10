#!/usr/bin/env python3
"""One-plate P5 event diagnostic for a theta/RG component pair.

The plate combines two complementary checks for one method-k setting:

1. Peri-event mean activity for the selected theta-like and
   ripple-gamma-like components, aligned to P2 theta/gamma/ripple event peaks.
2. Top activity windows for the same two components, with P2 event masks
   overlaid as rugs.

The goal is visual diagnosis, not statistical testing.
"""

from __future__ import annotations

import argparse
import math
import warnings
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)
warnings.filterwarnings("ignore", message="This figure includes Axes that are not compatible with tight_layout*")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BAND_ORDER = ["theta", "gamma", "ripple"]
BAND_COLORS = {
    "theta": "#2ca25f",
    "gamma": "#f0b429",
    "ripple": "#3182bd",
}
TRACE_COLORS = {
    "theta": "#238b45",
    "ripple_gamma": "#6a51a3",
}
DEFAULT_CONDITION = "complex_split_projected_vlambda_standardize"
DEFAULT_LABEL_SUFFIX = "p5_p2_band_event_response_v2_selective_envelope_20260528"


@dataclass
class BandEvents:
    label: str
    mask: np.ndarray
    centers: np.ndarray


@dataclass
class PairSpec:
    activity_transform: str
    method: str
    k: int
    theta_component: int
    rg_component: int
    theta_effect_z: float
    gamma_effect_z_for_theta: float
    ripple_effect_z_for_theta: float
    rg_theta_effect_z: float
    rg_gamma_effect_z: float
    rg_ripple_effect_z: float
    theta_source_file: str
    rg_source_file: str
    theta_window_samples: int
    rg_window_samples: int


def main() -> None:
    args = parse_args()
    workspace = Path(args.workspace)
    processed_root = Path(args.processed_root)
    dataset_root = processed_root / args.dataset
    out_dir = dataset_root / "pipeline5_event_diagnostic_plates" / args.condition / args.transform
    out_dir.mkdir(parents=True, exist_ok=True)

    bands, dt, session_start, session_end, n_time = load_p2_events(dataset_root, args.dataset)
    label_csv = workspace / "results" / f"{args.dataset}_{DEFAULT_LABEL_SUFFIX}" / "component_band_event_response_v2.csv"
    labels = pd.read_csv(label_csv)
    pairs = choose_pairs(labels, args.transform, args.method, args.k, args.all_method_k)

    manifest_rows: list[dict[str, object]] = []
    if not pairs:
        fig_file = write_no_pair_figure(args, out_dir)
        manifest_rows.append(
            {
                "dataset": args.dataset,
                "condition": args.condition,
                "activity_transform": args.transform,
                "method": args.method,
                "k": args.k,
                "theta_component": "",
                "ripple_gamma_component": "",
                "theta_effect_z": "",
                "ripple_gamma_effect_z": "",
                "ribbon": args.ribbon,
                "status": "no_strict_theta_rg_pair",
                "figure_file": str(fig_file),
            }
        )
        print(fig_file)
    for pair in pairs:
        theta_activity = load_activity(
            pair.theta_source_file,
            pair.theta_component,
            pair.activity_transform,
            pair.theta_window_samples,
            session_start,
            session_end,
        )
        rg_activity = load_activity(
            pair.rg_source_file,
            pair.rg_component,
            pair.activity_transform,
            pair.rg_window_samples,
            session_start,
            session_end,
        )

        fig_file = plot_plate(
            args,
            pair,
            theta_activity,
            rg_activity,
            bands,
            dt,
            session_start,
            session_end,
            n_time,
            out_dir,
        )
        manifest_rows.append(
            {
                "dataset": args.dataset,
                "condition": args.condition,
                "activity_transform": pair.activity_transform,
                "method": pair.method,
                "k": pair.k,
                "theta_component": pair.theta_component,
                "ripple_gamma_component": pair.rg_component,
                "theta_effect_z": pair.theta_effect_z,
                "ripple_gamma_effect_z": max(pair.rg_gamma_effect_z, pair.rg_ripple_effect_z),
                "ribbon": args.ribbon,
                "status": "ok",
                "figure_file": str(fig_file),
            }
        )
        print(fig_file)

    manifest = out_dir / f"diagnostic_plate_manifest_{args.ribbon}.csv"
    pd.DataFrame(manifest_rows).to_csv(manifest, index=False)
    print(f"Manifest: {manifest}")


def write_no_pair_figure(args: argparse.Namespace, out_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(14, 3.2))
    ax.axis("off")
    ax.text(
        0.02,
        0.82,
        f"{args.dataset}: no strict theta/RG pair | {args.transform}",
        fontsize=14,
        fontweight="bold",
        va="top",
    )
    ax.text(
        0.02,
        0.58,
        "P5 labels were generated, but this dataset/transform has no method-k with both\n"
        "theta_selective and ripple-gamma/gamma/ripple selective components under the current strict rule.",
        fontsize=11,
        va="top",
    )
    fig_file = out_dir / f"{args.dataset}_{args.transform}_no_strict_theta_rg_pair_{args.ribbon}_diagnostic_plate.png"
    fig.tight_layout()
    fig.savefig(fig_file, dpi=160)
    plt.close(fig)
    return fig_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot one P5 event diagnostic plate.")
    parser.add_argument("--dataset", default="e10gb1")
    parser.add_argument("--condition", default=DEFAULT_CONDITION)
    parser.add_argument("--workspace", default="/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
    parser.add_argument("--processed-root", default="/mnt/e/DataPons_processed")
    parser.add_argument("--transform", default="adaptive_envelope", choices=["abs", "adaptive_envelope"])
    parser.add_argument("--method", default="auto")
    parser.add_argument("--k", default="auto")
    parser.add_argument(
        "--all-method-k",
        action="store_true",
        help="Plot every method-k pair that has both theta-like and ripple-gamma-like candidates.",
    )
    parser.add_argument("--half-sec", type=float, default=1.0)
    parser.add_argument("--ribbon", default="sem", choices=["sem", "std", "iqr", "none"])
    parser.add_argument("--max-events", type=int, default=2000)
    parser.add_argument("--window-length-samples", type=int, default=6000)
    parser.add_argument("--n-source-windows", type=int, default=15)
    return parser.parse_args()


def load_p2_events(dataset_root: Path, dataset: str) -> tuple[dict[str, BandEvents], float, np.ndarray, np.ndarray, int]:
    p2_file = dataset_root / "pipeline2_event_detection" / f"{dataset}_bandpass_events_3bands.mat"
    with h5py.File(p2_file, "r") as f:
        detect = f["R/DetectResults"]
        n_time = int(np.asarray(f["R/session_end_idx"]).squeeze()[-1])
        dt = float(np.asarray(f["R/dx"]).squeeze())
        session_start = np.asarray(f["R/session_start_idx"]).squeeze().astype(int) - 1
        session_end = np.asarray(f["R/session_end_idx"]).squeeze().astype(int)
        bands: dict[str, BandEvents] = {}
        min_gap_samples = max(1, int(round(0.05 / dt)))
        for band_col in range(detect.shape[1]):
            mask = np.zeros(n_time, dtype=bool)
            centers_all: list[np.ndarray] = []
            label = ""
            for ch in range(detect.shape[0]):
                group = f[detect[ch, band_col]]
                label = decode_matlab_char(group["band_label"])
                wins = np.asarray(group["event_win"], dtype=float).T
                peaks = np.asarray(group["loc_peak"], dtype=float).squeeze()
                if wins.size:
                    starts = np.maximum(0, wins[:, 0].astype(int) - 1)
                    ends = np.minimum(n_time, wins[:, 1].astype(int))
                    for start, end in zip(starts, ends):
                        if end > start:
                            mask[start:end] = True
                if peaks.size:
                    centers_all.append(peaks.astype(int) - 1)
            centers = np.concatenate(centers_all) if centers_all else np.zeros(0, dtype=int)
            bands[label] = BandEvents(label=label, mask=mask, centers=unique_centers(centers, min_gap_samples))
    return bands, dt, session_start, session_end, n_time


def decode_matlab_char(dataset) -> str:
    arr = np.asarray(dataset).squeeze()
    return "".join(chr(int(x)) for x in arr if int(x) > 0)


def unique_centers(centers: np.ndarray, min_gap_samples: int) -> np.ndarray:
    if centers.size == 0:
        return centers.astype(int)
    centers = np.sort(centers.astype(int))
    kept = [int(centers[0])]
    for c in centers[1:]:
        if int(c) - kept[-1] >= min_gap_samples:
            kept.append(int(c))
    return np.asarray(kept, dtype=int)


def choose_pairs(labels: pd.DataFrame, transform: str, method: str, k: str, all_method_k: bool) -> list[PairSpec]:
    df = labels[labels["activity_transform"] == transform].copy()
    df["rg_effect_z"] = df[["gamma_effect_z", "ripple_effect_z"]].max(axis=1)
    if method != "auto":
        df = df[df["method"] == method]
    if k != "auto":
        df = df[df["k"] == int(str(k).replace("k", ""))]

    pairs: list[tuple[float, pd.Series, pd.Series]] = []
    for (_, _), sub in df.groupby(["method", "k"]):
        theta = sub[sub["strict_label"] == "theta_selective"].sort_values("theta_effect_z", ascending=False)
        rg = sub[
            sub["strict_label"].isin(["ripple_gamma_no_theta", "gamma_selective", "ripple_selective"])
        ].sort_values("rg_effect_z", ascending=False)
        if theta.empty or rg.empty:
            continue
        th = theta.iloc[0]
        rr = rg.iloc[0]
        pairs.append((float(th["theta_effect_z"]) + float(rr["rg_effect_z"]), th, rr))

    if not pairs:
        return []
    if not all_method_k and method == "auto" and k == "auto":
        pairs = [max(pairs, key=lambda item: item[0])]

    out: list[PairSpec] = []
    for _, th, rr in sorted(pairs, key=lambda item: (str(item[1]["method"]), int(item[1]["k"]))):
        out.append(
            PairSpec(
                activity_transform=transform,
                method=str(th["method"]),
                k=int(th["k"]),
                theta_component=int(th["component"]),
                rg_component=int(rr["component"]),
                theta_effect_z=float(th["theta_effect_z"]),
                gamma_effect_z_for_theta=float(th["gamma_effect_z"]),
                ripple_effect_z_for_theta=float(th["ripple_effect_z"]),
                rg_theta_effect_z=float(rr["theta_effect_z"]),
                rg_gamma_effect_z=float(rr["gamma_effect_z"]),
                rg_ripple_effect_z=float(rr["ripple_effect_z"]),
                theta_source_file=str(th["source_file"]),
                rg_source_file=str(rr["source_file"]),
                theta_window_samples=int(th["envelope_window_samples"]),
                rg_window_samples=int(rr["envelope_window_samples"]),
            )
        )
    return out


def load_activity(
    source_file: str,
    component: int,
    transform: str,
    window_samples: int,
    session_start: np.ndarray,
    session_end: np.ndarray,
) -> np.ndarray:
    with h5py.File(source_file, "r") as f:
        x = np.abs(
            np.asarray(
                f["result/core/temporal_components_time_by_comp"][component - 1, :],
                dtype=np.float32,
            )
        ).astype(np.float64)
    if transform == "adaptive_envelope":
        x = rms_envelope_sessionwise(x, session_start, session_end, window_samples)
    return x


def rms_envelope_sessionwise(x: np.ndarray, session_start: np.ndarray, session_end: np.ndarray, window_samples: int) -> np.ndarray:
    window_samples = max(1, int(window_samples))
    out = np.empty_like(x, dtype=np.float64)
    for start, end in zip(session_start, session_end):
        start = int(start)
        end = int(end)
        if end <= start:
            continue
        out[start:end] = moving_rms(x[start:end], window_samples)
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


def plot_plate(
    args: argparse.Namespace,
    pair: PairSpec,
    theta_activity: np.ndarray,
    rg_activity: np.ndarray,
    bands: dict[str, BandEvents],
    dt: float,
    session_start: np.ndarray,
    session_end: np.ndarray,
    n_time: int,
    out_dir: Path,
) -> Path:
    theta_windows = rank_windows(theta_activity, n_time, args.window_length_samples, args.n_source_windows, session_start, session_end)
    rg_windows = rank_windows(rg_activity, n_time, args.window_length_samples, args.n_source_windows, session_start, session_end)
    windows = [("theta-top", w) for w in theta_windows] + [("rg-top", w) for w in rg_windows]

    out_file = (
        out_dir
        / f"{args.dataset}_{pair.activity_transform}_{pair.method}_k{pair.k:02d}_"
        f"thetaC{pair.theta_component:02d}_rgC{pair.rg_component:02d}_"
        f"{args.ribbon}_diagnostic_plate.png"
    )

    fig = plt.figure(figsize=(30, 15.5))
    gs = fig.add_gridspec(
        6,
        7,
        width_ratios=[1.45, 1.45, 1, 1, 1, 1, 1],
        wspace=0.32,
        hspace=0.48,
    )
    ax_theta = fig.add_subplot(gs[0:3, 0:2])
    ax_rg = fig.add_subplot(gs[3:6, 0:2], sharex=ax_theta)

    x = np.arange(-int(round(args.half_sec / dt)), int(round(args.half_sec / dt)) + 1) * dt
    plot_perievent_axis(
        ax_theta,
        x,
        theta_activity,
        bands,
        dt,
        args.half_sec,
        args.ribbon,
        args.max_events,
        f"theta-like C{pair.theta_component} | T/G/R z="
        f"{pair.theta_effect_z:.2f}/{pair.gamma_effect_z_for_theta:.2f}/{pair.ripple_effect_z_for_theta:.2f}",
    )
    plot_perievent_axis(
        ax_rg,
        x,
        rg_activity,
        bands,
        dt,
        args.half_sec,
        args.ribbon,
        args.max_events,
        f"ripple-gamma-like C{pair.rg_component} | T/G/R z="
        f"{pair.rg_theta_effect_z:.2f}/{pair.rg_gamma_effect_z:.2f}/{pair.rg_ripple_effect_z:.2f}",
    )
    ax_rg.set_xlabel("time from P2 event peak (s)")
    ax_theta.legend(frameon=False, ncol=3, loc="upper right", fontsize=8)

    for idx in range(30):
        ax = fig.add_subplot(gs[idx // 5, 2 + idx % 5])
        if idx >= len(windows):
            ax.axis("off")
            continue
        source, (global_window_idx, start0, end0, score, session_id) = windows[idx]
        t = np.arange(end0 - start0) * dt
        y_theta = robust_scale(theta_activity[start0:end0])
        y_rg = robust_scale(rg_activity[start0:end0])
        ax.plot(t, 3.35 + y_theta, color=TRACE_COLORS["theta"], linewidth=0.75)
        ax.plot(t, 2.25 + y_rg, color=TRACE_COLORS["ripple_gamma"], linewidth=0.75)
        draw_event_raster(ax, start0, end0, bands, dt)
        ax.set_ylim(y_limits(y_theta, y_rg))
        ax.set_title(
            f"{source} #{(idx % args.n_source_windows) + 1} | gw{global_window_idx} s{session_id} | {score:.3g}",
            fontsize=7,
        )

    fig.suptitle(
        f"{args.dataset}: P5 event diagnostic plate | {pair.activity_transform} | {pair.method}_k{pair.k:02d}\n"
        f"left: all P2-event peri-event mean +/- {args.ribbon}; right: top activity windows with P2 event rugs",
        fontsize=16,
        weight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_file, dpi=170)
    plt.close(fig)
    return out_file


def plot_perievent_axis(
    ax,
    x: np.ndarray,
    activity: np.ndarray,
    bands: dict[str, BandEvents],
    dt: float,
    half_sec: float,
    ribbon: str,
    max_events: int,
    title: str,
) -> None:
    half = int(round(half_sec / dt))
    for band_name in BAND_ORDER:
        mat = perievent_matrix(activity, bands[band_name].centers, half, max_events)
        if mat.size == 0:
            continue
        mean = np.nanmean(mat, axis=0)
        lo, hi = ribbon_bounds(mat, mean, ribbon)
        ax.plot(x, mean, color=BAND_COLORS[band_name], linewidth=1.4, label=band_name)
        if lo is not None and hi is not None:
            ax.fill_between(x, lo, hi, color=BAND_COLORS[band_name], alpha=0.18, linewidth=0)
    ax.axvline(0, color="0.35", linestyle="--", linewidth=0.9)
    ax.set_title(title, fontsize=9.5)
    ax.set_ylabel("activity")
    ax.grid(True, color="0.9", linewidth=0.7)


def perievent_matrix(activity: np.ndarray, centers: np.ndarray, half: int, max_events: int) -> np.ndarray:
    valid = centers[(centers >= half) & (centers < activity.size - half)]
    if valid.size == 0:
        return np.empty((0, 2 * half + 1), dtype=np.float64)
    if valid.size > max_events:
        idx = np.linspace(0, valid.size - 1, max_events).astype(int)
        valid = valid[idx]
    return np.vstack([activity[c - half : c + half + 1] for c in valid])


def ribbon_bounds(mat: np.ndarray, mean: np.ndarray, ribbon: str) -> tuple[np.ndarray | None, np.ndarray | None]:
    if ribbon == "none":
        return None, None
    if ribbon == "std":
        spread = np.nanstd(mat, axis=0)
        return mean - spread, mean + spread
    if ribbon == "sem":
        spread = np.nanstd(mat, axis=0) / math.sqrt(max(1, mat.shape[0]))
        return mean - spread, mean + spread
    if ribbon == "iqr":
        return np.nanpercentile(mat, 25, axis=0), np.nanpercentile(mat, 75, axis=0)
    raise ValueError(f"Unknown ribbon: {ribbon}")


def rank_windows(
    activity: np.ndarray,
    n_time: int,
    window_len: int,
    n_top: int,
    session_start: np.ndarray,
    session_end: np.ndarray,
) -> list[tuple[int, int, int, float, int]]:
    n = min(int(n_time), int(activity.size))
    n_windows = n // window_len
    scores = activity[: n_windows * window_len].reshape(n_windows, window_len).mean(axis=1)
    ranked: list[tuple[int, int, int, float, int]] = []
    for idx in np.argsort(scores)[::-1]:
        start0 = int(idx) * window_len
        end0 = start0 + window_len
        session_id = resolve_session(start0, end0, session_start, session_end)
        if session_id < 1:
            continue
        ranked.append((int(idx) + 1, start0, end0, float(scores[idx]), session_id))
        if len(ranked) >= n_top:
            break
    return ranked


def resolve_session(start0: int, end0: int, session_start: np.ndarray, session_end: np.ndarray) -> int:
    for i, (s0, e0) in enumerate(zip(session_start, session_end), start=1):
        if int(s0) <= start0 and end0 <= int(e0):
            return i
    return -1


def robust_scale(x: np.ndarray) -> np.ndarray:
    lo, hi = np.nanpercentile(x, [2, 98])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(x, dtype=np.float64)
    return (x - lo) / (hi - lo)


def y_limits(y_theta: np.ndarray, y_rg: np.ndarray) -> tuple[float, float]:
    vals = np.concatenate([3.35 + y_theta[np.isfinite(y_theta)], 2.25 + y_rg[np.isfinite(y_rg)]])
    if vals.size == 0:
        return (-0.4, 4.7)
    ymin = min(-0.4, float(vals.min()))
    ymax = max(4.5, float(vals.max()))
    pad = 0.08 * max(1e-9, ymax - ymin)
    return ymin - pad, ymax + pad


def draw_event_raster(ax, start0: int, end0: int, bands: dict[str, BandEvents], dt: float) -> None:
    duration = (end0 - start0) * dt
    for y, band_name in enumerate(BAND_ORDER):
        mask = bands[band_name].mask[start0:end0]
        for seg_start, seg_end in mask_segments(mask):
            ax.axvspan(
                seg_start * dt,
                seg_end * dt,
                ymin=(y + 0.08) / 5.0,
                ymax=(y + 0.78) / 5.0,
                color=BAND_COLORS[band_name],
                alpha=0.75,
                lw=0,
            )
        ax.text(-0.03 * duration, y + 0.42, band_name, ha="right", va="center", fontsize=6.5, color=BAND_COLORS[band_name])
    ax.set_xlim(0, duration)
    ax.set_yticks([])
    ax.grid(axis="x", color="0.9", linewidth=0.5)


def mask_segments(mask: np.ndarray) -> list[tuple[int, int]]:
    if mask.size == 0 or not mask.any():
        return []
    edges = np.diff(np.r_[0, mask.astype(np.int8), 0])
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)
    return list(zip(starts, ends))


if __name__ == "__main__":
    main()
