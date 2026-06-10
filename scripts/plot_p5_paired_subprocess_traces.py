#!/usr/bin/env python3
"""Paired theta-like and ripple-gamma-like P5 dimred component traces.

For a given dataset/condition, choose two components from the same method-k:

- theta-like: highest theta event effect among theta_selective components
- ripple-gamma-like: highest max(gamma, ripple) effect among
  ripple_gamma_no_theta/gamma/ripple components

Then plot both component activity traces in the same windows, with P2
theta/gamma/ripple event masks overlaid.  The figure is meant for visual
validation of a two-subprocess interpretation, not for rank/overlap scoring.
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
    "theta_component": "#238b45",
    "ripple_gamma_component": "#6a51a3",
}
DEFAULT_CONDITION = "complex_split_projected_vlambda_standardize"
DEFAULT_LABEL_SUFFIX = "p5_p2_band_event_response_v2_selective_envelope_20260528"


@dataclass
class BandEvents:
    label: str
    mask: np.ndarray


@dataclass
class PairSpec:
    activity_transform: str
    method: str
    k: int
    theta_component: int
    rg_component: int
    theta_effect_z: float
    rg_effect_z: float
    theta_source_file: str
    rg_source_file: str
    theta_window_samples: int
    rg_window_samples: int


def main() -> None:
    args = parse_args()
    workspace = Path(args.workspace)
    processed_root = Path(args.processed_root)
    dataset_root = processed_root / args.dataset
    out_root = (
        dataset_root
        / "pipeline5_paired_subprocess_traces"
        / args.condition
    )
    out_root.mkdir(parents=True, exist_ok=True)

    bands, dt, session_start, session_end, n_time = load_p2_events(dataset_root, args.dataset)
    label_csv = workspace / "results" / f"{args.dataset}_{DEFAULT_LABEL_SUFFIX}" / "component_band_event_response_v2.csv"
    if not label_csv.exists():
        raise FileNotFoundError(f"P5/P2 component label CSV not found: {label_csv}")
    labels = pd.read_csv(label_csv)

    pair_rows = []
    for transform in args.transforms:
        pairs = choose_pairs(labels, transform, args.method, args.k, args.all_method_k)
        if not pairs:
            pair_rows.append(write_no_pair_figure(args, transform, out_root))
        for pair in pairs:
            pair_rows.append(
                plot_one_pair(
                    args,
                    pair,
                    bands,
                    dt,
                    session_start,
                    session_end,
                    n_time,
                    out_root,
                )
            )

    manifest = out_root / "paired_subprocess_trace_manifest.csv"
    pd.DataFrame(pair_rows).to_csv(manifest, index=False)
    print(f"Wrote paired subprocess manifest: {manifest}")
    for row in pair_rows:
        print(row["figure_file"])


def write_no_pair_figure(args: argparse.Namespace, transform: str, out_root: Path) -> dict[str, object]:
    out_dir = out_root / transform
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(14, 3.2))
    ax.axis("off")
    ax.text(
        0.02,
        0.82,
        f"{args.dataset}: no strict theta/RG pair | {transform}",
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
    fig_file = out_dir / f"{args.dataset}_{transform}_no_strict_theta_rg_pair_paired_top_windows.png"
    fig.tight_layout()
    fig.savefig(fig_file, dpi=160)
    plt.close(fig)
    return {
        "dataset": args.dataset,
        "condition": args.condition,
        "activity_transform": transform,
        "method": args.method,
        "k": args.k,
        "theta_component": "",
        "ripple_gamma_component": "",
        "theta_effect_z": "",
        "ripple_gamma_effect_z": "",
        "status": "no_strict_theta_rg_pair",
        "figure_file": str(fig_file),
    }


def plot_one_pair(
    args: argparse.Namespace,
    pair: PairSpec,
    bands: dict[str, BandEvents],
    dt: float,
    session_start: np.ndarray,
    session_end: np.ndarray,
    n_time: int,
    out_root: Path,
) -> dict[str, object]:
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
        theta_windows = rank_windows(theta_activity, n_time, args.window_length_samples, args.n_source_windows)
        rg_windows = rank_windows(rg_activity, n_time, args.window_length_samples, args.n_source_windows)
        windows = [("theta-top", w) for w in theta_windows] + [("rg-top", w) for w in rg_windows]
        fig_file = plot_pair_windows(
            args.dataset,
            pair,
            theta_activity,
            rg_activity,
            windows,
            bands,
            dt,
            out_root,
        )
        return {
            "dataset": args.dataset,
            "condition": args.condition,
            "activity_transform": pair.activity_transform,
            "method": pair.method,
            "k": pair.k,
            "theta_component": pair.theta_component,
            "ripple_gamma_component": pair.rg_component,
            "theta_effect_z": pair.theta_effect_z,
            "ripple_gamma_effect_z": pair.rg_effect_z,
            "theta_source_file": pair.theta_source_file,
            "ripple_gamma_source_file": pair.rg_source_file,
            "figure_file": str(fig_file),
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot paired theta/RG P5 dimred component traces.")
    parser.add_argument("--dataset", default="e10gb1")
    parser.add_argument("--condition", default=DEFAULT_CONDITION)
    parser.add_argument("--workspace", default="/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
    parser.add_argument("--processed-root", default="/mnt/e/DataPons_processed")
    parser.add_argument("--transforms", nargs="+", default=["abs", "adaptive_envelope"], choices=["abs", "adaptive_envelope"])
    parser.add_argument("--method", default="auto", help="Method to force, or auto.")
    parser.add_argument("--k", default="auto", help="Component count to force, or auto.")
    parser.add_argument(
        "--all-method-k",
        action="store_true",
        help="Plot every method-k pair that has both theta-like and ripple-gamma-like candidates.",
    )
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
        for band_col in range(detect.shape[1]):
            mask = np.zeros(n_time, dtype=bool)
            label = ""
            for ch in range(detect.shape[0]):
                group = f[detect[ch, band_col]]
                label = decode_matlab_char(group["band_label"])
                wins = np.asarray(group["event_win"], dtype=float).T
                if wins.size:
                    starts = np.maximum(0, wins[:, 0].astype(int) - 1)
                    ends = np.minimum(n_time, wins[:, 1].astype(int))
                    for start, end in zip(starts, ends):
                        if end > start:
                            mask[start:end] = True
            bands[label] = BandEvents(label=label, mask=mask)
    return bands, dt, session_start, session_end, n_time


def decode_matlab_char(dataset) -> str:
    arr = np.asarray(dataset).squeeze()
    return "".join(chr(int(x)) for x in arr if int(x) > 0)


def choose_pairs(labels: pd.DataFrame, transform: str, method: str, k: str, all_method_k: bool) -> list[PairSpec]:
    df = labels[labels["activity_transform"] == transform].copy()
    df["rg_effect_z"] = df[["gamma_effect_z", "ripple_effect_z"]].max(axis=1)
    if method != "auto":
        df = df[df["method"] == method]
    if k != "auto":
        k_int = int(str(k).replace("k", ""))
        df = df[df["k"] == k_int]

    pairs: list[PairSpec] = []
    for (m, kk), sub in df.groupby(["method", "k"]):
        theta = sub[sub["strict_label"] == "theta_selective"].sort_values("theta_effect_z", ascending=False)
        rg = sub[
            sub["strict_label"].isin(["ripple_gamma_no_theta", "gamma_selective", "ripple_selective"])
        ].sort_values("rg_effect_z", ascending=False)
        if theta.empty or rg.empty:
            continue
        th = theta.iloc[0]
        rr = rg.iloc[0]
        pairs.append(
            PairSpec(
                activity_transform=transform,
                method=str(m),
                k=int(kk),
                theta_component=int(th["component"]),
                rg_component=int(rr["component"]),
                theta_effect_z=float(th["theta_effect_z"]),
                rg_effect_z=float(rr["rg_effect_z"]),
                theta_source_file=str(th["source_file"]),
                rg_source_file=str(rr["source_file"]),
                theta_window_samples=int(th["envelope_window_samples"]),
                rg_window_samples=int(rr["envelope_window_samples"]),
            )
        )
    if not pairs:
        return []
    pairs = sorted(pairs, key=lambda p: (p.method, p.k))
    if all_method_k or method != "auto" or k != "auto":
        return pairs
    return [max(pairs, key=lambda p: p.theta_effect_z + p.rg_effect_z)]


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


def rank_windows(activity: np.ndarray, n_time: int, window_len: int, n_top: int) -> list[tuple[int, int, int, float]]:
    n = min(int(n_time), int(activity.size))
    n_windows = n // window_len
    scores = activity[: n_windows * window_len].reshape(n_windows, window_len).mean(axis=1)
    order = np.argsort(scores)[::-1][:n_top]
    return [(int(i) + 1, int(i) * window_len, int(i) * window_len + window_len, float(scores[i])) for i in order]


def plot_pair_windows(
    dataset: str,
    pair: PairSpec,
    theta_activity: np.ndarray,
    rg_activity: np.ndarray,
    windows: list[tuple[str, tuple[int, int, int, float]]],
    bands: dict[str, BandEvents],
    dt: float,
    out_root: Path,
) -> Path:
    out_dir = out_root / pair.activity_transform
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = (
        out_dir
        / f"{dataset}_{pair.activity_transform}_{pair.method}_k{pair.k:02d}_"
        f"thetaC{pair.theta_component:02d}_rgC{pair.rg_component:02d}_paired_top_windows.png"
    )
    n_cols = 10
    n_rows = int(math.ceil(len(windows) / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(30, 2.6 * n_rows), sharex=False, sharey=False)
    axes = np.atleast_1d(axes).ravel()
    for i, ax in enumerate(axes):
        if i >= len(windows):
            ax.axis("off")
            continue
        source, (global_window_idx, start0, end0, score) = windows[i]
        t = np.arange(end0 - start0) * dt
        y_theta = robust_scale(theta_activity[start0:end0])
        y_rg = robust_scale(rg_activity[start0:end0])
        ax.plot(t, 3.35 + y_theta, color=TRACE_COLORS["theta_component"], linewidth=0.85, label="theta-like")
        ax.plot(t, 2.25 + y_rg, color=TRACE_COLORS["ripple_gamma_component"], linewidth=0.85, label="ripple-gamma-like")
        draw_event_raster(ax, start0, end0, bands, dt)
        ymin, ymax = y_limits(y_theta, y_rg)
        ax.set_ylim(ymin, ymax)
        ax.set_title(f"{source} #{(i % 15) + 1} | gw{global_window_idx} | score={score:.3g}", fontsize=8)
        if i == 0:
            ax.legend(frameon=False, loc="upper right", fontsize=7)
    title = (
        f"{dataset}: paired P5 subprocess traces | {pair.activity_transform} | "
        f"{pair.method}_k{pair.k:02d}\n"
        f"theta-like C{pair.theta_component} (theta z={pair.theta_effect_z:.2f}) vs "
        f"ripple-gamma-like C{pair.rg_component} (RG z={pair.rg_effect_z:.2f})"
    )
    fig.suptitle(title, fontsize=15, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(out_file, dpi=175)
    plt.close(fig)
    return out_file


def robust_scale(x: np.ndarray) -> np.ndarray:
    lo, hi = np.nanpercentile(x, [2, 98])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return np.zeros_like(x)
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
        ax.text(-0.03 * duration, y + 0.42, band_name, ha="right", va="center", fontsize=7, color=BAND_COLORS[band_name])
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
