#!/usr/bin/env python3
"""Raw BLP saturation QC stratified by P2 event masks.

The main use case is checking whether a dataset's P2 theta events are driven by
raw signal clipping/saturation.  For each included session and selected BLP
channel, this script compares simple raw-level saturation proxies inside P2
theta/gamma/ripple event masks against non-event baseline samples.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BAND_ORDER = ["theta", "gamma", "ripple"]
BAND_COLORS = {
    "theta": "#2ca25f",
    "gamma": "#f0b429",
    "ripple": "#3182bd",
    "baseline": "#777777",
}


@dataclass
class BandMask:
    label: str
    mask: np.ndarray


def main() -> None:
    args = parse_args()
    workspace = Path(args.workspace)
    processed_root = Path(args.processed_root)
    raw_root = Path(args.raw_root) if args.raw_root else Path("/mnt/e/DataPons") / args.raw_dataset
    out_root = workspace / "results" / f"{args.dataset}_raw_blp_saturation_qc_20260528"
    out_root.mkdir(parents=True, exist_ok=True)

    bands, baseline_mask, session_start, session_end = load_p2_masks(processed_root, args.dataset)
    rows = []
    for i, sid in enumerate(args.sessions):
        session_mask_slice = slice(int(session_start[i]), int(session_end[i]))
        raw_file = raw_root / "blp" / f"{args.dataset}_{sid:04d}_blp.mat"
        if not raw_file.exists():
            raise FileNotFoundError(f"Raw BLP file not found: {raw_file}")
        x = load_raw_blp_first_page(raw_file, args.channels)
        if x.shape[0] != session_end[i] - session_start[i]:
            raise ValueError(
                f"Length mismatch for session {sid}: raw={x.shape[0]}, "
                f"P2={session_end[i] - session_start[i]}"
            )

        masks = {name: bands[name].mask[session_mask_slice] for name in BAND_ORDER}
        masks["baseline"] = baseline_mask[session_mask_slice]
        for local_ch, channel_id in enumerate(args.channels):
            rows.extend(
                summarize_channel(
                    x[:, local_ch],
                    masks,
                    dataset=args.dataset,
                    session_id=sid,
                    channel_id=channel_id,
                    edge_tol_scale=args.edge_tol_scale,
                    flat_tol_scale=args.flat_tol_scale,
                    robust_z=args.robust_z,
                )
            )

    out_csv = out_root / f"{args.dataset}_raw_blp_saturation_qc_by_session_channel.csv"
    with out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    df = pd.DataFrame(rows)
    summary_csv = out_root / f"{args.dataset}_raw_blp_saturation_qc_summary.csv"
    make_summary(df).to_csv(summary_csv, index=False)
    plot_summary(df, out_root, args.dataset)

    print(f"Wrote row-level QC: {out_csv}")
    print(f"Wrote summary QC: {summary_csv}")
    print(f"Wrote figures: {out_root}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Raw BLP saturation QC around P2 events.")
    parser.add_argument("--dataset", default="k13m17", help="Processed dataset stem, e.g. k13m17.")
    parser.add_argument("--raw-dataset", default="K13.m17", help="Raw dataset folder under /mnt/e/DataPons.")
    parser.add_argument("--raw-root", default="", help="Explicit raw dataset root. Should contain blp/*.mat.")
    parser.add_argument("--processed-root", default="/mnt/e/DataPons_processed")
    parser.add_argument("--workspace", default="/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
    parser.add_argument("--sessions", type=int, nargs="+", default=list(range(2, 19)))
    parser.add_argument("--channels", type=int, nargs="+", default=list(range(1, 9)))
    parser.add_argument("--edge-tol-scale", type=float, default=1e-6)
    parser.add_argument("--flat-tol-scale", type=float, default=1e-7)
    parser.add_argument("--robust-z", type=float, default=20.0)
    return parser.parse_args()


def load_p2_masks(processed_root: Path, dataset: str) -> tuple[dict[str, BandMask], np.ndarray, np.ndarray, np.ndarray]:
    p2_file = processed_root / dataset / "pipeline2_event_detection" / f"{dataset}_bandpass_events_3bands.mat"
    if not p2_file.exists():
        raise FileNotFoundError(f"P2 event file not found: {p2_file}")

    with h5py.File(p2_file, "r") as f:
        detect = f["R/DetectResults"]
        n_time = int(np.asarray(f["R/session_end_idx"]).squeeze()[-1])
        session_start = np.asarray(f["R/session_start_idx"]).squeeze().astype(int) - 1
        session_end = np.asarray(f["R/session_end_idx"]).squeeze().astype(int)
        bands: dict[str, BandMask] = {}

        for band_col in range(detect.shape[1]):
            label = ""
            mask = np.zeros(n_time, dtype=bool)
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
            bands[label] = BandMask(label=label, mask=mask)

    any_event = np.zeros_like(next(iter(bands.values())).mask)
    for band in bands.values():
        any_event |= band.mask
    baseline_mask = ~any_event
    return bands, baseline_mask, session_start, session_end


def decode_matlab_char(dataset) -> str:
    arr = np.asarray(dataset).squeeze()
    return "".join(chr(int(x)) for x in arr if int(x) > 0)


def load_raw_blp_first_page(raw_file: Path, channels: list[int]) -> np.ndarray:
    zero_based = [ch - 1 for ch in channels]
    with h5py.File(raw_file, "r") as f:
        dset = f["blp/dat"]
        if dset.ndim != 3:
            raise ValueError(f"Expected blp/dat to be 3D in {raw_file}, got shape {dset.shape}")
        # MATLAB shape is time x channel x page.  HDF5 stores it reversed:
        # page x channel x time.  MATLAB loader uses blp.dat(:, channels, 1).
        x = np.asarray(dset[0, zero_based, :], dtype=np.float64).T
    return x


def summarize_channel(
    x: np.ndarray,
    masks: dict[str, np.ndarray],
    *,
    dataset: str,
    session_id: int,
    channel_id: int,
    edge_tol_scale: float,
    flat_tol_scale: float,
    robust_z: float,
) -> list[dict[str, object]]:
    finite = np.isfinite(x)
    xf = x[finite]
    data_min = float(np.nanmin(xf))
    data_max = float(np.nanmax(xf))
    q001, q01, q50, q99, q999 = np.nanpercentile(xf, [0.1, 1, 50, 99, 99.9])
    robust_range = max(float(q999 - q001), np.finfo(float).eps)
    edge_tol = max(np.finfo(float).eps, edge_tol_scale * robust_range)
    flat_tol = max(np.finfo(float).eps, flat_tol_scale * robust_range)

    edge = (x <= data_min + edge_tol) | (x >= data_max - edge_tol)
    high_edge = x >= data_max - edge_tol
    low_edge = x <= data_min + edge_tol

    flat = np.zeros_like(x, dtype=bool)
    dx = np.abs(np.diff(x))
    same = dx <= flat_tol
    flat[:-1] |= same
    flat[1:] |= same

    baseline = masks["baseline"] & finite
    med = float(np.nanmedian(x[baseline])) if baseline.any() else float(q50)
    mad = float(np.nanmedian(np.abs(x[baseline] - med))) if baseline.any() else float(np.nanmedian(np.abs(xf - q50)))
    sigma = max(1.4826 * mad, np.finfo(float).eps)
    extreme = np.abs((x - med) / sigma) >= robust_z

    rows = []
    for mask_name in BAND_ORDER + ["baseline"]:
        mask = masks[mask_name] & finite
        n = int(mask.sum())
        row = {
            "dataset": dataset,
            "session_id": int(session_id),
            "channel_id": int(channel_id),
            "mask": mask_name,
            "n_samples": n,
            "sample_fraction": float(n / x.size),
            "raw_min": data_min,
            "raw_max": data_max,
            "q001": float(q001),
            "q01": float(q01),
            "q50": float(q50),
            "q99": float(q99),
            "q999": float(q999),
            "edge_tol": float(edge_tol),
            "flat_tol": float(flat_tol),
            "edge_fraction": safe_frac(edge, mask),
            "high_edge_fraction": safe_frac(high_edge, mask),
            "low_edge_fraction": safe_frac(low_edge, mask),
            "flat_fraction": safe_frac(flat, mask),
            "robust_extreme_fraction": safe_frac(extreme, mask),
            "mean": safe_mean(x, mask),
            "std": safe_std(x, mask),
        }
        rows.append(row)
    return rows


def safe_frac(flag: np.ndarray, mask: np.ndarray) -> float:
    n = int(mask.sum())
    if n == 0:
        return math.nan
    return float(np.mean(flag[mask]))


def safe_mean(x: np.ndarray, mask: np.ndarray) -> float:
    if not mask.any():
        return math.nan
    return float(np.nanmean(x[mask]))


def safe_std(x: np.ndarray, mask: np.ndarray) -> float:
    if not mask.any():
        return math.nan
    return float(np.nanstd(x[mask]))


def make_summary(df: pd.DataFrame) -> pd.DataFrame:
    metrics = ["edge_fraction", "high_edge_fraction", "low_edge_fraction", "flat_fraction", "robust_extreme_fraction"]
    rows = []
    for mask, sub in df.groupby("mask", sort=False):
        row = {"mask": mask, "n_session_channel": int(len(sub))}
        for metric in metrics:
            vals = sub[metric].astype(float)
            row[f"{metric}_mean"] = float(vals.mean())
            row[f"{metric}_median"] = float(vals.median())
            row[f"{metric}_p95"] = float(vals.quantile(0.95))
            row[f"{metric}_max"] = float(vals.max())
        rows.append(row)
    return pd.DataFrame(rows)


def plot_summary(df: pd.DataFrame, out_root: Path, dataset: str) -> None:
    metrics = [
        ("edge_fraction", "raw edge fraction"),
        ("flat_fraction", "flatline fraction"),
        ("robust_extreme_fraction", "robust extreme fraction"),
    ]
    summary = make_summary(df)
    x_labels = BAND_ORDER + ["baseline"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(15, 4.6), sharex=True)
    for ax, (metric, label) in zip(axes, metrics):
        vals = []
        p95 = []
        for mask in x_labels:
            row = summary[summary["mask"] == mask].iloc[0]
            vals.append(float(row[f"{metric}_median"]))
            p95.append(float(row[f"{metric}_p95"]))
        x = np.arange(len(x_labels))
        ax.bar(x, vals, color=[BAND_COLORS.get(m, "#777777") for m in x_labels], alpha=0.85, label="median")
        ax.scatter(x, p95, color="black", s=24, zorder=3, label="p95")
        ax.set_xticks(x, x_labels, rotation=20, ha="right")
        ax.set_title(label)
        ax.set_ylabel("fraction")
        ax.grid(axis="y", color="0.9")
    axes[0].legend(frameon=False)
    fig.suptitle(f"{dataset}: raw BLP saturation QC by P2 event mask", fontsize=14, weight="bold")
    fig.tight_layout()
    fig.savefig(out_root / f"{dataset}_raw_blp_saturation_qc_summary.png", dpi=180)
    plt.close(fig)

    pivot = (
        df[df["mask"].isin(["theta", "baseline"])]
        .pivot_table(index=["session_id", "channel_id"], columns="mask", values="edge_fraction")
        .reset_index()
    )
    if {"theta", "baseline"}.issubset(pivot.columns):
        pivot["theta_minus_baseline_edge_fraction"] = pivot["theta"] - pivot["baseline"]
        mat = pivot.pivot(index="channel_id", columns="session_id", values="theta_minus_baseline_edge_fraction")
        fig, ax = plt.subplots(figsize=(12, 4.8))
        vmax = max(1e-8, float(np.nanpercentile(np.abs(mat.to_numpy()), 98)))
        im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_xticks(range(mat.shape[1]), mat.columns.astype(str), rotation=45)
        ax.set_yticks(range(mat.shape[0]), mat.index.astype(str))
        ax.set_xlabel("session id")
        ax.set_ylabel("raw BLP channel")
        ax.set_title(f"{dataset}: theta minus baseline raw-edge fraction")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("theta edge fraction - baseline edge fraction")
        fig.tight_layout()
        fig.savefig(out_root / f"{dataset}_theta_minus_baseline_edge_fraction_heatmap.png", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
