#!/usr/bin/env python3
"""Band-level P5 component-label preview for E10gb1 standardized csplit.

This is intentionally separate from the earlier event-family preview.  The
primary label here comes from continuous theta/gamma/ripple band-power traces,
not from overlapping event families.
"""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


DATASET = "e10gb1"
CONDITION = "complex_split_projected_vlambda_standardize"

WORKSPACE = Path("/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
DATA_ROOT = Path("/mnt/e/DataPons_processed") / DATASET
P1_ABS_MAT = (
    DATA_ROOT
    / "pipeline1_spectrograms"
    / "e10gb1_mirrorpad_20s_regionmean_spectrograms_abs.mat"
)
P5_ROOT = DATA_ROOT / "pipeline5_eigenfunction_reduction" / CONDITION

OUT_ROOT = WORKSPACE / "results" / "e10gb1_p5_bandlevel_preview_20260528"
FIG_DIR = OUT_ROOT / "figures"
CACHE_DIR = OUT_ROOT / "cache"

BANDS = [
    ("theta", 2.0, 15.0),
    ("gamma", 30.0, 90.0),
    ("ripple", 90.0, 190.0),
]

METHODS = ["svd", "nmf", "mds", "umap"]
KS = list(range(3, 9))
BIN_SEC = 1.0
BORDER_SEC = 20.0
MIN_CORR = 0.08
RELATIVE_ACTIVE_FRAC = 0.60
RGBA = {
    "theta": "#2ca25f",
    "gamma": "#f0b429",
    "ripple": "#3182bd",
    "ripple_gamma": "#756bb1",
    "theta_gamma": "#74c476",
    "theta_ripple": "#9ecae1",
    "broad": "#a6611a",
    "weak": "#eeeeee",
}
LABEL_ORDER = [
    "theta",
    "gamma",
    "ripple",
    "ripple_gamma",
    "theta_gamma",
    "theta_ripple",
    "broad",
    "weak",
]


@dataclass
class ComponentRecord:
    method: str
    k: int
    component: int
    label: str
    theta_corr: float
    gamma_corr: float
    ripple_corr: float
    max_corr: float
    theta_score: float
    ripple_gamma_score: float
    source_file: str


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    band_cache = CACHE_DIR / "e10gb1_band_power_1s_bins.npz"
    if band_cache.exists():
        band_data = dict(np.load(band_cache, allow_pickle=True))
    else:
        band_data = build_band_power_cache(band_cache)

    comp_cache = CACHE_DIR / "e10gb1_p5_component_activity_1s_bins.npz"
    records_csv = OUT_ROOT / "component_band_scores.csv"
    if comp_cache.exists() and records_csv.exists():
        comp_data = dict(np.load(comp_cache, allow_pickle=True))
        records = read_records(records_csv)
    else:
        comp_data, records = build_component_cache_and_scores(band_data, comp_cache, records_csv)

    write_design_note()
    plot_band_score_heatmap(records)
    plot_label_grid(records)
    plot_label_composition(records)
    plot_theta_vs_rg_scatter(records)
    plot_two_subprocess_map(records)
    plot_pairing_diagram(records)
    plot_top_window_examples(records, band_data, comp_data)
    plot_band_trace_overview(band_data)

    print(f"Wrote {len(records)} component records to {OUT_ROOT}")
    print(f"Figures: {FIG_DIR}")


def build_valid_bins(f: h5py.File) -> tuple[np.ndarray, int, np.ndarray]:
    total_len = int(np.asarray(f["timesout"]).shape[0])
    session_start = np.asarray(f["session_start_idx"]).squeeze().astype(int) - 1
    session_end = np.asarray(f["session_end_idx"]).squeeze().astype(int)
    session_dx = np.asarray(f["session_dx"]).squeeze().astype(float)
    pad_samples = np.asarray(f["session_pad_samples"]).squeeze().astype(int)

    dx = float(np.nanmedian(session_dx))
    bin_samples = max(1, int(round(BIN_SEC / dx)))
    border_samples = max(0, int(round(BORDER_SEC / dx)))
    # Respect the saved mirror-padding metadata, but do not exceed the requested border.
    session_border = np.minimum(pad_samples, border_samples)

    bin_id = np.full(total_len, -1, dtype=np.int32)
    bin_session = []
    current = 0
    for sess_idx, (start, end, border) in enumerate(zip(session_start, session_end, session_border), start=1):
        valid_start = int(start + border)
        valid_end = int(end - border)
        if valid_end <= valid_start:
            continue
        n = valid_end - valid_start
        local_bins = np.arange(n, dtype=np.int64) // bin_samples
        local_bins = local_bins.astype(np.int32) + current
        bin_id[valid_start:valid_end] = local_bins
        bin_session.extend([sess_idx] * (int(local_bins[-1]) - current + 1))
        current = int(local_bins[-1]) + 1

    return bin_id, current, np.asarray(bin_session, dtype=np.int16)


def build_band_power_cache(cache_path: Path) -> dict[str, np.ndarray]:
    print("Building band-power cache from P1 spectrograms...")
    with h5py.File(P1_ABS_MAT, "r") as f:
        freqs = np.asarray(f["freqs"]).squeeze()
        dset = f["tmpall_mean_abs"]
        bin_id, n_bins, bin_session = build_valid_bins(f)
        band_indices = [
            np.where((freqs >= low) & (freqs <= high))[0] for _, low, high in BANDS
        ]
        if any(idx.size == 0 for idx in band_indices):
            raise RuntimeError("At least one band has no matching frequency bins.")

        sums = np.zeros((len(BANDS), n_bins), dtype=np.float64)
        counts = np.zeros(n_bins, dtype=np.float64)
        block = 32768
        total_len = dset.shape[1]
        for start in range(0, total_len, block):
            end = min(start + block, total_len)
            ids = bin_id[start:end]
            valid = ids >= 0
            if not np.any(valid):
                continue
            # Chunking is (region, time, frequency), so reading all frequencies once
            # per time block is faster than several narrow frequency slices.
            arr = np.asarray(dset[:, start:end, :], dtype=np.float32)
            ids_valid = ids[valid]
            counts += np.bincount(ids_valid, minlength=n_bins)
            for band_i, idx in enumerate(band_indices):
                values = np.nanmean(arr[:, valid, :][:, :, idx], axis=(0, 2))
                sums[band_i] += np.bincount(ids_valid, weights=values, minlength=n_bins)
            if start % (block * 20) == 0:
                print(f"  P1 block {start:,}/{total_len:,}")

        with np.errstate(invalid="ignore", divide="ignore"):
            band_power = sums / counts[None, :]
        band_power = np.log1p(np.maximum(band_power, 0.0))
        band_names = np.asarray([b[0] for b in BANDS], dtype=object)
        band_edges = np.asarray([(b[1], b[2]) for b in BANDS], dtype=float)

    np.savez_compressed(
        cache_path,
        band_power=band_power,
        band_names=band_names,
        band_edges=band_edges,
        bin_counts=counts,
        bin_session=bin_session,
        bin_sec=np.asarray([BIN_SEC], dtype=float),
        border_sec=np.asarray([BORDER_SEC], dtype=float),
    )
    return dict(np.load(cache_path, allow_pickle=True))


def build_component_cache_and_scores(
    band_data: dict[str, np.ndarray],
    comp_cache: Path,
    records_csv: Path,
) -> tuple[dict[str, np.ndarray], list[ComponentRecord]]:
    print("Building P5 component activity cache and band scores...")
    with h5py.File(P1_ABS_MAT, "r") as f:
        bin_id, n_bins, _ = build_valid_bins(f)

    band_power = np.asarray(band_data["band_power"], dtype=float)
    band_z = row_zscore(band_power)
    activity_rows = []
    keys = []
    records: list[ComponentRecord] = []

    for method in METHODS:
        for k in KS:
            mat_file = find_p5_mat(method, k)
            if mat_file is None:
                continue
            print(f"  {method}_k{k:02d}: {mat_file.name}")
            with h5py.File(mat_file, "r") as f:
                dset = f["result/core/temporal_components_time_by_comp"]
                n_comp = int(dset.shape[0])
                sums = np.zeros((n_comp, n_bins), dtype=np.float64)
                counts = np.zeros(n_bins, dtype=np.float64)
                block = 32768
                total_len = dset.shape[1]
                for start in range(0, total_len, block):
                    end = min(start + block, total_len)
                    ids = bin_id[start:end]
                    valid = ids >= 0
                    if not np.any(valid):
                        continue
                    x = np.asarray(dset[:, start:end], dtype=np.float64)
                    x = np.abs(x[:, valid])
                    ids_valid = ids[valid]
                    counts += np.bincount(ids_valid, minlength=n_bins)
                    for c in range(n_comp):
                        sums[c] += np.bincount(ids_valid, weights=x[c], minlength=n_bins)
                with np.errstate(invalid="ignore", divide="ignore"):
                    activity = sums / counts[None, :]

            activity_z = row_zscore(activity)
            corrs = np.matmul(activity_z, band_z.T) / max(1, activity_z.shape[1] - 1)
            corrs = np.nan_to_num(corrs, nan=0.0, posinf=0.0, neginf=0.0)

            for c in range(n_comp):
                theta_corr, gamma_corr, ripple_corr = [float(v) for v in corrs[c, :]]
                label = label_from_band_corr(theta_corr, gamma_corr, ripple_corr)
                rec = ComponentRecord(
                    method=method,
                    k=k,
                    component=c + 1,
                    label=label,
                    theta_corr=theta_corr,
                    gamma_corr=gamma_corr,
                    ripple_corr=ripple_corr,
                    max_corr=max(theta_corr, gamma_corr, ripple_corr),
                    theta_score=max(0.0, theta_corr),
                    ripple_gamma_score=max(0.0, gamma_corr, ripple_corr),
                    source_file=str(mat_file),
                )
                records.append(rec)
                activity_rows.append(activity[c].astype(np.float32))
                keys.append(f"{method}_k{k:02d}_c{c+1:02d}")

    activity_binned = np.vstack(activity_rows) if activity_rows else np.zeros((0, n_bins), dtype=np.float32)
    np.savez_compressed(
        comp_cache,
        activity_binned=activity_binned,
        component_keys=np.asarray(keys, dtype=object),
    )
    write_records(records_csv, records)
    return dict(np.load(comp_cache, allow_pickle=True)), records


def find_p5_mat(method: str, k: int) -> Path | None:
    mat_dir = P5_ROOT / f"{method}_k{k:02d}" / "mat"
    if not mat_dir.exists():
        return None
    candidates = sorted(mat_dir.glob("*.mat"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not candidates:
        return None
    return candidates[0]


def row_zscore(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    mu = np.nanmean(x, axis=1, keepdims=True)
    sd = np.nanstd(x, axis=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    z = (x - mu) / sd
    return np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)


def label_from_band_corr(theta: float, gamma: float, ripple: float) -> str:
    vals = np.asarray([theta, gamma, ripple], dtype=float)
    positive = np.maximum(vals, 0.0)
    max_score = float(np.max(positive))
    if max_score < MIN_CORR:
        return "weak"
    active = positive >= max(MIN_CORR, RELATIVE_ACTIVE_FRAC * max_score)
    t, g, r = [bool(v) for v in active]
    if t and not g and not r:
        return "theta"
    if g and not t and not r:
        return "gamma"
    if r and not t and not g:
        return "ripple"
    if g and r and not t:
        return "ripple_gamma"
    if t and g and not r:
        return "theta_gamma"
    if t and r and not g:
        return "theta_ripple"
    if t and g and r:
        return "broad"
    return "weak"


def write_records(path: Path, records: list[ComponentRecord]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(ComponentRecord.__dataclass_fields__.keys()))
        writer.writeheader()
        for rec in records:
            writer.writerow(rec.__dict__)


def read_records(path: Path) -> list[ComponentRecord]:
    records: list[ComponentRecord] = []
    with path.open("r", newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            records.append(
                ComponentRecord(
                    method=row["method"],
                    k=int(row["k"]),
                    component=int(row["component"]),
                    label=row["label"],
                    theta_corr=float(row["theta_corr"]),
                    gamma_corr=float(row["gamma_corr"]),
                    ripple_corr=float(row["ripple_corr"]),
                    max_corr=float(row["max_corr"]),
                    theta_score=float(row["theta_score"]),
                    ripple_gamma_score=float(row["ripple_gamma_score"]),
                    source_file=row["source_file"],
                )
            )
    return records


def record_key(rec: ComponentRecord) -> str:
    return f"{rec.method}_k{rec.k:02d}_c{rec.component:02d}"


def sorted_records(records: list[ComponentRecord]) -> list[ComponentRecord]:
    order = {m: i for i, m in enumerate(METHODS)}
    return sorted(records, key=lambda r: (order.get(r.method, 99), r.k, r.component))


def plot_band_score_heatmap(records: list[ComponentRecord]) -> None:
    rows = sorted_records(records)
    labels = [f"{r.method}_k{r.k:02d}_C{r.component}" for r in rows]
    mat = np.asarray([[r.theta_corr, r.gamma_corr, r.ripple_corr] for r in rows], dtype=float)
    fig_h = max(10, 0.13 * len(rows))
    fig, ax = plt.subplots(figsize=(8.5, fig_h))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-0.35, vmax=0.35)
    ax.set_xticks(range(3), ["theta", "gamma", "ripple"])
    ax.set_yticks(range(len(labels)), labels, fontsize=6)
    ax.set_title(f"{DATASET} P5 standardized csplit: component correlation with continuous band power")
    cbar = fig.colorbar(im, ax=ax, shrink=0.6)
    cbar.set_label("Pearson r, component abs activity vs log band power")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "01_component_band_correlation_heatmap.png", dpi=180)
    plt.close(fig)


def plot_label_grid(records: list[ComponentRecord]) -> None:
    by_mk = group_by_method_k(records)
    row_keys = [(m, k) for m in METHODS for k in KS if (m, k) in by_mk]
    max_k = max((len(by_mk[key]) for key in row_keys), default=8)
    fig, ax = plt.subplots(figsize=(1.15 * max_k + 3.0, 0.58 * len(row_keys) + 1.8))
    ax.set_xlim(-1.35, max_k)
    ax.set_ylim(-0.5, len(row_keys) - 0.5)
    ax.invert_yaxis()
    ax.axis("off")
    for y, key in enumerate(row_keys):
        method, k = key
        ax.text(-1.25, y, f"{method}_k{k:02d}", va="center", ha="left", fontsize=10, weight="bold")
        for x, rec in enumerate(sorted(by_mk[key], key=lambda r: r.component)):
            color = RGBA.get(rec.label, "#eeeeee")
            rect = plt.Rectangle((x, y - 0.32), 0.88, 0.64, facecolor=color, edgecolor="white", linewidth=0.8)
            ax.add_patch(rect)
            text = f"C{rec.component}\n{short_label(rec.label)}"
            ax.text(x + 0.44, y, text, ha="center", va="center", fontsize=7.5, color="black")
    add_label_legend(ax, x0=0.02, y0=-0.09)
    ax.set_title(f"{DATASET} P5 standardized csplit: band-level component label grid", fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "02_band_label_grid_by_method_k.png", dpi=180)
    plt.close(fig)


def plot_label_composition(records: list[ComponentRecord]) -> None:
    by_mk = group_by_method_k(records)
    row_keys = [(m, k) for m in METHODS for k in KS if (m, k) in by_mk]
    counts = np.zeros((len(row_keys), len(LABEL_ORDER)), dtype=float)
    for i, key in enumerate(row_keys):
        labels = [r.label for r in by_mk[key]]
        for j, label in enumerate(LABEL_ORDER):
            counts[i, j] = labels.count(label)
    frac = counts / np.maximum(1.0, counts.sum(axis=1, keepdims=True))
    fig, ax = plt.subplots(figsize=(12, max(5, 0.28 * len(row_keys) + 1.5)))
    left = np.zeros(len(row_keys))
    y = np.arange(len(row_keys))
    for j, label in enumerate(LABEL_ORDER):
        ax.barh(y, frac[:, j], left=left, color=RGBA[label], edgecolor="white", height=0.8, label=label)
        left += frac[:, j]
    ax.set_yticks(y, [f"{m}_k{k:02d}" for m, k in row_keys])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of components")
    ax.set_title(f"{DATASET} P5 standardized csplit: band-label composition by method-k")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "03_band_label_composition_by_method_k.png", dpi=180)
    plt.close(fig)


def plot_theta_vs_rg_scatter(records: list[ComponentRecord]) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    for label in LABEL_ORDER:
        xs = [r.theta_score for r in records if r.label == label]
        ys = [r.ripple_gamma_score for r in records if r.label == label]
        if xs:
            ax.scatter(xs, ys, s=42, alpha=0.78, color=RGBA[label], edgecolor="white", linewidth=0.5, label=label)
    ax.axvline(MIN_CORR, color="0.7", linestyle="--", linewidth=1)
    ax.axhline(MIN_CORR, color="0.7", linestyle="--", linewidth=1)
    ax.set_xlabel("theta score = max(0, corr(component activity, theta power))")
    ax.set_ylabel("ripple-gamma score = max(0, corr with gamma/ripple power)")
    ax.set_title(f"{DATASET} P5 standardized csplit: theta vs ripple-gamma band score")
    ax.grid(True, color="0.9", linewidth=0.8)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "04_theta_vs_ripple_gamma_band_score_scatter.png", dpi=180)
    plt.close(fig)


def plot_two_subprocess_map(records: list[ComponentRecord]) -> None:
    by_mk = group_by_method_k(records)
    matrix = np.full((len(METHODS), len(KS)), np.nan)
    annot = [["" for _ in KS] for _ in METHODS]
    for i, method in enumerate(METHODS):
        for j, k in enumerate(KS):
            rs = by_mk.get((method, k), [])
            if not rs:
                continue
            has_theta = any(r.label == "theta" for r in rs)
            has_rg = any(r.label in {"ripple_gamma", "gamma", "ripple"} for r in rs)
            value = (1 if has_theta else 0) + (2 if has_rg else 0)
            matrix[i, j] = value
            best_theta = max((r for r in rs if r.label == "theta"), key=lambda r: r.theta_score, default=None)
            best_rg = max(
                (r for r in rs if r.label in {"ripple_gamma", "gamma", "ripple"}),
                key=lambda r: r.ripple_gamma_score,
                default=None,
            )
            parts = []
            if best_theta:
                parts.append(f"T:C{best_theta.component}")
            if best_rg:
                parts.append(f"RG:C{best_rg.component}")
            annot[i][j] = "\n".join(parts) if parts else "-"

    cmap = plt.matplotlib.colors.ListedColormap(["#eeeeee", "#2ca25f", "#756bb1", "#225ea8"])
    norm = plt.matplotlib.colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    fig, ax = plt.subplots(figsize=(9, 4.4))
    im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(KS)), [f"k{k:02d}" for k in KS])
    ax.set_yticks(range(len(METHODS)), METHODS)
    ax.set_title(f"{DATASET} P5 standardized csplit: has theta-only and ripple/gamma no-theta candidates")
    for i in range(len(METHODS)):
        for j in range(len(KS)):
            ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=9, color="black")
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2, 3], shrink=0.8)
    cbar.ax.set_yticklabels(["none", "theta only", "RG only", "both"])
    fig.tight_layout()
    fig.savefig(FIG_DIR / "05_two_subprocess_candidate_map_bandlevel.png", dpi=180)
    plt.close(fig)


def plot_pairing_diagram(records: list[ComponentRecord]) -> None:
    by_mk = group_by_method_k(records)
    row_keys = [(m, k) for m in METHODS for k in KS if (m, k) in by_mk]
    fig, ax = plt.subplots(figsize=(15, max(7, 0.55 * len(row_keys) + 2)))
    ax.set_xlim(-0.5, 6.5)
    ax.set_ylim(-0.5, len(row_keys) - 0.5)
    ax.invert_yaxis()
    ax.axis("off")
    headers = ["method-k", "best theta", "best gamma", "best ripple", "best ripple-gamma/no-theta", "pair score"]
    xs = [0, 1.4, 2.6, 3.8, 5.0, 6.15]
    for x, h in zip(xs, headers):
        ax.text(x, -1.0, h, ha="center", va="center", fontsize=10, weight="bold")
    for y, key in enumerate(row_keys):
        rs = by_mk[key]
        method, k = key
        ax.text(xs[0], y, f"{method}_k{k:02d}", ha="center", va="center", fontsize=10, weight="bold")
        best_theta = max((r for r in rs if r.label == "theta"), key=lambda r: r.theta_score, default=None)
        best_gamma = max((r for r in rs if r.label == "gamma"), key=lambda r: r.gamma_corr, default=None)
        best_ripple = max((r for r in rs if r.label == "ripple"), key=lambda r: r.ripple_corr, default=None)
        best_rg = max(
            (r for r in rs if r.label in {"ripple_gamma", "gamma", "ripple"}),
            key=lambda r: r.ripple_gamma_score,
            default=None,
        )
        draw_candidate_box(ax, xs[1], y, best_theta, "theta", lambda r: r.theta_score)
        draw_candidate_box(ax, xs[2], y, best_gamma, "gamma", lambda r: max(0.0, r.gamma_corr))
        draw_candidate_box(ax, xs[3], y, best_ripple, "ripple", lambda r: max(0.0, r.ripple_corr))
        draw_candidate_box(ax, xs[4], y, best_rg, best_rg.label if best_rg else "weak", lambda r: r.ripple_gamma_score)
        if best_theta and best_rg:
            score = best_theta.theta_score + best_rg.ripple_gamma_score
            text = f"{score:.2f}"
        else:
            text = "-"
        ax.text(xs[5], y, text, ha="center", va="center", fontsize=9)
    ax.set_title(f"{DATASET} P5 standardized csplit: best band-level candidate pairing per method-k", fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "06_best_band_candidate_pairing_diagram.png", dpi=180)
    plt.close(fig)


def draw_candidate_box(ax, x: float, y: float, rec: ComponentRecord | None, label: str, score_fn) -> None:
    color = RGBA.get(label, "#eeeeee")
    rect = plt.Rectangle((x - 0.42, y - 0.28), 0.84, 0.56, facecolor=color, edgecolor="white", linewidth=0.8)
    ax.add_patch(rect)
    if rec is None:
        text = "-"
    else:
        text = f"C{rec.component}\n{score_fn(rec):.2f}"
    ax.text(x, y, text, ha="center", va="center", fontsize=8)


def plot_top_window_examples(
    records: list[ComponentRecord],
    band_data: dict[str, np.ndarray],
    comp_data: dict[str, np.ndarray],
) -> None:
    activity = np.asarray(comp_data["activity_binned"], dtype=float)
    keys = [str(k) for k in comp_data["component_keys"]]
    key_to_idx = {k: i for i, k in enumerate(keys)}
    band_power = np.asarray(band_data["band_power"], dtype=float)
    band_z = row_zscore(band_power)

    pairs = choose_best_pairs(records, limit=4)
    if not pairs:
        return

    fig, axes = plt.subplots(len(pairs), 2, figsize=(13, 2.8 * len(pairs)), sharex=False)
    if len(pairs) == 1:
        axes = np.asarray([axes])
    for row, (theta_rec, rg_rec) in enumerate(pairs):
        for col, rec in enumerate([theta_rec, rg_rec]):
            ax = axes[row, col]
            key = record_key(rec)
            idx = key_to_idx[key]
            comp = zscore_1d(activity[idx])
            center = int(np.nanargmax(comp))
            half = 60
            lo = max(0, center - half)
            hi = min(comp.size, center + half + 1)
            x = np.arange(lo, hi) * BIN_SEC
            ax.plot(x, comp[lo:hi], color="black", linewidth=1.4, label=f"{key} activity")
            for b_i, (band_name, _, _) in enumerate(BANDS):
                ax.plot(x, band_z[b_i, lo:hi], linewidth=1.0, alpha=0.85, label=band_name)
            ax.axvline(center * BIN_SEC, color="0.5", linestyle="--", linewidth=0.9)
            ax.set_title(
                f"{rec.method}_k{rec.k:02d} C{rec.component} | {rec.label} | "
                f"T/G/R={rec.theta_corr:.2f}/{rec.gamma_corr:.2f}/{rec.ripple_corr:.2f}",
                fontsize=9,
            )
            ax.set_xlabel("analysis time bin (s)")
            ax.set_ylabel("z")
            ax.grid(True, color="0.9", linewidth=0.8)
            if row == 0 and col == 1:
                ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=8)
    fig.suptitle(f"{DATASET} P5 standardized csplit: top binned windows for selected theta/RG candidates", fontsize=14, weight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "07_top_window_examples_for_band_candidates.png", dpi=180)
    plt.close(fig)


def plot_band_trace_overview(band_data: dict[str, np.ndarray]) -> None:
    band_power = np.asarray(band_data["band_power"], dtype=float)
    band_z = row_zscore(band_power)
    n = band_z.shape[1]
    # A compact downsampled overview helps catch obviously broken band traces.
    step = max(1, n // 3000)
    x = np.arange(0, n, step) * BIN_SEC
    fig, ax = plt.subplots(figsize=(13, 3.8))
    for b_i, (band_name, _, _) in enumerate(BANDS):
        ax.plot(x, band_z[b_i, ::step], linewidth=0.9, alpha=0.85, label=band_name)
    ax.set_title(f"{DATASET}: continuous band-power traces used for P5 band-level labels")
    ax.set_xlabel("analysis time bin (s), 1 s bins downsampled for display")
    ax.set_ylabel("z(log band power)")
    ax.grid(True, color="0.9", linewidth=0.8)
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "08_band_power_trace_overview.png", dpi=180)
    plt.close(fig)


def choose_best_pairs(records: list[ComponentRecord], limit: int = 4) -> list[tuple[ComponentRecord, ComponentRecord]]:
    pairs = []
    by_mk = group_by_method_k(records)
    for key, rs in by_mk.items():
        theta = max((r for r in rs if r.label == "theta"), key=lambda r: r.theta_score, default=None)
        rg = max(
            (r for r in rs if r.label in {"ripple_gamma", "gamma", "ripple"}),
            key=lambda r: r.ripple_gamma_score,
            default=None,
        )
        if theta is None or rg is None:
            continue
        pairs.append((theta.theta_score + rg.ripple_gamma_score, theta, rg))
    pairs.sort(reverse=True, key=lambda x: x[0])
    return [(t, rg) for _, t, rg in pairs[:limit]]


def zscore_1d(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    sd = float(np.nanstd(x))
    if not np.isfinite(sd) or sd == 0:
        sd = 1.0
    return np.nan_to_num((x - np.nanmean(x)) / sd)


def group_by_method_k(records: list[ComponentRecord]) -> dict[tuple[str, int], list[ComponentRecord]]:
    out: dict[tuple[str, int], list[ComponentRecord]] = {}
    for rec in records:
        out.setdefault((rec.method, rec.k), []).append(rec)
    return out


def short_label(label: str) -> str:
    return {
        "theta": "T",
        "gamma": "G",
        "ripple": "R",
        "ripple_gamma": "RG",
        "theta_gamma": "TG",
        "theta_ripple": "TR",
        "broad": "TGR",
        "weak": "-",
    }.get(label, label)


def add_label_legend(ax, x0: float, y0: float) -> None:
    handles = [
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=RGBA[label], markeredgecolor="white", markersize=10)
        for label in LABEL_ORDER
    ]
    ax.legend(
        handles,
        LABEL_ORDER,
        loc="lower left",
        bbox_to_anchor=(x0, y0),
        ncol=4,
        frameon=False,
        fontsize=8,
    )


def write_design_note() -> None:
    text = f"""# E10gb1 P5 Band-Level Component Label Preview

This preview intentionally uses continuous band power rather than overlapping
event-family labels as the primary component label.

Data:
- Dataset: `{DATASET}`
- P5 condition: `{CONDITION}`
- P1 band source: `{P1_ABS_MAT}`
- P5 reduction source root: `{P5_ROOT}`

Band definitions match the current P2 detector defaults:
- theta: 2-15 Hz
- gamma: 30-90 Hz
- ripple: 90-190 Hz

Method:
1. Exclude 20 s at each session edge using the saved P1 session metadata.
2. Bin P1 region-mean abs spectrogram band power into 1 s bins.
3. Bin each P5 dimred component absolute activity into the same 1 s bins.
4. Compute Pearson correlation between component activity and log band power.
5. Label each component by the positive band correlations.

The label is now a band profile:
- `theta`, `gamma`, `ripple`: one dominant band.
- `ripple_gamma`: gamma and ripple jointly active without theta.
- `theta_gamma`, `theta_ripple`, `broad`: mixed band profiles.
- `weak`: no band reaches the current positive-correlation threshold.

Event-family activity should be used later as validation, not as the first-pass
label definition.
"""
    (OUT_ROOT / "README_bandlevel_preview.md").write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
