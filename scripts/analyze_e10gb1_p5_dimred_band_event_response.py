#!/usr/bin/env python3
"""Preview P5 dimred component responses to P2 theta/gamma/ripple events.

This analysis uses P2 event_detection windows directly.  It does not rebuild
band power from P1 spectrograms, and it does not use overlapping event-family
labels as the primary component label.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


DATASET = "e10gb1"
CONDITION = "complex_split_projected_vlambda_standardize"

WORKSPACE = Path("/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
DATA_ROOT = Path("/mnt/e/DataPons_processed") / DATASET
P2_EVENT_MAT = DATA_ROOT / "pipeline2_event_detection" / f"{DATASET}_bandpass_events_3bands.mat"
P5_ROOT = DATA_ROOT / "pipeline5_eigenfunction_reduction" / CONDITION

OUT_ROOT = WORKSPACE / "results" / "e10gb1_p5_p2_band_event_response_20260528"
FIG_DIR = OUT_ROOT / "figures"

METHODS = ["svd", "nmf", "mds", "umap"]
KS = list(range(3, 9))
BAND_ORDER = ["theta", "gamma", "ripple"]

# These thresholds are deliberately soft for an exploratory preview.  The
# numeric CSV is the source of truth; labels are just a quick visual summary.
MIN_EFFECT_Z = 0.05
RELATIVE_ACTIVE_FRAC = 0.60

LABEL_ORDER = [
    "theta",
    "gamma",
    "ripple",
    "ripple_gamma",
    "theta_gamma",
    "theta_ripple",
    "broad",
    "inactive",
]
COLORS = {
    "theta": "#2ca25f",
    "gamma": "#f0b429",
    "ripple": "#3182bd",
    "ripple_gamma": "#756bb1",
    "theta_gamma": "#74c476",
    "theta_ripple": "#9ecae1",
    "broad": "#a6611a",
    "inactive": "#eeeeee",
}


@dataclass
class BandEvents:
    label: str
    passband_hz: tuple[float, float]
    mask: np.ndarray
    centers: np.ndarray
    n_raw_channel_events: int
    coverage_fraction: float


@dataclass
class ComponentRecord:
    method: str
    k: int
    component: int
    label: str
    theta_effect_z: float
    gamma_effect_z: float
    ripple_effect_z: float
    theta_event_mean: float
    gamma_event_mean: float
    ripple_event_mean: float
    baseline_mean: float
    baseline_sd: float
    theta_ratio: float
    gamma_ratio: float
    ripple_ratio: float
    source_file: str


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    bands, baseline_mask, dt = load_p2_band_events()
    records = compute_component_responses(bands, baseline_mask)
    write_records(OUT_ROOT / "component_band_event_response.csv", records)
    write_readme(bands, baseline_mask, dt)

    plot_event_coverage(bands, baseline_mask)
    plot_effect_heatmap(records)
    plot_label_grid(records)
    plot_label_composition(records)
    plot_theta_vs_ripple_gamma(records)
    plot_two_subprocess_map(records)
    plot_perievent_examples(records, bands, dt)

    print(f"Wrote {len(records)} component records")
    print(f"Output: {OUT_ROOT}")


def load_p2_band_events() -> tuple[dict[str, BandEvents], np.ndarray, float]:
    with h5py.File(P2_EVENT_MAT, "r") as f:
        detect = f["R/DetectResults"]
        n_time = int(np.asarray(f["R/session_end_idx"]).squeeze()[-1])
        dt = float(np.asarray(f["R/dx"]).squeeze())
        bands: dict[str, BandEvents] = {}

        for band_col in range(detect.shape[1]):
            mask = np.zeros(n_time, dtype=bool)
            centers_all: list[np.ndarray] = []
            n_raw = 0
            label = ""
            passband = (np.nan, np.nan)
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
            centers = unique_centers(centers, min_gap_samples=max(1, int(round(0.05 / dt))))
            bands[label] = BandEvents(
                label=label,
                passband_hz=passband,
                mask=mask,
                centers=centers,
                n_raw_channel_events=n_raw,
                coverage_fraction=float(np.mean(mask)),
            )

    any_event = np.zeros_like(next(iter(bands.values())).mask)
    for band in bands.values():
        any_event |= band.mask
    baseline_mask = ~any_event
    return bands, baseline_mask, dt


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


def compute_component_responses(
    bands: dict[str, BandEvents],
    baseline_mask: np.ndarray,
) -> list[ComponentRecord]:
    records: list[ComponentRecord] = []
    band_masks = {name: bands[name].mask for name in BAND_ORDER}

    for method in METHODS:
        for k in KS:
            mat_file = find_p5_mat(method, k)
            if mat_file is None:
                print(f"missing {method}_k{k:02d}")
                continue
            print(f"processing {method}_k{k:02d}: {mat_file.name}")
            with h5py.File(mat_file, "r") as f:
                dset = f["result/core/temporal_components_time_by_comp"]
                n_comp, n_time = dset.shape
                sums = {name: np.zeros(n_comp, dtype=np.float64) for name in BAND_ORDER}
                counts = {name: 0 for name in BAND_ORDER}
                base_sum = np.zeros(n_comp, dtype=np.float64)
                base_sumsq = np.zeros(n_comp, dtype=np.float64)
                base_count = 0

                block = 65536
                for start in range(0, n_time, block):
                    end = min(start + block, n_time)
                    x = np.abs(np.asarray(dset[:, start:end], dtype=np.float64))

                    bm = baseline_mask[start:end]
                    if np.any(bm):
                        xb = x[:, bm]
                        base_sum += np.sum(xb, axis=1)
                        base_sumsq += np.sum(xb * xb, axis=1)
                        base_count += int(np.sum(bm))

                    for name, mask in band_masks.items():
                        m = mask[start:end]
                        if np.any(m):
                            sums[name] += np.sum(x[:, m], axis=1)
                            counts[name] += int(np.sum(m))

            baseline_mean = base_sum / max(1, base_count)
            baseline_var = base_sumsq / max(1, base_count) - baseline_mean * baseline_mean
            baseline_sd = np.sqrt(np.maximum(baseline_var, 1e-12))

            band_means = {
                name: sums[name] / max(1, counts[name])
                for name in BAND_ORDER
            }
            band_effects = {
                name: (band_means[name] - baseline_mean) / baseline_sd
                for name in BAND_ORDER
            }
            band_ratios = {
                name: band_means[name] / np.maximum(baseline_mean, 1e-12)
                for name in BAND_ORDER
            }

            for comp in range(n_comp):
                theta = float(band_effects["theta"][comp])
                gamma = float(band_effects["gamma"][comp])
                ripple = float(band_effects["ripple"][comp])
                records.append(
                    ComponentRecord(
                        method=method,
                        k=k,
                        component=comp + 1,
                        label=label_from_effects(theta, gamma, ripple),
                        theta_effect_z=theta,
                        gamma_effect_z=gamma,
                        ripple_effect_z=ripple,
                        theta_event_mean=float(band_means["theta"][comp]),
                        gamma_event_mean=float(band_means["gamma"][comp]),
                        ripple_event_mean=float(band_means["ripple"][comp]),
                        baseline_mean=float(baseline_mean[comp]),
                        baseline_sd=float(baseline_sd[comp]),
                        theta_ratio=float(band_ratios["theta"][comp]),
                        gamma_ratio=float(band_ratios["gamma"][comp]),
                        ripple_ratio=float(band_ratios["ripple"][comp]),
                        source_file=str(mat_file),
                    )
                )
    return records


def find_p5_mat(method: str, k: int) -> Path | None:
    mat_dir = P5_ROOT / f"{method}_k{k:02d}" / "mat"
    if not mat_dir.exists():
        return None
    files = sorted(mat_dir.glob("*.mat"), key=lambda p: p.stat().st_mtime, reverse=True)
    return files[0] if files else None


def label_from_effects(theta: float, gamma: float, ripple: float) -> str:
    vals = np.asarray([theta, gamma, ripple], dtype=float)
    pos = np.maximum(vals, 0.0)
    max_effect = float(np.max(pos))
    if max_effect < MIN_EFFECT_Z:
        return "inactive"
    active = pos >= max(MIN_EFFECT_Z, RELATIVE_ACTIVE_FRAC * max_effect)
    t, g, r = [bool(x) for x in active]
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
    return "inactive"


def write_records(path: Path, records: list[ComponentRecord]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(ComponentRecord.__dataclass_fields__.keys()))
        writer.writeheader()
        for rec in records:
            writer.writerow(rec.__dict__)


def write_readme(bands: dict[str, BandEvents], baseline_mask: np.ndarray, dt: float) -> None:
    lines = [
        "# E10gb1 P5 Dimred Component Response to P2 Band Events",
        "",
        "Primary source for band labels is P2 event_detection, not P1 spectrogram power.",
        "",
        f"Dataset: `{DATASET}`",
        f"P5 condition: `{CONDITION}`",
        f"P2 event file: `{P2_EVENT_MAT}`",
        f"P5 root: `{P5_ROOT}`",
        f"Sample dt: `{dt}` sec",
        "",
        "Band event coverage:",
    ]
    for name in BAND_ORDER:
        b = bands[name]
        lines.append(
            f"- {name}: {b.passband_hz[0]:g}-{b.passband_hz[1]:g} Hz, "
            f"raw channel-events={b.n_raw_channel_events}, unique centers={b.centers.size}, "
            f"mask coverage={100*b.coverage_fraction:.3f}%"
        )
    lines.extend(
        [
            f"- non-event baseline coverage={100*np.mean(baseline_mask):.3f}%",
            "",
            "Score definition:",
            "",
            "`effect_z = (mean(abs(component) inside band event windows) - mean(abs(component) inside non-event baseline)) / sd(non-event baseline)`",
            "",
            "The visual label is a soft exploratory label from theta/gamma/ripple effect_z. The CSV keeps the numeric values.",
        ]
    )
    (OUT_ROOT / "README_p2_band_event_response.md").write_text("\n".join(lines), encoding="utf-8")


def grouped(records: list[ComponentRecord]) -> dict[tuple[str, int], list[ComponentRecord]]:
    out: dict[tuple[str, int], list[ComponentRecord]] = {}
    for rec in records:
        out.setdefault((rec.method, rec.k), []).append(rec)
    return out


def ordered_records(records: list[ComponentRecord]) -> list[ComponentRecord]:
    method_order = {m: i for i, m in enumerate(METHODS)}
    return sorted(records, key=lambda r: (method_order.get(r.method, 99), r.k, r.component))


def short_label(label: str) -> str:
    return {
        "theta": "T",
        "gamma": "G",
        "ripple": "R",
        "ripple_gamma": "RG",
        "theta_gamma": "TG",
        "theta_ripple": "TR",
        "broad": "TGR",
        "inactive": "-",
    }[label]


def plot_event_coverage(bands: dict[str, BandEvents], baseline_mask: np.ndarray) -> None:
    labels = BAND_ORDER + ["non-event baseline"]
    values = [100 * bands[name].coverage_fraction for name in BAND_ORDER] + [100 * np.mean(baseline_mask)]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(labels, values, color=[COLORS["theta"], COLORS["gamma"], COLORS["ripple"], "#777777"])
    ax.set_ylabel("sample coverage (%)")
    ax.set_title(f"{DATASET}: P2 event window coverage used for component response")
    ax.grid(axis="y", color="0.9")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "00_p2_band_event_window_coverage.png", dpi=180)
    plt.close(fig)


def plot_effect_heatmap(records: list[ComponentRecord]) -> None:
    rows = ordered_records(records)
    mat = np.asarray(
        [[r.theta_effect_z, r.gamma_effect_z, r.ripple_effect_z] for r in rows],
        dtype=float,
    )
    max_abs = float(np.nanpercentile(np.abs(mat), 98))
    vmax = max(0.15, min(1.5, max_abs))
    fig, ax = plt.subplots(figsize=(8.5, max(10, 0.13 * len(rows))))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(3), BAND_ORDER)
    ax.set_yticks(range(len(rows)), [f"{r.method}_k{r.k:02d}_C{r.component}" for r in rows], fontsize=6)
    ax.set_title(f"{DATASET} P5 standardized csplit: response to P2 band event windows")
    cbar = fig.colorbar(im, ax=ax, shrink=0.65)
    cbar.set_label("event effect z vs non-event baseline")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "01_component_p2_band_event_effect_heatmap.png", dpi=180)
    plt.close(fig)


def plot_label_grid(records: list[ComponentRecord]) -> None:
    by = grouped(records)
    row_keys = [(m, k) for m in METHODS for k in KS if (m, k) in by]
    max_k = max((len(by[key]) for key in row_keys), default=8)
    fig, ax = plt.subplots(figsize=(1.15 * max_k + 3.0, 0.58 * len(row_keys) + 1.8))
    ax.set_xlim(-1.35, max_k)
    ax.set_ylim(-0.5, len(row_keys) - 0.5)
    ax.invert_yaxis()
    ax.axis("off")
    for y, key in enumerate(row_keys):
        method, k = key
        ax.text(-1.25, y, f"{method}_k{k:02d}", va="center", ha="left", fontsize=10, weight="bold")
        for x, rec in enumerate(sorted(by[key], key=lambda r: r.component)):
            rect = plt.Rectangle(
                (x, y - 0.32),
                0.88,
                0.64,
                facecolor=COLORS[rec.label],
                edgecolor="white",
                linewidth=0.8,
            )
            ax.add_patch(rect)
            ax.text(x + 0.44, y, f"C{rec.component}\n{short_label(rec.label)}", ha="center", va="center", fontsize=7.5)
    add_legend(ax)
    ax.set_title(f"{DATASET} P5 standardized csplit: component labels from P2 band events", fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "02_p2_band_event_label_grid_by_method_k.png", dpi=180)
    plt.close(fig)


def plot_label_composition(records: list[ComponentRecord]) -> None:
    by = grouped(records)
    row_keys = [(m, k) for m in METHODS for k in KS if (m, k) in by]
    counts = np.zeros((len(row_keys), len(LABEL_ORDER)), dtype=float)
    for i, key in enumerate(row_keys):
        labels = [r.label for r in by[key]]
        for j, label in enumerate(LABEL_ORDER):
            counts[i, j] = labels.count(label)
    frac = counts / np.maximum(1, counts.sum(axis=1, keepdims=True))
    fig, ax = plt.subplots(figsize=(12, max(5, 0.28 * len(row_keys) + 1.5)))
    left = np.zeros(len(row_keys))
    y = np.arange(len(row_keys))
    for j, label in enumerate(LABEL_ORDER):
        ax.barh(y, frac[:, j], left=left, color=COLORS[label], edgecolor="white", height=0.8, label=label)
        left += frac[:, j]
    ax.set_yticks(y, [f"{m}_k{k:02d}" for m, k in row_keys])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of components")
    ax.set_title(f"{DATASET}: P2-band-event label composition by method-k")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "03_p2_band_event_label_composition_by_method_k.png", dpi=180)
    plt.close(fig)


def plot_theta_vs_ripple_gamma(records: list[ComponentRecord]) -> None:
    fig, ax = plt.subplots(figsize=(8.5, 7.2))
    for label in LABEL_ORDER:
        subset = [r for r in records if r.label == label]
        if not subset:
            continue
        x = [r.theta_effect_z for r in subset]
        y = [max(r.gamma_effect_z, r.ripple_effect_z) for r in subset]
        ax.scatter(x, y, s=42, alpha=0.78, color=COLORS[label], edgecolor="white", linewidth=0.5, label=label)
    ax.axvline(MIN_EFFECT_Z, color="0.7", linestyle="--", linewidth=1)
    ax.axhline(MIN_EFFECT_Z, color="0.7", linestyle="--", linewidth=1)
    ax.set_xlabel("theta event effect z")
    ax.set_ylabel("ripple/gamma event effect z = max(gamma, ripple)")
    ax.set_title(f"{DATASET}: P2 event-based theta vs ripple/gamma component response")
    ax.grid(True, color="0.9")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "04_theta_vs_ripple_gamma_p2_event_effect_scatter.png", dpi=180)
    plt.close(fig)


def plot_two_subprocess_map(records: list[ComponentRecord]) -> None:
    by = grouped(records)
    matrix = np.full((len(METHODS), len(KS)), np.nan)
    annot = [["" for _ in KS] for _ in METHODS]
    for i, method in enumerate(METHODS):
        for j, k in enumerate(KS):
            rs = by.get((method, k), [])
            if not rs:
                continue
            theta = max((r for r in rs if r.label == "theta"), key=lambda r: r.theta_effect_z, default=None)
            rg = max(
                (r for r in rs if r.label in {"ripple_gamma", "gamma", "ripple"}),
                key=lambda r: max(r.gamma_effect_z, r.ripple_effect_z),
                default=None,
            )
            matrix[i, j] = (1 if theta else 0) + (2 if rg else 0)
            parts = []
            if theta:
                parts.append(f"T:C{theta.component}")
            if rg:
                parts.append(f"RG:C{rg.component}")
            annot[i][j] = "\n".join(parts) if parts else "-"
    cmap = plt.matplotlib.colors.ListedColormap(["#eeeeee", COLORS["theta"], COLORS["ripple_gamma"], "#225ea8"])
    norm = plt.matplotlib.colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    fig, ax = plt.subplots(figsize=(9, 4.4))
    im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(KS)), [f"k{k:02d}" for k in KS])
    ax.set_yticks(range(len(METHODS)), METHODS)
    ax.set_title(f"{DATASET}: candidate theta + ripple/gamma subprocess pair from P2 events")
    for i in range(len(METHODS)):
        for j in range(len(KS)):
            ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=9)
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2, 3], shrink=0.8)
    cbar.ax.set_yticklabels(["none", "theta only", "RG only", "both"])
    fig.tight_layout()
    fig.savefig(FIG_DIR / "05_two_subprocess_candidate_map_p2_event_detection.png", dpi=180)
    plt.close(fig)


def plot_perievent_examples(records: list[ComponentRecord], bands: dict[str, BandEvents], dt: float) -> None:
    examples = choose_examples(records, limit=4)
    if not examples:
        return
    fig, axes = plt.subplots(len(examples), 1, figsize=(10, 2.6 * len(examples)), sharex=True)
    if len(examples) == 1:
        axes = [axes]
    half_sec = 1.0
    half = int(round(half_sec / dt))
    x = np.arange(-half, half + 1) * dt
    for ax, rec in zip(axes, examples):
        activity = load_component_activity(rec)
        for band_name in BAND_ORDER:
            centers = bands[band_name].centers
            avg = mean_perievent(activity, centers, half, max_events=1000)
            ax.plot(x, avg, linewidth=1.5, label=band_name, color=COLORS[band_name])
        ax.axvline(0, color="0.4", linestyle="--", linewidth=0.9)
        ax.set_ylabel("abs component")
        ax.set_title(
            f"{rec.method}_k{rec.k:02d} C{rec.component} | label={rec.label} | "
            f"T/G/R effect={rec.theta_effect_z:.2f}/{rec.gamma_effect_z:.2f}/{rec.ripple_effect_z:.2f}",
            fontsize=10,
        )
        ax.grid(True, color="0.9")
    axes[-1].set_xlabel("time from P2 event peak (s)")
    axes[0].legend(frameon=False, ncol=3, loc="upper right")
    fig.suptitle(f"{DATASET}: peri-event mean activity for selected P5 components", fontsize=14, weight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "06_selected_component_perievent_mean_by_band.png", dpi=180)
    plt.close(fig)


def choose_examples(records: list[ComponentRecord], limit: int) -> list[ComponentRecord]:
    selected: list[ComponentRecord] = []
    groups = [
        ("theta", lambda r: r.label == "theta", lambda r: r.theta_effect_z),
        ("ripple_gamma", lambda r: r.label == "ripple_gamma", lambda r: max(r.gamma_effect_z, r.ripple_effect_z)),
        ("gamma", lambda r: r.label == "gamma", lambda r: r.gamma_effect_z),
        ("ripple", lambda r: r.label == "ripple", lambda r: r.ripple_effect_z),
    ]
    for _, pred, score in groups:
        candidates = [r for r in records if pred(r)]
        if candidates:
            selected.append(max(candidates, key=score))
        if len(selected) >= limit:
            break
    return selected[:limit]


def load_component_activity(rec: ComponentRecord) -> np.ndarray:
    with h5py.File(rec.source_file, "r") as f:
        dset = f["result/core/temporal_components_time_by_comp"]
        return np.abs(np.asarray(dset[rec.component - 1, :], dtype=np.float64))


def mean_perievent(activity: np.ndarray, centers: np.ndarray, half: int, max_events: int) -> np.ndarray:
    valid = centers[(centers >= half) & (centers < activity.size - half)]
    if valid.size == 0:
        return np.full(2 * half + 1, np.nan)
    if valid.size > max_events:
        idx = np.linspace(0, valid.size - 1, max_events).astype(int)
        valid = valid[idx]
    out = np.zeros(2 * half + 1, dtype=np.float64)
    for c in valid:
        out += activity[c - half : c + half + 1]
    return out / max(1, valid.size)


def add_legend(ax) -> None:
    handles = [
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS[label], markeredgecolor="white", markersize=10)
        for label in LABEL_ORDER
    ]
    ax.legend(handles, LABEL_ORDER, loc="lower left", bbox_to_anchor=(0.0, -0.12), ncol=4, frameon=False, fontsize=8)


if __name__ == "__main__":
    main()
