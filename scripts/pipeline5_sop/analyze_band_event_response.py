#!/usr/bin/env python3
"""P2-event response labels for P5 dimred components, SOP version.

Adds two things missing from the first preview:

1. Strict selectivity labels for the formal P5 band-selectivity analysis.
2. The current SOP activity transform:
   - adaptive_envelope: session-aware RMS envelope of abs(component), with the
     window derived from component weighted source-mode timescale.

The legacy `abs` transform is still available through `--transforms abs`, but
it is not part of the default SOP.

Soft dominant labels remain as a compatibility/debug field, but default figures
and downstream interpretation are strict-label only.
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
warnings.filterwarnings("ignore", message="This figure includes Axes that are not compatible with tight_layout*")
import matplotlib.pyplot as plt
import numpy as np


DATASET = "e10gb1"
CONDITION = "complex_split_projected_vlambda_standardize"

WORKSPACE = Path(__file__).resolve().parents[2]
PROCESSED_ROOT = Path("/mnt/e/DataPons_processed") if str(WORKSPACE).startswith("/mnt/") else Path("E:/DataPons_processed")
DATA_ROOT = PROCESSED_ROOT / DATASET
P2_EVENT_MAT = DATA_ROOT / "pipeline2_event_detection" / f"{DATASET}_bandpass_events_3bands.mat"
P5_ROOT = DATA_ROOT / "pipeline5_eigenfunction_reduction" / CONDITION

OUT_ROOT = WORKSPACE / "results" / f"{DATASET}_p5_p2_band_event_response_v2_selective_envelope_20260528"
FIG_DIR = OUT_ROOT / "figures"

METHODS = ["svd", "nmf", "mds", "umap"]
KS = list(range(3, 17))
BAND_ORDER = ["theta", "gamma", "ripple"]
AVAILABLE_TRANSFORMS = ["abs", "adaptive_envelope"]
TRANSFORMS = ["adaptive_envelope"]

DOMINANT_MIN_EFFECT_Z = 0.05
DOMINANT_RELATIVE_ACTIVE_FRAC = 0.60

STRICT_ACTIVE_Z = 0.50
STRICT_OFF_FRACTION = 0.60
STRICT_MARGIN_Z = 0.50

ENVELOPE_ALPHA = 0.50
ENVELOPE_MIN_SEC = 0.03
ENVELOPE_MAX_SEC = 1.00
ENVELOPE_FALLBACK_SEC = 0.10

DOMINANT_LABELS = [
    "theta",
    "gamma",
    "ripple",
    "ripple_gamma",
    "theta_gamma",
    "theta_ripple",
    "broad",
    "inactive",
]
STRICT_LABELS = [
    "theta_selective",
    "gamma_selective",
    "ripple_selective",
    "ripple_gamma_no_theta",
    "mixed_or_partial",
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
    "theta_selective": "#2ca25f",
    "gamma_selective": "#f0b429",
    "ripple_selective": "#3182bd",
    "ripple_gamma_no_theta": "#756bb1",
    "mixed_or_partial": "#a6611a",
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
    activity_transform: str
    method: str
    k: int
    component: int
    dominant_label: str
    strict_label: str
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
    theta_selectivity_margin_z: float
    ripple_gamma_selectivity_margin_z: float
    component_timescale_sec: float
    envelope_window_sec: float
    envelope_window_samples: int
    top_source_modes: str
    source_file: str


def main() -> None:
    args = parse_args()
    configure_paths(args)

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    bands, baseline_mask, dt, session_start, session_end = load_p2_band_events()
    records = compute_all_records(bands, baseline_mask, dt, session_start, session_end)
    write_records(OUT_ROOT / "component_band_event_response_v2.csv", records)
    write_readme(bands, baseline_mask, dt)

    if len(TRANSFORMS) > 1:
        plot_label_count_comparison(records)
    if not args.summary_only:
        plot_event_coverage(bands, baseline_mask)
    for transform in TRANSFORMS:
        subset = [r for r in records if r.activity_transform == transform]
        plot_two_subprocess_map(subset, transform)
        if not args.summary_only:
            plot_effect_heatmap(subset, transform)
            if args.include_dominant:
                plot_label_grid(subset, transform, label_kind="dominant")
            plot_label_grid(subset, transform, label_kind="strict")
            plot_label_composition(subset, transform, label_kind="strict")
            plot_theta_vs_ripple_gamma(subset, transform)

    print(f"Wrote {len(records)} records")
    print(f"Output: {OUT_ROOT}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="P2-event strict selectivity labels for P5 standardized csplit dimred components."
    )
    parser.add_argument("--dataset", default=DATASET, help="Dataset folder name under processed root.")
    parser.add_argument("--condition", default=CONDITION, help="P5 condition folder name.")
    parser.add_argument("--processed-root", default=str(PROCESSED_ROOT), help="Processed data root.")
    parser.add_argument("--workspace", default=str(WORKSPACE), help="Repository/workspace root.")
    parser.add_argument("--output-root", default="", help="Optional explicit output directory.")
    parser.add_argument("--component-counts", nargs="+", type=int, default=KS)
    parser.add_argument("--transforms", nargs="+", default=TRANSFORMS, choices=AVAILABLE_TRANSFORMS)
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="Only save the strict two-subprocess candidate map.",
    )
    parser.add_argument(
        "--include-dominant",
        action="store_true",
        help="Also export the older soft dominant-label grid. Default analysis is strict-label only.",
    )
    return parser.parse_args()


def configure_paths(args: argparse.Namespace) -> None:
    global DATASET, CONDITION, WORKSPACE, PROCESSED_ROOT, DATA_ROOT, P2_EVENT_MAT, P5_ROOT, OUT_ROOT, FIG_DIR, KS, TRANSFORMS

    DATASET = str(args.dataset)
    CONDITION = str(args.condition)
    WORKSPACE = Path(args.workspace)
    PROCESSED_ROOT = Path(args.processed_root)
    DATA_ROOT = PROCESSED_ROOT / DATASET
    P2_EVENT_MAT = DATA_ROOT / "pipeline2_event_detection" / f"{DATASET}_bandpass_events_3bands.mat"
    P5_ROOT = DATA_ROOT / "pipeline5_eigenfunction_reduction" / CONDITION
    if args.output_root:
        OUT_ROOT = Path(args.output_root)
    else:
        OUT_ROOT = WORKSPACE / "results" / f"{DATASET}_p5_p2_band_event_response_v2_selective_envelope_20260528"
    FIG_DIR = OUT_ROOT / "figures"
    KS = sorted({int(k) for k in args.component_counts})
    TRANSFORMS = list(dict.fromkeys(str(t) for t in args.transforms))


def load_p2_band_events() -> tuple[dict[str, BandEvents], np.ndarray, float, np.ndarray, np.ndarray]:
    with h5py.File(P2_EVENT_MAT, "r") as f:
        detect = f["R/DetectResults"]
        n_time = int(np.asarray(f["R/session_end_idx"]).squeeze()[-1])
        dt = float(np.asarray(f["R/dx"]).squeeze())
        session_start = np.asarray(f["R/session_start_idx"]).squeeze().astype(int) - 1
        session_end = np.asarray(f["R/session_end_idx"]).squeeze().astype(int)
        bands: dict[str, BandEvents] = {}

        for band_col in range(detect.shape[1]):
            mask = np.zeros(n_time, dtype=bool)
            centers_all: list[np.ndarray] = []
            n_raw = 0
            label = ""
            passband = (math.nan, math.nan)

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
    return bands, baseline_mask, dt, session_start, session_end


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


def compute_all_records(
    bands: dict[str, BandEvents],
    baseline_mask: np.ndarray,
    dt: float,
    session_start: np.ndarray,
    session_end: np.ndarray,
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
                activity_abs = np.abs(np.asarray(dset, dtype=np.float32))
                meta = component_timescale_metadata(f, dt)

            transform_to_activity = {
                "abs": activity_abs,
                "adaptive_envelope": rms_envelope_sessionwise(
                    activity_abs,
                    session_start,
                    session_end,
                    meta["envelope_window_samples"],
                ),
            }

            for transform, activity in transform_to_activity.items():
                records.extend(
                    score_activity(
                        activity,
                        transform,
                        method,
                        k,
                        mat_file,
                        band_masks,
                        baseline_mask,
                        meta,
                    )
                )
    return records


def find_p5_mat(method: str, k: int) -> Path | None:
    mat_dir = P5_ROOT / f"{method}_k{k:02d}" / "mat"
    if not mat_dir.exists():
        return None
    files = sorted(mat_dir.glob("*.mat"), key=lambda p: p.stat().st_mtime, reverse=True)
    return files[0] if files else None


def component_timescale_metadata(f: h5py.File, dt: float) -> dict[str, np.ndarray | list[str]]:
    discrete = h5_complex_vector(f["result/data/evalues_discrete"])
    continuous = h5_complex_vector(f["result/data/evalues_bilinear"])
    tau_discrete = timescale_from_discrete(discrete, dt)
    tau_continuous = timescale_from_continuous(continuous)
    tau = tau_discrete.copy()
    missing = np.isnan(tau)
    tau[missing] = tau_continuous[missing]

    # MATLAB HDF5 stores the 275 x C matrix as C x 275 for h5py.
    weights = np.asarray(f["result/core/mode_weights_norm_mode_by_comp"], dtype=np.float64)
    if weights.shape[1] != tau.size and weights.shape[0] == tau.size:
        weights = weights.T
    n_comp = weights.shape[0]

    comp_tau = np.full(n_comp, np.nan, dtype=float)
    top_modes: list[str] = []
    for c in range(n_comp):
        w = np.abs(weights[c, :])
        valid = np.isfinite(w) & (w > 0) & np.isfinite(tau) & (tau > 0)
        if np.any(valid):
            comp_tau[c] = weighted_quantile(tau[valid], w[valid], 0.5)
        top = np.argsort(w)[::-1][:5] + 1
        top_modes.append(",".join(str(int(x)) for x in top))

    win_sec = ENVELOPE_ALPHA * comp_tau
    invalid = ~np.isfinite(win_sec) | (win_sec <= 0)
    win_sec[invalid] = ENVELOPE_FALLBACK_SEC
    win_sec = np.clip(win_sec, ENVELOPE_MIN_SEC, ENVELOPE_MAX_SEC)
    win_samples = np.maximum(1, np.round(win_sec / dt).astype(int))
    return {
        "component_timescale_sec": comp_tau,
        "envelope_window_sec": win_sec,
        "envelope_window_samples": win_samples,
        "top_source_modes": top_modes,
    }


def h5_complex_vector(dataset) -> np.ndarray:
    arr = np.asarray(dataset)
    if arr.dtype.fields and "real" in arr.dtype.fields and "imag" in arr.dtype.fields:
        return arr["real"].squeeze() + 1j * arr["imag"].squeeze()
    return np.asarray(arr).squeeze().astype(complex)


def timescale_from_discrete(lambda_d: np.ndarray, dt: float) -> np.ndarray:
    rho = np.abs(lambda_d)
    tau = np.full(rho.shape, np.nan, dtype=float)
    valid = np.isfinite(rho) & (rho > 0) & (rho < 1) & np.isfinite(dt) & (dt > 0)
    tau[valid] = -float(dt) / np.log(rho[valid])
    unstable = np.isfinite(rho) & (rho >= 1)
    tau[unstable] = np.inf
    return tau


def timescale_from_continuous(lambda_c: np.ndarray) -> np.ndarray:
    real = np.real(lambda_c)
    tau = np.full(real.shape, np.nan, dtype=float)
    valid = np.isfinite(real) & (real < 0)
    tau[valid] = -1.0 / real[valid]
    return tau


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    total = np.sum(w)
    if total <= 0:
        return math.nan
    cw = np.cumsum(w) / total
    return float(v[np.searchsorted(cw, q, side="left")])


def rms_envelope_sessionwise(
    abs_activity: np.ndarray,
    session_start: np.ndarray,
    session_end: np.ndarray,
    window_samples: np.ndarray,
) -> np.ndarray:
    n_comp, _ = abs_activity.shape
    out = np.empty_like(abs_activity, dtype=np.float32)
    for c in range(n_comp):
        w = int(max(1, window_samples[c]))
        left = (w - 1) // 2
        right = w // 2
        for start, end in zip(session_start, session_end):
            start_i = int(start)
            end_i = int(end)
            seg = abs_activity[c, start_i:end_i].astype(np.float64)
            n = seg.size
            if n == 0:
                continue
            sq = seg * seg
            cs = np.concatenate(([0.0], np.cumsum(sq)))
            idx = np.arange(n)
            lo = np.maximum(0, idx - left)
            hi = np.minimum(n, idx + right + 1)
            denom = np.maximum(1, hi - lo)
            out[c, start_i:end_i] = np.sqrt((cs[hi] - cs[lo]) / denom).astype(np.float32)
    return out


def score_activity(
    activity: np.ndarray,
    transform: str,
    method: str,
    k: int,
    mat_file: Path,
    band_masks: dict[str, np.ndarray],
    baseline_mask: np.ndarray,
    meta: dict[str, np.ndarray | list[str]],
) -> list[ComponentRecord]:
    n_comp = activity.shape[0]
    baseline = activity[:, baseline_mask]
    baseline_mean = np.nanmean(baseline, axis=1)
    baseline_sd = np.nanstd(baseline, axis=1)
    baseline_sd[~np.isfinite(baseline_sd) | (baseline_sd <= 0)] = 1.0

    band_mean = {}
    effects = {}
    ratios = {}
    for name in BAND_ORDER:
        vals = activity[:, band_masks[name]]
        band_mean[name] = np.nanmean(vals, axis=1)
        effects[name] = (band_mean[name] - baseline_mean) / baseline_sd
        ratios[name] = band_mean[name] / np.maximum(baseline_mean, 1e-12)

    records: list[ComponentRecord] = []
    comp_tau = np.asarray(meta["component_timescale_sec"], dtype=float)
    win_sec = np.asarray(meta["envelope_window_sec"], dtype=float)
    win_samples = np.asarray(meta["envelope_window_samples"], dtype=int)
    top_modes = list(meta["top_source_modes"])

    for comp in range(n_comp):
        theta = float(effects["theta"][comp])
        gamma = float(effects["gamma"][comp])
        ripple = float(effects["ripple"][comp])
        theta_margin = theta - max(gamma, ripple)
        rg = max(gamma, ripple)
        rg_margin = rg - theta
        records.append(
            ComponentRecord(
                activity_transform=transform,
                method=method,
                k=k,
                component=comp + 1,
                dominant_label=dominant_label(theta, gamma, ripple),
                strict_label=strict_label(theta, gamma, ripple),
                theta_effect_z=theta,
                gamma_effect_z=gamma,
                ripple_effect_z=ripple,
                theta_event_mean=float(band_mean["theta"][comp]),
                gamma_event_mean=float(band_mean["gamma"][comp]),
                ripple_event_mean=float(band_mean["ripple"][comp]),
                baseline_mean=float(baseline_mean[comp]),
                baseline_sd=float(baseline_sd[comp]),
                theta_ratio=float(ratios["theta"][comp]),
                gamma_ratio=float(ratios["gamma"][comp]),
                ripple_ratio=float(ratios["ripple"][comp]),
                theta_selectivity_margin_z=float(theta_margin),
                ripple_gamma_selectivity_margin_z=float(rg_margin),
                component_timescale_sec=float(comp_tau[comp]),
                envelope_window_sec=float(win_sec[comp]),
                envelope_window_samples=int(win_samples[comp]),
                top_source_modes=top_modes[comp],
                source_file=str(mat_file),
            )
        )
    return records


def dominant_label(theta: float, gamma: float, ripple: float) -> str:
    vals = np.asarray([theta, gamma, ripple], dtype=float)
    pos = np.maximum(vals, 0.0)
    max_effect = float(np.max(pos))
    if max_effect < DOMINANT_MIN_EFFECT_Z:
        return "inactive"
    active = pos >= max(DOMINANT_MIN_EFFECT_Z, DOMINANT_RELATIVE_ACTIVE_FRAC * max_effect)
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


def strict_label(theta: float, gamma: float, ripple: float) -> str:
    vals = np.asarray([theta, gamma, ripple], dtype=float)
    if np.nanmax(vals) < STRICT_ACTIVE_Z:
        return "inactive"
    active = vals >= STRICT_ACTIVE_Z
    theta_low_for_rg = theta <= STRICT_OFF_FRACTION * max(gamma, ripple)
    rg_low_for_theta = max(gamma, ripple) <= STRICT_OFF_FRACTION * theta

    if active[0] and rg_low_for_theta and (theta - max(gamma, ripple) >= STRICT_MARGIN_Z):
        return "theta_selective"

    gamma_ripple_joint = active[1] and active[2]
    gamma_only = active[1] and not active[2]
    ripple_only = active[2] and not active[1]
    rg_margin_ok = max(gamma, ripple) - theta >= STRICT_MARGIN_Z

    if gamma_ripple_joint and theta_low_for_rg and rg_margin_ok:
        return "ripple_gamma_no_theta"
    if gamma_only and theta_low_for_rg and (gamma - max(theta, ripple) >= STRICT_MARGIN_Z):
        return "gamma_selective"
    if ripple_only and theta_low_for_rg and (ripple - max(theta, gamma) >= STRICT_MARGIN_Z):
        return "ripple_selective"
    return "mixed_or_partial"


def write_records(path: Path, records: list[ComponentRecord]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(ComponentRecord.__dataclass_fields__.keys()))
        writer.writeheader()
        for rec in records:
            writer.writerow(rec.__dict__)


def write_readme(bands: dict[str, BandEvents], baseline_mask: np.ndarray, dt: float) -> None:
    lines = [
        "# E10gb1 P5 Dimred Component Response to P2 Band Events, v2",
        "",
        "This version adds strict selectivity labels and adaptive RMS-envelope activity.",
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
            "Activity transforms:",
            "",
            "- `abs`: `abs(component)`.",
            "- `adaptive_envelope`: session-aware `sqrt(movmean(abs(component)^2, window))`.",
            "",
            "Adaptive envelope window:",
            "",
            "`component_tau = weighted_median(source_mode_tau, abs(component_loading))`",
            "",
            "`window_sec = clamp(0.5 * component_tau, 0.03 sec, 1.0 sec)`",
            "",
            "Strict label thresholds:",
            "",
            f"- active effect: `{STRICT_ACTIVE_Z}` baseline SD.",
            f"- off-target must be no more than `{STRICT_OFF_FRACTION}` of target.",
            f"- target minus off-target margin must be at least `{STRICT_MARGIN_Z}` baseline SD.",
            "",
            "Formal figures and downstream interpretation use `strict_label` only.",
            "`dominant_label` is retained only as a compatibility/debug CSV field.",
        ]
    )
    (OUT_ROOT / "README_p2_band_event_response_v2.md").write_text("\n".join(lines), encoding="utf-8")


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
        "theta_selective": "T-sel",
        "gamma_selective": "G-sel",
        "ripple_selective": "R-sel",
        "ripple_gamma_no_theta": "RG-noT",
        "mixed_or_partial": "mixed",
    }[label]


def plot_event_coverage(bands: dict[str, BandEvents], baseline_mask: np.ndarray) -> None:
    labels = BAND_ORDER + ["non-event baseline"]
    values = [100 * bands[name].coverage_fraction for name in BAND_ORDER] + [100 * np.mean(baseline_mask)]
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(labels, values, color=[COLORS["theta"], COLORS["gamma"], COLORS["ripple"], "#777777"])
    ax.set_ylabel("sample coverage (%)")
    ax.set_title(f"{DATASET}: P2 event window coverage")
    ax.grid(axis="y", color="0.9")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "00_p2_band_event_window_coverage.png", dpi=180)
    plt.close(fig)


def plot_label_count_comparison(records: list[ComponentRecord]) -> None:
    labels = STRICT_LABELS
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(15.0, 4.8))
    for i, transform in enumerate(TRANSFORMS):
        subset = [r for r in records if r.activity_transform == transform]
        counts = [sum(r.strict_label == label for r in subset) for label in labels]
        ax.bar(x + (i - 0.5) * width, counts, width=width, label=transform)
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set_ylabel("component count")
    ax.set_title(f"{DATASET}: strict P2-band selectivity label counts, abs vs adaptive envelope")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="0.9")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "00_strict_label_counts_abs_vs_adaptive_envelope.png", dpi=180)
    plt.close(fig)


def plot_effect_heatmap(records: list[ComponentRecord], transform: str) -> None:
    rows = ordered_records(records)
    mat = np.asarray([[r.theta_effect_z, r.gamma_effect_z, r.ripple_effect_z] for r in rows], dtype=float)
    vmax = max(0.25, min(2.5, float(np.nanpercentile(np.abs(mat), 98))))
    fig, ax = plt.subplots(figsize=(15.5, max(10, 0.13 * len(rows))))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(3), BAND_ORDER)
    ax.set_yticks(range(len(rows)), [f"{r.method}_k{r.k:02d}_C{r.component}" for r in rows], fontsize=6)
    ax.set_title(f"{DATASET} P5 standardized csplit: P2 band event effect z | {transform}")
    cbar = fig.colorbar(im, ax=ax, shrink=0.65)
    cbar.set_label("event effect z vs non-event baseline")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"01_{transform}_component_p2_band_event_effect_heatmap.png", dpi=180)
    plt.close(fig)


def plot_label_grid(records: list[ComponentRecord], transform: str, label_kind: str) -> None:
    by = grouped(records)
    row_keys = [(m, k) for m in METHODS for k in KS if (m, k) in by]
    max_k = max((len(by[key]) for key in row_keys), default=8)
    fig, ax = plt.subplots(figsize=(max(15.5, 1.65 * max_k + 4.5), 0.58 * len(row_keys) + 1.9))
    ax.set_xlim(-1.35, max_k)
    ax.set_ylim(-0.5, len(row_keys) - 0.5)
    ax.invert_yaxis()
    ax.axis("off")
    labels = DOMINANT_LABELS if label_kind == "dominant" else STRICT_LABELS
    for y, key in enumerate(row_keys):
        method, k = key
        ax.text(-1.25, y, f"{method}_k{k:02d}", va="center", ha="left", fontsize=10, weight="bold")
        for x, rec in enumerate(sorted(by[key], key=lambda r: r.component)):
            label = rec.dominant_label if label_kind == "dominant" else rec.strict_label
            rect = plt.Rectangle((x, y - 0.32), 0.88, 0.64, facecolor=COLORS[label], edgecolor="white", linewidth=0.8)
            ax.add_patch(rect)
            ax.text(x + 0.44, y, f"C{rec.component}\n{short_label(label)}", ha="center", va="center", fontsize=7.2)
    add_legend(ax, labels)
    ax.set_title(f"{DATASET}: {label_kind} P2-band labels | {transform}", fontsize=15, weight="bold")
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"02_{transform}_{label_kind}_label_grid_by_method_k.png", dpi=180)
    plt.close(fig)


def plot_label_composition(records: list[ComponentRecord], transform: str, label_kind: str) -> None:
    by = grouped(records)
    row_keys = [(m, k) for m in METHODS for k in KS if (m, k) in by]
    labels = STRICT_LABELS if label_kind == "strict" else DOMINANT_LABELS
    counts = np.zeros((len(row_keys), len(labels)), dtype=float)
    for i, key in enumerate(row_keys):
        row_labels = [
            r.strict_label if label_kind == "strict" else r.dominant_label
            for r in by[key]
        ]
        for j, label in enumerate(labels):
            counts[i, j] = row_labels.count(label)
    frac = counts / np.maximum(1, counts.sum(axis=1, keepdims=True))
    fig, ax = plt.subplots(figsize=(16.5, max(5, 0.28 * len(row_keys) + 1.5)))
    left = np.zeros(len(row_keys))
    y = np.arange(len(row_keys))
    for j, label in enumerate(labels):
        ax.barh(y, frac[:, j], left=left, color=COLORS[label], edgecolor="white", height=0.8, label=label)
        left += frac[:, j]
    ax.set_yticks(y, [f"{m}_k{k:02d}" for m, k in row_keys])
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of components")
    ax.set_title(f"{DATASET}: {label_kind} label composition | {transform}")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"03_{transform}_{label_kind}_label_composition_by_method_k.png", dpi=180)
    plt.close(fig)


def plot_theta_vs_ripple_gamma(records: list[ComponentRecord], transform: str) -> None:
    fig, ax = plt.subplots(figsize=(12.8, 7.2))
    for label in STRICT_LABELS:
        subset = [r for r in records if r.strict_label == label]
        if not subset:
            continue
        x = [r.theta_effect_z for r in subset]
        y = [max(r.gamma_effect_z, r.ripple_effect_z) for r in subset]
        ax.scatter(x, y, s=46, alpha=0.8, color=COLORS[label], edgecolor="white", linewidth=0.5, label=label)
    ax.axvline(STRICT_ACTIVE_Z, color="0.7", linestyle="--", linewidth=1)
    ax.axhline(STRICT_ACTIVE_Z, color="0.7", linestyle="--", linewidth=1)
    ax.plot([-1, 3], [-1, 3], color="0.8", linewidth=0.8)
    ax.set_xlabel("theta event effect z")
    ax.set_ylabel("ripple/gamma event effect z = max(gamma, ripple)")
    ax.set_title(f"{DATASET}: strict selectivity plane | {transform}")
    ax.grid(True, color="0.9")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"04_{transform}_strict_theta_vs_ripple_gamma_effect_scatter.png", dpi=180)
    plt.close(fig)


def plot_two_subprocess_map(records: list[ComponentRecord], transform: str) -> None:
    by = grouped(records)
    matrix = np.full((len(METHODS), len(KS)), np.nan)
    annot = [["" for _ in KS] for _ in METHODS]
    for i, method in enumerate(METHODS):
        for j, k in enumerate(KS):
            rs = by.get((method, k), [])
            if not rs:
                continue
            theta = max(
                (r for r in rs if r.strict_label == "theta_selective"),
                key=lambda r: r.theta_effect_z,
                default=None,
            )
            rg = max(
                (r for r in rs if r.strict_label in {"ripple_gamma_no_theta", "gamma_selective", "ripple_selective"}),
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
    cmap = plt.matplotlib.colors.ListedColormap(["#eeeeee", COLORS["theta_selective"], COLORS["ripple_gamma_no_theta"], "#225ea8"])
    norm = plt.matplotlib.colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    fig, ax = plt.subplots(figsize=(14.5, 4.8))
    im = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto")
    ax.set_xticks(range(len(KS)), [f"k{k:02d}" for k in KS])
    ax.set_yticks(range(len(METHODS)), METHODS)
    ax.set_title(f"{DATASET}: strict theta + ripple/gamma candidates | {transform}")
    for i in range(len(METHODS)):
        for j in range(len(KS)):
            ax.text(j, i, annot[i][j], ha="center", va="center", fontsize=9)
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 1, 2, 3], shrink=0.8)
    cbar.ax.set_yticklabels(["none", "theta only", "RG only", "both"])
    fig.tight_layout()
    fig.savefig(FIG_DIR / f"05_{transform}_strict_two_subprocess_candidate_map.png", dpi=180)
    plt.close(fig)


def add_legend(ax, labels: list[str]) -> None:
    handles = [
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=COLORS[label], markeredgecolor="white", markersize=10)
        for label in labels
    ]
    ax.legend(handles, labels, loc="lower left", bbox_to_anchor=(0.0, -0.13), ncol=4, frameon=False, fontsize=8)


if __name__ == "__main__":
    main()
