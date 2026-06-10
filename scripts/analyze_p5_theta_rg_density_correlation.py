#!/usr/bin/env python3
"""Correlation between P5 theta-selective and RG-no-theta dimred densities.

This script answers a narrow SOP question:

    Are the P5 theta-selective and ripple-gamma-no-theta density traces
    independent, or do they co-vary within the same dataset/method-k?

It intentionally reads the P5 thresholded density MAT files used downstream by
P8/P10, not the continuous reduction components.  The current mainline scope is
standardized complex-split with adaptive RMS-envelope density.
"""

from __future__ import annotations

import math
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


WORKSPACE = Path("/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
PROCESSED_ROOT = Path("/mnt/e/DataPons_processed")
LABEL_SUFFIX = "p5_p2_band_event_response_v2_selective_envelope_20260528"
DENSITY_CONDITION = "complex_split_projected_vlambda_standardize_rmsenv_adaptive"
RESULT_ROOT = WORKSPACE / "results" / "standardized_csplit_k03_k16_all_current_20260607" / "p5_theta_rg_density_correlation"
FIG_ROOT = (
    PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "standardized_csplit_k03_k16_all_current_20260607"
    / "p5_theta_rg_density_correlation"
)

METHODS = ["svd", "nmf", "mds", "umap"]
THETA_LABEL = "theta_selective"
RG_LABEL = "ripple_gamma_no_theta"
TRANSFORM = "adaptive_envelope"


def main() -> None:
    RESULT_ROOT.mkdir(parents=True, exist_ok=True)
    FIG_ROOT.mkdir(parents=True, exist_ok=True)

    label_files = sorted(WORKSPACE.glob(f"results/*_{LABEL_SUFFIX}/component_band_event_response_v2.csv"))
    if not label_files:
        raise FileNotFoundError(f"No P5 label CSVs matching *_{LABEL_SUFFIX}")

    all_pair_rows: list[dict] = []
    canonical_rows: list[dict] = []

    for label_file in label_files:
        dataset = label_file.parent.name[: -len("_" + LABEL_SUFFIX)]
        print(f"processing {dataset}", flush=True)
        labels = pd.read_csv(label_file)
        labels = labels[labels["activity_transform"].eq(TRANSFORM)].copy()
        labels["k"] = labels["k"].astype(int)
        labels["component"] = labels["component"].astype(int)
        labels = labels[labels["k"].between(3, 16)]

        for method in METHODS:
            sub_m = labels[labels["method"].eq(method)]
            for k in sorted(sub_m["k"].unique()):
                sub = sub_m[sub_m["k"].eq(k)]
                theta = sub[sub["strict_label"].eq(THETA_LABEL)].copy()
                rg = sub[sub["strict_label"].eq(RG_LABEL)].copy()
                if theta.empty or rg.empty:
                    continue

                density_file = find_density_file(dataset, method, int(k))
                if density_file is None:
                    print(f"  missing density {dataset} {method}_k{int(k):02d}", flush=True)
                    continue
                density, session_id = load_density(density_file)

                pair_rows = []
                for _, trow in theta.iterrows():
                    t_idx = int(trow["component"]) - 1
                    if t_idx < 0 or t_idx >= density.shape[0]:
                        continue
                    for _, rrow in rg.iterrows():
                        r_idx = int(rrow["component"]) - 1
                        if r_idx < 0 or r_idx >= density.shape[0]:
                            continue
                        full = pearson(density[t_idx], density[r_idx])
                        sess = grouped_demeaned_corr(density[t_idx], density[r_idx], session_id)
                        row = {
                            "dataset": dataset,
                            "activity_transform": TRANSFORM,
                            "method": method,
                            "k": int(k),
                            "method_k": f"{method}_k{int(k):02d}",
                            "theta_component": int(trow["component"]),
                            "rg_component": int(rrow["component"]),
                            "theta_effect_z": float(trow["theta_effect_z"]),
                            "rg_effect_z": float(max(rrow["gamma_effect_z"], rrow["ripple_effect_z"])),
                            "gamma_effect_z": float(rrow["gamma_effect_z"]),
                            "ripple_effect_z": float(rrow["ripple_effect_z"]),
                            "density_corr": full,
                            "session_demeaned_density_corr": sess,
                            "n_density_windows": int(density.shape[1]),
                            "density_file": str(density_file),
                        }
                        all_pair_rows.append(row)
                        pair_rows.append(row)

                if not pair_rows:
                    continue

                # Canonical pair: strongest theta-effect component and strongest
                # RG-effect component within the same method-k.
                best_theta = theta.sort_values("theta_effect_z", ascending=False).iloc[0]
                rg = rg.assign(rg_effect_z=np.maximum(rg["gamma_effect_z"].astype(float), rg["ripple_effect_z"].astype(float)))
                best_rg = rg.sort_values("rg_effect_z", ascending=False).iloc[0]
                t_idx = int(best_theta["component"]) - 1
                r_idx = int(best_rg["component"]) - 1
                if 0 <= t_idx < density.shape[0] and 0 <= r_idx < density.shape[0]:
                    canonical_rows.append(
                        {
                            "dataset": dataset,
                            "activity_transform": TRANSFORM,
                            "method": method,
                            "k": int(k),
                            "method_k": f"{method}_k{int(k):02d}",
                            "theta_component": int(best_theta["component"]),
                            "rg_component": int(best_rg["component"]),
                            "theta_effect_z": float(best_theta["theta_effect_z"]),
                            "rg_effect_z": float(best_rg["rg_effect_z"]),
                            "gamma_effect_z": float(best_rg["gamma_effect_z"]),
                            "ripple_effect_z": float(best_rg["ripple_effect_z"]),
                            "density_corr": pearson(density[t_idx], density[r_idx]),
                            "session_demeaned_density_corr": grouped_demeaned_corr(density[t_idx], density[r_idx], session_id),
                            "n_theta_components": int(theta.shape[0]),
                            "n_rg_components": int(rg.shape[0]),
                            "n_density_windows": int(density.shape[1]),
                            "density_file": str(density_file),
                        }
                    )

    all_pairs = pd.DataFrame(all_pair_rows)
    canonical = pd.DataFrame(canonical_rows)
    all_pairs.to_csv(RESULT_ROOT / "theta_rg_density_all_pair_correlations.csv", index=False)
    canonical.to_csv(RESULT_ROOT / "theta_rg_density_canonical_pair_correlations.csv", index=False)

    summary = summarize(canonical, all_pairs)
    summary.to_csv(RESULT_ROOT / "theta_rg_density_correlation_dataset_summary.csv", index=False)

    plot_distribution(all_pairs)
    plot_dataset_canonical(canonical)
    plot_methodk_heatmap(canonical)
    write_readme(summary)

    print(f"Wrote {RESULT_ROOT}")
    print(f"Figures {FIG_ROOT}")


def find_density_file(dataset: str, method: str, k: int) -> Path | None:
    mat_dir = (
        PROCESSED_ROOT
        / dataset
        / "pipeline5_dimred_thresholded_density"
        / DENSITY_CONDITION
        / f"{method}_k{k:02d}"
        / "mat"
    )
    if not mat_dir.exists():
        return None
    files = sorted(mat_dir.glob("*.mat"), key=lambda p: p.stat().st_mtime, reverse=True)
    return files[0] if files else None


def load_density(path: Path) -> tuple[np.ndarray, np.ndarray | None]:
    with h5py.File(path, "r") as f:
        if "D/density_time_by_mode" in f:
            density = np.asarray(f["D/density_time_by_mode"], dtype=np.float32)
        elif "density_time_by_mode" in f:
            density = np.asarray(f["density_time_by_mode"], dtype=np.float32)
        elif "density_time_by_component" in f:
            density = np.asarray(f["density_time_by_component"], dtype=np.float32)
        elif "density" in f:
            density = np.asarray(f["density"], dtype=np.float32)
        else:
            raise KeyError(f"No density matrix in {path}")

        session_id = None
        for key in ["D/window_session_id", "D/window_session_idx", "window_session_id", "window_session_idx"]:
            if key in f:
                session_id = np.asarray(f[key]).squeeze()
                break
    if density.shape[0] > density.shape[1]:
        density = density.T
    if session_id is not None and session_id.size != density.shape[1]:
        session_id = None
    return density, session_id


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    if np.sum(ok) < 3:
        return math.nan
    x = x[ok]
    y = y[ok]
    x = x - np.mean(x)
    y = y - np.mean(y)
    sx = float(np.sum(x * x))
    sy = float(np.sum(y * y))
    if sx <= 0 or sy <= 0:
        return math.nan
    return float(np.sum(x * y) / math.sqrt(sx * sy))


def grouped_demeaned_corr(x: np.ndarray, y: np.ndarray, group: np.ndarray | None) -> float:
    if group is None:
        return math.nan
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    group = np.asarray(group).squeeze()
    xs = []
    ys = []
    for g in np.unique(group[np.isfinite(group)]):
        idx = group == g
        ok = idx & np.isfinite(x) & np.isfinite(y)
        if np.sum(ok) < 3:
            continue
        xs.append(x[ok] - np.mean(x[ok]))
        ys.append(y[ok] - np.mean(y[ok]))
    if not xs:
        return math.nan
    return pearson(np.concatenate(xs), np.concatenate(ys))


def summarize(canonical: pd.DataFrame, all_pairs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source_name, df in [("canonical_pair", canonical), ("all_pairs", all_pairs)]:
        if df.empty:
            continue
        for dataset, sub in df.groupby("dataset", sort=True):
            rows.append(
                {
                    "source": source_name,
                    "dataset": dataset,
                    "n_pairs_or_methodk": int(sub.shape[0]),
                    "density_corr_median": float(sub["density_corr"].median()),
                    "density_corr_mean": float(sub["density_corr"].mean()),
                    "session_demeaned_density_corr_median": float(sub["session_demeaned_density_corr"].median()),
                    "session_demeaned_density_corr_mean": float(sub["session_demeaned_density_corr"].mean()),
                    "frac_session_demeaned_abs_lt_0p2": float((sub["session_demeaned_density_corr"].abs() < 0.2).mean()),
                    "frac_session_demeaned_abs_lt_0p5": float((sub["session_demeaned_density_corr"].abs() < 0.5).mean()),
                }
            )
    return pd.DataFrame(rows)


def plot_distribution(all_pairs: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(10, 4))
    bins = np.linspace(-1, 1, 41)
    ax.hist(all_pairs["density_corr"].dropna(), bins=bins, alpha=0.5, label="full density", color="#4c78a8")
    ax.hist(
        all_pairs["session_demeaned_density_corr"].dropna(),
        bins=bins,
        alpha=0.5,
        label="session demeaned",
        color="#f58518",
    )
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("theta/RG density correlation")
    ax.set_ylabel("all theta/RG component pairs")
    ax.set_title("P5 theta-selective vs RG-no-theta density correlation | adaptive envelope")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG_ROOT / "01_theta_rg_density_corr_distribution_all_pairs.png", dpi=180)
    plt.close(fig)


def plot_dataset_canonical(canonical: pd.DataFrame) -> None:
    ds_order = sorted(canonical["dataset"].unique())
    values = [
        canonical[canonical["dataset"].eq(ds)]["session_demeaned_density_corr"].dropna().to_numpy()
        for ds in ds_order
    ]
    fig, ax = plt.subplots(figsize=(12, max(5, 0.35 * len(ds_order))))
    ax.boxplot(values, vert=False, labels=ds_order, showfliers=False, patch_artist=True)
    for i, vals in enumerate(values, start=1):
        if vals.size:
            y = np.full(vals.size, i) + np.linspace(-0.08, 0.08, vals.size)
            ax.scatter(vals, y, s=12, alpha=0.45, color="#4c78a8")
    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(0.2, color="#999999", lw=0.8, ls="--")
    ax.axvline(-0.2, color="#999999", lw=0.8, ls="--")
    ax.set_xlim(-1, 1)
    ax.set_xlabel("canonical theta/RG session-demeaned density corr")
    ax.set_ylabel("dataset")
    ax.set_title("P5 canonical theta/RG density correlation by dataset | adaptive envelope")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIG_ROOT / "02_theta_rg_density_corr_by_dataset_canonical.png", dpi=180)
    plt.close(fig)


def plot_methodk_heatmap(canonical: pd.DataFrame) -> None:
    if canonical.empty:
        return
    canonical = canonical.copy()
    canonical["method_k"] = pd.Categorical(
        canonical["method_k"],
        categories=[f"{m}_k{k:02d}" for m in METHODS for k in range(3, 17)],
        ordered=True,
    )
    pivot = canonical.pivot_table(
        index="dataset",
        columns="method_k",
        values="session_demeaned_density_corr",
        aggfunc="median",
        observed=False,
    ).sort_index()
    fig, ax = plt.subplots(figsize=(18, max(5, 0.35 * len(pivot))))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=90, fontsize=7)
    ax.set_title("Canonical theta/RG session-demeaned density corr by dataset and method-k")
    cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cbar.set_label("corr")
    fig.tight_layout()
    fig.savefig(FIG_ROOT / "03_theta_rg_density_corr_methodk_heatmap_adaptive_envelope.png", dpi=180)
    plt.close(fig)


def write_readme(summary: pd.DataFrame) -> None:
    primary = summary[summary["source"].eq("canonical_pair")].copy()
    lines = [
        "# P5 theta/RG density correlation",
        "",
        "Primary question: are strict `theta_selective` and `ripple_gamma_no_theta`",
        "dimred density traces independent within each dataset/method-k?",
        "",
        "Density source:",
        f"`pipeline5_dimred_thresholded_density/{DENSITY_CONDITION}`",
        "",
        "Primary metric: canonical pair session-demeaned Pearson correlation.",
        "Session demeaning removes shared per-session density offsets before computing correlation.",
        "",
        "Dataset summary:",
        "",
        "```text",
    ]
    for _, r in primary.sort_values("dataset").iterrows():
        lines.append(
            f"{r['dataset']:<8} n={int(r['n_pairs_or_methodk']):>3} "
            f"full_med={r['density_corr_median']:+.3f} "
            f"sess_med={r['session_demeaned_density_corr_median']:+.3f} "
            f"frac_abs<0.2={r['frac_session_demeaned_abs_lt_0p2']:.2f}"
        )
    lines.extend(
        [
            "```",
            "",
            "Files:",
            "- `theta_rg_density_all_pair_correlations.csv`",
            "- `theta_rg_density_canonical_pair_correlations.csv`",
            "- `theta_rg_density_correlation_dataset_summary.csv`",
        ]
    )
    (RESULT_ROOT / "README_theta_rg_density_correlation.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
