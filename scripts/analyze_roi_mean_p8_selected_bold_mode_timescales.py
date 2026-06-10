#!/usr/bin/env python3
"""P8 selected BOLD-mode eigenvalue/timescale checks for roi_mean SOP.

This asks whether P8 hits involving BOLD efun and deconv_efun select different
P7 BOLD modes in eigenvalue/timescale space.  P8 is the clean case because
`bold_mode_index` refers to a sorted P7 BOLD Koopman mode.  P10 selected units
are reduced P9 components and need a separate weighted-source-mode analysis.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Iterable, Sequence

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DATASETS = (
    "e10gb1",
    "e10fV1",
    "e10gh1",
    "e10gw1",
    "f12m01",
    "f12m02",
    "f12m03",
    "f12m05",
    "k13m17",
    "k13m18",
    "k13m19",
    "k13m20",
    "k13m21",
    "k13m23",
)
TARGETS = {
    "efun_theta": ("efun", "theta_selective"),
    "efun_rg": ("efun", "ripple_gamma_no_theta"),
    "deconv_theta": ("deconv_efun", "theta_selective"),
    "deconv_rg": ("deconv_efun", "ripple_gamma_no_theta"),
}
TARGET_ORDER = tuple(TARGETS)
TOP_NS = (3, 5, 10, 20)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--hits-csv",
        type=Path,
        default=Path("results")
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "p8_p10_strict_band_coupling"
        / "p8_p10_strict_band_hits_long.csv",
    )
    parser.add_argument(
        "--p7-roi-csv",
        type=Path,
        default=Path("results")
        / "pipeline_roi_profile_consistency_k03_k16_20260607"
        / "p7_roi_mean_profiles_direct_from_bold_post.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results")
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "roi_mean_p8_selected_bold_mode_timescales",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("/mnt/e/DataPons_processed")
        / "summary_figures"
        / "pipeline11_current_analysis_summary"
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "roi_mean_p8_selected_bold_mode_timescales",
    )
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS))
    parser.add_argument("--top-ns", nargs="+", type=int, default=list(TOP_NS))
    parser.add_argument(
        "--p5-background-csv",
        type=Path,
        default=Path("results")
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "p8_p10_strict_band_coupling"
        / "p5_strict_band_label_background.csv",
    )
    parser.add_argument(
        "--require-p5-pair-pass",
        action="store_true",
        help="Keep only dataset x density_method_k with both theta_selective and ripple_gamma_no_theta in P5.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    datasets = {x.lower() for x in args.datasets}

    p7_modes = read_p7_mode_metadata(args.p7_roi_csv, datasets)
    p5_pass = read_p5_pair_pass(args.p5_background_csv, datasets) if args.require_p5_pair_pass else None
    all_rows: list[dict[str, object]] = []
    for top_n in args.top_ns:
        hits = read_p8_hits(args.hits_csv, datasets, top_n, p5_pass)
        all_rows.extend(annotate_hits(hits, p7_modes, top_n))

    selected = pd.DataFrame(all_rows)
    write_csv(args.output_dir / "p8_selected_bold_mode_timescales_long.csv", all_rows)
    summary_rows = summarize(selected)
    write_csv(args.output_dir / "p8_selected_bold_mode_timescale_summary.csv", summary_rows)
    plot_all(selected, args.figure_dir)
    write_readme(args.output_dir, args.figure_dir, selected, pd.DataFrame(summary_rows), args.require_p5_pair_pass)
    print(f"Wrote {args.output_dir}")
    print(f"Figures {args.figure_dir}")


def read_p7_mode_metadata(path: Path, datasets: set[str]) -> pd.DataFrame:
    usecols = [
        "dataset",
        "bold_post_file",
        "sorted_position",
        "raw_index",
        "eigenvalue_real",
        "eigenvalue_imag",
        "eigenvalue_abs",
    ]
    df = pd.read_csv(path, usecols=usecols)
    df["dataset"] = df["dataset"].str.lower()
    df = df[df["dataset"].isin(datasets)].copy()
    df["sorted_position"] = pd.to_numeric(df["sorted_position"], errors="coerce").astype("Int64")
    df = df.drop_duplicates(["dataset", "sorted_position"]).copy()
    dt_by_file: dict[str, float] = {}
    dts = []
    for post_file in df["bold_post_file"].astype(str):
        if post_file not in dt_by_file:
            dt_by_file[post_file] = read_bold_post_dt(Path(post_file))
        dts.append(dt_by_file[post_file])
    df["dt_sec"] = dts
    return df


def read_p5_pair_pass(path: Path, datasets: set[str]) -> set[tuple[str, str]]:
    df = pd.read_csv(path, usecols=["dataset", "method_k", "strict_label"])
    df["dataset"] = df["dataset"].str.lower()
    df = df[df["dataset"].isin(datasets)].copy()
    labels = (
        df.groupby(["dataset", "method_k"])["strict_label"]
        .agg(lambda x: set(str(v) for v in x))
        .reset_index()
    )
    passed = labels[
        labels["strict_label"].apply(lambda s: "theta_selective" in s and "ripple_gamma_no_theta" in s)
    ]
    return set(zip(passed["dataset"].astype(str), passed["method_k"].astype(str)))


def read_bold_post_dt(path: Path) -> float:
    with h5py.File(path, "r") as handle:
        for key in (
            "BOLD_POST/dt",
            "BOLD_POST/EDMD_outputs/dt",
            "BOLD_POST/EDMD_outputs/sampling_period",
            "BOLD_POST/EDMD_outputs/sample_period",
            "BOLD_POST/params/default_dt",
        ):
            if key in handle:
                arr = np.asarray(handle[key][()])
                if arr.size:
                    value = float(arr.ravel(order="F")[0])
                    if math.isfinite(value) and value > 0:
                        return value
        if "BOLD_POST/EDMD_outputs/fs" in handle:
            fs = float(np.asarray(handle["BOLD_POST/EDMD_outputs/fs"][()]).ravel(order="F")[0])
            if math.isfinite(fs) and fs > 0:
                return 1.0 / fs
    return math.nan


def read_p8_hits(
    path: Path,
    datasets: set[str],
    top_n: int,
    p5_pass: set[tuple[str, str]] | None = None,
) -> pd.DataFrame:
    usecols = [
        "pipeline",
        "dataset",
        "run_tag",
        "bold_observable",
        "rank_in_source_csv",
        "density_class",
        "density_condition",
        "density_method_k",
        "density_index",
        "density_label",
        "strict_label",
        "component_timescale_sec",
        "bold_feature_family",
        "bold_mode_index",
        "peak_abs_corr",
        "peak_corr",
        "peak_lag_sec",
        "finite_peak",
    ]
    df = pd.read_csv(path, usecols=usecols, low_memory=False)
    df["dataset"] = df["dataset"].str.lower()
    df = df[
        df["dataset"].isin(datasets)
        & df["pipeline"].eq("P8")
        & df["bold_observable"].eq("roi_mean")
        & df["density_class"].eq("dimred_efun_density")
        & df["density_condition"].eq("csplit")
        & df["bold_feature_family"].isin(["efun", "deconv_efun"])
        & df["strict_label"].isin(["theta_selective", "ripple_gamma_no_theta"])
    ].copy()
    for col in ["rank_in_source_csv", "bold_mode_index", "peak_abs_corr", "peak_corr", "peak_lag_sec", "component_timescale_sec"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[(df["rank_in_source_csv"] <= top_n) & np.isfinite(df["peak_abs_corr"]) & (df["peak_abs_corr"] > 0)].copy()
    df["target"] = ""
    for target, (feature, label) in TARGETS.items():
        mask = df["bold_feature_family"].eq(feature) & df["strict_label"].eq(label)
        df.loc[mask, "target"] = target
    df = df[df["target"].ne("")].copy()
    if p5_pass is not None:
        keep = [
            (str(row.dataset).lower(), str(row.density_method_k)) in p5_pass
            for row in df.itertuples(index=False)
        ]
        df = df.loc[keep].copy()
        df["p5_pair_pass_filter"] = True
    else:
        df["p5_pair_pass_filter"] = False
    return df


def annotate_hits(hits: pd.DataFrame, p7_modes: pd.DataFrame, top_n: int) -> list[dict[str, object]]:
    if hits.empty:
        return []
    modes = p7_modes.rename(columns={"sorted_position": "bold_mode_index"}).copy()
    merged = hits.merge(modes, on=["dataset", "bold_mode_index"], how="left", suffixes=("", "_p7"))
    rows: list[dict[str, object]] = []
    for _, row in merged.iterrows():
        eig = complex(float_or_nan(row["eigenvalue_real"]), float_or_nan(row["eigenvalue_imag"]))
        dt = float_or_nan(row["dt_sec"])
        metrics = eigen_metrics(eig, dt)
        rows.append(
            {
                "top_n": top_n,
                "dataset": row["dataset"],
                "run_tag": row["run_tag"],
                "target": row["target"],
                "bold_feature_family": row["bold_feature_family"],
                "strict_label": row["strict_label"],
                "rank_in_source_csv": row["rank_in_source_csv"],
                "peak_abs_corr": row["peak_abs_corr"],
                "peak_corr": row["peak_corr"],
                "peak_lag_sec": row["peak_lag_sec"],
                "density_method_k": row.get("density_method_k", ""),
                "p5_pair_pass_filter": bool(row.get("p5_pair_pass_filter", False)),
                "density_index": row.get("density_index", ""),
                "density_label": row.get("density_label", ""),
                "component_timescale_sec": row.get("component_timescale_sec", math.nan),
                "bold_mode_index": row["bold_mode_index"],
                "bold_raw_index": row.get("raw_index", math.nan),
                "bold_post_file": row.get("bold_post_file", ""),
                "dt_sec": dt,
                **metrics,
            }
        )
    return rows


def eigen_metrics(eig: complex, dt: float) -> dict[str, float | str]:
    lam_abs = abs(eig)
    angle = math.atan2(eig.imag, eig.real)
    log_abs = math.log(lam_abs) if lam_abs > 0 and math.isfinite(lam_abs) else math.nan
    mag_tau = dt / abs(log_abs) if math.isfinite(dt) and dt > 0 and math.isfinite(log_abs) and abs(log_abs) > 0 else math.nan
    decay_tau = -dt / log_abs if math.isfinite(dt) and dt > 0 and math.isfinite(log_abs) and log_abs < 0 else math.nan
    growth_tau = dt / log_abs if math.isfinite(dt) and dt > 0 and math.isfinite(log_abs) and log_abs > 0 else math.nan
    period = 2 * math.pi * dt / abs(angle) if math.isfinite(dt) and dt > 0 and abs(angle) > 1e-12 else math.nan
    if math.isfinite(log_abs):
        stability = "growth" if log_abs > 0 else ("decay" if log_abs < 0 else "unit")
    else:
        stability = "unknown"
    return {
        "eigenvalue_real": eig.real,
        "eigenvalue_imag": eig.imag,
        "eigenvalue_abs": lam_abs,
        "eigenvalue_angle_rad": angle,
        "log_abs_lambda": log_abs,
        "magnitude_tau_sec": mag_tau,
        "decay_tau_sec": decay_tau,
        "growth_tau_sec": growth_tau,
        "oscillation_period_sec": period,
        "stability_class": stability,
    }


def summarize(df: pd.DataFrame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    if df.empty:
        return rows
    groupings = [
        ["top_n", "bold_feature_family"],
        ["top_n", "target"],
        ["top_n", "dataset", "bold_feature_family"],
        ["top_n", "dataset", "target"],
    ]
    for grouping in groupings:
        for keys, sub in df.groupby(grouping, dropna=False):
            if not isinstance(keys, tuple):
                keys = (keys,)
            row = dict(zip(grouping, keys))
            row.update(summary_stats(sub))
            rows.append(row)
    return rows


def summary_stats(sub: pd.DataFrame) -> dict[str, object]:
    return {
        "n_hit_rows": int(len(sub)),
        "n_datasets": int(sub["dataset"].nunique()),
        "n_unique_modes": int(sub[["dataset", "run_tag", "bold_mode_index"]].drop_duplicates().shape[0]),
        "median_bold_mode_index": median(sub["bold_mode_index"]),
        "median_eigenvalue_abs": median(sub["eigenvalue_abs"]),
        "median_log_abs_lambda": median(sub["log_abs_lambda"]),
        "median_magnitude_tau_sec": median_finite_clip(sub["magnitude_tau_sec"], 1e12),
        "median_component_timescale_sec": median(sub["component_timescale_sec"]),
        "growth_fraction": float(np.mean(sub["stability_class"].eq("growth"))) if len(sub) else math.nan,
        "median_peak_abs_corr": median(sub["peak_abs_corr"]),
        "median_abs_lag_sec": median(np.abs(pd.to_numeric(sub["peak_lag_sec"], errors="coerce"))),
    }


def median(values: Iterable[object]) -> float:
    arr = pd.to_numeric(pd.Series(list(values)), errors="coerce").dropna().to_numpy(dtype=float)
    if arr.size == 0:
        return math.nan
    return float(np.median(arr))


def median_finite_clip(values: Iterable[object], max_value: float) -> float:
    arr = pd.to_numeric(pd.Series(list(values)), errors="coerce").dropna().to_numpy(dtype=float)
    arr = arr[np.isfinite(arr) & (arr <= max_value)]
    if arr.size == 0:
        return math.nan
    return float(np.median(arr))


def float_or_nan(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def plot_all(df: pd.DataFrame, figure_dir: Path) -> None:
    if df.empty:
        return
    for top_n in sorted(df["top_n"].dropna().unique()):
        sub = df[df["top_n"].eq(top_n)].copy()
        plot_box(
            sub,
            "bold_feature_family",
            "magnitude_tau_sec",
            figure_dir / f"01_p8_selected_mode_magnitude_tau_by_feature__top{int(top_n)}.png",
            f"P8 selected BOLD mode magnitude timescale by feature | top{int(top_n)}",
            y_log=True,
        )
        plot_box(
            sub,
            "target",
            "magnitude_tau_sec",
            figure_dir / f"02_p8_selected_mode_magnitude_tau_by_target__top{int(top_n)}.png",
            f"P8 selected BOLD mode magnitude timescale by target | top{int(top_n)}",
            y_log=True,
        )
        plot_box(
            sub,
            "bold_feature_family",
            "bold_mode_index",
            figure_dir / f"03_p8_selected_mode_sorted_index_by_feature__top{int(top_n)}.png",
            f"P8 selected sorted BOLD mode index by feature | top{int(top_n)}",
            y_log=False,
        )
        plot_lambda_plane(
            sub,
            figure_dir / f"04_p8_selected_mode_lambda_plane__top{int(top_n)}.png",
            f"P8 selected BOLD eigenvalues | top{int(top_n)}",
        )
    plot_tau_by_topn(df, figure_dir / "05_p8_selected_mode_tau_by_topn_feature.png")


def plot_box(df: pd.DataFrame, group_col: str, value_col: str, path: Path, title: str, y_log: bool) -> None:
    order = [x for x in (TARGET_ORDER if group_col == "target" else ["efun", "deconv_efun"]) if x in set(df[group_col])]
    if not order:
        return
    data = []
    for label in order:
        vals = pd.to_numeric(df[df[group_col].eq(label)][value_col], errors="coerce").dropna().to_numpy(dtype=float)
        vals = vals[np.isfinite(vals) & (vals > 0 if y_log else np.isfinite(vals))]
        data.append(vals)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.boxplot(data, labels=order, showfliers=False)
    rng = np.random.default_rng(12345)
    for i, vals in enumerate(data, start=1):
        if vals.size:
            draw_vals = vals if vals.size <= 900 else rng.choice(vals, size=900, replace=False)
            ax.scatter(rng.normal(i, 0.055, size=draw_vals.size), draw_vals, s=7, alpha=0.18)
    if y_log:
        ax.set_yscale("log")
    ax.set_ylabel(value_col)
    ax.set_title(title)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_lambda_plane(df: pd.DataFrame, path: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7, 7))
    colors = {"efun": "#2b83ba", "deconv_efun": "#d95f02"}
    for feature, sub in df.groupby("bold_feature_family"):
        ax.scatter(sub["eigenvalue_real"], sub["eigenvalue_imag"], s=12, alpha=0.35, label=feature, color=colors.get(feature, "#777777"))
    theta = np.linspace(0, 2 * np.pi, 512)
    ax.plot(np.cos(theta), np.sin(theta), color="#999999", linewidth=1, linestyle="--")
    ax.axhline(0, color="#dddddd", linewidth=0.8)
    ax.axvline(0, color="#dddddd", linewidth=0.8)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("real(lambda)")
    ax.set_ylabel("imag(lambda)")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_tau_by_topn(df: pd.DataFrame, path: Path) -> None:
    rows = []
    for (top_n, feature), sub in df.groupby(["top_n", "bold_feature_family"]):
        rows.append({"top_n": top_n, "feature": feature, "tau": median_finite_clip(sub["magnitude_tau_sec"], 1e12)})
    summary = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(8, 5))
    for feature, sub in summary.groupby("feature"):
        sub = sub.sort_values("top_n")
        ax.plot(sub["top_n"], sub["tau"], marker="o", label=feature)
    ax.set_yscale("log")
    ax.set_xlabel("topN")
    ax.set_ylabel("median magnitude tau sec")
    ax.set_title("P8 selected BOLD mode tau by topN")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_readme(
    output_dir: Path,
    figure_dir: Path,
    selected: pd.DataFrame,
    summary: pd.DataFrame,
    require_p5_pair_pass: bool,
) -> None:
    lines = [
        "# P8 selected BOLD mode eigenvalue/timescale",
        "",
        "Scope: P8 roi_mean, dimred csplit density, strict theta_selective/RG-no-theta labels.",
        f"P5 pair-pass filter: `{require_p5_pair_pass}`.",
        "",
        "Timescale metric:",
        "- lambda is treated as a discrete Koopman eigenvalue from P7 BOLD_POST.",
        "- magnitude_tau_sec = dt / abs(log(abs(lambda))).",
        "- decay_tau_sec is only finite for abs(lambda)<1; growth_tau_sec is only finite for abs(lambda)>1.",
        "- Near-unit or slightly unstable modes should be interpreted as near-unit magnitude timescale, not pure decay.",
        "",
        f"Figures: `{figure_dir}`",
        "",
    ]
    if not summary.empty:
        lines.append("## Feature-level summary")
        cols = ["top_n", "bold_feature_family", "n_hit_rows", "n_unique_modes", "median_bold_mode_index", "median_magnitude_tau_sec", "growth_fraction", "median_peak_abs_corr"]
        feat = summary[[c for c in cols if c in summary.columns]].dropna(subset=["bold_feature_family"], how="any")
        for _, row in feat.sort_values(["top_n", "bold_feature_family"]).iterrows():
            lines.append(
                f"- top{int(row['top_n'])} {row['bold_feature_family']}: "
                f"n={int(row['n_hit_rows'])}, unique_modes={int(row['n_unique_modes'])}, "
                f"median_mode={float(row['median_bold_mode_index']):.1f}, "
                f"tau={float(row['median_magnitude_tau_sec']):.3g}s, "
                f"growth_frac={float(row['growth_fraction']):.2f}, "
                f"median|r|={float(row['median_peak_abs_corr']):.3f}"
            )
    output_dir.joinpath("README_p8_selected_bold_mode_timescales.md").write_text("\n".join(lines), encoding="utf-8")


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                fields.append(key)
                seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


if __name__ == "__main__":
    main()
