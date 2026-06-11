#!/usr/bin/env python3
"""Signed ROI tests for roi_mean efun-theta vs deconv-RG targets.

This script focuses on two checks requested after the mean_abs ROI summary was
found too coarse:

1. Selected BOLD unit overlap:
   Does efun x theta_selective select the same or nearby BOLD units as
   deconv_efun x ripple_gamma_no_theta?

2. Signed ROI / top ROI set:
   Do those two targets differ when ROI profiles are represented as signed
   real_mean?  mean_abs/positive_real/negative_real remain available as
   legacy/QC modes through --roi-value-modes, but the mainline default is
   real_mean only.

Scope is the current standardized complex-split roi_mean SOP.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
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
    "deconv_rg": ("deconv_efun", "ripple_gamma_no_theta"),
}
TARGET_PAIR = ("efun_theta", "deconv_rg")
ALL_ROI_VALUE_MODES = ("mean_abs", "real_mean", "positive_real", "negative_real")
ROI_VALUE_MODES = ("real_mean",)
TOP_NS = (5, 10, 20)
MAIN_TOP_N = 10
TOP_ROIS = 10


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
    parser.add_argument("--processed-root", type=Path, default=Path("/mnt/e/DataPons_processed"))
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results")
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "roi_mean_signed_target_mode_tests",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("/mnt/e/DataPons_processed")
        / "summary_figures"
        / "pipeline11_current_analysis_summary"
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "roi_mean_signed_target_mode_tests",
    )
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS))
    parser.add_argument("--top-ns", nargs="+", type=int, default=list(TOP_NS))
    parser.add_argument("--main-top-n", type=int, default=MAIN_TOP_N)
    parser.add_argument("--top-rois", type=int, default=TOP_ROIS)
    parser.add_argument(
        "--roi-value-modes",
        "--roi-value-mode",
        dest="roi_value_modes",
        nargs="+",
        choices=ALL_ROI_VALUE_MODES,
        default=list(ROI_VALUE_MODES),
        help=(
            "ROI profile value modes to compute. The current mainline default is "
            "real_mean; mean_abs/positive_real/negative_real are legacy/QC modes."
        ),
    )
    return parser.parse_args()


def main() -> None:
    global ROI_VALUE_MODES
    args = parse_args()
    ROI_VALUE_MODES = tuple(args.roi_value_modes)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    datasets = {d.lower() for d in args.datasets}

    roi_order, dataset_roi_order, bold_posts = read_p7_metadata(args.p7_roi_csv, datasets)
    hits_by_topn = {top_n: read_hits(args.hits_csv, datasets, top_n) for top_n in args.top_ns}

    profile_rows: list[dict[str, object]] = []
    overlap_rows: list[dict[str, object]] = []
    top_roi_rows: list[dict[str, object]] = []
    contrast_rows: list[dict[str, object]] = []

    p8_cache: dict[str, tuple[dict[int, np.ndarray], dict[int, np.ndarray], str]] = {}
    p10_cache: dict[tuple[str, str, str, str], tuple[dict[int, np.ndarray], str]] = {}

    for top_n, hits in hits_by_topn.items():
        unit_sets = collect_selected_units(hits)
        overlap_rows.extend(compute_unit_overlap(unit_sets, top_n))

        profiles = aggregate_target_profiles(
            hits,
            top_n,
            roi_order,
            dataset_roi_order,
            bold_posts,
            args.processed_root,
            p8_cache,
            p10_cache,
        )
        for key, payload in profiles.items():
            pipeline, dataset, target, mode = key
            vec = payload["profile"]
            for i, roi in enumerate(roi_order):
                value = vec[i]
                if math.isfinite(value):
                    profile_rows.append(
                        {
                            "top_n": top_n,
                            "pipeline": pipeline,
                            "dataset": dataset,
                            "target": target,
                            "feature_family": TARGETS[target][0],
                            "strict_label": TARGETS[target][1],
                            "roi_value_mode": mode,
                            "roi_index_global": i + 1,
                            "roi_label": roi,
                            "roi_value": value,
                            "weight_sum": payload["weight_sum"],
                            "n_hit_rows": payload["n_hit_rows"],
                        }
                    )

        top_roi_rows.extend(compute_top_roi_sets(profiles, roi_order, top_n, args.top_rois))
        contrast_rows.extend(compute_target_contrasts(profiles, top_n, args.top_rois))

    write_csv(args.output_dir / "selected_bold_unit_overlap_by_target.csv", overlap_rows)
    write_csv(args.output_dir / "signed_target_roi_profiles_long.csv", profile_rows)
    write_csv(args.output_dir / "signed_target_top_roi_sets.csv", top_roi_rows)
    write_csv(args.output_dir / "signed_target_roi_contrasts.csv", contrast_rows)

    plot_mode_overlap(pd.DataFrame(overlap_rows), args.figure_dir, args.main_top_n)
    plot_profile_corr(pd.DataFrame(contrast_rows), args.figure_dir, args.main_top_n)
    plot_top_roi_jaccard(pd.DataFrame(contrast_rows), args.figure_dir, args.main_top_n)
    plot_top_roi_heatmap(pd.DataFrame(top_roi_rows), args.figure_dir, args.main_top_n)
    write_readme(args.output_dir, args.figure_dir, overlap_rows, contrast_rows)

    print(f"Wrote {args.output_dir}")
    print(f"Figures {args.figure_dir}")


def read_p7_metadata(
    path: Path,
    datasets: set[str],
) -> tuple[list[str], dict[str, list[str]], dict[str, Path]]:
    usecols = ["dataset", "sorted_position", "roi_index", "roi_label", "bold_post_file"]
    df = pd.read_csv(path, usecols=usecols)
    df["dataset"] = df["dataset"].str.lower()
    df = df[df["dataset"].isin(datasets)].copy()
    df["sorted_position"] = pd.to_numeric(df["sorted_position"], errors="coerce")
    df["roi_index"] = pd.to_numeric(df["roi_index"], errors="coerce")

    roi_order: list[str] = []
    seen: set[str] = set()
    dataset_roi_order: dict[str, list[str]] = {}
    bold_posts: dict[str, Path] = {}
    for dataset in [d.lower() for d in DATASETS if d.lower() in set(df["dataset"])]:
        one = df[(df["dataset"].eq(dataset)) & (df["sorted_position"].eq(1))].sort_values("roi_index")
        labels = one["roi_label"].astype(str).tolist()
        dataset_roi_order[dataset] = labels
        for roi in labels:
            if roi not in seen:
                roi_order.append(roi)
                seen.add(roi)
        post = one["bold_post_file"].dropna().astype(str)
        if not post.empty:
            bold_posts[dataset] = to_wsl_path(Path(post.iloc[0]))
    if not roi_order:
        raise ValueError(f"No ROI labels read from {path}")
    return roi_order, dataset_roi_order, bold_posts


def to_wsl_path(path: Path) -> Path:
    text = str(path)
    if text.startswith("/mnt/"):
        return Path(text)
    if len(text) >= 3 and text[1:3] == ":\\":
        drive = text[0].lower()
        rest = text[3:].replace("\\", "/")
        return Path(f"/mnt/{drive}/{rest}")
    return Path(text.replace("\\", "/"))


def read_hits(path: Path, datasets: set[str], top_n: int) -> pd.DataFrame:
    usecols = [
        "pipeline",
        "dataset",
        "run_tag",
        "bold_observable",
        "p9_feature",
        "p9_method_k",
        "source_level",
        "rank_in_source_csv",
        "density_class",
        "density_condition",
        "strict_label",
        "bold_feature_family",
        "bold_mode_index",
        "bold_component_index",
        "peak_abs_corr",
        "finite_peak",
    ]
    df = pd.read_csv(path, usecols=usecols, low_memory=False)
    df["dataset"] = df["dataset"].str.lower()
    df = df[df["dataset"].isin(datasets)].copy()
    df = df[
        (df["bold_observable"] == "roi_mean")
        & (df["density_class"] == "dimred_efun_density")
        & (df["density_condition"] == "csplit")
        & (df["bold_feature_family"].isin({v[0] for v in TARGETS.values()}))
        & (df["strict_label"].isin({v[1] for v in TARGETS.values()}))
    ].copy()
    df["rank_in_source_csv"] = pd.to_numeric(df["rank_in_source_csv"], errors="coerce")
    df["peak_abs_corr"] = pd.to_numeric(df["peak_abs_corr"], errors="coerce")
    df["bold_mode_index"] = pd.to_numeric(df["bold_mode_index"], errors="coerce")
    df["bold_component_index"] = pd.to_numeric(df["bold_component_index"], errors="coerce")
    df = df[(df["rank_in_source_csv"] <= top_n) & np.isfinite(df["peak_abs_corr"]) & (df["peak_abs_corr"] > 0)].copy()
    df["target"] = ""
    for target, (feature, label) in TARGETS.items():
        mask = df["bold_feature_family"].eq(feature) & df["strict_label"].eq(label)
        df.loc[mask, "target"] = target
    return df[df["target"].ne("")].copy()


def collect_selected_units(hits: pd.DataFrame) -> dict[tuple[str, str, str], set[tuple]]:
    out: dict[tuple[str, str, str], set[tuple]] = defaultdict(set)
    for _, row in hits.iterrows():
        pipeline = str(row["pipeline"])
        dataset = str(row["dataset"]).lower()
        target = str(row["target"])
        if pipeline == "P8":
            mode = finite_int(row["bold_mode_index"])
            if mode is not None:
                out[(pipeline, dataset, target)].add(("p8", str(row["run_tag"]), mode))
        elif pipeline == "P10":
            comp = finite_int(row["bold_component_index"])
            if comp is not None:
                out[(pipeline, dataset, target)].add(("p10", str(row["run_tag"]), str(row["p9_feature"]), str(row["p9_method_k"]), comp))
    return out


def compute_unit_overlap(unit_sets: dict[tuple[str, str, str], set[tuple]], top_n: int) -> list[dict[str, object]]:
    rows = []
    keys = sorted({(pipeline, dataset) for pipeline, dataset, _ in unit_sets})
    for pipeline, dataset in keys:
        a = unit_sets.get((pipeline, dataset, TARGET_PAIR[0]), set())
        b = unit_sets.get((pipeline, dataset, TARGET_PAIR[1]), set())
        if not a or not b:
            continue
        inter = a & b
        union = a | b
        adjacent = 0
        if pipeline == "P8":
            b_by_run: dict[str, list[int]] = defaultdict(list)
            for _, run, mode in b:
                b_by_run[run].append(int(mode))
            for _, run, mode in a:
                if any(abs(int(mode) - int(other)) <= 1 for other in b_by_run.get(run, [])):
                    adjacent += 1
        else:
            b_by_context: dict[tuple[str, str, str], list[int]] = defaultdict(list)
            for _, run, feature, method_k, comp in b:
                b_by_context[(run, feature, method_k)].append(int(comp))
            for _, run, feature, method_k, comp in a:
                if any(abs(int(comp) - int(other)) <= 1 for other in b_by_context.get((run, feature, method_k), [])):
                    adjacent += 1
        rows.append(
            {
                "top_n": top_n,
                "pipeline": pipeline,
                "dataset": dataset,
                "target_a": TARGET_PAIR[0],
                "target_b": TARGET_PAIR[1],
                "n_units_a": len(a),
                "n_units_b": len(b),
                "n_exact_overlap": len(inter),
                "exact_jaccard": len(inter) / len(union) if union else math.nan,
                "adjacent_or_exact_fraction_of_a": adjacent / len(a) if a else math.nan,
            }
        )
    return rows


def aggregate_target_profiles(
    hits: pd.DataFrame,
    top_n: int,
    roi_order: Sequence[str],
    dataset_roi_order: dict[str, list[str]],
    bold_posts: dict[str, Path],
    processed_root: Path,
    p8_cache: dict[str, tuple[dict[int, np.ndarray], dict[int, np.ndarray], str]],
    p10_cache: dict[tuple[str, str, str, str], tuple[dict[int, np.ndarray], str]],
) -> dict[tuple[str, str, str, str], dict[str, object]]:
    n_roi = len(roi_order)
    sums: dict[tuple[str, str, str, str], np.ndarray] = defaultdict(lambda: np.zeros(n_roi, dtype=float))
    weights: dict[tuple[str, str, str, str], float] = defaultdict(float)
    n_hits: dict[tuple[str, str, str, str], int] = defaultdict(int)

    p8 = hits[hits["pipeline"].eq("P8")]
    for keys, g in p8.groupby(["dataset", "run_tag", "target", "bold_mode_index"], dropna=False):
        dataset, run_tag, target, mode_raw = keys
        mode = finite_int(mode_raw)
        if mode is None:
            continue
        modes_by_sorted, _, _ = load_p8_sorted_modes(str(dataset), roi_order, bold_posts, dataset_roi_order, p8_cache)
        complex_vec = modes_by_sorted.get(mode)
        if complex_vec is None:
            continue
        weight = float(pd.to_numeric(g["peak_abs_corr"], errors="coerce").sum())
        for roi_mode in ROI_VALUE_MODES:
            key = ("P8", str(dataset), str(target), roi_mode)
            sums[key] += transform_roi_vector(complex_vec, roi_mode) * weight
            weights[key] += weight
            n_hits[key] += int(g.shape[0])

    p10 = hits[hits["pipeline"].eq("P10")]
    for keys, g in p10.groupby(["dataset", "run_tag", "p9_feature", "p9_method_k", "target", "bold_component_index"], dropna=False):
        dataset, run_tag, p9_feature, p9_method_k, target, comp_raw = keys
        comp = finite_int(comp_raw)
        if comp is None:
            continue
        comp_vectors, _ = load_p10_complex_component_vectors(
            processed_root,
            roi_order,
            dataset_roi_order.get(str(dataset), []),
            str(dataset),
            str(run_tag),
            str(p9_feature),
            str(p9_method_k),
            p10_cache,
        )
        complex_vec = comp_vectors.get(comp)
        if complex_vec is None:
            continue
        weight = float(pd.to_numeric(g["peak_abs_corr"], errors="coerce").sum())
        for roi_mode in ROI_VALUE_MODES:
            key = ("P10", str(dataset), str(target), roi_mode)
            sums[key] += transform_roi_vector(complex_vec, roi_mode) * weight
            weights[key] += weight
            n_hits[key] += int(g.shape[0])

    out: dict[tuple[str, str, str, str], dict[str, object]] = {}
    for key, vec_sum in sums.items():
        w = weights[key]
        if w <= 0:
            continue
        out[key] = {"profile": vec_sum / w, "weight_sum": w, "n_hit_rows": n_hits[key], "top_n": top_n}
    return out


def load_p8_sorted_modes(
    dataset: str,
    roi_order: Sequence[str],
    bold_posts: dict[str, Path],
    dataset_roi_order: dict[str, list[str]],
    cache: dict[str, tuple[dict[int, np.ndarray], dict[int, np.ndarray], str]],
) -> tuple[dict[int, np.ndarray], dict[int, np.ndarray], str]:
    if dataset in cache:
        return cache[dataset]
    path = bold_posts.get(dataset)
    if path is None:
        return {}, {}, ""
    with h5py.File(path, "r") as handle:
        labels = [strip_roi_suffix(x) for x in decode_matlab_cellstr(handle, "BOLD_POST/observable/observable_labels")]
        evalues = read_complex_dataset(handle, "BOLD_POST/EDMD_outputs/evalues").reshape(-1)
        modes = read_complex_dataset(handle, "BOLD_POST/EDMD_outputs/kpm_modes")
    if modes.shape[0] != len(labels) and modes.shape[1] == len(labels):
        modes = modes.T
    order = np.argsort(np.abs(evalues))[::-1]
    roi_to_global = {roi: i for i, roi in enumerate(roi_order)}

    def to_global(local_vec: np.ndarray) -> np.ndarray:
        vec = np.full(len(roi_order), np.nan, dtype=complex)
        for local_i, roi in enumerate(labels):
            global_i = roi_to_global.get(roi)
            if global_i is not None:
                vec[global_i] = local_vec[local_i]
        return vec

    sorted_vectors = {i + 1: to_global(np.asarray(modes[:, raw_zero]).reshape(-1)) for i, raw_zero in enumerate(order)}
    raw_vectors = {int(raw_zero + 1): to_global(np.asarray(modes[:, raw_zero]).reshape(-1)) for raw_zero in range(modes.shape[1])}
    cache[dataset] = (sorted_vectors, raw_vectors, str(path))
    return cache[dataset]


def load_p10_complex_component_vectors(
    processed_root: Path,
    roi_order: Sequence[str],
    local_roi_order: Sequence[str],
    dataset: str,
    run_tag: str,
    p9_feature: str,
    p9_method_k: str,
    cache: dict[tuple[str, str, str, str], tuple[dict[int, np.ndarray], str]],
) -> tuple[dict[int, np.ndarray], str]:
    key = (dataset, run_tag, p9_feature, p9_method_k)
    if key in cache:
        return cache[key]
    p9_file = find_p9_file(processed_root, dataset, run_tag, p9_feature, p9_method_k)
    if p9_file is None or not local_roi_order:
        cache[key] = ({}, "")
        return cache[key]
    with h5py.File(p9_file, "r") as handle:
        weights = np.asarray(handle["result/core/mode_weights_mode_by_comp"][()], dtype=float)
        modes = read_complex_dataset(handle, "result/data/kpm_modes_mode_by_dict")
    if weights.shape[0] > weights.shape[1]:
        weights = weights.T
    if modes.shape[0] != len(local_roi_order) and modes.shape[1] == len(local_roi_order):
        modes = modes.T
    if modes.shape[0] != len(local_roi_order):
        cache[key] = ({}, str(p9_file))
        return cache[key]
    if modes.shape[1] != weights.shape[1] and modes.shape[1] == weights.shape[0]:
        weights = weights.T
    if modes.shape[1] != weights.shape[1]:
        cache[key] = ({}, str(p9_file))
        return cache[key]
    roi_to_global = {roi: i for i, roi in enumerate(roi_order)}
    out: dict[int, np.ndarray] = {}
    for comp_zero in range(weights.shape[0]):
        local_complex = np.asarray(modes @ weights[comp_zero, :].reshape(-1, 1)).reshape(-1)
        vec = np.full(len(roi_order), np.nan, dtype=complex)
        for local_i, roi in enumerate(local_roi_order):
            global_i = roi_to_global.get(roi)
            if global_i is not None:
                vec[global_i] = local_complex[local_i]
        out[comp_zero + 1] = vec
    cache[key] = (out, str(p9_file))
    return cache[key]


def transform_roi_vector(vec: np.ndarray, mode: str) -> np.ndarray:
    if mode == "mean_abs":
        return np.abs(vec).astype(float)
    real = np.real(vec).astype(float)
    if mode == "real_mean":
        return real
    if mode == "positive_real":
        return np.maximum(real, 0.0)
    if mode == "negative_real":
        return np.maximum(-real, 0.0)
    raise ValueError(mode)


def compute_top_roi_sets(
    profiles: dict[tuple[str, str, str, str], dict[str, object]],
    roi_order: Sequence[str],
    top_n: int,
    top_rois: int,
) -> list[dict[str, object]]:
    rows = []
    for (pipeline, dataset, target, mode), payload in sorted(profiles.items()):
        values = np.asarray(payload["profile"], dtype=float)
        if mode == "real_mean":
            order = np.argsort(np.abs(np.nan_to_num(values, nan=0.0)))[::-1]
        else:
            order = np.argsort(np.nan_to_num(values, nan=-np.inf))[::-1]
        for rank, idx in enumerate(order[:top_rois], start=1):
            rows.append(
                {
                    "top_n": top_n,
                    "pipeline": pipeline,
                    "dataset": dataset,
                    "target": target,
                    "roi_value_mode": mode,
                    "rank": rank,
                    "roi_label": roi_order[int(idx)],
                    "roi_value": float(values[int(idx)]),
                    "weight_sum": payload["weight_sum"],
                    "n_hit_rows": payload["n_hit_rows"],
                }
            )
    return rows


def compute_target_contrasts(
    profiles: dict[tuple[str, str, str, str], dict[str, object]],
    top_n: int,
    top_rois: int,
) -> list[dict[str, object]]:
    rows = []
    contexts = sorted({(pipeline, dataset, mode) for pipeline, dataset, _, mode in profiles})
    for pipeline, dataset, mode in contexts:
        a_key = (pipeline, dataset, TARGET_PAIR[0], mode)
        b_key = (pipeline, dataset, TARGET_PAIR[1], mode)
        if a_key not in profiles or b_key not in profiles:
            continue
        a = np.asarray(profiles[a_key]["profile"], dtype=float)
        b = np.asarray(profiles[b_key]["profile"], dtype=float)
        corr = pearson(a, b)
        set_a = top_roi_set(a, mode, top_rois)
        set_b = top_roi_set(b, mode, top_rois)
        rows.append(
            {
                "top_n": top_n,
                "pipeline": pipeline,
                "dataset": dataset,
                "roi_value_mode": mode,
                "target_a": TARGET_PAIR[0],
                "target_b": TARGET_PAIR[1],
                "profile_corr": corr,
                "top_roi_jaccard": jaccard(set_a, set_b),
                "n_top_roi_overlap": len(set_a & set_b),
                "top_rois": top_rois,
                "weight_sum_a": profiles[a_key]["weight_sum"],
                "weight_sum_b": profiles[b_key]["weight_sum"],
                "n_hit_rows_a": profiles[a_key]["n_hit_rows"],
                "n_hit_rows_b": profiles[b_key]["n_hit_rows"],
            }
        )
    return rows


def top_roi_set(values: np.ndarray, mode: str, n: int) -> set[int]:
    values = np.asarray(values, dtype=float)
    if mode == "real_mean":
        score = np.abs(np.nan_to_num(values, nan=0.0))
    else:
        score = np.nan_to_num(values, nan=-np.inf)
    order = np.argsort(score)[::-1]
    return {int(i) for i in order[:n]}


def plot_mode_overlap(df: pd.DataFrame, figure_dir: Path, top_n: int) -> None:
    sub = df[df["top_n"].eq(top_n)].copy()
    if sub.empty:
        return
    datasets = sorted(sub["dataset"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(15, max(4, 0.35 * len(datasets))), sharey=True)
    for ax, pipeline in zip(axes, ["P8", "P10"]):
        one = sub[sub["pipeline"].eq(pipeline)].set_index("dataset")
        vals = [float(one.loc[d, "exact_jaccard"]) if d in one.index else np.nan for d in datasets]
        adj = [float(one.loc[d, "adjacent_or_exact_fraction_of_a"]) if d in one.index else np.nan for d in datasets]
        y = np.arange(len(datasets))
        ax.barh(y - 0.16, vals, height=0.3, label="exact Jaccard", color="#4c78a8")
        ax.barh(y + 0.16, adj, height=0.3, label="same/adjacent frac", color="#f58518")
        ax.set_xlim(0, 1)
        ax.set_title(pipeline)
        ax.set_xlabel("selected BOLD unit overlap")
        ax.grid(axis="x", alpha=0.25)
        ax.set_yticks(y)
        ax.set_yticklabels(datasets)
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle(f"efun-theta vs deconv-RG selected BOLD unit overlap | top{top_n}")
    fig.tight_layout()
    fig.savefig(figure_dir / f"01_selected_bold_unit_overlap__top{top_n}.png", dpi=180)
    plt.close(fig)


def plot_profile_corr(df: pd.DataFrame, figure_dir: Path, top_n: int) -> None:
    sub = df[df["top_n"].eq(top_n)].copy()
    if sub.empty:
        return
    datasets = sorted(sub["dataset"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(16, max(4, 0.35 * len(datasets))), sharey=True)
    for ax, pipeline in zip(axes, ["P8", "P10"]):
        one = sub[sub["pipeline"].eq(pipeline)]
        x = np.arange(len(datasets))
        width = 0.18
        for j, mode in enumerate(ROI_VALUE_MODES):
            vals = []
            for d in datasets:
                hit = one[(one["dataset"].eq(d)) & (one["roi_value_mode"].eq(mode))]
                vals.append(float(hit["profile_corr"].iloc[0]) if not hit.empty else np.nan)
            ax.barh(x + (j - 1.5) * width, vals, height=width, label=mode)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlim(-1, 1)
        ax.set_title(pipeline)
        ax.set_xlabel("ROI profile corr")
        ax.grid(axis="x", alpha=0.25)
        ax.set_yticks(x)
        ax.set_yticklabels(datasets)
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle(f"efun-theta vs deconv-RG ROI profile correlation by value mode | top{top_n}")
    fig.tight_layout()
    fig.savefig(figure_dir / f"02_signed_roi_profile_corr__top{top_n}.png", dpi=180)
    plt.close(fig)


def plot_top_roi_jaccard(df: pd.DataFrame, figure_dir: Path, top_n: int) -> None:
    sub = df[df["top_n"].eq(top_n)].copy()
    if sub.empty:
        return
    datasets = sorted(sub["dataset"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(16, max(4, 0.35 * len(datasets))), sharey=True)
    for ax, pipeline in zip(axes, ["P8", "P10"]):
        one = sub[sub["pipeline"].eq(pipeline)]
        y = np.arange(len(datasets))
        width = 0.18
        for j, mode in enumerate(ROI_VALUE_MODES):
            vals = []
            for d in datasets:
                hit = one[(one["dataset"].eq(d)) & (one["roi_value_mode"].eq(mode))]
                vals.append(float(hit["top_roi_jaccard"].iloc[0]) if not hit.empty else np.nan)
            ax.barh(y + (j - 1.5) * width, vals, height=width, label=mode)
        ax.set_xlim(0, 1)
        ax.set_title(pipeline)
        ax.set_xlabel(f"top{TOP_ROIS} ROI-set Jaccard")
        ax.grid(axis="x", alpha=0.25)
        ax.set_yticks(y)
        ax.set_yticklabels(datasets)
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle(f"efun-theta vs deconv-RG top ROI set overlap by value mode | top{top_n}")
    fig.tight_layout()
    fig.savefig(figure_dir / f"03_signed_top_roi_set_jaccard__top{top_n}.png", dpi=180)
    plt.close(fig)


def plot_top_roi_heatmap(df: pd.DataFrame, figure_dir: Path, top_n: int) -> None:
    sub = df[(df["top_n"].eq(top_n)) & (df["roi_value_mode"].isin(["positive_real", "negative_real"]))].copy()
    sub = sub[sub["target"].isin(TARGET_PAIR)]
    if sub.empty:
        return
    sub["context"] = sub["pipeline"] + " " + sub["target"] + " " + sub["roi_value_mode"]
    # Only show rank 1-10 ROI labels as compact text table.
    contexts = ["P8 efun_theta positive_real", "P8 deconv_rg positive_real", "P8 efun_theta negative_real", "P8 deconv_rg negative_real",
                "P10 efun_theta positive_real", "P10 deconv_rg positive_real", "P10 efun_theta negative_real", "P10 deconv_rg negative_real"]
    datasets = sorted(sub["dataset"].unique())
    fig, axes = plt.subplots(len(contexts), 1, figsize=(20, max(8, 0.45 * len(contexts) * len(datasets))), sharex=True)
    if len(contexts) == 1:
        axes = [axes]
    for ax, context in zip(axes, contexts):
        one = sub[sub["context"].eq(context)]
        ax.set_axis_off()
        ax.set_title(context, loc="left", fontsize=10)
        lines = []
        for dataset in datasets:
            labels = (
                one[one["dataset"].eq(dataset)]
                .sort_values("rank")
                .head(TOP_ROIS)["roi_label"]
                .astype(str)
                .tolist()
            )
            if labels:
                lines.append(f"{dataset:<8} " + ", ".join(labels[:TOP_ROIS]))
        ax.text(0.0, 0.95, "\n".join(lines), ha="left", va="top", family="monospace", fontsize=7)
    fig.suptitle(f"Top signed ROI labels by target | top{top_n}")
    fig.tight_layout()
    fig.savefig(figure_dir / f"04_top_signed_roi_labels_by_target__top{top_n}.png", dpi=180)
    plt.close(fig)


def write_readme(output_dir: Path, figure_dir: Path, overlap_rows: Sequence[dict[str, object]], contrast_rows: Sequence[dict[str, object]]) -> None:
    overlap = pd.DataFrame(overlap_rows)
    contrast = pd.DataFrame(contrast_rows)
    lines = [
        "# roi_mean signed target mode tests",
        "",
        "Targets:",
        "- `efun_theta` = efun x theta_selective",
        "- `deconv_rg` = deconv_efun x ripple_gamma_no_theta",
        "",
        "ROI value modes:",
        "- `mean_abs`: magnitude baseline",
        "- `real_mean`: signed real part",
        "- `positive_real`: positive real part only",
        "- `negative_real`: negative real part magnitude",
        "",
    ]
    for top_n in sorted(contrast["top_n"].unique()) if not contrast.empty else []:
        sub = contrast[contrast["top_n"].eq(top_n)]
        lines.append(f"## top{top_n} summary")
        for pipeline in ["P8", "P10"]:
            one = sub[sub["pipeline"].eq(pipeline)]
            if one.empty:
                continue
            lines.append(f"{pipeline}:")
            for mode in ROI_VALUE_MODES:
                m = one[one["roi_value_mode"].eq(mode)]
                if m.empty:
                    continue
                lines.append(
                    f"  {mode:<13} profile_corr_median={m['profile_corr'].median():+.3f} "
                    f"topROI_jaccard_median={m['top_roi_jaccard'].median():.3f}"
                )
    lines.extend(
        [
            "",
            "Main figures:",
            f"- `{figure_dir / f'01_selected_bold_unit_overlap__top{MAIN_TOP_N}.png'}`",
            f"- `{figure_dir / f'02_signed_roi_profile_corr__top{MAIN_TOP_N}.png'}`",
            f"- `{figure_dir / f'03_signed_top_roi_set_jaccard__top{MAIN_TOP_N}.png'}`",
            f"- `{figure_dir / f'04_top_signed_roi_labels_by_target__top{MAIN_TOP_N}.png'}`",
        ]
    )
    (output_dir / "README_roi_mean_signed_target_mode_tests.md").write_text("\n".join(lines), encoding="utf-8")


def find_p9_file(processed_root: Path, dataset: str, run_tag: str, p9_feature: str, p9_method_k: str) -> Path | None:
    mat_dir = processed_root / dataset / "pipeline9_bold_eigenfunction_reduction" / run_tag / p9_feature / p9_method_k / "mat"
    if not mat_dir.is_dir():
        return None
    mats = sorted(mat_dir.glob("*.mat"), key=lambda p: p.stat().st_mtime, reverse=True)
    return mats[0] if mats else None


def strip_roi_suffix(label: str) -> str:
    return label[: -len("_mean")] if label.endswith("_mean") else label


def decode_matlab_string(handle: h5py.File, obj: object) -> str:
    if isinstance(obj, h5py.Reference):
        obj = handle[obj]
    if isinstance(obj, h5py.Dataset):
        data = obj[()]
    else:
        data = np.asarray(obj)
    if data.dtype.kind in "uif":
        vals = np.asarray(data).astype(np.uint16).ravel(order="F")
        return "".join(chr(int(v)) for v in vals if int(v) != 0)
    return str(data)


def decode_matlab_cellstr(handle: h5py.File, path: str) -> list[str]:
    data = handle[path][()]
    return [decode_matlab_string(handle, ref) for ref in np.asarray(data).ravel(order="F")]


def read_complex_dataset(handle: h5py.File, path: str) -> np.ndarray:
    data = handle[path][()]
    if data.dtype.fields and {"real", "imag"}.issubset(data.dtype.fields):
        return np.asarray(data["real"], dtype=float) + 1j * np.asarray(data["imag"], dtype=float)
    return np.asarray(data, dtype=float)


def finite_int(value: object) -> int | None:
    try:
        out = int(float(value))
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def pearson(a: Sequence[float], b: Sequence[float]) -> float:
    x = np.asarray(a, dtype=float)
    y = np.asarray(b, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if int(mask.sum()) < 3:
        return math.nan
    x = x[mask] - np.mean(x[mask])
    y = y[mask] - np.mean(y[mask])
    denom = float(np.sqrt(np.sum(x * x) * np.sum(y * y)))
    if denom <= 0:
        return math.nan
    return float(np.sum(x * y) / denom)


def jaccard(a: set[int], b: set[int]) -> float:
    return len(a & b) / len(a | b) if a or b else math.nan


def write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = []
        seen = set()
        for row in rows:
            for key in row:
                if key not in seen:
                    fields.append(key)
                    seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


if __name__ == "__main__":
    main()
