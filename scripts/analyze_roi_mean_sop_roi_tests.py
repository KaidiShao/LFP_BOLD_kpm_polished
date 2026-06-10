"""ROI tests for the current standardized-csplit roi_mean SOP.

This script is intentionally narrower than the older ROI consistency scripts:

    BLP condition: standardized complex-split, adaptive-envelope density
    BOLD observable: roi_mean
    Coupling source: current p8_p10_strict_band_hits_long.csv

It asks three questions:

1. What ROI profile is selected by RG-no-theta top xcorr hits in each dataset?
2. Are RG-no-theta selected ROI profiles more cross-dataset consistent than an
   all-mode BOLD ROI baseline?
3. Are efun-selected and deconv-selected RG-no-theta ROI profiles spatially
   separable, or do they select similar BOLD ROI subspaces?

P8 maps top hits directly to P7 roi_mean BOLD modes.  P10 maps top hits to P9
BOLD dimred components, reconstructed from the P9 component weights and P7
roi_mean Koopman modes stored in the P9 MAT files.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from itertools import combinations, product
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
LABELS = ("theta_selective", "ripple_gamma_no_theta", "mixed_or_partial", "inactive")
FEATURES = ("efun", "deconv_efun")
PIPELINES = ("P8", "P10")
TOP_NS = (5, 10, 20)
MAIN_TOP_N = 10
ROI_PANEL_COLORS = {
    "efun": "#2c7fb8",
    "deconv_efun": "#d95f0e",
}


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
        / "roi_mean_rg_no_theta_formal_roi_tests",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("/mnt/e/DataPons_processed")
        / "summary_figures"
        / "pipeline11_current_analysis_summary"
        / "standardized_csplit_k03_k16_all_current_20260607"
        / "roi_mean_rg_no_theta_formal_roi_tests",
    )
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS))
    parser.add_argument("--top-ns", nargs="+", type=int, default=list(TOP_NS))
    parser.add_argument("--main-top-n", type=int, default=MAIN_TOP_N)
    parser.add_argument("--max-baseline-pairs", type=int, default=0)
    return parser.parse_args()


def safe_name(text: object) -> str:
    out = []
    for ch in str(text):
        if ch.isalnum() or ch in ("-", "_"):
            out.append(ch)
        else:
            out.append("_")
    return "".join(out).strip("_") or "blank"


def finite_float(value: object, default: float = math.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) else default


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


def read_p7_roi_profiles(
    path: Path,
    datasets: set[str],
) -> tuple[
    list[str],
    dict[str, list[str]],
    dict[tuple[str, int], np.ndarray],
    dict[str, dict[int, np.ndarray]],
]:
    usecols = ["dataset", "sorted_position", "roi_index", "roi_label", "roi_value"]
    df = pd.read_csv(path, usecols=usecols)
    df["dataset"] = df["dataset"].str.lower()
    df = df[df["dataset"].isin(datasets)].copy()
    df["sorted_position"] = pd.to_numeric(df["sorted_position"], errors="coerce").astype("Int64")
    df["roi_index"] = pd.to_numeric(df["roi_index"], errors="coerce").astype("Int64")
    df["roi_value"] = pd.to_numeric(df["roi_value"], errors="coerce")

    present_datasets = [d.lower() for d in DATASETS if d.lower() in set(df["dataset"])]
    if not present_datasets:
        raise ValueError(f"No requested datasets found in {path}")
    dataset_roi_order: dict[str, list[str]] = {}
    roi_order: list[str] = []
    seen: set[str] = set()
    for dataset in present_datasets:
        one = (
            df[(df["dataset"] == dataset) & (df["sorted_position"] == 1)]
            .sort_values("roi_index")
            [["roi_label"]]
        )
        labels = one["roi_label"].astype(str).tolist()
        dataset_roi_order[dataset] = labels
        for roi in labels:
            if roi not in seen:
                roi_order.append(roi)
                seen.add(roi)
    roi_to_pos = {roi: i for i, roi in enumerate(roi_order)}

    profiles: dict[tuple[str, int], np.ndarray] = {}
    by_dataset: dict[str, dict[int, np.ndarray]] = defaultdict(dict)
    for (dataset, mode), g in df.groupby(["dataset", "sorted_position"], sort=False):
        if pd.isna(mode):
            continue
        vec = np.full(len(roi_order), np.nan, dtype=float)
        for _, row in g.iterrows():
            pos = roi_to_pos.get(str(row["roi_label"]))
            if pos is None:
                continue
            vec[pos] = finite_float(row["roi_value"])
        mode_int = int(mode)
        profiles[(str(dataset), mode_int)] = vec
        by_dataset[str(dataset)][mode_int] = vec
    return roi_order, dataset_roi_order, profiles, dict(by_dataset)


def read_complex_dataset(handle: h5py.File, path: str) -> np.ndarray:
    data = handle[path][()]
    if data.dtype.fields and {"real", "imag"}.issubset(data.dtype.fields):
        return np.asarray(data["real"], dtype=float) + 1j * np.asarray(data["imag"], dtype=float)
    return np.asarray(data, dtype=float)


def find_p9_file(processed_root: Path, dataset: str, run_tag: str, p9_feature: str, p9_method_k: str) -> Path | None:
    mat_dir = processed_root / dataset / "pipeline9_bold_eigenfunction_reduction" / run_tag / p9_feature / p9_method_k / "mat"
    if not mat_dir.is_dir():
        return None
    mats = sorted(mat_dir.glob("*.mat"), key=lambda p: p.stat().st_mtime, reverse=True)
    return mats[0] if mats else None


def load_p10_component_vectors(
    processed_root: Path,
    roi_order: Sequence[str],
    local_roi_order: Sequence[str],
    dataset: str,
    run_tag: str,
    p9_feature: str,
    p9_method_k: str,
) -> tuple[dict[int, np.ndarray], str]:
    p9_file = find_p9_file(processed_root, dataset, run_tag, p9_feature, p9_method_k)
    if p9_file is None:
        return {}, ""
    with h5py.File(p9_file, "r") as handle:
        weights = np.asarray(handle["result/core/mode_weights_mode_by_comp"][()], dtype=float)
        modes = read_complex_dataset(handle, "result/data/kpm_modes_mode_by_dict")
    # h5py exposes MATLAB arrays transposed for many v7.3 files.  The expected
    # orientation here is weights = comp x source_mode and modes = roi x source_mode.
    if weights.shape[0] > weights.shape[1]:
        weights = weights.T
    if modes.shape[0] != len(local_roi_order) and modes.shape[1] == len(local_roi_order):
        modes = modes.T
    if modes.shape[0] != len(local_roi_order):
        raise ValueError(
            f"Unexpected P9 mode/ROI shape in {p9_file}: "
            f"modes={modes.shape}, local_n_roi={len(local_roi_order)}, global_n_roi={len(roi_order)}"
        )
    if modes.shape[1] != weights.shape[1] and modes.shape[1] == weights.shape[0]:
        weights = weights.T
    if modes.shape[1] != weights.shape[1]:
        raise ValueError(f"Unexpected P9 weight shape in {p9_file}: modes={modes.shape}, weights={weights.shape}")
    out: dict[int, np.ndarray] = {}
    roi_to_global = {roi: i for i, roi in enumerate(roi_order)}
    for comp_zero in range(weights.shape[0]):
        comp_complex = modes @ weights[comp_zero, :].reshape(-1, 1)
        local_vec = np.abs(np.asarray(comp_complex).reshape(-1)).astype(float)
        vec = np.full(len(roi_order), np.nan, dtype=float)
        for local_i, roi in enumerate(local_roi_order):
            global_i = roi_to_global.get(roi)
            if global_i is not None:
                vec[global_i] = local_vec[local_i]
        out[comp_zero + 1] = vec
    return out, str(p9_file)


def read_hits(path: Path, datasets: set[str], top_n: int) -> pd.DataFrame:
    usecols = [
        "pipeline",
        "dataset",
        "run_tag",
        "bold_observable",
        "p9_feature",
        "p9_method",
        "p9_k",
        "p9_method_k",
        "source_level",
        "rank_in_source_csv",
        "density_class",
        "density_condition",
        "density_method",
        "density_k",
        "density_method_k",
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
        & (df["bold_feature_family"].isin(FEATURES))
        & (df["strict_label"].isin(LABELS))
    ].copy()
    df["rank_in_source_csv"] = pd.to_numeric(df["rank_in_source_csv"], errors="coerce")
    df["peak_abs_corr"] = pd.to_numeric(df["peak_abs_corr"], errors="coerce")
    df["bold_mode_index"] = pd.to_numeric(df["bold_mode_index"], errors="coerce")
    df["bold_component_index"] = pd.to_numeric(df["bold_component_index"], errors="coerce")
    df = df[(df["rank_in_source_csv"] <= top_n) & np.isfinite(df["peak_abs_corr"]) & (df["peak_abs_corr"] > 0)].copy()
    return df


def aggregate_vectors(
    hits: pd.DataFrame,
    roi_order: Sequence[str],
    dataset_roi_order: dict[str, list[str]],
    p7_profiles: dict[tuple[str, int], np.ndarray],
    processed_root: Path,
) -> tuple[dict[tuple[str, str, str, str], np.ndarray], list[dict[str, object]], list[dict[str, object]]]:
    n_roi = len(roi_order)
    sums: dict[tuple[str, str, str, str], np.ndarray] = defaultdict(lambda: np.zeros(n_roi, dtype=float))
    weights: dict[tuple[str, str, str, str], float] = defaultdict(float)
    n_hits: dict[tuple[str, str, str, str], int] = defaultdict(int)
    selected_mode_rows: list[dict[str, object]] = []
    missing_rows: list[dict[str, object]] = []
    p10_cache: dict[tuple[str, str, str, str], tuple[dict[int, np.ndarray], str]] = {}

    # Collapse repeated top-hit rows that select the same BOLD unit before ROI aggregation.
    p8 = hits[hits["pipeline"] == "P8"].copy()
    if not p8.empty:
        grouped = (
            p8.groupby(["dataset", "run_tag", "bold_feature_family", "strict_label", "bold_mode_index"], dropna=False)
            .agg(weight=("peak_abs_corr", "sum"), n_hit_rows=("peak_abs_corr", "size"))
            .reset_index()
        )
        for _, row in grouped.iterrows():
            dataset = str(row["dataset"]).lower()
            mode = int(row["bold_mode_index"])
            vec = p7_profiles.get((dataset, mode))
            if vec is None:
                missing_rows.append({"pipeline": "P8", "dataset": dataset, "bold_unit": mode, "reason": "missing_p7_roi_profile"})
                continue
            key = ("P8", dataset, str(row["bold_feature_family"]), str(row["strict_label"]))
            w = finite_float(row["weight"], 0.0)
            sums[key] += np.nan_to_num(vec, nan=0.0) * w
            weights[key] += w
            n_hits[key] += int(row["n_hit_rows"])
            selected_mode_rows.append(
                {
                    "pipeline": "P8",
                    "dataset": dataset,
                    "feature_family": row["bold_feature_family"],
                    "strict_label": row["strict_label"],
                    "bold_unit_type": "p7_sorted_mode",
                    "bold_unit_index": mode,
                    "weight_sum": w,
                    "n_hit_rows": int(row["n_hit_rows"]),
                }
            )

    p10 = hits[hits["pipeline"] == "P10"].copy()
    if not p10.empty:
        grouped = (
            p10.groupby(
                [
                    "dataset",
                    "run_tag",
                    "p9_feature",
                    "p9_method",
                    "p9_k",
                    "p9_method_k",
                    "bold_feature_family",
                    "strict_label",
                    "bold_component_index",
                ],
                dropna=False,
            )
            .agg(weight=("peak_abs_corr", "sum"), n_hit_rows=("peak_abs_corr", "size"))
            .reset_index()
        )
        for _, row in grouped.iterrows():
            dataset = str(row["dataset"]).lower()
            run_tag = str(row["run_tag"])
            p9_feature = str(row["p9_feature"])
            p9_method_k = str(row["p9_method_k"])
            comp = int(row["bold_component_index"])
            cache_key = (dataset, run_tag, p9_feature, p9_method_k)
            if cache_key not in p10_cache:
                try:
                    local_roi_order = dataset_roi_order.get(dataset, list(roi_order))
                    p10_cache[cache_key] = load_p10_component_vectors(
                        processed_root,
                        roi_order,
                        local_roi_order,
                        dataset,
                        run_tag,
                        p9_feature,
                        p9_method_k,
                    )
                except Exception as exc:  # keep going; record missing context
                    p10_cache[cache_key] = ({}, "")
                    missing_rows.append(
                        {
                            "pipeline": "P10",
                            "dataset": dataset,
                            "p9_feature": p9_feature,
                            "p9_method_k": p9_method_k,
                            "reason": f"p9_load_error:{exc}",
                        }
                    )
            comp_vectors, p9_file = p10_cache[cache_key]
            vec = comp_vectors.get(comp)
            if vec is None:
                missing_rows.append(
                    {
                        "pipeline": "P10",
                        "dataset": dataset,
                        "p9_feature": p9_feature,
                        "p9_method_k": p9_method_k,
                        "bold_unit": comp,
                        "reason": "missing_p10_component_vector",
                    }
                )
                continue
            key = ("P10", dataset, str(row["bold_feature_family"]), str(row["strict_label"]))
            w = finite_float(row["weight"], 0.0)
            sums[key] += np.nan_to_num(vec, nan=0.0) * w
            weights[key] += w
            n_hits[key] += int(row["n_hit_rows"])
            selected_mode_rows.append(
                {
                    "pipeline": "P10",
                    "dataset": dataset,
                    "feature_family": row["bold_feature_family"],
                    "strict_label": row["strict_label"],
                    "bold_unit_type": "p9_dimred_component",
                    "bold_unit_index": comp,
                    "p9_feature": p9_feature,
                    "p9_method": row["p9_method"],
                    "p9_k": row["p9_k"],
                    "p9_method_k": p9_method_k,
                    "p9_file": p9_file,
                    "weight_sum": w,
                    "n_hit_rows": int(row["n_hit_rows"]),
                }
            )

    vectors: dict[tuple[str, str, str, str], np.ndarray] = {}
    profile_rows: list[dict[str, object]] = []
    for key, vec_sum in sums.items():
        w = weights[key]
        if w <= 0:
            continue
        vec = vec_sum / w
        vectors[key] = vec
        pipeline, dataset, feature, label = key
        for roi_idx, (roi, value) in enumerate(zip(roi_order, vec), start=1):
            profile_rows.append(
                {
                    "pipeline": pipeline,
                    "dataset": dataset,
                    "feature_family": feature,
                    "strict_label": label,
                    "roi_index": roi_idx,
                    "roi_label": roi,
                    "roi_value_weighted": value,
                    "weight_sum": w,
                    "n_hit_rows": n_hits[key],
                }
            )
    return vectors, profile_rows, selected_mode_rows + missing_rows


def vector_pair_records(
    vectors: dict[tuple[str, str, str, str], np.ndarray],
    datasets: Sequence[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for pipeline in PIPELINES:
        for feature_a in FEATURES:
            for feature_b in FEATURES:
                for label_a in LABELS:
                    for label_b in LABELS:
                        vals = []
                        same_dataset_vals = []
                        cross_dataset_vals = []
                        for ds_a in datasets:
                            va = vectors.get((pipeline, ds_a.lower(), feature_a, label_a))
                            if va is None:
                                continue
                            for ds_b in datasets:
                                vb = vectors.get((pipeline, ds_b.lower(), feature_b, label_b))
                                if vb is None:
                                    continue
                                r = pearson(va, vb)
                                if not math.isfinite(r):
                                    continue
                                vals.append(r)
                                if ds_a.lower() == ds_b.lower():
                                    same_dataset_vals.append(r)
                                else:
                                    cross_dataset_vals.append(r)
                        for comparison, arr in (
                            ("all_dataset_pairs", vals),
                            ("same_dataset", same_dataset_vals),
                            ("cross_dataset", cross_dataset_vals),
                        ):
                            if not arr:
                                continue
                            rows.append(
                                {
                                    "pipeline": pipeline,
                                    "feature_a": feature_a,
                                    "label_a": label_a,
                                    "feature_b": feature_b,
                                    "label_b": label_b,
                                    "comparison": comparison,
                                    "n_corr": len(arr),
                                    "mean_corr": float(np.mean(arr)),
                                    "median_corr": float(np.median(arr)),
                                    "min_corr": float(np.min(arr)),
                                    "max_corr": float(np.max(arr)),
                                }
                            )
    return rows


def baseline_records(
    by_dataset: dict[str, dict[int, np.ndarray]],
    selected_vectors: dict[tuple[str, str, str, str], np.ndarray],
    datasets: Sequence[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    # Baseline: all P7 roi_mean mode pairs across dataset pairs.
    all_pair_corrs = []
    same_rank_corrs = []
    for ds_a, ds_b in combinations([d.lower() for d in datasets], 2):
        modes_a = by_dataset.get(ds_a, {})
        modes_b = by_dataset.get(ds_b, {})
        shared = sorted(set(modes_a) & set(modes_b))
        for mode in shared:
            r = pearson(modes_a[mode], modes_b[mode])
            if math.isfinite(r):
                same_rank_corrs.append(r)
        for ma in shared:
            va = modes_a[ma]
            for mb in shared:
                r = pearson(va, modes_b[mb])
                if math.isfinite(r):
                    all_pair_corrs.append(r)

    def add_dist(name: str, vals: Sequence[float]) -> None:
        if not vals:
            return
        rows.append(
            {
                "distribution": name,
                "n_corr": len(vals),
                "mean_corr": float(np.mean(vals)),
                "median_corr": float(np.median(vals)),
                "q25_corr": float(np.quantile(vals, 0.25)),
                "q75_corr": float(np.quantile(vals, 0.75)),
                "min_corr": float(np.min(vals)),
                "max_corr": float(np.max(vals)),
            }
        )

    add_dist("baseline_all_mode_pairs_cross_dataset", all_pair_corrs)
    add_dist("baseline_same_sorted_position_cross_dataset", same_rank_corrs)
    for pipeline in PIPELINES:
        for feature in FEATURES:
            vals = []
            for ds_a, ds_b in combinations([d.lower() for d in datasets], 2):
                va = selected_vectors.get((pipeline, ds_a, feature, "ripple_gamma_no_theta"))
                vb = selected_vectors.get((pipeline, ds_b, feature, "ripple_gamma_no_theta"))
                if va is None or vb is None:
                    continue
                r = pearson(va, vb)
                if math.isfinite(r):
                    vals.append(r)
            add_dist(f"{pipeline}_{feature}_rg_no_theta_selected_cross_dataset", vals)
    for pipeline in PIPELINES:
        vals = []
        for ds in [d.lower() for d in datasets]:
            ve = selected_vectors.get((pipeline, ds, "efun", "ripple_gamma_no_theta"))
            vd = selected_vectors.get((pipeline, ds, "deconv_efun", "ripple_gamma_no_theta"))
            if ve is None or vd is None:
                continue
            r = pearson(ve, vd)
            if math.isfinite(r):
                vals.append(r)
        add_dist(f"{pipeline}_efun_vs_deconv_rg_no_theta_same_dataset", vals)
    return rows


def write_roi_profile_matrix(profile_rows: Sequence[dict[str, object]], path: Path) -> None:
    fields = [
        "pipeline",
        "dataset",
        "feature_family",
        "strict_label",
        "roi_index",
        "roi_label",
        "roi_value_weighted",
        "weight_sum",
        "n_hit_rows",
    ]
    write_csv(path, profile_rows, fields)


def row_normalize(mat: np.ndarray) -> np.ndarray:
    out = np.asarray(mat, dtype=float).copy()
    for i in range(out.shape[0]):
        row = out[i]
        finite = np.isfinite(row)
        if not finite.any():
            continue
        lo = float(np.nanmin(row))
        hi = float(np.nanmax(row))
        if hi > lo:
            out[i, finite] = (row[finite] - lo) / (hi - lo)
        else:
            out[i, finite] = 0.0
    return out


def plot_all_roi_by_dataset(
    vectors: dict[tuple[str, str, str, str], np.ndarray],
    roi_order: Sequence[str],
    datasets: Sequence[str],
    figure_dir: Path,
    top_n: int,
) -> list[Path]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for pipeline in PIPELINES:
        for feature in FEATURES:
            label = "ripple_gamma_no_theta"
            rows = []
            present = []
            for ds in datasets:
                vec = vectors.get((pipeline, ds.lower(), feature, label))
                if vec is None:
                    vec = np.full(len(roi_order), np.nan)
                else:
                    present.append(ds)
                rows.append(vec)
            mat = np.vstack(rows)
            shown = row_normalize(mat)
            fig_w = max(18, len(roi_order) * 0.32)
            fig_h = max(5, len(datasets) * 0.36)
            fig, ax = plt.subplots(figsize=(fig_w, fig_h), constrained_layout=True)
            im = ax.imshow(shown, aspect="auto", cmap="viridis", vmin=0, vmax=1)
            ax.set_title(
                f"{pipeline} roi_mean RG-no-theta selected ROI profiles | {feature} | top{top_n}\n"
                "each dataset row normalized independently; ROI order preserved",
                weight="bold",
            )
            ax.set_yticks(range(len(datasets)))
            ax.set_yticklabels(datasets)
            ax.set_xticks(range(len(roi_order)))
            ax.set_xticklabels(roi_order, rotation=75, ha="right", fontsize=7)
            ax.set_xlabel("ROI")
            ax.set_ylabel("dataset")
            fig.colorbar(im, ax=ax, shrink=0.75, label="within-dataset normalized ROI value")
            path = figure_dir / f"01_all_roi_by_dataset__{pipeline.lower()}__{safe_name(feature)}__rg_no_theta__top{top_n}.png"
            fig.savefig(path, dpi=180)
            plt.close(fig)
            paths.append(path)
    return paths


def plot_confusion(pair_rows: Sequence[dict[str, object]], figure_dir: Path, top_n: int) -> list[Path]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    row_df = pd.DataFrame(pair_rows)
    if row_df.empty:
        return paths
    for pipeline in PIPELINES:
        for feature in FEATURES:
            sub = row_df[
                (row_df["pipeline"] == pipeline)
                & (row_df["feature_a"] == feature)
                & (row_df["feature_b"] == feature)
                & (row_df["comparison"] == "cross_dataset")
            ]
            if sub.empty:
                continue
            mat = np.full((len(LABELS), len(LABELS)), np.nan)
            for _, row in sub.iterrows():
                try:
                    i = LABELS.index(row["label_a"])
                    j = LABELS.index(row["label_b"])
                except ValueError:
                    continue
                mat[i, j] = float(row["mean_corr"])
            fig, ax = plt.subplots(figsize=(8.5, 7), constrained_layout=True)
            im = ax.imshow(mat, cmap="coolwarm", vmin=-0.2, vmax=1.0)
            ax.set_title(f"{pipeline} ROI cross-dataset label confusion | {feature} | top{top_n}", weight="bold")
            ax.set_xticks(range(len(LABELS)))
            ax.set_xticklabels(LABELS, rotation=35, ha="right")
            ax.set_yticks(range(len(LABELS)))
            ax.set_yticklabels(LABELS)
            for i in range(len(LABELS)):
                for j in range(len(LABELS)):
                    if math.isfinite(mat[i, j]):
                        ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=9)
            fig.colorbar(im, ax=ax, shrink=0.8, label="mean ROI profile correlation")
            path = figure_dir / f"02_label_confusion_cross_dataset__{pipeline.lower()}__{safe_name(feature)}__top{top_n}.png"
            fig.savefig(path, dpi=180)
            plt.close(fig)
            paths.append(path)

    for pipeline in PIPELINES:
        sub = row_df[
            (row_df["pipeline"] == pipeline)
            & (row_df["label_a"] == "ripple_gamma_no_theta")
            & (row_df["label_b"] == "ripple_gamma_no_theta")
            & (row_df["comparison"] == "all_dataset_pairs")
        ]
        if sub.empty:
            continue
        mat = np.full((len(FEATURES), len(FEATURES)), np.nan)
        for _, row in sub.iterrows():
            i = FEATURES.index(row["feature_a"])
            j = FEATURES.index(row["feature_b"])
            mat[i, j] = float(row["mean_corr"])
        fig, ax = plt.subplots(figsize=(5.5, 4.8), constrained_layout=True)
        im = ax.imshow(mat, cmap="coolwarm", vmin=-0.2, vmax=1.0)
        ax.set_title(f"{pipeline} RG-no-theta efun/deconv ROI similarity | top{top_n}", weight="bold")
        ax.set_xticks(range(len(FEATURES)))
        ax.set_xticklabels(FEATURES, rotation=25, ha="right")
        ax.set_yticks(range(len(FEATURES)))
        ax.set_yticklabels(FEATURES)
        for i in range(len(FEATURES)):
            for j in range(len(FEATURES)):
                if math.isfinite(mat[i, j]):
                    ax.text(j, i, f"{mat[i, j]:.2f}", ha="center", va="center", fontsize=10)
        fig.colorbar(im, ax=ax, shrink=0.8, label="mean ROI profile correlation")
        path = figure_dir / f"03_rg_no_theta_efun_deconv_similarity__{pipeline.lower()}__top{top_n}.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_baseline(baseline_rows: Sequence[dict[str, object]], figure_dir: Path, top_n: int) -> Path | None:
    df = pd.DataFrame(baseline_rows)
    if df.empty:
        return None
    wanted = [
        "baseline_all_mode_pairs_cross_dataset",
        "baseline_same_sorted_position_cross_dataset",
        "P8_efun_rg_no_theta_selected_cross_dataset",
        "P8_deconv_efun_rg_no_theta_selected_cross_dataset",
        "P10_efun_rg_no_theta_selected_cross_dataset",
        "P10_deconv_efun_rg_no_theta_selected_cross_dataset",
        "P8_efun_vs_deconv_rg_no_theta_same_dataset",
        "P10_efun_vs_deconv_rg_no_theta_same_dataset",
    ]
    df = df[df["distribution"].isin(wanted)].copy()
    if df.empty:
        return None
    df["distribution"] = pd.Categorical(df["distribution"], categories=wanted, ordered=True)
    df = df.sort_values("distribution")
    fig, ax = plt.subplots(figsize=(15, 5.5), constrained_layout=True)
    x = np.arange(len(df))
    med = df["median_corr"].astype(float).to_numpy()
    q25 = df["q25_corr"].astype(float).to_numpy()
    q75 = df["q75_corr"].astype(float).to_numpy()
    ax.bar(x, med, color="#6788c7")
    ax.errorbar(x, med, yerr=[med - q25, q75 - med], fmt="none", ecolor="black", capsize=4, lw=1)
    for i, row in enumerate(df.to_dict("records")):
        ax.text(i, med[i] + 0.025, f"n={int(row['n_corr'])}", ha="center", va="bottom", fontsize=8)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylim(-0.3, 1.05)
    ax.set_ylabel("ROI profile correlation median (IQR)")
    ax.set_xticks(x)
    ax.set_xticklabels([str(v).replace("_", "\n") for v in df["distribution"]], fontsize=8)
    ax.set_title(f"Selected RG-no-theta ROI consistency vs all-mode BOLD baseline | top{top_n}", weight="bold")
    path = figure_dir / f"04_selected_vs_allmode_baseline__top{top_n}.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main() -> None:
    args = parse_args()
    datasets = {d.lower() for d in args.datasets}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    roi_order, dataset_roi_order, p7_profiles, by_dataset = read_p7_roi_profiles(args.p7_roi_csv, datasets)
    all_summary_rows: list[dict[str, object]] = []
    all_figure_paths: list[Path] = []
    main_profile_rows: list[dict[str, object]] = []

    for top_n in args.top_ns:
        print(f"[ROI test] Loading hits top{top_n}")
        hits = read_hits(args.hits_csv, datasets, top_n)
        print(f"[ROI test] top{top_n} filtered rows: {len(hits)}")
        vectors, profile_rows, selected_rows = aggregate_vectors(
            hits,
            roi_order,
            dataset_roi_order,
            p7_profiles,
            args.processed_root,
        )
        pair_rows = vector_pair_records(vectors, args.datasets)
        baseline_rows = baseline_records(by_dataset, vectors, args.datasets)

        suffix = f"top{top_n}"
        write_roi_profile_matrix(profile_rows, args.output_dir / f"roi_mean_rg_no_theta_selected_roi_profiles_long__{suffix}.csv")
        write_csv(args.output_dir / f"roi_mean_rg_no_theta_selected_units_and_missing__{suffix}.csv", selected_rows)
        write_csv(args.output_dir / f"roi_pairwise_label_feature_confusion__{suffix}.csv", pair_rows)
        write_csv(args.output_dir / f"roi_selected_vs_allmode_baseline_summary__{suffix}.csv", baseline_rows)

        if top_n == args.main_top_n:
            main_profile_rows = profile_rows
            all_figure_paths.extend(plot_all_roi_by_dataset(vectors, roi_order, args.datasets, args.figure_dir, top_n))
            all_figure_paths.extend(plot_confusion(pair_rows, args.figure_dir, top_n))
            baseline_path = plot_baseline(baseline_rows, args.figure_dir, top_n)
            if baseline_path is not None:
                all_figure_paths.append(baseline_path)

        for row in baseline_rows:
            row = dict(row)
            row["top_n"] = top_n
            all_summary_rows.append(row)
        for row in pair_rows:
            if (
                row.get("comparison") in {"cross_dataset", "same_dataset"}
                and row.get("label_a") == "ripple_gamma_no_theta"
                and row.get("label_b") == "ripple_gamma_no_theta"
            ):
                row = dict(row)
                row["top_n"] = top_n
                all_summary_rows.append(row)

    write_csv(args.output_dir / "roi_test_compact_summary_all_topN.csv", all_summary_rows)
    with (args.output_dir / "README_roi_mean_rg_no_theta_formal_roi_tests.md").open("w", encoding="utf-8") as handle:
        handle.write("# roi_mean RG-no-theta formal ROI tests\n\n")
        handle.write(f"- Hits CSV: `{args.hits_csv}`\n")
        handle.write(f"- P7 ROI CSV: `{args.p7_roi_csv}`\n")
        handle.write(f"- Datasets: `{', '.join(args.datasets)}`\n")
        handle.write(f"- TopN values: `{', '.join(map(str, args.top_ns))}`\n")
        handle.write("- Main figures use topN = `%d`.\n" % args.main_top_n)
        handle.write("- All-ROI dataset panels preserve ROI order from P7 roi_mean export and normalize each dataset row independently.\n")
        handle.write("- P8 selected units map to P7 sorted BOLD modes; P10 selected units map to P9 BOLD dimred components reconstructed from P9 MAT files.\n\n")
        handle.write("## Main Figures\n\n")
        for path in all_figure_paths:
            handle.write(f"- `{path}`\n")
        handle.write("\n## Main Output CSVs\n\n")
        for name in [
            "roi_test_compact_summary_all_topN.csv",
            f"roi_mean_rg_no_theta_selected_roi_profiles_long__top{args.main_top_n}.csv",
            f"roi_pairwise_label_feature_confusion__top{args.main_top_n}.csv",
            f"roi_selected_vs_allmode_baseline_summary__top{args.main_top_n}.csv",
            f"roi_mean_rg_no_theta_selected_units_and_missing__top{args.main_top_n}.csv",
        ]:
            handle.write(f"- `{args.output_dir / name}`\n")
    print(f"[ROI test] Wrote {len(main_profile_rows)} main topN profile rows")
    print(f"[ROI test] Output: {args.output_dir}")
    print(f"[ROI test] Figures: {args.figure_dir}")


if __name__ == "__main__":
    main()
