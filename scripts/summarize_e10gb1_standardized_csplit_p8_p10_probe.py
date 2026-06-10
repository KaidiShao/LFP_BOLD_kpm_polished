"""Summarize the e10gb1 standardized complex-split P8/P10 probe.

This is intentionally narrow: it reads the custom xcorr tags generated for the
standardized complex-split branch and answers whether event/raw/dimred density
sources win, and what two-subprocess labels the winning dimred components have.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "e10gb1_standardized_csplit_p8_p10_probe"
DEFAULT_FIGURE_DIR = (
    DEFAULT_PROCESSED_ROOT
    / "summary_figures"
    / "pipeline11_current_analysis_summary"
    / "p8_p10_e10gb1_standardized_csplit_probe"
)
DEFAULT_LABEL_FILE = (
    REPO_ROOT
    / "results"
    / "e10gb1_p5_two_subprocess_candidates"
    / "component_subprocess_flags.csv"
)
DEFAULT_RAW_METADATA_FILE = (
    REPO_ROOT
    / "results"
    / "p8_p10_parameter_selection_current_rmsenv_adaptive_v2"
    / "p5_raw_density_mode_metadata.csv"
)

P8_TAG = "xcorr_csplit_standardize_rmsenv_adaptive"
P10_TAG = "dimred_xcorr_csplit_standardize_rmsenv_adaptive"
METHODS = ("svd", "nmf", "mds", "umap")
KS_BLP = tuple(range(3, 9))
SOURCE_KIND_ORDER = ("event", "raw", "dimred")
FAMILY_ORDER = ("efun", "deconv_efun")
RAW_TIMESCALE_FIELD = "timescale_sec_preferred"
RAW_FREQUENCY_FIELD = "frequency_hz_preferred"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def safe_float(value: object) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return math.nan
    return out


def safe_int(value: object) -> int | None:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def feature_family(name: str) -> str:
    text = name.lower()
    if "deconv" in text:
        return "deconv_efun"
    return "efun"


def parse_density_name(name: str) -> dict[str, object]:
    text = name.lower()
    if text.startswith("blp_evt"):
        return {"source_kind": "event", "blp_method": "", "blp_k": ""}
    if text.startswith("raw_csplit"):
        return {"source_kind": "raw", "blp_method": "", "blp_k": ""}
    m = re.search(r"dim_csplit_([a-z]+)(\d+)_q", text)
    if m:
        return {
            "source_kind": "dimred",
            "blp_method": m.group(1),
            "blp_k": int(m.group(2)),
        }
    return {"source_kind": "other", "blp_method": "", "blp_k": ""}


def load_subprocess_labels(label_file: Path) -> dict[tuple[str, str, int, int], dict[str, str]]:
    labels: dict[tuple[str, str, int, int], dict[str, str]] = {}
    if not label_file.is_file():
        return labels
    for row in read_csv(label_file):
        if row.get("branch") != "standardize" or row.get("condition") != "csplit":
            continue
        method = row.get("method", "").lower()
        k = safe_int(row.get("k"))
        comp = safe_int(row.get("component_idx"))
        if not method or k is None or comp is None:
            continue
        theta_strict = str(row.get("theta_strict", "")).lower() == "true"
        theta_relaxed = str(row.get("theta_relaxed", "")).lower() == "true"
        ripple_gamma = str(row.get("ripple_gamma_no_pure_theta", "")).lower() == "true"
        ripple_no_theta = str(row.get("ripple_no_pure_theta", "")).lower() == "true"
        if theta_strict and ripple_gamma:
            label = "theta+ripple_gamma"
        elif theta_strict:
            label = "theta_strict"
        elif theta_relaxed and ripple_gamma:
            label = "theta_relaxed+ripple_gamma"
        elif theta_relaxed:
            label = "theta_relaxed"
        elif ripple_gamma:
            label = "ripple_gamma_no_pure_theta"
        elif ripple_no_theta:
            label = "ripple_no_pure_theta"
        else:
            label = "other"
        labels[(method, "csplit", k, comp)] = {
            "two_subprocess_label": label,
            "theta_strict": str(theta_strict),
            "theta_relaxed": str(theta_relaxed),
            "ripple_gamma_no_pure_theta": str(ripple_gamma),
            "ripple_no_pure_theta": str(ripple_no_theta),
            "theta_active_events": row.get("theta_active_events", ""),
            "ripple_active_events": row.get("ripple_active_events", ""),
            "gamma_active_events": row.get("gamma_active_events", ""),
        }
    return labels


def norm_path(text: object) -> str:
    return str(text or "").replace("/", "\\").lower()


def load_raw_mode_metadata(metadata_file: Path) -> dict[tuple[str, int], dict[str, str]]:
    out: dict[tuple[str, int], dict[str, str]] = {}
    if not metadata_file.is_file():
        return out
    for row in read_csv(metadata_file):
        density_index = safe_int(row.get("density_index"))
        if density_index is None:
            continue
        out[(norm_path(row.get("density_file")), density_index)] = row
    return out


def add_common_fields(
    row: dict[str, str],
    analysis: str,
    observable: str,
    context_id: str,
    label_map: dict[tuple[str, str, int, int], dict[str, str]],
    p10_bold_feature: str = "",
    p10_bold_method_tag: str = "",
) -> dict[str, object]:
    density_name = row.get("density_name", "")
    parsed = parse_density_name(density_name)
    family = feature_family(row.get("bold_feature", "") or p10_bold_feature)
    density_index = safe_int(row.get("density_index"))
    blp_method = str(parsed["blp_method"])
    blp_k = parsed["blp_k"]
    label_info: dict[str, str] = {}
    if parsed["source_kind"] == "dimred" and density_index is not None and isinstance(blp_k, int):
        label_info = label_map.get((blp_method, "csplit", blp_k, density_index), {})
    out: dict[str, object] = {
        "analysis": analysis,
        "observable": observable,
        "context_id": context_id,
        "feature_family": family,
        "bold_feature": row.get("bold_feature", p10_bold_feature),
        "p10_bold_feature": p10_bold_feature,
        "p10_bold_method_tag": p10_bold_method_tag,
        "density_name": density_name,
        "density_type": row.get("density_type", ""),
        "density_file": row.get("density_file", ""),
        "source_kind": parsed["source_kind"],
        "blp_method": blp_method,
        "blp_k": blp_k,
        "density_index": density_index if density_index is not None else "",
        "density_label": row.get("density_label", ""),
        "bold_mode_index": row.get("bold_mode_index", ""),
        "peak_abs_corr": safe_float(row.get("peak_abs_corr")),
        "peak_corr": safe_float(row.get("peak_corr")),
        "peak_lag_sec": safe_float(row.get("peak_lag_sec")),
    }
    out.update(
        {
            "two_subprocess_label": label_info.get("two_subprocess_label", ""),
            "theta_strict": label_info.get("theta_strict", ""),
            "theta_relaxed": label_info.get("theta_relaxed", ""),
            "ripple_gamma_no_pure_theta": label_info.get("ripple_gamma_no_pure_theta", ""),
            "ripple_no_pure_theta": label_info.get("ripple_no_pure_theta", ""),
            "theta_active_events": label_info.get("theta_active_events", ""),
            "ripple_active_events": label_info.get("ripple_active_events", ""),
            "gamma_active_events": label_info.get("gamma_active_events", ""),
        }
    )
    return out


def enrich_raw_metadata(
    rows: list[dict[str, object]],
    raw_metadata: dict[tuple[str, int], dict[str, str]],
) -> None:
    for row in rows:
        if row.get("source_kind") != "raw":
            continue
        density_index = safe_int(row.get("density_index"))
        if density_index is None:
            continue
        meta = raw_metadata.get((norm_path(row.get("density_file")), density_index))
        if not meta:
            continue
        for key in (
            "raw_efun_index",
            "source_mode_index",
            "timescale_sec",
            "frequency_hz",
            "timescale_sec_preferred",
            "frequency_hz_preferred",
            "timescale_sec_discrete_log",
            "frequency_hz_discrete_angle",
            "timescale_sec_bilinear",
            "frequency_hz_bilinear",
            "timescale_source_preferred",
            "frequency_source_preferred",
            "envelope_window_sec",
            "envelope_timescale_sec",
            "envelope_window_status",
        ):
            row[key] = meta.get(key, "")


def collect_p8(processed_root: Path, label_map: dict[tuple[str, str, int, int], dict[str, str]]) -> list[dict[str, object]]:
    root = processed_root / "e10gb1" / "pipeline8_xcorr"
    rows: list[dict[str, object]] = []
    for file in sorted(root.glob(f"*/feature/*/{P8_TAG}_top__*.csv")):
        observable = file.parents[2].name
        for row in read_csv(file):
            family = feature_family(row.get("bold_feature", ""))
            context_id = f"P8::{observable}::{family}"
            rows.append(add_common_fields(row, "P8", observable, context_id, label_map))
    return rows


def collect_p10(processed_root: Path, label_map: dict[tuple[str, str, int, int], dict[str, str]]) -> list[dict[str, object]]:
    root = processed_root / "e10gb1" / "pipeline10_dimred_xcorr"
    rows: list[dict[str, object]] = []
    for file in sorted(root.glob(f"*/feature/*/{P10_TAG}_top__*.csv")):
        run_dir = file.parents[2].name
        parts = run_dir.split("__")
        observable = parts[0] if len(parts) >= 1 else run_dir
        p10_feature = parts[1] if len(parts) >= 2 else ""
        p10_method = parts[2] if len(parts) >= 3 else ""
        family = feature_family(p10_feature)
        context_id = f"P10::{observable}::{p10_feature}::{p10_method}"
        for row in read_csv(file):
            rows.append(
                add_common_fields(
                    row,
                    "P10",
                    observable,
                    context_id,
                    label_map,
                    p10_bold_feature=p10_feature,
                    p10_bold_method_tag=p10_method,
                )
            )
    return rows


def best_by(rows: list[dict[str, object]], keys: tuple[str, ...]) -> list[dict[str, object]]:
    best: dict[tuple[object, ...], dict[str, object]] = {}
    for row in rows:
        key = tuple(row.get(k, "") for k in keys)
        value = safe_float(row.get("peak_abs_corr"))
        if not math.isfinite(value):
            continue
        if key not in best or value > safe_float(best[key].get("peak_abs_corr")):
            best[key] = row
    return list(best.values())


def top_rows(rows: list[dict[str, object]], n: int) -> list[dict[str, object]]:
    return sorted(rows, key=lambda r: safe_float(r.get("peak_abs_corr")), reverse=True)[:n]


def mean(values: list[float]) -> float:
    clean = [v for v in values if math.isfinite(v)]
    if not clean:
        return math.nan
    return sum(clean) / len(clean)


def plot_source_kind_winners(best_context_rows: list[dict[str, object]], out: Path, title: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in best_context_rows:
        counts[str(row["feature_family"])][str(row["source_kind"])] += 1
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    x = list(range(len(FAMILY_ORDER)))
    width = 0.22
    for i, kind in enumerate(SOURCE_KIND_ORDER):
        vals = [counts[fam][kind] for fam in FAMILY_ORDER]
        ax.bar([xx + (i - 1) * width for xx in x], vals, width=width, label=kind)
    ax.set_xticks(x, FAMILY_ORDER)
    ax.set_ylabel("winner contexts")
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_source_kind_scores(best_kind_rows: list[dict[str, object]], out: Path, title: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    labels: list[str] = []
    data: list[list[float]] = []
    for fam in FAMILY_ORDER:
        for kind in SOURCE_KIND_ORDER:
            labels.append(f"{fam}\n{kind}")
            vals = [
                safe_float(r["peak_abs_corr"])
                for r in best_kind_rows
                if r.get("feature_family") == fam and r.get("source_kind") == kind
            ]
            data.append([v for v in vals if math.isfinite(v)])
    ax.boxplot(data, labels=labels, showfliers=False)
    for i, vals in enumerate(data, start=1):
        ax.scatter([i] * len(vals), vals, s=18, alpha=0.55)
    ax.set_ylabel("best |xcorr| per context/source kind")
    ax.set_title(title)
    ax.tick_params(axis="x", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_blp_method_k_heatmap(best_dimred_rows: list[dict[str, object]], out: Path, title: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = [f"{m}_k{k:02d}" for m in METHODS for k in KS_BLP]
    matrix: list[list[float]] = []
    for fam in FAMILY_ORDER:
        row_vals: list[float] = []
        for m in METHODS:
            for k in KS_BLP:
                vals = [
                    safe_float(r["peak_abs_corr"])
                    for r in best_dimred_rows
                    if r.get("feature_family") == fam
                    and r.get("blp_method") == m
                    and safe_int(r.get("blp_k")) == k
                ]
                row_vals.append(mean(vals))
        matrix.append(row_vals)
    fig, ax = plt.subplots(figsize=(14, 3.2))
    im = ax.imshow(matrix, aspect="auto", cmap="viridis")
    ax.set_yticks(range(len(FAMILY_ORDER)), FAMILY_ORDER)
    ax.set_xticks(range(len(cols)), cols, rotation=90, fontsize=7)
    ax.set_title(title)
    cbar = fig.colorbar(im, ax=ax, shrink=0.82)
    cbar.set_label("mean best |xcorr|")
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_dimred_label_counts(rows: list[dict[str, object]], out: Path, title: str, top_n: int) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    counters: dict[str, Counter[str]] = defaultdict(Counter)
    for analysis in ("P8", "P10"):
        subset = [r for r in rows if r.get("analysis") == analysis and r.get("source_kind") == "dimred"]
        for row in top_rows(subset, top_n):
            label = str(row.get("two_subprocess_label") or "unlabelled")
            counters[analysis][label] += 1
    labels = sorted(set(counters["P8"]) | set(counters["P10"]))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for ax, analysis in zip(axes, ("P8", "P10")):
        vals = [counters[analysis][lab] for lab in labels]
        ax.barh(labels, vals)
        ax.set_title(analysis)
        ax.set_xlabel(f"count in top{top_n} dimred hits")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_dimred_label_counts_by_family(rows: list[dict[str, object]], out: Path, title: str, top_n: int) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    labels = ("theta_strict", "ripple_gamma_no_pure_theta", "ripple_no_pure_theta", "theta_relaxed", "other", "unlabelled")
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 7.2), sharex=True)
    for ax, analysis in zip(axes[:, 0], ("P8", "P10")):
        ax.set_ylabel(analysis)
    for j, fam in enumerate(FAMILY_ORDER):
        axes[0, j].set_title(fam)
    for i, analysis in enumerate(("P8", "P10")):
        for j, fam in enumerate(FAMILY_ORDER):
            subset = [
                r
                for r in rows
                if r.get("analysis") == analysis
                and r.get("feature_family") == fam
                and r.get("source_kind") == "dimred"
            ]
            counts = Counter(str(r.get("two_subprocess_label") or "unlabelled") for r in top_rows(subset, top_n))
            vals = [counts.get(label, 0) for label in labels]
            ax = axes[i, j]
            ax.barh(labels, vals)
            ax.set_xlabel(f"count in top{top_n}")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def finite_values(rows: list[dict[str, object]], field: str) -> list[float]:
    values = [safe_float(r.get(field)) for r in rows]
    return [v for v in values if math.isfinite(v)]


def plot_ecdf(ax: plt.Axes, values: list[float], label: str) -> None:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return
    y = [(i + 1) / len(clean) for i in range(len(clean))]
    ax.step(clean, y, where="post", label=f"{label} (n={len(clean)})")


def plot_raw_index_ecdf_by_topn(
    rows: list[dict[str, object]],
    out: Path,
    top_ns: tuple[int, ...] = (10, 20, 50),
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, len(top_ns), figsize=(5.2 * len(top_ns), 7.2), sharey=True)
    for i, analysis in enumerate(("P8", "P10")):
        for j, top_n in enumerate(top_ns):
            ax = axes[i, j]
            for fam in FAMILY_ORDER:
                subset = [
                    r
                    for r in rows
                    if r.get("analysis") == analysis
                    and r.get("feature_family") == fam
                    and r.get("source_kind") == "raw"
                ]
                subset = top_rows(subset, top_n)
                values = finite_values(subset, "raw_efun_index")
                plot_ecdf(ax, values, fam)
            ax.set_title(f"{analysis} top{top_n}")
            ax.set_xlabel("raw BLP efun index")
            ax.grid(alpha=0.2)
            if j == 0:
                ax.set_ylabel("ECDF")
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Raw BLP efun index distribution: BOLD efun vs deconv_efun")
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def plot_raw_timescale_ecdf_by_topn(
    rows: list[dict[str, object]],
    out: Path,
    top_ns: tuple[int, ...] = (10, 20, 50),
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, len(top_ns), figsize=(5.2 * len(top_ns), 7.2), sharey=True)
    for i, analysis in enumerate(("P8", "P10")):
        for j, top_n in enumerate(top_ns):
            ax = axes[i, j]
            for fam in FAMILY_ORDER:
                subset = [
                    r
                    for r in rows
                    if r.get("analysis") == analysis
                    and r.get("feature_family") == fam
                    and r.get("source_kind") == "raw"
                ]
                subset = top_rows(subset, top_n)
                values = [v for v in finite_values(subset, RAW_TIMESCALE_FIELD) if v > 0]
                values = [math.log10(v) for v in values]
                plot_ecdf(ax, values, fam)
            ax.set_title(f"{analysis} top{top_n}")
            ax.set_xlabel("log10(raw BLP efun discrete-log timescale sec)")
            ax.grid(alpha=0.2)
            if j == 0:
                ax.set_ylabel("ECDF")
            ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Raw BLP efun preferred timescale distribution: BOLD efun vs deconv_efun")
    fig.tight_layout()
    fig.savefig(out, dpi=180)
    plt.close(fig)


def summarize_raw_distribution(rows: list[dict[str, object]], top_ns: tuple[int, ...]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for analysis in ("P8", "P10"):
        for fam in FAMILY_ORDER:
            subset_all = [
                r
                for r in rows
                if r.get("analysis") == analysis
                and r.get("feature_family") == fam
                and r.get("source_kind") == "raw"
            ]
            for top_n in top_ns:
                subset = top_rows(subset_all, top_n)
                idx_values = finite_values(subset, "raw_efun_index")
                tau_values = [v for v in finite_values(subset, RAW_TIMESCALE_FIELD) if v > 0]
                freq_values = [v for v in finite_values(subset, RAW_FREQUENCY_FIELD) if v >= 0]
                out.append(
                    {
                        "analysis": analysis,
                        "feature_family": fam,
                        "top_n": top_n,
                        "n_raw_hits": len(subset),
                        "n_with_timescale": len(tau_values),
                        "raw_efun_index_median": median(idx_values),
                        "raw_efun_index_min": min(idx_values) if idx_values else "",
                        "raw_efun_index_max": max(idx_values) if idx_values else "",
                        "timescale_field": RAW_TIMESCALE_FIELD,
                        "timescale_sec_median": median(tau_values),
                        "timescale_sec_min": min(tau_values) if tau_values else "",
                        "timescale_sec_max": max(tau_values) if tau_values else "",
                        "frequency_field": RAW_FREQUENCY_FIELD,
                        "frequency_hz_median": median(freq_values),
                    }
                )
    return out


def median(values: list[float]) -> float | str:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return ""
    n = len(clean)
    mid = n // 2
    if n % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def write_summary(
    output_dir: Path,
    p8_rows: list[dict[str, object]],
    p10_rows: list[dict[str, object]],
    best_context: list[dict[str, object]],
    top_dimred: list[dict[str, object]],
) -> None:
    lines: list[str] = [
        "# E10gb1 standardized complex-split P8/P10 probe",
        "",
        f"- P8 top rows read: `{len(p8_rows)}`",
        f"- P10 top rows read: `{len(p10_rows)}`",
        "- Branch: `complex_split_projected_vlambda_standardize_rmsenv_adaptive`",
        "- P8 density sources: event + raw csplit + dimred csplit SVD/NMF/MDS/UMAP k03:k08",
        "- P10 BOLD-side P9 grid: SVD/NMF/MDS/UMAP k05:k08 for available e10gb1 current P7/P9 outputs",
        "",
        "## Winner source kind by context",
    ]
    for analysis in ("P8", "P10"):
        subset = [r for r in best_context if r.get("analysis") == analysis]
        counts = Counter(str(r.get("source_kind")) for r in subset)
        lines.append(f"- {analysis}: " + ", ".join(f"{k}={counts.get(k, 0)}" for k in SOURCE_KIND_ORDER))
    lines.extend(["", "## Winner source kind by family"])
    for analysis in ("P8", "P10"):
        for fam in FAMILY_ORDER:
            subset = [
                r
                for r in best_context
                if r.get("analysis") == analysis and r.get("feature_family") == fam
            ]
            counts = Counter(str(r.get("source_kind")) for r in subset)
            lines.append(
                f"- {analysis} {fam}: "
                + ", ".join(f"{k}={counts.get(k, 0)}" for k in SOURCE_KIND_ORDER)
            )
    lines.extend(["", "## Top dimred hits"])
    for row in top_dimred[:15]:
        lines.append(
            "- "
            f"{row['analysis']} {row['observable']} {row['feature_family']} "
            f"{row.get('p10_bold_method_tag','')} | {row['density_name']} "
            f"comp {row['density_index']} | |r|={safe_float(row['peak_abs_corr']):.3f} "
            f"lag={safe_float(row['peak_lag_sec']):.1f}s | "
            f"{row.get('two_subprocess_label') or 'unlabelled'}"
        )
    lines.extend(
        [
            "",
            "## Notes",
            "",
            "This is an e10gb1-only probe.  It can test whether standardized",
            "complex-split produces interpretable P8/P10 hits, but it is not a",
            "cross-session conclusion yet.",
        ]
    )
    (output_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--label-file", type=Path, default=DEFAULT_LABEL_FILE)
    parser.add_argument("--raw-metadata-file", type=Path, default=DEFAULT_RAW_METADATA_FILE)
    parser.add_argument("--top-label-n", type=int, default=50)
    args = parser.parse_args()

    label_map = load_subprocess_labels(args.label_file)
    raw_metadata = load_raw_mode_metadata(args.raw_metadata_file)
    p8_rows = collect_p8(args.processed_root, label_map)
    p10_rows = collect_p10(args.processed_root, label_map)
    all_rows = p8_rows + p10_rows
    enrich_raw_metadata(all_rows, raw_metadata)

    best_context = best_by(all_rows, ("analysis", "context_id"))
    best_kind = best_by(all_rows, ("analysis", "context_id", "source_kind"))
    best_dimred_method_k = best_by(
        [r for r in all_rows if r.get("source_kind") == "dimred"],
        ("analysis", "context_id", "blp_method", "blp_k"),
    )
    top_dimred = top_rows([r for r in all_rows if r.get("source_kind") == "dimred"], 200)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figure_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "all_top_rows.csv", all_rows)
    write_csv(args.output_dir / "best_context_rows.csv", best_context)
    write_csv(args.output_dir / "best_context_source_kind_rows.csv", best_kind)
    write_csv(args.output_dir / "best_dimred_method_k_rows.csv", best_dimred_method_k)
    write_csv(args.output_dir / "top_dimred_rows.csv", top_dimred)
    raw_distribution = summarize_raw_distribution(all_rows, (10, 20, 50))
    write_csv(args.output_dir / "raw_blp_efun_distribution_summary_v3_discrete_log.csv", raw_distribution)

    for analysis in ("P8", "P10"):
        plot_source_kind_winners(
            [r for r in best_context if r.get("analysis") == analysis],
            args.figure_dir / f"{analysis.lower()}_source_kind_winners.png",
            f"{analysis} e10gb1 standardized csplit: winner source kind",
        )
        plot_source_kind_scores(
            [r for r in best_kind if r.get("analysis") == analysis],
            args.figure_dir / f"{analysis.lower()}_source_kind_score_distribution.png",
            f"{analysis} e10gb1 standardized csplit: event/raw/dimred score",
        )
        plot_blp_method_k_heatmap(
            [r for r in best_dimred_method_k if r.get("analysis") == analysis],
            args.figure_dir / f"{analysis.lower()}_blp_dimred_method_k_score.png",
            f"{analysis} BLP dimred method-k score | standardized csplit",
        )

    plot_dimred_label_counts(
        all_rows,
        args.figure_dir / "p8_p10_dimred_two_subprocess_label_counts_top50.png",
        "Top dimred hit labels | standardized csplit",
        args.top_label_n,
    )
    plot_dimred_label_counts_by_family(
        all_rows,
        args.figure_dir / "p8_p10_dimred_two_subprocess_label_counts_by_family_top50.png",
        "Top dimred hit labels by feature family | standardized csplit",
        args.top_label_n,
    )
    plot_raw_index_ecdf_by_topn(
        all_rows,
        args.figure_dir / "raw_blp_efun_index_ecdf_by_topN.png",
    )
    plot_raw_timescale_ecdf_by_topn(
        all_rows,
        args.figure_dir / "raw_blp_efun_timescale_ecdf_by_topN_v3_discrete_log.png",
    )
    write_summary(args.output_dir, p8_rows, p10_rows, best_context, top_dimred)
    print(f"P8 rows: {len(p8_rows)}")
    print(f"P10 rows: {len(p10_rows)}")
    print(f"Output: {args.output_dir}")
    print(f"Figures: {args.figure_dir}")


if __name__ == "__main__":
    main()
