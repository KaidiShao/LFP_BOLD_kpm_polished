#!/usr/bin/env python3
"""Plot E10gb1 P5-only component-label preview figures.

This script intentionally ignores P8/P10.  It visualizes whether the P5
standardized complex-split dimred components contain separable theta-like and
ripple-gamma-like subprocess candidates.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch, Rectangle


METHOD_ORDER = ("svd", "nmf", "mds", "umap")
K_VALUES = tuple(range(3, 9))
BASIC_LABELS = (
    "theta_selective",
    "ripple_selective",
    "gamma_selective",
    "mixed_selective",
    "nonselective_active",
    "other",
    "label_missing",
)
BASIC_COLORS = {
    "theta_selective": "#31a354",
    "ripple_selective": "#3182bd",
    "gamma_selective": "#f28e2b",
    "mixed_selective": "#9467bd",
    "nonselective_active": "#8c8c8c",
    "other": "#d9d9d9",
    "label_missing": "#ffffff",
}
ACTIVITY_LABELS = (
    "theta_selective",
    "ripple_gamma_no_pure_theta",
    "ripple_gamma_joint_with_theta",
    "pan_event",
    "partial_or_inactive",
    "unlabeled",
)
ACTIVITY_COLORS = {
    "theta_selective": "#2ca25f",
    "ripple_gamma_no_pure_theta": "#2171b5",
    "ripple_gamma_joint_with_theta": "#6baed6",
    "pan_event": "#a6611a",
    "partial_or_inactive": "#737373",
    "unlabeled": "#f0f0f0",
}
BASIC_CODES = {
    "theta_selective": "Tsel",
    "ripple_selective": "Rsel",
    "gamma_selective": "Gsel",
    "mixed_selective": "Mix",
    "nonselective_active": "Act",
    "other": "Oth",
    "label_missing": "-",
}


DEFAULT_LABEL_TABLE = (
    Path("results")
    / "pipeline5_dimred_component_process_labels_e10gb1_standardize_rmsenv_adaptive_all_components"
    / "dimred_efun_process_labels.csv"
)
DEFAULT_ACTIVITY_TABLE = (
    Path("results")
    / "peak_event_family_component_activity_e10gb1_standardize_rmsenv_adaptive"
    / "component_family_activity.csv"
)
DEFAULT_OUTPUT_DIR = Path("results") / "e10gb1_p5_label_preview_20260528"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label-table", type=Path, default=DEFAULT_LABEL_TABLE)
    parser.add_argument("--activity-table", type=Path, default=DEFAULT_ACTIVITY_TABLE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--dataset", default="e10gb1")
    parser.add_argument(
        "--condition",
        default="complex_split_projected_vlambda_standardize_rmsenv_adaptive",
    )
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


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
        writer.writerows(rows)


def truthy(value: object) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def split_items(value: object) -> set[str]:
    return {item.strip().lower() for item in str(value or "").split(";") if item.strip()}


def filter_labels(rows: Sequence[dict[str, str]], dataset: str, condition: str) -> list[dict[str, str]]:
    out = []
    for row in rows:
        if str(row.get("dataset", "")).lower() != dataset.lower():
            continue
        if str(row.get("condition", "")) != condition:
            continue
        out.append(row)
    return out


def build_activity_lookup(rows: Sequence[dict[str, str]], dataset: str, condition: str) -> dict[tuple[str, int, int], dict[str, bool]]:
    lookup: dict[tuple[str, int, int], dict[str, bool]] = defaultdict(dict)
    for row in rows:
        if str(row.get("dataset", "")).lower() != dataset.lower():
            continue
        if str(row.get("condition", "")) != condition:
            continue
        family = str(row.get("event_family", "")).lower()
        if family not in {"theta", "ripple", "gamma"}:
            continue
        key = (
            str(row.get("method", "")).lower(),
            as_int(row.get("component_count")),
            as_int(row.get("component_idx")),
        )
        active_events = split_items(row.get("active_events"))
        lookup[key][f"{family}_any_active"] = bool(active_events)
        lookup[key][f"{family}_all_family_active"] = truthy(row.get("all_family_events_active"))
        if family == "theta":
            lookup[key]["pure_theta_active"] = "theta" in active_events
    return lookup


def basic_label(row: dict[str, object]) -> str:
    selective = [
        fam
        for fam in ("theta", "ripple", "gamma")
        if bool(row.get(f"{fam}_selective"))
    ]
    active = [
        fam
        for fam in ("theta", "ripple", "gamma")
        if bool(row.get(f"{fam}_active"))
    ]
    if len(selective) == 1:
        return f"{selective[0]}_selective"
    if len(selective) > 1:
        return "mixed_selective"
    if active:
        return "nonselective_active"
    return "other"


def activity_label(row: dict[str, object]) -> str:
    if bool(row.get("ripple_gamma_no_pure_theta")):
        return "ripple_gamma_no_pure_theta"
    if bool(row.get("pan_event")):
        return "pan_event"
    if bool(row.get("ripple_gamma_joint")):
        return "ripple_gamma_joint_with_theta"
    if bool(row.get("theta_selective")):
        return "theta_selective"
    if bool(row.get("theta_active")) or bool(row.get("ripple_active")) or bool(row.get("gamma_active")):
        return "partial_or_inactive"
    return "unlabeled"


def build_component_rows(labels: Sequence[dict[str, str]], activity_lookup: dict[tuple[str, int, int], dict[str, bool]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in labels:
        method = str(row.get("method", "")).lower()
        k = as_int(row.get("k") or row.get("component_count"))
        comp = as_int(row.get("component_idx"))
        activity = activity_lookup.get((method, k, comp), {})
        out: dict[str, object] = {
            "dataset": row.get("dataset", ""),
            "condition": row.get("condition", ""),
            "method": method,
            "k": k,
            "method_k": f"{method}_k{k:02d}",
            "component_idx": comp,
            "primary_process_label": row.get("primary_process_label", ""),
            "theta_active": bool(activity.get("theta_any_active", truthy(row.get("theta_active")))),
            "theta_all_family_active": bool(activity.get("theta_all_family_active", truthy(row.get("theta_active")))),
            "theta_selective": truthy(row.get("theta_selective")),
            "ripple_active": bool(activity.get("ripple_any_active", truthy(row.get("ripple_active")))),
            "ripple_all_family_active": bool(activity.get("ripple_all_family_active", truthy(row.get("ripple_active")))),
            "ripple_selective": truthy(row.get("ripple_selective")),
            "gamma_active": bool(activity.get("gamma_any_active", truthy(row.get("gamma_active")))),
            "gamma_all_family_active": bool(activity.get("gamma_all_family_active", truthy(row.get("gamma_active")))),
            "gamma_selective": truthy(row.get("gamma_selective")),
        }
        out["pure_theta_active"] = bool(activity.get("pure_theta_active", False))
        out["ripple_gamma_joint"] = bool(out["ripple_active"]) and bool(out["gamma_active"])
        out["ripple_gamma_no_pure_theta"] = bool(out["ripple_gamma_joint"]) and not bool(out["pure_theta_active"])
        out["ripple_gamma_no_theta_selective"] = bool(out["ripple_gamma_joint"]) and not bool(out["theta_selective"])
        out["pan_event"] = bool(out["theta_active"]) and bool(out["ripple_active"]) and bool(out["gamma_active"])
        out["basic_selective_label"] = basic_label(out)
        out["activity_subprocess_label"] = activity_label(out)
        rows.append(out)
    return rows


def method_k_summary(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["method"]), int(row["k"]))].append(row)
    out: list[dict[str, object]] = []
    for method in METHOD_ORDER:
        for k in K_VALUES:
            group = grouped.get((method, k), [])
            item = {
                "method": method,
                "k": k,
                "method_k": f"{method}_k{k:02d}",
                "n_components": len(group),
                "n_theta_selective": sum(bool(r["theta_selective"]) for r in group),
                "n_ripple_selective": sum(bool(r["ripple_selective"]) for r in group),
                "n_gamma_selective": sum(bool(r["gamma_selective"]) for r in group),
                "n_mixed_selective": sum(str(r["basic_selective_label"]) == "mixed_selective" for r in group),
                "n_ripple_gamma_joint": sum(bool(r["ripple_gamma_joint"]) for r in group),
                "n_ripple_gamma_no_pure_theta": sum(bool(r["ripple_gamma_no_pure_theta"]) for r in group),
                "n_ripple_gamma_no_theta_selective": sum(bool(r["ripple_gamma_no_theta_selective"]) for r in group),
                "n_pan_event": sum(bool(r["pan_event"]) for r in group),
                "n_partial_or_inactive": sum(str(r["activity_subprocess_label"]) == "partial_or_inactive" for r in group),
            }
            item["has_theta_candidate"] = int(item["n_theta_selective"] > 0)
            item["has_ripple_gamma_candidate"] = int(item["n_ripple_gamma_no_pure_theta"] > 0)
            item["has_two_subprocess_candidate"] = int(item["has_theta_candidate"] and item["has_ripple_gamma_candidate"])
            if item["has_two_subprocess_candidate"]:
                item["candidate_state"] = "theta+ripple_gamma"
            elif item["has_theta_candidate"]:
                item["candidate_state"] = "theta_only"
            elif item["has_ripple_gamma_candidate"]:
                item["candidate_state"] = "ripple_gamma_only"
            else:
                item["candidate_state"] = "none"
            out.append(item)
    return out


def ordered_method_ks() -> list[str]:
    return [f"{method}_k{k:02d}" for method in METHOD_ORDER for k in K_VALUES]


def add_legend(fig, color_map: dict[str, str], labels: Iterable[str], loc: str = "upper center") -> None:
    handles = [Patch(facecolor=color_map[label], edgecolor="white", label=label) for label in labels]
    fig.legend(handles=handles, loc=loc, ncol=min(4, len(handles)), frameon=False)


def save_fig(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_basic_label_composition(rows: Sequence[dict[str, object]], out: Path) -> None:
    df = pd.DataFrame(rows)
    order = ordered_method_ks()
    counts = pd.crosstab(df["method_k"], df["basic_selective_label"]).reindex(order).fillna(0)
    for label in BASIC_LABELS:
        if label not in counts:
            counts[label] = 0
    counts = counts[list(BASIC_LABELS)]
    denom = counts.sum(axis=1).replace(0, np.nan)
    frac = counts.div(denom, axis=0).fillna(0)
    fig, ax = plt.subplots(figsize=(16, 5))
    bottom = np.zeros(len(frac))
    x = np.arange(len(frac))
    for label in BASIC_LABELS:
        values = frac[label].to_numpy()
        ax.bar(x, values, bottom=bottom, color=BASIC_COLORS[label], label=label, edgecolor="white", linewidth=0.4)
        bottom += values
    ax.set_xticks(x, order, rotation=50, ha="right")
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("fraction of components")
    ax.set_title("E10gb1 P5 standardized csplit: basic selective label composition")
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    save_fig(fig, out)


def heatmap_panel(ax, matrix: np.ndarray, title: str, vmax: float) -> None:
    im = ax.imshow(matrix, cmap="YlGnBu", vmin=0, vmax=vmax)
    ax.set_xticks(np.arange(len(K_VALUES)), [f"k{k:02d}" for k in K_VALUES])
    ax.set_yticks(np.arange(len(METHOD_ORDER)), METHOD_ORDER)
    ax.set_title(title, fontsize=10)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]
            ax.text(j, i, f"{int(val)}", ha="center", va="center", color="black", fontsize=9)
    return im


def plot_activity_heatmaps(summary: Sequence[dict[str, object]], out: Path) -> None:
    fields = [
        ("n_theta_selective", "theta selective"),
        ("n_ripple_selective", "ripple selective"),
        ("n_gamma_selective", "gamma selective"),
        ("n_ripple_gamma_joint", "ripple-gamma joint"),
        ("n_ripple_gamma_no_pure_theta", "ripple-gamma no pure theta"),
        ("n_pan_event", "pan-event"),
    ]
    df = pd.DataFrame(summary)
    vmax = max(float(df[field].max()) for field, _ in fields)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    for ax, (field, title) in zip(axes.ravel(), fields):
        mat = np.zeros((len(METHOD_ORDER), len(K_VALUES)))
        for i, method in enumerate(METHOD_ORDER):
            for j, k in enumerate(K_VALUES):
                match = df[(df.method == method) & (df.k == k)]
                mat[i, j] = float(match[field].iloc[0]) if not match.empty else 0
        im = heatmap_panel(ax, mat, title, vmax)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8)
    cbar.set_label("component count")
    fig.suptitle("E10gb1 P5 standardized csplit: subprocess/selectivity counts by method-k", fontsize=14, weight="bold")
    save_fig(fig, out)


def activity_code(row: dict[str, object]) -> str:
    label = str(row["activity_subprocess_label"])
    if label == "ripple_gamma_no_pure_theta":
        return "RG-noT"
    if label == "ripple_gamma_joint_with_theta":
        return "RG+T"
    if label == "pan_event":
        return "PAN"
    if label == "theta_selective":
        return "Tsel"
    if label == "partial_or_inactive":
        parts = []
        if row.get("theta_active"):
            parts.append("T")
        if row.get("ripple_active"):
            parts.append("R")
        if row.get("gamma_active"):
            parts.append("G")
        return "".join(parts) if parts else "part"
    return "-"


def basic_code(row: dict[str, object]) -> str:
    return BASIC_CODES.get(str(row.get("basic_selective_label", "")), "-")


def draw_component_tile(
    ax,
    x: float,
    y: float,
    width: float,
    height: float,
    row: dict[str, object],
    fontsize: float = 7,
) -> None:
    basic = str(row.get("basic_selective_label", "label_missing"))
    activity = str(row.get("activity_subprocess_label", "unlabeled"))
    face = BASIC_COLORS.get(basic, "#ffffff")
    stripe = ACTIVITY_COLORS.get(activity, "#eeeeee")
    ax.add_patch(Rectangle((x, y), width, height, facecolor=face, edgecolor="#333333", linewidth=0.7))
    ax.add_patch(Rectangle((x, y + height * 0.70), width, height * 0.30, facecolor=stripe, edgecolor="none"))
    ax.text(
        x + width / 2,
        y + height * 0.34,
        f"C{int(row['component_idx'])} {basic_code(row)}",
        ha="center",
        va="center",
        fontsize=fontsize,
        color="black",
    )
    ax.text(
        x + width / 2,
        y + height * 0.84,
        activity_code(row),
        ha="center",
        va="center",
        fontsize=fontsize,
        color="black",
    )


def add_two_layer_legends(fig) -> None:
    basic_handles = [
        Patch(facecolor=BASIC_COLORS[label], edgecolor="#333333", label=f"basic: {label}")
        for label in BASIC_LABELS
    ]
    activity_handles = [
        Patch(facecolor=ACTIVITY_COLORS[label], edgecolor="none", label=f"activity: {label}")
        for label in ACTIVITY_LABELS
    ]
    fig.legend(
        handles=basic_handles,
        loc="lower left",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.02, 0.005),
        fontsize=8,
    )
    fig.legend(
        handles=activity_handles,
        loc="lower right",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.98, 0.005),
        fontsize=8,
    )


def plot_component_grid(rows: Sequence[dict[str, object]], out: Path) -> None:
    row_order = ordered_method_ks()
    label_to_idx = {label: i for i, label in enumerate(BASIC_LABELS)}
    data = np.full((len(row_order), 8), label_to_idx["label_missing"], dtype=float)
    row_lookup = {(str(r["method_k"]), int(r["component_idx"])): r for r in rows}
    for i, mk in enumerate(row_order):
        for comp in range(1, 9):
            row = row_lookup.get((mk, comp))
            if row:
                data[i, comp - 1] = label_to_idx[str(row["basic_selective_label"])]
    cmap = ListedColormap([BASIC_COLORS[label] for label in BASIC_LABELS])
    fig, ax = plt.subplots(figsize=(10, 13))
    ax.imshow(data, cmap=cmap, vmin=-0.5, vmax=len(BASIC_LABELS) - 0.5, aspect="auto")
    ax.set_xticks(np.arange(8), [f"C{i}" for i in range(1, 9)])
    ax.set_yticks(np.arange(len(row_order)), row_order)
    ax.set_title("E10gb1 P5 standardized csplit: component label grid")
    for i, mk in enumerate(row_order):
        k = int(mk[-2:])
        for comp in range(1, 9):
            row = row_lookup.get((mk, comp))
            if comp > k:
                ax.add_patch(Rectangle((comp - 1.5, i - 0.5), 1, 1, facecolor="#f7f7f7", edgecolor="#eeeeee", hatch="//"))
                continue
            if row:
                ax.text(comp - 1, i, activity_code(row), ha="center", va="center", fontsize=7, color="black")
            ax.add_patch(Rectangle((comp - 1.5, i - 0.5), 1, 1, fill=False, edgecolor="white", linewidth=0.7))
    ax.set_xlabel("component index")
    add_legend(fig, BASIC_COLORS, BASIC_LABELS, loc="lower center")
    fig.subplots_adjust(bottom=0.12)
    save_fig(fig, out)


def plot_two_subprocess_map(summary: Sequence[dict[str, object]], out: Path) -> None:
    state_order = ("none", "theta_only", "ripple_gamma_only", "theta+ripple_gamma")
    state_colors = {
        "none": "#f0f0f0",
        "theta_only": "#a1d99b",
        "ripple_gamma_only": "#9ecae1",
        "theta+ripple_gamma": "#2b8cbe",
    }
    state_to_idx = {state: i for i, state in enumerate(state_order)}
    df = pd.DataFrame(summary)
    mat = np.zeros((len(METHOD_ORDER), len(K_VALUES)))
    label_mat: list[list[str]] = [["" for _ in K_VALUES] for _ in METHOD_ORDER]
    for i, method in enumerate(METHOD_ORDER):
        for j, k in enumerate(K_VALUES):
            row = df[(df.method == method) & (df.k == k)].iloc[0]
            state = str(row["candidate_state"])
            mat[i, j] = state_to_idx[state]
            if state == "theta+ripple_gamma":
                label_mat[i][j] = "T+RG"
            elif state == "theta_only":
                label_mat[i][j] = "T"
            elif state == "ripple_gamma_only":
                label_mat[i][j] = "RG"
            else:
                label_mat[i][j] = "-"
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    cmap = ListedColormap([state_colors[state] for state in state_order])
    ax.imshow(mat, cmap=cmap, vmin=-0.5, vmax=len(state_order) - 0.5)
    ax.set_xticks(np.arange(len(K_VALUES)), [f"k{k:02d}" for k in K_VALUES])
    ax.set_yticks(np.arange(len(METHOD_ORDER)), METHOD_ORDER)
    ax.set_title("E10gb1 P5 standardized csplit: two-subprocess candidate map")
    for i in range(len(METHOD_ORDER)):
        for j in range(len(K_VALUES)):
            ax.text(j, i, label_mat[i][j], ha="center", va="center", fontsize=10, weight="bold")
    add_legend(fig, state_colors, state_order, loc="lower center")
    fig.subplots_adjust(bottom=0.22)
    save_fig(fig, out)


def plot_count_scatter(summary: Sequence[dict[str, object]], out: Path) -> None:
    df = pd.DataFrame(summary)
    method_colors = {"svd": "#1f77b4", "nmf": "#ff7f0e", "mds": "#2ca02c", "umap": "#9467bd"}
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for method in METHOD_ORDER:
        sub = df[df.method == method]
        ax.scatter(
            sub["n_theta_selective"],
            sub["n_ripple_gamma_no_pure_theta"],
            s=35 + 18 * sub["k"],
            color=method_colors[method],
            alpha=0.75,
            edgecolor="white",
            linewidth=0.8,
            label=method,
        )
        for _, row in sub.iterrows():
            ax.text(row["n_theta_selective"] + 0.04, row["n_ripple_gamma_no_pure_theta"] + 0.04, f"k{int(row['k']):02d}", fontsize=8)
    ax.axvline(0.5, color="#cccccc", linewidth=1)
    ax.axhline(0.5, color="#cccccc", linewidth=1)
    ax.set_xlabel("theta selective component count")
    ax.set_ylabel("ripple-gamma no-pure-theta component count")
    ax.set_title("E10gb1 P5 standardized csplit: theta vs ripple-gamma candidate counts")
    ax.grid(color="#dddddd", linewidth=0.6)
    ax.legend(frameon=False, title="method")
    save_fig(fig, out)


def plot_pairing_diagram(rows: Sequence[dict[str, object]], summary: Sequence[dict[str, object]], out: Path) -> None:
    summary_df = pd.DataFrame(summary)
    candidate_mks = summary_df[summary_df.has_two_subprocess_candidate == 1]["method_k"].tolist()
    row_lookup = {(str(r["method_k"]), int(r["component_idx"])): r for r in rows}
    n_rows = len(candidate_mks)
    fig, ax = plt.subplots(figsize=(11, max(5, 0.42 * n_rows + 1.5)))
    ax.set_xlim(0, 9.8)
    ax.set_ylim(-0.5, n_rows - 0.5)
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title("E10gb1 P5 standardized csplit: candidate component pairing diagram", fontsize=13, weight="bold")
    for i, mk in enumerate(candidate_mks):
        ax.text(0, i, mk, va="center", ha="left", fontsize=9, weight="bold")
        k = int(mk[-2:])
        for comp in range(1, k + 1):
            row = row_lookup.get((mk, comp))
            if not row:
                continue
            x = 1.35 + (comp - 1) * 1.0
            draw_component_tile(ax, x, i - 0.32, 0.88, 0.64, row, fontsize=6.5)
    ax.text(
        1.35,
        -0.75,
        "tile: background = basic selective label; bottom stripe = activity/subprocess label",
        fontsize=8,
        ha="left",
        va="center",
        color="#333333",
    )
    add_two_layer_legends(fig)
    fig.subplots_adjust(bottom=0.20)
    save_fig(fig, out)


def choose_best_by_method(summary: Sequence[dict[str, object]]) -> list[str]:
    df = pd.DataFrame(summary)
    selected: list[str] = []
    for method in METHOD_ORDER:
        sub = df[(df.method == method) & (df.has_two_subprocess_candidate == 1)].copy()
        if sub.empty:
            sub = df[df.method == method].copy()
        sub["score"] = (
            100 * sub["has_two_subprocess_candidate"]
            + 5 * sub["n_ripple_gamma_no_pure_theta"]
            + 3 * sub["n_theta_selective"]
            + sub["n_ripple_gamma_joint"]
            - 0.01 * sub["k"]
        )
        best = sub.sort_values(["score", "k"], ascending=[False, True]).iloc[0]
        selected.append(str(best["method_k"]))
    return selected


def plot_best_contact_sheet(rows: Sequence[dict[str, object]], summary: Sequence[dict[str, object]], out: Path) -> None:
    selected_mks = choose_best_by_method(summary)
    row_lookup = {(str(r["method_k"]), int(r["component_idx"])): r for r in rows}
    summary_lookup = {str(r["method_k"]): r for r in summary}
    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.set_xlim(0, 9.9)
    ax.set_ylim(-0.5, len(selected_mks) - 0.5)
    ax.invert_yaxis()
    ax.axis("off")
    ax.set_title("E10gb1 P5 standardized csplit: best per-method candidate contact sheet", fontsize=13, weight="bold")
    for i, mk in enumerate(selected_mks):
        srow = summary_lookup[mk]
        subtitle = (
            f"{mk}  |  theta={srow['n_theta_selective']}  "
            f"RG_noT={srow['n_ripple_gamma_no_pure_theta']}  "
            f"RG_joint={srow['n_ripple_gamma_joint']}"
        )
        ax.text(0, i, subtitle, va="center", ha="left", fontsize=9, weight="bold")
        k = int(mk[-2:])
        for comp in range(1, k + 1):
            row = row_lookup.get((mk, comp))
            if not row:
                continue
            x = 3.0 + (comp - 1) * 0.82
            draw_component_tile(ax, x, i - 0.31, 0.74, 0.62, row, fontsize=6.2)
    ax.text(
        3.0,
        -0.72,
        "tile: background = basic selective label; bottom stripe = activity/subprocess label",
        fontsize=8,
        ha="left",
        va="center",
        color="#333333",
    )
    add_two_layer_legends(fig)
    fig.subplots_adjust(bottom=0.24)
    save_fig(fig, out)


def main() -> None:
    args = parse_args()
    labels = filter_labels(read_csv_rows(args.label_table), args.dataset, args.condition)
    activity = read_csv_rows(args.activity_table)
    activity_lookup = build_activity_lookup(activity, args.dataset, args.condition)
    component_rows = build_component_rows(labels, activity_lookup)
    summary_rows = method_k_summary(component_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = args.output_dir / "figures"
    write_csv(args.output_dir / "component_label_preview_rows.csv", component_rows)
    write_csv(args.output_dir / "method_k_label_summary.csv", summary_rows)

    plot_basic_label_composition(component_rows, fig_dir / "01_basic_selective_label_composition_by_method_k.png")
    plot_activity_heatmaps(summary_rows, fig_dir / "02_activity_subprocess_count_heatmaps_by_method_k.png")
    plot_component_grid(component_rows, fig_dir / "03_component_label_grid_by_method_k.png")
    plot_two_subprocess_map(summary_rows, fig_dir / "04_two_subprocess_candidate_map.png")
    plot_count_scatter(summary_rows, fig_dir / "05_theta_vs_ripple_gamma_count_scatter.png")
    plot_pairing_diagram(component_rows, summary_rows, fig_dir / "06_pairing_diagram_all_two_subprocess_candidates.png")
    plot_best_contact_sheet(component_rows, summary_rows, fig_dir / "07_best_candidate_contact_sheet_by_method.png")

    summary_df = pd.DataFrame(summary_rows)
    n_candidates = int(summary_df["has_two_subprocess_candidate"].sum())
    lines = [
        "# E10gb1 P5 label preview figures",
        "",
        f"- dataset: `{args.dataset}`",
        f"- condition: `{args.condition}`",
        f"- components: `{len(component_rows)}`",
        f"- method-k two-subprocess candidates: `{n_candidates}` / `{len(summary_rows)}`",
        "",
        "## Figures",
        "",
    ]
    for png in sorted(fig_dir.glob("*.png")):
        lines.append(f"- `{png.name}`")
    (args.output_dir / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote P5 label preview figures: {fig_dir}")


if __name__ == "__main__":
    main()
