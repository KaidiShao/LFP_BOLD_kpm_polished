"""Summarize P8/P10 coupling after restricting dimred densities to P5-stable method-k.

The purpose is to test whether the P8/P10 strict-label conclusion survives a
more conservative source pool:

    all method-k search
    vs
    P5-stable method-k only

By default the stable set is the 6/6 P5 subprocess-gate pass set after
excluding k13m17 from the current seven-dataset standardized complex-split run.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import itertools
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_HITS = (
    Path("results")
    / "pipeline8_10_strict_band_coupling_current"
    / "p8_p10_strict_band_hits_long.csv"
)
DEFAULT_RESULTS_DIR = (
    Path("results") / "pipeline8_10_strict_band_coupling_p5_stable_method_k_20260601"
)
DEFAULT_TOP_N = (1, 3, 5, 10, 20)
DEFAULT_EXCLUDE_DATASETS = ("k13m17",)
DEFAULT_STABLE_METHOD_K = (
    "mds_k04",
    "mds_k05",
    "mds_k06",
    "mds_k07",
    "mds_k08",
    "nmf_k04",
    "umap_k04",
    "umap_k05",
    "umap_k06",
    "umap_k08",
)

OBS_ORDER = ("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean")
PIPELINE_ORDER = ("P8", "P10")
FEATURE_ORDER = ("efun", "deconv_efun")
DENSITY_ORDER = ("event_density", "raw_efun_density", "dimred_efun_density")
LABEL_GROUP_ORDER = ("theta", "RG", "mixed", "inactive", "label_missing", "other")
LABEL_GROUP_COLORS = {
    "theta": "#2ca25f",
    "RG": "#756bb1",
    "mixed": "#b26a1b",
    "inactive": "#d9d9d9",
    "label_missing": "#111111",
    "other": "#8c8c8c",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits-csv", type=Path, default=DEFAULT_HITS)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--top-n-values", nargs="+", type=int, default=list(DEFAULT_TOP_N))
    parser.add_argument("--exclude-datasets", nargs="*", default=list(DEFAULT_EXCLUDE_DATASETS))
    parser.add_argument("--stable-method-k", nargs="+", default=list(DEFAULT_STABLE_METHOD_K))
    parser.add_argument("--max-rows", type=int, default=0)
    return parser.parse_args()


def as_float(value: object) -> float:
    try:
        out = float(str(value))
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


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


def label_group(strict_label: str) -> str:
    if strict_label == "theta_selective":
        return "theta"
    if strict_label in {"ripple_gamma_no_theta", "gamma_selective", "ripple_selective"}:
        return "RG"
    if strict_label == "mixed_or_partial":
        return "mixed"
    if strict_label == "inactive":
        return "inactive"
    if strict_label == "label_missing" or not strict_label:
        return "label_missing"
    return "other"


def context_key(row: dict[str, str]) -> tuple[str, str, str, str, str, str]:
    return (
        row.get("pipeline", ""),
        row.get("dataset", ""),
        row.get("run_tag", ""),
        row.get("bold_feature_family", ""),
        row.get("p9_feature", ""),
        row.get("p9_method_k", ""),
    )


def keep_for_scope(row: dict[str, str], method_scope: str, stable_method_k: set[str]) -> bool:
    if method_scope == "all_method_k":
        return True
    if row.get("density_class") != "dimred_efun_density":
        return True
    return row.get("density_method_k") in stable_method_k


def keep_for_dimred_scope(row: dict[str, str], method_scope: str, stable_method_k: set[str]) -> bool:
    if row.get("density_class") != "dimred_efun_density":
        return False
    if method_scope == "all_method_k":
        return True
    return row.get("density_method_k") in stable_method_k


def slim_row(row: dict[str, str]) -> dict[str, object]:
    return {
        "pipeline": row.get("pipeline", ""),
        "dataset": row.get("dataset", ""),
        "run_tag": row.get("run_tag", ""),
        "bold_observable": row.get("bold_observable", ""),
        "bold_feature_family": row.get("bold_feature_family", ""),
        "p9_feature": row.get("p9_feature", ""),
        "p9_method_k": row.get("p9_method_k", ""),
        "density_class": row.get("density_class", ""),
        "density_method_k": row.get("density_method_k", ""),
        "density_index": row.get("density_index", ""),
        "strict_label": row.get("strict_label", "") or "label_missing",
        "label_group": label_group(row.get("strict_label", "") or "label_missing"),
        "peak_abs_corr": as_float(row.get("peak_abs_corr")),
        "peak_lag_sec": as_float(row.get("peak_lag_sec")),
    }


def push_top(
    heaps: dict[tuple[str, str, tuple[str, str, str, str, str, str]], list[tuple[float, int, dict[str, object]]]],
    *,
    method_scope: str,
    selection_scope: str,
    key: tuple[str, str, str, str, str, str],
    row: dict[str, object],
    sequence_id: int,
    max_n: int,
) -> None:
    heap_key = (method_scope, selection_scope, key)
    heap = heaps.setdefault(heap_key, [])
    score = float(row["peak_abs_corr"])
    item = (score, sequence_id, row)
    if len(heap) < max_n:
        heapq.heappush(heap, item)
    elif score > heap[0][0]:
        heapq.heapreplace(heap, item)


def collect_contexts(args: argparse.Namespace) -> dict[
    tuple[str, str, tuple[str, str, str, str, str, str]],
    list[tuple[float, int, dict[str, object]]],
]:
    stable_method_k = set(args.stable_method_k)
    excluded = {name.lower() for name in args.exclude_datasets}
    max_n = max(args.top_n_values)
    heaps: dict[
        tuple[str, str, tuple[str, str, str, str, str, str]],
        list[tuple[float, int, dict[str, object]]],
    ] = {}
    seq = itertools.count()
    n_rows = 0
    n_kept = 0
    with args.hits_csv.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            n_rows += 1
            if args.max_rows and n_rows > args.max_rows:
                break
            if row.get("dataset", "").lower() in excluded:
                continue
            if row.get("density_condition") not in {"event", "csplit"}:
                continue
            if row.get("bold_feature_family") not in set(FEATURE_ORDER):
                continue
            peak = as_float(row.get("peak_abs_corr"))
            if not math.isfinite(peak):
                continue
            key = context_key(row)
            small = slim_row(row)
            n_kept += 1
            sequence_id = next(seq)
            for method_scope in ("all_method_k", "p5_stable_method_k"):
                if keep_for_scope(row, method_scope, stable_method_k):
                    push_top(
                        heaps,
                        method_scope=method_scope,
                        selection_scope="competed_topN",
                        key=key,
                        row=small,
                        sequence_id=sequence_id,
                        max_n=max_n,
                    )
                if keep_for_dimred_scope(row, method_scope, stable_method_k):
                    push_top(
                        heaps,
                        method_scope=method_scope,
                        selection_scope="dimred_only_topN",
                        key=key,
                        row=small,
                        sequence_id=sequence_id,
                        max_n=max_n,
                    )
    args.results_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        args.results_dir / "run_manifest.csv",
        [
            {
                "hits_csv": str(args.hits_csv),
                "results_dir": str(args.results_dir),
                "n_rows_read": n_rows,
                "n_candidate_rows": n_kept,
                "n_context_heaps": len(heaps),
                "exclude_datasets": " ".join(args.exclude_datasets),
                "stable_method_k": " ".join(args.stable_method_k),
                "top_n_values": " ".join(str(v) for v in args.top_n_values),
            }
        ],
    )
    return heaps


def ranked_rows(
    heaps: dict[tuple[str, str, tuple[str, str, str, str, str, str]], list[tuple[float, int, dict[str, object]]]],
    top_n: int,
) -> Iterable[tuple[str, str, tuple[str, str, str, str, str, str], list[dict[str, object]]]]:
    for (method_scope, selection_scope, key), heap in heaps.items():
        ordered = [item[2] for item in sorted(heap, key=lambda item: (item[0], item[1]), reverse=True)]
        yield method_scope, selection_scope, key, ordered[:top_n]


def summarize_density(
    heaps: dict[tuple[str, str, tuple[str, str, str, str, str, str]], list[tuple[float, int, dict[str, object]]]],
    top_n_values: Sequence[int],
) -> list[dict[str, object]]:
    rows_out: list[dict[str, object]] = []
    for top_n in top_n_values:
        counts: dict[tuple[str, str, str, str], Counter[str]] = defaultdict(Counter)
        totals: Counter[tuple[str, str, str, str]] = Counter()
        for method_scope, selection_scope, _key, rows in ranked_rows(heaps, top_n):
            if selection_scope != "competed_topN":
                continue
            for row in rows:
                group = (
                    method_scope,
                    str(row["pipeline"]),
                    str(row["bold_observable"]),
                    str(row["bold_feature_family"]),
                )
                dclass = str(row["density_class"])
                if dclass not in DENSITY_ORDER:
                    continue
                counts[group][dclass] += 1
                totals[group] += 1
        for (method_scope, pipeline, observable, family), counter in sorted(counts.items()):
            total = totals[(method_scope, pipeline, observable, family)]
            for dclass in DENSITY_ORDER:
                rows_out.append(
                    {
                        "method_scope": method_scope,
                        "selection_scope": "competed_topN",
                        "pipeline": pipeline,
                        "bold_observable": observable,
                        "bold_feature_family": family,
                        "top_n": top_n,
                        "density_class": dclass,
                        "count": counter[dclass],
                        "fraction": counter[dclass] / total if total else math.nan,
                        "n_top_hits": total,
                    }
                )
    return rows_out


def summarize_labels(
    heaps: dict[tuple[str, str, tuple[str, str, str, str, str, str]], list[tuple[float, int, dict[str, object]]]],
    top_n_values: Sequence[int],
) -> list[dict[str, object]]:
    rows_out: list[dict[str, object]] = []
    for top_n in top_n_values:
        counts: dict[tuple[str, str, str, str, str], Counter[str]] = defaultdict(Counter)
        totals: Counter[tuple[str, str, str, str, str]] = Counter()
        for method_scope, selection_scope, _key, rows in ranked_rows(heaps, top_n):
            for row in rows:
                if row["density_class"] != "dimred_efun_density":
                    continue
                group = (
                    method_scope,
                    selection_scope,
                    str(row["pipeline"]),
                    str(row["bold_observable"]),
                    str(row["bold_feature_family"]),
                )
                counts[group][str(row["label_group"])] += 1
                totals[group] += 1
        for (method_scope, selection_scope, pipeline, observable, family), counter in sorted(counts.items()):
            total = totals[(method_scope, selection_scope, pipeline, observable, family)]
            for group_name in LABEL_GROUP_ORDER:
                rows_out.append(
                    {
                        "method_scope": method_scope,
                        "selection_scope": selection_scope,
                        "pipeline": pipeline,
                        "bold_observable": observable,
                        "bold_feature_family": family,
                        "top_n": top_n,
                        "label_group": group_name,
                        "count": counter[group_name],
                        "fraction": counter[group_name] / total if total else math.nan,
                        "n_dimred_hits": total,
                    }
                )
    return rows_out


def add_all_observable_summary(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    out = list(rows)
    grouped: dict[tuple[object, ...], Counter[str]] = defaultdict(Counter)
    totals: Counter[tuple[object, ...]] = Counter()
    value_field = "label_group" if "label_group" in rows[0] else "density_class"
    count_field = "n_dimred_hits" if value_field == "label_group" else "n_top_hits"
    for row in rows:
        key = tuple(
            row[k]
            for k in (
                "method_scope",
                "selection_scope",
                "pipeline",
                "bold_feature_family",
                "top_n",
            )
        )
        grouped[key][str(row[value_field])] += int(row["count"])
    for key, counter in grouped.items():
        total = sum(counter.values())
        for category in (LABEL_GROUP_ORDER if value_field == "label_group" else DENSITY_ORDER):
            method_scope, selection_scope, pipeline, family, top_n = key
            out.append(
                {
                    "method_scope": method_scope,
                    "selection_scope": selection_scope,
                    "pipeline": pipeline,
                    "bold_observable": "all",
                    "bold_feature_family": family,
                    "top_n": top_n,
                    value_field: category,
                    "count": counter[category],
                    "fraction": counter[category] / total if total else math.nan,
                    count_field: total,
                }
            )
    return out


def get_fraction(
    rows: Sequence[dict[str, object]],
    *,
    method_scope: str,
    selection_scope: str,
    pipeline: str,
    observable: str,
    family: str,
    top_n: int,
    label: str,
) -> float:
    for row in rows:
        if (
            row.get("method_scope") == method_scope
            and row.get("selection_scope") == selection_scope
            and row.get("pipeline") == pipeline
            and row.get("bold_observable") == observable
            and row.get("bold_feature_family") == family
            and int(row.get("top_n", -1)) == int(top_n)
            and row.get("label_group") == label
        ):
            return as_float(row.get("fraction"))
    return math.nan


def plot_rg_topn(label_rows: Sequence[dict[str, object]], out_dir: Path, selection_scope: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2), sharex=True, sharey=True)
    for ax, (pipeline, family) in zip(axes.ravel(), [(p, f) for p in PIPELINE_ORDER for f in FEATURE_ORDER]):
        for method_scope, style in [("all_method_k", "--"), ("p5_stable_method_k", "-")]:
            xs = sorted({int(row["top_n"]) for row in label_rows})
            ys = [
                get_fraction(
                    label_rows,
                    method_scope=method_scope,
                    selection_scope=selection_scope,
                    pipeline=pipeline,
                    observable="all",
                    family=family,
                    top_n=x,
                    label="RG",
                )
                for x in xs
            ]
            ax.plot(xs, ys, style, marker="o", label=method_scope.replace("_", " "))
        ax.set_title(f"{pipeline} {family}")
        ax.set_ylim(0, 1)
        ax.grid(True, alpha=0.25)
        ax.set_ylabel("RG fraction")
        ax.set_xlabel("topN")
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"RG label fraction by topN | {selection_scope}")
    fig.tight_layout()
    fig.savefig(out_dir / f"rg_fraction_topN__{selection_scope}.png", dpi=180)
    plt.close(fig)


def plot_stable_top10_by_observable(label_rows: Sequence[dict[str, object]], out_dir: Path, selection_scope: str) -> None:
    subset = [
        row
        for row in label_rows
        if row.get("method_scope") == "p5_stable_method_k"
        and row.get("selection_scope") == selection_scope
        and int(row.get("top_n", -1)) == 10
        and row.get("bold_observable") in OBS_ORDER
    ]
    fig, axes = plt.subplots(4, 2, figsize=(11.5, 10.5), sharex=True)
    for i, observable in enumerate(OBS_ORDER):
        for j, family in enumerate(FEATURE_ORDER):
            ax = axes[i, j]
            y_labels = list(PIPELINE_ORDER)
            left = [0.0 for _ in y_labels]
            for group_name in LABEL_GROUP_ORDER:
                vals = []
                for pipeline in PIPELINE_ORDER:
                    val = get_fraction(
                        subset,
                        method_scope="p5_stable_method_k",
                        selection_scope=selection_scope,
                        pipeline=pipeline,
                        observable=observable,
                        family=family,
                        top_n=10,
                        label=group_name,
                    )
                    vals.append(0.0 if not math.isfinite(val) else val)
                ax.barh(
                    y_labels,
                    vals,
                    left=left,
                    color=LABEL_GROUP_COLORS[group_name],
                    label=group_name if i == 0 and j == 0 else None,
                )
                left = [a + b for a, b in zip(left, vals)]
            ax.set_xlim(0, 1)
            ax.set_title(f"{observable} | {family}", fontsize=9)
            ax.grid(True, axis="x", alpha=0.25)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=6, frameon=False, fontsize=8)
    fig.suptitle(f"P5-stable method-k only: top10 dimred label groups | {selection_scope}")
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    fig.savefig(out_dir / f"stable_method_k_top10_label_groups_by_observable__{selection_scope}.png", dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        args.results_dir / "stable_method_k_list.csv",
        [{"method_k": value} for value in args.stable_method_k],
    )
    heaps = collect_contexts(args)
    density_rows = add_all_observable_summary(summarize_density(heaps, args.top_n_values))
    label_rows = add_all_observable_summary(summarize_labels(heaps, args.top_n_values))
    write_csv(args.results_dir / "density_class_fraction_by_observable.csv", density_rows)
    write_csv(args.results_dir / "dimred_label_group_fraction_by_observable.csv", label_rows)
    plot_rg_topn(label_rows, args.results_dir, "competed_topN")
    plot_rg_topn(label_rows, args.results_dir, "dimred_only_topN")
    plot_stable_top10_by_observable(label_rows, args.results_dir, "competed_topN")
    plot_stable_top10_by_observable(label_rows, args.results_dir, "dimred_only_topN")


if __name__ == "__main__":
    main()
