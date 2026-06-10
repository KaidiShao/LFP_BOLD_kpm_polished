"""Compact report of strongest P8/P10 xcorr density hits.

The goal is deliberately simple: read the numeric top CSVs, keep the density
identity visible, and sort by absolute peak cross-correlation.  This is meant
as the first look after a new density version such as ``rmsenv_adaptive``.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Iterable


DEFAULT_PROCESSED_ROOT = Path(r"E:\DataPons_processed")
DEFAULT_DATASETS = ("e10gb1", "e10fV1", "e10gh1", "e10gw1", "f12m01")
DEFAULT_RUN_TAGS = ("pv_gsvd100", "pv_gsvd100_ds", "pv_hp100", "pv_roi")
DEFAULT_OUTPUT_DIR = Path("results") / "p8_p10_top_xcorr_density_current_rmsenv_adaptive"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--datasets", nargs="+", default=list(DEFAULT_DATASETS))
    parser.add_argument("--run-tags", nargs="+", default=list(DEFAULT_RUN_TAGS))
    parser.add_argument("--p8-save-tag", default="xcorr_rmsenv_adaptive")
    parser.add_argument("--p10-save-tag", default="dimred_xcorr_rmsenv_adaptive")
    parser.add_argument("--p10-features", nargs="*", default=[])
    parser.add_argument("--p10-method-tags", nargs="*", default=[])
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--top-n", type=int, default=200)
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open("r", newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def as_float(value: object) -> float:
    try:
        out = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return math.nan
    return out if math.isfinite(out) else math.nan


def parse_p10_dir_name(name: str) -> dict[str, object]:
    m = re.match(r"^(pv_.+)__(.+)__([a-z]+)_k(\d+)$", name)
    if not m:
        return {"p9_feature": "", "p9_method": "", "p9_k": ""}
    _run_tag, feature, method, k = m.groups()
    return {"p9_feature": feature, "p9_method": method, "p9_k": int(k)}


def source_level_from_path(path: Path) -> str:
    parts = {p.lower() for p in path.parts}
    if "feature" in parts:
        return "per_density_feature"
    if "density" in parts:
        return "per_density"
    return "combined"


def add_common_fields(
    out: dict[str, object],
    row: dict[str, str],
    source_csv: Path,
) -> dict[str, object]:
    peak_abs = as_float(row.get("peak_abs_corr"))
    peak = as_float(row.get("peak_corr"))
    lag = as_float(row.get("peak_lag_sec"))
    out.update(
        source_level=source_level_from_path(source_csv),
        density_name=row.get("density_name", ""),
        density_type=row.get("density_type", ""),
        density_index=row.get("density_index", ""),
        density_label=row.get("density_label", ""),
        density_file=row.get("density_file", ""),
        density_field=row.get("density_field", ""),
        bold_feature=row.get("bold_feature", ""),
        bold_mode_index=row.get("bold_mode_index", ""),
        bold_component_index=row.get("bold_component_index", ""),
        peak_abs_corr=peak_abs,
        peak_corr=peak,
        peak_lag_sec=lag,
        zero_corr=as_float(row.get("zero_corr")),
        source_csv=str(source_csv),
    )
    return out


def p8_csvs(root: Path, dataset: str, run_tag: str, save_tag: str) -> Iterable[Path]:
    run_dir = root / dataset / "pipeline8_xcorr" / run_tag
    yield run_dir / f"{save_tag}_top.csv"
    yield from sorted((run_dir / "density").glob(f"{save_tag}_top__*.csv"))
    yield from sorted((run_dir / "feature").glob(f"*/{save_tag}_top__*.csv"))


def p10_csvs(root: Path, dataset: str, run_tag: str, save_tag: str) -> Iterable[Path]:
    p10_root = root / dataset / "pipeline10_dimred_xcorr"
    if not p10_root.is_dir():
        return
    for p10_dir in sorted(p10_root.glob(f"{run_tag}__*__*_k*")):
        if not p10_dir.is_dir():
            continue
        meta = parse_p10_dir_name(p10_dir.name)
        if getattr(p10_csvs, "features", None):
            if meta.get("p9_feature") not in getattr(p10_csvs, "features"):
                continue
        if getattr(p10_csvs, "method_tags", None):
            method_tag = f"{meta.get('p9_method')}_k{int(meta.get('p9_k')):02d}" if meta.get("p9_k") != "" else ""
            if method_tag not in getattr(p10_csvs, "method_tags"):
                continue
        yield p10_dir / f"{save_tag}_top.csv"
        yield from sorted((p10_dir / "density").glob(f"{save_tag}_top__*.csv"))
        yield from sorted((p10_dir / "feature").glob(f"*/{save_tag}_top__*.csv"))


def collect_p8(args: argparse.Namespace) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for dataset in args.datasets:
        for run_tag in args.run_tags:
            for csv_path in p8_csvs(args.processed_root, dataset, run_tag, args.p8_save_tag):
                for rank, row in enumerate(read_rows(csv_path), start=1):
                    rows.append(
                        add_common_fields(
                            {
                                "pipeline": "P8",
                                "dataset": dataset,
                                "run_tag": run_tag,
                                "csv_rank": rank,
                                "p9_feature": "",
                                "p9_method": "",
                                "p9_k": "",
                            },
                            row,
                            csv_path,
                        )
                    )
    return rows


def collect_p10(args: argparse.Namespace) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    p10_csvs.features = set(args.p10_features) if args.p10_features else set()  # type: ignore[attr-defined]
    p10_csvs.method_tags = set(args.p10_method_tags) if args.p10_method_tags else set()  # type: ignore[attr-defined]
    for dataset in args.datasets:
        for run_tag in args.run_tags:
            for csv_path in p10_csvs(args.processed_root, dataset, run_tag, args.p10_save_tag):
                p9_meta = parse_p10_dir_name(csv_path.parents[1].name if csv_path.parent.name in {"density", "feature"} else csv_path.parent.name)
                if csv_path.parent.name not in {"density", "feature"} and csv_path.parent.parent.name == "feature":
                    p9_meta = parse_p10_dir_name(csv_path.parents[2].name)
                for rank, row in enumerate(read_rows(csv_path), start=1):
                    base = {
                        "pipeline": "P10",
                        "dataset": dataset,
                        "run_tag": run_tag,
                        "csv_rank": rank,
                    }
                    base.update(p9_meta)
                    rows.append(add_common_fields(base, row, csv_path))
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                keys.append(key)
                seen.add(key)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def top_sorted(rows: list[dict[str, object]], n: int) -> list[dict[str, object]]:
    valid = [r for r in rows if math.isfinite(as_float(r.get("peak_abs_corr")))]
    valid.sort(key=lambda r: as_float(r.get("peak_abs_corr")), reverse=True)
    return valid[:n]


def main() -> None:
    args = parse_args()
    p8_rows = top_sorted(collect_p8(args), args.top_n)
    p10_rows = top_sorted(collect_p10(args), args.top_n)
    all_rows = top_sorted(p8_rows + p10_rows, args.top_n)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "p8_top_xcorr_density_hits.csv", p8_rows)
    write_csv(args.output_dir / "p10_top_xcorr_density_hits.csv", p10_rows)
    write_csv(args.output_dir / "p8_p10_top_xcorr_density_hits.csv", all_rows)

    print(f"P8 rows: {len(p8_rows)}")
    print(f"P10 rows: {len(p10_rows)}")
    print(f"Output dir: {args.output_dir}")
    for row in all_rows[:20]:
        print(
            "{pipeline} {dataset} {run_tag} {source_level} "
            "| r={peak_abs_corr:.4g} lag={peak_lag_sec:.4g} "
            "| density={density_name}[{density_index}:{density_label}] "
            "| bold={bold_feature} m={bold_mode_index} c={bold_component_index} "
            "| p9={p9_feature}/{p9_method}/k{p9_k}".format(**row)
        )


if __name__ == "__main__":
    main()
