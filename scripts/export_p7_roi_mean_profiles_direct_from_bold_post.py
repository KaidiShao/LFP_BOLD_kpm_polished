"""Export P7 roi_mean ROI profiles directly from BOLD_POST kpm_modes.

This is a narrow fallback for roi_mean observable runs where the generic
voxel-mapping ROI exporter fails because one ROI label cannot be resolved into
the voxel-space atlas.  For roi_mean runs, the BOLD source space is already one
feature per ROI, so the ROI profile can be read directly from the Koopman mode
matrix without a voxel back-projection step.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from typing import Sequence

import h5py
import numpy as np


HEADER = [
    "dataset",
    "dataset_id",
    "observable",
    "residual_form",
    "run_name",
    "output_dir",
    "autodl_root",
    "best_val_metric",
    "best_outer_epoch",
    "completed_outer_epochs",
    "training_state_path",
    "selection_metric_source",
    "bold_post_file",
    "observable_file",
    "roi_ts_file",
    "sorted_position",
    "mode_label",
    "raw_index",
    "eigenvalue_real",
    "eigenvalue_imag",
    "eigenvalue_abs",
    "n_export_modes",
    "n_rois",
    "roi_index",
    "roi_label",
    "roi_value",
    "roi_value_mode",
    "feature_reduce",
    "source_space_kind",
    "map_coverage_fraction",
    "map_n_covered_voxels",
    "map_n_uncovered_voxels",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=Path("/mnt/e/DataPons_processed"))
    parser.add_argument("--datasets", nargs="+", default=["k13m23"])
    parser.add_argument("--observable", default="roi_mean")
    parser.add_argument("--residual-form", default="projected_vlambda")
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=Path("results") / "pipeline_roi_profile_consistency_current" / "p7_roi_mean_profiles_direct_from_bold_post.csv",
    )
    parser.add_argument(
        "--roi-value-mode",
        choices=["mean_abs", "abs_mean", "real_mean", "imag_mean", "positive_real", "negative_real"],
        default="mean_abs",
    )
    parser.add_argument("--max-modes", type=int, default=0, help="0 means all modes")
    return parser.parse_args()


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
    labels: list[str] = []
    for ref in np.asarray(data).ravel(order="F"):
        labels.append(decode_matlab_string(handle, ref))
    return labels


def read_complex_dataset(handle: h5py.File, path: str) -> np.ndarray:
    data = handle[path][()]
    if data.dtype.fields and {"real", "imag"}.issubset(data.dtype.fields):
        return np.asarray(data["real"], dtype=float) + 1j * np.asarray(data["imag"], dtype=float)
    return np.asarray(data, dtype=float)


def read_scalar(handle: h5py.File, path: str) -> float:
    try:
        data = np.asarray(handle[path][()])
    except KeyError:
        return math.nan
    if data.size == 0:
        return math.nan
    try:
        return float(data.ravel(order="F")[0])
    except (TypeError, ValueError):
        return math.nan


def read_string_path(handle: h5py.File, path: str) -> str:
    try:
        return decode_matlab_string(handle, handle[path])
    except KeyError:
        return ""


def strip_roi_suffix(label: str) -> str:
    if label.endswith("_mean"):
        return label[: -len("_mean")]
    return label


def roi_value(value: complex, mode: str) -> float:
    if mode == "mean_abs":
        return float(abs(value))
    if mode == "abs_mean":
        return float(abs(value))
    if mode == "real_mean":
        return float(np.real(value))
    if mode == "imag_mean":
        return float(np.imag(value))
    if mode == "positive_real":
        return float(max(np.real(value), 0.0))
    if mode == "negative_real":
        return float(max(-np.real(value), 0.0))
    raise ValueError(f"Unsupported roi value mode: {mode}")


def find_bold_post(processed_root: Path, dataset: str, observable: str, residual_form: str) -> Path | None:
    root = processed_root / dataset / "pipeline7_bold_reskoopnet_postprocessing"
    if not root.exists():
        return None
    candidates: list[Path] = []
    for run_dir in root.iterdir():
        if not run_dir.is_dir():
            continue
        name = run_dir.name.lower()
        if observable.lower() not in name or residual_form.lower() not in name:
            continue
        candidates.extend((run_dir / "mat").glob("*_bold_post.mat"))
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def rows_for_bold_post(path: Path, dataset: str, args: argparse.Namespace) -> list[dict[str, object]]:
    with h5py.File(path, "r") as handle:
        labels = [strip_roi_suffix(x) for x in decode_matlab_cellstr(handle, "BOLD_POST/observable/observable_labels")]
        evalues = read_complex_dataset(handle, "BOLD_POST/EDMD_outputs/evalues").reshape(-1)
        modes = read_complex_dataset(handle, "BOLD_POST/EDMD_outputs/kpm_modes")
        if modes.shape[0] != len(labels) and modes.shape[1] == len(labels):
            modes = modes.T
        if modes.shape[0] != len(labels):
            raise ValueError(f"Mode/label shape mismatch for {path}: modes={modes.shape}, labels={len(labels)}")
        order = np.argsort(np.abs(evalues))[::-1]
        if args.max_modes and args.max_modes > 0:
            order = order[: args.max_modes]

        run_name = read_string_path(handle, "BOLD_POST/run_info/run_name") or path.name.replace("_bold_post.mat", "")
        dataset_id = read_string_path(handle, "BOLD_POST/run_info/dataset_id") or dataset
        residual_form = read_string_path(handle, "BOLD_POST/run_info/residual_form") or args.residual_form
        observable = read_string_path(handle, "BOLD_POST/run_info/observable_mode") or args.observable
        output_dir = read_string_path(handle, "BOLD_POST/run_info/output_dir")
        autodl_root = read_string_path(handle, "BOLD_POST/run_info/autodl_root") or "processed_p7_direct_bold_post"
        best_val_metric = read_scalar(handle, "BOLD_POST/run_info/best_val_metric")
        best_outer_epoch = read_scalar(handle, "BOLD_POST/run_info/best_outer_epoch")
        completed_outer_epochs = read_scalar(handle, "BOLD_POST/run_info/completed_outer_epochs")
        selection_metric_source = read_string_path(handle, "BOLD_POST/run_info/selection_metric_source") or "direct_bold_post"
        observable_file = read_string_path(handle, "BOLD_POST/observable_file")

        rows: list[dict[str, object]] = []
        for sorted_pos, raw_zero in enumerate(order, start=1):
            eig = evalues[raw_zero]
            for roi_idx, label in enumerate(labels, start=1):
                value = roi_value(modes[roi_idx - 1, raw_zero], args.roi_value_mode)
                if not math.isfinite(value):
                    continue
                rows.append(
                    {
                        "dataset": dataset,
                        "dataset_id": dataset_id,
                        "observable": observable,
                        "residual_form": residual_form,
                        "run_name": run_name,
                        "output_dir": output_dir,
                        "autodl_root": autodl_root,
                        "best_val_metric": best_val_metric,
                        "best_outer_epoch": best_outer_epoch,
                        "completed_outer_epochs": completed_outer_epochs,
                        "training_state_path": "",
                        "selection_metric_source": selection_metric_source,
                        "bold_post_file": str(path),
                        "observable_file": observable_file,
                        "roi_ts_file": "",
                        "sorted_position": sorted_pos,
                        "mode_label": f"sorted{sorted_pos:03d}",
                        "raw_index": int(raw_zero + 1),
                        "eigenvalue_real": float(np.real(eig)),
                        "eigenvalue_imag": float(np.imag(eig)),
                        "eigenvalue_abs": float(abs(eig)),
                        "n_export_modes": int(len(order)),
                        "n_rois": int(len(labels)),
                        "roi_index": roi_idx,
                        "roi_label": label,
                        "roi_value": value,
                        "roi_value_mode": args.roi_value_mode,
                        "feature_reduce": "direct_roi_mean",
                        "source_space_kind": "roi_mean_direct_from_bold_post",
                        "map_coverage_fraction": "",
                        "map_n_covered_voxels": "",
                        "map_n_uncovered_voxels": "",
                    }
                )
        return rows


def write_csv(path: Path, rows: Sequence[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in HEADER})


def main() -> None:
    args = parse_args()
    rows: list[dict[str, object]] = []
    for dataset in args.datasets:
        bold_post = find_bold_post(args.processed_root, dataset, args.observable, args.residual_form)
        if bold_post is None:
            print(f"MISSING {dataset}: no BOLD_POST found")
            continue
        dataset_rows = rows_for_bold_post(bold_post, dataset, args)
        rows.extend(dataset_rows)
        print(f"{dataset}: {len(dataset_rows)} rows from {bold_post}")
    write_csv(args.output_csv, rows)
    print(f"Wrote {len(rows)} rows to {args.output_csv}")


if __name__ == "__main__":
    main()
