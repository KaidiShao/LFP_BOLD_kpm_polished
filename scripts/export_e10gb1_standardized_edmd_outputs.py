#!/usr/bin/env python3
"""Export full e10gb1 standardized BLP EDMD outputs from best checkpoints.

The 2026-05-23 standardized BLP runs were exported as compact summaries only.
This script restores the best checkpoints, streams the full observable matrix in
chunks, computes EDMD eigenfunctions, and writes normal chunked EDMD output MAT
files into a separate derived autodl-style root.
"""

from __future__ import annotations

import argparse
import gc
import json
import shutil
import sys
from pathlib import Path

import h5py
import numpy as np
import scipy.io


REPO_ROOT = Path(__file__).resolve().parents[1]
TMP_DIR = REPO_ROOT / "tmp"
if str(TMP_DIR) not in sys.path:
    sys.path.insert(0, str(TMP_DIR))

import plot_e10gb1_blp_raw_vs_std_spectra_loss_efun_20260523 as raw_vs_std


DEFAULT_OUTPUT_ROOT = Path(
    "/mnt/e/DataPons_processed/derived_autodl_results_standardize"
)
RUN_KEYS = ("abs_std", "complex_std")

CHUNK_EXCLUDED_SHARED_FIELDS = {
    "snapshot_valid_idx_x",
    "snapshot_valid_idx_y",
    "snapshot_session_idx",
    "snapshot_train_pair_idx",
    "snapshot_valid_pair_idx",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stream-export full standardized e10gb1 BLP EDMD outputs."
    )
    parser.add_argument(
        "--run",
        choices=("both", *RUN_KEYS),
        default="both",
        help="Which standardized run to export.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="Derived autodl-style root; outputs go under <root>/e10gb1/mlp/outputs/<run>.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=5000,
        help="Number of time points per exported MAT chunk.",
    )
    parser.add_argument(
        "--export-precision",
        choices=("single", "double"),
        default="single",
        help="Precision for saved efuns. Single keeps disk use manageable for P5.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Remove existing output chunks for the selected run before export.",
    )
    parser.add_argument(
        "--max-chunks",
        type=int,
        default=None,
        help="Debug option: export only the first N chunks.",
    )
    return parser.parse_args()


def wsl_to_win(path: Path | str) -> str:
    text = str(path)
    if text.startswith("/mnt/") and len(text) >= 7:
        drive = text[5].upper()
        rest = text[7:].replace("/", "\\")
        return f"{drive}:\\{rest}"
    return text


def run_name_from_spec(spec: raw_vs_std.RunSpec) -> str:
    return spec.output_summary.parent.name


def dataset_from_spec(spec: raw_vs_std.RunSpec) -> str:
    parts = spec.data_file.parts
    for i, part in enumerate(parts[:-1]):
        if part == "DataPons_processed" and i + 1 < len(parts):
            return parts[i + 1]
    return "e10gb1"


def load_summary_payload(summary_path: Path) -> tuple[dict, dict | None]:
    mat = scipy.io.loadmat(summary_path, simplify_cells=True)
    if "EDMD_outputs" not in mat:
        raise KeyError(f"EDMD_outputs missing from {summary_path}")
    shared = dict(mat["EDMD_outputs"])
    shared.pop("efuns", None)
    normalization = mat.get("Data_normalization")
    if normalization is not None:
        normalization = dict(normalization)
    return shared, normalization


def n_samples_in_file(spec: raw_vs_std.RunSpec) -> int:
    with h5py.File(spec.data_file, "r") as handle:
        return int(handle[raw_vs_std.DATASET_KEY].shape[1])


def clear_old_outputs(output_dir: Path) -> None:
    if not output_dir.exists():
        return
    for path in output_dir.glob("*_outputs_*.mat"):
        path.unlink()
    for path in output_dir.glob("*_summary.mat"):
        path.unlink()
    for name in (
        "training_loss_diagnostics.png",
        "training_inner_loss_metrics.png",
        "export_manifest.json",
    ):
        path = output_dir / name
        if path.exists():
            path.unlink()


def export_one(
    run_key: str,
    output_root: Path,
    chunk_size: int,
    export_precision: str,
    overwrite: bool,
    max_chunks: int | None,
) -> dict:
    spec = raw_vs_std.RUNS[run_key]
    if not spec.standardize:
        raise ValueError(f"{run_key} is not a standardized run")
    if spec.output_summary is None or not spec.output_summary.exists():
        raise FileNotFoundError(f"Missing summary MAT for {run_key}: {spec.output_summary}")

    run_name = run_name_from_spec(spec)
    dataset = dataset_from_spec(spec)
    output_dir = output_root / dataset / "mlp" / "outputs" / run_name
    output_dir.mkdir(parents=True, exist_ok=True)
    if overwrite:
        clear_old_outputs(output_dir)

    shared, normalization = load_summary_payload(spec.output_summary)
    if normalization is None:
        raise RuntimeError(f"Data_normalization missing from {spec.output_summary}")

    # Keep summary and training diagnostics next to the derived chunks so the
    # derived output directory looks like a normal P4 output directory.
    summary_copy = output_dir / spec.output_summary.name
    if not summary_copy.exists() or overwrite:
        scipy.io.savemat(
            summary_copy,
            {"EDMD_outputs": shared, "Data_normalization": normalization},
            long_field_names=True,
        )
    for png in ("training_loss_diagnostics.png", "training_inner_loss_metrics.png"):
        source_png = spec.output_summary.parent / png
        if source_png.exists():
            target_png = output_dir / png
            if not target_png.exists() or overwrite:
                shutil.copy2(source_png, target_png)

    existing_chunks = sorted(output_dir.glob("*_outputs_*.mat"))
    manifest_path = output_dir / "export_manifest.json"
    if existing_chunks and not overwrite:
        old_manifest = {}
        if manifest_path.exists():
            try:
                old_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                old_manifest = {}
        if old_manifest.get("complete") is True or old_manifest.get("status") == "complete":
            print(
                f"[{run_key}] complete existing chunks found in {output_dir}; "
                "skipping.",
                flush=True,
            )
            return {
                "run_key": run_key,
                "run_name": run_name,
                "output_dir": str(output_dir),
                "status": "exists",
                "n_existing_chunks": len(existing_chunks),
            }
        print(
            f"[{run_key}] found {len(existing_chunks)} existing chunk(s); "
            "resuming missing chunks.",
            flush=True,
        )

    print(f"[{run_key}] restoring checkpoint", flush=True)
    solver, checkpoint_path = raw_vs_std.restore_solver(spec)
    normalizer = raw_vs_std.load_norm(spec.output_summary)
    if normalizer is None:
        raise RuntimeError(f"Could not load normalizer for {run_key}")

    n_samples = n_samples_in_file(spec)
    n_chunks = int(np.ceil(n_samples / float(chunk_size)))
    if max_chunks is not None:
        n_chunks_to_export = min(n_chunks, max_chunks)
    else:
        n_chunks_to_export = n_chunks

    out_base = spec.data_file.stem
    n_dict = int(np.asarray(shared["N_dict"]).reshape(-1)[0])
    complex_dtype = np.complex64 if export_precision == "single" else np.complex128
    chunk_shared = {
        key: value
        for key, value in shared.items()
        if key not in CHUNK_EXCLUDED_SHARED_FIELDS
    }

    manifest = {
        "run_key": run_key,
        "dataset": dataset,
        "run_name": run_name,
        "checkpoint_path": checkpoint_path,
        "source_summary": str(spec.output_summary),
        "data_file": str(spec.data_file),
        "output_dir": str(output_dir),
        "output_dir_windows": wsl_to_win(output_dir),
        "chunk_size": int(chunk_size),
        "export_precision": export_precision,
        "n_samples": int(n_samples),
        "n_chunks_total": int(n_chunks),
        "n_chunks_exported": 0,
        "complete": False,
        "status": "running",
    }
    (output_dir / "export_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

    for chunk_idx in range(n_chunks_to_export):
        start = chunk_idx * chunk_size
        end = min((chunk_idx + 1) * chunk_size, n_samples)
        out_file = output_dir / (
            f"{out_base}_Python_resdmd_Layer_100_Ndict_{n_dict}_outputs_{chunk_idx + 1}.mat"
        )
        if out_file.exists() and not overwrite:
            print(f"[{run_key}] skip existing chunk {chunk_idx + 1}/{n_chunks}", flush=True)
            continue
        rows = raw_vs_std.load_rows_by_index(
            spec,
            np.arange(start, end, dtype=np.int64),
            normalization=normalizer,
        )
        efuns = raw_vs_std.compute_efuns(solver, rows).astype(complex_dtype, copy=False)
        scipy.io.savemat(
            out_file,
            {
                "EDMD_outputs": {
                    "efuns": efuns,
                    **chunk_shared,
                },
                "Data_normalization": normalization,
            },
            long_field_names=True,
        )
        print(
            f"[{run_key}] saved chunk {chunk_idx + 1}/{n_chunks} "
            f"samples {start + 1}:{end} -> {out_file.name}",
            flush=True,
        )
        del rows, efuns
        gc.collect()

    n_chunks_exported = len(list(output_dir.glob("*_outputs_*.mat")))
    manifest["complete"] = bool(n_chunks_to_export == n_chunks and n_chunks_exported >= n_chunks)
    manifest["n_chunks_exported"] = int(n_chunks_exported)
    manifest["status"] = "complete" if manifest["complete"] else "partial"
    (output_dir / "export_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return manifest


def main() -> None:
    args = parse_args()
    run_keys = RUN_KEYS if args.run == "both" else (args.run,)
    manifests = []
    for run_key in run_keys:
        manifests.append(
            export_one(
                run_key=run_key,
                output_root=args.output_root,
                chunk_size=args.chunk_size,
                export_precision=args.export_precision,
                overwrite=args.overwrite,
                max_chunks=args.max_chunks,
            )
        )
    print(json.dumps(manifests, indent=2), flush=True)


if __name__ == "__main__":
    main()
