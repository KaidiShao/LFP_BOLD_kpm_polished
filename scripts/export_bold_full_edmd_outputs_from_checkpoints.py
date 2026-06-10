#!/usr/bin/env python3
"""Export full BOLD EDMD output chunks from summary-only BOLD checkpoints.

Some current K13 BOLD P4 runs were saved as compact ``*_summary.mat`` files
without ``EDMD_outputs.efuns``.  Pipeline 7/9/8/10 need the normal chunked
``*_outputs_*.mat`` files.  This script restores the existing best checkpoints
and writes full chunks into a derived AutoDL-style root, leaving the original
AutoDL output folders untouched.
"""

from __future__ import annotations

import argparse
import gc
import json
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import scipy.io


REPO_ROOT = Path(__file__).resolve().parents[1]
TMP_DIR = REPO_ROOT / "tmp"
if str(TMP_DIR) not in sys.path:
    sys.path.insert(0, str(TMP_DIR))

import plot_e10gb1_blp_raw_vs_std_spectra_loss_efun_20260523 as restore_helpers


DEFAULT_SOURCE_ROOT = Path("/mnt/e/autodl_results_local/bold_wsl")
DEFAULT_PROCESSED_ROOT = Path("/mnt/e/DataPons_processed")
DEFAULT_OUTPUT_ROOT = Path("/mnt/e/DataPons_processed/derived_autodl_results_bold_full_export")
DATASET_KEY = "/obs"
N_PSI_TRAIN = 100

CHUNK_EXCLUDED_SHARED_FIELDS = {
    "snapshot_valid_idx_x",
    "snapshot_valid_idx_y",
    "snapshot_session_idx",
    "snapshot_train_pair_idx",
    "snapshot_valid_pair_idx",
}


@dataclass(frozen=True)
class RunSpec:
    dataset: str
    run_name: str
    observable_mode: str
    data_file: Path
    checkpoint_dir: Path
    output_summary: Path
    target_dim: int
    default_reg: float
    standardize: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="+", required=True)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--run-names", nargs="*", default=None)
    parser.add_argument("--only-current-best", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=5000)
    parser.add_argument("--export-precision", choices=("single", "double"), default="single")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def read_json(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


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


def infer_n_dict(summary_path: Path) -> int:
    shared, _ = load_summary_payload(summary_path)
    return int(np.asarray(shared["N_dict"]).reshape(-1)[0])


def observable_file(processed_root: Path, dataset: str, observable_mode: str) -> Path:
    obs_dir = processed_root / dataset / "pipeline3_bold_observables"
    files = sorted(obs_dir.glob(f"*_bold_observables_{observable_mode}.mat"))
    if not files:
        raise FileNotFoundError(f"No P3 BOLD observable file for {dataset} {observable_mode}: {obs_dir}")
    if len(files) > 1:
        raise RuntimeError(f"Ambiguous P3 BOLD observable files for {dataset} {observable_mode}: {files}")
    return files[0]


def latest_summary(output_dir: Path) -> Path:
    files = sorted(output_dir.glob("*_summary.mat"), key=lambda p: p.stat().st_mtime)
    if not files:
        raise FileNotFoundError(f"No summary MAT found in {output_dir}")
    return files[-1]


def discover_specs(args: argparse.Namespace) -> list[RunSpec]:
    specs: list[RunSpec] = []
    run_filter = {name.lower() for name in args.run_names or []}
    for dataset in args.datasets:
        output_parent = args.source_root / dataset / "mlp" / "outputs"
        checkpoint_parent = args.source_root / dataset / "mlp" / "checkpoints"
        if not output_parent.is_dir():
            continue
        for run_dir in sorted(p for p in output_parent.iterdir() if p.is_dir()):
            run_name = run_dir.name
            if run_filter and run_name.lower() not in run_filter:
                continue
            checkpoint_dir = checkpoint_parent / run_name
            state_file = checkpoint_dir / "final" / "training_state.json"
            best_dir = checkpoint_dir / "best"
            if not state_file.is_file() or not best_dir.is_dir():
                continue
            state = read_json(state_file)
            meta = state.get("run_metadata") or {}
            observable_mode = str(meta.get("observable_mode") or "")
            if not observable_mode:
                continue
            summary = latest_summary(run_dir)
            n_dict = infer_n_dict(summary)
            target_dim = n_dict - N_PSI_TRAIN - 1
            if target_dim <= 0:
                raise ValueError(f"Invalid target_dim={target_dim} for {summary}")
            specs.append(
                RunSpec(
                    dataset=dataset,
                    run_name=run_name,
                    observable_mode=observable_mode,
                    data_file=observable_file(args.processed_root, dataset, observable_mode),
                    checkpoint_dir=checkpoint_dir,
                    output_summary=summary,
                    target_dim=target_dim,
                    default_reg=float(state.get("best_reg", state.get("current_reg", 0.001))),
                    standardize=bool(meta.get("standardize_data", False)),
                )
            )
    if args.only_current_best:
        specs = keep_current_best(specs)
    return specs


def keep_current_best(specs: list[RunSpec]) -> list[RunSpec]:
    grouped: dict[tuple[str, str], list[RunSpec]] = {}
    for spec in specs:
        grouped.setdefault((spec.dataset.lower(), spec.observable_mode.lower()), []).append(spec)
    keep: list[RunSpec] = []
    for _, group in sorted(grouped.items()):
        best = min(group, key=lambda s: best_val_metric(s))
        keep.append(best)
    return sorted(keep, key=lambda s: (s.dataset.lower(), s.observable_mode.lower(), s.run_name.lower()))


def best_val_metric(spec: RunSpec) -> float:
    state = read_json(spec.checkpoint_dir / "final" / "training_state.json")
    try:
        return float(state.get("best_val_metric", np.inf))
    except (TypeError, ValueError):
        return float("inf")


def n_samples_in_file(path: Path) -> int:
    with h5py.File(path, "r") as handle:
        return int(handle[DATASET_KEY].shape[1])


def load_rows(path: Path, row_idx: np.ndarray, normalization: dict | None) -> np.ndarray:
    row_idx = np.asarray(row_idx, dtype=np.int64)
    order = np.argsort(row_idx)
    inv = np.empty_like(order)
    inv[order] = np.arange(order.size)
    with h5py.File(path, "r") as handle:
        rows = np.asarray(handle[DATASET_KEY][:, row_idx[order]], dtype=np.float64).T
    rows = rows[inv]
    if normalization is not None:
        mean = np.asarray(normalization["mean"], dtype=np.float64).reshape(1, -1)
        scale = np.asarray(normalization["scale"], dtype=np.float64).reshape(1, -1)
        rows = (rows - mean) / scale
    return rows


def make_restore_spec(spec: RunSpec) -> restore_helpers.RunSpec:
    return restore_helpers.RunSpec(
        key=spec.run_name,
        observable=spec.observable_mode,
        condition="standardized" if spec.standardize else "raw",
        label=f"{spec.dataset} {spec.observable_mode}",
        data_file=spec.data_file,
        checkpoint_dir=spec.checkpoint_dir,
        output_summary=spec.output_summary,
        solver_name="resdmd_batch4",
        training_policy="float64",
        analysis_dtype="float64",
        gram_dtype="float64",
        spectral_dtype="float64",
        target_dim=spec.target_dim,
        default_reg=spec.default_reg,
        standardize=spec.standardize,
    )


def clear_existing(output_dir: Path) -> None:
    if not output_dir.exists():
        return
    for path in output_dir.glob("*_outputs_*.mat"):
        path.unlink()
    for path in output_dir.glob("*_summary.mat"):
        path.unlink()
    for name in ("training_loss_diagnostics.png", "training_inner_loss_metrics.png", "export_manifest.json"):
        path = output_dir / name
        if path.exists():
            path.unlink()


def export_one(spec: RunSpec, args: argparse.Namespace) -> dict:
    output_dir = args.output_root / spec.dataset / "mlp" / "outputs" / spec.run_name
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.overwrite:
        clear_existing(output_dir)
    existing = sorted(output_dir.glob("*_outputs_*.mat"))
    manifest_file = output_dir / "export_manifest.json"
    if existing and not args.overwrite:
        old_manifest = {}
        if manifest_file.exists():
            try:
                old_manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                old_manifest = {}
        if old_manifest.get("status") == "complete":
            return {
                "dataset": spec.dataset,
                "run_name": spec.run_name,
                "status": "exists",
                "n_existing_chunks": len(existing),
                "output_dir": str(output_dir),
            }
        print(
            f"[{spec.dataset}] {spec.run_name}: found {len(existing)} existing chunk(s); "
            "resuming missing chunks.",
            flush=True,
        )

    shared, normalization = load_summary_payload(spec.output_summary)
    if spec.standardize and normalization is None:
        raise RuntimeError(f"Standardized run is missing Data_normalization: {spec.output_summary}")

    summary_copy = output_dir / spec.output_summary.name
    scipy.io.savemat(
        summary_copy,
        {"EDMD_outputs": shared, **({"Data_normalization": normalization} if normalization is not None else {})},
        long_field_names=True,
    )
    for png in ("training_loss_diagnostics.png", "training_inner_loss_metrics.png"):
        source = spec.output_summary.parent / png
        if source.exists():
            shutil.copy2(source, output_dir / png)

    restore_spec = make_restore_spec(spec)
    solver, checkpoint_path = restore_helpers.restore_solver(restore_spec)
    n_samples = n_samples_in_file(spec.data_file)
    n_chunks = int(np.ceil(n_samples / float(args.chunk_size)))
    complex_dtype = np.complex64 if args.export_precision == "single" else np.complex128
    chunk_shared = {key: value for key, value in shared.items() if key not in CHUNK_EXCLUDED_SHARED_FIELDS}
    base = spec.data_file.stem
    n_dict = int(np.asarray(shared["N_dict"]).reshape(-1)[0])

    manifest = {
        "dataset": spec.dataset,
        "run_name": spec.run_name,
        "observable_mode": spec.observable_mode,
        "source_summary": str(spec.output_summary),
        "checkpoint_path": checkpoint_path,
        "data_file": str(spec.data_file),
        "output_dir": str(output_dir),
        "chunk_size": int(args.chunk_size),
        "export_precision": args.export_precision,
        "n_samples": int(n_samples),
        "n_chunks": int(n_chunks),
        "status": "running",
    }
    (output_dir / "export_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    for i in range(n_chunks):
        start = i * args.chunk_size
        end = min((i + 1) * args.chunk_size, n_samples)
        out_file = output_dir / f"{base}_Python_resdmd_Layer_100_Ndict_{n_dict}_outputs_{i + 1}.mat"
        if out_file.exists() and not args.overwrite:
            print(f"[{spec.dataset}] {spec.run_name}: skip existing chunk {i + 1}/{n_chunks}", flush=True)
            continue
        rows = load_rows(spec.data_file, np.arange(start, end, dtype=np.int64), normalization)
        efuns = restore_helpers.compute_efuns(solver, rows).astype(complex_dtype, copy=False)
        payload = {"EDMD_outputs": {"efuns": efuns, **chunk_shared}}
        if normalization is not None:
            payload["Data_normalization"] = normalization
        scipy.io.savemat(out_file, payload, long_field_names=True)
        print(f"[{spec.dataset}] {spec.run_name}: saved chunk {i + 1}/{n_chunks}", flush=True)
        del rows, efuns, payload
        gc.collect()

    manifest["n_chunks_exported"] = len(list(output_dir.glob("*_outputs_*.mat")))
    manifest["status"] = "complete" if manifest["n_chunks_exported"] >= n_chunks else "partial"
    (output_dir / "export_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def main() -> int:
    args = parse_args()
    specs = discover_specs(args)
    if not specs:
        print("No matching runs found.", flush=True)
        return 1
    print(f"Runs to export: {len(specs)}", flush=True)
    for spec in specs:
        print(f"  {spec.dataset} | {spec.observable_mode} | {spec.run_name}", flush=True)
    manifests = []
    for spec in specs:
        manifests.append(export_one(spec, args))
    print(json.dumps(manifests, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
