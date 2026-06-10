#!/usr/bin/env python3
"""Export full standardized complex-split BLP EDMD chunks for P5.

The current standardized complex-split P4 runs save compact summary MAT files
without ``EDMD_outputs.efuns``.  P5 needs normal ``*_outputs_*.mat`` chunks, so
this wrapper reuses the existing e10gb1 exporter machinery and supplies
dataset-specific paths inferred from the mainline stdComplexPair run name.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import scipy.io


REPO_ROOT = Path(__file__).resolve().parents[1]
TMP_DIR = REPO_ROOT / "tmp"
if str(TMP_DIR) not in sys.path:
    sys.path.insert(0, str(TMP_DIR))

import export_e10gb1_standardized_edmd_outputs as exporter
import plot_e10gb1_blp_raw_vs_std_spectra_loss_efun_20260523 as template


DEFAULT_OUTPUT_ROOT = Path("/mnt/e/DataPons_processed/derived_autodl_results_standardize")
DEFAULT_AUTODL_ROOT = Path("/mnt/e/autodl_results_new")
DEFAULT_PROCESSED_ROOT = Path("/mnt/e/DataPons_processed")
RUN_STEM_TEMPLATE = (
    "mlp_obs_blp_vlambda_complex_split_stdComplexPair_l1e4_r1e3_b2000_i2_pat40_"
    "20260522_pat40_allblp_{dataset}_seed1234_projected_vlambda_complex_split"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="+", required=True)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--autodl-root", type=Path, default=DEFAULT_AUTODL_ROOT)
    parser.add_argument("--processed-root", type=Path, default=DEFAULT_PROCESSED_ROOT)
    parser.add_argument("--chunk-size", type=int, default=5000)
    parser.add_argument("--export-precision", choices=("single", "double"), default="single")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-chunks", type=int, default=None)
    return parser.parse_args()


def latest_summary(output_dir: Path) -> Path:
    files = sorted(output_dir.glob("*_summary.mat"), key=lambda p: p.stat().st_mtime)
    if not files:
        raise FileNotFoundError(f"No summary MAT found in {output_dir}")
    return files[-1]


def infer_n_dict(summary_file: Path) -> int:
    mat = scipy.io.loadmat(summary_file, simplify_cells=True, variable_names=["EDMD_outputs"])
    edmd = mat["EDMD_outputs"]
    if "N_dict" not in edmd:
        raise KeyError(f"N_dict missing from {summary_file}")
    return int(edmd["N_dict"])


def build_spec(dataset: str, autodl_root: Path, processed_root: Path):
    run_name = RUN_STEM_TEMPLATE.format(dataset=dataset)
    output_dir = autodl_root / dataset / "mlp" / "outputs" / run_name
    summary_file = latest_summary(output_dir)
    checkpoint_dir = autodl_root / dataset / "mlp" / "checkpoints" / run_name
    state_file = checkpoint_dir / "final" / "training_state.json"
    best_dir = checkpoint_dir / "best"
    if not state_file.is_file():
        raise FileNotFoundError(f"Missing training state: {state_file}")
    if not best_dir.is_dir():
        raise FileNotFoundError(f"Missing best checkpoint dir: {best_dir}")

    n_dict = infer_n_dict(summary_file)
    target_dim = n_dict - template.N_PSI_TRAIN - 1
    if target_dim <= 0:
        raise ValueError(f"Invalid target_dim={target_dim} inferred from N_dict={n_dict}")

    data_file = (
        processed_root
        / dataset
        / "pipeline1_reskoopnet_dictionary"
        / f"{dataset}_low50_high250_g2_complex_split_single.mat"
    )
    if not data_file.is_file():
        raise FileNotFoundError(f"Missing P1 complex-split data file: {data_file}")

    return template.RunSpec(
        key="complex_std",
        observable="complex_split",
        condition="std_complex_pair",
        label=f"{dataset} complex_split std_complex_pair",
        data_file=data_file,
        checkpoint_dir=checkpoint_dir,
        output_summary=summary_file,
        solver_name="resdmd_batch_mixedgpu",
        training_policy="float32",
        analysis_dtype="float64",
        gram_dtype="float64",
        spectral_dtype="float64",
        target_dim=target_dim,
        default_reg=0.001,
        standardize=True,
    )


def main() -> int:
    args = parse_args()
    manifests = []
    for dataset in args.datasets:
        spec = build_spec(dataset, args.autodl_root, args.processed_root)
        template.RUNS = {"complex_std": spec}
        exporter.raw_vs_std = template
        manifest = exporter.export_one(
            run_key="complex_std",
            output_root=args.output_root,
            chunk_size=args.chunk_size,
            export_precision=args.export_precision,
            overwrite=args.overwrite,
            max_chunks=args.max_chunks,
        )
        manifests.append(manifest)
    print(json.dumps(manifests, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
