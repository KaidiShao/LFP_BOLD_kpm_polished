#!/usr/bin/env python3
"""Export full e10gh1 standardized BLP EDMD outputs from best checkpoints.

This is the E10gH1 counterpart of ``export_e10gb1_standardized_edmd_outputs``.
It reuses the same streaming export implementation but swaps in the E10gH1
raw-vs-standardized monitor module, where the accepted standardized complex
split run is defined.
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TMP_DIR = REPO_ROOT / "tmp"
if str(TMP_DIR) not in sys.path:
    sys.path.insert(0, str(TMP_DIR))

import export_e10gb1_standardized_edmd_outputs as exporter
import plot_e10gh1_blp_raw_vs_std_spectra_loss_efun_20260524 as e10gh1_specs


DEFAULT_OUTPUT_ROOT = Path("/mnt/e/DataPons_processed/derived_autodl_results_standardize")
RUN_KEYS = ("abs_std", "complex_std")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", choices=RUN_KEYS, default="complex_std")
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--chunk-size", type=int, default=5000)
    parser.add_argument("--export-precision", choices=("single", "double"), default="single")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-chunks", type=int, default=None)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    exporter.raw_vs_std = e10gh1_specs
    manifest = exporter.export_one(
        run_key=args.run,
        output_root=args.output_root,
        chunk_size=args.chunk_size,
        export_precision=args.export_precision,
        overwrite=args.overwrite,
        max_chunks=args.max_chunks,
    )
    wrong_dir = Path(manifest["output_dir"])
    correct_dir = args.output_root / "e10gh1" / "mlp" / "outputs" / wrong_dir.name
    if wrong_dir != correct_dir:
        correct_dir.parent.mkdir(parents=True, exist_ok=True)
        if correct_dir.exists():
            if args.overwrite:
                shutil.rmtree(correct_dir)
            else:
                raise FileExistsError(
                    f"Correct E10gH1 output directory already exists: {correct_dir}"
                )
        shutil.move(str(wrong_dir), str(correct_dir))
        manifest["output_dir"] = str(correct_dir)
        manifest["output_dir_windows"] = exporter.wsl_to_win(correct_dir)
        (correct_dir / "export_manifest.json").write_text(
            json.dumps(manifest, indent=2), encoding="utf-8"
        )
    print(json.dumps([manifest], indent=2), flush=True)


if __name__ == "__main__":
    main()
