import argparse
import shutil
from pathlib import Path

import numpy as np
import scipy.io


DEFAULT_EXCLUDED_FIELDS = (
    "snapshot_valid_idx_x",
    "snapshot_valid_idx_y",
    "snapshot_session_idx",
    "snapshot_train_pair_idx",
    "snapshot_valid_pair_idx",
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Rewrite EDMD output chunk .mat files without large per-run snapshot "
            "metadata that should live only in the summary file."
        )
    )
    parser.add_argument("--source-dir", required=True, help="Directory containing original EDMD .mat outputs.")
    parser.add_argument("--target-dir", required=True, help="Directory to write compacted outputs into.")
    parser.add_argument(
        "--excluded-field",
        action="append",
        default=list(DEFAULT_EXCLUDED_FIELDS),
        help="EDMD_outputs struct field to omit from chunk files. Can be repeated.",
    )
    parser.add_argument("--overwrite-target", action="store_true", help="Delete target dir before writing.")
    parser.add_argument("--progress-interval", type=int, default=25)
    return parser.parse_args()


def load_edmd_struct(path):
    mat = scipy.io.loadmat(path, squeeze_me=False, struct_as_record=False)
    if "EDMD_outputs" not in mat:
        raise KeyError(f"Missing EDMD_outputs in {path}")
    return mat["EDMD_outputs"][0, 0]


def struct_to_dict(edmd_struct, excluded_fields=frozenset()):
    fields = getattr(edmd_struct, "_fieldnames", None)
    if not fields:
        raise ValueError("EDMD_outputs is not a MATLAB struct with fields")
    return {
        field: getattr(edmd_struct, field)
        for field in fields
        if field not in excluded_fields
    }


def is_output_chunk(path):
    name = path.name
    return "_outputs_" in name and "_outputs_Psi_" not in name and name.endswith(".mat")


def verify_no_excluded_fields(path, excluded_fields):
    edmd_struct = load_edmd_struct(path)
    fields = set(getattr(edmd_struct, "_fieldnames", []) or [])
    present = sorted(fields.intersection(excluded_fields))
    if present:
        raise ValueError(f"{path} still contains excluded fields: {present}")
    if "efuns" not in fields:
        raise ValueError(f"{path} missing efuns after compaction")
    efuns = np.asarray(getattr(edmd_struct, "efuns"))
    if efuns.size == 0:
        raise ValueError(f"{path} has empty efuns after compaction")


def main():
    args = parse_args()
    source_dir = Path(args.source_dir)
    target_dir = Path(args.target_dir)
    excluded_fields = set(args.excluded_field)

    if not source_dir.is_dir():
        raise FileNotFoundError(f"Source directory not found: {source_dir}")

    if target_dir.exists() and args.overwrite_target:
        shutil.rmtree(target_dir)
    target_dir.mkdir(parents=True, exist_ok=True)

    all_files = sorted(path for path in source_dir.iterdir() if path.is_file())
    output_files = [path for path in all_files if is_output_chunk(path)]
    passthrough_files = [path for path in all_files if not is_output_chunk(path)]

    print(f"source_dir={source_dir}")
    print(f"target_dir={target_dir}")
    print(f"output_chunks={len(output_files)} passthrough_files={len(passthrough_files)}")
    print(f"excluded_fields={sorted(excluded_fields)}")

    for path in passthrough_files:
        shutil.copy2(path, target_dir / path.name)

    for idx, path in enumerate(output_files, start=1):
        target_path = target_dir / path.name
        temp_path = target_dir / f"{path.name}.tmp"
        edmd_struct = load_edmd_struct(path)
        compact_dict = struct_to_dict(edmd_struct, excluded_fields)
        scipy.io.savemat(temp_path, {"EDMD_outputs": compact_dict})
        verify_no_excluded_fields(temp_path, excluded_fields)
        temp_path.replace(target_path)

        if idx == 1 or idx == len(output_files) or idx % args.progress_interval == 0:
            src_mib = path.stat().st_size / (1024 ** 2)
            dst_mib = target_path.stat().st_size / (1024 ** 2)
            print(f"compacted {idx}/{len(output_files)}: {path.name} {src_mib:.1f} MiB -> {dst_mib:.1f} MiB")

    summary_count = len(list(target_dir.glob("*_summary.mat")))
    output_count = len([path for path in target_dir.iterdir() if path.is_file() and is_output_chunk(path)])
    png_count = len(list(target_dir.glob("training_loss_diagnostics.png")))
    pdf_count = len(list(target_dir.glob("training_loss_diagnostics.pdf")))

    if summary_count != len(list(source_dir.glob("*_summary.mat"))):
        raise RuntimeError("Summary file count mismatch after compaction")
    if output_count != len(output_files):
        raise RuntimeError("Output chunk count mismatch after compaction")

    total_gib = sum(path.stat().st_size for path in target_dir.iterdir() if path.is_file()) / (1024 ** 3)
    print(
        "COMPACT_OK "
        f"summary={summary_count} outputs={output_count} png={png_count} pdf={pdf_count} "
        f"target_gib={total_gib:.3f}"
    )


if __name__ == "__main__":
    main()
