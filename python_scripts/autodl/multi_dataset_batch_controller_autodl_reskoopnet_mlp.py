import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


DATASET_BATCH_PATH = Path(__file__).with_name("dataset_batch_controller_autodl_reskoopnet_mlp.py")


def sanitize_tag(value):
    return str(value).replace("/", "-").replace(" ", "_")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run multiple datasets sequentially. Each dataset delegates to the "
            "dataset-level batch controller, which then runs observable batches and residual-form batches."
        )
    )
    parser.add_argument("--ssh-host", required=True)
    parser.add_argument("--ssh-user", default="root")
    parser.add_argument("--ssh-port", type=int, default=22)
    parser.add_argument("--ssh-key", default=None)

    parser.add_argument("--remote-code-dir", default="/root/autodl-tmp/code/reskoopnet_mlp")
    parser.add_argument("--remote-data-root", default="/root/autodl-tmp/data")
    parser.add_argument("--remote-output-parent", default="/root/autodl-tmp/outputs")
    parser.add_argument("--remote-checkpoint-parent", default="/root/autodl-tmp/checkpoints")
    parser.add_argument("--remote-log-parent", default="/root/autodl-tmp/logs")
    parser.add_argument("--remote-env-setup", default="")
    parser.add_argument("--remote-python", default="python")
    parser.add_argument("--remote-script-name", default="run_autodl_reskoopnet_mlp.py")

    parser.add_argument("--dataset-stems", nargs="+", required=True)
    parser.add_argument("--experiment-name", required=True)
    parser.add_argument(
        "--observable-modes",
        nargs="+",
        default=["abs", "complex_split"],
        choices=[
            "abs", "complex", "complex_split", "eleHP", "HP", "identity",
            "roi_mean", "slow_band_power", "svd", "HP_svd100", "global_svd100",
        ],
    )
    parser.add_argument("--residual-forms", nargs="+", default=["projected_kv", "projected_vlambda"])
    parser.add_argument(
        "--loss-mode",
        default="squared",
        choices=["squared", "per_dim", "relative_target"],
        help="Legacy compatibility flag; optimization now always uses squared residual loss.",
    )
    parser.add_argument(
        "--loss-epsilon",
        type=float,
        default=1e-12,
        help="Epsilon used only for relative_target diagnostic metrics.",
    )
    parser.add_argument(
        "--remote-data-subdir-prefix",
        default=None,
        help=(
            "Optional prefix used to build a per-dataset remote input directory. "
            "If omitted, each dataset uses '<dataset_stem>_batch'."
        ),
    )

    parser.add_argument("--run-name-base", default="mlp_obs")
    parser.add_argument("--selected-device", default="gpu", choices=["cpu", "gpu"])
    parser.add_argument("--solver-name", default="resdmd_batch")
    parser.add_argument("--file-type", default=".h5", choices=[".h5", ".mat"])
    parser.add_argument("--field-name", default="obs")
    parser.add_argument("--layer-sizes", type=int, nargs="+", default=[100, 100, 100])
    parser.add_argument("--n-psi-train", type=int, default=100)
    parser.add_argument("--train-ratio", type=float, default=0.7)
    parser.add_argument("--reg", type=float, default=0.1)
    parser.add_argument("--rounds", type=int, default=1)
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--batch-size", type=int, default=2000)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--log-interval", type=int, default=1)
    parser.add_argument("--lr-decay-factor", type=float, default=0.8)
    parser.add_argument("--inner-epochs", type=int, default=5)
    parser.add_argument("--end-condition", type=float, default=1e-9)
    parser.add_argument("--chunk-size", type=int, default=5000)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--fresh-checkpoints", action="store_true")
    parser.set_defaults(export_psi=False)
    parser.add_argument("--export-psi", dest="export_psi", action="store_true")
    parser.add_argument("--use-smoke-inputs", action="store_true")
    parser.add_argument("--smoke-input-root", default=str(Path(__file__).with_name("test_inputs")))
    parser.add_argument("--smoke-pattern-tag", default="smoke20k_v73")
    parser.add_argument(
        "--local-download-root-base",
        default=None,
        help=(
            "Optional base directory for downloads. Each dataset will use "
            "<base>/<dataset_stem>/mlp as its local download root."
        ),
    )
    parser.add_argument("--download-checkpoints", action="store_true")
    parser.add_argument("--skip-data-upload", action="store_true")
    parser.add_argument("--delete-remote-run-after-download", action="store_true")
    parser.add_argument("--delete-remote-input-after-observable", action="store_true")
    parser.add_argument("--poll-interval", type=int, default=30)
    parser.add_argument("--tail-lines", type=int, default=20)
    parser.add_argument("--launch-timeout-sec", type=int, default=30)
    parser.add_argument("--continue-on-error", action="store_true")
    return parser.parse_args()


def build_dataset_batch_cmd(args, dataset_stem):
    cmd = [
        sys.executable,
        str(DATASET_BATCH_PATH),
        "--ssh-host",
        args.ssh_host,
        "--ssh-user",
        args.ssh_user,
        "--ssh-port",
        str(args.ssh_port),
        "--remote-code-dir",
        args.remote_code_dir,
        "--remote-data-root",
        args.remote_data_root,
        "--remote-output-parent",
        args.remote_output_parent,
        "--remote-checkpoint-parent",
        args.remote_checkpoint_parent,
        "--remote-log-parent",
        args.remote_log_parent,
        "--remote-python",
        args.remote_python,
        "--remote-script-name",
        args.remote_script_name,
        "--dataset-stem",
        dataset_stem,
        "--experiment-name",
        args.experiment_name,
        "--run-name-base",
        args.run_name_base,
        "--selected-device",
        args.selected_device,
        "--solver-name",
        args.solver_name,
        "--loss-mode",
        args.loss_mode,
        "--loss-epsilon",
        str(args.loss_epsilon),
        "--file-type",
        args.file_type,
        "--field-name",
        args.field_name,
        "--n-psi-train",
        str(args.n_psi_train),
        "--train-ratio",
        str(args.train_ratio),
        "--reg",
        str(args.reg),
        "--rounds",
        str(args.rounds),
        "--epochs",
        str(args.epochs),
        "--batch-size",
        str(args.batch_size),
        "--lr",
        str(args.lr),
        "--log-interval",
        str(args.log_interval),
        "--lr-decay-factor",
        str(args.lr_decay_factor),
        "--inner-epochs",
        str(args.inner_epochs),
        "--end-condition",
        str(args.end_condition),
        "--chunk-size",
        str(args.chunk_size),
        "--poll-interval",
        str(args.poll_interval),
        "--tail-lines",
        str(args.tail_lines),
        "--launch-timeout-sec",
        str(args.launch_timeout_sec),
        "--observable-modes",
        *args.observable_modes,
        "--residual-forms",
        *args.residual_forms,
    ]

    if args.ssh_key:
        cmd.extend(["--ssh-key", args.ssh_key])
    if args.remote_env_setup:
        cmd.extend(["--remote-env-setup", args.remote_env_setup])
    if args.use_smoke_inputs:
        cmd.append("--use-smoke-inputs")
        cmd.extend(["--smoke-input-root", args.smoke_input_root])
        cmd.extend(["--smoke-pattern-tag", args.smoke_pattern_tag])
    if args.local_download_root_base:
        local_root = Path(args.local_download_root_base).expanduser().resolve() / dataset_stem / "mlp"
        cmd.extend(["--local-download-root", str(local_root)])

    remote_subdir = (
        f"{args.remote_data_subdir_prefix}_{dataset_stem}"
        if args.remote_data_subdir_prefix
        else f"{dataset_stem}_batch"
    )
    cmd.extend(["--remote-data-subdir", remote_subdir])

    for layer_size in args.layer_sizes:
        if "--layer-sizes" not in cmd:
            cmd.extend(["--layer-sizes", str(layer_size)])
        else:
            cmd.append(str(layer_size))

    if args.resume:
        cmd.append("--resume")
    if args.fresh_checkpoints:
        cmd.append("--fresh-checkpoints")
    if args.export_psi:
        cmd.append("--export-psi")
    if args.download_checkpoints:
        cmd.append("--download-checkpoints")
    if args.skip_data_upload:
        cmd.append("--skip-data-upload")
    if args.delete_remote_run_after_download:
        cmd.append("--delete-remote-run-after-download")
    if args.delete_remote_input_after_observable:
        cmd.append("--delete-remote-input-after-observable")
    if args.continue_on_error:
        cmd.append("--continue-on-error")

    return cmd


def main():
    args = parse_args()

    manifest_dir = (
        Path(args.local_download_root_base).expanduser().resolve()
        if args.local_download_root_base
        else Path(__file__).resolve().parent
    )
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = manifest_dir / f"{args.experiment_name}_multi_dataset_manifest.json"

    manifest = {
        "experiment_name": args.experiment_name,
        "dataset_stems": list(args.dataset_stems),
        "observable_modes": list(args.observable_modes),
        "residual_forms": list(args.residual_forms),
        "loss_mode": "squared",
        "requested_loss_mode": args.loss_mode,
        "loss_epsilon": args.loss_epsilon,
        "use_smoke_inputs": bool(args.use_smoke_inputs),
        "smoke_pattern_tag": args.smoke_pattern_tag,
        "started_at_unix": time.time(),
        "datasets": [],
    }

    for dataset_stem in args.dataset_stems:
        cmd = build_dataset_batch_cmd(args, dataset_stem)

        print("")
        print("=" * 80)
        print(f"[multi-dataset] starting dataset_stem={dataset_stem}")
        print("=" * 80)

        started_at = time.time()
        result = subprocess.run(cmd, check=False)
        duration_sec = time.time() - started_at

        manifest["datasets"].append(
            {
                "dataset_stem": dataset_stem,
                "loss_mode": "squared",
                "requested_loss_mode": args.loss_mode,
                "return_code": int(result.returncode),
                "duration_sec": duration_sec,
            }
        )
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

        if result.returncode != 0:
            print(f"[multi-dataset] dataset_stem={dataset_stem} failed with return code {result.returncode}")
            if not args.continue_on_error:
                raise SystemExit(result.returncode)
        else:
            print(f"[multi-dataset] dataset_stem={dataset_stem} finished successfully in {duration_sec:.1f}s")

    manifest["finished_at_unix"] = time.time()
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("")
    print(f"[multi-dataset] finished. Manifest saved to: {manifest_path}")


if __name__ == "__main__":
    main()
