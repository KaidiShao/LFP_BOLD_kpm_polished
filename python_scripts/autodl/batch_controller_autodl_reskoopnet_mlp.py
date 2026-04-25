import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


CONTROLLER_PATH = Path(__file__).with_name("controller_autodl_reskoopnet_mlp.py")
OBSERVABLE_TAGS = {
    "abs": "abs",
    "complex": "complex_split",
    "complex_split": "complex_split",
    "eleHP": "eleHP",
    "HP": "HP",
    "identity": "identity",
    "roi_mean": "roi_mean",
    "slow_band_power": "slow_band_power",
    "svd": "svd",
    "HP_svd100": "HP_svd100",
    "global_svd100": "global_svd100",
}


def sanitize_tag(value):
    return str(value).replace("/", "-").replace(" ", "_")


def build_run_label(run_name_base, experiment_name, residual_form, observable_mode, loss_mode="squared"):
    observable_tag = OBSERVABLE_TAGS[observable_mode]
    effective_experiment_name = f"{experiment_name}_{sanitize_tag(residual_form)}"
    return f"{run_name_base}_{effective_experiment_name}_{observable_tag}"


def local_output_complete(
    local_download_root,
    run_name_base,
    experiment_name,
    residual_form,
    observable_mode,
    loss_mode="squared",
):
    if not local_download_root:
        return False, None

    run_label = build_run_label(run_name_base, experiment_name, residual_form, observable_mode, loss_mode)
    output_dir = Path(local_download_root).resolve() / "outputs" / run_label
    if not output_dir.exists():
        return False, output_dir

    summary_files = sorted(output_dir.glob("*_summary.mat"))
    output_chunks = sorted(
        path for path in output_dir.glob("*_outputs_*.mat") if "_outputs_Psi_" not in path.name
    )
    return bool(summary_files and output_chunks), output_dir


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run multiple residual_form conditions sequentially for the same observable package. "
            "The first run uploads the input unless --skip-data-upload is given; later runs reuse the same remote input."
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

    parser.add_argument("--local-data-file", required=True)
    parser.add_argument("--local-obs-info-file", default=None)
    parser.add_argument("--remote-data-subdir", default=None)

    parser.add_argument("--run-name-base", default="mlp_obs")
    parser.add_argument("--experiment-name", required=True)
    parser.add_argument("--dataset-stem", default=None)
    parser.add_argument("--selected-device", default="gpu", choices=["cpu", "gpu"])
    parser.add_argument("--solver-name", default="resdmd_batch")
    parser.add_argument(
        "--residual-forms",
        nargs="+",
        default=["projected_kv", "projected_vlambda"],
        choices=["projected_kv", "projected_vlambda"],
    )
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
        "--observable-mode",
        default="abs",
        choices=[
            "abs", "complex", "complex_split", "eleHP", "HP", "identity",
            "roi_mean", "slow_band_power", "svd", "HP_svd100", "global_svd100",
        ],
    )
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
    parser.add_argument("--local-download-root", default=None)
    parser.add_argument("--download-checkpoints", action="store_true")
    parser.add_argument("--skip-data-upload", action="store_true")
    parser.add_argument("--recover-completed-remote-runs", action="store_true")
    parser.add_argument(
        "--skip-completed-local-runs",
        action="store_true",
        help=(
            "Skip a residual_form if its local output directory already contains a summary "
            "and at least one output chunk. Useful when resuming after the local controller crashed."
        ),
    )
    parser.add_argument("--delete-remote-run-after-download", action="store_true")
    parser.add_argument("--poll-interval", type=int, default=15)
    parser.add_argument("--tail-lines", type=int, default=80)
    parser.add_argument("--launch-timeout-sec", type=int, default=30)
    parser.add_argument(
        "--continue-on-error",
        action="store_true",
        help="Continue to the next residual_form even if one run fails.",
    )
    return parser.parse_args()


def build_base_controller_cmd(args):
    cmd = [
        sys.executable,
        str(CONTROLLER_PATH),
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
        "--local-data-file",
        args.local_data_file,
        "--run-name-base",
        args.run_name_base,
        "--experiment-name",
        args.experiment_name,
        "--dataset-stem",
        args.dataset_stem,
        "--selected-device",
        args.selected_device,
        "--solver-name",
        args.solver_name,
        "--loss-mode",
        args.loss_mode,
        "--loss-epsilon",
        str(args.loss_epsilon),
        "--observable-mode",
        args.observable_mode,
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
    ]

    if args.ssh_key:
        cmd.extend(["--ssh-key", args.ssh_key])
    if args.remote_env_setup:
        cmd.extend(["--remote-env-setup", args.remote_env_setup])
    if args.local_obs_info_file:
        cmd.extend(["--local-obs-info-file", args.local_obs_info_file])
    if args.remote_data_subdir:
        cmd.extend(["--remote-data-subdir", args.remote_data_subdir])
    if args.local_download_root:
        cmd.extend(["--local-download-root", args.local_download_root])

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
    if args.recover_completed_remote_runs:
        cmd.append("--recover-completed-remote-run")
    if args.delete_remote_run_after_download:
        cmd.append("--delete-remote-run-after-download")

    return cmd


def main():
    args = parse_args()

    if args.dataset_stem is None:
        raise ValueError("--dataset-stem must be provided for the batch driver.")

    manifest = {
        "experiment_name": args.experiment_name,
        "dataset_stem": args.dataset_stem,
        "observable_mode": args.observable_mode,
        "loss_mode": "squared",
        "requested_loss_mode": args.loss_mode,
        "loss_epsilon": args.loss_epsilon,
        "local_data_file": args.local_data_file,
        "residual_forms": list(args.residual_forms),
        "skip_completed_local_runs": bool(args.skip_completed_local_runs),
        "started_at_unix": time.time(),
        "runs": [],
    }

    manifest_dir = Path(args.local_download_root).resolve() if args.local_download_root else Path(args.local_data_file).resolve().parent
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = (
        manifest_dir
        / f"{args.dataset_stem}_{args.experiment_name}_{args.observable_mode}_batch_manifest.json"
    )

    base_cmd = build_base_controller_cmd(args)
    remote_input_available = bool(args.skip_data_upload)

    for index, residual_form in enumerate(args.residual_forms):
        if args.skip_completed_local_runs:
            is_complete, output_dir = local_output_complete(
                args.local_download_root,
                args.run_name_base,
                args.experiment_name,
                residual_form,
                args.observable_mode,
                args.loss_mode,
            )
            if is_complete:
                print("")
                print("=" * 80)
                print(f"[batch] skipping residual_form={residual_form}; local outputs already complete")
                print(f"[batch] local output dir: {output_dir}")
                print("=" * 80)
                manifest["runs"].append(
                    {
                        "residual_form": residual_form,
                        "loss_mode": "squared",
                        "requested_loss_mode": args.loss_mode,
                        "upload_input": False,
                        "return_code": 0,
                        "duration_sec": 0.0,
                        "skipped_existing_local_output": True,
                        "local_output_dir": str(output_dir),
                    }
                )
                manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
                continue

        cmd = list(base_cmd)
        cmd.extend(["--residual-form", residual_form])

        should_upload_this_run = not remote_input_available
        if not should_upload_this_run:
            cmd.append("--skip-data-upload")

        print("")
        print("=" * 80)
        print(f"[batch] starting residual_form={residual_form} (upload_input={should_upload_this_run})")
        print("=" * 80)

        started_at = time.time()
        result = subprocess.run(cmd, check=False)
        duration_sec = time.time() - started_at
        if result.returncode == 0 and should_upload_this_run:
            remote_input_available = True

        manifest["runs"].append(
            {
                "residual_form": residual_form,
                "loss_mode": "squared",
                "requested_loss_mode": args.loss_mode,
                "upload_input": should_upload_this_run,
                "return_code": int(result.returncode),
                "duration_sec": duration_sec,
                "skipped_existing_local_output": False,
            }
        )
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

        if result.returncode != 0:
            print(f"[batch] residual_form={residual_form} failed with return code {result.returncode}")
            if not args.continue_on_error:
                raise SystemExit(result.returncode)
        else:
            print(f"[batch] residual_form={residual_form} finished successfully in {duration_sec:.1f}s")

    manifest["finished_at_unix"] = time.time()
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("")
    print(f"[batch] finished. Manifest saved to: {manifest_path}")


if __name__ == "__main__":
    main()
