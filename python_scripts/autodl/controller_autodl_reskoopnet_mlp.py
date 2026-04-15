import argparse
import json
import os
import posixpath
import shlex
import subprocess
import sys
import time
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
AUTODL_DIR = Path(__file__).resolve().parent


def resolve_local_results_root():
    env_root = os.environ.get("LFP_BOLD_KPM_RESULTS_ROOT")
    if env_root:
        return Path(env_root)
    if os.name == "nt":
        return Path("E:/autodl_results")
    return Path("/mnt/e/autodl_results")


def sanitize_tag(value):
    return str(value).replace("/", "-").replace(" ", "_")


def infer_dataset_stem_from_observable_file(path_value):
    filename = Path(path_value).name
    for marker in ["_low", "_Python_", "_single"]:
        if marker in filename:
            return filename.split(marker, 1)[0]
    return Path(path_value).stem


def default_local_download_root(dataset_stem, model_family):
    return resolve_local_results_root() / sanitize_tag(dataset_stem) / model_family


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Upload AutoDL MLP pipeline code and input data, launch a remote run, "
            "poll logs, and download the resulting artifacts."
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
    parser.add_argument(
        "--remote-env-setup",
        default="",
        help="Optional shell fragment run before launching Python, e.g. 'source ~/.bashrc && conda activate myenv'",
    )
    parser.add_argument(
        "--remote-python",
        default="python",
        help="Remote Python executable used to launch the training script.",
    )

    parser.add_argument("--local-data-file", required=True)
    parser.add_argument("--local-obs-info-file", default=None)
    parser.add_argument(
        "--remote-data-subdir",
        default=None,
        help="Remote data subdirectory under remote-data-root. Defaults to the local data file parent directory name.",
    )

    parser.add_argument("--run-name-base", default="mlp_obs")
    parser.add_argument("--experiment-name", default="f12m01_base")
    parser.add_argument("--dataset-stem", default=None)
    parser.add_argument("--selected-device", default="gpu", choices=["cpu", "gpu"])
    parser.add_argument("--solver-name", default="resdmd_batch")
    parser.add_argument(
        "--residual-form",
        default="projected_kv",
        choices=["projected_kv", "projected_vlambda"],
    )
    parser.add_argument(
        "--observable-mode",
        default="abs",
        choices=["abs", "complex", "complex_split"],
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
    parser.add_argument("--skip-psi-export", action="store_true")

    parser.add_argument("--local-download-root", default=None)
    parser.add_argument("--download-checkpoints", action="store_true")
    parser.add_argument("--skip-code-sync", action="store_true")
    parser.add_argument("--skip-data-upload", action="store_true")
    parser.add_argument("--delete-remote-after-download", action="store_true")
    parser.add_argument("--poll-interval", type=int, default=30)
    parser.add_argument("--tail-lines", type=int, default=20)
    return parser.parse_args()
def resolve_observable_tag(observable_mode):
    observable_mode_map = {
        "abs": "abs",
        "complex": "complex_split",
        "complex_split": "complex_split",
    }
    return observable_mode_map[observable_mode]


def build_run_label(run_name_base, experiment_name, residual_form, observable_mode):
    observable_tag = resolve_observable_tag(observable_mode)
    residual_form_tag = sanitize_tag(residual_form)
    effective_experiment_name = f"{experiment_name}_{residual_form_tag}"
    run_label = f"{run_name_base}_{effective_experiment_name}_{observable_tag}"
    return run_label, effective_experiment_name, observable_tag


def ssh_target(args):
    return f"{args.ssh_user}@{args.ssh_host}"


def build_ssh_base(args):
    cmd = ["ssh", "-p", str(args.ssh_port)]
    if args.ssh_key:
        cmd.extend(["-i", args.ssh_key])
    cmd.append(ssh_target(args))
    return cmd


def build_scp_base(args):
    cmd = ["scp", "-P", str(args.ssh_port)]
    if args.ssh_key:
        cmd.extend(["-i", args.ssh_key])
    return cmd


def run_command(cmd, check=True, capture_output=False, text=True):
    result = subprocess.run(
        cmd,
        check=check,
        capture_output=capture_output,
        text=text,
    )
    return result


def run_ssh(args, remote_command, capture_output=False):
    cmd = build_ssh_base(args) + [remote_command]
    return run_command(cmd, check=True, capture_output=capture_output, text=True)


def scp_upload(args, local_path, remote_path):
    cmd = build_scp_base(args) + [str(local_path), f"{ssh_target(args)}:{remote_path}"]
    run_command(cmd)


def scp_download_dir(args, remote_dir, local_parent):
    local_parent.mkdir(parents=True, exist_ok=True)
    cmd = build_scp_base(args) + ["-r", f"{ssh_target(args)}:{remote_dir}", str(local_parent)]
    run_command(cmd)


def quote_remote(path_value):
    return shlex.quote(path_value)


def ensure_remote_dirs(args, *remote_dirs):
    mkdir_targets = " ".join(quote_remote(path_value) for path_value in remote_dirs)
    run_ssh(args, f"mkdir -p {mkdir_targets}")


def build_remote_python_command(args, remote_script_path, remote_data_filename, remote_data_subdir):
    parts = [
        args.remote_python,
        remote_script_path,
        "--project-root",
        args.remote_code_dir,
        "--solver-dir",
        args.remote_code_dir,
        "--data-root",
        args.remote_data_root,
        "--output-parent",
        args.remote_output_parent,
        "--checkpoint-parent",
        args.remote_checkpoint_parent,
        "--log-parent",
        args.remote_log_parent,
        "--run-name-base",
        args.run_name_base,
        "--experiment-name",
        args.experiment_name,
        "--selected-device",
        args.selected_device,
        "--solver-name",
        args.solver_name,
        "--residual-form",
        args.residual_form,
        "--data-subdir",
        remote_data_subdir,
        "--observable-mode",
        args.observable_mode,
        "--data-filename",
        remote_data_filename,
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
    ]

    for layer_size in args.layer_sizes:
        if "--layer-sizes" not in parts:
            parts.extend(["--layer-sizes", str(layer_size)])
        else:
            parts.append(str(layer_size))

    if args.resume:
        parts.append("--resume")
    if args.fresh_checkpoints:
        parts.append("--fresh-checkpoints")
    if args.skip_psi_export:
        parts.append("--skip-psi-export")

    return " ".join(shlex.quote(part) for part in parts)


def write_local_controller_manifest(local_manifest_path, payload):
    local_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    local_manifest_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main():
    args = parse_args()

    local_data_file = Path(args.local_data_file).resolve()
    if not local_data_file.exists():
        raise FileNotFoundError(f"Local data file not found: {local_data_file}")

    local_obs_info_file = None
    if args.local_obs_info_file:
        local_obs_info_file = Path(args.local_obs_info_file).resolve()
        if not local_obs_info_file.exists():
            raise FileNotFoundError(f"Local obs info file not found: {local_obs_info_file}")

    dataset_stem = args.dataset_stem or infer_dataset_stem_from_observable_file(local_data_file.name)
    if args.local_download_root is None:
        args.local_download_root = str(default_local_download_root(dataset_stem, "mlp"))

    remote_data_subdir = args.remote_data_subdir or local_data_file.parent.name
    remote_data_dir = posixpath.join(args.remote_data_root, remote_data_subdir)

    run_label, effective_experiment_name, observable_tag = build_run_label(
        args.run_name_base,
        args.experiment_name,
        args.residual_form,
        args.observable_mode,
    )
    remote_output_dir = posixpath.join(args.remote_output_parent, run_label)
    remote_checkpoint_dir = posixpath.join(args.remote_checkpoint_parent, run_label)
    remote_log_dir = posixpath.join(args.remote_log_parent, run_label)
    remote_log_file = posixpath.join(remote_log_dir, "controller_run.log")
    remote_pid_file = posixpath.join(remote_log_dir, "controller_run.pid")
    remote_exit_file = posixpath.join(remote_log_dir, "controller_run.exit")

    remote_script_path = posixpath.join(args.remote_code_dir, "run_autodl_reskoopnet_mlp.py")

    local_download_root = Path(args.local_download_root).resolve()
    local_outputs_parent = local_download_root / "outputs"
    local_logs_parent = local_download_root / "logs"
    local_checkpoints_parent = local_download_root / "checkpoints"
    local_manifest_path = local_download_root / f"{run_label}_controller_manifest.json"

    manifest = {
        "run_label": run_label,
        "dataset_stem": dataset_stem,
        "effective_experiment_name": effective_experiment_name,
        "observable_tag": observable_tag,
        "remote_code_dir": args.remote_code_dir,
        "remote_data_dir": remote_data_dir,
        "remote_output_dir": remote_output_dir,
        "remote_checkpoint_dir": remote_checkpoint_dir,
        "remote_log_dir": remote_log_dir,
        "remote_log_file": remote_log_file,
        "local_data_file": str(local_data_file),
        "local_obs_info_file": str(local_obs_info_file) if local_obs_info_file else None,
    }
    write_local_controller_manifest(local_manifest_path, manifest)

    print(f"run_label: {run_label}")
    print(f"effective_experiment_name: {effective_experiment_name}")
    print(f"remote_data_dir: {remote_data_dir}")
    print(f"remote_output_dir: {remote_output_dir}")
    print(f"remote_log_file: {remote_log_file}")
    print(f"local_download_root: {local_download_root}")
    print(f"controller_manifest: {local_manifest_path}")

    ensure_remote_dirs(
        args,
        args.remote_code_dir,
        remote_data_dir,
        remote_log_dir,
        remote_output_dir,
        remote_checkpoint_dir,
    )

    if not args.skip_code_sync:
        print("Uploading AutoDL code...")
        code_files = [
            AUTODL_DIR / "edmd_utils.py",
            AUTODL_DIR / "solver_resdmd_batch2.py",
            AUTODL_DIR / "solver_resdmd_batch3.py",
            AUTODL_DIR / "run_autodl_reskoopnet_mlp.py",
        ]
        for code_file in code_files:
            scp_upload(args, code_file, posixpath.join(args.remote_code_dir, code_file.name))
    else:
        print("Skipping code sync.")

    if not args.skip_data_upload:
        print("Uploading data package...")
        scp_upload(args, local_data_file, posixpath.join(remote_data_dir, local_data_file.name))
        if local_obs_info_file:
            scp_upload(
                args,
                local_obs_info_file,
                posixpath.join(remote_data_dir, local_obs_info_file.name),
            )
    else:
        print("Skipping data upload.")

    remote_python_command = build_remote_python_command(
        args,
        remote_script_path=remote_script_path,
        remote_data_filename=local_data_file.name,
        remote_data_subdir=remote_data_subdir,
    )
    if args.remote_env_setup:
        remote_python_command = f"{args.remote_env_setup} && {remote_python_command}"

    remote_body = (
        f"cd {quote_remote(args.remote_code_dir)} && "
        f"{remote_python_command}; "
        f"status=$?; "
        f"echo $status > {quote_remote(remote_exit_file)}; "
        f"exit $status"
    )
    launch_command = (
        f"mkdir -p {quote_remote(remote_log_dir)} && "
        f"rm -f {quote_remote(remote_pid_file)} {quote_remote(remote_exit_file)} && "
        f"nohup bash -lc {quote_remote(remote_body)} "
        f"> {quote_remote(remote_log_file)} 2>&1 < /dev/null & "
        f"echo $! > {quote_remote(remote_pid_file)} && "
        f"echo STARTED"
    )

    print("Launching remote run...")
    launch_result = run_ssh(args, launch_command, capture_output=True)
    launch_stdout = launch_result.stdout.strip()
    if launch_stdout:
        print(launch_stdout)

    print("Polling remote run...")
    previous_tail = None
    exit_status = None
    while True:
        status_command = (
            f"if [ -f {quote_remote(remote_pid_file)} ]; then "
            f"PID=$(cat {quote_remote(remote_pid_file)}); "
            f"if kill -0 \"$PID\" 2>/dev/null; then "
            f"echo RUNNING:$PID; "
            f"else "
            f"if [ -f {quote_remote(remote_exit_file)} ]; then "
            f"echo EXITED:$(cat {quote_remote(remote_exit_file)}); "
            f"else "
            f"echo STOPPED_NO_EXIT; "
            f"fi; "
            f"fi; "
            f"else "
            f"if [ -f {quote_remote(remote_exit_file)} ]; then "
            f"echo EXITED:$(cat {quote_remote(remote_exit_file)}); "
            f"else "
            f"echo PID_MISSING; "
            f"fi; "
            f"fi"
        )
        status_result = run_ssh(args, status_command, capture_output=True)
        status_text = status_result.stdout.strip()

        tail_command = f"tail -n {int(args.tail_lines)} {quote_remote(remote_log_file)} 2>/dev/null || true"
        tail_result = run_ssh(args, tail_command, capture_output=True)
        tail_text = tail_result.stdout.strip()
        if tail_text and tail_text != previous_tail:
            print("\n--- remote log tail ---")
            print(tail_text)
            print("--- end remote log tail ---\n")
            previous_tail = tail_text

        if status_text.startswith("EXITED:"):
            exit_status = int(status_text.split(":", 1)[1])
            print(f"Remote run finished with exit status {exit_status}.")
            break

        print(f"Remote status: {status_text}. Waiting {args.poll_interval}s...")
        time.sleep(args.poll_interval)

    print("Downloading remote logs...")
    scp_download_dir(args, remote_log_dir, local_logs_parent)

    if exit_status == 0:
        print("Downloading remote outputs...")
        scp_download_dir(args, remote_output_dir, local_outputs_parent)
        if args.download_checkpoints:
            print("Downloading remote checkpoints...")
            scp_download_dir(args, remote_checkpoint_dir, local_checkpoints_parent)
    else:
        print("Skipping output download because the remote run did not exit cleanly.")

    if args.delete_remote_after_download and exit_status == 0:
        print("Deleting remote run artifacts...")
        delete_command = (
            f"rm -rf "
            f"{quote_remote(remote_output_dir)} "
            f"{quote_remote(remote_checkpoint_dir)} "
            f"{quote_remote(remote_log_dir)}"
        )
        run_ssh(args, delete_command)

    if exit_status != 0:
        raise SystemExit(exit_status)

    print("Controller finished successfully.")


if __name__ == "__main__":
    main()
