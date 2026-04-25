import argparse
import json
import shlex
import subprocess
import sys
import time
from pathlib import Path


OBSERVABLE_BATCH_PATH = Path(__file__).with_name("batch_controller_autodl_reskoopnet_mlp.py")


def sanitize_tag(value):
    return str(value).replace("/", "-").replace(" ", "_")


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run multiple observable packages for one dataset. Each observable package "
            "runs its residual_form batch sequentially, reusing the existing per-observable batch controller."
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

    parser.add_argument("--dataset-stem", required=True)
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
    parser.add_argument("--remote-data-subdir", default=None)

    parser.add_argument("--local-data-file-abs", default=None)
    parser.add_argument("--local-obs-info-file-abs", default=None)
    parser.add_argument("--local-data-file-complex-split", default=None)
    parser.add_argument("--local-obs-info-file-complex-split", default=None)
    parser.add_argument("--local-data-file-complex", default=None)
    parser.add_argument("--local-obs-info-file-complex", default=None)
    parser.add_argument("--local-data-file-eleHP", default=None)
    parser.add_argument("--local-obs-info-file-eleHP", default=None)
    parser.add_argument("--local-data-file-HP", default=None)
    parser.add_argument("--local-obs-info-file-HP", default=None)
    parser.add_argument("--local-data-file-identity", default=None)
    parser.add_argument("--local-obs-info-file-identity", default=None)
    parser.add_argument("--local-data-file-roi-mean", default=None)
    parser.add_argument("--local-obs-info-file-roi-mean", default=None)
    parser.add_argument("--local-data-file-slow-band-power", default=None)
    parser.add_argument("--local-obs-info-file-slow-band-power", default=None)
    parser.add_argument("--local-data-file-svd", default=None)
    parser.add_argument("--local-obs-info-file-svd", default=None)
    parser.add_argument("--local-data-file-HP-svd100", default=None)
    parser.add_argument("--local-obs-info-file-HP-svd100", default=None)
    parser.add_argument("--local-data-file-global-svd100", default=None)
    parser.add_argument("--local-obs-info-file-global-svd100", default=None)
    parser.add_argument(
        "--use-smoke-inputs",
        action="store_true",
        help=(
            "Resolve local observable inputs automatically from "
            "python_scripts/autodl/test_inputs/<dataset_stem>/."
        ),
    )
    parser.add_argument(
        "--smoke-input-root",
        default=str(Path(__file__).with_name("test_inputs")),
        help="Root directory containing per-dataset smoke inputs.",
    )
    parser.add_argument(
        "--smoke-pattern-tag",
        default="smoke20k_v73",
        help="Suffix tag used to locate smoke input files.",
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
    parser.add_argument("--local-download-root", default=None)
    parser.add_argument("--download-checkpoints", action="store_true")
    parser.add_argument("--skip-data-upload", action="store_true")
    parser.add_argument("--recover-completed-remote-runs", action="store_true")
    parser.add_argument("--skip-completed-local-runs", action="store_true")
    parser.add_argument(
        "--wait-for-remote-input",
        action="store_true",
        help=(
            "When used with --skip-data-upload, poll the expected remote input path until the "
            "observable file appears before starting that observable batch."
        ),
    )
    parser.add_argument(
        "--remote-input-timeout-sec",
        type=int,
        default=72 * 3600,
        help="Maximum time to wait for a manually uploaded remote input file.",
    )
    parser.add_argument(
        "--remote-input-poll-sec",
        type=int,
        default=60,
        help="Polling interval while waiting for a manually uploaded remote input file.",
    )
    parser.add_argument("--delete-remote-run-after-download", action="store_true")
    parser.add_argument(
        "--delete-remote-input-after-observable",
        action="store_true",
        help="After one observable batch succeeds locally, delete its remote input file(s). Disabled by default.",
    )
    parser.add_argument("--poll-interval", type=int, default=15)
    parser.add_argument("--tail-lines", type=int, default=80)
    parser.add_argument("--launch-timeout-sec", type=int, default=30)
    parser.add_argument("--continue-on-error", action="store_true")
    return parser.parse_args()


def observable_file_bundle(args, observable_mode):
    mapping = {
        "abs": {
            "data": args.local_data_file_abs,
            "obs_info": args.local_obs_info_file_abs,
        },
        "complex_split": {
            "data": args.local_data_file_complex_split,
            "obs_info": args.local_obs_info_file_complex_split,
        },
        "complex": {
            "data": args.local_data_file_complex or args.local_data_file_complex_split,
            "obs_info": args.local_obs_info_file_complex or args.local_obs_info_file_complex_split,
        },
        "eleHP": {
            "data": args.local_data_file_eleHP,
            "obs_info": args.local_obs_info_file_eleHP,
        },
        "HP": {
            "data": args.local_data_file_HP,
            "obs_info": args.local_obs_info_file_HP,
        },
        "identity": {
            "data": args.local_data_file_identity,
            "obs_info": args.local_obs_info_file_identity,
        },
        "roi_mean": {
            "data": args.local_data_file_roi_mean,
            "obs_info": args.local_obs_info_file_roi_mean,
        },
        "slow_band_power": {
            "data": args.local_data_file_slow_band_power,
            "obs_info": args.local_obs_info_file_slow_band_power,
        },
        "svd": {
            "data": args.local_data_file_svd,
            "obs_info": args.local_obs_info_file_svd,
        },
        "HP_svd100": {
            "data": args.local_data_file_HP_svd100,
            "obs_info": args.local_obs_info_file_HP_svd100,
        },
        "global_svd100": {
            "data": args.local_data_file_global_svd100,
            "obs_info": args.local_obs_info_file_global_svd100,
        },
    }
    bundle = dict(mapping[observable_mode])

    if args.use_smoke_inputs:
        smoke_data_file = resolve_smoke_data_file(args, observable_mode)
        bundle["data"] = str(smoke_data_file)
        smoke_obs_info = smoke_data_file.with_name(f"{smoke_data_file.stem}_obs_info.csv")
        bundle["obs_info"] = str(smoke_obs_info) if smoke_obs_info.exists() else None

    if not bundle["data"]:
        raise ValueError(
            f"No local data file configured for observable_mode={observable_mode}. "
            "Provide the matching --local-data-file-* argument."
        )
    return bundle


def resolve_smoke_data_file(args, observable_mode):
    smoke_dataset_dir = Path(args.smoke_input_root).expanduser().resolve() / args.dataset_stem
    if not smoke_dataset_dir.exists():
        raise FileNotFoundError(
            f"Smoke input directory not found for dataset {args.dataset_stem}: {smoke_dataset_dir}"
        )

    patterns = [f"{args.dataset_stem}_*_{observable_mode}_single_{args.smoke_pattern_tag}.mat"]
    if observable_mode == "complex":
        patterns.append(f"{args.dataset_stem}_*_complex_split_single_{args.smoke_pattern_tag}.mat")

    matches = []
    for pattern in patterns:
        matches.extend(sorted(smoke_dataset_dir.glob(pattern)))

    unique_matches = list(dict.fromkeys(path.resolve() for path in matches))
    if len(unique_matches) == 1:
        return unique_matches[0]
    if not unique_matches:
        raise FileNotFoundError(
            "No smoke input matched "
            f"observable_mode={observable_mode} under {smoke_dataset_dir} "
            f"with tag {args.smoke_pattern_tag}."
        )
    raise ValueError(
        "Multiple smoke inputs matched "
        f"observable_mode={observable_mode}: {', '.join(str(path) for path in unique_matches)}"
    )


def build_batch_cmd(args, observable_mode, bundle):
    cmd = [
        sys.executable,
        str(OBSERVABLE_BATCH_PATH),
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
        bundle["data"],
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
        observable_mode,
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
        "--residual-forms",
        *args.residual_forms,
    ]

    if args.ssh_key:
        cmd.extend(["--ssh-key", args.ssh_key])
    if args.remote_env_setup:
        cmd.extend(["--remote-env-setup", args.remote_env_setup])
    if bundle["obs_info"]:
        cmd.extend(["--local-obs-info-file", bundle["obs_info"]])
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
        cmd.append("--recover-completed-remote-runs")
    if args.skip_completed_local_runs:
        cmd.append("--skip-completed-local-runs")
    if args.skip_data_upload:
        cmd.append("--skip-data-upload")
    if args.delete_remote_run_after_download:
        cmd.append("--delete-remote-run-after-download")
    if args.continue_on_error:
        cmd.append("--continue-on-error")

    return cmd


def resolve_remote_input_paths(args, bundle):
    remote_data_subdir = args.remote_data_subdir or Path(bundle["data"]).resolve().parent.name
    remote_base = args.remote_data_root.rstrip("/")
    remote_paths = [f"{remote_base}/{remote_data_subdir}/{Path(bundle['data']).name}"]
    if bundle["obs_info"]:
        remote_paths.append(f"{remote_base}/{remote_data_subdir}/{Path(bundle['obs_info']).name}")
    return remote_paths


def build_ssh_probe_command(args, remote_path):
    cmd = ["ssh", "-p", str(args.ssh_port)]
    if args.ssh_key:
        cmd.extend(["-i", args.ssh_key])
    cmd.append(f"{args.ssh_user}@{args.ssh_host}")
    cmd.append(f"test -f {shlex.quote(remote_path)}")
    return cmd


def wait_for_remote_input(args, remote_path, observable_mode):
    deadline = time.time() + args.remote_input_timeout_sec
    print(
        f"[dataset-batch] waiting for remote input for observable_mode={observable_mode}: {remote_path}"
    )
    while True:
        probe = subprocess.run(
            build_ssh_probe_command(args, remote_path),
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if probe.returncode == 0:
            print(
                f"[dataset-batch] remote input detected for observable_mode={observable_mode}: {remote_path}"
            )
            return True

        now = time.time()
        if now >= deadline:
            print(
                f"[dataset-batch] timed out while waiting for remote input for "
                f"observable_mode={observable_mode}: {remote_path}"
            )
            return False

        remaining_sec = int(deadline - now)
        print(
            f"[dataset-batch] remote input not found yet for observable_mode={observable_mode}. "
            f"Rechecking in {args.remote_input_poll_sec}s (remaining ~{remaining_sec}s)."
        )
        time.sleep(args.remote_input_poll_sec)


def build_ssh_delete_command(args, remote_paths):
    cmd = ["ssh", "-p", str(args.ssh_port)]
    if args.ssh_key:
        cmd.extend(["-i", args.ssh_key])
    cmd.append(f"{args.ssh_user}@{args.ssh_host}")
    quoted_paths = " ".join(shlex.quote(path) for path in remote_paths)
    remote_data_dir = str(Path(remote_paths[0]).parent).replace("\\", "/")
    cmd.append(f"rm -f {quoted_paths} && rmdir {shlex.quote(remote_data_dir)} 2>/dev/null || true")
    return cmd


def main():
    args = parse_args()

    manifest = {
        "dataset_stem": args.dataset_stem,
        "experiment_name": args.experiment_name,
        "observable_modes": list(args.observable_modes),
        "residual_forms": list(args.residual_forms),
        "loss_mode": "squared",
        "requested_loss_mode": args.loss_mode,
        "loss_epsilon": args.loss_epsilon,
        "use_smoke_inputs": bool(args.use_smoke_inputs),
        "smoke_input_root": args.smoke_input_root,
        "smoke_pattern_tag": args.smoke_pattern_tag,
        "started_at_unix": time.time(),
        "observables": [],
    }

    manifest_dir = Path(args.local_download_root).resolve() if args.local_download_root else Path(__file__).resolve().parent
    manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = (
        manifest_dir / f"{args.dataset_stem}_{args.experiment_name}_dataset_batch_manifest.json"
    )

    for observable_mode in args.observable_modes:
        bundle = observable_file_bundle(args, observable_mode)
        remote_input_paths = resolve_remote_input_paths(args, bundle)
        batch_cmd = build_batch_cmd(args, observable_mode, bundle)

        print("")
        print("=" * 80)
        print(f"[dataset-batch] starting observable_mode={observable_mode}")
        print("=" * 80)
        print(f"[dataset-batch] local data file: {bundle['data']}")
        if bundle["obs_info"]:
            print(f"[dataset-batch] local obs_info file: {bundle['obs_info']}")
        if args.skip_data_upload:
            print(f"[dataset-batch] expected remote input file: {remote_input_paths[0]}")

        if args.skip_data_upload and args.wait_for_remote_input:
            wait_started_at = time.time()
            found_remote_input = wait_for_remote_input(args, remote_input_paths[0], observable_mode)
            wait_duration_sec = time.time() - wait_started_at
            if not found_remote_input:
                observable_record = {
                    "observable_mode": observable_mode,
                    "local_data_file": bundle["data"],
                    "local_obs_info_file": bundle["obs_info"],
                    "loss_mode": "squared",
                    "requested_loss_mode": args.loss_mode,
                    "return_code": 1,
                    "duration_sec": 0.0,
                    "waited_for_remote_input": True,
                    "remote_input_wait_duration_sec": wait_duration_sec,
                    "remote_input_deleted": False,
                    "remote_input_paths": remote_input_paths,
                }
                manifest["observables"].append(observable_record)
                manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
                if not args.continue_on_error:
                    raise SystemExit(1)
                print(
                    f"[dataset-batch] skipping observable_mode={observable_mode} after remote input wait timed out"
                )
                continue

        started_at = time.time()
        result = subprocess.run(batch_cmd, check=False)
        duration_sec = time.time() - started_at

        observable_record = {
            "observable_mode": observable_mode,
            "local_data_file": bundle["data"],
            "local_obs_info_file": bundle["obs_info"],
            "loss_mode": "squared",
            "requested_loss_mode": args.loss_mode,
            "return_code": int(result.returncode),
            "duration_sec": duration_sec,
            "waited_for_remote_input": bool(args.skip_data_upload and args.wait_for_remote_input),
            "remote_input_deleted": False,
            "remote_input_paths": remote_input_paths,
        }

        if result.returncode == 0 and args.delete_remote_input_after_observable:
            print(f"[dataset-batch] deleting remote input for observable_mode={observable_mode}:")
            for remote_path in remote_input_paths:
                print(f"  {remote_path}")
            delete_result = subprocess.run(build_ssh_delete_command(args, remote_input_paths), check=False)
            observable_record["remote_input_deleted"] = delete_result.returncode == 0
            observable_record["remote_input_delete_return_code"] = int(delete_result.returncode)
            if delete_result.returncode == 0:
                print(f"[dataset-batch] remote input deleted for observable_mode={observable_mode}")
            else:
                print(
                    f"[dataset-batch] remote input deletion failed for observable_mode={observable_mode} "
                    f"with return code {delete_result.returncode}"
                )

        manifest["observables"].append(observable_record)
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

        if result.returncode != 0:
            print(f"[dataset-batch] observable_mode={observable_mode} failed with return code {result.returncode}")
            if not args.continue_on_error:
                raise SystemExit(result.returncode)
        else:
            print(f"[dataset-batch] observable_mode={observable_mode} finished successfully in {duration_sec:.1f}s")

    manifest["finished_at_unix"] = time.time()
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("")
    print(f"[dataset-batch] finished. Manifest saved to: {manifest_path}")


if __name__ == "__main__":
    main()
