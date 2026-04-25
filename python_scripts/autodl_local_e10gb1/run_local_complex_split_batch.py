import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


DATASET_CONFIGS = {
    "e10gb1": {
        "dataset_stem": "e10gb1",
        "data_subdir": "e10gb1/reskoopnet_dictionary",
        "data_filename": "e10gb1_low50_high250_g2_complex_split_single.mat",
        "results_subdir": "e10gb1",
    },
    "e10fV1": {
        "dataset_stem": "e10fV1",
        "data_subdir": "E10fV1/reskoopnet_dictionary",
        "data_filename": "e10fV1_low50_high250_g2_complex_split_single.mat",
        "results_subdir": "e10fV1",
    },
    "e10gh1": {
        "dataset_stem": "e10gh1",
        "data_subdir": "E10gH1/reskoopnet_dictionary",
        "data_filename": "e10gh1_low50_high250_g2_complex_split_single.mat",
        "results_subdir": "e10gh1",
    },
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run local complex_split ResKoopNet MLP jobs sequentially for multiple datasets "
            "using the local WSL-ready solver copy."
        )
    )
    parser.add_argument(
        "--dataset-stems",
        nargs="+",
        default=["e10gb1", "e10fV1", "e10gh1"],
        choices=sorted(DATASET_CONFIGS.keys()),
    )
    parser.add_argument("--experiment-name", default="local_complexsplit_60_20260418")
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
    parser.add_argument("--epochs", type=int, default=60)
    parser.add_argument("--python-exe", default=sys.executable)
    parser.add_argument("--data-root", default="/mnt/e/DataPons_processed")
    parser.add_argument("--results-root", default="/mnt/e/autodl_results")
    parser.add_argument("--selected-device", default="gpu", choices=["cpu", "gpu"])
    parser.add_argument("--solver-name", default="resdmd_batch")
    parser.add_argument("--file-type", default=".h5", choices=[".h5", ".mat"])
    parser.add_argument("--field-name", default="obs")
    parser.add_argument("--layer-sizes", type=int, nargs="+", default=[100, 100, 100])
    parser.add_argument("--n-psi-train", type=int, default=100)
    parser.add_argument("--train-ratio", type=float, default=0.7)
    parser.add_argument("--reg", type=float, default=0.1)
    parser.add_argument("--rounds", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=2000)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--log-interval", type=int, default=1)
    parser.add_argument("--lr-decay-factor", type=float, default=0.8)
    parser.add_argument("--inner-epochs", type=int, default=5)
    parser.add_argument("--end-condition", type=float, default=1e-9)
    parser.add_argument("--chunk-size", type=int, default=5000)
    parser.add_argument("--export-psi", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--fresh-checkpoints", action="store_true")
    parser.add_argument("--continue-on-error", action="store_true")
    return parser.parse_args()


def build_run_script_path():
    repo_root = Path(__file__).resolve().parents[2]
    return repo_root / "python_scripts" / "autodl" / "run_autodl_reskoopnet_mlp.py"


def build_local_solver_root():
    return Path(__file__).resolve().parent


def build_results_dirs(results_root, results_subdir):
    base = Path(results_root) / results_subdir / "mlp"
    return {
        "base": base,
        "outputs": base / "outputs",
        "checkpoints": base / "checkpoints",
        "logs": base / "logs",
    }


def sanitize_tag(value):
    return str(value).replace("/", "-").replace(" ", "_")


def build_run_label(experiment_name, residual_form, loss_mode="squared"):
    effective_experiment_name = f"{experiment_name}_{sanitize_tag(residual_form)}"
    return f"mlp_obs_{effective_experiment_name}_complex_split"


def generate_loss_diagnostics(output_dir, title):
    output_dir = Path(output_dir)
    summary_files = sorted(output_dir.glob("*_summary.mat"))
    if not summary_files:
        return {
            "ok": False,
            "reason": f"No summary file found in {output_dir}",
            "png": None,
            "pdf": None,
        }

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import scipy.io
    except Exception as exc:
        return {
            "ok": False,
            "reason": f"Could not import plotting dependencies: {type(exc).__name__}: {exc}",
            "png": None,
            "pdf": None,
        }

    summary_path = summary_files[0]
    try:
        edmd_outputs = scipy.io.loadmat(
            summary_path,
            squeeze_me=True,
            struct_as_record=False,
        )["EDMD_outputs"]
    except Exception as exc:
        return {
            "ok": False,
            "reason": f"Could not read summary file {summary_path}: {type(exc).__name__}: {exc}",
            "png": None,
            "pdf": None,
        }

    def field_array(name):
        value = getattr(edmd_outputs, name, np.array([]))
        try:
            return np.ravel(np.asarray(value, dtype=float))
        except Exception:
            return np.array([], dtype=float)

    loss = field_array("loss")
    val_loss = field_array("val_loss")
    outer_epoch = field_array("outer_epoch_history")
    outer_train = field_array("outer_train_metric_history")
    outer_val = field_array("outer_val_metric_history")
    outer_cond = field_array("outer_eigvec_cond_history")
    inner_metric_specs = [
        (
            "squared",
            field_array("inner_train_metric_squared_history"),
            field_array("inner_val_metric_squared_history"),
        ),
        (
            "per_dim",
            field_array("inner_train_metric_per_dim_history"),
            field_array("inner_val_metric_per_dim_history"),
        ),
        (
            "relative_target",
            field_array("inner_train_metric_relative_target_history"),
            field_array("inner_val_metric_relative_target_history"),
        ),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    fig.suptitle(title, fontsize=14)

    ax = axes[0, 0]
    if loss.size:
        ax.plot(range(1, loss.size + 1), loss, label="inner train loss", linewidth=1.2)
    if val_loss.size:
        ax.plot(range(1, val_loss.size + 1), val_loss, label="inner val loss", linewidth=1.2)
    ax.set_xlabel("inner fit epoch")
    ax.set_ylabel("loss")
    ax.set_title("Inner loss history")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[0, 1]
    if outer_epoch.size and outer_train.size:
        ax.plot(outer_epoch, outer_train, marker="o", label="outer train metric", linewidth=1.2)
    if outer_epoch.size and outer_val.size:
        ax.plot(outer_epoch, outer_val, marker="o", label="outer val metric", linewidth=1.2)
    ax.set_xlabel("outer epoch")
    ax.set_ylabel("metric")
    ax.set_title("Outer residual metric")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[1, 0]
    if outer_epoch.size and outer_cond.size:
        ax.plot(
            outer_epoch,
            np.log10(np.maximum(outer_cond, np.finfo(float).tiny)),
            marker="o",
            linewidth=1.2,
        )
    ax.set_xlabel("outer epoch")
    ax.set_ylabel("log10(eigvec cond)")
    ax.set_title("Spectral conditioning")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    finite_val = np.isfinite(outer_val)
    if outer_epoch.size and outer_val.size and np.any(finite_val):
        ax.scatter(outer_epoch, outer_val, s=24, label="outer val metric")
        best_idx = int(np.nanargmin(outer_val))
        ax.scatter(
            [outer_epoch[best_idx]],
            [outer_val[best_idx]],
            s=80,
            label=f"best epoch {int(outer_epoch[best_idx])}",
        )
    ax.set_xlabel("outer epoch")
    ax.set_ylabel("val metric")
    ax.set_title("Best validation epoch")
    ax.grid(True, alpha=0.3)
    ax.legend()

    png_path = output_dir / "training_loss_diagnostics.png"
    pdf_path = output_dir / "training_loss_diagnostics.pdf"
    fig.savefig(png_path, dpi=180)
    fig.savefig(pdf_path)
    plt.close(fig)

    inner_png_path = None
    inner_pdf_path = None
    has_inner_metric_history = any(train.size or val.size for _, train, val in inner_metric_specs)
    if has_inner_metric_history:
        fig2, axes2 = plt.subplots(3, 1, figsize=(11, 9), constrained_layout=True)
        fig2.suptitle(f"{title} - inner metrics", fontsize=14)
        for ax, (metric_name, train_metric, val_metric) in zip(axes2, inner_metric_specs):
            if train_metric.size:
                ax.plot(
                    range(1, train_metric.size + 1),
                    train_metric,
                    label=f"train {metric_name}",
                    linewidth=1.2,
                )
            if val_metric.size:
                ax.plot(
                    range(1, val_metric.size + 1),
                    val_metric,
                    label=f"val {metric_name}",
                    linewidth=1.2,
                )
            finite_values = np.concatenate(
                [
                    np.ravel(train_metric[np.isfinite(train_metric)]),
                    np.ravel(val_metric[np.isfinite(val_metric)]),
                ]
            )
            if finite_values.size and np.all(finite_values > 0):
                ax.set_yscale("log")
            ax.set_xlabel("inner fit epoch")
            ax.set_ylabel(metric_name)
            ax.grid(True, alpha=0.3)
            ax.legend()

        inner_png_path = output_dir / "training_inner_loss_metrics.png"
        inner_pdf_path = output_dir / "training_inner_loss_metrics.pdf"
        fig2.savefig(inner_png_path, dpi=180)
        fig2.savefig(inner_pdf_path)
        plt.close(fig2)

    return {
        "ok": True,
        "reason": "Generated training loss diagnostics.",
        "png": str(png_path),
        "pdf": str(pdf_path),
        "inner_png": str(inner_png_path) if inner_png_path else None,
        "inner_pdf": str(inner_pdf_path) if inner_pdf_path else None,
        "loss_len": int(loss.size),
        "val_loss_len": int(val_loss.size),
        "outer_epochs": int(outer_epoch.size),
    }


def build_run_cmd(args, dataset_key, residual_form):
    cfg = DATASET_CONFIGS[dataset_key]
    run_script = build_run_script_path()
    local_solver_root = build_local_solver_root()
    results_dirs = build_results_dirs(args.results_root, cfg["results_subdir"])

    cmd = [
        args.python_exe,
        str(run_script),
        "--project-root",
        str(local_solver_root),
        "--solver-dir",
        str(local_solver_root),
        "--data-root",
        args.data_root,
        "--output-parent",
        str(results_dirs["outputs"]),
        "--checkpoint-parent",
        str(results_dirs["checkpoints"]),
        "--log-parent",
        str(results_dirs["logs"]),
        "--experiment-name",
        args.experiment_name,
        "--selected-device",
        args.selected_device,
        "--solver-name",
        args.solver_name,
        "--residual-form",
        residual_form,
        "--loss-mode",
        args.loss_mode,
        "--loss-epsilon",
        str(args.loss_epsilon),
        "--data-subdir",
        cfg["data_subdir"],
        "--dataset-stem",
        cfg["dataset_stem"],
        "--observable-mode",
        "complex_split",
        "--data-filename",
        cfg["data_filename"],
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
        if "--layer-sizes" not in cmd:
            cmd.extend(["--layer-sizes", str(layer_size)])
        else:
            cmd.append(str(layer_size))

    if args.export_psi:
        cmd.append("--export-psi")

    if args.resume:
        cmd.append("--resume")
    elif args.fresh_checkpoints or not args.resume:
        cmd.append("--fresh-checkpoints")

    return cmd, results_dirs


def main():
    args = parse_args()

    local_solver_root = build_local_solver_root()
    run_script = build_run_script_path()
    manifest_path = local_solver_root / (
        f"{args.experiment_name}_local_complex_split_batch_manifest.json"
    )

    manifest = {
        "experiment_name": args.experiment_name,
        "dataset_stems": list(args.dataset_stems),
        "residual_forms": list(args.residual_forms),
        "loss_mode": "squared",
        "requested_loss_mode": args.loss_mode,
        "loss_epsilon": args.loss_epsilon,
        "observable_mode": "complex_split",
        "epochs": int(args.epochs),
        "run_script": str(run_script),
        "local_solver_root": str(local_solver_root),
        "started_at_unix": time.time(),
        "runs": [],
    }

    for dataset_key in args.dataset_stems:
        for residual_form in args.residual_forms:
            cmd, results_dirs = build_run_cmd(args, dataset_key, residual_form)

            print("")
            print("=" * 80)
            print(f"[local-batch] dataset={dataset_key} residual_form={residual_form}")
            print(f"[local-batch] results_root={results_dirs['base']}")
            print("=" * 80)

            started_at = time.time()
            result = subprocess.run(cmd, check=False)
            duration_sec = time.time() - started_at

            manifest["runs"].append(
                {
                    "dataset_stem": dataset_key,
                    "residual_form": residual_form,
                    "loss_mode": "squared",
                    "requested_loss_mode": args.loss_mode,
                    "return_code": int(result.returncode),
                    "duration_sec": duration_sec,
                    "results_root": str(results_dirs["base"]),
                }
            )
            manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

            if result.returncode != 0:
                print(
                    f"[local-batch] dataset={dataset_key} residual_form={residual_form} "
                    f"failed with return code {result.returncode}"
                )
                if not args.continue_on_error:
                    raise SystemExit(result.returncode)
            else:
                run_label = build_run_label(args.experiment_name, residual_form, args.loss_mode)
                output_dir = results_dirs["outputs"] / run_label
                plot_result = generate_loss_diagnostics(
                    output_dir,
                    f"{dataset_key} complex_split {residual_form} squared training diagnostics",
                )
                manifest["runs"][-1]["loss_plot"] = plot_result
                manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
                if plot_result["ok"]:
                    print(f"[local-batch] loss diagnostics: {plot_result['png']}")
                    if plot_result.get("inner_png"):
                        print(f"[local-batch] inner loss metrics: {plot_result['inner_png']}")
                else:
                    print(f"[local-batch] loss diagnostics skipped: {plot_result['reason']}")
                print(
                    f"[local-batch] dataset={dataset_key} residual_form={residual_form} "
                    f"finished successfully in {duration_sec:.1f}s"
                )

    manifest["finished_at_unix"] = time.time()
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("")
    print(f"[local-batch] finished. Manifest saved to: {manifest_path}")


if __name__ == "__main__":
    main()
