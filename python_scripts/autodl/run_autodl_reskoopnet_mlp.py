import argparse
import gc
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import scipy.io
import tensorflow as tf
from koopmanlib.dictionary import PsiNN

import edmd_utils


CHUNK_EXCLUDED_SHARED_FIELDS = {
    "snapshot_valid_idx_x",
    "snapshot_valid_idx_y",
    "snapshot_session_idx",
    "snapshot_train_pair_idx",
    "snapshot_valid_pair_idx",
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run the AutoDL MLP ResKoopNet pipeline and export outputs from the best checkpoint."
    )
    parser.add_argument("--project-root", default="/root/projects")
    parser.add_argument("--solver-dir", default=None)
    parser.add_argument("--data-root", default="/root/autodl-tmp/data")
    parser.add_argument("--output-parent", default="/root/autodl-tmp/outputs")
    parser.add_argument("--checkpoint-parent", default="/root/autodl-tmp/checkpoints")
    parser.add_argument("--log-parent", default="/root/autodl-tmp/logs")
    parser.add_argument("--run-name-base", default="mlp_obs")
    parser.add_argument("--experiment-name", default="f12m01_base")
    parser.add_argument("--selected-device", default="gpu", choices=["cpu", "gpu"])
    parser.add_argument(
        "--require-gpu",
        action="store_true",
        help="Fail immediately if --selected-device gpu is requested but TensorFlow cannot see a GPU.",
    )
    parser.add_argument("--solver-name", default="resdmd_batch")
    parser.add_argument(
        "--residual-form",
        default="projected_kv",
        choices=["projected_kv", "projected_vlambda"],
    )
    parser.add_argument(
        "--loss-mode",
        default="squared",
        choices=["squared", "per_dim", "relative_target"],
        help=(
            "Legacy compatibility flag. Optimization now always uses the "
            "historical squared residual loss; per_dim and relative_target "
            "are still exported as diagnostic metrics."
        ),
    )
    parser.add_argument(
        "--loss-epsilon",
        type=float,
        default=1e-12,
        help="Epsilon used only when computing relative_target diagnostic metrics.",
    )
    parser.add_argument("--data-subdir", default="f12m01")
    parser.add_argument("--dataset-stem", default="f12m01")
    parser.add_argument(
        "--observable-mode",
        default="abs",
        choices=[
            "abs", "complex", "complex_split", "eleHP", "HP", "identity",
            "roi_mean", "slow_band_power", "svd", "HP_svd100", "global_svd100",
        ],
    )
    parser.add_argument("--data-filename", default=None)
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
    parser.add_argument(
        "--resume-mode",
        default=None,
        choices=["fresh", "final", "best"],
        help="Checkpoint/state restore mode. Defaults to final when --resume is set, otherwise fresh.",
    )
    parser.add_argument("--fresh-checkpoints", action="store_true")
    parser.set_defaults(export_psi=False)
    parser.add_argument(
        "--export-psi",
        dest="export_psi",
        action="store_true",
        help="Opt in to exporting Psi_X/Psi_Y chunks. Disabled by default to save remote disk space.",
    )
    parser.add_argument(
        "--skip-psi-export",
        dest="export_psi",
        action="store_false",
        help=argparse.SUPPRESS,
    )
    return parser.parse_args()


def sanitize_tag(value):
    return str(value).replace("/", "-").replace(" ", "_")


def resolve_observable_tag(observable_mode):
    observable_mode_map = {
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
    return observable_mode_map[observable_mode]


def ensure_sys_path(project_root, solver_dir):
    for path in [project_root, solver_dir]:
        if path not in sys.path:
            sys.path.append(path)


def ensure_roots(*paths):
    for path in paths:
        os.makedirs(path, exist_ok=True)


def restore_best_checkpoint_for_export(solver, best_checkpoint_path):
    checkpoint_dir = os.path.dirname(best_checkpoint_path)
    latest_checkpoint = tf.train.latest_checkpoint(checkpoint_dir)
    if latest_checkpoint is None:
        raise FileNotFoundError(
            f"No best checkpoint found under {checkpoint_dir}. "
            "Export is configured to require the best checkpoint."
        )

    k_var = tf.Variable(solver.K, trainable=False, name="K")
    reg_var = tf.Variable(solver.reg, trainable=False, name="reg")
    checkpoint = tf.train.Checkpoint(
        model=solver.model,
        K=k_var,
        eigenvectors=solver.eigenvectors_var,
        eigenvalues=solver.eigenvalues_var,
        reg=reg_var,
    )
    checkpoint.restore(latest_checkpoint)

    solver.K = tf.identity(k_var.read_value())
    solver.eigenvectors = solver.eigenvectors_var.numpy().astype(np.complex128, copy=False)
    solver.eigenvalues = solver.eigenvalues_var.numpy().astype(np.complex128, copy=False)
    solver.reg = float(reg_var.numpy())
    solver.eigvec_cond = float(np.linalg.cond(solver.eigenvectors))
    solver.compute_mode()
    solver._sync_training_spectral_state()
    print(f"Restored best checkpoint for export: {latest_checkpoint}")


def as_scalar_or_nan(value):
    try:
        if value is None:
            return np.nan
        arr = np.asarray(value).squeeze()
        if arr.size == 0:
            return np.nan
        if np.iscomplexobj(arr):
            arr = np.real(arr)
        return float(arr.item())
    except Exception:
        return np.nan


def history_array(history, key):
    return np.asarray(
        [as_scalar_or_nan(item.get(key, np.nan)) for item in history],
        dtype=np.float32,
    )


def solver_attr_array(solver, key):
    return np.ravel(np.asarray(getattr(solver, key, []), dtype=np.float32))


def build_shared_outputs(solver, loss_all, val_loss_all, export_metadata=None, snapshot_metadata=None):
    evalues = np.asarray(solver.eigenvalues.T, dtype=np.complex64)
    kpm_modes = np.asarray(solver.compute_mode().T, dtype=np.complex64)
    loss_array = np.asarray(loss_all, dtype=np.float32)
    val_loss_array = np.asarray(val_loss_all, dtype=np.float32)
    n_dict = int(np.shape(evalues)[0])
    export_metadata = export_metadata or {}
    snapshot_metadata = snapshot_metadata or {}

    outer_history = getattr(solver, "outer_history", None) or []
    residual_form = str(getattr(solver, "residual_form", ""))
    loss_mode = str(getattr(solver, "loss_mode", export_metadata.get("loss_mode", "squared")))
    requested_loss_mode = str(
        getattr(
            solver,
            "requested_loss_mode",
            export_metadata.get("requested_loss_mode", loss_mode),
        )
    )
    loss_epsilon = np.float64(
        as_scalar_or_nan(getattr(solver, "loss_epsilon", export_metadata.get("loss_epsilon", np.nan)))
    )
    time_metadata = export_metadata.get("time_metadata", {}) if export_metadata else {}
    final_eigvec_cond = np.float32(as_scalar_or_nan(getattr(solver, "eigvec_cond", np.nan)))
    best_checkpoint_path = str(getattr(solver, "best_checkpoint_path", "") or "")
    final_checkpoint_path = str(getattr(solver, "final_checkpoint_path", "") or "")
    num_outer_epochs = int(len(outer_history))
    inner_train_metric_squared_history = solver_attr_array(solver, "inner_train_metric_squared_history")
    inner_val_metric_squared_history = solver_attr_array(solver, "inner_val_metric_squared_history")
    inner_train_metric_per_dim_history = solver_attr_array(solver, "inner_train_metric_per_dim_history")
    inner_val_metric_per_dim_history = solver_attr_array(solver, "inner_val_metric_per_dim_history")
    inner_train_metric_relative_target_history = solver_attr_array(
        solver, "inner_train_metric_relative_target_history"
    )
    inner_val_metric_relative_target_history = solver_attr_array(
        solver, "inner_val_metric_relative_target_history"
    )

    if num_outer_epochs > 0:
        outer_epoch_history = history_array(outer_history, "outer_epoch")
        outer_train_metric_history = history_array(outer_history, "train_metric")
        outer_val_metric_history = history_array(outer_history, "val_metric")
        outer_eigvec_cond_history = history_array(outer_history, "eigvec_cond")
        outer_lr_history = history_array(outer_history, "lr")
        outer_reg_history = history_array(outer_history, "reg")
        outer_inner_train_last_history = history_array(outer_history, "inner_train_last")
        outer_inner_val_last_history = history_array(outer_history, "inner_val_last")
        outer_train_metric_squared_history = history_array(outer_history, "train_metric_squared")
        outer_val_metric_squared_history = history_array(outer_history, "val_metric_squared")
        outer_train_metric_per_dim_history = history_array(outer_history, "train_metric_per_dim")
        outer_val_metric_per_dim_history = history_array(outer_history, "val_metric_per_dim")
        outer_train_metric_relative_target_history = history_array(
            outer_history, "train_metric_relative_target"
        )
        outer_val_metric_relative_target_history = history_array(
            outer_history, "val_metric_relative_target"
        )

        finite_val_mask = np.isfinite(outer_val_metric_history)
        if np.any(finite_val_mask):
            val_candidates = np.where(finite_val_mask, outer_val_metric_history, np.inf)
            best_idx = int(np.argmin(val_candidates))
            best_outer_index = best_idx
            best_outer_epoch = np.float32(outer_epoch_history[best_idx])
            best_val_metric = np.float32(outer_val_metric_history[best_idx])
            best_train_metric_at_best = np.float32(outer_train_metric_history[best_idx])
            best_eigvec_cond_at_best = np.float32(outer_eigvec_cond_history[best_idx])
            best_lr_at_best = np.float32(outer_lr_history[best_idx])
            best_reg_at_best = np.float32(outer_reg_history[best_idx])
        else:
            best_outer_index = np.nan
            best_outer_epoch = np.nan
            best_val_metric = np.nan
            best_train_metric_at_best = np.nan
            best_eigvec_cond_at_best = np.nan
            best_lr_at_best = np.nan
            best_reg_at_best = np.nan
    else:
        empty_history = np.array([], dtype=np.float32)
        outer_epoch_history = empty_history
        outer_train_metric_history = empty_history
        outer_val_metric_history = empty_history
        outer_eigvec_cond_history = empty_history
        outer_lr_history = empty_history
        outer_reg_history = empty_history
        outer_inner_train_last_history = empty_history
        outer_inner_val_last_history = empty_history
        outer_train_metric_squared_history = empty_history
        outer_val_metric_squared_history = empty_history
        outer_train_metric_per_dim_history = empty_history
        outer_val_metric_per_dim_history = empty_history
        outer_train_metric_relative_target_history = empty_history
        outer_val_metric_relative_target_history = empty_history
        best_outer_index = np.nan
        best_outer_epoch = np.nan
        best_val_metric = np.nan
        best_train_metric_at_best = np.nan
        best_eigvec_cond_at_best = np.nan
        best_lr_at_best = np.nan
        best_reg_at_best = np.nan

    shared_outputs = {
        "evalues": evalues,
        "kpm_modes": kpm_modes,
        "N_dict": n_dict,
        "loss": loss_array,
        "val_loss": val_loss_array,
        "residual_form": residual_form,
        "loss_mode": loss_mode,
        "requested_loss_mode": requested_loss_mode,
        "loss_epsilon": loss_epsilon,
        "inner_train_metric_squared_history": inner_train_metric_squared_history,
        "inner_val_metric_squared_history": inner_val_metric_squared_history,
        "inner_train_metric_per_dim_history": inner_train_metric_per_dim_history,
        "inner_val_metric_per_dim_history": inner_val_metric_per_dim_history,
        "inner_train_metric_relative_target_history": inner_train_metric_relative_target_history,
        "inner_val_metric_relative_target_history": inner_val_metric_relative_target_history,
        "num_outer_epochs": num_outer_epochs,
        "outer_epoch_history": outer_epoch_history,
        "outer_train_metric_history": outer_train_metric_history,
        "outer_val_metric_history": outer_val_metric_history,
        "outer_train_metric_squared_history": outer_train_metric_squared_history,
        "outer_val_metric_squared_history": outer_val_metric_squared_history,
        "outer_train_metric_per_dim_history": outer_train_metric_per_dim_history,
        "outer_val_metric_per_dim_history": outer_val_metric_per_dim_history,
        "outer_train_metric_relative_target_history": outer_train_metric_relative_target_history,
        "outer_val_metric_relative_target_history": outer_val_metric_relative_target_history,
        "outer_eigvec_cond_history": outer_eigvec_cond_history,
        "outer_lr_history": outer_lr_history,
        "outer_reg_history": outer_reg_history,
        "outer_inner_train_last_history": outer_inner_train_last_history,
        "outer_inner_val_last_history": outer_inner_val_last_history,
        "best_outer_index": best_outer_index,
        "best_outer_epoch": best_outer_epoch,
        "best_val_metric": best_val_metric,
        "best_train_metric_at_best": best_train_metric_at_best,
        "best_eigvec_cond_at_best": best_eigvec_cond_at_best,
        "best_lr_at_best": best_lr_at_best,
        "best_reg_at_best": best_reg_at_best,
        "final_eigvec_cond": final_eigvec_cond,
        "best_checkpoint_path": best_checkpoint_path,
        "final_checkpoint_path": final_checkpoint_path,
        "observable_tag": str(export_metadata.get("observable_tag", "")),
        "observable_mode": str(export_metadata.get("observable_mode", "")),
        "data_filename": str(export_metadata.get("data_filename", "")),
        "data_root": str(export_metadata.get("data_root", "")),
        "data_subdir": str(export_metadata.get("data_subdir", "")),
        "data_full_path": str(export_metadata.get("data_full_path", "")),
        "observable_file": str(time_metadata.get("observable_file", "")),
        "time_metadata_source": str(time_metadata.get("time_metadata_source", "")),
        "dataset_stem": str(export_metadata.get("dataset_stem", "")),
        "run_label": str(export_metadata.get("run_label", "")),
    }

    for key in ["dx", "dt", "fs", "sampling_period", "sample_period", "sampling_frequency"]:
        value = time_metadata.get(key, None)
        scalar_value = as_scalar_or_nan(value)
        if np.isfinite(scalar_value) and scalar_value > 0:
            shared_outputs[key] = np.float64(scalar_value)

    for key in [
        "session_dx",
        "session_fs",
        "session_ids",
        "session_lengths",
        "session_start_idx",
        "session_end_idx",
        "border_idx",
    ]:
        value = time_metadata.get(key, None)
        if value is not None:
            arr = np.asarray(value).squeeze()
            if arr.size:
                shared_outputs[key] = arr

    for key in ["valid_idx_x", "valid_idx_y", "session_idx", "train_pair_idx", "valid_pair_idx"]:
        value = snapshot_metadata.get(key, None)
        if value is not None:
            arr = np.asarray(value).squeeze()
            if arr.size:
                shared_outputs[f"snapshot_{key}"] = arr.astype(np.int64, copy=False) + 1

    if "lag" in snapshot_metadata:
        shared_outputs["snapshot_lag"] = int(snapshot_metadata["lag"])

    return shared_outputs, n_dict


def save_training_loss_diagnostics(output_root, title, loss_all, val_loss_all, solver):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Skipping training loss diagnostics plot: matplotlib import failed ({exc})")
        return None

    loss = np.ravel(np.asarray(loss_all, dtype=float))
    val_loss = np.ravel(np.asarray(val_loss_all, dtype=float))
    outer_history = getattr(solver, "outer_history", []) or []
    loss_mode = str(getattr(solver, "loss_mode", "squared"))
    outer_epoch = np.asarray([row.get("outer_epoch", np.nan) for row in outer_history], dtype=float)
    outer_train = np.asarray([row.get("train_metric", np.nan) for row in outer_history], dtype=float)
    outer_val = np.asarray([row.get("val_metric", np.nan) for row in outer_history], dtype=float)
    outer_cond = np.asarray([row.get("eigvec_cond", np.nan) for row in outer_history], dtype=float)
    inner_metric_specs = [
        (
            "squared",
            solver_attr_array(solver, "inner_train_metric_squared_history"),
            solver_attr_array(solver, "inner_val_metric_squared_history"),
        ),
        (
            "per_dim",
            solver_attr_array(solver, "inner_train_metric_per_dim_history"),
            solver_attr_array(solver, "inner_val_metric_per_dim_history"),
        ),
        (
            "relative_target",
            solver_attr_array(solver, "inner_train_metric_relative_target_history"),
            solver_attr_array(solver, "inner_val_metric_relative_target_history"),
        ),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), constrained_layout=True)
    fig.suptitle(title, fontsize=14)

    ax = axes[0, 0]
    if loss.size:
        ax.plot(np.arange(1, loss.size + 1), loss, label="inner train loss", linewidth=1.2)
    if val_loss.size:
        ax.plot(np.arange(1, val_loss.size + 1), val_loss, label="inner val loss", linewidth=1.2)
    ax.set_xlabel("inner fit epoch")
    ax.set_ylabel(f"{loss_mode} loss")
    ax.set_title("Inner loss history")
    ax.grid(True, alpha=0.3)
    ax.legend()

    ax = axes[0, 1]
    if outer_epoch.size and outer_train.size:
        ax.plot(outer_epoch, outer_train, marker="o", label="outer train metric", linewidth=1.2)
    if outer_epoch.size and outer_val.size:
        ax.plot(outer_epoch, outer_val, marker="o", label="outer val metric", linewidth=1.2)
    ax.set_xlabel("outer epoch")
    ax.set_ylabel(f"{loss_mode} metric")
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
    if outer_epoch.size and outer_val.size and np.any(np.isfinite(outer_val)):
        ax.scatter(outer_epoch, outer_val, s=24, label="outer val metric")
        best_idx = int(np.nanargmin(outer_val))
        ax.scatter(
            [outer_epoch[best_idx]],
            [outer_val[best_idx]],
            s=80,
            label=f"best epoch {int(outer_epoch[best_idx])}",
        )
    ax.set_xlabel("outer epoch")
    ax.set_ylabel(f"val {loss_mode} metric")
    ax.set_title("Best validation epoch")
    ax.grid(True, alpha=0.3)
    ax.legend()

    png_path = os.path.join(output_root, "training_loss_diagnostics.png")
    pdf_path = os.path.join(output_root, "training_loss_diagnostics.pdf")
    fig.savefig(png_path, dpi=180)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"saved training diagnostics: {png_path}")
    print(f"saved training diagnostics: {pdf_path}")

    has_inner_metric_history = any(train.size or val.size for _, train, val in inner_metric_specs)
    if has_inner_metric_history:
        fig2, axes2 = plt.subplots(3, 1, figsize=(11, 9), constrained_layout=True)
        fig2.suptitle(f"{title} - inner metrics", fontsize=14)
        for ax, (metric_name, train_metric, val_metric) in zip(axes2, inner_metric_specs):
            if train_metric.size:
                ax.plot(
                    np.arange(1, train_metric.size + 1),
                    train_metric,
                    label=f"train {metric_name}",
                    linewidth=1.2,
                )
            if val_metric.size:
                ax.plot(
                    np.arange(1, val_metric.size + 1),
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

        inner_png_path = os.path.join(output_root, "training_inner_loss_metrics.png")
        inner_pdf_path = os.path.join(output_root, "training_inner_loss_metrics.pdf")
        fig2.savefig(inner_png_path, dpi=180)
        fig2.savefig(inner_pdf_path)
        plt.close(fig2)
        print(f"saved inner loss metrics: {inner_png_path}")
        print(f"saved inner loss metrics: {inner_pdf_path}")
    return png_path


def get_num_samples(data_full_path, dataset_key):
    with h5py.File(data_full_path, "r") as f:
        return int(f[dataset_key].shape[1])


def _numeric_array_or_none(value):
    if value is None:
        return None
    arr = np.asarray(value).squeeze()
    if arr.size == 0 or not np.issubdtype(arr.dtype, np.number):
        return None
    if np.iscomplexobj(arr):
        arr = np.real(arr)
    return np.asarray(arr, dtype=np.float64)


def _positive_scalar_or_none(value):
    arr = _numeric_array_or_none(value)
    if arr is None or arr.size != 1:
        return None
    scalar = float(arr.item())
    if np.isfinite(scalar) and scalar > 0:
        return scalar
    return None


def _uniform_positive_value_or_none(value, abs_tol=1e-12, rel_tol=1e-5):
    arr = _numeric_array_or_none(value)
    if arr is None:
        return None
    arr = arr[np.isfinite(arr) & (arr > 0)]
    if arr.size == 0:
        return None
    center = float(np.median(arr))
    tol = max(abs_tol, rel_tol * abs(center))
    if np.all(np.abs(arr - center) <= tol):
        return center
    return None


def _read_observable_metadata_values(data_full_path, field_names):
    values = {}
    if h5py.is_hdf5(data_full_path):
        with h5py.File(data_full_path, "r") as f:
            for field_name in field_names:
                key = field_name if field_name.startswith("/") else f"/{field_name}"
                if key in f:
                    values[field_name] = np.asarray(f[key])
        return values

    mat = scipy.io.loadmat(data_full_path, variable_names=field_names)
    for field_name in field_names:
        if field_name in mat:
            values[field_name] = mat[field_name]
    return values


def read_observable_time_metadata(data_full_path):
    field_names = [
        "dx",
        "dt",
        "fs",
        "sampling_period",
        "sample_period",
        "sampling_frequency",
        "session_dx",
        "session_fs",
        "session_ids",
        "session_lengths",
        "session_start_idx",
        "session_end_idx",
        "border_idx",
    ]
    metadata = {
        "observable_file": str(data_full_path),
        "time_metadata_source": "missing",
    }
    values = _read_observable_metadata_values(data_full_path, field_names)

    for key, value in values.items():
        arr = _numeric_array_or_none(value)
        if arr is not None and arr.size:
            metadata[key] = arr

    dx = None
    source = ""
    for key in ["dx", "dt", "sampling_period", "sample_period"]:
        dx = _positive_scalar_or_none(values.get(key, None))
        if dx is not None:
            source = f"observable.{key}"
            break

    if dx is None:
        dx = _uniform_positive_value_or_none(values.get("session_dx", None))
        if dx is not None:
            source = "observable.session_dx"

    fs = None
    for key in ["fs", "sampling_frequency"]:
        fs = _positive_scalar_or_none(values.get(key, None))
        if fs is not None:
            break

    if fs is None:
        fs = _uniform_positive_value_or_none(values.get("session_fs", None))

    if dx is None and fs is not None:
        dx = 1.0 / fs
        source = "observable.fs"

    if dx is not None:
        fs = fs if fs is not None else 1.0 / dx
        metadata["dx"] = np.float64(dx)
        metadata["dt"] = np.float64(dx)
        metadata["sampling_period"] = np.float64(dx)
        metadata["sample_period"] = np.float64(dx)
        metadata["fs"] = np.float64(fs)
        metadata["sampling_frequency"] = np.float64(fs)
        metadata["time_metadata_source"] = source

    return metadata


def load_observable_chunk(data_full_path, dataset_key, start_idx, end_idx):
    with h5py.File(data_full_path, "r") as f:
        chunk = np.array(f[dataset_key][:, start_idx:end_idx]).T
    return chunk.astype(np.float32, copy=False)


def load_observable_rows(data_full_path, dataset_key, row_idx):
    row_idx = np.asarray(row_idx, dtype=np.int64).reshape(-1)
    with h5py.File(data_full_path, "r") as f:
        chunk = np.array(f[dataset_key][:, row_idx]).T
    return chunk.astype(np.float32, copy=False)


def export_outputs(
    solver,
    data_full_path,
    dataset_key,
    output_root,
    data_filename,
    chunk_size,
    loss_all,
    val_loss_all,
    export_metadata=None,
    snapshot_metadata=None,
):
    num_samples = get_num_samples(data_full_path, dataset_key)
    num_chunks = int(np.ceil(num_samples / float(chunk_size)))
    out_base = Path(data_filename).stem
    shared_outputs, n_dict = build_shared_outputs(
        solver,
        loss_all,
        val_loss_all,
        export_metadata=export_metadata,
        snapshot_metadata=snapshot_metadata,
    )

    summary_file = os.path.join(
        output_root,
        f"{out_base}_Python_resdmd_Layer_100_Ndict_{n_dict}_summary.mat",
    )
    scipy.io.savemat(summary_file, {"EDMD_outputs": shared_outputs}, long_field_names=True)
    print(f"saved summary: {summary_file}")

    # These snapshot arrays describe the full dataset split and can be tens of
    # MiB each. Keep them once in the summary instead of repeating them in every
    # output chunk.
    chunk_shared_outputs = {
        key: value
        for key, value in shared_outputs.items()
        if key not in CHUNK_EXCLUDED_SHARED_FIELDS
    }

    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, num_samples)
        chunk_data = load_observable_chunk(data_full_path, dataset_key, start_idx, end_idx)
        efuns = np.asarray(solver.eigenfunctions(chunk_data), dtype=np.complex64)

        edmd_outputs = {
            "EDMD_outputs": {
                "efuns": efuns,
                **chunk_shared_outputs,
            }
        }

        out_file = os.path.join(
            output_root,
            f"{out_base}_Python_resdmd_Layer_100_Ndict_{n_dict}_outputs_{i + 1}.mat",
        )
        scipy.io.savemat(out_file, edmd_outputs, long_field_names=True)
        print(f"saved: {out_file}")

        del chunk_data
        del efuns
        del edmd_outputs
        gc.collect()

    return n_dict


def export_psi_outputs(
    solver,
    data_full_path,
    dataset_key,
    output_root,
    data_filename,
    chunk_size,
    n_dict,
    snapshot_metadata=None,
):
    out_base = Path(data_filename).stem
    snapshot_metadata = snapshot_metadata or {}
    valid_idx_x = snapshot_metadata.get("valid_idx_x", None)
    valid_idx_y = snapshot_metadata.get("valid_idx_y", None)

    if valid_idx_x is not None and valid_idx_y is not None:
        valid_idx_x = np.asarray(valid_idx_x, dtype=np.int64).reshape(-1)
        valid_idx_y = np.asarray(valid_idx_y, dtype=np.int64).reshape(-1)
        if valid_idx_x.size != valid_idx_y.size:
            raise ValueError("snapshot valid_idx_x and valid_idx_y must have the same length")

        num_pairs = valid_idx_x.size
        num_chunks = int(np.ceil(num_pairs / float(chunk_size)))
        for i in range(num_chunks):
            start_idx = i * chunk_size
            end_idx = min((i + 1) * chunk_size, num_pairs)
            pair_idx = np.arange(start_idx, end_idx, dtype=np.int64)
            chunk_x = load_observable_rows(data_full_path, dataset_key, valid_idx_x[pair_idx])
            chunk_y = load_observable_rows(data_full_path, dataset_key, valid_idx_y[pair_idx])
            psi_x = solver.dic.call(chunk_x).numpy().T.astype(np.float32, copy=False)
            psi_y = solver.dic.call(chunk_y).numpy().T.astype(np.float32, copy=False)
            edmd_outputs = {
                "EDMD_outputs": {
                    "Psi_X": psi_x,
                    "Psi_Y": psi_y,
                    "snapshot_pair_idx": pair_idx + 1,
                    "snapshot_valid_idx_x": valid_idx_x[pair_idx] + 1,
                    "snapshot_valid_idx_y": valid_idx_y[pair_idx] + 1,
                }
            }

            out_file = os.path.join(
                output_root,
                f"{out_base}_Python_resdmd_Layer_100_Ndict_{n_dict}_outputs_Psi_{i + 1}.mat",
            )
            scipy.io.savemat(out_file, edmd_outputs, long_field_names=True)
            print(f"saved: {out_file}")

            del chunk_x
            del chunk_y
            del psi_x
            del psi_y
            del edmd_outputs
            gc.collect()
        return

    num_samples = get_num_samples(data_full_path, dataset_key)
    num_chunks = int(np.ceil(num_samples / float(chunk_size)))

    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, num_samples)
        chunk_data = load_observable_chunk(data_full_path, dataset_key, start_idx, end_idx)
        psi = solver.dic.call(chunk_data).numpy().T.astype(np.float32, copy=False)
        psi_x = psi[:, 0:-1]
        psi_y = psi[:, 1:]
        edmd_outputs = {
            "EDMD_outputs": {
                "Psi_X": psi_x,
                "Psi_Y": psi_y,
            }
        }

        out_file = os.path.join(
            output_root,
            f"{out_base}_Python_resdmd_Layer_100_Ndict_{n_dict}_outputs_Psi_{i + 1}.mat",
        )
        scipy.io.savemat(out_file, edmd_outputs, long_field_names=True)
        print(f"saved: {out_file}")

        del chunk_data
        del psi
        del psi_x
        del psi_y
        del edmd_outputs
        gc.collect()


def main():
    args = parse_args()

    project_root = args.project_root
    solver_dir = args.solver_dir or os.path.join(project_root, "solver")
    ensure_sys_path(project_root, solver_dir)
    ensure_roots(args.data_root, args.output_parent, args.checkpoint_parent, args.log_parent)

    print(f"project_root: {project_root}")
    print(f"solver_dir: {solver_dir}")
    print(f"data_root: {args.data_root}")
    print(f"output_parent: {args.output_parent}")
    print(f"checkpoint_parent: {args.checkpoint_parent}")
    print(f"log_parent: {args.log_parent}")
    print(f"run_name_base: {args.run_name_base}")
    print(f"experiment_name: {args.experiment_name}")

    actual_device = edmd_utils.set_device(args.selected_device, require_gpu=args.require_gpu)
    print(f"actual_device: {actual_device}")
    koopman_solver_cls = edmd_utils.import_solver(
        args.solver_name,
        search_dirs=[solver_dir, project_root],
    )

    basis_function = PsiNN(
        layer_sizes=list(args.layer_sizes),
        n_psi_train=args.n_psi_train,
    )

    observable_tag = resolve_observable_tag(args.observable_mode)
    data_filename = args.data_filename or (
        f"{args.dataset_stem}_low50_high250_g2_{observable_tag}_single.mat"
    )
    residual_form_tag = sanitize_tag(args.residual_form)
    effective_loss_mode = "squared"
    requested_loss_mode = args.loss_mode
    effective_experiment_name = f"{args.experiment_name}_{residual_form_tag}"
    run_label = f"{args.run_name_base}_{effective_experiment_name}_{observable_tag}"

    output_root = os.path.join(args.output_parent, run_label)
    checkpoint_root = os.path.join(args.checkpoint_parent, run_label)
    log_root = os.path.join(args.log_parent, run_label)
    checkpoint_path = os.path.join(checkpoint_root, "final", "model.ckpt")
    best_checkpoint_path = os.path.join(checkpoint_root, "best", "model.ckpt")

    ensure_roots(
        output_root,
        checkpoint_root,
        log_root,
        os.path.dirname(checkpoint_path),
        os.path.dirname(best_checkpoint_path),
    )

    file_type = args.file_type
    field_name = args.field_name
    data_input_path = os.path.join(args.data_root, args.data_subdir)
    data_full_path = os.path.join(data_input_path, data_filename)
    dataset_key = field_name if field_name.startswith("/") else f"/{field_name}"
    time_metadata = read_observable_time_metadata(data_full_path)

    print(f"observable_mode: {args.observable_mode} -> file tag: {observable_tag}")
    print(f"residual_form: {args.residual_form}")
    print(
        f"loss_mode: {effective_loss_mode} (requested={requested_loss_mode}), "
        f"loss_epsilon: {args.loss_epsilon:.3e}"
    )
    print(f"resume: {args.resume}, resume_mode: {args.resume_mode or ('final' if args.resume else 'fresh')}")
    print(f"effective_experiment_name: {effective_experiment_name}")
    print(f"run_label: {run_label}")
    print(f"output_root: {output_root}")
    print(f"checkpoint_root: {checkpoint_root}")
    print(f"log_root: {log_root}")
    print(f"checkpoint_path: {checkpoint_path}")
    print(f"best_checkpoint_path: {best_checkpoint_path}")
    print(f"data_full_path: {data_full_path}")
    if "dx" in time_metadata:
        print(
            "observable time metadata: "
            f"dx={float(time_metadata['dx']):.12g}, "
            f"fs={float(time_metadata['fs']):.12g}, "
            f"source={time_metadata['time_metadata_source']}"
        )
    else:
        print(
            "observable time metadata: missing dx/dt/session_dx; "
            "EDMD output chunks will not contain sampling interval metadata."
        )

    data = edmd_utils.load_data(
        data_filename,
        data_input_path,
        file_type,
        field_name,
    ).T.astype(np.float32, copy=False)
    print("Loaded DATA shape:", data.shape)
    print("Loaded DATA dtype:", data.dtype)
    print(f"Loaded DATA size: {data.nbytes / (1024 ** 3):.2f} GiB")

    snapshot_metadata = {}
    has_session_bounds = (
        "session_start_idx" in time_metadata
        and "session_end_idx" in time_metadata
        and np.asarray(time_metadata["session_start_idx"]).size > 0
        and np.asarray(time_metadata["session_end_idx"]).size > 0
    )
    if has_session_bounds:
        data_train, data_valid, snapshot_metadata = edmd_utils.transfer_data_format_session_aware(
            data,
            time_metadata["session_start_idx"],
            time_metadata["session_end_idx"],
            train_ratio=args.train_ratio,
            lag=1,
            seed=100,
            matlab_indexing=True,
        )
        print(
            "Using session-aware snapshot pairs: "
            f"{snapshot_metadata['valid_idx_x'].size} pairs, "
            f"{np.asarray(time_metadata['session_start_idx']).size} sessions"
        )
    else:
        data_train, data_valid = edmd_utils.transfer_data_format(data, train_ratio=args.train_ratio)
        print("Using continuous snapshot pairs because session metadata was not found.")
    print("Data shape: ", data_train[1].shape)
    del data
    gc.collect()

    solver = koopman_solver_cls(
        dic=basis_function,
        target_dim=data_train[0].shape[-1],
        reg=args.reg,
        residual_form=args.residual_form,
        loss_mode=effective_loss_mode,
        loss_epsilon=args.loss_epsilon,
    )

    checkpoint_cleanup_root = checkpoint_root
    print(f"Active experiment: {args.experiment_name}")
    print(f"Effective experiment: {effective_experiment_name}")
    print(f"Experiment checkpoint root: {checkpoint_cleanup_root}")
    if args.fresh_checkpoints:
        edmd_utils.remove_checkpoint(checkpoint_cleanup_root)
        ensure_roots(
            output_root,
            checkpoint_root,
            log_root,
            os.path.dirname(checkpoint_path),
            os.path.dirname(best_checkpoint_path),
        )

    loss_all = []
    val_loss_all = []
    completed_rounds = 0
    stop_flag = False

    for n_round in range(args.rounds):
        print("Round number: ", n_round)
        temp_loss, temp_val_loss, stop_flag, _ = solver.build(
            data_train=data_train,
            data_valid=data_valid,
            epochs=args.epochs,
            batch_size=args.batch_size,
            lr=args.lr,
            log_interval=args.log_interval,
            lr_decay_factor=args.lr_decay_factor,
            Nepoch=args.inner_epochs,
            end_condition=args.end_condition,
            checkpoint_path=checkpoint_path,
            best_checkpoint_path=best_checkpoint_path,
            resume=args.resume,
            resume_mode=args.resume_mode,
        )
        loss_all.extend(temp_loss)
        val_loss_all.extend(temp_val_loss)
        completed_rounds = n_round + 1
        if stop_flag:
            break

    print(f"Completed rounds: {completed_rounds}, stop_flag={stop_flag}")
    print(f"Collected losses: train={len(loss_all)}, valid={len(val_loss_all)}")

    restore_best_checkpoint_for_export(solver, best_checkpoint_path)
    export_metadata = {
        "observable_tag": observable_tag,
        "observable_mode": args.observable_mode,
        "data_filename": data_filename,
        "data_root": args.data_root,
        "data_subdir": args.data_subdir,
        "data_full_path": data_full_path,
        "dataset_stem": args.dataset_stem,
        "run_label": run_label,
        "loss_mode": effective_loss_mode,
        "requested_loss_mode": requested_loss_mode,
        "loss_epsilon": args.loss_epsilon,
        "time_metadata": time_metadata,
    }
    n_dict = export_outputs(
        solver=solver,
        data_full_path=data_full_path,
        dataset_key=dataset_key,
        output_root=output_root,
        data_filename=data_filename,
        chunk_size=args.chunk_size,
        loss_all=loss_all,
        val_loss_all=val_loss_all,
        export_metadata=export_metadata,
        snapshot_metadata=snapshot_metadata,
    )
    save_training_loss_diagnostics(
        output_root=output_root,
        title=(
            f"{args.dataset_stem} {observable_tag} {args.residual_form} "
            f"{effective_loss_mode} training diagnostics"
        ),
        loss_all=loss_all,
        val_loss_all=val_loss_all,
        solver=solver,
    )

    if args.export_psi:
        print("Exporting Psi_X/Psi_Y chunk-by-chunk.")
        export_psi_outputs(
            solver=solver,
            data_full_path=data_full_path,
            dataset_key=dataset_key,
            output_root=output_root,
            data_filename=data_filename,
            chunk_size=args.chunk_size,
            n_dict=n_dict,
            snapshot_metadata=snapshot_metadata,
        )
    else:
        print("Skipping Psi_X/Psi_Y export.")


if __name__ == "__main__":
    main()
