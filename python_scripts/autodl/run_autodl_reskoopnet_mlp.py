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
    parser.add_argument("--solver-name", default="resdmd_batch")
    parser.add_argument(
        "--residual-form",
        default="projected_kv",
        choices=["projected_kv", "projected_vlambda"],
    )
    parser.add_argument("--data-subdir", default="f12m01")
    parser.add_argument("--dataset-stem", default="f12m01")
    parser.add_argument(
        "--observable-mode",
        default="abs",
        choices=["abs", "complex", "complex_split"],
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
    parser.add_argument("--fresh-checkpoints", action="store_true")
    parser.add_argument("--skip-psi-export", action="store_true")
    return parser.parse_args()


def sanitize_tag(value):
    return str(value).replace("/", "-").replace(" ", "_")


def resolve_observable_tag(observable_mode):
    observable_mode_map = {
        "abs": "abs",
        "complex": "complex_split",
        "complex_split": "complex_split",
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


def build_shared_outputs(solver, loss_all, val_loss_all):
    evalues = np.asarray(solver.eigenvalues.T, dtype=np.complex64)
    kpm_modes = np.asarray(solver.compute_mode().T, dtype=np.complex64)
    loss_array = np.asarray(loss_all, dtype=np.float32)
    val_loss_array = np.asarray(val_loss_all, dtype=np.float32)
    n_dict = int(np.shape(evalues)[0])

    outer_history = getattr(solver, "outer_history", None) or []
    residual_form = str(getattr(solver, "residual_form", ""))
    final_eigvec_cond = np.float32(as_scalar_or_nan(getattr(solver, "eigvec_cond", np.nan)))
    best_checkpoint_path = str(getattr(solver, "best_checkpoint_path", "") or "")
    final_checkpoint_path = str(getattr(solver, "final_checkpoint_path", "") or "")
    num_outer_epochs = int(len(outer_history))

    if num_outer_epochs > 0:
        outer_epoch_history = history_array(outer_history, "outer_epoch")
        outer_train_metric_history = history_array(outer_history, "train_metric")
        outer_val_metric_history = history_array(outer_history, "val_metric")
        outer_eigvec_cond_history = history_array(outer_history, "eigvec_cond")
        outer_lr_history = history_array(outer_history, "lr")
        outer_reg_history = history_array(outer_history, "reg")
        outer_inner_train_last_history = history_array(outer_history, "inner_train_last")
        outer_inner_val_last_history = history_array(outer_history, "inner_val_last")

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
        "num_outer_epochs": num_outer_epochs,
        "outer_epoch_history": outer_epoch_history,
        "outer_train_metric_history": outer_train_metric_history,
        "outer_val_metric_history": outer_val_metric_history,
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
    }
    return shared_outputs, n_dict


def get_num_samples(data_full_path, dataset_key):
    with h5py.File(data_full_path, "r") as f:
        return int(f[dataset_key].shape[1])


def load_observable_chunk(data_full_path, dataset_key, start_idx, end_idx):
    with h5py.File(data_full_path, "r") as f:
        chunk = np.array(f[dataset_key][:, start_idx:end_idx]).T
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
):
    num_samples = get_num_samples(data_full_path, dataset_key)
    num_chunks = int(np.ceil(num_samples / float(chunk_size)))
    out_base = Path(data_filename).stem
    shared_outputs, n_dict = build_shared_outputs(solver, loss_all, val_loss_all)

    summary_file = os.path.join(
        output_root,
        f"{out_base}_Python_resdmd_Layer_100_Ndict_{n_dict}_summary.mat",
    )
    scipy.io.savemat(summary_file, {"EDMD_outputs": shared_outputs})
    print(f"saved summary: {summary_file}")

    for i in range(num_chunks):
        start_idx = i * chunk_size
        end_idx = min((i + 1) * chunk_size, num_samples)
        chunk_data = load_observable_chunk(data_full_path, dataset_key, start_idx, end_idx)
        efuns = np.asarray(solver.eigenfunctions(chunk_data), dtype=np.complex64)

        edmd_outputs = {
            "EDMD_outputs": {
                "efuns": efuns,
                **shared_outputs,
            }
        }

        out_file = os.path.join(
            output_root,
            f"{out_base}_Python_resdmd_Layer_100_Ndict_{n_dict}_outputs_{i + 1}.mat",
        )
        scipy.io.savemat(out_file, edmd_outputs)
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
):
    num_samples = get_num_samples(data_full_path, dataset_key)
    num_chunks = int(np.ceil(num_samples / float(chunk_size)))
    out_base = Path(data_filename).stem

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
        scipy.io.savemat(out_file, edmd_outputs)
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

    edmd_utils.set_device(args.selected_device)
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

    print(f"observable_mode: {args.observable_mode} -> file tag: {observable_tag}")
    print(f"residual_form: {args.residual_form}")
    print(f"effective_experiment_name: {effective_experiment_name}")
    print(f"run_label: {run_label}")
    print(f"output_root: {output_root}")
    print(f"checkpoint_root: {checkpoint_root}")
    print(f"log_root: {log_root}")
    print(f"checkpoint_path: {checkpoint_path}")
    print(f"best_checkpoint_path: {best_checkpoint_path}")
    print(f"data_full_path: {data_full_path}")

    data = edmd_utils.load_data(
        data_filename,
        data_input_path,
        file_type,
        field_name,
    ).T.astype(np.float32, copy=False)
    print("Loaded DATA shape:", data.shape)
    print("Loaded DATA dtype:", data.dtype)
    print(f"Loaded DATA size: {data.nbytes / (1024 ** 3):.2f} GiB")

    data_train, data_valid = edmd_utils.transfer_data_format(data, train_ratio=args.train_ratio)
    print("Data shape: ", data_train[1].shape)
    del data
    gc.collect()

    solver = koopman_solver_cls(
        dic=basis_function,
        target_dim=data_train[0].shape[-1],
        reg=args.reg,
        residual_form=args.residual_form,
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
        )
        loss_all.extend(temp_loss)
        val_loss_all.extend(temp_val_loss)
        completed_rounds = n_round + 1
        if stop_flag:
            break

    print(f"Completed rounds: {completed_rounds}, stop_flag={stop_flag}")
    print(f"Collected losses: train={len(loss_all)}, valid={len(val_loss_all)}")

    restore_best_checkpoint_for_export(solver, best_checkpoint_path)
    n_dict = export_outputs(
        solver=solver,
        data_full_path=data_full_path,
        dataset_key=dataset_key,
        output_root=output_root,
        data_filename=data_filename,
        chunk_size=args.chunk_size,
        loss_all=loss_all,
        val_loss_all=val_loss_all,
    )

    if args.skip_psi_export:
        print("Skipping Psi_X/Psi_Y export.")
    else:
        print("Exporting Psi_X/Psi_Y chunk-by-chunk.")
        export_psi_outputs(
            solver=solver,
            data_full_path=data_full_path,
            dataset_key=dataset_key,
            output_root=output_root,
            data_filename=data_filename,
            chunk_size=args.chunk_size,
            n_dict=n_dict,
        )


if __name__ == "__main__":
    main()
