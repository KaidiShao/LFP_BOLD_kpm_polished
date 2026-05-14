import argparse
import importlib
import json
import math
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.integrate import solve_ivp
import tensorflow as tf


REPO_ROOT = Path(__file__).resolve().parents[2]
SOLVER_DIR = REPO_ROOT / "python_scripts" / "autodl"
if str(SOLVER_DIR) not in sys.path:
    sys.path.insert(0, str(SOLVER_DIR))


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Reproduce the ResKoopNet pendulum TensorFlow experiment using a "
            "selected local solver implementation."
        )
    )
    parser.add_argument(
        "--output-dir",
        default=str(REPO_ROOT / "results" / "pendulum_solver3_test"),
        help="Directory for generated data, checkpoints, plots, and summary files.",
    )
    parser.add_argument(
        "--solver-module",
        default="solver_resdmd_batch3",
        help="Python module name under python_scripts/autodl that provides KoopmanNN and KoopmanSolver.",
    )
    parser.add_argument(
        "--experiment-tag",
        default="pendulum_solver3",
        help="Prefix used in exported filenames and summary metadata.",
    )
    parser.add_argument(
        "--device",
        default="gpu",
        choices=["cpu", "gpu"],
        help="TensorFlow execution device preference.",
    )
    parser.add_argument(
        "--require-gpu",
        action="store_true",
        help="Fail if --device gpu is requested but TensorFlow cannot see a GPU.",
    )
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--delta-t", type=float, default=0.5)
    parser.add_argument(
        "--grid-angle-count",
        type=int,
        default=10,
        help="Corresponds to M1 in pendulum.m; 10 reproduces the 90-trajectory variant.",
    )
    parser.add_argument(
        "--grid-velocity-count",
        type=int,
        default=10,
        help="Corresponds to M2 in pendulum.m; 10 reproduces the 90-trajectory variant.",
    )
    parser.add_argument("--velocity-limit", type=float, default=15.0)
    parser.add_argument(
        "--integration-intervals",
        type=int,
        default=1000,
        help="Number of one-step intervals used to unfold each trajectory.",
    )
    parser.add_argument("--ode-rtol", type=float, default=1e-12)
    parser.add_argument("--ode-atol", type=float, default=1e-12)
    parser.add_argument("--train-ratio", type=float, default=0.7)
    parser.add_argument("--dic-size", type=int, default=300)
    parser.add_argument(
        "--layer-sizes",
        type=int,
        nargs="+",
        default=[300, 300, 300],
    )
    parser.add_argument("--reg", type=float, default=0.1)
    parser.add_argument("--outer-epochs", type=int, default=30)
    parser.add_argument(
        "--inner-epochs",
        type=int,
        default=2,
        help="Matches the fixed inner train_psi epochs in the original notebook.",
    )
    parser.add_argument("--batch-size", type=int, default=50000)
    parser.add_argument("--lr", type=float, default=1e-5)
    parser.add_argument("--log-interval", type=int, default=10)
    parser.add_argument("--lr-decay-factor", type=float, default=0.8)
    parser.add_argument("--end-condition", type=float, default=1e-5)
    parser.add_argument(
        "--residual-form",
        default="projected_vlambda",
        choices=["projected_kv", "projected_vlambda"],
        help="projected_vlambda best matches the original pendulum TF notebook.",
    )
    parser.add_argument(
        "--training-policy",
        default="float32",
        choices=["float32", "float64", "mixed_float16"],
        help="Use float32 by default for stability on the small pendulum test.",
    )
    parser.add_argument("--analysis-dtype", default="float64")
    parser.add_argument("--gram-dtype", default="float64")
    parser.add_argument("--spectral-dtype", default="float64")
    parser.add_argument("--loss-mode", default="squared")
    parser.add_argument("--loss-epsilon", type=float, default=1e-12)
    parser.add_argument("--train-shuffle", action="store_true")
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip saving the eigenvalue scatter plot.",
    )
    return parser.parse_args()


def configure_tensorflow(device_choice, require_gpu):
    if device_choice == "cpu":
        tf.config.set_visible_devices([], "GPU")
        return "cpu"

    gpus = tf.config.list_physical_devices("GPU")
    if not gpus:
        if require_gpu:
            raise RuntimeError("GPU was requested but TensorFlow did not detect one.")
        print("No GPU detected, falling back to CPU.")
        return "cpu"

    for gpu in gpus:
        tf.config.experimental.set_memory_growth(gpu, True)
    return "gpu"


def load_solver_module(module_name):
    solver_module = importlib.import_module(module_name)
    return solver_module.KoopmanNN, solver_module.KoopmanSolver


def wrap_angle(theta):
    return (theta + np.pi) % (2.0 * np.pi) - np.pi


def pendulum_rhs(_t, state):
    return np.array([state[1], -math.sin(state[0])], dtype=np.float64)


def build_pendulum_dataset(
    delta_t,
    grid_angle_count,
    grid_velocity_count,
    velocity_limit,
    integration_intervals,
    ode_rtol,
    ode_atol,
):
    x1 = np.linspace(-np.pi, np.pi, grid_angle_count + 1, dtype=np.float64)
    x1 = x1 + (x1[1] - x1[0]) / 2.0
    x1 = x1[:-1]
    x2 = np.linspace(-velocity_limit, velocity_limit, grid_velocity_count, dtype=np.float64)

    grid_theta, grid_omega = np.meshgrid(x1[:-1], x2)
    initial_conditions = np.column_stack([grid_theta.ravel(), grid_omega.ravel()])
    t_eval = np.linspace(1e-6, delta_t, integration_intervals + 1, dtype=np.float64)

    data_all = np.zeros((initial_conditions.shape[0], t_eval.size, 2), dtype=np.float64)
    for idx, y0 in enumerate(initial_conditions):
        solution = solve_ivp(
            pendulum_rhs,
            (float(t_eval[0]), float(t_eval[-1])),
            y0.astype(np.float64, copy=False),
            t_eval=t_eval,
            rtol=ode_rtol,
            atol=ode_atol,
        )
        if not solution.success:
            raise RuntimeError(f"solve_ivp failed for trajectory {idx}: {solution.message}")
        trajectory = solution.y.T
        trajectory[:, 0] = wrap_angle(trajectory[:, 0])
        data_all[idx] = trajectory

    temp_x = data_all[:, :-1, :]
    temp_y = data_all[:, 1:, :]
    data_x = temp_x.reshape(-1, 2)
    data_y = temp_y.reshape(-1, 2)

    metadata = {
        "num_initial_conditions": int(initial_conditions.shape[0]),
        "integration_points": int(t_eval.size),
        "num_snapshot_pairs": int(data_x.shape[0]),
        "grid_angle_count": int(grid_angle_count),
        "grid_velocity_count": int(grid_velocity_count),
        "velocity_limit": float(velocity_limit),
        "delta_t": float(delta_t),
        "integration_intervals": int(integration_intervals),
    }
    return data_x, data_y, metadata


def split_train_valid(data_x, data_y, train_ratio):
    if not (0.0 < train_ratio < 1.0):
        raise ValueError("train_ratio must be between 0 and 1.")

    len_all = data_x.shape[0]
    split_idx = int(train_ratio * len_all)
    if split_idx < 1 or split_idx >= len_all - 1:
        raise ValueError("train_ratio results in an empty train or validation split.")

    data_x_train = data_x[:split_idx]
    data_x_valid = data_x[split_idx + 1 :]
    data_y_train = data_y[:split_idx]
    data_y_valid = data_y[split_idx + 1 :]
    return [data_x_train, data_y_train], [data_x_valid, data_y_valid]


def maybe_numpy(value):
    if hasattr(value, "numpy"):
        return value.numpy()
    return np.asarray(value)


def save_eigenvalue_plot(evalues, output_path):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(evalues.real, evalues.imag, s=12, alpha=0.7)
    phi = np.linspace(0.0, 2.0 * np.pi, 1000)
    ax.plot(np.cos(phi), np.sin(phi), color="red", linewidth=1.0)
    ax.set_xlabel("Real")
    ax.set_ylabel("Imag")
    ax.set_title("Pendulum Eigenvalues")
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    data_dir = output_dir / "data"
    checkpoint_dir = output_dir / "checkpoints"
    plot_dir = output_dir / "plots"
    for path in [output_dir, data_dir, checkpoint_dir, plot_dir]:
        path.mkdir(parents=True, exist_ok=True)

    np.random.seed(args.seed)
    tf.random.set_seed(args.seed)
    actual_device = configure_tensorflow(args.device, args.require_gpu)
    KoopmanNN, KoopmanSolver = load_solver_module(args.solver_module)

    start_time = time.time()
    data_x, data_y, dataset_meta = build_pendulum_dataset(
        delta_t=args.delta_t,
        grid_angle_count=args.grid_angle_count,
        grid_velocity_count=args.grid_velocity_count,
        velocity_limit=args.velocity_limit,
        integration_intervals=args.integration_intervals,
        ode_rtol=args.ode_rtol,
        ode_atol=args.ode_atol,
    )

    dataset_name = (
        f"data_pendulum_"
        f"{dataset_meta['num_initial_conditions']}.mat"
    )
    dataset_path = data_dir / dataset_name
    sio.savemat(dataset_path, {"DATA_X": data_x, "DATA_Y": data_y})
    data_train, data_valid = split_train_valid(data_x, data_y, args.train_ratio)

    basis_function = KoopmanNN(
        layer_sizes=args.layer_sizes,
        n_psi_train=args.dic_size - 3,
    )
    solver = KoopmanSolver(
        dic=basis_function,
        target_dim=int(data_x.shape[1]),
        reg=args.reg,
        training_policy=args.training_policy,
        analysis_dtype=args.analysis_dtype,
        gram_dtype=args.gram_dtype,
        spectral_dtype=args.spectral_dtype,
        residual_form=args.residual_form,
        loss_mode=args.loss_mode,
        loss_epsilon=args.loss_epsilon,
    )

    final_ckpt_prefix = checkpoint_dir / "final" / "ckpt"
    best_ckpt_prefix = checkpoint_dir / "best" / "ckpt"
    losses, val_losses, stop_flag, lr_history = solver.build(
        data_train=data_train,
        data_valid=data_valid,
        epochs=args.outer_epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        log_interval=args.log_interval,
        lr_decay_factor=args.lr_decay_factor,
        Nepoch=args.inner_epochs,
        end_condition=args.end_condition,
        train_shuffle=args.train_shuffle,
        checkpoint_path=str(final_ckpt_prefix),
        best_checkpoint_path=str(best_ckpt_prefix),
        save_best_only=True,
        resume=False,
        run_metadata={
            "experiment": args.experiment_tag,
            "solver_module": args.solver_module,
            "actual_device": actual_device,
        },
    )

    evalues = np.asarray(solver.eigenvalues.T)
    efuns = np.asarray(solver.eigenfunctions(data_x))
    kpm_modes = np.asarray(solver.compute_mode().T)
    psi_x = np.asarray(solver.get_Psi_X())
    psi_y = np.asarray(solver.get_Psi_Y())
    koopman_matrix = maybe_numpy(solver.K)

    outputs = {
        "efuns": efuns.astype(np.complex128),
        "evalues": evalues.astype(np.complex128),
        "kpm_modes": kpm_modes.astype(np.complex128),
        "N_dict": np.array([evalues.shape[0]], dtype=np.float64),
        "K": np.asarray(koopman_matrix, dtype=np.float64),
        "Psi_X": np.asarray(psi_x, dtype=np.float64),
        "Psi_Y": np.asarray(psi_y, dtype=np.float64),
        "loss": np.asarray(losses, dtype=np.float64),
        "val_loss": np.asarray(val_losses, dtype=np.float64),
        "lr_history": np.asarray(lr_history, dtype=np.float64),
    }
    outputs_path = output_dir / (
        f"{args.experiment_tag}_{dataset_meta['num_initial_conditions']}"
        f"traj_{args.residual_form}_{args.dic_size}basis.mat"
    )
    sio.savemat(outputs_path, outputs)

    if not args.skip_plots:
        save_eigenvalue_plot(evalues, plot_dir / "eigenvalues.png")

    summary = {
        "status": "ok",
        "solver_module": args.solver_module,
        "experiment_tag": args.experiment_tag,
        "actual_device": actual_device,
        "dataset_path": str(dataset_path),
        "outputs_path": str(outputs_path),
        "residual_form": args.residual_form,
        "dic_size": int(args.dic_size),
        "layer_sizes": [int(v) for v in args.layer_sizes],
        "batch_size": int(args.batch_size),
        "outer_epochs_requested": int(args.outer_epochs),
        "inner_epochs_requested": int(args.inner_epochs),
        "outer_epochs_completed": int(len(getattr(solver, "outer_history", []))),
        "stop_flag": int(stop_flag),
        "best_val_metric": float(getattr(solver, "best_val_metric", np.nan)),
        "best_outer_epoch": (
            None if getattr(solver, "best_outer_epoch", None) is None
            else int(solver.best_outer_epoch)
        ),
        "eigvec_cond": float(getattr(solver, "eigvec_cond", np.nan)),
        "final_loss_tail": [float(x) for x in np.asarray(losses, dtype=np.float64)[-5:]],
        "final_val_loss_tail": [float(x) for x in np.asarray(val_losses, dtype=np.float64)[-5:]],
        "max_abs_eigenvalue": float(np.max(np.abs(evalues))),
        "num_eigenvalues": int(evalues.shape[0]),
        "train_pairs": int(data_train[0].shape[0]),
        "valid_pairs": int(data_valid[0].shape[0]),
        "runtime_seconds": float(time.time() - start_time),
        "dataset": dataset_meta,
    }
    summary_path = output_dir / "summary.json"
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
