import argparse
import json
import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from numpy import linalg as la
from scipy import linalg as sla
from tqdm import tqdm

_Q_MAT = None
_A_MAT = None
_L_MAT = None
_G_MAT = None


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot the pendulum pseudospectrum from exported ResKoopNet solver outputs."
    )
    parser.add_argument(
        "--input-mat",
        required=True,
        help="MAT file exported by run_pendulum_reskoopnet_solver3.py",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for pseudospectrum outputs.",
    )
    parser.add_argument("--x-min", type=float, default=-1.5)
    parser.add_argument("--x-max", type=float, default=1.5)
    parser.add_argument("--step", type=float, default=0.05)
    parser.add_argument(
        "--epsilon-threshold",
        type=float,
        default=0.25,
        help="Matches the notebook threshold used for continuous-spectrum visualization.",
    )
    parser.add_argument(
        "--lower-bound",
        type=float,
        default=1e-16,
        help="Clamp exact zeros before log scaling.",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, (os.cpu_count() or 2) - 1),
        help="Number of worker processes for pseudospectrum evaluation.",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=96,
        help="Number of z-points per worker task.",
    )
    return parser.parse_args()


def _init_worker(q_mat, a_mat, l_mat, g_mat):
    global _Q_MAT, _A_MAT, _L_MAT, _G_MAT
    _Q_MAT = q_mat
    _A_MAT = a_mat
    _L_MAT = l_mat
    _G_MAT = g_mat


def _eval_chunk(z_chunk):
    res_list = []
    for z_val in z_chunk:
        l_adjusted = (
            _L_MAT
            - z_val * np.conj(_A_MAT.T)
            - np.conj(z_val) * _A_MAT
            + (np.abs(z_val) ** 2) * _G_MAT
        )
        curr_ev = sla.eigvals(l_adjusted, overwrite_a=True, check_finite=False)
        min_abs_idx = int(np.argmin(np.abs(curr_ev)))
        curr_res = np.emath.sqrt(np.real(curr_ev[min_abs_idx]))
        res_list.append(curr_res)
    return np.asarray(res_list)


def koop_pseudo_spec_qr(px, py, w, z_pts, workers=1, chunk_size=96):
    q_mat, r_mat = la.qr(np.sqrt(w) * px)
    c1 = (np.sqrt(w) * py) @ la.inv(r_mat)
    l_mat = np.conj(c1.T) @ c1
    g_mat = np.eye(px.shape[1], dtype=np.float64)
    a_mat = np.conj(q_mat.T) @ c1

    if workers <= 1:
        _init_worker(q_mat, a_mat, l_mat, g_mat)
        return _eval_chunk(z_pts)

    z_chunks = [
        z_pts[start : start + chunk_size]
        for start in range(0, len(z_pts), chunk_size)
    ]
    results = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=_init_worker,
        initargs=(q_mat, a_mat, l_mat, g_mat),
    ) as executor:
        for chunk_res in tqdm(
            executor.map(_eval_chunk, z_chunks),
            total=len(z_chunks),
            desc="pseudospectrum",
            unit="chunk",
        ):
            results.append(chunk_res)
    return np.concatenate(results, axis=0)


def plot_notebook_style(x_mesh, y_mesh, res_grid, n_dict, output_path):
    fig, ax = plt.subplots(figsize=(4, 4))
    v = np.array([0.25, 1e-64], dtype=np.float64)
    levels = np.log10(1.0 / v)
    field = np.log10(1.0 / np.real(res_grid))
    ax.contourf(x_mesh, y_mesh, field, levels=levels)
    ax.contour(x_mesh, y_mesh, field, levels=levels)
    circle_phi = np.arange(0.0, 2.0 * np.pi, 2.0 * np.pi / 1000.0)
    ax.plot(np.cos(circle_phi), np.sin(circle_phi), color="red")
    ax.set_title(f"N_K = {n_dict}")
    ax.grid(True)
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_threshold_style(x_mesh, y_mesh, res_grid, n_dict, threshold, output_path):
    fig, ax = plt.subplots(figsize=(4, 4))
    mask = np.real(res_grid) <= threshold
    img = ax.imshow(
        mask.astype(float),
        origin="lower",
        extent=[x_mesh.min(), x_mesh.max(), y_mesh.min(), y_mesh.max()],
        cmap=matplotlib.colors.ListedColormap(["white", "lightgreen"]),
        interpolation="nearest",
        aspect="equal",
    )
    del img
    circle_phi = np.arange(0.0, 2.0 * np.pi, 2.0 * np.pi / 1000.0)
    ax.plot(np.cos(circle_phi), np.sin(circle_phi), color="red", linewidth=1.2)
    ax.set_title(f"Continuous Spectrum, N_K = {n_dict}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    input_mat = Path(args.input_mat).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    data = sio.loadmat(input_mat)
    psi_x = np.asarray(data["Psi_X"], dtype=np.float64)
    psi_y = np.asarray(data["Psi_Y"], dtype=np.float64)
    n_dict = int(np.asarray(data["N_dict"]).ravel()[0])

    x_pts = np.arange(args.x_min, args.x_max + args.step, args.step, dtype=np.float64)
    y_pts = x_pts.copy()
    x_mesh, y_mesh = np.meshgrid(x_pts, y_pts)
    z_pts = (x_mesh + 1j * y_mesh).ravel()
    weights = np.ones((psi_x.shape[0], 1), dtype=np.float64)

    res = koop_pseudo_spec_qr(
        psi_x,
        psi_y,
        weights,
        z_pts,
        workers=max(1, args.workers),
        chunk_size=max(1, args.chunk_size),
    ).reshape(x_mesh.shape)
    zero_real_indices = np.nonzero(np.real(res) == 0.0)
    if len(zero_real_indices[0]) > 0:
        res[zero_real_indices] = args.lower_bound + 1j * np.imag(res[zero_real_indices])

    notebook_plot = output_dir / "pseudospectrum_notebook_style.png"
    threshold_plot = output_dir / "pseudospectrum_threshold_style.png"
    plot_notebook_style(x_mesh, y_mesh, res, n_dict, notebook_plot)
    plot_threshold_style(
        x_mesh,
        y_mesh,
        res,
        n_dict,
        args.epsilon_threshold,
        threshold_plot,
    )

    summary = {
        "input_mat": str(input_mat),
        "n_dict": n_dict,
        "psi_x_shape": list(psi_x.shape),
        "psi_y_shape": list(psi_y.shape),
        "grid_size": [int(x_mesh.shape[0]), int(x_mesh.shape[1])],
        "epsilon_threshold": float(args.epsilon_threshold),
        "res_min_real": float(np.min(np.real(res))),
        "res_max_real": float(np.max(np.real(res))),
        "notebook_plot": str(notebook_plot),
        "threshold_plot": str(threshold_plot),
    }
    summary_path = output_dir / "pseudospectrum_summary.json"
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
