#!/usr/bin/env python3
"""Build missing K13m23 global/gsvd BOLD SVD observables without MATLAB.

This is a targeted fallback for the 2026-06-01 K13m23 backfill: MATLAB batch
startup was failing before the normal P3 builder could run, while P4 only needs
the canonical BOLD observable MAT fields, especially ``obs`` and session
metadata.  The computation mirrors ``script_build_one_cfg_bold_observables.m``
for ``global_svd100`` and ``gsvd100_ds`` closely enough for downstream P4/P7.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import h5py
import hdf5storage
import numpy as np


@dataclass(frozen=True)
class DatasetSpec:
    stem: str
    dataset_id: str
    raw_root: Path
    sessions: tuple[int, ...]


K13M23 = DatasetSpec(
    stem="k13m23",
    dataset_id="K13.m23",
    raw_root=Path("/mnt/e/DataPons/K13.m23"),
    sessions=tuple(range(27, 37)),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processed-root", type=Path, default=Path("/mnt/e/DataPons_processed"))
    parser.add_argument("--modes", nargs="+", default=["global_svd100", "gsvd100_ds"])
    parser.add_argument("--n-components", type=int, default=100)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def decode_matlab_char(dataset: h5py.Dataset) -> str:
    arr = np.asarray(dataset)
    return "".join(chr(int(x)) for x in arr.ravel(order="F") if int(x) != 0)


def read_roi_session(path: Path) -> tuple[np.ndarray, np.ndarray, list[str], float]:
    with h5py.File(path, "r") as handle:
        refs = np.asarray(handle["roiTs"]).ravel()
        data_parts: list[np.ndarray] = []
        coord_parts: list[np.ndarray] = []
        labels: list[str] = []
        dx_values: list[float] = []
        for ref in refs:
            roi = handle[ref]
            name = decode_matlab_char(roi["name"])
            # MATLAB stores these transposed in HDF5.  The normal MATLAB loader
            # sees dat as time x voxel and coords as voxel x xyz.
            data = np.asarray(roi["dat"], dtype=np.float64).T
            coords = np.asarray(roi["coords"], dtype=np.float64).T
            dx = float(np.asarray(roi["dx"]).reshape(-1)[0])
            if coords.shape[0] != data.shape[1]:
                raise ValueError(f"{path}: ROI {name} coords/data mismatch {coords.shape} vs {data.shape}")
            data_parts.append(data)
            coord_parts.append(coords)
            labels.extend(f"{name}_v{idx:04d}" for idx in range(1, data.shape[1] + 1))
            dx_values.append(dx)

    if not dx_values or not np.allclose(dx_values, dx_values[0], atol=1e-12, rtol=0):
        raise ValueError(f"{path}: inconsistent ROI dx values")
    return np.concatenate(data_parts, axis=1), np.concatenate(coord_parts, axis=0), labels, dx_values[0]


def load_dataset(spec: DatasetSpec) -> dict[str, object]:
    data_parts: list[np.ndarray] = []
    session_lengths: list[int] = []
    session_dx: list[float] = []
    labels_ref: list[str] | None = None
    coords_ref: np.ndarray | None = None
    source_files: list[str] = []

    for session_id in spec.sessions:
        path = spec.raw_root / "roits" / f"{spec.stem}_{session_id:04d}_roits.mat"
        if not path.is_file():
            raise FileNotFoundError(path)
        print(f"Loading {path}")
        data, coords, labels, dx = read_roi_session(path)
        if labels_ref is None:
            labels_ref = labels
            coords_ref = coords
        else:
            if labels != labels_ref:
                raise ValueError(f"{path}: ROI/voxel labels changed across sessions")
            if coords_ref is None or coords.shape != coords_ref.shape or not np.allclose(coords, coords_ref, equal_nan=True):
                raise ValueError(f"{path}: voxel coordinates changed across sessions")
        data_parts.append(data)
        session_lengths.append(data.shape[0])
        session_dx.append(dx)
        source_files.append(str(path))
        print(f"  session {session_id}: data={data.shape}, dx={dx:g}")

    data_all = np.concatenate(data_parts, axis=0)
    starts = np.cumsum([1] + session_lengths[:-1])
    ends = np.cumsum(session_lengths)
    return {
        "data": data_all,
        "coords": coords_ref,
        "labels": labels_ref or [],
        "source_files": source_files,
        "session_ids": np.asarray(spec.sessions, dtype=np.float64).reshape(-1, 1),
        "session_lengths": np.asarray(session_lengths, dtype=np.float64).reshape(-1, 1),
        "session_dx": np.asarray(session_dx, dtype=np.float64).reshape(-1, 1),
        "session_start_idx": starts.astype(np.float64).reshape(-1, 1),
        "session_end_idx": ends.astype(np.float64).reshape(-1, 1),
        "border_idx": ends[:-1].astype(np.float64).reshape(-1, 1),
        "dx": float(session_dx[0]),
    }


def stable_unique_rows(coords: np.ndarray) -> np.ndarray:
    _, first_idx = np.unique(coords, axis=0, return_index=True)
    return np.sort(first_idx)


def svd_scores_via_time_covariance(x: np.ndarray, n_components: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=np.float64)
    mu = np.mean(x, axis=0)
    xc = x - mu
    gram = xc @ xc.T
    vals, vecs = np.linalg.eigh(gram)
    order = np.argsort(vals)[::-1]
    vals = np.maximum(vals[order], 0.0)
    vecs = vecs[:, order]
    n_comp = min(n_components, x.shape[0] - 1, x.shape[1], vals.size)
    singular = np.sqrt(vals[:n_comp])
    score = vecs[:, :n_comp] * singular.reshape(1, -1)
    latent = vals[:n_comp] / max(x.shape[0] - 1, 1)
    return score, latent, mu


def build_mode(dataset: dict[str, object], mode: str, n_components: int) -> dict[str, object]:
    x = np.asarray(dataset["data"], dtype=np.float64)
    coords = np.asarray(dataset["coords"], dtype=np.float64)
    labels = list(dataset["labels"])
    pre_svd: dict[str, object] = {}

    if mode == "gsvd100_ds":
        keep = stable_unique_rows(coords)
        x = x[:, keep]
        labels = [labels[i] for i in keep]
        pre_svd = {
            "original_n_variables": float(coords.shape[0]),
            "keep_variable_idx": (keep + 1).astype(np.float64).reshape(-1, 1),
            "dedup_by_coords": True,
            "session_center": True,
            "session_detrend": False,
            "variable_zscore": False,
            "dedup_removed_variables": float(coords.shape[0] - keep.size),
        }
        for start, end in zip(np.ravel(dataset["session_start_idx"]).astype(int), np.ravel(dataset["session_end_idx"]).astype(int)):
            sl = slice(start - 1, end)
            x[sl, :] -= np.mean(x[sl, :], axis=0, keepdims=True)
    elif mode != "global_svd100":
        raise ValueError(f"Unsupported mode: {mode}")

    print(f"Computing {mode}: matrix={x.shape}")
    score, latent, mu = svd_scores_via_time_covariance(x, n_components)
    if mode == "gsvd100_ds":
        total_var_x = x - np.mean(x, axis=0, keepdims=True)
        total_latent = float(np.sum(total_var_x * total_var_x) / max(x.shape[0] - 1, 1))
        explained = 100.0 * latent / total_latent if total_latent > 0 else np.zeros_like(latent)
    else:
        explained = 100.0 * latent / np.sum(latent) if np.sum(latent) > 0 else np.zeros_like(latent)

    observable_labels = [f"{mode}{idx:03d}" for idx in range(1, score.shape[1] + 1)]
    observable_info = {
        "observable_idx": np.arange(1, score.shape[1] + 1, dtype=np.float64).reshape(-1, 1),
        "source": np.asarray([mode] * score.shape[1], dtype=object).reshape(-1, 1),
        "observable_label": np.asarray(observable_labels, dtype=object).reshape(-1, 1),
        "latent": latent.astype(np.float64).reshape(-1, 1),
        "explained": explained.astype(np.float64).reshape(-1, 1),
    }
    params = {
        "source": "data",
        "observable_branch": mode,
        "mode": "svd",
        "precision": "single",
        "n_components": float(n_components),
    }
    if mode == "gsvd100_ds":
        params.update(
            {
                "dedup_by_coords": True,
                "session_center": True,
                "session_detrend": False,
                "variable_zscore": False,
                "explained_use_total_variance": True,
            }
        )
    model = {
        "type": "svd",
        "center": True,
        "mu": mu.astype(np.float32),
        "latent": latent.astype(np.float64).reshape(-1, 1),
        "explained": explained.astype(np.float64).reshape(-1, 1),
        "pre_svd": pre_svd,
        "input_variable_labels": np.asarray(labels, dtype=object).reshape(-1, 1),
    }
    return {
        "obs": score.astype(np.float32),
        "obs_info": observable_info,
        "observable_labels": np.asarray(observable_labels, dtype=object).reshape(-1, 1),
        "params": params,
        "model": model,
    }


def save_artifact(path: Path, spec: DatasetSpec, dataset: dict[str, object], mode_payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dx = float(dataset["dx"])
    o_struct = {
        "data": mode_payload["obs"],
        "observable_info": mode_payload["obs_info"],
        "observable_labels": mode_payload["observable_labels"],
        "params": mode_payload["params"],
        "model": mode_payload["model"],
        "source_data_kind": "data",
        "session_ids": dataset["session_ids"],
        "session_lengths": dataset["session_lengths"],
        "session_dx": dataset["session_dx"],
        "session_start_idx": dataset["session_start_idx"],
        "session_end_idx": dataset["session_end_idx"],
        "border_idx": dataset["border_idx"],
        "dx": dx,
        "fs": 1.0 / dx,
        "dataset_id": spec.dataset_id,
        "file_stem": spec.stem,
    }
    valid_x: list[np.ndarray] = []
    valid_y: list[np.ndarray] = []
    sess_idx: list[np.ndarray] = []
    sess_id: list[np.ndarray] = []
    for idx, (sid, start, end) in enumerate(
        zip(
            np.ravel(dataset["session_ids"]),
            np.ravel(dataset["session_start_idx"]).astype(int),
            np.ravel(dataset["session_end_idx"]).astype(int),
        ),
        start=1,
    ):
        ix = np.arange(start, end, dtype=np.float64)
        valid_x.append(ix)
        valid_y.append(ix + 1)
        sess_idx.append(np.full(ix.shape, idx, dtype=np.float64))
        sess_id.append(np.full(ix.shape, sid, dtype=np.float64))
    snapshot_valid_idx_x = np.concatenate(valid_x).reshape(-1, 1)
    snapshot_valid_idx_y = np.concatenate(valid_y).reshape(-1, 1)
    snapshot_session_idx = np.concatenate(sess_idx).reshape(-1, 1)
    snapshot_session_id = np.concatenate(sess_id).reshape(-1, 1)
    snap = {
        "X": mode_payload["obs"][snapshot_valid_idx_x.astype(int).ravel() - 1, :],
        "Y": mode_payload["obs"][snapshot_valid_idx_y.astype(int).ravel() - 1, :],
        "valid_idx_x": snapshot_valid_idx_x,
        "valid_idx_y": snapshot_valid_idx_y,
        "session_idx": snapshot_session_idx,
        "session_id": snapshot_session_id,
        "lag": 1.0,
        "drop_initial": 0.0,
        "n_pairs": float(snapshot_valid_idx_x.size),
    }
    cfg_saved = {
        "dataset_id": spec.dataset_id,
        "file_stem": spec.stem,
        "raw_data_root": str(spec.raw_root),
    }
    payload = {
        "obs": mode_payload["obs"],
        "obs_info": mode_payload["obs_info"],
        "dx": dx,
        "dt": dx,
        "fs": 1.0 / dx,
        "sampling_period": dx,
        "sample_period": dx,
        "sampling_frequency": 1.0 / dx,
        "session_ids": dataset["session_ids"],
        "session_lengths": dataset["session_lengths"],
        "session_dx": dataset["session_dx"],
        "session_start_idx": dataset["session_start_idx"],
        "session_end_idx": dataset["session_end_idx"],
        "border_idx": dataset["border_idx"],
        "snapshot_lag_saved": 1.0,
        "snapshot_valid_idx_x": snapshot_valid_idx_x,
        "snapshot_valid_idx_y": snapshot_valid_idx_y,
        "snapshot_session_idx": snapshot_session_idx,
        "snapshot_session_id": snapshot_session_id,
        "O": o_struct,
        "snap": snap,
        "params": mode_payload["params"],
        "pre": {
            "demean": False,
            "zscore": False,
            "detrend": False,
            "notch_hz": np.empty((0, 0), dtype=np.float64),
            "bandpass_hz": np.empty((0, 0), dtype=np.float64),
            "save_filtered": False,
        },
        "cfg_saved": cfg_saved,
    }
    print(f"Saving {path}")
    hdf5storage.savemat(
        path,
        payload,
        format="7.3",
        matlab_compatible=True,
        store_python_metadata=False,
    )


def main() -> None:
    args = parse_args()
    spec = K13M23
    out_dir = args.processed_root / spec.stem / "pipeline3_bold_observables"
    missing_modes = [
        mode for mode in args.modes
        if args.overwrite or not (out_dir / f"{spec.dataset_id}_bold_observables_{mode}.mat").is_file()
    ]
    if not missing_modes:
        print("All requested K13m23 BOLD observable files already exist.")
        return
    dataset = load_dataset(spec)
    for mode in missing_modes:
        payload = build_mode(dataset, mode, args.n_components)
        save_artifact(out_dir / f"{spec.dataset_id}_bold_observables_{mode}.mat", spec, dataset, payload)
    print("Done.")


if __name__ == "__main__":
    main()
