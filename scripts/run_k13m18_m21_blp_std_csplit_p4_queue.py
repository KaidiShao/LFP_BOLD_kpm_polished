#!/usr/bin/env python3
"""Run BLP P4 standardized complex-split training for K13m18-K13m21.

This is the current standardized complex-split BLP P4 branch:

  complex_split observable
  projected_vlambda residual
  std_complex_pair normalization
  block-wise GPU training until stale/material plateau
  compact summary export only

It mirrors the accepted 2026-05-22 standardized branch used by P5/P8/P10,
but restricts the target set to the four new K13 datasets.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from contextlib import contextmanager
from pathlib import Path

try:
    import fcntl
except ImportError:  # pragma: no cover - intended for WSL/Linux.
    fcntl = None


REPO = Path("/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished")
PYTHON = Path("/home/kdshao/anaconda3/bin/python")
RUNNER = REPO / "python_scripts" / "autodl" / "run_autodl_reskoopnet_mlp.py"
SOLVER_DIR = REPO / "python_scripts" / "autodl"
RUN_LOG_DIR = REPO / "tmp" / "run_logs"
MAIN_LOG = RUN_LOG_DIR / "blp_std_csplit_k13m18_m21_p4_queue_20260601.log"
LOCK_PATH = RUN_LOG_DIR / "blp_std_csplit_k13m18_m21_p4_queue_20260601.lock"

DATA_ROOT = Path("/mnt/e/DataPons_processed")
RESULT_ROOT = Path("/mnt/e/autodl_results_new")
RUN_SUFFIX = "20260522_pat40_allblp"

DATASETS = ("k13m18", "k13m19", "k13m20", "k13m21")
OBSERVABLE = "complex_split"

PARAM = {
    "lr": "1e-4",
    "reg": "0.001",
    "batch_size": "2000",
    "inner_epochs": "2",
    "tag": "stdComplexPair_l1e4_r1e3_b2000_i2_pat40",
}

BLOCK_EPOCHS = 40
MIN_EPOCHS = 80
STALE_EPOCHS = 40
MAX_TOTAL_EPOCHS = 1000
MATERIAL_MIN_REL_IMPROVEMENT = 0.01
MATERIAL_MIN_ABS_IMPROVEMENT = 1e-8
MATERIAL_PATIENCE_EPOCHS = 40
DONE_STATUSES = {"stale_done", "delta_done", "max_epoch_reached"}

THREAD_ENV = {
    "CUDA_VISIBLE_DEVICES": "0",
    "TF_FORCE_GPU_ALLOW_GROWTH": "true",
    "OMP_NUM_THREADS": "2",
    "OPENBLAS_NUM_THREADS": "2",
    "MKL_NUM_THREADS": "2",
    "TF_NUM_INTRAOP_THREADS": "2",
    "TF_NUM_INTEROP_THREADS": "1",
}


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(message: str) -> None:
    RUN_LOG_DIR.mkdir(parents=True, exist_ok=True)
    line = f"[{now()}] {message}"
    print(line, flush=True)
    with MAIN_LOG.open("a", encoding="utf-8") as handle:
        handle.write(line + "\n")


@contextmanager
def training_lock(stem: str):
    RUN_LOG_DIR.mkdir(parents=True, exist_ok=True)
    with LOCK_PATH.open("a+", encoding="utf-8") as handle:
        if fcntl is not None:
            log(f"waiting BLP std-csplit lock: {stem}")
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
            handle.seek(0)
            handle.truncate()
            handle.write(f"{now()} {stem}\n")
            handle.flush()
        try:
            yield
        finally:
            if fcntl is not None:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
                log(f"released BLP std-csplit lock: {stem}")


def output_parent(stem: str) -> Path:
    return RESULT_ROOT / stem / "mlp" / "outputs"


def checkpoint_parent(stem: str) -> Path:
    return RESULT_ROOT / stem / "mlp" / "checkpoints"


def log_parent(stem: str) -> Path:
    return RESULT_ROOT / stem / "mlp" / "logs"


def experiment_name(stem: str) -> str:
    return f"blp_vlambda_{OBSERVABLE}_{PARAM['tag']}_{RUN_SUFFIX}_{stem}_seed1234"


def run_label(stem: str) -> str:
    return f"mlp_obs_{experiment_name(stem)}_projected_vlambda_{OBSERVABLE}"


def state_path(stem: str) -> Path:
    return checkpoint_parent(stem) / run_label(stem) / "final" / "training_state.json"


def summary_path() -> Path:
    return RUN_LOG_DIR / "blp_std_csplit_k13m18_m21_p4_queue_20260601_summary.json"


def find_observable_file(stem: str) -> Path | None:
    obs_dir = DATA_ROOT / stem / "pipeline1_reskoopnet_dictionary"
    if not obs_dir.exists():
        return None
    matches = sorted(obs_dir.glob(f"{stem}_*_g2_{OBSERVABLE}_single.mat"))
    matches = [
        path
        for path in matches
        if ".corrupt_" not in path.name and path.stat().st_size > 0
    ]
    return matches[0] if matches else None


def read_state(stem: str) -> dict | None:
    path = state_path(stem)
    if not path.exists():
        return None
    state = json.loads(path.read_text(encoding="utf-8"))
    outer = state.get("outer_history") or []
    vals = [
        float(row["val_metric"])
        for row in outer
        if row.get("val_metric") is not None
    ]
    state["_first_val_metric"] = vals[0] if vals else None
    state["_final_val_metric"] = vals[-1] if vals else None
    state["_material_best"] = material_best_summary(outer)
    return state


def material_best_summary(outer: list[dict]) -> dict:
    material_best_val = None
    material_best_epoch = None
    best_val = None
    best_epoch = None
    for idx, row in enumerate(outer):
        val = row.get("val_metric")
        if val is None:
            continue
        val = float(val)
        epoch = int(row.get("outer_epoch") or row.get("epoch") or (idx + 1))
        if best_val is None or val < best_val:
            best_val = val
            best_epoch = epoch
        if material_best_val is None:
            material_best_val = val
            material_best_epoch = epoch
            continue
        rel_drop = (material_best_val - val) / max(abs(material_best_val), 1e-30)
        abs_drop = material_best_val - val
        if rel_drop >= MATERIAL_MIN_REL_IMPROVEMENT and abs_drop >= MATERIAL_MIN_ABS_IMPROVEMENT:
            material_best_val = val
            material_best_epoch = epoch
    last_epoch = 0
    if outer:
        last = outer[-1]
        last_epoch = int(last.get("outer_epoch") or last.get("epoch") or len(outer))
    return {
        "best_val": best_val,
        "best_epoch": best_epoch,
        "material_best_val": material_best_val,
        "material_best_epoch": material_best_epoch,
        "material_stale_epochs": (
            last_epoch - material_best_epoch
            if last_epoch and material_best_epoch is not None
            else 0
        ),
    }


def summarize_state(state: dict | None) -> dict:
    if not state:
        return {
            "epochs": 0,
            "best_val": None,
            "best_epoch": None,
            "final_val": None,
            "stale_epochs": 0,
            "final_over_first": None,
            "best_over_first": None,
        }
    first = state.get("_first_val_metric")
    final = state.get("_final_val_metric")
    best = state.get("best_val_metric")
    material = state.get("_material_best") or {}
    return {
        "epochs": int(state.get("completed_outer_epochs") or 0),
        "best_val": best,
        "best_epoch": state.get("best_outer_epoch"),
        "final_val": final,
        "stale_epochs": int(state.get("epochs_since_best") or 0),
        "material_best_val": material.get("material_best_val"),
        "material_best_epoch": material.get("material_best_epoch"),
        "material_stale_epochs": int(material.get("material_stale_epochs") or 0),
        "material_min_rel_improvement": MATERIAL_MIN_REL_IMPROVEMENT,
        "material_min_abs_improvement": MATERIAL_MIN_ABS_IMPROVEMENT,
        "final_over_first": (
            final / first
            if first not in (None, 0) and final is not None
            else None
        ),
        "best_over_first": (
            best / first
            if first not in (None, 0) and best is not None
            else None
        ),
    }


def task_status(stem: str) -> tuple[str, dict]:
    summary = summarize_state(read_state(stem))
    epochs = int(summary["epochs"] or 0)
    stale = int(summary["stale_epochs"] or 0)
    if epochs <= 0:
        return "not_started", summary
    if (
        epochs >= MIN_EPOCHS
        and int(summary.get("material_stale_epochs") or 0) >= MATERIAL_PATIENCE_EPOCHS
    ):
        return "delta_done", summary
    if epochs >= MIN_EPOCHS and stale >= STALE_EPOCHS:
        return "stale_done", summary
    if epochs >= MAX_TOTAL_EPOCHS:
        return "max_epoch_reached", summary
    return "running", summary


def write_summary(records: list[dict]) -> None:
    path = summary_path()
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
    tmp.replace(path)


def build_command(stem: str, obs_file: Path, epochs: int, resume: bool) -> list[str]:
    command = [
        str(PYTHON),
        str(RUNNER),
        "--project-root",
        str(REPO),
        "--solver-dir",
        str(SOLVER_DIR),
        "--data-root",
        str(DATA_ROOT / stem),
        "--data-subdir",
        "pipeline1_reskoopnet_dictionary",
        "--output-parent",
        str(output_parent(stem)),
        "--checkpoint-parent",
        str(checkpoint_parent(stem)),
        "--log-parent",
        str(log_parent(stem)),
        "--run-name-base",
        "mlp_obs",
        "--experiment-name",
        experiment_name(stem),
        "--selected-device",
        "gpu",
        "--require-gpu",
        "--solver-name",
        "resdmd_batch_mixedgpu",
        "--residual-form",
        "projected_vlambda",
        "--training-policy",
        "float32",
        "--analysis-dtype",
        "float64",
        "--gram-dtype",
        "float64",
        "--spectral-dtype",
        "float64",
        "--dataset-stem",
        stem,
        "--observable-mode",
        OBSERVABLE,
        "--data-filename",
        obs_file.name,
        "--file-type",
        ".mat",
        "--field-name",
        "obs",
        "--n-psi-train",
        "100",
        "--train-ratio",
        "0.7",
        "--reg",
        PARAM["reg"],
        "--rounds",
        "1",
        "--epochs",
        str(epochs),
        "--batch-size",
        PARAM["batch_size"],
        "--lr",
        PARAM["lr"],
        "--log-interval",
        "1",
        "--lr-decay-factor",
        "0.8",
        "--inner-epochs",
        PARAM["inner_epochs"],
        "--end-condition",
        "1e-9",
        "--chunk-size",
        "5000",
        "--layer-sizes",
        "100",
        "100",
        "100",
        "--train-shuffle",
        "--spectral-sync-mode",
        "pre_only",
        "--seed",
        "1234",
        "--export-mode",
        "summary_only",
        "--skip-diagnostic-pdf",
        "--standardize-data",
        "--standardize-eps",
        "1e-8",
        "--standardize-mode",
        "std_complex_pair",
    ]
    if resume:
        command.extend(["--resume", "--resume-mode", "final"])
    else:
        command.append("--fresh-checkpoints")
    return command


def latest_logged_block_index(stem: str) -> int:
    label = run_label(stem)
    max_index = 0
    for path in RUN_LOG_DIR.glob(f"{label}_block*.stdout.log"):
        suffix = path.name.removeprefix(f"{label}_block").removesuffix(".stdout.log")
        try:
            max_index = max(max_index, int(suffix))
        except ValueError:
            continue
    return max_index


def run_block(stem: str, obs_file: Path, block_index: int, resume: bool) -> int:
    label = run_label(stem)
    stdout_path = RUN_LOG_DIR / f"{label}_block{block_index:03d}.stdout.log"
    stderr_path = RUN_LOG_DIR / f"{label}_block{block_index:03d}.stderr.log"
    command = build_command(stem, obs_file, BLOCK_EPOCHS, resume=resume)
    env = {**os.environ, **THREAD_ENV}
    with training_lock(stem):
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
            "w", encoding="utf-8"
        ) as stderr:
            result = subprocess.run(
                command, cwd=str(REPO), env=env, stdout=stdout, stderr=stderr
            )
    log(
        f"block done: {stem} block={block_index} rc={result.returncode} "
        f"stdout={stdout_path.name}"
    )
    return int(result.returncode)


def make_tasks(selected: list[str] | None = None) -> list[dict]:
    stems = selected or list(DATASETS)
    tasks = []
    for stem in stems:
        obs_file = find_observable_file(stem)
        tasks.append({"stem": stem, "obs_file": str(obs_file) if obs_file else None})
    return tasks


def run_task(task: dict, records: list[dict]) -> None:
    stem = task["stem"]
    obs_file_text = task.get("obs_file")
    if not obs_file_text:
        record = {
            "stem": stem,
            "observable": OBSERVABLE,
            "status": "missing_observable",
            "param": PARAM,
            "label": run_label(stem),
        }
        records.append(record)
        write_summary(records)
        log(f"missing observable: {stem} {OBSERVABLE}")
        return

    obs_file = Path(obs_file_text)
    block_index = latest_logged_block_index(stem)
    while True:
        status, summary = task_status(stem)
        record = {
            "stem": stem,
            "observable": OBSERVABLE,
            "status": status,
            "param": PARAM,
            "label": run_label(stem),
            "state_path": str(state_path(stem)),
            "obs_file": str(obs_file),
            **summary,
        }
        records.append(record)
        write_summary(records)
        if status in DONE_STATUSES:
            log(f"task done: {stem} status={status} summary={summary}")
            return

        resume = status != "not_started"
        block_index += 1
        log(f"run start: {stem} block={block_index} resume={resume} before={summary}")
        rc = run_block(stem, obs_file, block_index, resume=resume)
        after_status, after_summary = task_status(stem)
        records.append(
            {
                "stem": stem,
                "observable": OBSERVABLE,
                "status": after_status if rc == 0 else "failed",
                "returncode": rc,
                "param": PARAM,
                "label": run_label(stem),
                "state_path": str(state_path(stem)),
                "obs_file": str(obs_file),
                **after_summary,
            }
        )
        write_summary(records)
        log(f"run done: {stem} rc={rc} after={after_summary}")
        if rc != 0:
            return


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", nargs="*", default=list(DATASETS))
    args = parser.parse_args()

    RUN_LOG_DIR.mkdir(parents=True, exist_ok=True)
    tasks = make_tasks(args.datasets)
    records: list[dict] = [
        {
            "event": "queue_start",
            "datasets": args.datasets,
            "observable": OBSERVABLE,
            "run_suffix": RUN_SUFFIX,
            "block_epochs": BLOCK_EPOCHS,
            "min_epochs": MIN_EPOCHS,
            "stale_epochs": STALE_EPOCHS,
            "material_patience_epochs": MATERIAL_PATIENCE_EPOCHS,
            "material_min_rel_improvement": MATERIAL_MIN_REL_IMPROVEMENT,
            "material_min_abs_improvement": MATERIAL_MIN_ABS_IMPROVEMENT,
            "max_total_epochs": MAX_TOTAL_EPOCHS,
            "param": PARAM,
            "tasks": tasks,
            "started_at": now(),
        }
    ]
    write_summary(records)
    log(f"BLP std-csplit P4 queue started; tasks={len(tasks)}")
    for task in tasks:
        run_task(task, records)
    records.append({"event": "queue_finished", "finished_at": now()})
    write_summary(records)
    log("BLP std-csplit P4 queue finished")


if __name__ == "__main__":
    main()
