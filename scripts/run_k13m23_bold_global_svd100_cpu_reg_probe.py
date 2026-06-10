#!/usr/bin/env python3
"""CPU parameter probe for K13m23 BOLD global_svd100 P4.

This run does not overwrite the current pat40 mainline.  It keeps the current
standardized BOLD P4 setup but increases EDMD regularization from 0.001 to
0.01, so we can test whether the early K13m23 global_svd100 best checkpoint is
a weak-regularization artifact.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tmp.run_bold_stdT_roi_global_queue_20260520 as queue


MIN_EPOCHS = 80

queue.DATASETS = [{"stem": "k13m23", "group": "b"}]
queue.OBSERVABLES = ("global_svd100",)
queue.PARAM = {
    "lr": "1e-4",
    "reg": "0.01",
    "batch_size": "2000",
    "inner_epochs": "2",
    "tag": "stdT_l1e4_r1e2_b2000_i2_regprobe20260601",
}
queue.RUN_SUFFIX = "k13m23_global_regprobe"
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_cpu_k13m23_global_svd100_regprobe_20260601.log"
queue.BLOCK_EPOCHS = 20
queue.STALE_EPOCHS = 40
queue.MAX_TOTAL_EPOCHS = 400
queue.THREAD_ENV = {
    **queue.THREAD_ENV,
    "OMP_NUM_THREADS": "4",
    "OPENBLAS_NUM_THREADS": "4",
    "MKL_NUM_THREADS": "4",
    "NUMEXPR_NUM_THREADS": "4",
    "TF_NUM_INTRAOP_THREADS": "4",
    "TF_NUM_INTEROP_THREADS": "1",
    "CUDA_VISIBLE_DEVICES": "-1",
}


_original_build_command = queue.build_command


def _summary_path(group: str) -> Path:
    return queue.RUN_LOG_DIR / "bold_cpu_k13m23_global_svd100_regprobe_20260601_summary.json"


def _task_status(stem: str, observable: str) -> tuple[str, dict]:
    summary = queue.summarize_state(queue.read_state(stem, observable))
    epochs = int(summary["epochs"] or 0)
    stale = int(summary["stale_epochs"] or 0)
    if epochs <= 0:
        return "not_started", summary
    if epochs >= MIN_EPOCHS and stale >= queue.STALE_EPOCHS:
        return "stale_done", summary
    if epochs >= queue.MAX_TOTAL_EPOCHS:
        return "max_epoch_reached", summary
    return "running", summary


def _cpu_build_command(stem: str, observable: str, obs_file: Path, epochs: int, resume: bool) -> list[str]:
    command = _original_build_command(stem, observable, obs_file, epochs, resume=resume)
    cpu_command: list[str] = []
    skip_next = False
    for token in command:
        if skip_next:
            skip_next = False
            continue
        if token == "--selected-device":
            cpu_command.extend(["--selected-device", "cpu"])
            skip_next = True
            continue
        if token == "--require-gpu":
            continue
        cpu_command.append(token)
    return cpu_command


queue.summary_path = _summary_path
queue.task_status = _task_status
queue.build_command = _cpu_build_command


if __name__ == "__main__":
    if "--group" not in sys.argv:
        sys.argv.extend(["--group", "b"])
    queue.main()
