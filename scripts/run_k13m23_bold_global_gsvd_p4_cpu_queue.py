#!/usr/bin/env python3
"""Run K13m23 BOLD global/gsvd P4 on CPU.

This keeps the same run labels/parameters as the standardized pat40 BOLD P4
queue, but removes the GPU requirement so it can run while BLP training uses
the GPU.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tmp.run_bold_stdT_roi_global_pat40_queue_20260522 as wrapped


queue = wrapped.queue
queue.DATASETS = [{"stem": "k13m23", "group": "b"}]
queue.OBSERVABLES = ("global_svd100", "gsvd100_ds")
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_cpu_k13m23_global_gsvd_pat40_20260601.log"
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
    return queue.RUN_LOG_DIR / "bold_cpu_k13m23_global_gsvd_pat40_20260601_summary.json"


def _cpu_build_command(stem: str, observable: str, obs_file: Path, epochs: int, resume: bool) -> list[str]:
    command = _original_build_command(stem, observable, obs_file, epochs, resume=resume)
    cpu_command: list[str] = []
    skip_next = False
    for idx, token in enumerate(command):
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
queue.build_command = _cpu_build_command


if __name__ == "__main__":
    if "--group" not in sys.argv:
        sys.argv.extend(["--group", "b"])
    queue.main()
