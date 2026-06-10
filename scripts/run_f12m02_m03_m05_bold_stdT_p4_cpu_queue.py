#!/usr/bin/env python3
"""Run BOLD P4 standardized training for F12m02/F12m03/F12m05 on CPU.

This keeps the current BOLD stdT pat40 run labels/parameters, but forces CPU
execution so BLP P4 can own the GPU.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tmp.run_bold_stdT_roi_global_pat40_queue_20260522 as wrapped


queue = wrapped.queue

queue.DATASETS = [
    {"stem": "f12m02", "group": "b"},
    {"stem": "f12m03", "group": "b"},
    {"stem": "f12m05", "group": "b"},
]
queue.OBSERVABLES = ("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean")
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_cpu_f12m02_m03_m05_stdT_pat40_20260605.log"
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

DONE_STATUSES = {"stale_done", "delta_done", "max_epoch_reached"}
_original_build_command = queue.build_command


def _summary_path(group: str) -> Path:
    return queue.RUN_LOG_DIR / "bold_cpu_f12m02_m03_m05_stdT_pat40_20260605_summary.json"


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


def _latest_logged_block_index(stem: str, observable: str) -> int:
    label = queue.run_label(stem, observable)
    max_index = 0
    for path in queue.RUN_LOG_DIR.glob(f"{label}_block*.stdout.log"):
        suffix = path.name.removeprefix(f"{label}_block").removesuffix(".stdout.log")
        try:
            max_index = max(max_index, int(suffix))
        except ValueError:
            continue
    return max_index


def _run_task(group: str, task: dict, records: list[dict]) -> None:
    stem = task["stem"]
    observable = task["observable"]
    obs_file_text = task.get("obs_file")
    if not obs_file_text:
        record = {
            "stem": stem,
            "observable": observable,
            "status": "missing_observable",
            "param": queue.PARAM,
            "label": queue.run_label(stem, observable),
        }
        records.append(record)
        queue.write_summary(group, records)
        queue.log(f"group={group} missing observable: {stem} {observable}")
        return

    obs_file = Path(obs_file_text)
    block_index = _latest_logged_block_index(stem, observable)
    while True:
        status, summary = queue.task_status(stem, observable)
        record = {
            "stem": stem,
            "observable": observable,
            "status": status,
            "param": queue.PARAM,
            "label": queue.run_label(stem, observable),
            "state_path": str(queue.state_path(stem, observable)),
            "obs_file": str(obs_file),
            **summary,
        }
        records.append(record)
        queue.write_summary(group, records)
        if status in DONE_STATUSES:
            queue.log(f"group={group} task done: {stem} {observable} status={status} summary={summary}")
            return

        resume = status != "not_started"
        block_index += 1
        queue.log(
            f"group={group} run start: {stem} {observable} block={block_index} "
            f"resume={resume} before={summary}"
        )
        rc = queue.run_block(group, stem, observable, obs_file, block_index, resume=resume)
        after_status, after_summary = queue.task_status(stem, observable)
        records.append(
            {
                "stem": stem,
                "observable": observable,
                "status": after_status if rc == 0 else "failed",
                "returncode": rc,
                "param": queue.PARAM,
                "label": queue.run_label(stem, observable),
                "state_path": str(queue.state_path(stem, observable)),
                "obs_file": str(obs_file),
                **after_summary,
            }
        )
        queue.write_summary(group, records)
        queue.log(f"group={group} run done: {stem} {observable} rc={rc} after={after_summary}")
        if rc != 0:
            return


queue.summary_path = _summary_path
queue.build_command = _cpu_build_command
queue.run_task = _run_task


if __name__ == "__main__":
    if "--group" not in sys.argv:
        sys.argv.extend(["--group", "b"])
    queue.main()
