#!/usr/bin/env python3
"""Run all-current BOLD P4 standardized roi_mean training on CPU, npsi=200.

This creates a separate run family from the existing npsi=100 stdT pat40
results. Only --n-psi-train changes to 200; hidden layers remain 100 100 100.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tmp.run_bold_stdT_roi_global_pat40_queue_20260522 as wrapped


queue = wrapped.queue

queue.DATASETS = [
    {"stem": "a13pf1", "group": "b"},
    {"stem": "a13pt1", "group": "b"},
    {"stem": "a13qk1", "group": "b"},
    {"stem": "a13qo1", "group": "b"},
    {"stem": "a13s11", "group": "b"},
    {"stem": "e10aw1", "group": "b"},
    {"stem": "e10bv1", "group": "b"},
    {"stem": "e10fV1", "group": "b"},
    {"stem": "e10gb1", "group": "b"},
    {"stem": "e10gh1", "group": "b"},
    {"stem": "e10gw1", "group": "b"},
    {"stem": "f12m01", "group": "b"},
    {"stem": "f12m02", "group": "b"},
    {"stem": "f12m03", "group": "b"},
    {"stem": "f12m05", "group": "b"},
    {"stem": "k13m17", "group": "b"},
    {"stem": "k13m18", "group": "b"},
    {"stem": "k13m19", "group": "b"},
    {"stem": "k13m20", "group": "b"},
    {"stem": "k13m21", "group": "b"},
    {"stem": "k13m23", "group": "b"},
]
queue.OBSERVABLES = ("roi_mean",)
queue.PARAM = {
    "lr": "1e-4",
    "reg": "0.001",
    "batch_size": "2000",
    "inner_epochs": "2",
    "tag": "stdT_npsi200_l1e4_r1e3_b2000_i2_pat40",
}
queue.RUN_SUFFIX = "20260611_npsi200_roi_mean_allbold"
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_cpu_all_current_roi_mean_stdT_npsi200_p4_queue_20260611.log"
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
    return queue.RUN_LOG_DIR / "bold_cpu_all_current_roi_mean_stdT_npsi200_p4_queue_20260611_summary.json"


def _cpu_npsi200_build_command(
    stem: str, observable: str, obs_file: Path, epochs: int, resume: bool
) -> list[str]:
    command = _original_build_command(stem, observable, obs_file, epochs, resume=resume)
    next_replacement: str | None = None
    rewritten: list[str] = []
    skip_next = False
    for token in command:
        if skip_next:
            skip_next = False
            continue
        if next_replacement is not None:
            rewritten.append(next_replacement)
            next_replacement = None
            continue
        if token == "--selected-device":
            rewritten.extend(["--selected-device", "cpu"])
            skip_next = True
            continue
        if token == "--require-gpu":
            continue
        if token == "--n-psi-train":
            rewritten.append(token)
            next_replacement = "200"
            continue
        rewritten.append(token)
    return rewritten


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
queue.build_command = _cpu_npsi200_build_command
queue.run_task = _run_task


if __name__ == "__main__":
    if "--group" not in sys.argv:
        sys.argv.extend(["--group", "b"])
    queue.main()
