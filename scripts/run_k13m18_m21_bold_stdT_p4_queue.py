#!/usr/bin/env python3
"""Run current standardized BOLD P4 for K13m18-K13m21.

This wrapper reuses the accepted 2026-05-22 BOLD stdT P4 queue parameters
and restricts the target set to the four new K13 datasets. It runs the current
four BOLD observables used downstream by P7/P8/P9/P10:

  global_svd100, gsvd100_ds, HP_svd100, roi_mean

The wrapper keeps its own summary/log names so it does not affect older queues.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tmp.run_bold_stdT_roi_global_pat40_queue_20260522 as wrapped


queue = wrapped.queue

queue.DATASETS = [
    {"stem": "k13m18", "group": "b"},
    {"stem": "k13m19", "group": "b"},
    {"stem": "k13m20", "group": "b"},
    {"stem": "k13m21", "group": "b"},
]
queue.OBSERVABLES = ("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean")
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_stdT_k13m18_m21_p4_queue_20260602.log"

DONE_STATUSES = {"stale_done", "delta_done", "max_epoch_reached"}


def _summary_path(group: str) -> Path:
    return queue.RUN_LOG_DIR / "bold_stdT_k13m18_m21_p4_queue_20260602_summary.json"


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
            queue.log(
                f"group={group} task done: {stem} {observable} "
                f"status={status} summary={summary}"
            )
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
queue.run_task = _run_task


if __name__ == "__main__":
    queue.main()
