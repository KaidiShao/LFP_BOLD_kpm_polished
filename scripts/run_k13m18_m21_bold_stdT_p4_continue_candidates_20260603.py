#!/usr/bin/env python3
"""Continue selected K13m18-K13m21 BOLD P4 stdT runs.

The 2026-06-02 queue stopped all runs with a short patience. A few curves were
still trending down near the end, so this wrapper keeps the accepted labels and
parameters, but restricts continuation to those candidate dataset/observable
pairs and uses a wider stale window.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scripts.run_k13m18_m21_bold_stdT_p4_queue as previous_queue


queue = previous_queue.queue

CANDIDATES = (
    ("k13m21", "HP_svd100"),
    ("k13m20", "HP_svd100"),
    ("k13m20", "roi_mean"),
    ("k13m18", "HP_svd100"),
    ("k13m20", "gsvd100_ds"),
)

queue.BLOCK_EPOCHS = 20
queue.STALE_EPOCHS = 100
queue.MAX_TOTAL_EPOCHS = 500
queue.OBSERVABLES = tuple(dict.fromkeys(observable for _, observable in CANDIDATES))
queue.MAIN_LOG = (
    queue.RUN_LOG_DIR / "bold_stdT_k13m18_m21_p4_continue_candidates_20260603.log"
)


def _summary_path(group: str) -> Path:
    return (
        queue.RUN_LOG_DIR
        / "bold_stdT_k13m18_m21_p4_continue_candidates_20260603_summary.json"
    )


def _task_status(stem: str, observable: str) -> tuple[str, dict]:
    summary = queue.summarize_state(queue.read_state(stem, observable))
    epochs = int(summary["epochs"] or 0)
    stale = int(summary["stale_epochs"] or 0)
    if epochs <= 0:
        return "not_started", summary
    if epochs >= queue.MAX_TOTAL_EPOCHS:
        return "max_epoch_reached", summary
    if stale >= queue.STALE_EPOCHS:
        return "stale_done", summary
    return "running", summary


def _make_tasks(group: str) -> list[dict]:
    if group != "b":
        return []
    tasks: list[dict] = []
    for stem, observable in CANDIDATES:
        obs_file = queue.find_observable_file(stem, observable)
        tasks.append(
            {
                "stem": stem,
                "observable": observable,
                "obs_file": str(obs_file) if obs_file else None,
            }
        )
    return tasks


queue.summary_path = _summary_path
queue.task_status = _task_status
queue.make_tasks = _make_tasks


if __name__ == "__main__":
    queue.main()
