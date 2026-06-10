#!/usr/bin/env python3
"""Continue selected K13 BOLD P4 stdT runs with patience 250."""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scripts.run_k13m18_m21_bold_stdT_p4_queue as previous_queue


queue = previous_queue.queue

SELECTED = (
    ("k13m20", "HP_svd100"),
    ("k13m21", "global_svd100"),
    ("k13m21", "roi_mean"),
)

queue.BLOCK_EPOCHS = 20
queue.STALE_EPOCHS = 250
queue.MAX_TOTAL_EPOCHS = 1200
queue.OBSERVABLES = ("HP_svd100", "global_svd100", "roi_mean")
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_stdT_k13m20_m21_p4_continue_selected_pat250_20260603.log"


def _summary_path(group: str) -> Path:
    return queue.RUN_LOG_DIR / "bold_stdT_k13m20_m21_p4_continue_selected_pat250_20260603_summary.json"


def _make_tasks(group: str) -> list[dict]:
    tasks: list[dict] = []
    if group != "b":
        return tasks
    for stem, observable in SELECTED:
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
queue.make_tasks = _make_tasks


if __name__ == "__main__":
    queue.main()
