#!/usr/bin/env python3
"""Continue only f12m02 BOLD stdT roi_mean after the first pat150 pass."""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scripts.run_f12m02_m03_m05_bold_stdT_selected_continue_pat150_20260606 as base


queue = base.queue

queue.DATASETS = [{"stem": "f12m02", "group": "b"}]
queue.OBSERVABLES = ("roi_mean",)
queue.STALE_EPOCHS = 150
queue.MAX_TOTAL_EPOCHS = 1600
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_cpu_f12m02_roi_mean_stdT_continue_extra_pat150_20260606.log"


def _summary_path(group: str) -> Path:
    return queue.RUN_LOG_DIR / "bold_cpu_f12m02_roi_mean_stdT_continue_extra_pat150_20260606_summary.json"


def _make_tasks(group: str) -> list[dict]:
    if group != "b":
        return []
    obs_file = queue.find_observable_file("f12m02", "roi_mean")
    return [
        {
            "stem": "f12m02",
            "observable": "roi_mean",
            "resume_mode": "final",
            "obs_file": str(obs_file) if obs_file else None,
        }
    ]


queue.summary_path = _summary_path
queue.make_tasks = _make_tasks


if __name__ == "__main__":
    if "--group" not in sys.argv:
        sys.argv.extend(["--group", "b"])
    queue.main()
