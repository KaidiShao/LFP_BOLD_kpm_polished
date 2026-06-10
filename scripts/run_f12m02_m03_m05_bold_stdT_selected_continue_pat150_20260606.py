#!/usr/bin/env python3
"""Continue selected F12 BOLD stdT P4 runs with patience 150 on CPU."""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scripts.run_f12m02_m03_m05_bold_stdT_p4_cpu_queue as base


queue = base.queue

SELECTED_TASKS = (
    # The first continuation attempt reset f12m02/global_svd100 from best.
    # Continue from that new final state now, otherwise repeated best-resume
    # blocks can replay the same short segment if no new best is found.
    ("f12m02", "global_svd100", "final"),
    ("f12m02", "HP_svd100", "final"),
    ("f12m02", "roi_mean", "final"),
    ("f12m03", "HP_svd100", "final"),
    ("f12m03", "roi_mean", "final"),
    ("f12m05", "gsvd100_ds", "final"),
    ("f12m05", "roi_mean", "final"),
)

queue.DATASETS = [
    {"stem": "f12m02", "group": "b"},
    {"stem": "f12m03", "group": "b"},
    {"stem": "f12m05", "group": "b"},
]
queue.OBSERVABLES = tuple(sorted({observable for _, observable, _ in SELECTED_TASKS}))
queue.STALE_EPOCHS = 150
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_cpu_f12m02_m03_m05_stdT_selected_continue_pat150_20260606.log"

_original_summary_path = queue.summary_path
_original_build_command = queue.build_command


def _summary_path(group: str) -> Path:
    return queue.RUN_LOG_DIR / "bold_cpu_f12m02_m03_m05_stdT_selected_continue_pat150_20260606_summary.json"


def _task_resume_mode(stem: str, observable: str) -> str:
    for task_stem, task_observable, resume_mode in SELECTED_TASKS:
        if stem == task_stem and observable == task_observable:
            return resume_mode
    return "final"


def _build_command(stem: str, observable: str, obs_file: Path, epochs: int, resume: bool) -> list[str]:
    command = _original_build_command(stem, observable, obs_file, epochs, resume=resume)
    if not resume:
        return command
    resume_mode = _task_resume_mode(stem, observable)
    for idx, token in enumerate(command):
        if token == "--resume-mode" and idx + 1 < len(command):
            command[idx + 1] = resume_mode
            break
    return command


def _make_tasks(group: str) -> list[dict]:
    tasks: list[dict] = []
    if group != "b":
        return tasks
    for stem, observable, resume_mode in SELECTED_TASKS:
        obs_file = queue.find_observable_file(stem, observable)
        tasks.append(
            {
                "stem": stem,
                "observable": observable,
                "resume_mode": resume_mode,
                "obs_file": str(obs_file) if obs_file else None,
            }
        )
    return tasks


queue.summary_path = _summary_path
queue.build_command = _build_command
queue.make_tasks = _make_tasks


if __name__ == "__main__":
    if "--group" not in sys.argv:
        sys.argv.extend(["--group", "b"])
    queue.main()
