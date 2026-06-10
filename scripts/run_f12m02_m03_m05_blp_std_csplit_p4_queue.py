#!/usr/bin/env python3
"""Run BLP P4 standardized complex-split training for F12m02/F12m03/F12m05.

This is a thin F12 wrapper around the current standardized complex-split BLP
P4 queue. It keeps the accepted training parameters and only changes the
dataset list plus log/summary names.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scripts.run_k13m18_m21_blp_std_csplit_p4_queue as queue


queue.DATASETS = ("f12m02", "f12m03", "f12m05")
queue.RUN_SUFFIX = "20260522_pat40_allblp"
queue.MAIN_LOG = queue.RUN_LOG_DIR / "blp_std_csplit_f12m02_m03_m05_p4_queue_20260605.log"
queue.LOCK_PATH = queue.RUN_LOG_DIR / "blp_std_csplit_f12m02_m03_m05_p4_queue_20260605.lock"


def _summary_path() -> Path:
    return queue.RUN_LOG_DIR / "blp_std_csplit_f12m02_m03_m05_p4_queue_20260605_summary.json"


queue.summary_path = _summary_path


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.extend(["--datasets", *queue.DATASETS])
    queue.main()
