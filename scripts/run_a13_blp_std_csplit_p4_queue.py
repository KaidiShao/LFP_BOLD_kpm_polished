#!/usr/bin/env python3
"""Run A13 BLP P4 standardized complex-split queue.

This is an A13-only wrapper around the accepted standardized complex-split
BLP P4 branch. It keeps the run labels/parameters used by downstream current
standardized-csplit analyses and only changes the dataset list plus log names.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scripts.run_k13m18_m21_blp_std_csplit_p4_queue as queue


queue.DATASETS = ("a13pf1", "a13pt1", "a13qk1", "a13qo1", "a13s11")
queue.RUN_SUFFIX = "20260522_pat40_allblp"
queue.MAIN_LOG = queue.RUN_LOG_DIR / "blp_std_csplit_a13_p4_queue_20260611.log"
queue.LOCK_PATH = queue.RUN_LOG_DIR / "blp_std_csplit_a13_p4_queue_20260611.lock"


def _summary_path() -> Path:
    return queue.RUN_LOG_DIR / "blp_std_csplit_a13_p4_queue_20260611_summary.json"


queue.summary_path = _summary_path


if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.extend(["--datasets", *queue.DATASETS])
    queue.main()
