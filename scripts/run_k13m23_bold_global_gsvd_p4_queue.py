#!/usr/bin/env python3
"""Run standardized BOLD P4 for K13m23 global_svd100 and gsvd100_ds.

This is the targeted continuation after
``build_k13m23_bold_svd_observables_python.py`` fills the missing P3 BOLD
observable files.
"""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import tmp.run_bold_stdT_roi_global_pat40_queue_20260522 as wrapped


queue = wrapped.queue
queue.DATASETS = [{"stem": "k13m23", "group": "b"}]
queue.OBSERVABLES = ("global_svd100", "gsvd100_ds")
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_stdT_k13m23_global_gsvd_pat40_20260601.log"


def _summary_path(group: str) -> Path:
    return queue.RUN_LOG_DIR / "bold_stdT_k13m23_global_gsvd_pat40_20260601_summary.json"


queue.summary_path = _summary_path


if __name__ == "__main__":
    queue.main()
