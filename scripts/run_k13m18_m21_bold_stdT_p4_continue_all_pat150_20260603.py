#!/usr/bin/env python3
"""Continue all K13m18-K13m21 BOLD P4 stdT runs with patience 150."""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import scripts.run_k13m18_m21_bold_stdT_p4_queue as previous_queue


queue = previous_queue.queue

queue.BLOCK_EPOCHS = 20
queue.STALE_EPOCHS = 150
queue.MAX_TOTAL_EPOCHS = 1000
queue.OBSERVABLES = ("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean")
queue.MAIN_LOG = queue.RUN_LOG_DIR / "bold_stdT_k13m18_m21_p4_continue_all_pat150_20260603.log"


def _summary_path(group: str) -> Path:
    return queue.RUN_LOG_DIR / "bold_stdT_k13m18_m21_p4_continue_all_pat150_20260603_summary.json"


queue.summary_path = _summary_path


if __name__ == "__main__":
    queue.main()
