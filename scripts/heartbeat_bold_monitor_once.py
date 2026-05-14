from __future__ import annotations

import re
from pathlib import Path


REPORT = Path("tmp/monitor_reports/bold_global_svd_monitor_latest.md")
POSTED = Path("tmp/monitor_reports/bold_global_svd_last_posted.txt")


def main() -> int:
    if not REPORT.exists():
        print("DECISION=DONT_NOTIFY")
        print("MESSAGE=BOLD monitor report is not available yet.")
        return 0

    text = REPORT.read_text(encoding="utf-8")
    match = re.search(r"\*\*BOLD Monitor ([^*]+)\*\*", text)
    stamp = match.group(1).strip() if match else REPORT.stat().st_mtime_ns.__str__()
    previous = POSTED.read_text(encoding="utf-8").strip() if POSTED.exists() else ""

    if stamp == previous:
        print("DECISION=DONT_NOTIFY")
        print(f"MESSAGE=No new BOLD monitor report since {stamp}.")
        return 0

    POSTED.parent.mkdir(parents=True, exist_ok=True)
    POSTED.write_text(stamp, encoding="utf-8")
    print("DECISION=NOTIFY")
    print(f"MESSAGE=New BOLD monitor report {stamp}.")
    print("REPORT_START")
    print(text)
    print("REPORT_END")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
