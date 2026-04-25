#!/usr/bin/env bash
set -euo pipefail

EXPERIMENT_NAME="${1:-local_complexsplit_60_20260418}"
RESUME_FLAG="${2:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_EXE="${PYTHON_EXE:-python3}"

CMD=(
  "$PYTHON_EXE"
  "$SCRIPT_DIR/run_local_complex_split_batch.py"
  "--experiment-name" "$EXPERIMENT_NAME"
  "--epochs" "60"
)

if [[ "$RESUME_FLAG" == "--resume" ]]; then
  CMD+=("--resume")
else
  CMD+=("--fresh-checkpoints")
fi

exec "${CMD[@]}"
