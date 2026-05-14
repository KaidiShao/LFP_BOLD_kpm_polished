#!/usr/bin/env bash
set -euo pipefail

# LEGACY local complex-split batch wrapper.
# Prefer python_scripts/local/run_bold_observables_mlp_reskoopnet_wsl.sh.

EXPERIMENT_NAME="${1:-local_complexsplit_60_20260418}"
RESUME_FLAG="${2:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_EXE="${PYTHON_EXE:-python3}"

echo "LEGACY launcher warning: run_E10gb1_E10fV1_E10gH1_complex_split_60epoch.sh is kept only for archived local batch reproduction." >&2

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
