#!/usr/bin/env bash
set -euo pipefail

# LEGACY one-off local launcher for archived E10fV1 four-condition runs.
# Prefer python_scripts/local/run_bold_observables_mlp_reskoopnet_wsl.sh.

EXPERIMENT_NAME="${EXPERIMENT_NAME:-e10fV1_local4cond_run01}"
EPOCHS="${EPOCHS:-60}"
RESUME="${RESUME:-0}"
FORCE_RERUN="${FORCE_RERUN:-0}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-0}"
SELECTED_DEVICE="${SELECTED_DEVICE:-gpu}"
PYTHON_EXE="${PYTHON_EXE:-python3}"
RESULTS_ROOT="${RESULTS_ROOT:-/mnt/e/autodl_results}"
DATA_ROOT="${DATA_ROOT:-/mnt/e/DataPons_processed}"
BATCH_SIZE="${BATCH_SIZE:-2000}"
CHUNK_SIZE="${CHUNK_SIZE:-5000}"
N_PSI_TRAIN="${N_PSI_TRAIN:-100}"
TRAIN_RATIO="${TRAIN_RATIO:-0.7}"
REG="${REG:-0.1}"
LR="${LR:-1e-4}"
INNER_EPOCHS="${INNER_EPOCHS:-5}"
ONLY_OBSERVABLE_MODES="${ONLY_OBSERVABLE_MODES:-}"
ONLY_RESIDUAL_FORMS="${ONLY_RESIDUAL_FORMS:-}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
RUNNER="${REPO_ROOT}/python_scripts/autodl/run_autodl_reskoopnet_mlp.py"
SOLVER_DIR="${SCRIPT_DIR}"

DATASET_STEM="e10fV1"
DATA_SUBDIR="${DATASET_STEM}/pipeline1_reskoopnet_dictionary"
OUTPUT_BASE="${RESULTS_ROOT}/${DATASET_STEM}/mlp"
OUTPUT_PARENT="${OUTPUT_BASE}/outputs"
CHECKPOINT_PARENT="${OUTPUT_BASE}/checkpoints"
LOG_PARENT="${OUTPUT_BASE}/logs"
CONSOLE_LOG_PARENT="${OUTPUT_BASE}/console_logs"

OBSERVABLE_MODES=("abs" "complex_split")
RESIDUAL_FORMS=("projected_kv" "projected_vlambda")

declare -A DATA_FILENAMES=(
  ["abs"]="e10fV1_low50_high250_g2_abs_single.mat"
  ["complex_split"]="e10fV1_low50_high250_g2_complex_split_single.mat"
)

split_csv_into_array() {
  local csv="$1"
  local -n output_ref="$2"
  output_ref=()
  [[ -z "$csv" ]] && return 0
  IFS=',' read -r -a output_ref <<<"$csv"
}

array_contains() {
  local needle="$1"
  shift
  local candidate
  for candidate in "$@"; do
    [[ "$candidate" == "$needle" ]] && return 0
  done
  return 1
}

FILTER_OBSERVABLE_MODES=()
FILTER_RESIDUAL_FORMS=()
split_csv_into_array "$ONLY_OBSERVABLE_MODES" FILTER_OBSERVABLE_MODES
split_csv_into_array "$ONLY_RESIDUAL_FORMS" FILTER_RESIDUAL_FORMS

run_label() {
  local experiment_name="$1"
  local residual_form="$2"
  local observable_mode="$3"
  printf 'mlp_obs_%s_%s_%s' "$experiment_name" "$residual_form" "$observable_mode"
}

is_complete() {
  local output_parent="$1"
  local label="$2"
  local output_dir="${output_parent}/${label}"
  compgen -G "${output_dir}/*_summary.mat" >/dev/null || return 1
  compgen -G "${output_dir}/*_outputs_*.mat" >/dev/null || return 1
  return 0
}

require_file() {
  local file_path="$1"
  if [[ ! -f "$file_path" ]]; then
    echo "Required file not found: $file_path" >&2
    exit 1
  fi
}

require_file "$RUNNER"
for mode in "${OBSERVABLE_MODES[@]}"; do
  require_file "${DATA_ROOT}/${DATA_SUBDIR}/${DATA_FILENAMES[$mode]}"
done

if ! "$PYTHON_EXE" - <<'PY'
import importlib.util
import sys

required = ["h5py", "numpy", "scipy", "tensorflow", "koopmanlib"]
missing = [name for name in required if importlib.util.find_spec(name) is None]
if missing:
    print("Missing required Python modules in the selected WSL Python:", file=sys.stderr)
    print("  " + ", ".join(missing), file=sys.stderr)
    raise SystemExit(3)
PY
then
  echo >&2
  echo "Use a WSL Python environment that already has the ResKoopNet dependencies." >&2
  echo "Recommended install flow for this repo:" >&2
  echo "  1. Install Anaconda/Miniconda inside WSL" >&2
  echo "  2. cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/python_scripts/settings" >&2
  echo "  3. conda env create -f reskoopnet.yml -n reskoopnet" >&2
  echo "  4. conda activate reskoopnet" >&2
  echo "  5. sudo apt install -y build-essential" >&2
  echo "  6. git clone https://github.com/MLDS-NUS/KoopmanDL.git ~/KoopmanDL" >&2
  echo "  7. cd ~/KoopmanDL && pip install -r requirements.txt && pip install -e ." >&2
  echo >&2
  echo "Then rerun this launcher with something like:" >&2
  echo "  PYTHON_EXE=~/anaconda3/envs/reskoopnet/bin/python bash run_E10fV1_4conditions_60epoch_local.sh" >&2
  exit 3
fi

mkdir -p "$OUTPUT_PARENT" "$CHECKPOINT_PARENT" "$LOG_PARENT" "$CONSOLE_LOG_PARENT"

echo
printf '=%.0s' {1..80}
echo
echo "LEGACY launcher warning: run_E10fV1_4conditions_60epoch_local.sh is kept only for archived four-condition reproduction."
echo "E10fV1 local 4-condition BLP ResKoopNet launcher"
echo "Experiment name: ${EXPERIMENT_NAME}"
echo "Epochs: ${EPOCHS}"
echo "Selected device: ${SELECTED_DEVICE}"
echo "Observable modes: ${ONLY_OBSERVABLE_MODES:-all}"
echo "Residual forms: ${ONLY_RESIDUAL_FORMS:-all}"
echo "Results root: ${OUTPUT_BASE}"
echo "Resume: ${RESUME}"
echo "Force rerun: ${FORCE_RERUN}"
printf '=%.0s' {1..80}
echo

for observable_mode in "${OBSERVABLE_MODES[@]}"; do
  if ((${#FILTER_OBSERVABLE_MODES[@]} > 0)) && ! array_contains "$observable_mode" "${FILTER_OBSERVABLE_MODES[@]}"; then
    continue
  fi
  for residual_form in "${RESIDUAL_FORMS[@]}"; do
    if ((${#FILTER_RESIDUAL_FORMS[@]} > 0)) && ! array_contains "$residual_form" "${FILTER_RESIDUAL_FORMS[@]}"; then
      continue
    fi
    label="$(run_label "${EXPERIMENT_NAME}" "${residual_form}" "${observable_mode}")"
    data_file="${DATA_FILENAMES[$observable_mode]}"
    console_log="${CONSOLE_LOG_PARENT}/${label}.log"

    if [[ "${FORCE_RERUN}" != "1" ]] && is_complete "$OUTPUT_PARENT" "$label"; then
      echo
      echo "Skipping complete run: ${label}"
      continue
    fi

    cmd=(
      "$PYTHON_EXE" "$RUNNER"
      "--project-root" "$SOLVER_DIR"
      "--solver-dir" "$SOLVER_DIR"
      "--data-root" "$DATA_ROOT"
      "--output-parent" "$OUTPUT_PARENT"
      "--checkpoint-parent" "$CHECKPOINT_PARENT"
      "--log-parent" "$LOG_PARENT"
      "--run-name-base" "mlp_obs"
      "--experiment-name" "$EXPERIMENT_NAME"
      "--selected-device" "$SELECTED_DEVICE"
      "--solver-name" "resdmd_batch"
      "--residual-form" "$residual_form"
      "--data-subdir" "$DATA_SUBDIR"
      "--dataset-stem" "$DATASET_STEM"
      "--observable-mode" "$observable_mode"
      "--data-filename" "$data_file"
      "--file-type" ".mat"
      "--field-name" "obs"
      "--n-psi-train" "$N_PSI_TRAIN"
      "--train-ratio" "$TRAIN_RATIO"
      "--reg" "$REG"
      "--rounds" "1"
      "--epochs" "$EPOCHS"
      "--batch-size" "$BATCH_SIZE"
      "--lr" "$LR"
      "--log-interval" "1"
      "--lr-decay-factor" "0.8"
      "--inner-epochs" "$INNER_EPOCHS"
      "--end-condition" "1e-9"
      "--chunk-size" "$CHUNK_SIZE"
      "--layer-sizes" "100" "100" "100"
    )

    if [[ "${RESUME}" == "1" ]]; then
      cmd+=("--resume")
    else
      cmd+=("--fresh-checkpoints")
    fi

    echo
    printf '%s\n' "--------------------------------------------------------------------------------"
    echo "Starting local BLP ResKoopNet: ${label}"
    echo "Observable mode: ${observable_mode}"
    echo "Residual form: ${residual_form}"
    echo "Input file: ${DATA_ROOT}/${DATA_SUBDIR}/${data_file}"
    echo "Console log: ${console_log}"
    printf '%s\n' "--------------------------------------------------------------------------------"

    set +e
    "${cmd[@]}" 2>&1 | tee "$console_log"
    exit_code="${PIPESTATUS[0]}"
    set -e

    if [[ "$exit_code" -ne 0 ]]; then
      echo "Run failed: ${label} (exit code ${exit_code})" >&2
      if [[ "${CONTINUE_ON_ERROR}" == "1" ]]; then
        continue
      fi
      exit "$exit_code"
    fi
  done
done

echo
echo "All requested E10fV1 local BLP ResKoopNet runs finished successfully."
