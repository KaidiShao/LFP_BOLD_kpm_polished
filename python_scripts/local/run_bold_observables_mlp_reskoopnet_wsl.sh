#!/usr/bin/env bash
set -euo pipefail

EXPERIMENT_NAME="${EXPERIMENT_NAME:-bold_wsl_20260423}"
LOCAL_OUTPUT_ROOT_BASE="${LOCAL_OUTPUT_ROOT_BASE:-/mnt/e/autodl_results_local/bold_wsl}"
OBSERVABLE_MODES="${OBSERVABLE_MODES:-eleHP HP roi_mean slow_band_power svd HP_svd100 global_svd100}"
RESIDUAL_FORMS="${RESIDUAL_FORMS:-projected_kv projected_vlambda}"
START_WITH_DATASET="${START_WITH_DATASET:-e10gb1}"
ONLY_START_DATASET="${ONLY_START_DATASET:-0}"
RESUME="${RESUME:-0}"
FORCE_RERUN="${FORCE_RERUN:-0}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-0}"
EXPORT_PSI="${EXPORT_PSI:-0}"
SELECTED_DEVICE="${SELECTED_DEVICE:-gpu}"
EPOCHS="${EPOCHS:-30}"
BATCH_SIZE="${BATCH_SIZE:-2000}"
CHUNK_SIZE="${CHUNK_SIZE:-5000}"
N_PSI_TRAIN="${N_PSI_TRAIN:-100}"
TRAIN_RATIO="${TRAIN_RATIO:-0.7}"
REG="${REG:-0.1}"
LR="${LR:-1e-4}"
LAYER_SIZES="${LAYER_SIZES:-100 100 100}"
PYTHON_EXE="${PYTHON_EXE:-python3}"
DRY_RUN="${DRY_RUN:-0}"
MAX_PARALLEL="${MAX_PARALLEL:-1}"
BLAS_THREADS="${BLAS_THREADS:-2}"
TF_INTRAOP_THREADS="${TF_INTRAOP_THREADS:-2}"
TF_INTEROP_THREADS="${TF_INTEROP_THREADS:-1}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$BLAS_THREADS}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-$BLAS_THREADS}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-$BLAS_THREADS}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-$BLAS_THREADS}"
export TF_NUM_INTRAOP_THREADS="${TF_NUM_INTRAOP_THREADS:-$TF_INTRAOP_THREADS}"
export TF_NUM_INTEROP_THREADS="${TF_NUM_INTEROP_THREADS:-$TF_INTEROP_THREADS}"

REPO_ROOT="/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished"
RUNNER="${REPO_ROOT}/python_scripts/autodl/run_autodl_reskoopnet_mlp.py"
SOLVER_DIR="${REPO_ROOT}/python_scripts/autodl"
BOLD_ROOT="/mnt/e/DataPons_processed/bold_observables"

dataset_stems=(e10gb1 e10gh1 e10fV1 f12m01)
dataset_ids=("E10.gb1" "E10.gH1" "E10.fV1" "F12.m01")

observable_file() {
  local dataset_id="$1"
  local mode="$2"
  printf '%s/%s/%s_bold_observables_%s.mat' "$BOLD_ROOT" "$dataset_id" "$dataset_id" "$mode"
}

run_label() {
  local experiment_name_for_dataset="$1"
  local residual_form="$2"
  local observable_mode="$3"
  printf 'mlp_obs_%s_%s_%s' "$experiment_name_for_dataset" "$residual_form" "$observable_mode"
}

is_complete() {
  local output_parent="$1"
  local label="$2"
  local output_dir="${output_parent}/${label}"
  compgen -G "${output_dir}/*_summary.mat" >/dev/null || return 1
  compgen -G "${output_dir}/*_outputs_*.mat" >/dev/null || return 1
  return 0
}

if [[ ! -f "$RUNNER" ]]; then
  echo "Local runner not found: $RUNNER" >&2
  exit 1
fi

if [[ "$DRY_RUN" != "1" ]]; then
  "$PYTHON_EXE" - <<'PY'
import importlib.util
import sys

required = ["h5py", "numpy", "scipy", "tensorflow", "koopmanlib"]
missing = [name for name in required if importlib.util.find_spec(name) is None]
if missing:
    print("Missing required Python modules in the selected WSL Python:", file=sys.stderr)
    print("  " + ", ".join(missing), file=sys.stderr)
    print("", file=sys.stderr)
    print("Use a Python/conda environment that already has the ResKoopNet dependencies,", file=sys.stderr)
    print("or create one in WSL, for example:", file=sys.stderr)
    print("  conda create -n reskoopnet python=3.10 -y", file=sys.stderr)
    print("  conda activate reskoopnet", file=sys.stderr)
    print("  pip install numpy scipy h5py matplotlib tqdm memory_profiler tensorflow koopmanlib", file=sys.stderr)
    print("", file=sys.stderr)
    print("Then rerun this launcher with:", file=sys.stderr)
    print("  -PythonExe /home/<user>/miniconda3/envs/reskoopnet/bin/python", file=sys.stderr)
    raise SystemExit(3)
PY
fi

start_index=-1
for i in "${!dataset_stems[@]}"; do
  if [[ "${dataset_stems[$i]}" == "$START_WITH_DATASET" ]]; then
    start_index="$i"
    break
  fi
done
if [[ "$start_index" -lt 0 ]]; then
  echo "Unknown START_WITH_DATASET: $START_WITH_DATASET" >&2
  exit 1
fi

ordered_indices=()
if [[ "$ONLY_START_DATASET" == "1" ]]; then
  ordered_indices=("$start_index")
else
  for ((offset = 0; offset < ${#dataset_stems[@]}; offset++)); do
    ordered_indices+=($(((start_index + offset) % ${#dataset_stems[@]})))
  done
fi

read -r -a mode_array <<< "$OBSERVABLE_MODES"
read -r -a residual_array <<< "$RESIDUAL_FORMS"
read -r -a layer_array <<< "$LAYER_SIZES"

if ! [[ "$MAX_PARALLEL" =~ ^[0-9]+$ ]] || [[ "$MAX_PARALLEL" -lt 1 ]]; then
  echo "MAX_PARALLEL must be a positive integer, got: $MAX_PARALLEL" >&2
  exit 1
fi

active_pids=()
active_labels=()
failed_count=0

wait_for_available_slot() {
  while [[ "$(jobs -rp | wc -l)" -ge "$MAX_PARALLEL" ]]; do
    sleep 5
  done
}

wait_for_all_parallel_jobs() {
  local pid
  local label
  local i

  for i in "${!active_pids[@]}"; do
    pid="${active_pids[$i]}"
    label="${active_labels[$i]}"
    if wait "$pid"; then
      echo "Parallel run finished successfully: $label"
    else
      echo "Parallel run failed: $label" >&2
      failed_count=$((failed_count + 1))
    fi
  done
}

for idx in "${ordered_indices[@]}"; do
  dataset_id="${dataset_ids[$idx]}"
  for mode in "${mode_array[@]}"; do
    file="$(observable_file "$dataset_id" "$mode")"
    if [[ ! -f "$file" ]]; then
      echo "Missing BOLD observable file for ${dataset_id} / ${mode}: ${file}" >&2
      exit 1
    fi
  done
done

for idx in "${ordered_indices[@]}"; do
  dataset_stem="${dataset_stems[$idx]}"
  dataset_id="${dataset_ids[$idx]}"
  dataset_experiment_name="${EXPERIMENT_NAME}_${dataset_stem}"
  dataset_result_root="${LOCAL_OUTPUT_ROOT_BASE}/${dataset_stem}/mlp"
  output_parent="${dataset_result_root}/outputs"
  checkpoint_parent="${dataset_result_root}/checkpoints"
  log_parent="${dataset_result_root}/logs"
  console_log_parent="${dataset_result_root}/console_logs"

  if [[ "$DRY_RUN" != "1" ]]; then
    mkdir -p "$output_parent" "$checkpoint_parent" "$log_parent" "$console_log_parent"
  fi

  for mode in "${mode_array[@]}"; do
    for residual_form in "${residual_array[@]}"; do
      label="$(run_label "$dataset_experiment_name" "$residual_form" "$mode")"

      if [[ "$FORCE_RERUN" != "1" ]] && is_complete "$output_parent" "$label"; then
        echo
        echo "Skipping complete local run: $label"
        continue
      fi

      data_file="$(observable_file "$dataset_id" "$mode")"
      console_log="${console_log_parent}/${label}.log"

      cmd=(
        "$PYTHON_EXE" "$RUNNER"
        "--project-root" "$REPO_ROOT"
        "--solver-dir" "$SOLVER_DIR"
        "--data-root" "$BOLD_ROOT"
        "--output-parent" "$output_parent"
        "--checkpoint-parent" "$checkpoint_parent"
        "--log-parent" "$log_parent"
        "--run-name-base" "mlp_obs"
        "--experiment-name" "$dataset_experiment_name"
        "--selected-device" "$SELECTED_DEVICE"
        "--solver-name" "resdmd_batch"
        "--residual-form" "$residual_form"
        "--data-subdir" "$dataset_id"
        "--dataset-stem" "$dataset_stem"
        "--observable-mode" "$mode"
        "--data-filename" "$(basename "$data_file")"
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
        "--inner-epochs" "5"
        "--end-condition" "1e-9"
        "--chunk-size" "$CHUNK_SIZE"
        "--layer-sizes" "${layer_array[@]}"
      )

      if [[ "$RESUME" == "1" ]]; then
        cmd+=("--resume")
      else
        cmd+=("--fresh-checkpoints")
      fi
      if [[ "$EXPORT_PSI" == "1" ]]; then
        cmd+=("--export-psi")
      fi

      echo
      printf '=%.0s' {1..80}
      echo
      echo "Starting WSL BOLD MLP ResKoopNet: $label"
      echo "Data file: $data_file"
      echo "Console log: $console_log"
      printf '=%.0s' {1..80}
      echo

      if [[ "$DRY_RUN" == "1" ]]; then
        printf 'Dry run command:\n  '
        printf '%q ' "${cmd[@]}"
        echo
        continue
      fi

      if [[ "$MAX_PARALLEL" -gt 1 ]]; then
        wait_for_available_slot
        (
          set +e
          "${cmd[@]}" 2>&1 | tee "$console_log"
          exit_code="${PIPESTATUS[0]}"
          echo "$exit_code" > "${console_log}.exit"
          if [[ "$exit_code" -ne 0 ]]; then
            echo "WSL BOLD MLP ResKoopNet failed for ${label} with exit code ${exit_code}" >&2
          fi
          exit "$exit_code"
        ) &
        active_pids+=("$!")
        active_labels+=("$label")
      else
        set +e
        "${cmd[@]}" 2>&1 | tee "$console_log"
        exit_code="${PIPESTATUS[0]}"
        set -e

        if [[ "$exit_code" -ne 0 ]]; then
          echo "WSL BOLD MLP ResKoopNet failed for ${label} with exit code ${exit_code}" >&2
          if [[ "$CONTINUE_ON_ERROR" == "1" ]]; then
            continue
          fi
          exit "$exit_code"
        fi
      fi
    done
  done
done

if [[ "$MAX_PARALLEL" -gt 1 && "$DRY_RUN" != "1" ]]; then
  wait_for_all_parallel_jobs
  if [[ "$failed_count" -gt 0 ]]; then
    echo "${failed_count} parallel run(s) failed." >&2
    if [[ "$CONTINUE_ON_ERROR" != "1" ]]; then
      exit 1
    fi
  fi
fi

if [[ "$DRY_RUN" == "1" ]]; then
  echo
  echo "Dry run finished successfully. No WSL jobs were launched."
else
  echo
  echo "All requested WSL BOLD MLP ResKoopNet runs finished successfully."
fi
