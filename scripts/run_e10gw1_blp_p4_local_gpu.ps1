param(
    [string]$WslDistro = "Ubuntu-22.04",
    [string]$PythonExe = "/home/kdshao/anaconda3/bin/python",
    [string]$RepoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished",
    [string]$DataRoot = "/mnt/e/DataPons_processed/e10gw1",
    [string]$OutputRoot = "/mnt/e/autodl_results_new/e10gw1/mlp",
    [string]$ExperimentName = "blp_vlambda_mainline_torchlike_20260517_e10gw1_seed1234",
    [string[]]$Modes = @("abs", "complex_split"),
    [int]$Epochs = 50,
    [int]$BatchSize = 2000,
    [int]$ChunkSize = 5000,
    [string]$CudaVisibleDevices = "0",
    [switch]$Resume,
    [switch]$DryRun,
    [switch]$NoRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function ConvertTo-WslPath {
    param([string]$Path)

    if ($Path -match '^([A-Za-z]):\\(.*)$') {
        $drive = $Matches[1].ToLowerInvariant()
        $rest = $Matches[2] -replace '\\', '/'
        return "/mnt/$drive/$rest"
    }

    return ($Path -replace '\\', '/')
}

function Quote-Bash {
    param([string]$Value)
    return "'" + ($Value -replace "'", "'`"`"'`"`"'") + "'"
}

$repoWsl = ConvertTo-WslPath $RepoRoot
$tmpDir = Join-Path $RepoRoot "tmp\local_gpu_e10gw1_p4"
New-Item -ItemType Directory -Force -Path $tmpDir | Out-Null

$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$bashPathWin = Join-Path $tmpDir "run_e10gw1_blp_p4_local_gpu_$stamp.sh"
$bashPathWsl = ConvertTo-WslPath $bashPathWin
$modeExpr = ($Modes | ForEach-Object { Quote-Bash $_ }) -join " "
$resumeFlag = if ($Resume) { "1" } else { "0" }
$dryRunFlag = if ($DryRun) { "1" } else { "0" }

$template = @'
#!/usr/bin/env bash
set -euo pipefail

export CUDA_VISIBLE_DEVICES=__CUDA_VISIBLE_DEVICES__
export TF_FORCE_GPU_ALLOW_GROWTH=true
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-2}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-2}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-2}"
export TF_NUM_INTRAOP_THREADS="${TF_NUM_INTRAOP_THREADS:-2}"
export TF_NUM_INTEROP_THREADS="${TF_NUM_INTEROP_THREADS:-1}"

PY=__PYTHON_EXE__
REPO=__REPO_WSL__
DATA_ROOT=__DATA_ROOT__
OUT_ROOT=__OUTPUT_ROOT__
EXP=__EXPERIMENT_NAME__
EPOCHS=__EPOCHS__
BATCH_SIZE=__BATCH_SIZE__
CHUNK_SIZE=__CHUNK_SIZE__
RESUME=__RESUME_FLAG__
DRY_RUN=__DRY_RUN_FLAG__
modes=(__MODES__)

mkdir -p "$OUT_ROOT/outputs" "$OUT_ROOT/checkpoints" "$OUT_ROOT/logs" "$OUT_ROOT/console_logs"

if [[ "$DRY_RUN" != "1" ]]; then
  "$PY" - <<'PYGPU'
import tensorflow as tf
gpus = tf.config.list_physical_devices("GPU")
print("GPUS:", gpus)
if not gpus:
    raise SystemExit("No TensorFlow GPU visible")
PYGPU
fi

for MODE in "${modes[@]}"; do
  LOG="$OUT_ROOT/console_logs/local_gpu_fresh_${MODE}_$(date +%Y%m%d_%H%M%S).log"
  echo
  echo "================================================================================"
  echo "E10.gW1 BLP P4 local GPU: $MODE"
  echo "Experiment: $EXP"
  echo "Log: $LOG"
  echo "================================================================================"

  cmd=(
    "$PY" "$REPO/python_scripts/autodl/run_autodl_reskoopnet_mlp.py"
    --project-root "$REPO"
    --solver-dir "$REPO/python_scripts/autodl"
    --data-root "$DATA_ROOT"
    --data-subdir pipeline1_reskoopnet_dictionary
    --output-parent "$OUT_ROOT/outputs"
    --checkpoint-parent "$OUT_ROOT/checkpoints"
    --log-parent "$OUT_ROOT/logs"
    --run-name-base mlp_obs
    --experiment-name "$EXP"
    --selected-device gpu
    --require-gpu
    --solver-name resdmd_batch
    --residual-form projected_vlambda
    --dataset-stem e10gw1
    --observable-mode "$MODE"
    --file-type .mat
    --field-name obs
    --n-psi-train 100
    --train-ratio 0.7
    --reg 0.1
    --rounds 1
    --epochs "$EPOCHS"
    --batch-size "$BATCH_SIZE"
    --lr 1e-4
    --log-interval 1
    --lr-decay-factor 0.8
    --inner-epochs 5
    --end-condition 1e-9
    --chunk-size "$CHUNK_SIZE"
    --layer-sizes 100 100 100
    --train-shuffle
    --spectral-sync-mode pre_only
    --seed 1234
    --export-mode full
    --skip-diagnostic-pdf
  )

  if [[ "$RESUME" == "1" ]]; then
    cmd+=(--resume --resume-mode best)
  else
    cmd+=(--fresh-checkpoints)
  fi

  if [[ "$DRY_RUN" == "1" ]]; then
    printf 'Dry run command:\n  '
    printf '%q ' "${cmd[@]}"
    echo
  else
    "${cmd[@]}" 2>&1 | tee "$LOG"
  fi
done

echo
echo "E10.gW1 local GPU P4 runner finished."
echo "Output root: $OUT_ROOT/outputs"
'@

$bashContent = $template.
    Replace("__CUDA_VISIBLE_DEVICES__", (Quote-Bash $CudaVisibleDevices)).
    Replace("__PYTHON_EXE__", (Quote-Bash $PythonExe)).
    Replace("__REPO_WSL__", (Quote-Bash $repoWsl)).
    Replace("__DATA_ROOT__", (Quote-Bash $DataRoot)).
    Replace("__OUTPUT_ROOT__", (Quote-Bash $OutputRoot)).
    Replace("__EXPERIMENT_NAME__", (Quote-Bash $ExperimentName)).
    Replace("__EPOCHS__", [string]$Epochs).
    Replace("__BATCH_SIZE__", [string]$BatchSize).
    Replace("__CHUNK_SIZE__", [string]$ChunkSize).
    Replace("__RESUME_FLAG__", $resumeFlag).
    Replace("__DRY_RUN_FLAG__", $dryRunFlag).
    Replace("__MODES__", $modeExpr)

Set-Content -LiteralPath $bashPathWin -Value $bashContent -Encoding ASCII

Write-Host "Generated WSL runner:"
Write-Host "  $bashPathWin"
Write-Host "WSL path:"
Write-Host "  $bashPathWsl"

if ($NoRun) {
    Write-Host "NoRun requested; not launching WSL."
    exit 0
}

& wsl.exe -d $WslDistro bash $bashPathWsl
$code = $LASTEXITCODE
if ($code -ne 0) {
    throw "WSL local GPU runner failed with exit code $code"
}
