<# 
LEGACY one-off launcher for the archived E10fV1 four-condition BLP run.

This is not the canonical pipeline 4 entry.
Prefer:
  python_scripts/autodl/run_blp_observables_mlp_reskoopnet.ps1
  python_scripts/local/run_blp_observables_mlp_reskoopnet_local.ps1
#>

param(
    [string]$ExperimentName = "e10fV1_4cond_run01",
    [int]$Epochs = 60,
    [switch]$Resume,
    [switch]$OnlyAbs,
    [switch]$OnlyComplexSplit,
    [switch]$ContinueOnError,
    [string]$RemoteDataSubdirAbs = "e10fV1_abs",
    [string]$LocalDownloadRoot = "E:\autodl_results\e10fV1\mlp",
    [string]$WslPython = "python3"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

if ($OnlyAbs -and $OnlyComplexSplit) {
    throw "OnlyAbs and OnlyComplexSplit cannot both be set."
}

$repoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished"
$pythonExe = "python"
$processedRoot = "E:\DataPons_processed"

$remoteDatasetController = Join-Path $repoRoot "python_scripts\autodl\dataset_batch_controller_autodl_reskoopnet_mlp.py"
$localComplexSplitRunner = Join-Path $repoRoot "python_scripts\autodl_local_e10gb1\run_local_complex_split_batch.py"

$absObservableFile = Join-Path $processedRoot "e10fV1\pipeline1_reskoopnet_dictionary\e10fV1_low50_high250_g2_abs_single.mat"
$absObsInfoFile = Join-Path $processedRoot "e10fV1\pipeline1_reskoopnet_dictionary\e10fV1_low50_high250_g2_abs_single_obs_info.csv"
$complexObservableFile = Join-Path $processedRoot "e10fV1\pipeline1_reskoopnet_dictionary\e10fV1_low50_high250_g2_complex_split_single.mat"

function Assert-PathExists {
    param(
        [string]$PathValue,
        [string]$Label
    )

    if (-not (Test-Path -LiteralPath $PathValue)) {
        throw "$Label not found: $PathValue"
    }
}

function Convert-ToWslPath {
    param([string]$WindowsPath)

    $normalized = $WindowsPath -replace "\\", "/"
    if ($normalized -match "^([A-Za-z]):/(.*)$") {
        $drive = $matches[1].ToLowerInvariant()
        $rest = $matches[2]
        return "/mnt/$drive/$rest"
    }

    throw "Could not convert Windows path to WSL path: $WindowsPath"
}

Assert-PathExists -PathValue $remoteDatasetController -Label "Remote dataset controller"
Assert-PathExists -PathValue $localComplexSplitRunner -Label "Local complex_split runner"
Assert-PathExists -PathValue $absObservableFile -Label "ABS observable file"
Assert-PathExists -PathValue $absObsInfoFile -Label "ABS obs_info CSV"
Assert-PathExists -PathValue $complexObservableFile -Label "complex_split observable file"

if (-not (Test-Path -LiteralPath $LocalDownloadRoot)) {
    New-Item -ItemType Directory -Path $LocalDownloadRoot -Force | Out-Null
}

$commonRemoteArgs = @(
    $remoteDatasetController,
    "--ssh-host", "connect.westb.seetacloud.com",
    "--ssh-user", "root",
    "--ssh-port", "19241",
    "--ssh-key", "$HOME\.ssh\id_ed25519_autodl",
    "--dataset-stem", "e10fV1",
    "--experiment-name", $ExperimentName,
    "--observable-modes", "abs",
    "--residual-forms", "projected_kv", "projected_vlambda",
    "--file-type", ".mat",
    "--remote-python", "/root/miniconda3/envs/reskoopnet/bin/python",
    "--epochs", [string]$Epochs,
    "--remote-data-subdir", $RemoteDataSubdirAbs,
    "--local-download-root", $LocalDownloadRoot,
    "--local-data-file-abs", $absObservableFile,
    "--local-obs-info-file-abs", $absObsInfoFile,
    "--recover-completed-remote-runs",
    "--skip-completed-local-runs",
    "--delete-remote-run-after-download",
    "--delete-remote-input-after-observable"
)

if ($Resume) {
    $commonRemoteArgs += "--resume"
} else {
    $commonRemoteArgs += "--fresh-checkpoints"
}

if ($ContinueOnError) {
    $commonRemoteArgs += "--continue-on-error"
}

$wslRepoRoot = Convert-ToWslPath -WindowsPath $repoRoot
$wslResultsRoot = "/mnt/e/autodl_results"
$localComplexArgs = @(
    "--dataset-stems", "e10fV1",
    "--experiment-name", $ExperimentName,
    "--residual-forms", "projected_kv", "projected_vlambda",
    "--epochs", [string]$Epochs,
    "--file-type", ".mat",
    "--data-root", "/mnt/e/DataPons_processed",
    "--results-root", $wslResultsRoot,
    "--selected-device", "gpu"
)

if ($Resume) {
    $localComplexArgs += "--resume"
} else {
    $localComplexArgs += "--fresh-checkpoints"
}

if ($ContinueOnError) {
    $localComplexArgs += "--continue-on-error"
}

$quotedLocalComplexArgs = ($localComplexArgs | ForEach-Object {
    "'" + ($_ -replace "'", "'\\''") + "'"
}) -join " "
$wslCommand = "cd '$wslRepoRoot/python_scripts/autodl_local_e10gb1' && $WslPython run_local_complex_split_batch.py $quotedLocalComplexArgs"

Write-Host ""
Write-Host ("=" * 80)
Write-Warning "Legacy launcher: run_E10fV1_4conditions_60epoch.ps1 is kept only for archived four-condition reproduction."
Write-Host "E10fV1 four-condition BLP ResKoopNet launcher"
Write-Host ("ExperimentName: {0}" -f $ExperimentName)
Write-Host ("Epochs: {0}" -f $Epochs)
Write-Host ("Local download root: {0}" -f $LocalDownloadRoot)
Write-Host ("Resume mode: {0}" -f ([bool]$Resume))
Write-Host ("ContinueOnError: {0}" -f ([bool]$ContinueOnError))
Write-Host ("=" * 80)

if (-not $OnlyComplexSplit) {
    Write-Host ""
    Write-Host ("-" * 80)
    Write-Host "Starting AutoDL ABS branch:"
    Write-Host "  abs + projected_kv"
    Write-Host "  abs + projected_vlambda"
    Write-Host ("Remote data subdir: {0}" -f $RemoteDataSubdirAbs)
    Write-Host ("-" * 80)

    & $pythonExe @commonRemoteArgs
    if ($LASTEXITCODE -ne 0) {
        throw "AutoDL ABS branch failed with exit code $LASTEXITCODE"
    }
}

if (-not $OnlyAbs) {
    Write-Host ""
    Write-Host ("-" * 80)
    Write-Host "Starting local WSL complex_split branch:"
    Write-Host "  complex_split + projected_kv"
    Write-Host "  complex_split + projected_vlambda"
    Write-Host ("WSL command: {0}" -f $wslCommand)
    Write-Host ("-" * 80)

    & wsl.exe bash -lc $wslCommand
    if ($LASTEXITCODE -ne 0) {
        throw "Local WSL complex_split branch failed with exit code $LASTEXITCODE"
    }
}

Write-Host ""
Write-Host "All requested E10fV1 conditions finished successfully."
