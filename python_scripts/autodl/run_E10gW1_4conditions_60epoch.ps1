<# 
LEGACY dataset-specific launcher for the archived E10gW1 four-condition run.

This is not a canonical pipeline 4 entry.
#>

param(
    [string]$ExperimentName = "e10gw1_batch60_20260428",
    [switch]$Resume,
    [switch]$UseLocalUpload,
    [int]$RemoteInputWaitHours = 72,
    [int]$RemoteInputPollSeconds = 60,
    [int]$Epochs = 60,
    [string]$RemoteDataSubdir = "e10gw1_4cond",
    [string]$LocalDownloadRoot = "E:\autodl_results\e10gw1\mlp"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished"
$controller = Join-Path $repoRoot "python_scripts\autodl\dataset_batch_controller_autodl_reskoopnet_mlp.py"
$pythonExe = "python"
$processedRoot = "E:\DataPons_processed"
$datasetStem = "e10gw1"
$dictDir = Join-Path (Join-Path $processedRoot $datasetStem) "pipeline1_reskoopnet_dictionary"

function Get-ObservableFile {
    param(
        [string]$Mode
    )

    return Join-Path $dictDir ("{0}_low50_high250_g2_{1}_single.mat" -f $datasetStem, $Mode)
}

function Get-ObservableInfoFile {
    param(
        [string]$Mode
    )

    return Join-Path $dictDir ("{0}_low50_high250_g2_{1}_single_obs_info.csv" -f $datasetStem, $Mode)
}

$absFile = Get-ObservableFile -Mode "abs"
$absInfoFile = Get-ObservableInfoFile -Mode "abs"
$complexFile = Get-ObservableFile -Mode "complex_split"
$complexInfoFile = Get-ObservableInfoFile -Mode "complex_split"

foreach ($path in @($controller, $absFile, $absInfoFile, $complexFile, $complexInfoFile)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required file not found: $path"
    }
}

$commonArgs = @(
    $controller,
    "--ssh-host", "connect.westb.seetacloud.com",
    "--ssh-user", "root",
    "--ssh-port", "19241",
    "--ssh-key", "$HOME\.ssh\id_ed25519_autodl",
    "--dataset-stem", $datasetStem,
    "--experiment-name", $ExperimentName,
    "--observable-modes", "abs", "complex_split",
    "--residual-forms", "projected_kv", "projected_vlambda",
    "--file-type", ".mat",
    "--field-name", "obs",
    "--remote-python", "/root/miniconda3/envs/reskoopnet/bin/python",
    "--epochs", [string]$Epochs,
    "--remote-data-subdir", $RemoteDataSubdir,
    "--local-download-root", $LocalDownloadRoot,
    "--local-data-file-abs", $absFile,
    "--local-obs-info-file-abs", $absInfoFile,
    "--local-data-file-complex-split", $complexFile,
    "--local-obs-info-file-complex-split", $complexInfoFile,
    "--recover-completed-remote-runs",
    "--skip-completed-local-runs",
    "--delete-remote-run-after-download",
    "--delete-remote-input-after-observable"
)

if (-not $UseLocalUpload) {
    $commonArgs += @(
        "--skip-data-upload",
        "--wait-for-remote-input",
        "--remote-input-timeout-sec", [string]($RemoteInputWaitHours * 3600),
        "--remote-input-poll-sec", [string]$RemoteInputPollSeconds
    )
}

if ($Resume) {
    $commonArgs += "--resume"
} else {
    $commonArgs += "--fresh-checkpoints"
}

Write-Host ""
Write-Host ("=" * 80)
Write-Warning "Legacy launcher: run_E10gW1_4conditions_60epoch.ps1 is kept only for archived dataset-specific reproduction."
Write-Host "Starting E10gW1 AutoDL batch"
Write-Host ("ExperimentName: {0}" -f $ExperimentName)
Write-Host ("Epochs: {0}" -f $Epochs)
Write-Host ("RemoteDataSubdir: {0}" -f $RemoteDataSubdir)
Write-Host ("LocalDownloadRoot: {0}" -f $LocalDownloadRoot)
if ($UseLocalUpload) {
    Write-Host "Input mode: upload observables from local via controller"
} else {
    Write-Host "Input mode: reuse pre-uploaded remote inputs and wait if needed"
}
Write-Host "Observable modes: abs, complex_split"
Write-Host "Residual forms: projected_kv, projected_vlambda"
Write-Host ("=" * 80)

& $pythonExe @commonArgs
if ($LASTEXITCODE -ne 0) {
    throw "E10gW1 AutoDL batch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "E10gW1 AutoDL batch finished successfully."
