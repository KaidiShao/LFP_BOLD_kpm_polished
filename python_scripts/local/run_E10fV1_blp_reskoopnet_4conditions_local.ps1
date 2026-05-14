<# 
LEGACY local one-off launcher for the archived E10fV1 four-condition BLP run.

This is not the canonical local pipeline 4 entry.
Prefer:
  python_scripts/local/run_blp_observables_mlp_reskoopnet_local.ps1
  python_scripts/local/run_blp_observables_mlp_reskoopnet_wsl.sh
#>

param(
    [string]$ExperimentName = "e10fV1_local4cond_run01",
    [int]$Epochs = 60,
    [switch]$Resume,
    [switch]$ForceRerun,
    [switch]$ContinueOnError,
    [ValidateSet("cpu", "gpu")]
    [string]$SelectedDevice = "gpu",
    [ValidateSet("abs", "complex_split")]
    [string[]]$ObservableMode,
    [ValidateSet("projected_kv", "projected_vlambda")]
    [string[]]$ResidualForm,
    [string]$WslDistro = "Ubuntu-22.04",
    [string]$WslPython = "/home/kdshao/anaconda3/bin/python",
    [string]$ResultsRoot = "/mnt/e/autodl_results",
    [string]$DataRoot = "/mnt/e/DataPons_processed"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished"
$wslScript = "/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/python_scripts/autodl_local_e10gb1/run_E10fV1_4conditions_60epoch_local.sh"

function Quote-BashValue {
    param([string]$Value)
    return "'" + ($Value -replace "'", "'\''") + "'"
}

$envParts = @(
    "EXPERIMENT_NAME=$(Quote-BashValue $ExperimentName)",
    "EPOCHS=$(Quote-BashValue ([string]$Epochs))",
    "RESUME=$(if ($Resume) { '1' } else { '0' })",
    "FORCE_RERUN=$(if ($ForceRerun) { '1' } else { '0' })",
    "CONTINUE_ON_ERROR=$(if ($ContinueOnError) { '1' } else { '0' })",
    "SELECTED_DEVICE=$(Quote-BashValue $SelectedDevice)",
    "PYTHON_EXE=$(Quote-BashValue $WslPython)",
    "RESULTS_ROOT=$(Quote-BashValue $ResultsRoot)",
    "DATA_ROOT=$(Quote-BashValue $DataRoot)"
)

if ($ObservableMode) {
    $envParts += "ONLY_OBSERVABLE_MODES=$(Quote-BashValue ($ObservableMode -join ','))"
}

if ($ResidualForm) {
    $envParts += "ONLY_RESIDUAL_FORMS=$(Quote-BashValue ($ResidualForm -join ','))"
}

$bashCommand = ($envParts -join " ") + " bash " + (Quote-BashValue $wslScript)

Write-Host ""
Write-Host ("=" * 80)
Write-Warning "Legacy launcher: run_E10fV1_blp_reskoopnet_4conditions_local.ps1 is kept only for archived four-condition reproduction."
Write-Host "Launching local E10fV1 BLP ResKoopNet 4-condition batch through WSL"
Write-Host ("ExperimentName: {0}" -f $ExperimentName)
Write-Host ("Epochs: {0}" -f $Epochs)
Write-Host ("SelectedDevice: {0}" -f $SelectedDevice)
Write-Host ("ObservableMode: {0}" -f $(if ($ObservableMode) { $ObservableMode -join ", " } else { "all" }))
Write-Host ("ResidualForm: {0}" -f $(if ($ResidualForm) { $ResidualForm -join ", " } else { "all" }))
Write-Host ("WslDistro: {0}" -f $WslDistro)
Write-Host ("WSL Python: {0}" -f $WslPython)
Write-Host ("=" * 80)

& wsl.exe -d $WslDistro bash -lc $bashCommand
if ($LASTEXITCODE -ne 0) {
    throw "WSL local E10fV1 BLP ResKoopNet batch failed with exit code $LASTEXITCODE"
}

Write-Host ""
Write-Host "Local E10fV1 BLP ResKoopNet 4-condition batch finished successfully."
