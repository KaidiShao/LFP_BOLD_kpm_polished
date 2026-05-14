param(
    [string]$DatasetStem = "e10gb1",
    [string]$ExperimentName = "blp_vlambda_seed_sweep_20260506",
    [string]$DataRootBase = "E:\DataPons_processed",
    [string]$DataRoot = "",
    [string]$DataSubdir = "pipeline1_reskoopnet_dictionary",
    [string]$AbsFilename = "",
    [string]$ComplexSplitFilename = "",
    [string]$OutputRootBase = "E:\autodl_results",
    [string]$RepoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished",
    [string]$WslPython = "/home/kdshao/anaconda3/bin/python",
    [string]$Distro = "Ubuntu-22.04",
    [ValidateSet("cpu", "gpu")]
    [string]$SelectedDevice = "gpu",
    [int]$Epochs = 12,
    [int]$BatchSize = 2000,
    [int]$ChunkSize = 5000,
    [int]$NPsiTrain = 100,
    [double]$TrainRatio = 0.7,
    [double]$Reg = 0.3,
    [double]$LearningRate = 1e-5,
    [int]$InnerEpochs = 2,
    [int[]]$LayerSizes = @(100, 100, 100),
    [int[]]$Seeds = @(1234, 2026, 3407, 7777, 9999),
    [string[]]$ObservableModes = @("abs"),
    [switch]$Resume,
    [switch]$ForceRerun,
    [switch]$ExportPsi,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$singleRunner = Join-Path $RepoRoot "python_scripts\local\run_one_blp_reskoopnet_vlambda_wsl.ps1"
if (-not (Test-Path -LiteralPath $singleRunner -PathType Leaf)) {
    throw "Required script not found: $singleRunner"
}

foreach ($seed in $Seeds) {
    Write-Host ""
    Write-Host ("#" * 80)
    Write-Host ("Starting seed sweep run for seed {0}" -f $seed)
    Write-Host ("#" * 80)

    $invokeArgs = @{
        DatasetStem = $DatasetStem
        ExperimentName = $ExperimentName
        DataRootBase = $DataRootBase
        DataSubdir = $DataSubdir
        OutputRootBase = $OutputRootBase
        RepoRoot = $RepoRoot
        WslPython = $WslPython
        Distro = $Distro
        SelectedDevice = $SelectedDevice
        Epochs = $Epochs
        BatchSize = $BatchSize
        ChunkSize = $ChunkSize
        NPsiTrain = $NPsiTrain
        TrainRatio = $TrainRatio
        Reg = $Reg
        LearningRate = $LearningRate
        InnerEpochs = $InnerEpochs
        LayerSizes = $LayerSizes
        Seed = $seed
        ObservableModes = $ObservableModes
    }

    if (-not [string]::IsNullOrWhiteSpace($DataRoot)) {
        $invokeArgs.DataRoot = $DataRoot
    }
    if (-not [string]::IsNullOrWhiteSpace($AbsFilename)) {
        $invokeArgs.AbsFilename = $AbsFilename
    }
    if (-not [string]::IsNullOrWhiteSpace($ComplexSplitFilename)) {
        $invokeArgs.ComplexSplitFilename = $ComplexSplitFilename
    }
    if ($Resume) {
        $invokeArgs.Resume = $true
    }
    if ($ForceRerun) {
        $invokeArgs.ForceRerun = $true
    }
    if ($ExportPsi) {
        $invokeArgs.ExportPsi = $true
    }
    if ($DryRun) {
        $invokeArgs.DryRun = $true
    }

    & $singleRunner @invokeArgs
}

Write-Host ""
Write-Host "Seed sweep finished."
