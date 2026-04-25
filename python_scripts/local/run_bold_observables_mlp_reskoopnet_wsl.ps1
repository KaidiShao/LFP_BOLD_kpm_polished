param(
    [string]$ExperimentName = "bold_wsl_20260423",
    [string[]]$ObservableModes = @("eleHP", "HP", "roi_mean", "slow_band_power", "svd", "HP_svd100", "global_svd100"),
    [string[]]$ResidualForms = @("projected_kv", "projected_vlambda"),
    [ValidateSet("e10gb1", "e10gh1", "e10fV1", "f12m01")]
    [string]$StartWithDataset = "e10gb1",
    [switch]$OnlyStartDataset,
    [switch]$Resume,
    [switch]$ForceRerun,
    [switch]$ContinueOnError,
    [switch]$ExportPsi,
    [ValidateSet("cpu", "gpu")]
    [string]$SelectedDevice = "gpu",
    [int]$Epochs = 30,
    [int]$BatchSize = 2000,
    [int]$ChunkSize = 5000,
    [int]$NPsiTrain = 100,
    [double]$TrainRatio = 0.7,
    [double]$Reg = 0.1,
    [double]$LearningRate = 1e-4,
    [int[]]$LayerSizes = @(100, 100, 100),
    [string]$PythonExe = "python3",
    [int]$MaxParallel = 1,
    [int]$BlasThreads = 2,
    [int]$TfIntraopThreads = 2,
    [int]$TfInteropThreads = 1,
    [string]$Distro = "",
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$scriptPath = "/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/python_scripts/local/run_bold_observables_mlp_reskoopnet_wsl.sh"

$envParts = @(
    "EXPERIMENT_NAME='$ExperimentName'",
    "OBSERVABLE_MODES='$($ObservableModes -join ' ')'",
    "RESIDUAL_FORMS='$($ResidualForms -join ' ')'",
    "START_WITH_DATASET='$StartWithDataset'",
    "ONLY_START_DATASET='$([int][bool]$OnlyStartDataset)'",
    "RESUME='$([int][bool]$Resume)'",
    "FORCE_RERUN='$([int][bool]$ForceRerun)'",
    "CONTINUE_ON_ERROR='$([int][bool]$ContinueOnError)'",
    "EXPORT_PSI='$([int][bool]$ExportPsi)'",
    "SELECTED_DEVICE='$SelectedDevice'",
    "EPOCHS='$Epochs'",
    "BATCH_SIZE='$BatchSize'",
    "CHUNK_SIZE='$ChunkSize'",
    "N_PSI_TRAIN='$NPsiTrain'",
    "TRAIN_RATIO='$TrainRatio'",
    "REG='$Reg'",
    "LR='$LearningRate'",
    "LAYER_SIZES='$($LayerSizes -join ' ')'",
    "PYTHON_EXE='$PythonExe'",
    "MAX_PARALLEL='$MaxParallel'",
    "BLAS_THREADS='$BlasThreads'",
    "TF_INTRAOP_THREADS='$TfIntraopThreads'",
    "TF_INTEROP_THREADS='$TfInteropThreads'",
    "DRY_RUN='$([int][bool]$DryRun)'"
)

$bashCommand = "cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished && " +
    ($envParts -join " ") + " bash '$scriptPath'"

if ([string]::IsNullOrWhiteSpace($Distro)) {
    & wsl.exe bash -lc $bashCommand
} else {
    & wsl.exe -d $Distro bash -lc $bashCommand
}

if ($LASTEXITCODE -ne 0) {
    throw "WSL BOLD MLP ResKoopNet launcher failed with exit code $LASTEXITCODE"
}
