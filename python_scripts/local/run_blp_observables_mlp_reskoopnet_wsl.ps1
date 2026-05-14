param(
    [string]$ExperimentName = "blp_vlambda_mainline_torchlike_20260509",
    [string]$OutputRootBase = "E:\autodl_results_new",
    [string[]]$ObservableModes = @("abs", "complex_split"),
    [string[]]$ResidualForms = @("projected_vlambda"),
    [ValidateSet("e10gh1", "e10fV1", "e10gw1", "f12m01")]
    [string]$StartWithDataset = "e10gh1",
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
    [double]$Reg = 0.3,
    [double]$LearningRate = 1e-5,
    [int]$InnerEpochs = 2,
    [int[]]$LayerSizes = @(100, 100, 100),
    [int]$Seed = 1234,
    [switch]$TrainShuffle,
    [switch]$TorchLikeSync,
    [string]$PythonExe = "/home/kdshao/anaconda3/bin/python",
    [int]$MaxParallel = 1,
    [int]$BlasThreads = 2,
    [int]$TfIntraopThreads = 2,
    [int]$TfInteropThreads = 1,
    [string]$Distro = "Ubuntu-22.04",
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$scriptPath = "/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/python_scripts/local/run_blp_observables_mlp_reskoopnet_wsl.sh"

function Convert-ToWslPath {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Path
    )

    $fullPath = [System.IO.Path]::GetFullPath($Path)
    if ($fullPath -match '^(?<drive>[A-Za-z]):\\(?<rest>.*)$') {
        $drive = $Matches.drive.ToLowerInvariant()
        $rest = ($Matches.rest -replace '\\', '/')
        if ([string]::IsNullOrWhiteSpace($rest)) {
            return "/mnt/$drive"
        }
        return "/mnt/$drive/$rest"
    }

    throw "Cannot convert path to WSL form: $Path"
}

$outputRootBaseWsl = Convert-ToWslPath $OutputRootBase
$effectiveTrainShuffle = if ($PSBoundParameters.ContainsKey("TrainShuffle")) {
    [bool]($TrainShuffle.IsPresent)
} else {
    $true
}
$effectiveTorchLikeSync = if ($PSBoundParameters.ContainsKey("TorchLikeSync")) {
    [bool]($TorchLikeSync.IsPresent)
} else {
    $true
}

$envParts = @(
    "EXPERIMENT_NAME='$ExperimentName'",
    "LOCAL_OUTPUT_ROOT_BASE='$outputRootBaseWsl'",
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
    "INNER_EPOCHS='$InnerEpochs'",
    "LAYER_SIZES='$($LayerSizes -join ' ')'",
    "SEED='$Seed'",
    "TRAIN_SHUFFLE='$([int]$effectiveTrainShuffle)'",
    "SPECTRAL_SYNC_MODE='$(if ($effectiveTorchLikeSync) { "pre_only" } else { "dual" })'",
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
    throw "WSL BLP MLP ResKoopNet launcher failed with exit code $LASTEXITCODE"
}
