param(
    [string]$DatasetStem = "e10gb1",
    [string]$ExperimentName = "blp_vlambda_hparam_sweep_20260506",
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
    [int]$Epochs = 6,
    [int]$BatchSize = 2000,
    [int]$ChunkSize = 5000,
    [int]$NPsiTrain = 100,
    [double]$TrainRatio = 0.7,
    [int[]]$Seeds = @(1234, 2026),
    [double[]]$LearningRates = @(1e-5),
    [double[]]$Regs = @(0.2, 0.3, 0.5),
    [int[]]$InnerEpochsList = @(2),
    [int[]]$LayerSizes = @(100, 100, 100),
    [string[]]$ObservableModes = @("abs"),
    [ValidateSet("full", "summary_only")]
    [string]$ExportMode = "summary_only",
    [ValidateSet("single", "double")]
    [string]$ExportPrecision = "single",
    [int]$DiagnosticDpi = 96,
    [switch]$KeepDiagnosticPdf,
    [switch]$ExportPsi,
    [switch]$ForceRerun,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

function Format-TagValue {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Prefix,
        [Parameter(Mandatory = $true)]
        [string]$Value
    )

    $tag = $Value.ToLowerInvariant()
    $tag = $tag.Replace("+", "")
    $tag = $tag.Replace(".", "p")
    $tag = $tag.Replace("-", "m")
    return "{0}{1}" -f $Prefix, $tag
}

function Format-DoubleTag {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Prefix,
        [Parameter(Mandatory = $true)]
        [double]$Value
    )

    $raw = $Value.ToString("G", [System.Globalization.CultureInfo]::InvariantCulture)
    return Format-TagValue -Prefix $Prefix -Value $raw
}

$singleRunner = Join-Path $RepoRoot "python_scripts\local\run_one_blp_reskoopnet_vlambda_wsl.ps1"
if (-not (Test-Path -LiteralPath $singleRunner -PathType Leaf)) {
    throw "Required script not found: $singleRunner"
}

$totalRuns = $Seeds.Count * $LearningRates.Count * $Regs.Count * $InnerEpochsList.Count * $ObservableModes.Count
Write-Host ("Planned sweep: {0} total runs" -f $totalRuns)
Write-Host ("Seeds: {0}" -f (($Seeds | ForEach-Object { $_.ToString() }) -join ", "))
Write-Host ("LearningRates: {0}" -f (($LearningRates | ForEach-Object { $_.ToString("G", [System.Globalization.CultureInfo]::InvariantCulture) }) -join ", "))
Write-Host ("Regs: {0}" -f (($Regs | ForEach-Object { $_.ToString("G", [System.Globalization.CultureInfo]::InvariantCulture) }) -join ", "))
Write-Host ("InnerEpochs: {0}" -f (($InnerEpochsList | ForEach-Object { $_.ToString() }) -join ", "))
Write-Host ("ObservableModes: {0}" -f ($ObservableModes -join ", "))

foreach ($learningRate in $LearningRates) {
    foreach ($reg in $Regs) {
        foreach ($innerEpochs in $InnerEpochsList) {
            $comboExperimentName = "{0}_{1}_{2}_{3}" -f `
                $ExperimentName, `
                (Format-DoubleTag -Prefix "lr" -Value $learningRate), `
                (Format-DoubleTag -Prefix "reg" -Value $reg), `
                (Format-TagValue -Prefix "ie" -Value $innerEpochs.ToString())

            Write-Host ""
            Write-Host ("#" * 80)
            Write-Host ("Starting hyperparameter combo: {0}" -f $comboExperimentName)
            Write-Host ("lr={0}, reg={1}, inner_epochs={2}" -f $learningRate, $reg, $innerEpochs)
            Write-Host ("#" * 80)

            foreach ($seed in $Seeds) {
                $invokeArgs = @{
                    DatasetStem = $DatasetStem
                    ExperimentName = $comboExperimentName
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
                    Reg = $reg
                    LearningRate = $learningRate
                    InnerEpochs = $innerEpochs
                    LayerSizes = $LayerSizes
                    Seed = $seed
                    ObservableModes = $ObservableModes
                    ExportMode = $ExportMode
                    ExportPrecision = $ExportPrecision
                    DiagnosticDpi = $DiagnosticDpi
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
                if (-not $KeepDiagnosticPdf) {
                    $invokeArgs.SkipDiagnosticPdf = $true
                }
                if ($ExportPsi) {
                    $invokeArgs.ExportPsi = $true
                }
                if ($ForceRerun) {
                    $invokeArgs.ForceRerun = $true
                }
                if ($DryRun) {
                    $invokeArgs.DryRun = $true
                }

                & $singleRunner @invokeArgs
            }
        }
    }
}

Write-Host ""
Write-Host "Hyperparameter sweep finished."
