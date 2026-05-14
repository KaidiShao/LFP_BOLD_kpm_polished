param(
    [string]$DatasetStem = "e10gb1",
    [string]$ExperimentName = "blp_vlambda_mainline_torchlike_20260509",
    [string]$DataRootBase = "E:\DataPons_processed",
    [string]$DataRoot = "",
    [string]$DataSubdir = "pipeline1_reskoopnet_dictionary",
    [string]$AbsFilename = "",
    [string]$ComplexSplitFilename = "",
    [string]$OutputRootBase = "E:\autodl_results_new",
    [string]$RepoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished",
    [string]$WslPython = "/home/kdshao/anaconda3/bin/python",
    [string]$Distro = "Ubuntu-22.04",
    [ValidateSet("cpu", "gpu")]
    [string]$SelectedDevice = "gpu",
    [string]$SolverName = "resdmd_batch",
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
    [string[]]$ObservableModes = @("abs", "complex_split"),
    [ValidateSet("full", "summary_only")]
    [string]$ExportMode = "full",
    [ValidateSet("single", "double")]
    [string]$ExportPrecision = "double",
    [int]$DiagnosticDpi = 120,
    [switch]$SkipDiagnosticPdf,
    [switch]$Resume,
    [switch]$ForceRerun,
    [switch]$ExportPsi,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

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

function Test-RunComplete {
    param(
        [string]$OutputParent,
        [string]$RunLabel,
        [string]$ExportMode
    )

    $outputDir = Join-Path $OutputParent $RunLabel
    if (-not (Test-Path -LiteralPath $outputDir)) {
        return $false
    }

    $summaryCount = (Get-ChildItem -LiteralPath $outputDir -Filter "*_summary.mat" -ErrorAction SilentlyContinue | Measure-Object).Count
    $chunkCount = (Get-ChildItem -LiteralPath $outputDir -Filter "*_outputs_*.mat" -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -notlike "*_outputs_Psi_*" } | Measure-Object).Count

    if ($ExportMode -eq "summary_only") {
        return ($summaryCount -gt 0)
    }

    return ($summaryCount -gt 0 -and $chunkCount -gt 0)
}

if ([string]::IsNullOrWhiteSpace($DataRoot)) {
    $DataRoot = Join-Path $DataRootBase $DatasetStem
}

if ([string]::IsNullOrWhiteSpace($AbsFilename)) {
    $AbsFilename = "{0}_low50_high250_g2_abs_single.mat" -f $DatasetStem
}
if ([string]::IsNullOrWhiteSpace($ComplexSplitFilename)) {
    $ComplexSplitFilename = "{0}_low50_high250_g2_complex_split_single.mat" -f $DatasetStem
}

$runner = Join-Path $RepoRoot "python_scripts\autodl\run_autodl_reskoopnet_mlp.py"
$solverDir = Join-Path $RepoRoot "python_scripts\autodl"
$dataDir = Join-Path $DataRoot $DataSubdir
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

$datasetExperimentName = "{0}_{1}_seed{2}" -f $ExperimentName, $DatasetStem, $Seed
$datasetResultRoot = Join-Path (Join-Path $OutputRootBase $DatasetStem) "mlp"
$outputParent = Join-Path $datasetResultRoot "outputs"
$checkpointParent = Join-Path $datasetResultRoot "checkpoints"
$logParent = Join-Path $datasetResultRoot "logs"
$consoleLogParent = Join-Path $datasetResultRoot "console_logs"

foreach ($path in @($runner, $solverDir, $dataDir)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required path not found: $path"
    }
}

$allowedObservableModes = @("abs", "complex_split")
foreach ($observableMode in $ObservableModes) {
    if ($allowedObservableModes -notcontains $observableMode) {
        throw "Unsupported observable mode: $observableMode"
    }
}

$allRuns = @(
    @{
        ObservableMode = "abs"
        FileName = $AbsFilename
    },
    @{
        ObservableMode = "complex_split"
        FileName = $ComplexSplitFilename
    }
)

$runs = @($allRuns | Where-Object { $ObservableModes -contains $_.ObservableMode })

foreach ($run in $runs) {
    $filePath = Join-Path $dataDir $run.FileName
    if (-not (Test-Path -LiteralPath $filePath -PathType Leaf)) {
        throw "Missing observable file: $filePath"
    }
}

if (-not $DryRun) {
    foreach ($path in @($outputParent, $checkpointParent, $logParent, $consoleLogParent)) {
        if (-not (Test-Path -LiteralPath $path)) {
            New-Item -ItemType Directory -Path $path | Out-Null
        }
    }
}

$repoRootWsl = Convert-ToWslPath $RepoRoot
$runnerWsl = Convert-ToWslPath $runner
$solverDirWsl = Convert-ToWslPath $solverDir
$dataRootWsl = Convert-ToWslPath $DataRoot
$outputParentWsl = Convert-ToWslPath $outputParent
$checkpointParentWsl = Convert-ToWslPath $checkpointParent
$logParentWsl = Convert-ToWslPath $logParent

foreach ($run in $runs) {
    $observableMode = $run.ObservableMode
    $fileName = $run.FileName
    $runLabel = "mlp_obs_{0}_projected_vlambda_{1}" -f $datasetExperimentName, $observableMode
    $effectiveSpectralSyncMode = if ($effectiveTorchLikeSync -or $SolverName -eq "resdmd_batch_static") {
        "pre_only"
    } else {
        "dual"
    }

    if ((-not $ForceRerun) -and (Test-RunComplete -OutputParent $outputParent -RunLabel $runLabel -ExportMode $ExportMode)) {
        Write-Host ""
        Write-Host ("Skipping complete run: {0}" -f $runLabel)
        continue
    }

    $consoleLog = Join-Path $consoleLogParent ("{0}.log" -f $runLabel)
    $wslArgs = @(
        "-d", $Distro,
        $WslPython,
        $runnerWsl,
        "--project-root", $repoRootWsl,
        "--solver-dir", $solverDirWsl,
        "--data-root", $dataRootWsl,
        "--output-parent", $outputParentWsl,
        "--checkpoint-parent", $checkpointParentWsl,
        "--log-parent", $logParentWsl,
        "--run-name-base", "mlp_obs",
        "--experiment-name", $datasetExperimentName,
        "--selected-device", $SelectedDevice,
        "--solver-name", $SolverName,
        "--residual-form", "projected_vlambda",
        "--training-policy", "float64",
        "--analysis-dtype", "float64",
        "--gram-dtype", "float64",
        "--spectral-dtype", "float64",
        "--data-subdir", $DataSubdir,
        "--dataset-stem", $DatasetStem,
        "--observable-mode", $observableMode,
        "--data-filename", $fileName,
        "--file-type", ".mat",
        "--field-name", "obs",
        "--n-psi-train", [string]$NPsiTrain,
        "--train-ratio", [string]$TrainRatio,
        "--reg", [string]$Reg,
        "--rounds", "1",
        "--epochs", [string]$Epochs,
        "--batch-size", [string]$BatchSize,
        "--lr", [string]$LearningRate,
        "--log-interval", "1",
        "--lr-decay-factor", "0.8",
        "--inner-epochs", [string]$InnerEpochs,
        "--end-condition", "1e-9",
        "--seed", [string]$Seed,
        "--export-mode", $ExportMode,
        "--export-precision", $ExportPrecision,
        "--diagnostic-dpi", [string]$DiagnosticDpi,
        "--chunk-size", [string]$ChunkSize,
        "--layer-sizes"
    )
    $wslArgs += [string[]]$LayerSizes

    if ($Resume) {
        $wslArgs += "--resume"
    } else {
        $wslArgs += "--fresh-checkpoints"
    }
    if ($effectiveTrainShuffle) {
        $wslArgs += "--train-shuffle"
    }
    if ($effectiveSpectralSyncMode -eq "pre_only") {
        $wslArgs += @("--spectral-sync-mode", "pre_only")
    }
    if ($ExportPsi) {
        $wslArgs += "--export-psi"
    }
    if ($SkipDiagnosticPdf) {
        $wslArgs += "--skip-diagnostic-pdf"
    }

    Write-Host ""
    Write-Host ("=" * 80)
    Write-Host ("Starting BLP projected_vlambda run: {0}" -f $runLabel)
    Write-Host ("Observable mode: {0}" -f $observableMode)
    Write-Host ("Seed: {0}" -f $Seed)
    Write-Host ("Solver: {0}" -f $SolverName)
    Write-Host ("Hyperparameters: lr={0}, reg={1}, inner_epochs={2}" -f $LearningRate, $Reg, $InnerEpochs)
    Write-Host ("Training schedule: shuffle={0}, spectral_sync_mode={1}" -f $effectiveTrainShuffle, $effectiveSpectralSyncMode)
    Write-Host ("Export policy: mode={0}, precision={1}, dpi={2}, skip_pdf={3}" -f $ExportMode, $ExportPrecision, $DiagnosticDpi, $SkipDiagnosticPdf.IsPresent)
    Write-Host ("Data file: {0}" -f (Join-Path $dataDir $fileName))
    Write-Host ("Console log: {0}" -f $consoleLog)
    Write-Host ("=" * 80)

    if ($DryRun) {
        Write-Host "Dry run command:"
        Write-Host ("  wsl.exe {0}" -f ($wslArgs -join " "))
        continue
    }

    $stdoutLog = "{0}.stdout.tmp" -f $consoleLog
    $stderrLog = "{0}.stderr.tmp" -f $consoleLog
    foreach ($tempLog in @($stdoutLog, $stderrLog, $consoleLog)) {
        if (Test-Path -LiteralPath $tempLog) {
            Remove-Item -LiteralPath $tempLog -Force
        }
    }

    $process = Start-Process -FilePath "wsl.exe" `
        -ArgumentList $wslArgs `
        -NoNewWindow `
        -Wait `
        -PassThru `
        -RedirectStandardOutput $stdoutLog `
        -RedirectStandardError $stderrLog
    $exitCode = $process.ExitCode

    $combinedOutput = @()
    foreach ($tempLog in @($stdoutLog, $stderrLog)) {
        if (Test-Path -LiteralPath $tempLog) {
            $combinedOutput += Get-Content -LiteralPath $tempLog
        }
    }
    if ($combinedOutput.Count -gt 0) {
        $combinedOutput | Tee-Object -FilePath $consoleLog
    } else {
        New-Item -ItemType File -Path $consoleLog -Force | Out-Null
    }

    foreach ($tempLog in @($stdoutLog, $stderrLog)) {
        if (Test-Path -LiteralPath $tempLog) {
            Remove-Item -LiteralPath $tempLog -Force
        }
    }

    if ($exitCode -ne 0) {
        throw "Run failed for ${runLabel} with exit code ${exitCode}"
    }
}

if ($DryRun) {
    Write-Host ""
    Write-Host "Dry run finished successfully. No jobs were launched."
} else {
    Write-Host ""
    Write-Host "All requested projected_vlambda BLP runs finished successfully."
}
