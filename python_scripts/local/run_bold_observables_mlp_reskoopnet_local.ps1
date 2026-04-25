param(
    [string]$ExperimentName = "bold_local_20260423",
    [string]$LocalOutputRootBase = "E:\autodl_results_local\bold",
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
    [string]$PythonExe = "python",
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished"
$runner = Join-Path $repoRoot "python_scripts\autodl\run_autodl_reskoopnet_mlp.py"
$solverDir = Join-Path $repoRoot "python_scripts\autodl"
$boldRoot = "E:\DataPons_processed\bold_observables"

function Get-ObservableFile {
    param(
        [string]$DatasetId,
        [string]$Mode
    )

    return Join-Path (Join-Path $boldRoot $DatasetId) ("{0}_bold_observables_{1}.mat" -f $DatasetId, $Mode)
}

function Test-LocalRunComplete {
    param(
        [string]$OutputParent,
        [string]$RunLabel
    )

    $outputDir = Join-Path $OutputParent $RunLabel
    if (-not (Test-Path -LiteralPath $outputDir)) {
        return $false
    }

    $summaryCount = (Get-ChildItem -LiteralPath $outputDir -Filter "*_summary.mat" -ErrorAction SilentlyContinue | Measure-Object).Count
    $chunkCount = (Get-ChildItem -LiteralPath $outputDir -Filter "*_outputs_*.mat" -ErrorAction SilentlyContinue |
        Where-Object { $_.Name -notlike "*_outputs_Psi_*" } | Measure-Object).Count

    return ($summaryCount -gt 0 -and $chunkCount -gt 0)
}

function Get-RunLabel {
    param(
        [string]$ExperimentNameForDataset,
        [string]$ResidualForm,
        [string]$ObservableMode
    )

    return "mlp_obs_{0}_{1}_{2}" -f $ExperimentNameForDataset, $ResidualForm, $ObservableMode
}

$datasets = @(
    @{
        DatasetStem = "e10gb1"
        DatasetId = "E10.gb1"
    },
    @{
        DatasetStem = "e10gh1"
        DatasetId = "E10.gH1"
    },
    @{
        DatasetStem = "e10fV1"
        DatasetId = "E10.fV1"
    },
    @{
        DatasetStem = "f12m01"
        DatasetId = "F12.m01"
    }
)

$startIndex = -1
for ($i = 0; $i -lt $datasets.Count; $i++) {
    if ($datasets[$i].DatasetStem -eq $StartWithDataset) {
        $startIndex = $i
        break
    }
}
if ($startIndex -lt 0) {
    throw "Unknown StartWithDataset: $StartWithDataset"
}
if ($startIndex -gt 0) {
    $datasets = @($datasets[$startIndex..($datasets.Count - 1)] + $datasets[0..($startIndex - 1)])
}

if ($OnlyStartDataset) {
    $datasets = @($datasets | Where-Object { $_.DatasetStem -eq $StartWithDataset })
}

if (Test-Path -LiteralPath $runner -PathType Leaf) {
    Write-Host "Using local runner: $runner"
} else {
    throw "Local runner not found: $runner"
}

foreach ($dataset in $datasets) {
    foreach ($mode in $ObservableModes) {
        $file = Get-ObservableFile -DatasetId $dataset.DatasetId -Mode $mode
        if (-not (Test-Path -LiteralPath $file -PathType Leaf)) {
            throw "Missing BOLD observable file for $($dataset.DatasetId) / ${mode}: $file"
        }
    }
}

foreach ($dataset in $datasets) {
    $datasetStem = $dataset.DatasetStem
    $datasetId = $dataset.DatasetId
    $datasetExperimentName = "{0}_{1}" -f $ExperimentName, $datasetStem
    $datasetResultRoot = Join-Path $LocalOutputRootBase $datasetStem
    $outputParent = Join-Path (Join-Path $datasetResultRoot "mlp") "outputs"
    $checkpointParent = Join-Path (Join-Path $datasetResultRoot "mlp") "checkpoints"
    $logParent = Join-Path (Join-Path $datasetResultRoot "mlp") "logs"
    $consoleLogParent = Join-Path (Join-Path $datasetResultRoot "mlp") "console_logs"

    if (-not $DryRun) {
        foreach ($path in @($outputParent, $checkpointParent, $logParent, $consoleLogParent)) {
            if (-not (Test-Path -LiteralPath $path)) {
                New-Item -ItemType Directory -Path $path | Out-Null
            }
        }
    }

    foreach ($mode in $ObservableModes) {
        foreach ($residualForm in $ResidualForms) {
            $runLabel = Get-RunLabel -ExperimentNameForDataset $datasetExperimentName -ResidualForm $residualForm -ObservableMode $mode
            if ((-not $ForceRerun) -and (Test-LocalRunComplete -OutputParent $outputParent -RunLabel $runLabel)) {
                Write-Host ""
                Write-Host ("Skipping complete local run: {0}" -f $runLabel)
                continue
            }

            $dataFile = Get-ObservableFile -DatasetId $datasetId -Mode $mode
            $consoleLog = Join-Path $consoleLogParent ("{0}.log" -f $runLabel)

            $cmdArgs = @(
                $runner,
                "--project-root", $repoRoot,
                "--solver-dir", $solverDir,
                "--data-root", $boldRoot,
                "--output-parent", $outputParent,
                "--checkpoint-parent", $checkpointParent,
                "--log-parent", $logParent,
                "--run-name-base", "mlp_obs",
                "--experiment-name", $datasetExperimentName,
                "--selected-device", $SelectedDevice,
                "--solver-name", "resdmd_batch",
                "--residual-form", $residualForm,
                "--data-subdir", $datasetId,
                "--dataset-stem", $datasetStem,
                "--observable-mode", $mode,
                "--data-filename", (Split-Path -Leaf $dataFile),
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
                "--inner-epochs", "5",
                "--end-condition", "1e-9",
                "--chunk-size", [string]$ChunkSize,
                "--layer-sizes"
            )
            $cmdArgs += [string[]]$LayerSizes

            if ($Resume) {
                $cmdArgs += "--resume"
            } else {
                $cmdArgs += "--fresh-checkpoints"
            }
            if ($ExportPsi) {
                $cmdArgs += "--export-psi"
            }

            Write-Host ""
            Write-Host ("=" * 80)
            Write-Host ("Starting local BOLD MLP ResKoopNet: {0}" -f $runLabel)
            Write-Host ("Data file: {0}" -f $dataFile)
            Write-Host ("Console log: {0}" -f $consoleLog)
            Write-Host ("=" * 80)

            if ($DryRun) {
                Write-Host ("Dry run command:")
                Write-Host ("  {0} {1}" -f $PythonExe, ($cmdArgs -join " "))
                continue
            }

            & $PythonExe @cmdArgs 2>&1 | Tee-Object -FilePath $consoleLog
            $exitCode = $LASTEXITCODE

            if ($exitCode -ne 0) {
                $message = "Local BOLD MLP ResKoopNet failed for ${runLabel} with exit code ${exitCode}"
                if ($ContinueOnError) {
                    Write-Warning $message
                    continue
                }
                throw $message
            }
        }
    }
}

if ($DryRun) {
    Write-Host ""
    Write-Host "Dry run finished successfully. No local jobs were launched."
} else {
    Write-Host ""
    Write-Host "All requested local BOLD MLP ResKoopNet runs finished successfully."
}
