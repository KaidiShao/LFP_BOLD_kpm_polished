param(
    [string]$ExperimentName = "blp_mlp_20260430",
    [string]$SshHost = "connect.westb.seetacloud.com",
    [string]$SshUser = "root",
    [int]$SshPort = 19241,
    [string]$SshKey = "$HOME\.ssh\id_ed25519_autodl",
    [string]$RemoteCodeDir = "/root/autodl-tmp/code/reskoopnet_mlp",
    [string]$RemotePython = "/root/miniconda3/envs/reskoopnet/bin/python",
    [string]$RemoteDataSubdirPrefix = "blp",
    [string]$LocalDownloadRootBase = "E:\autodl_results\blp",
    [string[]]$ObservableModes = @("abs", "complex_split"),
    [string[]]$ResidualForms = @("projected_kv", "projected_vlambda"),
    [ValidateSet("e10gb1", "e10gh1", "e10fV1", "f12m01")]
    [string]$StartWithDataset = "e10gb1",
    [switch]$OnlyStartDataset,
    [switch]$Resume,
    [switch]$SkipCodeSync,
    [switch]$SkipDataUpload,
    [int]$RemoteInputWaitHours = 72,
    [int]$RemoteInputPollSeconds = 60,
    [int]$Epochs = 60,
    [int]$BatchSize = 2000,
    [int]$ChunkSize = 5000,
    [int]$NPsiTrain = 100,
    [double]$TrainRatio = 0.7,
    [double]$Reg = 0.1,
    [double]$LearningRate = 1e-4,
    [int[]]$LayerSizes = @(100, 100, 100),
    [switch]$ContinueOnError,
    [switch]$DownloadCheckpoints,
    [switch]$ExportPsi,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished"
$autodlDir = Join-Path $repoRoot "python_scripts\autodl"
$controller = Join-Path $autodlDir "dataset_batch_controller_autodl_reskoopnet_mlp.py"
$pythonExe = "python"
$blpRoot = "E:\DataPons_processed"

function Sync-RemoteCode {
    Write-Host "Syncing local MLP ResKoopNet code to ${SshUser}@${SshHost}:${RemoteCodeDir}"
    & ssh -p $SshPort -i $SshKey "${SshUser}@${SshHost}" "mkdir -p '$RemoteCodeDir'"
    if ($LASTEXITCODE -ne 0) {
        throw "Failed to create remote code directory: $RemoteCodeDir"
    }

    $codeFiles = @(
        "run_autodl_reskoopnet_mlp.py",
        "controller_autodl_reskoopnet_mlp.py",
        "batch_controller_autodl_reskoopnet_mlp.py",
        "dataset_batch_controller_autodl_reskoopnet_mlp.py",
        "multi_dataset_batch_controller_autodl_reskoopnet_mlp.py",
        "edmd_utils.py",
        "solver_resdmd_batch3.py",
        "solver_resdmd_batch2.py"
    )

    foreach ($fileName in $codeFiles) {
        $localFile = Join-Path $autodlDir $fileName
        if (Test-Path -LiteralPath $localFile) {
            & scp -P $SshPort -i $SshKey $localFile "${SshUser}@${SshHost}:${RemoteCodeDir}/$fileName"
            if ($LASTEXITCODE -ne 0) {
                throw "Failed to upload code file: $localFile"
            }
        }
    }
}

function Get-ObservableFile {
    param(
        [string]$DatasetStem,
        [string]$Mode
    )

    return Join-Path (Join-Path (Join-Path $blpRoot $DatasetStem) "pipeline1_reskoopnet_dictionary") ("{0}_low50_high250_g2_{1}_single.mat" -f $DatasetStem, $Mode)
}

function Get-ObservableInfoFile {
    param(
        [string]$DatasetStem,
        [string]$Mode
    )

    return Join-Path (Join-Path (Join-Path $blpRoot $DatasetStem) "pipeline1_reskoopnet_dictionary") ("{0}_low50_high250_g2_{1}_single_obs_info.csv" -f $DatasetStem, $Mode)
}

$datasets = @(
    @{ DatasetStem = "e10gb1" },
    @{ DatasetStem = "e10gh1" },
    @{ DatasetStem = "e10fV1" },
    @{ DatasetStem = "f12m01" }
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

foreach ($dataset in $datasets) {
    foreach ($mode in $ObservableModes) {
        $file = Get-ObservableFile -DatasetStem $dataset.DatasetStem -Mode $mode
        if (-not (Test-Path -LiteralPath $file)) {
            throw "Missing BLP observable file for $($dataset.DatasetStem) / ${mode}: $file"
        }
    }
}

if ((-not $SkipCodeSync) -and (-not $DryRun)) {
    Sync-RemoteCode
} elseif ($DryRun) {
    Write-Host "Dry run: skipping remote code sync."
}

$commonArgs = @(
    $controller,
    "--ssh-host", $SshHost,
    "--ssh-user", $SshUser,
    "--ssh-port", [string]$SshPort,
    "--ssh-key", $SshKey,
    "--remote-code-dir", $RemoteCodeDir,
    "--remote-python", $RemotePython,
    "--observable-modes"
)
$commonArgs += $ObservableModes
$commonArgs += @(
    "--residual-forms"
)
$commonArgs += $ResidualForms
$commonArgs += @(
    "--file-type", ".mat",
    "--field-name", "obs",
    "--epochs", [string]$Epochs,
    "--batch-size", [string]$BatchSize,
    "--chunk-size", [string]$ChunkSize,
    "--n-psi-train", [string]$NPsiTrain,
    "--train-ratio", [string]$TrainRatio,
    "--reg", [string]$Reg,
    "--lr", [string]$LearningRate,
    "--recover-completed-remote-runs",
    "--skip-completed-local-runs",
    "--delete-remote-run-after-download",
    "--delete-remote-input-after-observable",
    "--poll-interval", "30",
    "--tail-lines", "40",
    "--launch-timeout-sec", "30",
    "--layer-sizes"
)
$commonArgs += [string[]]$LayerSizes

if ($SkipDataUpload) {
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
if ($ContinueOnError) {
    $commonArgs += "--continue-on-error"
}
if ($DownloadCheckpoints) {
    $commonArgs += "--download-checkpoints"
}
if ($ExportPsi) {
    $commonArgs += "--export-psi"
}

foreach ($dataset in $datasets) {
    $datasetExperimentName = "{0}_{1}" -f $ExperimentName, $dataset.DatasetStem
    $remoteDataSubdir = "{0}_{1}" -f $RemoteDataSubdirPrefix, $dataset.DatasetStem
    $localDownloadRoot = Join-Path (Join-Path $LocalDownloadRootBase $dataset.DatasetStem) "mlp"

    $datasetArgs = @(
        "--dataset-stem", $dataset.DatasetStem,
        "--experiment-name", $datasetExperimentName,
        "--remote-data-subdir", $remoteDataSubdir,
        "--local-download-root", $localDownloadRoot,
        "--local-data-file-abs", (Get-ObservableFile -DatasetStem $dataset.DatasetStem -Mode "abs"),
        "--local-obs-info-file-abs", (Get-ObservableInfoFile -DatasetStem $dataset.DatasetStem -Mode "abs"),
        "--local-data-file-complex-split", (Get-ObservableFile -DatasetStem $dataset.DatasetStem -Mode "complex_split"),
        "--local-obs-info-file-complex-split", (Get-ObservableInfoFile -DatasetStem $dataset.DatasetStem -Mode "complex_split")
    )

    Write-Host ""
    Write-Host ("=" * 80)
    Write-Host ("Starting BLP MLP ResKoopNet dataset: {0}" -f $dataset.DatasetStem)
    Write-Host ("Experiment: {0}" -f $datasetExperimentName)
    Write-Host ("Observable modes: {0}" -f ($ObservableModes -join ", "))
    Write-Host ("Residual forms: {0}" -f ($ResidualForms -join ", "))
    Write-Host ("Local download root: {0}" -f $localDownloadRoot)
    Write-Host ("Remote data subdir: {0}" -f $remoteDataSubdir)
    Write-Host ("=" * 80)

    if ($DryRun) {
        Write-Host "Dry run command:"
        Write-Host ("  {0} {1}" -f $pythonExe, (($commonArgs + $datasetArgs) -join " "))
    } else {
        & $pythonExe @commonArgs @datasetArgs
        if ($LASTEXITCODE -ne 0) {
            throw "BLP MLP ResKoopNet batch failed for $($dataset.DatasetStem) with exit code $LASTEXITCODE"
        }
    }
}

Write-Host ""
if ($DryRun) {
    Write-Host "Dry run finished successfully. No remote jobs were launched."
} else {
    Write-Host "All requested BLP MLP ResKoopNet batches finished successfully."
}
