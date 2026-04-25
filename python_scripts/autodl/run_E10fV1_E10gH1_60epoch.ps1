param(
    [string]$ExperimentName = "real_batch60_20260418",
    [switch]$Resume,
    [switch]$UseLocalUpload,
    [int]$RemoteInputWaitHours = 72,
    [int]$RemoteInputPollSeconds = 60,
    [ValidateSet("e10fV1", "e10gh1")]
    [string]$StartWithDataset = "e10fV1",
    [switch]$OnlyStartDataset,
    [string]$RemoteDataSubdirE10fV1 = "e10fV1_abs",
    [string]$RemoteDataSubdirE10gH1 = "e10gH1_abs"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished"
$controller = Join-Path $repoRoot "python_scripts\autodl\dataset_batch_controller_autodl_reskoopnet_mlp.py"
$pythonExe = "python"

$commonArgs = @(
    $controller,
    "--ssh-host", "connect.westb.seetacloud.com",
    "--ssh-user", "root",
    "--ssh-port", "19241",
    "--ssh-key", "$HOME\.ssh\id_ed25519_autodl",
    "--observable-modes", "abs",
    "--residual-forms", "projected_kv", "projected_vlambda",
    "--file-type", ".h5",
    "--remote-python", "/root/miniconda3/envs/reskoopnet/bin/python",
    "--epochs", "60",
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

$datasets = @(
    @{
        DatasetStem = "e10fV1"
        RemoteDataSubdir = $RemoteDataSubdirE10fV1
        LocalDownloadRoot = "E:\autodl_results\e10fV1\mlp"
        LocalDataFileAbs = "E:\DataPons_processed\E10fV1\reskoopnet_dictionary\e10fV1_low50_high250_g2_abs_single.mat"
        LocalObsInfoFileAbs = "E:\DataPons_processed\E10fV1\reskoopnet_dictionary\e10fV1_low50_high250_g2_abs_single_obs_info.csv"
        LocalDataFileComplexSplit = "E:\DataPons_processed\E10fV1\reskoopnet_dictionary\e10fV1_low50_high250_g2_complex_split_single.mat"
        LocalObsInfoFileComplexSplit = "E:\DataPons_processed\E10fV1\reskoopnet_dictionary\e10fV1_low50_high250_g2_complex_split_single_obs_info.csv"
    },
    @{
        DatasetStem = "e10gh1"
        RemoteDataSubdir = $RemoteDataSubdirE10gH1
        LocalDownloadRoot = "E:\autodl_results\e10gh1\mlp"
        LocalDataFileAbs = "E:\DataPons_processed\E10gH1\reskoopnet_dictionary\e10gh1_low50_high250_g2_abs_single.mat"
        LocalObsInfoFileAbs = "E:\DataPons_processed\E10gH1\reskoopnet_dictionary\e10gh1_low50_high250_g2_abs_single_obs_info.csv"
        LocalDataFileComplexSplit = "E:\DataPons_processed\E10gH1\reskoopnet_dictionary\e10gh1_low50_high250_g2_complex_split_single.mat"
        LocalObsInfoFileComplexSplit = "E:\DataPons_processed\E10gH1\reskoopnet_dictionary\e10gh1_low50_high250_g2_complex_split_single_obs_info.csv"
    }
)

if ($StartWithDataset -eq "e10gh1") {
    $datasets = @($datasets[1], $datasets[0])
}

if ($OnlyStartDataset) {
    $datasets = @($datasets | Where-Object { $_.DatasetStem -eq $StartWithDataset })
}

foreach ($dataset in $datasets) {
    $datasetExperimentName = $ExperimentName
    if ($dataset.DatasetStem -ne "e10fV1") {
        $datasetExperimentName = "{0}_{1}" -f $ExperimentName, $dataset.DatasetStem
    }

    Write-Host ""
    Write-Host ("=" * 80)
    Write-Host ("Starting dataset: {0}" -f $dataset.DatasetStem)
    Write-Host ("Base experiment name: {0}" -f $ExperimentName)
    Write-Host ("Dataset experiment name: {0}" -f $datasetExperimentName)
    Write-Host ("Remote data subdir: {0}" -f $dataset.RemoteDataSubdir)
    Write-Host ("Local download root: {0}" -f $dataset.LocalDownloadRoot)
    if ($UseLocalUpload) {
        Write-Host "Input mode: upload from local via controller"
    } else {
        Write-Host "Input mode: reuse pre-uploaded remote input and wait if it is not there yet"
        Write-Host ("Expected remote data file: /root/autodl-tmp/data/{0}/{1}" -f $dataset.RemoteDataSubdir, [System.IO.Path]::GetFileName($dataset.LocalDataFileAbs))
    }
    Write-Host ("=" * 80)

    $datasetArgs = @(
        "--dataset-stem", $dataset.DatasetStem,
        "--experiment-name", $datasetExperimentName,
        "--remote-data-subdir", $dataset.RemoteDataSubdir,
        "--local-download-root", $dataset.LocalDownloadRoot,
        "--local-data-file-abs", $dataset.LocalDataFileAbs,
        "--local-obs-info-file-abs", $dataset.LocalObsInfoFileAbs,
        "--local-data-file-complex-split", $dataset.LocalDataFileComplexSplit,
        "--local-obs-info-file-complex-split", $dataset.LocalObsInfoFileComplexSplit
    )

    & $pythonExe @commonArgs @datasetArgs
    if ($LASTEXITCODE -ne 0) {
        throw "Dataset batch failed for $($dataset.DatasetStem) with exit code $LASTEXITCODE"
    }
}

Write-Host ""
Write-Host "All dataset batches finished successfully."
