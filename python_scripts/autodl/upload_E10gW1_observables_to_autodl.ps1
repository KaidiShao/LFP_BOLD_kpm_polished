param(
    [string]$RemoteDataSubdir = "e10gw1_4cond",
    [string]$SshHost = "connect.westb.seetacloud.com",
    [int]$SshPort = 19241,
    [string]$SshUser = "root",
    [string]$SshKey = "$HOME\.ssh\id_ed25519_autodl",
    [string]$ProcessedRoot = "E:\DataPons_processed"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$datasetStem = "e10gw1"
$dictDir = Join-Path (Join-Path $ProcessedRoot $datasetStem) "pipeline1_reskoopnet_dictionary"
$remoteDataDir = "/root/autodl-tmp/data/$RemoteDataSubdir"

$files = @(
    (Join-Path $dictDir ("{0}_low50_high250_g2_abs_single.mat" -f $datasetStem))
    (Join-Path $dictDir ("{0}_low50_high250_g2_abs_single_obs_info.csv" -f $datasetStem))
    (Join-Path $dictDir ("{0}_low50_high250_g2_complex_split_single.mat" -f $datasetStem))
    (Join-Path $dictDir ("{0}_low50_high250_g2_complex_split_single_obs_info.csv" -f $datasetStem))
)

foreach ($path in $files) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required file not found: $path"
    }
}

$sshArgs = @(
    "-p", [string]$SshPort,
    "-i", $SshKey,
    ("{0}@{1}" -f $SshUser, $SshHost)
)

$scpBaseArgs = @(
    "-P", [string]$SshPort,
    "-i", $SshKey
)

Write-Host ""
Write-Host ("=" * 80)
Write-Host "Preparing remote input directory"
Write-Host ("RemoteDataDir: {0}" -f $remoteDataDir)
Write-Host ("=" * 80)

& ssh @sshArgs "mkdir -p $remoteDataDir"
if ($LASTEXITCODE -ne 0) {
    throw "Failed to create remote directory: $remoteDataDir"
}

foreach ($path in $files) {
    Write-Host ("Uploading: {0}" -f $path)
    & scp @scpBaseArgs $path ("{0}@{1}:{2}/" -f $SshUser, $SshHost, $remoteDataDir)
    if ($LASTEXITCODE -ne 0) {
        throw "Upload failed: $path"
    }
}

Write-Host ""
Write-Host ("=" * 80)
Write-Host "Remote directory listing"
Write-Host ("=" * 80)

& ssh @sshArgs "ls -lh $remoteDataDir"
if ($LASTEXITCODE -ne 0) {
    throw "Failed to list remote directory: $remoteDataDir"
}

Write-Host ""
Write-Host "E10gW1 observable upload finished successfully."
