param(
    [string[]] $Datasets = @("e10fV1", "e10gw1"),
    [string] $LogPath = "",
    [string] $RunReduction = "true",
    [string] $RunDensity = "true",
    [string] $RunPeakStats = "true"
)

$ErrorActionPreference = "Stop"

function Convert-ToBoolString {
    param([string] $Value)
    $v = $Value.Trim().ToLowerInvariant()
    if ($v -in @("1", "true", "yes", "y")) {
        return "true"
    }
    if ($v -in @("0", "false", "no", "n")) {
        return "false"
    }
    throw "Expected a boolean-like value, got: $Value"
}

$repoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $repoRoot

if ([string]::IsNullOrWhiteSpace($LogPath)) {
    $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
    $LogPath = Join-Path $repoRoot "tmp\logs\p5_standardized_csplit_e10fv1_e10gw1_$stamp.log"
}

$logDir = Split-Path -Parent $LogPath
if (-not (Test-Path -LiteralPath $logDir)) {
    New-Item -ItemType Directory -Path $logDir | Out-Null
}

$datasetLiteral = ($Datasets | ForEach-Object { "'" + ($_ -replace "'", "''") + "'" }) -join ","
$matlabRunReduction = Convert-ToBoolString $RunReduction
$matlabRunDensity = Convert-ToBoolString $RunDensity
$matlabRunPeakStats = Convert-ToBoolString $RunPeakStats
$batch = "cd('$($repoRoot -replace "'", "''")'); target_dataset_stems={$datasetLiteral}; run_reduction=$matlabRunReduction; run_density=$matlabRunDensity; run_peak_stats=$matlabRunPeakStats; run('scripts/script_run_standardized_csplit_p5_for_p12.m');"
$matlab = "C:\Program Files\MATLAB\R2025b\bin\matlab.exe"

"[$(Get-Date -Format s)] Starting MATLAB P5 standardized csplit run" | Tee-Object -FilePath $LogPath
"Datasets: $($Datasets -join ', ')" | Tee-Object -FilePath $LogPath -Append
"RunReduction: $matlabRunReduction | RunDensity: $matlabRunDensity | RunPeakStats: $matlabRunPeakStats" | Tee-Object -FilePath $LogPath -Append
"LogPath : $LogPath" | Tee-Object -FilePath $LogPath -Append

& $matlab -batch $batch *>&1 | Tee-Object -FilePath $LogPath -Append
$exitCode = $LASTEXITCODE

"[$(Get-Date -Format s)] MATLAB exited with code $exitCode" | Tee-Object -FilePath $LogPath -Append
exit $exitCode
