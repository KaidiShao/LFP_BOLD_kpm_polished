param(
    [string]$WaitForPids = "",
    [int]$PollSeconds = 60,
    [string]$MatlabExe = "C:\Program Files\MATLAB\R2023a\bin\matlab.exe",
    [string]$Datasets = "e10gb1,e10fV1,e10gh1,f12m01,e10gw1",
    [int]$MaxLagBins = 80,
    [switch]$ForceRecompute,
    [switch]$ContinueOnError,
    [switch]$DryRun
)

$ErrorActionPreference = "Stop"

$repoPath = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$logDir = Join-Path $repoPath "tmp\p6_spkt_lagged\$stamp"
$logFile = Join-Path $logDir "lagged.log"
$latestFile = Join-Path $repoPath "tmp\p6_spkt_lagged\latest_logdir.txt"

New-Item -ItemType Directory -Force -Path $logDir | Out-Null
Set-Content -Path $latestFile -Value $logDir

function Write-Log {
    param([string]$Message)
    $line = "[{0}] {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $Message
    $line | Tee-Object -FilePath $logFile -Append
}

function Convert-ToMatlabCell {
    param([string[]]$Values)
    $quoted = @()
    foreach ($value in $Values) {
        $quoted += "'" + ($value -replace "'", "''") + "'"
    }
    return "{" + ($quoted -join ", ") + "}"
}

$pidList = @()
if (-not [string]::IsNullOrWhiteSpace($WaitForPids)) {
    $pidList = $WaitForPids -split "[,\s]+" |
        Where-Object { -not [string]::IsNullOrWhiteSpace($_) } |
        ForEach-Object { [int]$_ }
}

$datasetList = $Datasets -split "[,\s]+" |
    Where-Object { -not [string]::IsNullOrWhiteSpace($_) }

Write-Log "Current-five SPKT residual lagged runner"
Write-Log "Repo: $repoPath"
Write-Log "Datasets: $($datasetList -join ', ')"
Write-Log "Max lag bins: $MaxLagBins"
Write-Log "MATLAB: $MatlabExe"
Write-Log "DryRun: $([bool]$DryRun)"
Write-Log "ForceRecompute: $([bool]$ForceRecompute)"
Write-Log "ContinueOnError: $([bool]$ContinueOnError)"
Write-Log "Log dir: $logDir"

if ($pidList.Count -gt 0) {
    Write-Log "Waiting for PID(s): $($pidList -join ', ')"
    while ($true) {
        $alive = @()
        foreach ($pidValue in $pidList) {
            $proc = Get-Process -Id $pidValue -ErrorAction SilentlyContinue
            if ($proc) {
                $alive += $proc
            }
        }

        if ($alive.Count -eq 0) {
            Write-Log "All waited PID(s) have exited."
            break
        }

        $aliveDesc = $alive | ForEach-Object {
            "{0}:{1}" -f $_.Id, $_.ProcessName
        }
        Write-Log "Still waiting: $($aliveDesc -join ', ')"
        Start-Sleep -Seconds $PollSeconds
    }
}

$repoForMatlab = $repoPath -replace "'", "''"
$datasetsForMatlab = Convert-ToMatlabCell -Values $datasetList
$forceForMatlab = if ($ForceRecompute) { "true" } else { "false" }
$continueForMatlab = if ($ContinueOnError) { "true" } else { "false" }
$dryForMatlab = if ($DryRun) { "true" } else { "false" }

$batch = @(
    "cd('$repoForMatlab')",
    "dataset_stems=$datasetsForMatlab",
    "max_lag_bins=$MaxLagBins",
    "lag_features={'abs_rms', 'abs_mean'}",
    "force_recompute=$forceForMatlab",
    "continue_on_error=$continueForMatlab",
    "dry_run=$dryForMatlab",
    "run('scripts/script_run_current5_spkt_residual_lagged.m')"
) -join "; "

Write-Log "Starting MATLAB lagged run."
Write-Log "Command: $MatlabExe -batch ""$batch"""

& $MatlabExe -batch $batch 2>&1 | Tee-Object -FilePath $logFile -Append
$exitCode = $LASTEXITCODE

Write-Log "MATLAB exit code: $exitCode"
exit $exitCode
