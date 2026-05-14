param(
    [int]$PollSeconds = 60,
    [int]$SettleSeconds = 20,
    [int]$Rounds = 2,
    [switch]$SkipWaitExisting
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$driver = Join-Path $PSScriptRoot "backfill_current5_minimal.ps1"
$monitor = Join-Path $PSScriptRoot "monitor_current5_backfill.ps1"
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$supervisorRoot = Join-Path $repoRoot "tmp/backfill_current5_supervisor/$stamp"
$supervisorLog = Join-Path $supervisorRoot "supervisor.log"

New-Item -ItemType Directory -Force -Path $supervisorRoot | Out-Null

function Write-SupervisorLog {
    param([string]$Message)
    $line = "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] $Message"
    $line | Tee-Object -FilePath $supervisorLog -Append
}

function Get-LatestBackfillDir {
    $root = Join-Path $repoRoot "tmp/backfill_current5"
    if (-not (Test-Path $root)) {
        return $null
    }
    return Get-ChildItem -Path $root -Directory |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
}

function Test-BackfillDirFinished {
    param([System.IO.DirectoryInfo]$Dir)

    if ($null -eq $Dir -or -not (Test-Path $Dir.FullName)) {
        return $true
    }

    $logs = @(Get-ChildItem -Path $Dir.FullName -Filter "*.log" -File)
    if ($logs.Count -eq 0) {
        return $true
    }

    foreach ($log in $logs) {
        if (-not (Select-String -Path $log.FullName -Pattern "\] EXIT \d+" -Quiet)) {
            return $false
        }
    }
    return $true
}

function Invoke-LoggedPowerShell {
    param([string[]]$Arguments)

    $psExe = Join-Path $PSHOME "powershell.exe"
    Write-SupervisorLog "RUN $psExe $($Arguments -join ' ')"
    & $psExe @Arguments 2>&1 | Tee-Object -FilePath $supervisorLog -Append
    $code = $LASTEXITCODE
    Write-SupervisorLog "EXIT $code"
    return $code
}

Write-SupervisorLog "Current-five unattended restart supervisor started."
Write-SupervisorLog "Repo: $repoRoot"

if (-not $SkipWaitExisting) {
    $latest = Get-LatestBackfillDir
    if ($null -ne $latest) {
        Write-SupervisorLog "Waiting for existing backfill dir to finish: $($latest.FullName)"
        while (-not (Test-BackfillDirFinished $latest)) {
            $openLogs = Get-ChildItem -Path $latest.FullName -Filter "*.log" -File |
                Where-Object { -not (Select-String -Path $_.FullName -Pattern "\] EXIT \d+" -Quiet) } |
                Select-Object -ExpandProperty Name
            Write-SupervisorLog "Still running/waiting: $($openLogs -join ', ')"
            Start-Sleep -Seconds $PollSeconds
            $latest.Refresh()
        }
        Write-SupervisorLog "Existing backfill dir finished."
        if ($SettleSeconds -gt 0) {
            Write-SupervisorLog "Settling for $SettleSeconds seconds before restart."
            Start-Sleep -Seconds $SettleSeconds
        }
    }
}

$driverCode = 0
for ($round = 1; $round -le $Rounds; $round++) {
    Write-SupervisorLog "Starting independent backfill round $round of $Rounds."
    $driverArgs = @(
        "-NoProfile",
        "-ExecutionPolicy", "Bypass",
        "-File", $driver,
        "-IndependentGroups",
        "-Execute",
        "-ContinueOnError"
    )
    $driverCode = Invoke-LoggedPowerShell $driverArgs
    if ($driverCode -ne 0) {
        Write-SupervisorLog "Driver round $round failed with exit code $driverCode; stopping remaining rounds."
        break
    }
    if (($round -lt $Rounds) -and ($SettleSeconds -gt 0)) {
        Write-SupervisorLog "Settling for $SettleSeconds seconds before next independent round."
        Start-Sleep -Seconds $SettleSeconds
    }
}

Write-SupervisorLog "Final monitor snapshot follows."
$monitorArgs = @(
    "-NoProfile",
    "-ExecutionPolicy", "Bypass",
    "-File", $monitor,
    "-Once",
    "-NoClear",
    "-ShowLogs"
)
[void](Invoke-LoggedPowerShell $monitorArgs)

Write-SupervisorLog "Supervisor finished. Driver exit code: $driverCode"
exit $driverCode
