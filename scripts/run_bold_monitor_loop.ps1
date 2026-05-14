param(
    [int]$IntervalSeconds = 300,
    [string]$PythonExe = "python"
)

$ErrorActionPreference = "Continue"
$RepoRoot = Resolve-Path (Join-Path $PSScriptRoot "..")
$ReportDir = Join-Path $RepoRoot "tmp\monitor_reports"
$LogPath = Join-Path $ReportDir "bold_global_svd_monitor_loop.log"
$PidPath = Join-Path $ReportDir "bold_global_svd_monitor_loop.pid"

New-Item -ItemType Directory -Force -Path $ReportDir | Out-Null
Set-Content -Path $PidPath -Value $PID -Encoding UTF8

while ($true) {
    $stamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss K"
    Add-Content -Path $LogPath -Value "[$stamp] monitor tick pid=$PID" -Encoding UTF8
    try {
        Push-Location $RepoRoot
        & $PythonExe "scripts\monitor_bold_global_svd_runs.py" *>&1 |
            Add-Content -Path $LogPath -Encoding UTF8
        Pop-Location
    } catch {
        Add-Content -Path $LogPath -Value "[$stamp] ERROR: $_" -Encoding UTF8
        try { Pop-Location } catch {}
    }
    Start-Sleep -Seconds $IntervalSeconds
}
