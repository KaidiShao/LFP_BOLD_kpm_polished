param(
    [int]$WaitForPid = 0,
    [int]$PollSeconds = 60,
    [string]$Distro = "Ubuntu-22.04"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repo = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished"
$logDir = Join-Path $repo "tmp\run_logs"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null
$logPath = Join-Path $logDir "bold_stdT_k13m23_global_gsvd_after_wait_20260601_launcher.log"

function Write-RunLog {
    param([string]$Message)
    $line = "[{0}] {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $Message
    $line | Tee-Object -FilePath $logPath -Append
}

Write-RunLog "K13m23 BOLD global/gsvd P4 launcher started."
Write-RunLog "Repo: $repo"
Write-RunLog "WaitForPid: $WaitForPid"

if ($WaitForPid -gt 0) {
    while ($true) {
        & wsl.exe -d $Distro bash -lc "ps -p $WaitForPid -o pid= >/dev/null 2>&1"
        if ($LASTEXITCODE -ne 0) {
            Write-RunLog "Wait WSL PID $WaitForPid finished; launching BOLD P4 queue."
            break
        }
        $desc = & wsl.exe -d $Distro bash -lc "ps -p $WaitForPid -o pid=,etime=,pcpu=,pmem=,cmd= | head -n 1"
        Write-RunLog ("Waiting for WSL PID {0}: {1}" -f $WaitForPid, $desc)
        Start-Sleep -Seconds $PollSeconds
    }
}

$wslRepo = "/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished"
$cmd = "cd '$wslRepo' && /home/kdshao/anaconda3/bin/python -u scripts/run_k13m23_bold_global_gsvd_p4_queue.py"
Write-RunLog "WSL command: $cmd"

& wsl.exe -d $Distro bash -lc $cmd 2>&1 | Tee-Object -FilePath $logPath -Append
$code = $LASTEXITCODE
Write-RunLog "K13m23 BOLD global/gsvd P4 launcher finished with exit code $code."
exit $code
