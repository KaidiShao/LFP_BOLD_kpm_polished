param(
    [string[]] $Datasets = @(
        'e10gb1',
        'e10fV1',
        'e10gh1',
        'e10gw1',
        'f12m01',
        'k13m17',
        'k13m23'
    ),
    [string] $RepoWsl = '/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished',
    [string] $PythonWsl = '/home/kdshao/anaconda3/bin/python',
    [string] $LogDir = 'tmp/pipeline5_band_selectivity_current7_logs'
)

$ErrorActionPreference = 'Stop'

$repoRoot = Split-Path -Parent $PSScriptRoot
$logRoot = Join-Path $repoRoot $LogDir
New-Item -ItemType Directory -Force -Path $logRoot | Out-Null

$summaryRows = @()

foreach ($ds in $Datasets) {
    $stamp = Get-Date -Format 'yyyyMMdd_HHmmss'
    $logFile = Join-Path $logRoot ("{0}_{1}.log" -f $ds, $stamp)
    $start = Get-Date

    "=== START $ds $($start.ToString('s')) ===" | Tee-Object -FilePath $logFile

    $wslCommand = @(
        $PythonWsl,
        'scripts/run_pipeline5_band_selectivity.py',
        '--dataset', $ds,
        '--ribbons', 'sem',
        '2>&1'
    ) -join ' '

    & wsl.exe -d Ubuntu-22.04 --cd $RepoWsl bash -lc $wslCommand |
        Tee-Object -FilePath $logFile -Append

    $exitCode = $LASTEXITCODE
    $finish = Get-Date
    $elapsed = [Math]::Round(($finish - $start).TotalMinutes, 2)
    "=== FINISH $ds $($finish.ToString('s')) exit=$exitCode elapsed_min=$elapsed ===" | Tee-Object -FilePath $logFile -Append

    $summaryRows += [pscustomobject]@{
        dataset = $ds
        exit_code = $exitCode
        started = $start.ToString('s')
        finished = $finish.ToString('s')
        elapsed_min = $elapsed
        log_file = $logFile
    }
}

$summaryFile = Join-Path $logRoot ('current7_summary_{0}.csv' -f (Get-Date -Format 'yyyyMMdd_HHmmss'))
$summaryRows | Export-Csv -NoTypeInformation -Path $summaryFile
"Summary: $summaryFile"
