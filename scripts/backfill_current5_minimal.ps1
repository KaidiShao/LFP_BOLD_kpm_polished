param(
    [string[]]$Datasets = @("e10gb1", "e10fV1", "e10gh1", "f12m01", "e10gw1"),
    [string[]]$P3Modes = @("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean"),
    [string[]]$BoldModes = @("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean"),
    [int[]]$ConsistencyComponentCounts = @(3, 4, 5, 6, 7, 8),
    [string]$MatlabExe = "matlab",
    [string]$PythonExe = "python",
    [switch]$Execute,
    [switch]$ParallelGroups,
    [switch]$IndependentGroups,
    [switch]$SkipPreBold,
    [switch]$SkipBlp,
    [switch]$SkipP5Derived,
    [switch]$SkipPostBold,
    [switch]$ContinueOnError
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$matlabBackfillScript = "scripts/script_backfill_current5_minimal_pipeline_set.m"
$dryRun = -not [bool]$Execute

$datasetMap = @{
    "e10gb1" = @{ Cfg = "E10gb1"; Stem = "e10gb1"; DatasetId = "E10.gb1" }
    "e10fV1" = @{ Cfg = "E10fV1"; Stem = "e10fV1"; DatasetId = "E10.fV1" }
    "e10gh1" = @{ Cfg = "E10gH1"; Stem = "e10gh1"; DatasetId = "E10.gH1" }
    "f12m01" = @{ Cfg = "F12m01"; Stem = "f12m01"; DatasetId = "F12.m01" }
    "e10gw1" = @{ Cfg = "E10gW1"; Stem = "e10gw1"; DatasetId = "E10.gW1" }
}

function ConvertTo-MatlabCellLiteral {
    param([string[]]$Values)
    $quoted = foreach ($value in $Values) {
        "'" + ($value -replace "'", "''") + "'"
    }
    return "{" + ($quoted -join ", ") + "}"
}

function ConvertTo-MatlabDatasetStructLiteral {
    param([string[]]$DatasetNames)

    $cfgNames = @()
    $stems = @()
    $datasetIds = @()
    foreach ($name in $DatasetNames) {
        if (-not $datasetMap.ContainsKey($name)) {
            throw "Unknown dataset '$name'. Known datasets: $($datasetMap.Keys -join ', ')"
        }
        $cfgNames += $datasetMap[$name].Cfg
        $stems += $datasetMap[$name].Stem
        $datasetIds += $datasetMap[$name].DatasetId
    }

    $cfgCell = ConvertTo-MatlabCellLiteral $cfgNames
    $stemCell = ConvertTo-MatlabCellLiteral $stems
    $idCell = ConvertTo-MatlabCellLiteral $datasetIds
    return "struct('cfg_name', $cfgCell, 'stem', $stemCell, 'dataset_id', $idCell)"
}

function Invoke-MatlabBackfillPhase {
    param([string]$Phase)

    $cmd = Get-MatlabBackfillCommand $Phase
    Write-Host ""
    Write-Host ("=" * 80)
    Write-Host "MATLAB backfill phase: $Phase"
    Write-Host "dry_run: $dryRun"
    Write-Host ("=" * 80)
    & $MatlabExe -batch $cmd
    if ($LASTEXITCODE -ne 0) {
        throw "MATLAB backfill phase '$Phase' failed with exit code $LASTEXITCODE"
    }
}

function Invoke-ExternalChecked {
    param(
        [string]$Label,
        [string]$FilePath,
        [string[]]$Arguments
    )

    Write-Host ""
    Write-Host $Label
    Write-Host "  $FilePath $($Arguments -join ' ')"
    & $FilePath @Arguments
    if ($LASTEXITCODE -ne 0) {
        $message = "$Label failed with exit code $LASTEXITCODE"
        if ($ContinueOnError) {
            Write-Warning $message
        } else {
            throw $message
        }
    }
}

function Get-MatlabBackfillCommand {
    param([string]$Phase)

    $datasetLiteral = ConvertTo-MatlabDatasetStructLiteral $Datasets
    $p3ModeLiteral = ConvertTo-MatlabCellLiteral $P3Modes
    $modeLiteral = ConvertTo-MatlabCellLiteral $BoldModes
    $dry = if ($dryRun) { "true" } else { "false" }
    $continue = if ($ContinueOnError) { "true" } else { "false" }
    $escapedRepo = $repoRoot -replace "'", "''"
    return "cd('$escapedRepo'); dry_run=$dry; backfill_phase='$Phase'; continue_on_error=$continue; p3_bold_modes=$p3ModeLiteral; minimal_bold_modes=$modeLiteral; backfill_datasets=$datasetLiteral; run('$matlabBackfillScript');"
}

function ConvertTo-PowerShellStringLiteral {
    param([string]$Value)
    return "'" + ($Value -replace "'", "''") + "'"
}

function Start-LoggedExternalJob {
    param(
        [string]$Name,
        [string]$FilePath,
        [string[]]$Arguments,
        [string]$LogPath
    )

    $logDir = Split-Path -Parent $LogPath
    New-Item -ItemType Directory -Force -Path $logDir | Out-Null
    $wrapperPath = Join-Path $logDir "$Name.run.ps1"
    $argLiteral = "@(" + (($Arguments | ForEach-Object { ConvertTo-PowerShellStringLiteral $_ }) -join ", ") + ")"
    $scriptLines = @(
        '$ErrorActionPreference = "Continue"',
        '$exe = ' + (ConvertTo-PowerShellStringLiteral $FilePath),
        '$argList = ' + $argLiteral,
        '$log = ' + (ConvertTo-PowerShellStringLiteral $LogPath),
        'New-Item -ItemType Directory -Force -Path (Split-Path -Parent $log) | Out-Null',
        '"[$(Get-Date -Format ''yyyy-MM-dd HH:mm:ss'')] START $exe $($argList -join '' '')" | Tee-Object -FilePath $log',
        '& $exe @argList 2>&1 | Tee-Object -FilePath $log -Append',
        '$code = $LASTEXITCODE',
        '"[$(Get-Date -Format ''yyyy-MM-dd HH:mm:ss'')] EXIT $code" | Tee-Object -FilePath $log -Append',
        'exit $code'
    )
    Set-Content -Path $wrapperPath -Value $scriptLines -Encoding ASCII

    $psExe = Join-Path $PSHOME "powershell.exe"
    $process = Start-Process -FilePath $psExe `
        -ArgumentList @("-NoProfile", "-ExecutionPolicy", "Bypass", "-File", $wrapperPath) `
        -PassThru -WindowStyle Hidden
    Write-Host "Launched $Name (pid=$($process.Id)); log: $LogPath"
    return [pscustomobject]@{ Name = $Name; Process = $process; Log = $LogPath; Wrapper = $wrapperPath }
}

function Wait-BackfillJobs {
    param([System.Collections.ArrayList]$Jobs)

    $failed = 0
    while ($Jobs.Count -gt 0) {
        foreach ($job in @($Jobs.ToArray())) {
            $job.Process.Refresh()
            if (-not $job.Process.HasExited) {
                continue
            }

            $code = $job.Process.ExitCode
            if ($code -eq 0) {
                Write-Host "Finished $($job.Name) (exit 0); log: $($job.Log)"
            } else {
                $failed += 1
                Write-Warning "Group job failed: $($job.Name) (exit $code); log: $($job.Log)"
                if (Test-Path $job.Log) {
                    Write-Host "Last log lines for $($job.Name):"
                    Get-Content -Path $job.Log -Tail 30
                }
            }

            [void]$Jobs.Remove($job)

            if (($code -ne 0) -and (-not $ContinueOnError)) {
                throw "Group job failed: $($job.Name)"
            }
        }

        if ($Jobs.Count -gt 0) {
            Start-Sleep -Seconds 5
        }
    }

    if ($failed -gt 0) {
        Write-Warning "$failed group job(s) failed. Check logs above."
    }
}

function Invoke-ParallelBackfillGroups {
    $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
    $logRoot = Join-Path $repoRoot "tmp/backfill_current5/$stamp"
    New-Item -ItemType Directory -Force -Path $logRoot | Out-Null
    $launchedAny = $false

    Write-Host ""
    Write-Host "Parallel group mode enabled."
    Write-Host "Each group runs as one worker. Dependency gates are still respected:"
    Write-Host "  wave 1: P3 and P5+P6 can run together"
    Write-Host "  wave 2: P7+P8 runs after upstream MATLAB workers finish"
    Write-Host "  note: P4 BOLD training is manual and is never launched by this driver"
    Write-Host "Logs: $logRoot"

    $upstreamJobs = [System.Collections.ArrayList]::new()

    if (-not $SkipPreBold) {
        $cmd = Get-MatlabBackfillCommand "pre_bold"
        $log = Join-Path $logRoot "P3_pre_bold.log"
        [void]$upstreamJobs.Add((Start-LoggedExternalJob -Name "P3_pre_bold" -FilePath $MatlabExe -Arguments @("-batch", $cmd) -LogPath $log))
        $launchedAny = $true
    }
    if (-not $SkipBlp) {
        $cmd = Get-MatlabBackfillCommand "blp"
        $log = Join-Path $logRoot "P5_P6_blp.log"
        [void]$upstreamJobs.Add((Start-LoggedExternalJob -Name "P5_P6_blp" -FilePath $MatlabExe -Arguments @("-batch", $cmd) -LogPath $log))
        $launchedAny = $true
    }

    if ($upstreamJobs.Count -gt 0) {
        Wait-BackfillJobs $upstreamJobs
    }

    if (-not $SkipPostBold) {
        $cmd = Get-MatlabBackfillCommand "post_bold"
        $log = Join-Path $logRoot "P7_P8_post_bold.log"
        $postJobs = [System.Collections.ArrayList]::new()
        [void]$postJobs.Add((Start-LoggedExternalJob -Name "P7_P8_post_bold" -FilePath $MatlabExe -Arguments @("-batch", $cmd) -LogPath $log))
        $launchedAny = $true
        Wait-BackfillJobs $postJobs
    }

    if (-not $launchedAny) {
        Write-Host "No groups were selected."
        return
    }
}

function Invoke-IndependentBackfillGroups {
    $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
    $logRoot = Join-Path $repoRoot "tmp/backfill_current5/$stamp"
    New-Item -ItemType Directory -Force -Path $logRoot | Out-Null
    $jobs = [System.Collections.ArrayList]::new()

    Write-Host ""
    Write-Host "Independent group mode enabled."
    Write-Host "Each group runs as one worker and independently detects what can be filled now:"
    Write-Host "  P3 worker: BOLD observables/QC"
    Write-Host "  P5+P6 worker: BLP reduction/postprocessing"
    Write-Host "  P7+P8 worker: BOLD postprocessing/coupling for currently available upstream data"
    Write-Host "  note: P4 BOLD training is manual and is never launched by this driver"
    Write-Host "Logs: $logRoot"

    if (-not $SkipPreBold) {
        $cmd = Get-MatlabBackfillCommand "pre_bold"
        $log = Join-Path $logRoot "P3_pre_bold.log"
        [void]$jobs.Add((Start-LoggedExternalJob -Name "P3_pre_bold" -FilePath $MatlabExe -Arguments @("-batch", $cmd) -LogPath $log))
    }
    if (-not $SkipBlp) {
        $cmd = Get-MatlabBackfillCommand "blp"
        $log = Join-Path $logRoot "P5_P6_blp.log"
        [void]$jobs.Add((Start-LoggedExternalJob -Name "P5_P6_blp" -FilePath $MatlabExe -Arguments @("-batch", $cmd) -LogPath $log))
    }
    if (-not $SkipPostBold) {
        $cmd = Get-MatlabBackfillCommand "post_bold"
        $log = Join-Path $logRoot "P7_P8_post_bold.log"
        [void]$jobs.Add((Start-LoggedExternalJob -Name "P7_P8_post_bold" -FilePath $MatlabExe -Arguments @("-batch", $cmd) -LogPath $log))
    }

    if ($jobs.Count -eq 0) {
        Write-Host "No groups were selected."
        return
    }

    Wait-BackfillJobs $jobs
}

function Invoke-P5DerivedReports {
    Write-Host ""
    Write-Host ("=" * 80)
    Write-Host "P5 derived inspection outputs"
    Write-Host "dry_run: $dryRun"
    Write-Host ("=" * 80)

    if ($dryRun) {
        Write-Host "Would generate:"
        Write-Host "  peak-state statistics with figures"
        Write-Host "  flat peak-statistic summary figures"
        Write-Host "  per-dataset method consistency figures for k=$($ConsistencyComponentCounts -join ', ')"
        Write-Host "  cross-dataset consistency figures for k=$($ConsistencyComponentCounts -join ', ')"
        return
    }

    $datasetLiteral = ConvertTo-MatlabCellLiteral $Datasets
    $continue = if ($ContinueOnError) { "true" } else { "false" }
    $escapedRepo = $repoRoot -replace "'", "''"

    $peakCmd = "cd('$escapedRepo'); dataset_stems=$datasetLiteral; peak_mode='max_abs'; save_figures=true; continue_on_error=$continue; dry_run=false; run('scripts/script_run_current_pipeline5_peak_state_all_datasets.m');"
    Invoke-ExternalChecked -Label "P5 peak-state statistics + figures" -FilePath $MatlabExe -Arguments @("-batch", $peakCmd)

    $summaryCmd = "cd('$escapedRepo'); dataset_stems=$datasetLiteral; peak_stage='pipeline5_eigenfunction_peaks_by_state_maxabs'; run('scripts/plot_pipeline5_peak_statistics_summary_figures.m');"
    Invoke-ExternalChecked -Label "P5 flat peak-statistic summary figures" -FilePath $MatlabExe -Arguments @("-batch", $summaryCmd)

    Push-Location $repoRoot
    try {
        $summarizeScript = Join-Path $repoRoot "scripts/summarize_peak_state_method_consistency.py"
        $reportScript = Join-Path $repoRoot "scripts/report_peak_state_method_consistency_all_datasets.py"
        $plotPerDatasetScript = Join-Path $repoRoot "scripts/plot_peak_state_method_consistency_per_dataset.py"

        foreach ($k in $ConsistencyComponentCounts) {
            $kk = "{0:D2}" -f $k

            foreach ($ds in $Datasets) {
                $root = "E:\DataPons_processed\$ds\pipeline5_eigenfunction_peaks_by_state_maxabs"
                $outDir = "results\peak_state_method_consistency_current_maxabs_k$kk`_$ds"
                Invoke-ExternalChecked `
                    -Label "P5 method consistency summary k=$kk dataset=$ds" `
                    -FilePath $PythonExe `
                    -Arguments @($summarizeScript, "--root", $root, "--output-dir", $outDir, "--component-counts", "$k")
            }

            Invoke-ExternalChecked `
                -Label "P5 cross-dataset consistency report k=$kk" `
                -FilePath $PythonExe `
                -Arguments @(
                    $reportScript,
                    "--input-prefix", "peak_state_method_consistency_current_maxabs_k$kk`_",
                    "--output-dirname", "peak_state_method_consistency_current_all_datasets_maxabs_k$kk"
                )
        }

        Invoke-ExternalChecked `
            -Label "P5 per-dataset consistency figures" `
            -FilePath $PythonExe `
            -Arguments @($plotPerDatasetScript, "--input-glob", "peak_state_method_consistency_current_maxabs_k*_*")
    } finally {
        Pop-Location
    }
}

Write-Host "Current-five minimal backfill driver"
Write-Host "Repo: $repoRoot"
Write-Host "Datasets: $($Datasets -join ', ')"
Write-Host "P3 modes: $($P3Modes -join ', ')"
Write-Host "BOLD modes: $($BoldModes -join ', ')"
Write-Host "P5 consistency k: $($ConsistencyComponentCounts -join ', ')"
Write-Host "Execute: $([bool]$Execute)"
Write-Host "P4 BOLD training: manual only; not launched by this driver"
Write-Host "ParallelGroups: $([bool]$ParallelGroups)"
Write-Host "IndependentGroups: $([bool]$IndependentGroups)"
Write-Host "P5 derived reports: $(-not [bool]$SkipP5Derived)"

if ($ParallelGroups -and $IndependentGroups) {
    throw "Use either -ParallelGroups or -IndependentGroups, not both."
}

if ($IndependentGroups) {
    Invoke-IndependentBackfillGroups
    if (-not $SkipP5Derived -and -not $SkipBlp) {
        Invoke-P5DerivedReports
    }
    Write-Host ""
    if ($dryRun) {
        Write-Host "Independent dry-run finished. Pass -Execute to launch the missing work."
    } else {
        Write-Host "Independent backfill groups finished."
    }
    return
}

if ($ParallelGroups) {
    Invoke-ParallelBackfillGroups
    if (-not $SkipP5Derived -and -not $SkipBlp) {
        Invoke-P5DerivedReports
    }
    Write-Host ""
    if ($dryRun) {
        Write-Host "Parallel dry-run finished. Pass -Execute to launch the missing work."
    } else {
        Write-Host "Parallel backfill groups finished."
    }
    return
}

if (-not $SkipPreBold) {
    Invoke-MatlabBackfillPhase "pre_bold"
}
if (-not $SkipBlp) {
    Invoke-MatlabBackfillPhase "blp"
}
if (-not $SkipP5Derived -and -not $SkipBlp) {
    Invoke-P5DerivedReports
}
Write-Host ""
Write-Host "Skipping P4 BOLD training by design. Run the WSL launcher manually when you want to inspect training."
if (-not $SkipPostBold) {
    Invoke-MatlabBackfillPhase "post_bold"
}

Write-Host ""
if ($dryRun) {
    Write-Host "Dry-run finished. Pass -Execute to launch the missing work."
} else {
    Write-Host "Backfill driver finished."
}
