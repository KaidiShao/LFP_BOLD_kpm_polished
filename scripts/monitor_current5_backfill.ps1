param(
    [string[]]$Datasets = @("e10gb1", "e10fV1", "e10gh1", "f12m01", "e10gw1"),
    [string[]]$P3Modes = @("global_svd100", "gsvd100_ds", "HP", "HP_svd100", "roi_mean"),
    [string[]]$BoldModes = @("global_svd100", "gsvd100_ds", "HP", "HP_svd100", "roi_mean"),
    [string]$ProcessedRoot = "E:\DataPons_processed",
    [string]$BoldResultRoot = "E:\autodl_results_local\bold_wsl",
    [string]$BackfillLogRoot = "",
    [int]$RefreshSeconds = 30,
    [int]$TailLines = 8,
    [switch]$Once,
    [switch]$ShowLogs,
    [switch]$ShowProcesses,
    [switch]$NoClear
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
if ([string]::IsNullOrWhiteSpace($BackfillLogRoot)) {
    $BackfillLogRoot = Join-Path $repoRoot "tmp\backfill_current5"
}

$datasetMap = @{
    "e10gb1" = @{ Label = "E10.gb1"; Stem = "e10gb1"; DatasetId = "E10.gb1" }
    "e10fV1" = @{ Label = "E10.fV1"; Stem = "e10fV1"; DatasetId = "E10.fV1" }
    "e10gh1" = @{ Label = "E10.gH1"; Stem = "e10gh1"; DatasetId = "E10.gH1" }
    "f12m01" = @{ Label = "F12.m01"; Stem = "f12m01"; DatasetId = "F12.m01" }
    "e10gw1" = @{ Label = "E10.gW1"; Stem = "e10gw1"; DatasetId = "E10.gW1" }
}

function Count-Files {
    param([string]$Path, [string]$Filter = "*")
    if (-not (Test-Path -LiteralPath $Path)) { return 0 }
    return @(Get-ChildItem -LiteralPath $Path -Recurse -File -Filter $Filter -ErrorAction SilentlyContinue).Count
}

function Format-Count {
    param([int]$Done, [int]$Total)
    if ($Total -le 0) { return "$Done" }
    if ($Done -ge $Total) { return "$Done/$Total ok" }
    if ($Done -gt 0) { return "$Done/$Total part" }
    return "0/$Total miss"
}

function Get-P8Tag {
    param([string]$Mode)
    switch ($Mode.ToLowerInvariant()) {
        "svd" { return "pv_svd" }
        "global_svd100" { return "pv_gsvd100" }
        "hp_svd100" { return "pv_hp100" }
        "hp" { return "pv_hp" }
        "elehp" { return "pv_elehp" }
        "roi_mean_slow_band_power" { return "pv_roi_sbp" }
        "slow_band_power_svd" { return "pv_sbp_svd" }
        "gsvd100_ds" { return "pv_gsvd100_ds" }
        default { return "pv_" + (($Mode.ToLowerInvariant()) -replace "[^a-z0-9]+", "_") }
    }
}

function Has-P4BoldMode {
    param([string]$Stem, [string]$Mode)
    $root = Join-Path $BoldResultRoot "$Stem\mlp\outputs"
    if (-not (Test-Path -LiteralPath $root)) { return $false }
    $dirs = Get-ChildItem -LiteralPath $root -Directory -ErrorAction SilentlyContinue |
        Where-Object { Test-RunNameMode -RunName $_.Name -Mode $Mode }
    foreach ($dir in $dirs) {
        if (@(Get-ChildItem -LiteralPath $dir.FullName -File -Filter "*summary.mat" -ErrorAction SilentlyContinue).Count -gt 0) {
            return $true
        }
    }
    return $false
}

function Has-P7Mode {
    param([string]$Stem, [string]$Mode)
    $root = Join-Path $ProcessedRoot "$Stem\pipeline7_bold_reskoopnet_postprocessing"
    if (-not (Test-Path -LiteralPath $root)) { return $false }
    $dirs = Get-ChildItem -LiteralPath $root -Directory -ErrorAction SilentlyContinue |
        Where-Object { Test-RunNameMode -RunName $_.Name -Mode $Mode }
    foreach ($dir in $dirs) {
        if ((Count-Files $dir.FullName "*.mat") -gt 0 -and (Count-Files $dir.FullName "*.png") -gt 0) {
            return $true
        }
    }
    return $false
}

function Test-RunNameMode {
    param([string]$RunName, [string]$Mode)
    $suffix = "projected_vlambda_$Mode"
    return $RunName.EndsWith($suffix, [System.StringComparison]::Ordinal)
}

function Has-P8Mode {
    param([string]$Stem, [string]$Mode)
    $tag = Get-P8Tag $Mode
    $xcorr = Join-Path $ProcessedRoot "$Stem\pipeline8_xcorr\$tag"
    $maps = Join-Path $ProcessedRoot "$Stem\pipeline8_top_maps\$tag"
    $hasXcorr = (Count-Files $xcorr "*.csv") -gt 0 -and (Count-Files $xcorr "*.mat") -gt 0
    $hasMaps = (Count-Files $maps "*.png") -gt 0 -and (Count-Files $maps "*.mat") -gt 0
    return ($hasXcorr -and $hasMaps)
}

function Get-LatestWrite {
    param([string[]]$Paths)
    $latest = $null
    foreach ($path in $Paths) {
        if (-not (Test-Path -LiteralPath $path)) { continue }
        $item = Get-Item -LiteralPath $path -ErrorAction SilentlyContinue
        if ($item -and (($null -eq $latest) -or ($item.LastWriteTime -gt $latest))) {
            $latest = $item.LastWriteTime
        }
        $file = Get-ChildItem -LiteralPath $path -Recurse -File -ErrorAction SilentlyContinue |
            Sort-Object LastWriteTime -Descending |
            Select-Object -First 1
        if ($file -and (($null -eq $latest) -or ($file.LastWriteTime -gt $latest))) {
            $latest = $file.LastWriteTime
        }
    }
    if ($latest) { return $latest.ToString("MM-dd HH:mm") }
    return ""
}

function Get-StatusRows {
    $rows = foreach ($name in $Datasets) {
        if (-not $datasetMap.ContainsKey($name)) {
            Write-Warning "Unknown dataset skipped: $name"
            continue
        }
        $d = $datasetMap[$name]
        $stem = $d.Stem
        $datasetId = $d.DatasetId
        $root = Join-Path $ProcessedRoot $stem

        $p3ObsDone = 0
        $p3QcDone = 0
        foreach ($mode in $P3Modes) {
            $obs = Join-Path $root "pipeline3_bold_observables\$($datasetId)_bold_observables_$mode.mat"
            if (Test-Path -LiteralPath $obs) { $p3ObsDone++ }
            $qc = Join-Path $root "pipeline3_figures_bold_pre_reskoopnet_qc\$mode"
            if ((Count-Files $qc "*.png") -gt 0) { $p3QcDone++ }
        }

        $p4Done = 0
        $p7Done = 0
        $p8Done = 0
        foreach ($mode in $BoldModes) {
            if (Has-P4BoldMode $stem $mode) { $p4Done++ }
            if (Has-P7Mode $stem $mode) { $p7Done++ }
            if (Has-P8Mode $stem $mode) { $p8Done++ }
        }

        $p5Consensus = 0
        foreach ($method in @("svd", "mds", "umap")) {
            $consensus = Join-Path $root "pipeline5_summary_figures\abs_projected_vlambda\eigenfunction_reduction\consensus_state_space"
            if (@(Get-ChildItem -LiteralPath $consensus -File -Filter "$($method)_k08__efun_ssc__*.png" -ErrorAction SilentlyContinue).Count -gt 0) {
                $p5Consensus++
            }
        }
        $p5Density = Count-Files (Join-Path $root "pipeline5_dimred_thresholded_density\abs_projected_vlambda\umap_k08\mat\q070") "*.mat"
        $p5 = "$(Format-Count $p5Consensus 3); dens=$p5Density"

        $p6Top = Count-Files (Join-Path $root "pipeline6_top_state_diversity_postprocessing") "*.png"
        $p6Time = Count-Files (Join-Path $root "pipeline6_figures_timescale_diagnostics") "*.png"
        $p6Cross = (Count-Files (Join-Path $root "pipeline6_spkt_residual_cross_correlation") "*.mat") +
            (Count-Files (Join-Path $root "pipeline6_mua_residual_cross_correlation") "*.mat")
        $p6 = if ($p6Top -gt 0 -and $p6Time -gt 0 -and $p6Cross -gt 0) {
            "ok"
        } elseif ($p6Top -gt 0 -or $p6Time -gt 0 -or $p6Cross -gt 0) {
            "part"
        } else {
            "miss"
        }

        [pscustomobject]@{
            Dataset = $d.Label
            P3_Obs = Format-Count $p3ObsDone $P3Modes.Count
            P3_QC = Format-Count $p3QcDone $P3Modes.Count
            P4_BOLD = Format-Count $p4Done $BoldModes.Count
            P5 = $p5
            P6 = $p6
            P7 = Format-Count $p7Done $BoldModes.Count
            P8 = Format-Count $p8Done $BoldModes.Count
            Latest = Get-LatestWrite @($root, (Join-Path $BoldResultRoot $stem))
        }
    }
    return $rows
}

function Show-ActiveProcesses {
    $procs = Get-Process -ErrorAction SilentlyContinue |
        Where-Object { $_.ProcessName -match "MATLAB|matlab|wsl|python|powershell" } |
        Select-Object ProcessName, Id, StartTime, CPU |
        Sort-Object ProcessName, StartTime
    if ($procs) {
        Write-Host ""
        Write-Host "Active related processes:"
        $procs | Format-Table -AutoSize
    }
}

function Show-LatestLogs {
    if (-not (Test-Path -LiteralPath $BackfillLogRoot)) {
        Write-Host ""
        Write-Host "No backfill log root yet: $BackfillLogRoot"
        return
    }
    $latestRun = Get-ChildItem -LiteralPath $BackfillLogRoot -Directory -ErrorAction SilentlyContinue |
        Sort-Object LastWriteTime -Descending |
        Select-Object -First 1
    if (-not $latestRun) {
        Write-Host ""
        Write-Host "No backfill logs yet."
        return
    }

    Write-Host ""
    Write-Host "Latest backfill log dir: $($latestRun.FullName)"
    $logs = Get-ChildItem -LiteralPath $latestRun.FullName -File -Filter "*.log" -ErrorAction SilentlyContinue |
        Sort-Object Name
    foreach ($log in $logs) {
        Write-Host ""
        Write-Host "--- $($log.Name)  [$($log.LastWriteTime.ToString('MM-dd HH:mm:ss'))]"
        Get-Content -LiteralPath $log.FullName -Tail $TailLines -ErrorAction SilentlyContinue
    }
}

do {
    if (-not $NoClear) { Clear-Host }
    Write-Host "Current-five backfill monitor  $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
    Write-Host "P3 modes: $($P3Modes -join ', ')"
    Write-Host "P4/P7/P8 modes: $($BoldModes -join ', ')"
    Write-Host ""
    Get-StatusRows | Format-Table -AutoSize
    if ($ShowProcesses) { Show-ActiveProcesses }
    if ($ShowLogs) { Show-LatestLogs }

    if ($Once) { break }
    Start-Sleep -Seconds $RefreshSeconds
} while ($true)
