param(
    [switch]$Execute,
    [switch]$NoPeakStats,
    [string]$MatlabExe = ""
)

$ErrorActionPreference = "Stop"

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$repoRoot = Split-Path -Parent $scriptDir
$repoMatlab = $repoRoot.Replace("'", "''")
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$logDir = Join-Path $repoRoot "tmp\p5_missing_dimred_sources\$stamp"
New-Item -ItemType Directory -Force -Path $logDir | Out-Null

if ([string]::IsNullOrWhiteSpace($MatlabExe)) {
    $candidate = "C:\Program Files\MATLAB\R2025b\bin\matlab.exe"
    if (Test-Path -LiteralPath $candidate) {
        $MatlabExe = $candidate
    } else {
        $MatlabExe = "matlab.exe"
    }
}

$commonP5 = @"
autodl_root='E:\autodl_results_new';
condition_tag_mode='condition';
make_overview_plot=false;
make_state_space_plot=false;
make_consensus_state_space_plot=false;
make_spectrum_diagnostics=false;
make_top30_window_plots=false;
make_thresholded_density=false;
make_thresholded_events=false;
make_dimred_thresholded_density=false;
make_dimred_thresholded_events=false;
continue_on_error=true;
run('scripts/script_run_one_cfg_blp_eigenfunction_reduction.m');
"@ -replace "\r?\n", " "

$steps = @(
    @{
        Name = "e10fV1_complex_split_mds_k03_k08"
        Code = "cfg_name='E10fV1'; condition_key_filter={'e10fV1|complex_split|projected_vlambda'}; method_filter={'mds'}; component_count_sweep=3:8; spectrum_component_count_sweep=3:8; $commonP5"
    },
    @{
        Name = "e10fV1_abs_umap_k05_k07"
        Code = "cfg_name='E10fV1'; condition_key_filter={'e10fV1|abs|projected_vlambda'}; method_filter={'umap'}; component_count_sweep=5:7; spectrum_component_count_sweep=5:7; $commonP5"
    },
    @{
        Name = "e10fV1_complex_split_umap_k03_k07"
        Code = "cfg_name='E10fV1'; condition_key_filter={'e10fV1|complex_split|projected_vlambda'}; method_filter={'umap'}; component_count_sweep=3:7; spectrum_component_count_sweep=3:7; $commonP5"
    },
    @{
        Name = "f12m01_abs_complex_mds_umap_k03_k07"
        Code = "cfg_name='F12m01'; condition_key_filter={'f12m01|abs|projected_vlambda','f12m01|complex_split|projected_vlambda'}; method_filter={'mds','umap'}; component_count_sweep=3:7; spectrum_component_count_sweep=3:7; $commonP5"
    }
)

if (-not $NoPeakStats) {
    $peakTags = @(
        "mds_k03", "mds_k04", "mds_k05", "mds_k06", "mds_k07", "mds_k08",
        "umap_k03", "umap_k04", "umap_k05", "umap_k06", "umap_k07", "umap_k08"
    )
    $peakTagLiteral = ($peakTags | ForEach-Object { "'$_'" }) -join ","
    $steps += @{
        Name = "peak_state_stats_e10fV1_f12m01"
        Code = "dataset_stems={'e10fV1','f12m01'}; variant_filter={'abs_projected_vlambda','complex_split_projected_vlambda'}; method_tag_filter={$peakTagLiteral}; force_recompute=false; save_figures=false; continue_on_error=true; run('scripts/script_run_current_pipeline5_peak_state_all_datasets.m');"
    }
}

Write-Host "P5 missing dimred-source backfill"
Write-Host "Repo: $repoRoot"
Write-Host "MATLAB: $MatlabExe"
Write-Host "Execute: $Execute"
Write-Host "Peak stats: $(-not $NoPeakStats)"
Write-Host "Logs: $logDir"
Write-Host ""

foreach ($step in $steps) {
    $batch = "cd('$repoMatlab'); $($step.Code)"
    $logFile = Join-Path $logDir "$($step.Name).log"
    if (-not $Execute) {
        Write-Host "[dry-run] $($step.Name)"
        Write-Host "  log: $logFile"
        Write-Host "  matlab -batch `"$batch`""
        continue
    }

    Write-Host "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] START $($step.Name)"
    $oldErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    try {
        & $MatlabExe -batch $batch *> $logFile
        $exitCode = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $oldErrorActionPreference
    }
    Write-Host "[$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')] EXIT $exitCode $($step.Name)"
    if ($exitCode -ne 0) {
        Write-Host "Failed log: $logFile"
        throw "Step failed: $($step.Name)"
    }
}

Write-Host ""
if ($Execute) {
    Write-Host "Done. Logs: $logDir"
} else {
    Write-Host "Dry-run only. Re-run with -Execute to launch the backfill."
}
