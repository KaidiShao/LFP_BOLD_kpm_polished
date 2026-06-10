param(
    [string]$MatlabExe = "C:\Program Files\MATLAB\R2025b\bin\matlab.exe",
    [string]$PythonExe = "python",
    [string]$SshHost = "connect.westb.seetacloud.com",
    [string]$SshUser = "root",
    [int]$SshPort = 19241,
    [string]$SshKey = "$HOME\.ssh\id_ed25519_autodl",
    [string]$RemotePython = "/root/miniconda3/envs/reskoopnet/bin/python",
    [string]$ExperimentName = "blp_vlambda_mainline_torchlike_20260515_e10gw1_seed1234",
    [string]$RemoteDataSubdir = "e10gw1_blp_rebuild_20260515",
    [string]$LocalDownloadRoot = "E:\autodl_results_new\e10gw1\mlp",
    [int]$Epochs = 50,
    [switch]$ResumeP4,
    [switch]$SkipP1P2,
    [switch]$SkipP4,
    [switch]$SkipP5P6,
    [switch]$SkipP7,
    [switch]$SkipP8,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = Split-Path -Parent $PSScriptRoot
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$logRoot = Join-Path $repoRoot "tmp\e10gw1_rebuild_sequence\$stamp"
New-Item -ItemType Directory -Force -Path $logRoot | Out-Null

function ConvertTo-MatlabStringLiteral {
    param([string]$Value)
    return "'" + ($Value -replace "'", "''") + "'"
}

function Invoke-LoggedExternal {
    param(
        [string]$Name,
        [string]$FilePath,
        [string[]]$Arguments
    )

    $logPath = Join-Path $logRoot "$Name.log"
    Write-Host ""
    Write-Host ("=" * 80)
    Write-Host "STEP: $Name"
    Write-Host "LOG : $logPath"
    Write-Host ("=" * 80)
    Write-Host "$FilePath $($Arguments -join ' ')"

    if ($DryRun) {
        "DRY RUN: $FilePath $($Arguments -join ' ')" | Set-Content -Path $logPath -Encoding ASCII
        return
    }

    & $FilePath @Arguments 2>&1 | Tee-Object -FilePath $logPath
    $code = $LASTEXITCODE
    if ($code -ne 0) {
        throw "Step '$Name' failed with exit code $code. See log: $logPath"
    }
}

function Invoke-MatlabBatch {
    param(
        [string]$Name,
        [string]$Code
    )
    Invoke-LoggedExternal -Name $Name -FilePath $MatlabExe -Arguments @("-batch", $Code)
}

$repoLit = ConvertTo-MatlabStringLiteral $repoRoot

Write-Host "E10.gW1 sequential rebuild"
Write-Host "Repo: $repoRoot"
Write-Host "Logs: $logRoot"
Write-Host "P4 experiment: $ExperimentName"
Write-Host "P4 local download root: $LocalDownloadRoot"

if (-not $SkipP1P2) {
    $p1Code = "cd($repoLit); set(groot,'defaultFigureVisible','off'); cfg_name='E10gW1'; force_recompute=true; save_precision='single'; chunk_size=200000; dict_modes={'abs','complex_split'}; load_metadata_only=true; cache_raw_to_disk=false; run('scripts/script_preprocess_one_cfg_to_observables_streamed.m');"
    Invoke-MatlabBatch -Name "01_P1_spectrogram_dictionary" -Code $p1Code

    $p2Code = "cd($repoLit); set(groot,'defaultFigureVisible','off'); cfg_name='E10gW1'; force_event_recompute=true; force_density_recompute=true; force_consensus_recompute=true; force_summary_recompute=true; force_window_recompute=true; make_top_window_plots=false; run('scripts/script_run_one_cfg_to_consensus_state_top_windows.m');"
    Invoke-MatlabBatch -Name "02_P2_consensus_states" -Code $p2Code
}

if (-not $SkipP4) {
    $dictRoot = "E:\DataPons_processed\e10gw1\pipeline1_reskoopnet_dictionary"
    $absFile = Join-Path $dictRoot "e10gw1_low50_high250_g2_abs_single.mat"
    $absCsv = Join-Path $dictRoot "e10gw1_low50_high250_g2_abs_single_obs_info.csv"
    $complexFile = Join-Path $dictRoot "e10gw1_low50_high250_g2_complex_split_single.mat"
    $complexCsv = Join-Path $dictRoot "e10gw1_low50_high250_g2_complex_split_single_obs_info.csv"

    foreach ($path in @($absFile, $absCsv, $complexFile, $complexCsv)) {
        if (-not $DryRun -and -not (Test-Path -LiteralPath $path)) {
            throw "Required P1 dictionary file is missing before P4: $path"
        }
    }

    $controller = Join-Path $repoRoot "python_scripts\autodl\dataset_batch_controller_autodl_reskoopnet_mlp.py"
    $p4Args = @(
        $controller,
        "--ssh-host", $SshHost,
        "--ssh-user", $SshUser,
        "--ssh-port", [string]$SshPort,
        "--ssh-key", $SshKey,
        "--dataset-stem", "e10gw1",
        "--experiment-name", $ExperimentName,
        "--observable-modes", "abs", "complex_split",
        "--residual-forms", "projected_vlambda",
        "--remote-data-subdir", $RemoteDataSubdir,
        "--local-download-root", $LocalDownloadRoot,
        "--local-data-file-abs", $absFile,
        "--local-obs-info-file-abs", $absCsv,
        "--local-data-file-complex-split", $complexFile,
        "--local-obs-info-file-complex-split", $complexCsv,
        "--file-type", ".mat",
        "--field-name", "obs",
        "--remote-python", $RemotePython,
        "--epochs", [string]$Epochs,
        "--recover-completed-remote-runs",
        "--skip-completed-local-runs",
        "--delete-remote-run-after-download",
        "--delete-remote-input-after-observable"
    )
    if ($ResumeP4) {
        $p4Args += "--resume"
    } else {
        $p4Args += "--fresh-checkpoints"
    }
    Invoke-LoggedExternal -Name "03_P4_autodl_blp_abs_complex_split" -FilePath $PythonExe -Arguments $p4Args
}

$psExe = Join-Path $PSHOME "powershell.exe"

if (-not $SkipP5P6) {
    $p56Args = @(
        "-NoProfile",
        "-ExecutionPolicy", "Bypass",
        "-File", (Join-Path $repoRoot "scripts\backfill_current5_minimal.ps1"),
        "-Execute",
        "-Datasets", "e10gw1",
        "-SkipPreBold",
        "-SkipPostBold",
        "-ContinueOnError",
        "-MatlabExe", $MatlabExe
    )
    Invoke-LoggedExternal -Name "04_P5_P6_blp" -FilePath $psExe -Arguments $p56Args
}

if (-not $SkipP7) {
    $p7Code = "cd($repoLit); set(groot,'defaultFigureVisible','off'); cfg_name='E10gW1'; observable_modes={'global_svd100','gsvd100_ds','HP_svd100','roi_mean'}; autodl_roots={'E:\autodl_results_local\bold_wsl','E:\autodl_results\bold'}; force_recompute=true; make_main_plot=true; compute_deconv=true; make_deconv_plot=true; make_timescale_plot=true; make_intrinsic_activation_maps=true; make_intrinsic_roi_summary=true; close_figures_after_each_run=true; run('scripts/script_run_one_cfg_bold_reskoopnet_postprocessing.m');"
    Invoke-MatlabBatch -Name "05_P7_bold_postprocessing_refresh" -Code $p7Code
}

if (-not $SkipP8) {
    $p8Args = @(
        "-NoProfile",
        "-ExecutionPolicy", "Bypass",
        "-File", (Join-Path $repoRoot "scripts\backfill_current5_minimal.ps1"),
        "-Execute",
        "-Datasets", "e10gw1",
        "-SkipPreBold",
        "-SkipBlp",
        "-ContinueOnError",
        "-MatlabExe", $MatlabExe
    )
    Invoke-LoggedExternal -Name "06_P8_coupling" -FilePath $psExe -Arguments $p8Args
}

Write-Host ""
Write-Host "E10.gW1 sequential rebuild finished."
Write-Host "Logs: $logRoot"
