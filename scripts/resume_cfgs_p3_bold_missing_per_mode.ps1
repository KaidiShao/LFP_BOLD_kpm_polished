param(
    [string[]]$CfgNames = @("K13m17", "K13m23"),
    [string[]]$ObservableModes = @("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean"),
    [string]$MatlabExe = "C:\Program Files\MATLAB\R2023a\bin\matlab.exe",
    [int]$WaitForPid = 0,
    [int]$PollSeconds = 60,
    [switch]$ForceRecompute,
    [switch]$SkipQC,
    [switch]$ContinueOnError,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$CfgNames = @(
    $CfgNames |
        ForEach-Object { $_ -split "," } |
        ForEach-Object { $_.Trim() } |
        Where-Object { $_.Length -gt 0 }
)
$ObservableModes = @(
    $ObservableModes |
        ForEach-Object { $_ -split "," } |
        ForEach-Object { $_.Trim() } |
        Where-Object { $_.Length -gt 0 }
)

$repoRoot = Split-Path -Parent $PSScriptRoot
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$logRoot = Join-Path $repoRoot "tmp\p3_bold_per_mode\$stamp"
New-Item -ItemType Directory -Force -Path $logRoot | Out-Null
$supervisorLog = Join-Path $logRoot "supervisor.log"

$DatasetMap = @{
    "E10aw1" = @{ Stem = "e10aw1"; DatasetId = "E10.aW1" }
    "E10bv1" = @{ Stem = "e10bv1"; DatasetId = "E10.bv1" }
    "K13m17" = @{ Stem = "k13m17"; DatasetId = "K13.m17" }
    "K13m23" = @{ Stem = "k13m23"; DatasetId = "K13.m23" }
}

function Write-Log {
    param([string]$Message)
    $line = "[{0}] {1}" -f (Get-Date -Format "yyyy-MM-dd HH:mm:ss"), $Message
    Write-Host $line
    Add-Content -Path $supervisorLog -Value $line -Encoding ASCII
}

function ConvertTo-MatlabStringLiteral {
    param([string]$Value)
    return "'" + ($Value -replace "'", "''") + "'"
}

function ConvertTo-MatlabCellLiteral {
    param([string[]]$Values)
    $items = @()
    foreach ($value in $Values) {
        $items += (ConvertTo-MatlabStringLiteral $value)
    }
    return "{" + ($items -join ",") + "}"
}

function ConvertTo-MatlabLogical {
    param([bool]$Value)
    if ($Value) {
        return "true"
    }
    return "false"
}

function Test-AnyQcOutput {
    param([string]$QcDir)
    if (Test-Path -LiteralPath $QcDir -PathType Container) {
        $png = Get-ChildItem -LiteralPath $QcDir -Filter "*.png" -File -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($null -ne $png) {
            return $true
        }
    }
    return $false
}

function Invoke-LoggedExternal {
    param(
        [string]$Name,
        [string]$FilePath,
        [string[]]$Arguments
    )

    $safeName = $Name -replace "[^A-Za-z0-9_.-]", "_"
    $logPath = Join-Path $logRoot "$safeName.log"
    Write-Log ("STEP {0}" -f $Name)
    Write-Log ("LOG  {0}" -f $logPath)
    Write-Log ("CMD  {0} {1}" -f $FilePath, ($Arguments -join " "))

    if ($DryRun) {
        "DRY RUN: $FilePath $($Arguments -join ' ')" | Set-Content -Path $logPath -Encoding ASCII
        return 0
    }

    & $FilePath @Arguments 2>&1 | Tee-Object -FilePath $logPath | ForEach-Object {
        Write-Host $_
    }
    $exitCode = $LASTEXITCODE
    return $exitCode
}

function Invoke-MatlabBatch {
    param(
        [string]$Name,
        [string]$Code
    )

    $code = Invoke-LoggedExternal -Name $Name -FilePath $MatlabExe -Arguments @("-batch", $Code)
    if ($code -ne 0) {
        throw "Step '$Name' failed with exit code $code."
    }
}

try {
    Write-Log "P3 BOLD per-mode resume runner"
    Write-Log ("Repo: {0}" -f $repoRoot)
    Write-Log ("CfgNames: {0}" -f ($CfgNames -join ", "))
    Write-Log ("ObservableModes: {0}" -f ($ObservableModes -join ", "))
    Write-Log ("MATLAB: {0}" -f $MatlabExe)
    Write-Log ("Logs: {0}" -f $logRoot)
    Write-Log ("ForceRecompute: {0}; SkipQC: {1}; ContinueOnError: {2}; DryRun: {3}" -f `
        $ForceRecompute.IsPresent, $SkipQC.IsPresent, $ContinueOnError.IsPresent, $DryRun.IsPresent)

    if (-not (Test-Path -LiteralPath $MatlabExe)) {
        throw "MATLAB executable not found: $MatlabExe"
    }

    if ($WaitForPid -gt 0) {
        while ($true) {
            $proc = Get-Process -Id $WaitForPid -ErrorAction SilentlyContinue
            if ($null -eq $proc) {
                Write-Log ("Wait PID {0}: no longer running; starting P3 resume." -f $WaitForPid)
                break
            }
            Write-Log ("Waiting for PID {0} ({1}); CPU {2}." -f $WaitForPid, $proc.ProcessName, $proc.TotalProcessorTime)
            Start-Sleep -Seconds $PollSeconds
        }
    }

    $repoLit = ConvertTo-MatlabStringLiteral $repoRoot
    $forceLit = ConvertTo-MatlabLogical $ForceRecompute.IsPresent
    $processedRoot = "E:\DataPons_processed"
    $attempted = 0
    $skipped = 0
    $failed = 0

    foreach ($cfg in $CfgNames) {
        if (-not $DatasetMap.ContainsKey($cfg)) {
            throw "Unknown cfg mapping: $cfg. Add it to DatasetMap in this script."
        }

        $stem = $DatasetMap[$cfg].Stem
        $datasetId = $DatasetMap[$cfg].DatasetId
        $obsDir = Join-Path $processedRoot "$stem\pipeline3_bold_observables"
        $qcRoot = Join-Path $processedRoot "$stem\pipeline3_figures_bold_pre_reskoopnet_qc"

        foreach ($mode in $ObservableModes) {
            $obsFile = Join-Path $obsDir ("{0}_bold_observables_{1}.mat" -f $datasetId, $mode)
            $qcDir = Join-Path $qcRoot $mode
            $hasObs = (Test-Path -LiteralPath $obsFile -PathType Leaf)
            $hasQc = Test-AnyQcOutput -QcDir $qcDir

            if ($hasObs -and $hasQc -and -not $ForceRecompute.IsPresent) {
                Write-Log ("SKIP {0} {1}: observable and QC already present." -f $cfg, $mode)
                $skipped++
                continue
            }

            $cfgLit = ConvertTo-MatlabStringLiteral $cfg
            $modeLit = ConvertTo-MatlabCellLiteral @($mode)
            $stepPrefix = "{0}_{1}" -f $cfg, $mode

            if (-not $hasObs -or $ForceRecompute.IsPresent) {
                $attempted++
                try {
                    $buildCode = "cd($repoLit); set(groot,'defaultFigureVisible','off'); close all force; cfg_name=$cfgLit; observable_modes=$modeLit; force_recompute=$forceLit; save_precision='single'; run('scripts/script_build_one_cfg_bold_observables.m');"
                    Invoke-MatlabBatch -Name "${stepPrefix}_P3_observables" -Code $buildCode
                } catch {
                    $failed++
                    Write-Log ("FAILED {0} {1} observables: {2}" -f $cfg, $mode, $_.Exception.Message)
                    if (-not $ContinueOnError.IsPresent) {
                        throw
                    }
                    continue
                }
            } else {
                Write-Log ("SKIP {0} {1}: observable already present." -f $cfg, $mode)
                $skipped++
            }

            if ($SkipQC.IsPresent) {
                continue
            }

            $hasObsAfter = (Test-Path -LiteralPath $obsFile -PathType Leaf)
            $hasQcAfter = Test-AnyQcOutput -QcDir $qcDir
            if ($hasObsAfter -and (-not $hasQcAfter -or $ForceRecompute.IsPresent)) {
                $attempted++
                try {
                    $qcCode = "cd($repoLit); set(groot,'defaultFigureVisible','off'); close all force; cfg_name=$cfgLit; observable_modes=$modeLit; save_qc=true; qc_visible='off'; run('scripts/script_plot_one_cfg_bold_pre_reskoopnet_qc.m');"
                    Invoke-MatlabBatch -Name "${stepPrefix}_P3_QC" -Code $qcCode
                } catch {
                    $failed++
                    Write-Log ("FAILED {0} {1} QC: {2}" -f $cfg, $mode, $_.Exception.Message)
                    if (-not $ContinueOnError.IsPresent) {
                        throw
                    }
                }
            } elseif ($hasQcAfter) {
                Write-Log ("SKIP {0} {1}: QC already present." -f $cfg, $mode)
                $skipped++
            } else {
                Write-Log ("SKIP {0} {1}: QC waiting for observable file." -f $cfg, $mode)
                $skipped++
            }
        }
    }

    Write-Log ("P3 BOLD per-mode resume finished. attempted={0}; skipped={1}; failed={2}" -f $attempted, $skipped, $failed)
    if ($failed -gt 0) {
        exit 1
    }
} catch {
    Write-Log ("ERROR: {0}" -f $_.Exception.Message)
    throw
}
