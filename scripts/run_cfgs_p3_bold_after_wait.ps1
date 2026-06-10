param(
    [string[]]$CfgNames = @("E10aw1", "E10bv1", "K13m17", "K13m23"),
    [string[]]$ObservableModes = @("global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean"),
    [int]$WaitForPid = 0,
    [int]$PollSeconds = 60,
    [string]$MatlabExe = "C:\Program Files\MATLAB\R2023a\bin\matlab.exe",
    [switch]$ForceRecompute,
    [switch]$SkipQC,
    [switch]$NoWait,
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
$logRoot = Join-Path $repoRoot "tmp\p3_bold_cfgs\$stamp"
New-Item -ItemType Directory -Force -Path $logRoot | Out-Null
$supervisorLog = Join-Path $logRoot "supervisor.log"

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

function Invoke-LoggedExternal {
    param(
        [string]$Name,
        [string]$FilePath,
        [string[]]$Arguments
    )

    $logPath = Join-Path $logRoot "$Name.log"
    Write-Log ("STEP {0}" -f $Name)
    Write-Log ("LOG  {0}" -f $logPath)
    Write-Log ("CMD  {0} {1}" -f $FilePath, ($Arguments -join " "))

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

try {
    Write-Log "P3 BOLD cfg runner"
    Write-Log ("Repo: {0}" -f $repoRoot)
    Write-Log ("CfgNames: {0}" -f ($CfgNames -join ", "))
    Write-Log ("ObservableModes: {0}" -f ($ObservableModes -join ", "))
    Write-Log ("MATLAB: {0}" -f $MatlabExe)
    Write-Log ("Logs: {0}" -f $logRoot)

    if (-not (Test-Path -LiteralPath $MatlabExe)) {
        throw "MATLAB executable not found: $MatlabExe"
    }

    if ($WaitForPid -gt 0 -and -not $NoWait) {
        while ($true) {
            $proc = Get-Process -Id $WaitForPid -ErrorAction SilentlyContinue
            if ($null -eq $proc) {
                Write-Log ("Wait PID {0}: no longer running; starting P3." -f $WaitForPid)
                break
            }
            Write-Log ("Waiting for PID {0} ({1}); elapsed {2}." -f $WaitForPid, $proc.ProcessName, $proc.TotalProcessorTime)
            Start-Sleep -Seconds $PollSeconds
        }
    }

    $repoLit = ConvertTo-MatlabStringLiteral $repoRoot
    $modesLit = ConvertTo-MatlabCellLiteral $ObservableModes
    $forceLit = ConvertTo-MatlabLogical $ForceRecompute.IsPresent

    for ($i = 0; $i -lt $CfgNames.Count; $i++) {
        $cfg = $CfgNames[$i]
        $cfgLit = ConvertTo-MatlabStringLiteral $cfg
        $prefix = "{0:D2}_{1}" -f ($i + 1), $cfg

        $buildCode = "cd($repoLit); set(groot,'defaultFigureVisible','off'); close all force; cfg_name=$cfgLit; observable_modes=$modesLit; force_recompute=$forceLit; save_precision='single'; run('scripts/script_build_one_cfg_bold_observables.m');"
        Invoke-MatlabBatch -Name "${prefix}_P3_observables" -Code $buildCode

        if (-not $SkipQC) {
            $qcCode = "cd($repoLit); set(groot,'defaultFigureVisible','off'); close all force; cfg_name=$cfgLit; observable_modes=$modesLit; save_qc=true; qc_visible='off'; run('scripts/script_plot_one_cfg_bold_pre_reskoopnet_qc.m');"
            Invoke-MatlabBatch -Name "${prefix}_P3_QC" -Code $qcCode
        }
    }

    Write-Log "P3 BOLD cfg runner finished."
} catch {
    Write-Log ("ERROR: {0}" -f $_.Exception.Message)
    throw
}
