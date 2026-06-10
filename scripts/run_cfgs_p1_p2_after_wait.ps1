param(
    [string[]]$CfgNames = @("E10aw1", "E10bv1", "K13m17", "K13m23"),
    [int]$WaitForPid = 0,
    [int]$PollSeconds = 60,
    [string]$MatlabExe = "C:\Program Files\MATLAB\R2023a\bin\matlab.exe",
    [int]$ChunkSize = 200000,
    [switch]$ForceRecompute,
    [switch]$MakeTopWindowPlots,
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

$repoRoot = Split-Path -Parent $PSScriptRoot
$stamp = Get-Date -Format "yyyyMMdd_HHmmss"
$logRoot = Join-Path $repoRoot "tmp\p1_p2_cfgs\$stamp"
New-Item -ItemType Directory -Force -Path $logRoot | Out-Null
$supervisorLog = Join-Path $logRoot "supervisor.log"
$latestFile = Join-Path $repoRoot "tmp\p1_p2_cfgs\latest_logdir.txt"
Set-Content -Path $latestFile -Value $logRoot -Encoding ASCII

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

    "" | Set-Content -Path $logPath -Encoding ASCII
    $oldErrorActionPreference = $ErrorActionPreference
    $ErrorActionPreference = "Continue"
    try {
        & $FilePath @Arguments 2>&1 | ForEach-Object {
            $text = ($_ | Out-String).TrimEnd()
            if ($text.Length -gt 0) {
                Write-Host $text
                Add-Content -Path $logPath -Value $text -Encoding ASCII
            }
        }
        $code = $LASTEXITCODE
    } finally {
        $ErrorActionPreference = $oldErrorActionPreference
    }

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
    Write-Log "P1+P2 cfg runner"
    Write-Log ("Repo: {0}" -f $repoRoot)
    Write-Log ("CfgNames: {0}" -f ($CfgNames -join ", "))
    Write-Log ("ChunkSize: {0}" -f $ChunkSize)
    Write-Log ("Logs: {0}" -f $logRoot)

    if (-not (Test-Path -LiteralPath $MatlabExe)) {
        throw "MATLAB executable not found: $MatlabExe"
    }

    if ($WaitForPid -gt 0 -and -not $NoWait) {
        while ($true) {
            $proc = Get-Process -Id $WaitForPid -ErrorAction SilentlyContinue
            if ($null -eq $proc) {
                Write-Log ("Wait PID {0}: no longer running; starting P1+P2." -f $WaitForPid)
                break
            }
            Write-Log ("Waiting for PID {0} ({1}); elapsed {2}." -f $WaitForPid, $proc.ProcessName, $proc.TotalProcessorTime)
            Start-Sleep -Seconds $PollSeconds
        }
    }

    $repoLit = ConvertTo-MatlabStringLiteral $repoRoot
    $cfgCellLit = ConvertTo-MatlabCellLiteral $CfgNames
    $forceLit = ConvertTo-MatlabLogical $ForceRecompute.IsPresent
    $topPlotLit = ConvertTo-MatlabLogical $MakeTopWindowPlots.IsPresent
    $p1InternalLog = ConvertTo-MatlabStringLiteral (Join-Path $logRoot "01_P1_internal.log")
    $p2InternalLog = ConvertTo-MatlabStringLiteral (Join-Path $logRoot "02_P2_internal.log")

    $p1Code = "cd($repoLit); set(groot,'defaultFigureVisible','off'); close all force; cfg_names=$cfgCellLit; force_recompute=$forceLit; save_precision='single'; chunk_size=$ChunkSize; dict_modes={'abs','complex_split'}; load_metadata_only=true; cache_raw_to_disk=false; run_log_file=$p1InternalLog; run('scripts/script_preprocess_cfgs_to_observables_streamed.m');"
    Invoke-MatlabBatch -Name "01_P1_spectrogram_dictionary" -Code $p1Code

    $p2Code = "cd($repoLit); set(groot,'defaultFigureVisible','off'); close all force; cfg_names=$cfgCellLit; force_event_recompute=$forceLit; force_density_recompute=$forceLit; force_consensus_recompute=$forceLit; force_summary_recompute=$forceLit; force_window_recompute=$forceLit; make_top_window_plots=$topPlotLit; run_log_file=$p2InternalLog; run('scripts/script_run_cfgs_to_consensus_state_top_windows.m');"
    Invoke-MatlabBatch -Name "02_P2_consensus_states" -Code $p2Code

    Write-Log "P1+P2 cfg runner finished."
} catch {
    Write-Log ("ERROR: {0}" -f $_.Exception.Message)
    throw
}
