param(
    [string]$ExperimentName = "e10gw1_batch60_tonight_fix1",
    [string]$RunNameBase = "mlp_obs",
    [string]$LocalDownloadRoot = "E:\autodl_results\e10gw1\mlp",
    [switch]$LogsOnly,
    [switch]$SkipCheckpoints,
    [string]$SshHost = "connect.westb.seetacloud.com",
    [int]$SshPort = 19241,
    [string]$SshUser = "root",
    [string]$SshKey = "$HOME\.ssh\id_ed25519_autodl"
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$remoteOutputParent = "/root/autodl-tmp/outputs"
$remoteCheckpointParent = "/root/autodl-tmp/checkpoints"
$remoteLogParent = "/root/autodl-tmp/logs"
$runPrefix = "{0}_{1}_" -f $RunNameBase, $ExperimentName

function Resolve-OpenSshExecutable {
    param(
        [string]$LeafName,
        [string]$FallbackCommand
    )

    $windir = if ($env:WINDIR) { $env:WINDIR } else { "C:\Windows" }
    $candidate = Join-Path (Join-Path $windir "System32\OpenSSH") $LeafName
    if (Test-Path -LiteralPath $candidate) {
        return $candidate
    }
    return $FallbackCommand
}

function Quote-BashSingle {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Value
    )

    $replacement = "'" + '"' + "'" + '"' + "'"
    return "'" + $Value.Replace("'", $replacement) + "'"
}

function Invoke-RemoteBash {
    param(
        [Parameter(Mandatory = $true)]
        [string]$ScriptText
    )

    $normalizedScript = $ScriptText -replace "`r`n", "`n"
    $normalizedScript = $normalizedScript -replace "`r", "`n"
    $lines = $normalizedScript | & $script:sshExe @script:sshArgs "tr -d '\r' | bash -s"
    if ($LASTEXITCODE -ne 0) {
        throw "Remote command failed."
    }
    return @($lines)
}

function Get-RemoteRunDirectories {
    param(
        [Parameter(Mandatory = $true)]
        [string]$RemoteParent
    )

    $remoteParentQuoted = Quote-BashSingle -Value $RemoteParent
    $runPrefixQuoted = Quote-BashSingle -Value $script:runPrefix
    $scriptText = @(
        "set -e"
        "parent=$remoteParentQuoted"
        "prefix=$runPrefixQuoted"
        'if [ ! -d "$parent" ]; then'
        '  exit 0'
        'fi'
        'find "$parent" -mindepth 1 -maxdepth 1 -type d -name "${prefix}*" -print | sort'
    ) -join "`n"

    $lines = Invoke-RemoteBash -ScriptText $scriptText
    return @($lines | Where-Object { $_ -and $_.Trim() })
}

function Get-RemoteMatchingProcesses {
    $experimentQuoted = Quote-BashSingle -Value $script:ExperimentName
    $scriptText = @"
pgrep -af $experimentQuoted || true
"@

    $lines = Invoke-RemoteBash -ScriptText $scriptText
    return @($lines | Where-Object { $_ -and $_.Trim() })
}

function Copy-RemoteDirectories {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Kind,
        [Parameter(Mandatory = $true)]
        [string]$LocalParent,
        [Parameter(Mandatory = $true)]
        [string[]]$RemoteDirs
    )

    foreach ($remoteDir in $RemoteDirs) {
        Write-Host ("Downloading {0}: {1}" -f $Kind, $remoteDir)
        & $script:scpExe @script:scpBaseArgs -r ("{0}@{1}:{2}" -f $script:SshUser, $script:SshHost, $remoteDir) $LocalParent
        if ($LASTEXITCODE -ne 0) {
            throw "Download failed for ${Kind}: $remoteDir"
        }
    }
}

$sshExe = Resolve-OpenSshExecutable -LeafName "ssh.exe" -FallbackCommand "ssh"
$scpExe = Resolve-OpenSshExecutable -LeafName "scp.exe" -FallbackCommand "scp"

$sshArgs = @(
    "-p", [string]$SshPort,
    "-i", $SshKey,
    ("{0}@{1}" -f $SshUser, $SshHost)
)

$scpBaseArgs = @(
    "-P", [string]$SshPort,
    "-i", $SshKey
)

$localOutputsDir = Join-Path $LocalDownloadRoot "outputs"
$localLogsDir = Join-Path $LocalDownloadRoot "logs"
$localCheckpointsDir = Join-Path $LocalDownloadRoot "checkpoints"

foreach ($dir in @($LocalDownloadRoot, $localOutputsDir, $localLogsDir, $localCheckpointsDir)) {
    New-Item -ItemType Directory -Force -Path $dir | Out-Null
}

Write-Host ""
Write-Host ("=" * 80)
Write-Host "Downloading E10gW1 AutoDL results"
Write-Host ("ExperimentName: {0}" -f $ExperimentName)
Write-Host ("RunPrefix: {0}" -f $runPrefix)
Write-Host ("LocalDownloadRoot: {0}" -f $LocalDownloadRoot)
Write-Host ("SSH: {0}@{1}:{2}" -f $SshUser, $SshHost, $SshPort)
Write-Host ("=" * 80)

$remoteProcesses = @(Get-RemoteMatchingProcesses)
if ($remoteProcesses.Count -gt 0) {
    Write-Warning "The matching experiment still appears to be running remotely. Downloaded files may be incomplete."
    $remoteProcesses | ForEach-Object { Write-Host ("ACTIVE: {0}" -f $_) }
}

$remoteOutputs = @()
if (-not $LogsOnly) {
    $remoteOutputs = @(Get-RemoteRunDirectories -RemoteParent $remoteOutputParent)
}
$remoteLogs = @(Get-RemoteRunDirectories -RemoteParent $remoteLogParent)
$remoteCheckpoints = @()
if ((-not $LogsOnly) -and (-not $SkipCheckpoints)) {
    $remoteCheckpoints = @(Get-RemoteRunDirectories -RemoteParent $remoteCheckpointParent)
}

Write-Host ""
Write-Host ("Remote output dirs found: {0}" -f $remoteOutputs.Count)
$remoteOutputs | ForEach-Object { Write-Host ("  {0}" -f $_) }
Write-Host ("Remote log dirs found: {0}" -f $remoteLogs.Count)
$remoteLogs | ForEach-Object { Write-Host ("  {0}" -f $_) }
if ((-not $LogsOnly) -and (-not $SkipCheckpoints)) {
    Write-Host ("Remote checkpoint dirs found: {0}" -f $remoteCheckpoints.Count)
    $remoteCheckpoints | ForEach-Object { Write-Host ("  {0}" -f $_) }
}

if (($remoteOutputs.Count -eq 0) -and ($remoteLogs.Count -eq 0) -and ($SkipCheckpoints -or $remoteCheckpoints.Count -eq 0)) {
    throw "No remote run directories matched prefix: $runPrefix"
}

if ($remoteOutputs.Count -gt 0) {
    Copy-RemoteDirectories -Kind "outputs" -LocalParent $localOutputsDir -RemoteDirs $remoteOutputs
}

if ($remoteLogs.Count -gt 0) {
    Copy-RemoteDirectories -Kind "logs" -LocalParent $localLogsDir -RemoteDirs $remoteLogs
}

if ((-not $LogsOnly) -and (-not $SkipCheckpoints) -and ($remoteCheckpoints.Count -gt 0)) {
    Copy-RemoteDirectories -Kind "checkpoints" -LocalParent $localCheckpointsDir -RemoteDirs $remoteCheckpoints
}

Write-Host ""
Write-Host ("=" * 80)
Write-Host "Local directory listing"
Write-Host ("=" * 80)
Get-ChildItem -Path $LocalDownloadRoot | Select-Object Name, LastWriteTime

Write-Host ""
Write-Host "E10gW1 remote result download finished successfully."
