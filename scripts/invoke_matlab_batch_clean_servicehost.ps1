param(
    [Parameter(Mandatory = $true)]
    [string]$BatchCommand,

    [string]$MatlabExe = "C:\Program Files\MATLAB\R2025b\bin\matlab.exe"
)

$ErrorActionPreference = "Stop"

function Reset-MathWorksEndpoint {
    Get-Process |
        Where-Object { $_.ProcessName -match '^MathWorksServiceHost|^MathWorksServiceHost-Monitor|^MATLAB$' } |
        Stop-Process -Force -ErrorAction SilentlyContinue

    Start-Sleep -Seconds 2

    $endpoint = Join-Path $env:LOCALAPPDATA "MathWorks\mwEndpointRegistry"
    if (Test-Path $endpoint) {
        $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
        $parent = Split-Path $endpoint -Parent
        Move-Item -LiteralPath $endpoint -Destination (Join-Path $parent "mwEndpointRegistry_backup_codex_$stamp") -Force
    }
}

Reset-MathWorksEndpoint
try {
    & $MatlabExe -batch $BatchCommand
    $exitCode = $LASTEXITCODE
}
finally {
    Reset-MathWorksEndpoint
}

exit $exitCode
