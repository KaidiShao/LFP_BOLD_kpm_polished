param(
    [string]$OutputName = "turbulence_precision_safe_fp32",
    [string]$TrainingPolicy = "float32",
    [string]$AnalysisDtype = "float64",
    [string]$GramDtype = "float64",
    [string]$SpectralDtype = "float64",
    [double]$FocusCorrThreshold = 0.995,
    [string]$PythonExe = "/home/kdshao/anaconda3/bin/python",
    [string]$Distro = "Ubuntu-22.04",
    [int]$Seed = 1234,
    [int]$OuterEpochs = 100,
    [int]$InnerEpochs = 2,
    [int]$BatchSize = 500,
    [double]$LearningRate = 1e-4,
    [double]$Reg = 0.1,
    [switch]$DryRun
)

Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$repoRoot = "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished"
$dataMat = Join-Path $repoRoot "tmp\turbulence_lab\original_data\pressure_data.mat"
$officialMat = Join-Path $repoRoot "tmp\turbulence_lab\original_data\turbulence_resdmd_250basis_official.mat"
$runScript = Join-Path $repoRoot "tmp\active_validation\run_active_turbulence.py"
$compareScript = Join-Path $repoRoot "tmp\active_validation\compare_turbulence_vs_official.py"
$outputDir = Join-Path $repoRoot "tmp\active_validation\outputs\$OutputName"

function Convert-ToWslPath {
    param(
        [Parameter(Mandatory = $true)]
        [string]$Path
    )

    $fullPath = [System.IO.Path]::GetFullPath($Path)
    if ($fullPath -match '^(?<drive>[A-Za-z]):\\(?<rest>.*)$') {
        $drive = $Matches.drive.ToLowerInvariant()
        $rest = ($Matches.rest -replace '\\', '/')
        if ([string]::IsNullOrWhiteSpace($rest)) {
            return "/mnt/$drive"
        }
        return "/mnt/$drive/$rest"
    }

    throw "Cannot convert path to WSL form: $Path"
}

foreach ($path in @($dataMat, $officialMat, $runScript, $compareScript)) {
    if (-not (Test-Path -LiteralPath $path)) {
        throw "Required file not found: $path"
    }
}

if (-not $DryRun) {
    New-Item -ItemType Directory -Force -Path $outputDir | Out-Null
}

$repoRootWsl = Convert-ToWslPath $repoRoot
$dataMatWsl = Convert-ToWslPath $dataMat
$officialMatWsl = Convert-ToWslPath $officialMat
$runScriptWsl = Convert-ToWslPath $runScript
$compareScriptWsl = Convert-ToWslPath $compareScript
$outputDirWsl = Convert-ToWslPath $outputDir
$comparisonJsonWsl = "$outputDirWsl/comparison_vs_official.json"

$runCmd = @(
    [string]$PythonExe,
    $runScriptWsl,
    "--data-mat", $dataMatWsl,
    "--output-dir", $outputDirWsl,
    "--active-solver-file", "$repoRootWsl/python_scripts/autodl/solver_resdmd_batch3.py",
    "--residual-form", "projected_vlambda",
    "--seed", [string]$Seed,
    "--outer-epochs", [string]$OuterEpochs,
    "--inner-epochs", [string]$InnerEpochs,
    "--batch-size", [string]$BatchSize,
    "--lr", [string]$LearningRate,
    "--reg", [string]$Reg,
    "--training-policy", $TrainingPolicy,
    "--analysis-dtype", $AnalysisDtype,
    "--gram-dtype", $GramDtype,
    "--spectral-dtype", $SpectralDtype
)

$compareCmd = @(
    [string]$PythonExe,
    $compareScriptWsl,
    "--our-mat", "$outputDirWsl/turbulence_active_projected_vlambda_250basis.mat",
    "--official-mat", $officialMatWsl,
    "--output-json", $comparisonJsonWsl,
    "--focus-our-mode", "1",
    "--focus-official-mode", "1",
    "--focus-corr-threshold", [string]$FocusCorrThreshold
)

$runCmdRendered = ($runCmd | ForEach-Object { "'$_'" }) -join " "
$compareCmdRendered = ($compareCmd | ForEach-Object { "'$_'" }) -join " "
$bashCommand = "cd '$repoRootWsl' && $runCmdRendered && $compareCmdRendered"

if ($DryRun) {
    Write-Host "Dry run command:"
    Write-Host "wsl.exe -d $Distro bash -lc `"$bashCommand`""
    return
}

& wsl.exe -d $Distro bash -lc $bashCommand

if ($LASTEXITCODE -ne 0) {
    throw "Turbulence precision regression failed with exit code $LASTEXITCODE"
}

Write-Host "Finished turbulence precision regression."
Write-Host "Output directory: $outputDir"
Write-Host "Comparison JSON: $($outputDir)\comparison_vs_official.json"
