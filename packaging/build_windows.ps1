<#
.SYNOPSIS
    Build the Melanopsin Model one-click Windows executable.

.DESCRIPTION
    Creates (or reuses) an isolated virtual environment, installs the pinned
    build dependencies, and runs PyInstaller against
    packaging/melanopsin_gui.spec. The resulting windowed executable is written
    to dist/MelanopsinModel-v<version>.exe, where <version> comes from
    myutils/_version.py.

.PARAMETER Clean
    Remove the existing build/ and dist/ folders before building.

.PARAMETER SkipVenv
    Build using the current Python environment instead of creating .venv-build.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File packaging\build_windows.ps1 -Clean
#>
[CmdletBinding()]
param(
    [switch]$Clean,
    [switch]$SkipVenv
)

$ErrorActionPreference = "Stop"

# Resolve the repository root (the parent of this script's packaging/ folder).
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$RepoRoot = Split-Path -Parent $ScriptDir
Set-Location $RepoRoot
Write-Host "Repository root: $RepoRoot"

if ($Clean) {
    Write-Host "Cleaning build/ and dist/ ..."
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue (Join-Path $RepoRoot "build")
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue (Join-Path $RepoRoot "dist")
}

$Python = "python"

if (-not $SkipVenv) {
    $VenvDir = Join-Path $RepoRoot ".venv-build"
    if (-not (Test-Path $VenvDir)) {
        Write-Host "Creating build virtual environment at $VenvDir ..."
        & $Python -m venv $VenvDir
    }
    $Python = Join-Path $VenvDir "Scripts\python.exe"
}

Write-Host "Using Python interpreter: $Python"
& $Python --version

Write-Host "Upgrading pip and installing build dependencies ..."
& $Python -m pip install --upgrade pip
& $Python -m pip install -r (Join-Path $ScriptDir "requirements-build.txt")

Write-Host "Running PyInstaller ..."
& $Python -m PyInstaller --noconfirm --clean (Join-Path $ScriptDir "melanopsin_gui.spec")

$Exe = Get-ChildItem -Path (Join-Path $RepoRoot "dist") -Filter "MelanopsinModel-v*.exe" -ErrorAction SilentlyContinue |
    Select-Object -First 1

if ($Exe) {
    Write-Host ""
    Write-Host "Build complete: $($Exe.FullName)"
    Write-Host "Distribute this .exe inside a cloned copy of the repository so it"
    Write-Host "can find the data/ and myutils/ folders next to it."
}
else {
    throw "Build finished but no executable was found in dist/."
}
