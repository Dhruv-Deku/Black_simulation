
# Builds the simulator and runs the physics tests.
#   .\build.ps1            Release build + tests
#   .\build.ps1 -Run       ...then start the simulator
#   .\build.ps1 -Debug     Debug build
param(
    [switch]$Run,
    [switch]$Debug
)
$ErrorActionPreference = 'Stop'
Set-Location $PSScriptRoot

# Find every cmake (PATH and Visual Studio installs) and use the newest; the presets
# need CMake 3.25 or later, and older Visual Studio versions bundle older ones.
$candidates = @()
$onPath = Get-Command cmake -ErrorAction SilentlyContinue
if ($onPath) { $candidates += $onPath.Source }
$vsRoots = @("$env:ProgramFiles\Microsoft Visual Studio", "${env:ProgramFiles(x86)}\Microsoft Visual Studio")
foreach ($root in $vsRoots) {
    $candidates += Get-ChildItem "$root\*\*\Common7\IDE\CommonExtensions\Microsoft\CMake\CMake\bin\cmake.exe" `
        -ErrorAction SilentlyContinue | ForEach-Object FullName
}
$cmake = $null
$best = [version]'0.0'
foreach ($c in $candidates | Select-Object -Unique) {
    if ((& $c --version | Select-Object -First 1) -match '(\d+\.\d+\.\d+)' -and [version]$Matches[1] -gt $best) {
        $best = [version]$Matches[1]
        $cmake = $c
    }
}
if (-not $cmake -or $best -lt [version]'3.25') {
    throw 'CMake 3.25+ not found. Install Visual Studio 2026 (or 2022) with the "Desktop development with C++" workload.'
}
Write-Host "Using CMake $best ($cmake)"

# Use Visual Studio 2026 when it is installed, otherwise 2022 (both build into build\).
$vswhere = "${env:ProgramFiles(x86)}\Microsoft Visual Studio\Installer\vswhere.exe"
$preset = 'vs2026'
if (Test-Path $vswhere) {
    $vs18 = & $vswhere -products * -version '[18.0,19.0)' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    $vs17 = & $vswhere -products * -version '[17.0,18.0)' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $vs18 -and $vs17) { $preset = 'vs2022' }
}
Write-Host "Using preset $preset"

$config = if ($Debug) { 'Debug' } else { 'Release' }
& $cmake --preset $preset
if ($LASTEXITCODE) { exit $LASTEXITCODE }
$buildPreset = $config.ToLower()
if ($preset -ne 'vs2026') { $buildPreset += "-$preset" }
& $cmake --build --preset $buildPreset --parallel
if ($LASTEXITCODE) { exit $LASTEXITCODE }

& ".\build\$config\kerr_tests.exe"
if ($LASTEXITCODE) { Write-Warning 'Physics tests failed.'; exit $LASTEXITCODE }

if ($Run) { & ".\build\$config\blackhole.exe" }
