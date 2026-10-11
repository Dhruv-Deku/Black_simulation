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
    throw 'CMake 3.25+ not found. Install Visual Studio 2026 with the "Desktop development with C++" workload.'
}
Write-Host "Using CMake $best ($cmake)"

$config = if ($Debug) { 'Debug' } else { 'Release' }
& $cmake --preset vs2026
if ($LASTEXITCODE) { exit $LASTEXITCODE }
& $cmake --build --preset ($config.ToLower()) --parallel
if ($LASTEXITCODE) { exit $LASTEXITCODE }

& ".\build\$config\kerr_tests.exe"
if ($LASTEXITCODE) { Write-Warning 'Physics tests failed.'; exit $LASTEXITCODE }

if ($Run) { & ".\build\$config\blackhole.exe" }
