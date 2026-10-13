# sc-audit 一键安装（Windows PowerShell）
#   .\install.ps1            联网安装
#   .\install.ps1 -Offline   用同目录 wheels\ 离线安装
param([switch]$Offline)
$ErrorActionPreference = "Stop"
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
$Py = if ($env:PYTHON) { $env:PYTHON } else { "python" }
& $Py -c "import sys; sys.exit(0 if sys.version_info[:2] >= (3,9) else 1)"
if ($LASTEXITCODE -ne 0) { Write-Host "x Python 需 >=3.9"; exit 1 }
$venv = if ($env:SC_AUDIT_VENV) { $env:SC_AUDIT_VENV } else { "$HOME\.sc-audit-venv" }
Write-Host "== sc-audit 一键安装 =="
Write-Host "   虚拟环境: $venv"
& $Py -m venv $venv
$vpy = "$venv\Scripts\python.exe"
& $vpy -m pip install --upgrade pip -q
$wheel = Get-ChildItem "$Here\sc_audit-*.whl" -ErrorAction SilentlyContinue | Select-Object -First 1
$target = if ($wheel) { $wheel.FullName } else { $Here }
if ($Offline) {
  if (-not (Test-Path "$Here\wheels")) { Write-Host "x -Offline 需要同目录 wheels\"; exit 1 }
  & $vpy -m pip install --no-index --find-links "$Here\wheels" $target
} else {
  & $vpy -m pip install $target
}
Write-Host "-- 验证"
& "$venv\Scripts\sc-audit.exe" --version
Write-Host "OK 安装完成。用法： $venv\Scripts\sc-audit.exe --help"
