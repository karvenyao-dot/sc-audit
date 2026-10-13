#!/usr/bin/env bash
# sc-audit 一键安装（macOS / Linux）
#   bash install.sh            联网安装（依赖从 PyPI 取）
#   bash install.sh --offline  用同目录 wheels/ 离线安装
# 可用环境变量：PYTHON / SC_AUDIT_VENV / SC_AUDIT_VENV_ARGS / SC_AUDIT_PIP_ARGS
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OFFLINE=0; [ "${1:-}" = "--offline" ] && OFFLINE=1

PY="${PYTHON:-}"
if [ -z "$PY" ]; then
  for c in python3 python; do command -v "$c" >/dev/null 2>&1 && { PY="$c"; break; }; done
fi
[ -n "$PY" ] || { echo "✗ 未找到 Python（需 ≥3.9）"; exit 1; }
"$PY" -c 'import sys; sys.exit(0 if sys.version_info[:2] >= (3,9) else 1)' \
  || { echo "✗ Python 需 ≥3.9，当前：$("$PY" -V 2>&1)"; exit 1; }

VENV="${SC_AUDIT_VENV:-$HOME/.sc-audit-venv}"
echo "== sc-audit 一键安装 =="
echo "   Python : $("$PY" -V 2>&1)"
echo "   虚拟环境: $VENV"

"$PY" -m venv ${SC_AUDIT_VENV_ARGS:-} "$VENV"
PIP=("$VENV/bin/python" -m pip)
"${PIP[@]}" install --upgrade pip -q >/dev/null 2>&1 || true

WHEEL="$(ls "$HERE"/sc_audit-*.whl 2>/dev/null | head -1 || true)"
PIPARGS=${SC_AUDIT_PIP_ARGS:-}
if [ "$OFFLINE" = "1" ]; then
  [ -d "$HERE/wheels" ] || { echo "✗ --offline 需要同目录 wheels/（离线依赖包）"; exit 1; }
  TARGET="$HERE"; [ -n "$WHEEL" ] && TARGET="$WHEEL"
  "${PIP[@]}" install --no-index --find-links "$HERE/wheels" $PIPARGS "$TARGET"
else
  TARGET="$HERE"; [ -n "$WHEEL" ] && TARGET="$WHEEL"
  "${PIP[@]}" install $PIPARGS "$TARGET"
fi

echo "-- 验证"
"$VENV/bin/sc-audit" --version
"$VENV/bin/sc-audit" check --help >/dev/null && echo "✓ sc-audit 可用"
echo
echo "✓ 安装完成。用法："
echo "   source \"$VENV/bin/activate\" && sc-audit --help"
echo "   或直接   \"$VENV/bin/sc-audit\" --help"
