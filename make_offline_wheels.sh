#!/usr/bin/env bash
# 生成离线安装依赖包：把 sc-audit 的 wheel 与全部依赖下到 ./wheels/
# 之后把 wheels/ 与 install.sh 一起拷走，即可在无网机器上 bash install.sh --offline
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-python3}"
mkdir -p "$HERE/wheels"
"$PY" -m pip download --dest "$HERE/wheels" \
  numpy pandas scipy anndata scanpy matplotlib pyyaml
WHEEL="$(ls "$HERE"/sc_audit-*.whl 2>/dev/null | head -1 || true)"
[ -n "$WHEEL" ] && cp "$WHEEL" "$HERE/wheels/"
echo "✓ wheels/ 就绪：$(ls "$HERE/wheels" | wc -l | tr -d ' ') 个文件"
echo "  分发方法：把 install.sh + wheels/ + (可选)本目录一起拷到离线机，运行 bash install.sh --offline"
