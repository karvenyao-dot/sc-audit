#!/usr/bin/env bash
# sc-audit 演示：跑一遍四项审计 + v0.3 可识别性层，并生成 HTML 报告。
#
# 用法:
#   bash examples/run_demo.sh                          # 用**合成**演示数据（自包含，推荐）
#   bash examples/run_demo.sh <input.h5ad> [ref.csv] [outdir]
#
# 说明：本脚本不依赖任何外部/私有数据集；未给输入时会调用 examples/make_demo_data.py 生成合成数据。
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$HERE")"
INPUT="${1:-}"
REF="${2:-}"
OUT="${3:-./sc_audit_demo_out}"

if [ -z "$INPUT" ]; then
  if [ -n "${SC_AUDIT_DEMO_H5AD:-}" ] && [ -f "${SC_AUDIT_DEMO_H5AD}" ]; then
    INPUT="${SC_AUDIT_DEMO_H5AD}"          # 也可以用环境变量指向自己的数据
  else
    INPUT="$HERE/demo/demo.h5ad"
    [ -f "$INPUT" ] || python3 "$HERE/make_demo_data.py" "$INPUT"
  fi
fi

ARGS=(--input "$INPUT" --condition-col condition --donor-col donor
      --target-genes MYL3,PLN,TNNT2,ACTB,GAPDH
      --outdir "$OUT")
[ -n "$REF" ] && ARGS+=(--reference "$REF")

echo "== 预检（只报告，不采纳）=="
python3 -m sc_audit.cli check "${ARGS[@]}" || true

echo
echo "== 四项审计 + 可识别性层 =="
python3 -m sc_audit.cli run "${ARGS[@]}"

echo
echo "产物目录：$OUT"
ls "$OUT"
echo "打开 $OUT/report.html 查看可读报告。"
