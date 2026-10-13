#!/usr/bin/env bash
# macOS 双击运行（Finder 里双击本文件即安装）
cd "$(dirname "$0")" && bash install.sh
echo; read -n 1 -s -r -p "按任意键关闭…"
