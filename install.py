#!/usr/bin/env python3
"""sc-audit 一键安装（跨平台，纯标准库）。

用法:
    python3 install.py            # 联网安装
    python3 install.py --offline  # 用同目录 wheels/ 离线安装
环境变量: SC_AUDIT_VENV（默认 ~/.sc-audit-venv）、PYTHON
"""
import argparse, glob, os, subprocess, sys, venv

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    a = ap.parse_args()
    here = os.path.dirname(os.path.abspath(__file__))
    if sys.version_info[:2] < (3, 9):
        print("x Python 需 >=3.9，当前", sys.version.split()[0]); return 1
    venv_dir = os.environ.get("SC_AUDIT_VENV", os.path.expanduser("~/.sc-audit-venv"))
    print("== sc-audit 一键安装 ==")
    print("   Python:", sys.version.split()[0], "-> venv:", venv_dir)
    _sysites = "--system-site-packages" in os.environ.get("SC_AUDIT_VENV_ARGS", "")
    venv.EnvBuilder(with_pip=True, system_site_packages=_sysites).create(venv_dir)
    vpy = os.path.join(venv_dir, "Scripts" if os.name == "nt" else "bin",
                       "python.exe" if os.name == "nt" else "python")
    subprocess.run([vpy, "-m", "pip", "install", "--upgrade", "pip", "-q"], check=False)
    whls = sorted(glob.glob(os.path.join(here, "sc_audit-*.whl")))
    target = whls[0] if whls else here
    cmd = [vpy, "-m", "pip", "install"]
    cmd += [a for a in os.environ.get("SC_AUDIT_PIP_ARGS", "").split() if a]
    if a.offline:
        wheels = os.path.join(here, "wheels")
        if not os.path.isdir(wheels):
            print("x --offline 需要同目录 wheels/"); return 1
        cmd += ["--no-index", "--find-links", wheels]
    cmd.append(target)
    if subprocess.run(cmd).returncode != 0:
        print("x 安装失败"); return 1
    exe = os.path.join(venv_dir, "Scripts" if os.name == "nt" else "bin",
                       "sc-audit.exe" if os.name == "nt" else "sc-audit")
    subprocess.run([exe, "--version"], check=False)
    print("\nOK 安装完成。用法：", exe, "--help")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
