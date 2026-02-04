import shutil
import subprocess
from pathlib import Path
import sys
import os

def find_sigma2p():
    """Find sigma2p_at_hPa executable"""

    # 1) Check system PATH
    exe = shutil.which("sigma2p_at_hPa")
    if exe:
        return exe

    # 2) Check inside installed package
    pkg_dir = Path(__file__).resolve().parent
    exe_pkg = pkg_dir / "bin" / "sigma2p_at_hPa"

    if exe_pkg.exists():
        return str(exe_pkg)

    raise RuntimeError(
        "sigma2p_at_hPa not found.\n"
        "Run build_sigma/build_sigma.sh first."
    )

def run_sigma2p_entry():
    """Entry point ให้เรียก sigma2p_at_hPa จาก PATH ของ package"""
    script_dir = os.path.dirname(__file__)
    exe_path = os.path.join(script_dir, "bin", "sigma2p_at_hPa")

    if not os.path.exists(exe_path):
        raise FileNotFoundError(f"sigma2p_at_hPa not found at {exe_path}")

    cmd = [exe_path] + sys.argv[1:]
    subprocess.check_call(cmd)
