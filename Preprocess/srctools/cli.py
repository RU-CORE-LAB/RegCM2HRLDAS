import subprocess
import pathlib
import sys


SCRIPT_DIR = pathlib.Path(__file__).parent / "scripts"


def _run_script(script_name):
    script_path = SCRIPT_DIR / script_name

    if not script_path.exists():
        raise FileNotFoundError(f"Script not found: {script_path}")

    cmd = ["bash", str(script_path)] + sys.argv[1:]
    #print("Running:", " ".join(cmd))
    subprocess.check_call(cmd)


def atm_main():
    #print("- REGCM -> HRLDAS ATM")
    _run_script("regcm2hrldas_atm.sh")


def srf_main():
    #print("- REGCM -> HRLDAS SRF")
    _run_script("regcm2hrldas_srf.sh")


def setup_main():
    #print("- REGCM -> HRLDAS SETUP")
    _run_script("regcm2hrldas_setup.sh")

