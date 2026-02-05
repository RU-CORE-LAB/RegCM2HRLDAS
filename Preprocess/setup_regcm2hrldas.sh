#!/usr/bin/env bash
#
# ==========================================================================
#
# RegCM → HRLDAS Toolchain Installer
# ----------------------------------
# Developed by RU-CORE | Author: Ratchanan Srisawadwong (Nick)
#
# This script builds and installs the full RegCM-to-HRLDAS preprocessing
# toolchain, including:
#
#   1) Compilation of RegCM sigma conversion utilities
#      • sigma2p_at_hPa
#
#   2) Creation of a dedicated Conda environment
#
#   3) Installation of the Python package:
#      regcm2hrldas_tools (editable/development mode)
#
# Requirements:
#   • Pre-compiled RegCM installation
#   • Conda (Anaconda or Miniconda)
#   • gfortran and NetCDF libraries (for sigma tools)
#
# Usage:
#   Edit configuration file : ./config/regcm2hrldas.cfg
#   then run:
#
#       bash setup_regcm2hrldas.sh
#
# After installation:
#   conda activate regcm2hrldas
#
# ============================================================================
#
# --------------------------------------------------
# 01 compile sigm2p and sigm2z from RegCM installed.
# --------------------------------------------------
source ./config/regcm2hrldas.cfg
curr_home=$PWD



cd build_sigma
bash build_sigma.sh
cd ${curr_home}

# ------------------------------
# Checking Anaconda or miniconda
# ------------------------------
if command -v conda > /dev/null 2>&1; then
	CONDA_EXE=$(which conda 2> /dev/null)
	CONDA_EXE=$(readlink -f "$CONDA_EXE")
	CONDA_BASE=$(dirname "$(dirname "$CONDA_EXE")")
        conda --version
else
        echo "Not Found 'conda' Please install Anaconda/miniconda."
        exit 1
fi

CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
source ${CONDA_SH}
echo "source ${CONDA_SH}"
conda deactivate
#CONDA_ENV="regcm2hrldas"

if conda env list | awk '{print $1}' | grep -qx "$CONDA_ENV"; then
	echo ${CONDA_ENV} is already.
else 	
	conda create -n ${CONDA_ENV} python=3.10 -y
fi
echo "conda activate ${CONDA_ENV}"
source ${CONDA_SH}
conda activate ${CONDA_ENV}

# ----------------------------------------------
# 02 Install regcm2hrldas_tools python package.
# ----------------------------------------------
cd ${curr_home}
pip install -e .


echo install completed.
echo "# To activate this environment, use"
echo "#"
echo "#     $ conda activate ${CONDA_ENV}"
echo "#"
echo "# To deactivate an active environment, use"
echo "#"
echo "#     $ conda deactivate"
