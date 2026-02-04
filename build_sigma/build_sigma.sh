#!/bin/bash
set -e

ROOT_DIR=$(cd "$(dirname "$0")/.." && pwd)
CONFIG_FILE="$ROOT_DIR/config/regcm2hrldas.cfg"
BIN_DIR="$ROOT_DIR/srctools/bin"  

mkdir -p "$BIN_DIR"

echo "=============================="
echo " Building RegCM Sigma Tools"
echo "=============================="

# -----------------------------
# Load config
# -----------------------------
if [ -f "$CONFIG_FILE" ]; then
    echo "Loading config: $CONFIG_FILE"
    source "$CONFIG_FILE"
else
    echo "Config file not found!"
    exit 1
fi

FC=${FC:-gfortran}
USE_OPENMP=${USE_OPENMP:-1}

SHARE_DIR="$REGCM_DIR/Share"
EXT_DIR="$REGCM_DIR/external"

if [ ! -f "$SHARE_DIR/librcmlib.a" ]; then
    echo "librcmlib.a not found in $SHARE_DIR"
    exit 1
fi

# -----------------------------
# NetCDF Fortran
# -----------------------------
if ! command -v nf-config &> /dev/null; then
    echo "nf-config not found. Load netcdf-fortran first."
    exit 1
fi

NETCDF_INC=$(nf-config --includedir)
NETCDF_FLIBS=$(nf-config --flibs)

# -----------------------------
# Compiler flags
# -----------------------------
if [ "$FC" = "ifort" ]; then
    OMPFLAG=""
    [ "$USE_OPENMP" = "1" ] && OMPFLAG="-qopenmp"
    FFLAGS="-O3 $OMPFLAG -convert big_endian"
else
    OMPFLAG=""
    [ "$USE_OPENMP" = "1" ] && OMPFLAG="-fopenmp"
    FFLAGS="-O3 $OMPFLAG -fconvert=big-endian -fno-range-check"
fi

INCLUDES="-I$NETCDF_INC -I$SHARE_DIR -I$EXT_DIR"
LIBS="$NETCDF_FLIBS -L$SHARE_DIR -lrcmlib -lhdf5_hl -lhdf5 -lz -lm"

echo "Compiler : $FC"
echo "RegCM dir: $REGCM_DIR"
echo "Building sigma2p_at_hPa..."

cd "$ROOT_DIR/build_sigma"
make clean
make FC="$FC" FFLAGS="$FFLAGS" INCLUDES="$INCLUDES" LIBS="$LIBS"

mv sigma2p_at_hPa "$BIN_DIR/"

echo "Build complete!"
echo "Binary installed to: $BIN_DIR/sigma2p_at_hPa"

