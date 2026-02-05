#!/usr/bin/env bash
#
# ===========================================================================
# RegCM → HRLDAS Soil Setup File Generator
# ----------------------------------------
# Developed by RU-CORE | Author: Ratchanan Srisawadwong (Nick)
#
# Extracts land surface initial state variables from RegCM SRF output
# and converts them into an HRLDAS-compatible soil setup file.
#
# Main processes:
#   • Spatial subsetting to target domain
#   • Extraction of soil & surface state variables
#   • Conversion of RegCM soil layers to ERA5/HRLDAS soil layer depths
#     using thickness-weighted averaging
#
# Input variables:
#   ts    : Surface Skin Temperature
#   tsoil : Soil Temperature
#   mrsos : Soil Moisture
#   snw   : Snow Water Equivalent
#
# Output:
#   <YYYYMM>0100_setup.nc
#
# Requirements:
#   - CDO
#   - regcm2setupfile (custom RU-CORE soil layer conversion tool)
#
# ===========================================================================

set -euo pipefail

usage() {
cat << EOF
Usage:
  regcm2hrldas-setup -m <gcm_name> -p <scenario> -s <year> \
                     -d lon1,lon2,lat1,lat2 -i <input_dir> -o <output_dir>
EOF
exit 1
}

while getopts ":m:p:d:s:i:o:h" opt; do
  case $opt in
    m) model="$OPTARG" ;;
    p) sno="$OPTARG" ;;
    s) sYM="$OPTARG" ;;
    d) dom="$OPTARG" ;;
    i) input="$OPTARG" ;;
    o) output="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "${model:-}" ]] && echo "Missing -m GCM Name" && exit 1
[[ -z "${sno:-}" ]] && echo "Missing -p Scenario" && exit 1
[[ -z "${sYM:-}" ]] && echo "Missing -s start year" && exit 1
[[ -z "${input:-}" ]] && echo "Missing -i input path" && exit 1
[[ -z "${output:-}" ]] && echo "Missing -o output path" && exit 1
[[ -z "${dom:-}" ]] && echo "Missing -d lon1,lon2,lat1,lat2" && exit 1

if [[ ! "$dom" =~ ^-?[0-9]+([.][0-9]+)?,-?[0-9]+([.][0-9]+)?,-?[0-9]+([.][0-9]+)?,-?[0-9]+([.][0-9]+)?$ ]]; then
  echo "Error: -d must be lon1,lon2,lat1,lat2"
  exit 1
fi

command -v cdo >/dev/null || { echo "CDO not found"; exit 1; }
command -v regcm2setupfile >/dev/null || { echo "regcm2setupfile not found"; exit 1; }


IFS=',' read -r lo1 lo2 la1 la2 <<< "$dom"

echo "Looking for SRF file: *_SRF.${sYM}*.nc"

fl=$(ls "$input"/*_SRF.${sYM}*.nc 2>/dev/null | head -n 1 || true)

if [[ -z "$fl" ]]; then
  echo "No SRF file found for ${sYM} in ${input}"
  exit 1
fi

select_vars="ts,tsoil,mrsos,snw"

mkdir -p "$output"

tmp_srf=$(mktemp -d "${TMPDIR:-/tmp}/regcm2hrldas_setup_XXXX")
trap 'rm -rf "$tmp_srf"' EXIT

echo "--------------------------------------------------"
echo "Temp dir: $tmp_srf"

fname=$(basename "$fl")
ymon=$(echo "$fname" | cut -d '.' -f2 | cut -c1-6)

echo "--------------------------------------------------"
echo "Processing : $fname"

printf "Step 1: Spatial subset & vars (%s) ... " "$select_vars"
if cdo -s -O -setgridtype,lonlat \
    -seltimestep,1 \
    -selvar,$select_vars \
    -sellonlatbox,$lo1,$lo2,$la1,$la2 \
    "$fl" "$tmp_srf/$fname"; then
    printf "[OK]\n"
else
    printf "[FAIL]\n"
    exit 1
fi

printf "Step 2: Soil layer to ERA5 levels ... \n"
if regcm2setupfile "$tmp_srf/$fname" "$output/${ymon}0100_setup.nc"; then
    printf "[OK]\n"
else
    printf "[FAIL]\n"
    exit 1
fi

echo "Setup file created: ${ymon}0100_setup.nc"

