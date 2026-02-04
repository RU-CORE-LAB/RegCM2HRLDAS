#!/usr/bin/env bash
# ===========================================================================
# RegCM → HRLDAS Surface Forcing Converter
# ----------------------------------------
# Developed by RU-CORE | Author: Ratchanan Srisawadwong (Nick)
#
# Converts RegCM surface (SRF) output into HRLDAS single-layer forcing files.
#
# Main processes:
#   • Spatial subsetting to target domain
#   • Unit conversion:
#       - Radiation: W m⁻² → J m⁻² (hourly accumulation)
#       - Precipitation: kg m⁻² s⁻¹ → m hr⁻¹
#   • Variable renaming to HRLDAS convention
#   • Dimension renaming (time, lat, lon)
#   • Removal of time bounds (if present)
#
# Input variables:
#   ps   → sp    (Surface Pressure)
#   rlds → strd  (Longwave Downward Radiation)
#   rsds → ssrd  (Shortwave Downward Radiation)
#   pr   → tp    (Total Precipitation)
#
# Output:
#   <YYYYMM>_ssrd_strd_sp_tp_<GCM>_<SCENARIO>_single_layer.nc
#
# Requirements: CDO, NCO (ncap2, ncrename, ncatted, ncks, ncwa)
#
# ===========================================================================
set -euo pipefail

usage() {
cat << EOF
Usage:
  regcm2hrldas-srf -m <gcm_name> -p <scenario> -s <start_year> -e <end_year> \
                  -d lon1,lon2,lat1,lat2 -i <input_dir> -o <output_dir>
EOF
exit 1
}

while getopts ":m:p:s:e:d:i:o:h" opt; do
  case $opt in
    m) model="$OPTARG" ;;
    p) sno="$OPTARG" ;;
    s) sYR="$OPTARG" ;;
    e) eYR="$OPTARG" ;;
    d) dom="$OPTARG" ;;
    i) input="$OPTARG" ;;
    o) output="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "${model:-}" || -z "${sno:-}" || -z "${sYR:-}" || -z "${eYR:-}" || -z "${dom:-}" || -z "${input:-}" || -z "${output:-}" ]] && usage

if [[ ! "$dom" =~ ^-?[0-9]+([.][0-9]+)?,-?[0-9]+([.][0-9]+)?,-?[0-9]+([.][0-9]+)?,-?[0-9]+([.][0-9]+)?$ ]]; then
  echo "Error: -d must be lon1,lon2,lat1,lat2"
  exit 1
fi

command -v cdo >/dev/null || { echo "CDO not found"; exit 1; }
command -v ncap2 >/dev/null || { echo "NCO not found"; exit 1; }

mkdir -p "$output"
tmp_srf=$(mktemp -d -t regcm2hrldas_srf_XXXX)
echo "--------------------------------------------------"
echo "Temp dir: $tmp_srf"

mapfile -t ls_file < <(
  find "$input" -maxdepth 1 -name "*_SRF.*.nc" | awk -v s=$sYR -v e=$eYR '
    match($0, /_SRF\.([0-9]{4})/, a) { if (a[1] >= s && a[1] <= e) print }'
)

select_vars="ps,rlds,rsds,pr"

# ================== MAIN LOOP ==================
for fl in "${ls_file[@]}"; do
    fname=$(basename "$fl")
    ymon=$(echo "$fname" | cut -d '.' -f2 | cut -c1-6)

    echo "--------------------------------------------------"
    echo "Processing: $fname"

    # Step 1
    printf "Step 1: Spatial subset ... "
    if cdo -s -O -setgridtype,lonlat -selvar,$select_vars -sellonlatbox,$dom "$fl" "$tmp_srf/$fname"; then
        printf "[OK]\n"
    else
        printf "[FAIL]\n"; continue
    fi

    # Step 2
    printf "Step 2: Unit conversion ... "
    if ncap2 -O -s 'rlds=rlds*3600; rsds=rsds*3600; pr=pr*3.6' "$tmp_srf/$fname" "$tmp_srf/$fname"; then
        ncatted -a units,rlds,o,c,"J m-2" -a units,rsds,o,c,"J m-2" -a units,pr,o,c,"m" "$tmp_srf/$fname"
        printf "[OK]\n"
    else
        printf "[FAIL]\n"; continue
    fi

    # Step 3
    printf "Step 3: Rename variables/dims ... "
    if ncrename -O \
      -d time,valid_time -v time,valid_time \
      -d lon,longitude   -v lon,longitude \
      -d lat,latitude    -v lat,latitude \
      -v pr,tp -v ps,sp -v rlds,strd -v rsds,ssrd \
      "$tmp_srf/$fname" "$tmp_srf/$fname.tmp"; then
        printf "[OK]\n"
    else
        printf "[FAIL]\n"; continue
    fi

    infile="$tmp_srf/$fname.tmp"
    outfile="${output}/${ymon}_ssrd_strd_sp_tp_${model}_${sno}_single_layer.nc"

    # Step 4
    printf "Step 4: Remove time bounds (if any) ... "
    if ! ncks -m "$infile" | grep -q "time_bnds"; then
        cp "$infile" "$outfile"
        printf "[SKIP]\n"
    else
        ncks -O -C -x -v time_bnds "$infile" "$outfile"
        ncks -m "$outfile" | grep -q "bnds =" && ncwa -O -a bnds "$outfile" "$outfile"
        printf "[OK]\n"
    fi

    echo "Done: $outfile"
done

# ------------------ Cleanup ------------------
rm -rf "$tmp_srf"
echo "All SRF processing completed."
