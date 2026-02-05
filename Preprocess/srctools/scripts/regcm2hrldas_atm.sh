#!/usr/bin/env bash
#
# ======================================================================================
# RegCM → HRLDAS Atmospheric Forcing Converter
# --------------------------------------------
# Developed by RU-CORE | Author: Ratchanan Srisawadwong (Nick)
#
# Converts 6-hourly RegCM atmospheric output into hourly HRLDAS pressure-level forcing.
#
# Main processes:
#   • Spatial subsetting to target domain
#   • Sigma → pressure level conversion
#   • Variable extraction (T, U, V, Q, HGT, TOPO)
#   • Temporal interpolation (6-hourly → hourly)
#   • Geopotential height calculation (z = hgt + topo)
#   • NetCDF reformatting for HRLDAS
#
# Output:
#   <YYYYMM>_t_u_v_q_z_<GCM>_<SCENARIO>_pressure_level.nc
#
# Requirements: CDO, NCO, interp1h, sigma2p_at_hPa
#
# =======================================================================================
set -euo pipefail

usage() {
cat << EOF
Usage:
  regcm2hrldas-atm -m <gcm_name> -p <scenario> -s <start_year> -e <end_year> \
                  -d lon1,lon2,lat1,lat2 -i <input_dir> -o <output_dir>

Example:
  regcm2hrldas-atm -m MPI -p ssp245 -s 2030 -e 2035 \
                  -d 96,107,4,22 -i ./regcm_output -o ./hrldas_forcing
EOF
exit 1
}

# ------------------ Parse arguments ------------------
while getopts ":m:p:s:e:d:i:o:h" opt; do
  case $opt in
    m) model="$OPTARG" ;;
    p) scenario="$OPTARG" ;;
    s) sYR="$OPTARG" ;;
    e) eYR="$OPTARG" ;;
    d) dom="$OPTARG" ;;
    i) input="$OPTARG" ;;
    o) output="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

[[ -z "${model:-}" || -z "${scenario:-}" || -z "${sYR:-}" || -z "${eYR:-}" || -z "${dom:-}" || -z "${input:-}" || -z "${output:-}" ]] && usage

# ------------------ Validate domain ------------------
if [[ ! "$dom" =~ ^-?[0-9]+([.][0-9]+)?,-?[0-9]+([.][0-9]+)?,-?[0-9]+([.][0-9]+)?,-?[0-9]+([.][0-9]+)?$ ]]; then
  echo "Error: -d must be lon1,lon2,lat1,lat2"
  exit 1
fi

# ------------------ Dependencies ------------------
command -v cdo >/dev/null || { echo "CDO not found"; exit 1; }
command -v ncrename >/dev/null || { echo "NCO not found"; exit 1; }
command -v interp1h >/dev/null || { echo "interp1h not found"; exit 1; }

sigma2hPa=$(command -v sigma2p_at_hPa) || { echo "sigma2p_at_hPa not found"; exit 1; }

# ------------------ Prepare dirs ------------------
mkdir -p "$output"
tmp_atm=$(mktemp -d -t regcm2hrldas_atm_XXXX)
echo "--------------------------------------------------"
echo "Temp dir: $tmp_atm"

# ------------------ Find input files ------------------
mapfile -t files < <(
  find "$input" -maxdepth 1 -name "*_ATM.*.nc" | awk -v s=$sYR -v e=$eYR '
  match($0, /_ATM\.([0-9]{4})/, a) { if (a[1] >= s && a[1] <= e) print }'
)

[[ ${#files[@]} -eq 0 ]] && { echo "No files found"; exit 1; }

# ================== MAIN LOOP ==================
for fl in "${files[@]}"; do
  	fname=`basename ${fl}`
        ymon=`echo ${fname} | cut -d '.' -f2 | cut -c1-6`
        fhPa=`echo ${fname} | cut -d '.' -f1-2`"_pressure.nc"
        echo "--------------------------------------------------"
        echo "Processing : $fname"

        # *-------------------------------------------------
        printf "Step 1: Spatial subset ... "

        if cdo -s -O -sellonlatbox,${dom} "$fl" "${tmp_atm}/${fname}"; then
                printf "[OK]\n"
        else
                printf "[SKIP]\n"
                continue
        fi
        # *-------------------------------------------------
        printf "Step 2: sigma2p ... "

        if ${sigma2hPa} "${tmp_atm}/${fname}" "${tmp_atm}"; then
                printf "[OK]\n"
                rm ${tmp_atm}/${fname}
        else
                printf "[SKIP]\n"
                continue
        fi

        # *------------------------------------------------
        printf "Step 3: Select : ta ua va hus hgt topo ..."

        cdo -s -O -setgridtype,lonlat -sellevel,1010.8487 -selvar,ta,ua,va,hus,hgt ${tmp_atm}/${fhPa} ${tmp_atm}/ta_ua_va_hus_hgt.nc
        cdo -s -O -setgridtype,lonlat -selvar,topo ${tmp_atm}/${fhPa} ${tmp_atm}/topo.nc

        # ---------- topo 2d -> topo 4d
        cdo -s -O -selvar,hgt ${tmp_atm}/ta_ua_va_hus_hgt.nc ${tmp_atm}/hgt.nc
        cdo -s -O merge ${tmp_atm}/topo.nc ${tmp_atm}/hgt.nc ${tmp_atm}/hgt_topo.nc

        #
        # hgt(time, plev, lat, lon) = 0
        # topo(time, plev, lat, lon) = 0
        # topo(time, plev, lat, lon) = topo(1)
        cdo -w -s -O add -chname,hgt,topo -mulc,0 -selname,hgt ${tmp_atm}/hgt_topo.nc -selname,topo ${tmp_atm}/hgt_topo.nc  ${tmp_atm}/topo_4d.nc


        printf " [OK]\n"

        printf "Step 5: merge and formating ${tmp_atm}/topo_4d.nc  ..."
        #
        # t u v q merge z
        #
        cdo -s -w -O merge ${tmp_atm}/ta_ua_va_hus_hgt.nc ${tmp_atm}/topo_4d.nc ${tmp_atm}/ta_ua_va_hus_hgt_topo.nc
        ncrename -O \
        -d time,valid_time   -v time,valid_time \
        -d lon,longitude     -v lon,longitude \
        -d lat,latitude      -v lat,latitude \
        -d plev,model_level  -v plev,model_level \
        -v ta,t \
        -v ua,u \
        -v va,v \
        -v hus,q\
        ${tmp_atm}/ta_ua_va_hus_hgt_topo.nc ${tmp_atm}/t_u_v_q_hgt_topo.nc

        ncap2 -O -s 'model_level[$model_level]=136.0' ${tmp_atm}/t_u_v_q_hgt_topo.nc ${tmp_atm}/t_u_v_q_hgt_topo.nc

        printf " [OK]\n"

        printf "Step 6: Temporal interpolation (6hr) -> 1(hr) ... \n"
	# ----------------------
	# before interpolation : 
	# ----------------------
	curr_day(){ cdo -s showtimestamp -seltimestep,1 "$1" | head -n1 | awk -F'T' '{print $1}' | tr -d '-' | tr -d '[:space:]'; }
	prev_day(){ cdo -s showtimestamp -shifttime,-1day -seltimestep,1 "$1" | head -n1 | awk -F'T' '{print $1}' | tr -d '-' | tr -d '[:space:]'; }
	next_day(){ cdo -s showtimestamp -seltimestep,-1 "$1" | head -n1 | awk -F'T' '{print $1}' | tr -d '-' | tr -d '[:space:]'; }

	
	mkdir -p .sav_t_step
	# -----------------------
	# Check savefile prev_day.
	# -----------------------
	nt_day=$(next_day "${tmp_atm}/t_u_v_q_hgt_topo.nc")
	pv_day=$(prev_day "${tmp_atm}/t_u_v_q_hgt_topo.nc")
	cr_day=$(curr_day "${tmp_atm}/t_u_v_q_hgt_topo.nc")
	used_sav=0
	f_time_of_day=".sav_t_step/${model}_${scenario}_${cr_day}00.tmp"
	if [[ -f $f_time_of_day ]]; then
		cdo -w -s -O mergetime ${f_time_of_day} ${tmp_atm}/t_u_v_q_hgt_topo.nc ${tmp_atm}/t_u_v_q_hgt_topo_add_frist_t_step.nc 
		rm ${f_time_of_day}
		used_sav=1
	else
		mv ${tmp_atm}/t_u_v_q_hgt_topo.nc ${tmp_atm}/t_u_v_q_hgt_topo_add_frist_t_step.nc
	fi
	# --------------------------------
	# keep last timestep for next file
	# --------------------------------
	l_time_of_day=".sav_t_step/${model}_${scenario}_${nt_day}00.tmp"
	echo ${l_time_of_day}
	cdo -s -O -seltimestep,-1 ${tmp_atm}/t_u_v_q_hgt_topo_add_frist_t_step.nc ${l_time_of_day}
	
	# ----------------------
        # interpolarion : interp1h
	# ----------------------
	tmp_outfile="${tmp_atm}/${ymon}_t_u_v_q_hgt_topo.nc"
        interp1h ${tmp_atm}/t_u_v_q_hgt_topo_add_frist_t_step.nc ${tmp_outfile}.tmp

	# ----------------------
	# After interp1h remove frist timestep when merge last timestep of prev file.
	# ----------------------
	#
	if [[ ${used_sav} -eq 1 ]];then
		cdo -w -s -O delete,timestep=1 ${tmp_outfile}.tmp ${tmp_outfile}
	else
		mv ${tmp_outfile}.tmp ${tmp_outfile}
	fi


        printf "Step 7: Calculate Geopontential(z): hgt + topo .."
        # z = hgt + topo
        cdo -w -s -O chname,hgt,z -add -selname,hgt ${tmp_outfile} ${tmp_atm}/topo_4d.nc ${tmp_atm}/z_msl.nc
        cdo -w -s -O setattribute,z@standard_name=geopotential_height,z@long_name="Geopotential Height",z@units="m" ${tmp_atm}/z_msl.nc ${tmp_atm}/z_msl_fixed.nc

        outfile="${output}/${ymon}_t_u_v_q_z_${model}_${scenario}_pressure_level.nc"
        cdo -w -s -O merge ${tmp_outfile} ${tmp_atm}/z_msl_fixed.nc ${outfile}
        printf "[OK]\n"
        printf "\tDone. ${outfile}\n"

        rm -rf ${tmp_atm}/*
done

# ------------------ Cleanup ------------------
rm -rf "$tmp_atm"
echo "All ATM processing completed."

