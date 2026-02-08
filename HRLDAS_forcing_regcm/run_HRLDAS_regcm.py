"""
This script is modified from the original work by

    DOAN Quang Van and XUE Lingbo
    Center for Computational Sciences (CCS), University of Tsukuba, Japan

Original purpose:
   Automated downloading and preprocessing of ERA5 data,
   computation of geopotential on model levels, and generation
   of HRLDAS setup and LDASIN forcing files. 

Modifications by:
    Ratchanan Srisawadwong
    Ramkhamkaeng University, Center of Regional Climate Change and Renewable Energy (RU-CORE).
    2026

Main modifications include:

This version is adapted for research and educational use.
"""

#from download_era5_hrldas import *
from create_forcing_prl import *

import os
import shutil
import subprocess


gcm_name = "CMCC-ESM2"
scenario = "ssp245"
start_year = 2020
end_year = 2020
#months = ['01','02']
months = ['01','02','03','04','05','06','07','08','09','10','11','12']
loop_start_date = '01-01'
loop_end_date = '12-31'
area = [14.25, 99.9, 13.45, 101]
levelists = ['136']
dir_raw = '/home/admin/WRF/RUN/LSP-DS/Data_114/run_regcm/RegCM_Data/CMCC-ESM2/ssp245/'
dir_hrldas = '/home/admin/WRF/RUN/LSP-DS/Data_114/run_regcm/hrldas_simulation/'
geo_em_file = '/home/admin/WRF/RUN/LSP-DS/Data_114/run_regcm/geo_em.d01.nc'
#urbanParamTable = '../test/ERA5/tables/URBPARM.TBL'
#nameList = '../test/ERA5/namelists/namelist.hrldas'
#exe_directory = '../hrldas/hrldas/run/'
levelist = '136'
ZLVL = 30
calendar="gregorian"


create_lai_vegfra(geo_em_file, dir_hrldas)
print("LAI : OK")
create_setup_file(f'{str(start_year)}-{loop_start_date}',dir_raw, dir_hrldas,geo_em_file)
print("create_setup_file : OK" )
num_cores = 20

for year in range(start_year, end_year+1):
    create_LDASIN_files_parallel(f'{str(year)}-{loop_start_date}', f'{str(year)}-{loop_end_date}', \
                        dir_raw, dir_hrldas, \
                        geo_em_file, gcm_name,scenario ,levelist, ZLVL, calendar, num_cores)

