##############################################################################
# History:
#   2023.03.31  Original version created by
#               DOAN Quang Van and XUE Lingbo (CCS, Tsukuba, Japan)
#
#   2026.02.03  Modified and extended by Ratchanan Srisawadwong
#
#   Major updates in this version:
#     - Added support for GCM/Scenario-based forcing (RegCM Output datasets)
#     - Added parallel processing for LDASIN generation (multiprocessing)
#     - Implemented calendar-safe time handling using cftime
#     - Unified dataset reader (nco_dataset) for mixed calendars and domains
#     - Replaced external geopotential file with height field from pressure-level data
#     - Added geopotential height (ZGH) as forcing variable
#     - Improved latitude ordering and spatial subsetting logic
#     - Generalized filenames and workflow for multi-year batch processing
##############################################################################
import numpy as np
import os
import pandas as pd
import xarray as xr
#import mpl_toolkits.basemap
#from mpl_toolkits.basemap import interp
from scipy.interpolate import RectBivariateSpline
from multiprocessing import Pool
from functools import partial

# N: add cftime
import cftime
# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
# N: Method for open dataset that different calendar.
def nco_dataset(filename,geo_lat,geo_lon,time=None):
    raw_data_file = xr.open_dataset(filename, engine="netcdf4", decode_times=True, use_cftime=True )
    if raw_data_file.latitude[0] > raw_data_file.latitude[-1]:
        raw_data_file = raw_data_file.sortby('latitude')    
    lat_min = float(geo_lat.min())
    lat_max = float(geo_lat.max())
    lon_min = float(geo_lon.min())
    lon_max = float(geo_lon.max())    
    data_var_raw = raw_data_file.sel(latitude=slice(lat_min, lat_max),longitude=slice(lon_min, lon_max))
    raw_lat = data_var_raw.latitude.values
    raw_lon = data_var_raw.longitude.values
    
    # --- Match calendar only if time is provided ---
    if time is not None:

        # Get time coordinate (time / valid_time / Times )
        time_coord_name = [c for c in data_var_raw.coords if 'time' in c.lower()][0]
        sample_time = data_var_raw[time_coord_name].values[0]
        CFTimeType = type(sample_time)
        cal_name = CFTimeType.__name__

        y, m, d, h = time.year, time.month, time.day, time.hour

        # --- NOLEAP ---
        if cal_name == "DatetimeNoLeap" and m == 2 and d == 29:
            d = 28

        # --- 360_DAY ---
        if cal_name == "Datetime360Day" and d > 30:
            d = 30

        time = CFTimeType(y, m, d, h)

    
    return data_var_raw, raw_lat, raw_lon, time

def floor_to_day_cf(t):
    return type(t)(t.year, t.month, t.day)

# *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*

def process_single_hour(time, gcm_name, scenario,  variables, geo_em_file, raw_data_dir, output_dir, levelist, ZLVL):
    import warnings
    warnings.filterwarnings("ignore")
    
    date = floor_to_day_cf(time)  # date = time.floor('D')
    geo_em = xr.open_dataset(geo_em_file)
    geo_lat, geo_lon = geo_em.XLAT_M.values[0], geo_em.XLONG_M.values[0]
    geo_lat_flat, geo_lon_flat = geo_lat.ravel(), geo_lon.ravel()
    # N: remove z file and use z in regcm file. 
    # z_file = xr.open_dataset(z_file_path)
    # ZGH = variables['z']


    LDASIN_file = xr.Dataset()

    for var in variables:
        if variables[var]['name'] in ['LWDOWN', 'SWDOWN','RAINRATE',"PSFC"]: 
            # N:
            # rename support GCM and scenario names. 
            filename = os.path.join(raw_data_dir, f"{date.strftime('%Y%m')}_ssrd_strd_sp_tp_{gcm_name}_{scenario}_single_layer.nc")
            # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
            # call def nco_dataset()
            data_vars_raw, raw_lat, raw_lon, time = nco_dataset(filename,geo_lat,geo_lon,time) 
            # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
            val = data_vars_raw[var].sel(valid_time=time, method='nearest')[::-1].values
            if variables[var]['name'] in ['LWDOWN', 'SWDOWN']:
                data_var = val / 3600
            elif variables[var]['name'] in ['RAINRATE']:
                data_var = val / 3600 * 1000
            else:
                data_var = val     

        elif variables[var]['name'] in ['T2D']:
            # N:
            # rename support GCM and scenario names. 
            filename = os.path.join(raw_data_dir, f"{date.strftime('%Y%m')}_t_u_v_q_z_{levelist}_{gcm_name}_{scenario}_pressure_layer.nc")
            # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
            # call def nco_dataset()
            data_vars_raw, raw_lat, raw_lon, time = nco_dataset(filename,geo_lat,geo_lon, time) 
            # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
            t_raw = data_vars_raw[var].sel(valid_time=time, model_level=levelist, method='nearest').values[::-1, :]
            #
            # Use 'z' in pressure_layer file.
            # 'z' is already Geopotential Height (meters), not geopotential (m^2 s^-2),
            # so we must NOT divide by gravity (9.80665) again. 
            z_raw = data_vars_raw['z'].sel(valid_time=time, model_level=levelist, method='nearest').values[::-1, :]
            data_var_correct_to_msl = t_raw - (-0.0065 * z_raw )    # lapse rate ~ 6.5 K/km    
            # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*         

        elif variables[var]['name'] in ['LAI', 'VEGFRA']:
            if date.is_leap_year:
                raw_data_file = xr.open_dataset(os.path.join(output_dir,'LDASIN', f'{var}_leap.nc'))
                data_var = [raw_data_file[var].sel(date='2020'+str(date.date())[-6:]).values]
            else:
                raw_data_file = xr.open_dataset(os.path.join(output_dir,'LDASIN', f'{var}.nc'))
                data_var = [raw_data_file[var].sel(date='2021'+str(date.date())[-6:]).values]

        else:
            # N:
            # rename support GCM and scenario names. 
            filename = os.path.join(raw_data_dir, f"{date.strftime('%Y%m')}_t_u_v_q_z_{levelist}_{gcm_name}_{scenario}_pressure_layer.nc")
            # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
            # call def nco_dataset()
            data_vars_raw, raw_lat, raw_lon, time = nco_dataset(filename,geo_lat,geo_lon, time) 
            # *=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
            data_var = data_vars_raw[var].sel(valid_time=time, model_level=levelist, method='nearest')[::-1].values

        if variables[var]['name'] in ['T2D']:
            # data_var_interpolated = mpl_toolkits.basemap.interp(data_var_correct_to_msl, 
            #                                                 raw_lon, raw_lat, 
            #                                                 geo_lon, geo_lat, 
            #                                                 checkbounds=False, masked=False, order=1)
            interp_spline = RectBivariateSpline(raw_lat, raw_lon, data_var_correct_to_msl, kx=1, ky=1)
            data_var_interpolated = interp_spline.ev(geo_lat_flat, geo_lon_flat).reshape(geo_lat.shape)
            data_var_correct_to_HGT_M = data_var_interpolated + ( -0.0065 * (geo_em['HGT_M'].values.squeeze()+ZLVL))
            LDASIN_file[variables[var]['name']] = (('Time','south_north','west_east'), [data_var_correct_to_HGT_M])
            LDASIN_file[variables[var]['name']].attrs['units'] = variables[var]['attrs']['units']
        elif variables[var]['name'] in ['LAI', 'VEGFRA']:
            LDASIN_file[variables[var]['name']] = (('Time','south_north','west_east'), data_var)
            LDASIN_file[variables[var]['name']].attrs['units'] = variables[var]['attrs']['units']
        else:
            # data_var_interpolated = mpl_toolkits.basemap.interp(data_var, 
            #                                                 raw_lon, raw_lat, 
            #                                                 geo_lon, geo_lat, 
            #                                                 checkbounds=False, masked=False, order=1)
            interp_spline = RectBivariateSpline(raw_lat, raw_lon, data_var, kx=1, ky=1)
            data_var_interpolated = interp_spline.ev(geo_lat_flat, geo_lon_flat).reshape(geo_lat.shape)
            LDASIN_file[variables[var]['name']] = (('Time','south_north','west_east'), [data_var_interpolated])
            LDASIN_file[variables[var]['name']].attrs['units'] = variables[var]['attrs']['units']


    encoding=[{var: {'_FillValue': None}} for var in LDASIN_file.variables]    
    output_filename = f"{time.strftime('%Y%m%d%H')}.LDASIN_DOMAIN{geo_em_file[-4]}"
    LDASIN_file.to_netcdf(os.path.join(output_dir, 'LDASIN', output_filename), encoding=encoding[0])
    print(output_filename)

def create_LDASIN_files_parallel(start_date, end_date, raw_data_dir, output_dir, geo_em_file, gcm_name, scenario, levelist, ZLVL, calendar="gregorian", num_cores=1):

    if not os.path.exists(output_dir+"/LDASIN"):
        os.makedirs(output_dir+"/LDASIN")

    variables = {
        't': {'name':'T2D', 'attrs':{'units':'K'}}, 
        'q': {'name':'Q2D', 'attrs':{'units':'kg/kg'}},
        'u': {'name':'U2D', 'attrs':{'units':'m/s'}},
        'v': {'name':'V2D', 'attrs':{'units':'m/s'}},
        'z': {'name':'ZGH', 'attrs': {'units' : 'm' }},  # N: add z in pressure layers file.
        'sp': {'name':'PSFC', 'attrs':{'units':'Pa'}},
        'strd': {'name':'LWDOWN', 'attrs':{'units':'W/m^2'}},
        'ssrd': {'name':'SWDOWN', 'attrs':{'units':'W/m^2'}},
        'tp': {'name':'RAINRATE', 'attrs':{'units':'kg/m^2/s'}},
        'LAI12M':{'name':'LAI', 'attrs':{'units':'m^2/m^2'}},
        'GREENFRAC':{'name':'VEGFRA', 'attrs':{'units':'%'}},
    }
    # N: remove z file and use z in RegCM file. 
    # z_file_path = os.path.join(raw_data_dir, 'z_out.grib')
    #
   
    time_list = [t for date in pd.date_range(start_date, end_date, freq='d')
                    for t in pd.date_range(date, periods=24, freq='h')]
    
   
    process_func = partial(
        process_single_hour,
        variables=variables,
        geo_em_file=geo_em_file,
        raw_data_dir=raw_data_dir,
        output_dir=output_dir,
        gcm_name=gcm_name,
        scenario=scenario,
        levelist=levelist,
        ZLVL=ZLVL
    )

    with Pool(processes=num_cores) as pool:
        pool.map(process_func, time_list)

def create_setup_file(start_date, raw_data_dir, output_dir, geo_em_file):

    if not os.path.exists(output_dir+"/LDASIN"):
        os.makedirs(output_dir+"/LDASIN")

    variables = {
    
        "Times": {'units':''} ,

        # from geo_em file
        "XLAT": {'units': 'degree_north', 'geoname': 'XLAT_M'} , 
        "XLONG": {'units': 'degree_east', 'geoname': 'XLONG_M'} , 
        "HGT": {'units': 'm', 'geoname': 'HGT_M'} , 
        "MAPFAC_MX": {'units': '', 'geoname': 'MAPFAC_MX'} , 
        "MAPFAC_MY": {'units': '', 'geoname': 'MAPFAC_MY'} ,  
        "URBLANDUSEF": {'units': '', 'geoname': 'LANDUSEF'} ,      # for 2D urban fraction

        # edit from geo_em file
        "TMN": {'units': 'K', 'geoname': 'SOILTEMP'} ,  # adjust to elevation
        "SHDMAX": {'units': '%'} ,  # max(100*GREENFRAC) 
        "SHDMIN": {'units': '%'} ,  # min(100*GREENFRAC) 
        "LAI": {'units': 'm^2/m^2'} , # LAI12M after interpolated
        "XLAND": {'units': '', 'geoname': 'LU_INDEX'} ,  # if LU_INDEX==iswater or islake,2; else 1.  
        "ISLTYP": {'units': '', 'geoname': 'SOILCTOP'} ,
        "IVGTYP": {'units': '', 'geoname': 'LU_INDEX'} ,

        # from raw data file
        "TSK": {'units': 'K'} ,  # skin temperature from ERA5 
        "TSLB": {'units': 'K'} ,  # soil layer temp 
        "SMOIS": {'units': 'm^3/m^3'},  # layer volumetric total water content [m3/m3] !!!
        "DZS": {'units': 'm'} ,  # each soil layer depth
        "ZS": {'units': 'm'} ,   # soil layer 
        "SNOW": {'units': 'kg/m^2'} , #snow depth

        # add
        "SEAICE": {'units': ''} ,       # sea ice fraction (=0 for a land point)
        "CANWAT": {'units': 'kg/m^2'} , # set CANWAT = 0
    }

    geo_em = xr.open_dataset(geo_em_file)
    geo_lat, geo_lon = geo_em.XLAT_M.values[0], geo_em.XLONG_M.values[0]
    geo_lat_flat, geo_lon_flat = geo_lat.ravel(), geo_lon.ravel()

    iswater = int(geo_em.attrs['ISWATER'])
    islake = int(geo_em.attrs['ISLAKE'])
    issoilwater = int(geo_em.attrs['ISOILWATER'])
    water_mask = (geo_em.LU_INDEX.values[0] == iswater) | (geo_em.LU_INDEX.values[0] == islake)

    if pd.Timestamp(start_date).is_leap_year:
        LAI = xr.open_dataset(os.path.join(output_dir, 'LDASIN', 'LAI12M_leap.nc'))
    else:
        LAI = xr.open_dataset(os.path.join(output_dir, 'LDASIN', 'LAI12M.nc'))

    # N: 
    fname_setup = os.path.join(raw_data_dir, pd.to_datetime(start_date).strftime('%Y%m%d')+'00_setup.nc')
    data_vars_raw, raw_lat, raw_lon, time = nco_dataset(fname_setup,geo_lat,geo_lon)

    soil_data = []
    for var in ['skt', 'swvl1', 'swvl2', 'swvl3', 'swvl4', 'stl1', 'stl2', 'stl3', 'stl4', 'sd']:
        # data_var = data_vars_raw[var].rio.write_crs("epsg:4326",inplace=True).rio.interpolate_na()[0][::-1].values
        data_var = data_vars_raw[var].rio.write_crs("epsg:4326").rio.interpolate_na().squeeze().values[::-1, :]
        interp_spline = RectBivariateSpline(raw_lat, raw_lon, data_var, kx=1, ky=1)
        data_var_interpolated = interp_spline.ev(geo_lat_flat, geo_lon_flat).reshape(geo_lat.shape)
        # data_var_interpolated = xr.where((geo_em.LU_INDEX==iswater)|(geo_em.LU_INDEX==islake), np.nan, data_var_interpolated)
        data_var_interpolated_masked = np.where(water_mask, np.nan, data_var_interpolated)
        soil_data.append(data_var_interpolated_masked) 

    setup_file = xr.Dataset()

    for var in variables:

        dims = ('Time', 'south_north', 'west_east')
        dim4 = ('Time', 'soil_layers_stag', 'south_north', 'west_east')
        dim2 = ('Time', 'soil_layers_stag')

        if var == 'Times':
            data_var, dims =  [pd.to_datetime(start_date).strftime('%Y-%m-%d_%H:%M:%S')], ( 'Time' )

        ##################
        # from geo_em file
        ##################
        elif var in ['XLAT', 'XLONG', 'HGT', "MAPFAC_MX", "MAPFAC_MY", "IVGTYP"]:
            data_var = geo_em[variables[var]['geoname']].values
        
        elif var == 'URBLANDUSEF':
            data_var = geo_em[variables[var]['geoname']].sel(land_cat=12).values       # urban
            # if unresonable values appear, set them to 0
            data_var[data_var < 0] = 0              
            data_var[data_var > 1] = 0

        #######################
        # edit from geo_em file
        #######################
        # adjust to elevation
        elif var == 'TMN':
            data_var = geo_em[variables[var]['geoname']].values[0] - 0.0065 * geo_em['HGT_M'].values[0]
            # data_var = xr.where((geo_em.LU_INDEX==iswater)|(geo_em.LU_INDEX==islake), -1.e36, data_var).values
            data_var = [np.where(water_mask, -1.e36, data_var)]

        # gvfmax%field(:,:) = maxval(geo_em%veg,3)
        elif var == 'SHDMAX':
            data_var = geo_em.GREENFRAC.max(axis=1).values*100

        # gvfmin%field(:,:) = minval(geo_em%veg,3)
        elif var == 'SHDMIN':
            data_var = geo_em.GREENFRAC.min(axis=1).values*100

        elif var == 'LAI':
            if pd.Timestamp(start_date).is_leap_year:
                data_var = [LAI.sel(date='2020'+start_date[-6:])['LAI12M'].values]
            else:
                data_var = [LAI.sel(date='2021'+start_date[-6:])['LAI12M'].values]

        # if LU_INDEX==iswater or islake,2; else 1.
        elif var == 'XLAND':
            LU_data = geo_em[variables[var]['geoname']]
            # data_var = xr.where((LU_data==iswater)|(LU_data==islake), 2, 1).values
            data_var = [np.where(water_mask, 2, 1)]

        elif var == 'ISLTYP':
            dominant_index = geo_em['SOILCTOP'].argmax(dim='soil_cat') + 1
            dominant_value = geo_em['SOILCTOP'].max(dim='soil_cat')
            dominant_index_corrected = xr.where(
                (dominant_value < 0.01) | (dominant_value > 1.0), 8, dominant_index)
            data_var = xr.where(setup_file['XLAND']==2, issoilwater, dominant_index_corrected)
            data_var = xr.where((setup_file['XLAND']!=2)&(data_var==14), 8, data_var).values
        # if urbanlandusef exist and resonable, set IVGTYP=13(urban)
        elif var == 'IVGTYP':
            LU = xr.where(((geo_em.LANDUSEF).sel(land_cat=12)<=1)&((geo_em.LANDUSEF).sel(land_cat=12)>0), 13, geo_em.LU_INDEX)  
            data_var = LU.values

        ####################
        # from raw data file
        ####################
        elif var == 'TSK':
            data_var = [soil_data[0]]

        elif var == 'TSLB':
            data_var, dims = [ np.array(soil_data[5:9]) ], dim4

        elif var == 'SMOIS':
            data_var, dims = [ np.array(soil_data[1:5]) ], dim4

        elif var == 'ZS':
            data_var, dims = [ [0.035, 0.175, 0.64, 1.945] ], dim2

        elif var == 'DZS':
            data_var, dims = [ [0.07 , 0.21 , 0.72, 1.89 ] ], dim2

        elif var == 'SNOW':
            data_var =  [ soil_data[-1] * 1000]

        ########################
        # add SEAICE and CANWAT
        ########################
        elif var == 'SEAICE':
            data_var = np.zeros(geo_em.LU_INDEX.shape)

        elif var == 'CANWAT':
            data_var = np.zeros(geo_em.LU_INDEX.shape)

        print(var)
        print(dims, np.array(data_var).shape)

        setup_file[var] = ( dims, data_var )

        setup_file[var].attrs['units'] = variables[var]['units']

    setup_file = setup_file.fillna({'SNOW': -999})

    setup_file.attrs = geo_em.attrs

    output_filename = f"HRLDAS_setup_{pd.to_datetime(start_date).strftime('%Y%m%d')}00_d{geo_em_file[-4]}"

    setup_file.to_netcdf(os.path.join(output_dir, 'LDASIN', output_filename))

def create_lai_vegfra(geo_em_file, output_dir):

    if not os.path.exists(output_dir+"/LDASIN"):
        os.makedirs(output_dir+"/LDASIN")

    for var in ('LAI12M', 'GREENFRAC'):

        geo = xr.open_dataset(geo_em_file)
        LAI_geo = geo[var].sel(Time=0)

        LAI = xr.concat([LAI_geo, LAI_geo, LAI_geo, LAI_geo], dim="month")
        month = pd.date_range('2019-01-01', periods=48, freq='MS') + pd.DateOffset(days=14)
        LAI["month"] = ("month", month)

        date = pd.date_range('2019-01-15', '2022-12-15')
        LAI=LAI.rename({'month': 'date'})
        LAI=LAI.interp(date=date).to_dataset()

        # vegfra calibration
        if var=='GREENFRAC':
            LAI[var] = xr.where(LAI[var]<=0 ,0.01, LAI[var])
            LAI = LAI * 100

        # (iswater || islake ) == 0
        iswater = int(geo.attrs['ISWATER'])
        islake = int(geo.attrs['ISLAKE'])
        LU_geo = geo['LU_INDEX'].sel(Time=0)
        mask = ((LU_geo==iswater)|(LU_geo==islake)).expand_dims(dim={"date": date}, axis=0)
        LAI[var] = xr.where(mask, 0, LAI[var])

        if var=='LAI12M':
            LAI.sel(date=slice('2020-01-01','2020-12-31')).to_netcdf(os.path.join(output_dir, 'LDASIN', 'LAI12M_leap.nc'))
            LAI.sel(date=slice('2021-01-01','2021-12-31')).to_netcdf(os.path.join(output_dir, 'LDASIN', 'LAI12M.nc'))
        else:
            LAI.sel(date=slice('2020-01-01','2020-12-31')).to_netcdf(os.path.join(output_dir, 'LDASIN', 'GREENFRAC_leap.nc'))
            LAI.sel(date=slice('2021-01-01','2021-12-31')).to_netcdf(os.path.join(output_dir, 'LDASIN', 'GREENFRAC.nc'))


if __name__ == '__main__':
    gcm_name = "CMCC-ESM2"
    scenario = "ssp245"
    start_year = 2000
    end_year = 2022
    months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    loop_start_date = '01-01'
    loop_end_date = '12-31'
    area = [14.25, 99.9, 13.45, 101]
    levelists = '136'
    raw_data_dir = './ERA5_data/'
    output_dir = './hrldas_simulation/'
    geo_em_file = './geo_em.d01.nc'
    calendar = "gregorian"    


    #os.system(f"python3 compute_geopotential_on_ml.py {os.path.join(raw_data_dir, 'tq_ml.grib')} {os.path.join(raw_data_dir, 'zlnsp_ml.grib')} -o {os.path.join(raw_data_dir, 'z_out.grib')}")



    ZLVL = 30

    create_lai_vegfra(geo_em_file, output_dir)
    print("LAI : OK")
    create_setup_file(f'{str(start_year)}-{loop_start_date}',raw_data_dir, output_dir, geo_em_file)
    print("create_setup_file : OK" )
    num_cores = 30

    for year in range(start_year, end_year+1):
        create_LDASIN_files_parallel(f'{year}-{loop_start_date}', f'{year}-{loop_end_date}',
                                     raw_data_dir, output_dir, geo_em_file, gcm_name, scenario, levelists, ZLVL, calendar ,num_cores)


