#!/usr/bin/env python3
"""
RegCM → HRLDAS Soil Initialization Converter
--------------------------------------------
Developed by RU-CORE | Author: Ratchanan Srisawadwong (Nick)

Converts RegCM land surface output into an HRLDAS-compatible soil
initialization (setup) file by remapping RegCM soil layers to
ERA5/HRLDAS standard soil depths using thickness-weighted averaging.

Main processes:
  • Convert RegCM 10 soil layers → ERA5 4 soil layers
  • Thickness-weighted averaging for soil temperature
  • Water-mass-conserving conversion for soil moisture
  • Unit conversion:
      - Snow water equivalent → snow depth (m)
  • Output formatted for HRLDAS initialization

Input variables:
  ts     : Surface skin temperature → skt
  snw    : Snow water equivalent → sd
  tsoil  : Soil temperature profile → stl1–stl4
  mrsos  : Soil moisture → swvl1–swvl4

Output:
  HRLDAS soil setup NetCDF file

Usage:
  python regcm2setupfile.py <input.nc> <output.nc>

Requirements:
  xarray, numpy, netCDF4
"""

def run(infile, outfile):
    import xarray as xr
    import numpy as np
    # =====================================================
    #  ERA5 soil layer bounds (m)
    # =====================================================
    era5_bounds = np.array([
        [0.00, 0.07],
        [0.07, 0.28],
        [0.28, 1.00],
        [1.00, 2.89]
    ])
    era5_thickness = era5_bounds[:, 1] - era5_bounds[:, 0]
    
    # =====================================================
    #  RegCM soil layer bounds from midpoint
    # =====================================================
    regcm_mid_cm = np.array([0.7, 2.8, 6.2, 11.9, 21.2, 36.6, 62.0, 103.8, 172.8, 286.5])
    regcm_mid = regcm_mid_cm / 100.0
    
    interfaces = np.zeros(len(regcm_mid) + 1)
    interfaces[1:-1] = (regcm_mid[:-1] + regcm_mid[1:]) / 2
    interfaces[0] = 0.0
    interfaces[-1] = regcm_mid[-1] + (regcm_mid[-1] - regcm_mid[-2]) / 2
    
    regcm_bounds = np.column_stack([interfaces[:-1], interfaces[1:]])
    regcm_thickness = regcm_bounds[:, 1] - regcm_bounds[:, 0]
    
    # =====================================================
    #  Overlap matrix
    # =====================================================
    def calc_overlap(a_top, a_bot, b_top, b_bot):
        return max(0.0, min(a_bot, b_bot) - max(a_top, b_top))
    
    overlap = np.zeros((4, 10))
    for k in range(4):
        for i in range(10):
            overlap[k, i] = calc_overlap(
                era5_bounds[k, 0], era5_bounds[k, 1],
                regcm_bounds[i, 0], regcm_bounds[i, 1]
            )
    
    # =====================================================
    #  Load data (support all CF calendars)
    # =====================================================
    import cftime
    ds = xr.open_dataset(infile, decode_times=True, use_cftime=True)

    ds0 = ds.isel(time=0)
    ds0 = ds0.rename({
        "time": "valid_time",
        "lat": "latitude",
        "lon": "longitude"
    })

    time_val = ds0["valid_time"].item()

    
    # =====================================================
    #  Surface vars
    # =====================================================
    skt = ds0["ts"].astype("float32").reset_coords(drop=True)
    skt.name = "skt"
    skt.attrs = ds["ts"].attrs.copy()
    
    sd = (ds0["snw"] / 1000.0).astype("float32").reset_coords(drop=True)
    sd.name = "sd"
    sd.attrs = ds["snw"].attrs.copy()
    sd.attrs["units"] = "m"
    
    # =====================================================
    #  Soil layers
    # =====================================================
    rho_water = 1000.0
    stl_vars = {}
    swvl_vars = {}
    
    for k in range(4):
        weights = overlap[k, :]
        den = np.sum(weights)
    
        # ---- Soil temperature
        num = 0
        for i in range(10):
            if weights[i] > 0:
                num += ds0["tsoil"].isel(soil_layer=i) * weights[i]
    
        stl = (num / den).astype("float32").reset_coords(drop=True)
        stl.name = f"stl{k+1}"
        stl.attrs = ds["tsoil"].attrs.copy()
        stl_vars[stl.name] = stl
    
        # ---- Soil moisture
        sm_sum = 0
        for i in range(10):
            if weights[i] > 0:
                frac = weights[i] / regcm_thickness[i]
                sm_sum += ds0["mrsos"].isel(soil_layer=i) * frac
    
        swvl = (sm_sum / (rho_water * era5_thickness[k])).astype("float32").reset_coords(drop=True)
        swvl.name = f"swvl{k+1}"
        swvl.attrs = ds["mrsos"].attrs.copy()
        swvl.attrs["units"] = "m3 m-3"
        swvl_vars[swvl.name] = swvl
    
    # =====================================================
    #  Combine dataset
    # =====================================================
    ds_out = xr.Dataset({
        "skt": skt,
        "sd": sd,
        **stl_vars,
        **swvl_vars
    })
    
    ds_out = ds_out.expand_dims("valid_time")
    ds_out = ds_out.transpose("valid_time", "latitude", "longitude")
    ds_out = ds_out.assign_coords(valid_time=("valid_time", [time_val]))

    ds_out.attrs = ds.attrs.copy()
    
    # =====================================================
    #  Save
    # =====================================================
    ds_out.to_netcdf(outfile)
    #print(ds_out)
    print("Done:", outfile)
 


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Convert RegCM soil data to setup format")
    parser.add_argument("infile")
    parser.add_argument("outfile")
    args = parser.parse_args()
    run(args.infile, args.outfile)

if __name__ == "__main__":
    main()

