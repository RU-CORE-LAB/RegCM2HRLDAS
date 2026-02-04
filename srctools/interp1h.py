#!/usr/bin/env python3
"""
RegCM → HRLDAS Temporal Interpolation Tool
------------------------------------------
Developed by RU-CORE | Author: Ratchanan Srisawadwong (Nick)

Interpolates RegCM atmospheric forcing data from 6-hourly to hourly resolution
while preserving CF-compliant time handling for all calendar types
(standard, noleap, 360_day, etc.).

Interpolation methods:
  • Cubic  : u, v, t, hgt, topo
  • Linear : q (specific humidity)

Key features:
  • Supports CFTime calendars (noleap, 360_day)
  • Preserves metadata and variable attributes
  • Applies NetCDF4 compression and chunking
  • Keeps HRLDAS-compatible dimension structure

Input:
  Pressure-level atmospheric forcing file (6-hourly)

Output:
  Hourly interpolated forcing file (HRLDAS-ready)

Usage:
  python interp1h.py <input.nc> <output.nc>

Requirements:
  xarray, numpy, pandas, cftime, netCDF4
"""

def run(infile, outfile):
    import xarray as xr
    import numpy as np
    import pandas as pd
    import cftime

    time_dim = "valid_time"

    vars_cubic  = ["u", "v", "t", "hgt", "topo"]
    vars_linear = ["q"]

    ds = xr.open_dataset(
        infile,
        engine="netcdf4",
        decode_times=True,
        use_cftime=True,   # force cftime for noleap / 360_day
        chunks={
            time_dim: 10,
            "latitude": 45,
            "longitude": 65
        }
    )

    # =====================
    # FIX model_level
    # =====================
    if "model_level" in ds.variables and ds["model_level"].ndim == 0:
        level_val = int(ds["model_level"].values)
        ds = ds.drop_vars("model_level")
        ds = ds.expand_dims({"model_level": [level_val]})

    # =====================
    # CREATE NEW TIME AXIS (ALL CALENDARS)
    # =====================
    time_var = ds[time_dim]
    t0 = time_var.values[0]
    t1 = time_var.values[-1]

    calendar = time_var.attrs.get("calendar", "standard")

    # ---- Standard / Gregorian (numpy datetime64) ----
    if np.issubdtype(time_var.dtype, np.datetime64):
        t_new = pd.date_range(
            start=pd.to_datetime(t0),
            end=pd.to_datetime(t1),
            freq="1h"
        )
    # ---- CFTime calendars (noleap / 360_day / etc.) ----
    elif isinstance(t0, cftime.datetime):
        t_new = xr.cftime_range(
            start=t0,
            end=t1,
            freq="1h",
            calendar=calendar
        )

    # ---- Numeric time (rare but kept safe) ----
    elif np.issubdtype(time_var.dtype, np.number):
        t_new = np.arange(t0, t1 + 1, 1)

    else:
        raise RuntimeError(f"Unsupported time coordinate type: {type(t0)}")

    # =====================
    # INTERPOLATION
    # =====================
    ds_cubic_i = ds[vars_cubic].interp(
        **{time_dim: t_new},
        method="cubic",
        kwargs={"fill_value": "extrapolate"}
    )

    ds_linear_i = ds[vars_linear].interp(
        **{time_dim: t_new},
        method="linear",
        kwargs={"fill_value": "extrapolate"}
    )

    ds_out = xr.merge([ds_cubic_i, ds_linear_i])

    # =====================
    # COPY ATTRIBUTES
    # =====================
    ds_out.attrs = ds.attrs
    ds_out[time_dim].attrs = ds[time_dim].attrs
    for v in ds_out.data_vars:
        ds_out[v].attrs = ds[v].attrs

    # =====================
    # ENCODING
    # =====================
    encoding = {
        v: {
            "zlib": True,
            "complevel": 4,
            "shuffle": True,
            "chunksizes": (
                1,
                ds_out.sizes["model_level"],
                ds_out.sizes["latitude"],
                ds_out.sizes["longitude"]
            )
        }
        for v in ds_out.data_vars
    }

    ds_out.to_netcdf(
        outfile,
        format="NETCDF4",
        engine="netcdf4",
        encoding=encoding
    )

    print(f"\tDone: {outfile} [OK]")

def main(): 
    import argparse 
    parser = argparse.ArgumentParser(description="Interpolate ATM forcing to 1-hour interval") 
    parser.add_argument("infile") 
    parser.add_argument("outfile") 
    args = parser.parse_args() 
    run(args.infile, args.outfile) 

if __name__ == "__main__": 
    main()
