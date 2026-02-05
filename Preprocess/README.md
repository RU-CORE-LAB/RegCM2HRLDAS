# RegCM2HRLDAS forcing Toolchain Installer

Developed by: RU-CORE  
Author: Ratchanan Srisawadwong (Nick)

This toolchain is a preprocessing system designed to convert output from the RegCM climate model into forcing and initialization files that can be used by HRLDAS.

The system automates key steps including:

- Program compilation  
- Environment configuration  
- Data conversion and formatting

## Toolchain Capabilities

The installer builds and configures all required components for converting RegCM outputs into HRLDAS input data as follows:

### 1. RegCM Sigma-Level Conversion Tool

The installer compiles:

sigma2p_at_hPa.F90

This code has been adapted from the original RegCM post-processing utilities to support vertical level conversion consistent with ERA5 pressure levels:  
https://confluence.ecmwf.int/plugins/servlet/mobile?contentId=108117123#content/view/108117123

This tool is used to:

- Convert RegCM sigma levels → pressure levels  
- Align vertical atmospheric structure with ERA5 pressure levels  
- Serve as a critical step in preparing atmospheric forcing data for HRLDAS

### 2. Python Environment Setup

The system creates a dedicated Conda environment for processing and installs all required Python libraries for the data transformation workflow

### 3. Command-Line Utilities for RegCM Output Processing

After installation, the following command-line tools will be available:

- regcm2hrldas-atm  
- regcm2hrldas-srf  
- regcm2hrldas-setup  
- regcm2setupfile  
- interp1h  

These tools cover atmospheric forcing, temporal processing, surface fluxes, and soil initialization data preparation for HRLDAS.

## System Requirements

(Required for both Toolchain installation and data processing)

### For Installation

| Software | Purpose |
|----------|---------|
| RegCM (compiled) | Provides required libraries and files for compiling the atmospheric vertical conversion tool |
| Conda (Anaconda or Miniconda) | Manages the Python environment used by the processing tools |
| gfortran or ifort | Compiles Fortran programs |
| NetCDF libraries | Required to compile Fortran programs that read/write NetCDF files |

### For Data Processing

| Tool | Purpose |
|------|---------|
| CDO (Climate Data Operators) | Used for metadata handling and climate data processing such as subsetting, calculations, and format adjustments |
| NCO (NetCDF Operators) | Used for editing NetCDF files and managing metadata |

## Setup

### Download the Toolchain

Clone the repository from GitHub:
```bash
git clone git@github.com:RU-CORE-LAB/RegCM2HRLDAS.git
cd RegCM2HRLDAS/Preprocess
```


### Configure Before Installation

Edit the user configuration file:
```bash
vim ./config/regcm2hrldas.cfg
```

Example configuration:
```bash
REGCM_DIR=/home/regcm/RegCM_2025/RegCM-5.0.0
FC=gfortran
CONDA_ENV="regcm2hrldas"
```
### Installation

Run the following command to start the installation process:

```bash
bash setup_regcm2hrldas.sh
```

The installer will automatically perform the following steps:

- Compile `sigma2p_at_hPa.F90`
- Check whether Conda is installed
- Create a Conda environment
- Install required Python libraries
- Install the processing toolchain package

### After Installation

Activate the environment with:

```bash
conda activate regcm2hrldas
```

The system will then be ready to use all tools included in this toolchain.

## Command-Line Tools

### 1. regcm2hrldas-atm

Converts 6-hourly RegCM atmospheric output into hourly HRLDAS forcing data on pressure levels.

**Process**

- Read ATM (6-hourly) files from RegCM output
- Spatial subsetting to target domain
- Convert sigma → pressure levels
- Extract variables (ta, ua, va, hus, hgt, topo)
- Temporal interpolation (6-hourly → hourly) using `interp1h`
- Compute geopotential height
- Format output to HRLDAS-compatible structure

**Output file**

```
<YYYYMM>_t_u_v_q_z_<GCM>_<SCENARIO>_pressure_level.nc
```

### 2. regcm2hrldas-srf

Converts RegCM surface flux data into HRLDAS surface forcing.

**Process**

- Read SRF (hourly) files from RegCM output
- Spatial subsetting to target domain
- Extract variables (ps, rlds, rsds, pr)
- Convert variable units as required
- Format output to HRLDAS-compatible structure

**Output file**

```
<YYYYMM>_ssrd_strd_sp_tp_<GCM>_<SCENARIO>_single_layer.nc
```

### 3. regcm2hrldas-setup

Generates the HRLDAS soil initialization file from RegCM SRF output.

**Process**

- Read SRF (hourly) files from RegCM output
- Spatial subsetting to target domain
- Extract soil and surface state variables
- Convert RegCM soil layers → ERA5/HRLDAS soil layers  
  (thickness-weighted averaging) using `regcm2setupfile`

**Output file**

```
<YYYYMM>0100_setup.nc
```
### 4. interp1h

Performs temporal interpolation of RegCM atmospheric forcing data  
from 6-hourly to hourly resolution while supporting all CF-compliant calendars  
(e.g., standard, noleap, 360_day).

**Features**

- Temporal interpolation from 6-hourly → hourly
- Supports CFTime calendars (Gregorian, noleap, 360_day)
- Applies NetCDF4 compression and chunking
- Preserves metadata and variable attributes
- Maintains HRLDAS-compatible dimension structure

### 5. regcm2setupfile

Converts RegCM land surface data into HRLDAS soil layer structure.

**Main Process**

- Spatial subsetting to match the target domain
- Extract soil and surface state variables from RegCM data
- Convert RegCM soil layer depths into standard ERA5/HRLDAS soil layer depths  
  using thickness-weighted averaging

## Example

An example workflow is provided in:
```bash
vim ./example_convert/run2hrldas.sh
```

### Example Script

```bash
inpath="/home/regcm/_STORAE_/_SAS_/CMIP6_25km_RAW/RCPs/EC-Earth-Veg/output/"

model="EC-Earth3-Veg"
scenario="ssp245"
sYR=2020
eYR=2021
domain=98,102,12.5,15

outpath="./output"


# model layer
regcm2hrldas-atm -m ${model} -p ${scenario} -s ${sYR} -e ${eYR} -d ${domain} -i ${inpath} -o ${outpath}

# single layer
regcm2hrldas-atm -m ${model} -p ${scenario} -s ${sYR} -e ${eYR} -d ${domain} -i ${inpath} -o ${outpath}

# setupfile
regcm2hrldas-srf -m ${model} -p ${scenario} -s ${sYR} -e ${eYR} -d ${domain} -i ${inpath} -o ${outpath}

month=01
regcm2hrldas-setup -m ${model} -p ${scenario} -s ${sYR}${month} -d ${domain} -i ${inpath} -o ${outpath}




