RUN: RegCM2HRLDAS Tools
```bash
#!/bin/bash

inpath="/home/regcm/_STORAE_/_SAS_/CMIP6_25km_RAW/RCPs/CMCC_ssp245/output/"

model="CMCC-ESM2"
scenario="ssp245"
sYR=2020
eYR=2021
domain=98,102,12.5,15

outpath="./output"

regcm2hrldas-atm -m ${model} -p ${scenario} -s ${sYR} -e ${eYR} -d ${domain} -i ${inpath} -o ${outpath}
regcm2hrldas-srf -m ${model} -p ${scenario} -s ${sYR} -e ${eYR} -d ${domain} -i ${inpath} -o ${outpath}

month=01
regcm2hrldas-setup -m ${model} -p ${scenario} -s ${sYR}${month} -d ${domain} -i ${inpath} -o ${outpath}
```
