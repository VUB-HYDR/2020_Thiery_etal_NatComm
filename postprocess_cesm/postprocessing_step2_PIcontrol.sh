#!/bin/bash -l



#==============================================================================
# main script which computes percentile values in the 20C reference simulations
#==============================================================================


# time information
START=$(date +%s.%N)


# Notes:
# - in contrast to earlier version of this script, it is now operating on the pixel level (instead of computing the srex-average time series)



#==============================================================================
# initialisation
#==============================================================================


# set case name ...  (comment out 'done' at the end of the file)
CASE=f.e122.F1850PDC5.f09_g16.20cirr-io192
# CASE=f.e122.F1850PDC5.f09_g16.20cc-io192


# run settings
flag_var=3;    # 1: TREFHTMX (daily maximum air temperature)
               # 2: AT       (apparent temperature,        i.e. heat stress)
               # 3: WBGT     (wet bulb global temperature, i.e. heat stress)
flag_block=2   # 2: atm data
flag_stream=1  # 1: h1 output block


# define all possible variables
vars=('empty' 'TREFHTMX' 'AT' 'WBGT')


# set run settings based on flags
var=${vars[$flag_var]}


# set output directory name
outdir=/net/cfc/landclim1/wthiery/cesm_output


# define start and end year
# startyear=1861 # exclude spinup years
# endyear=1890   # last year of the simulation
startyear=1901 # exclude spinup years
endyear=1930   # last year of the simulation


# define percentages to compute percentiles from
percentages=('50.000' '90.000' '95.000' '97.500' '99.000' '99.500' '99.750' '99.900' '99.950' '99.975' '99.990')



#==============================================================================
# manipulations
#==============================================================================


# define all possible blocks
blocks=('empty' 'lnd' 'atm')


# define all possible output streams
streams=('empty' 'h1')


# security
if   [ $flag_block -eq 1 ];  then  # lnd data
    echo "$stream data from the $block model is not processed"; exit 1
fi


# set block value based on flag_block
block=${blocks[$flag_block]}
stream=${streams[$flag_stream]}
echo; echo "This script will process $stream output from the $block model for case $CASE"



#==============================================================================
# manipulations: concatenate ensemble members
#==============================================================================


# create directory where postprocessed data will be stored
postdir=$outdir/ensmean/$block/postprocessed
if [ ! -d "$postdir" ]; then  # check if directory exists
  mkdir -p $postdir
fi


# concatenate ensemble members
echo; echo; echo "concatenate ensemble members: ${var}" 
cdo -O cat $outdir/${CASE}.???/$block/postprocessed/${CASE}.???_${block}_${stream}_${var}.nc  ${postdir}/${CASE}.merged_${block}_${stream}_${var}.nc



#==============================================================================
# manipulations: compute percentiles
#==============================================================================


# change to output directory
cd ${postdir}
echo; echo; echo "changing to postprocessing directory:" 
pwd
echo; echo;


# loop over percetiles and compute the maps
for percentage in "${percentages[@]}"; do

  echo 'computing' ${percentage}'th percentile in pic - '${var}
  cdo -s setname,${var}_P${percentage} -timpctl,${percentage} ${CASE}.merged_${block}_${stream}_${var}.nc -timmin ${CASE}.merged_${block}_${stream}_${var}.nc -timmax ${CASE}.merged_${block}_${stream}_${var}.nc ${CASE}.merged_${block}_${stream}_${var}_P${percentage}.nc
  
done
cdo -O merge ${CASE}.merged_${block}_${stream}_${var}_P*.nc ${CASE}.merged_${block}_${stream}_${var}_pctl.nc


# clean up
rm -f ${CASE}.merged_${block}_${stream}_${var}_P*.nc


# time information
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo; echo "Elapsed time is" $DIFF "seconds."; echo


