#!/bin/bash -l



#==============================================================================
# main script which postprocesses CESM present output to generate ensemble means
#==============================================================================


# time information
START=$(date +%s.%N)



#==============================================================================
# initialisation
#==============================================================================


# control or irrigation experiment
# run=control
# run=irrigation
for run in control       \
	       irrigation; do


# set case name
CASE=f.e122.F1850PDC5.f09_g16.${run}-io192


# run settings
flag_var=3;    # 1: TREFHTMX (daily maximum air temperature)
               # 2: AT       (apparent temperature,        i.e. heat stress)
               # 3: WBGT     (wet bulb global temperature, i.e. heat stress)
flag_ref=2;    # 1: 20cirr - including irrigation (default)
               # 2: 20cc   - excluding irrigation (alternative)
flag_block=2   # 2: atm data
flag_stream=1  # 1: h1 output block


# define all possible variables
vars=('empty' 'TREFHTMX' 'AT' 'WBGT')
refs=('empty' '20cirr' '20cc')


# set run settings based on flags
var=${vars[$flag_var]}


# define whether 20th Century reference run ('pic' hereafter) includes or excludes irrigation
ref=${refs[$flag_ref]}


# define percentile values file from 20cirr
ref_pctl=/net/cfc/landclim1/wthiery/cesm_output/ensmean/atm/postprocessed/f.e122.F1850PDC5.f09_g16.${ref}-io192.merged_atm_h1_${var}_pctl.nc


# set output directory name
outdir=/net/cfc/landclim1/wthiery/cesm_output


# define percentages to compute percentiles from
percentages=('50.000' '90.000' '95.000' '97.500' '99.000' '99.500' '99.750' '99.900' '99.950' '99.975' '99.990')



#==============================================================================
# manipulations: checks
#==============================================================================


# define all possible blocks
blocks=('empty' 'lnd' 'atm')


# define all possible output streams
streams=('empty' 'h1')


# set block value based on flag_block
block=${blocks[$flag_block]}
stream=${streams[$flag_stream]}
echo; echo "This script will process $stream output from the $block model for case $CASE"


# create directory where postprocessed data will be stored
postdir=$outdir/ensmean/$block/postprocessed
if [ ! -d "$postdir" ]; then  # check if directory exists
  mkdir -p $postdir
fi


# change to postprocessing directory
cd ${postdir}
pwd


# count number of ensemble members
nensmembers=($outdir/${CASE}.???/)
nensmembers=${#nensmembers[@]}
echo; echo "counting $nensmembers ensemble members:"


# print them to screen
ensmembers=$outdir/${CASE}.???/
for ensmember in $ensmembers; do
  echo "$ensmember"
done


# get all the files $block/postprocessed directory of the first member
echo
for ensmember in $ensmembers; do
  nfiles=($ensmember/$block/postprocessed/*)
  nfiles=${#nfiles[@]}
  echo "counting $nfiles files for ensemble member $ensmember"
done



#==============================================================================
# manipulations: concatenate ensemble members
#==============================================================================


# concatenate ensemble members
echo; echo; echo "concatenate ensemble members: ${var}" 
cdo -O cat $outdir/${CASE}.???/$block/postprocessed/${CASE}.???_${block}_${stream}_${var}.nc  ${postdir}/${CASE}.merged_${block}_${stream}_${var}.nc


# change to output directory
cd ${postdir}
echo; echo; echo "changing to postprocessing directory:" 
pwd
echo; echo



#==============================================================================
# manipulations: percentiles wrt predefined probabilities
# PDcontrol only
#==============================================================================


if    [[ $run = "control" ]];  then  # control run, so compute percentiles needed for ctl2irr

  # loop over percentages and compute the percentile maps
  for percentage in "${percentages[@]}"; do
    echo 'computing' ${percentage}'th percentile in ctl - '${var}
    cdo -s setname,${var}_P${percentage} -timpctl,${percentage} ${CASE}.merged_${block}_${stream}_${var}.nc -timmin ${CASE}.merged_${block}_${stream}_${var}.nc -timmax ${CASE}.merged_${block}_${stream}_${var}.nc pctl_ctl_${var}_P${percentage}.nc
  done
  cdo -O merge pctl_ctl_${var}_P*.nc ${CASE}.merged_${block}_${stream}_${var}_pctl.nc
  rm -f pctl_ctl_${var}_P*.nc


fi



#==============================================================================
# manipulations: probabilities wrt percentiles from PIcontrol and PDcontrol
#==============================================================================


if    [[ $run = "control" ]];  then  # control run, so compute pic2ctl



  #================================
  # pic2ctl: GHG effect
  #================================

  
  # print message
  echo; echo; echo 'pic2ctl: GHG effect'
  
  
  # loop over percentages
  for percentage in "${percentages[@]}"; do
    echo 'computing percentage in ctl lying below' ${percentage}'th percentile in pic'
    cdo -s setname,prob_pic2ctl_P${percentage} -mulc,100 -timmean -lt ${CASE}.merged_${block}_${stream}_${var}.nc -selname,${var}_P${percentage} ${ref_pctl} prob_pic2ctl_${var}_P${percentage}.nc
  done
# # #   cdo -O merge prob_pic2ctl_P*.nc ${CASE}.merged_${block}_${stream}_${var}_prob_pic2ctl.nc
  cdo -O merge prob_pic2ctl_${var}_P*.nc ${CASE}.merged_${block}_${stream}_${var}_prob_${ref}2ctl.nc
  rm -f prob_pic2ctl_${var}_P*.nc


elif    [[ $run = "irrigation" ]];  then  # irrigation run, so compute ctl2irr and pic2irr



  #================================
  # ctl2irr: irrigation effect
  #================================

  
  # print message
  echo; echo; echo 'ctl2irr: irrigation effect'
  
  
  # loop over percentages
  for percentage in "${percentages[@]}"; do
    echo 'computing percentage in irr lying below' ${percentage}'th percentile in ctl'
    cdo -s setname,prob_ctl2irr_P${percentage} -mulc,100 -timmean -lt ${CASE}.merged_${block}_${stream}_${var}.nc -selname,${var}_P${percentage} f.e122.F1850PDC5.f09_g16.control-io192.merged_${block}_${stream}_${var}_pctl.nc prob_ctl2irr_${var}_P${percentage}.nc
  done
  cdo -O merge prob_ctl2irr_${var}_P*.nc ${CASE}.merged_${block}_${stream}_${var}_prob_ctl2irr.nc
  rm -f prob_ctl2irr_${var}_P*.nc
  
  
  
  #================================
  # pic2irr: GHG + irrigation effect
  #================================

  
  # print message
  echo; echo; echo 'pic2irr: GHG + irrigation effect'
  
  
  # loop over percentages
  for percentage in "${percentages[@]}"; do
    echo 'computing percentage in irr lying below' ${percentage}'th percentile in pic'
    cdo -s setname,prob_pic2irr_P${percentage} -mulc,100 -timmean -lt ${CASE}.merged_${block}_${stream}_${var}.nc -selname,${var}_P${percentage} ${ref_pctl} prob_pic2irr_${var}_P${percentage}.nc
  done
# # #   cdo -O merge prob_pic2irr_P*.nc ${CASE}.merged_${block}_${stream}_${var}_prob_pic2irr.nc
  cdo -O merge prob_pic2irr_${var}_P*.nc ${CASE}.merged_${block}_${stream}_${var}_prob_${ref}2irr.nc
  rm -f prob_pic2irr_${var}_P*.nc

  
fi



done # loop over runs


# time information
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo; echo "Elapsed time is" $DIFF "seconds."; echo
