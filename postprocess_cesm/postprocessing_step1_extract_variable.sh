#!/bin/bash -l



#==============================================================================
# main script which postprocesses raw CESM present output
# extract variable of interest and/or compute heat stress indices per ensmember
#==============================================================================


# time information
START=$(date +%s.%N)



#==============================================================================
# initialisation
#==============================================================================


# set case name ...  (comment out 'done' at the end of the file)
# CASE=f.e122.F1850PDC5.f09_g16.control-io192.001
# CASE=f.e122.F1850PDC5.f09_g16.irrigation-io192.001


# .. or loop over ensemble members
for CASE in f.e122.F1850PDC5.f09_g16.control-io192.001    \
            f.e122.F1850PDC5.f09_g16.control-io192.002    \
            f.e122.F1850PDC5.f09_g16.control-io192.003    \
            f.e122.F1850PDC5.f09_g16.control-io192.004    \
            f.e122.F1850PDC5.f09_g16.control-io192.005    \
            f.e122.F1850PDC5.f09_g16.irrigation-io192.001 \
            f.e122.F1850PDC5.f09_g16.irrigation-io192.002 \
            f.e122.F1850PDC5.f09_g16.irrigation-io192.003 \
            f.e122.F1850PDC5.f09_g16.irrigation-io192.004 \
            f.e122.F1850PDC5.f09_g16.irrigation-io192.005; do
# for CASE in f.e122.F1850PDC5.f09_g16.pic-io192.001    \
# 	        f.e122.F1850PDC5.f09_g16.pic-io192.002    \
# 	        f.e122.F1850PDC5.f09_g16.pic-io192.003    \
# 	        f.e122.F1850PDC5.f09_g16.pic-io192.004    \
# 	        f.e122.F1850PDC5.f09_g16.pic-io192.005; do
# # for CASE in f.e122.F1850PDC5.f09_g16.20cirr-io192.001    \
# # 	        f.e122.F1850PDC5.f09_g16.20cirr-io192.002    \
# # 	        f.e122.F1850PDC5.f09_g16.20cirr-io192.003    \
# # 	        f.e122.F1850PDC5.f09_g16.20cirr-io192.004    \
# # 	        f.e122.F1850PDC5.f09_g16.20cirr-io192.005    \
# #             f.e122.F1850PDC5.f09_g16.20cc-io192.001      \
# # 	        f.e122.F1850PDC5.f09_g16.20cc-io192.002      \
# # 	        f.e122.F1850PDC5.f09_g16.20cc-io192.003      \
# # 	        f.e122.F1850PDC5.f09_g16.20cc-io192.004      \
# # 	        f.e122.F1850PDC5.f09_g16.20cc-io192.005; do


# run settings
flag_var=1;    # 1: TREFHTMX (daily maximum air temperature)
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
startyear=1981 # exclude spinup years
endyear=2010   # last year of the simulation
# startyear=1861 # exclude spinup years
# endyear=1890   # last year of the simulation
# # startyear=1901 # exclude spinup years
# # endyear=1930   # last year of the simulation



#==============================================================================
# manipulations
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
postdir=$outdir/${CASE}/$block/postprocessed
if [ ! -d "$postdir" ]; then  # check if directory exists
  mkdir -p $postdir
fi


# loop over years
for YEAR in $(seq -w $startyear $endyear); do


  # change to output directory
  cd $outdir/${CASE}/$block/hist/year$YEAR
  pwd

  
  # send message to screen
  echo; echo; echo "processing $block data for year $YEAR"

  
  
  #================================
  # lnd data
  #================================

    
  # merge over all time steps
  if   [ $flag_block -eq 1 ];  then  # lnd data


 
    echo "$stream data from the $block model is not processed"; exit 1
    

    
  
  #================================
  # atm data
  #================================

  
  elif [ $flag_block -eq 2 ];  then  # atm data


    # merge over all time steps
    if   [ $flag_stream -gt 1 ];  then
      echo "$stream data from the $block model does not exist"; exit 1      
    fi
    
    
    if   [ $flag_var -eq 1 ];  then  # TREFHTMX

    
      # select variable TREFHTMX only
      cdo selname,TREFHTMX ${CASE}.cam.${stream}.${YEAR}-01-01-00000.nc $postdir/TREFHTMX_${stream}_${YEAR}.nc

      
    elif   [ $flag_var -eq 2 ];  then  # AT
    
    
      # compute Water vapour pressure
      # intermediate step, not used in final calculation
      #cdo -setunit,'hPa' -setname,'e' -expr,'e = RHREFHT / 100 * 6.105 * exp ( 17.27 * (TREFHTMX-273.15) / ( 237.7 + (TREFHTMX-273.15) ) )' ${CASE}.cam.${stream}.${YEAR}-01-01-00000.nc $postdir/e_${stream}_${YEAR}.nc

    
      # compute Apparent temperature - Version including the effects of temperature, humidity, and wind
      # http://www.bom.gov.au/info/thermal_stress/?cid=003bl08
      # !!! uses 24-h average T and RH !!!
      cdo -setunit,'degC' -setname,'AT' -expr,'AT = (TREFHT-273.15) + 0.33 * (RHREFHT / 100 * 6.105 * exp ( 17.27 * (TREFHT-273.15) / ( 237.7 + (TREFHT-273.15) ) ) ) - 0.70 * U10 - 4.00' ${CASE}.cam.${stream}.${YEAR}-01-01-00000.nc $postdir/${var}_${stream}_${YEAR}.nc

      
    elif   [ $flag_var -eq 3 ];  then  # WBGT
    
    
      # compute Wet Bulb Global Temperature - approximation to the WBGT used by the Bureau of Meteorology
      # http://www.bom.gov.au/info/thermal_stress/?cid=003bl08
      # !!! uses 24-h average T and RH !!!
      cdo -setunit,'degC' -setname,'WBGT' -expr,'WBGT = 0.567 * (TREFHT-273.15) + 0.393 * (RHREFHT / 100 * 6.105 * exp ( 17.27 * (TREFHTMX-273.15) / ( 237.7 + (TREFHTMX-273.15) ) )) + 3.94' ${CASE}.cam.${stream}.${YEAR}-01-01-00000.nc $postdir/${var}_${stream}_${YEAR}.nc

      
    fi
    
    
  fi

  
done


# change to output directory
cd ${postdir}
echo; echo; echo "changing to postprocessing directory:" 
pwd
echo; echo;


# merge data in time - daily data
cdo -O cat ${var}_${stream}_????.nc ${CASE}_${block}_${stream}_${var}.nc


# clean up
rm -f ${var}_${stream}_????.nc


# time information
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo; echo "Elapsed time is" $DIFF "seconds."; echo


done
