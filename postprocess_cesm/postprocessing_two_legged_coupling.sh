#!/bin/bash -l



#==============================================================================
# main script which postprocesses CESM present output
#==============================================================================



#==============================================================================
# initialisation
#==============================================================================


# set case name ...  (comment out 'done' at the end of the file)
# CASE=f.e122.F1850PDC5.f09_g16.control-io192.001
# CASE=f.e122.F1850PDC5.f09_g16.irrigation-io192.001
# CASE=f.e122.F1850PDC5.f09_g16.20cirr-io192.005

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
#for CASE in f.e122.F1850PDC5.f09_g16.pic-io192.001    \
#            f.e122.F1850PDC5.f09_g16.pic-io192.002    \
#            f.e122.F1850PDC5.f09_g16.pic-io192.003    \
#            f.e122.F1850PDC5.f09_g16.pic-io192.004    \
#            f.e122.F1850PDC5.f09_g16.pic-io192.005; do
# # # for CASE in f.e122.F1850PDC5.f09_g16.20cc-io192.001    \
# # #             f.e122.F1850PDC5.f09_g16.20cc-io192.002    \
# # #             f.e122.F1850PDC5.f09_g16.20cc-io192.003    \
# # #             f.e122.F1850PDC5.f09_g16.20cc-io192.004    \
# # #             f.e122.F1850PDC5.f09_g16.20cc-io192.005    \
# # #             f.e122.F1850PDC5.f09_g16.20cirr-io192.001  \
# # #             f.e122.F1850PDC5.f09_g16.20cirr-io192.002  \
# # #             f.e122.F1850PDC5.f09_g16.20cirr-io192.003  \
# # #             f.e122.F1850PDC5.f09_g16.20cirr-io192.004  \
# # #             f.e122.F1850PDC5.f09_g16.20cirr-io192.005; do



# run settings
flag_block=1  # 1: lnd data
              # 2: atm data - does not work
flag_stream=1 # 1: h1 output block
              # 2: h2 output block - does not work
flag_srex=0   # 0: do not compute averages for srex region
              # 1: compute averages for srex region (works only on h1 data) - does not work


# set output directory name
outdir=/net/cfc/landclim1/wthiery/cesm_output


# define start and end year
startyear=1981 # exclude spinup years
endyear=2010   # last year of the simulation
#startyear=1861 # exclude spinup years
#endyear=1890   # last year of the simulation
#startyear=1901 # exclude spinup years
#endyear=1930   # last year of the simulation

              
# define srex mask netcdf file              
srex_masks=/home/wthiery/documents/research/cesm_present/observations/srex/srex_masks.nc              
srex_masks_land=/home/wthiery/documents/research/cesm_present/observations/srex/srex_masks_land.nc              


# srex regions for which we want spatial averages
srex_names=('WNA' 'CNA' 'MED' 'WAS' 'SAS' 'SEA' 'EAS')
# srex_names=('ENA' 'SSA' 'CEU' 'SAF' 'SAF' 'SAU')  # For Annette Hirsch



#==============================================================================
# manipulations
#==============================================================================


# time information
START=$(date +%s.%N)


# define all possible blocks
blocks=('empty' 'lnd' 'atm')


# define all possible output streams
streams=('empty' 'h1' 'h2')


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
  echo; echo; echo "processing ${block} data for year $YEAR"

  
  
  #================================
  # lnd data
  #================================

    
  # merge over all time steps
  if   [ $flag_block -eq 1 ];  then  # lnd data


  
    #================================
    # grid cell data (h1)
    #================================

    
    # merge over all time steps
    if   [ $flag_stream -eq 1 ];  then  # h1 data


      # select MAMJJA, warmest months in most of the irrigation hotspots
      cdo selseason,MAMJJA -selname,FSH                                       ${CASE}.clm2.${stream}.${YEAR}-01-01-00000.nc FSH_${stream}_${YEAR}_MAMJJA.nc
      cdo selseason,MAMJJA -selname,SOILWATER_10CM                            ${CASE}.clm2.${stream}.${YEAR}-01-01-00000.nc SM10_${stream}_${YEAR}_MAMJJA.nc
      cdo selseason,MAMJJA -selname,TREFHT $outdir/${CASE}/atm/hist/year$YEAR/${CASE}.cam.${stream}.${YEAR}-01-01-00000.nc  TREFHT_${stream}_${YEAR}_MAMJJA.nc
      
      
      # Calculate the standard deviation over all timesteps for QH and T2
      cdo -s timstd FSH_${stream}_${YEAR}_MAMJJA.nc    FSH_${stream}_${YEAR}_MAMJJA_timstd.nc
      cdo -s timstd TREFHT_${stream}_${YEAR}_MAMJJA.nc TREFHT_${stream}_${YEAR}_MAMJJA_timstd.nc
   
      
      # Calculate the correlation between variables
      # Soil moisture and QH
      cdo timcor SM10_${stream}_${YEAR}_MAMJJA.nc FSH_${stream}_${YEAR}_MAMJJA.nc cor_SM10_FSH_${YEAR}.nc
      # QH and T2
      cdo timcor FSH_${stream}_${YEAR}_MAMJJA.nc TREFHT_${stream}_${YEAR}_MAMJJA.nc cor_FSH_TREFHT_${YEAR}.nc
      # Soil moisture and T2
      cdo timcor SM10_${stream}_${YEAR}_MAMJJA.nc TREFHT_${stream}_${YEAR}_MAMJJA.nc cor_SM10_TREFHT_${YEAR}.nc
 
      # Calculate the terrestrial-leg coupling
      cdo mul FSH_${stream}_${YEAR}_MAMJJA_timstd.nc cor_SM10_FSH_${YEAR}.nc IL_${YEAR}_MAMJJA.nc    
      ncrename -v FSH,IL IL_${YEAR}_MAMJJA.nc
      ncatted -a units,IL,o,c,"" IL_${YEAR}_MAMJJA.nc
      ncatted -a long_name,IL,o,c,"Thermal Terrestrial Coupling" IL_${YEAR}_MAMJJA.nc
      ncatted -a standard_name,IL,o,c,"Thermal Terrestrial Coupling" IL_${YEAR}_MAMJJA.nc

      # Calculate the atmospheric-leg coupling
      cdo mul TREFHT_${stream}_${YEAR}_MAMJJA_timstd.nc cor_FSH_TREFHT_${YEAR}.nc IA_${YEAR}_MAMJJA.nc
      ncrename -v TREFHT,IA IA_${YEAR}_MAMJJA.nc
      ncatted -a units,IA,o,c,"" IA_${YEAR}_MAMJJA.nc
      ncatted -a long_name,IA,o,c,"Thermal Atmospheric Coupling" IA_${YEAR}_MAMJJA.nc
      ncatted -a standard_name,IA,o,c,"Thermal Atmospheric Coupling" IA_${YEAR}_MAMJJA.nc

      # Calculate the full land-atmosphere coupling
      cdo mul TREFHT_${stream}_${YEAR}_MAMJJA_timstd.nc cor_SM10_TREFHT_${YEAR}.nc ILA_${YEAR}_MAMJJA.nc
      ncrename -v TREFHT,ILA ILA_${YEAR}_MAMJJA.nc
      ncatted -a units,ILA,o,c,"" ILA_${YEAR}_MAMJJA.nc
      ncatted -a long_name,ILA,o,c,"Thermal Land-Atmospheric Coupling" ILA_${YEAR}_MAMJJA.nc
      ncatted -a standard_name,ILA,o,c,"Thermal Land-Atmospheric Coupling" ILA_${YEAR}_MAMJJA.nc

      # Merge into one netcdf file
      cdo merge IL_${YEAR}_MAMJJA.nc IA_${YEAR}_MAMJJA.nc ILA_${YEAR}_MAMJJA.nc ${postdir}/coupling_${YEAR}_MAMJJA.nc

      # Clean-up
      rm FSH_${stream}_${YEAR}_MAMJJA.nc
      rm SM10_${stream}_${YEAR}_MAMJJA.nc
      rm TREFHT_${stream}_${YEAR}_MAMJJA.nc
      rm FSH_${stream}_${YEAR}_MAMJJA_timstd.nc
      rm TREFHT_${stream}_${YEAR}_MAMJJA_timstd.nc
      rm cor_SM10_FSH_${YEAR}.nc
      rm cor_FSH_TREFHT_${YEAR}.nc
      rm cor_SM10_TREFHT_${YEAR}.nc
      rm IL_${YEAR}_MAMJJA.nc
      rm IA_${YEAR}_MAMJJA.nc
      rm ILA_${YEAR}_MAMJJA.nc

     
    
    #================================
    # PFT data (h2)
    #================================

  
    elif [ $flag_stream -eq 2 ];  then  # h2 data
    
    
      # do nothing
      echo 'do nothing'
      
    
    fi

    
  
  #================================
  # atm data
  #================================

  
  elif [ $flag_block -eq 2 ];  then  # atm data


      # do nothing
      echo 'do nothing'
 
    
  fi

  
done


# change to output directory
cd ${postdir}
echo; echo; echo "changing to postprocessing directory:" 
pwd


# merge variables in time - yearly coupling data
cdo -O mergetime coupling_????_MAMJJA.nc ${CASE}_${block}_${stream}_coupling_MAMJJA_year.nc
cdo timmean ${CASE}_${block}_${stream}_coupling_MAMJJA_year.nc ${CASE}_${block}_${stream}_coupling_MAMJJA_timmean.nc


# clean up
rm -f coupling_????_MAMJJA.nc


# time information
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo; echo "Elapsed time is" $DIFF "seconds."; echo


done
