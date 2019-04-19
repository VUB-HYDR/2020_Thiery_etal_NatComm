#!/bin/bash -l



#==============================================================================
# main script which postprocesses Fluxnet-MTE reanalysis
#==============================================================================


# time information
START=$(date +%s.%N)



#==============================================================================
# initialisation
#==============================================================================


# set file name
file_LHF=gleam_v3a.e.025deg.1d.????.nc


# set output directory name
outdir=/net/exo/landclim/data/dataset/GLEAM/v3a/0.25deg_lat-lon_1d/processed/netcdf


# set refgrid file name
refgrid=../f.e122.F1850PDC5.f09_g16.control-io192.001_atm_h1_timmean_refgrid.nc


# define start and end year - exclude spinup years
startyear=1981
endyear=2010



#==============================================================================
# manipulations
#==============================================================================


# merge files for different years
cdo mergetime $outdir/$file_LHF GLEAM_merged.nc


# get temporal means - all interpolated to CESM grid
cdo setname,ET_mean -setunit,mm/yr -timmean -remapcon2,$refgrid -mulc,365.25 -selyear,${startyear}/${endyear} -selname,e GLEAM_merged.nc GLEAM_${startyear}-${endyear}_timmean.nc


# get multiyear monthly - all interpolated to CESM grid
cdo setname,ET_mean -setunit,mm/d -ymonmean -remapcon2,$refgrid -selyear,${startyear}/${endyear} -selname,e GLEAM_merged.nc GLEAM_${startyear}-${endyear}_ymonmean.nc


# clean up
rm -f GLEAM_merged.nc


# time information
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo; echo "Elapsed time is" $DIFF "seconds."; echo
