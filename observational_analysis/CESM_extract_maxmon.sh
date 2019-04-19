#!/bin/bash

#============================================================
# CESM_extract_maxmon.sh
# Script which calculates the annual maximum monthly mean
# daily maximum temperature from CESM model output.
#
# Author: Auke Visser
# Date: 5.1.2017
#============================================================

#set input file name
# infile1=f.e122.F1850PDC5.f09_g16.irrigation-io192.001_atm_h1_TREFHTMX.nc
# infile2=f.e122.F1850PDC5.f09_g16.irrigation-io192.002_atm_h1_TREFHTMX.nc
# infile3=f.e122.F1850PDC5.f09_g16.irrigation-io192.003_atm_h1_TREFHTMX.nc
# infile4=f.e122.F1850PDC5.f09_g16.irrigation-io192.004_atm_h1_TREFHTMX.nc
# infile5=f.e122.F1850PDC5.f09_g16.irrigation-io192.005_atm_h1_TREFHTMX.nc

infile_irr=f.e122.F1850PDC5.f09_g16.irrigation-io192.ensmean_atm_h1_monmean.nc
infile_ctl=f.e122.F1850PDC5.f09_g16.control-io192.ensmean_atm_h1_monmean.nc
infile_pic=f.e122.F1850PDC5.f09_g16.pic-io192.ensmean_atm_h1_monmean.nc
infile_20cc=f.e122.F1850PDC5.f09_g16.20cc-io192.ensmean_atm_h1_monmean.nc
#set output file name
outfile_irr=CESM.WT.1981-2010.TREFHT_irr_ensmonmean_yearmax.nc
outfile_ctl=CESM.WT.1981-2010.TREFHT_ctl_ensmonmean_yearmax.nc
outfile_pic=CESM.WT.1981-2010.TREFHT_pic_ensmonmean_yearmax.nc
outfile_20cc=CESM.WT.1981-2010.TREFHT_20cc_ensmonmean_yearmax.nc

#set input directory name
# indir1=/net/cfc/landclim1/wthiery/cesm_output/f.e122.F1850PDC5.f09_g16.irrigation-io192.001/atm/postprocessed
# indir2=/net/cfc/landclim1/wthiery/cesm_output/f.e122.F1850PDC5.f09_g16.irrigation-io192.002/atm/postprocessed
# indir3=/net/cfc/landclim1/wthiery/cesm_output/f.e122.F1850PDC5.f09_g16.irrigation-io192.003/atm/postprocessed
# indir4=/net/cfc/landclim1/wthiery/cesm_output/f.e122.F1850PDC5.f09_g16.irrigation-io192.004/atm/postprocessed
# indir5=/net/cfc/landclim1/wthiery/cesm_output/f.e122.F1850PDC5.f09_g16.irrigation-io192.005/atm/postprocessed

# indir=/net/cfc/landclim1/wthiery/cesm_output/ensmean/atm/postprocessed
indir=/net/cfc/landclim1/wthiery/cesm_output/ensmean/atm/postprocessed

#set output directory name
outdir=/net/exo/landclim/wthiery/observational_analysis/Data/CESM_output

#============================================================
# manipulations
#============================================================

# cdo ensmean $indir1/$infile1 $indir2/$infile2 $indir3/$infile3 $indir4/$infile4 $indir5/$infile5 $outdir/ensmean.nc
# cdo yearmax -monmean $outdir/ensmean.nc $outdir/$outfile
# rm $outdir/ensmean.nc

cdo yearmax -selname,TREFHT $indir/$infile_irr $outdir/$outfile_irr
# cdo yearmax -selname,TREFHTMX $indir/$infile_ctl $outdir/$outfile_ctl
# cdo yearmax -selname,TREFHTMX $indir/$infile_pic $outdir/$outfile_pic
# cdo yearmax -selname,TREFHT $indir/$infile_20cc $outdir/$outfile_20cc
