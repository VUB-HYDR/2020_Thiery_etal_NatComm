#!/bin/bash


#============================================================
# regrid.sh
# Script which regrids CRU data to the CESM domain using
# CDO second-order conservative remapping after
# calculating the hottest month of the year (mean monthly 
# daily max from CRU tmx)
#
# Author: Auke Visser
# Last modified: 06.01.2017
#============================================================

#set input file name
infile=cru_ts3.22.1901.2013.tmp.dat.nc

#set output file name
# outfile_HID=hid_v1.0.f_irr_hyde_final_ir.0.5deg.10y.1900-2005.nc
# outfile_HID=hid_v1.0.f_irr_hyde_final_ir.2.5x3.75deg.10y.1900-2005.nc
outfile=cru_ts3.22.1901.2013.0.9375x1.25.tmp_max.dat.nc

#set input directory name
indir=/net/exo/landclim/data/dataset/CRUTS/v3.22/0.5deg_lat-lon_1m/original

#set output directory name
outdir=/net/exo/landclim/wthiery/observational_analysis/Data

#set regrid file
# refgrid=/net/exo/landclim/data/dataset/CRUTS/v3.22/0.5deg_lat-lon_1m/original/cru_ts3.22.1901.2013.tmn.dat.nc
# refgrid=/net/exo/landclim/data/dataset/E-OBS/v14/0.25deg_lat-lon_1d/original/tg_0.25deg_reg_v14.0.nc
# refgrid=/net/exo/landclim/data/dataset/HadEX2/20150106/2.5x3.75deg_lat-lon_1y/original/H2_TXx_1901-2010_RegularGrid_global_3.75x2.5deg_LSmask.nc
# refgrid=/net/cfc/landclim1/wthiery/cesm_output/ensmean/atm/postprocessed/f.e122.F1850PDC5.f09_g16.control-io192.ensmean_atm_h1_monmean.nc
# refgrid=/home/wthiery/documents/research/cesm_present/postprocessing/constants/f.e122.F1850PDC5.f09_g16.irrigation-io192.001_constants.nc
refgrid=/net/cfc/landclim1/wthiery/cesm_output/ensmean/atm/postprocessed/f.e122.F1850PDC5.f09_g16.control-io192.ensmean_atm_h1_TREFHTMX.nc


#============================================================
# manipulations
#============================================================

#regrid 
# if [ -f $outdir/$outfile_HID ]; then
#     echo "File $outfile_HID already exists in the directory $outdir"
# else
#     cdo gridarea $indir/$infile_HID $outdir/HID_area.nc
#     cdo mulc,10000 $indir/$infile_HID $outdir/HID_m2.nc
#     cdo div $outdir/HID_m2.nc $outdir/HID_area.nc $outdir/HID_irrfrac_nativeres.nc
#     cdo remapcon2,$refgrid $outdir/HID_irrfrac_nativeres.nc $outdir/$outfile_HID
#     #echo "File $outfile_HID does not exit in the directory $outdir"
# fi

# if [ -f $outdir/$outfile_HID ]; then
#     echo "File $outfile_HID already exists in the directory $outdir"
# else
#     cdo -f nc remapcon2,$refgrid -div -mulc,10000 $indir/$infile_HID -gridarea $indir/$infile_HID $outdir/$outfile_HID
# fi

cdo remapcon2,$refgrid -yearmax $indir/$infile $outdir/$outfile
