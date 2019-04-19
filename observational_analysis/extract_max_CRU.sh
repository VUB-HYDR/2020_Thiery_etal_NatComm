#!/bin/bash


#============================================================
# extract_max_CRU.sh
# Script which calculates annual minimum/maximum of an 
# input data set
#
# Author: Auke Visser
# Date: 8.11.2016
#============================================================

#set input file name
infile=cru_ts3.22.1901.2013.tmn.dat.nc

#set output file name
outfile=cru_ts3.22.1901.2013.tmn_max.dat.nc

#set input directory name
indir=/net/exo/landclim/data/dataset/CRUTS/v3.22/0.5deg_lat-lon_1m/original

#set output directory name
outdir=/net/exo/landclim/wthiery/observational_analysis/Data

#============================================================
# manipulations
#============================================================

# cdo yearmin $indir/$infile $outdir/$outfile
cdo yearmax $indir/$infile $outdir/$outfile
