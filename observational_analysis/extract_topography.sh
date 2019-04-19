#!/bin/bash


#============================================================
# extract_topography.sh
# Script which extracts topography and land fraction from 
# an input file downloaded from the CLM community website.
#
# Author: Auke Visser
# Date: 25.10.2016
#============================================================

#set input file name
#File was obtained from http://www.clm-community.eu/index.php?menuid=221&reporeid=260
infile=domain2016102710109.nc

#set output file name
outfile=CLMdata_topography.0.5deg.nc


#set input directory name
dir=/net/exo/landclim/wthiery/observational_analysis/Data

#set output directory name
dir=/net/exo/landclim/wthiery/observational_analysis/Data

#============================================================
# manipulations
#============================================================

#copy infile to outfile
cdo select,name=FR_LAND,HSURF $dir/$infile $dir/$outfile
rm $dir/$infile
