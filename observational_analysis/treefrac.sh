#!/bin/bash


#============================================================
# treefrac.sh
# Script which calculates the total tree fraction out of seven 
# categories.
#
# Author: Auke Visser
# Date: 24.11.2016
#============================================================

#set input file name
indir_tf=/net/exo/landclim/wthiery/observational_analysis/Data/LUC/HYDE_AREAVEG

#============================================================
# manipulations
#============================================================

# add 7 fields from every file, loop over 110 files

for i in {1900..1902};
do
#     cdo expr,'tf=TrpEBF+TrpDBF+TmpEBF+TmpENF+TmpDBF+BorENF+BorDNF;' $indir_tf/land-cover_hyde_landcover_yr$i.nc $indir_tf/land-cover_hyde_totaltreefrac_yr$i.nc
done
