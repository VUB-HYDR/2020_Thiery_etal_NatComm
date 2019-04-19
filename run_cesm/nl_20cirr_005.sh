#! /bin/bash
#Author: Wim Thiery (Original: Mathias Hauser and Micah Wilhelm)


#use this file to change the settings for the run
#i.e. source (. ./) it in the simulation-creation file


# note:
# - we changed surfdata_0.9x1.25_simyr2000_irrcr_c100916.nc to surfdata_0.9x1.25_simyr2000_irrcr_c100916_1915.nc (the original file was modified to impose the irrigation extend map corresponding to the year 1915)
# - we changed aerocom_mam3_dms_surf_1849-2100_c111017.nc to aerocom_mam3_dms_surf_1849-2006_c090804.nc (needed to run prior to 1944)
# - by setting hist_dov2xy to .false., you generate output in vector format, which allows you to get variables per pft through postprocessing.


#CAM
#add history fields
cat > user_nl_cam << EOF
fincl1 = 'Z500', 'T', 'FSNTOA', 'FSNRTOAS'
fincl2 = 'FLNS', 'FSNS', 'LHFLX', 'PRECT', 'PRECTMX', 'PRECC', 'PRECL', 'PS', 'PSL', 'RHREFHT', 'SHFLX', 'TREFHT', 'TREFHTMN', 'TREFHTMX', 'U10', 'U850', 'V850', 'OMEGA850', 'Q850', 'T850', 'TMQ', 'TS', 'TSMN', 'TSMX'
nhtfrq = 0, -24
mfilt  = 1, 365
pertlim=5.e-14
srf_emis_specifier		= 'DMS       -> /cluster/work/climate/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/aerocom_mam3_dms_surf_1849-2006_c090804.nc',
         'SO2       -> /cluster/work/climate/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/RCP45_mam3_so2_surf_1850-2020_c20120308.nc',
         'SOAG      -> /cluster/work/climate/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/RCP45_mam3_soag_1.5_surf_1850-2020_c20120308.nc',
         'bc_a1     -> /cluster/work/climate/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/RCP45_mam3_bc_surf_1850-2020_c20120308.nc',
         'num_a1    -> /cluster/work/climate/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/RCP45_mam3_num_a1_surf_1850-2020_c20120308.nc',
         'num_a2    -> /cluster/work/climate/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/RCP45_mam3_num_a2_surf_1850-2020_c20120308.nc',
         'pom_a1    -> /cluster/work/climate/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/RCP45_mam3_oc_surf_1850-2020_c20120308.nc',
         'so4_a1    -> /cluster/work/climate/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/RCP45_mam3_so4_a1_surf_1850-2020_c20120308.nc',
         'so4_a2    -> /cluster/work/climate/cesm/inputdata/atm/cam/chem/trop_mozart_aero/emis/RCP45_mam3_so4_a2_surf_1850-2020_c20120308.nc'
EOF


#CLM
#add history fields
cat > user_nl_clm << EOF
hist_fincl2 = 'TREFMNAV', 'TREFMXAV', 'TSA', 'TV', 'TG', 'TSOI_10CM', 'Rnet', 'FSA', 'FIRA', 'SWdown', 'SWup', 'LWdown', 'LWup', 'Qle', 'FCEV', 'FCTR', 'FGEV', 'FSH', 'FSH_G', 'FSH_V', 'FGR', 'FGR12', 'WASTEHEAT', 'HEAT_FROM_AC', 'FSM', 'HCSOI', 'SOILWATER_10CM', 'QDRAI', 'QINFL', 'QOVER', 'QRUNOFF', 'QSOIL', 'QVEGE', 'QVEGT', 'QIRRIG', 'QCHARGE', 'QRGWL', 'QSNWCPICE', 'QFLOOD', 'RAIN', 'SNOW', 'SOILLIQ:I', 'SOILICE:I', 'TLAKE', 'TBOT', 'THBOT', 'WT', 'WA:I', 'H2OCAN:I', 'H2OSNO:I', 'ELAI', 'TLAI', 'TSAI', 'RHO', 'RAH', 'TAF', 'EMV', 'EMG'
hist_fincl3 = 'TREFMNAV', 'TREFMXAV', 'TSA', 'TV',                    'Rnet', 'FSA', 'FIRA',           'SWup',           'LWup', 'Qle', 'FCEV', 'FCTR', 'FGEV',        'FSH_G', 'FSH_V', 'FGR',          'WASTEHEAT', 'HEAT_FROM_AC',                 'SOILWATER_10CM',                                       'QSOIL', 'QVEGE', 'QVEGT',                                                                      'SOILLIQ:I',                                                                                                                    'TAF'                          
hist_nhtfrq = 0, -24, -24
hist_mfilt = 1, 365, 365
hist_dov2xy = .true., .true., .false.
fsurdat = '/cluster/home/wthiery/cesm/setup_newrun/SurfdataMods/surfdata_0.9x1.25_simyr2000_irrcr_c100916_irr1915.nc'
EOF


#CICE
#store nothing
cat > user_nl_cice << EOF
histfreq = 'x','x','x','x','x'
EOF


