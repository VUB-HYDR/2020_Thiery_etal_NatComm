#! /bin/bash
#Author: Wim Thiery (Original: Mathias Hauser and Micah Wilhelm)


#use this file to change the settings for the run
#i.e. source (. ./) it in the simulation-creation file


# note:
# - by setting hist_dov2xy to .false., you generate output in vector format, which allows you to get variables per pft through postprocessing.


#CAM
#add history fields
cat > user_nl_cam << EOF
fincl1 = 'Z500', 'T', 'FSNTOA', 'FSNRTOAS'
fincl2 = 'FLNS', 'FSNS', 'LHFLX', 'PRECT', 'PRECTMX', 'PRECC', 'PRECL', 'PS', 'PSL', 'RHREFHT', 'SHFLX', 'TREFHT', 'TREFHTMN', 'TREFHTMX', 'U10', 'U850', 'V850', 'OMEGA850', 'Q850', 'T850', 'TMQ', 'TS', 'TSMN', 'TSMX'
nhtfrq = 0, -24
mfilt  = 1, 365
pertlim=1.e-14
EOF


#CLM
#add history fields
cat > user_nl_clm << EOF
hist_fincl2 = 'TREFMNAV', 'TREFMXAV', 'TSA', 'TV', 'TG', 'TSOI_10CM', 'Rnet', 'FSA', 'FIRA', 'SWdown', 'SWup', 'LWdown', 'LWup', 'Qle', 'FCEV', 'FCTR', 'FGEV', 'FSH', 'FSH_G', 'FSH_V', 'FGR', 'FGR12', 'WASTEHEAT', 'HEAT_FROM_AC', 'FSM', 'HCSOI', 'SOILWATER_10CM', 'QDRAI', 'QINFL', 'QOVER', 'QRUNOFF', 'QSOIL', 'QVEGE', 'QVEGT', 'QIRRIG', 'QCHARGE', 'QRGWL', 'QSNWCPICE', 'QFLOOD', 'RAIN', 'SNOW', 'SOILLIQ:I', 'SOILICE:I', 'TLAKE', 'TBOT', 'THBOT', 'WT', 'WA:I', 'H2OCAN:I', 'H2OSNO:I', 'ELAI', 'TLAI', 'TSAI', 'RHO', 'RAH', 'TAF', 'EMV', 'EMG'
hist_fincl3 = 'TREFMNAV', 'TREFMXAV', 'TSA', 'TV',                    'Rnet', 'FSA', 'FIRA',           'SWup',           'LWup', 'Qle', 'FCEV', 'FCTR', 'FGEV',        'FSH_G', 'FSH_V', 'FGR',          'WASTEHEAT', 'HEAT_FROM_AC',                 'SOILWATER_10CM',                                       'QSOIL', 'QVEGE', 'QVEGT',                                                                      'SOILLIQ:I',                                                                                                                    'TAF'                          
hist_nhtfrq = 0, -24, -24
hist_mfilt = 1, 365, 365
hist_dov2xy = .true., .true., .false.
EOF


#CICE
#store nothing
cat > user_nl_cice << EOF
histfreq = 'x','x','x','x','x'
EOF


