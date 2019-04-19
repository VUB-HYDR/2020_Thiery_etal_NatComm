#!/bin/bash


# bash script to set up a CESM1.2.2 transient historical run with year-2000 land-use
# information taken from:
# https://wiki.iac.ethz.ch/Climphys/ProjectCESM122SetupControlRun



#=====================================================
# initialisation
#=====================================================


# define path where cases are stored
CASEDIR=~/cesm/cases


# define path where CESM scripts are stored
SCRIPTSDIR=~/cesm/cesm122-trunk/scripts


# define machine
MACH=euler           # new ETH Cluster


# define compiler
COMPILER=intel;  MPILIB=openmpi      # io


# define compset - more information below
COMP=F1850PDC5       # PIPD_CAM5_CLM40%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV !!! TRANSIENT LAND COVER DATASET MISSES PFT16 !!! - stick to CLM40 for now because CLM4.5 requires you to run the crop model and we don't want that


# define resolution
RES=f09_g16          # 0.9x1.25_gx1v6 (1°) resolution (default NTASKS=128)


# define number of cores
NTASKS=192           # number of tasks (= cores used) - on euler  192 = 8 nodes each has 24 cores


# define other information
ver=e122                              # version, used in CASE name, e122 = CESM 1.2.2
flavor="${COMPILER:0:1}${MPILIB:0:1}" # either im, io, pm, po
type=20cc                             # type of run
desc=$type-$flavor$NTASKS             # description, used in CASE name
nnn=001                               # unique number, used in CASE name
comp=$( echo ${COMP:0:1} | tr '[:upper:]' '[:lower:]' )


# define run settings
start_year=1896  # start year of the simulation (1856)
end_year=2012    # end year of the SST dataset
simul_length=12  # months


# define namelist script
nl_file=nl_pic_${nnn}.sh


# set whether you are in production mode or test mode
production=true # true=final production runs; false=testing
nresubmit=34    # if simul_length=12, this is number of years -1; use for final production runs including spinup



#=====================================================
# 1. create new case
#=====================================================


# get path where setup scripts are
SETUPDIR=$(pwd)


# Change into the new case directory 
cd $SCRIPTSDIR


# generate casename
CASE=$comp.$ver.$COMP.$RES.$desc.$nnn


# create a new case in the directory 'cases'
./create_newcase -case $CASEDIR/$CASE  -res $RES -compset $COMP -mach $MACH -compiler $COMPILER -mpilib $MPILIB



#=====================================================
# 2. invoke cesm_setup
#=====================================================


# Change into the new case directory 
cd $CASEDIR/$CASE


# copy this script to the case directory for reference
cp $SETUPDIR/`basename "$0"` .


# Modify env_mach_pes.xml
# Note: If you change env_mach_pes.xml later, you have to do ./cesm_setup -clean; and again ./cesm_setup; ./$CASE.$MACH.build ! 
./xmlchange -file env_mach_pes.xml -id NTASKS_ATM -val $NTASKS
./xmlchange -file env_mach_pes.xml -id NTASKS_LND -val $NTASKS
./xmlchange -file env_mach_pes.xml -id NTASKS_ICE -val $NTASKS
./xmlchange -file env_mach_pes.xml -id NTASKS_OCN -val $NTASKS
./xmlchange -file env_mach_pes.xml -id NTASKS_CPL -val $NTASKS
./xmlchange -file env_mach_pes.xml -id NTASKS_GLC -val $NTASKS
./xmlchange -file env_mach_pes.xml -id NTASKS_ROF -val $NTASKS
./xmlchange -file env_mach_pes.xml -id NTASKS_WAV -val $NTASKS


# modify env_run.xml
./xmlchange -file env_run.xml -id RUN_STARTDATE     -val ${start_year}"-01-01"  # year in which CESM starts running (1st January)
./xmlchange -file env_run.xml -id SSTICE_YEAR_ALIGN -val ${start_year}          # model year when you want the data to start
./xmlchange -file env_run.xml -id SSTICE_YEAR_START -val ${start_year}          # date on the file when the data starts
./xmlchange -file env_run.xml -id SSTICE_YEAR_END   -val ${end_year}            # date on the file when the data ends
./xmlchange -file env_run.xml -id RUN_TYPE          -val startup                # Set to run type to startup (is the default)


# modify env_run.xml - MAY BE CHANGED ANYTIME during a run
./xmlchange -file env_run.xml -id STOP_OPTION       -val nmonths
### ./xmlchange -file env_run.xml -id STOP_OPTION       -val ndays
./xmlchange -file env_run.xml -id STOP_N            -val ${simul_length}


# force 2000_control configuration in a transient simulation (e.g. static 2000 land use)
# see $SCRIPTSDIR/models/lnd/clm/bld/build-namelist for the information
# and $SCRIPTSDIR/models/lnd/clm/bld/namelist_files/use_cases/ for the cases
./xmlchange CLM_BLDNML_OPTS="-use_case 2000_control" -append


# introduce source code modifications
cp $SETUPDIR/SourceMods/src.clm/* ./SourceMods/src.clm


# introduce changes to CLM namelist
cp $SETUPDIR/$nl_file .
. ./$nl_file


# Configure case
./cesm_setup



#=====================================================
# 3. Build/Compile the model
#    note: restart from here when making source mods
#=====================================================


# Build/Compile the model
./$CASE.build



#=====================================================
# 4. Submit run to the batch queue
#=====================================================


# change walltime from 23:59h to 03:59h - for 2degree run
### sed -i 's/#BSUB -W 23:59/#BSUB -W 03:59/g' $CASE.run


# run the model
if [ "$production" = true ]; then  # production runs

  ./xmlchange -file env_run.xml -id CONTINUE_RUN -val FALSE
  ./xmlchange -file env_run.xml -id RESUBMIT -val ${nresubmit}
  ./$CASE.submit

  # Urs his way, never got short term archiver to work properly
  #./xmlchange -file env_run.xml -id CONTINUE_RUN -val FALSE
  #cp $SCRIPTSDIR/resubmit_cesm .
  #./resubmit_cesm -case $CASE -start 1 -end ${nresubmit}

else                               # single year run

  ./$CASE.submit

fi


# Check the submitted job with the command bjobs -w
bjobs -w



#=====================================================
# Notes
#=====================================================


# for creating a user_compset, do:
#  note: USER_COMP has to be long name, COMP has to be short name if you define a case yourself. else just use short name and adapt create_case
#  note: 20TR_CAM5_CLM40%CN_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV is what Urs uses as 'hist' - I want to change this to CLM45 and hopefully CN to SP
#USER_COMP=some_long_name   # for defining a user_compset
#COMP=F20C5TR      # 20TR_CAM5_CLM40%SP_CICE%PRES_DOCN%DOM_RTM_SGLC_SWAV - stick to CLM40 for now because CLM4.5 requires you to run the crop model and we don't want that


#  then create the case by invoking
#./create_newcase -case $CASEDIR/$CASE  -res $RES -user_compset $USER_COMP -mach $MACH -compiler $COMPILER -mpilib $MPILIB


















