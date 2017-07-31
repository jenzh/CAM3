#! /bin/tcsh -f
#
#=======================================================================
#
# build_isomodel.csh
#
# build isotopic cam3 on olympus
# David Noone <dcn@colorado.edu> - Thu Jul  1 10:06:38 MDT 2004
#
module load impi/4.0.3.008

set isovers = 'cam3_0_iso0_20'
set MYPWD=${PWD}  # current working directory

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: initial, restart, or branch.
## $nelapse is the number of timesteps to integrate, or number of days if negative.

set case         = tst5		# < ------ *** CHANGE THIS ***
set runtype      = initial
#set nelapse      = -123			# -123 to stop end year1
#set nelapse      = -488
#set nelapse      = -1218
set nelapse      = -365


## Additional namelist settings
set user_namelist = "trace_water=.true., wisotope=.true., readtrace=.false., reset_csim_iceprops=.true., trace_gas=.false."
set user_namelist = "${user_namelist}, mfilt=1"
set user_namelist = "${user_namelist}, bndtag='${MYPWD}/input_data/tag.20161010.nc'"
#set user_namelist = "${user_namelist}, nhtfrq=-24"
# Single column continuous output - Boulder, Darwin, Casey, Ascension Island
#set user_namelist = "${user_namelist}, fincl1lonlat='255e_40n','130e_0n','111e_66s','346e_8s'"

## Set resolution of model
set dycore       = eul
set resolution   = 64x128
set nlev         = 26

## Set number of tracers
set nadv = 102   # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + 10 tagged regions
#set nadv = 48   # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + 4 tagged regions
#set nadv = 21   # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + 1 taged region
#set nadv = 12   # Qx3 + H2Ox3 + HDOx3 + H218Ox3
#set nadv = 18   # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + H217Ox3 + HTOx3
#set nadv = 16    # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + (CH4, N2O, CFC11, CFC12)
#set nadv = 22   # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + H217Ox3 + HTOx3 + (CH4, N2O, CFC11, CFC12)

#######################################################################
## ENVIRONMENTAL VARIABLES NEEDED TO CONFIGURE CAM
##
## This part of the script that is set up for OLYMPUS
## It very likely needs to be changed for another 
## system. Indeed it does for RANGER with regard to 
## the file locations 
set OS = `uname -s`;
#setenv PGI_FC TRUE
#setenv MPIF90_FC TRUE
setenv USER_CC mpiicc
setenv USER_FC mpiifort
#setenv USER_FC /usr/local/mvapich/bin/mpif90
setenv INC_MPI /software/intel/impi/4.0.3.008/include64
echo "${0}:  Set INC_MPI to $INC_MPI";
setenv LIB_MPI /software/intel/impi/4.0.3.008/lib64
echo "${0}:  Set LIB_MPI to $LIB_MPI";
#echo "${0}:  USE MPI PF90 WRAPPER";
#echo "${0}:  USING pgf90 COMPILER WITH pgcc"
setenv INC_NETCDF /software/apps/netcdf/4.2/i1214-hdf5-1.8.9/include
echo "${0}:  Set INC_NETCDF to $INC_NETCDF";
setenv MOD_NETCDF $INC_NETCDF
echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF";
setenv LIB_NETCDF /software/apps/netcdf/4.2/i1214-hdf5-1.8.9/lib
echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF";
#######################################################################

## Do our best to get sufficient stack memory
limit stacksize unlimited

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
#set camroot      = /data2/dcn/models/CAM3/isotope/${isovers}/cam1
set camroot      = ${MYPWD}/${isovers}/cam1

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
#setenv CSMDATA     /fs/cgd/csm/inputdata
setenv CSMDATA    ${MYPWD}/input_data

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
#set wrkdir       = /ptmp/$LOGNAME
set wrkdir       = ${MYPWD}/camruns
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case/rundir
set cfgdir       = $camroot/models/atm/cam/bld
# User modificatoins
set mymods       = ${MYPWD}/SourceMods

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

#=======================================================================
## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir || echo "cd $blddir failed" && exit 1

    $cfgdir/configure -spmd \
                      -nadv $nadv \
                      -dyn $dycore \
                      -nlev $nlev \
                      -usr_src $mymods \
                      -cppdefs '-DNO_R16' -cppdefs '-DTAGGING' -fopt '-O2' \
                      -res $resolution     || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
#    rm -f Depends
    gmake -j 128 >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

#=======================================================================
## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -case $case -runtype $runtype -o $rundir/namelist \
 -namelist "&camexp nelapse=$nelapse mss_irt=0 $user_namelist /"  || echo "build-namelist failed" && exit 1
#
#=======================================================================
## create run script (in rundir)

set jobName	= $rundir/run_isomodel.csh
echo $jobName
rm $jobName
cat > $jobName << EOF
#!/bin/csh -f
#SBATCH -t 06:35:00
#SBATCH -N 1
#SBATCH --job-name=CAMiso
#SBATCH --output=slurm.out

echo "mhmmm"
cd $rundir
echo "mi sto incazzando"
rm -rf cam.out
touch cam.out

setenv OMP_NUM_THREADS 1
echo "running CAM in $rundir"
mpprun ../bld/cam < namelist >> cam.out || echo "CAM run failed" && exit 1
echo "come cazzo si fa"
EOF
chmod 755 $jobName

exit 0

