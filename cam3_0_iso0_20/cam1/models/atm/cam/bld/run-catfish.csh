#! /bin/tcsh -f
#
#=======================================================================
#
# run-catfish.csh
#
# Run cam3 on catfish, which does not use PBS.
# David Noone <dcn@colorado.edu> - Thu Jul  1 10:06:38 MDT 2004
#
#
#=======================================================================
nice +10
#=======================================================================

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: initial, restart, or branch.
## $nelapse is the number of timesteps to integrate, or number of days if negative.

set case         = camrun		# < ------ *** CHANGE THIS ***
set runtype      = initial
set nelapse      = -1			# -123 to stop end year1
                                        # -1218 from sept year 0 to Dec31 year 3

## Additional namelist settings
set user_namelist = ' '

## Set resolution of model
#set res = "256x512"            # T159 (ERA40)
#set res = "128x256"            # T85 (IPCC)
set res = "64x128"              # T42 (DEFAULT)
#set res = '48x96'              # T31 (Paleoclimate)
#set res = "32x64"              # T21
#set res = "8x16"               # T5 (REALLY fast for testing)


## Set parallel stratergy
set nthreads = 2
set smp  = FALSE		# DOES NOT WORK!
set spmd = TRUE

## Local paths and compiler settings
##
## Home of netcdf, fortran compiler (pgi by default, but let's be explict)
##
setenv USER_FC pgf90
setenv USER_CC pgcc
##
set arch = x86          # 32 bit compile
#set arch = x86_64      # 64 bit compile

if ($arch == "x86_64") then
#  setenv NCROOT  /usr/local/netcdf-3.5.1
  setenv NCROOT  /usr/local/netcdf-3.5.0_pgi64
  setenv MPIROOT /usr/local/mpich-1.2.5.2_pgi64
else
  setenv NCROOT  /usr/local/netcdf-3.5.0_pgi32
  setenv MPIROOT /usr/local/mpich-1.2.5.2_pgi32
endif

setenv LIB_NETCDF $NCROOT/lib
setenv INC_NETCDF $NCROOT/include
setenv LIB_MPI    $MPIROOT/lib
setenv INC_MPI    $MPIROOT/include
alias mpirun $MPIROOT/bin/mpirun


#=======================================================================
## Do our best to get sufficient stack memory
limit stacksize unlimited

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
#set camroot      = /fs/cgd/data0/$LOGNAME/cam2_0_2_devNN/cam1
set camroot      = /data2/csm/cam3_0/cam1

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
#setenv CSMDATA     /fs/cgd/csm/inputdata
setenv CSMDATA     /data2/csm/inputdata

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
#set wrkdir       = /ptmp/$LOGNAME
set wrkdir       = /scratch/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld

## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

#=======================================================================
## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1

    set parops = ' ' 
    if ($smp  == "TRUE") set parops = "$parops -smp"
    if ($spmd == "TRUE") set parops = "$parops -spmd"

    $cfgdir/configure $parops -res $res     || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j $nthreads >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

#=======================================================================
## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -case $case -runtype $runtype -o $rundir/namelist \
 -namelist "&camexp nelapse=$nelapse mss_irt=0 $user_namelist /"  || echo "build-namelist failed" && exit 1

#=======================================================================
## Run CAM
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CAM in $rundir"
setenv MPSTKZ 400M
if ($spmd == 'TRUE') then
  echo "Running in SPMD mode (npe = $nthreads)"
  setenv OMP_NUM_THREADS 1
  mpirun -np $nthreads $blddir/cam < namelist  || echo "CAM run failed" && exit 1
else
  echo "Running in SMP (or serial) mode (threads = $nthreads)"
  setenv OMP_NUM_THREADS $nthreads
  $blddir/cam < namelist  || echo "CAM run failed" && exit 1
endif

exit 0
