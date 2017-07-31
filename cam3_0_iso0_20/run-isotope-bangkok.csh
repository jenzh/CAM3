#! /bin/tcsh -f
#
#=======================================================================
#
#  run-isotope-bangkok.csh
#
#  Generic batch submission script for PC-linux using PBS.  
#  MODIFIED for bangkok.cgd.ucar.edu by Andrew Gettelman
#  MODIFIED for isotope model by David Noone
#
#-----------------------------------------------------------------------
# Batch options for machine with PBS batch system (bangkok)
#
#
# Usage for Lahey compiler:
#   env OPT_BLD=lahey qsub run-pc.csh
#
# Usage for pgf90 compiler with pgcc: 
#   env OPT_BLD=pgf90-pgcc qsub run-pc.csh
#
# Usage for pgf90 compiler with gcc: 
#   env OPT_BLD=pgf90-gcc qsub run-pc.csh
#
# Usage for interactive queue:
#   qsub -I -l nodes=$nproc run-pc.csh
#
#-----------------------------------------------------------------------
#
### Run name (for sterr and stout files)
#PBS -N camiso
#PBS -o pbs_camiso.out
#PBS -e pbs_camiso.err
### Name of the queue (bangkok: small medium long verylong monster)
#PBS -q long
### Maximum number of nodes/processors
###PBS -l nodes=4:ppn=2
###PBS -l nodes=2:ppn=2VV
#PBS -l nodes=8:ppn=2
### Export all Environment variables
#PBS -V
# End of options
#=======================================================================
## Do our best to get sufficient stack memory
limit stacksize unlimited
#
# Set version of isotope code
#
set isovers = 'cam3_0_iso0_19'
#
#
set OS = `uname -s`;
switch ( $OS )
  case Linux:
     if ( ! $?PBS_JOBID ) then
       echo "${0}: ERROR::  This batch script must be submitted via PBS";
       echo "${0}:          on a Linux machine\!";
       exit;
     else
       echo "${0}: Running CAM on Linux using PBS";
     endif
     set job_id = `echo ${PBS_JOBID} | cut -f1 -d'.'`
     echo "${0}:  Set job_id to $job_id";
#
# Set a default build option
# pgf90-gcc, lahey, any other string defaults to pgf90-pgcc
#
     if ( ! $?OPT_BLD ) then
       setenv OPT_BLD pgf90-pgcc
     endif

     if ($OPT_BLD == "lahey") then
       setenv PATH .:/usr/local/lf95/bin:$PATH
       setenv LD_LIBRARY_PATH /usr/local/lf95/lib
  
       echo "${0}:  USING lf95 COMPILER"
       echo "${0}:  which lf95 = `which lf95`"
       setenv USER_FC lf95

       setenv NETCDF_ROOT /usr/local/netcdf-gcc-g++-lf95
       setenv MPI_ROOT    /usr/local/mpich-gcc-g++-lf95

     else		# assume pgf90... but which c compiler?
       if ( $OPT_BLD == "pgf90-gcc" ) then
         echo "${0}:  USING pgf90 COMPILER WITH gcc"
         setenv USER_CC gcc
         echo "${0}:  Set USER_CC to $USER_CC";
       else
         echo "${0}:  USING pgf90 COMPILER WITH pgcc"
       endif

       setenv PATH .:/usr/local/pgi-hpf-cc-5.1-6/linux86/bin:$PATH
       setenv LD_LIBRARY_PATH /usr/local/pgi/linux86/lib

       setenv NETCDF_ROOT /usr/local/netcdf-pgi-hpf-cc/
       setenv MPI_ROOT    /usr/local/mpich-1.2.5-pgi-hpf-cc-5.1-3
 
     endif
     breaksw;
  default:
    echo "${0}:  Use default values for number of nodes and shared memory CPUs";    exit;
endsw

#
# Set the environment
#
setenv LIB_NETCDF $NETCDF_ROOT/lib
setenv INC_NETCDF $NETCDF_ROOT/include
setenv MOD_NETCDF $INC_NETCDF
setenv LIB_MPI    $MPI_ROOT/lib
setenv INC_MPI    $MPI_ROOT/include
set    mpirun =   $MPI_ROOT/bin/mpirun

#
# Report environment
#
echo "${0}:  Set LIB_NETCDF to $LIB_NETCDF";
echo "${0}:  Set INC_NETCDF to $INC_NETCDF";
echo "${0}:  Set MOD_NETCDF to $MOD_NETCDF";
echo "${0}:  Set LIB_MPI to $LIB_MPI";
echo "${0}:  Set INC_MPI to $INC_MPI";
echo "${0}:  Set mpirun to $mpirun";

#
# Number of SPMD processes (let openmp defaults set SMP threads)
# *** match this with nodes in PBS at top
#
cd ${PBS_O_WORKDIR};
set nproc = `cat $PBS_NODEFILE | wc -l`
echo NUMBER OF NODES: $nproc
echo 'NODES'
cat  "$PBS_NODEFILE" 
echo

#
# Set needed namelist additions
#
set user_namelist = "trace_water=.true., wisotope=.true., readtrace=.false., reset_csim_iceprops=.true."
set user_namelist = "${user_namelist}, mfilt=1"

## Set resolution of model
set dycore = 'eul'
#set res = "256x512"            # T159 (ERA40)
#set res = "128x256"            # T85 (IPCC)
#set res = "64x128"             # T42 (DEFAULT)
#set res = '48x96'              # T31 (Paleoclimate/BGC)
set res = "32x64"              # T21
#set res = "8x16"               # T5 (REALLY fast for testing)

#set dycore = 'fv'
#set res = "2x2.5"              # normal (about T42)
#set res = "4x5"                # slower fastish  (about T21)
#set res = "10x15"              # fastish (25 x 27 point, about T9)

## Set number of tracers
#set nadv = 12   # Qx3 + H2Ox3 + HDOx3 + H218Ox3
set nadv = 18   # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + H217Ox3 + HTOx3

## ROOT OF CAM DISTRIBUTION - probably needs to be customized.
## Contains the source code for the CAM distribution.
## (the root directory contains the subdirectory "models")
set camroot      = /fs/cgd/data0/dcn/isotope/cam3_0/${isovers}/cam1

## ROOT OF CAM DATA DISTRIBUTION - needs to be customized unless running at NCAR.
## Contains the initial and boundary data for the CAM distribution.
## (the root directory contains the subdirectories "atm" and "lnd")
setenv CSMDATA     /fs/cgd/csm/inputdata

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: initial, restart, or branch.
## $nelapse is the number of timesteps to integrate, or number of days if negative.
set case         = camiso
set runtype      = initial
#set nelapse      = -488
#set nelapse      = -63
#set nelapse      = -1218
set nelapse      = -1
set histf        = -24

## $wrkdir is a working directory where the model will be built and run.
## $blddir is the directory where model will be compiled.
## $rundir is the directory where the model will be run.
## $cfgdir is the directory containing the CAM configuration scripts.
set wrkdir       = /scratch/cluster/$LOGNAME
set blddir       = $wrkdir/$case/bld
set rundir       = $wrkdir/$case
set cfgdir       = $camroot/models/atm/cam/bld


## Ensure that run and build directories exist
mkdir -p $rundir                || echo "cannot create $rundir" && exit 1
mkdir -p $blddir                || echo "cannot create $blddir" && exit 1

## If an executable doesn't exist, build one.
if ( ! -x $blddir/cam ) then
    cd $blddir                  || echo "cd $blddir failed" && exit 1
## hack for USER SOURCE CODE
    $cfgdir/configure -spmd  \
                      -nadv $nadv \
                      -dyn $dycore \
                      -res $res     || echo "configure failed" && exit 1

    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j $nproc >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -case $case -runtype $runtype -o $rundir/namelist \
 -namelist "&camexp nelapse=$nelapse nhtfrq=$histf mss_irt=0 $user_namelist /"  || echo "build-namelist failed" && exit 1

## Run CAM
cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CAM in $rundir"
$mpirun -np $nproc $blddir/cam < namelist >& $case.log  || echo "CAM run failed" && exit 1

exit 0
