#! /bin/tcsh -f
#
#=======================================================================
#
# run-isotope-catfish.csh
#
# Run isotopic cam3 on catfish, which does not use PBS.
# David Noone <dcn@colorado.edu> - Thu Jul  1 10:06:38 MDT 2004
#
#
#=======================================================================
##nice +10
#=======================================================================
set norun = 0			# compile one

set isovers = 'cam3_0_iso0_20'

#setenv DEBUG TRUE	# fails esmf build with defailt settings

## Default namelist settings:
## $case is the case identifier for this run. It will be placed in the namelist.
## $runtype is the run type: initial, restart, or branch.
## $nelapse is the number of timesteps to integrate, or number of days if negative.

set case         = cami20		# < ------ *** CHANGE THIS ***
set runtype      = initial
#set nelapse      = -123			# -123 to stop end year1
#set nelapse      = -488
#set nelapse      = -1218
set nelapse      = -1	


## Additional namelist settings
set user_namelist = "trace_water=.true., wisotope=.true., readtrace=.false., reset_csim_iceprops=.true., trace_gas=.true."
set user_namelist = "${user_namelist}, mfilt=1"
#set user_namelist = "${user_namelist}, nhtfrq=-24"
# Single column continuous output - Boulder, Darwin, Casey, Ascension Island
#set user_namelist = "${user_namelist}, fincl1lonlat='255e_40n','130e_0n','111e_66s','346e_8s'"

## Set resolution of model
set dycore = 'eul'
#set res = "256x512"            # T159 (ERA40)
#set res = "128x256"            # T85 (IPCC)
set res = "64x128"             # T42 (DEFAULT)
#set res = '48x96'              # T31 (Paleoclimate/BGC)
#set res = "32x64"              # T21
#set res = "8x16"               # T5 (REALLY fast for testing)

#set dycore = 'fv'
#set res = "10x15"               # fastish (about T10 - like the old GISS 8x10)
#set res = "4x5"                # fastish (about T21)
#set res = "2x2.5"              # default (about T42)

## Set number of tracers
#set nadv = 12   # Qx3 + H2Ox3 + HDOx3 + H218Ox3
#set nadv = 18   # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + H217Ox3 + HTOx3
set nadv = 16    # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + (CH4, N2O, CFC11, CFC12)
#set nadv = 22   # Qx3 + H2Ox3 + HDOx3 + H218Ox3 + H217Ox3 + HTOx3 + (CH4, N2O, CFC11, CFC12)

## Set parallel stratergy
set smp  = FALSE
#set spmd = FALSE
set spmd = TRUE
#
# Set number of nodes (SPMD), or threads (SMP) (not used with serial code)
#
set nodes = 4		# default 2 proc per node, so nodes=2 means 8 cpus
set nthreads = 1	# crashes with OMP

## Local paths and compiler settings
##
## Home of netcdf, fortran compiler (pgi by default, but let's be explict)
##
setenv USER_FC pgf90
setenv USER_CC pgcc
##
#set arch = x86          # 32 bit compile
set arch = x86_64      # 64 bit compile

if ($arch == "x86_64") then
#  setenv NCROOT  /usr/local/netcdf-3.5.1
  setenv NCROOT  /usr/local/netcdf-3.5.0_pgi64
  setenv MPIROOT /usr/local/mpich-1.2.6_pgi64
else
  setenv NCROOT  /usr/local/netcdf-3.5.0_pgi32
  setenv MPIROOT /usr/local/mpich-1.2.6_pgi32
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
#set camroot      = /data2/csm/cam3_0/cam1
set camroot      = /data2/dcn/models/CAM3/isotope/${isovers}/cam1

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

    $cfgdir/configure $parops \
                      -ocn som \
                      -nadv $nadv \
                      -dyn $dycore \
                      -res $res     || echo "configure failed" && exit 1
    echo "building CAM in $blddir ..."
    rm -f Depends
    gmake -j 2 >&! MAKE.out      || echo "CAM build failed: see $blddir/MAKE.out" && exit 1
endif

#=======================================================================
## Create the namelist
cd $blddir                      || echo "cd $blddir failed" && exit 1
$cfgdir/build-namelist -s -case $case -runtype $runtype -o $rundir/namelist \
 -namelist "&camexp nelapse=$nelapse mss_irt=0 $user_namelist /"  || echo "build-namelist failed" && exit 1

#=======================================================================
## Run CAM
if ($norun > 0) exit		# early bailout

cd $rundir                      || echo "cd $rundir failed" && exit 1
echo "running CAM in $rundir"

if ($spmd == 'TRUE') then
  echo "Prepaging for cam run with the pbs for spmd mode ($nodes nodes, 2 proc per node)"
  cat - > $case.job << EOF
#\!/bin/csh
#-----------------------------------------------------------------------
# Runs CAM vi the pbs. Type: qsub $case.job
#-----------------------------------------------------------------------
#PBS -q workq
#PBS -l nodes=${nodes}:ppn=2
#PBS -N $case
#PBS -o $case.out
#PBS -j oe
#PBS -m e
#PBS -V
#-----------------------------------------------------------------------
setenv OMP_NUM_THREADS $nthreads
cd $rundir
set nproc = \`cat \$PBS_NODEFILE | wc -l\`
mpirun -nolocal -np \$nproc -machinefile \$PBS_NODEFILE bld/cam < namelist > runlog
EOF

  echo PBS script created. Now YOU need to submit it.
  echo "TYPE: qsub $rundir/$case.job"
else
  echo "Running in SMP (or serial) mode (threads = $nthreads)"
  setenv OMP_NUM_THREADS $nthreads
  $blddir/cam < namelist  || echo "CAM run failed" && exit 1
endif

exit 0
