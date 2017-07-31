#include <misc.h>
#include <params.h>

module runtime_opts

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for reading CAM namelist camexp 
!          and broadcasting namelist values if needed.  
! 
! Author:
!   Original routines:  CMS
!   Module:             T. Henderson, September 2003
!
! $Id$
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use history
   use pspect
   use shr_orb_mod
   use units
   use constituents, only: pcnst, readtrace
   use chemistry,    only: trace_gas
   use soxbnd, only: scenario_prognostic_sulfur, rampyear_prognostic_sulfur
   use ghg_surfvals, only: scenario_ghg, rampYear_ghg, &
      ch4vmr, n2ovmr, f11vmr, f12vmr, co2vmr, ramp_co2_start_ymd, &
      ramp_co2_annual_rate, ramp_co2_cap
   use tracers, only: tracers_flag
   use time_manager, only: calendar, dtime, nestep, nelapse,      &
      start_ymd, start_tod, stop_ymd, stop_tod, ref_ymd, ref_tod, &
      perpetual_run, perpetual_ymd, tm_aqua_planet
   use filenames, only: nrevsn, ncdata, bndtvs, bndtvo, bndtvaer, &
      absems_data, bndtvg, aeroptics, bndtvvolc, bndtvcarbonscale,&
      bndtvsf6,&
      mss_wpass, rest_pfile, mss_irt, caseid, init_filepaths,     &
      get_archivedir, isccpdata,                                  &
      co_emis, bndtvdms, soil_erod, bndtvoxid, bndtvsox,               &
      brnch_retain_casename
   use restart, only: set_restart_filepath
#if ( ! defined COUP_CSM )
   use ice_dh, only: prognostic_icesnow,reset_csim_iceprops, icemodel_is
#endif
   use prescribed_aerosols, only: radforce, strat_volcanic,       &
      sulscl_rf, carscl_rf, ssltscl_rf, dustscl_rf, volcscl_rf,   &
      sulscl, carscl, ssltscl, dustscl, volcscl,                  &
      bgscl_rf, tauback, scenario_carbon_scale,                   &
      scenario_prescribed_sulfur, rampyear_prescribed_sulfur,     &
      prescribed_sulfur
   use cloudsimulator, only: doisccp
   use dycore, only: dycore_is
   use abortutils, only: endrun
   use ramp_scon, only: bndtvscon
   use ghg_surfvals, only: bndtvghg
   use water_tracers, only: trace_water
   use water_isotopes, only: wisotope

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make the default access private
   save


!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public runtime_options    ! Set and/or get all runtime options


!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! SOMEWHAT ALPHABETICAL listing of variables in the camexp namelist:
!
! variable                description
! --------             -----------------
!
! calendar             Calendar to use in date calculations.  'no_leap' (default) or 'gregorian'
!
! ctitle               Case title for header.
! 
! bndtvs               Path and filename of time-variant boundary
!                      dataset for sst's.
! 
! bndtvg               Path and filename of time-variant boundary 
!                      dataset for greenhouse loss rates.
!                      (required if trace_gas is set to true)
!
! bndtvsf6               Path and filename of time-variant boundary
!                      dataset for trtest emissions.
!		       (required if tracers_flag is set to true)
!
! bndtvo               Path and filename of time-variant boundary 
!                      dataset for ozone.
!
! bndtvaer             Path and filename of time-variant boundary
!                      dataset for aerosols.
!
! bndtvcarbonscale     Path and filename of time-variant boundary
!                      data of carbon scaling
!
! bndtvvolc            Path and filename of time-variant boundary
!                      dataset for stratospheric aerosol masses
!
! bndtvscon            Path and filename of time-variant boundary
!                      dataset for solar constant.
!
! bndtvghg             Path and filename of time-variant boundary
!                      dataset for greenhouse gas surface values.
!
! aeroptics            Path and filename of time-invariant 
!                      aerosol optics.
!
! co_emis              Path and filename of time-variant boundary 
!                      data set for fossil fuel carbon surface emissions.  
!
! bndtvdms             Path and filename of time-variant boundary 
!                      data set for DMS surface emissions.  
!
! soil_erod            Path and filename of time-variant boundary 
!                      data set for soil erodibility factors.  
!
! bndtvoxid                 Path and filename of time-variant boundary 
!                      data set for oxidants.  
!
! bndtvsox             Path and filename of time-variant boundary 
!                      data set for SOx surface emissions.  
!
! scenario_prognostic_sulfur 
!                      values can be 'FIXED' or 'RAMPED'
!                      sets so2,so4 surface flux
!                      FIXED  =>  not implemented (ends run)
!                      RAMPED =>  uses boundary data set bndtvsox
!                      Default: RAMPED
!
! rampyear_prognostic_sulfur
!                      Set to YYYY in order to cycle that year of sox emissions
!                      Default: not set ( does not cycle )
!
! prescribed_sulfur   'off', 'passive' or 'direct'
!                     default: 'direct'
!                     off is not implemented
!                     passive is an implicit method when prognostic is on
!                     direct means interacts with radiation code.
!
! scenario_prescribed_sulfur 
!                      values can be 'FIXED' or 'RAMPED'
!                      FIXED  =>  uses climatology
!                      RAMPED =>  not implemented
!
! rampyear_prescribed_sulfur
!                      Default: not set ( does not cycle )
!                      no other option is valid
!
! absems_data          Dataset with absorption and emissivity factors.
!
! aero_carbon          Set to .TRUE. to turn on carbon prognostic aerosols.  
   logical :: aero_carbon
! 
! aero_feedback_carbon     Set to .TRUE. to enable feedback of carbon
!                          prognostic aerosols.  
   logical :: aero_feedback_carbon
! 
! aero_sea_salt        Set to .TRUE. to turn on sea salt prognostic aerosols.  
   logical :: aero_sea_salt
! 
! aero_feedback_sea_salt   Set to .TRUE. to enable feedback of sea salt
!                          prognostic aerosols.  
   logical :: aero_feedback_sea_salt
! 
! prognostic_sulfur    "off", "passive", "direct"
!                      off = no prognostic sulfur (default)
!                      passive = prognostic sulfur, no radiative interaction
!                      direct = prognostic sulfur drive radiative interaction
   character(len=16) :: prognostic_sulfur
! 
! caseid               Case name for model run.  32 characters max.
!                      Included in mass store path name for history and
!                      restart files.
! 
! TBH:  Move declarations of dif2 and dif4 here from read_namelist() 
! TBH:  once dif2 and dif4 are made private in module comhd.  
!   real(r8) :: dif2
! dif2 = nnn.n,        del2 horizontal diffusion coeff. Default value 
!                      defined in module comhd.  
! 
!   real(r8) :: dif4
! dif4 = nnn.n,        del4 horizontal diffusion coeff. Default value 
!                      defined in module comhd.  
! 
! divdampn = 0.        Number of days (from nstep 0) to run divergence
!                      damper
!
! dtime = nnnn,        Model time step in seconds. Default is dycore dependent.
! 
! eccen                The eccentricity of the earths orbit to use (1.e36 to
!                      use the default -- defined as SHR_ORB_UNDEF_REAL).
!                      (Unitless typically 0 - 0.1)
! 
! eps = nnn.n,         time filter coefficient. Defaults to 0.06.
! 
! fincl1 = 'field1', 'field2',...
!                      List of fields to add to the primary history file.
! fincl1lonlat = 'longitude by latitude','longitude by latitude',...
!                      List of columns ('longitude_latitude') or contiguous 
!                      columns ('longitude:longitude_latitude:latitude') at 
!                      which the fincl1 fields will be output. Individual 
!                      columns are specified as a string using a longitude
!                      degree (greater or equal to 0.) followed by a single 
!                      character (e)ast/(w)est identifer, an
!                      underscore '_' , and a latitude degree followed by a 
!                      single character (n)orth/(s)outh identifier.
!                      example '10e_20n' would pick the model column closest
!                      to 10 degrees east longitude by 20 degrees north 
!                      latitude.  A group of contiguous columns can be 
!                      specified by using lon lat ranges with their single
!                      character east/west or north/south identifiers
!                      example '10e:20e_15n:20n'.  Would outfield all 
!                      fincl1 fields at the model columns which fall
!                      with in the longitude range from 10 east to 20 east
!                      and the latitude range from 15 north to 20 north
!
! fincl[2..6] = 'field1', 'field2',...
!                      List of fields to add to the auxiliary history file.
!
! fincl2..6]lonlat = 'longitude by latitude','longitude by latitude',...
!                      List of columns ('longitude_latitude') or contiguous 
!                      columns ('longitude:longitude_latitude:latitude') at 
!                      which the fincl[2..6] fields will be output. Individual 
!                      columns are specified as a string using a longitude
!                      degree (greater or equal to 0.) followed by a single 
!                      character (e)ast/(w)est identifer, an
!                      underscore '_' , and a latitude degree followed by a 
!                      singel character (n)orth/(s)outh identifier.
!                      example '10e_20n' would pick the model column closest
!                      to 10 degrees east longitude by 20 degrees north 
!                      latitude.  A group of contiguous columns can be 
!                      specified by using lon lat ranges with their single
!                      character east/west or north/south identifiers
!                      example '10e:20e_15n:20n'.  Would outfield all 
!                      fincl[2..6] fields at the model columns which fall
!                      with in the longitude range from 10 east to 20 east
!                      and the latitude range from 15 north to 20 north
!
! fexcl1 = 'field1','field2',... 
!                      List of field names to exclude from default
!                      primary history file (default fields on the 
!                      Master Field List).
! 
! fexcl[2..6] = 'field1','field2',... 
!                      List of field names to exclude from
!                      auxiliary history files.
! 
! fhstpr1 = 'field1', 'field2',...
!                      List of fields to change buffer size in
!                      primary history file
!
! fhstpr[2..6] = 'field1', 'field2',...
!                      List of fields to change buffer size in auxiliary files
!
! fwrtpr1 = 'field1', 'field2',...
!                      List of fields to change output data type in
!                      primary history file
!
! fwrtpr[2..6] = 'field1', 'field2',...
!                      List of fields to change output data type in
!                      auxiliary files
!
! iradae = nnn,        frequency of absorp/emis calc in time steps
!                      (positive) or hours (negative).
! 
! iradlw = nnn,        frequency of longwave rad. calc. in time steps
!                      (positive) or hours (negative).
! 
! iradsw = nnn,        freq. of shortwave radiation calc in time steps
!                      (positive) or hours (negative).
! 
! mss_irt              Mass Store retention time for history files
!                      in days.
! 
! itsst = nnn,         frequency of SST update in time steps
! 
! TBH:  Move declaration of kmxhdc here from read_namelist() 
! TBH:  once kmxhdc is made private in module comhd.  
! kmxhdc = nn          number of levels (starting from model top) to
!                      apply Courant limiter.  Default value defined 
!                      in module comhd.  
! 
! mfilt = nn,nn,nn     Array containing the maximum number of time 
!                      samples per disk history file. Defaults to 5.
!                      The first value applies to the primary hist. file,
!                      the second to the first aux. hist. file, etc.
! 
! mvelp                The longitude of vernal equinox of the earths orbit to 
!                      use (1.e36 to use the default -- defined as 
!                      SHR_ORB_UNDEF_REAL).  (0-360 degrees')
! 
! ncdata               Path and filename of initial condition dataset.
! 
! nelapse = nnn,       Specify the ending time for the run as an interval
!                      starting at the current time in either timesteps
!                      (if positive) or days (if negative).
!                      Either nestep or (stop_ymd,stop_tod) take precedence.
! 
! nestep = nnnn,       Specify the ending time for the run as an interval
!                      starting at (start_ymd,start_tod) in either timesteps
!                      (if positive) or days (if negative).
!                      (stop_ymd,stop_tod) takes precedence if set.
! 
! nhtfrq = nn,nn,nn,.. Output history frequency for each tape
!
!                      If = 0 : monthly average
!                      If > 0 : output every nhtfrq time steps.
!                      If < 0 : output every abs(nhtfrq) hours.
! 
! nlvdry = nn,         Number of layers over which to do dry
!                      adjustment. Defaults to 3.
! 
! nrefrq = nn,         Frequency of restart dataset writes. 
!                      For non-flux coupled runs, restart files are 
!                      written and disposed for every dispose of the 
!                      primary history file. If this variable is 0, then 
!                      no restart are written.
!                      NOTE: NOW DUE TO NEW LSM: THIS VARIABLE CAN 
!                      ONLY BE 1 or 0. 
!                      For flux coupled runs, insist that restart files
!                      are written
! 
! nrevsn               Filename of dataset to branch from (nsrest=3)
!                      Full pathname of dataset required.
! 
!------------------------------------------------------------------
! The following 5 are specific to f-v dynamics (see dynpkg for info)
!------------------------------------------------------------------
! nsplit               Lagrangian time splits for Lin-Rood.
! iord                 scheme to be used for E-W transport (default: 4)
! jord                 scheme to be used for N-S transport (default: 4)
! kord                 scheme to be used for vertical mapping (default: 4)
! use_eta              flag to use ETA values from dynamics/lr/set_eta.F90
!                      Default is .false. (use eta values from IC)
!------------------------------------------------------------------
! 
!------------------------------------------------------------------
! The following 7 are specific to f-v decomposition and transposes 
! (see spmd_dyn for info)
!------------------------------------------------------------------
! npr_yz(4)            yz and xy decompositions
   integer :: npr_yz(4)
! geopktrans           geopotential method (routine geopk)
   integer :: geopktrans
! tracertrans          number of simultaneously transposed tracers
   integer :: tracertrans
! ompnest              option for nested openmp
   integer :: ompnest
! force_2d             option to force transpose computation for 1D decomp.
   integer :: force_2d
! modcomm_transpose    mod_comm transpose method (varies with mpi/mpi2 choice)
   integer :: modcomm_transpose
! modcomm_geopk        mod_comm geopk method (varies with mpi/mpi2 choice)
   integer :: modcomm_geopk
!------------------------------------------------------------------
! The following 2 are specific to eul/sld communication algorithms
! (see { eul | sld }/spmd_dyn for info)
!------------------------------------------------------------------
! dyn_alltoall         dynamics transpose option.
   integer :: dyn_alltoall
! dyn_allgather        dynamics gather option.
   integer :: dyn_allgather
!------------------------------------------------------------------
! The following 2 are specific to the swap communication module, used
! in the point-to-point implementations of eul/sld and physics
! communication algorithms (see swap_comm for info)
!------------------------------------------------------------------
! swap_comm_order      Performance tuning option for swap communication.
   integer :: swap_comm_order
! swap_comm_protocol   Performance tuning option for swap communication.
   integer :: swap_comm_protocol
!------------------------------------------------------------------
!
! nsrest               Code for type of run: 0=initial, 1=restart,
!                      or 3=branch
! 
! archive_dir          Archive directory name
!
! hfilename_spec       Flexible filename specifier for history files
!
! rest_pfile           Name of Restart Pointer file
! 
! mss_wpass            Write password for model output files.
! 
! ozncyc = .T.,        If false, do not cycle ozone dataset(assume
!                      multiyear)
!
! obliq                The obliquity of the earths orbit to use (1.e36 to
!                      use the default -- defined as SHR_ORB_UNDEF_REAL). 
!                      (Degree's)
!
! perpetual_run = .F.  Set to .true. to specify that the run will use a perpetual
!                      calendar.  If perpetual_ymd is not set then read the perpetual
!                      date from the initial file.
!
! perpetual_ymd        Perpetual date specified as (year*1000 + month*100 + day).
!                      This date overrides the date from the initial file.
!                      If aqua_planet=.true. then perpetual_ymd is ignored and the
!                      perpetual date is set to 321.
! 
! pertlim = n.n        Max size of perturbation to apply to initial
!                      temperature field.
!
! phys_alltoall        Dynamics/physics transpose option. See phys_grid module.
!
   integer :: phys_alltoall
! 
! phys_loadbalance     Load balance option for performance tuning of 
!                      physics chunks.  See phys_grid module.  
   integer :: phys_loadbalance
! 
! phys_chnk_per_thd    Performance tuning option for physics chunks.  See 
!                      phys_grid module.  
   integer :: phys_chnk_per_thd
! 
! ref_ymd              Reference date for time coordinate encoded in yearmmdd format.
!                      Default value is start_ymd.
!
! ref_tod              Reference time of day for time coordinate in seconds since 0Z.
!                      Default value is start_tod.
!
! sstcyc = .T.,        If false, do not cycle sst dataset(assume
!                      multiyear)
! 
! logical reset_csim_iceprops = .F.,
!
!                    ! if true => resets the csim ice properties to base state
!                    ! No Snow Cover, TSICE and TS1-4 are all set to
!                    ! freezing. Default is false.
!                    ! The csim is sensitive to imbalances between the
!                    ! surface temperature and ice temperatures. When
!                    ! using an initial conditions dataset interpolated
!                    ! from a different resolution you may have to set this
!                    ! to true to get csim to run.  If set to true you will
!                    ! have to allow time for the ice to "spin-up".
!
! start_ymd            Starting date for run encoded in yearmmdd format.  Default value
!                      is read from initial conditions file.
!
! start_tod            Starting time of day for run in seconds since 0Z.  Default value
!                      is read from initial conditions file.
!
! stop_ymd             Stopping date for run encoded in yearmmdd format.  No default.
!
! stop_tod             Stopping time of day for run in seconds since 0Z.  Default: 0.
!
! adiabatic = .F.      Don't call physics
!
! ideal_phys = .F.     Only run the "idealized" dynamical core
!                      (dynamics + specified physics) of the model.
!
! aqua_planet = .F.    Run in "aqua_planet" mode.  Physics remains on but is run for
!                      perpetual vernal equinox conditions; phis = 0; ocean
!                      everywhere - no land and no sea-ice; SST's specified analytically
!
! flxave = .T.         If true, only send data to the flux coupler on
!                      radiation time steps. This namelist variable is
!                      only used when running through the flux coupler.
!
! precc_thresh         Precipitation threshold to use for PRECCINT and PRECCFRQ (mm/hr)
!                      Defaults to 0.1.
!
! precl_thresh         Precipitation threshold to use for PRECLINT and PRECLFRQ (mm/hr)
!                      Defaults to 0.05.
!
! trace_gas = .F.      If true, turn on greenhouse gas code for
!                      CH4, N2O, CFC11 and CFC12 . (Must add 4 to pcnst)
!
! trace_water = .F.    If true, turn on nwater tracers (use with wisotope)
! wisotope = .F.       If true, define water tracers to be water isotopes
!
! tracers_flag = .F.    If true, implement tracer test code. Number of tracers determined
!                      in tracers_suite.F90 must agree with PCNST in params.h
!
! readtrace = .T.      If true, tracer initial conditions obtained from 
!                      initial file. 
!
! co2vmr               global       co2 volume mixing ratio
! ch4vmr               tropospheric ch4 volume mixing ratio
! n2ovmr               tropospheric n2o volume mixing ratio
! f11vmr               tropospheric f11 volume mixing ratio
! f12vmr               tropospheric f12 volume mixing ratio
!
! iyear_AD             The year AD to calculate the orbital parameters for.  
!                      By default this is set to 2000000000 (defined to SHR_ORB_UNDEF_INT) 
!                      which means use the input values o: eccen, obliq and mvelp.
!
! inithist             Generate initial dataset as auxillary history file
!                      can be set to '6-HOURLY', 'DAILY', 'MONTHLY', 'YEARLY' or 'NONE'. 
!                      default: 'MONTHLY '
!
! prognostic_icesnow = .T,  prognostic snow over ice, currently limited to
!                      0.5m.  If this is false then a snow climatology
!                      is used (default .T.)
!
! linebuf              true => force buffer flush of stdout with each 
!                      newline generated (useful for debugging)
!
! empty_htapes         true => no fields by default on history tapes
!
! print_step_cost      true => print per timestep cost info
!
! avgflag_pertape      A, I, X, or M means avg, instantaneous, max or min for all fields on
!                      that tape
!
! scenario_ghg         values can be 'FIXED' or 'RAMPED' or 'RAMP_CO2_ONLY'
!                      sets co2,ch4,n2o,cfcf11,cfc12 volume mixing ratios
!                      FIXED => volume mixing ratios are fixed and are
!                      either have preset or namelist input values
!                      RAMPED => volume mixing ratios are ramped
!                      RAMP_CO2_ONLY => only co2 mixing ratios are ramped
!                      DEFAULT: FIXED 
!
! ramp_co2_start_ymd     date on which ramping of co2 begins; REQUIRED to be set 
!                        for scenario_ghg='RAMP_CO2_ONLY'
! ramp_co2_annual_rate   percentage amount of co2 ramping per yr; default is 1.0 
! ramp_co2_cap           co2 ramp cap if rate>0, floor otherwise; 
!                        specified as multiple or fraction of inital value;
!                        ex. 4.0 => will cap at 4x initial co2 setting;
!                        default is boundless if rate>0, zero otherwise
!
! doisccp              whether to do ISCCP calcs and history output (default false)
!
   character*16 scenario_scon
!                    ! values can be 'FIXED' or 'RAMPED'
!                    ! FIXED => scon is fixed and can either have preset or
!                    ! namelist value
!                    ! RAMPED => scon is ramped
!                    ! DEFAULT => FIXED
!
! rampYear_ghg         ramped gases fixed at this year if set to a value
!                      greater than zero.  Default value is 0.
!
   integer rampYear_scon
!                    ! ramped scon fixed at this year if set to a value
!                    ! greater than zero.  Default value is 0.
!
!   logical indirect     
!                    ! true => include indirect radiative effects of
!                    ! sulfate aerosols.  Default is false.
!
! radforce             Compute forcing from aerosols (Default is false)
!
! strat_volcanic       Use stratospheric volcanic aerosols masses and
!                      couple with radiative forcing computations
!
! scenario_carbon_scale
!                     'FIXED' or 'RAMPED'
!                      FIXED means use carscl
!                      RAMPED means use data from file bndtvcarbonscale
!
! sulscl_rf, carscl_rf, ssltscl_rf, dustscl_rf, bgscl_rf, volcscl_rf
!                      Set corresponding aerosols to 0.0 mmr
!                      for radiative forcing.  These do not affect
!                      mmr's used for climate integration.
!
! tauback              Optical depth of (rh = .8, sulfate-like) 
!                      background aerosol
! 
! sulscl, carscl, ssltscl, dustscl, volcscl
!                      Scale corresponding aerosols in 
!                      climatology by this amount for the
!                      purpose of the climate integration
!
! Define the camexp namelist
!
! TBH:  NOTE that the definition of camexp SHOULD APPEAR here, not 
! TBH:  inside read_namelist().  If it did, then we could easily 
! TBH:  write other methods (like a proposed method to dump the 
! TBH:  namlist to a log file) that use camexp.  However, before 
! TBH:  the definition can be moved outside of read_namelist(), 
! TBH:  common blocks in comctl.h, comtfc.h, comsol.h, 
! TBH:  comadj.h, and perturb.h must be converted to modules.  


!-----------------------------------------------------------------------
! Subroutines and functions --------------------------------------------
!-----------------------------------------------------------------------
contains


subroutine read_namelist

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read data from namelist camexp to define the run. Process some of the
! namelist variables to determine history and restart/branch file path 
! names.  Check input namelist variables for validity and print them
! to standard output. 
! 
! Method: 
! Important Note for running on SUN systems: "implicit automatic (a-z)"
! will not work because namelist data must be static.
!
! Author: 
! Original version:  CCM1
! Standardized:      L. Bath, June 1992
!                    T. Acker, March 1996
!     
!-----------------------------------------------------------------------
!
! $Id$
!
!-----------------------------------------------------------------------

   use infnan,       only: inf
   use string_utils, only: to_upper
   ! Note that the following interfaces are prototypes proposed by Henderson 
   ! and Eaton.  They minimize coupling with other modules.  Design of these 
   ! interfaces should be refined via review by other CAM developers.  
   ! Interface *_defaultopts() gets default values from the responsible 
   ! module (Expert) prior to namelist read.  
   ! Interface *_setopts() sends values to the responsible module (Expert) 
   ! after namelist read.  Erroneous values are handled by Experts.  
   ! TBH  9/8/03 
   use phys_grid, only: phys_grid_defaultopts, phys_grid_setopts
#if ( defined SPMD )
   use swap_comm, only: swap_comm_defaultopts, swap_comm_setopts
   use spmd_dyn, only: spmd_dyn_defaultopts, spmd_dyn_setopts
#endif
   use aerosol_intr, only: aerosol_defaultopts, aerosol_setopts
   use comhd, only: comhd_defaultopts, comhd_setopts
   use ramp_scon, only: rampnl_scon


#include <comadj.h>
#include <comctl.h>
#include <comtfc.h>
#include <perturb.h>
#include <comsol.h>

!-----------------------------------------------------------------------
   include 'netcdf.inc'
!
!---------------------------Local variables-----------------------------
! 
! TBH:  Move declarations of dif2, dif4, and kmxhdc to module data 
! TBH:  once dif2 and dif4 are made private in module comhd.  
   real(r8) :: dif2, dif4
   integer :: kmxhdc

   logical linebuf
   character(len=256) :: archive_dir = ''
!
#if ( defined SUNOS )
!
! Namelist variables may not be on the stack on SUN
!
   save linebuf, archive_dir
#endif
   data linebuf/.false./ ! Default: allow system to buffer stdout
!
   character ctemp*8      ! Temporary character strings
   integer ntspdy         ! number of timesteps per day
   integer t              ! history tape index
   integer lastchar       ! index to last char of a char variable
   integer ierr           ! error code

#if ( defined COUP_CSM )
   logical prognostic_icesnow,reset_csim_iceprops
#endif

!
! Define the camexp namelist
!
! TBH:  NOTE:  Move the definition of camexp outside of this routine 
! TBH:  as soon as common blocks in comctl.h, comtfc.h, 
! TBH:  comsol.h, comadj.h, and perturb.h have been converted to 
! TBH:  modules.  
!        
!
#if ( ! defined T3D )
!
! Disclaimer: The namelist items, nhstpr, fhstpr1-fhstpr6, fhstwrtpr1-fwrtpr6,
! ideal_phys, trace_gas, bndtvg, scenario_ghg, tracers_flag, bndtvsf6,
! scenario_scon, rampYear_ghg, and rampYear_scon
! are considered unsuported features. The code may not even run with
! these options and has NOT been verified to create correct science.
! As such these options should only be used with caution.
!
! If a namelist option is not mentioned in the CAM Users Guide, it may 
! not be supported.  
!       
  namelist /camexp/ ctitle  ,ncdata  ,bndtvs  ,bndtvo  , bndtvg , &
                    bndtvaer, bndtvvolc, aeroptics, bndtvcarbonscale,&
                    co_emis, bndtvdms, soil_erod, bndtvoxid, bndtvsox, &
                    scenario_prognostic_sulfur, rampyear_prognostic_sulfur, &
                    rest_pfile,mss_wpass,nsrest  ,mss_irt , archive_dir, &
                    nrevsn  ,nhstpr  ,ndens   ,nhtfrq  , &
                    nrefrq  ,mfilt   ,absems_data , &
                    fincl1  ,fincl2  ,fincl3  ,fincl4  ,fincl5  , &
                    fincl1lonlat,fincl2lonlat,fincl3lonlat, &
                    fincl4lonlat  ,fincl5lonlat  , &
                    fincl6  ,fexcl1  ,fexcl2  ,fexcl3  ,fexcl4  , &
                    fexcl5  ,fexcl6  ,hfilename_spec, &
                    fhstpr1 ,fhstpr2 ,fhstpr3 ,fhstpr4 ,fhstpr5 ,fhstpr6 , &
                    fwrtpr1 ,fwrtpr2 ,fwrtpr3, fwrtpr4 ,fwrtpr5 ,fwrtpr6 , &
                    calendar, dtime, nelapse, nestep, start_ymd, start_tod,  &
                    stop_ymd, stop_tod, ref_ymd, ref_tod, perpetual_run, &
                    perpetual_ymd,   precc_thresh, precl_thresh, &
                    eps     ,dif2    ,dif4    ,kmxhdc  ,iradsw  , &
                    iradlw  ,iradae  ,itsst   ,nlvdry  ,sstcyc  , &
                    ozncyc  ,pertlim ,divdampn,caseid  ,adiabatic,flxave , &
                    trace_gas, readtrace, &
                    trace_water, wisotope, &
                    tracers_flag, bndtvsf6, &
                    co2vmr  ,ch4vmr  ,n2ovmr  ,f11vmr  ,f12vmr  , &
                    obliq   ,eccen   ,mvelp   ,iyear_AD,scon    , &
                    inithist, linebuf ,ideal_phys, &
                    aqua_planet, indirect, nsplit, &
                    iord, jord, kord, use_eta, &
                    npr_yz, geopktrans, tracertrans, ompnest, &
                    force_2d, modcomm_transpose, modcomm_geopk, &
                    dyn_alltoall, dyn_allgather, &
                    swap_comm_order, swap_comm_protocol, &
                    scenario_ghg, scenario_scon, &
                    rampYear_ghg, rampYear_scon, empty_htapes, &
                    print_step_cost, avgflag_pertape,prognostic_icesnow, &
                    reset_csim_iceprops, som_conschk_frq, ice_conschk_frq, &
                    doisccp, isccpdata, radforce, &
                    strat_volcanic, scenario_carbon_scale, &
                    sulscl_rf, carscl_rf, ssltscl_rf, dustscl_rf, &
                    bgscl_rf, volcscl_rf, &
                    tauback, sulscl, carscl, ssltscl, dustscl, volcscl,&
                    scenario_prescribed_sulfur, rampyear_prescribed_sulfur, &
                    phys_alltoall, phys_loadbalance, phys_chnk_per_thd, &
                    prognostic_sulfur, &
                    prescribed_sulfur, &
                    aero_carbon, aero_feedback_carbon, &
                    aero_sea_salt, aero_feedback_sea_salt, &
                    brnch_retain_casename, bndtvscon, bndtvghg, &
                    ramp_co2_start_ymd, ramp_co2_annual_rate, ramp_co2_cap

#endif

! 
!-----------------------------------------------------------------------
!
! Preset scenario variables and ramping year
!
   scenario_scon = 'FIXED'
   rampYear_scon = 0
!
! Finite volume code only: Set Lagrangian time splits.  A default of zero indicates the number
! should be automatically computed unless the user enters something.
!
   nsplit = 0
   iord = 4
   jord = 4
   kord = 4
   use_eta = .false.        ! Use a's and b's from the initial file
!
! Preset sulfate aerosol related variables

   indirect  = .false.
! 
! Set anncyc true, no longer in namelist
! 
   anncyc = .true.

! 
! Get default values of runtime options for spmd_dyn
!
#if ( defined SPMD )
   if ( dycore_is ('LR') ) then
      call spmd_dyn_defaultopts(                 &
             npr_yz_out         =npr_yz,         &
             geopktrans_out     =geopktrans,     &
             tracertrans_out    =tracertrans,    &
             ompnest_out        =ompnest,        &
             force_2d_out       =force_2d,       &
             modcomm_transpose_out =modcomm_transpose, &
             modcomm_geopk_out     =modcomm_geopk)
   endif
   if ( dycore_is ('EUL') .or. dycore_is ('SLD') ) then
      call spmd_dyn_defaultopts(                 &
             dyn_alltoall_out   =dyn_alltoall,   &
             dyn_allgather_out  =dyn_allgather   )
   endif
! 
! Get default values of runtime options for swap module.
!
   call swap_comm_defaultopts(                       &
          swap_comm_order_out=swap_comm_order,       &
          swap_comm_protocol_out=swap_comm_protocol)
#endif
! 
! Get default values of runtime options for physics chunking.
!
   call phys_grid_defaultopts(                    &
          phys_loadbalance_out =phys_loadbalance, &
          phys_alltoall_out    =phys_alltoall,   &
          phys_chnk_per_thd_out=phys_chnk_per_thd)
! 
! Get default values of runtime options for prognostic aerosols
!
   call aerosol_defaultopts(                               &
          prognostic_sulfur_out     =prognostic_sulfur,    &
          aero_carbon_out           =aero_carbon,          &
          aero_feedback_carbon_out  =aero_feedback_carbon, &
          aero_sea_salt_out         =aero_sea_salt,        &
          aero_feedback_sea_salt_out=aero_feedback_sea_salt)
! 
! Get default values of runtime options for comhd
!
   call comhd_defaultopts(dif2_out  =dif2, &
                          dif4_out  =dif4, &
                          kmxhdc_out=kmxhdc)

   if (masterproc) then
!
! Read in the camexp namelist from standard input
!
      read (5,camexp,iostat=ierr)
      if (ierr /= 0) then
         write(6,*)'READ_NAMELIST: Namelist read returns ',ierr
         call endrun
      end if
! 
! Check CASE namelist variable
!
      if (caseid==' ') then
         call endrun ('READ_NAMELIST: Namelist variable CASEID must be set')
      end if

      lastchar = len(caseid)
      if (caseid(lastchar:lastchar) /= ' ') then
         write(6,*)'READ_NAMELIST: CASEID must not exceed ', len(caseid)-1, ' characters'
         call endrun
      end if
      icecyc = sstcyc    ! ice-cycling is tied to the sst-dataset
#ifndef COUP_CSM
!
! Data ice-model can not use prognostic snow-depth or reset the ice properties
!
      if ( icemodel_is('data') )then
         if ( .not. prognostic_icesnow ) &
            write(6,*) 'Warning: prognostic_icesnow for data-ice-model is always false'
         prognostic_icesnow = .false.
         if ( .not. reset_csim_iceprops ) &
            write(6,*) 'Warning: reset_csim_iceprops for data-ice-model is always false'
         reset_csim_iceprops = .false.
      end if
#endif
   end if
!
! Line buffer stdout if requested
!
   if (linebuf) then
!        call flush(6)
      call linebuf_stdout ()
   end if
!
! Precipitation thresholds (check range and convert to mm/hr)
!
   if ( precc_thresh < 0.0_r8 ) then
      call endrun ('READ_NAMELIST: PRECC threshold needs to be >= 0.0.')
   endif
   if ( precc_thresh > 9.99_r8 ) then
      call endrun ('READ_NAMELIST: PRECC threshold needs to be <= 9.99 mm/hr.')
   endif
   if ( precl_thresh < 0.0_r8 ) then
      call endrun ('READ_NAMELIST: PRECL threshold needs to be >= 0.0.')
   endif
   if ( precl_thresh > 9.99_r8 ) then
      call endrun ('READ_NAMELIST: PRECL threshold needs to be <= 9.99 mm/hr.')
   endif
   precc_thresh = precc_thresh/(1000.0*3600.0) ! convert to m/sec
   precl_thresh = precl_thresh/(1000.0*3600.0) ! convert to m/sec
#if ( defined SPMD )
   call distnl ( )
#endif

! Communicate to time manager (there should be a method for this).
   tm_aqua_planet = aqua_planet

! 
! Set continuation run flags
! 
   if (nsrest>0) then
      nlres  = .true.
   endif
   if (nsrest==2) then
      call endrun ('READ_NAMELIST: The regeneration option is no longer available')
   end if
   if (nsrest==3) then
      nlhst  = .true.
      lbrnch = .true.
   endif

#if ( defined COUP_CSM )
!
! Check that flxave occurs only if iradsw is gt 1
!
   if (flxave .and. iradsw==1 ) then
      call endrun ('READ_NAMELIST: iradsw must be greater that one if flux averaging option is enabled')
   endif
#endif
!++mv
!
! Determine ramping logic
!
   if (scenario_scon == 'FIXED') then
      doRamp_scon = .false.
   else if (scenario_scon == 'RAMPED') then
      doRamp_scon = .true.
   else
      call endrun ('READ_NAMELIST: SCENARIO_SCON must be set to either FIXED or RAMPED')
   endif
!       
! Initialize namelist related scon info
!
   if (doRamp_scon) then
      call rampnl_scon( rampYear_scon )
      if (masterproc) write(6,*)'scon set by ramp code'
   else
      if (masterproc) write(6,*)'scon set to fixed value of ',scon 
   endif
!
! Auxiliary history files:
! Store input auxf values in array aux (from common block /comhst/).
!
! If generate an initial conditions history file as an auxillary tape:
!
   ctemp = to_upper(inithist) 
   inithist = trim(ctemp)
   if (inithist /= '6-HOURLY' .and. inithist /= 'DAILY' .and. &
       inithist /= 'MONTHLY'  .and. inithist /= 'YEARLY') then
      inithist = 'NONE'
   endif
!
! Ensure that monthly averages have not been specified for aux. tapes
!
   do t=2,ptapes
      if (nhtfrq(t) == 0) then
         call endrun ('READ_NAMELIST: Only the primary history file may be monthly averaged')
      end if
   end do
! 
! History file write up times
! Convert write freq. of hist files from hours to timesteps if necessary.
! 
   do t=1,ptapes
      if (nhtfrq(t) < 0) then
         nhtfrq(t) = nint((-nhtfrq(t)*3600.)/dtime)
      end if
   end do
!
! Initialize the filename specifier if not already set
! This is the format for the history filenames:
! %c= caseid, %t=tape no., %y=year, %m=month, %d=day, %s=second, %%=%
! See the filenames module for more information
!
   do t = 1, ptapes
      if ( len_trim(hfilename_spec(t)) == 0 )then
         if ( nhtfrq(t) == 0 )then
            hfilename_spec(t) = '%c.cam2.h%t.%y-%m.nc'        ! Monthly files
         else
            hfilename_spec(t) = '%c.cam2.h%t.%y-%m-%d-%s.nc'
         end if
      end if
      if ( masterproc ) then
         write(6,*) 'Filename specifier for tape ', t, ' = ', &
                    trim(hfilename_spec(t))
      end if
   end do
!
! Only one time sample allowed per monthly average file
! 
   if (nhtfrq(1) == 0) mfilt(1) = 1
!
! Check validity of per-tape averaging flag
!
   do t=1,ptapes
      if (avgflag_pertape(t) /= ' ') then
         if (avgflag_pertape(t) == 'A' .or. avgflag_pertape(t) == 'I' .or. &
             avgflag_pertape(t) == 'X' .or. avgflag_pertape(t) == 'M') then
            write(6,*)'Unless overridden by namelist input on a per-field basis (FINCL),'
            write(6,*)'All fields on history file ',t,' will have averaging flag ',avgflag_pertape(t)
         else
            write(6,*)'Invalid per-tape averaging flag specified:', avgflag_pertape(t)
            call endrun ('READ_NAMELIST')
         end if
      end if
   end do
! 
! Convert iradsw and iradlw from hours to timesteps if necessary
! 
   if (iradsw < 0) iradsw = nint((-iradsw*3600.)/dtime)
   if (iradlw < 0) iradlw = nint((-iradlw*3600.)/dtime)
! 
! Convert iradae from hours to timesteps if necessary and check that
! iradae must be an even multiple of iradlw
! 
   if (iradae < 0) iradae = nint((-iradae*3600.)/dtime)
   if (mod(iradae,iradlw)/=0) then
      write(6,*)'READ_NAMELIST:iradae must be an even multiple of iradlw.'
      write(6,*)'     iradae = ',iradae,', iradlw = ',iradlw
      call endrun
   end if
! 
! Do absorptivities/emissivities have to go on a restart dataset?
! 
   ntspdy = nint(86400./dtime) ! no. timesteps per day
   if (nhtfrq(1) /= 0) then
      if (masterproc .and. mod(nhtfrq(1),iradae)/=0) then
         write(6,*)'READ_NAMELIST: *** NOTE: Extra overhead invoked putting',  &
            ' a/e numbers on restart dataset. ***   ',         &
            ' To avoid, make mod(nhtfrq,iradae) = 0'
      end if
   else
      if (masterproc) then
         if (mod(ntspdy,iradae) /= 0 .or. iradae > ntspdy) then
            write(6,*)'READ_NAMELIST: *** NOTE: Extra overhead invoked',  &
                      ' putting a/e numbers on restart dataset. ***'
            write(6,*)' To avoid, make mod(timesteps per day,iradae)= 0'
         end if
      end if
   end if
! 
! Build MSS pathname for restart file for branch run.
! Note that full (absolute) pathname must be input as nrevsn.
! 
   if (lbrnch .and. (nrevsn(1:1) /= '/') ) then
      call endrun ('READ_NAMELIST: NREVSN must be a full pathname for BRANCH run.')
   endif
!
! Restart files write frequency (on or off)
!
#if ( defined COUP_CSM )
   nrefrq = 1
#else
   if (nrefrq /= 0) then
      if ((nrefrq /= 1)) then
         call endrun ('READ_NAMELIST: the value of NREFRQ must be 1 or 0')
      endif
   end if
#endif

#if ( defined SPMD )
! 
! Set runtime options for spmd_dyn
!
   if ( dycore_is ('LR') ) then
      call spmd_dyn_setopts(                    &
             npr_yz_in         =npr_yz,         &
             geopktrans_in     =geopktrans,     &
             tracertrans_in    =tracertrans,    &
             ompnest_in        =ompnest,        &
             force_2d_in       =force_2d,       &
             modcomm_transpose_in =modcomm_transpose, &
             modcomm_geopk_in     =modcomm_geopk)
   endif
   if ( dycore_is ('EUL') .or. dycore_is ('SLD') ) then
      call spmd_dyn_setopts(                    &
             dyn_alltoall_in   =dyn_alltoall,   &
             dyn_allgather_in  =dyn_allgather   )
   endif
! 
! Set runtime options for swap communications.
!
   call swap_comm_setopts(                          &
          swap_comm_order_in=swap_comm_order,       &
          swap_comm_protocol_in=swap_comm_protocol)
#endif
! 
! Set runtime options for physics chunking.
!
   call phys_grid_setopts(                       &
          phys_loadbalance_in =phys_loadbalance, &
          phys_alltoall_in    =phys_alltoall,   &
          phys_chnk_per_thd_in=phys_chnk_per_thd)
! 
! exit if conflicts between prognostics and prescribed
!
   if(.not.( prescribed_sulfur == 'direct' .or. prognostic_sulfur == 'direct' )) then
     write(6,*)'either prescribed_sulfur or prognostic_sulfur must be direct'
     call endrun
   endif

   if(prescribed_sulfur == prognostic_sulfur ) then
     write(6,*)'prescribed_sulfur and prognostic_sulfur cannot be the same'
     call endrun
   endif

! 
! Set values of runtime options for prognostic aerosols
!
   call aerosol_setopts(                                  &
          prognostic_sulfur_in     =prognostic_sulfur,    &
          aero_carbon_in           =aero_carbon,          &
          aero_feedback_carbon_in  =aero_feedback_carbon, &
          aero_sea_salt_in         =aero_sea_salt,        &
          aero_feedback_sea_salt_in=aero_feedback_sea_salt)
! 
! Set runtime options for comhd
!
   call comhd_setopts( dif2_in  =dif2, &
                       dif4_in  =dif4, &
                       kmxhdc_in=kmxhdc)

!
! Initialize file paths module
!
   call init_filepaths( archivedirname=archive_dir )
!
! If branch set restart filepath to path given on namelist
!
   if ( lbrnch ) call set_restart_filepath( nrevsn )
! 
! Print camexp input variables to standard output
!
! TBH:  Need to prepend standard CCSM text...  
! 
   if (masterproc) then
      write(6,*)'READ_NAMELIST:rest_pfile= ',rest_pfile
      write(6,*)' ------------------------------------------'
      write(6,*)'     *** INPUT VARIABLES (CAMEXP) ***'
      write(6,*)' ------------------------------------------'
      if (nlres) then
         write(6,*) '  Continuation of an earlier run'
      else
         write(6,*) '         Initial run'
      end if
      write(6,*) ' ********** CASE = ',trim(caseid),' **********'
      write(6,'(1x,a)') ctitle
      if (len_trim(ncdata) > 0) then
         write(6,*) 'Initial dataset is: ',trim(ncdata)
      end if
      write(6,*) ' History-file archive directory = ', trim(get_archivedir('hist'))
      write(6,*) ' Restart-file archive directory = ', trim(get_archivedir('rest'))
      write(6,*) ' Initial-file archive directory = ', trim(get_archivedir('init'))
#if ( ! defined COUP_CSM )
      write(6,*)'Time-variant boundary dataset (sst) is: ', trim(bndtvs)
#endif
      write(6,*)'Time-variant boundary dataset (ozone) is: ', trim(bndtvo)
      write(6,*)'Time-invariant (absorption/emissivity) factor dataset is: ', trim(absems_data)

      write(6,*)'Time-variant boundary dataset (aerosols) is: ', trim(bndtvaer)
      write(6,*)'Time-variant boundary dataset (carbonscale) is: ', trim(bndtvcarbonscale)
      write(6,*)'Time-variant boundary dataset (solar constant) is: ', trim(bndtvscon)
      write(6,*)'Time-variant boundary dataset (greenhouse gas surface values) is: ', trim(bndtvghg)
      write(6,*)'Time-variant boundary dataset (volcanics) is: ', trim(bndtvvolc)
      write(6,*)'Aerosol Optics dataset is: ', trim(aeroptics)

      write(6,*)'Time-variant boundary dataset (carbon emissions) is: ', trim(co_emis)
      write(6,*)'Time-variant boundary dataset (DMS emissions) is: ', trim(bndtvdms)
      write(6,*)'Time-variant boundary dataset (soil erodibility) is: ', trim(soil_erod)
      write(6,*)'Time-variant boundary dataset (oxidants) is: ', trim(bndtvoxid)
      write(6,*)'Time-variant boundary dataset (SOx emissions) is: ', trim(bndtvsox)

!
! Restart files info
!
      if (nrefrq == 1) then
         write(6,*)'READ_NAMELIST3:rest_pfile=',rest_pfile
         write(6,*)'Restart pointer file is: ',trim(rest_pfile)
      else if (nrefrq==0) then 
         write(6,*) 'NO RESTART DATASET will be written'
      endif
#if ( defined COUP_CSM )
      write(6,*)'Restart files will be written only when specified by the flux coupler'
#endif
!
! Write password
!
      if (mss_wpass /='        ') then
         write(6,*)'Write passwd for output tapes (MSS_WPASS) is ', mss_wpass
      end if
!
! Type of run
!
      write(6,*)'Restart flag (NSREST) 0=no,1=yes,3=branch ',nsrest
   end if
!
! Print retention period for mass store
!
   if (mss_irt > 0) then
      if (mss_irt > 4096) then
         mss_irt = 4096
      end if
      if (masterproc) then
         write(6,*) 'Retention time for output files = ',mss_irt,' days'
      end if
   else
      if (masterproc) write(6,*) 'Output files will NOT be disposed to Mass Store'
   end if
!
! History file info 
!
   if (masterproc) then
      if (inithist == '6-HOURLY' ) then
         write(6,*)'Initial conditions history files will be written 6-hourly.'
      else if (inithist == 'DAILY' ) then
         write(6,*)'Initial conditions history files will be written daily.'
      else if (inithist == 'MONTHLY' ) then
         write(6,*)'Initial conditions history files will be written monthly.'
      else if (inithist == 'YEARLY' ) then
         write(6,*)'Initial conditions history files will be written yearly.'
      else
         write(6,*)'Initial conditions history files will not be created'
      end if
   end if
!
! Write physics variables from namelist camexp to std. output
!
   if (masterproc) then

#if ( defined COUP_CSM )
      write(6,*)'Ending time step determined by flux coupler'
#endif

      write(6,9108) eps,dif2,dif4,kmxhdc,nlvdry
      write(6,9110) iradsw,iradlw,iradae,itsst

9108 format(' Time filter coefficient (EPS)                 ',f10.3,/,&
            ' DEL2 Horizontal diffusion coefficient (DIF2)  ',e10.3/, &
            ' DEL4 Horizontal diffusion coefficient (DIF4)  ',e10.3/, &
            ' Number of levels Courant limiter applied      ',i10/,   &
            ' Lowest level for dry adiabatic adjust (NLVDRY)',i10)

9110 format(' Frequency of Shortwave Radiation calc. (IRADSW)     ',i5/, &
            ' Frequency of Longwave Radiation calc. (IRADLW)      ',i5/,  &
            ' Frequency of Absorptivity/Emissivity calc. (IRADAE) ',i5/, &
            ' Frequency of SST Initialization calc. (ITSST)       ',i5)

      if (sstcyc) then
         write(6,*)'SST dataset will be reused for each model year'
      else
         write(6,*)'SST dataset will not be cycled'
      end if

#ifndef COUP_CSM
      if ( icemodel_is('csim') .and. reset_csim_iceprops) then
         write(6,*)'CSIM ICE properties being reset to a new base state'
      end if
      if (prognostic_icesnow) then
         write(6,*)'Snow will accumulate to a maximum over sea-ice'
      else
         write(6,*)'Snow over sea-ice will be set to a climatology'
      end if
#endif

      if (icecyc) then
         write(6,*)'ICE dataset will be reused for each model year'
      else
         write(6,*)'ICE dataset will not be cycled'
      end if

      if (ozncyc) then
         write(6,*)'OZONE dataset will be reused for each model year'
      else
         write(6,*)'OZONE dataset will not be cycled'
      end if

      write (6,*)'Output files will be disposed ASYNCHRONOUSLY'

      if (divdampn > 0.) then
         write(6,*) 'Divergence damper invoked for days 0. to ',divdampn,' of this case'
      elseif (divdampn < 0.) then
         call endrun ('READ_NAMELIST: divdampn must be a positive number')
      else
         write(6,*) 'divergence damper NOT invoked'
      endif

      if ( (adiabatic .and. ideal_phys) .or. (adiabatic .and. aqua_planet) .or. &
           (ideal_phys .and. aqua_planet) ) then
         call endrun ('READ_NAMELIST: Only one of ADIABATIC, IDEAL_PHYS, or AQUA_PLANET can be .true.')
      end if

      if (adiabatic)   write(6,*) 'Model will run ADIABATICALLY (i.e. no physics)'
      if (ideal_phys)  write(6,*) 'Run ONLY the "idealized" dynamical core of the ', &
                                  'model  (dynamics + Held&Suarez-specified physics)'
      if (aqua_planet) write(6,*) 'Run model in "AQUA_PLANET" mode'
   end if

#ifdef PERGRO
   if (masterproc) then
      write(6,*)'pergro for cloud water is true'
   end if
#endif

   if (masterproc) then
      write(6,*) 'Visible optical depth (tauback) = ',tauback

#if ( defined COUP_CSM )
!
! Write coupled model input
!
      if (flxave) then
         write (6,*) 'Data will be sent to the flux coupler ', &
              'only on solar radiation time steps and ', &
              'the precipitation fluxes will be averaged ', &
              'on steps where communication with the flux ', &
              'coupler does not occur'
      else
         write (6,*) 'Data will be sent and received to/from ', &
              'the flux coupler at every time step except for ', &
              'nstep=1'
      endif

      if (        (iyear_AD /= SHR_ORB_UNDEF_INT )  &
           .or. (eccen    /= SHR_ORB_UNDEF_REAL)  &
           .or. (obliq    /= SHR_ORB_UNDEF_REAL)  &
           .or. (mvelp    /= SHR_ORB_UNDEF_REAL) )then
         write(6,*)' WARNING: Orbital parameters set from namelist'
         write(6,*)' will be overwritten by those obtained from coupler'
      end if
      write(6,*)' ------------------------------------------'
#else
      call shr_orb_print( iyear_AD, eccen, obliq, mvelp )
      write(6,*)' ------------------------------------------'
#endif
   end if

#if ( ! defined COUP_CSM )
#ifdef COUP_SOM
   if (som_conschk_frq < 0) then
      som_conschk_frq = -som_conschk_frq*ntspdy
   end if
   if (masterproc) then
      write(6,*)'SOM option is ENABLED'
      if (som_conschk_frq > 0) then
         write(6,*)'SOM global energy checking will be done every ',som_conschk_frq,' timesteps'
      end if
   end if
#endif

   if (ice_conschk_frq < 0) then
      ice_conschk_frq = -ice_conschk_frq*ntspdy
   end if
   if (masterproc .and. ice_conschk_frq > 0) then
      write(6,*)'ICE global energy checking will be done every ',ice_conschk_frq,' timesteps'
   end if
#endif

   if (masterproc) then
      if (doisccp) then
         write(6,*)'ISCCP calcs and history IO will be done'
      else
         write(6,*)'ISCCP calcs and history IO will NOT be done'
      end if
   end if

   return
end subroutine read_namelist


!=======================================================================

#ifdef SPMD
subroutine distnl
!-----------------------------------------------------------------------
!     
! Purpose:     
! Distribute namelist data to all processors.
!
! The cpp SPMD definition provides for the funnelling of all program i/o
! through the master processor. Processor 0 either reads restart/history
! data from the disk and distributes it to all processors, or collects
! data from all processors and writes it to disk.
!     
!---------------------------Code history-------------------------------
!
! Original version:  CCM2
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   use mpishorthand
   use comhd, only: dif2, dif4, kmxhdc
!-----------------------------------------------------------------------

#include <comadj.h>
#include <comctl.h>
#include <comsol.h>
#include <comtfc.h>
!
!-----------------------------------------------------------------------
! 
   call mpibcast (calendar,   32,mpichar,0,mpicom)
   call mpibcast (dtime,       1,mpiint,0,mpicom)
   call mpibcast (nestep,      1,mpiint,0,mpicom)
   call mpibcast (nelapse,     1,mpiint,0,mpicom)
   call mpibcast (start_ymd,   1,mpiint,0,mpicom)
   call mpibcast (start_tod,   1,mpiint,0,mpicom)
   call mpibcast (stop_ymd,    1,mpiint,0,mpicom)
   call mpibcast (stop_tod,    1,mpiint,0,mpicom)
   call mpibcast (ref_ymd,     1,mpiint,0,mpicom)
   call mpibcast (ref_tod,     1,mpiint,0,mpicom)
   call mpibcast (perpetual_run, 1,mpilog,0,mpicom)
   call mpibcast (perpetual_ymd, 1,mpiint,0,mpicom)

   call mpibcast (nhstpr  ,ptapes,mpiint,0,mpicom)
   call mpibcast (ndens   ,ptapes,mpiint,0,mpicom)
   call mpibcast (nhtfrq  ,ptapes,mpiint,0,mpicom)
   call mpibcast (mfilt   ,ptapes,mpiint,0,mpicom)
   call mpibcast (nsrest  ,1,mpiint,0,mpicom)
   call mpibcast (mss_irt ,1,mpiint,0,mpicom)
   call mpibcast (nrefrq  ,1,mpiint,0,mpicom)
   call mpibcast (kmxhdc  ,1,mpiint,0,mpicom)
   call mpibcast (iradsw  ,1,mpiint,0,mpicom)
   call mpibcast (iradlw  ,1,mpiint,0,mpicom)
   call mpibcast (iradae  ,1,mpiint,0,mpicom)
   call mpibcast (itsst   ,1,mpiint,0,mpicom)
   call mpibcast (nlvdry  ,1,mpiint,0,mpicom)
#if ( ! defined COUP_CSM )
   call mpibcast (reset_csim_iceprops,1,mpilog,0,mpicom)
   call mpibcast (prognostic_icesnow,1,mpilog,0,mpicom)
#endif
   call mpibcast (som_conschk_frq,1,mpiint,0,mpicom)
   call mpibcast (ice_conschk_frq,1,mpiint,0,mpicom)
! f-v dynamics specific
   call mpibcast (nsplit  ,1,mpiint,0,mpicom)
   call mpibcast (iord    ,1,mpiint,0,mpicom)
   call mpibcast (jord    ,1,mpiint,0,mpicom)
   call mpibcast (kord    ,1,mpiint,0,mpicom)
   call mpibcast (use_eta ,1,mpilog,0,mpicom)

   call mpibcast (divdampn,1,mpir8,0,mpicom)
   call mpibcast (co2vmr  ,1,mpir8,0,mpicom)
   call mpibcast (ch4vmr  ,1,mpir8,0,mpicom)
   call mpibcast (n2ovmr  ,1,mpir8,0,mpicom)
   call mpibcast (f11vmr  ,1,mpir8,0,mpicom)
   call mpibcast (f12vmr  ,1,mpir8,0,mpicom)
   call mpibcast (eps     ,1,mpir8,0,mpicom)
   call mpibcast (dif2    ,1,mpir8,0,mpicom)
   call mpibcast (dif4    ,1,mpir8,0,mpicom)

   call mpibcast (precc_thresh,1,mpir8,0,mpicom)
   call mpibcast (precl_thresh,1,mpir8,0,mpicom)

   call mpibcast (flxave      ,1,mpilog,0,mpicom)
   call mpibcast (adiabatic   ,1,mpilog,0,mpicom)
   call mpibcast (trace_gas   ,1,mpilog,0,mpicom)
   call mpibcast (trace_water ,1,mpilog,0,mpicom)
   call mpibcast (wisotope    ,1,mpilog,0,mpicom)
   call mpibcast (tracers_flag  ,1,mpilog,0,mpicom)
   call mpibcast (readtrace   ,1,mpilog,0,mpicom)
   call mpibcast (sstcyc      ,1,mpilog,0,mpicom)
   call mpibcast (icecyc      ,1,mpilog,0,mpicom)
   call mpibcast (ozncyc      ,1,mpilog,0,mpicom)
   call mpibcast (ideal_phys  ,1,mpilog,0,mpicom)
   call mpibcast (aqua_planet ,1,mpilog,0,mpicom)
   call mpibcast (empty_htapes,1,mpilog,0,mpicom)
   call mpibcast (print_step_cost,1,mpilog,0,mpicom)
   call mpibcast (doisccp     ,1,mpilog,0,mpicom)

   call mpibcast (caseid  ,len(caseid) ,mpichar,0,mpicom)
   call mpibcast (avgflag_pertape, ptapes, mpichar,0,mpicom)
   call mpibcast (ctitle  ,len(ctitle),mpichar,0,mpicom)
   call mpibcast (ncdata  ,len(ncdata) ,mpichar,0,mpicom)
   call mpibcast (bndtvs  ,len(bndtvs) ,mpichar,0,mpicom)
   call mpibcast (bndtvo  ,len(bndtvo) ,mpichar,0,mpicom)
   call mpibcast (bndtvaer  ,len(bndtvaer) ,mpichar,0,mpicom)
   call mpibcast (bndtvcarbonscale  ,len(bndtvcarbonscale) ,mpichar,0,mpicom)
   call mpibcast (bndtvscon  ,len(bndtvscon) ,mpichar,0,mpicom)
   call mpibcast (bndtvghg  ,len(bndtvghg) ,mpichar,0,mpicom)
   call mpibcast (bndtvvolc ,len(bndtvvolc) ,mpichar,0,mpicom)
   call mpibcast (co_emis  ,len(co_emis) ,mpichar,0,mpicom)
   call mpibcast (bndtvdms  ,len(bndtvdms) ,mpichar,0,mpicom)
   call mpibcast (soil_erod  ,len(soil_erod) ,mpichar,0,mpicom)
   call mpibcast (bndtvoxid  ,len(bndtvoxid) ,mpichar,0,mpicom)
   call mpibcast (bndtvsox  ,len(bndtvsox) ,mpichar,0,mpicom)
   call mpibcast (aeroptics  ,len(aeroptics) ,mpichar,0,mpicom)
   call mpibcast (absems_data,len(absems_data),mpichar,0,mpicom)
   call mpibcast (bndtvg  ,len(bndtvg),mpichar,0,mpicom)
   call mpibcast (bndtvsf6  ,len(bndtvsf6),mpichar,0,mpicom)
   call mpibcast (mss_wpass,len(mss_wpass)  ,mpichar,0,mpicom)
   call mpibcast (nrevsn  ,len(nrevsn) ,mpichar,0,mpicom)
   call mpibcast (inithist,len(inithist)  ,mpichar,0,mpicom)
   call mpibcast (fincl   ,10*pflds*ptapes,mpichar,0,mpicom)
   call mpibcast (fexcl   , 8*pflds*ptapes,mpichar,0,mpicom)
   call mpibcast (fhstpr  ,10*pflds*ptapes,mpichar,0,mpicom)
   call mpibcast (fwrtpr  ,10*pflds*ptapes,mpichar,0,mpicom)
!
! Orbital stuff
!
   call mpibcast (scon    ,1  ,mpir8 ,0,mpicom)
   call mpibcast (eccen   ,1  ,mpir8 ,0,mpicom)
   call mpibcast (obliq   ,1  ,mpir8 ,0,mpicom)
   call mpibcast (mvelp   ,1  ,mpir8 ,0,mpicom)
   call mpibcast (iyear_ad,1  ,mpiint,0,mpicom)
!
   call mpibcast (scenario_ghg ,16 ,mpichar,0,mpicom)
   call mpibcast (scenario_prognostic_sulfur ,16 ,mpichar,0,mpicom)
   call mpibcast (scenario_scon,16 ,mpichar,0,mpicom)
   call mpibcast (rampYear_ghg , 1 ,mpiint, 0,mpicom)
   call mpibcast (rampyear_prognostic_sulfur , 1 ,mpiint, 0,mpicom)
   call mpibcast (rampYear_scon, 1 ,mpiint, 0,mpicom)
   call mpibcast (indirect     , 1 ,mpilog, 0,mpicom)

   call mpibcast (ramp_co2_start_ymd     , 1 ,mpiint, 0,mpicom)
   call mpibcast (ramp_co2_annual_rate     , 1 ,mpir8, 0,mpicom)
   call mpibcast (ramp_co2_cap     , 1 ,mpir8, 0,mpicom)

!
!  Aerosol stuff
!
   call mpibcast (radforce,        1, mpilog, 0,mpicom)
   call mpibcast (scenario_carbon_scale, 16, mpichar, 0,mpicom)
   call mpibcast (scenario_prescribed_sulfur ,16 ,mpichar,0,mpicom)
   call mpibcast (rampyear_prescribed_sulfur ,1 ,mpiint,0,mpicom)
   call mpibcast (prescribed_sulfur ,16 ,mpichar,0,mpicom)

   call mpibcast (bgscl_rf,   1, mpir8, 0,mpicom)
   call mpibcast (tauback,    1, mpir8,0,mpicom)

   call mpibcast (sulscl_rf,  1, mpir8, 0,mpicom)
   call mpibcast (carscl_rf,  1, mpir8, 0,mpicom)
   call mpibcast (ssltscl_rf, 1, mpir8, 0,mpicom)
   call mpibcast (dustscl_rf, 1, mpir8, 0,mpicom)

   call mpibcast (sulscl,     1, mpir8 ,0,mpicom)
   call mpibcast (carscl,     1, mpir8 ,0,mpicom)
   call mpibcast (ssltscl,    1, mpir8 ,0,mpicom)
   call mpibcast (dustscl,    1, mpir8 ,0,mpicom)

   call mpibcast (strat_volcanic,  1, mpilog, 0,mpicom)
   call mpibcast (volcscl_rf, 1, mpir8, 0,mpicom)
   call mpibcast (volcscl,    1, mpir8 ,0,mpicom)

!
!  spmd_dyn stuff
!
   if ( dycore_is ('LR') ) then
      call mpibcast (npr_yz        ,4,mpiint,0,mpicom)
      call mpibcast (geopktrans    ,1,mpiint,0,mpicom)
      call mpibcast (tracertrans   ,1,mpiint,0,mpicom)
      call mpibcast (ompnest       ,1,mpiint,0,mpicom)
      call mpibcast (force_2d      ,1,mpiint,0,mpicom)
      call mpibcast (modcomm_transpose,1,mpiint,0,mpicom)
      call mpibcast (modcomm_geopk    ,1,mpiint,0,mpicom)
   endif
   if ( dycore_is ('EUL') .or. dycore_is ('SLD') ) then
      call mpibcast (dyn_alltoall  ,1,mpiint,0,mpicom)
      call mpibcast (dyn_allgather ,1,mpiint,0,mpicom)
   endif

!
!  Physics chunk tuning stuff
!
   call mpibcast (phys_loadbalance ,1,mpiint,0,mpicom)
   call mpibcast (phys_alltoall    ,1,mpiint,0,mpicom)
   call mpibcast (phys_chnk_per_thd,1,mpiint,0,mpicom)

!
!  Interprocessor communication tuning stuff
!
   call mpibcast (swap_comm_order ,1,mpiint,0,mpicom)
   call mpibcast (swap_comm_protocol,1,mpiint,0,mpicom)

!
!  Prognostic aerosol stuff
!
   call mpibcast (prognostic_sulfur,     16 ,mpichar, 0,mpicom)
   call mpibcast (aero_carbon,           1 ,mpilog, 0,mpicom)
   call mpibcast (aero_feedback_carbon,  1 ,mpilog, 0,mpicom)
   call mpibcast (aero_sea_salt,         1 ,mpilog, 0,mpicom)
   call mpibcast (aero_feedback_sea_salt,1 ,mpilog, 0,mpicom)

   return
end subroutine distnl
#endif



subroutine preset
!----------------------------------------------------------------------- 
! 
! Purpose: Preset namelist CAMEXP input variables and initialize some other variables
! 
! Method: Hardwire the values
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
   use history,      only: fincl, fexcl, fhstpr, fwrtpr,fincllonlat
   use rgrid
#if ( ! defined COUP_CSM )
   use ice_dh, only: prognostic_icesnow,reset_csim_iceprops
#endif
   use time_manager, only: timemgr_preset
!------------------------------Commons----------------------------------
#include <comadj.h>
!-----------------------------------------------------------------------
#include <comctl.h>
!-----------------------------------------------------------------------
#include <comlun.h>
!-----------------------------------------------------------------------
#include <comtfc.h>
!-----------------------------------------------------------------------
#include <comsol.h>
!-----------------------------------------------------------------------
#include <perturb.h>
!-----------------------------------------------------------------------
   include 'netcdf.inc'
!-----------------------------------------------------------------------
!
! Preset character history variables here because module initialization of character arrays
! does not work on all machines
! $$$ TBH:  is this still true?  12/14/03
!
   fincl(:,:)  = ' '
   fincllonlat(:,:)  = ' '
   fexcl(:,:)  = ' '
   fhstpr(:,:) = ' '
   fwrtpr(:,:) = ' '
!
! Flags
!
   nlend       = .false.       ! end of run flag
   nlres       = .false.       ! continuation run flag
   nlhst       = .false.       ! regen run or branch run flag
   lbrnch      = .false.       ! branch run flag
   adiabatic   = .false.       ! no physics
   ideal_phys  = .false.       ! "idealized" model configuration
   aqua_planet = .false.       ! global oceans/analytical SST's
   print_step_cost = .false.   ! print per timestep cost info
!
! Ice flags
!
#if ( ! defined COUP_CSM )
   prognostic_icesnow = .true.    ! snow falls on ice by default but
                                  ! it is limited to 0.5 meter.
   reset_csim_iceprops= .false.   ! use initial condition info unless
                                  ! need to reset ice properties in csim
   ice_conschk_frq = 0
   som_conschk_frq = 0
#endif
!
! Default run type is initialization
!
   nsrest = 0
!
! Default value for writing restart files
!
   nrefrq = 1      ! normal run, dispose with full history file
!
! Cycling flags for input boundary data files
!
   sstcyc = .true.
   icecyc = .true.
   ozncyc = .true.
!
! Model time defaults
!
   call timemgr_preset()
!
! Frequency in iterations of absorptivity/emissivity calc (negative
! values in model hours)
!
   iradae = -12
!
! Frequency of annual cycle sst update
!
   itsst  =  1
!
! Default frequency of shortwave and longwave radiation computations: 
! once per hour (negative value implies model hours)
!
   iradsw = -1
   iradlw = -1
!
! Numerical scheme default values
!
   eps    = 0.06
   nlvdry = 3
!
! No divergence damping
!
   divdampn = 0.
!
! Precipitation threshold for PRECCINT, PRECLINT, PRECCFRQ, and PRECLFRQ output fields
! (mm/hr)
!
   precc_thresh = 0.1
   precl_thresh = 0.05
!
! Orbital parameters.
! NOTE: if iyear_AD is set to SHR_ORB_UNDEF_INT after namelist input
! then namelist values of obliq,eccen,and mvelp are used otherwise
! obliq,eccen and mvelp are calculated based on iyear_AD
!
   iyear_ad = shr_orb_undef_int  
   obliq    = shr_orb_undef_real
   eccen    = shr_orb_undef_real
   mvelp    = shr_orb_undef_real
!
! Solar constant
!
   scon       = 1.367e6

#if ( defined COUP_CSM )
!
! Communications with the flux coupler
!
   flxave = .true.
#endif
!
! rgrid: set default to full grid
!
   nlon(:) = plon
!
! Unit numbers: set to invalid
!
   nsds     = -1
   nrg      = -1
   nrg2     = -1
   ncid_ini = -1
   ncid_oz  = -1
   ncid_sst = -1
   ncid_trc = -1
   luhrest  = -1
!
! /perturb/
!
  pertlim = 0.0

   return
end subroutine preset


!=======================================================================
  subroutine runtime_options( )
!----------------------------------------------------------------------- 
!
! Purpose:  Set default values of runtime options 
!           before namelist camexp is read, then
!           read namelist (and broadcast, if SPMD).  
!
! Method:   Calls preset() and read_namelist (which 
!           used to be called parse_namelist()).  
!
! Author:  Tom Henderson
!
!-----------------------------------------------------------------------

    !
    ! Set defaults then override with user-specified input
    !
    call preset ()
#if (!defined SCAM) 
    call read_namelist ()    ! used to be called parse_namelist()
#endif 
  end subroutine runtime_options


end module runtime_opts
