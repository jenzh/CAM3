#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: driver
!
! !INTERFACE:
subroutine driver (doalb, eccen, obliqr, lambm0, mvelpp)
!
! !DESCRIPTION:
! This subroutine provides the main CLM driver calling sequence.  Most
! computations occurs over ``clumps'' of gridcells (and associated subgrid
! scale entities) assigned to each MPI process.  Computation is further
! parallelized by looping over clumps on each process using shared memory
! OpenMP or Cray Streaming Directives.
!
! The main CLM driver calling sequence is as follows:
! \begin{verbatim}
! * Communicate with flux coupler [COUP_CSM]
! + interpMonthlyVeg      interpolate monthly vegetation data [!DGVM]
!   + readMonthlyVegetation read vegetation data for two months [!DGVM]
! ==== Begin Loop 1 over clumps ====
!  -> DriverInit          save of variables from previous time step
!  -> Hydrology1          canopy interception and precip on ground
!     -> FracWet          fraction of wet vegetated surface and dry elai
!  -> SurfaceRadiation    surface solar radiation
!  -> Biogeophysics1      leaf temperature and surface fluxes
!  -> BareGroundFluxes    surface fluxes for bare soil or snow-covered
!                         vegetation patches
!     -> MoninObukIni     first-guess Monin-Obukhov length and wind speed
!     -> FrictionVelocity friction velocity and potential temperature and
!                         humidity profiles
!  -> CanopyFluxes        leaf temperature and surface fluxes for vegetated
!                         patches
!     -> QSat             saturated vapor pressure, specific humidity, &
!                         derivatives at leaf surface
!     -> MoninObukIni     first-guess Monin-Obukhov length and wind speed
!     -> FrictionVelocity friction velocity and potential temperature and
!                         humidity profiles
!     -> Stomata          stomatal resistance and photosynthesis for
!                         sunlit leaves
!     -> Stomata          stomatal resistance and photosynthesis for
!                         shaded leaves
!     -> QSat             recalculation of saturated vapor pressure,
!                         specific humidity, & derivatives at leaf surface
!  -> Biogeophysics_Lake  lake temperature and surface fluxes
!   + VOCEmission         compute VOC emission [VOC]
!   + DGVMRespiration     CO2 respriation and plant production [DGVM]
!   + DGVMEcosystemDyn    DGVM ecosystem dynamics: vegetation phenology [!DGVM]
!  -> EcosystemDyn        "static" ecosystem dynamics: vegetation phenology
!                         and soil carbon [!DGVM]
!  -> SurfaceAlbedo       albedos for next time step
!  -> Biogeophysics2      soil/snow & ground temp and update surface fluxes
!  -> pft2col             Average from PFT level to column level
!  ====  End Loop 1 over clumps  ====
! * Average fluxes over time interval and send to flux coupler [COUP_CSM]
!  ==== Begin Loop 2 over clumps ====
!  -> Hydrology2          surface and soil hydrology
!  -> Hydrology_Lake      lake hydrology
!  -> SnowAge             update snow age for surface albedo calcualtion
!  -> BalanceCheck        check for errors in energy and water balances
!  ====  End Loop 2 over clumps  ====
!  -> write_diagnostic    output diagnostic if appropriate
!   + Rtmriverflux        calls RTM river routing model [RTM]
!  -> updateAccFlds       update accumulated fields
!  -> update_hbuf         accumulate history fields for time interval
!  Begin DGVM calculations at end of model year [DGVM]
!    ==== Begin Loop over clumps ====
!     + lpj                 LPJ ecosystem dynamics: reproduction, turnover,
!                           kill, allocation, light, mortality, fire
!     + lpjreset1           reset variables & initialize for next year
!    ====  End Loop over clumps  ====
!  End DGVM calculations at end of model year [DGVM]
!  -> htapes_wrapup       write history tapes if appropriate
!  Begin DGVM calculations at end of model year [DGVM]
!    ==== Begin Loop over clumps ====
!     + lpjreset2           reset variables and patch weights
!    ====  End Loop over clumps  ====
!  End DGVM calculations at end of model year [DGVM]
!  -> restart             write restart file if appropriate
!  -> inicfile            write initial file if appropriate
! \end{verbatim}
! Optional subroutines are denoted by an plus (+) with the associated
! CPP variable in brackets at the end of the line.  Coupler communication
! when coupled with CCSM components is denoted by an asterisk (*).
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
#if (defined COUP_CSM)
  use clm_varctl          , only : wrtdia, fsurdat, csm_doflxave
#else
  use clm_varctl          , only : wrtdia, fsurdat, nsrest
#endif
  use spmdMod             , only : masterproc
  use decompMod           , only : get_proc_clumps, get_clump_bounds
  use filterMod           , only : filter
  use clm_varcon          , only : zlnd
#if (defined COUP_CAM)
  use time_manager        , only : get_step_size, get_curr_calday, &
                                   get_curr_date, get_ref_date, get_nstep, &
                                   is_perpetual
#else
  use time_manager        , only : get_step_size, get_curr_calday, &
                                   get_curr_date, get_ref_date, get_nstep
#endif
  use histFileMod         , only : update_hbuf, htapes_wrapup
  use restFileMod         , only : restart
#if (defined COUP_CAM)
  use inicFileMod         , only : inicfile, do_inicwrite, inicperp
#else
  use inicFileMod         , only : inicfile, do_inicwrite
#endif
  use DriverInitMod       , only : DriverInit
  use BalanceCheckMod     , only : BalanceCheck
  use SurfaceRadiationMod , only : SurfaceRadiation
  use Hydrology1Mod       , only : Hydrology1
  use Hydrology2Mod       , only : Hydrology2
  use HydrologyLakeMod    , only : HydrologyLake
  use Biogeophysics1Mod   , only : Biogeophysics1
  use BareGroundFluxesMod , only : BareGroundFluxes
  use CanopyFluxesMod     , only : CanopyFluxes
  use Biogeophysics2Mod   , only : Biogeophysics2
  use BiogeophysicsLakeMod, only : BiogeophysicsLake
  use SurfaceAlbedoMod    , only : SurfaceAlbedo, Snowage
  use pft2colMod          , only : pft2col
  use accFldsMod          , only : updateAccFlds
#if (defined DGVM)
  use DGVMEcosystemDynMod , only : DGVMEcosystemDyn, DGVMRespiration
  use DGVMMod             , only : lpj, lpjreset1, lpjreset2, &
                                   gatherWeightsDGVM, histDGVM
#else
  use STATICEcosysDynMod  , only : EcosystemDyn, interpMonthlyVeg
#endif
#if (defined VOC)
  use VOCEmissionMod      , only : VOCEmission
#endif
#if (defined RTM)
  use RtmMod              , only : Rtmriverflux
#endif
#if (defined COUP_CSM)
  use clm_csmMod          , only : csm_dosndrcv, csm_recv, csm_send, &
                                   csm_flxave, dorecv, dosend, csmstop_now
#endif
  use lnd2atmMod          , only : lnd2atm
  use abortutils          , only : endrun
!
! !ARGUMENTS:
  implicit none
  logical , intent(in) :: doalb  !true if time for surface albedo
                                 !calculation
  real(r8), intent(in) :: eccen  !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0 !Mean longitude of perihelion at the
                                 !vernal equinox (radians)
  real(r8), intent(in) :: mvelpp !Earth's moving vernal equinox longitude
                                 !of perihelion + pi (radians)
!
! !CALLED FROM:
! program program_off (if COUP_OFFLINE cpp variable is defined)
! program program_csm (if COUP_CSM cpp variable is defined)
! subroutine atm_lnddrv in module atm_lndMod (if COUP_CAM cpp variable
!   is defined)
!
! !REVISION HISTORY:
! 2002.10.01  Mariana Vertenstein latest update to new data structures
!
!EOP
!
! !LOCAL VARIABLES:
  integer  :: nc, c         ! indices
  integer  :: nclumps       ! number of clumps on this processor
  integer  :: yrp1          ! year (0, ...) for nstep+1
  integer  :: monp1         ! month (1, ..., 12) for nstep+1
  integer  :: dayp1         ! day of month (1, ..., 31) for nstep+1
  integer  :: secp1         ! seconds into current date for nstep+1
  integer  :: yr            ! year (0, ...)
  integer  :: mon           ! month (1, ..., 12)
  integer  :: day           ! day of month (1, ..., 31)
  integer  :: sec           ! seconds of the day
  integer  :: ncdate        ! current date
  integer  :: nbdate        ! base date (reference date)
  integer  :: kyr           ! thousand years, equals 2 at end of first year
  real(r8) :: caldayp1      ! calendar day for nstep+1
  real(r8) :: dtime         ! land model time step (sec)
  integer  :: nstep         ! time step number
  integer  :: begp, endp    ! clump beginning and ending pft indices
  integer  :: begc, endc    ! clump beginning and ending column indices
  integer  :: begl, endl    ! clump beginning and ending landunit indices
  integer  :: begg, endg    ! clump beginning and ending gridcell indices
  type(column_type)  , pointer :: cptr    ! pointer to column derived subtype
  logical, external :: do_restwrite ! determine if time to write restart
!-----------------------------------------------------------------------

  ! Set pointers into derived type

  cptr => clm3%g%l%c

  call t_startf('clm_driver')

#if (defined COUP_CSM)
  ! ============================================================================
  ! Coupler receive
  ! ============================================================================

  ! Determine if information should be sent/received to/from flux coupler
  call csm_dosndrcv(doalb)

  ! Get atmospheric state and fluxes from flux coupler
  if (dorecv) then
     call csm_recv()
     if (csmstop_now) then
        call t_stopf('clm_driver')
        RETURN
     endif
  endif
#endif

  ! ============================================================================
  ! Calendar information for next time step
  ! o caldayp1 = calendar day (1.00 -> 365.99) for cosine solar zenith angle
  !   calday is based on Greenwich time
  ! o get_curr_calday in the cam time manager know about perpetual mode
  !   and perpetual model is only used within cam
  ! ============================================================================

  dtime = get_step_size()
  caldayp1 = get_curr_calday(offset=int(dtime))

#if (!defined DGVM)
  ! ============================================================================
  ! Determine weights for time interpolation of monthly vegetation data.
  ! This also determines whether it is time to read new monthly vegetation and
  ! obtain updated leaf area index [mlai1,mlai2], stem area index [msai1,msai2],
  ! vegetation top [mhvt1,mhvt2] and vegetation bottom [mhvb1,mhvb2]. The
  ! weights obtained here are used in subroutine ecosystemdyn to obtain time
  ! interpolated values.
  ! ============================================================================

  if (doalb) call interpMonthlyVeg()
#endif

  ! ============================================================================
  ! Loop1
  ! ============================================================================

  nclumps = get_proc_clumps()

  call t_startf('loop1')

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
  do nc = 1,nclumps

     ! ============================================================================
     ! Determine clump boundaries
     ! ============================================================================

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! ============================================================================
     ! Initialize variables from previous time step and
     ! Determine canopy interception and precipitation onto ground surface.
     ! Determine the fraction of foliage covered by water and the fraction
     ! of foliage that is dry and transpiring. Initialize snow layer if the
     ! snow accumulation exceeds 10 mm.
     ! ============================================================================

     call t_startf('drvinit')
     call DriverInit(begc, endc, begp, endp, &
          filter(nc)%num_nolakec, filter(nc)%nolakec, filter(nc)%num_lakec, filter(nc)%lakec)
     call t_stopf('drvinit')

     ! ============================================================================
     ! Hydrology1
     ! ============================================================================

     call t_startf('hydro1')
     call Hydrology1(begc, endc, begp, endp, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('hydro1')

     ! ============================================================================
     ! Surface Radiation
     ! ============================================================================

     call t_startf('surfrad')
     call SurfaceRadiation(begp, endp)
     call t_stopf('surfrad')

     ! ============================================================================
     ! Determine leaf temperature and surface fluxes based on ground
     ! temperature from previous time step.
     ! ============================================================================

     call t_startf('bgp1')
     call Biogeophysics1(begg, endg, begc, endc, begp, endp, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp1')

     ! ============================================================================
     ! Determine bare soil or snow-covered vegetation surface temperature and fluxes
     ! Calculate Ground fluxes (frac_veg_nosno is either 1 or 0)
     ! ============================================================================

     call t_startf('bgflux')
     call BareGroundFluxes(begp, endp, &
                           filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgflux')

     ! ============================================================================
     ! Determine non snow-covered vegetation surface temperature and fluxes
     ! Calculate canopy temperature, latent and sensible fluxes from the canopy,
     ! and leaf water change by evapotranspiration
     ! ============================================================================

     call t_startf('canflux')
     call CanopyFluxes(begg, endg, begc, endc, begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('canflux')

     ! ============================================================================
     ! Determine lake temperature and surface fluxes
     ! ============================================================================

     call t_startf('bgplake')
     call BiogeophysicsLake(begc, endc, begp, endp, &
                            filter(nc)%num_lakec, filter(nc)%lakec, &
                            filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('bgplake')

     ! ============================================================================
     ! Determine VOC and DGVM Respiration if appropriate
     ! ============================================================================

     call t_startf('bgc')

#if (defined VOC)
     ! VOC emission (A. Guenther's model)
     call VOCEmission(begp, endp, &
                      filter(nc)%num_nolakep, filter(nc)%nolakep)
#endif

     call t_stopf('bgc')

     ! ============================================================================
     ! Ecosystem dynamics: phenology, vegetation, soil carbon, snow fraction
     ! ============================================================================

     call t_startf('ecosysdyn')
#if (defined DGVM)
     ! Surface biogeochemical fluxes: co2 respiration and plant production
     call DGVMRespiration(begc, endc, begp, endp, &
                          filter(nc)%num_nolakec, filter(nc)%nolakec, &
                          filter(nc)%num_nolakep, filter(nc)%nolakep)

     call DGVMEcosystemDyn(begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep, &
                       doalb, endofyr=.false.)
#else
     call EcosystemDyn(begp, endp, &
                       filter(nc)%num_nolakep, filter(nc)%nolakep, &
                       doalb)
#endif
     call t_stopf('ecosysdyn')

     ! ============================================================================
     ! Determine albedos for next time step
     ! ============================================================================

     if (doalb) then
        call t_startf('surfalb')
        call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, &
                           caldayp1, eccen, obliqr, lambm0, mvelpp)
        call t_stopf('surfalb')
     end if

     ! ============================================================================
     ! Determine soil/snow temperatures including ground temperature and
     ! update surface fluxes for new ground temperature.
     ! ============================================================================

     call t_startf('bgp2')
     call Biogeophysics2(begc, endc, begp, endp, &
                         filter(nc)%num_nolakec, filter(nc)%nolakec, &
                         filter(nc)%num_nolakep, filter(nc)%nolakep)
     call t_stopf('bgp2')

     ! ============================================================================
     ! Perform averaging from PFT level to column level
     ! ============================================================================

     call t_startf('pft2col')
     call pft2col(begc, endc, filter(nc)%num_nolakec, filter(nc)%nolakec)
     call t_stopf('pft2col')

  end do
!$OMP END PARALLEL DO
!CSD$ END PARALLEL DO

  call t_stopf('loop1')

#if (defined COUP_CSM)
  ! ============================================================================
  ! Coupler send
  ! ============================================================================

  ! Average fluxes over interval if appropriate
  ! Surface states sent to the flux coupler states are not time averaged
  if (csm_doflxave) call csm_flxave()

  ! Send fields to flux coupler
  ! Send states[n] (except for snow[n-1]), time averaged fluxes for [n,n-1,n-2],
  ! albedos[n+1], and ocnrof_vec[n]
  if (dosend) call csm_send()
#endif

  ! ============================================================================
  ! Loop2
  ! ============================================================================

  call t_startf('loop2')

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
  do nc = 1,nclumps

     ! ============================================================================
     ! Determine clump boundaries
     ! ============================================================================

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! ============================================================================
     ! Vertical (column) soil and surface hydrology
     ! ============================================================================

     call t_startf('hydro2')
     call Hydrology2(begc, endc, &
                     filter(nc)%num_nolakec, filter(nc)%nolakec, &
                     filter(nc)%num_soilc, filter(nc)%soilc, &
                     filter(nc)%num_snowc, filter(nc)%snowc, &
                     filter(nc)%num_nosnowc, filter(nc)%nosnowc)
     call t_stopf('hydro2')

     ! ============================================================================
     ! Lake hydrology
     ! ============================================================================

     call t_startf('hylake')
     call HydrologyLake(begp, endp, &
                        filter(nc)%num_lakep, filter(nc)%lakep)
     call t_stopf('hylake')

     ! ============================================================================
     ! Update Snow Age (needed for surface albedo calculation
     ! ============================================================================

     call t_startf('snowage')
     call SnowAge(begc, endc)
     call t_stopf('snowage')

     ! ============================================================================
     ! ! Fraction of soil covered by snow (Z.-L. Yang U. Texas)
     ! ============================================================================

!dir$ concurrent
!cdir nodep
     do c = begc,endc
        cptr%cps%frac_sno(c) = cptr%cps%snowdp(c) / (10.*zlnd + cptr%cps%snowdp(c))
     end do

     ! ============================================================================
     ! Check the energy and water balance
     ! ============================================================================

     call t_startf('balchk')
     call BalanceCheck(begp, endp, begc, endc)
     call t_stopf('balchk')

  end do
!$OMP END PARALLEL DO
!CSD$ END PARALLEL DO

  call t_stopf('loop2')

#if (defined OFFLINE)
  ! ============================================================================
  ! Determine fields that would be sent to atm for diagnostic purposes
  ! When not in offline mode, this is called from clm_csmMod and lp_coupling
  ! ============================================================================

  call lnd2atm()
#endif

  ! ============================================================================
  ! Write global average diagnostics to standard output
  ! ============================================================================

  nstep = get_nstep()
  call write_diagnostic(wrtdia, nstep)

#if (defined RTM)
  ! ============================================================================
  ! Route surface and subsurface runoff into rivers
  ! ============================================================================

  call t_startf('clmrtm')
  call Rtmriverflux()
  call t_stopf('clmrtm')
#endif

#if (defined COUP_CAM)
  ! ============================================================================
  ! Read initial snow and soil moisture data at each time step
  ! ============================================================================

  call t_startf('inicperp')
  if (is_perpetual()) call inicperp()
  call t_stopf('inicperp')
#endif

  ! ============================================================================
  ! Update accumulators
  ! ============================================================================

  call t_startf('accum')
  call updateAccFlds()
  call t_stopf('accum')

  ! ============================================================================
  ! Update history buffer
  ! ============================================================================

  call t_startf('hbuf')
  call update_hbuf()
  call t_stopf('hbuf')

#if (defined DGVM)
  ! ============================================================================
  ! Call DGVM (Dynamic Global Vegetation Model)
  ! LPJ is called at last time step of year. Then reset vegetation distribution
  ! and some counters for the next year.
  ! NOTE: monp1, dayp1, and secp1 correspond to nstep+1
  ! NOTE: lpjreset1 must be called after update_accum and update_hbuf
  ! in order to have the correct values of the accumulated variables
  ! NOTE: lpjreset2 is called after htape_wrapup call in order to use the old
  ! old weights in December's history calculations.
  ! ============================================================================

  call get_curr_date(yrp1, monp1, dayp1, secp1, offset=int(dtime))
  if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then

     ! Get date info.  kyr is used in lpj().  At end of first year, kyr = 2.
     call get_curr_date(yr, mon, day, sec)
     ncdate = yr*10000 + mon*100 + day
     call get_ref_date(yr, mon, day, sec)
     nbdate = yr*10000 + mon*100 + day
     kyr = ncdate/10000 - nbdate/10000 + 1

     if (masterproc) write(6,*) 'End of year. DGVM called now: ncdate=', &
                     ncdate,' nbdate=',nbdate,' kyr=',kyr,' nstep=', nstep

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
     do nc = 1,nclumps
        call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
        call lpj(begg, endg, begp, endp, filter(nc)%num_natvegp, filter(nc)%natvegp, kyr)
        call lpjreset1(begg, endg, begc, endc, begp, endp, filter(nc)%num_nolakep, filter(nc)%nolakep, &
                       caldayp1, eccen, obliqr, lambm0, mvelpp)
    end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO
  end if
#endif

  call t_stopf('clm_driver')
  call t_startf('clm_driver_io')

  ! ============================================================================
  ! Create history and write history tapes if appropriate
  ! ============================================================================

#if ( !defined SCAM )

  call t_startf('wrapup')
  call htapes_wrapup()
  call t_stopf('wrapup')

#endif

#if (defined DGVM)
  ! ============================================================================
  ! Finish DGVM calculation
  ! ============================================================================

  if (monp1==1 .and. dayp1==1 .and. secp1==dtime .and. nstep>0)  then
!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
     do nc = 1,nclumps
        call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)
        call lpjreset2(begg, endg, begl, endl, begc, endc, begp, endp)
     end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO
#ifdef SPMD
     call gatherWeightsDGVM()
#endif
     call histDGVM()
     if (masterproc) write(6,*) 'Annual DGVM calculations are complete'
  end if
#endif

#if ( !defined SCAM )
  ! ============================================================================
  ! Write restart files if appropriate
  ! ============================================================================

  if (do_restwrite()) call restart('write')

  ! ============================================================================
  ! Write intial files if appropriate
  ! ============================================================================

  if (do_inicwrite()) call inicfile(flag='write')

#endif 

  call t_stopf('clm_driver_io')

end subroutine driver

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_diagnostic
!
! !INTERFACE:
subroutine write_diagnostic (wrtdia, nstep)
!
! !DESCRIPTION:
! Write diagnostic surface temperature output each timestep.  Written to
! be fast but not bit-for-bit because order of summations can change each
! timestep.
!
! !USES:
  use clmtype
  use decompMod  , only : get_proc_bounds, get_proc_global
#if (defined SPMD)
  use spmdMod    , only : masterproc, npes, MPI_REAL8, MPI_ANY_SOURCE, &
                          MPI_STATUS_SIZE, mpicom
#else
  use spmdMod    , only : masterproc
#endif
  use shr_sys_mod, only : shr_sys_flush
  use abortutils , only : endrun
!
! !ARGUMENTS:
  implicit none
  logical, intent(in) :: wrtdia     !true => write diagnostic
  integer, intent(in) :: nstep      !model time step
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
  integer :: p                       ! loop index
  integer :: begp, endp              ! per-proc beginning and ending pft indices
  integer :: begc, endc              ! per-proc beginning and ending column indices
  integer :: begl, endl              ! per-proc beginning and ending landunit indices
  integer :: begg, endg              ! per-proc gridcell ending gridcell indices
  integer :: numg                    ! total number of gridcells across all processors
  integer :: numl                    ! total number of landunits across all processors
  integer :: numc                    ! total number of columns across all processors
  integer :: nump                    ! total number of pfts across all processors
  integer :: ier                     ! error status
  real(r8):: psum                    ! partial sum of ts
  real(r8):: tsum                    ! sum of ts
  real(r8):: tsxyav                  ! average ts for diagnostic output
#if (defined SPMD)
  integer :: status(MPI_STATUS_SIZE) ! mpi status
#endif
!------------------------------------------------------------------------

#if (!defined COUP_CAM)

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)

  if (wrtdia) then

#if (defined TIMING_BARRIERS)
     call t_startf ('sync_write_diag')
     call mpi_barrier (mpicom, ier)
     call t_stopf ('sync_write_diag')
#endif
     psum = sum(clm3%g%l2as%t_rad(begg:endg))
#if (defined SPMD)
     if (masterproc) then
        tsum = psum
        do p = 1, npes-1
           call mpi_recv(psum, 1, MPI_REAL8, p, 999, mpicom, status, ier)
           if (ier/=0) then
              write(6,*) 'write_diagnostic: Error in mpi_recv()'
              call endrun
           end if
           tsum = tsum + psum
        end do
     else
        call mpi_send(psum, 1, MPI_REAL8, 0, 999, mpicom, ier)
        if (ier/=0) then
           write(6,*) 'write_diagnostic: Error in mpi_send()'
           call endrun
        end if
     end if
#else
     tsum = psum
#endif
     if (masterproc) then
        tsxyav = tsum / numg
        write (6,1000) nstep, tsxyav
        call shr_sys_flush(6)
     end if

  else

     if (masterproc) then
        write(6,*)'clm2: completed timestep ',nstep
        call shr_sys_flush(6)
     end if

  endif

1000 format (1x,'nstep = ',i10,'   TS = ',e21.15)

#endif

end subroutine write_diagnostic
