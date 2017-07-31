#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: iniTimeVar
!
! !INTERFACE:
subroutine iniTimeVar(readini, eccen, obliqr, lambm0 , mvelpp)
!
! !DESCRIPTION:
! Initializes the following time varying variables:
! water : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
! snow : snowdp, snowage, snl, dz, z, zi
! temperature: t_soisno, t_veg, t_grnd
! The variable, h2osoi_vol, is needed by the soil albedo routine - this is not needed
! on restart since it is computed before the soil albedo computation is called.
! The remaining variables are initialized by calls to ecosystem dynamics and
! albedo subroutines.
!
! !USES:
  use shr_kind_mod         , only : r8 => shr_kind_r8
  use clmtype
  use spmdMod              , only : masterproc
  use decompMod            , only : get_proc_clumps, get_clump_bounds
  use filterMod            , only : filter
  use clm_varpar           , only : nlevsoi, nlevsno, nlevlak
  use clm_varcon           , only : denice, denh2o, zlnd
  use time_manager         , only : get_curr_calday
  use inicFileMod          , only : inicfile
  use FracWetMod           , only : FracWet
  use SurfaceAlbedoMod     , only : SurfaceAlbedo
#ifdef DGVM
  use DGVMMod              , only : resetTimeConstDGVM
  use DGVMEcosystemDynMod  , only : DGVMEcosystemDyn
#else
  use STATICEcosysDynMod, only : EcosystemDyn, interpMonthlyVeg
#endif
!
! !ARGUMENTS:
  implicit none
  logical , intent(in) :: readini  !true if read in initial data set
  real(r8), intent(in) :: eccen    !Earth's orbital eccentricity
  real(r8), intent(in) :: obliqr   !Earth's obliquity in radians
  real(r8), intent(in) :: lambm0   !Mean longitude of perihelion at the vernal equinox (radians)
  real(r8), intent(in) :: mvelpp   !Earth's moving vernal equinox long. of perihelion + pi (radians)
!
! !CALLED FROM:
! subroutine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: plandunit(:)      ! landunit index associated with each pft
  logical , pointer :: lakpoi(:)         ! true => landunit is a lake point
  real(r8), pointer :: dz(:,:)           ! layer thickness depth (m)
  real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
  real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
  integer , pointer :: frac_veg_nosno_alb(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
!
! local pointers to implicit out arguments
!
  real(r8), pointer :: h2osoi_vol(:,:)   ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
  real(r8), pointer :: snowdp(:)         ! snow height (m)
  real(r8), pointer :: frac_sno(:)       ! fraction of ground covered by snow (0 to 1)
  integer , pointer :: frac_veg_nosno(:) ! fraction of vegetation not covered by snow (0 OR 1) [-]
  real(r8), pointer :: fwet(:)           ! fraction of canopy that is wet (0 to 1) (pft-level)
!
! local pointers to implicit out arguments (lake points only)
!
  real(r8), pointer :: fdry(:)     ! fraction of foliage that is green and dry [-] (new)
  real(r8), pointer :: tlai(:)     ! one-sided leaf area index, no burying by snow
  real(r8), pointer :: tsai(:)     ! one-sided stem area index, no burying by snow
  real(r8), pointer :: htop(:)     ! canopy top (m)
  real(r8), pointer :: hbot(:)     ! canopy bottom (m)
  real(r8), pointer :: elai(:)     ! one-sided leaf area index with burying by snow
  real(r8), pointer :: esai(:)     ! one-sided stem area index with burying by snow
!
!EOP
!
! !OTHER LOCAL VARIABLES:
  integer :: nc,j,l,c,p   ! indices
  integer :: nclumps      ! number of clumps on this processor
  integer :: begp, endp   ! per-clump beginning and ending pft indices
  integer :: begc, endc   ! per-clump beginning and ending column indices
  integer :: begl, endl   ! per-clump beginning and ending landunit indices
  integer :: begg, endg   ! per-clump gridcell ending gridcell indices
  real(r8):: calday       ! calendar day
!-----------------------------------------------------------------------

  ! Assign local pointers to derived subtypes components (landunit-level)

  lakpoi              => clm3%g%l%lakpoi

  ! Assign local pointers to derived subtypes components (column-level)

  dz                  => clm3%g%l%c%cps%dz
  h2osoi_ice          => clm3%g%l%c%cws%h2osoi_ice
  h2osoi_liq          => clm3%g%l%c%cws%h2osoi_liq
  h2osoi_vol          => clm3%g%l%c%cws%h2osoi_vol
  snowdp              => clm3%g%l%c%cps%snowdp
  frac_sno            => clm3%g%l%c%cps%frac_sno

  ! Assign local pointers to derived subtypes components (pft-level)

  plandunit          => clm3%g%l%c%p%landunit
  frac_veg_nosno_alb => clm3%g%l%c%p%pps%frac_veg_nosno_alb
  frac_veg_nosno     => clm3%g%l%c%p%pps%frac_veg_nosno
  fwet               => clm3%g%l%c%p%pps%fwet

  ! Assign local pointers to derived subtypes components (pft-level)
  ! The folowing pointers will only be used for lake points in this routine

  htop               => clm3%g%l%c%p%pps%htop
  hbot               => clm3%g%l%c%p%pps%hbot
  tlai               => clm3%g%l%c%p%pps%tlai
  tsai               => clm3%g%l%c%p%pps%tsai
  elai               => clm3%g%l%c%p%pps%elai
  esai               => clm3%g%l%c%p%pps%esai
  fdry               => clm3%g%l%c%p%pps%fdry

  ! ========================================================================
  ! Initialize water and temperature based on:
  ! readini = true : read initial data set -- requires netCDF codes
  ! readini = false: arbitrary initialization
  ! ========================================================================

  if (readini) then
     if ( masterproc ) write (6,*) 'Reading initial data from initial dataset'
     call inicfile(flag='read')
  else
     if ( masterproc ) write (6,*) 'Setting initial data to non-spun up values'
     call mkarbinit()
  end if

  ! ========================================================================
  ! Remaining variables are initialized by calls to ecosystem dynamics and
  ! albedo subroutines.
  ! Note: elai, esai, frac_veg_nosno_alb are computed in
  ! Ecosysdyn and needed by routines FracWet and SurfaceAlbedo
  ! frac_veg_nosno is needed by FracWet
  ! fwet is needed in routine TwoStream (called by SurfaceAlbedo)
  ! frac_sno is needed by SoilAlbedo (called by SurfaceAlbedo)
  ! ========================================================================

  ! Determine current calday

  calday = get_curr_calday()

#if (!defined DGVM)
  ! Read monthly vegetation data for interpolation to daily values

  call interpMonthlyVeg()
#endif

  ! Determine clump bounds for this processor

  nclumps = get_proc_clumps()

  ! Loop over clumps on this processor

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp,p,c,l)
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp,p,c,l)
  do nc = 1,nclumps

     ! Determine clump bounds

     call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

     ! Determine variables needed by SurfaceAlbedo for lake points

!dir$ concurrent
!cdir nodep
     do p = begp,endp
        l = plandunit(p)
        if (lakpoi(l)) then
           fwet(p) = 0.
           fdry(p) = 0.
           elai(p) = 0.
           esai(p) = 0.
           htop(p) = 0.
           hbot(p) = 0.
           tlai(p) = 0.
           tsai(p) = 0.
           frac_veg_nosno_alb(p) = 0.
           frac_veg_nosno(p) = 0.
        end if
     end do

     ! Determine variables needed for SurfaceAlbedo for non-lake points

#ifdef DGVM
     call resetTimeConstDGVM(begp, endp)
     call DGVMEcosystemDyn(begp, endp, filter(nc)%num_nolakep, filter(nc)%nolakep, &
                           doalb=.true., endofyr=.false.)
#else
     call EcosystemDyn(begp, endp, filter(nc)%num_nolakep, filter(nc)%nolakep, &
                       doalb=.true.)
#endif

!dir$ concurrent
!cdir nodep
     do p = begp, endp
        l = plandunit(p)
        if (.not. lakpoi(l)) then
           frac_veg_nosno(p) = frac_veg_nosno_alb(p)
           fwet(p) = 0.
        end if
     end do

     call FracWet(filter(nc)%num_nolakep, filter(nc)%nolakep)

     ! Compute Surface Albedo - all land points (including lake)
     ! Needs as input fracion of soil covered by snow (Z.-L. Yang U. Texas)

!dir$ concurrent
!cdir nodep
     do c = begc, endc
        frac_sno(c) = snowdp(c)/(10.*zlnd + snowdp(c))
     end do

     call SurfaceAlbedo(begg, endg, begc, endc, begp, endp, &
                        calday, eccen, obliqr, lambm0, mvelpp)

  end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

end subroutine iniTimeVar

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: mkarbinit
!
! !INTERFACE:
subroutine mkarbinit()
!
! !DESCRIPTION:
! Initializes the following time varying variables:
! water      : h2osno, h2ocan, h2osoi_liq, h2osoi_ice, h2osoi_vol
! snow       : snowdp, snowage, snl, dz, z, zi
! temperature: t_soisno, t_veg, t_grnd
! The variable, h2osoi_vol, is needed by clm_soilalb -this is not needed on
! restart since it is computed before the soil albedo computation is called.
! The remaining variables are initialized by calls to ecosystem dynamics
! and albedo subroutines.
!
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use clmtype
  use decompMod    , only : get_proc_bounds
  use clm_varpar   , only : nlevsoi, nlevsno, nlevlak
  use clm_varcon   , only : bdsno, istice, istwet, istsoil, denice, denh2o, spval, sb
  use shr_const_mod, only : SHR_CONST_TKFRZ
#ifdef DGVM
  use DGVMMod      , only : resetWeightsDGVM, gatherWeightsDGVM
#endif
!
! !ARGUMENTS:
  implicit none
!
! !CALLED FROM:
! subroutine iniTimeVar
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: pcolumn(:)        ! column index associated with each pft
  integer , pointer :: clandunit(:)      ! landunit index associated with each column
  integer , pointer :: ltype(:)          ! landunit type
  logical , pointer :: lakpoi(:)         ! true => landunit is a lake point
  real(r8), pointer :: dz(:,:)           ! layer thickness depth (m)
  real(r8), pointer :: watsat(:,:)       ! volumetric soil water at saturation (porosity) (nlevsoi)
  real(r8), pointer :: h2osoi_ice(:,:)   ! ice lens (kg/m2)
  real(r8), pointer :: h2osoi_liq(:,:)   ! liquid water (kg/m2)
!
! local pointers to implicit out arguments
!
  integer , pointer :: snl(:)            ! number of snow layers
  real(r8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
  real(r8), pointer :: t_lake(:,:)       ! lake temperature (Kelvin)  (1:nlevlak)
  real(r8), pointer :: t_grnd(:)         ! ground temperature (Kelvin)
  real(r8), pointer :: t_veg(:)          ! vegetation temperature (Kelvin)
  real(r8), pointer :: h2osoi_vol(:,:)   ! volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
  real(r8), pointer :: h2ocan_col(:)     ! canopy water (mm H2O) (column-level)
  real(r8), pointer :: h2ocan_pft(:)     ! canopy water (mm H2O) (pft-level)
  real(r8), pointer :: h2osno(:)         ! snow water (mm H2O)
  real(r8), pointer :: snowdp(:)         ! snow height (m)
  real(r8), pointer :: snowage(:)        ! non dimensional snow age [-] (new)
  real(r8), pointer :: eflx_lwrad_out(:) ! emitted infrared (longwave) radiation (W/m**2)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
  integer :: j,l,c,p      ! indices
  integer :: begp, endp   ! per-proc beginning and ending pft indices
  integer :: begc, endc   ! per-proc beginning and ending column indices
  integer :: begl, endl   ! per-proc beginning and ending landunit indices
  integer :: begg, endg   ! per-proc gridcell ending gridcell indices
!-----------------------------------------------------------------------

  ! Assign local pointers to derived subtypes components (landunit-level)

  ltype => clm3%g%l%itype
  lakpoi => clm3%g%l%lakpoi

  ! Assign local pointers to derived subtypes components (column-level)

  clandunit  => clm3%g%l%c%landunit
  snl        => clm3%g%l%c%cps%snl
  dz         => clm3%g%l%c%cps%dz
  watsat     => clm3%g%l%c%cps%watsat
  h2osoi_ice => clm3%g%l%c%cws%h2osoi_ice
  h2osoi_liq => clm3%g%l%c%cws%h2osoi_liq
  h2osoi_vol => clm3%g%l%c%cws%h2osoi_vol
  h2ocan_col => clm3%g%l%c%cws%pws_a%h2ocan
  snowage    => clm3%g%l%c%cps%snowage
  snowdp     => clm3%g%l%c%cps%snowdp
  h2osno     => clm3%g%l%c%cws%h2osno
  t_soisno   => clm3%g%l%c%ces%t_soisno
  t_lake     => clm3%g%l%c%ces%t_lake
  t_grnd     => clm3%g%l%c%ces%t_grnd

  ! Assign local pointers to derived subtypes components (pft-level)

  pcolumn => clm3%g%l%c%p%column
  h2ocan_pft => clm3%g%l%c%p%pws%h2ocan
  t_veg => clm3%g%l%c%p%pes%t_veg
  eflx_lwrad_out => clm3%g%l%c%p%pef%eflx_lwrad_out

  ! Determine subgrid bounds on this processor

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

  ! ========================================================================
  ! Set snow water
  ! ========================================================================

  ! NOTE: h2ocan, h2osno, snowdp and snowage has valid values everywhere

  ! canopy water (pft level)

  do p = begp, endp
     h2ocan_pft(p) = 0._r8
  end do

!dir$ concurrent
!cdir nodep
  do c = begc,endc

     ! canopy water (column level)

     h2ocan_col(c) = 0._r8

     ! snow water

     l = clandunit(c)
     if (ltype(l) == istice) then
        h2osno(c) = 1000._r8
     else
        h2osno(c) = 0._r8
     endif

     ! snow depth

     snowdp(c)  = h2osno(c) / bdsno

     ! snow age

     snowage(c) = 0.

  end do

  ! ========================================================================
  ! Set snow layer number, depth and thickiness
  ! ========================================================================

  call snowdp2lev(begc, endc)

  ! ========================================================================
  ! Set snow/soil temperature
  ! ========================================================================

  ! NOTE:
  ! t_soisno only has valid values over non-lake
  ! t_lake   only has valid values over lake
  ! t_grnd has valid values over all land
  ! t_veg  has valid values over all land

!dir$ concurrent
!cdir nodep
  do c = begc,endc

     t_soisno(c,-nlevsno+1:nlevsoi) = spval
     t_lake(c,1:nlevlak) = spval

     l = clandunit(c)
     if (.not. lakpoi(l)) then  !not lake
        t_soisno(c,-nlevsno+1:0) = spval
        if (snl(c) < 0) then    !snow layer temperatures
           do j = snl(c)+1, 0
              t_soisno(c,j) = 250._r8
           enddo
        endif
        if (ltype(l) == istice) then
           do j = 1, nlevsoi
              t_soisno(c,j) = 250._r8
           end do
        else if (ltype(l) == istwet) then
           do j = 1, nlevsoi
              t_soisno(c,j) = 277._r8
           end do
        else
           do j = 1, nlevsoi
              t_soisno(c,j) = 283._r8
           end do
        endif
        t_grnd(c) = t_soisno(c,snl(c)+1)
     else                     !lake
        t_lake(c,1:nlevlak) = 277._r8
        t_grnd(c) = t_lake(c,1)
     endif

  end do

!dir$ concurrent
!cdir nodep
  do p = begp, endp
     c = pcolumn(p)
     t_veg(p) = 283._r8
     eflx_lwrad_out(p) = sb * (t_grnd(c))**4
  end do

  ! ========================================================================
  ! Set snow/soil ice and liquid mass
  ! ========================================================================

  ! volumetric water is set first and liquid content and ice lens are obtained
  ! NOTE: h2osoi_vol, h2osoi_liq and h2osoi_ice only have valid values over soil

  h2osoi_vol(begc:endc,         1:nlevsoi) = spval
  h2osoi_liq(begc:endc,-nlevsno+1:nlevsoi) = spval
  h2osoi_ice(begc:endc,-nlevsno+1:nlevsoi) = spval

  do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
     do c = begc,endc
        l = clandunit(c)
        if (.not. lakpoi(l)) then  !not lake
           ! volumetric water
           if (ltype(l) == istsoil) then
              h2osoi_vol(c,j) = 0.3_r8
           else
              h2osoi_vol(c,j) = 1.0_r8
           endif
           h2osoi_vol(c,j) = min(h2osoi_vol(c,j),watsat(c,j))
           ! soil layers
           if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
              h2osoi_ice(c,j)  = dz(c,j)*denice*h2osoi_vol(c,j)
              h2osoi_liq(c,j) = 0._r8
           else
              h2osoi_ice(c,j) = 0._r8
              h2osoi_liq(c,j) = dz(c,j)*denh2o*h2osoi_vol(c,j)
           endif
        end if
     end do
  end do

  do j = -nlevsno+1, 0
!dir$ concurrent
!cdir nodep
     do c = begc,endc
        l = clandunit(c)
        if (.not. lakpoi(l)) then  !not lake
           ! snow
           if (j > snl(c)) then
              h2osoi_ice(c,j) = dz(c,j)*250.
              h2osoi_liq(c,j) = 0._r8
           end if
        end if
     end do
  end do

#ifdef DGVM
  ! Determine new subgrid weights and areas (obtained
  ! from new value of fpcgrid read in above) - this is needed
  ! here to avoid round off level errors on restart before
  ! lpj is called the first time

  call resetWeightsDGVM(begg, endg, begc, endc, begp, endp)
#ifdef SPMD
  call gatherWeightsDGVM()
#endif
#endif

end subroutine mkarbinit
