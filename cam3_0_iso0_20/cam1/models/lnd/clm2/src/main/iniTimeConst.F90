#include <misc.h>
#include <preproc.h>

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: iniTimeConst
!
! !INTERFACE:
subroutine iniTimeConst
!
! !DESCRIPTION:
! Initialize time invariant clm variables
! 1) removed references to shallow lake - since it is not used
! 2) ***Make c%z, c%zi and c%dz allocatable depending on if you
!    have lake or soil
! 3) rootfr only initialized for soil points
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use nanMod
  use clmtype
  use decompMod , only : get_proc_bounds, get_proc_global
  use clm_varpar, only : nlevsoi, nlevlak, lsmlon, lsmlat, numpft
  use clm_varsur, only : soic2d, sand3d, clay3d
  use clm_varcon, only : istice, istdlak, istwet, isturb, &
                         zlak, dzlak, zsoi, dzsoi, zisoi, spval
  use clm_varctl, only : nsrest
  use pftvarcon , only : ncorn, nwheat, noveg, ntree, roota_par, rootb_par,  &
                         z0mr, displar, dleaf, rhol, rhos, taul, taus, xl, &
                         qe25, vcmx25, mp, c3psn, &
                         pftpar , tree   , summergreen, raingreen  , sla     , &
                         lm_sapl, sm_sapl, hm_sapl    , rm_sapl    , latosa  , &
                         allom1 , allom2 , allom3     , reinickerp , wooddens
  use time_manager, only : get_step_size
  use abortutils, only : endrun
!
! !ARGUMENTS:
  implicit none
!
! !CALLED FROM:
! subroutine initialize in module initializeMod.
!
! !REVISION HISTORY:
! Created by Gordon Bonan.
! Updated to clm2.1 data structrues by Mariana Vertenstein
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
  integer , pointer :: ivt(:)             !  vegetation type index
  integer , pointer :: ixy(:)             ! xy lon index (column-level)
  integer , pointer :: jxy(:)             ! xy lat index (column_level)
  integer , pointer :: pcolumn(:)         ! column index of corresponding pft
  integer , pointer :: clandunit(:)       ! landunit index of corresponding column
  integer , pointer :: ltype(:)           ! landunit type index
!
! local pointers to implicit out arguments
!
  real(r8), pointer :: z(:,:)             ! layer depth (m)
  real(r8), pointer :: zi(:,:)            ! interface level below a "z" level (m)
  real(r8), pointer :: dz(:,:)            ! layer thickness depth (m)
  real(r8), pointer :: rootfr(:,:)        ! fraction of roots in each soil layer
  real(r8), pointer :: dewmx(:)           ! maximum allowed dew [mm]
  real(r8), pointer :: bsw(:,:)           ! Clapp and Hornberger "b" (nlevsoi)
  real(r8), pointer :: watsat(:,:)        ! volumetric soil water at saturation (porosity) (nlevsoi)
  real(r8), pointer :: hksat(:,:)         ! hydraulic conductivity at saturation (mm H2O /s) (nlevsoi)
  real(r8), pointer :: sucsat(:,:)        ! minimum soil suction (mm) (nlevsoi)
  real(r8), pointer :: csol(:,:)          ! heat capacity, soil solids (J/m**3/Kelvin) (nlevsoi)
  real(r8), pointer :: tkmg(:,:)          ! thermal conductivity, soil minerals  [W/m-K] (new) (nlevsoi)
  real(r8), pointer :: tkdry(:,:)         ! thermal conductivity, dry soil (W/m/Kelvin) (nlevsoi)
  real(r8), pointer :: tksatu(:,:)        ! thermal conductivity, saturated soil [W/m-K] (new) (nlevsoi)
  real(r8), pointer :: wtfact(:)          ! Fraction of model area with high water table
  real(r8), pointer :: smpmin(:)          ! restriction for min of soil potential (mm) (new)
  integer , pointer :: isoicol(:)         ! soil color class
  real(r8), pointer :: gwc_thr(:)         ! threshold soil moisture based on clay content
  real(r8), pointer :: mss_frc_cly_vld(:) ! [frc] Mass fraction clay limited to 0.20
!
!EOP
!
! !OTHER LOCAL VARIABLES:
  integer :: i,j,ib,lev       ! indices
  integer :: g,l,c,p          ! indices
  integer :: m                ! vegetation type index
  real(r8):: bd               ! bulk density of dry soil material [kg/m^3]
  real(r8):: tkm              ! mineral conductivity
  real(r8):: xksat            ! maximum hydraulic conductivity of soil [mm/s]
  real(r8):: scalez = 0.025   ! Soil layer thickness discretization (m)
  real(r8):: hkdepth = 0.5    ! Length scale for Ksat decrease (m)
  real(r8):: clay,sand        ! temporaries
  integer :: begp, endp       ! per-proc beginning and ending pft indices
  integer :: begc, endc       ! per-proc beginning and ending column indices
  integer :: begl, endl       ! per-proc beginning and ending landunit indices
  integer :: begg, endg       ! per-proc gridcell ending gridcell indices
  integer :: numg             ! total number of gridcells across all processors
  integer :: numl             ! total number of landunits across all processors
  integer :: numc             ! total number of columns across all processors
  integer :: nump             ! total number of pfts across all processors
!------------------------------------------------------------------------

  ! Assign local pointers to derived subtypes components (gridcell-level)

  wtfact          => clm3%g%gps%wtfact

  ! Assign local pointers to derived subtypes components (landunit-level)

  ltype           => clm3%g%l%itype

  ! Assign local pointers to derived subtypes components (column-level)

  ixy             => clm3%g%l%c%ixy
  jxy             => clm3%g%l%c%jxy
  clandunit       => clm3%g%l%c%landunit
  z               => clm3%g%l%c%cps%z
  dz              => clm3%g%l%c%cps%dz
  zi              => clm3%g%l%c%cps%zi
  bsw             => clm3%g%l%c%cps%bsw
  watsat          => clm3%g%l%c%cps%watsat
  hksat           => clm3%g%l%c%cps%hksat
  sucsat          => clm3%g%l%c%cps%sucsat
  tkmg            => clm3%g%l%c%cps%tkmg
  tksatu          => clm3%g%l%c%cps%tksatu
  tkdry           => clm3%g%l%c%cps%tkdry
  csol            => clm3%g%l%c%cps%csol
  smpmin          => clm3%g%l%c%cps%smpmin
  isoicol         => clm3%g%l%c%cps%isoicol
  gwc_thr         => clm3%g%l%c%cps%gwc_thr
  mss_frc_cly_vld => clm3%g%l%c%cps%mss_frc_cly_vld

  ! Assign local pointers to derived subtypes components (pft-level)

  ivt             => clm3%g%l%c%p%itype
  pcolumn         => clm3%g%l%c%p%column
  dewmx           => clm3%g%l%c%p%pps%dewmx
  rootfr          => clm3%g%l%c%p%pps%rootfr

  ! Determine necessary subgrid bounds

  call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
  call get_proc_global(numg, numl, numc, nump)

  ! --------------------------------------------------------------------
  ! Initialize time constant arrays of ecophysiological constants and
  ! arrays of dgvm ecophysiological constants
  ! --------------------------------------------------------------------

!dir$ concurrent
!cdir nodep
   do m = 0,numpft
      pftcon%ncorn(m) = ncorn
      pftcon%nwheat(m) = nwheat
      pftcon%noveg(m) = noveg
      pftcon%ntree(m) = ntree
      pftcon%z0mr(m) = z0mr(m)
      pftcon%displar(m) = displar(m)
      pftcon%dleaf(m) = dleaf(m)
      pftcon%xl(m) = xl(m)
      do ib = 1,numrad
         pftcon%rhol(m,ib) = rhol(m,ib)
         pftcon%rhos(m,ib) = rhos(m,ib)
         pftcon%taul(m,ib) = taul(m,ib)
         pftcon%taus(m,ib) = taus(m,ib)
      end do
      pftcon%qe25(m) = qe25(m)
      pftcon%vcmx25(m) = vcmx25(m)
      pftcon%mp(m) = mp(m)
      pftcon%c3psn(m) = c3psn(m)
      pftcon%sla(m) = sla(m)
   end do

!dir$ concurrent
!cdir nodep
   do m = 0,numpft
      dgv_pftcon%respcoeff(m) = pftpar(m,5)
      dgv_pftcon%flam(m) = pftpar(m,6)
      dgv_pftcon%resist(m) = pftpar(m,8)
      dgv_pftcon%l_turn(m) = pftpar(m,9)
      dgv_pftcon%l_long(m) = pftpar(m,10)
      dgv_pftcon%s_turn(m) = pftpar(m,11)
      dgv_pftcon%r_turn(m) = pftpar(m,12)
      dgv_pftcon%l_cton(m) = pftpar(m,13)
      dgv_pftcon%s_cton(m) = pftpar(m,14)
      dgv_pftcon%r_cton(m) = pftpar(m,15)
      dgv_pftcon%l_morph(m) = pftpar(m,16)
      dgv_pftcon%l_phen(m) = pftpar(m,17)
      dgv_pftcon%lmtorm(m) = pftpar(m,18)
      dgv_pftcon%crownarea_max(m) = pftpar(m,20)
      dgv_pftcon%init_lai(m) = pftpar(m,21)
      dgv_pftcon%x(m) = pftpar(m,22)
      dgv_pftcon%tcmin(m) = pftpar(m,28)
      dgv_pftcon%tcmax(m) = pftpar(m,29)
      dgv_pftcon%gddmin(m) = pftpar(m,30)
      dgv_pftcon%twmax(m) = pftpar(m,31)
      dgv_pftcon%lm_sapl(m) = lm_sapl(m)
      dgv_pftcon%sm_sapl(m) = sm_sapl(m)
      dgv_pftcon%hm_sapl(m) = hm_sapl(m)
      dgv_pftcon%rm_sapl(m) = rm_sapl(m)
      dgv_pftcon%tree(m) = tree(m)
      dgv_pftcon%summergreen(m) = summergreen(m)
      dgv_pftcon%raingreen(m) = raingreen(m)
      dgv_pftcon%reinickerp(m) = reinickerp
      dgv_pftcon%wooddens(m) = wooddens
      dgv_pftcon%latosa(m) = latosa
      dgv_pftcon%allom1(m) = allom1
      dgv_pftcon%allom2(m) = allom2
      dgv_pftcon%allom3(m) = allom3
   end do

   ! --------------------------------------------------------------------
   ! Define layer structure for soil and lakes
   ! Vertical profile of snow is initialized in routine iniTimeVar
   ! --------------------------------------------------------------------

   ! check that lake and soil levels are the same for now

   if (nlevlak /= nlevsoi) then
      write(6,*)'number of soil levels and number of lake levels must be the same'
      write(6,*)'nlevsoi= ',nlevsoi,' nlevlak= ',nlevlak
      call endrun
   endif

   ! Lake layers (assumed same for all lake patches)

   dzlak(1) = 0.1
   dzlak(2) = 1.
   dzlak(3) = 2.
   dzlak(4) = 3.
   dzlak(5) = 4.
   dzlak(6) = 5.
   dzlak(7) = 7.
   dzlak(8) = 7.
   dzlak(9) = 10.45
   dzlak(10)= 10.45

   zlak(1) =  0.05
   zlak(2) =  0.6
   zlak(3) =  2.1
   zlak(4) =  4.6
   zlak(5) =  8.1
   zlak(6) = 12.6
   zlak(7) = 18.6
   zlak(8) = 25.6
   zlak(9) = 34.325
   zlak(10)= 44.775

   ! Soil layers and interfaces (assumed same for all non-lake patches)
   ! "0" refers to soil surface and "nlevsoi" refers to the bottom of model soil

   do j = 1, nlevsoi
      zsoi(j) = scalez*(exp(0.5*(j-0.5))-1.)    !node depths
   enddo

   dzsoi(1) = 0.5*(zsoi(1)+zsoi(2))             !thickness b/n two interfaces
   do j = 2,nlevsoi-1
      dzsoi(j)= 0.5*(zsoi(j+1)-zsoi(j-1))
   enddo
   dzsoi(nlevsoi) = zsoi(nlevsoi)-zsoi(nlevsoi-1)

   zisoi(0) = 0.
   do j = 1, nlevsoi-1
      zisoi(j) = 0.5*(zsoi(j)+zsoi(j+1))         !interface depths
   enddo
   zisoi(nlevsoi) = zsoi(nlevsoi) + 0.5*dzsoi(nlevsoi)

   ! --------------------------------------------------------------------
   ! Initialize soil and lake levels
   ! Initialize soil color, thermal and hydraulic properties
   ! --------------------------------------------------------------------

   ! Grid level initialization
   do g = begg, endg

      ! Initialize fraction of model area with high water table
      wtfact(g) = 0.3

   end do

   ! Column level initialization
!dir$ concurrent
!cdir nodep
   do c = begc, endc

      ! Set gridcell and landunit indices
      i = ixy(c)
      j = jxy(c)
      l = clandunit(c)

      ! Initialize restriction for min of soil potential (mm)
      smpmin(c) = -1.e8

      ! Soil color
      isoicol(c) = soic2d(i,j)

      ! Soil hydraulic and thermal properties
      if (ltype(l)==istdlak .or. ltype(l)==istwet .or. &
          ltype(l)==istice .or. ltype(l)==isturb ) then
         do lev = 1,nlevsoi
            bsw(c,lev) = spval
            watsat(c,lev) = spval
            hksat(c,lev) = spval
            sucsat(c,lev) = spval
            tkmg(c,lev) = spval
            tksatu(c,lev) = spval
            tkdry(c,lev) = spval
            csol(c,lev) = spval
         end do
      else
         do lev = 1,nlevsoi
            clay = clay3d(i,j,lev)
            sand = sand3d(i,j,lev)
            watsat(c,lev) = 0.489 - 0.00126*sand
            bd = (1.-watsat(c,lev))*2.7e3
            xksat = 0.0070556 *( 10.**(-0.884+0.0153*sand) ) ! mm/s
            tkm = (8.80*sand+2.92*clay)/(sand+clay)          ! W/(m K)

            bsw(c,lev) = 2.91 + 0.159*clay
            hksat(c,lev) = xksat * exp(-zisoi(lev)/hkdepth)
            sucsat(c,lev) = 10. * ( 10.**(1.88-0.0131*sand) )
            tkmg(c,lev) = tkm ** (1.- watsat(c,lev))
            tksatu(c,lev) = tkmg(c,lev)*0.57**watsat(c,lev)
            tkdry(c,lev) = (0.135*bd + 64.7) / (2.7e3 - 0.947*bd)
            csol(c,lev) = (2.128*sand+2.385*clay) / (sand+clay)*1.e6  ! J/(m3 K)
         end do
      endif

      ! Define lake or non-lake levels layers
      if (ltype(l) == istdlak) then
         z(c,1:nlevlak) = zlak(1:nlevlak)
         dz(c,1:nlevlak) = dzlak(1:nlevlak)
      else
         z(c,1:nlevsoi) = zsoi(1:nlevsoi)
         zi(c,0:nlevsoi) = zisoi(0:nlevsoi)
         dz(c,1:nlevsoi) = dzsoi(1:nlevsoi)
      end if

      ! Initialize terms needed for dust model
      clay = clay3d(i,j,1)
      gwc_thr(c) = 0.17 + 0.14*clay*0.01
      mss_frc_cly_vld(c) = min(clay*0.01_r8, 0.20_r8)

   end do

   ! pft level initialization
!dir$ concurrent
!cdir nodep
   do p = begp, endp

      ! Initialize maximum allowed dew

      dewmx(p)  = 0.1

      ! Initialize root fraction (computing from surface, d is depth in meter):
      ! Y = 1 -1/2 (exp(-ad)+exp(-bd) under the constraint that
      ! Y(d =0.1m) = 1-beta^(10 cm) and Y(d=d_obs)=0.99 with
      ! beta & d_obs given in Zeng et al. (1998).

      c = pcolumn(p)
      if (ivt(p) /= noveg) then
         do lev = 1, nlevsoi-1
            rootfr(p,lev) = .5*( exp(-roota_par(ivt(p)) * zi(c,lev-1))  &
                               + exp(-rootb_par(ivt(p)) * zi(c,lev-1))  &
                               - exp(-roota_par(ivt(p)) * zi(c,lev  ))  &
                               - exp(-rootb_par(ivt(p)) * zi(c,lev  )) )
         end do
         rootfr(p,nlevsoi) = .5*( exp(-roota_par(ivt(p)) * zi(c,nlevsoi-1))  &
                                + exp(-rootb_par(ivt(p)) * zi(c,nlevsoi-1)) )
      else
         rootfr(p,1:nlevsoi) = spval
      endif

   end do ! end pft level initialization

end subroutine iniTimeConst
