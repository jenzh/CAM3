#include <misc.h>
#include <preproc.h>

module clmtypeInitMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clmtypeInitMod
!
! !DESCRIPTION:
! Allocate clmtype components and initialize them to signaling NaN.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use nanMod
  use clmtype
  use clm_varpar, only: maxpatch_pft
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: initClmtype
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: init_pft_type
  private :: init_column_type
  private :: init_landunit_type
  private :: init_gridcell_type
  private :: init_energy_balance_type
  private :: init_water_balance_type
  private :: init_pft_pstate_type
  private :: init_pft_ecophys_constants
  private :: init_pft_DGVMecophys_constants
  private :: init_pft_estate_type
  private :: init_pft_wstate_type
  private :: init_pft_pdgvstate_type
  private :: init_pft_eflux_type
  private :: init_pft_mflux_type
  private :: init_pft_wflux_type
  private :: init_pft_cflux_type
  private :: init_pft_vflux_type
  private :: init_pft_dflux_type
  private :: init_column_pstate_type
  private :: init_column_estate_type
  private :: init_column_wstate_type
  private :: init_column_cstate_type
  private :: init_column_eflux_type
  private :: init_column_wflux_type
  private :: init_landunit_pstate_type
  private :: init_gridcell_pstate_type
  private :: init_atm2lnd_state_type
  private :: init_lnd2atm_state_type
  private :: init_atm2lnd_flux_type
  private :: init_lnd2atm_flux_type
  private :: init_gridcell_wflux_type
!EOP
!----------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initClmtype
!
! !INTERFACE:
  subroutine initClmtype()
!
! !DESCRIPTION:
! Initialize clmtype components to signaling nan
! The following clmtype components should NOT be initialized here
! since they are set in routine clm_map which is called before this
! routine is invoked
!    *%area, *%wt, *%wtlnd, *%wtxy, *%ixy, *%jxy, *%mxy, %snindex
!    *%ifspecial, *%ityplun, *%itype
!    *%pfti, *%pftf, *%pftn
!    *%coli, *%colf, *%coln
!    *%luni, *%lunf, *%lunn
!
! !USES:
    use decompMod, only : get_proc_bounds, get_proc_global
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
!------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! ***NOTE: all topological info has global extent***

    call init_pft_type     (1, nump, clm3%g%l%c%p)
    call init_column_type  (1, numc, clm3%g%l%c)
    call init_landunit_type(1, numl, clm3%g%l)
    call init_gridcell_type(1, numg, clm3%g)

    ! pft ecophysiological constants

    call init_pft_ecophys_constants()

    ! pft DGVM-specific ecophysiological constants

    call init_pft_DGVMecophys_constants()

    ! energy balance structures (all levels)

    call init_energy_balance_type(begp, endp, clm3%g%l%c%p%pebal)
    call init_energy_balance_type(begc, endc, clm3%g%l%c%cebal)
    call init_energy_balance_type(begl, endl, clm3%g%l%lebal)
    call init_energy_balance_type(begg, endg, clm3%g%gebal)
    call init_energy_balance_type(1,       1, clm3%mebal)

    ! water balance structures (all levels)

    call init_water_balance_type(begp, endp, clm3%g%l%c%p%pwbal)
    call init_water_balance_type(begc, endc, clm3%g%l%c%cwbal)
    call init_water_balance_type(begl, endl, clm3%g%l%lwbal)
    call init_water_balance_type(begg, endg, clm3%g%gwbal)
    call init_water_balance_type(1,       1, clm3%mwbal)

    ! pft physical state variables at pft level and averaged to the column

    call init_pft_pstate_type(begp, endp, clm3%g%l%c%p%pps)
    call init_pft_pstate_type(begc, endc, clm3%g%l%c%cps%pps_a)

    ! pft energy state variables at the pft level and averaged to the column

    call init_pft_estate_type(begp, endp, clm3%g%l%c%p%pes)
    call init_pft_estate_type(begc, endc, clm3%g%l%c%ces%pes_a)

    ! pft water state variables at the pft level and averaged to the column

    call init_pft_wstate_type(begp, endp, clm3%g%l%c%p%pws)
    call init_pft_wstate_type(begc, endc, clm3%g%l%c%cws%pws_a)

    ! pft DGVM state variables at pft level and averaged to column

    call init_pft_pdgvstate_type(begp, endp, clm3%g%l%c%p%pdgvs)
    call init_pft_pdgvstate_type(begc, endc, clm3%g%l%c%cdgvs%pdgvs_a)

    ! pft energy flux variables at pft level and averaged to column

    call init_pft_eflux_type(begp, endp, clm3%g%l%c%p%pef)
    call init_pft_eflux_type(begc, endc, clm3%g%l%c%cef%pef_a)

    ! pft momentum flux variables at pft level and averaged to the column

    call init_pft_mflux_type(begp, endp, clm3%g%l%c%p%pmf)
    call init_pft_mflux_type(begc, endc, clm3%g%l%c%cmf%pmf_a)

    ! pft water flux variables

    call init_pft_wflux_type(begp, endp, clm3%g%l%c%p%pwf)
    call init_pft_wflux_type(begc, endc, clm3%g%l%c%cwf%pwf_a)

    ! pft carbon flux variables at pft level and averaged to column

    call init_pft_cflux_type(begp, endp, clm3%g%l%c%p%pcf)
    call init_pft_cflux_type(begc, endc, clm3%g%l%c%ccf%pcf_a)

    ! pft VOC flux variables at pft level and averaged to column

    call init_pft_vflux_type(begp, endp, clm3%g%l%c%p%pvf)
    call init_pft_vflux_type(begc, endc, clm3%g%l%c%cvf%pvf_a)

    ! pft dust flux variables at pft level and averaged to column

    call init_pft_dflux_type(begp, endp, clm3%g%l%c%p%pdf)
    call init_pft_dflux_type(begc, endc, clm3%g%l%c%cdf%pdf_a)

    ! column physical state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_pstate_type(begc, endc, clm3%g%l%c%cps)
    call init_column_pstate_type(begl, endl, clm3%g%l%lps%cps_a)
    call init_column_pstate_type(begg, endg, clm3%g%gps%cps_a)
    call init_column_pstate_type(1,       1, clm3%mps%cps_a)

    ! column energy state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_estate_type(begc, endc, clm3%g%l%c%ces)
    call init_column_estate_type(begl, endl, clm3%g%l%les%ces_a)
    call init_column_estate_type(begg, endg, clm3%g%ges%ces_a)
    call init_column_estate_type(1,       1, clm3%mes%ces_a)

    ! column water state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_wstate_type(begc, endc, clm3%g%l%c%cws)
    call init_column_wstate_type(begl, endl, clm3%g%l%lws%cws_a)
    call init_column_wstate_type(begg, endg, clm3%g%gws%cws_a)
    call init_column_wstate_type(1,       1, clm3%mws%cws_a)

    ! column carbon state variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_cstate_type(begc, endc, clm3%g%l%c%ccs)
    call init_column_cstate_type(begl, endl, clm3%g%l%lcs%ccs_a)
    call init_column_cstate_type(begg, endg, clm3%g%gcs%ccs_a)
    call init_column_cstate_type(1,       1, clm3%mcs%ccs_a)

    ! column energy flux variables at column level and averaged to
    ! the landunit and gridcell and model

    call init_column_eflux_type(begc, endc, clm3%g%l%c%cef)
    call init_column_eflux_type(begl, endl, clm3%g%l%lef%cef_a)
    call init_column_eflux_type(begg, endg, clm3%g%gef%cef_a)
    call init_column_eflux_type(1,       1, clm3%mef%cef_a)

    ! column water flux variables at column level and averaged to
    ! landunit, gridcell and model level

    call init_column_wflux_type(begc, endc, clm3%g%l%c%cwf)
    call init_column_wflux_type(begl, endl, clm3%g%l%lwf%cwf_a)
    call init_column_wflux_type(begg, endg, clm3%g%gwf%cwf_a)
    call init_column_wflux_type(1,       1, clm3%mwf%cwf_a)

    ! land unit physical state variables

    call init_landunit_pstate_type(begl, endl, clm3%g%l%lps)

    ! gridcell DGVM variables
    call init_gridcell_dgvstate_type(begg, endg, clm3%g%gdgvs)

    ! gridcell physical state variables

    call init_gridcell_pstate_type(begg, endg, clm3%g%gps)

    ! gridcell: water flux variables

    call init_gridcell_wflux_type(begg, endg, clm3%g%gwf)

    ! gridcell atmosphere->land state and flux variables

    call init_atm2lnd_state_type(begg, endg, clm3%g%a2ls)
    call init_atm2lnd_flux_type (begg, endg, clm3%g%a2lf)

    ! gridcell land->atmosphere state and flux variables

    call init_lnd2atm_state_type(begg, endg, clm3%g%l2as)
    call init_lnd2atm_flux_type (begg, endg, clm3%g%l2af)

  end subroutine initClmtype

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_type
!
! !INTERFACE:
  subroutine init_pft_type (beg, end, p)
!
! !DESCRIPTION:
! Initialize components of pft_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(pft_type), intent(inout):: p
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(p%column(beg:end))
    allocate(p%landunit(beg:end))
    allocate(p%gridcell(beg:end))
    allocate(p%itype(beg:end))
    allocate(p%area(beg:end))
    allocate(p%wtcol(beg:end))
    allocate(p%wtlunit(beg:end))
    allocate(p%wtgcell(beg:end))
    allocate(p%latdeg(beg:end))
    allocate(p%londeg(beg:end))
    allocate(p%ixy(beg:end))
    allocate(p%jxy(beg:end))
    allocate(p%mxy(beg:end))
    allocate(p%snindex(beg:end))

  end subroutine init_pft_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_type
!
! !INTERFACE:
  subroutine init_column_type (beg, end, c)
!
! !DESCRIPTION:
! Initialize components of column_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(column_type), intent(inout):: c
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(c%pfti(beg:end))
   allocate(c%pftf(beg:end))
   allocate(c%npfts(beg:end))
   allocate(c%landunit(beg:end))
   allocate(c%gridcell(beg:end))
   allocate(c%itype(beg:end))
   allocate(c%area(beg:end))
   allocate(c%wtgcell(beg:end))
   allocate(c%wtlunit(beg:end))
   allocate(c%ixy(beg:end))
   allocate(c%jxy(beg:end))
   allocate(c%latdeg(beg:end))
   allocate(c%londeg(beg:end))
   allocate(c%snindex(beg:end))

  end subroutine init_column_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_type
!
! !INTERFACE:
  subroutine init_landunit_type (beg, end,l)
!
! !DESCRIPTION:
! Initialize components of landunit_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(landunit_type), intent(inout):: l
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(l%coli(beg:end))
   allocate(l%colf(beg:end))
   allocate(l%ncolumns(beg:end))
   allocate(l%pfti(beg:end))
   allocate(l%pftf(beg:end))
   allocate(l%npfts(beg:end))
   allocate(l%gridcell(beg:end))
   allocate(l%itype(beg:end))
   allocate(l%area(beg:end))
   allocate(l%wtgcell(beg:end))
   allocate(l%ixy(beg:end))
   allocate(l%jxy(beg:end))
   allocate(l%latdeg(beg:end))
   allocate(l%londeg(beg:end))
   allocate(l%ifspecial(beg:end))
   allocate(l%lakpoi(beg:end))
   allocate(l%snindex(beg:end))

  end subroutine init_landunit_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_type
!
! !INTERFACE:
  subroutine init_gridcell_type (beg, end,g)
!
! !DESCRIPTION:
! Initialize components of gridcell_type structure
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(gridcell_type), intent(inout):: g
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

   allocate(g%luni(beg:end))
   allocate(g%lunf(beg:end))
   allocate(g%nlandunits(beg:end))
   allocate(g%coli(beg:end))
   allocate(g%colf(beg:end))
   allocate(g%ncolumns(beg:end))
   allocate(g%pfti(beg:end))
   allocate(g%pftf(beg:end))
   allocate(g%npfts(beg:end))
   allocate(g%itype(beg:end))
   allocate(g%area(beg:end))
   allocate(g%wtglob(beg:end))
   allocate(g%ixy(beg:end))
   allocate(g%jxy(beg:end))
   allocate(g%lat(beg:end))
   allocate(g%lon(beg:end))
   allocate(g%latdeg(beg:end))
   allocate(g%londeg(beg:end))
   allocate(g%landfrac(beg:end))
   allocate(g%snindex(beg:end))

  end subroutine init_gridcell_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_energy_balance_type
!
! !INTERFACE:
  subroutine init_energy_balance_type(beg, end, ebal)
!
! !DESCRIPTION:
! Initialize energy balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(energy_balance_type), intent(inout):: ebal
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ebal%errsoi(beg:end))
    allocate(ebal%errseb(beg:end))
    allocate(ebal%errsol(beg:end))
    allocate(ebal%errlon(beg:end))

    ebal%errsoi(beg:end) = nan
    ebal%errseb(beg:end) = nan
    ebal%errsol(beg:end) = nan
    ebal%errlon(beg:end) = nan

  end subroutine init_energy_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_water_balance_type
!
! !INTERFACE:
  subroutine init_water_balance_type(beg, end, wbal)
!
! !DESCRIPTION:
! Initialize water balance variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type(water_balance_type), intent(inout):: wbal
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(wbal%begwb(beg:end))
    allocate(wbal%endwb(beg:end))
    allocate(wbal%errh2o(beg:end))

    wbal%begwb(beg:end) = nan
    wbal%endwb(beg:end) = nan
    wbal%errh2o(beg:end) = nan

  end subroutine init_water_balance_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pstate_type
!
! !INTERFACE:
  subroutine init_pft_pstate_type(beg, end, pps)
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_pstate_type), intent(inout):: pps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pps%frac_veg_nosno(beg:end))
    allocate(pps%frac_veg_nosno_alb(beg:end))
    allocate(pps%rsw(beg:end))
    allocate(pps%emv(beg:end))
    allocate(pps%z0mv(beg:end))
    allocate(pps%z0hv(beg:end))
    allocate(pps%z0qv(beg:end))
    allocate(pps%rootfr(beg:end,1:nlevsoi))
    allocate(pps%rootr(beg:end,1:nlevsoi))
    allocate(pps%dewmx(beg:end))
    allocate(pps%rssun(beg:end))
    allocate(pps%rssha(beg:end))
    allocate(pps%laisun(beg:end))
    allocate(pps%laisha(beg:end))
    allocate(pps%btran(beg:end))
    allocate(pps%fsun(beg:end))
    allocate(pps%tlai(beg:end))
    allocate(pps%tsai(beg:end))
    allocate(pps%elai(beg:end))
    allocate(pps%esai(beg:end))
    allocate(pps%igs(beg:end))
    allocate(pps%stembio(beg:end))
    allocate(pps%rootbio(beg:end))
    allocate(pps%fwet(beg:end))
    allocate(pps%fdry(beg:end))
    allocate(pps%dt_veg(beg:end))
    allocate(pps%htop(beg:end))
    allocate(pps%hbot(beg:end))
    allocate(pps%z0m(beg:end))
    allocate(pps%displa(beg:end))
    allocate(pps%albd(beg:end,1:numrad))
    allocate(pps%albi(beg:end,1:numrad))
    allocate(pps%fabd(beg:end,1:numrad))
    allocate(pps%fabi(beg:end,1:numrad))
    allocate(pps%ftdd(beg:end,1:numrad))
    allocate(pps%ftid(beg:end,1:numrad))
    allocate(pps%ftii(beg:end,1:numrad))
    allocate(pps%u10(beg:end))
    allocate(pps%fv(beg:end))
    allocate(pps%ram1(beg:end))

    pps%frac_veg_nosno(beg:end) = bigint
    pps%frac_veg_nosno_alb(beg:end) = bigint
    pps%rsw(beg:end) = nan
    pps%emv(beg:end) = nan
    pps%z0mv(beg:end) = nan
    pps%z0hv(beg:end) = nan
    pps%z0qv(beg:end) = nan
    pps%rootfr(beg:end,:nlevsoi) = nan
    pps%rootr (beg:end,:nlevsoi) = nan
    pps%dewmx(beg:end) = nan
    pps%rssun(beg:end) = nan
    pps%rssha(beg:end) = nan
    pps%laisun(beg:end) = nan
    pps%laisha(beg:end) = nan
    pps%btran(beg:end) = nan
    pps%fsun(beg:end) = nan
    pps%tlai(beg:end) = nan
    pps%tsai(beg:end) = nan
    pps%elai(beg:end) = nan
    pps%esai(beg:end) = nan
    pps%igs(beg:end) = nan
    pps%stembio(beg:end) = nan
    pps%rootbio(beg:end) = nan
    pps%fwet(beg:end) = nan
    pps%fdry(beg:end) = nan
    pps%dt_veg(beg:end) = nan
    pps%htop(beg:end) = nan
    pps%hbot(beg:end) = nan
    pps%z0m(beg:end) = nan
    pps%displa(beg:end) = nan
    pps%albd(beg:end,:numrad) = nan
    pps%albi(beg:end,:numrad) = nan
    pps%fabd(beg:end,:numrad) = nan
    pps%fabi(beg:end,:numrad) = nan
    pps%ftdd(beg:end,:numrad) = nan
    pps%ftid(beg:end,:numrad) = nan
    pps%ftii(beg:end,:numrad) = nan
    pps%u10(beg:end) = nan
    pps%fv(beg:end) = nan
    pps%ram1(beg:end) = nan
  end subroutine init_pft_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_ecophys_constants
!
! !INTERFACE:
  subroutine init_pft_ecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pftcon%ncorn(0:numpft))
    allocate(pftcon%nwheat(0:numpft))
    allocate(pftcon%noveg(0:numpft))
    allocate(pftcon%ntree(0:numpft))
    allocate(pftcon%smpmax(0:numpft))
    allocate(pftcon%foln(0:numpft))
    allocate(pftcon%dleaf(0:numpft))
    allocate(pftcon%c3psn(0:numpft))
    allocate(pftcon%vcmx25(0:numpft))
    allocate(pftcon%mp(0:numpft))
    allocate(pftcon%qe25(0:numpft))
    allocate(pftcon%xl(0:numpft))
    allocate(pftcon%rhol(0:numpft,numrad))
    allocate(pftcon%rhos(0:numpft,numrad))
    allocate(pftcon%taul(0:numpft,numrad))
    allocate(pftcon%taus(0:numpft,numrad))
    allocate(pftcon%z0mr(0:numpft))
    allocate(pftcon%displar(0:numpft))
    allocate(pftcon%roota_par(0:numpft))
    allocate(pftcon%rootb_par(0:numpft))
    allocate(pftcon%sla(0:numpft))

    pftcon%ncorn(:) = bigint
    pftcon%nwheat(:) = bigint
    pftcon%noveg(:) = bigint
    pftcon%ntree(:) = bigint
    pftcon%smpmax(:) = -1.5e5
    pftcon%foln(:) = nan
    pftcon%dleaf(:) = nan
    pftcon%c3psn(:) = nan
    pftcon%vcmx25(:) = nan
    pftcon%mp(:) = nan
    pftcon%qe25(:) = nan
    pftcon%xl(:) = nan
    pftcon%rhol(:,:numrad) = nan
    pftcon%rhos(:,:numrad) = nan
    pftcon%taul(:,:numrad) = nan
    pftcon%taus(:,:numrad) = nan
    pftcon%z0mr(:) = nan
    pftcon%displar(:) = nan
    pftcon%roota_par(:) = nan
    pftcon%rootb_par(:) = nan
    pftcon%sla(:) = nan

  end subroutine init_pft_ecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_DGVMecophys_constants
!
! !INTERFACE:
  subroutine init_pft_DGVMecophys_constants()
!
! !DESCRIPTION:
! Initialize pft physical state
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(dgv_pftcon%respcoeff(0:numpft))
    allocate(dgv_pftcon%flam(0:numpft))
    allocate(dgv_pftcon%resist(0:numpft))
    allocate(dgv_pftcon%l_turn(0:numpft))
    allocate(dgv_pftcon%l_long(0:numpft))
    allocate(dgv_pftcon%s_turn(0:numpft))
    allocate(dgv_pftcon%r_turn(0:numpft))
    allocate(dgv_pftcon%l_cton(0:numpft))
    allocate(dgv_pftcon%s_cton(0:numpft))
    allocate(dgv_pftcon%r_cton(0:numpft))
    allocate(dgv_pftcon%l_morph(0:numpft))
    allocate(dgv_pftcon%l_phen(0:numpft))
    allocate(dgv_pftcon%lmtorm(0:numpft))
    allocate(dgv_pftcon%crownarea_max(0:numpft))
    allocate(dgv_pftcon%init_lai(0:numpft))
    allocate(dgv_pftcon%x(0:numpft))
    allocate(dgv_pftcon%tcmin(0:numpft))
    allocate(dgv_pftcon%tcmax(0:numpft))
    allocate(dgv_pftcon%gddmin(0:numpft))
    allocate(dgv_pftcon%twmax(0:numpft))
    allocate(dgv_pftcon%lm_sapl(0:numpft))
    allocate(dgv_pftcon%sm_sapl(0:numpft))
    allocate(dgv_pftcon%hm_sapl(0:numpft))
    allocate(dgv_pftcon%rm_sapl(0:numpft))
    allocate(dgv_pftcon%tree(0:numpft))
    allocate(dgv_pftcon%summergreen(0:numpft))
    allocate(dgv_pftcon%raingreen(0:numpft))
    allocate(dgv_pftcon%reinickerp(0:numpft))
    allocate(dgv_pftcon%wooddens(0:numpft))
    allocate(dgv_pftcon%latosa(0:numpft))
    allocate(dgv_pftcon%allom1(0:numpft))
    allocate(dgv_pftcon%allom2(0:numpft))
    allocate(dgv_pftcon%allom3(0:numpft))

    dgv_pftcon%respcoeff(:) = nan
    dgv_pftcon%flam(:) = nan
    dgv_pftcon%resist(:) = nan
    dgv_pftcon%l_turn(:) = nan
    dgv_pftcon%l_long(:) = nan
    dgv_pftcon%s_turn(:) = nan
    dgv_pftcon%r_turn(:) = nan
    dgv_pftcon%l_cton(:) = nan
    dgv_pftcon%s_cton(:) = nan
    dgv_pftcon%r_cton(:) = nan
    dgv_pftcon%l_morph(:) = nan
    dgv_pftcon%l_phen(:) = nan
    dgv_pftcon%lmtorm(:) = nan
    dgv_pftcon%crownarea_max(:) = nan
    dgv_pftcon%init_lai(:) = nan
    dgv_pftcon%x(:) = nan
    dgv_pftcon%tcmin(:) = nan
    dgv_pftcon%tcmax(:) = nan
    dgv_pftcon%gddmin(:) = nan
    dgv_pftcon%twmax(:) = nan
    dgv_pftcon%lm_sapl(:) = nan
    dgv_pftcon%sm_sapl(:) = nan
    dgv_pftcon%hm_sapl(:) = nan
    dgv_pftcon%rm_sapl(:) = nan
    dgv_pftcon%tree(:) = .false.
    dgv_pftcon%summergreen(:) = .false.
    dgv_pftcon%raingreen(:) = .false.
    dgv_pftcon%reinickerp(:) = nan
    dgv_pftcon%wooddens(:) = nan
    dgv_pftcon%latosa(:) = nan
    dgv_pftcon%allom1(:) = nan
    dgv_pftcon%allom2(:) = nan
    dgv_pftcon%allom3(:) = nan

  end subroutine init_pft_DGVMecophys_constants

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_energy_type
!
! !INTERFACE:
  subroutine init_pft_estate_type(beg, end, pes)
!
! !DESCRIPTION:
! Initialize pft energy state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_estate_type), intent(inout):: pes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pes%t_ref2m(beg:end))
    allocate(pes%t_ref2m_min(beg:end))
    allocate(pes%t_ref2m_max(beg:end))
    allocate(pes%t_ref2m_min_inst(beg:end))
    allocate(pes%t_ref2m_max_inst(beg:end))
    allocate(pes%q_ref2m(beg:end))
    allocate(pes%t_veg(beg:end))

    pes%t_ref2m(beg:end) = nan
    pes%t_ref2m_min(beg:end) = nan
    pes%t_ref2m_max(beg:end) = nan
    pes%t_ref2m_min_inst(beg:end) = nan
    pes%t_ref2m_max_inst(beg:end) = nan
    pes%q_ref2m(beg:end) = nan
    pes%t_veg(beg:end) = nan

  end subroutine init_pft_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wstate_type
!
! !INTERFACE:
  subroutine init_pft_wstate_type(beg, end, pws)
!
! !DESCRIPTION:
! Initialize pft water state
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wstate_type), intent(inout):: pws !pft water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pws%h2ocan(beg:end))
    pws%h2ocan(beg:end) = nan

  end subroutine init_pft_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_pdgvstate_type
!
! !INTERFACE:
  subroutine init_pft_pdgvstate_type(beg, end, pdgvs)
!
! !DESCRIPTION:
! Initialize pft DGVM state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dgvstate_type), intent(inout):: pdgvs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdgvs%agdd0(beg:end))
    allocate(pdgvs%agdd5(beg:end))
    allocate(pdgvs%agddtw(beg:end))
    allocate(pdgvs%agdd(beg:end))
    allocate(pdgvs%t10(beg:end))
    allocate(pdgvs%t_mo(beg:end))
    allocate(pdgvs%t_mo_min(beg:end))
    allocate(pdgvs%fnpsn10(beg:end))
    allocate(pdgvs%prec365(beg:end))
    allocate(pdgvs%agdd20(beg:end))
    allocate(pdgvs%tmomin20(beg:end))
    allocate(pdgvs%t10min(beg:end))
    allocate(pdgvs%tsoi25(beg:end))
    allocate(pdgvs%annpsn(beg:end))
    allocate(pdgvs%annpsnpot(beg:end))
    allocate(pdgvs%present(beg:end))
    allocate(pdgvs%dphen(beg:end))
    allocate(pdgvs%leafon(beg:end))
    allocate(pdgvs%leafof(beg:end))
    allocate(pdgvs%nind(beg:end))
    allocate(pdgvs%lm_ind(beg:end))
    allocate(pdgvs%sm_ind(beg:end))
    allocate(pdgvs%hm_ind(beg:end))
    allocate(pdgvs%rm_ind(beg:end))
    allocate(pdgvs%lai_ind(beg:end))
    allocate(pdgvs%fpcinc(beg:end))
    allocate(pdgvs%fpcgrid(beg:end))
    allocate(pdgvs%crownarea(beg:end))
    allocate(pdgvs%bm_inc(beg:end))
    allocate(pdgvs%afmicr(beg:end))
    allocate(pdgvs%firelength (beg:end))
    allocate(pdgvs%litterag(beg:end))
    allocate(pdgvs%litterbg(beg:end))
    allocate(pdgvs%cpool_fast(beg:end))
    allocate(pdgvs%cpool_slow(beg:end))
    allocate(pdgvs%k_fast_ave(beg:end))
    allocate(pdgvs%k_slow_ave(beg:end))
    allocate(pdgvs%litter_decom_ave(beg:end))
    allocate(pdgvs%turnover_ind(beg:end))

    pdgvs%agdd0(beg:end) = nan
    pdgvs%agdd5(beg:end) = nan
    pdgvs%agddtw(beg:end) = nan
    pdgvs%agdd(beg:end) = nan
    pdgvs%t10(beg:end) = nan
    pdgvs%t_mo(beg:end) = nan
    pdgvs%t_mo_min(beg:end) = nan
    pdgvs%fnpsn10(beg:end) = nan
    pdgvs%prec365(beg:end) = nan
    pdgvs%agdd20(beg:end) = nan
    pdgvs%tmomin20(beg:end) = nan
    pdgvs%t10min(beg:end) = nan
    pdgvs%tsoi25(beg:end) = nan
    pdgvs%annpsn(beg:end) = nan
    pdgvs%annpsnpot(beg:end) = nan
    pdgvs%present(beg:end) = .false.
    pdgvs%dphen(beg:end) = nan
    pdgvs%leafon(beg:end) = nan
    pdgvs%leafof(beg:end) = nan
    pdgvs%nind(beg:end) = nan
    pdgvs%lm_ind(beg:end) = nan
    pdgvs%sm_ind(beg:end) = nan
    pdgvs%hm_ind(beg:end) = nan
    pdgvs%rm_ind(beg:end) = nan
    pdgvs%lai_ind(beg:end) = nan
    pdgvs%fpcinc(beg:end) = nan
    pdgvs%fpcgrid(beg:end) = nan
    pdgvs%crownarea(beg:end) = nan
    pdgvs%bm_inc(beg:end) = nan
    pdgvs%afmicr(beg:end) = nan
    pdgvs%firelength (beg:end) = nan
    pdgvs%litterag(beg:end) = nan
    pdgvs%litterbg(beg:end) = nan
    pdgvs%cpool_fast(beg:end) = nan
    pdgvs%cpool_slow(beg:end) = nan
    pdgvs%k_fast_ave(beg:end) = nan
    pdgvs%k_slow_ave(beg:end) = nan
    pdgvs%litter_decom_ave(beg:end) = nan
    pdgvs%turnover_ind(beg:end) = nan

  end subroutine init_pft_pdgvstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_eflux_type
!
! !INTERFACE:
  subroutine init_pft_eflux_type(beg, end, pef)
!
! !DESCRIPTION:
! Initialize pft energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_eflux_type), intent(inout):: pef
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pef%sabg(beg:end))
    allocate(pef%sabv(beg:end))
    allocate(pef%fsa(beg:end))
    allocate(pef%fsr(beg:end))
    allocate(pef%parsun(beg:end))
    allocate(pef%parsha(beg:end))
    allocate(pef%dlrad(beg:end))
    allocate(pef%ulrad(beg:end))
    allocate(pef%eflx_lh_tot(beg:end))
    allocate(pef%eflx_lh_grnd(beg:end))
    allocate(pef%eflx_soil_grnd(beg:end))
    allocate(pef%eflx_sh_tot(beg:end))
    allocate(pef%eflx_sh_grnd(beg:end))
    allocate(pef%eflx_sh_veg(beg:end))
    allocate(pef%eflx_lh_vege(beg:end))
    allocate(pef%eflx_lh_vegt(beg:end))
    allocate(pef%cgrnd(beg:end))
    allocate(pef%cgrndl(beg:end))
    allocate(pef%cgrnds(beg:end))
    allocate(pef%eflx_gnet(beg:end))
    allocate(pef%dgnetdT(beg:end))
    allocate(pef%eflx_lwrad_out(beg:end))
    allocate(pef%eflx_lwrad_net(beg:end))
    allocate(pef%fsds_vis_d(beg:end))
    allocate(pef%fsds_nir_d(beg:end))
    allocate(pef%fsds_vis_i(beg:end))
    allocate(pef%fsds_nir_i(beg:end))
    allocate(pef%fsr_vis_d(beg:end))
    allocate(pef%fsr_nir_d(beg:end))
    allocate(pef%fsr_vis_i(beg:end))
    allocate(pef%fsr_nir_i(beg:end))
    allocate(pef%fsds_vis_d_ln(beg:end))
    allocate(pef%fsds_nir_d_ln(beg:end))
    allocate(pef%fsr_vis_d_ln(beg:end))
    allocate(pef%fsr_nir_d_ln(beg:end))

    pef%sabg(beg:end) = nan
    pef%sabv(beg:end) = nan
    pef%fsa(beg:end) = nan
    pef%fsr(beg:end) = nan
    pef%parsun(beg:end) = nan
    pef%parsha(beg:end) = nan
    pef%dlrad(beg:end) = nan
    pef%ulrad(beg:end) = nan
    pef%eflx_lh_tot(beg:end) = nan
    pef%eflx_lh_grnd(beg:end) = nan
    pef%eflx_soil_grnd(beg:end) = nan
    pef%eflx_sh_tot(beg:end) = nan
    pef%eflx_sh_grnd(beg:end) = nan
    pef%eflx_sh_veg(beg:end) = nan
    pef%eflx_lh_vege(beg:end) = nan
    pef%eflx_lh_vegt(beg:end) = nan
    pef%cgrnd(beg:end) = nan
    pef%cgrndl(beg:end) = nan
    pef%cgrnds(beg:end) = nan
    pef%eflx_gnet(beg:end) = nan
    pef%dgnetdT(beg:end) = nan
    pef%eflx_lwrad_out(beg:end) = nan
    pef%eflx_lwrad_net(beg:end) = nan
    pef%fsds_vis_d(beg:end) = nan
    pef%fsds_nir_d(beg:end) = nan
    pef%fsds_vis_i(beg:end) = nan
    pef%fsds_nir_i(beg:end) = nan
    pef%fsr_vis_d(beg:end) = nan
    pef%fsr_nir_d(beg:end) = nan
    pef%fsr_vis_i(beg:end) = nan
    pef%fsr_nir_i(beg:end) = nan
    pef%fsds_vis_d_ln(beg:end) = nan
    pef%fsds_nir_d_ln(beg:end) = nan
    pef%fsr_vis_d_ln(beg:end) = nan
    pef%fsr_nir_d_ln(beg:end) = nan

  end subroutine init_pft_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_mflux_type
!
! !INTERFACE:
  subroutine init_pft_mflux_type(beg, end, pmf)
!
! !DESCRIPTION:
! Initialize pft momentum flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_mflux_type), intent(inout) :: pmf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pmf%taux(beg:end))
    allocate(pmf%tauy(beg:end))

    pmf%taux(beg:end) = nan
    pmf%tauy(beg:end) = nan

  end subroutine init_pft_mflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_wflux_type
!
! !INTERFACE:
  subroutine init_pft_wflux_type(beg, end, pwf)
!
! !DESCRIPTION:
! Initialize pft water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_wflux_type), intent(inout) :: pwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pwf%qflx_prec_intr(beg:end))
    allocate(pwf%qflx_prec_grnd(beg:end))
    allocate(pwf%qflx_rain_grnd(beg:end))
    allocate(pwf%qflx_snow_grnd(beg:end))
    allocate(pwf%qflx_snowcap(beg:end))
    allocate(pwf%qflx_evap_veg(beg:end))
    allocate(pwf%qflx_tran_veg(beg:end))
    allocate(pwf%qflx_evap_can(beg:end))
    allocate(pwf%qflx_evap_soi(beg:end))
    allocate(pwf%qflx_evap_tot(beg:end))
    allocate(pwf%qflx_evap_grnd(beg:end))
    allocate(pwf%qflx_dew_grnd(beg:end))
    allocate(pwf%qflx_sub_snow(beg:end))
    allocate(pwf%qflx_dew_snow(beg:end))

    pwf%qflx_prec_intr(beg:end) = nan
    pwf%qflx_prec_grnd(beg:end) = nan
    pwf%qflx_rain_grnd(beg:end) = nan
    pwf%qflx_snow_grnd(beg:end) = nan
    pwf%qflx_snowcap(beg:end) = nan
    pwf%qflx_evap_veg(beg:end) = nan
    pwf%qflx_tran_veg(beg:end) = nan
    pwf%qflx_evap_can(beg:end) = nan
    pwf%qflx_evap_soi(beg:end) = nan
    pwf%qflx_evap_tot(beg:end) = nan
    pwf%qflx_evap_grnd(beg:end) = nan
    pwf%qflx_dew_grnd(beg:end) = nan
    pwf%qflx_sub_snow(beg:end) = nan
    pwf%qflx_dew_snow(beg:end) = nan

  end subroutine init_pft_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_cflux_type
!
! !INTERFACE:
  subroutine init_pft_cflux_type(beg, end, pcf)
!
! !DESCRIPTION:
! Initialize pft carbon flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_cflux_type), intent(inout) :: pcf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pcf%psnsun(beg:end))
    allocate(pcf%psnsha(beg:end))
    allocate(pcf%fpsn(beg:end))
    allocate(pcf%frm(beg:end))
    allocate(pcf%frmf(beg:end))
    allocate(pcf%frms(beg:end))
    allocate(pcf%frmr(beg:end))
    allocate(pcf%frg(beg:end))
    allocate(pcf%dmi(beg:end))
    allocate(pcf%fco2(beg:end))
    allocate(pcf%fmicr(beg:end))

    pcf%psnsun(beg:end) = nan
    pcf%psnsha(beg:end) = nan
    pcf%fpsn(beg:end) = nan
    pcf%frm(beg:end) = nan
    pcf%frmf(beg:end) = nan
    pcf%frms(beg:end) = nan
    pcf%frmr(beg:end) = nan
    pcf%frg(beg:end) = nan
    pcf%dmi(beg:end) = nan
    pcf%fco2(beg:end) = nan
    pcf%fmicr(beg:end) = nan

  end subroutine init_pft_cflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_vflux_type
!
! !INTERFACE:
  subroutine init_pft_vflux_type(beg, end, pvf)
!
! !DESCRIPTION:
! Initialize pft VOC flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_vflux_type), intent(inout) :: pvf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pvf%vocflx_tot(beg:end))
    allocate(pvf%vocflx(beg:end,1:nvoc))
    allocate(pvf%vocflx_1(beg:end))
    allocate(pvf%vocflx_2(beg:end))
    allocate(pvf%vocflx_3(beg:end))
    allocate(pvf%vocflx_4(beg:end))
    allocate(pvf%vocflx_5(beg:end))

    pvf%vocflx_tot(beg:end) = nan
    pvf%vocflx(beg:end,1:nvoc) = nan
    pvf%vocflx_1(beg:end) = nan
    pvf%vocflx_2(beg:end) = nan
    pvf%vocflx_3(beg:end) = nan
    pvf%vocflx_4(beg:end) = nan
    pvf%vocflx_5(beg:end) = nan

  end subroutine init_pft_vflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_pft_dflux_type
!
! !INTERFACE:
  subroutine init_pft_dflux_type(beg, end, pdf)
!
! !DESCRIPTION:
! Initialize pft dust flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (pft_dflux_type), intent(inout):: pdf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(pdf%flx_mss_vrt_dst(beg:end,1:ndst))
    allocate(pdf%flx_mss_vrt_dst_tot(beg:end))
    allocate(pdf%vlc_trb(beg:end,1:ndst))
    allocate(pdf%vlc_trb_1(beg:end))
    allocate(pdf%vlc_trb_2(beg:end))
    allocate(pdf%vlc_trb_3(beg:end))
    allocate(pdf%vlc_trb_4(beg:end))

    pdf%flx_mss_vrt_dst(beg:end,1:ndst) = nan
    pdf%flx_mss_vrt_dst_tot(beg:end) = nan
    pdf%vlc_trb(beg:end,1:ndst) = nan
    pdf%vlc_trb_1(beg:end) = nan
    pdf%vlc_trb_2(beg:end) = nan
    pdf%vlc_trb_3(beg:end) = nan
    pdf%vlc_trb_4(beg:end) = nan

  end subroutine init_pft_dflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_pstate_type
!
! !INTERFACE:
  subroutine init_column_pstate_type(beg, end, cps)
!
! !DESCRIPTION:
! Initialize column physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_pstate_type), intent(inout):: cps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cps%snl(beg:end))      !* cannot be averaged up
    allocate(cps%isoicol(beg:end))  !* cannot be averaged up
    allocate(cps%bsw(beg:end,nlevsoi))
    allocate(cps%watsat(beg:end,nlevsoi))
    allocate(cps%hksat(beg:end,nlevsoi))
    allocate(cps%sucsat(beg:end,nlevsoi))
    allocate(cps%csol(beg:end,nlevsoi))
    allocate(cps%tkmg(beg:end,nlevsoi))
    allocate(cps%tkdry(beg:end,nlevsoi))
    allocate(cps%tksatu(beg:end,nlevsoi))
    allocate(cps%smpmin(beg:end))
    allocate(cps%gwc_thr(beg:end))
    allocate(cps%mss_frc_cly_vld(beg:end))
    allocate(cps%mbl_bsn_fct(beg:end))
    allocate(cps%do_capsnow(beg:end))
    allocate(cps%snowdp(beg:end))
    allocate(cps%snowage(beg:end))
    allocate(cps%frac_sno (beg:end))
    allocate(cps%zi(beg:end,-nlevsno+0:nlevsoi))
    allocate(cps%dz(beg:end,-nlevsno+1:nlevsoi))
    allocate(cps%z (beg:end,-nlevsno+1:nlevsoi))
    allocate(cps%frac_iceold(beg:end,-nlevsno+1:nlevsoi))
    allocate(cps%imelt(beg:end,-nlevsno+1:nlevsoi))
    allocate(cps%eff_porosity(beg:end,nlevsoi))
    allocate(cps%sfact(beg:end))
    allocate(cps%sfactmax(beg:end))
    allocate(cps%emg(beg:end))
    allocate(cps%z0mg(beg:end))
    allocate(cps%z0hg(beg:end))
    allocate(cps%z0qg(beg:end))
    allocate(cps%htvp(beg:end))
    allocate(cps%beta(beg:end))
    allocate(cps%zii(beg:end))
    allocate(cps%albgrd(beg:end,numrad))
    allocate(cps%albgri(beg:end,numrad))
    allocate(cps%rootr_column(beg:end,nlevsoi))
    allocate(cps%wf(beg:end))

    cps%isoicol(beg:end) = bigint
    cps%bsw(beg:end,1:nlevsoi) = nan
    cps%watsat(beg:end,1:nlevsoi) = nan
    cps%hksat(beg:end,1:nlevsoi) = nan
    cps%sucsat(beg:end,1:nlevsoi) = nan
    cps%csol(beg:end,1:nlevsoi) = nan
    cps%tkmg(beg:end,1:nlevsoi) = nan
    cps%tkdry(beg:end,1:nlevsoi) = nan
    cps%tksatu(beg:end,1:nlevsoi) = nan
    cps%smpmin(beg:end) = nan
    cps%gwc_thr(beg:end) = nan
    cps%mss_frc_cly_vld(beg:end) = nan
    cps%mbl_bsn_fct(beg:end) = nan
    cps%do_capsnow (beg:end)= .false.
    cps%snowdp(beg:end) = nan
    cps%snowage(beg:end) = nan
    cps%frac_sno(beg:end) = nan
    cps%zi(beg:end,-nlevsno+0:nlevsoi) = nan
    cps%dz(beg:end,-nlevsno+1:nlevsoi) = nan
    cps%z (beg:end,-nlevsno+1:nlevsoi) = nan
    cps%frac_iceold(beg:end,-nlevsno+1:nlevsoi) = nan
    cps%imelt(beg:end,-nlevsno+1:nlevsoi) = bigint
    cps%eff_porosity(beg:end,1:nlevsoi) = nan
    cps%sfact(beg:end) = nan
    cps%sfactmax(beg:end) = nan
    cps%emg(beg:end) = nan
    cps%z0mg(beg:end) = nan
    cps%z0hg(beg:end) = nan
    cps%z0qg(beg:end) = nan
    cps%htvp(beg:end) = nan
    cps%beta(beg:end) = nan
    cps%zii(beg:end) = nan
    cps%albgrd(beg:end,:numrad) = nan
    cps%albgri(beg:end,:numrad) = nan
    cps%rootr_column(beg:end,1:nlevsoi) = nan
    cps%wf(beg:end) = nan

  end subroutine init_column_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_estate_type
!
! !INTERFACE:
  subroutine init_column_estate_type(beg, end, ces)
!
! !DESCRIPTION:
! Initialize column energy state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_estate_type), intent(inout):: ces
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    allocate(ces%t_grnd(beg:end))
    allocate(ces%dt_grnd(beg:end))
    allocate(ces%t_soisno(beg:end,-nlevsno+1:nlevsoi))
    allocate(ces%t_lake(beg:end,1:nlevlak))
    allocate(ces%tssbef(beg:end,-nlevsno+1:nlevsoi))
    allocate(ces%t_snow(beg:end))
    allocate(ces%thv(beg:end))
    allocate(ces%thm(beg:end))

    ces%t_grnd(beg:end) = nan
    ces%dt_grnd(beg:end) = nan
    ces%t_soisno(beg:end,-nlevsno+1:nlevsoi) = nan
    ces%t_lake(beg:end,1:nlevlak)= nan
    ces%tssbef(beg:end,-nlevsno+1:nlevsoi) = nan
    ces%t_snow(beg:end) = nan
    ces%thv(beg:end) = nan
    ces%thm(beg:end) = nan

  end subroutine init_column_estate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wstate_type
!
! !INTERFACE:
  subroutine init_column_wstate_type(beg, end, cws)
!
! !DESCRIPTION:
! Initialize column water state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wstate_type), intent(inout):: cws !column water state
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cws%h2osno(beg:end))
    allocate(cws%h2osoi_liq(beg:end,-nlevsno+1:nlevsoi))
    allocate(cws%h2osoi_ice(beg:end,-nlevsno+1:nlevsoi))
    allocate(cws%h2osoi_vol(beg:end,1:nlevsoi))
    allocate(cws%h2osno_old(beg:end))
    allocate(cws%qg(beg:end))
    allocate(cws%dqgdT(beg:end))
    allocate(cws%snowice(beg:end))
    allocate(cws%snowliq(beg:end))

    cws%h2osno(beg:end) = nan
    cws%h2osoi_liq(beg:end,-nlevsno+1:nlevsoi)= nan
    cws%h2osoi_ice(beg:end,-nlevsno+1:nlevsoi) = nan
    cws%h2osoi_vol(beg:end,1:nlevsoi) = nan
    cws%h2osno_old(beg:end) = nan
    cws%qg(beg:end) = nan
    cws%dqgdT(beg:end) = nan
    cws%snowice(beg:end) = nan
    cws%snowliq(beg:end) = nan

  end subroutine init_column_wstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_cstate_type
!
! !INTERFACE:
  subroutine init_column_cstate_type(beg, end, ccs)
!
! !DESCRIPTION:
! Initialize column carbon state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_cstate_type), intent(inout):: ccs
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(ccs%soilc(beg:end))
    ccs%soilc(beg:end) = nan

  end subroutine init_column_cstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_eflux_type
!
! !INTERFACE:
  subroutine init_column_eflux_type(beg, end, cef)
!
! !DESCRIPTION:
! Initialize column energy flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_eflux_type), intent(inout):: cef
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cef%eflx_snomelt(beg:end))
    allocate(cef%eflx_impsoil(beg:end))

    cef%eflx_snomelt(beg:end) = nan
    cef%eflx_impsoil(beg:end) = nan

  end subroutine init_column_eflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_column_wflux_type
!
! !INTERFACE:
  subroutine init_column_wflux_type(beg, end, cwf)
!
! !DESCRIPTION:
! Initialize column water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (column_wflux_type), intent(inout):: cwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(cwf%qflx_infl(beg:end))
    allocate(cwf%qflx_surf(beg:end))
    allocate(cwf%qflx_drain(beg:end))
    allocate(cwf%qflx_top_soil(beg:end))
    allocate(cwf%qflx_snomelt(beg:end))
    allocate(cwf%qflx_qrgwl(beg:end))
    allocate(cwf%qmelt(beg:end))

    cwf%qflx_infl(beg:end) = nan
    cwf%qflx_surf(beg:end) = nan
    cwf%qflx_drain(beg:end) = nan
    cwf%qflx_top_soil(beg:end) = nan
    cwf%qflx_snomelt(beg:end) = nan
    cwf%qflx_qrgwl(beg:end) = nan
    cwf%qmelt(beg:end) = nan

  end subroutine init_column_wflux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_landunit_pstate_type
!
! !INTERFACE:
  subroutine init_landunit_pstate_type(beg, end, lps)
!
! !DESCRIPTION:
! Initialize landunit physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (landunit_pstate_type), intent(inout):: lps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

  ! currently nothing is here - just a place holder

  end subroutine init_landunit_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_dgvstate_type
!
! !INTERFACE:
  subroutine init_gridcell_dgvstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell DGVM variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_dgvstate_type), intent(inout):: gps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gps%afirefrac(beg:end))
    allocate(gps%acfluxfire(beg:end))
    allocate(gps%bmfm(beg:end,maxpatch_pft))
    allocate(gps%afmicr(beg:end,maxpatch_pft))
    allocate(gps%begwater(beg:end))
    allocate(gps%endwater(beg:end))
    allocate(gps%begenergy(beg:end))
    allocate(gps%endenergy(beg:end))
    gps%afirefrac(beg:end) = nan
    gps%acfluxfire(beg:end) = nan
    gps%bmfm(beg:end,1:maxpatch_pft) = nan
    gps%afmicr(beg:end,1:maxpatch_pft) = nan
    gps%begwater(beg:end) = nan
    gps%endwater(beg:end) = nan
    gps%begenergy(beg:end) = nan
    gps%endenergy(beg:end) = nan

  end subroutine init_gridcell_dgvstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_pstate_type
!
! !INTERFACE:
  subroutine init_gridcell_pstate_type(beg, end, gps)
!
! !DESCRIPTION:
! Initialize gridcell physical state variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_pstate_type), intent(inout):: gps
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(gps%wtfact(beg:end))
    gps%wtfact(beg:end) = nan

  end subroutine init_gridcell_pstate_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_atm2lnd_state_type
!
! !INTERFACE:
  subroutine init_atm2lnd_state_type(beg, end, a2ls)
!
! !DESCRIPTION:
! Initialize atmospheric state variables required by the land
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (atm2lnd_state_type), intent(inout):: a2ls
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

#if (defined OFFLINE)
    allocate(a2ls%flfall(beg:end))
#endif
    allocate(a2ls%forc_t(beg:end))
    allocate(a2ls%forc_u(beg:end))
    allocate(a2ls%forc_v(beg:end))
    allocate(a2ls%forc_wind(beg:end))
    allocate(a2ls%forc_q(beg:end))
    allocate(a2ls%forc_hgt(beg:end))
    allocate(a2ls%forc_hgt_u(beg:end))
    allocate(a2ls%forc_hgt_t(beg:end))
    allocate(a2ls%forc_hgt_q(beg:end))
    allocate(a2ls%forc_pbot(beg:end))
    allocate(a2ls%forc_th(beg:end))
    allocate(a2ls%forc_vp(beg:end))
    allocate(a2ls%forc_rho(beg:end))
    allocate(a2ls%forc_co2(beg:end))
    allocate(a2ls%forc_o2(beg:end))
    allocate(a2ls%forc_psrf(beg:end))

#if (defined OFFLINE)
    a2ls%flfall(beg:end) = nan
#endif
    a2ls%forc_t(beg:end) = nan
    a2ls%forc_u(beg:end) = nan
    a2ls%forc_v(beg:end) = nan
    a2ls%forc_wind(beg:end) = nan
    a2ls%forc_q(beg:end) = nan
    a2ls%forc_hgt(beg:end) = nan
    a2ls%forc_hgt_u(beg:end) = nan
    a2ls%forc_hgt_t(beg:end) = nan
    a2ls%forc_hgt_q(beg:end) = nan
    a2ls%forc_pbot(beg:end) = nan
    a2ls%forc_th(beg:end) = nan
    a2ls%forc_vp(beg:end) = nan
    a2ls%forc_rho(beg:end) = nan
    a2ls%forc_co2(beg:end) = nan
    a2ls%forc_o2(beg:end) = nan
    a2ls%forc_psrf(beg:end) = nan

  end subroutine init_atm2lnd_state_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2atm_state_type
!
! !INTERFACE:
  subroutine init_lnd2atm_state_type(beg, end, l2as)
!
! !DESCRIPTION:
! Initialize land state variables required by the atmosphere
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (lnd2atm_state_type), intent(inout):: l2as
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(l2as%t_rad(beg:end))
    allocate(l2as%t_ref2m(beg:end))
    allocate(l2as%q_ref2m(beg:end))
    allocate(l2as%h2osno(beg:end))
    allocate(l2as%albd(beg:end,1:numrad))
    allocate(l2as%albi(beg:end,1:numrad))

    l2as%t_rad(beg:end) = nan
    l2as%t_ref2m(beg:end) = nan
    l2as%q_ref2m(beg:end) = nan
    l2as%h2osno(beg:end) = nan
    l2as%albd(beg:end,1:numrad) = nan
    l2as%albi(beg:end,1:numrad) = nan

  end subroutine init_lnd2atm_state_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_atm2lnd_flux_type
!
! !INTERFACE:
  subroutine init_atm2lnd_flux_type(beg, end, a2lf)
!
! !DESCRIPTION:
! Initialize atmospheric fluxes required by the land
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (atm2lnd_flux_type), intent(inout):: a2lf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(a2lf%forc_lwrad(beg:end))
    allocate(a2lf%forc_solad(beg:end,numrad))
    allocate(a2lf%forc_solai(beg:end,numrad))
    allocate(a2lf%forc_solar(beg:end))
    allocate(a2lf%forc_rain(beg:end))
    allocate(a2lf%forc_snow(beg:end))

    a2lf%forc_lwrad(beg:end) = nan
    a2lf%forc_solad(beg:end,1:numrad) = nan
    a2lf%forc_solai(beg:end,1:numrad) = nan
    a2lf%forc_solar(beg:end) = nan
    a2lf%forc_rain(beg:end) = nan
    a2lf%forc_snow(beg:end) = nan

  end subroutine init_atm2lnd_flux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_lnd2atm_flux_type
!
! !INTERFACE:
  subroutine init_lnd2atm_flux_type(beg, end, l2af)
!
! !DESCRIPTION:
! Initialize land fluxes required by the atmosphere
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (lnd2atm_flux_type), intent(inout):: l2af
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------

    allocate(l2af%taux(beg:end))
    allocate(l2af%tauy(beg:end))
    allocate(l2af%eflx_lwrad_out(beg:end))
    allocate(l2af%eflx_sh_tot(beg:end))
    allocate(l2af%eflx_lh_tot(beg:end))
    allocate(l2af%qflx_evap_tot(beg:end))
    allocate(l2af%fsa(beg:end))

    l2af%taux(beg:end) = nan
    l2af%tauy(beg:end) = nan
    l2af%eflx_lwrad_out(beg:end) = nan
    l2af%eflx_sh_tot(beg:end) = nan
    l2af%eflx_lh_tot(beg:end) = nan
    l2af%qflx_evap_tot(beg:end) = nan
    l2af%fsa(beg:end) = nan

  end subroutine init_lnd2atm_flux_type

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_gridcell_wflux_type
!
! !INTERFACE:
  subroutine init_gridcell_wflux_type(beg, end, gwf)
!
! !DESCRIPTION:
! Initialize gridcell water flux variables
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: beg, end
    type (gridcell_wflux_type), intent(inout):: gwf
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!------------------------------------------------------------------------
    allocate(gwf%qchan2(beg:end))
    allocate(gwf%qchocn2(beg:end))

    gwf%qchan2(beg:end) = 0.
    gwf%qchocn2(beg:end) = 0.

  end subroutine init_gridcell_wflux_type

end module clmtypeInitMod


