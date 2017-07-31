#include <misc.h>
#include <preproc.h>

module BalanceCheckMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: BalanceCheckMod
!
! !DESCRIPTION:
! Water and energy balance check.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: BalanceCheck ! Water and energy balance check
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: BalanceCheck
!
! !INTERFACE:
  subroutine BalanceCheck(lbp, ubp, lbc, ubc)
!
! !DESCRIPTION:
! This subroutine accumulates the numerical truncation errors of the water
! and energy balance calculation. It is helpful to see the performance of
! the process of integration.
!
! The error for energy balance:
!
! error = abs(Net radiation - change of internal energy - Sensible heat
!             - Latent heat)
!
! The error should be less than 0.02 W/m$^2$ in each time integration
! interval.
!
! The error for water balance:
!
! error = abs(precipitation - change of water storage - evaporation - runoff)
!
! The error should be less than 0.001 mm in each time integration interval.
!
! !USES:
    use clmtype
    use subgridAveMod
    use time_manager , only: get_step_size, get_nstep
!
! !ARGUMENTS:
    implicit none
    integer :: lbp, ubp ! pft-index bounds
    integer :: lbc, ubc ! column-index bounds
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! 15 September 1999: Yongjiu Dai; Initial code
! 15 December 1999:  Paul Houser and Jon Radakovich; F90 Revision
! 10 November 2000: Mariana Vertenstein
! Migrated to new data structures by Mariana Vertenstein and
! Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to original implicit in arrays
    integer , pointer :: pgridcell(:)       ! pft's gridcell index
    integer , pointer :: cgridcell(:)       ! column's gridcell index
    real(r8), pointer :: forc_rain(:)       ! rain rate [mm/s]
    real(r8), pointer :: forc_snow(:)       ! snow rate [mm/s]
    real(r8), pointer :: forc_lwrad(:)      ! downward infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: endwb(:)           ! water mass end of the time step
    real(r8), pointer :: begwb(:)           ! water mass begining of the time step
    real(r8), pointer :: fsa(:)             ! solar radiation absorbed (total) (W/m**2)
    real(r8), pointer :: fsr(:)             ! solar radiation reflected (W/m**2)
    real(r8), pointer :: eflx_lwrad_out(:)  ! emitted infrared (longwave) radiation (W/m**2)
    real(r8), pointer :: eflx_lwrad_net(:)  ! net infrared (longwave) rad (W/m**2) [+ = to atm]
    real(r8), pointer :: sabv(:)            ! solar radiation absorbed by vegetation (W/m**2)
    real(r8), pointer :: sabg(:)            ! solar radiation absorbed by ground (W/m**2)
    real(r8), pointer :: eflx_sh_tot(:)     ! total sensible heat flux (W/m**2) [+ to atm]
    real(r8), pointer :: eflx_lh_tot(:)     ! total latent heat flux (W/m8*2)  [+ to atm]
    real(r8), pointer :: eflx_soil_grnd(:)  ! soil heat flux (W/m**2) [+ = into soil]
    real(r8), pointer :: qflx_evap_tot(:)   ! qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
    real(r8), pointer :: qflx_surf(:)       ! surface runoff (mm H2O /s)
    real(r8), pointer :: qflx_qrgwl(:)      ! qflx_surf at glaciers, wetlands, lakes
    real(r8), pointer :: qflx_drain(:)      ! sub-surface runoff (mm H2O /s)
!
! local pointers to original implicit out arrays
!
    real(r8), pointer :: errh2o(:)          ! water conservation error (mm H2O)
    real(r8), pointer :: errsol(:)          ! solar radiation conservation error (W/m**2)
    real(r8), pointer :: errlon(:)          ! longwave radiation conservation error (W/m**2)
    real(r8), pointer :: errseb(:)          ! surface energy conservation error (W/m**2)
!
! local pointers to original implicit in multi-level arrays
!
    real(r8), pointer :: forc_solad(:,:)    ! direct beam radiation (vis=forc_sols , nir=forc_soll )
    real(r8), pointer :: forc_solai(:,:)    ! diffuse radiation     (vis=forc_solsd, nir=forc_solld)
!
!EOP
!
! !OTHER LOCAL VARIABLES:
    integer  :: p,c,g                      ! indices
    real(r8) :: dtime                      ! land model time step (sec)
    integer  :: nstep                      ! time step number
    logical  :: found                      ! flag in search loop
    integer  :: index                      ! index of first found in search loop
!-----------------------------------------------------------------------

    ! Assign local pointers to derived type scalar members (gridcell-level)

    forc_rain         => clm3%g%a2lf%forc_rain
    forc_snow         => clm3%g%a2lf%forc_snow
    forc_lwrad        => clm3%g%a2lf%forc_lwrad
    forc_solad        => clm3%g%a2lf%forc_solad
    forc_solai        => clm3%g%a2lf%forc_solai

    ! Assign local pointers to derived type scalar members (column-level)

    cgridcell         => clm3%g%l%c%gridcell
    endwb             => clm3%g%l%c%cwbal%endwb
    begwb             => clm3%g%l%c%cwbal%begwb
    qflx_surf         => clm3%g%l%c%cwf%qflx_surf
    qflx_qrgwl        => clm3%g%l%c%cwf%qflx_qrgwl
    qflx_drain        => clm3%g%l%c%cwf%qflx_drain
    qflx_evap_tot     => clm3%g%l%c%cwf%pwf_a%qflx_evap_tot
    errh2o            => clm3%g%l%c%cwbal%errh2o

    ! Assign local pointers to derived type scalar members (pft-level)

    pgridcell         => clm3%g%l%c%p%gridcell
    fsa               => clm3%g%l%c%p%pef%fsa
    fsr               => clm3%g%l%c%p%pef%fsr
    eflx_lwrad_out    => clm3%g%l%c%p%pef%eflx_lwrad_out
    eflx_lwrad_net    => clm3%g%l%c%p%pef%eflx_lwrad_net
    sabv              => clm3%g%l%c%p%pef%sabv
    sabg              => clm3%g%l%c%p%pef%sabg
    eflx_sh_tot       => clm3%g%l%c%p%pef%eflx_sh_tot
    eflx_lh_tot       => clm3%g%l%c%p%pef%eflx_lh_tot
    eflx_soil_grnd    => clm3%g%l%c%p%pef%eflx_soil_grnd
    errsol            => clm3%g%l%c%p%pebal%errsol
    errseb            => clm3%g%l%c%p%pebal%errseb
    errlon            => clm3%g%l%c%p%pebal%errlon

    ! Get step size and time step

    nstep = get_nstep()
    dtime = get_step_size()

    ! Water balance check

!dir$ concurrent
!cdir nodep
    do c = lbc, ubc
       g = cgridcell(c)

       errh2o(c) = endwb(c) - begwb(c) &
            - (forc_rain(g) + forc_snow(g) - qflx_evap_tot(c) - qflx_surf(c) &
            - qflx_qrgwl(c) - qflx_drain(c)) * dtime
    end do

    found = .false.
    do c = lbc, ubc
       if (abs(errh2o(c)) > .10) then
          found = .true.
          index = c
          exit
       end if
    end do
    if ( found ) then
       write(6,200)'water balance error',nstep,index,errh2o(index)
       write(6,*)'clm model is stopping'
       call endrun()
    end if

    ! Energy balance checks

!dir$ concurrent
!cdir nodep
    do p = lbp, ubp
       g = pgridcell(p)

       ! Solar radiation energy balance
       errsol(p) = fsa(p) + fsr(p) - (forc_solad(g,1) + forc_solad(g,2) &
            + forc_solai(g,1) + forc_solai(g,2))

       ! Longwave radiation energy balance
       errlon(p) = eflx_lwrad_out(p) - eflx_lwrad_net(p) - forc_lwrad(g)

       ! Surface energy balance
       errseb(p) = sabv(p) + sabg(p) + forc_lwrad(g) - eflx_lwrad_out(p) &
            - eflx_sh_tot(p) - eflx_lh_tot(p) - eflx_soil_grnd(p)
    end do

    ! Solar radiation energy balance check

    found = .false.
    do p = lbp, ubp
       if (abs(errsol(p)) > .10 ) then
          found = .true.
          index = p
          exit
       end if
    end do
    if ( found ) then
       write(6,100)'solar radiation balance error',nstep,index,errsol(index)
       write(6,*)'clm model is stopping'
       call endrun()
    end if

    ! Longwave radiation energy balance check

    found = .false.
    do p = lbp, ubp
       if (abs(errlon(p)) > .10 ) then
          found = .true.
          index = p
          exit
       end if
    end do
    if ( found ) then
       write(6,100)'longwave enery balance error',nstep,index,errlon(index)
       write(6,*)'clm model is stopping'
       call endrun()
    end if

    ! Surface energy balance check

    found = .false.
    do p = lbp, ubp
       if (abs(errseb(p)) > .10 ) then
          found = .true.
          index = p
          exit
       end if
    end do
    if ( found ) then
       write(6,100)'surface flux energy balance error',nstep,index,errseb(index)
       write(6,*)'clm model is stopping'
       call endrun()
    end if

100 format (1x,a14,' nstep =',i10,' point =',i6,' imbalance =',f8.2,' W/m2')
200 format (1x,a14,' nstep =',i10,' point =',i6,' imbalance =',f8.2,' mm')

  end subroutine BalanceCheck

end module BalanceCheckMod
