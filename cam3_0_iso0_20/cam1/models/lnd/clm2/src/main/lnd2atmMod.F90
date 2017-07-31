#include <misc.h>
#include <preproc.h>

module lnd2atmMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: lnd2atmMod
!
! !DESCRIPTION:
! Compute l2a component of gridcell derived type
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: lnd2atm
!
! !REVISION HISTORY:
! 03-04-27 : Created by Mariana Vertenstein
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: makel2a
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: lnd2atm
!
! !INTERFACE: subroutine lnd2atm(init)
  subroutine lnd2atm(init)
!
! !DESCRIPTION:
! Compute l2a component of gridcell derived type
!
! !USES:
    use decompMod, only : get_proc_clumps, get_clump_bounds
!
! !ARGUMENTS:
    implicit none
    logical, optional, intent(in) :: init  ! if true=>only set a subset of arguments
!
! !REVISION HISTORY:
! Mariana Vertenstein: created 03/10-25
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: nc              ! clump index
    integer :: nclumps         ! number of clumps on this processor
    integer :: begp, endp      ! per-proc beginning and ending pft indices
    integer :: begc, endc      ! per-proc beginning and ending column indices
    integer :: begl, endl      ! per-proc beginning and ending landunit indices
    integer :: begg, endg      ! per-proc gridcell ending gridcell indices
!------------------------------------------------------------------------

    ! Determine clump bounds for this processor

    nclumps = get_proc_clumps()

    ! Loop over clumps on this processor

!$OMP PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
!CSD$ PARALLEL DO PRIVATE (nc,begg,endg,begl,endl,begc,endc,begp,endp)
    do nc = 1,nclumps

       call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

       if (present(init)) then
          call makel2a(begp, endp, begc, endc, begg, endg, init=.true.)
       else
          call makel2a(begp, endp, begc, endc, begg, endg)
       end if

    end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

  end subroutine lnd2atm

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: makel2a
!
! !INTERFACE: subroutine lnd2atm(init)
  subroutine makel2a(lbp, ubp, lbc, ubc, lbg, ubg, init)
!
! !DESCRIPTION:
! Compute l2a component of gridcell derived type
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clmtype
    use clm_varcon  , only : sb
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbp, ubp   ! per-proc beginning and ending pft indices
    integer, intent(in) :: lbc, ubc   ! per-proc beginning and ending column indices
    integer, intent(in) :: lbg, ubg   ! per-proc gridcell ending gridcell indices
    logical, optional, intent(in) :: init  ! if true=>only set a subset of arguments
!
! !REVISION HISTORY:
! 03-04-27 : Created by Mariana Vertenstein
! 03-08-25 : Updated to vector data structure (Mariana Vertenstein)
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: pi,p,c,l,g   ! indices
    real(r8):: wt           ! temporary weight
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Initialize gridcell land->atm components before grid averaging

    gptr%l2as%albd(lbg:ubg,:) = 0.
    gptr%l2as%albi(lbg:ubg,:) = 0.
    gptr%l2as%t_ref2m(lbg:ubg) = 0.
    gptr%l2as%q_ref2m(lbg:ubg) = 0.
    gptr%l2as%h2osno(lbg:ubg) = 0.
    gptr%l2af%taux(lbg:ubg) = 0.
    gptr%l2af%tauy(lbg:ubg) = 0.
    gptr%l2af%eflx_lh_tot(lbg:ubg) = 0.
    gptr%l2af%eflx_sh_tot(lbg:ubg) = 0.
    gptr%l2af%eflx_lwrad_out(lbg:ubg) = 0.
    gptr%l2af%qflx_evap_tot(lbg:ubg) = 0.
    gptr%l2af%fsa(lbg:ubg) = 0.
    gptr%l2as%t_rad(lbg:ubg) = 0.

    ! Compute gridcell averages. Note that gridcell value for the
    ! radiative temperature (l2as%t_rad) must be computed after the
    ! gridcell average of eflx_lwrad_out is computed.

   if (present(init)) then

      if (init) then

         do c = lbc,ubc
            g = cptr%gridcell(c)
            wt = cptr%wtgcell(c)
            gptr%l2as%h2osno(g) = gptr%l2as%h2osno(g) + cptr%cws%h2osno(c)/1000. * wt
         end do

         do pi = 1,maxpatch
!dir$ concurrent
!cdir nodep
            do g = lbg,ubg
               if ( pi <=  gptr%npfts(g) ) then
                  p = gptr%pfti(g) + pi - 1
                  wt = pptr%wtgcell(p)
                  gptr%l2as%albd(g,1)         = gptr%l2as%albd(g,1)         + pptr%pps%albd(p,1) * wt
                  gptr%l2as%albd(g,2)         = gptr%l2as%albd(g,2)         + pptr%pps%albd(p,2) * wt
                  gptr%l2as%albi(g,1)         = gptr%l2as%albi(g,1)         + pptr%pps%albi(p,1) * wt
                  gptr%l2as%albi(g,2)         = gptr%l2as%albi(g,2)         + pptr%pps%albi(p,2) * wt
                  gptr%l2af%eflx_lwrad_out(g) = gptr%l2af%eflx_lwrad_out(g) + pptr%pef%eflx_lwrad_out(p) * wt
               end if
            end do
         end do

!dir$ concurrent
!cdir nodep
         do g = lbg,ubg
            gptr%l2as%t_rad(g) = sqrt(sqrt(gptr%l2af%eflx_lwrad_out(g)/sb))
         end do
      end if

   else

      do c = lbc,ubc
         g = cptr%gridcell(c)
         wt = cptr%wtgcell(c)
         gptr%l2as%h2osno(g) = gptr%l2as%h2osno(g) + cptr%cws%h2osno(c)/1000. * wt
      end do

      do pi = 1,maxpatch
!dir$ concurrent
!cdir nodep
         do g = lbg,ubg
            if ( pi <=  gptr%npfts(g) ) then
               p = gptr%pfti(g) + pi - 1
               wt = pptr%wtgcell(p)
               gptr%l2as%albd(g,1)         = gptr%l2as%albd(g,1)         + pptr%pps%albd(p,1) * wt
               gptr%l2as%albd(g,2)         = gptr%l2as%albd(g,2)         + pptr%pps%albd(p,2) * wt
               gptr%l2as%albi(g,1)         = gptr%l2as%albi(g,1)         + pptr%pps%albi(p,1) * wt
               gptr%l2as%albi(g,2)         = gptr%l2as%albi(g,2)         + pptr%pps%albi(p,2) * wt
               gptr%l2as%t_ref2m(g)        = gptr%l2as%t_ref2m(g)        + pptr%pes%t_ref2m(p) * wt
               gptr%l2as%q_ref2m(g)        = gptr%l2as%q_ref2m(g)        + pptr%pes%q_ref2m(p) * wt
               gptr%l2af%taux(g)           = gptr%l2af%taux(g)           + pptr%pmf%taux(p) * wt
               gptr%l2af%tauy(g)           = gptr%l2af%tauy(g)           + pptr%pmf%tauy(p) * wt
               gptr%l2af%eflx_lh_tot(g)    = gptr%l2af%eflx_lh_tot(g)    + pptr%pef%eflx_lh_tot(p) * wt
               gptr%l2af%eflx_sh_tot(g)    = gptr%l2af%eflx_sh_tot(g)    + pptr%pef%eflx_sh_tot(p) * wt
               gptr%l2af%eflx_lwrad_out(g) = gptr%l2af%eflx_lwrad_out(g) + pptr%pef%eflx_lwrad_out(p) * wt
               gptr%l2af%qflx_evap_tot(g)  = gptr%l2af%qflx_evap_tot(g)  + pptr%pwf%qflx_evap_tot(p) * wt
               gptr%l2af%fsa(g)            = gptr%l2af%fsa(g)            + pptr%pef%fsa(p) * wt
            end if
         end do
      end do

!dir$ concurrent
!cdir nodep
      do g = lbg,ubg
         gptr%l2as%t_rad(g) = sqrt(sqrt(gptr%l2af%eflx_lwrad_out(g)/sb))
      end do

   end if

 end subroutine makel2a

end module lnd2atmMod
