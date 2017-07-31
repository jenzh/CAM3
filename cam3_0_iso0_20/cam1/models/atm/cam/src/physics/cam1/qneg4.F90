#include <misc.h>
#include <params.h>

subroutine qneg4 (subnam  ,lchnk   ,ncol    ,ztodt   ,        &
                  qbot    ,srfrpdel,shflx   ,lhflx   ,qflx    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check if moisture flux into the ground is exceeding the total
! moisture content of the lowest model layer (creating negative moisture
! values).  If so, then subtract the excess from the moisture and
! latent heat fluxes and add it to the sensible heat flux.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Olson
! Modified for water tracers: David Noone <dcn@colorado.edu> - Tue Jul 13 17:40:44 MDT 2004
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,    only: get_lat_p, get_lon_p
   use physconst,    only: gravit, latvap
   use constituents, only: qmin, pcnst, pnats

   use water_tracers, only: trace_water, ixwti, ixwtx, wtrc_is_vap, iwspec, wtrc_ratio


   implicit none

!
! Input arguments
!
   character*8, intent(in) :: subnam         ! name of calling routine
!
   integer, intent(in) :: lchnk              ! chunk index
   integer, intent(in) :: ncol               ! number of atmospheric columns
!
   real(r8), intent(in) :: ztodt             ! two times model timestep (2 delta-t)
   real(r8), intent(in) :: qbot(pcols,pcnst+pnats)   ! moisture at lowest model level
   real(r8), intent(in) :: srfrpdel(pcols)   ! 1./(pint(K+1)-pint(K))
!
! Input/Output arguments
!
   real(r8), intent(inout) :: shflx(pcols)   ! Surface sensible heat flux (J/m2/s)
   real(r8), intent(inout) :: lhflx(pcols)   ! Surface latent   heat flux (J/m2/s)
   real(r8), intent(inout) :: qflx (pcols,pcnst+pnats)   ! surface water flux (kg/m^2/s)
!
!---------------------------Local workspace-----------------------------
!
   integer :: i,ii              ! longitude indices
   integer :: iw                ! i index of worst violator
   integer :: indxexc(pcols)    ! index array of points with excess flux
   integer :: nptsexc           ! number of points with excess flux
   integer :: m                 ! tracer index
!
   real(r8):: worst             ! biggest violator
   real(r8):: excess(pcols)     ! Excess downward sfc latent heat flux
   real(r8):: qfxo(pcols,pcnst+pnats)   ! initial tracer flux
   real(r8):: rat               ! tracer ratio

!
!-----------------------------------------------------------------------
!
! Store old value to input for water tracers
!
!CDIR$ IVDEP
   do m = 1, pcnst+pnats
     do i = 1, ncol
       qfxo(i,m) = qflx(i,m)
     end do
   end do
!
! Compute excess downward (negative) q flux compared to a theoretical
! maximum downward q flux.  The theoretical max is based upon the
! given moisture content of lowest level of the model atmosphere.
!
   nptsexc = 0
   do i = 1,ncol
      excess(i) = qflx(i,1) - (qmin(1) - qbot(i,1))/(ztodt*gravit*srfrpdel(i))
!
! If there is an excess downward (negative) q flux, then subtract
! excess from "qflx" and "lhflx" and add to "shflx".
!
      if (excess(i) < 0.) then
         nptsexc = nptsexc + 1
         indxexc(nptsexc) = i
         qflx (i,1) = qflx (i,1) - excess(i)
         lhflx(i) = lhflx(i) - excess(i)*latvap
         shflx(i) = shflx(i) + excess(i)*latvap
      end if
   end do
!
! Write out worst value if excess
!
   if (nptsexc.gt.0) then
      worst = 0.
      do ii=1,nptsexc
         i = indxexc(ii)
         if (excess(i) < worst) then
            worst = excess(i)
            iw = i
         end if
      end do
      write(6,9000) subnam,get_lat_p(lchnk,iw),nptsexc,worst,get_lon_p(lchnk,iw)
   end if
!
! Water tracers: where total has change, modify tracers to conserve ! ratios
!
   if (trace_water) then
     do m = ixwti, ixwtx
       if (wtrc_is_vap(m)) then
         do ii = 1, nptsexc
           i = indxexc(ii)
           rat = wtrc_ratio(qfxo(i,m),qfxo(i,1))
           qflx(i,m) = qflx(i,m) - rat*excess(i)
         end do
       end if
     end do
   end if
!
   return
9000 format(' QNEG4 WARNING from ',a8,', lchnk = ',i3,';', &
            ' Max possible LH flx exceeded at ',i4,' points. ', &
            ' Worst excess = ',1pe12.4,' at i = ',i4)
end subroutine qneg4
