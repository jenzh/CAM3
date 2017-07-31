#include <misc.h>
#include <params.h>

subroutine wtrc_ice_checkscale(npts,indx,snowh,wtsnowh)
!-----------------------------------------------------------------------
! Corrects for any numerical imprecision by scaling tracers to prognostic
! while preserving tracer ratios
!-----------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8

  use ppgrid, only: pcols
  use constituents, only: pcnst, pnats
  use water_tracers, only: ixwti,ixwtx,wtrc_is_vap, wtrc_ratio,ixh2oq
  use abortutils,   only: endrun

  implicit none
!---------------------------- Arguments --------------------------------
  integer , intent(in)    :: npts              ! number of ice points
  integer , intent(in)    :: indx(pcols)       ! indicies to ice points
  real(r8), intent(in)    :: snowh(pcols)      ! snow height
  real(r8), intent(inout) :: wtsnowh(pcols,pcnst+pnats) ! tracer snow height

  real(r8) :: htol = 1.e-14           ! meters error allowed

!------------------------- Local Variables -----------------------------
  integer i,ii		! column indicies
  integer m		! tracer index

  real(r8) dice(pcols)	! error
  real(r8) rat		! tracer ratio to preserve
!-----------------------------------------------------------------------

  dice(:) = 0.
  do ii = 1, npts
    i = indx(ii)
    dice(i) = snowh(i) - wtsnowh(i,ixh2oq)
  end do
!
! Check is difference is too big. 
! This can happen when zm convection produces negatives in updrafts.
! Thus differences in rain, then evap, the pflx, then zm_evap
! difference, and ultimately precc, which is accumulated into snowfall
! and added to the snow pack...
!
  if (count(abs(dice)>htol) > 0) then
           write(*,*) '(wtrc_ice_checkscale) ice snowh differences. '
       write(*,*) count(abs(dice)>htol),' points'
       write(*,*) 'Worst:',maxval(abs(dice))
       call endrun('wtrc_ice_checkscale: aborted')
  end if
!
! All OK to required precision 
!  do a rescale based on trace ratio applied to prognostic mass
!
  do m = ixwti, ixwtx
    if (wtrc_is_vap(m)) then
      do ii = 1, npts
        i = indx(ii)
        rat = wtrc_ratio(wtsnowh(i,m), wtsnowh(i,ixh2oq))
        wtsnowh(i,m) = rat*wtsnowh(i,1)
      end do
    end if
  end do
!
  return
end subroutine wtrc_ice_checkscale

