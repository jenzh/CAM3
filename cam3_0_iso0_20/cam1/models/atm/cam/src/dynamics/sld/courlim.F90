#include <misc.h>
#include <params.h>

subroutine courlim (vmax2d, vmax2dt, vcour)
!-----------------------------------------------------------------------
!
! Purpose:
! Find out whether Courant limiter needs to be applied
!
! Method: 
! 
! Author: 
! Original version:  CCM2
! 
!-----------------------------------------------------------------------
!
! $Id: courlim.F90,v 1.5.2.4 2003/12/15 18:53:04 hender Exp $
! $Author: hender $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use physconst, only: rga
  use time_manager, only: get_nstep, is_first_step
  use comhd
#ifdef TIMING_BARRIERS
  use mpishorthand
#endif

  implicit none

#include <comsta.h>

!
! Arguments
!
  real(r8), intent(inout)   :: vmax2d (plev,plat)  ! Max. wind at each level, latitude
  real(r8), intent(inout)   :: vmax2dt(plev,plat)  ! Max. truncated wind at each lvl,lat
  real(r8), intent(inout)   :: vcour  (plev,plat)  ! Maximum Courant number in slice
!
!--------------------------Local Variables------------------------------
!
  integer  k,lat               ! Indices
  integer :: nstep             ! Current timestep number

  real(r8) vcourmax            ! Max courant number in the vertical wind field
  real(r8) vmaxt               ! Temporary var. for max. wind speed
  real(r8) vmax1dt(plev)       ! Sqrt of max wind speed
  real(r8) cn                  ! Estimate of truncated Courant number 
  real(r8) cnmax               ! Max. courant no. horiz. wind field
  real(r8) psurfsum            ! Summing variable - global mass
  real(r8) stqsum              ! Summing variable - global moisture
  real(r8) rmsdsum             ! Summing variable - global divergence
  real(r8) rmstsum             ! Summing variable - global temperature
  real(r8) stps                ! Global Mass integral
  real(r8) stqf                ! Global Moisture integral
  real(r8) rmsdf               ! Global RMS Divergence
  real(r8) rmstf               ! Global RMS Temperature
!
!-----------------------------------------------------------------------
!
#if ( defined SPMD )
#if ( defined TIMING_BARRIERS )
  call t_startf ('sync_realloc7')
  call mpibarrier (mpicom)
  call t_stopf ('sync_realloc7')
#endif
  call t_startf ('realloc7')
  call realloc7 (vmax2d  ,vmax2dt ,vcour   )
  call t_stopf ('realloc7')
#endif

  nstep = get_nstep()
!
! Compute maximum wind speed for each level
!
  do k=1,plev
     vmax1dt(k) = sqrt (maxval (vmax2dt(k,:)))
  end do
!
! Compute max. vertical Courant number (k is index to Model interfaces)
!
  vcourmax = maxval (vcour(2:,:))
!
! Compute Max Courant # for whole atmosphere for diagnostic output
!
  cnmax = 0.
  do k=1,plev-1
     cn = vmax1dt(k)*cnfac   ! estimate of Courant number
     cnmax = max(cnmax,cn)
  end do
!
! Write out statisitics to standard output
!
  psurfsum = 0.
  stqsum = 0.
  rmsdsum = 0.
  rmstsum = 0.

  do lat=1,plat
     psurfsum = psurfsum + psurf(lat)
     stqsum = stqsum + stq(lat)
     rmsdsum = rmsdsum + rmsd(lat)
     rmstsum = rmstsum + rmst(lat)
  end do

  stps = 0.5*psurfsum
  stqf = 0.5*rga*stqsum
  rmsdf = sqrt(0.5*rmsdsum)
  rmstf = sqrt(0.5*rmstsum)
  if (masterproc) then
     if (is_first_step()) write(6,810)
     write (6,820) nstep, rmsdf, rmstf, stps, stqf, cnmax, vcourmax
  end if
!
  return
!
! Formats
!
810 format(/101x,'COURANT'/10x,'NSTEP',4x,'RMSD',19x, &
                 'RMST',19x,'STPS',9x,'STQ',19x,'HOR   VERT')
820 format(' NSTEP =',i8,1x,1p2e23.15,1p1e13.5,e23.15,0p1f5.2,f6.2)

end subroutine courlim

