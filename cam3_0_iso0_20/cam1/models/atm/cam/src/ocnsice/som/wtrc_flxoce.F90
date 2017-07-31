#include <misc.h>
#include <params.h>

subroutine wtrc_flxoce( npts   ,indx   ,rbot   ,zbot   ,qbot   , &
                        ubot   ,vbot   ,ts     ,ustar  ,re     , &
                        ssq    ,qflx )

!-----------------------------------------------------------------------
!
! Purpose: compute water tracer exchange from ocean
!
! Method:
!   Used diagnostics output from (./dom/)flxoce to ensure
!   quantities are exactly equal for constituent number 1.
!   Isotopic fractionation (equilibrium and kinetci) is applied, 
!   when needed.
!
!   Assume surface latent heating is governed by form:
!     E = fac (q - qs(ts))
!
!   where fac is some exchange efficiency and qs is the saturation
!   vapour mixing rati at the surface temperature. These are needed
!   from calling routine to solve isotopic equivilent.
!
!     Ei = fac (1-kmol) (qi - qs(ts)*Rocn/alpha)
!
!   To compute the kinetic drag modifneed also
!
!
!
! Author:
!   David Noone <dcn@caltech.edu> - Mon Jun 30 10:24:49 MDT 2003
! This version is for the SOM, but is the same as for DOM
!   David Noone <dcn@colorado.edu> - Thu May 26 13:39:51 MDT 2005
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols
  use constituents, only: pcnst, pnats
  use water_tracers, only: trace_water, wtrc_is_vap, iwspec, ixwti, ixwtx
  use water_isotopes, only: wisotope, wiso_kmol, wiso_alpl,wiso_get_roce, &
                               wiso_alpi
  use abortutils, only: endrun

  implicit none
  
!---------------------------- Arguments --------------------------------
!
  integer , intent(in)  :: npts         ! number of point to compute
  integer , intent(in)  :: indx(pcols)  ! indicies to points
  real(r8), intent(in)  :: rbot(pcols)  ! density of lowest layer (kg/m3)
  real(r8), intent(in)  :: zbot(pcols)  ! height of lowest level (m)
  real(r8), intent(in)  :: qbot(pcols,pcnst+pnats)  ! constituents at lowest level (m)
  real(r8), intent(in)  :: ubot(pcols)  ! U wind (m/s)
  real(r8), intent(in)  :: vbot(pcols)  ! V wind (m/s)
  real(r8), intent(in)  :: ts(pcols)    ! (sea) surface temperature K
  real(r8), intent(in)  :: ustar(pcols) ! friction velocity (m/s)
  real(r8), intent(in)  :: re(pcols)    ! Reynolds number ?
  real(r8), intent(in)  :: ssq(pcols)   ! s.hum. saturation at Ts
!
  real(r8), intent(inout) :: qflx(pcols,pcnst+pnats) ! constituentflux (kg/kg/s)
!
!------------------------- Local Variables -----------------------------
  real(r8) alpkn(pcols,pcnst+pnats)     ! kinetic fractionation efficiency (m)
  real(r8) tau                          ! stress
  real(r8) delq                         ! spec. hum. difference
  real(r8) qstar                        ! spec. hum,. mixing scale
  real(r8) Roce                         ! water tracer ratio of ocean surface
  real(r8) alpha                        ! fractionation factor
  integer ii,i                          ! column indicies
  integer m                             ! constituent index
!-----------------------------------------------------------------------
!
  if (.not. trace_water) then
    write(6,*) 'WTRC_FLXOCE: called, but not needed.'
    call endrun
  end if
  if (.not. wisotope) then
    write(6,*) 'WTRC_FLXOCE: isotopes needed for water tracers.',wisotope
    call endrun
  end if
!
! Compute the isotopic kinetic enrichment
!
  alpkn(:,:) = 0.
  do m = ixwti, ixwtx
    if (wtrc_is_vap(m)) then
!!      call wiso_kmolv10(pcols,npts,indx,iwspec(m),zbot,ubot,vbot,alpkn(1:pcols,m))
      call wiso_kmol(pcols,npts,indx,iwspec(m),rbot,zbot,ustar,alpkn(1:pcols,m))
    end if
  end do
!
! Compute the vapour deficit then, get the fluxes
!
  do m = ixwti, ixwtx
    if (wtrc_is_vap(m)) then
      do ii = 1,npts
        i = indx(ii)
!
! Do some tests:
!
        Roce  = wiso_get_roce(iwspec(m))
        alpha = wiso_alpl(iwspec(m),ts(i))
        delq  = qbot(i,m) - ssq(i)*Roce/alpha
!
!        write(*,*) 'ALPHA:',alpha,ts(i),alpkn(i,m)
!
        qstar = re(i)*delq
        tau   = rbot(i) * ustar(i) * ustar(i)
        qflx(i,m) = -tau*alpkn(i,m)*qstar/ustar(i)
!
      end do
    end if
  end do	!m, constituents
!
  return
end subroutine wtrc_flxoce


