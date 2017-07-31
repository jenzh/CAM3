module drydep_mod
#include <params.h>

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid
  use pmgrid,  only: plat
  use constituents, only: pcnst, pnats

      ! Shared Data for dry deposition calculation.

      real(r8) rair                ! Gas constant for dry air (J/K/kg)
      real(r8) gravit              ! Gravitational acceleration
      real(r8) phi(plat)           ! grid latitudes (radians)11

contains

!##############################################################################

! $Id: drydep_mod.F90,v 1.1.6.3 2004/01/27 17:10:56 rosinski Exp $
#include <params.h>

      subroutine inidrydep( xrair, xgravit, xphi )

! Initialize dry deposition parameterization.

      implicit none

! Input arguments:
      real(r8), intent(in) :: xrair                ! Gas constant for dry air
      real(r8), intent(in) :: xgravit              ! Gravitational acceleration
      real(r8), intent(in) :: xphi(plat)           ! grid latitudes (radians)

! Local variables:
      integer i, j, ncid, vid
!-----------------------------------------------------------------------

      rair = xrair
      gravit = xgravit
      do j = 1, plat
         phi(j) = xphi(j)
      end do

      return
      end subroutine inidrydep

!##############################################################################

      subroutine setdvel( ncol, landfrac, icefrac, ocnfrac, vgl, vgo, vgsi, vg )

! Set the deposition velocity depending on whether we are over
! land, ocean, and snow/ice


      implicit none

! Input arguments:

      integer, intent(in) :: ncol
      real (r8), intent(in) :: landfrac(pcols)       ! land fraction
      real (r8), intent(in) :: icefrac(pcols)       ! ice fraction
      real (r8), intent(in) :: ocnfrac(pcols)       ! ocean fraction

      real(r8), intent(in) :: vgl                  ! dry deposition velocity in m/s (land)
      real(r8), intent(in) :: vgo                  ! dry deposition velocity in m/s (ocean)
      real(r8), intent(in) :: vgsi                 ! dry deposition velocity in m/s (snow/ice)

! Output arguments:
      real(r8), intent(out) ::  vg(pcols) ! dry deposition velocity in m/s

! Local variables:

      integer i
      real(r8) a


      do i = 1, ncol
         vg(i) = landfrac(i)*vgl + ocnfrac(i)*vgo + icefrac(i)*vgsi
!         if (ioro(i).eq.0) then
!            vg(i) = vgo
!         else if (ioro(i).eq.1) then
!            vg(i) = vgl
!         else
!            vg(i) = vgsi
!         endif
      end do

      return
      end subroutine setdvel

!##############################################################################

      subroutine ddflux( ncol, vg, q, p, tv, flux )

! Compute surface flux due to dry deposition processes.


      implicit none

! Input arguments:
      integer , intent(in) :: ncol
      real(r8), intent(in) ::    vg(pcols)  ! dry deposition velocity in m/s
      real(r8), intent(in) ::    q(pcols)   ! tracer conc. in surface layer (kg tracer/kg moist air)
      real(r8), intent(in) ::    p(pcols)   ! midpoint pressure in surface layer (Pa)
      real(r8), intent(in) ::    tv(pcols)  ! midpoint virtual temperature in surface layer (K)

! Output arguments:

      real(r8), intent(out) ::    flux(pcols) ! flux due to dry deposition in kg/m^s/sec

! Local variables:

      integer i

      do i = 1, ncol
         flux(i) = -vg(i) * q(i) * p(i) /(tv(i) * rair)
      end do

      return
      end subroutine ddflux

!##############################################################################
end module drydep_mod
