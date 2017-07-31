#include <params.h>
#include <misc.h>

subroutine wiso_cldliq_equilibrate(state, ptend, ztodt)
!-----------------------------------------------------------------------
!
! Purpose: Performs isotopic equilibration between (cloud) liquid and
!    large scale vapour. This is called after non-fractionatin
!    sedimentation, and is approximatedly valid as the cloud droplets are
!    very small. Assume ice have tiny amount of equilibration during
!    a single timestep. As this is so small, we can in fact to this as
!    two parts, rather than a single 3-body problem. 
!
!    For liq., it is assumed that it reaches the equilibrium  
!    fractionation at a time scale given by Jouzel (1986) pg. 67 
!
!    (could also equilibrate surface, precipitation flux if needed)
!
! I could possibly be convinced that we should equilibrate only wit the
! cloudy part of the grid cell
!
! THIS EQUILIBRATION NEEDS TO USE AN ADJUSTMENT TIME SCALE.
!
!
! Code history: Was hacked into cld_sediment_tend, and uses Sun Wong's
!               diffusion limitation for equilibration.
!
! Author: David Noone <dcn@colorado.edu> - Wed Apr  7 19:57:37 MDT 2004
!
!
!-----------------------------------------------------------------------

  use shr_kind_mod,     only: r8=>shr_kind_r8
  use ppgrid,           only: pcols, pver
  use constituents,     only: cnst_get_ind
  use physics_types,    only: physics_state, physics_ptend
  use water_tracers,    only: wtrc_is_vap, ixwti, ixwtx, iwspec
  use water_isotopes,   only: wiso_alpl, wiso_alpi, wiso_liqvap_equil
  use wtrc_camtune,     only: feq_liq, feq_ice

!-----------------------------------------------------------------------
  implicit none
!------------------------- Input Arguments -----------------------------
  type(physics_state), intent(in)  :: state   ! state variables
  real(r8), intent(in)             :: ztodt	  ! timestep
!------------------------ Output Arguments -----------------------------
  type(physics_ptend), intent(out) :: ptend   ! package tendencies
!------------------------ Local Variables -----------------------------
  real(r8) :: frain(pcols,pver) ! fraction of vapour ramaining
  real(r8) :: qtot		! total isotopic water
  real(r8) :: qvap		! final isotopic vapour
  real(r8) :: vapiso 		! isotopic vapour
  real(r8) :: liqiso 		! isotopic cloud liquid
  real(r8) :: iceiso 		! isotopic cloud ice

  real(r8) :: alpliq		! fractionation factor
  real(r8) :: alpice		! fractionation factor

  real(r8) :: dliqvap		! change in liquid from vapour
  real(r8) :: dicevap		! change in ice from vapour

  real(r8) :: efac		! equilibration factor for vapour

  integer ncol			! number of columns
  integer ixcldliq		! prognostic cloud liquid index
  integer ixcldice		! prognostic cloud ice index
  integer mvap, mliq, mice	! triplet pointers
  integer i,k,m
!-----------------------------------------------------------------------
  real(r8) :: qtiny = 1.e-22		! small q for denominator
  real(r8) :: flim  = 1.e-7		! limit for useful frain 
!-----------------------------------------------------------------------
  ncol = state%ncol
  call cnst_get_ind('CLDLIQ', ixcldliq)
  call cnst_get_ind('CLDICE', ixcldice)
!
! Initialize parameterization tendency (m=1,ixcldliq, ixcldice, = 0)
!
  ptend%q(:,:,:) = 0.
  ptend%lq(:) = .false.
!
! Precompute liquid/vapour fraction. This form has better numerical
! paroperties than l/v
!
  frain(:ncol,:) = 0.
  do k = 1, pver
     do i = 1, ncol
         qtot = state%q(i,k,1) + state%q(i,k,ixcldliq)
         frain(i,k) = state%q(i,k,1) / max(qtot,qtiny)
         frain(i,k) = max(frain(i,k), 0.)
         frain(i,k) = min(frain(i,k), 1.)
    end do
  end do
!
! Loop over all water tracer tripplets
!
  do m = ixwti, ixwtx
    if (wtrc_is_vap(m)) then
      mvap = m
      mliq = m+1
      mice = m+2		! do nothing to the ice
!
      ptend%lq(mvap) = .true.
      ptend%lq(mliq) = .true.
      ptend%lq(mice) = .true.
!
! Solve isotopic budgets to infer tendencies
!
     do k = 1, pver
       do i = 1, ncol
!
! Construct isotopic tendency (notice dv = -dl for conservation)
! Only do this if liquid present (nominally affect 1 ten thousands of a permil)
!
          if (frain(i,k) < 1.0-flim .and. frain(i,k) > flim) then ! saves numerics
!
            vapiso = state%q(i,k,mvap)
            liqiso = state%q(i,k,mliq)
            iceiso = state%q(i,k,mice)

            alpliq = wiso_alpl(iwspec(m), state%t(i,k))
            alpice = wiso_alpi(iwspec(m), state%t(i,k))
!
            call wiso_liqvap_equil(alpliq,feq_liq, state%q(i,k,1),state%q(i,k,ixcldliq), &
                         vapiso,liqiso,dliqvap)

            if (feq_ice > 0.0) then	! note vapiso has been updated
              call wiso_liqvap_equil(alpice,feq_ice, state%q(i,k,1),state%q(i,k,ixcldice), &
                         vapiso,iceiso,dicevap)
            end if
!
            ptend%q(i,k,mliq) = dliqvap / ztodt
            ptend%q(i,k,mice) = dicevap / ztodt
            ptend%q(i,k,mvap) = -(ptend%q(i,k,mliq) + ptend%q(i,k,mice))
          end if	! some liquid available 
!
        end do	! i columns
      end do	! k levels
!
    end if	! m is a water vapour tracer
  end do	! m tracers
!
  return
end subroutine wiso_cldliq_equilibrate

