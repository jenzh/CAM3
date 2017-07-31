#include <misc.h>
#include <params.h>

subroutine hordif1(rearth,phi)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Horizontal diffusion of z,d,t,q
! 
! Method: 
! 1. implicit del**2 form above level kmnhd4
! 2. implicit del**4 form at level kmnhd4 and below
! 3. courant number based truncation at level kmxhdc and above
! 4. increased del**2 coefficient at level kmxhd2 and above
!
! Computational note: this routine is multitasked by level, hence it 
! is called once for each k
!
! Author: 
! 
!-----------------------------------------------------------------------
!
! $Id: hordif1.F90,v 1.1.2.3 2004/04/28 23:43:47 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  implicit none

!------------------------------Arguments--------------------------------
  real(r8), intent(in)    :: rearth     ! radius of earth
  real(r8), intent(inout) :: phi(psp)   ! used in spectral truncation of phis
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer ir,ii            ! spectral indices       
  integer mr,mc            ! spectral indices
  real(r8) k42             ! Nominal  Del^4 diffusion coeff at T42
  real(r8) k63             ! Nominal  Del^4 diffusion coeff at T63
  real(r8) knn             ! Computed Del^4 diffusion coeff at TNN
  real(r8) tmp             ! temp space
  real(r8) hdfst4(pnmax)
  integer  expon
  integer  m,n             ! spectral indices
!-----------------------------------------------------------------------
!
! Compute Del^4 diffusion coefficient
!
  k42   = 1.e+16
  k63   = 5.e+15
  expon = 25

  if(pmax-1 <= 42) then
     knn = k42
  elseif(pmax-1 == 63) then
     knn = k63
  else
     if(pmax-1 < 63) then
        tmp = log(k42/k63)/log(63._r8*64._r8/42._r8/43._r8)
     else
        tmp = 2.
     endif
     knn = k63*(63._r8*64._r8/float(pmax)/float(pmax-1))**tmp
  endif
!
! Set the Del^4 diffusion coefficients for each wavenumber
!
  hdfst4(1) = 0.
  do n=2,pnmax
     hdfst4(n) = knn * (n*(n-1)*n*(n-1)  ) / rearth**4
  end do
!
! Set the horizontal diffusion factors for each wavenumer at this level
! del^4 diffusion is to be applied and compute time-split implicit
! factors.
!
  do m=1,pmmax
     mr = nstart(m)
     mc = 2*mr
     do n=1,nlen(m)
        ir = mc + 2*n - 1
        ii = ir + 1
        phi(ir)  = phi(ir)/(1. + 3600.*hdfst4(n+m-1))**expon
        phi(ii)  = phi(ii)/(1. + 3600.*hdfst4(n+m-1))**expon
     end do
  end do

  return
end subroutine hordif1
