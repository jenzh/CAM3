#include <misc.h>
#include <params.h>

subroutine qmassa(cwava   ,w       ,q3      ,pdel    ,hw1lat  , &
                  nlon    ,q0, lat)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate contribution of current latitude to mass of constituents
! being advected by slt.
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id: qmassa.F90,v 1.1.2.3 2004/03/10 23:02:09 bundy Exp $
! $Author: bundy $
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use constituents, only: pcnst, pnats, cnst_get_type_byind

  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)  :: nlon                 ! longitude dimension  
  real(r8), intent(in)  :: cwava                ! normalization factor    l/(g*plon)
  real(r8), intent(in)  :: w                    ! gaussian weight this latitude
  real(r8), intent(in)  :: q3(plond,plev,pcnst) ! constituents
  real(r8), intent(in)  :: q0(plond,plev,pcnst) ! constituents at begining of time step
  real(r8), intent(in)  :: pdel(plond,plev)     ! pressure diff between interfaces
  real(r8), intent(out) :: hw1lat(pcnst)        ! accumulator
  integer lat
!-----------------------------------------------------------------------
!
!---------------------------Local variables-----------------------------
  integer i,k,m             ! longitude, level, constituent indices
  real(r8) const            ! temporary constant
!-----------------------------------------------------------------------
!
! Integration factor (the 0.5 factor arises because gaussian weights sum to 2)
!
  const = cwava*w*0.5
  do m=1,pcnst
     hw1lat(m) = 0.
  end do
!
! Compute mass integral for water ONLY
!
  do k=1,plev
     do i=1,nlon
        hw1lat(1) = hw1lat(1) + q3(i,k,1)*pdel(i,k)
     end do
  end do
!
! Compute mass integral for non-water constituents (on either WET or DRY basis)
!
  do m=2,pcnst
     if (cnst_get_type_byind(m).eq.'dry' ) then
        do k=1,plev
           do i=1,nlon
              hw1lat(m) = hw1lat(m) + q3(i,k,m)*(1. - q0(i,k,1))*pdel(i,k)
           end do
        end do
     else 
        do k=1,plev
           do i=1,nlon
              hw1lat(m) = hw1lat(m) + q3(i,k,m)*(1. - q3(i,k,1))*pdel(i,k)
           end do
        end do
     end if
  end do  

  do m = 1,pcnst
     hw1lat(m) = hw1lat(m)*const
  end do

  return
end subroutine qmassa

