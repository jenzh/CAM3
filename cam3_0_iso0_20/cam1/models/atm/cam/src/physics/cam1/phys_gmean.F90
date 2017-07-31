#include <misc.h>
module phys_gmean
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Perform mixed layer global calculations for energy conservation checks.
!
! Method: 
! Gather to a master processor who does all the work
!
! Author: Byron Boville from SOM code by Jim Rosinski/Bruce Briegleb
! Modified: P. Worley to aggregate calculations (4/04)
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,       only: pcols, begchunk, endchunk
   use rgrid,        only: nlon
   use commap,       only: w
   use phys_grid,    only: gather_chunk_to_field
   use pmgrid,       only: plon, plat, masterproc
#if ( defined SPMD )
   use mpishorthand
#endif

   implicit none

   CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine gmean (arr, arr_gmean, nflds)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the global mean of each field in "arr" in the physics 
! chunked decomposition
!
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in)  :: nflds                         ! number of fields
      real(r8), intent(in) :: arr(pcols,begchunk:endchunk,nflds) ! Input array, chunked
      real(r8), intent(out):: arr_gmean(nflds)              ! global means
!
! Local workspace
!
      real(r8) :: arr_field(plon,plat,nflds)   ! rectangular version of arr
      real(r8) :: zmean                        ! zonal mean value
      real(r8) :: tmean                        ! temp global mean value
      integer :: i, j, ifld                    ! longitude, latitude, field indices

      call t_startf ('gmean')
      call gather_chunk_to_field (1, 1, nflds, plon, arr, arr_field)

      if (masterproc) then
         do ifld=1,nflds
            tmean = 0.
            do j=1,plat
               zmean = 0.
               do i=1,nlon(j)
                  zmean = zmean + arr_field(i,j,ifld)
               end do
               tmean = tmean + zmean * 0.5*w(j)/nlon(j)
            end do
            arr_gmean(ifld) = tmean
         end do
      end if
#if ( defined SPMD )
      call mpibcast (arr_gmean, 3, mpir8, 0, mpicom)
#endif
      call t_stopf ('gmean')

      return
    end subroutine gmean
end module phys_gmean
