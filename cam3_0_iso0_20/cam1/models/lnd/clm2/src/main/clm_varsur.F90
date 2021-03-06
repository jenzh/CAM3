#include <misc.h>
#include <preproc.h>

module clm_varsur

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varsur
!
! !DESCRIPTION:
! Module containing 2-d surface boundary data information
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar, only : lsmlon, lsmlat, nlevsoi
  use abortutils, only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! land model grid
!
  integer :: numlon(lsmlat)             !longitude points for each latitude strip
  real(r8):: latixy(lsmlon,lsmlat)      !latitude of grid cell (degrees)
  real(r8):: longxy(lsmlon,lsmlat)      !longitude of grid cell (degrees)
  real(r8):: area(lsmlon,lsmlat)        !grid cell area (km**2)
  real(r8):: landarea                   !total land area for all gridcells (km^2)
  real(r8):: lats(lsmlat+1)             !grid cell latitude, southern edge (degrees)
  real(r8):: lonw(lsmlon+1,lsmlat)      !grid cell longitude, western edge (degrees)
  real(r8):: lsmedge(4)                 !North,East,South,West edges of grid (deg)
  logical :: pole_points                !true => grid has pole points
  logical :: fullgrid  = .true.         !true => no grid reduction towards poles
  logical :: offline_rdgrid             !true => read offline grid rather than creating it
!
! fractional land and mask
!
  integer  landmask(lsmlon,lsmlat)      !land mask: 1 = land. 0 = ocean
  real(r8) landfrac(lsmlon,lsmlat)      !fractional land
!
! surface boundary data
!
  integer , allocatable :: soic2d(:,:)   !soil color
  real(r8), allocatable :: sand3d(:,:,:) !soil texture: percent sand
  real(r8), allocatable :: clay3d(:,:,:) !soil texture: percent clay
  real(r8), allocatable :: pctgla(:,:)   !percent of grid cell that is glacier
  real(r8), allocatable :: pctlak(:,:)   !percent of grid cell that is lake
  real(r8), allocatable :: pctwet(:,:)   !percent of grid cell that is wetland
  real(r8), allocatable :: pcturb(:,:)   !percent of grid cell that is urbanized
  logical, public :: all_pfts_on_srfdat = .false.  ! old format dataset is used by default
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: varsur_alloc    !allocates 2d surface data needed for initialization
  public :: varsur_dealloc  !deallocates 2d surface data needed for initialization
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: varsur_alloc
!
! !INTERFACE:
  subroutine varsur_alloc
!
! !DESCRIPTION:
! Allocate dynamic memory for module variables
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! LOCAL VARIABLES:
    integer :: ier                    !error code
!-----------------------------------------------------------------------

    allocate (soic2d(lsmlon,lsmlat), &
              sand3d(lsmlon,lsmlat,nlevsoi), &
              clay3d(lsmlon,lsmlat,nlevsoi), &
              pctgla(lsmlon,lsmlat), &
              pctlak(lsmlon,lsmlat), &
              pctwet(lsmlon,lsmlat), &
              pcturb(lsmlon,lsmlat), stat=ier)
    if (ier /= 0) then
       write(6,*)'varsur_alloc(): allocation error'
       call endrun
    endif

  end subroutine varsur_alloc

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: varsur_dealloc
!
! !INTERFACE:
  subroutine varsur_dealloc
!
! !DESCRIPTION:
! Deallocate dynamic memory for module variables
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    deallocate (soic2d, &
                sand3d, &
                clay3d, &
                pctgla, &
                pctlak, &
                pctwet, &
                pcturb)

  end subroutine varsur_dealloc

end module clm_varsur
