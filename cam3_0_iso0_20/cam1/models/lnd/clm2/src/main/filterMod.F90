#include <misc.h>
#include <preproc.h>

module filterMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: filterMod
!
! !DESCRIPTION:
! Module of filters used for processing columns and pfts of particular
! types, including lake, non-lake, soil, snow, non-snow, and
! naturally-vegetated patches.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils, only : endrun
!
! !PUBLIC TYPES:
  implicit none
  save

  type clumpfilter
#ifdef DGVM
     integer, pointer :: natvegp(:) ! DGVM naturally-vegetated (present)
                                    ! filter (pfts)
     integer :: num_natvegp         ! number of pfts in naturally-vegetated
                                    ! filter
#endif
     integer, pointer :: lakep(:)   ! lake filter (pfts)
     integer :: num_lakep           ! number of pfts in lake filter
     integer, pointer :: nolakep(:) ! non-lake filter (pfts)
     integer :: num_nolakep         ! number of pfts in non-lake filter
     integer, pointer :: lakec(:)   ! lake filter (columns)
     integer :: num_lakec           ! number of columns in lake filter
     integer, pointer :: nolakec(:) ! non-lake filter (columns)
     integer :: num_nolakec         ! number of columns in non-lake filter
     integer, pointer :: soilc(:)   ! soil filter (columns)
     integer :: num_soilc           ! number of columns in soil filter
     integer, pointer :: snowc(:)   ! snow filter (columns)
     integer :: num_snowc           ! number of columns in snow filter
     integer, pointer :: nosnowc(:) ! non-snow filter (columns)
     integer :: num_nosnowc         ! number of columns in non-snow filter
  end type clumpfilter
  type(clumpfilter), allocatable, public :: filter(:)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2004.04.27 DGVM naturally-vegetated filter added by Forrest Hoffman
!
!EOP
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initFilters
!
! !INTERFACE:
  subroutine initFilters()
!
! !DESCRIPTION:
! Initialize CLM filters.
!
! !USES:
    use clmtype
    use decompMod , only : get_nclumps, get_proc_clumps, get_clump_bounds
    use clm_varcon, only : istsoil
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2004.04.27 DGVM naturally-vegetated filter added by Forrest Hoffman
!
!EOP
!
! LOCAL VARAIBLES:
    integer :: nc          ! clump index
    integer :: c,l,p       ! column, landunit, pft indices
    integer :: nclumps     ! total number of clumps on this processor
    integer :: fl          ! lake filter index
    integer :: fnl         ! non-lake filter index
    integer :: fs          ! soil filter index
    integer :: begp, endp  ! per-clump beginning and ending pft indices
    integer :: begc, endc  ! per-clump beginning and ending column indices
    integer :: begl, endl  ! per-clump beginning and ending landunit indices
    integer :: begg, endg  ! per-clump beginning and ending gridcell indices
    integer :: ier         ! error status
!------------------------------------------------------------------------

    ! Determine clump variables for this processor

    nclumps = get_proc_clumps()
    ier = 0
    if( .not. allocated(filter))allocate(filter(nclumps), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initFilters(): allocation error for clumpsfilters'
       call endrun
    end if

    ! Loop over clumps on this processor

    do nc = 1, nclumps

       ! Determine clump boundaries

       call get_clump_bounds(nc, begg, endg, begl, endl, begc, endc, begp, endp)

       ! Create lake and non-lake filters at column-level
       ! Allocate enough space for all columns

       allocate (filter(nc)%lakec(endc-begc+1))
       allocate (filter(nc)%nolakec(endc-begc+1))

       fl = 0
       fnl = 0
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          l = clm3%g%l%c%landunit(c)
          if (clm3%g%l%lakpoi(l)) then
             fl = fl + 1
             filter(nc)%lakec(fl) = c
          else
             fnl = fnl + 1
             filter(nc)%nolakec(fnl) = c
          end if
       end do
       filter(nc)%num_lakec = fl
       filter(nc)%num_nolakec = fnl

       ! Create lake and non-lake filters at pft-level
       ! Allocate enough space for all pfts

       allocate (filter(nc)%lakep(endp-begp+1))
       allocate (filter(nc)%nolakep(endp-begp+1))

       fl = 0
       fnl = 0
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          l = clm3%g%l%c%p%landunit(p)
          if (clm3%g%l%lakpoi(l)) then
             fl = fl + 1
             filter(nc)%lakep(fl) = p
          else
             fnl = fnl + 1
             filter(nc)%nolakep(fnl) = p
          end if
       end do
       filter(nc)%num_lakep = fl
       filter(nc)%num_nolakep = fnl

       ! Create column-level soil filter
       ! Allocate enough space for all columns

       allocate (filter(nc)%soilc(endc-begc+1))

       fs = 0
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          l = clm3%g%l%c%landunit(c)
          if (clm3%g%l%itype(l) == istsoil) then
             fs = fs + 1
             filter(nc)%soilc(fs) = c
          end if
       end do
       filter(nc)%num_soilc = fs

       ! Snow filters are reconstructed each time step in Hydrology2
       ! Just allocate enough space for all columns once here.

       allocate(filter(nc)%snowc(endc-begc+1))
       allocate(filter(nc)%nosnowc(endc-begc+1))

#ifdef DGVM
       ! DGVM present vegetated filter is reconstructed each time DGVM
       ! is run.  Just allocate enough space for all PFTs once here.
       allocate (filter(nc)%natvegp(endp-begp+1))
#endif

    end do

  end subroutine initFilters

end module FilterMod
