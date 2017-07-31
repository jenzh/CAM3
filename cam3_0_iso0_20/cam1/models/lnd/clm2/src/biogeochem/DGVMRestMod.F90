#include <misc.h>
#include <preproc.h>

module DGVMRestMod

#if (defined DGVM)
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: restDGVMMod
!
! !DESCRIPTION:
! Read/Write to/from DGVM info to CLM restart file.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils  , only: endrun
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: restart_dgvm
!
! !REVISION HISTORY:
! Module created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_dgvm
!
! !INTERFACE:
  subroutine restart_dgvm (nio, flag)
!
! !DESCRIPTION:
! Read/write DGVM restart data
!
! !USES:
    use clmtype
    use iobinary
    use decompMod, only : get_proc_bounds, get_proc_global
    use DGVMMod  , only : resetWeightsDGVM, gatherWeightsDGVM
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             ! restart unit
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,p          ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    integer :: ier          ! error status
    integer , pointer :: ibuf1dp(:) ! pointer to memory to be allocated
    real(r8), pointer :: rbuf1dp(:) ! pointer to memory to be allocated
    real(r8), pointer :: rbuf1dc(:) ! pointer to memory to be allocated
!-----------------------------------------------------------------------

    ! Determine necessary subgrid bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

   ! Allocate necessary 1d buffers

    allocate (ibuf1dp(nump), rbuf1dp(nump), stat=ier)
    if (ier /= 0) then
       write (6,*) 'restart_dgvm: allocation error '
       call endrun
    end if

    !column type dgvm physical state variable - wf
    allocate (rbuf1dc(numc))
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=namec)
!dir$ concurrent
!cdir nodep
    do c = begc,endc
       if (flag == 'read' ) clm3%g%l%c%cps%wf(c) = rbuf1dc(c)
       if (flag == 'write') rbuf1dc(c) = clm3%g%l%c%cps%wf(c)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=namec)
    deallocate (rbuf1dc)

    ! pft type dgvm physical state - t_mo_min
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%t_mo_min(p)  = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%t_mo_min(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - annpsn
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%annpsn(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%annpsn(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - annpsnpot
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%annpsnpot(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%annpsnpot(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft cflux tye - fmicr
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pcf%fmicr(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pcf%fmicr(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - bm_inc
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%bm_inc(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%bm_inc(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - afmicr
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%afmicr(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%afmicr(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - t10min
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%t10min(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%t10min(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - tmomin20
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%tmomin20(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%tmomin20(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - agdd20
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%agdd20(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%agdd20(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - itypveg
    if (flag == 'read') call readin (nio, ibuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%itype(p) = ibuf1dp(p)
       if (flag == 'write') ibuf1dp(p) = clm3%g%l%c%p%itype(p)
    end do
    if (flag == 'write') call wrtout (nio, ibuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - fpcgrid
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%fpcgrid(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%fpcgrid(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - lai_ind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%lai_ind(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%lai_ind(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - crownarea
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%crownarea(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%crownarea(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - dphen
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%dphen(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%dphen(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - leafon
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%leafon(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%leafon(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - leafof
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%leafof(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%leafof(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - firelength
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%firelength(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%firelength(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - litterag
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%litterag(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%litterag(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - litterbg
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%litterbg(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%litterbg(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - cpool_fast
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%cpool_fast(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%cpool_fast(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - cpool_slow
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%cpool_slow(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%cpool_slow(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - k_fast_ave
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%k_fast_ave(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%k_fast_ave(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - k_slow_ave
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read') clm3%g%l%c%p%pdgvs%k_slow_ave(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%k_slow_ave(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - litter_decom_ave
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%litter_decom_ave(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%litter_decom_ave(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - nind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%nind(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%nind(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - lm_ind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%lm_ind(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%lm_ind(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - sm_ind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%sm_ind(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%sm_ind(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - hm_ind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%hm_ind(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%hm_ind(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - rm_ind
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pdgvs%rm_ind(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pdgvs%rm_ind(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type dgvm physical state - present
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read') then
          clm3%g%l%c%p%pdgvs%present(p) = .false.
          if (rbuf1dp(p) == 1.0) clm3%g%l%c%p%pdgvs%present(p) = .true.
       end if
       if (flag == 'write') then
          rbuf1dp(p) = 0
          if (clm3%g%l%c%p%pdgvs%present(p)) rbuf1dp(p) = 1.0
       end if
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! Determine new pft, column and land properties that result from
    ! restart data input

    if (flag == 'read') then
       call resetWeightsDGVM(begg, endg, begc, endc, begp, endp)
#ifdef SPMD
       call gatherWeightsDGVM()
#endif
    end if

    deallocate (ibuf1dp, rbuf1dp)

  end subroutine restart_dgvm

#endif

end module DGVMRestMod
