#include <misc.h>
#include <params.h>

module spmd_dyn

!----------------------------------------------------------------------- 
! 
! Purpose: SPMD implementation of CAM spectral SLD dynamics.
! 
! Author: CCM Core Group
! Modified: P. Worley, November 2003, December 2003
! 
!-----------------------------------------------------------------------

#if (defined SPMD)

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plat, masterproc, iam, beglatex, endlatex, numbnd, numlats, numlatsex, &
                           beglat, endlat, begirow, endirow, plev
   use constituents, only: pcnst
   use mpishorthand, only: mpir8, mpicom
   use infnan, only: inf
   use abortutils, only: endrun
 
   implicit none
   save

   private
   public spmdinit_dyn, compute_gsfactors, spmdbuf_resize
   public spmd_dyn_defaultopts, spmd_dyn_setopts

   integer, public :: npes                 ! Total number of MPI tasks
   integer, public :: nsmps                ! Total number of SMP nodes
   integer, public :: cut(2,0:plat-1)      ! partition for MPI tasks
   integer, public :: cutex(2,0:plat-1)    ! extended partition 
   integer, public :: proc(plat)           ! MPI task id associated with a given lat.
   integer, public :: neighs               ! number of south neighbors to comm guardcells
   integer, public, allocatable :: neighs_proc(:)    ! sorted south process neighbors
   integer, public :: neighn               ! number of north neighbors to comm guardcells
   integer, public, allocatable :: neighn_proc(:)    ! sorted north process neighbors
   integer, public :: npessp               ! number of MPI tasks in spectral space
   integer, public :: bsiz                 ! buffer size (in bytes)
   integer, public :: spmdbuf_siz          ! buffer size (in r8s)
   integer, public :: maxlats              ! max number of lats on any MPI task
   integer, public, allocatable :: nlat_p(:)    ! number of latitudes per MPI task
   integer, public, allocatable :: proc_smp_map(:) ! map of process/SMP node assignments
   integer, public :: realloc4_steps       ! number of swaps in realloc4 algorithms
   integer, public, allocatable :: realloc4_proc(:)
                                           ! swap partner in each step of 
                                           ! realloc4 algorithms
   integer, public, allocatable :: realloc4_step(:)
                                           ! step in realloc4 algorithms
                                           ! in which communicate with a given
                                           ! process
   integer, public, allocatable :: allgather_proc(:)
                                           ! swap partner in each step of 
                                           ! allgather (realloc5/7) algorithm
   integer, public, allocatable :: allgather_step(:)
                                           ! step in allgather (realloc5/7) algorithm
                                           ! in which communicate with a given
                                           ! process
!
   logical, private, parameter :: def_mirror = .false.          ! default
   logical, private :: mirror = def_mirror ! flag indicating whether latitudes and their
                                           ! reflections across the equator should assigned 
                                           ! to consecutive processes
!
! Dynamics communication transpose algorithm option:
!  0: use mpi_alltoallv
!  1: use point-to-point implementation
   integer, private, parameter :: min_alltoall = 0
   integer, private, parameter :: max_alltoall = 1
   integer, private, parameter :: def_alltoall = 0         ! default
   integer, public :: dyn_alltoall  = def_alltoall
!
! Dynamics communication allgather (realloc5/7) algorithm option:
!  0: use mpi_allgatherv
!  1: use point-to-point implementation
   integer, private, parameter :: min_allgather = 0
   integer, private, parameter :: max_allgather = 1
   integer, private, parameter :: def_allgather = 0         ! default
   integer, public :: dyn_allgather = def_allgather
!
   real(r8), public, allocatable :: buf1(:),buf2(:) ! buffers for packing MPI msgs

CONTAINS

!========================================================================

  subroutine spmd_dyn_defaultopts(npr_yz_out, geopktrans_out,    &
               tracertrans_out, ompnest_out, force_2d_out,       &
               modcomm_transpose_out, modcomm_geopk_out,         &
               dyn_alltoall_out, dyn_allgather_out  )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: Art Mirin / Pat Worley
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
! FV-only arguments
     integer, intent(out), optional :: npr_yz_out(4)
     integer, intent(out), optional :: geopktrans_out
     integer, intent(out), optional :: tracertrans_out
     integer, intent(out), optional :: ompnest_out
     integer, intent(out), optional :: force_2d_out
     integer, intent(out), optional :: modcomm_transpose_out
     integer, intent(out), optional :: modcomm_geopk_out
! EUL/SLD arguments
     integer, intent(out), optional :: dyn_alltoall_out
     integer, intent(out), optional :: dyn_allgather_out
!-----------------------------------------------------------------------
     if ( present(dyn_alltoall_out) ) then
       dyn_alltoall_out = def_alltoall
     endif
     if ( present(dyn_allgather_out) ) then
       dyn_allgather_out = def_allgather
     endif
  end subroutine spmd_dyn_defaultopts

!========================================================================

  subroutine spmd_dyn_setopts(npr_yz_in, geopktrans_in,       &
               tracertrans_in, ompnest_in, force_2d_in,       &
               modcomm_transpose_in, modcomm_geopk_in,        &
               dyn_alltoall_in, dyn_allgather_in  )
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: Art Mirin / Pat Worley
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
! FV-only arguments
     integer, intent(in), optional :: npr_yz_in(4)
     integer, intent(in), optional :: geopktrans_in
     integer, intent(in), optional :: tracertrans_in
     integer, intent(in), optional :: ompnest_in
     integer, intent(in), optional :: force_2d_in
     integer, intent(in), optional :: modcomm_transpose_in
     integer, intent(in), optional :: modcomm_geopk_in
! EUL/SLD arguments
     integer, intent(in), optional :: dyn_alltoall_in
     integer, intent(in), optional :: dyn_allgather_in
!-----------------------------------------------------------------------
     if ( present(dyn_alltoall_in) ) then
       dyn_alltoall = dyn_alltoall_in
       if ((dyn_alltoall.lt.min_alltoall).or. &
           (dyn_alltoall.gt.max_alltoall)) then
         write(6,*)                                          &
           'SPMD_DYN_SETOPTS:  ERROR:  dyn_alltoall=', &
           dyn_alltoall_in,                              &
           '  is out of range.  It must be between ',        &
           min_alltoall,' and ',max_alltoall
         call endrun
       endif
     endif
     if ( present(dyn_allgather_in) ) then
       dyn_allgather = dyn_allgather_in
       if ((dyn_allgather.lt.min_allgather).or. &
           (dyn_allgather.gt.max_allgather)) then
         write(6,*)                                          &
           'SPMD_DYN_SETOPTS:  ERROR:  dyn_allgather=', &
           dyn_allgather_in,                              &
           '  is out of range.  It must be between ',        &
           min_allgather,' and ',max_allgather
         call endrun
       endif
     endif
  end subroutine spmd_dyn_setopts

!========================================================================

   subroutine spmdinit_dyn ()
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute latitudes among available processors
! 
! Method: Distribution is S->N for processors 0->npes
! 
! Author: CCM Core Group
! Modified: P. Worley, November 2003 to improve SMP load balance, and to
!           change distribution to 
!             S->E for processors 0,2,..,npes-2
!           and 
!             N->E for processors 1,3,..,npes-1
!           when mirror flag is set (at request of physics)
! 
!-----------------------------------------------------------------------
      use comspe, only: numm
      use spmd_utils
#if (defined MODCM_DP_TRANSPOSE)
      use parutilitiesmodule, only : parinit
#endif
!-----------------------------------------------------------------------
!
! Local workspace
!
      integer procj     ! process offset loop index
      integer procid    ! process id
      integer procids   ! process id SH
      integer procidn   ! process id NH
      integer smpid     ! SMP id
      integer smpids    ! SMP id for SH process
      integer smpidn    ! SMP id for NH process
      integer nlat_base ! minimum number of latitudes per proc
      integer workleft  ! amount of work still to be parcelled out
      integer lat       ! latitude index
      integer iend      ! ending latitude band of work for a given proc
      integer smostlat  ! southern-most latitude index
      integer nmostlat  ! northern-most latitude index
      integer m2,m3,m5  ! 2, 3, 5 prime factors for problem decomposition
      integer nlat_smp(0:npes-1)  ! number of latitudes per SMP
      integer nproc_smp(0:npes-1) ! number of MPI processes per SMP
      real(r8) avgnlat_proc(0:npes-1) ! average number of latitudes per 
                                      ! MPI process in a given SMP node
      real(r8) minavgnlat_proc        ! minimum average number of latitudes per 
                                      ! MPI process over SMP nodes
      integer neighn_minlat(plat)    ! minimum latitude in north neighbor
      integer neighs_maxlat(plat)    ! maximum latitude in south neighbor
      integer allgather_steps         ! number of swaps in realloc7 algorithm
!
!-----------------------------------------------------------------------
!
! Initialize Pilgrim library
!
#if (defined MODCM_DP_TRANSPOSE)
      call parinit(mpicom)
#endif
!
! Initialize mirror flag
!
      mirror = phys_mirror_decomp_req
!
! Allocate memory for number of lats per proc
!
      allocate (nlat_p (0:npes-1))
      nlat_p(0:npes-1) = 0
!
! Make sure number of PEs and number of latitudes are kosher
!
      call factor (plat, m2, m3, m5)

      if (m2 < 1) then
         call endrun ('SPMDINIT_DYN: Problem size is not divisible by 2')
      end if

      if (masterproc) then
         write (6,*) 'Problem factors: 2**',m2,' * 3**',m3,' * 5**',m5
      end if
      call factor (npes, m2, m3, m5)
      
      if (mod(npes,2) /= 0) then
         write(6,*)'SPMDINIT_DYN: nprocs(',npes,') must be a multiple of 2'
         call endrun
      end if

!
! Determine minimum number of latitudes for each process
!
      nlat_base = plat/npes
      do procids=0,npes/2-1
         procidn = npes - procids - 1
         nlat_p(procids) = nlat_base
         nlat_p(procidn) = nlat_p(procids)
      enddo
      maxlats = nlat_base
!
! Calculate initial distribution of latitudes and 
! distribution of processes by SMP
!
      nlat_smp(0:npes-1) = 0
      nproc_smp(0:npes-1) = 0
      do procid=0,npes-1
         smpid = proc_smp_map(procid)
         nproc_smp(smpid) = nproc_smp(smpid) + 1
      enddo
!
      do smpid=0,nsmps-1
         nlat_smp(smpid)     = nlat_base*nlat_smp(smpid)
         avgnlat_proc(smpid) = float(nlat_base)
      enddo
!
! Equi-distribute remaining latitudes across SMPs
! without increasing per process imbalance beyond minimum
!
      workleft = plat - npes*nlat_base
      if (workleft > 0) maxlats = maxlats + 1
      do while (workleft > 0)
!
! (a) Find minimun number of latitudes assigned to an SMP
!
         minavgnlat_proc = avgnlat_proc(0)
         do smpid=1,nsmps-1
            if (minavgnlat_proc > avgnlat_proc(smpid)) then
               minavgnlat_proc = avgnlat_proc(smpid)
            endif
         enddo
!
! (b) Assign an additional latitude to processes with nlat_base
!     latitudes in SMPs with the minimum average number of 
!     latitudes
!
         do procid=npes/2-1,0,-1
            if (mirror) then
               procids = 2*procid
               procidn = procids + 1
            else
               procids = procid
               procidn = npes - procids - 1
            endif
!
            smpids = proc_smp_map(procids)
            smpidn = proc_smp_map(procidn)
            if ((nlat_p(procids) .eq. nlat_base)  .and. &
                ((avgnlat_proc(smpids) .eq. minavgnlat_proc) .or. &
                 (avgnlat_proc(smpidn) .eq. minavgnlat_proc)) .and. &
                (workleft > 0)) then
!
               nlat_p(procids) = nlat_p(procids) + 1
               nlat_smp(smpids) = nlat_smp(smpids) + 1
               avgnlat_proc(smpids) = &
                  float(nlat_smp(smpids))/float(nproc_smp(smpids))
!
               nlat_p(procidn) = nlat_p(procids)
               nlat_smp(smpidn) = nlat_smp(smpidn) + 1
               avgnlat_proc(smpidn) = &
                  float(nlat_smp(smpidn))/float(nproc_smp(smpidn))
!
               workleft = workleft - 2
            endif
         enddo
      end do
!
! Determine latitude assignments
!
      iend = 0
      do procid=0,npes/2-1
         if (mirror) then
            procids = 2*procid
            procidn = procids + 1
         else
            procids = procid
            procidn = npes - procids - 1
         endif
!
         cut(1,procids) = iend + 1
         cut(2,procids) = iend + nlat_p(procids)
         iend = iend + nlat_p(procids)
!
! Assign mirror latitudes
!
         cut(1,procidn) = plat - cut(2,procids) + 1
         cut(2,procidn) = plat - cut(1,procids) + 1
!
! Save local information
!
         if (iam == procids .or. iam == procidn) then
            beglat = cut(1,iam)
            endlat = cut(2,iam)
            numlats = nlat_p(iam)
            begirow = cut(1,procids)
            endirow = cut(2,procids)
         end if
!
      enddo
!
      do procid=0,npes-1
         if (masterproc) then
            write(6,*)'procid ',procid,' assigned ', &
                      cut(2,procid)-cut(1,procid)+1,' latitude values from', &
                      cut(1,procid),' through ',cut(2,procid)
         end if
!
! Determine which processor is responsible for the defined latitudes
!
         do lat=cut(1,procid),cut(2,procid)
            proc(lat) = procid
         end do
!
! The extended regions are simply "numbnd" wider at each
! side. The extended region do not go beyond 1 and plat, though
!
         cutex(1,procid) = cut(1,procid) - numbnd
         cutex(2,procid) = cut(2,procid) + numbnd
         if (iam == procid) then
            beglatex = cutex(1,procid) + numbnd
            endlatex = cutex(2,procid) + numbnd
            numlatsex = endlatex - beglatex + 1
         end if
      end do
!
! Determine neighbor processors needed for boundary communication.  
! North first.
!
      neighn = 0
      neighn_minlat(:) = -1
      do procid=0,npes-1
         if (procid /= iam) then
            if ((cut(1,procid) > cut(2,iam)) .and. &
                (cut(1,procid) <= cut(2,iam)+numbnd)) then
               neighn_minlat(cut(1,procid)) = procid
               neighn = neighn + 1
            endif
         endif
      enddo
!
! Sort north processes by increasing latitude
!
      allocate (neighn_proc (neighn))
      neighn = 0
      do lat=1,plat
         if (neighn_minlat(lat) /= -1) then
            neighn = neighn + 1
            neighn_proc(neighn) = neighn_minlat(lat)
         endif
      enddo
!
! South next.
!
      neighs = 0
      neighs_maxlat(:) = -1
      do procid=0,npes-1
         if (procid /= iam) then
            if ((cut(2,procid) < cut(1,iam)) .and. &
                (cut(2,procid) >= cut(1,iam)-numbnd)) then
               neighs_maxlat(cut(2,procid)) = procid
               neighs = neighs + 1
            endif
         endif
      enddo
!
! Sort south processes by decreasing latitude
!
      allocate (neighs_proc (neighs))
      neighs = 0
      do lat=plat,1,-1
         if (neighs_maxlat(lat) /= -1) then
            neighs = neighs + 1
            neighs_proc(neighs) = neighs_maxlat(lat)
         endif
      enddo
!
      if (masterproc) then
         write(6,*)'-----------------------------------------'
         write(6,*)'Number of lats passed north & south = ',numbnd
         write(6,*)'Node  Partition  Extended Partition'
         write(6,*)'-----------------------------------------'
         do procid=0,npes-1
            write(6,200) procid,cut(1,procid),cut(2,procid) ,cutex(1,procid), cutex(2,procid)
200         format(i3,4x,i3,'-',i3,7x,i3,'-',i3)
         end do
      end if
!      write(6,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!      write(6,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

      call decomp_wavenumbers ()
      call spmdbuf ()
!
! Precompute swap partners and number of steps in realloc4 alltoall algorithm.
! First, determine number of swaps.
!
      realloc4_steps = 0
      do procj=1,ceil2(npes)-1
         procid = pair(npes,procj,iam)
         if (procid >= 0) then
            if ((numm(iam) > 0 .or. numm(procid) > 0)) then
               realloc4_steps = realloc4_steps + 1
            end if
         end if
      end do
!
! Second, determine swap partners.
!
      allocate( realloc4_proc(realloc4_steps) )
      allocate( realloc4_step(0:npes-1) )
      realloc4_step(:) = -1
      realloc4_steps = 0
      do procj=1,ceil2(npes)-1
         procid = pair(npes,procj,iam)
         if (procid >= 0) then
            if ((numm(iam) > 0 .or. numm(procid) > 0)) then
               realloc4_steps = realloc4_steps + 1
               realloc4_proc(realloc4_steps) = procid
               realloc4_step(procid) = realloc4_steps
            end if
         end if
      end do
!
! Precompute swap partners in realloc5/7 allgather algorithm.
      allocate( allgather_proc(npes-1) )
      allocate( allgather_step(0:npes-1) )
      allgather_step(:) = -1
      allgather_steps = 0
      do procj=1,ceil2(npes)-1
         procid = pair(npes,procj,iam)
         if (procid >= 0) then
            allgather_steps = allgather_steps + 1
            allgather_proc(allgather_steps) = procid
            allgather_step(procid) = allgather_steps
         end if
      end do
!
      return
   end subroutine spmdinit_dyn

!========================================================================

   subroutine factor (nitems, m2, m3, m5)
!----------------------------------------------------------------------- 
! 
! Purpose: Factor a given number into powers of 2,3,5
! 
! Method: Brute force application of "mod" function
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Arguments
!
      integer, intent(in) :: nitems      ! Number to be factored into powers of 2,3,5
      integer, intent(out) :: m2,m3,m5   ! Powers of 2, 3, and 5 respectively
!
! Local workspace
!
      integer num                        ! current number to be factored
!
!-----------------------------------------------------------------------
!
      num = nitems
      m2 = 0
      m3 = 0
      m5 = 0
      
2     if (mod(num,2) == 0) then
         m2 = m2 + 1
         num = num/2
         goto 2
      end if
      
3     if (mod(num,3) == 0) then
         m3 = m3 + 1
         num = num/3
         goto 3
      end if
      
5     if (mod(num,5) == 0) then
         m5 = m5 + 1
         num = num/5
         goto 5
      end if
      
      if (num /= 1) then
         write(6,*) 'FACTOR: ',nitems,' has a prime factor other than 2, 3, or 5.  Aborting...'
         call endrun
      end if
      
      return
   end subroutine factor

!========================================================================

   subroutine decomp_wavenumbers
!----------------------------------------------------------------------- 
! 
! Purpose: partition the spectral work among the given number of processors
! 
! Method: Approximately equidistribute both the number of spectral 
!         coefficients and the number of wavenumbers assigned to each 
!         MPI task using a modified version of the mapping due to
!         Barros and Kauranne.
! 
! Author: P. Worley, September 2002
! 
!-----------------------------------------------------------------------
      use pspect, only: pmmax
      use comspe, only: numm, maxm, locm, nlen
      use infnan, only: bigint
!
! Local workspace
!
      integer procid      ! processor id
      integer m, lm       ! global and local fourier wavenumber indices
      integer speccount(0:npes-1)
                          ! number of spectral coefficients assigned to
                          ! each MPI task
      integer mstride     ! Stride over wavenumbers used in decomposition
      integer begm1       ! Starting Fourier wavenumbers owned by an MPI task
      integer begm2       !  when using Barros & Kauranne decomposition
!-----------------------------------------------------------------------
!
      if (mod(pmmax,npes) .eq. 0) then
         maxm = pmmax/npes
      else
         maxm = (pmmax/npes) + 1
      endif
      allocate ( locm(1:maxm, 0:npes-1) )
!
      mstride = 2*npes
      npessp = 0
      do procid = 0,npes-1
         numm(procid) = 0
         speccount(procid) = 0
         begm1 = procid + 1
         begm2 = mstride - procid
         do m=begm1,pmmax,mstride
            numm(procid) = numm(procid) + 1
            locm(numm(procid),procid) = m
            speccount(procid) = speccount(procid) + nlen(m)
         enddo
         do m=begm2,pmmax,mstride
            numm(procid) = numm(procid) + 1
            locm(numm(procid),procid) = m
            speccount(procid) = speccount(procid) + nlen(m)
         enddo
!
         if (numm(procid) .gt. 0) then
            npessp = npessp + 1
         endif
!
      enddo
!
      do procid = 0,npes-1
         if (masterproc) then
            write(6,*)'procid ',procid,' assigned ', speccount(procid), &
                      ' spectral coefficients and ', numm(procid), &
                      ' m values: ', (locm(lm,procid),lm=1,numm(procid))
         end if
         do lm=numm(procid)+1,maxm
            locm(lm,procid) = bigint
         enddo
      enddo
!   
      return
   end subroutine decomp_wavenumbers

!========================================================================

  subroutine spmdbuf
!----------------------------------------------------------------------- 
! 
! Purpose: allocate spmd pack buffers used in pairwise all-all exchanges
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
     use error_messages, only: alloc_err
     use comspe, only: nlen, maxm
!-----------------------------------------------------------------------
!
! Local workspace
!
     integer maxcount(4),m
     integer length,i,lm,istat
!
! realloc4a max: 4  2 plev*numm*numlats (e.g. tdyn)
!                1  2     *numm*numlats (bpstr)
!
     maxcount(1) = (npes-1)*maxlats*(2*maxm*(plev*4 + 1))
!
! realloc4b max: 11 2 plev*numm*numlats (e.g. vort)
!                4  2     *numm*numlats (e.g. dps)
!
     maxcount(2) = (npes-1)*maxlats*(2*maxm*(plev*11 + 4))
!
! realloc5 max: 6 numlats         (tmass)
!               5 numlats  *pcnst (e.g. hw1lat)
!               2 4*numlats*pcnst (e.g.hw2al)
!
     maxcount(3) = npes*maxlats*(6 + (5 + 2*4)*pcnst)
!
! realloc7 max: 3 plev *numlats    (e.g. vmax2d)
!               4      *numlats    (e.g. psurf)
!
     maxcount(4) = npes*maxlats*(3*plev + 4)
!
     m = maxval(maxcount)
     call mpipack_size (m, mpir8, mpicom, bsiz)
     write(6,*) 'SPMDBUF: Allocating SPMD buffers of size ',bsiz
     spmdbuf_siz = bsiz/8 + 1
     allocate(buf1(spmdbuf_siz), stat=istat)
     call alloc_err( istat, 'spmdbuf', 'buf1', spmdbuf_siz )
     allocate(buf2(spmdbuf_siz), stat=istat)
     call alloc_err( istat, 'spmdbuf', 'buf2', spmdbuf_siz )
     return
  end subroutine spmdbuf

!========================================================================

  subroutine spmdbuf_resize(bufsiz)
!----------------------------------------------------------------------- 
! 
! Purpose: reallocate spmd pack buffers when current size is too small
! 
! Author: P. Worley
! 
!-----------------------------------------------------------------------
     use error_messages, only: alloc_err
!-----------------------------------------------------------------------
!
! Arguments
!
     integer, intent(in) :: bufsiz      ! required buffer size (in r8s)
!-----------------------------------------------------------------------
!
! Local workspace
!
     integer istat
!
!-----------------------------------------------------------------------
!
     if (bufsiz > spmdbuf_siz) then
!
        deallocate(buf1)
        deallocate(buf2)
        call mpipack_size (bufsiz, mpir8, mpicom, bsiz)
        write(6,*) 'SPMDBUF_RESIZE: Allocating SPMD buffers of size ',bsiz
        spmdbuf_siz = bsiz/8 + 1
        allocate(buf1(spmdbuf_siz), stat=istat)
        call alloc_err( istat, 'spmdbuf', 'buf1', spmdbuf_siz )
        allocate(buf2(spmdbuf_siz), stat=istat)
        call alloc_err( istat, 'spmdbuf', 'buf2', spmdbuf_siz )
!
     endif
!
     return
  end subroutine spmdbuf_resize

!========================================================================

  subroutine compute_gsfactors (numperlat, numtot, numperproc, displs)
!----------------------------------------------------------------------- 
! 
! Purpose: Compute arguments for gatherv, scatterv
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
!
! Input arguments
!
     integer, intent(in) :: numperlat    ! number of elements per latitude
!
! Output arguments
!
     integer, intent(out) :: numtot               ! total number of elements (to send or recv)
     integer, intent(out) :: numperproc(0:npes-1) ! per-PE number of items to receive
     integer, intent(out) :: displs(0:npes-1)     ! per-PE displacements
!
! Local variables
!
     integer :: p                    ! index
   
     numtot = numperlat*numlats
   
     do p=0,npes-1
        numperproc(p) = numperlat*nlat_p(p)
     end do
     
     displs(0) = 0
     do p=1,npes-1
        displs(p) = numperlat*(cut(1,p)-1)
     end do
     
  end subroutine compute_gsfactors

#endif

end module spmd_dyn
