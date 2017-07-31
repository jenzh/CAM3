#include <misc.h>
#include <preproc.h>

module lp_coupling

#if (defined COUP_CAM)

!------------------------------------------------------------------------------
!BOP
!
! !MODULE: lp_coupling
!
! !DESCRIPTION:
! Module provides coupling between the atmosphere physics (decomposed into
! chunks) and the land (decomposed into clumps).
!
! !USES:
   use shr_kind_mod, only : r8 => shr_kind_r8

#if (defined SPMD)
   use mpishorthand, only : mpir8, mpicom
   use spmd_dyn    , only : npes
   use pmgrid      , only : iam
#else
   use spmdMod     , only : npes, iam
#endif
   use decompMod   , only : get_nclumps, get_clump_owner_id, &
                            get_clump_ncells_id, get_clump_coord_id, &
                            get_clump_gcell_info
   use phys_grid   , only : get_chunk_coord_owner_p
   use abortutils  , only : endrun
!
! !PUBLIC TYPES:
   implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
   public lp_coupling_init             ! initialize clump<-->chunk mapping
   public lp_coupling_finalize         ! destroy clump<-->chunk mapping
   public alltoall_clump_to_chunk_init ! communicate fluxes from lnd to atm
   public alltoall_clump_to_chunk      ! communicate fluxes from lnd to atm
   public alltoall_chunk_to_clump      ! communicate fluxes from atm to lnd
   SAVE
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!
! !PRIVATE TYPES:
   private

   type clump2chunk
      integer :: lchunk
      integer :: col
   end type clump2chunk
   type(clump2chunk), dimension(:,:), allocatable, private :: clump2chunks

   type chunk2clump
      integer :: clumpid
      integer :: cell
   end type chunk2clump
   type(chunk2clump), dimension(:,:), allocatable, private :: chunk2clumps

   real(r8), dimension(:), target, allocatable :: lp_sendbuf ! lnd->phys send buf
   real(r8), dimension(:), target, allocatable :: lp_recvbuf ! lnd->phys receive buf
   real(r8), dimension(:), target, allocatable :: pl_sendbuf ! phys->lnd send buf
   real(r8), dimension(:), target, allocatable :: pl_recvbuf ! phys->lnd receive buf

   integer, dimension(:), allocatable :: lp_blkcnts ! l->p send/p->l recv blocks
   integer, dimension(:), allocatable :: lp_sndcnts ! lnd->phys send counts
   integer, dimension(:), allocatable :: lp_rcvcnts ! lnd->phys receive counts
   integer, dimension(:), allocatable :: lp_sdispls ! lnd->phys send dsplsmnt
   integer, dimension(:), allocatable :: lp_rdispls ! lnd->phys receive dsplsmnt

   integer, dimension(:), allocatable :: pl_blkcnts ! p->l send/l->p recv blocks
   integer, dimension(:), allocatable :: pl_sndcnts ! phys->lnd send counts
   integer, dimension(:), allocatable :: pl_rcvcnts ! phys->lnd receive counts
   integer, dimension(:), allocatable :: pl_sdispls ! phys->lnd send dsplsmnt
   integer, dimension(:), allocatable :: pl_rdispls ! phys->lnd receive dsplsmnt

   integer, parameter :: pl_nval = 16        ! phys->lnd flux values
   integer, parameter :: lp_nval = 13        ! lnd->phys flux values

   logical :: lpc_init_flag = .false.        ! set if initialized
!------------------------------------------------------------------------------
! lp_coupling.F90,v 1.1.6.11 2004/04/27 04:11:46 forrest Exp
! forrest
!------------------------------------------------------------------------------

   contains

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lp_coupling_init
!
! !INTERFACE:
   subroutine lp_coupling_init()
!
! !DESCRIPTION:
! This subroutine initializes the mapping between the atmosphere physics
! chunks and the land clumps.  It may (and must) be called repeatedly to
! re-initialize the mapping if the decomposition of either the atmosphere
! physics or the land changes.  It allocates communication buffers
! constructs vectors of counts and displacements used for subsequent
! communication between MPI processes.

! !ARGUMENTS:
   implicit none
!
! !LOCAL VARIABLES:
   integer :: p, c, g                            ! loop indices
   integer :: nclumps                            ! number of clumps defined
   integer :: ncells                             ! number of clump cells
   integer :: clump_owner                        ! clump owner
   integer, dimension(:), allocatable :: lons    ! clump longitudes
   integer, dimension(:), allocatable :: lats    ! clump latitudes
   integer, dimension(:), allocatable :: lchnks  ! chunk ids
   integer, dimension(:), allocatable :: cols    ! chunk columns
   integer, dimension(:), allocatable :: chunk_owners  ! chunk owners
   integer :: max_gpc = 0                        ! max cells per clump
   integer :: ier                                ! error codes
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!------------------------------------------------------------------------------
! lp_coupling.F90,v 1.1.6.11 2004/04/27 04:11:46 forrest Exp
! forrest
!------------------------------------------------------------------------------

   ! If already initialized, then deallocate buffers and re-initialize everything

   call lp_coupling_finalize()

   allocate(lp_blkcnts(0:npes-1), lp_sndcnts(0:npes-1), lp_rcvcnts(0:npes-1), &
            pl_blkcnts(0:npes-1), pl_sndcnts(0:npes-1), pl_rcvcnts(0:npes-1), &
            stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for lp_blkcnts, ', &
         'lp_sndcnts, pl_blkcnts, and pl_sndcnts'
      call endrun
   end if
   lp_blkcnts(:) = 0
   lp_sndcnts(:) = 0
   lp_rcvcnts(:) = 0
   pl_blkcnts(:) = 0
   pl_sndcnts(:) = 0
   pl_rcvcnts(:) = 0

   ! Determine max_gpc and allocate dynamic memory

   nclumps = get_nclumps()
   do c = 1,nclumps
      ncells = get_clump_ncells_id(c)
      if (ncells > max_gpc) max_gpc = ncells
   end do
   allocate(lons(max_gpc), lats(max_gpc), lchnks(max_gpc), cols(max_gpc), &
        chunk_owners(max_gpc), stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for local lons, ', &
         'lats, lchnks, cols, and chunk_owners variables'
      call endrun
   end if

   ! Found above already
   ! nclumps = get_nclumps()
   ! lp_blkcnts: for each cam pid, determine the total number of sends that
   ! will be sent from clm
   ! pl_blkcnts: for each clm pid, determine the total number of sends that
   ! will be sent from cam

   do c = 1,nclumps
      clump_owner = get_clump_owner_id(c)
      ncells = get_clump_ncells_id(c)
      call get_clump_coord_id(c, ncells, lons, lats)
      call get_chunk_coord_owner_p(ncells, lons, lats, lchnks, cols, chunk_owners)
      do g = 1,ncells
         if (clump_owner == iam) then
            p = chunk_owners(g)
            lp_blkcnts(p) = lp_blkcnts(p) + 1
         endif
         if (chunk_owners(g) == iam) then
            p = clump_owner
            pl_blkcnts(p) = pl_blkcnts(p) + 1
         endif
      end do
   end do

   allocate(clump2chunks(0:npes-1, 1:maxval(pl_blkcnts)), &
            chunk2clumps(0:npes-1, 1:maxval(lp_blkcnts)), stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for clump2chunks ', &
         'and chunk2clumps'
      call endrun
   end if
   clump2chunks(:,:)%lchunk = 0
   clump2chunks(:,:)%col = 0
   chunk2clumps(:,:)%clumpid = 0
   chunk2clumps(:,:)%cell = 0


   ! Found above already
   ! nclumps = get_nclumps()

   lp_blkcnts(:) = 0
   pl_blkcnts(:) = 0
   do c = 1,nclumps
      clump_owner = get_clump_owner_id(c)
      ncells = get_clump_ncells_id(c)
      call get_clump_coord_id(c, ncells, lons, lats)
      call get_chunk_coord_owner_p(ncells,lons,lats,lchnks,cols,chunk_owners)
      do g = 1,ncells
         if (clump_owner == iam) then
            p = chunk_owners(g)
            lp_blkcnts(p) = lp_blkcnts(p) + 1
            chunk2clumps(p,lp_blkcnts(p))%clumpid=c
            chunk2clumps(p,lp_blkcnts(p))%cell = g
         end if
         if (chunk_owners(g) == iam) then
            p = clump_owner
            pl_blkcnts(p) = pl_blkcnts(p) + 1
            clump2chunks(p,pl_blkcnts(p))%lchunk = lchnks(g)
            clump2chunks(p,pl_blkcnts(p))%col = cols(g)
         end if
      end do
   end do

   deallocate(lons, lats, lchnks, cols, chunk_owners)

   pl_sndcnts(:) = pl_blkcnts(:) * pl_nval
   lp_rcvcnts(:) = pl_blkcnts(:) * lp_nval
   lp_sndcnts(:) = lp_blkcnts(:) * lp_nval
   pl_rcvcnts(:) = lp_blkcnts(:) * pl_nval

   allocate(pl_sendbuf(0:sum(pl_sndcnts)-1), lp_recvbuf(0:sum(lp_rcvcnts)-1), &
      stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for pl_sendbuf and ', &
         'lp_recvbuf'
      call endrun
   end if
   allocate(lp_sendbuf(0:sum(lp_sndcnts)-1), pl_recvbuf(0:sum(pl_rcvcnts)-1), &
      stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for lp_sendbuf and ', &
         'pl_recvbuf'
      call endrun
   end if

   allocate(lp_sdispls(0:npes-1), lp_rdispls(0:npes-1), pl_sdispls(0:npes-1), &
      pl_rdispls(0:npes-1), stat=ier)
   if (ier /= 0) then
      write (6,*) 'lp_coupling_init(): allocation error for lp_sdispls, ', &
         'lp_rdispls, pl_sdispls, and pl_rdispls'
      call endrun
   end if

   lp_sdispls(0) = 0
   lp_rdispls(0) = 0
   pl_sdispls(0) = 0
   pl_rdispls(0) = 0
   do p = 1,npes-1
      lp_sdispls(p) = lp_sdispls(p-1) + lp_sndcnts(p-1)
      lp_rdispls(p) = lp_rdispls(p-1) + lp_rcvcnts(p-1)
      pl_sdispls(p) = pl_sdispls(p-1) + pl_sndcnts(p-1)
      pl_rdispls(p) = pl_rdispls(p-1) + pl_rcvcnts(p-1)
   end do

   lpc_init_flag = .true.

   end subroutine lp_coupling_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: lp_coupling_finalize
!
! !INTERFACE:
   subroutine lp_coupling_finalize()
!
! !ARGUMENTS:
   implicit none
!
! !DESCRIPTION:
! This subroutine destroys the mapping between the atmsphere physics
! chunks and the land clumps if the lpc\_init\_flag flag is set.  It is
! called from lp\_coupling\_init() to ensure memory is recycled when a
! new mapping is to be created.
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!EOP
!------------------------------------------------------------------------------
! lp_coupling.F90,v 1.1.6.11 2004/04/27 04:11:46 forrest Exp
! forrest
!------------------------------------------------------------------------------

   if (lpc_init_flag) then
      deallocate(clump2chunks, chunk2clumps)
      deallocate(lp_sendbuf, lp_recvbuf, pl_sendbuf, pl_recvbuf)
      deallocate(lp_blkcnts, pl_blkcnts)
      deallocate(lp_sndcnts, lp_rcvcnts, pl_sndcnts, pl_rcvcnts)
      deallocate(lp_sdispls, lp_rdispls, pl_sdispls, pl_rdispls)
      lpc_init_flag = .false.
   endif

   end subroutine lp_coupling_finalize

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: alltoall_clump_to_chunk_init
!
! !INTERFACE:
   subroutine alltoall_clump_to_chunk_init (srfflx2d)
!
! !DESCRIPTION:
! This subroutine performs the initial communication from the land model
! to the atmosphere physics (from clumps to chunks) based on the mapping
! constructed in lp\_coupling\_init().
!
! !USES:
   use ppgrid    , only : begchunk, endchunk
   use comsrf    , only : srfflx_parm, snowhland
   use clm_varcon, only : sb
   use clmtype
   use lnd2atmMod, only : lnd2atm
!
! !ARGUMENTS:
   implicit none
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx2d
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
! 2003.05.01  Mariana Vertenstein Updated to l2as data structures
!
!EOP
!
! !LOCAL VARIABLES:
   integer  :: p, n, k, m, g              ! indices
   integer  :: is                         ! send buffer index
   integer  :: lchnk, i                   ! local chunk and column
   integer  :: ier                        ! returned error code
   real(r8), pointer :: lp_rbufp(:)       ! recv buffer pointer
   type(gridcell_type), pointer :: gptr   ! pointer to gridcell derived subtype
!------------------------------------------------------------------------------
! lp_coupling.F90,v 1.1.6.11 2004/04/27 04:11:46 forrest Exp
! forrest
!------------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

   ! Determine gridcell averaged properties to send to atm

   call lnd2atm(init=.true.)

   ! Fill lnd->phys send buffer

   lp_sendbuf(:) = 0.0_r8

!$OMP PARALLEL DO PRIVATE(p,n,g,is)
!CSD$ PARALLEL DO PRIVATE(p,n,g,is)
   do p = 0, npes-1
      do n = 1, lp_blkcnts(p)

         ! Determine clump 1d gridcell index
         call get_clump_gcell_info (chunk2clumps(p,n)%clumpid, chunk2clumps(p,n)%cell, g)

         ! Fill in send buffer
         is = lp_sdispls(p)+(n-1)*lp_nval + 0
         lp_sendbuf(is) = sqrt(sqrt(gptr%l2af%eflx_lwrad_out(g)/sb)) ! tsxy
         is = lp_sdispls(p)+(n-1)*lp_nval + 1
         lp_sendbuf(is) = gptr%l2as%albd(g,1)                ! asdir
         is = lp_sdispls(p)+(n-1)*lp_nval + 2
         lp_sendbuf(is) = gptr%l2as%albd(g,2)                ! aldir
         is = lp_sdispls(p)+(n-1)*lp_nval + 3
         lp_sendbuf(is) = gptr%l2as%albi(g,1)                ! asdif
         is = lp_sdispls(p)+(n-1)*lp_nval + 4
         lp_sendbuf(is) = gptr%l2as%albi(g,2)                ! aldif
         is = lp_sdispls(p)+(n-1)*lp_nval + 5
         lp_sendbuf(is) = gptr%l2as%h2osno(g)                ! snow [mm]->[m]
         is = lp_sdispls(p)+(n-1)*lp_nval + 6
         lp_sendbuf(is) = 1.e36                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 7
         lp_sendbuf(is) = 1.e36                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 8
         lp_sendbuf(is) = 1.e36                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 9
         lp_sendbuf(is) = 1.e36                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 10
         lp_sendbuf(is) = gptr%l2af%eflx_lwrad_out(g)        ! lwrad
         is = lp_sdispls(p)+(n-1)*lp_nval + 11
         lp_sendbuf(is) = 1.e36                              ! spval
         is = lp_sdispls(p)+(n-1)*lp_nval + 12
         lp_sendbuf(is) = 1.e36                              ! spval

      end do
   end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

#ifdef SPMD
   call mpi_alltoallv (lp_sendbuf, lp_sndcnts, lp_sdispls, mpir8, &
                       lp_recvbuf, lp_rcvcnts, lp_rdispls, mpir8, mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'alltoall_clump_to_chunk_init(): MPI error ', ier, &
         ' from mpi_alltoallv()'
      call endrun
   endif
   lp_rbufp => lp_recvbuf
#else
   lp_rbufp => lp_sendbuf
#endif

   ! Extract lnd->phys receive buffer

!$OMP PARALLEL DO PRIVATE(p,n,lchnk,i)
!CSD$ PARALLEL DO PRIVATE(p,n,lchnk,i)
   do p = 0,npes-1
      do n = 1,pl_blkcnts(p)
         lchnk = clump2chunks(p,n)%lchunk
         i     = clump2chunks(p,n)%col
         srfflx2d(lchnk)%ts(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+0)
         srfflx2d(lchnk)%asdir(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+1)
         srfflx2d(lchnk)%aldir(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+2)
         srfflx2d(lchnk)%asdif(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+3)
         srfflx2d(lchnk)%aldif(i) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+4)
         snowhland(i,lchnk)       = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+5)
         srfflx2d(lchnk)%lwup(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+10)
      end do
   end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

   end subroutine alltoall_clump_to_chunk_init

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: alltoall_clump_to_chunk
!
! !INTERFACE:
   subroutine alltoall_clump_to_chunk (srfflx2d)
!
! !DESCRIPTION:
! This subroutine performs the communication from the land model to
! the atmosphere physics (from clumps to chunks) based on the mapping
! constructed in lp\_coupling\_init().
!
! !USES:
   use ppgrid      , only : begchunk, endchunk
   use comsrf      , only : srfflx_parm, snowhland
   use constituents, only : pcnst, pnats
   use clmtype
   use lnd2atmMod  , only : lnd2atm
!
! !ARGUMENTS:
   implicit none
   type(srfflx_parm), intent(inout), dimension(begchunk:endchunk) :: srfflx2d
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
! 2003.05.01  Mariana Vertenstein Updated to l2as data structures
!
!EOP
!
! !LOCAL VARIABLES:
   integer  :: p, n, k, m, g              ! loop indices
   integer  :: bpatch, npatch             ! patch id and count
   integer  :: lchnk, i                   ! local chunk and column
   integer  :: ier                        ! returned error code
   integer  :: is                         ! send buffer index
   real(r8) :: wt                         ! patch wt
   real(r8), pointer :: lp_rbufp(:)       ! recv buffer pointer
   type(gridcell_type), pointer :: gptr   ! pointer to gridcell derived subtype
!------------------------------------------------------------------------------
! lp_coupling.F90,v 1.1.6.11 2004/04/27 04:11:46 forrest Exp
! forrest
!------------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

   ! Determine gridcell averaged properties to send to atm

   call lnd2atm()

   ! Fill lnd->phys send buffer

   lp_sendbuf(:) = 0.0_r8

!$OMP PARALLEL DO PRIVATE(p,n,g,is)
!CSD$ PARALLEL DO PRIVATE(p,n,g,is)
   do p = 0,npes-1
      do n = 1,lp_blkcnts(p)

         ! Determine clump gridcell info
         call get_clump_gcell_info(chunk2clumps(p,n)%clumpid, chunk2clumps(p,n)%cell, g)

         ! Fill in send buffer
         is = lp_sdispls(p)+(n-1)*lp_nval+0   ! tsxy
         lp_sendbuf(is) = gptr%l2as%t_rad(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+5   ! snow [mm]->[m]
         lp_sendbuf(is) = gptr%l2as%h2osno(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+1   ! asdir
         lp_sendbuf(is) = gptr%l2as%albd(g,1)
         is = lp_sdispls(p)+(n-1)*lp_nval+2   ! aldir
         lp_sendbuf(is) = gptr%l2as%albd(g,2)
         is = lp_sdispls(p)+(n-1)*lp_nval+3   ! asdif
         lp_sendbuf(is) = gptr%l2as%albi(g,1)
         is = lp_sdispls(p)+(n-1)*lp_nval+4   ! aldif
         lp_sendbuf(is) = gptr%l2as%albi(g,2)
         is = lp_sdispls(p)+(n-1)*lp_nval+6   ! taux
         lp_sendbuf(is) = gptr%l2af%taux(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+7   ! tauy
         lp_sendbuf(is) = gptr%l2af%tauy(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+8   ! lhflx
         lp_sendbuf(is) = gptr%l2af%eflx_lh_tot(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+9   ! shflx
         lp_sendbuf(is) = gptr%l2af%eflx_sh_tot(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+10  ! lwrad
         lp_sendbuf(is) = gptr%l2af%eflx_lwrad_out(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+11  ! qflx
         lp_sendbuf(is) = gptr%l2af%qflx_evap_tot(g)
         is = lp_sdispls(p)+(n-1)*lp_nval+12  ! tref
         lp_sendbuf(is) = gptr%l2as%t_ref2m(g)

      end do   ! end loop over destination processes
   end do   ! end loop over source processes
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

#ifdef SPMD
   call mpi_alltoallv (lp_sendbuf, lp_sndcnts, lp_sdispls, mpir8, &
                       lp_recvbuf, lp_rcvcnts, lp_rdispls, mpir8, mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'alltoall_clump_to_chunk(): MPI error ', ier, &
         ' from mpi_alltoallv()'
      call endrun
   endif
   lp_rbufp => lp_recvbuf
#else
   lp_rbufp => lp_sendbuf
#endif

   ! Extract lnd->phys receive buffer

!$OMP PARALLEL DO PRIVATE(p,n,lchnk,i,m)
!CSD$ PARALLEL DO PRIVATE(p,n,lchnk,i,m)
   do p = 0,npes-1
      do n = 1,pl_blkcnts(p)
         lchnk = clump2chunks(p,n)%lchunk
         i     = clump2chunks(p,n)%col
         srfflx2d(lchnk)%ts(i)     = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+0)
         srfflx2d(lchnk)%asdir(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+1)
         srfflx2d(lchnk)%aldir(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+2)
         srfflx2d(lchnk)%asdif(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+3)
         srfflx2d(lchnk)%aldif(i)  = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+4)
         snowhland(i,lchnk)        = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+5)
         srfflx2d(lchnk)%wsx(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+6)
         srfflx2d(lchnk)%wsy(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+7)
         srfflx2d(lchnk)%lhf(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+8)
         srfflx2d(lchnk)%shf(i)    = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+9)
         srfflx2d(lchnk)%lwup(i)   = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+10)
         srfflx2d(lchnk)%cflx(i,1) = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+11)
         srfflx2d(lchnk)%tref(i)   = lp_rbufp(lp_rdispls(p)+(n-1)*lp_nval+12)

         ! Reset all other constituent surface fluxes to zero over land
         do m = 2,pcnst+pnats
            srfflx2d(lchnk)%cflx(i,m) = 0.0_r8
         end do
      end do
   end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

   end subroutine alltoall_clump_to_chunk

!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: alltoall_chunk_to_clump
!
! !INTERFACE:
   subroutine alltoall_chunk_to_clump (srf_state)
!
! !DESCRIPTION:
! This subroutine performs communication from the atmosphere physics to
! the land model (from chunks to clumps) based on the mapping constructed
! in lp\_coupling\_init().
!
! !USES:
   use ppgrid    , only: begchunk, endchunk
   use comsrf    , only: surface_state
   use clmtype
   use clm_varcon, only: rair, po2, pco2
!
! !ARGUMENTS:
   implicit none
   type(surface_state), intent(in), dimension(begchunk:endchunk) :: srf_state
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
! 2002.11.18  Mariana Vertenstein and Forrest Hoffman Updated to clm2.1.
!
!
!EOP
!
! !LOCAL VARIABLES:
   integer  :: g, p, n, k               ! indices
   integer  :: lchnk, i                 ! local chunk and column
   integer  :: ier                      ! returned error code
   real(r8) :: forc_rainc, forc_rainl   ! rainxy [mm/s]
   real(r8) :: forc_snowc, forc_snowl   ! snowfxy [mm/s]
   real(r8), pointer :: pl_rbufp(:)     ! recv buffer pointer
   type(gridcell_type), pointer :: gptr ! pointer to gridcell derived subtype
!------------------------------------------------------------------------------
! lp_coupling.F90,v 1.1.6.11 2004/04/27 04:11:46 forrest Exp
! forrest
!------------------------------------------------------------------------------

    ! Set pointers into derived type

    gptr => clm3%g

   ! Fill phys->lnd send buffer

!$OMP PARALLEL DO PRIVATE(p,n,lchnk,i)
!CSD$ PARALLEL DO PRIVATE(p,n,lchnk,i)
   do p = 0,npes-1
      do n = 1,pl_blkcnts(p)
         lchnk = clump2chunks(p,n)%lchunk
         i = clump2chunks(p,n)%col

         ! Atmoshperic state variable [m]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+0) = srf_state(lchnk)%zbot(i)

         ! Atmoshperic state variable [m/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+1) = srf_state(lchnk)%ubot(i)

         ! Atmoshperic state variable [m/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+2) = srf_state(lchnk)%vbot(i)

         ! Atmoshperic state variable [K]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+3) = srf_state(lchnk)%thbot(i)

         ! Atmoshperic state variable [kg/kg]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+4) = srf_state(lchnk)%qbot(i,1)

         ! Atmoshperic state variable [Pa]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+5) = srf_state(lchnk)%pbot(i)

         ! Atmoshperic state variable [K]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+6) = srf_state(lchnk)%tbot(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+7) = srf_state(lchnk)%flwds(i)

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+8) = &
              srf_state(lchnk)%precsc(i) * 1000.

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+9) = &
            srf_state(lchnk)%precsl(i) * 1000.

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+10) = &
            (srf_state(lchnk)%precc(i) - srf_state(lchnk)%precsc(i)) * 1000.

         ! Atmoshperic flux [m/s] -> [mm/s]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+11) = &
            (srf_state(lchnk)%precl(i) - srf_state(lchnk)%precsl(i)) * 1000.

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+12) = srf_state(lchnk)%soll(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+13) = srf_state(lchnk)%sols(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+14) = srf_state(lchnk)%solld(i)

         ! Atmoshperic flux [W/m^2]
         pl_sendbuf(pl_sdispls(p)+(n-1)*pl_nval+15) = srf_state(lchnk)%solsd(i)

      end do
   end do
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

#ifdef SPMD
   call mpi_alltoallv (pl_sendbuf, pl_sndcnts, pl_sdispls, mpir8, &
                       pl_recvbuf, pl_rcvcnts, pl_rdispls, mpir8, mpicom, ier)
   if (ier /= 0) then
      write (6,*) 'alltoall_chunk_to_clump(): MPI error ', ier, &
         ' from mpi_alltoallv()'
      call endrun
   endif
   pl_rbufp => pl_recvbuf
#else
   pl_rbufp => pl_sendbuf
#endif

   ! Extract phys->lnd receive buffer

!$OMP PARALLEL DO PRIVATE(p,n,g,forc_rainc,forc_rainl,forc_snowc,forc_snowl)
!CSD$ PARALLEL DO PRIVATE(p,n,g,forc_rainc,forc_rainl,forc_snowc,forc_snowl)
   do p = 0,npes-1
      do n = 1,lp_blkcnts(p)

         ! Determine clump gridcell info
         call get_clump_gcell_info(chunk2clumps(p,n)%clumpid, chunk2clumps(p,n)%cell, g)

         ! Fill in clmtype gridcell info
         gptr%a2ls%forc_hgt(g)     = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+0)
         gptr%a2ls%forc_u(g)       = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+1)
         gptr%a2ls%forc_v(g)       = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+2)
         gptr%a2ls%forc_th(g)      = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+3)
         gptr%a2ls%forc_q(g)       = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+4)
         gptr%a2ls%forc_pbot(g)    = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+5)
         gptr%a2ls%forc_t(g)       = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+6)
         gptr%a2lf%forc_lwrad(g)   = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+7)
         forc_snowc                = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+8)
         forc_snowl                = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+9)
         forc_rainc                = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+10)
         forc_rainl                = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+11)
         gptr%a2lf%forc_solad(g,2) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+12)
         gptr%a2lf%forc_solad(g,1) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+13)
         gptr%a2lf%forc_solai(g,2) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+14)
         gptr%a2lf%forc_solai(g,1) = pl_rbufp(pl_rdispls(p)+(n-1)*pl_nval+15)

         ! Determine derived quantities
         gptr%a2ls%forc_hgt_u(g) = gptr%a2ls%forc_hgt(g)  ! obs height of wind [m]
         gptr%a2ls%forc_hgt_t(g) = gptr%a2ls%forc_hgt(g)  ! obs height of temperature [m]
         gptr%a2ls%forc_hgt_q(g) = gptr%a2ls%forc_hgt(g)  ! obs height of humidity [m]
         gptr%a2ls%forc_vp(g)    = gptr%a2ls%forc_q(g) * gptr%a2ls%forc_pbot(g) &
                                   / (0.622 + 0.378 * gptr%a2ls%forc_q(g))
         gptr%a2ls%forc_rho(g)   = (gptr%a2ls%forc_pbot(g) - 0.378 * gptr%a2ls%forc_vp(g)) &
                                   / (rair * gptr%a2ls%forc_t(g))
         gptr%a2ls%forc_co2(g)   = pco2*gptr%a2ls%forc_pbot(g)
         gptr%a2ls%forc_o2(g)    = po2*gptr%a2ls%forc_pbot(g)
         gptr%a2ls%forc_wind(g)  = sqrt(gptr%a2ls%forc_u(g)**2 + gptr%a2ls%forc_v(g)**2)
         gptr%a2lf%forc_solar(g) = gptr%a2lf%forc_solad(g,1) + gptr%a2lf%forc_solai(g,1) + &
                                   gptr%a2lf%forc_solad(g,2) + gptr%a2lf%forc_solai(g,2)

         ! Determine precipitation needed by clm
#ifdef PERGRO
         ! For error growth only, allow rain, not snowfall
         forc_rainc           = forc_rainc + forc_snowc
         forc_rainl           = forc_rainl + forc_snowl
         forc_snowc           = 0.0_r8
         forc_snowl           = 0.0_r8
#endif
         gptr%a2lf%forc_rain(g) = forc_rainc + forc_rainl
         gptr%a2lf%forc_snow(g) = forc_snowc + forc_snowl

      end do   ! end loop over destination processes
   end do   ! end loop over source processes
!CSD$ END PARALLEL DO
!$OMP END PARALLEL DO

   end subroutine alltoall_chunk_to_clump

#endif

end module lp_coupling
