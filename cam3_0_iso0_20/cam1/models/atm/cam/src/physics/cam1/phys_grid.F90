#include <misc.h>
module phys_grid
!----------------------------------------------------------------------- 
! 
! Purpose: Definition of physics computational horizontal grid.
!
! Method: Variables are private; interface routines used to extract
!         information for use in user code.
! 
! Entry points:
!      phy_grid_init       initialize chunk'ed data structure
!
!      phys_grid_defaultopts   get default runtime options
!      phys_grid_setopts       set runtime options
!
!      get_chunk_indices_p get local chunk index range
!      get_ncols_p         get number of columns for a given chunk
!      get_xxx_all_p       get global indices or coordinates for a given
!                          chunk
!      get_xxx_vec_p       get global indices or coordinates for a subset
!                          of the columns in a chunk
!      get_xxx_p           get global indices or coordinates for a single
!                          column
!      where xxx is
!       lat                for global latitude index
!       lon                for global longitude index
!       rlat               for latitude coordinate (in radians)
!       rlon               for longitude coordinate (in radians)
!
!      get_chunk_owner_p   get owner of chunk
!                          for given (lon,lat) coordinate
!
!      get_chunk_coord_p   get local chunk and column indices
!                          for given (lon,lat) coordinates
!
!      get_chunk_coord_owner_p
!                          get owner, local chunk, and column indices
!                          for given (lon,lat) coordinates
!
!      scatter_field_to_chunk
!                          distribute longitude/latitude field
!                          to decomposed chunk data structure
!      gather_chunk_to_field
!                          reconstruct longitude/latitude field
!                          from decomposed chunk data structure
!
!      read_chunk_from_field
!                          read and distribute longitude/latitude field
!                          to decomposed chunk data structure
!      write_field_from_chunk
!                          write longitude/latitude field
!                          from decomposed chunk data structure
!
!      block_to_chunk_send_pters
!                          return pointers into send buffer where data
!                          from decomposed longitude/latitude fields should
!                          be copied to
!      block_to_chunk_recv_pters
!                          return pointers into receive buffer where data
!                          for decomposed chunk data structures should
!                          be copied from
!      transpose_block_to_chunk
!                          transpose buffer containing decomposed 
!                          longitude/latitude fields to buffer
!                          containing decomposed chunk data structures
!
!      chunk_to_block_send_pters
!                          return pointers into send buffer where data
!                          from decomposed chunk data structures should
!                          be copied to
!      chunk_to_block_recv_pters
!                          return pointers into receive buffer where data
!                          for decomposed longitude/latitude fields should
!                          be copied from
!      transpose_chunk_to_block
!                          transpose buffer containing decomposed
!                          chunk data structures to buffer
!                          containing decomposed longitude/latitude fields
!
!      chunk_index         identify whether index is for a latitude or
!                          a chunk
!
! Author: Patrick Worley and John Drake
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8, r4 => shr_kind_r4
   use ppgrid, only: pcols, pver, begchunk, endchunk
   use pmgrid, only: plon, plat, beglat, endlat
   use abortutils, only: endrun
#if ( defined SPMD )
   use spmd_dyn, only: proc, npes, nsmps, proc_smp_map
   use mpishorthand
#endif

   implicit none

   save

#if ( ! defined SPMD )
   integer :: npes = 1
   integer :: nsmps = 1
   integer :: proc_smp_map(0:0)
#endif
   integer, private :: dp_coup_steps     ! number of swaps in transpose algorithm
   integer, dimension(:), private, allocatable :: dp_coup_proc
                                         ! swap partner in each step of 
                                         !  transpose algorithm

! chunk data structures
   type chunk
     integer  :: ncols                 ! number of vertical columns
     integer  :: lon(pcols)            ! global longitude indices
     integer  :: lat(pcols)            ! global latitude indices
     integer  :: owner                 ! id of process where chunk assigned
     integer  :: lchunk                ! local chunk index
   end type chunk

   integer :: nchunks                  ! global chunk count
   type (chunk), dimension(:), allocatable, public :: chunks  
                                       ! global computational grid

   integer, private :: nlchunks        ! local chunk count
   integer, dimension(:), allocatable, private :: lchunks 
                                       ! local chunks

   type knuhc
     integer  :: chunkid               ! chunk id
     integer  :: col                   ! column index in chunk
   end type knuhc

   type (knuhc), dimension(:,:), allocatable, public :: knuhcs
                                       ! map from global (lon,lat) coordinates
                                       ! to chunk'ed grid

! column mapping data structures
   type column_map
     integer  :: chunk                 ! global chunk index
     integer  :: ccol                  ! column ordering in chunk
   end type column_map

   integer :: ngcols                   ! global column count
   integer :: nlcols                   ! local column count
   type (column_map), dimension(:), allocatable, private :: pgcols
                                       ! ordered list of columns (for use in gather/scatter)
                                       ! NOTE: consistent with local ordering

! column remap data structures
   integer, dimension(:), allocatable, private :: gs_col_num
                                       ! number of columns scattered to each process in
                                       ! field_to_chunk scatter
   integer, dimension(:), allocatable, private :: gs_col_offset
                                       ! offset of columns (-1) in pgcols scattered to
                                       ! each process in field_to_chunk scatter

   integer, dimension(:), allocatable, private :: btofc_blk_num
                                       ! number of grid points scattered to each process in
                                       ! block_to_chunk alltoallv, and gathered from each
                                       ! process in chunk_to_block alltoallv

   integer, dimension(:), allocatable, private :: btofc_chk_num
                                       ! number of grid points gathered from each process in
                                       ! block_to_chunk alltoallv, and scattered to each
                                       ! process in chunk_to_block alltoallv

   type btofc_pters
     integer :: ncols                  ! number of columns in block
     integer :: nlvls                  ! number of levels in columns
     integer, dimension(:,:), pointer :: pter 
   end type btofc_pters
   type (btofc_pters), dimension(:), allocatable, private :: btofc_blk_offset
                                       ! offset in btoc send array (-1) where 
                                       ! (blockid, bcid, k) column should be packed in
                                       ! block_to_chunk alltoallv, AND
                                       ! offset in ctob receive array (-1) from which
                                       ! (blockid, bcid, k) column should be unpacked in
                                       ! chunk_to_block alltoallv

   type (btofc_pters), dimension(:), allocatable, private :: btofc_chk_offset
                                       ! offset in btoc receive array (-1) from which
                                       ! (lchnk, i, k) data should be unpacked in
                                       ! block_to_chunk alltoallv, AND
                                       ! offset in ctob send array (-1) where
                                       ! (lchnk, i, k) data should be packed in
                                       ! chunk_to_block alltoallv

   integer :: block_buf_nrecs          ! number of local grid points (lon,lat,lev)
                                       ! in dynamics decomposition (including level 0)
   integer :: chunk_buf_nrecs          ! number of local grid points (lon,lat,lev)
                                       ! in physics decomposition (including level 0)

! miscellaneous phys_grid data
   real(r8) :: clat_p(plat)            ! physics grid latitudes (radians)
   integer  :: nlon_p(plat)            ! num longitudes per latitude
   real(r8) :: clon_p(plon,plat)       ! physics grid longitudes (radians)
   logical :: physgrid_set = .false.   ! flag indicates physics grid has been set
   logical :: local_dp_map = .false.   ! flag indicates that mapping between dynamics 
                                       ! and physics decompositions does not require 
                                       ! interprocessor communication

! Physics grid decomposition options:  
! -1: each chunk is a latitude line
!  0: chunk definitions and assignments do not require interprocess comm.
!  1: chunk definitions and assignments do not require internode comm.
!  2: optimal diurnal, seasonal, and latitude load-balanced chunk definition and assignments
!  3: chunk definitions and assignments only require communication with one other process
!  4: concatenated blocks, no load balancing, no interprocessor communication
   integer, private, parameter :: min_lbal_opt = -1
   integer, private, parameter :: max_lbal_opt = 4
   integer, private, parameter :: def_lbal_opt = 0                ! default
   integer, private :: lbal_opt = def_lbal_opt

! target number of chunks per thread
   integer, private, parameter :: min_chunks_per_thread = 1
   integer, private, parameter :: def_chunks_per_thread = &
                                    min_chunks_per_thread         ! default
   integer, private :: chunks_per_thread = def_chunks_per_thread

! Dynamics/physics transpose method for nonlocal load-balance:
!  0: use mpi_alltoallv
!  1: use point-to-point implementation
!  11-14: use mod_comm, choosing any of several methods internal to mod_comm.
!      The method within mod_comm (denoted mod_method) has possible values 0,1,2,3 and
!      is set according to mod_method = phys_alltoall - modmin_alltoall, where
!      modmin_alltoall is 11.
   integer, private, parameter :: min_alltoall = 0
# if defined(MODCM_DP_TRANSPOSE)
   integer, private, parameter :: max_alltoall = 14
   integer, private, parameter :: modmin_alltoall = 11
# else
   integer, private, parameter :: max_alltoall = 1
# endif
   integer, private, parameter :: def_alltoall = 0                ! default
   integer, private :: phys_alltoall = def_alltoall

contains
!========================================================================

   subroutine phys_grid_init( )
!----------------------------------------------------------------------- 
! 
! Purpose: Physics mapping initialization routine:  
! 
! Method: 
! 
! Author: John Drake and Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, plev, plond, platd, masterproc
   use pspect, only: pmmax, pnmax
   use rgrid, only: nlon
   use commap, only: clat, clon
   use dyn_grid, only: get_block_coord_cnt_d, get_block_coord_d, &
                       get_block_col_cnt_d, get_block_lvl_cnt_d, &
                       get_lon_d, get_lat_d, get_block_bounds_d, &
                       get_block_owner_d, get_block_levels_d
   use spmd_utils, only: pair, ceil2
!
!------------------------------Arguments--------------------------------
!
!
!---------------------------Local workspace-----------------------------
!
   integer :: i, j, jb, k, lchnk, p      ! loop indices
   integer :: tchunks                    ! target number of chunks per thread
   integer :: cid                        ! chunk id
   integer :: pchunkid                   ! chunk global ordering
   integer :: begpchunk, endpchunk       ! segment of chunk global ordering on 
                                         !  a given process
   integer :: plchunks                   ! number of chunks for a given process
   integer :: curgcol                    ! current global column index
   integer :: firstblock, lastblock      ! global block indices
   integer :: blksiz                     ! current block size
   integer :: glbcnt, curcnt             ! running grid point counts
   integer :: curp                       ! current process id
   integer :: block_cnt                  ! number of blocks containing data
                                         ! for a given vertical column
   integer :: numlvl                     ! number of vertical levels in block 
                                         ! column
   integer :: levels(plev+1)             ! vertical level indices
   integer :: owner_d                    ! processor owning given block column
   integer :: owner_p                    ! processor owning given chunk column
   integer :: ncol                       ! number of columns in current chunk
   integer :: blockids(plev+1)           ! block indices
   integer :: bcids(plev+1)              ! block column indices
   integer :: glon, glat                 ! global (lon,lat) indices
   integer :: ntmp1, ntmp2               ! work variables
!-----------------------------------------------------------------------
!
! Initialize physics grid, using dynamics grid
!
   do j=1,plat
      clat_p(j) = clat(j)
      nlon_p(j) = nlon(j)
      do i=1,nlon(j)
         clon_p(i,j) = clon(i,j)
      enddo
   enddo
!
! Determine total number of columns and block index bounds
!
   ngcols = 0
   do j=1,plat
      ngcols = ngcols + nlon_p(j)
   enddo
   call get_block_bounds_d(firstblock,lastblock)
!
! Option -1: each latitude line is a single chunk, same as 1D dynamics decompositions.
!            
   if (lbal_opt == -1) then
!
! Check that pcols == plon
!
      if (pcols /= plon) then
         call endrun ('PHYS_GRID_INIT error: phys_loadbalance -1 specified but PCOLS /= PLON')
      endif
!
! Determine total number of chunks
!
      nchunks = plat
!
! Allocate and initialize chunks and knuhcs data structures
!
      allocate ( chunks(1:nchunks) )
      allocate ( knuhcs(1:plond, 1:platd) )

      cid = 0
      do j=1,plat
         chunks(j)%ncols = nlon_p(j)
         do i=1,chunks(j)%ncols
            chunks(j)%lon(i) = i
            chunks(j)%lat(i) = j
            knuhcs(i,j)%chunkid = j
            knuhcs(i,j)%col = i
         enddo
      enddo
!
! Determine parallel decomposition (assuming 1D latitude decomposition in dynamics)
!
      do j=1,plat
#if (defined SPMD)
         chunks(j)%owner = proc(j)
#else
         chunks(j)%owner = 0
#endif
      enddo
!
! (including allocating and initializing data structures for gather/scatter)
!  
      allocate ( pgcols(1:ngcols) )
      allocate ( gs_col_num(0:npes-1) )
      allocate ( gs_col_offset(0:npes) )

      pchunkid = 0
      endpchunk = 0
      curgcol = 0
      do p=0,npes-1
         gs_col_offset(p) = curgcol + 1
         begpchunk = endpchunk + 1
         plchunks = 0
         gs_col_num(p) = 0
         do cid=1,nchunks
            if (chunks(cid)%owner == p) then
               pchunkid = pchunkid + 1
               plchunks = plchunks + 1

               do i=1,chunks(cid)%ncols
                  curgcol = curgcol + 1
                  pgcols(curgcol)%chunk = cid
                  pgcols(curgcol)%ccol = i
                  gs_col_num(p) = gs_col_num(p) + 1
               enddo

            endif
         enddo
         endpchunk = begpchunk + plchunks - 1
      enddo
      gs_col_offset(npes) = curgcol + 1

      do j=1,plat
         chunks(j)%lchunk = j
      enddo
      nlchunks = endlat-beglat+1
      nlcols = gs_col_num(iam)
!
! Local chunk indices are identical to global latitudes {beglat,...,endlat}
!
      begchunk = beglat
      endchunk = endlat
      allocate ( lchunks(begchunk:endchunk) )
      do j=begchunk,endchunk
         lchunks(j) = j
      enddo
!
! Set flag indicating columns in physics and dynamics 
! decompositions reside on the same processors
!
      local_dp_map = .true. 
!
   else
!
! Option == 0: split local longitude/latitude blocks into chunks,
!               while attempting to create load-balanced chunks.
!               Does not work with vertically decomposed blocks.
!               (default)
! Option == 1: split SMP-local longitude/latitude blocks into chunks,
!               while attempting to create load-balanced chunks.
!               Does not work with vertically decomposed blocks.
! Option == 2: load balance chunks with respect to diurnal and
!               seaonsal cycles and wth respect to latitude, 
!               and assign chunks to processor 
!               in a way that attempts to minimize communication costs
! Option == 3: divide processes into pairs and split 
!               longitude/latitude blocks assigned to these pairs into 
!               chunks, attempting to create load-balanced chunks.
!               The process pairs are chosen to maximize load balancing
!               opportunities.
!               Does not work with vertically decomposed blocks.
! Option == 4: concatenate local longitude/latitude blocks, then
!               divide into chunks.
!               Does not work with vertically decomposed blocks.
! Option == 5: split indiviudal longitude/latitude blocks into chunks,
!               assigning columns using block ordering
!
! Allocate and initialize chunks and knuhcs data structures, then
! assign chunks to processes.
!
      call create_chunks(lbal_opt, chunks_per_thread)
!
! Determine whether dynamics and physics decompositions
! are colocated, not requiring any interprocessor communication
! in the coupling.
      local_dp_map = .true.   
      do cid=1,nchunks
         do i=1,chunks(cid)%ncols
            glon = chunks(cid)%lon(i)
            glat = chunks(cid)%lat(i) 
            block_cnt = get_block_coord_cnt_d(glon,glat)
            call get_block_coord_d(glon,glat,block_cnt,blockids,bcids)
            do jb=1,block_cnt
               owner_d = get_block_owner_d(blockids(jb)) 
               if (owner_d .ne. chunks(cid)%owner) then
                  local_dp_map = .false.   
               endif
            enddo
         enddo
      enddo
!
! Allocate and initialize data structures for gather/scatter
!  
      allocate ( pgcols(1:ngcols) )
      allocate ( gs_col_num(0:npes-1) )
      allocate ( gs_col_offset(0:npes) )

      pchunkid = 0
      endpchunk = 0
      curgcol = 0
      do p=0,npes-1
         gs_col_offset(p) = curgcol + 1
         begpchunk = endpchunk + 1
         plchunks = 0
         gs_col_num(p) = 0
         do cid=1,nchunks
            if (chunks(cid)%owner == p) then
               pchunkid = pchunkid + 1

               plchunks = plchunks + 1
               chunks(cid)%lchunk = pchunkid + lastblock

               do i=1,chunks(cid)%ncols
                  curgcol = curgcol + 1
                  pgcols(curgcol)%chunk = cid
                  pgcols(curgcol)%ccol = i
                  gs_col_num(p) = gs_col_num(p) + 1
               enddo

            endif
         enddo
         endpchunk = begpchunk + plchunks - 1
         if (iam == p) then
!
! Local chunk index range chosen so that it does not overlap 
! {begblock,...,endblock}
! 
            nlchunks = plchunks
            begchunk = begpchunk + lastblock
            endchunk = endpchunk + lastblock
         endif
      enddo
      gs_col_offset(npes) = curgcol + 1
      nlcols = gs_col_num(iam)
!
      allocate ( lchunks(begchunk:endchunk) )
      do cid=1,nchunks
         if (chunks(cid)%owner == iam) then
            lchunks(chunks(cid)%lchunk) = cid
         endif
      enddo
!
   endif
!
   if (.not. local_dp_map) then
!
! allocate and initialize data structures for transposes
!  
      allocate ( btofc_blk_num(0:npes-1) )
      allocate ( btofc_blk_offset(firstblock:lastblock) )
      do jb = firstblock,lastblock
         nullify( btofc_blk_offset(jb)%pter )
      enddo
!
      glbcnt = 0
      curcnt = 0
      curp = 0
      do curgcol=1,ngcols
         cid = pgcols(curgcol)%chunk
         i   = pgcols(curgcol)%ccol
         owner_p   = chunks(cid)%owner
         do while (curp < owner_p)
            btofc_blk_num(curp) = curcnt
            curcnt = 0
            curp = curp + 1
         enddo
         glon = chunks(cid)%lon(i)
         glat = chunks(cid)%lat(i)
         block_cnt = get_block_coord_cnt_d(glon,glat)
         call get_block_coord_d(glon,glat,block_cnt,blockids,bcids)
         do jb = 1,block_cnt
            owner_d = get_block_owner_d(blockids(jb))
            if (iam == owner_d) then
               if (.not. associated(btofc_blk_offset(blockids(jb))%pter)) then
                  blksiz = get_block_col_cnt_d(blockids(jb))
                  numlvl = get_block_lvl_cnt_d(blockids(jb),bcids(jb))
                  btofc_blk_offset(blockids(jb))%ncols = blksiz
                  btofc_blk_offset(blockids(jb))%nlvls = numlvl
                  allocate ( btofc_blk_offset(blockids(jb))%pter(blksiz,numlvl) )
               endif
               do k=1,btofc_blk_offset(blockids(jb))%nlvls
                  btofc_blk_offset(blockids(jb))%pter(bcids(jb),k) = glbcnt
                  curcnt = curcnt + 1
                  glbcnt = glbcnt + 1
               enddo
            endif
         enddo
      enddo
      btofc_blk_num(curp) = curcnt
      block_buf_nrecs = glbcnt
!  
      allocate ( btofc_chk_num(0:npes-1) )
      allocate ( btofc_chk_offset(begchunk:endchunk) )
      do lchnk=begchunk,endchunk
         ncol = chunks(lchunks(lchnk))%ncols
         btofc_chk_offset(lchnk)%ncols = ncol
         btofc_chk_offset(lchnk)%nlvls = pver+1
         allocate ( btofc_chk_offset(lchnk)%pter(ncol,pver+1) )
      enddo
!
      curcnt = 0
      glbcnt = 0
      do p=0,npes-1
         do curgcol=gs_col_offset(iam),gs_col_offset(iam+1)-1
            cid  = pgcols(curgcol)%chunk
            owner_p  = chunks(cid)%owner
            if (iam == owner_p) then
               i    = pgcols(curgcol)%ccol
               lchnk = chunks(cid)%lchunk
               glon   = chunks(cid)%lon(i)
               glat   = chunks(cid)%lat(i)
               block_cnt = get_block_coord_cnt_d(glon,glat)
               call get_block_coord_d(glon,glat,block_cnt,blockids,bcids)
               do jb = 1,block_cnt
                  owner_d = get_block_owner_d(blockids(jb))
                  if (p == owner_d) then
                     numlvl = get_block_lvl_cnt_d(blockids(jb),bcids(jb))
                     call get_block_levels_d(blockids(jb),bcids(jb),numlvl,levels)
                     do k=1,numlvl
                        btofc_chk_offset(lchnk)%pter(i,levels(k)+1) = glbcnt
                        curcnt = curcnt + 1
                        glbcnt = glbcnt + 1
                     enddo
                  endif
               enddo
            endif
         enddo
         btofc_chk_num(p) = curcnt
         curcnt = 0
      enddo
      chunk_buf_nrecs = glbcnt
!
! Precompute swap partners and number of steps in point-to-point
! implementations of alltoall algorithm.
! First, determine number of swaps.
!
      dp_coup_steps = 0
      do i=1,ceil2(npes)-1
         p = pair(npes,i,iam)
         if (p >= 0) then
            if ((btofc_blk_num(p) > 0 .or. btofc_chk_num(p) > 0)) then
               dp_coup_steps = dp_coup_steps + 1
            end if
         end if
      end do
!
! Second, determine swap partners.
!
      allocate( dp_coup_proc(dp_coup_steps) )
      dp_coup_steps = 0
      do i=1,ceil2(npes)-1
         p = pair(npes,i,iam)
         if (p >= 0) then
            if ((btofc_blk_num(p) > 0 .or. btofc_chk_num(p) > 0)) then
               dp_coup_steps = dp_coup_steps + 1
               dp_coup_proc(dp_coup_steps) = p
            end if
         end if
      end do
!
   endif
!
   physgrid_set = .true.   ! Set flag indicating physics grid is now set
!
   if (masterproc) then
      write(6,*) 'PHYS_GRID_INIT:  Using PCOLS=',pcols,     &
                 '  phys_loadbalance=',lbal_opt,            &
                 '  phys_alltoall=',phys_alltoall,          &
                 '  chunks_per_thread=',chunks_per_thread
   endif
!
   return
   end subroutine phys_grid_init
!
!========================================================================
!
   subroutine phys_grid_defaultopts(phys_loadbalance_out, &
                                    phys_alltoall_out, &
                                    phys_chnk_per_thd_out )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: Tom Henderson
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
     ! grid optimization option
     integer, intent(out), optional :: phys_loadbalance_out
     ! alltoall option
     integer, intent(out), optional :: phys_alltoall_out
     ! number of chunks per thread
     integer, intent(out), optional :: phys_chnk_per_thd_out
!-----------------------------------------------------------------------
     if ( present(phys_loadbalance_out) ) then
       phys_loadbalance_out = def_lbal_opt
     endif
     if ( present(phys_alltoall_out) ) then
       phys_alltoall_out = def_alltoall
     endif
     if ( present(phys_chnk_per_thd_out) ) then
       phys_chnk_per_thd_out = def_chunks_per_thread
     endif
   end subroutine phys_grid_defaultopts
!
!========================================================================
!
   subroutine phys_grid_setopts(phys_loadbalance_in, &
                                phys_alltoall_in,    &
                                phys_chnk_per_thd_in )
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: Tom Henderson
!-----------------------------------------------------------------------
   use pmgrid, only: masterproc
   use spmd_utils, only: phys_mirror_decomp_req
   use dycore, only: dycore_is
#if defined(MODCM_DP_TRANSPOSE)
   use mod_comm, only: phys_transpose_mod
#endif
!------------------------------Arguments--------------------------------
     ! grid optimization option
     integer, intent(in), optional :: phys_loadbalance_in
     ! alltoall option
     integer, intent(in), optional :: phys_alltoall_in
     ! number of chunks per thread
     integer, intent(in), optional :: phys_chnk_per_thd_in
!-----------------------------------------------------------------------
     if ( present(phys_loadbalance_in) ) then
        lbal_opt = phys_loadbalance_in
        if ((lbal_opt < min_lbal_opt).or.(lbal_opt > max_lbal_opt)) then
           if (masterproc) then
              write(6,*)                                          &
                 'PHYS_GRID_SETOPTS:  ERROR:  phys_loadbalance=', &
                 phys_loadbalance_in,                             &
                 '  is out of range.  It must be between ',       &
                 min_lbal_opt,' and ',max_lbal_opt
           endif
           call endrun
        endif
        if (lbal_opt .eq. 3) then
           phys_mirror_decomp_req = .true.
        else
           phys_mirror_decomp_req = .false.
        endif
     endif
!
     if ( present(phys_alltoall_in) ) then
        phys_alltoall = phys_alltoall_in
        if ((phys_alltoall .lt. min_alltoall) .or. &
           (phys_alltoall .gt. max_alltoall)) then
           if (masterproc) then
              write(6,*)                                          &
                 'PHYS_GRID_SET_OPTS:  ERROR:  phys_alltoall=',   &
                  phys_alltoall_in,                               &
                  '  is out of range.  It must be between ',      &
                  min_alltoall,' and ',max_alltoall
           endif
           call endrun
        endif
#if defined(SPMD)
# if defined(MODCM_DP_TRANSPOSE)
        phys_transpose_mod = phys_alltoall
# endif
#endif
     endif
!
     if ( present(phys_chnk_per_thd_in) ) then
        chunks_per_thread = phys_chnk_per_thd_in
        if (chunks_per_thread < min_chunks_per_thread) then
           if (masterproc) then
              write(6,*)                                          &
                 'PHYS_GRID_SETOPTS:  ERROR:  phys_chnk_per_thd=',&
                 phys_chnk_per_thd_in,                            &
                 ' is too small.  It must not be smaller than ',  &
                 min_chunks_per_thread
           endif
           call endrun
        endif
     endif
   end subroutine phys_grid_setopts
!
!========================================================================
!
   subroutine get_chunk_indices_p(index_beg, index_end)
!----------------------------------------------------------------------- 
! 
! Purpose: Return range of indices for local chunks
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(out) :: index_beg  ! first index used for local chunks
   integer, intent(out) :: index_end  ! last index used for local chunks
!-----------------------------------------------------------------------

   index_beg = begchunk
   index_end = endchunk

   return
   end subroutine get_chunk_indices_p
!
!========================================================================
!
   integer function get_ncols_p(lchunkid)
!----------------------------------------------------------------------- 
! 
! Purpose: Return number of columns in chunk given the local chunk id.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id

!---------------------------Local workspace-----------------------------
   integer              :: chunkid       ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   get_ncols_p = chunks(chunkid)%ncols

   return
   end function get_ncols_p
!
!========================================================================
!
   subroutine get_lat_all_p(lchunkid, latdim, lats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all global latitude indices for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: latdim        ! declared size of output array

   integer, intent(out) :: lats(latdim)  ! array of global latitude indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   do i=1,chunks(chunkid)%ncols
     lats(i) = chunks(chunkid)%lat(i)
   enddo

   return
   end subroutine get_lat_all_p
!
!========================================================================

   subroutine get_lat_vec_p(lchunkid, lth, cols, lats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global latitude indices for set of chunk columns
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid

!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: lth           ! number of column indices
   integer, intent(in)  :: cols(lth)     ! column indices

   integer, intent(out) :: lats(lth)     ! array of global latitude indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   do i=1,lth
     lats(i) = chunks(chunkid)%lat(cols(i))
   enddo

   return
   end subroutine get_lat_vec_p
!
!========================================================================

   integer function get_lat_p(lchunkid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global latitude index for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   get_lat_p = chunks(chunkid)%lat(col)

   return
   end function get_lat_p
!
!========================================================================
!
   subroutine get_lon_all_p(lchunkid, londim, lons)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all global longitude indices for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: londim        ! declared size of output array

   integer, intent(out) :: lons(londim)  ! array of global longitude indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   do i=1,chunks(chunkid)%ncols
     lons(i) = chunks(chunkid)%lon(i)
   enddo

   return
   end subroutine get_lon_all_p
!
!========================================================================

   subroutine get_lon_vec_p(lchunkid, lth, cols, lons)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global longitude indices for set of chunk columns
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: lth           ! number of column indices
   integer, intent(in)  :: cols(lth)     ! column indices

   integer, intent(out) :: lons(lth)     ! array of global longitude indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   do i=1,lth
     lons(i) = chunks(chunkid)%lon(cols(i))
   enddo

   return
   end subroutine get_lon_vec_p
!
!========================================================================

   integer function get_lon_p(lchunkid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return global longitude index for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   get_lon_p = chunks(chunkid)%lon(col)

   return
   end function get_lon_p
!
!========================================================================
!
   subroutine get_rlat_all_p(lchunkid, rlatdim, rlats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all latitudes (in radians) for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: rlatdim        ! declared size of output array

   real(r8), intent(out) :: rlats(rlatdim)! array of latitudes

!---------------------------Local workspace-----------------------------
   integer :: i                           ! loop index
   integer :: chunkid                     ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   do i=1,chunks(chunkid)%ncols
     rlats(i) = clat_p(chunks(chunkid)%lat(i))
   enddo

   return
   end subroutine get_rlat_all_p
!
!========================================================================

   subroutine get_rlat_vec_p(lchunkid, lth, cols, rlats)
!----------------------------------------------------------------------- 
! 
! Purpose: Return latitudes (in radians) for set of chunk columns
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: lth           ! number of column indices
   integer, intent(in)  :: cols(lth)     ! column indices

   real(r8), intent(out) :: rlats(lth)   ! array of latitudes

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   do i=1,lth
     rlats(i) = clat_p(chunks(chunkid)%lat(cols(i)))
   enddo

   return
   end subroutine get_rlat_vec_p
!
!========================================================================

   real(r8) function get_rlat_p(lchunkid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return latitude (in radians) for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   get_rlat_p = clat_p(chunks(chunkid)%lat(col))

   return
   end function get_rlat_p
!
!
!========================================================================
!
   subroutine get_rlon_all_p(lchunkid, rlondim, rlons)
!----------------------------------------------------------------------- 
! 
! Purpose: Return all longitudes (in radians) for chunk
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: rlondim        ! declared size of output array

   real(r8), intent(out) :: rlons(rlondim)! array of longitudes

!---------------------------Local workspace-----------------------------
   integer :: i                           ! loop index
   integer :: chunkid                     ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   do i=1,chunks(chunkid)%ncols
     rlons(i) = clon_p(chunks(chunkid)%lon(i),chunks(chunkid)%lat(i))
   enddo

   return
   end subroutine get_rlon_all_p
!
!========================================================================

   subroutine get_rlon_vec_p(lchunkid, lth, cols, rlons)
!----------------------------------------------------------------------- 
! 
! Purpose: Return longitudes (in radians) for set of chunk columns
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: lth           ! number of column indices
   integer, intent(in)  :: cols(lth)     ! column indices

   real(r8), intent(out) :: rlons(lth)   ! array of longitudes

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   do i=1,lth
     rlons(i) = clon_p(chunks(chunkid)%lon(cols(i)), &
                       chunks(chunkid)%lat(cols(i)))
   enddo

   return
   end subroutine get_rlon_vec_p
!
!========================================================================

   real(r8) function get_rlon_p(lchunkid, col)
!----------------------------------------------------------------------- 
! 
! Purpose: Return longitude (in radians) for chunk column
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use ppgrid
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lchunkid      ! local chunk id
   integer, intent(in)  :: col           ! column index

!---------------------------Local workspace-----------------------------
   integer :: chunkid                    ! global chunk id

!-----------------------------------------------------------------------
   chunkid = lchunks(lchunkid)
   get_rlon_p = clon_p(chunks(chunkid)%lon(col),chunks(chunkid)%lat(col))

   return
   end function get_rlon_p
!
!========================================================================

logical function chunk_index (idx)
!----------------------------------------------------------------------- 
! 
! Purpose: Identify whether index is for a latitude or a chunk
! 
! Method: Quick hack, using convention that local chunk indices do not
!         overlap latitude index range
! 
! Author: Pat Worley
! 
!-----------------------------------------------------------------------
   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: idx              ! latitude or chunk index
!
!-----------------------------------------------------------------------
!
   if ((idx >= begchunk) .and. (idx <= endchunk)) then
      chunk_index = .true.
   else
      chunk_index = .false.
   endif
!
   return
   end function chunk_index

!
!========================================================================

   integer function get_chunk_owner_p(loni,latj)
!----------------------------------------------------------------------- 
! 
! Purpose: Return owner of chunk at location loni, latj
! 
! Method: 
! 
! Author: R. Jacob
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: loni     ! longitude index
   integer, intent(in)  :: latj     ! latitude index

!-----------------------------------------------------------------------
   
   get_chunk_owner_p = chunks(knuhcs(loni,latj)%chunkid)%owner

   return
   end function get_chunk_owner_p
!
!========================================================================

   subroutine get_chunk_coord_p(lth, xylons, xylats, ckcols, ckcids)
!----------------------------------------------------------------------- 
! 
! Purpose: Return local chunk coordinates for corresponding global 
!          (lon,lat) coordinates
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lth           ! number of coordinates
   integer, intent(in)  :: xylons(lth)   ! longitude indices
   integer, intent(in)  :: xylats(lth)   ! latitude indices

   integer, intent(out) :: ckcols(lth)   ! column indices
   integer, intent(out) :: ckcids(lth)   ! local chunk indices

!---------------------------Local workspace-----------------------------
   integer :: i                          ! loop index

!-----------------------------------------------------------------------
   do i=1,lth
      if (chunks(knuhcs(xylons(i),xylats(i))%chunkid)%owner .eq. iam) then
         ckcols(i) = knuhcs(xylons(i),xylats(i))%col
         ckcids(i) = chunks(knuhcs(xylons(i),xylats(i))%chunkid)%lchunk
      else
         ckcols(i) = -1
         ckcids(i) = -1
      endif
   enddo

   return
   end subroutine get_chunk_coord_p
!
!========================================================================


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: get_chunk_coord_owner_p
!
! !INTERFACE:
   subroutine get_chunk_coord_owner_p(lth, lons, lats, lchunks, cols, owners)
!
! !PARAMETERS:
   implicit none
   integer, intent(in)  :: lth               ! number of column indices
   integer, intent(in)  :: lons(lth)         ! longitude vector
   integer, intent(in)  :: lats(lth)         ! latitude vector
!
! !RETURN VALUE:
   integer, intent(out) :: lchunks(lth)      ! local chunk index vector
   integer, intent(out) :: cols(lth)         ! column vector
   integer, intent(out) :: owners(lth)       ! column owner vector
!
! !LOCAL VARIABLES:
   integer :: i                              ! loop index
!
! !CALLED FROM:
! subroutine lp_coupling_init() in module lp_coupling (lp_coupling.F90)
!
! !DESCRIPTION:
! Fill vectors of lchunks, cols, and owners for each longitude/latitude
! index pair.
!
! !REVISION HISTORY:
! 2002.09.11  Forrest Hoffman  Creation.
!
!EOP
!-----------------------------------------------------------------------
! $Id: phys_grid.F90,v 1.7.4.21.4.1 2004/05/20 18:36:08 mvr Exp $
! $Author: mvr $
!-----------------------------------------------------------------------
   do i = 1,lth
      lchunks(i) = chunks(knuhcs(lons(i),lats(i))%chunkid)%lchunk
      cols(i) = knuhcs(lons(i),lats(i))%col
      owners(i) = chunks(knuhcs(lons(i),lats(i))%chunkid)%owner
   end do

   end subroutine get_chunk_coord_owner_p
!
!========================================================================

   subroutine buff_to_chunk(mdim,nlond,lbuff, localchunks)
!-----------------------------------------------------------------------
!
! Purpose: copy local long/lat buffer to local chunk data structure.
!          Needed for cpl6.
!
! Method:
!
! Author: Pat Worley and Robert Jacob
!
!-----------------------------------------------------------------------
   use pmgrid, only: iam
   use rgrid, only: nlon
!------------------------------Arguments--------------------------------
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: nlond      ! declared length of middle dimension
   real(r8), intent(in) :: lbuff(nlcols,mdim) ! local buff

   real(r8), intent(out):: localchunks(pcols,mdim,begchunk:endchunk) ! local chunks


!---------------------------Local workspace-----------------------------
   integer :: i,j,m,n                      ! loop indices
!-----------------------------------------------------------------------

   n = 1
   do j = 1, plat
     do i=1,nlon(j)
       if(chunks(knuhcs(i,j)%chunkid)%owner .eq. iam) then
         do m=1,mdim
           localchunks(knuhcs(i,j)%col,m,chunks(knuhcs(i,j)%chunkid)%lchunk) = lbuff(n,m)
         end do
         n = n + 1
       endif
     enddo
   end do

   return
   end subroutine buff_to_chunk



   subroutine scatter_field_to_chunk(fdim,mdim,ldim, &
                                     nlond,globalfield,localchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute longitude/latitude field
!          to decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, masterproc
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: nlond     ! declared number of longitudes
   real(r8), intent(in) :: globalfield(fdim,nlond,mdim,plat,ldim) 
                                    ! global field

   real(r8), intent(out):: localchunks(fdim,pcols,mdim, &
                                       begchunk:endchunk,ldim) 
                                    ! local chunks

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local longitude index

#if ( defined SPMD )
   real(r8) gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be scattered
   real(r8) lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of scattered
                                         !  vector
   integer :: displs(0:npes-1)           ! scatter displacements
   integer :: sndcnts(0:npes-1)          ! scatter send counts
   integer :: recvcnt                    ! scatter receive count
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
#if ( defined SPMD )
   displs(0) = 0
   sndcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + sndcnts(p-1)
     sndcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   recvcnt = fdim*mdim*ldim*nlcols

   if (masterproc) then

! copy field into global (process-ordered) chunked data structure

      do l=1,ldim
         do m=1,mdim
            do f=1,fdim
!DIR$ PREFERVECTOR
               do i=1,ngcols
                  cid = pgcols(i)%chunk
                  lid = pgcols(i)%ccol
                  gfield_p(f,m,l,i) = &
                     globalfield(f,chunks(cid)%lon(lid), m, &
                                 chunks(cid)%lat(lid),l)
               end do
            end do
         end do
      end do
   endif

! scatter to other processes
! (pgcols ordering consistent with begchunk:endchunk 
! local ordering)

#if ( defined TIMING_BARRIERS )
   call t_startf ('sync_scat_ftoc')
   call mpibarrier (mpicom)
   call t_stopf ('sync_scat_ftoc')
#endif

   call mpiscatterv(gfield_p, sndcnts, displs, mpir8, &
                    lfield_p, recvcnt, mpir8, 0, mpicom)

! copy into local chunked data structure

   do l=1,ldim
      do m=1,mdim
         do f=1,fdim
!DIR$ PREFERVECTOR
            do i=1,nlcols
               cid = pgcols(beglcol+i)%chunk
               lcid = chunks(cid)%lchunk
               lid = pgcols(beglcol+i)%ccol
               localchunks(f,lid,m,lcid,l) = &
                 lfield_p(f, m, l, i)
            end do
         end do
      end do
   end do
#else

! copy field into chunked data structure
! (pgcol ordering chosen to reflect begchunk:endchunk 
!  local ordering)

   do l=1,ldim
      do m=1,mdim
         do f=1,fdim
!DIR$ PREFERVECTOR
            do i=1,ngcols
               cid = pgcols(i)%chunk
               lcid = chunks(cid)%lchunk
               lid = pgcols(i)%ccol
               localchunks(f,lid,m,lcid,l) = &
                  globalfield(f,chunks(cid)%lon(lid), m, &
                              chunks(cid)%lat(lid),l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine scatter_field_to_chunk
!========================================================================

   subroutine scatter_field_to_chunk4(fdim,mdim,ldim, &
                                      nlond,globalfield,localchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute longitude/latitude field
!          to decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, masterproc

   implicit none
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: nlond     ! declared number of longitudes
   real(r4), intent(in) :: globalfield(fdim,nlond,mdim,plat,ldim) 
                                    ! global field

   real(r4), intent(out):: localchunks(fdim,pcols,mdim, &
                                       begchunk:endchunk,ldim) 
                                    ! local chunks

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local longitude index

#if ( defined SPMD )
   real(r4) gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be scattered
   real(r4) lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of scattered
                                         !  vector
   integer :: displs(0:npes-1)           ! scatter displacements
   integer :: sndcnts(0:npes-1)          ! scatter send counts
   integer :: recvcnt                    ! scatter receive count
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
#if ( defined SPMD )
   displs(0) = 0
   sndcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + sndcnts(p-1)
     sndcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   recvcnt = fdim*mdim*ldim*nlcols

   if (masterproc) then
      ! copy field into global (process-ordered) chunked data structure
      do i=1,ngcols
         cid = pgcols(i)%chunk
         lid = pgcols(i)%ccol
         do l=1,ldim
            do m=1,mdim
               do f=1,fdim
                  gfield_p(f,m,l,i) = &
                     globalfield(f,chunks(cid)%lon(lid), m, &
                                 chunks(cid)%lat(lid),l)
               end do
            end do
         end do
      end do
   endif

! scatter to other processes
! (pgcols ordering consistent with begchunk:endchunk 
!  local ordering)

#if ( defined TIMING_BARRIERS )
   call t_startf ('sync_scat_ftoc')
   call mpibarrier (mpicom)
   call t_stopf ('sync_scat_ftoc')
#endif

   call mpiscatterv(gfield_p, sndcnts, displs, mpir4, &
                    lfield_p, recvcnt, mpir4, 0, mpicom)

! copy into local chunked data structure

   do i=1,nlcols
      cid = pgcols(beglcol+i)%chunk
      lcid = chunks(cid)%lchunk
      lid = pgcols(beglcol+i)%ccol
      do l=1,ldim
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                 lfield_p(f, m, l, i)
            end do
         end do
      end do
   end do
#else

   ! copy field into chunked data structure
   ! (pgcol ordering chosen to reflect begchunk:endchunk 
   !  local ordering)
   do l=1,ldim
      do i=1,ngcols
         cid = pgcols(i)%chunk
         lcid = chunks(cid)%lchunk
         lid = pgcols(i)%ccol
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                  globalfield(f,chunks(cid)%lon(lid), m, &
                              chunks(cid)%lat(lid),l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine scatter_field_to_chunk4
!========================================================================

   subroutine scatter_field_to_chunk_int(fdim,mdim,ldim, &
                                         nlond,globalfield,localchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Distribute longitude/latitude field
!          to decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, masterproc
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: nlond     ! declared number of longitudes
   integer, intent(in) :: globalfield(fdim,nlond,mdim,plat,ldim) 
                                    ! global field

   integer, intent(out):: localchunks(fdim,pcols,mdim, &
                                       begchunk:endchunk,ldim) 
                                    ! local chunks

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local longitude index

#if ( defined SPMD )
   integer gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be scattered
   integer lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of scattered
                                         !  vector
   integer :: displs(0:npes-1)           ! scatter displacements
   integer :: sndcnts(0:npes-1)          ! scatter send counts
   integer :: recvcnt                    ! scatter receive count
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
#if ( defined SPMD )
   displs(0) = 0
   sndcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + sndcnts(p-1)
     sndcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   recvcnt = fdim*mdim*ldim*nlcols

   if (masterproc) then

! copy field into global (process-ordered) chunked data structure

      do i=1,ngcols
         cid = pgcols(i)%chunk
         lid = pgcols(i)%ccol
         do l=1,ldim
            do m=1,mdim
               do f=1,fdim
                  gfield_p(f,m,l,i) = &
                     globalfield(f,chunks(cid)%lon(lid), m, &
                                 chunks(cid)%lat(lid),l)
               end do
            end do
         end do
      end do
   endif

! scatter to other processes
! (pgcols ordering consistent with begchunk:endchunk 
!  local ordering)

#if ( defined TIMING_BARRIERS )
   call t_startf ('sync_scat_ftoc')
   call mpibarrier (mpicom)
   call t_stopf ('sync_scat_ftoc')
#endif

   call mpiscatterv(gfield_p, sndcnts, displs, mpiint, &
                    lfield_p, recvcnt, mpiint, 0, mpicom)

! copy into local chunked data structure

   do i=1,nlcols
      cid = pgcols(beglcol+i)%chunk
      lcid = chunks(cid)%lchunk
      lid = pgcols(beglcol+i)%ccol
      do l=1,ldim
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                 lfield_p(f, m, l, i)
            end do
         end do
      end do
   end do
#else

! copy field into chunked data structure
! (pgcol ordering chosen to reflect begchunk:endchunk 
!  local ordering)

   do l=1,ldim
      do i=1,ngcols
         cid = pgcols(i)%chunk
         lcid = chunks(cid)%lchunk
         lid = pgcols(i)%ccol
         do m=1,mdim
            do f=1,fdim
               localchunks(f,lid,m,lcid,l) = &
                  globalfield(f,chunks(cid)%lon(lid), m, &
                              chunks(cid)%lat(lid),l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine scatter_field_to_chunk_int
!
!========================================================================
!
   subroutine chunk_to_buff(mdim,nlond,localchunks,lbuff)

!-----------------------------------------------------------------------
!
! Purpose: Copy from local chunk data structure
!          to local longitude/latitude buffer.  Needed for cpl6
!          (local = on processor)
!
! Method:
!
! Author: Pat Worley and Robert Jacob
!-----------------------------------------------------------------------
   use pmgrid, only: iam
   use rgrid, only: nlon
!------------------------------Arguments--------------------------------
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: nlond     ! declared number of longitudes
   real(r8), intent(in):: localchunks(pcols,mdim, begchunk:endchunk) ! local chunks

   real(r8), intent(out) :: lbuff(nlcols,mdim) ! local buff

!---------------------------Local workspace-----------------------------
   integer :: i,j,m,n                  ! loop indices

!-----------------------------------------------------------------------
   n = 1
   do j = 1, plat
     do i=1,nlon(j)
       if(chunks(knuhcs(i,j)%chunkid)%owner .eq. iam) then
         do m=1,mdim
           lbuff(n,m)=localchunks(knuhcs(i,j)%col,m,chunks(knuhcs(i,j)%chunkid)%lchunk)
         end do
         n = n + 1
       endif
     enddo
   end do

   return
   end subroutine chunk_to_buff

!
!========================================================================
!
   subroutine gather_chunk_to_field(fdim,mdim,ldim, &
                                     nlond,localchunks,globalfield)

!----------------------------------------------------------------------- 
! 
! Purpose: Reconstruct longitude/latitude field
!          from decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, masterproc
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: nlond     ! declared number of longitudes
   real(r8), intent(in):: localchunks(fdim,pcols,mdim, &
                                      begchunk:endchunk,ldim) 
                                    ! local chunks

   real(r8), intent(out) :: globalfield(fdim,nlond,mdim,plat,ldim) 
                                    ! global field

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local longitude index

#if ( defined SPMD )
   real(r8) gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be gathered
   real(r8) lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of gather
                                         !  vector
   integer :: displs(0:npes-1)           ! gather displacements
   integer :: rcvcnts(0:npes-1)          ! gather receive count
   integer :: sendcnt                    ! gather send counts
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
#if ( defined SPMD )
   displs(0) = 0
   rcvcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + rcvcnts(p-1)
     rcvcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   sendcnt = fdim*mdim*ldim*nlcols

! copy into local gather data structure

   do l=1,ldim
      do m=1,mdim
         do f=1,fdim
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
            do i=1,nlcols
               cid = pgcols(beglcol+i)%chunk
               lcid = chunks(cid)%lchunk
               lid = pgcols(beglcol+i)%ccol
               lfield_p(f, m, l, i) = &
                  localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

! gather from other processes

#if ( defined TIMING_BARRIERS )
   call t_startf ('sync_gath_ctof')
   call mpibarrier (mpicom)
   call t_stopf ('sync_gath_ctof')
#endif

   call mpigatherv(lfield_p, sendcnt, mpir8, &
                   gfield_p, rcvcnts, displs, mpir8, 0, mpicom)

   if (masterproc) then

! copy gathered columns into lon/lat field

      do l=1,ldim
         do m=1,mdim
            do f=1,fdim
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
               do i=1,ngcols
                  cid = pgcols(i)%chunk
                  lid = pgcols(i)%ccol
                  globalfield(f,chunks(cid)%lon(lid), m, &
                              chunks(cid)%lat(lid),l)    &
                  = gfield_p(f,m,l,i)
               end do
            end do
         end do
      end do
   endif

#else

   ! copy chunked data structure into lon/lat field
   ! (pgcol ordering chosen to reflect begchunk:endchunk 
   !  local ordering)
   do l=1,ldim
      do m=1,mdim
         do f=1,fdim
!DIR$ PREFERVECTOR, PREFERSTREAM
!DIR$ CONCURRENT
            do i=1,ngcols
               cid = pgcols(i)%chunk
               lcid = chunks(cid)%lchunk
               lid = pgcols(i)%ccol
               globalfield(f,chunks(cid)%lon(lid), m, &
                           chunks(cid)%lat(lid),l)    &
               = localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine gather_chunk_to_field

!
!========================================================================
!
   subroutine gather_chunk_to_field4 (fdim,mdim,ldim, &
                                      nlond,localchunks,globalfield)

!----------------------------------------------------------------------- 
! 
! Purpose: Reconstruct longitude/latitude field
!          from decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, masterproc
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: nlond     ! declared number of longitudes
   real(r4), intent(in):: localchunks(fdim,pcols,mdim, &
                                      begchunk:endchunk,ldim) 
                                    ! local chunks

   real(r4), intent(out) :: globalfield(fdim,nlond,mdim,plat,ldim) 
                                    ! global field

!---------------------------Local workspace-----------------------------
   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local longitude index

#if ( defined SPMD )
   real(r4) gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be gathered
   real(r4) lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of gather
                                         !  vector
   integer :: displs(0:npes-1)           ! gather displacements
   integer :: rcvcnts(0:npes-1)          ! gather receive count
   integer :: sendcnt                    ! gather send counts
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
#if ( defined SPMD )
   displs(0) = 0
   rcvcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + rcvcnts(p-1)
     rcvcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   sendcnt = fdim*mdim*ldim*nlcols

! copy into local gather data structure

   do i=1,nlcols
      cid = pgcols(beglcol+i)%chunk
      lcid = chunks(cid)%lchunk
      lid = pgcols(beglcol+i)%ccol
      do l=1,ldim
         do m=1,mdim
            do f=1,fdim
               lfield_p(f, m, l, i) = &
                  localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

! gather from other processes

#if ( defined TIMING_BARRIERS )
   call t_startf ('sync_gath_ctof')
   call mpibarrier (mpicom)
   call t_stopf ('sync_gath_ctof')
#endif

   call mpigatherv(lfield_p, sendcnt, mpir4, &
                   gfield_p, rcvcnts, displs, mpir4, 0, mpicom)

   if (masterproc) then

! copy gathered columns into lon/lat field

      do i=1,ngcols
         cid = pgcols(i)%chunk
         lid = pgcols(i)%ccol
         do l=1,ldim
            do m=1,mdim
               do f=1,fdim
                  globalfield(f,chunks(cid)%lon(lid), m, &
                              chunks(cid)%lat(lid),l)    &
                  = gfield_p(f,m,l,i)
               end do
            end do
         end do
      end do
   endif

#else

! copy chunked data structure into lon/lat field
! (pgcol ordering chosen to reflect begchunk:endchunk 
!  local ordering)

   do l=1,ldim
      do i=1,ngcols
         cid = pgcols(i)%chunk
         lcid = chunks(cid)%lchunk
         lid = pgcols(i)%ccol
         do m=1,mdim
            do f=1,fdim
               globalfield(f,chunks(cid)%lon(lid), m, &
                           chunks(cid)%lat(lid),l)    &
               = localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine gather_chunk_to_field4

!
!========================================================================
!
   subroutine gather_chunk_to_field_int (fdim,mdim,ldim, &
                                         nlond,localchunks,globalfield)

!----------------------------------------------------------------------- 
! 
! Purpose: Reconstruct longitude/latitude field
!          from decomposed chunk data structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, masterproc
!------------------------------Arguments--------------------------------
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   integer, intent(in) :: nlond     ! declared number of longitudes
   integer, intent(in):: localchunks(fdim,pcols,mdim,begchunk:endchunk,ldim) ! local chunks
!JR Changed globalfield to inout because slaves under lf95 pass a bogus argument, which will result
!JR in trash being written to useful memory if intent(out) is specified.  THIS SHOULD BE FIXED!!!
   integer, intent(inout) :: globalfield(fdim,nlond,mdim,plat,ldim) ! global field

!---------------------------Local workspace-----------------------------

   integer :: f,i,m,l,p                  ! loop indices
   integer :: cid                        ! global chunk id
   integer :: lcid                       ! local chunk id
   integer :: lid                        ! local longitude index

#if ( defined SPMD )
   integer gfield_p(fdim,mdim,ldim,ngcols) 
                                         ! vector to be gathered
   integer lfield_p(fdim,mdim,ldim,nlcols) 
                                         ! local component of gather
                                         !  vector
   integer :: displs(0:npes-1)           ! gather displacements
   integer :: rcvcnts(0:npes-1)          ! gather receive count
   integer :: sendcnt                    ! gather send counts
   integer :: beglcol                    ! beginning index for local columns
                                         !  in global column ordering
#endif

!-----------------------------------------------------------------------
#if ( defined SPMD )
   displs(0) = 0
   rcvcnts(0) = fdim*mdim*ldim*gs_col_num(0)
   beglcol = 0
   do p=1,npes-1
     displs(p) = displs(p-1) + rcvcnts(p-1)
     rcvcnts(p) = fdim*mdim*ldim*gs_col_num(p)
     if (p <= iam) then
        beglcol = beglcol + gs_col_num(p-1)
     endif
   enddo
   sendcnt = fdim*mdim*ldim*nlcols

! copy into local gather data structure

   do i=1,nlcols
      cid = pgcols(beglcol+i)%chunk
      lcid = chunks(cid)%lchunk
      lid = pgcols(beglcol+i)%ccol
      do l=1,ldim
         do m=1,mdim
            do f=1,fdim
               lfield_p(f, m, l, i) = &
                  localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

! gather from other processes

#if ( defined TIMING_BARRIERS )
   call t_startf ('sync_gath_ctof')
   call mpibarrier (mpicom)
   call t_stopf ('sync_gath_ctof')
#endif

   call mpigatherv(lfield_p, sendcnt, mpiint, &
                   gfield_p, rcvcnts, displs, mpiint, 0, mpicom)

   if (masterproc) then

! copy gathered columns into lon/lat field

      do i=1,ngcols
         cid = pgcols(i)%chunk
         lid = pgcols(i)%ccol
         do l=1,ldim
            do m=1,mdim
               do f=1,fdim
                  globalfield(f,chunks(cid)%lon(lid), m, &
                              chunks(cid)%lat(lid),l)    &
                  = gfield_p(f,m,l,i)
               end do
            end do
         end do
      end do
   endif

#else

   ! copy chunked data structure into lon/lat field
   ! (pgcol ordering chosen to reflect begchunk:endchunk 
   !  local ordering)
   do l=1,ldim
      do i=1,ngcols
         cid = pgcols(i)%chunk
         lcid = chunks(cid)%lchunk
         lid = pgcols(i)%ccol
         do m=1,mdim
            do f=1,fdim
               globalfield(f,chunks(cid)%lon(lid), m, &
                           chunks(cid)%lat(lid),l)    &
               = localchunks(f,lid,m,lcid,l)
            end do
         end do
      end do
   end do

#endif

   return
   end subroutine gather_chunk_to_field_int

!
!========================================================================
!
   subroutine write_field_from_chunk(iu,fdim,mdim,ldim,localchunks)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Write longitude/latitude field from decomposed chunk data 
!          structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: masterproc
!------------------------------Arguments--------------------------------
   integer, intent(in) :: iu        ! logical unit
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension
   real(r8), intent(in):: localchunks(fdim,pcols,mdim,begchunk:endchunk,ldim) ! local chunks

!---------------------------Local workspace-----------------------------

   integer :: ioerr                 ! error return

   real(r8), allocatable :: globalfield(:,:,:,:,:)
                                    ! global field
!-----------------------------------------------------------------------

   allocate(globalfield(fdim,plon,mdim,plat,ldim))

   call gather_chunk_to_field (fdim,mdim,ldim,plon,localchunks,globalfield)
                               
   if (masterproc) then
      write (iu,iostat=ioerr) globalfield
      if (ioerr /= 0 ) then
         write (6,*) 'WRITE_FIELD_FROM_CHUNK ioerror ', ioerr,' on i/o unit = ',iu
         call endrun
      end if
   endif

   deallocate(globalfield)

   return
   end subroutine write_field_from_chunk

!
!========================================================================
!
   subroutine read_chunk_from_field(iu,fdim,mdim,ldim,localchunks)

!----------------------------------------------------------------------- 
! 
!                          
! Purpose: Write longitude/latitude field from decomposed chunk data 
!          structure
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: masterproc
!------------------------------Arguments--------------------------------
   integer, intent(in) :: iu        ! logical unit
   integer, intent(in) :: fdim      ! declared length of first dimension
   integer, intent(in) :: mdim      ! declared length of middle dimension
   integer, intent(in) :: ldim      ! declared length of last dimension

   real(r8), intent(out):: localchunks(fdim,pcols,mdim,begchunk:endchunk,ldim) ! local chunks

!---------------------------Local workspace-----------------------------

   integer :: ioerr                 ! error return

   real(r8), allocatable :: globalfield(:,:,:,:,:)
                                    ! global field
!-----------------------------------------------------------------------

   allocate(globalfield(fdim,plon,mdim,plat,ldim))

   if (masterproc) then
      read (iu,iostat=ioerr) globalfield
      if (ioerr /= 0 ) then
         write (6,*) 'READ_CHUNK_FROM_FIELD ioerror ', ioerr,' on i/o unit = ',iu
         call endrun
      end if
   endif

   call scatter_field_to_chunk (fdim,mdim,ldim,plon,globalfield,localchunks)

   deallocate(globalfield)

   return
   end subroutine read_chunk_from_field
!
!========================================================================

   subroutine transpose_block_to_chunk(record_size, &
                                       block_buffer, chunk_buffer)
                                       
!----------------------------------------------------------------------- 
! 
! Purpose: Transpose buffer containing decomposed 
!          longitude/latitude fields to buffer
!          containing decomposed chunk data structures
! 
! Method: 
! 
! Author: Patrick Worley
! Modified: Art Mirin, Jan 04, to add support for mod_comm
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, masterproc
#if ( defined SPMD )
   use swap_comm
# if defined(MODCM_DP_TRANSPOSE)
   use mod_comm, only: numpro, blockdescriptor, mp_sendirr, mp_recvirr,  &
                       get_partneroffset
# endif
#endif
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtag  = 1000
!------------------------------Arguments--------------------------------
   integer, intent(in) :: record_size  ! per column amount of data 
   real(r8), intent(in) :: block_buffer(record_size*block_buf_nrecs)
                                       ! buffer of block data to be
                                       ! transposed

   real(r8), intent(out):: chunk_buffer(record_size*chunk_buf_nrecs)
                                       ! buffer of chunk data 
                                       ! transposed into

!---------------------------Local workspace-----------------------------
#if ( defined SPMD )
   integer :: sdispls(0:npes-1)        ! send displacements
   integer :: sndcnts(0:npes-1)        ! send counts
   integer :: rdispls(0:npes-1)        ! receive displacements
   integer :: rcvcnts(0:npes-1)        ! receive counts
   integer :: offset_s                 ! send displacement + 1
   integer :: offset_r                 ! receive displacement + 1
   integer :: i, p                     ! loop indices
   integer :: procid                   ! swap processor id
   integer :: step                     ! step in alltoall alg.
   integer :: sndids(dp_coup_steps)    ! nonblocking MPI send request ids
   integer :: rcvids(dp_coup_steps)    ! nonblocking MPI recv request ids
# if defined(MODCM_DP_TRANSPOSE)
   type (blockdescriptor), allocatable, save :: sendbl(:), recvbl(:)
# endif
   integer ione, ierror, mod_method
   integer first_time_through
   data first_time_through / 0 /
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
   sdispls(0) = 0
   sndcnts(0) = record_size*btofc_blk_num(0)
   do p=1,npes-1
     sdispls(p) = sdispls(p-1) + sndcnts(p-1)
     sndcnts(p) = record_size*btofc_blk_num(p)
   enddo
!
   rdispls(0) = 0
   rcvcnts(0) = record_size*btofc_chk_num(0)
   do p=1,npes-1
     rdispls(p) = rdispls(p-1) + rcvcnts(p-1)
     rcvcnts(p) = record_size*btofc_chk_num(p)
   enddo

#if ( defined TIMING_BARRIERS )
   call t_startf ('sync_tran_btoc')
   call mpibarrier (mpicom)
   call t_stopf ('sync_tran_btoc')
#endif

   if (phys_alltoall .eq. 0) then
!
      call mpialltoallv(block_buffer, sndcnts, sdispls, mpir8, &
                        chunk_buffer, rcvcnts, rdispls, mpir8, &
                        mpicom)
!
   elseif (phys_alltoall .eq. 1) then
!
! Post receive requests.
      call swap1m(dp_coup_steps, msgtag, dp_coup_proc, rcvcnts, & 
                  rdispls, record_size*chunk_buf_nrecs, chunk_buffer, rcvids)
!     
! Copy local data to new location.
      if (sndcnts(iam) > 0) then
         do i=1,sndcnts(iam)
            chunk_buffer(rdispls(iam)+i) = block_buffer(sdispls(iam)+i)
         enddo
      end if
!
! Post send requests and wait for receive requests to complete.
      do step=1,dp_coup_steps
         procid = dp_coup_proc(step)
         offset_s = sdispls(procid)+1
         offset_r = rdispls(procid)+1
         call swap2(msgtag, procid, &
                    sndcnts(procid), block_buffer(offset_s), &
                    sndids(step), rcvcnts(procid), &
                    chunk_buffer(offset_r), rcvids(step))
      enddo
!
! Wait for send requests to complete.
      call swap3m(dp_coup_steps, msgtag, dp_coup_proc, sndids, rcvcnts, & 
                  rdispls, record_size*chunk_buf_nrecs, chunk_buffer, rcvids)
!
# if defined(MODCM_DP_TRANSPOSE)
   elseif (phys_alltoall .ge. modmin_alltoall) then
!
! This branch uses mod_comm. Admissable values of phys_alltoall are 11,12,13 and 14. Each value corresponds
!   to a differerent option within mod_comm of implementing the communication. That option is expressed
!   internally to mod_comm using the variable mod_method defined below; mod_method will have values 0,1,2
!   or 3 and is defined as phys_alltoall - modmin_alltoall, where modmin_alltoall equals 11.
!
      if (first_time_through .eq. 0) then
         first_time_through = 1
         mod_method = phys_alltoall - modmin_alltoall
         ione = 1
         allocate( sendbl(0:numpro-1) )
         allocate( recvbl(0:numpro-1) )

         do p = 0,numpro-1

            sendbl(p)%method = mod_method
            recvbl(p)%method = mod_method

            allocate( sendbl(p)%blocksizes(1) )
            allocate( sendbl(p)%displacements(1) )
            allocate( recvbl(p)%blocksizes(1) )
            allocate( recvbl(p)%displacements(1) )

            if ( sndcnts(p) .ne. 0 ) then

               call MPI_TYPE_INDEXED(ione, sndcnts(p),   &
                    sdispls(p), mpir8, &
                    sendbl(p)%type, ierror)
               call MPI_TYPE_COMMIT(sendbl(p)%type, ierror)

               sendbl(p)%blocksizes(1) = sndcnts(p)
               sendbl(p)%displacements(1) = sdispls(p)
               sendbl(p)%partneroffset = 0 ! not properly initialized - do not use Mpi2

            else

               sendbl(p)%type = MPI_DATATYPE_NULL

               sendbl(p)%blocksizes(1) = 0
               sendbl(p)%displacements(1) = 0
               sendbl(p)%partneroffset = 0

            endif

            if ( rcvcnts(p) .ne. 0) then

               call MPI_TYPE_INDEXED(ione, rcvcnts(p),   &
                    rdispls(p), mpir8, &
                    recvbl(p)%type, ierror)
               call MPI_TYPE_COMMIT(recvbl(p)%type, ierror)

               recvbl(p)%blocksizes(1) = rcvcnts(p)
               recvbl(p)%displacements(1) = rdispls(p)
               recvbl(p)%partneroffset = 0 ! not properly initialized - do not use Mpi2
            else

               recvbl(p)%type = MPI_DATATYPE_NULL

               recvbl(p)%blocksizes(1) = 0
               recvbl(p)%displacements(1) = 0
               recvbl(p)%partneroffset = 0

            endif

         enddo

         call get_partneroffset(sendbl, recvbl)

      endif

      call mp_sendirr(block_buffer, sendbl, recvbl, chunk_buffer)
      call mp_recvirr(chunk_buffer, recvbl)
# endif

   endif
!
#endif
   return
   end subroutine transpose_block_to_chunk
!
!========================================================================

   subroutine block_to_chunk_send_pters(blockid, fdim, ldim, &
                                        record_size, pter)
!----------------------------------------------------------------------- 
! 
! Purpose: Return pointers into send buffer where column from decomposed 
!          longitude/latitude fields should be copied to
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid      ! block index
   integer, intent(in) :: fdim         ! first dimension of pter array
   integer, intent(in) :: ldim         ! last dimension of pter array
   integer, intent(in) :: record_size  ! per coordinate amount of data 

   integer, intent(out) :: pter(fdim,ldim)  ! buffer offsets
!---------------------------Local workspace-----------------------------
   integer :: i, k                     ! loop indices
!-----------------------------------------------------------------------
   if ((btofc_blk_offset(blockid)%ncols > fdim) .or. &
       (btofc_blk_offset(blockid)%nlvls > ldim)) then
      write(6,*) "BLOCK_TO_CHUNK_SEND_PTERS: pter array dimensions ", &
                 "not large enough: (",fdim,",",ldim,") not >= (", &
                  btofc_blk_offset(blockid)%ncols,",", &
                  btofc_blk_offset(blockid)%nlvls,")"
      call endrun()
   endif
!
   do k=1,btofc_blk_offset(blockid)%nlvls
      do i=1,btofc_blk_offset(blockid)%ncols
         pter(i,k) = 1 + record_size* &
                     (btofc_blk_offset(blockid)%pter(i,k))
      enddo
      do i=btofc_blk_offset(blockid)%ncols+1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   do k=btofc_blk_offset(blockid)%nlvls+1,ldim
      do i=1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   return
   end subroutine block_to_chunk_send_pters
!
!========================================================================

   subroutine block_to_chunk_recv_pters(lchunkid, fdim, ldim, &
                                        record_size, pter)
!----------------------------------------------------------------------- 
! 
! Purpose: Return pointers into receive buffer where data for
!          decomposed chunk data strctures should be copied from
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in) :: lchunkid     ! local chunk id
   integer, intent(in) :: fdim         ! first dimension of pter array
   integer, intent(in) :: ldim         ! last dimension of pter array
   integer, intent(in) :: record_size  ! per coordinate amount of data 

   integer, intent(out) :: pter(fdim,ldim)  ! buffer offset
!---------------------------Local workspace-----------------------------
   integer :: i, k                     ! loop indices
!-----------------------------------------------------------------------
   if ((btofc_chk_offset(lchunkid)%ncols > fdim) .or. &
       (btofc_chk_offset(lchunkid)%nlvls > ldim)) then
      write(6,*) "BLOCK_TO_CHUNK_RECV_PTERS: pter array dimensions ", &
                 "not large enough: (",fdim,",",ldim,") not >= (", &
                  btofc_chk_offset(lchunkid)%ncols,",", &
                  btofc_chk_offset(lchunkid)%nlvls,")"
      call endrun()
   endif
!
   do k=1,btofc_chk_offset(lchunkid)%nlvls
      do i=1,btofc_chk_offset(lchunkid)%ncols
         pter(i,k) = 1 + record_size* &
                     (btofc_chk_offset(lchunkid)%pter(i,k))
      enddo
      do i=btofc_chk_offset(lchunkid)%ncols+1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   do k=btofc_chk_offset(lchunkid)%nlvls+1,ldim
      do i=1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   return
   end subroutine block_to_chunk_recv_pters
!
!========================================================================

   subroutine transpose_chunk_to_block(record_size, &
                                       chunk_buffer, block_buffer)
!----------------------------------------------------------------------- 
! 
! Purpose: Transpose buffer containing decomposed 
!          chunk data structures to buffer
!          containing decomposed longitude/latitude fields 
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, masterproc
#if ( defined SPMD )
   use swap_comm
# if defined(MODCM_DP_TRANSPOSE)
   use mod_comm, only: numpro, blockdescriptor, mp_sendirr, mp_recvirr,  &
                       get_partneroffset
# endif
#endif
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtag  = 1000
!------------------------------Arguments--------------------------------
   integer, intent(in) :: record_size  ! per column amount of data 
   real(r8), intent(out):: chunk_buffer(record_size*chunk_buf_nrecs)
                                       ! buffer of chunk data to be
                                       ! transposed

   real(r8), intent(out) :: block_buffer(record_size*block_buf_nrecs)
                                       ! buffer of block data to
                                       ! transpose into

!---------------------------Local workspace-----------------------------
#if ( defined SPMD )
   integer :: sdispls(0:npes-1)        ! send displacements
   integer :: sndcnts(0:npes-1)        ! send counts
   integer :: rdispls(0:npes-1)        ! receive displacements
   integer :: rcvcnts(0:npes-1)        ! receive counts
   integer :: offset_s                 ! send displacement + 1
   integer :: offset_r                 ! receive displacement + 1
   integer :: i, p                     ! loop indices
   integer :: procid                   ! swap processor id
   integer :: step                     ! step in alltoall alg.
   integer :: sndids(dp_coup_steps)    ! nonblocking MPI send request ids
   integer :: rcvids(dp_coup_steps)    ! nonblocking MPI recv request ids
# if defined(MODCM_DP_TRANSPOSE)
   type (blockdescriptor), allocatable, save :: sendbl(:), recvbl(:)
# endif
   integer ione, ierror, mod_method
   integer first_time_through
   data first_time_through / 0 /
#endif
!-----------------------------------------------------------------------
#if ( defined SPMD )
   sdispls(0) = 0
   sndcnts(0) = record_size*btofc_chk_num(0)
   do p=1,npes-1
     sdispls(p) = sdispls(p-1) + sndcnts(p-1)
     sndcnts(p) = record_size*btofc_chk_num(p)
   enddo
!
   rdispls(0) = 0
   rcvcnts(0) = record_size*btofc_blk_num(0)
   do p=1,npes-1
     rdispls(p) = rdispls(p-1) + rcvcnts(p-1)
     rcvcnts(p) = record_size*btofc_blk_num(p)
   enddo

#if ( defined TIMING_BARRIERS )
   call t_startf ('sync_tran_ctob')
   call mpibarrier (mpicom)
   call t_stopf ('sync_tran_ctob')
#endif

   if (phys_alltoall .eq. 0) then
!
      call mpialltoallv(chunk_buffer, sndcnts, sdispls, mpir8, &
                        block_buffer, rcvcnts, rdispls, mpir8, &
                        mpicom)
!
   elseif (phys_alltoall .eq. 1) then
!
! Post receive requests.
      call swap1m(dp_coup_steps, msgtag, dp_coup_proc, rcvcnts, & 
                  rdispls, record_size*block_buf_nrecs, block_buffer, rcvids)
!     
! Copy local data to new location
      if (sndcnts(iam) > 0) then
         do i=1,sndcnts(iam)
            block_buffer(rdispls(iam)+i) = chunk_buffer(sdispls(iam)+i)
         enddo
      end if
!
! Post send requests and wait for receive requests to complete.
      do step=1,dp_coup_steps
         procid = dp_coup_proc(step)
         offset_s = sdispls(procid)+1
         offset_r = rdispls(procid)+1
         call swap2(msgtag, procid, &
                    sndcnts(procid), chunk_buffer(offset_s), &
                    sndids(step), rcvcnts(procid), &
                    block_buffer(offset_r), rcvids(step))
      enddo
!
! Wait for send requests to complete.
      call swap3m(dp_coup_steps, msgtag, dp_coup_proc, sndids, rcvcnts, & 
                  rdispls, record_size*block_buf_nrecs, block_buffer, rcvids)
!
# if defined(MODCM_DP_TRANSPOSE)
   elseif (phys_alltoall .ge. modmin_alltoall) then
!
! This branch uses mod_comm. Admissable values of phys_alltoall are 11,12,13 and 14. Each value corresponds
!   to a differerent option within mod_comm of implementing the communication. That option is expressed
!   internally to mod_comm using the variable mod_method defined below; mod_method will have values 0,1,2
!   or 3 and is defined as phys_alltoall - modmin_alltoall, where modmin_alltoall equals 11.
!
      if (first_time_through .eq. 0) then
         first_time_through = 1
         mod_method = phys_alltoall - modmin_alltoall
         ione = 1
         allocate( sendbl(0:numpro-1) )
         allocate( recvbl(0:numpro-1) )

         do p = 0,numpro-1

            sendbl(p)%method = mod_method
            recvbl(p)%method = mod_method

            allocate( sendbl(p)%blocksizes(1) )
            allocate( sendbl(p)%displacements(1) )
            allocate( recvbl(p)%blocksizes(1) )
            allocate( recvbl(p)%displacements(1) )

            if ( sndcnts(p) .ne. 0 ) then

               call MPI_TYPE_INDEXED(ione, sndcnts(p),   &
                   sdispls(p), mpir8, &
                   sendbl(p)%type, ierror)
               call MPI_TYPE_COMMIT(sendbl(p)%type, ierror)

               sendbl(p)%blocksizes(1) = sndcnts(p)
               sendbl(p)%displacements(1) = sdispls(p)
               sendbl(p)%partneroffset = 0 ! not properly initialized - do not use Mpi2
            else

               sendbl(p)%type = MPI_DATATYPE_NULL

               sendbl(p)%blocksizes(1) = 0
               sendbl(p)%displacements(1) = 0
               sendbl(p)%partneroffset = 0

            endif

            if ( rcvcnts(p) .ne. 0) then

               call MPI_TYPE_INDEXED(ione, rcvcnts(p),   &
                    rdispls(p), mpir8, &
                    recvbl(p)%type, ierror)
               call MPI_TYPE_COMMIT(recvbl(p)%type, ierror)

               recvbl(p)%blocksizes(1) = rcvcnts(p)
               recvbl(p)%displacements(1) = rdispls(p)
               recvbl(p)%partneroffset = 0 ! not properly initialized - do not use Mpi2
            else

               recvbl(p)%type = MPI_DATATYPE_NULL

               recvbl(p)%blocksizes(1) = 0
               recvbl(p)%displacements(1) = 0
               recvbl(p)%partneroffset = 0

            endif

         enddo

         call get_partneroffset(sendbl, recvbl)

      endif
!
      call mp_sendirr(chunk_buffer, sendbl, recvbl, block_buffer)
      call mp_recvirr(block_buffer, recvbl)
# endif

   endif
!
#endif

   return
   end subroutine transpose_chunk_to_block
!
!========================================================================

   subroutine chunk_to_block_send_pters(lchunkid, fdim, ldim, &
                                        record_size, pter)
!----------------------------------------------------------------------- 
! 
! Purpose: Return pointers into send buffer where data for
!          decomposed chunk data strctures should be copied to
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in) :: lchunkid     ! local chunk id
   integer, intent(in) :: fdim         ! first dimension of pter array
   integer, intent(in) :: ldim         ! last dimension of pter array
   integer, intent(in) :: record_size  ! per coordinate amount of data 

   integer, intent(out) :: pter(fdim,ldim)  ! buffer offset
!---------------------------Local workspace-----------------------------
   integer :: i, k                     ! loop indices
!-----------------------------------------------------------------------
   if ((btofc_chk_offset(lchunkid)%ncols > fdim) .or. &
       (btofc_chk_offset(lchunkid)%nlvls > ldim)) then
      write(6,*) "CHUNK_TO_BLOCK_SEND_PTERS: pter array dimensions ", &
                 "not large enough: (",fdim,",",ldim,") not >= (", &
                  btofc_chk_offset(lchunkid)%ncols,",", &
                  btofc_chk_offset(lchunkid)%nlvls,")"
      call endrun()
   endif
!
   do k=1,btofc_chk_offset(lchunkid)%nlvls
      do i=1,btofc_chk_offset(lchunkid)%ncols
         pter(i,k) = 1 + record_size* &
                     (btofc_chk_offset(lchunkid)%pter(i,k))
      enddo
      do i=btofc_chk_offset(lchunkid)%ncols+1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   do k=btofc_chk_offset(lchunkid)%nlvls+1,ldim
      do i=1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   return
   end subroutine chunk_to_block_send_pters
!
!========================================================================

   subroutine chunk_to_block_recv_pters(blockid, fdim, ldim, &
                                        record_size, pter)
!----------------------------------------------------------------------- 
! 
! Purpose: Return pointers into receive buffer where column from decomposed 
!          longitude/latitude fields should be copied from
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   integer, intent(in) :: blockid      ! block index
   integer, intent(in) :: fdim         ! first dimension of pter array
   integer, intent(in) :: ldim         ! last dimension of pter array
   integer, intent(in) :: record_size  ! per coordinate amount of data 

   integer, intent(out) :: pter(fdim,ldim)  ! buffer offsets
!---------------------------Local workspace-----------------------------
   integer :: i, k                     ! loop indices
!-----------------------------------------------------------------------
   if ((btofc_blk_offset(blockid)%ncols > fdim) .or. &
       (btofc_blk_offset(blockid)%nlvls > ldim)) then
      write(6,*) "CHUNK_TO_BLOCK_RECV_PTERS: pter array dimensions ", &
                 "not large enough: (",fdim,",",ldim,") not >= (", &
                  btofc_blk_offset(blockid)%ncols,",", &
                  btofc_blk_offset(blockid)%nlvls,")"
      call endrun()
   endif
!
   do k=1,btofc_blk_offset(blockid)%nlvls
      do i=1,btofc_blk_offset(blockid)%ncols
         pter(i,k) = 1 + record_size* &
                     (btofc_blk_offset(blockid)%pter(i,k))
      enddo
      do i=btofc_blk_offset(blockid)%ncols+1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   do k=btofc_blk_offset(blockid)%nlvls+1,ldim
      do i=1,fdim
         pter(i,k) = -1
      enddo
   enddo
!
   return
   end subroutine chunk_to_block_recv_pters
!
!========================================================================

   subroutine create_chunks(opt, chunks_per_thread)
!----------------------------------------------------------------------- 
! 
! Purpose: Decompose physics computational grid into chunks, for
!          improved serial efficiency and parallel load balance.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, plev
   use dyn_grid, only: get_block_coord_cnt_d, get_block_coord_d, &
                       get_block_col_cnt_d, &
                       get_lon_d, get_lat_d, get_block_bounds_d, &
                       get_block_owner_d
   use pmgrid, only: plond, platd
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: opt           ! chunking option
      !  0: chunks may cross block boundaries, but retain same
      !     processor mapping as blocks. If possible, processors assigned
      !     as day/night pairs. Columns (or pairs) are wrap-mapped.
      !     May not work with vertically decomposed blocks. (default)
      !  1: chunks may cross block boundaries, but retain same
      !     SMP-node mapping as blocks.  If possible, processors assigned
      !     as day/night pairs.  Columns (or pairs) are wrap-mapped.
      !     May not work with vertically decomposed blocks.
      !  2: 2-column day/night and season column pairs wrap-mapped
      !     to chunks to also balance assignment of polar, mid-latitude, 
      !     and equatorial columns across  chunks.
      !  3: same as 1 except that SMP defined to be pairs of consecutive
      !     processors
      !  4: chunks may cross block boundaries, but retain same
      !     processor mapping as blocks. Columns assigned to chunks
      !     in block ordering.
      !     May not work with vertically decomposed blocks.
      !  5: Chunks do not cross  latitude boundaries, and are block-mapped.
   integer, intent(in)  :: chunks_per_thread 
                                         ! target number of chunks per
                                         !  thread
!---------------------------Local workspace-----------------------------
   integer :: i, j, p                    ! loop indices
   integer :: nlthreads                  ! number of local OpenMP threads
   integer :: npthreads(0:npes-1)        ! number of OpenMP threads per process
   integer :: proc_smp_mapx(0:npes-1)    ! process/virtual SMP map
   integer :: block_cnt                  ! number of blocks containing data
                                         ! for a given vertical column
   integer :: blockids(plev+1)           ! block indices
   integer :: bcids(plev+1)              ! block column indices
   integer :: col_smp_mapx(plond,platd)  ! column/virtual SMP map
   integer :: nsmpx                      ! virtual SMP count 
   integer :: nsmpthreads(0:npes-1)      ! number of OpenMP threads per SMP
   integer :: nsmpcolumns(0:npes-1)      ! number of columns assigned to
                                         !  a given SMP
   integer :: nsmpchunks(0:npes-1)       ! number of chunks assigned to 
                                         !  a given process
   integer :: maxcol_chk(0:npes-1)       ! maximum number of columns assigned 
                                         !  to a chunk in a given SMP
   integer :: smp                        ! SMP index
   integer :: cid                        ! chunk id
   integer :: cid_offset(0:npes-1)       ! chunk id processor offset
   integer :: local_cid(0:npes-1)        ! processor-local chunk id
   integer :: jb, ib                     ! global block and columns indices
   integer :: blksiz                     ! current block size
   integer :: ntmp1, ntmp2               ! work variables
   integer :: firstblock, lastblock      ! global block index bounds
   integer :: lon, lat, twinlon, twinlat ! longitude and latitude indices
   integer :: cbeg                       ! beginning longitude index for 
                                         !  current chunk

#if ( defined _OPENMP )
   integer omp_get_max_threads
   external omp_get_max_threads
#endif
!-----------------------------------------------------------------------
!
! determine number of threads per MPI process
!
   nlthreads = 1
#if ( defined _OPENMP )
   nlthreads = OMP_GET_MAX_THREADS()
#endif
!
#if ( defined SPMD )
   call mpiallgatherint(nlthreads, 1, npthreads, 1, mpicom)
#else
   npthreads(0) = nlthreads
   proc_smp_map(0) = 0
#endif
!
! determine index range for dynamics blocks
!
   call get_block_bounds_d(firstblock,lastblock)
!
! Determine virtual SMP count and process/virtual SMP map.
!  If option 0 or >3, pretend that each SMP has only one processor. 
!  If option 1, use SMP information.
!  If option 2, pretend that all processors are in one SMP node. 
!  If option 3, pretend that each SMP node is made up of two 
!     processes, chosen to maximize load-balancing opportunities.
!
   if (opt == 0) then
      nsmpx = npes
      do p=0,npes-1
         proc_smp_mapx(p) = p
      enddo
   elseif (opt == 1) then
      nsmpx = nsmps
      do p=0,npes-1
         proc_smp_mapx(p) = proc_smp_map(p)
      enddo
   elseif (opt == 2) then
      nsmpx = 1
      do p=0,npes-1
         proc_smp_mapx(p) = 0
      enddo
   elseif (opt == 3) then
      call find_partners(opt,nsmpx,proc_smp_mapx)
   else
      nsmpx = npes
      do p=0,npes-1
         proc_smp_mapx(p) = p
      enddo
   endif
!
! Determine number of columns assigned to each
! SMP in block decomposition
!
   do j=1,platd
      do i=1,plond
         col_smp_mapx(i,j) = -1
      enddo
   enddo
!
   do j=1,plat
      do i=1,nlon_p(j)
         block_cnt = get_block_coord_cnt_d(i,j)
         call get_block_coord_d(i,j,block_cnt,blockids,bcids)
         do jb=1,block_cnt
            p = get_block_owner_d(blockids(jb)) 
            if (col_smp_mapx(i,j) .eq. -1) then
               col_smp_mapx(i,j) = proc_smp_mapx(p)
            elseif (col_smp_mapx(i,j) .ne. proc_smp_mapx(p)) then
               write(6,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
                  "but vertical decomposition not limited to virtual SMP"
               call endrun()
            endif
         enddo
      end do
   end do
!
   nsmpcolumns(:) = 0
   do j=1,plat
      do i=1,nlon_p(j)
         smp = col_smp_mapx(i,j)
         nsmpcolumns(smp) = nsmpcolumns(smp) + 1
      end do
   end do
!
! Options 0-3: split local longitude/latitude blocks into chunks,
!              using wrap-map assignment of columns and
!              day/night and north/south column pairs
!              to chunks to improve load balance
!  Option 0: local is per process
!  Option 1: local is subset of`processes assigned to same SMP node
!  Option 2: local is global
!  Option 3: local is pair of processes chosen to maximize load-balance
!            wrt restriction that only communicate with one other
!            process.
! Option 4: split local longitude/latitude blocks into chunks,
!           using block-map assignment of columns
!             
   if ((opt >= 0) .and. (opt <= 4)) then
!
! Calculate number of threads available in each SMP node. 
!
      nsmpthreads(:) = 0
      do p=0,npes-1
         smp = proc_smp_mapx(p)
         nsmpthreads(smp) = nsmpthreads(smp) + npthreads(p)
      enddo
!
! Determine number of chunks to keep all threads busy
!
      nchunks = 0
      do smp=0,nsmpx-1
         nsmpchunks(smp) = nsmpcolumns(smp)/pcols
         if (mod(nsmpcolumns(smp), pcols) .ne. 0) then
            nsmpchunks(smp) = nsmpchunks(smp) + 1
         endif
         if (nsmpchunks(smp) < chunks_per_thread*nsmpthreads(smp)) then
            nsmpchunks(smp) = chunks_per_thread*nsmpthreads(smp)
         endif
         do while (mod(nsmpchunks(smp), nsmpthreads(smp)) .ne. 0)
            nsmpchunks(smp) = nsmpchunks(smp) + 1
         enddo
         if (nsmpchunks(smp) > nsmpcolumns(smp)) then
            nsmpchunks(smp) = nsmpcolumns(smp)
         endif
         nchunks = nchunks + nsmpchunks(smp)
      enddo
!
! Determine maximum number of columns to assign to chunks
! in a given SMP
!
      do smp=0,nsmpx-1
         ntmp1 = nsmpcolumns(smp)/nsmpchunks(smp)
         ntmp2 = mod(nsmpcolumns(smp),nsmpchunks(smp))
         if (ntmp2 > 0) then
            maxcol_chk(smp) = ntmp1 + 1
         else
            maxcol_chk(smp) = ntmp1
         endif
      enddo
!
! Allocate chunks and knuhcs data structures
!
      allocate ( chunks(1:nchunks) )
      allocate ( knuhcs(1:plond, 1:platd) )
!
! Initialize chunks and knuhcs data structures
!
      do cid=1,nchunks
         chunks(cid)%ncols = 0
      enddo
!
      do j=1,platd
         do i=1,plond
            knuhcs(i,j)%chunkid = -1
         enddo
      enddo
!
! Determine chunk id ranges for each SMP
!
      cid_offset(0) = 1
      local_cid(0) = 0
      do smp=1,nsmpx-1
         cid_offset(smp) = cid_offset(smp-1) + nsmpchunks(smp-1)
         local_cid(smp) = 0
      enddo
!
! Assign columns to chunks
!
      do jb=firstblock,lastblock
         p = get_block_owner_d(jb)
         smp = proc_smp_mapx(p)
         blksiz = get_block_col_cnt_d(jb)
         do ib = 1,blksiz
!
! Assign column to a chunk if not already assigned
            lon = get_lon_d(jb,ib)
            lat = get_lat_d(jb,ib)
            if (knuhcs(lon,lat)%chunkid == -1) then
!
! Find next chunk with space
               cid = cid_offset(smp) + local_cid(smp)
               do while (chunks(cid)%ncols >=  maxcol_chk(smp))
                  local_cid(smp) = mod(local_cid(smp)+1,nsmpchunks(smp))
                  cid = cid_offset(smp) + local_cid(smp)
               enddo
               chunks(cid)%ncols = chunks(cid)%ncols + 1
!
               i = chunks(cid)%ncols
               chunks(cid)%lon(i) = lon
               chunks(cid)%lat(i) = lat
               knuhcs(chunks(cid)%lon(i),chunks(cid)%lat(i))%chunkid = cid
               knuhcs(chunks(cid)%lon(i),chunks(cid)%lat(i))%col = i
!
               if (opt < 4) then
!
! If space available, look to assign a load-balancing "twin" to same chunk
                  if (chunks(cid)%ncols <  maxcol_chk(smp)) then

                     call find_twin(lon, lat, smp, &
                                    proc_smp_mapx, twinlon, twinlat)

                     if ((twinlon > 0) .and. (twinlat > 0)) then
                        chunks(cid)%ncols = chunks(cid)%ncols + 1
!
                        i = chunks(cid)%ncols
                        chunks(cid)%lon(i) = twinlon
                        chunks(cid)%lat(i) = twinlat
                        knuhcs(chunks(cid)%lon(i),chunks(cid)%lat(i))%chunkid = cid
                        knuhcs(chunks(cid)%lon(i),chunks(cid)%lat(i))%col = i
                     endif
!
                  endif
!
! Move on to next chunk (wrap map)
                  local_cid(smp) = mod(local_cid(smp)+1,nsmpchunks(smp))
!
               endif
!
            endif
         enddo
      enddo
!
   else
!
! Option 5: split individual longitude/latitude blocks into chunks,
!            assigning consecutive columns to the same chunk
!
! Determine total number of chunks and maximum block size
!  (assuming no vertical decomposition)
      nchunks = 0
      do j=firstblock,lastblock
         blksiz = get_block_col_cnt_d(j)
         nchunks = nchunks + blksiz/pcols
         if (pcols*(blksiz/pcols) /= blksiz) then
            nchunks = nchunks + 1
         endif
      enddo
!
! Allocate chunks and knuhcs data structures
!
      allocate ( chunks(1:nchunks) )
      allocate ( knuhcs(1:plond, 1:platd) )
!
! Initialize chunks and knuhcs data structures
!
      cid = 0
      do j=firstblock,lastblock
         cbeg = 1
         blksiz = get_block_col_cnt_d(j)
         do while (cbeg <= blksiz)
            cid = cid + 1
            chunks(cid)%ncols = min(pcols,blksiz-(cbeg-1))
            do i=1,chunks(cid)%ncols
               chunks(cid)%lon(i) = get_lon_d(j,i+(cbeg-1))
               chunks(cid)%lat(i) = get_lat_d(j,i+(cbeg-1))
               knuhcs(chunks(cid)%lon(i),chunks(cid)%lat(i))%chunkid = cid
               knuhcs(chunks(cid)%lon(i),chunks(cid)%lat(i))%col = i
            enddo
            cbeg = cbeg + chunks(cid)%ncols
         enddo
      enddo
!
   endif
!
   if ((opt >= 0) .and. (opt <= 3)) then
!
! Arrange columns in chunks in approximate lat, lon order
!
      call sort_chunks()
   endif
!
! Assign chunks to processes.
!
   call assign_chunks(opt, npthreads, &
                      nsmpx, proc_smp_mapx, &
                      nsmpthreads, nsmpchunks)
!
   return
   end subroutine create_chunks
!
!========================================================================

   subroutine find_partners(opt, nsmpx, proc_smp_mapx)
!----------------------------------------------------------------------- 
! 
! Purpose: Divide processes into pairs, attempting to maximize the
!          the number of columns in one process whose twins are in the 
!          other process.
! 
! Method: The day/night and north/south hemisphere complement is defined
!         to be the column twin.
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use dyn_grid, only: get_block_coord_cnt_d, get_block_coord_d, &
                       get_block_owner_d
   use pmgrid, only: plev, masterproc

!------------------------------Arguments--------------------------------
   integer, intent(in)  :: opt           ! chunking option
   integer, intent(out) :: nsmpx         ! calculated number of SMP nodes
   integer, intent(out) :: proc_smp_mapx(0:npes-1)
                                         ! proc/virtual smp map
!---------------------------Local workspace-----------------------------
   integer :: i, j, twini, twinj         ! longitude and latitude indices
   integer :: block_cnt                  ! number of blocks containing data
                                         ! for a given vertical column
   integer :: blockids(plev+1)           ! block indices
   integer :: bcids(plev+1)              ! block column indices
   integer :: jb                         ! block index
   integer :: p, twp                     ! process indices
   integer :: col_proc_mapx(plon,plat)   ! location of columns in 
                                         !  dynamics decomposition
   integer :: twin_proc_mapx(plon,plat)  ! location of column twins in 
                                         !  dynamics decomposition
   integer :: twin_cnt(0:npes-1)         ! for each process, number of twins 
                                         !  in each of the other processes
   logical :: assigned(0:npes-1)         ! flag indicating whether process
                                         !  assigned to an SMP node yet
   integer :: maxpartner, maxcnt         ! process with maximum number of 
                                         !  twins and this count
!-----------------------------------------------------------------------
!
! Determine process location of column and its twin in dynamics decomposition
!
   col_proc_mapx(:,:) = -1
   twin_proc_mapx(:,:) = -1
   do j=1,plat
      do i=1,nlon_p(j)
!
         twinj = (plat+1-j)
         twini = mod((i-1)+(nlon_p(j)/2), nlon_p(j)) + 1
         twini = (nlon_p(twinj)*twini)/nlon_p(j)
         if (twini < 1) twini = 1
         if (twini > nlon_p(twinj)) twini = nlon_p(twinj)
!
         block_cnt = get_block_coord_cnt_d(i,j)
         call get_block_coord_d(i,j,block_cnt,blockids,bcids)

         do jb=1,block_cnt
            p = get_block_owner_d(blockids(jb)) 
            if (col_proc_mapx(i,j) .eq. -1) then
               col_proc_mapx(i,j) = p
            elseif (col_proc_mapx(i,j) .ne. p) then
               if (masterproc) then
                  write(6,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
                     "but vertical decomposition not limited to single process"
               endif
               call endrun()
            endif
         enddo

         block_cnt = get_block_coord_cnt_d(twini,twinj)
         call get_block_coord_d(twini,twinj,block_cnt,blockids,bcids)

         do jb=1,block_cnt
            p = get_block_owner_d(blockids(jb)) 
            if (twin_proc_mapx(i,j) .eq. -1) then
               twin_proc_mapx(i,j) = p
            elseif (twin_proc_mapx(i,j) .ne. p) then
               if (masterproc) then
                  write(6,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
                     "but vertical decomposition not limited to single process"
               endif
               call endrun()
            endif
         enddo

      end do
   end do
!
! Assign process pairs to SMPs, attempting to maximize the number of column,twin
! pairs in same SMP.
!
   assigned(:) = .false.
   twin_cnt(:) = 0
   nsmpx = 0
   do p=0,npes-1
      if (.not. assigned(p)) then
!
! For each process, determine number of twins in each of the other processes
! (running over all columns multiple times to minimize memory requirements).
!
         do j=1,plat
            do i=1,nlon_p(j)
               if (col_proc_mapx(i,j) .eq. p) then
                  twin_cnt(twin_proc_mapx(i,j)) = &
                     twin_cnt(twin_proc_mapx(i,j)) + 1
               endif
            enddo
         enddo
!
! Find process with maximum number of twins which has not yet been designated
! a partner.
!
         maxpartner = -1
         maxcnt = 0
         do twp=0,npes-1
            if ((.not. assigned(twp)) .and. (twp .ne. p)) then
               if (twin_cnt(twp) >= maxcnt) then
                  maxcnt = twin_cnt(twp)
                  maxpartner = twp
               endif
            endif
         enddo
!
! Assign p and twp to the same SMP node
!
         if (maxpartner .ne. -1) then
            assigned(p) = .true.
            assigned(maxpartner) = .true.
            proc_smp_mapx(p) = nsmpx
            proc_smp_mapx(maxpartner) = nsmpx
            nsmpx = nsmpx + 1
         else
            if (masterproc) then
               write(6,*) "PHYS_GRID_INIT error: opt", opt, "specified, ", &
                  "but could not divide processes into pairs."
            endif
            call endrun()
         endif
!
      endif
!      
   enddo
!
   return
   end subroutine find_partners
!
!========================================================================

   subroutine find_twin(lon, lat, smp, &
                        proc_smp_mapx, twinlon_f, twinlat_f)
!----------------------------------------------------------------------- 
! 
! Purpose: Find column that when paired with (lon,lat) in a chunk
!          balances the load. A column is a candidate to be paired with
!          (lon,lat) if it is in the same SMP node as (lon,lat) as defined
!          by proc_smp_mapx.
! 
! Method: The day/night and north/south hemisphere complement is
!         tried first. If it is not a candidate or if it has already been
!         assigned, then the day/night complement is tried next. If that
!         also is not available, then nothing is returned.
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use dyn_grid, only: get_block_coord_d, get_block_owner_d

!------------------------------Arguments--------------------------------
   integer, intent(in)  :: lon, lat      ! global indices for column
                                         ! seeking a twin for
   integer, intent(in)  :: smp           ! index of SMP node (lon,lat)
                                         ! currently assigned to
   integer, intent(in)  :: proc_smp_mapx(0:npes-1)
                                         ! proc/virtual smp map
   integer, intent(out) :: twinlon_f, twinlat_f
                                         ! global indices for twin
!---------------------------Local workspace-----------------------------
   logical :: found                      ! found an acceptable twin
   integer :: twinlon, twinlat           ! indices of twin candidate
   integer :: jbtwin(npes)               ! global block indices
   integer :: ibtwin(npes)               ! global column indices
   integer :: twinproc, twinsmp          ! process and smp ids
!-----------------------------------------------------------------------
   twinlon_f = -1
   twinlat_f = -1
   found = .false.
!
! Try day/night and north/south hemisphere complement first
   twinlon = mod((lon-1)+(nlon_p(lat)/2), nlon_p(lat)) + 1
   twinlat = (plat+1-lat)

   call get_block_coord_d(twinlon,twinlat,npes,jbtwin,ibtwin)
   twinproc = get_block_owner_d(jbtwin(1))
   twinsmp  = proc_smp_mapx(twinproc)
!
   if ((twinsmp .eq. smp) .and. &
       (knuhcs(twinlon,twinlat)%chunkid == -1)) then
      found = .true.
      twinlon_f = twinlon
      twinlat_f = twinlat
   endif
!
! Try day/night complement next
   if (.not. found) then
      twinlon = mod((lon-1)+(nlon_p(lat)/2), nlon_p(lat)) + 1
      twinlat = lat
!
      call get_block_coord_d(twinlon,twinlat,npes,jbtwin,ibtwin)
      twinproc = get_block_owner_d(jbtwin(1))
      twinsmp  = proc_smp_mapx(twinproc)
!
      if ((twinsmp .eq. smp) .and. &
          (knuhcs(twinlon,twinlat)%chunkid == -1)) then
         found = .true.
         twinlon_f = twinlon
         twinlat_f = twinlat
      endif
!
   endif
!
   return
   end subroutine find_twin
!
!========================================================================

   subroutine sort_chunks()
!----------------------------------------------------------------------- 
! 
! Purpose: Arrange columns within each chunk into lat,lon ordering
!          for a range of latitudes, in an attempt to put day columns
!          together.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------

!---------------------------Local workspace-----------------------------
   integer :: i, j, col                  ! loop indices
   integer :: jblksiz                    ! number of latitudes to group
                                         ! together in lat,lon ordering
   logical :: done                       ! reassignment completion flag
   integer :: cid                        ! chunk id
   integer :: ncols                      ! number of columns in chunk
   integer :: cid_map(plon,plat)         ! mapping of global longitude/
                                         !  latitude coordinates within a 
                                         !  chunk
   integer :: col_map(plon,plat)         ! mapping of global longitude/
                                         !  latitude coordinates within a 
                                         !  chunk
   integer :: colmax(nchunks)            ! number of columns in current 
                                         !  chunk
   integer :: colstart(nchunks)          ! latitude index for first column
                                         !  in chunk during a given pass
                                         !  through the columns
   integer :: holdlon, holdlat           ! work variables
   integer :: newlon, newlat             !  ""
   integer :: holdf, holdi               !  ""
   integer :: forward(pcols,nchunks)     ! mapping from current column index
                                         !  to new column index
   integer :: inverse(pcols,nchunks)     ! mapping from new column index
                                         !  to current column index
!-----------------------------------------------------------------------
!
! Initialize column map array
!
   do j=1,plat
      do i=1,plon
         col_map(i,j) = -1
      enddo
   enddo
!
! Fill column map array
!
   do cid=1,nchunks
      ncols = chunks(cid)%ncols
      do col=1,ncols
         i = chunks(cid)%lon(col)
         j = chunks(cid)%lat(col) 
         col_map(i,j) = col
         cid_map(i,j) = cid
      enddo
   enddo
! 
! Determine where columns should move to implement lat,lon ordering
!
#if (defined SCAM)
   jblksiz = 1
#else
   jblksiz = min(plat/npes, plat/4)
#endif
   colmax(:) = 0
   done = .false.
   do while (.not. done)
      colstart(:) = 0
      done = .true.
      do i=1,plon
         do j=1,plat
            if (col_map(i,j) .ne. -1) then
               cid = cid_map(i,j)
               if (colstart(cid) .eq. 0) then
                  colstart(cid) = j
               endif
               if (j .le. colstart(cid)+jblksiz-1) then
                  col = col_map(i,j)
                  colmax(cid) = colmax(cid) + 1
                  inverse(colmax(cid),cid) = col
                  forward(col,cid) = colmax(cid)
                  col_map(i,j) = -1
               else
                  done = .false.
               endif
            endif
         enddo
      enddo
   enddo
!
! Implement rearrangement
!
   do cid=1,nchunks
      do col=1,colmax(cid)
!
! Save current contents of chunks(cid)%xxx(col)
!
         holdlon = chunks(cid)%lon(col)
         holdlat = chunks(cid)%lat(col)
!
! Move column into new location
!
         newlon = chunks(cid)%lon(inverse(col,cid))
         newlat = chunks(cid)%lat(inverse(col,cid))
         chunks(cid)%lon(col) = newlon
         chunks(cid)%lat(col) = newlat
         knuhcs(newlon,newlat)%col = col
!
! Move saved column into old location
!
         chunks(cid)%lon(inverse(col,cid)) = holdlon
         chunks(cid)%lat(inverse(col,cid)) = holdlat
         knuhcs(holdlon,holdlat)%col = inverse(col,cid)
!
! Update forward and inverse functions
!            
         holdi = forward(col,cid)
         holdf = inverse(col,cid)
         forward(holdf,cid) = holdi
         inverse(holdi,cid) = holdf
!
      enddo
!
   enddo
!
   return
   end subroutine sort_chunks
!
!========================================================================

   subroutine assign_chunks(opt, npthreads, &
                            nsmpx, proc_smp_mapx, &
                            nsmpthreads, nsmpchunks)
!----------------------------------------------------------------------- 
! 
! Purpose: Assign chunks to processes.
! 
! Method: 
! 
! Author: Patrick Worley
! 
!-----------------------------------------------------------------------
   use pmgrid, only: iam, plev
   use dyn_grid, only: get_block_coord_cnt_d, get_block_coord_d,&
                       get_block_owner_d 
!------------------------------Arguments--------------------------------
   integer, intent(in)  :: opt           ! mapping option
                                         !  0: keep columns on same processors
                                         !     for both blocks and chunks
                                         !  1: keep columns on same SMP node
                                         !     for both blocks and chunks
                                         !  0-3: 
                                         !     load balance chunks and minimize
                                         !     communication costs in dp_coupling
   integer, intent(in)  :: npthreads(0:npes-1)
                                         ! number of OpenMP threads per process
   integer, intent(in)  :: nsmpx         ! virtual smp count
   integer, intent(in)  :: proc_smp_mapx(0:npes-1)
                                         ! proc/virtual smp map
   integer, intent(in)  :: nsmpthreads(0:npes-1)
                                         ! number of OpenMP threads 
                                         ! per virtual SMP
   integer, intent(in)  :: nsmpchunks(0:npes-1)
                                         ! number of chunks assigned 
                                         ! to a given virtual SMP
!---------------------------Local workspace-----------------------------
   integer :: i, jb, p                   ! loop indices
   integer :: cid                        ! chunk id
   integer :: smp                        ! SMP index
   integer :: glon, glat                 ! global (lon,lat) indices
   integer :: block_cnt                  ! number of blocks containing data
                                         ! for a given vertical column
   integer :: blockids(plev+1)           ! block indices
   integer :: bcids(plev+1)              ! block column indices
   integer :: ntmp1, ntmp2               ! work variables
   integer :: ntmp1_smp(0:npes-1)
   integer :: ntmp2_smp(0:npes-1)
   integer :: npchunks(0:npes-1)         ! number of chunks to be assigned to
                                         !  a given process
   integer :: cur_npchunks(0:npes-1)     ! current number of chunks assigned 
                                         !  to a given process
   integer :: column_count(0:npes-1)     ! number of columns from current chunk
                                         !  assigned to each process in dynamics
                                         !  decomposition
!-----------------------------------------------------------------------
!
! Options 0,4: keep chunks on same processors as corresponding blocks
! Option 1: assign same number of chunks to each thread in the
!           same SMP while minimizing communication costs
! Option 2:
! Option 3: assign same number of chunks to each thread
!           while minimizing communication costs
!
! Determine number of chunks to assign to each process
!
   do smp=0,nsmpx-1
      ntmp1_smp(smp) = nsmpchunks(smp)/nsmpthreads(smp)
      ntmp2_smp(smp) = mod(nsmpchunks(smp),nsmpthreads(smp))
   enddo
!
   do p=0,npes-1
      smp = proc_smp_mapx(p)
      if (ntmp2_smp(smp) > 0) then
         npchunks(p) = ntmp1_smp(smp)*npthreads(p) + 1
         ntmp2_smp(smp) = ntmp2_smp(smp) - 1
      else
         npchunks(p) = ntmp1_smp(smp)*npthreads(p)
      endif
   enddo
!
! Assign chunks to processors: 
!
   do p=0,npes-1
      cur_npchunks(p) = 0
   enddo
!
   do cid=1,nchunks
!
      do p=0,npes-1
         column_count(p) = 0
      enddo
!
!  For each chunk, determine number of columns in each
!  process within the dynamics.
      do i=1,chunks(cid)%ncols
         glon = chunks(cid)%lon(i)
         glat = chunks(cid)%lat(i) 
         block_cnt = get_block_coord_cnt_d(glon,glat)
         call get_block_coord_d(glon,glat,block_cnt,blockids,bcids)
         do jb=1,block_cnt
            p = get_block_owner_d(blockids(jb)) 
            column_count(p) = column_count(p) + 1
         enddo
      enddo
!
!  Eliminate processes that already have their quota of chunks
      do p=0,npes-1
         if (cur_npchunks(p) == npchunks(p)) then
            column_count(p) = -1
         endif
      enddo
!
!  Assign chunk to process with most
!  columns from chunk, from among those still available
      ntmp1 = -1
      ntmp2 = -1
      do p=0,npes-1
         if (column_count(p) > ntmp1) then
            ntmp1 = column_count(p)
            ntmp2 = p
         endif
      enddo
      cur_npchunks(ntmp2) = cur_npchunks(ntmp2) + 1
      chunks(cid)%owner   = ntmp2
!
   enddo
!
   return
   end subroutine assign_chunks
!
!========================================================================

!#######################################################################

end module phys_grid
