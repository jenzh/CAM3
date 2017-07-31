!BOP
!
! !MODULE: mod_comm --- SPMD parallel decompostion/communication module
      module mod_comm
!
! !DESCRIPTION:
!
!  \paragraph{Overview}
!
!    This module contains SPMD parallelism decomposition and
!    communication routines.  This library was originally written by
!    W. Putman and S.-J. Lin for simple gridded communications in the
!    Finite-Volume General Circulation Model (FVGCM).  Most of the
!    member functions are specific to the type of gridded data, 
!    ghost communication and decompositions used in FVGCM (which 
!    are, however, very common in atmospheric models).
!
!    The module was extended for irregular communication by W. Sawyer
!    and A. Mirin, and subsequently optimized in many ways.  It is now
!    a more general tool and has been incorporated into the Parallel
!    Library for Grid Manipulations (PILGRIM) which is used in the
!    Community Atmospheric Model (CAM) and The Physical-space
!    Statistical Analysis System (PSAS).
!
!    Irregular communication is based on the {\tt blockdescriptor}
!    derived type, which defines a set of chunks which are to be
!    send to (or received from) another PE.  The irregular 
!    communication routines operate on arrays of block descriptors
!    whose length is equal to number of PEs involved in the
!    communication.  This means the irregular communication primitives
!    are merely non-blocking all-to-all primitives.
! 
!    Primitives using MPI1 and MPI2 have been implemented, SHMEM is
!    partially implemented.  The module can also make use of OpenMP,
!    and certain communications (MPI2) can be multitasked.  
!
!
!  \paragraph{Use of Global Arrays}
!
!    The module uses the concept of global arrays (coined from former
!    usage of shared memory arenas in the "multi-level parallelism"
!    (MLP) paradigm).  Global arrays are merely buffers into which
!    data are packed for the transfer to other PEs.
!
!    All global arrays are 1-dimensional, they are
!    accessed as needed inside the Ga_Put/Ga_Get routines with offset
!    vars.  (Ga_Put/Ga_Get routines are all 4d with openmp on the 3rd
!    (k) dim.)
!
!  \paragraph{Use of MPI-2 Windows}
!
!       All implementations use real*8, real*4, and integer*4 windows
!       which are used with global arrays as follows:
!
!       \begin{itemize}
!         \item   r8\_win -> ga\_r8 - for use with real*8 types
!         \item   i4\_win -> ga\_i4 - for use with integer*4 types
!       \end{itemize}
!
!       note: MPI routines need 2 buffers per GA, ga\_<type>\_s & ga\_<type>\_r
!             ga\_<type>\_r is used for the MPI2 windows
!
!  \paragraph{Compilation}
!
!    This module contains numerous optimizations for various platforms
!    and underlying communication primitives.   To take advantage of
!    these, the various CPP tokens can be defined.
!
!    \begin{itemize}
!      \item {\tt STAND_ALONE}:  Use as stand-alone library (if
!                                defined) or as part of CAM (if 
!                                undefined)
!      \item {\tt MODCM_TIMING}: Turn on CAM timing routines (only
!                                available if compiled in CAM framework)
!      \item {\tt USE\_MPI2}:    Use MPI-2 (one-sided communication)
!                                for underlying communication.
!      \item {\tt USE\_SHMEM}:   Use SHMEM (one-sided communication)
!                                for underlying communication.
!      \item {\tt USE\_VT}:      Initialize and finalize MPI
!                                internally (when used in absense
!                                of other libraries)
!      \item {\tt MODCM\_ALLOC}: Dynamic allocation of buffers in
!                                module (if defined), static 
!                                allocation at run-time (when undefined);
!                                for MPI2 uses mpi_alloc_mem.
!      \item {\tt MT\_OFF}:      Explicitly turn off multitasking
!      \item {\tt SET\_CPUS}:    Set number of threads (for multitasking)
!                                from the environment variable
!                                \verb$AGCM_N_THREADS_PER_PROCESS$
!      \item {\tt PIN_CPUS}:     (SGI specific) map tasks to hardware
!                                cpus for performance
!      \item {\tt IRIX64}:       Compilation specific to SGI IRIX64
!      \item {\tt OSF1}:         Compilation specific to Compaq (HP) OSF1
!      \item {\tt AIX}:          Compilation specific to IBM AIX
!      \item {\tt Linux}:        Compilation specific to Linux
!      \item {\tt _OPENMP}:      Implicit token (controlled by
!                                compiler) to enable OpenMP
!    \end{itemize}
!
!    
!  \paragraph{Usage}
!
!    NOTE - must call PILGRIM routine parinit to initialize before
!    making any other calls (unless {\tt USE\_VT} defined).   
!    SHMEM (IRIX64) supported for border communications but 
!    not irregular communications. IRIX64 -specific
!    commands not recently tested.
!
!    The public members of this module are:
!
!      \begin{itemize}
!         \item {\tt mp\_init}:          Initialize module
!         \item {\tt mp\_exit}:          Exit module
!         \item {\tt mp\_send4d\_ns}:    Ghost 4D array on north/south
!         \item {\tt mp\_recv4d\_ns}:    Complete 4D N/S ghost operation
!         \item {\tt mp\_send2\_ns}:     Ghost 2 3D arrays on north/south
!         \item {\tt mp\_recv2\_ns}:     Complete 2x3D N/S ghost operation
!         \item {\tt mp\_send3d\_2}:     Send 2x3D general ghost regions
!         \item {\tt mp\_recv3d\_2}:     Complete 2x3D general ghost operation
!         \item {\tt mp\_sendirr}:       Initiate all-to-all send of chunks
!         \item {\tt mp\_recvirr}:       Complete all-to-all chunk commun.
!         \item {\tt mp\_sendirr\_i4}:   Initiate all-to-all send of
!                                        chunks (integer*4)
!         \item {\tt mp\_recvirr\_i4}:   Complete all-to-all chunk
!                                        commun. (integer*4)
!         \item {\tt y\_decomp}:         Generate YZ decomposition (CAM)
!         \item {\tt set\_decomp}:       Set YZ decomposition (CAM)
!       \end{itemize}
!
!     There are other public routines, but these are only used internally
!     in PILGRIM, and they should not be called by user applications.
!
! !REVISION HISTORY:
!    2001.09.01   Lin
!    2002.04.16   Putman  Modified for Global Array code
!    2002.04.16   Putman  Added ProTeX documentation
!    2002.05.28   Putman  Added use of precision module
!    2003.06.24   Sawyer  Minor additions for use with mod_irreg
!    2004.01.08   Sawyer  Removed older functionality, no longer needed
!    2004.02.10   Mirin   Major restructuring and simplification. Documentation
!    2004.03.06   Sawyer  Additional documentation; cosmetics
! !USES:
!
#include "pilgrim.h"
!
! Mod_comm has option for stand-alone use as well as within CAM
!
#if !defined( STAND_ALONE )
#include "misc.h"
#endif

#if defined ( SPMD )

#if defined( STAND_ALONE )
#define r8 selected_real_kind(12)
#define r4 selected_real_kind( 6)
#define i8 selected_int_kind(13)
#define i4 selected_int_kind( 6)
#else
      use shr_kind_mod, only : r8 => shr_kind_r8, r4 => shr_kind_r4,  &
                               i8 => shr_kind_i8, i4 => shr_kind_i4
#endif

      implicit none

!
! Shmem option presently implemented for regular communications only
!
#if defined ( USE_SHMEM )
#if defined(USE_VT)
#include "mpif.h"
#endif
#if defined(IRIX64)
#include "mpp/shmem.fh"
#else
#include "shmem.fh"
#endif
#else
#include "mpif.h"
#endif

#if defined(STAND_ALONE)
#  define PLON        144
#  define PLAT         91
#  define PLEV         55
#  define PCNST         1
#  define PNATS         0
#else
#include "params.h"
#endif
 
! !PUBLIC MEMBER FUNCTIONS:
      public mp_init, mp_exit, y_decomp, set_decomp,               &
             mp_send4d_ns, mp_recv4d_ns, mp_send2_ns, mp_recv2_ns, &
             mp_send3d, mp_recv3d, mp_send3d_2, mp_recv3d_2,       &
             mp_sendirr, mp_recvirr, mp_sendirr_i4, mp_recvirr_i4, &
             mp_barrier, mp_r8, mp_r4, mp_i4
#if !defined(USE_MPI2)
      public sqest, rqest, nsend, nread
#endif

!------------------------------------------------------------------------------
!  type declaration for describing an arbitrary number of contiguous chunks
!  this is for irregular communications
!------------------------------------------------------------------------------
      type blockdescriptor
         integer              :: method             ! transpose method
         integer              :: type               ! Ptr to MPI derived type
         integer, pointer     :: displacements(:)   ! Offsets in local segment
         integer, pointer     :: blocksizes(:)      ! Block sizes to transfer
         integer              :: partneroffset      ! Aggregated partner offset
         integer              :: partnertype        ! Ptr to partner's MPI derived type
      end type blockdescriptor
      integer maxwin
      parameter(maxwin=1)   !   Do not change without reading below (AAM, 12/08/03)
!
!   WARNING     WARNING     WARNING     WARNING     WARNING     WARNING
!
! Normal operation is to use dedicated target windows for Mpi2; however, one can use
!    local windows. For that, see immediately below.
!
! Maxwin refers to the number of stored local target windows with Mpi2. For each
!    window, one stores the target array location and window size. If the window entry
!    is already in the stored table, a new window need not be created. This approach
!    assumes that a given variable name is always stored at the same location. That is
!    often not the case and is definitely not the case with Cam/Ccsm, as temporary
!    variables are reallocated each timestep. Hence, one will run out of windows
!    quickly and need to create new windows anyway. What's worse, is that for one
!    MPI task a target variable might be at the same location as before, and for
!    another task at a different location. That will cause some tasks to attempt
!    to create a new window while others do not. This leads to deadlock and is fatal.
!    If you wish to run with maxwin greater than 1, be sure to negate the following
!    line in mp_sendirr, which forces nwinhit to 0:
!       nwinhit = 0  ! force window recreation to avoid partial hit deadlock
!
      integer :: winlocal(maxwin)
      integer*8 :: wintable(maxwin,2)
      integer*8 winaddr
#if defined (USE_MPI2)
      integer (kind=MPI_ADDRESS_KIND) winsize
#endif
      integer :: nwininit=0
      integer nwinused

! transpose methods (method)
!      0 for contiguous temporary buffer
!      (MPI-1) 1 for direct communication (derived types)
!      (MPI-2) >1 or 2 for direct put into contiguous temporary target window
!                 1 for OpenMP over contiguous segments destined for target task
!                 2 for OpenMP over target tasks, using derived types for each target
!               3 for direct put into local target window coincident with target array 
!                 uses derived types for source and target

      INTEGER, SAVE :: InHandle(MAX_PAX, MAX_TRF)
      INTEGER, SAVE :: OutHandle(MAX_PAX, MAX_TRF)
      INTEGER, SAVE :: BegTrf = 0  ! Ongoing overlapped begintransfer #
      INTEGER, SAVE :: EndTrf = 0  ! Ongoing overlapped endtransfer #

! !PUBLIC DATA MEMBERS:
      integer, SAVE:: gid                         ! PE id
      integer(i4), SAVE:: numpro                  ! Permanent No. of PEs
      integer(i4), SAVE:: numcpu                  ! No. of threads
      integer, SAVE:: commglobal                  ! Global Communicator

!------------------------------------------------------------------------------
!  Local parameters for use with MPI-2 and MPI-1
!------------------------------------------------------------------------------
      integer, parameter:: nbuf = 2               ! Max No. of sends per call
      integer, parameter:: nghost = 3             ! No. of ghost indices
      integer, parameter:: max_nq = PCNST + PNATS ! No. of tracers
      integer, parameter:: max_call = 2           ! Max No. of back-to-back...
                                                  ! ...mp_send calls

#if defined(USE_SHMEM)
      integer, parameter:: idimsize = PLON*nghost*(PLEV+1)*max_nq
      integer, parameter :: mp_r4 = r4
      integer, parameter :: mp_r8 = r8
      integer, parameter :: mp_i4 = i4
      integer, SAVE::  reduce_sync(SHMEM_REDUCE_SYNC_SIZE)
      integer, SAVE:: collect_sync(SHMEM_COLLECT_SYNC_SIZE)
      integer, SAVE::   bcast_sync(SHMEM_BCAST_SYNC_SIZE)
#else
      integer, parameter:: idimsize = PLON*nghost*(PLEV+1)*max_nq
                                                  ! Size of MPI buffer region
                                                  ! in mp_send/mp_recv calls, used
                                                  ! to determine offset in GA
      integer, parameter :: mp_r4 = MPI_REAL
      integer, parameter :: mp_r8 = MPI_DOUBLE_PRECISION
      integer, parameter :: mp_i4 = MPI_INTEGER
#endif

!------------------------------------------------------------------------------
!  Local variables for use with MPI-2, MPI-1 and SHMEM
!------------------------------------------------------------------------------
      integer, SAVE:: gsize                       ! No. of PEs
      integer, SAVE:: nowpro                      ! Temp. PE id
      integer, SAVE:: np_loop                     ! No. of sends for bcasts 
                                                  !     np_loop = numpro for MPI
      integer, allocatable, SAVE:: yfirst(:)      ! First latitude
      integer, allocatable, SAVE:: ylast(:)       ! Last latitude
      integer, allocatable, SAVE:: zfirst(:)      ! First level
      integer, allocatable, SAVE:: zlast(:)       ! Last level
      integer, SAVE:: ncall_r, ncall_s

      integer, SAVE:: sizet1, sizer8, sizei4, tracmax, dpvarmax, totvar
      integer, SAVE:: tracertrans_mod, phys_transpose_mod, phys_transpose_modmin
      integer, SAVE:: phys_transpose_vars, idimsizz
      data tracertrans_mod / -1 /
      data phys_transpose_mod / -1 /
      data phys_transpose_modmin / 11 /
      data phys_transpose_vars / 7 /
!
! tracertrans_mod is maximum number of tracers to be simultaneously transposed; it is communicated
!     from CAM.
! phys_transpose_mod is the communication method for dynamics/physics transposes; admissable values
!     are >= phys_transpose_modmin; it is communicated from CAM.
! phys_transpose_vars is the number of non-tracer variables transposed between dynamics and
!     physics instanciations in CAM.

!------------------------------------------------------------------------------
!  Variables to control global array locations and window synchronization
!------------------------------------------------------------------------------
      integer win_count                           ! Counts No. of windows in use
      integer lastwin                             ! ID of last synch'd window
      integer pkgs_per_pro                        ! No. of MPI packages per PE
      integer igosouth, igonorth                  ! Index of send direction
      integer ifromsouth, ifromnorth              ! Index of recv direction

!------------------------------------------------------------------------------
!  Local type declaration for mp_windows
!------------------------------------------------------------------------------
      type window
         integer :: id            ! Window id
         integer :: size          ! Size of global window (point based)
         integer :: ncall_s       ! Count send calls on window
         integer :: ncall_r       ! Count recv calls on window
#if defined(USE_MPI2)
#if defined(LINUX)
#define MPI_OFFSET_KIND 8
#define MPI_ADDRESS_KIND 8
#endif
         integer(kind=MPI_ADDRESS_KIND) :: offset_s
         integer(kind=MPI_ADDRESS_KIND) :: offset_r
#else
         integer :: offset_s      ! Starting position in GA send
         integer :: offset_r      ! Starting position in GA recv
#endif
         integer :: dest          ! For use with send calls
         integer :: src           ! For use with recv calls
         integer :: size_r        ! Size of incoming message
     end type window

!------------------------------------------------------------------------------
! Beginning Global Array variable declaration:
!------------------------------------------------------------------------------

      type (window) :: r8_win
      type (window) :: i4_win


!
!  SHMEM variable declarations
!

#if defined(USE_SHMEM)
      integer ga_ptr

#define ga_t1_r ga_r8_r
#define ga_t1_s ga_r8_s
#define t1_win  r8_win

      real(r8), SAVE, TARGET:: ga_r8a_r(MAX( PLON*PLAT*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))
      real(r8), SAVE, TARGET:: ga_r8a_s(MAX( PLON*PLAT*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))

      real(r8), SAVE, TARGET:: ga_r8b_r(MAX( PLON*PLAT*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))
      real(r8), SAVE, TARGET:: ga_r8b_s(MAX( PLON*PLAT*(PLEV+1)*max_nq, idimsize*nbuf*max_call ))

      real(r8), DIMENSION(:), SAVE, POINTER:: ga_r8_r
      real(r8), DIMENSION(:), SAVE, POINTER:: ga_r8_s

      integer(i4), SAVE, TARGET:: ga_i4a_r(PLON*PLAT*PLEV)
      integer(i4), SAVE, TARGET:: ga_i4a_s(PLON*PLAT*PLEV)

      integer(i4), SAVE:: ga_i4c_r(PLON*PLAT*PLEV)
      integer(i4), SAVE:: ga_i4c_s(PLON*PLAT*PLEV)

      integer(i4), SAVE, TARGET:: ga_i4b_r(PLON*PLAT*PLEV)
      integer(i4), SAVE, TARGET:: ga_i4b_s(PLON*PLAT*PLEV)

      integer(i4), DIMENSION(:), SAVE, POINTER:: ga_i4_r
      integer(i4), DIMENSION(:), SAVE, POINTER:: ga_i4_s

      integer ierror
#endif

!
!  MPI-1 and MPI-2 window variable declarations
!

#if !defined(USE_SHMEM)
      type (window) :: t1_win
# if defined(MODCM_ALLOC)
#  if defined(USE_MPI2)
      real(r8) ga_t1_r
      real(r8) ga_t1_s
      real(r8) ga_r8_r
      real(r8) ga_r8_s
      integer(i4) ga_i4_r
      integer(i4) ga_i4_s
!
! Cray pointers required for mpi_alloc_mem
!
      pointer (pa_t1_r, ga_t1_r(1))
      pointer (pa_t1_s, ga_t1_s(1))
      pointer (pa_r8_r, ga_r8_r(1))
      pointer (pa_r8_s, ga_r8_s(1))
      pointer (pa_i4_r, ga_i4_r(1))
      pointer (pa_i4_s, ga_i4_s(1))
      save pa_t1_r
      save pa_t1_s
      save pa_r8_r
      save pa_r8_s
      save pa_i4_r
      save pa_i4_s
#  else
      real(r8),allocatable,    SAVE:: ga_t1_r(:)
      real(r8),allocatable,    SAVE:: ga_t1_s(:)
      real(r8),allocatable,    SAVE:: ga_r8_r(:)
      real(r8),allocatable,    SAVE:: ga_r8_s(:)
      integer(i4),allocatable, SAVE:: ga_i4_r(:)
      integer(i4),allocatable, SAVE:: ga_i4_s(:)
#  endif
# else
      real(r8),    SAVE:: ga_t1_r(idimsize*nbuf*max_call)
      real(r8),    SAVE:: ga_t1_s(idimsize*nbuf*max_call)
      real(r8),    SAVE:: ga_r8_r(PLON*PLAT*(PLEV+1)*max_nq)
      real(r8),    SAVE:: ga_r8_s(PLON*PLAT*(PLEV+1)*max_nq)
      integer(i4), SAVE:: ga_i4_r(PLON*PLAT*PLEV)
      integer(i4), SAVE:: ga_i4_s(PLON*PLAT*PLEV)
# endif
#endif

!
!  MPI-2 variable declarations
!

#if defined(USE_MPI2)

!------------------------------------------------------------------------------
!  The lines immediately below assume that MPI-2 is not implemented
!  on the specified architectures. This faulty assumption is historical
!  and must be modified for machines that do support MPI-2.
!
#if defined(LINUX) || defined(OSF1)
#if defined(LINUX)
#define MPI_ADDRESS_KIND 8
#endif
      integer, parameter:: MPI_MODE_NOCHECK    = 0 
      integer, parameter:: MPI_MODE_NOSTORE    = 0    
      integer, parameter:: MPI_MODE_NOPUT      = 0 
      integer, parameter:: MPI_MODE_NOPRECEDE  = 0 
      integer, parameter:: MPI_MODE_NOSUCCEED  = 0 
#endif
!------------------------------------------------------------------------------

      integer(kind=MPI_ADDRESS_KIND) bsize
      integer, SAVE:: Status(MPI_STATUS_SIZE)
      integer ierror
#endif

!
!  MPI-1 variable declarations
!

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      integer, SAVE:: nsend                   ! Number of messages out-going
      integer, SAVE:: nrecv                   ! Number of messages in-coming
      integer, SAVE:: nread                   ! Number of messages read
      integer, allocatable, SAVE:: sqest(:)   ! Request handler for sends
      integer, allocatable, SAVE:: rqest(:)   ! Request handler for recvs
      integer, SAVE:: bsize             
      integer, SAVE:: Status(MPI_STATUS_SIZE)
      integer, allocatable, SAVE:: Stats(:)
      integer ierror
#endif

#if defined (SEMA)
      integer semid
#endif

!
!EOP
!------------------------------------------------------------------------------
      contains
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_init --- Initialize SPMD parallel communication
!
! !INTERFACE:
      subroutine mp_init( comm )
!
! !INPUT PARAMETERS:
      integer, optional :: comm
! !DESCRIPTION:
!
!     Initialize SPMD parallel communication
!
! !REVISION HISTORY: 
!    2001.09.01   Lin
!    2002.02.15   Putman	Modified for Global Array code
!    2002.04.09   Putman	Added ProTeX documentation
!    2002.08.06   Sawyer        Added optional communicator input argument
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mysize
        integer firsttime
        data firsttime / 0 /
#if defined(USE_SHMEM)
        integer num_pes, my_pe
#endif
#if defined(USE_MPI2)
        integer (KIND=MPI_ADDRESS_KIND) sizet18
        integer (KIND=MPI_ADDRESS_KIND) sizer88
        integer (KIND=MPI_ADDRESS_KIND) sizei48
#endif

        win_count = 0

#if defined(USE_SHMEM)
# if defined(IRIX64)
        CALL START_PES(0)
# else
        call SHMEM_INIT()
# endif
# if defined(USE_VT)
        call MPI_INIT(ierror)
# endif
        numpro = num_pes()
        gid = my_pe()
	ga_ptr = 1
	ga_r8_r => ga_r8a_r
	ga_r8_s => ga_r8a_s
	ga_i4_r => ga_i4a_r
	ga_i4_s => ga_i4a_s
	reduce_sync = SHMEM_SYNC_VALUE
	collect_sync = SHMEM_SYNC_VALUE
	bcast_sync = SHMEM_SYNC_VALUE
        call omp_start
#else
        if ( present(comm) ) then
          call mpi_start( comm )
        else
          call mpi_start( MPI_COMM_WORLD )
        endif
        call omp_start
# if !defined(USE_MPI2)
        allocate( sqest(MAX(nbuf,numpro)*max_call) )
        allocate( rqest(MAX(nbuf,numpro)*max_call) )
        allocate( Stats(MAX(nbuf,numpro)*max_call*MPI_STATUS_SIZE) )
# endif
#endif

      idimsizz = idimsize
#if defined(MODCM_ALLOC)
!
! Dynamically allocate target global arrays
!
      if (firsttime .eq. 0) then
         firsttime = 1
! NOTE: tracertrans_mod computed in spmd_dyn_setopts; it is the maximum number of tracers
!    that can be simultaneously transposed
         if (tracertrans_mod .eq. -1) then
            if (gid .eq. 0) print *,       &
              '(MOD_COMM) - tracertrans_mod not properly initialized - ignoring'
            tracmax = max_nq
         else
            tracmax = min(max_nq, tracertrans_mod)
         endif

         totvar = tracmax
! NOTE: phys_transpose_mod computed in phys_grid_setopts
         if (phys_transpose_mod .eq. -1) then
            if (gid .eq. 0) print *,       &
              '(MOD_COMM) - phys_transpose_mod not properly initialized - ignoring'
            dpvarmax = 0
!
! If phys_transpose_mod is >= phys_transpose_modmin, that is a signal that mod_comm is to be used
!    for dynamics/physics transposes in cam/ccsm. In that case, one must allocate enough window
!    storage for those transposes.
!    Presently, the number of transposed variables equals phys_transpose_vars plus the number of
!    constituents.
!
         elseif (phys_transpose_mod .ge. phys_transpose_modmin) then
            dpvarmax = phys_transpose_vars + max_nq
         else
            dpvarmax = 0
         endif
         totvar = max(dpvarmax, totvar)
!
! totvar is the maximum of the number of tracers to be simultaneously transposed and the number of
!    variables to be transposed between dynamics and physics instanciations in CAM
         sizet1 = nghost*nbuf*max_call*PLON*(PLEV+1)*tracmax
         idimsizz = (idimsize/max_nq)*tracmax
! Compute local storage requirement by dividing global requirement by the number of tasks.
!    Allow factor-of-2 slack to account for nonuniformity of decomposition.
         sizer8 = 2*ceiling( real(PLON*PLAT*(PLEV+1)*totvar)/real(numpro) )
         sizei4 = 1    !   i4 operations temporarily disabled
# if defined(USE_MPI2)
         sizet18 = 8*sizet1
         sizer88 = 8*sizer8
         sizei48 = 8*sizei4
         call mpi_alloc_mem(sizet18, mpi_info_null, pa_t1_r, ierror)
         if (ierror .ne. 0) print *, 'MPI_ALLOC_MEM error'
         call mpi_alloc_mem(sizet18, mpi_info_null, pa_t1_s, ierror)
         if (ierror .ne. 0) print *, 'MPI_ALLOC_MEM error'
         call mpi_alloc_mem(sizer88, mpi_info_null, pa_r8_r, ierror)
         if (ierror .ne. 0) print *, 'MPI_ALLOC_MEM error'
         call mpi_alloc_mem(sizer88, mpi_info_null, pa_r8_s, ierror)
         if (ierror .ne. 0) print *, 'MPI_ALLOC_MEM error'
         call mpi_alloc_mem(sizei48, mpi_info_null, pa_i4_r, ierror)
         if (ierror .ne. 0) print *, 'MPI_ALLOC_MEM error'
         call mpi_alloc_mem(sizei48, mpi_info_null, pa_i4_s, ierror)
         if (ierror .ne. 0) print *, 'MPI_ALLOC_MEM error'
# else
         allocate( ga_t1_r(sizet1) )
         allocate( ga_t1_s(sizet1) )
         allocate( ga_r8_r(sizer8) )
         allocate( ga_r8_s(sizer8) )
         allocate( ga_i4_r(sizei4) )
         allocate( ga_i4_s(sizei4) )
# endif
      endif
#endif

#if !defined(USE_SHMEM)
# if defined(MODCM_ALLOC)
        mysize = sizet1
# else
        mysize = idimsize*nbuf*max_call
# endif
        call win_init_r8(t1_win, ga_t1_r, mysize)
#endif

#if defined(MODCM_ALLOC)
        mysize = sizer8
#else
        mysize = PLON*PLAT*(PLEV+1)*max_nq
#endif
        call win_init_r8(r8_win, ga_r8_r, mysize)

#if defined(MODCM_ALLOC)
        mysize = sizei4
#else
        mysize = PLON*PLAT*PLEV
#endif
        call win_init_i4(i4_win, ga_i4_r, mysize)

        igosouth   = 0
        igonorth   = 1
        ifromsouth = 1
        ifromnorth = 0

        ncall_s = 0
        ncall_r = 0
        lastwin = r8_win%id

        np_loop = numpro

        allocate( yfirst( numpro ) )
        allocate( ylast( numpro ) )
        allocate( zfirst( numpro ) )
        allocate( zlast( numpro ) )

#if defined(USE_SHMEM)
        print *, 'Using Shmem communications in mod_comm'
#elif defined(USE_MPI2)
        if (gid .eq. 0) print *, 'Using Mpi2 communications in mod_comm'
#else
        if (gid .eq. 0) print *, 'Using Mpi1 communications in mod_comm'
#endif
!EOC
      end subroutine mp_init
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_exit --- End SPMD parallel communication
!
! !INTERFACE:
      subroutine mp_exit
! !DESCRIPTION:
!
!     End SPMD parallel communication
!
! !REVISION HISTORY: 
!    2001.09.01   Lin
!    2002.02.15   Putman	Modified for Global Array code
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

        deallocate( yfirst )
        deallocate( ylast )
        deallocate( zfirst )
        deallocate( zlast )
#if !defined(USE_SHMEM)
#if defined(USE_MPI2)
        call MPI_WIN_FREE( t1_win%id, ierror )
        call MPI_WIN_FREE( r8_win%id, ierror )
        call MPI_WIN_FREE( i4_win%id, ierror )
#endif
        call MPI_FINALIZE (ierror)
#elif defined(USE_SHMEM)
        call mp_barrier()
#if defined(USE_VT)
    	call MPI_FINALIZE (ierror)
#endif
#endif
        return
!EOC
      end subroutine mp_exit
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: omp_start --- Start openMP parallelism
!
! !INTERFACE:
      subroutine omp_start
! !DESCRIPTION:
!
!     Start openMP parallelism
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer ios, n, nowpro, nowcpu

!  Accommodate different means of specifying the number of OpenMP threads
!  See comments at beginning of module
!
#if defined(SET_CPUS)

        character*80 evalue
        call getenv('AGCM_N_THREADS_PER_PROCESS',evalue)
        if (gid == 0) then
          read(evalue,*,iostat=ios) numcpu
          if ( ios .ne. 0 ) then
              print *, 'ERROR: cannot read AGCM_N_THREADS_PER_PROCESS', &
                       trim(evalue)
               call exit(1)
          end if
        endif
        call mp_bcst_int(numcpu)

# if defined(_OPENMP)

#  if defined (IRIX64)
       call mp_set_numthreads(numcpu)  !keep it for a while, :)
#  else
       call omp_set_num_threads(numcpu)
#  endif

# else
       print *, 'ERROR: must turn on OPENMP to set numcpu in mod_comm'
       call exit(1)
# endif

#else 

# if defined(_OPENMP)

#  if defined (IRIX64)
        integer mp_suggested_numthreads
        numcpu = mp_suggested_numthreads(0)
#  else
        integer omp_get_num_threads
!$omp parallel
        numcpu = omp_get_num_threads()
!$omp end parallel
#  endif

# else
        numcpu = 1
# endif

#endif

#if defined(USE_MPI2) || defined(USE_SHMEM)
# if defined(MT_OFF)
        pkgs_per_pro = 1
# else
        pkgs_per_pro = numcpu
# endif
#endif

#if defined(PIN_CPUS)
!$omp parallel do private(n,nowcpu)
        nowpro = gid
        do n=1,numcpu
          nowcpu = n + (nowpro) * numcpu-1
          call mp_assign_to_cpu(nowcpu)
        enddo
#endif
!EOC
      end subroutine omp_start
!------------------------------------------------------------------------------

#if !defined(USE_SHMEM)
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mpi_start --- Start MPI parallelism
!
! !INTERFACE:
      subroutine mpi_start( comm )
! !INPUT PARAMETERS:
      integer :: comm   ! Global communicator (may be MPI_COMM_WORLD)
! !DESCRIPTION:
!
!     Start MPI parallelism
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!    02.08.06   Sawyer  Added communicator input arguments
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        logical flag
        integer npthreads

#if defined(USE_MPI2) && !defined(AIX) && (!defined LINUX) && (!defined OSF1)
        call MPI_INITIALIZED( flag, ierror )
        if ( .not. flag ) then
          call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE, npthreads, ierror)
        endif
        call MPI_QUERY_THREAD(npthreads, ierror)
#if !defined(MT_OFF)
        if (npthreads /= MPI_THREAD_MULTIPLE) then
          write(*,*) gid, 'did not provide MPI_THREAD_MULTIPLE. ', &
                'Change to MPI_THREAD_MULTIPLE with MPI_INIT_THREAD ', &
                'for multi-threading MPI2'
          call MPI_FINALIZE(ierror)
          call exit(1)
        endif
#endif
#else
        call MPI_INITIALIZED( flag, ierror )
        if ( .not. flag ) then
          call MPI_INIT( ierror )
        endif
#endif
        call MPI_COMM_RANK (comm, gid, ierror)
        call MPI_COMM_SIZE (comm, numpro, ierror)
        call MPI_COMM_DUP  (comm, commglobal, ierror)
!EOC
      end subroutine mpi_start
!------------------------------------------------------------------------------
#endif

!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_r8 --- Initialize real*8 communication window
!
! !INTERFACE:
      subroutine win_init_r8(win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: isize
        real(r8), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize real*8 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mp_size, info

        win_count = win_count + 1

#if defined(USE_MPI2)
        call MPI_INFO_CREATE(info, ierror)
        call MPI_INFO_SET(info, "no_locks", "true", ierror)
#if defined(AIX)
        info = MPI_INFO_NULL
#endif
        call MPI_TYPE_SIZE(mp_r8, mp_size, ierror)
        bsize = isize*mp_size
        call MPI_WIN_CREATE(ga, bsize, mp_size, info, commglobal, &
                            win%id, ierror)
#if !defined(AIX)
        call MPI_INFO_FREE(info, ierror)
#endif
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
#else
        win%id = win_count
#endif
        win%size = isize
        win%ncall_s = 0
        win%ncall_r = 0
!EOC
      end subroutine win_init_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: win_init_i4 --- Initialize integer*4 communication window
!
! !INTERFACE:
      subroutine win_init_i4(win, ga, isize)
! !INPUT PARAMETERS:
        integer, intent(in) :: isize
        integer(i4), intent(in) :: ga(isize)
! !OUTPUT PARAMETERS:
        type (window), intent(inout) :: win
! !DESCRIPTION:
!
!     Initialize integer*4 communication window
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
        integer mp_size, info

        win_count = win_count + 1

#if defined(USE_MPI2)
        call MPI_INFO_CREATE(info, ierror)
        call MPI_INFO_SET(info, "no_locks", "true", ierror)
#if defined(AIX)
        info = MPI_INFO_NULL
#endif
        call MPI_TYPE_SIZE(mp_i4, mp_size, ierror)
        bsize = isize*mp_size
        call MPI_WIN_CREATE(ga, bsize, mp_size, info, commglobal, &
                            win%id, ierror)
#if !defined(AIX)
        call MPI_INFO_FREE(info, ierror)
#endif
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
#else
        win%id = win_count
#endif
        win%size = isize
        win%ncall_s = 0
        win%ncall_r = 0
!EOC
      end subroutine win_init_i4
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!BOP
! !ROUTINE: y_decomp --- Decompose the Latitude & Level direction
!
! !INTERFACE:
      subroutine y_decomp(jm, km, jfirst, jlast, kfirst, klast, myid)
! !INPUT PARAMETERS:
      integer, intent(in):: jm     ! Dimensions
      integer, intent(in):: km     ! Levels
      integer, intent(in):: myid
! !OUTPUT PARAMETERS:
      integer, intent(inout):: jfirst, jlast, kfirst, klast
! !DESCRIPTION:
!
!     Decompose the Latitude & Level direction for SPMD parallelism
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer p, p1, p2, lats, pleft
      integer, allocatable:: ydist(:)

      if (myid == 0) print *, "numpro", numpro, "numcpu", numcpu
      allocate( ydist( numpro ) )

      lats = jm / numpro
      pleft = jm - lats * numpro

      if( lats < 3 ) then
         write(*,*) 'Number of Proc is too large for jm=',jm
!BMP         call exit(1)
      endif

      do p=1,numpro
         ydist(p) = lats
      enddo

      if ( pleft .ne. 0 ) then
          p1 = (numpro+1) / 2 
          p2 = p1 + 1
        do while ( pleft .ne. 0 )
           if( p1 .eq. 1 ) p1 = numpro
               ydist(p1) = ydist(p1) + 1
               pleft = pleft - 1
               if ( pleft .ne. 0 ) then
                    ydist(p2) = ydist(p2) + 1
                    pleft = pleft - 1
               endif
               p2 = p2 + 1
               p1 = p1 - 1
        enddo
      endif

! Safety check:
      lats = 0
      do p = 1, numpro
         lats = lats + ydist(p)
      enddo

      if ( lats .ne. jm ) then
         print *, "Decomp: big trouble sum(ydist) = ", lats, "!=", jm
      endif
 
      jfirst = 1
      jlast  = ydist(1)
      yfirst(1) = jfirst
      ylast(1) = jlast
      kfirst = 1
      klast = km
      zfirst(1) = kfirst
      zlast(1) = klast

      do p = 1,numpro-1
         yfirst(p+1) = ylast(p) + 1
         ylast(p+1) = ylast(p) + ydist(p+1) 
         if( p == myid ) then
            jfirst = yfirst(p+1)
            jlast  = ylast (p+1)
         endif
         zfirst(p+1) = kfirst
         zlast(p+1) = klast
      enddo

      deallocate (ydist)
!EOC
      end subroutine y_decomp
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: set_decomp --- Set the Latitude & Level decomposition
!
! !INTERFACE:
      subroutine set_decomp(nprocs, jm, km, ydist, zdist)
! !INPUT PARAMETERS:
      integer, intent(in):: nprocs
      integer, intent(in):: jm, km
      integer, intent(in):: ydist(nprocs)
      integer, intent(in):: zdist(nprocs)   ! Currently not used
!
! !DESCRIPTION:
!
!     Set the Latitude & Level decomposition:
!             if it is defined external to mod_comm
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer lats, p

! Safety check:
      lats = 0
      do p = 1, nprocs
         lats = lats + ydist(p)
      enddo

      if ( lats .ne. jm ) then
         print *, "Decomp: big trouble sum(ydist) = ", lats, "!=", jm
      endif

      yfirst(1) = 1
      ylast(1) = ydist(1)
      zfirst(1) = 1
      zlast(1) = km

      do p = 1,nprocs-1
         yfirst(p+1) = ylast(p) + 1
         ylast(p+1) = ylast(p) + ydist(p+1)
         zfirst(p+1) = 1
         zlast(p+1) = km
      enddo
!EOC
      end subroutine set_decomp
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send4d_ns --- Send 4d north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_send4d_ns(im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
      real(r8), intent(in):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Send 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call Win_Open(t1_win)

! Send to south
      if ( jfirst > 1 ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        t1_win%src = gid - 1
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        t1_win%size_r = im*ng_s*(klast-kfirst+1)*nq
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
#endif
        t1_win%dest = gid - 1
        t1_win%offset_s = igosouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8(q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst, jfirst+ng_n-1, kfirst, klast, 1, nq,   &
                         ga_t1_s, ga_t1_r )
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        t1_win%src = gid + 1
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        t1_win%size_r = im*ng_n*(klast-kfirst+1)*nq
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
#endif
        t1_win%dest = gid + 1
        t1_win%offset_s = igonorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8(q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast-ng_s+1, jlast, kfirst, klast, 1, nq,     &
                         ga_t1_s, ga_t1_r )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send4d_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv4d_ns --- Receive 4d north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_recv4d_ns(im, jm, km, nq, jfirst, jlast, kfirst, klast, &
                              ng_s, ng_n, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km, nq
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast
      integer, intent(in):: ng_s      ! southern zones to ghost 
      integer, intent(in):: ng_n      ! northern zones to ghost 
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: q(im,jfirst-ng_s:jlast+ng_n,kfirst:klast,nq)
!
! !DESCRIPTION:
!
!     Receive 4d north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call Win_Close(t1_win)


! Recv from south
      if ( jfirst > 1 ) then
        t1_win%src  = gid-1
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8(q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jfirst-ng_s, jfirst-1,   kfirst, klast, 1, nq, &
                         ga_t1_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        t1_win%src  = gid+1
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8(q, t1_win, im, jm, km, nq, &
                         1, im, jfirst-ng_s, jlast+ng_n, kfirst, klast, 1, nq, &
                         1, im, jlast+1,     jlast+ng_n, kfirst, klast, 1, nq, &
                         ga_t1_r  )
      endif

      call Win_Finalize(t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv4d_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send2_ns --- Send 2 variables north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_send2_ns(im, jm, km, jfirst, jlast, kfirst, klast, &
                             nd, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
      real(r8), intent(in):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real(r8), intent(in):: q2(im,jfirst-nd:jlast+nd,kfirst:klast) 
!
! !DESCRIPTION:
!
!     Send 2 variables north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call Win_Open(t1_win)

! Send to south
      if ( jfirst > 1 ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        t1_win%src  = gid - 1
        t1_win%size_r = im*(klast-kfirst+1)
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
#endif
        t1_win%dest = gid - 1
        t1_win%offset_s = igosouth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, jfirst,    jfirst,   kfirst, klast, 1, 1, &
                          ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + im*(klast-kfirst+1)
        call Ga_Put4d_r8( q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, jfirst,    jfirst,   kfirst, klast, 2, 2, &
                          ga_t1_s, ga_t1_r  )
      endif
! Send to north
      if ( jlast < jm ) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        t1_win%src  = gid + 1
        t1_win%size_r = im*(klast-kfirst+1)
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
#endif
        t1_win%dest = gid + 1
        t1_win%offset_s = igonorth*idimsizz + (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, jlast,     jlast,    kfirst, klast, 1, 1, &
                          ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + im*(klast-kfirst+1)
        call Ga_Put4d_r8( q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, jlast,     jlast,    kfirst, klast, 2, 2, &
                          ga_t1_s, ga_t1_r  )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send2_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv2_ns --- Receive 2 variables north/south ghost latitudes
!
! !INTERFACE:
      subroutine mp_recv2_ns(im, jm, km, jfirst, jlast, kfirst, klast, &
                             nd, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: im, jm, km
      integer, intent(in):: jfirst, jlast
      integer, intent(in):: kfirst, klast !careful: klast might be klast+1
      integer, intent(in):: nd
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: q1(im,jfirst-nd:jlast+nd,kfirst:klast) 
      real(r8), intent(inout):: q2(im,jfirst-nd:jlast+nd,kfirst:klast)
!
! !DESCRIPTION:
!
!     Receive 2 variables north/south ghost latitudes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.02.15   Putman	Modified for Global Arrays code   
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer j

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call Win_Close(t1_win)

! Recv from south
      if ( jfirst > 1 ) then
        j = jfirst - 1
        t1_win%src  = gid - 1
        t1_win%offset_r = ifromsouth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( q1, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, j,         j,        kfirst, klast, 1, 1, &
                          ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_Get4d_r8( q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, j,         j,        kfirst, klast, 2, 2, &
                          ga_t1_r  )
      endif
! Recv from north
      if ( jlast < jm ) then
        j = jlast + 1
        t1_win%src  = gid + 1
        t1_win%offset_r = ifromnorth*idimsizz + (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( q1, t1_win, im, jm, km, 2, & 
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 1, 1, &
                          1, im, j,         j,        kfirst, klast, 1, 1, &
                          ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + im*(klast-kfirst+1)
        call Ga_Get4d_r8( q2, t1_win, im, jm, km, 2, &
                          1, im, jfirst-nd, jlast+nd, kfirst, klast, 2, 2, &
                          1, im, j,         j,        kfirst, klast, 2, 2, &
                          ga_t1_r  )
      endif

      call Win_Finalize(t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv2_ns
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send3d --- Send ghost region
!
! !INTERFACE:
      subroutine mp_send3d(dest, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                                  i1, i2, j1, j2, k1, k2, q)
!
! !INPUT PARAMETERS:
      integer, intent(in):: dest, src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl 
      integer, intent(in):: i1, i2, j1, j2, k1, k2
      real(r8), intent(in):: q(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Send a general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman  
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call Win_Open(t1_win)

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
! Init Recv src
      if ( src >= 0 .and. src < numpro ) then     ! is PE in valid range?
        t1_win%src = src
        t1_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)  ! chunk size
        t1_win%offset_r = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
      endif
#endif
! Send ghost region
      if ( dest >= 0 .and. dest < numpro ) then
        t1_win%dest = dest
        t1_win%offset_s = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( q, t1_win, im, jm, km, 1, &
                          if, il, jf, jl, kf, kl, 1, 1,  &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_s, ga_t1_r  )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv3d --- Recv ghost region
!
! !INTERFACE:
      subroutine mp_recv3d(src, im, jm, km, if, il, jf, jl, kf, kl, &
                                                i1, i2, j1, j2, k1, k2, qout)
!
! !INPUT PARAMETERS:
      integer, intent(in):: src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl 
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qout(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Recv a general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman  
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call Win_Close(t1_win)

! Recv from src
      if ( src >= 0 .and. src < numpro ) then     ! is PE in valid range?
        t1_win%src  = src
        t1_win%offset_r = (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( qout, t1_win, im, jm, km, 1, &
                          if, il, jf, jl, kf, kl, 1, 1,   &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_r  )
      endif

      call Win_Finalize(t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv3d
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_send3d_2 --- Send 2 ghost regions
!
! !INTERFACE:
      subroutine mp_send3d_2(dest, src, im, jm, km, if, il, jf, jl, kf, kl, &
                                        i1, i2, j1, j2, k1, k2, q1, q2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: dest, src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl
      integer, intent(in):: i1, i2, j1, j2, k1, k2
      real(r8), intent(in):: q1(if:il, jf:jl, kf:kl)
      real(r8), intent(in):: q2(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Send two general 3d real*8 ghost region
!
! !REVISION HISTORY:
!    02.04.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call Win_Open(t1_win)

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
! Init Recv src
      if ( src >= 0 .and. src < numpro ) then     ! is PE in valid range?
        t1_win%src = src
        t1_win%size_r = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)  ! chunk size
        t1_win%offset_r = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
        t1_win%offset_r = t1_win%offset_r + t1_win%size_r 
        call Ga_RecvInit_r8(t1_win, ga_t1_r)
      endif
#endif
! Send ghost region
      if ( dest >= 0 .and. dest < numpro ) then
        t1_win%dest = dest
        t1_win%offset_s = (t1_win%ncall_s-1)*idimsizz*nbuf
        call Ga_Put4d_r8( q1, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 1, 1,  &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_s, ga_t1_r  )
        t1_win%offset_s = t1_win%offset_s + (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        call Ga_Put4d_r8( q2, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 2, 2,  &
                          i1, i2, j1, j2, k1, k2, 2, 2, ga_t1_s, ga_t1_r  )
      endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_send3d_2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recv3d_2 --- Recv 2 ghost regions
!
! !INTERFACE:
      subroutine mp_recv3d_2(src, im, jm, km, if, il, jf, jl, kf, kl, &
                                  i1, i2, j1, j2, k1, k2, qout1, qout2)
!
! !INPUT PARAMETERS:
      integer, intent(in):: src
      integer, intent(in):: im, jm, km
      integer, intent(in):: if, il, jf, jl, kf, kl
      integer, intent(in):: i1, i2, j1, j2, k1, k2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: qout1(if:il, jf:jl, kf:kl)
      real(r8), intent(inout):: qout2(if:il, jf:jl, kf:kl)
!
! !DESCRIPTION:
!
!     Recv two general 3d real*8 ghost regions
!
! !REVISION HISTORY:
!    02.04.15   Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      call Win_Close(t1_win)

! Recv from src
      if ( src >= 0 .and. src < numpro ) then     ! is PE in valid range?
        t1_win%src  = src
        t1_win%offset_r = (t1_win%ncall_r-1)*idimsizz*nbuf
        call Ga_Get4d_r8( qout1, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 1, 1,   &
                          i1, i2, j1, j2, k1, k2, 1, 1, ga_t1_r  )
        t1_win%offset_r = t1_win%offset_r + (i2-i1+1)*(j2-j1+1)*(k2-k1+1)
        call Ga_Get4d_r8( qout2, t1_win, im, jm, km, 2, &
                          if, il, jf, jl, kf, kl, 2, 2,   &
                          i1, i2, j1, j2, k1, k2, 2, 2, ga_t1_r  )
      endif

      call Win_Finalize(t1_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recv3d_2
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_barrier --- Synchronize all SPMD processes
!
! !INTERFACE:
      subroutine mp_barrier
!
! !DESCRIPTION:
!
!     Synchronize all SPMD processes
!
! !REVISION HISTORY: 
!    2001.09.01   Lin    
!    2002.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

#if defined (USE_SHMEM)
        call SHMEM_BARRIER_ALL()
#else
	call MPI_BARRIER(commglobal, ierror)
#endif

!EOC
      end subroutine mp_barrier
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Open --- Open a communication window
!
! !INTERFACE:
      subroutine Win_Open(win)
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     Begin a communication epoch, by opening a comm window.
!     Update number of send calls on the window (win%ncall_s).
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      win%ncall_s = win%ncall_s + 1
      ncall_s = ncall_s + 1

#if defined(USE_SHMEM)
      if (ncall_s == 1) then
         if (ga_ptr == 1) then 
            ga_ptr = 2
            ga_r8_r => ga_r8b_r
            ga_r8_s => ga_r8b_s
            ga_r8_r => ga_r8b_r
            ga_r8_s => ga_r8b_s
            ga_i4_r => ga_i4b_r
            ga_i4_s => ga_i4b_s
         else
            ga_ptr = 1
            ga_r8_r => ga_r8a_r
            ga_r8_s => ga_r8a_s
            ga_r8_r => ga_r8a_r
            ga_r8_s => ga_r8a_s
            ga_i4_r => ga_i4a_r
            ga_i4_s => ga_i4a_s
         endif
      endif
#endif

!EOC
      end subroutine Win_Open
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Close --- Close a communication window
!
! !INTERFACE:
      subroutine Win_Close(win)
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     End a communication epoch, by closing a comm window.
!     Update number of receive calls on the window (win%ncall_r).
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      win%ncall_r = win%ncall_r + 1
      ncall_r = ncall_r + 1
#if defined(USE_SHMEM)
      if (ncall_r == 1) then
          call mp_barrier()
      endif
#endif
#if defined(USE_MPI2)
      if (win%ncall_r == 1) then
          call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED, &
                             win%id, ierror)
      endif
#endif

!EOC
      end subroutine Win_Close
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Win_Finalize --- Reset a communication window after a comm epoch.
!
! !INTERFACE:
      subroutine Win_Finalize(win)
! !OUTPUT PARAMETERS:
      type(window), intent(inout):: win
!
! !DESCRIPTION:
!
!     Complete a communication epoch and reset a comm window.
!     Update global lastwin with win%id.
!     Barrier synchronzize if necessary.
!
! !REVISION HISTORY:
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC

      if (ncall_s == ncall_r) then
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
        call MPI_WAITALL(nsend, sqest, Stats, ierror)
        nsend = 0
        nrecv = 0
        nread = 0
#endif
        ncall_s = 0
        ncall_r = 0
      endif

      if (win%ncall_s == win%ncall_r) then
#if defined(USE_MPI2)
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, win%id, ierror)
#endif
        lastwin = win%id
        win%ncall_s = 0
        win%ncall_r = 0
      endif

!EOC
      end subroutine Win_Finalize
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_r8 --- Write to real*8 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_r8 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r8), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: ga_s(win%size)
      real(r8), intent(inout):: ga_r(win%size)
!
! !DESCRIPTION:
!
!     Write to real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      integer send_tag, qsize
#if defined(USE_MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#if defined(USE_SHMEM)
#if defined(OSF1)
      integer(kind=8) p, tmpsize, mysize, mydisp
#else
      integer p, tmpsize, mysize, mydisp
#endif
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length 
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends ( MPI-2,MPI-1 )
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2

                ga_s(inc+i) = q(i,j,k,iq)

              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
#if defined(USE_MPI2) || defined(USE_SHMEM)
      tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
!$omp parallel do private(p,mysize,mydisp,ierror)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = win%offset_s + (p-1)*tmpsize
#if defined(USE_MPI2)
        call MPI_PUT(ga_s(mydisp+1), mysize, mp_r8, &
                     win%dest, mydisp, mysize, mp_r8, &
                     win%id, ierror)
#else
           call SHMEM_PUT64(ga_r(mydisp+1), ga_s(mydisp+1), &
                                 mysize, win%dest)
#endif
      enddo
#else
      send_tag = gid
      nsend = nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_r8, win%dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
#endif

!EOC
      end subroutine Ga_Put4d_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_r8 --- Initiate real*8 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_r8( win, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      real(r8), intent(inout):: ga(win%size)
!
! !DESCRIPTION:
!
!     Initiate real*8 Non-Blocking receive ( MPI-1 ).
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!    03.06.06   Sawyer  Added else clause
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      if (win%size >= win%offset_r + win%size_r) then
        recv_tag = win%src
        qsize    = win%size_r
        nrecv    = nrecv + 1
        call MPI_IRECV(ga(win%offset_r+1), qsize, mp_r8, win%src, &
                       recv_tag, commglobal, rqest(nrecv), ierror)
      else
        print *, "Fatal ga_recvinit_r8: receive window out of space - exiting"
        call exit(1)
      endif
#endif

!EOC
      end subroutine Ga_RecvInit_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_r8 --- Read from real*8 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_r8 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      real(r8), intent(in)  :: ga(win%size)
! !OUTPUT PARAMETERS:
      real(r8), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from real*8 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Put4d_i4 --- Write to integer*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Put4d_i4 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga_s, ga_r )
! !INPUT PARAMETERS:
      type(window), intent(in)  :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      integer(i4), intent(in)  :: q(ifrom:ito,jfrom:jto,kfrom:kto,nqfrom:nqto)
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: ga_s(win%size)
      integer(i4), intent(inout):: ga_r(win%size)
!
! !DESCRIPTION:
!
!     Write to integer*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length

      integer send_tag, qsize
#if defined(USE_MPI2)
      integer p, tmpsize, mysize
      integer(kind=MPI_ADDRESS_KIND) mydisp
#endif
#if defined(USE_SHMEM)
#if defined(OSF1)
      integer(kind=8) p, tmpsize, mysize, mydisp
#else
      integer p, tmpsize, mysize, mydisp
#endif
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1

      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

! Begin Non-Blocking Sends ( MPI-2,MPI-1 )
        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_s) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                ga_s(inc+i) = q(i,j,k,iq)
              enddo
            enddo
          enddo
        enddo

      qsize = (i2-i1+1)*(j2-j1+1)*(k2-k1+1)*(nq2-nq1+1)
#if defined(USE_MPI2) || defined(USE_SHMEM)
      tmpsize = ceiling(real(qsize)/real(pkgs_per_pro))
!$omp parallel do private(p,mysize,mydisp,ierror)
      do p=1,MIN(pkgs_per_pro,qsize)
        mysize = MIN(tmpsize, MAX(qsize-(tmpsize*(p-1)),0))
        mydisp = win%offset_s + (p-1)*tmpsize
#if defined(USE_MPI2)
        call MPI_PUT(ga_s(mydisp+1), mysize, mp_i4, &
                     win%dest, mydisp, mysize, mp_i4, &
                     win%id, ierror)
#else
           call SHMEM_PUT32(ga_r(mydisp+1), ga_s(mydisp+1), &
                                 mysize, win%dest)
#endif
      enddo
#else
      send_tag = gid
      nsend = nsend + 1
      call MPI_ISEND(ga_s(win%offset_s+1), qsize, mp_i4, win%dest, &
                     send_tag, commglobal, sqest(nsend), ierror)
#endif

!EOC
      end subroutine Ga_Put4d_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_RecvInit_i4 --- Initiate integer*4 Non-Blocking receive
!
! !INTERFACE:
      subroutine Ga_RecvInit_i4( win, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout):: ga(win%size)
!
! !DESCRIPTION:
!
!     Initiate integer*4 Non-Blocking receive ( MPI-1 ).
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer qsize, recv_tag

#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      recv_tag = win%src
      qsize    = win%size_r
      nrecv    = nrecv + 1
      call MPI_IRECV(ga(win%offset_r+1), qsize, mp_i4, win%src, &
                     recv_tag, commglobal, rqest(nrecv), ierror)
#endif
!EOC
      end subroutine Ga_RecvInit_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Get4d_i4 --- Read from integer*4 4d global array
!
! !INTERFACE:
      subroutine Ga_Get4d_i4 ( q, win, im, jm, km, nq, &
                                  ifrom, ito, jfrom, jto, kfrom, kto, &
                                  nqfrom, nqto, i1, i2, j1, j2, k1, k2, &
                                  nq1, nq2, ga )
! !INPUT PARAMETERS:
      type(window), intent(in)    :: win           ! Global Array Window
      integer, intent(in)  :: im, jm, km, nq
      integer, intent(in)  :: i1, i2, j1, j2, k1, k2, nq1, nq2
      integer, intent(in)  :: ifrom, ito, jfrom, jto, kfrom, kto, nqfrom, nqto
      integer(i4), intent(in)  :: ga(win%size)
! !OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: q(ifrom:ito, jfrom:jto, kfrom:kto, nqfrom:nqto)
!
! !DESCRIPTION:
!
!     Read from integer*4 4d global array.
!
! !REVISION HISTORY: 
!    02.02.15   Putman
!    02.04.09   Putman	Added ProTeX documentation
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer i, j, k, iq, inc, inc1
      integer i_length, j_length, k_length, ij_length, ijk_length
#if !defined(USE_SHMEM) && !defined(USE_MPI2)
      nread = nread + 1
      call MPI_WAIT(rqest(nread), Status, ierror)
#endif

      i_length   = i2-i1+1
      j_length   = j2-j1+1
      k_length   = k2-k1+1
      ij_length  = i_length*j_length
      ijk_length = i_length*j_length*k_length

        do iq = nq1, nq2
!$omp parallel do private(i,j,k,inc,inc1)
          do k = k1, k2
                inc1 = (win%offset_r) + ((iq-nq1)*ijk_length) &
                       + ((k-k1)*ij_length) -i1+1
            do j = j1, j2
                inc = inc1 + (j-j1)*i_length
              do i = i1, i2
                q(i,j,k,iq) = ga(inc+i)
              enddo
            enddo
          enddo
        enddo

!EOC
      end subroutine Ga_Get4d_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_r8 --- Broadcast an real*8 1d global array
!
! !INTERFACE:
	subroutine Ga_Broadcast_r8 ( q, isize )
! !INPUT PARAMETERS:
	integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
	real(r8), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!	Broadcast an real*8 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
	integer i
#if defined(USE_SHMEM)
#if defined(OSF1)
        integer(kind=8) p, tmpsize, mysize, mydisp
#else
        integer p, tmpsize, mysize, mydisp
#endif
#endif

#if defined(USE_SHMEM)
        mysize = isize
	call mp_barrier
        do i=1,mysize
           ga_r8_s(i) = q(i)
        enddo
	call SHMEM_BROADCAST64(ga_r8_r, ga_r8_s, mysize, 0, 0, 0, numpro, bcast_sync)
	if (gid /= 0) then
	do i=1,isize
           q(i) = ga_r8_r(i)
	enddo
	endif
#endif

#if !defined(USE_SHMEM)
	call MPI_BCAST(q, isize, mp_r8, 0, commglobal, ierror)
#endif

!EOC
      end subroutine Ga_Broadcast_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_Broadcast_i4 --- Broadcast an integer*4 1d global array
!
! !INTERFACE:
	subroutine Ga_Broadcast_i4 ( q, isize )
! !INPUT PARAMETERS:
	integer, intent(in)  :: isize
! !OUTPUT PARAMETERS:
	integer(i4), intent(inout) :: q(isize)
!
! !DESCRIPTION:
!
!	Broadcast an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
	integer i
#if defined(USE_SHMEM)
#if defined(OSF1)
        integer(kind=8) mysize
#else
        integer mysize
#endif
#endif

#if defined(USE_SHMEM)
	mysize = isize
	call mp_barrier
        do i=1,mysize
           ga_i4_s(i) = q(i)
        enddo
	call SHMEM_BROADCAST32(ga_i4_r, ga_i4_s, mysize, 0, 0, 0, numpro, bcast_sync)
	if (gid /= 0) then
	do i=1,isize
           q(i) = ga_i4_r(i)
	enddo
	endif
#endif

#if !defined(USE_SHMEM)
	call MPI_BCAST(q, isize, mp_i4, 0, commglobal, ierror)
#endif

!EOC
      end subroutine Ga_Broadcast_i4
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_r8 --- All to All of an integer*4 1d global array
!
! !INTERFACE:
	subroutine Ga_AllToAll_r8 ( q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
	integer, intent(in)  :: Gsize	! Global size of array
	integer, intent(in)  :: Lsize	! size of Local portion
	integer, intent(in)  :: istart	! starting point
! !OUTPUT PARAMETERS:
	real(r8), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!	All to All of an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
	 integer n, i
#if defined(USE_SHMEM)
#if defined(OSF1)
        integer(kind=8) my_Lsize
#else
        integer my_Lsize
#endif
#endif

#if defined(USE_SHMEM)
	my_Lsize = Lsize
	call mp_barrier
        do i=1,Lsize
           ga_r8_s(i) = q(i+istart-1)
        enddo
	call SHMEM_COLLECT8(ga_r8_r, ga_r8_s, my_Lsize, 0, 0, numpro, collect_sync)
	do i=1,Gsize
	   q(i) = ga_r8_r(i)
	enddo
#endif
#if !defined(USE_SHMEM)
        call MPI_ALLGATHER(q(istart), Lsize, mp_r8, q, Lsize, mp_r8, commglobal, ierror)
#endif

!EOC
	end subroutine Ga_AllToAll_r8
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: Ga_AllToAll_i4 --- All to All of an integer*4 1d global array
!
! !INTERFACE:
	subroutine Ga_AllToAll_i4 ( q, Gsize, Lsize, istart )
! !INPUT PARAMETERS:
	integer, intent(in)  :: Gsize   ! Global size of array
        integer, intent(in)  :: Lsize   ! size of Local portion
        integer, intent(in)  :: istart  ! starting point
! !OUTPUT PARAMETERS:
	integer(i4), intent(inout) :: q(Gsize)
!
! !DESCRIPTION:
!
!	All to All of an integer*4 1d global array.
!
! !REVISION HISTORY:
!    03.04.02	Putman
!
!EOP
!------------------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
	integer n, i
#if defined(USE_SHMEM)
#if defined(OSF1)
        integer(kind=8) my_Lsize
#else
        integer my_Lsize
#endif
#endif

#if defined(USE_SHMEM)
	my_Lsize = Lsize
	call mp_barrier
        do i=1,Lsize
           ga_i4_s(i) = q(i+istart-1)
        enddo
	call SHMEM_COLLECT4(ga_i4_r, ga_i4_s, my_Lsize, 0, 0, numpro, collect_sync)
	do i=1,Gsize
	   q(i) = ga_i4_r(i)
	enddo
#endif
#if !defined(USE_SHMEM)
        call MPI_ALLGATHER(q(istart), Lsize, mp_i4, q, Lsize, mp_i4, commglobal, ierror)
#endif

!EOC
	end subroutine Ga_AllToAll_i4
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: get_partneroffset --- Computes partneroffset/type from descriptor
!
! !INTERFACE:
      subroutine get_partneroffset ( send_bl, recv_bl )

! !INPUT/OUTPUT PARAMETERS:
      type(blockdescriptor), intent(inout)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(inout)  :: recv_bl(:) ! receive blocks

!
! !DESCRIPTION:
!     Compute partneroffsets/types from other blockdescriptor
!     information.  Used exclusively for irregular communication 
!     in PILGRIM.
!
! !REVISION HISTORY: 
!    03.10.31   Mirin       Creation
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

      integer :: i, j, k, ns, pos, por, numpsq, ierror
      integer :: ami(numpro,numpro), am(numpro,numpro)
      integer :: nrblcks(numpro)
      type sizedisp
         integer size
         integer, pointer :: data(:,:)
      end type sizedisp
      type (sizedisp) nrblckdata(numpro)
      integer nstot, nrtot, nn, ind, mod_method, mpi2_flag
      integer, pointer :: sendtot(:,:), recvtot(:,:)
      integer sendcounts(numpro), sdispls(numpro)
      integer recvcounts(numpro), rdispls(numpro)

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif

      do j = 1, numpro
         send_bl(j)%partneroffset = 0
         recv_bl(j)%partneroffset = 0
         send_bl(j)%partnertype = MPI_DATATYPE_NULL
         recv_bl(j)%partnertype = MPI_DATATYPE_NULL
      enddo

    if (mpi2_flag .eq. 1) then

      mod_method = recv_bl(1)%method

     if (mod_method .lt. 3) then

! partner offset (dedicated window)
! am(i,j) is number of words sent from i to j (1-based indices)
! compute global table by using reduction on local parts
      ami(:,:) = 0
      am(:,:) = 0
      i = gid + 1
      numpsq = numpro*numpro
      do j = 1, numpro
         ns = size(send_bl(j)%blocksizes)
         do k = 1, ns
            ami(i,j) = ami(i,j) + send_bl(j)%blocksizes(k)
         enddo
      enddo

      call mpi_allreduce(ami, am, numpsq, mpi_integer, MPI_MAX, commglobal, ierror)

      do j = 1, numpro
         send_bl(j)%partneroffset = 0
         recv_bl(j)%partneroffset = 0
         pos = 0
         por = 0
         do k = 1, i-1
            pos = pos + am(k,j)
            por = por + am(j,k)
         enddo
         if (ami(i,j) .ne. 0) send_bl(j)%partneroffset = pos
         if (ami(j,i) .ne. 0) recv_bl(j)%partneroffset = por
      enddo

     endif  ! mod_method .lt. 3

     if (mod_method .ge. 3) then

! partner receive type (for mpi_put) (local window)
! am(i,j) is the number of blocks sending from i to j (1-based indices)
! compute global table by using reduction on local parts
! prepare to distribute receive block information
      ami(:,:) = 0
      am(:,:) = 0
      j = gid + 1
      nstot = 0
      do i = 1, numpro
         ns = size(recv_bl(i)%blocksizes)
         ami(i,j) = ns
         nstot = nstot + ns
      enddo
! nstot is total number of received contiguous blocks
      sendcounts(:) = 0
      sdispls(:) = 0
      if (nstot .gt. 0) then
         allocate ( sendtot(2,nstot) )
         ind = 0
         do i = 1, numpro
            ns = size(recv_bl(i)%blocksizes)
            do nn = 1, ns
               ind = ind + 1
               sendtot(1,ind) = recv_bl(i)%blocksizes(nn)
               sendtot(2,ind) = recv_bl(i)%displacements(nn)
            enddo
            sendcounts(i) = 2*ns
            if (i .gt. 1) sdispls(i) = sdispls(i-1) + sendcounts(i-1)
         enddo
      else
         allocate ( sendtot(2,1) )
      endif

      numpsq = numpro*numpro

! globalize number of blocks sending from i to j
      call mpi_allreduce(ami, am, numpsq, mpi_integer, MPI_MAX, commglobal, ierror)

! nrblcks is the number of blocks to be received by target "i" from source "j"
      nrtot = 0
      do i = 1, numpro
         nrblcks(i) = am(j,i)
         if (nrblcks(i) .gt. 0) then
            nrtot = nrtot + nrblcks(i)
            allocate ( nrblckdata(i)%data(nrblcks(i),2) )
            nrblckdata(i)%size = nrblcks(i)
         else
            allocate ( nrblckdata(i)%data(1,2) )
            nrblckdata(i)%size = 0
         endif
      enddo
! nrtot is the total number of blocks to be received from source "j"
      if (nrtot .gt. 0) then
         allocate ( recvtot(2,nrtot) )
         recvcounts(:) = 0
         rdispls(:) = 0
         do i = 1, numpro
            recvcounts(i) = 2*nrblcks(i)
            if (i .gt. 1) rdispls(i) = rdispls(i-1) + recvcounts(i-1)
         enddo
      else
         allocate ( recvtot(2,1) )
      endif

! distribute receive block information
      call mpi_alltoallv(sendtot, sendcounts, sdispls, mpi_integer, recvtot,    &
                         recvcounts, rdispls, mpi_integer, commglobal, ierror)

      if (nrtot .gt. 0) then
         ind = 0
         do i = 1, numpro
            do nn = 1, nrblcks(i)
               ind = ind + 1
               nrblckdata(i)%data(nn,1) = recvtot(1,ind)
               nrblckdata(i)%data(nn,2) = recvtot(2,ind)
            enddo
         enddo
      endif

! compute derived types for receive blocks
      do i = 1, numpro
         if (nrblckdata(i)%size .gt. 0) then
            call mpi_type_indexed( nrblckdata(i)%size, nrblckdata(i)%data(1,1),   &
                                   nrblckdata(i)%data(1,2), mpi_double_precision, &
                                   send_bl(i)%partnertype, ierror)
            call mpi_type_commit( send_bl(i)%partnertype, ierror)
         else
            send_bl(i)%partnertype = MPI_DATATYPE_NULL
         endif
         recv_bl(i)%partnertype = MPI_DATATYPE_NULL
      enddo

      deallocate (sendtot)
      deallocate (recvtot)
      do i = 1, numpro
         deallocate (nrblckdata(i)%data)
      enddo

     endif  ! mod_method .ge. 3

    endif  ! mpi2_flag .eq. 1

      end subroutine get_partneroffset
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr --- Write r8 contiguous chunks to global array
!
! !INTERFACE:
      subroutine mp_sendirr ( q, send_bl, recv_bl, qout )
 
! !INPUT PARAMETERS:
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      real(r8), intent(in) :: q(*)                     ! local data segment

! !OUTPUT PARAMETERS:
      real(r8), intent(out) :: qout(*)                 ! local output segment
!
! !DESCRIPTION:
!     Send a number of contiguous chunks destined to an arbitrary set 
!     of PEs.  This is basically the non-blocking start of an
!     all-to-all communcation primitive.  This fundamental
!     routine forms a basis for higher level primitives.
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!
! !BUGS:
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nchunks, offset_s, offset_r, ierr, mod_method, mpi2_flag
      integer p, maxset, mysize, nthpc, minsize, nthrd, pn, pt, tmpsize
#if defined (USE_MPI2)
      integer (kind=MPI_ADDRESS_KIND) target_disp, winsize8
      integer (kind=MPI_ADDRESS_KIND) swin
#endif
      integer, allocatable :: addset (:)
      integer nwin, nwinhit, ns, i, j, qsiz, send_tag
      logical locwin

!
!     initialize window

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif

      mod_method = recv_bl(1)%method

! locwin is flag for local window
      locwin = mod_method .gt. 2 .and. mpi2_flag .eq. 1
      if (.not. locwin) call Win_Open(r8_win)
    if (mod_method .gt. 0) then

#if defined (USE_MPI2)
     if (mod_method .eq. 1) then
!  directly put contiguous segments into global target window
      maxset = 0
      do ipe = 1, numpro
         nchunks = size( send_bl(ipe)%displacements )
         maxset = max (maxset, nchunks)
      enddo
      allocate (addset(maxset))
      do ipe = 1, numpro
         qsize = SUM( send_bl(ipe)%blocksizes )
         if (qsize .ne. 0) then
            nchunks = size( send_bl(ipe)%displacements )
            minsize = send_bl(ipe)%blocksizes(1)
            addset(1) = 0
            do p = 2, nchunks
                minsize = min(minsize, send_bl(ipe)%blocksizes(p))
                addset(p) = addset(p-1) + send_bl(ipe)%blocksizes(p-1)
            enddo
            r8_win%dest = ipe-1
            nthpc = ceiling(real(pkgs_per_pro)/real(nchunks))
            nthrd = min(nthpc,minsize)      !  number of threads per block
            r8_win%offset_r = send_bl(ipe)%partneroffset
!$omp parallel do private(pn,p,pt,tmpsize,mysize,offset_s,ierr)
            do pn = 1, nchunks*nthrd
               p = (pn-1)/nthrd + 1        ! block number
               pt = mod(pn-1,nthrd) + 1    ! thread number
               tmpsize = ceiling(real(send_bl(ipe)%blocksizes(p))/real(nthrd))
               mysize = min(tmpsize, max(send_bl(ipe)%blocksizes(p)-tmpsize*(pt-1),0))
               offset_s = send_bl(ipe)%displacements(p)
               call MPI_PUT(q(offset_s+(pt-1)*tmpsize+1), mysize,       &
                            mp_r8, r8_win%dest,                         &
                            r8_win%offset_r+addset(p)+(pt-1)*tmpsize,   &
                            mysize, mp_r8, r8_win%id, ierr)
            enddo
         endif
      enddo
      deallocate (addset)
     elseif (mod_method .eq. 2) then
! directly put derived types into global target window
!$omp parallel do private(ipe,mysize,ierr,target_disp)
        do ipe = 1, numpro
           if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
               mysize = SUM( send_bl(ipe)%blocksizes )
               target_disp = send_bl(ipe)%partneroffset
               call MPI_PUT(q, 1, send_bl(ipe)%type, ipe-1, target_disp,   &
                            mysize, mp_r8, r8_win%id, ierr)
           endif
        enddo
     elseif (mod_method .ge. 3) then
! directly put derived types into derived types of local target window
        if (nwininit .eq. 0) then
           nwininit = 1
           wintable(:,:) = 0
           nwinused = 0
        endif
        winaddr = loc(qout)
        winsize = 1
        do ipe = 1, numpro
           ns = size(recv_bl(ipe)%blocksizes)
           if (ns .gt. 0) then
              swin = recv_bl(ipe)%blocksizes(ns) + recv_bl(ipe)%displacements(ns)
              winsize = max(winsize, swin)
           endif
        enddo
! check to see if window is already active
        nwinhit = 0
        do nwin = 1, nwinused
           if (winaddr .eq. wintable(nwin,1) .and. winsize .eq. wintable(nwin,2))   &
                  nwinhit = nwin
        enddo

! See comments at beginning of module regarding validity of this approach
        nwinhit = 0  ! force window recreation to avoid partial hit deadlock
! See comments at beginning of module regarding validity of this approach

        if (nwinhit .eq. 0) then
           nwinused = nwinused + 1
           if (nwinused .gt. maxwin) then
              do nwin = 1, maxwin
                 call MPI_WIN_FREE(winlocal(nwin), ierr)
              enddo
              wintable(:,:) = 0
              nwinused = 1
           endif
           nwinhit = nwinused
           wintable(nwinhit,1) = winaddr
           wintable(nwinhit,2) = winsize
           winsize8 = 8*winsize
           call MPI_WIN_CREATE(qout, winsize8, 8, MPI_INFO_NULL,         &
                               commglobal, winlocal(nwinhit), ierr)
        endif
        target_disp = 0
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, winlocal(nwinhit), ierr)
!$omp parallel do private(ipe, ierr)
        do ipe = 1, numpro
           if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
              call MPI_PUT(q, 1, send_bl(ipe)%type, ipe-1, target_disp, 1,   &
                           send_bl(ipe)%partnertype, winlocal(nwinhit),   &
                           ierr)
           endif
        enddo
     endif
#else
!
! mpi-1 with derived types
! Increment the ongoing transfer number
      BegTrf = MOD(BegTrf,MAX_TRF) + 1
!
! MPI: Irecv over all processes
!
      do ipe=1, numpro

!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        OutHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_irecv( qout, 1, recv_bl(ipe)%type, ipe-1, ipe-1,     &
                          commglobal, OutHandle(ipe,BegTrf), ierr )
        endif
      enddo

!
! MPI: Isend over all processes
!
      do ipe=1, numpro

!
! Send the individual buffers with non-blocking sends
!
        InHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_isend( q, 1, send_bl(ipe)%type, ipe-1, gid,        &
                          commglobal, InHandle(ipe,BegTrf), ierr )
        endif
      enddo
#endif
    else

! temporary contiguous buffers
      blocksize = r8_win%size / numpro   ! size alloted to this PE
      offset = gid*blocksize
      offset_s = 0
      offset_r = 0

      do ipe=1, numpro
#if !defined(USE_MPI2)
         r8_win%size_r = SUM(recv_bl(ipe)%blocksizes)
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            if (offset_r .gt. r8_win%size) then
              print *, "Fatal mp_sendirr: receive window out of space - exiting"
              call exit(1)
            endif
            r8_win%src = ipe-1
            call Ga_RecvInit_r8(r8_win, ga_r8_r)
         endif
#endif
         qsize = SUM( send_bl(ipe)%blocksizes )
         if (qsize .ne. 0) then
            r8_win%dest = ipe-1
            r8_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s .gt. r8_win%size) then
              print *, "Fatal mp_sendirr: send window out of space - exiting"
              call exit(1)
            endif
            nchunks = size( send_bl(ipe)%displacements )
            qsiz = 0
            do j = 1, nchunks
               do i = send_bl(ipe)%displacements(j)+1, send_bl(ipe)%displacements(j)+send_bl(ipe)%blocksizes(j)
                  qsiz = qsiz+1
                  ga_r8_s(r8_win%offset_s+qsiz) = q(i)
               enddo
            enddo

#if defined(USE_MPI2)
            r8_win%offset_r = send_bl(ipe)%partneroffset
            if ( r8_win%offset_r+qsiz > r8_win%size ) then
               print *, "Fatal mp_sendirr: target window out of space - exiting"
               print *, 'gid r8_win%offset_r qsiz r8_win%size = ', gid,      &
                         r8_win%offset_r, qsiz, r8_win%size
               call exit(1)
            endif
            tmpsize = ceiling(real(qsiz)/real(pkgs_per_pro))
!$omp parallel do private(p,mysize,ierr)
            do p=1,MIN(pkgs_per_pro,qsiz)
              mysize = MIN(tmpsize, MAX(qsiz-(tmpsize*(p-1)),0))
              call MPI_PUT(ga_r8_s(r8_win%offset_s+(p-1)*tmpsize+1), mysize, mp_r8,     &
                           r8_win%dest, r8_win%offset_r+(p-1)*tmpsize, mysize, mp_r8, &
                           r8_win%id, ierr)
            enddo
#else
            send_tag = gid
            nsend = nsend + 1
            call MPI_ISEND(ga_r8_s(r8_win%offset_s+1), qsiz, mp_r8, r8_win%dest, &
                           send_tag, commglobal, sqest(nsend), ierr)
#endif
         endif
         offset = offset + qsize
      enddo

    endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

      end subroutine mp_sendirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr --- Read r8 contiguous chunks to global array
!
! !INTERFACE:
      subroutine mp_recvirr ( qout, recv_bl )
 
! !INPUT PARAMETERS:
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! Global Array Window
! !INPUT/OUTPUT PARAMETERS:
      real(r8), intent(inout) :: qout(*)               ! local data segment
!
! !DESCRIPTION:
!
!     Complete transfer of a generalized region initiated by {\tt mp\_sendirr}.
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, nchunks, offset_r, mod_method, mpi2_flag
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
      integer nwin, nwinhit, ns
#if defined (USE_MPI2)
      integer (kind=MPI_ADDRESS_KIND) swin
#endif
      integer qsiz, i, j
      logical locwin


#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif
      mod_method = recv_bl(1)%method

! locwin is flag for local window
    locwin = mod_method .gt. 2 .and. mpi2_flag .eq. 1
    if (.not. locwin) call Win_Close(r8_win)
    if (mod_method .gt. 0 .and. mpi2_flag .eq. 0) then

! mpi-1 derived types
      EndTrf = MOD(EndTrf,MAX_TRF) + 1
      CALL MPI_WAITALL( numpro, InHandle(:,EndTrf), InStats, Ierr )
      CALL MPI_WAITALL( numpro, OutHandle(:,EndTrf), OutStats, Ierr )

    elseif (mod_method .lt. 3 .or. mpi2_flag .eq. 0) then

! temporary contiguous buffer / global window
      blocksize = r8_win%size / numpro   ! size alloted to this PE
      offset_r = 0
      do ipe=1, numpro
         r8_win%size_r = SUM(recv_bl(ipe)%blocksizes)
         if (r8_win%size_r .ne. 0) then
            r8_win%offset_r = offset_r
            offset_r = offset_r + r8_win%size_r
            if (offset_r .gt. r8_win%size) then
              print *, "Fatal mp_recvirr: receive window out of space - exiting"
              call exit(1)
            endif
            nchunks = SIZE( recv_bl(ipe)%displacements )

#if !defined(USE_MPI2)
            nread = nread + 1
            call MPI_WAIT(rqest(nread), Status, ierr)
#endif
            qsiz = 0
            do j = 1, nchunks
               do i = recv_bl(ipe)%displacements(j)+1, recv_bl(ipe)%displacements(j)+recv_bl(ipe)%blocksizes(j)
                  qsiz = qsiz+1
                  qout(i) = ga_r8_r(r8_win%offset_r+qsiz)
               enddo
            enddo
         endif
      enddo

    else

! local window
#if defined (USE_MPI2)
       if (nwininit .eq. 0) then
          print *, "Fatal mp_recvirr: window initialization should have already occurred - exiting"
          call exit(1)
       endif
       winaddr = loc(qout)
       winsize = 1
       do ipe = 1, numpro
          ns = size(recv_bl(ipe)%blocksizes)
          if (ns .gt. 0) then
             swin = recv_bl(ipe)%blocksizes(ns) + recv_bl(ipe)%displacements(ns)
             winsize = max(winsize, swin)
          endif
       enddo
       nwinhit = 0
       do nwin = 1, nwinused
          if (winaddr .eq. wintable(nwin,1) .and. winsize .eq. wintable(nwin,2))   &
                 nwinhit = nwin
       enddo
! See comments at beginning of module regarding validity of this approach
       if (nwinhit .eq. 0) then
          print *, "Fatal mp_recvirr: window should have been initialized - exiting"
          call exit(1)
       endif

       call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED,         &
                          winlocal(nwinhit), ierr)
#endif

    endif
    if (.not. locwin) call Win_Finalize(r8_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recvirr
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_sendirr_i4 --- Write i4 contiguous chunks to global array
!
! !INTERFACE:
      subroutine mp_sendirr_i4 ( q, send_bl, recv_bl, qout )
 
! !INPUT PARAMETERS:
      type(blockdescriptor), intent(in)  :: send_bl(:) ! send blocks
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! receive blocks
      integer(i4), intent(in) :: q(*)                     ! local data segment

! !OUTPUT PARAMETERS:
      integer(i4), intent(out) :: qout(*)                 ! local output segment
!
! !DESCRIPTION:  Integer I4 version
!     Send a number of contiguous chunks destined to an arbitrary set 
!     of PEs.  This is basically the non-blocking start of an
!     all-to-all communcation primitive.  This fundamental
!     routine forms a basis for higher level primitives.  
!
! !REVISION HISTORY: 
!    02.08.13   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Use partneroffset
!    03.06.24   Sawyer      Integrated Use_Mpi_Types; added qout
!    04.02.24   Mirin       Various mpi2 options
!
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
      integer ipe, qsize, offset, blocksize, nchunks, offset_s, offset_r, ierr, mod_method, mpi2_flag
      integer p, maxset, mysize, nthpc, minsize, nthrd, pn, pt, tmpsize
#if defined (USE_MPI2)
      integer (kind=MPI_ADDRESS_KIND) target_disp, winsize4
      integer (kind=MPI_ADDRESS_KIND) swin
#endif
      integer, allocatable :: addset (:)
      integer nwin, nwinhit, ns, i, j, qsiz, send_tag
      logical locwin

!
!     initialize window

#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif

      mod_method = recv_bl(1)%method

! locwin is flag for local window
      locwin = mod_method .gt. 2 .and. mpi2_flag .eq. 1
      if (.not. locwin) call Win_Open(i4_win)
    if (mod_method .gt. 0) then

#if defined (USE_MPI2)
     if (mod_method .eq. 1) then
!  directly put contiguous segments into global target window
      maxset = 0
      do ipe = 1, numpro
         nchunks = size( send_bl(ipe)%displacements )
         maxset = max (maxset, nchunks)
      enddo
      allocate (addset(maxset))
      do ipe = 1, numpro
         qsize = SUM( send_bl(ipe)%blocksizes )
         if (qsize .ne. 0) then
            nchunks = size( send_bl(ipe)%displacements )
            minsize = send_bl(ipe)%blocksizes(1)
            addset(1) = 0
            do p = 2, nchunks
                minsize = min(minsize, send_bl(ipe)%blocksizes(p))
                addset(p) = addset(p-1) + send_bl(ipe)%blocksizes(p-1)
            enddo
            i4_win%dest = ipe-1
            nthpc = ceiling(real(pkgs_per_pro)/real(nchunks))
            nthrd = min(nthpc,minsize)      !  number of threads per block
            i4_win%offset_r = send_bl(ipe)%partneroffset
!$omp parallel do private(pn,p,pt,tmpsize,mysize,offset_s,ierr)
            do pn = 1, nchunks*nthrd
               p = (pn-1)/nthrd + 1        ! block number
               pt = mod(pn-1,nthrd) + 1    ! thread number
               tmpsize = ceiling(real(send_bl(ipe)%blocksizes(p))/real(nthrd))
               mysize = min(tmpsize, max(send_bl(ipe)%blocksizes(p)-tmpsize*(pt-1),0))
               offset_s = send_bl(ipe)%displacements(p)
               call MPI_PUT(q(offset_s+(pt-1)*tmpsize+1), mysize,       &
                            mp_i4, i4_win%dest,                         &
                            i4_win%offset_r+addset(p)+(pt-1)*tmpsize,   &
                            mysize, mp_i4, i4_win%id, ierr)
            enddo
         endif
      enddo
      deallocate (addset)
     elseif (mod_method .eq. 2) then
! directly put derived types into global target window
!$omp parallel do private(ipe,mysize,ierr,target_disp)
        do ipe = 1, numpro
           if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
               mysize = SUM( send_bl(ipe)%blocksizes )
               target_disp = send_bl(ipe)%partneroffset
               call MPI_PUT(q, 1, send_bl(ipe)%type, ipe-1, target_disp,   &
                            mysize, mp_i4, i4_win%id, ierr)
           endif
        enddo
     elseif (mod_method .ge. 3) then
! directly put derived types into derived types of local target window
        if (nwininit .eq. 0) then
           nwininit = 1
           wintable(:,:) = 0
           nwinused = 0
        endif
        winaddr = loc(qout)
        winsize = 1
        do ipe = 1, numpro
           ns = size(recv_bl(ipe)%blocksizes)
           if (ns .gt. 0) then
              swin = recv_bl(ipe)%blocksizes(ns) + recv_bl(ipe)%displacements(ns)
              winsize = max(winsize, swin)
           endif
        enddo
! check to see if window is already active
        nwinhit = 0
        do nwin = 1, nwinused
           if (winaddr .eq. wintable(nwin,1) .and. winsize .eq. wintable(nwin,2))   &
                  nwinhit = nwin
        enddo

! See comments at beginning of module regarding validity of this approach
        nwinhit = 0  ! force window recreation to avoid partial hit deadlock
! See comments at beginning of module regarding validity of this approach

        if (nwinhit .eq. 0) then
           nwinused = nwinused + 1
           if (nwinused .gt. maxwin) then
              do nwin = 1, maxwin
                 call MPI_WIN_FREE(winlocal(nwin), ierr)
              enddo
              wintable(:,:) = 0
              nwinused = 1
           endif
           nwinhit = nwinused
           wintable(nwinhit,1) = winaddr
           wintable(nwinhit,2) = winsize
           winsize4 = 4*winsize
           call MPI_WIN_CREATE(qout, winsize4, 4, MPI_INFO_NULL,         &
                               commglobal, winlocal(nwinhit), ierr)
        endif
        target_disp = 0
        call MPI_WIN_FENCE(MPI_MODE_NOPRECEDE, winlocal(nwinhit), ierr)
!$omp parallel do private(ipe, ierr)
        do ipe = 1, numpro
           if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
              call MPI_PUT(q, 1, send_bl(ipe)%type, ipe-1, target_disp, 1,   &
                           send_bl(ipe)%partnertype, winlocal(nwinhit),   &
                           ierr)
           endif
        enddo
     endif
#else
!
! mpi-1 with derived types
! Increment the ongoing transfer number
      BegTrf = MOD(BegTrf,MAX_TRF) + 1
!
! MPI: Irecv over all processes
!
      do ipe=1, numpro

!
! Receive the buffers with MPI_Irecv. Non-blocking
!
        OutHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( recv_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_irecv( qout, 1, recv_bl(ipe)%type, ipe-1, ipe-1,     &
                          commglobal, OutHandle(ipe,BegTrf), ierr )
        endif
      enddo

!
! MPI: Isend over all processes
!
      do ipe=1, numpro

!
! Send the individual buffers with non-blocking sends
!
        InHandle(ipe,BegTrf) = MPI_REQUEST_NULL
        if ( send_bl(ipe)%type /= MPI_DATATYPE_NULL ) then
          call mpi_isend( q, 1, send_bl(ipe)%type, ipe-1, gid,        &
                          commglobal, InHandle(ipe,BegTrf), ierr )
        endif
      enddo
#endif
    else

! temporary contiguous buffers
      blocksize = i4_win%size / numpro   ! size alloted to this PE
      offset = gid*blocksize
      offset_s = 0
      offset_r = 0

      do ipe=1, numpro
#if !defined(USE_MPI2)
         i4_win%size_r = SUM(recv_bl(ipe)%blocksizes)
         if (i4_win%size_r .ne. 0) then
            i4_win%offset_r = offset_r
            offset_r = offset_r + i4_win%size_r
            if (offset_r .gt. i4_win%size) then
              print *, "Fatal mp_sendirr_i4: receive window out of space - exiting"
              call exit(1)
            endif
            i4_win%src = ipe-1
            call Ga_RecvInit_i4(i4_win, ga_i4_r)
         endif
#endif
         qsize = SUM( send_bl(ipe)%blocksizes )
         if (qsize .ne. 0) then
            i4_win%dest = ipe-1
            i4_win%offset_s = offset_s
            offset_s = offset_s + qsize
            if (offset_s .gt. i4_win%size) then
              print *, "Fatal mp_sendirr_i4: send window out of space - exiting"
              call exit(1)
            endif
            nchunks = size( send_bl(ipe)%displacements )
            qsiz = 0
            do j = 1, nchunks
               do i = send_bl(ipe)%displacements(j)+1, send_bl(ipe)%displacements(j)+send_bl(ipe)%blocksizes(j)
                  qsiz = qsiz+1
                  ga_i4_s(i4_win%offset_s+qsiz) = q(i)
               enddo
            enddo

#if defined(USE_MPI2)
            i4_win%offset_r = send_bl(ipe)%partneroffset
            if ( i4_win%offset_r+qsiz > i4_win%size ) then
               print *, "Fatal mp_sendirr_i4: target window out of space - exiting"
               print *, 'gid i4_win%offset_r qsiz i4_win%size = ', gid,      &
                         i4_win%offset_r, qsiz, i4_win%size
               call exit(1)
            endif
            tmpsize = ceiling(real(qsiz)/real(pkgs_per_pro))
!$omp parallel do private(p,mysize,ierr)
            do p=1,MIN(pkgs_per_pro,qsiz)
              mysize = MIN(tmpsize, MAX(qsiz-(tmpsize*(p-1)),0))
              call MPI_PUT(ga_i4_s(i4_win%offset_s+(p-1)*tmpsize+1), mysize, mp_i4,     &
                           i4_win%dest, i4_win%offset_r+(p-1)*tmpsize, mysize, mp_i4, &
                           i4_win%id, ierr)
            enddo
#else
            send_tag = gid
            nsend = nsend + 1
            call MPI_ISEND(ga_i4_s(i4_win%offset_s+1), qsiz, mp_i4, i4_win%dest, &
                           send_tag, commglobal, sqest(nsend), ierr)
#endif
         endif
         offset = offset + qsize
      enddo

    endif

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

      end subroutine mp_sendirr_i4
!------------------------------------------------------------------------------
!
!------------------------------------------------------------------------------
!BOP
! !ROUTINE: mp_recvirr_i4 --- Read i4 contiguous chunks to global array
!
! !INTERFACE:
      subroutine mp_recvirr_i4 ( qout, recv_bl )
 
! !INPUT PARAMETERS:
      type(blockdescriptor), intent(in)  :: recv_bl(:) ! Global Array Window
! !INPUT/OUTPUT PARAMETERS:
      integer(i4), intent(inout) :: qout(*)               ! local data segment
!
! !DESCRIPTION:
!
!     Complete transfer of a generalized region initiated by {\tt
!     mp\_sendirr\_i4}
!
! !REVISION HISTORY:
!    02.08.15   Sawyer      Creation
!    02.11.06   Mirin       Optimizations
!    03.03.03   Sawyer      Now using packed arrays for MPI2
!    04.02.24   Mirin       Various mpi2 options
!
!EOP
!------------------------------------------------------------------------------
!BOC
      integer :: ipe, blocksize, nchunks, offset_r, mod_method, mpi2_flag
      integer Ierr
      integer InStats(numpro*MPI_STATUS_SIZE)
      integer OutStats(numpro*MPI_STATUS_SIZE)
      integer nwin, nwinhit, ns
#if defined (USE_MPI2)
      integer (kind=MPI_ADDRESS_KIND) swin
#endif
      integer qsiz, i, j
      logical locwin


#if defined( MODCM_TIMING )
      call t_startf('mod_comm communication')
#endif

      mpi2_flag = 0
#if defined (USE_MPI2)
      mpi2_flag = 1
#endif
      mod_method = recv_bl(1)%method

! locwin is flag for local window
    locwin = mod_method .gt. 2 .and. mpi2_flag .eq. 1
    if (.not. locwin) call Win_Close(i4_win)
    if (mod_method .gt. 0 .and. mpi2_flag .eq. 0) then

! mpi-1 derived types
      EndTrf = MOD(EndTrf,MAX_TRF) + 1
      CALL MPI_WAITALL( numpro, InHandle(:,EndTrf), InStats, Ierr )
      CALL MPI_WAITALL( numpro, OutHandle(:,EndTrf), OutStats, Ierr )

    elseif (mod_method .lt. 3 .or. mpi2_flag .eq. 0) then

! temporary contiguous buffer / global window
      blocksize = i4_win%size / numpro   ! size alloted to this PE
      offset_r = 0
      do ipe=1, numpro
         i4_win%size_r = SUM(recv_bl(ipe)%blocksizes)
         if (i4_win%size_r .ne. 0) then
            i4_win%offset_r = offset_r
            offset_r = offset_r + i4_win%size_r
            if (offset_r .gt. i4_win%size) then
              print *, "Fatal mp_recvirr_i4: receive window out of space - exiting"
              call exit(1)
            endif
            nchunks = SIZE( recv_bl(ipe)%displacements )

#if !defined(USE_MPI2)
            nread = nread + 1
            call MPI_WAIT(rqest(nread), Status, ierr)
#endif
            qsiz = 0
            do j = 1, nchunks
               do i = recv_bl(ipe)%displacements(j)+1, recv_bl(ipe)%displacements(j)+recv_bl(ipe)%blocksizes(j)
                  qsiz = qsiz+1
                  qout(i) = ga_i4_r(i4_win%offset_r+qsiz)
               enddo
            enddo
         endif
      enddo

    else

! local window
#if defined (USE_MPI2)
       if (nwininit .eq. 0) then
          print *, "Fatal mp_recvirr_i4: window initialization should have already occurred - exiting"
          call exit(1)
       endif
       winaddr = loc(qout)
       winsize = 1
       do ipe = 1, numpro
          ns = size(recv_bl(ipe)%blocksizes)
          if (ns .gt. 0) then
             swin = recv_bl(ipe)%blocksizes(ns) + recv_bl(ipe)%displacements(ns)
             winsize = max(winsize, swin)
          endif
       enddo
       nwinhit = 0
       do nwin = 1, nwinused
          if (winaddr .eq. wintable(nwin,1) .and. winsize .eq. wintable(nwin,2))   &
                 nwinhit = nwin
       enddo
! See comments at beginning of module regarding validity of this approach
       if (nwinhit .eq. 0) then
          print *, "Fatal mp_recvirr_i4: window should have been initialized - exiting"
          call exit(1)
       endif

       call MPI_WIN_FENCE(MPI_MODE_NOSTORE + MPI_MODE_NOSUCCEED,         &
                          winlocal(nwinhit), ierr)
#endif

    endif
    if (.not. locwin) call Win_Finalize(i4_win)

#if defined( MODCM_TIMING )
      call t_stopf('mod_comm communication')
#endif

!EOC
      end subroutine mp_recvirr_i4
!------------------------------------------------------------------------------

#endif
      end module mod_comm

