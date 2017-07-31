#include <misc.h>
#include <params.h>

module spmd_dyn
!BOP
!
! !MODULE: Subroutines to initialize SPMD implementation of CAM
!

#if (defined SPMD)

!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plat, plon, masterproc, iam, numbnd, &
                     numlats, beglat, endlat, &
                     plev, beglev, endlev, endlevp1, &
                     endlevp, myid_y, myid_z, npr_y, npr_z, plevp, &
                     myidxy_x, myidxy_y, nprxy_x, nprxy_y, &
                     beglonxy, endlonxy, beglatxy, endlatxy, &
                     twod_decomp, spmd_on, mod_transpose, mod_geopk
   use constituents, only: ppcnst
   use mpishorthand, only: mpir8, mpicom, mpiint
   use decompmodule, only: decomptype, decompcreate
   use ghostmodule, only:  ghosttype
   use parutilitiesmodule, only: parpatterntype
   use mod_comm, only: tracertrans_mod
   use infnan, only: inf
   use abortutils, only: endrun

   implicit none

! !PUBLIC MEMBER FUNCTIONS:

   public spmdinit_dyn, decomp_wavenumbers
   public compute_gsfactors, spmdbuf

! !PUBLIC DATA MEMBERS:

   integer :: npes                    ! Total number of MPI tasks
   integer :: nsmps                   ! Total number of SMP nodes
   integer :: proc(plat)              ! processor id associated with a given lat.
   integer, allocatable :: cut(:,:)   ! partition for MPI tasks
   integer, allocatable :: nlat_p(:)  ! number of latitudes per subdomain
   integer, allocatable :: proc_smp_map(:) ! map of process/SMP node assignments
   integer, allocatable :: kextent(:) ! number of levels per subdomain

   integer comm_y            ! communicator in latitude
   integer comm_z            ! communicator in vertical
   integer commxy_x          ! communicator in longitude (xy second. decomp.)
   integer commxy_y          ! communicator in latitude (xy second. decomp.)
   integer, allocatable :: lonrangexy(:,:)   ! global xy-longitude subdomain index
   integer, allocatable :: latrangexy(:,:)   ! global xy-latitude subdomain index
   logical geopk16byte ! method for geopotential calculation with 2D decomp.

   integer m_ttrans, q_ttrans, r_ttrans

   type (ghosttype), save  :: ghostpe_yz, ghostpe1_yz
   type (parpatterntype)   :: ikj_xy_to_yz, ijk_yz_to_xy, &
                              ijk_xy_to_yz, q_to_qxy, qxy_to_q,         &
                              pexy_to_pe, pkxy_to_pkc, r_to_rxy, rxy_to_r
!
! !DESCRIPTION: 
!   {\bf Purpose:} Subroutines to initialize SPMD implementation of CAM
!
! !REVISION HISTORY:
!   ??.??.??  CCM Core Group     Creation
!   00.09.30  Sawyer             Alterations for LR SPMD mode
!   01.05.09  Mirin              2-D yz decomposition
!   01.06.27  Mirin              Secondary 2-D xy decomposition
!   01.12.20  Sawyer             Changed index order of Q3 decomposition
!   02.12.11  Sawyer             Use parbegin/endtransfer for transposes
!   03.05.07  Sawyer             Removed unneeded decompositions
!
!EOP
!-----------------------------------------------------------------------

contains

  subroutine spmd_dyn_defaultopts(npr_yz_out, geopktrans_out,    &
               tracertrans_out, ompnest_out, force_2d_out,       &
               modcomm_transpose_out, modcomm_geopk_out,         &
               dyn_alltoall_out, dyn_allgather_out  )
!----------------------------------------------------------------------
! Purpose: Return default runtime options
! Author: Art Mirin
!----------------------------------------------------------------------
!------------------------------Arguments-------------------------------
     ! yz and xy decompositions
     integer, intent(out), optional :: npr_yz_out(4)
     ! geopotential method (routine geopk)
     integer, intent(out), optional :: geopktrans_out
     ! number of simultaneously transposed tracers
     integer, intent(out), optional :: tracertrans_out
     ! option for nested openmp
     integer, intent(out), optional :: ompnest_out
     ! option to force transpose computation for 1D decomp.
     integer, intent(out), optional :: force_2d_out
     ! mod_comm transpose method (varies with mpi flag)
     integer, intent(out), optional :: modcomm_transpose_out
     ! mod_comm geopk method (varies with mpi flag)
     integer, intent(out), optional :: modcomm_geopk_out
! EUL/SLD-only arguments
     integer, intent(out), optional :: dyn_alltoall_out
     integer, intent(out), optional :: dyn_allgather_out
!----------------------------------------------------------------------
     if (present(npr_yz_out) ) then
        npr_yz_out(1) = npes
        npr_yz_out(2) = 1
        npr_yz_out(3) = 1
        npr_yz_out(4) = npes
     endif
     if (present(geopktrans_out) ) then
        geopktrans_out = 0
     endif
     if (present(tracertrans_out) ) then
        tracertrans_out = 6
     endif
     if (present(ompnest_out) ) then
        ompnest_out = 1
     endif
     if (present(force_2d_out) ) then
        force_2d_out = 0
     endif
     if (present(modcomm_transpose_out) ) then
        modcomm_transpose_out = 1
     endif
     if (present(modcomm_geopk_out) ) then
        modcomm_geopk_out = 1
     endif

     return

  end subroutine spmd_dyn_defaultopts

  subroutine spmd_dyn_setopts(npr_yz_in, geopktrans_in,       &
               tracertrans_in, ompnest_in, force_2d_in,       &
               modcomm_transpose_in, modcomm_geopk_in,        &
               dyn_alltoall_in, dyn_allgather_in  )
  use pmgrid, only: omptotal, ompouter, ompinner
!----------------------------------------------------------------------
! Purpose: Set runtime options
! Author: Art Mirin
!----------------------------------------------------------------------
!------------------------------Arguments-------------------------------
! yz and xy decompositions (npr_y, npr_z, nprxy_x, nprxy_y)
! must satisfy npes = npr_y*npr_z = nprxy_x*nprxy_y
! beneficial to satisfy npr_y=nprxy_y and npr_z=nprxy_x
     integer, intent(in), optional :: npr_yz_in(4)
! geopotential method (routine geopk)
! 0 for transpose method, 1 for method using semi-global z communication with optional
!   16-byte arithmetic; the transpose and semi-global z communication (with 16-byte
!   arithmetic) methods are both bit-for-bit across decompositions; the former scales
!   better with npr_z, and the latter is superior for small npr_z; optimum speed is
!   attained using the semi-global z communication method with 8-byte arithmetic (standard
!   for geopk16); see geopk16 within geopk.F90.
     integer, intent(in), optional :: geopktrans_in
! number of simultaneously transposed tracers
! this is a tradeoff between latency and communication/computation overlap
     integer, intent(in), optional :: tracertrans_in
! option for nested openmp (ibm only)
     integer, intent(in), optional :: ompnest_in
! option to force transpose computation for 1D decomp.
! the only purpose for invoking this option is debugging
     integer, intent(in), optional :: force_2d_in
! mod_comm transpose/geopk method (varies with mpi flag)
! with standard mpi:
!   0 for temporary contiguous buffers
!   1 for mpi derived types (default)
! with use_mpi2:
!   0 for temporary contiguous buffers
!   1 for direct mpi_put's of contiguous segments into temporary contiguous window, with
!     threading over the segments (default)
!   2 for mpi derived types into temporary contiguous window, with threading over the target
!   3 for derived types at source and target, with threading over the target
     integer, intent(in), optional :: modcomm_transpose_in, modcomm_geopk_in
! EUL/SLD-only arguments
     integer, intent(in), optional :: dyn_alltoall_in
     integer, intent(in), optional :: dyn_allgather_in
!----------------------------------------------------------------------
     integer omp_get_num_threads
!----------------------------------------------------------------------
     if (present(npr_yz_in) ) then
        npr_y   = npr_yz_in(1)
        npr_z   = npr_yz_in(2)
        nprxy_x = npr_yz_in(3)
        nprxy_y = npr_yz_in(4)
        if (masterproc) then
           write (6,*) 'npr_y = ', npr_y, '  npr_z = ', npr_z
           write (6,*) 'nprxy_x= ', nprxy_x, '  nprxy_y = ', nprxy_y
        endif
        if (npr_y*npr_z /= npes .or. nprxy_x*nprxy_y /= npes) then
           call endrun ('SPMD_DYN_SET : incorrect domain decomposition - aborting')
        endif
     else
        npr_y   = npes
        npr_z   = 1
        nprxy_x = 1
        nprxy_y = npes
        if (masterproc) then
           write (6,*) 'WARNING : npr_yz not present - using 1-D domain decomposition'
        endif
     endif
     myid_z   = iam/npr_y
     myid_y   = iam - myid_z*npr_y
     myidxy_y = iam/nprxy_x
     myidxy_x = iam - myidxy_y*nprxy_x

     geopk16byte = .false.
     if (present(geopktrans_in) ) then
        if (geopktrans_in .ne. 0) geopk16byte = .true.
        if (masterproc) then
           write (6,*) 'non-transpose geopk communication method = ', geopk16byte
        endif
     else
        if (masterproc) then
           write (6,*) 'WARNING : geopktrans not present - using transpose method'
        endif
     endif

     m_ttrans = 1
     if (present(tracertrans_in) ) then
        m_ttrans = tracertrans_in
        if (masterproc) then
           write (6,*) 'number of tracers simultaneously transposed = ', m_ttrans
        endif
     else
        if (masterproc) then
           write (6,*) 'WARNING : tracertrans not present - transposing one at a time'
        endif
     endif
     tracertrans_mod = m_ttrans

#if defined( _OPENMP )
!$omp parallel
     omptotal = omp_get_num_threads()
!$omp end parallel
#else
     omptotal = 1
#endif

     if (present(ompnest_in) ) then
        if (mod(omptotal,ompnest_in) .ne. 0) then
           call endrun ('SPMD_DYN_SET : bad value of ompnest - aborting')
        endif
        ompinner = ompnest_in
        ompouter = omptotal/ompinner
        if (masterproc) then
           write (6,*) 'omptotal ompouter ompinner = ', omptotal, ompouter, ompinner
        endif
     else
        ompinner = 1
        ompouter = omptotal
        if (masterproc) then
           write (6,*) 'WARNING : ompnest not present - no openmp nesting'
        endif
     endif

     twod_decomp = 1

     if (present(force_2d_in) ) then
        if (npr_z .eq. 1 .and. nprxy_x .eq. 1 .and. force_2d_in .eq. 0) then
           twod_decomp = 0
           if (masterproc) then
              print *, 'decomposition is effectively 1D - skipping transposes'
           endif
        else
           if (masterproc) then
              print *, 'using multi-2d decomposition methodology'
           endif
        endif
     else
        if (npr_z .eq. 1 .and. nprxy_x .eq. 1 ) twod_decomp = 0
        if (masterproc) then
           write (6,*) 'WARNING : force_2d not present - ignoring'
        endif
     endif

     if (present(modcomm_transpose_in) ) then
        mod_transpose = modcomm_transpose_in
        if (masterproc) then
           write (6,*) 'modcomm transpose method = ', mod_transpose
        endif
     else
        mod_transpose = 1
        if (masterproc) then
           write (6,*) 'WARNING : modcomm_transpose not present - ignoring'
        endif
     endif

     if (present(modcomm_geopk_in) ) then
        mod_geopk = modcomm_geopk_in
        if (masterproc) then
           write (6,*) 'modcomm geopk method = ', mod_geopk
        endif
     else
        mod_geopk = 1
        if (masterproc) then
           write (6,*) 'WARNING : modcomm_geopk not present - ignoring'
        endif
     endif

     return

  end subroutine spmd_dyn_setopts

!-----------------------------------------------------------------------
!BOP
! !IROUTINE: spmdinit_dyn --- SPMD initialization for dynamics
!
! !INTERFACE:

   subroutine spmdinit_dyn ()

! !USES:
      use parutilitiesmodule, only : parinit, parsplit
      use decompmodule, only : decompcreate
      use pmgrid, only: strip2d, strip3dxyz, strip3dxzy,                     &
                        strip3dxyzp, strip3dxzyp, strip3zaty,                &
                        strip3zatypt, strip3yatz, strip3yatzp,               &
                        strip3kxyz, strip3kxzy, strip3kxyzp, strip3kxzyp

! !DESCRIPTION:
!
!   SPMD initialization routine: get number of cpus, processes, tids, etc
!
! !REVISION HISTORY:
!   ??.??.??  CCM Core Group     Creation
!   00.09.30  Sawyer             Added LR-specific initialization
!   01.03.26  Sawyer             Added ProTeX documentation
!   01.06.27  Mirin              Secondary 2-D xy decomposition
!   01.10.16  Sawyer             Added Y at each Z decompositions
!   03.07.22  Sawyer             Removed decomps used by highp2
!
!EOP
!-----------------------------------------------------------------------
!BOC
! !LOCAL VARIABLES:
      integer procid    ! processor id
      integer procids   ! processor id SH
      integer procidn   ! processor id NH
      integer lat       ! latitude index
      integer iend      ! ending latitude band of work for a given proc
      integer workleft  ! amount of work still to be parcelled out
      integer actual    ! actual amount of work parcelled out
      integer ideal     ! ideal amt of work to parcel out
      integer pesleft   ! number of procs still to be given work
      integer isum      ! running total of work parcelled out
      integer smostlat  ! southern-most latitude index
      integer nmostlat  ! northern-most latitude index
      integer m2,m3,m5  ! 2, 3, 5 prime factors for problem decomposition
      integer xdist(1)  ! number of lons per subdomain
      integer, allocatable :: ydist(:) ! number of lats per subdomain
      integer, allocatable :: zdist(:) ! number of levels per subdomain
      integer, allocatable :: zdistq(:) ! number of levels per subdomain for Q3
      integer ier       ! error flag
      integer rank_y, size_y   !  rank and size wrt y-communicator
      integer rank_z, size_z   !  rank and size wrt z-communicator
      integer rankxy_x, sizexy_x   !  rank and size wrt xy x-communicator
      integer rankxy_y, sizexy_y   !  rank and size wrt xy y-communicator
      integer zdist1(1) ! used for misc. decomposition definitions
      integer, allocatable :: xdistxy(:) ! number of xy-longs per subdomain
      integer, allocatable :: ydistxy(:) ! number of xy-lats per subdomain
      integer, allocatable :: ydistqxy(:) ! number of xy tracer/lats per subdomain
      integer zdistxy(1)  ! number of xy-verts per subdomain
      integer j, k, vert, lonn
      integer ydistk(1)

      spmd_on = 1

! Default 2D decomposition
      beglev = 1
      endlev = plev
      endlevp1 = plev + 1
      endlevp = plev + 1
!
! Addition for LR dynamical core to initialize PILGRIM library
!
      call parinit(mpicom)
!
! Form separate communicators
!
      call parsplit(mpicom, myid_z, iam, comm_y, rank_y, size_y)
      call parsplit(mpicom, myid_y, iam, comm_z, rank_z, size_z)
      call parsplit(mpicom, myidxy_y, iam, commxy_x, rankxy_x, sizexy_x)
      call parsplit(mpicom, myidxy_x, iam, commxy_y, rankxy_y, sizexy_y)

!
!-----------------------------------------------------------------------
!
! Compute y decomposition
!
      allocate (ydist  (npr_y))
      allocate (nlat_p (0:npes-1))
      allocate (cut    (2,0:npes-1))

      ydist(:) = 0
      nlat_p(0:npes-1) = 0

      lat = plat / npr_y
      workleft = plat - lat * npr_y
      if ( lat < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 latitudes per subdomain')
      endif
!
! Be careful:  ydist is 1-based.  NCARs arrays, e.g., cut,  are 0-based
!
      do procid=1,npr_y
         ydist(procid) = lat
      enddo

      if ( workleft /= 0 ) then
         procids = (npr_y+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = npr_y
            ydist(procids) = ydist(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               ydist(procidn) = ydist(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(ydist) /= plat ) then
         write(6,*)'SPMDINIT_DYN:', ydist,' does not add up to ', plat
         call endrun
      endif

      if (workleft/=0) then
         write(6,*)'SPMDINIT_DYN: Workleft(y) not zero.  Value is ',workleft
         call endrun
      end if

! Set the NCAR data structures

      lat  = 0
      do procid=0,npr_y-1
         cut(1,procid) = lat+1
         lat = lat + ydist(procid+1)
         cut(2,procid) = lat
         nlat_p(procid) = ydist(procid+1)

         if (masterproc) then
            write(6,*) 'nlat_p(',procid,') = ', nlat_p(procid)
         end if

         if (myid_y == procid) then
            beglat  = cut(1,myid_y)
            endlat  = cut(2,myid_y)
            numlats = ydist(procid+1)
         end if
      enddo

      do k = 1, npr_z-1
         do j = 0, npr_y-1
            procid = j + k*npr_y
            cut(1,procid) = cut(1,j)
            cut(2,procid) = cut(2,j)
            nlat_p(procid) = nlat_p(j)
         enddo
      enddo
!
! Compute z decomposition
!
      allocate (zdist (npr_z))
      allocate (zdistq(npr_z))
      allocate (kextent(npr_z))

      zdist(:) = 0

      vert = plev / npr_z
      workleft = plev - vert * npr_z
      if ( vert < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 verticals per subdomain')
      endif

      do procid=1,npr_z
         zdist(procid) = vert
      enddo

      if ( workleft /= 0 ) then
         procids = (npr_z+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = npr_z
            zdist(procids) = zdist(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               zdist(procidn) = zdist(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(zdist) /= plev ) then
         write(6,*)'SPMDINIT_DYN:', zdist,' does not add up to ', plev
         call endrun
      endif

      if (workleft/=0) then
         write(6,*)'SPMDINIT_DYN: Workleft(z) not zero.  Value is ',workleft
         call endrun
      end if

! kextent is global, zdist is local to this module
      kextent(:) = zdist(:)

! Compute local limits

      beglev = 1
      endlev = zdist(1)
      do procid = 1, myid_z
         beglev = endlev + 1
         endlev = beglev + zdist(procid+1) - 1
      enddo
      endlevp1 = endlev + 1
      endlevp = endlev
      if (myid_z == npr_z-1) endlevp = endlev + 1

!
! Compute x secondary decomposition
!
      allocate (xdistxy (nprxy_x))

      xdistxy(:) = 0

      lonn = plon / nprxy_x
      workleft = plon - lonn * nprxy_x
      if ( lonn < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 xy-longitudes per subdomain')
      endif

      do procid=1,nprxy_x
         xdistxy(procid) = lonn
      enddo

      if ( workleft /= 0 ) then
         procids = (nprxy_x+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = nprxy_x
            xdistxy(procids) = xdistxy(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               xdistxy(procidn) = xdistxy(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(xdistxy) /= plon ) then
         write(6,*)'SPMDINIT_DYN:', xdistxy,' does not add up to ', plon
         call endrun
      endif

      if (workleft/=0) then
         write(6,*)'SPMDINIT_DYN: Workleft(xy-x) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglonxy = 1
      endlonxy = xdistxy(1)
      do procid = 1, myidxy_x
         beglonxy = endlonxy + 1
         endlonxy = beglonxy + xdistxy(procid+1) - 1
      enddo

! Compute global table

      allocate (lonrangexy(2,nprxy_x))
      lonrangexy(1,1) = 1
      lonrangexy(2,1) = xdistxy(1)
      do procid = 2, nprxy_x
         lonrangexy(1,procid) = lonrangexy(2,procid-1) + 1
         lonrangexy(2,procid) = lonrangexy(1,procid) + xdistxy(procid) - 1
      enddo
!
! Compute y secondary decomposition
!
      allocate (ydistxy (nprxy_y))

      ydistxy(:) = 0

      lat = plat / nprxy_y
      workleft = plat - lat * nprxy_y
      if ( lat < 3 ) then
         call endrun ('SPMDINIT_DYN: less than 3 xy-latitudes per subdomain')
      endif

      do procid=1,nprxy_y
         ydistxy(procid) = lat
      enddo

      if ( workleft /= 0 ) then
         procids = (nprxy_y+1) / 2
         procidn = procids + 1
         do while ( workleft /= 0 )
            if ( procids == 1 ) procids = nprxy_y
            ydistxy(procids) = ydistxy(procids) + 1
            workleft = workleft - 1  
            if ( workleft /= 0 ) then
               ydistxy(procidn) = ydistxy(procidn) + 1
               workleft = workleft - 1
            endif
            procidn = procidn + 1
            procids = procids - 1
         enddo
      endif

! Safety check:

      if ( sum(ydistxy) /= plat ) then
         write(6,*)'SPMDINIT_DYN:', ydistxy,' does not add up to ', plat
         call endrun
      endif

      if (workleft/=0) then
         write(6,*)'SPMDINIT_DYN: Workleft(xy-y) not zero.  Value is ',workleft
         call endrun
      end if

! Compute local limits

      beglatxy = 1
      endlatxy = ydistxy(1)
      do procid = 1, myidxy_y
         beglatxy = endlatxy + 1
         endlatxy = beglatxy + ydistxy(procid+1) - 1
      enddo

! Compute global table

      allocate (latrangexy(2,nprxy_y))
      latrangexy(1,1) = 1
      latrangexy(2,1) = ydistxy(1)
      do procid = 2, nprxy_y
         latrangexy(1,procid) = latrangexy(2,procid-1) + 1
         latrangexy(2,procid) = latrangexy(1,procid) + ydistxy(procid) - 1
      enddo

! Define variables for grouping tracer transposes (set group size)

      q_ttrans = ppcnst/m_ttrans
      r_ttrans = ppcnst - m_ttrans*q_ttrans
!
! WS: create decompositions for NCAR data structures
!
      xdist(1) = plon
!
! Create PILGRIM decompositions (see decompmodule)
!
      call decompcreate( 1, npr_y, xdist, ydist, strip2d )
      call decompcreate( 1, npr_y, npr_z, xdist, ydist, zdist, strip3dxyz )
      call decompcreate( "xzy", 1, npr_z, npr_y, xdist, zdist, ydist, strip3dxzy )


! For y communication within z subdomain (klast version)
      zdist1(1) = endlev-beglev+1
      call decompcreate( 1, npr_y, 1, xdist, ydist, zdist1, strip3yatz )

! For z communication within y subdomain

      ydistk(1) = endlat-beglat+1
      call decompcreate( 1, 1, npr_z, xdist, ydistk, zdist, strip3zaty )

! Arrays dimensioned plev+1

      zdist(npr_z) = zdist(npr_z) + 1
      call decompcreate( 1, npr_y, npr_z, xdist, ydist, zdist, strip3dxyzp )
      call decompcreate( "xzy", 1, npr_z, npr_y, xdist, zdist, ydist, strip3dxzyp )

! Arrays dimensioned plev+1, within y subdomain

      ydistk(1) = endlat-beglat+1
      call decompcreate( "xzy", 1, npr_z, 1, xdist, zdist, ydistk, strip3zatypt )

! For y communication within z subdomain (klast+1 version)
      zdist1(1) = endlev-beglev+2
      call decompcreate( 1, npr_y, 1, xdist, ydist, zdist1, strip3yatzp )

! Secondary xy decomposition
!
      if (twod_decomp == 1) then
        zdistxy(1) = plev
        call decompcreate( nprxy_x, nprxy_y, 1, xdistxy, ydistxy, zdistxy, strip3kxyz )
        call decompcreate( "xzy", nprxy_x, 1, nprxy_y, xdistxy, zdistxy, ydistxy, strip3kxzy )

        zdistxy(1) = zdistxy(1) + 1
        call decompcreate( nprxy_x, nprxy_y, 1, xdistxy, ydistxy, zdistxy, strip3kxyzp )
        call decompcreate( "xzy", nprxy_x, 1, nprxy_y, xdistxy, zdistxy, ydistxy, strip3kxzyp )

! Initialize transposes
!
      endif

!
! Do generic NCAR decomposition
!
      do procid=0,npes-1
         if (iam == 0) then
            write(6,*)'procid ',procid,' assigned ', &
                 cut(2,procid)-cut(1,procid)+1,' latitude values from', &
                 cut(1,procid),' through ',cut(2,procid)
         endif
!
! Determine which processor is responsible for the defined latitudes
!
         do lat=cut(1,procid),cut(2,procid)
            proc(lat) = procid
         end do
      end do
!
! Number of neighbor processors needed for boundary communication.  North
! first.
!
      isum = 0
      do procid=myid_y+1,npr_y-1
         nmostlat = cut(2,procid)
         isum = isum + cut(2,procid) - cut(1,procid) + 1
         if (isum >= numbnd) goto 20
      end do
20    if (myid_y /= npr_y-1 .and. isum < numbnd .and. nmostlat /= plat)then
         call endrun ('SPMDINIT_DYN: Something wrong in computation of northern neighbors')
      end if

      isum = 0
      do procid=myid_y-1,0,-1
         smostlat = cut(1,procid)
         isum = isum + cut(2,procid) - cut(1,procid) + 1
         if (isum >= numbnd) goto 30
      end do
30    if (myid_y /= 0 .and. isum < numbnd .and. smostlat /= 1) then
         call endrun ('SPMDINIT_DYN: Something wrong in computation of southern neighbors')
      end if

!      write(6,*)'-----------------------------------------'
!      write(6,*)'Number of lats passed north & south = ',numbnd
!      write(6,*)'Node  Partition'
!      write(6,*)'-----------------------------------------'
!      do procid=0,npes-1
!         write(6,200) procid,cut(1,procid),cut(2,procid)
!      end do
!      write(6,*)'iam=',iam,'Number of south neighbors needed for bndry exchange = ',neighs
!      write(6,*)'iam=',iam,'Number of north neighbors needed for bndry exchange = ',neighn

      deallocate (ydist)
      deallocate (zdist)

      return
!
! Formats
!
200   format(i3,4x,i3,'-',i3,7x,i3,'-',i3)

!EOC
   end subroutine spmdinit_dyn

!========================================================================

   subroutine decomp_wavenumbers
!----------------------------------------------------------------------- 
! 
! Purpose: partition the spectral work among the given number of processors
! 
! Method: Make the labor division as equal as possible given loop lengths
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      implicit none
      
      call endrun ('decomp_wavenumbers() should never be called in LR dynamics')

   end subroutine decomp_wavenumbers

   subroutine spmdbuf
!----------------------------------------------------------------------- 
! 
! Purpose: placeholder for buffer allocation routine 
! 
! Method: Make the labor division as equal as possible given loop lengths
! 
! Author: CCM Core Group
! 
!-----------------------------------------------------------------------
      implicit none
      
      call endrun ('spmdbuf() should never be called in LR dynamics')

   end subroutine spmdbuf

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
      integer, intent(in) :: numperlat            ! number of elements per latitude
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
         displs(p) = displs(p-1) + numperproc(p-1)
      end do
     
   end subroutine compute_gsfactors

#endif

end module spmd_dyn

