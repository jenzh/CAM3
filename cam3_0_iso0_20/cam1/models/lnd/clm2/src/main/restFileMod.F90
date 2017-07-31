#include <misc.h>
#include <preproc.h>

module restFileMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: restFileMod
!
! !DESCRIPTION:
! Reads from or  writes to/ the CLM restart file.
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use spmdMod     , only : masterproc
  use abortutils,   only : endrun
!
! !PUBLIC TYPES:
  implicit none
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: restart
!
! !PRIVATE MEMBER FUNCTIONS:
  private :: restart_setup               !Setup restart files, perform consistency checks
  private :: restart_time                !Read/write restart time manager data
  private :: restart_biogeophys          !Read/write restart biogeophysics data
  private :: restart_wrapup              !Close restart file and write restart pointer file
  private :: write_rest_pfile            !Writes restart pointer file
  private :: set_restart_filename        !Sets restart filename
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !PRIVATE TYPES:
  private
  character(len=256) :: locfnr           !restart file name
  integer, parameter :: rest_id = 7      !restart id
  integer            :: rest_id_input    !restart id of input restart file
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart
!
! !INTERFACE:
  subroutine restart (flag)
!
! !DESCRIPTION:
! Read CLM restart file.
!
! !USES:
    use clm_varctl , only : nsrest, rpntdir, rpntfil, nrevsn, &
                            archive_dir, mss_irt, mss_wpass, &
                            csm_doflxave, caseid
    use accumulMod , only : restart_accum
#if (defined RTM)
    use RtmMod     , only : restart_rtm
#endif
    use histFileMod, only : restart_history
#if (defined DGVM)
    use DGVMRestMod, only : restart_dgvm
#endif
#if (defined COUP_CSM)
    use clm_csmMod , only : restart_coupler
#endif
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine driver
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: nio                         !Fortran unit number
!-----------------------------------------------------------------------

    if (masterproc) write(6,*)'in restart'

    if (flag /= 'read' .and. flag  /= 'write') then
       write(6,*)'ERROR: flag must be set to read or write'
       call endrun
    endif

    call restart_setup(nio, flag)

    call restart_time(nio, flag)

    call restart_biogeophys(nio, flag)

#if (defined COUP_CSM)
    call restart_coupler(nio, flag)
#endif

#if (defined RTM)
    call restart_rtm(nio, flag, rest_id_input)
#endif

#if (defined DGVM)
    call restart_dgvm(nio, flag)
#endif

    call restart_accum(nio, flag)

    call restart_history(nio, flag)

    call restart_wrapup(nio, flag)

  end subroutine restart

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_setup
!
! !INTERFACE:
  subroutine restart_setup (nio, flag)
!
! !DESCRIPTION:
! Setup restart file and perform necessary consistency checks
!
! !USES:
    use fileutils   , only : opnfil, getfil, getavu, relavu
    use time_manager, only : get_step_size, get_nstep
    use clm_varctl  , only : nsrest, nrevsn, rpntdir, rpntfil, caseid, &
                             brnch_retain_casename
!
! !ARGUMENTS:
    implicit none
    include 'netcdf.inc'
    integer, intent(out) :: nio             !restart unit
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer i,n                            !indices
    character(len=256) :: fnamer           !full name of restart file
    character(len= 32) :: casename         !case name read in from restart
!-----------------------------------------------------------------------

    if (masterproc) then

       if (flag == 'write') then
          ! Open restart file for writing
          write (6,*) 'Writout out restart data at nstep = ',get_nstep()
          write (6,'(72a1)') ("-",i=1,60)
          locfnr = set_restart_filename()
          nio = getavu()
          call opnfil (locfnr, nio, 'u')
       else
          ! Obtain the restart file for reading
          ! For restart runs, the restart pointer file contains the
          ! full pathname of the restart file. For branch runs, the
          ! namelist variable [nrevsn] contains the full pathname of the
          ! restart file.  New history files are created for branch runs.
          write (6,*) 'Reading in restart data .....'
          write (6,'(72a1)') ("-",i=1,60)
          if (nsrest == 1) then     !restart
             nio = getavu()
             locfnr = trim(rpntdir) //'/'// trim(rpntfil)
             call opnfil (locfnr, nio, 'f')
             read (nio,'(a256)') fnamer
             call relavu (nio)
          else                      !branch run
             fnamer = nrevsn
          end if
          nio = getavu()
          call getfil (fnamer, locfnr, 0)
          call opnfil (locfnr, nio, 'u')
       endif

       if (flag == 'write') then
          ! Write restart id
          write(nio) rest_id
       else
          ! Check restart id (restart file id should match restart code version id)
          read(nio) rest_id_input
          if (rest_id_input == rest_id) then
             write(6,*)'RESTRD: using restart file id ',rest_id
          else
             write(6,*)'RESTRD WARNING: input restart file id ',rest_id_input, &
                  ' does not match required restart id ',rest_id
          end if
       end if

       if (flag == 'write') then
          ! Write case name
          write(nio) caseid
       else
          ! Check case name consistency (case name must be different for branch run)
          read (nio) casename
          if (nsrest == 3 .and. casename==caseid .and. .not.(brnch_retain_casename)) then
             write(6,*) 'RESTRD ERROR: Must change case name on branch run'
             write(6,*) 'Prev case = ',trim(casename),' current case = ',trim(caseid)
             call endrun
          end if
       endif

    endif  ! masterproc if-block

  end subroutine restart_setup

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_time
!
! !INTERFACE:
  subroutine restart_time (nio, flag)
!
! !DESCRIPTION:
! Read/write time manager information to/from restart
!
! !USES:
#if (defined COUP_CSM)
    use controlMod, only : csm_dtime       !dtime from input namelist
#endif
    use time_manager, only : get_step_size, timemgr_write_restart, &
                             timemgr_read_restart, timemgr_restart, get_nstep
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             !restart unit
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
#if (defined COUP_CAM)
    integer :: clm_nstep                     !nstep from restart file
#endif
!-----------------------------------------------------------------------

#if (defined OFFLINE)

    if (flag == 'write') then
       if (masterproc) call timemgr_write_restart(nio)
    else
       if (masterproc) call timemgr_read_restart(nio)
       call timemgr_restart()
    end if

#elif (defined COUP_CAM)

    if (masterproc) then
       if (flag == 'write') then
          ! Write out time step

          write(nio) get_nstep()
       else
          ! Read in time step - check that clm restart time step is the
          ! same as cam restart time step. In cam mode, the time manager
          ! restart variables have already been read in by calls to
          ! routines timemgr_read_restart and timemgr_restart.

          read(nio) clm_nstep
          if ((clm_nstep+1) /= get_nstep()) then
             write(6,*)'RESTART_TIME: incompatibility in clm and cam restart dates'
             write(6,*)'  restart step from cam = ',get_nstep()
             write(6,*)'  restart step from clm = ',clm_nstep + 1
             call endrun
          endif
       endif
    endif

#elif (defined COUP_CSM)

    if (flag == 'write') then
       if (masterproc) call timemgr_write_restart(nio)
    else
       if (masterproc) call timemgr_read_restart(nio)
       call timemgr_restart()
    end if
    if (masterproc) then
       if (csm_dtime /= get_step_size()) then
          write(6,*)'(RESTRD): error '
          write(6,*)'namelist dtime on restart does not match input dtime'
          call endrun
       endif
    endif

#endif

  end subroutine restart_time

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_biogeophys
!
! !INTERFACE:
  subroutine restart_biogeophys (nio, flag)
!
! !DESCRIPTION:
! Read/Write biogeophysics information to/from restart file.
!
! !USES:
    use clmtype
    use iobinary
    use decompMod , only : get_proc_bounds, get_proc_global
    use clm_varpar, only : nlevsoi, numrad, lsmlon, lsmlat
    use clm_varcon, only : denice, denh2o
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: nio             !restart unit
    character(len=*), intent(in) :: flag  !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: g,l,c,p,j    ! indices
    real(r8):: pftsum       ! temporary used for pft averaging for columns
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    integer , pointer :: ibuf1dg(:)    ! temporary
    integer , pointer :: ibuf1dl(:)    ! temporary
    integer , pointer :: ibuf1dc(:)    ! temporary
    integer , pointer :: ibuf1dp(:)    ! temporary
    real(r8), pointer :: rbuf1dg(:)    ! temporary
    real(r8), pointer :: rbuf1dl(:)    ! temporary
    real(r8), pointer :: rbuf1dc(:)    ! temporary
    real(r8), pointer :: rbuf1dp(:)    ! temporary
    real(r8), pointer :: rbuf2dc(:,:)  ! temporary
    real(r8), pointer :: rbuf2dp(:,:)  ! temporary
!-----------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

   ! Allocate necessary 1d buffers

    allocate (ibuf1dg(numg))
    allocate (ibuf1dl(numl))
    allocate (ibuf1dc(numc))
    allocate (ibuf1dp(nump))

    allocate (rbuf1dg(numg))
    allocate (rbuf1dl(numl))
    allocate (rbuf1dc(numc))
    allocate (rbuf1dp(nump))

    ! Column physical state - snl
    if (flag == 'read') call readin (nio, ibuf1dc, clmlevel=namec)
!dir$ concurrent
!cdir nodep
    do c = begc,endc
       if (flag == 'read' ) clm3%g%l%c%cps%snl(c) = ibuf1dc(c)
       if (flag == 'write') ibuf1dc(c) = clm3%g%l%c%cps%snl(c)
    end do
    if (flag == 'write') call wrtout (nio, ibuf1dc, clmlevel=namec)

    ! Column physical state - snowdp
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=namec)
!dir$ concurrent
!cdir nodep
    do c = begc,endc
       if (flag == 'read' ) clm3%g%l%c%cps%snowdp(c) = rbuf1dc(c)
       if (flag == 'write') rbuf1dc(c) = clm3%g%l%c%cps%snowdp(c)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=namec)

    ! Column physical state - snowage
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=namec)
!dir$ concurrent
!cdir nodep
    do c = begc,endc
       if (flag == 'read' ) clm3%g%l%c%cps%snowage(c) = rbuf1dc(c)
       if (flag == 'write') rbuf1dc(c) = clm3%g%l%c%cps%snowage(c)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=namec)

    ! Column physical state - frac_sno
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=namec)
!dir$ concurrent
!cdir nodep
    do c = begc,endc
       if (flag == 'read' ) clm3%g%l%c%cps%frac_sno(c) = rbuf1dc(c)
       if (flag == 'write') rbuf1dc(c) = clm3%g%l%c%cps%frac_sno(c)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=namec)

    ! Column physical state - dz (snow)
    allocate (rbuf2dc(-nlevsno+1:0,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=namec)
    do j = -nlevsno+1,0
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (flag == 'read' ) clm3%g%l%c%cps%dz(c,j) = rbuf2dc(j,c)
          if (flag == 'write') rbuf2dc(j,c) = clm3%g%l%c%cps%dz(c,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=namec)
    deallocate(rbuf2dc)

    ! Column physical state - z (snow)
    allocate (rbuf2dc(-nlevsno+1:0,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=namec)
    do j = -nlevsno+1,0
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (flag == 'read' ) clm3%g%l%c%cps%z(c,j) = rbuf2dc(j,c)
          if (flag == 'write') rbuf2dc(j,c) = clm3%g%l%c%cps%z(c,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=namec)
    deallocate(rbuf2dc)

    allocate (rbuf2dc(-nlevsno:0,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=namec)
    do j = -nlevsno,0
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (flag == 'read' ) clm3%g%l%c%cps%zi(c,j) = rbuf2dc(j,c)
          if (flag == 'write') rbuf2dc(j,c) = clm3%g%l%c%cps%zi(c,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=namec)
    deallocate(rbuf2dc)

    !pft type physical state variable - albd
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=namep)
    do j = 1,numrad
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          if (flag == 'read' ) clm3%g%l%c%p%pps%albd(p,j) = rbuf2dp(j,p)
          if (flag == 'write') rbuf2dp(j,p) = clm3%g%l%c%p%pps%albd(p,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=namep)
    deallocate(rbuf2dp)

    !pft type physical state variable - albi
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=namep)
    do j = 1,numrad
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          if (flag == 'read' ) clm3%g%l%c%p%pps%albi(p,j) = rbuf2dp(j,p)
          if (flag == 'write') rbuf2dp(j,p) = clm3%g%l%c%p%pps%albi(p,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=namep)
    deallocate(rbuf2dp)

    !column type physical state variable - albgrd
    allocate (rbuf2dc(numrad,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=namec)
    do j = 1,numrad
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (flag == 'read' ) clm3%g%l%c%cps%albgrd(c,j) = rbuf2dc(j,c)
          if (flag == 'write') rbuf2dc(j,c) = clm3%g%l%c%cps%albgrd(c,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=namec)
    deallocate(rbuf2dc)

    !column type physical state variable - albgri
    allocate (rbuf2dc(numrad,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=namec)
    do j = 1,numrad
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (flag == 'read' ) clm3%g%l%c%cps%albgri(c,j) = rbuf2dc(j,c)
          if (flag == 'write') rbuf2dc(j,c) = clm3%g%l%c%cps%albgri(c,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=namec)
    deallocate(rbuf2dc)

   ! column water state variable - h2osno
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=namec)
!dir$ concurrent
!cdir nodep
    do c = begc,endc
       if (flag == 'read' ) clm3%g%l%c%cws%h2osno(c) = rbuf1dc(c)
       if (flag == 'write') rbuf1dc(c) = clm3%g%l%c%cws%h2osno(c)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=namec)

   ! column water state variable - h2osoi_liq
    allocate (rbuf2dc(-nlevsno+1:nlevsoi,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=namec)
    do j = -nlevsno+1,nlevsoi
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (flag == 'read' ) clm3%g%l%c%cws%h2osoi_liq(c,j) = rbuf2dc(j,c)
          if (flag == 'write') rbuf2dc(j,c) = clm3%g%l%c%cws%h2osoi_liq(c,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=namec)
    deallocate(rbuf2dc)

   ! column water state variable - h2osoi_ice
    allocate (rbuf2dc(-nlevsno+1:nlevsoi,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=namec)
    do j = -nlevsno+1,nlevsoi
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (flag == 'read' ) clm3%g%l%c%cws%h2osoi_ice(c,j) = rbuf2dc(j,c)
          if (flag == 'write') rbuf2dc(j,c) = clm3%g%l%c%cws%h2osoi_ice(c,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=namec)
    deallocate(rbuf2dc)

   ! column energy state variable - t_grnd
    if (flag == 'read') call readin (nio, rbuf1dc, clmlevel=namec)
!dir$ concurrent
!cdir nodep
    do c = begc,endc
       if (flag == 'read' ) clm3%g%l%c%ces%t_grnd(c) = rbuf1dc(c)
       if (flag == 'write') rbuf1dc(c) = clm3%g%l%c%ces%t_grnd(c)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dc, clmlevel=namec)

   ! pft energy state variable - t_ref2m_min
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pes%t_ref2m_min(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pes%t_ref2m_min(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

   ! pft energy state variable - t_ref2m_max
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pes%t_ref2m_max(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pes%t_ref2m_max(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

   ! pft energy state variable - t_ref2m_min_inst
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pes%t_ref2m_min_inst(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pes%t_ref2m_min_inst(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

   ! pft energy state variable - t_ref2m_max_inst
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pes%t_ref2m_max_inst(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pes%t_ref2m_max_inst(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

   ! column energy state variable - t_soisno
    allocate (rbuf2dc(-nlevsno+1:nlevsoi,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=namec)
    do j = -nlevsno+1,nlevsoi
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (flag == 'read' ) clm3%g%l%c%ces%t_soisno(c,j) = rbuf2dc(j,c)
          if (flag == 'write') rbuf2dc(j,c) = clm3%g%l%c%ces%t_soisno(c,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=namec)
    deallocate(rbuf2dc)

    !column type energy state variable - t_lake
    allocate (rbuf2dc(1:nlevlak,numc))
    if (flag == 'read') call readin (nio, rbuf2dc, clmlevel=namec)
    do j = 1,nlevlak
!dir$ concurrent
!cdir nodep
       do c = begc,endc
          if (flag == 'read' ) clm3%g%l%c%ces%t_lake(c,j) = rbuf2dc(j,c)
          if (flag == 'write') rbuf2dc(j,c) = clm3%g%l%c%ces%t_lake(c,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dc, clmlevel=namec)
    deallocate(rbuf2dc)

    ! pft type physical state variable - frac_veg_nosno_alb
    if (flag == 'read') call readin (nio, ibuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pps%frac_veg_nosno_alb(p) = ibuf1dp(p)
       if (flag == 'write') ibuf1dp(p) = clm3%g%l%c%p%pps%frac_veg_nosno_alb(p)
    end do
    if (flag == 'write') call wrtout (nio, ibuf1dp, clmlevel=namep)

    ! pft type physical state variable - fwet
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pps%fwet(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pps%fwet(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type physical state variable - tlai
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pps%tlai(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pps%tlai(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type physical state variable - tsai
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pps%tsai(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pps%tsai(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type physical state variable - elai
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pps%elai(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pps%elai(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type physical state variable - esai
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pps%esai(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p)= clm3%g%l%c%p%pps%esai(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type physical state variable - fsun
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read') clm3%g%l%c%p%pps%fsun(p) = rbuf1dp(p)
       if (flag == 'write')rbuf1dp(p)= clm3%g%l%c%p%pps%fsun(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type physical state variable - htop
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pps%htop(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p)= clm3%g%l%c%p%pps%htop(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type physical state variable - hbot
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pps%hbot(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p)= clm3%g%l%c%p%pps%hbot(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type physical state variable - fabd
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=namep)
    do j = 1,numrad
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          if (flag == 'read' ) clm3%g%l%c%p%pps%fabd(p,j) = rbuf2dp(j,p)
          if (flag == 'write') rbuf2dp(j,p) = clm3%g%l%c%p%pps%fabd(p,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=namep)
    deallocate(rbuf2dp)

    ! pft type physical state variable - fabi
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=namep)
    do j = 1,numrad
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          if (flag == 'read' ) clm3%g%l%c%p%pps%fabi(p,j) = rbuf2dp(j,p)
          if (flag == 'write') rbuf2dp(j,p) = clm3%g%l%c%p%pps%fabi(p,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=namep)
    deallocate(rbuf2dp)

    ! pft type physical state variable - ftdd
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=namep)
    do j = 1,numrad
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          if (flag == 'read' ) clm3%g%l%c%p%pps%ftdd(p,j) = rbuf2dp(j,p)
          if (flag == 'write') rbuf2dp(j,p) = clm3%g%l%c%p%pps%ftdd(p,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=namep)
    deallocate(rbuf2dp)

    ! pft type physical state variable - ftid
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=namep)
    do j = 1,numrad
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          if (flag == 'read' ) clm3%g%l%c%p%pps%ftid(p,j) = rbuf2dp(j,p)
          if (flag == 'write') rbuf2dp(j,p) = clm3%g%l%c%p%pps%ftid(p,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=namep)
    deallocate(rbuf2dp)

    ! pft type physical state variable - ftii
    allocate (rbuf2dp(numrad,nump))
    if (flag == 'read') call readin (nio, rbuf2dp, clmlevel=namep)
    do j = 1,numrad
!dir$ concurrent
!cdir nodep
       do p = begp,endp
          if (flag == 'read' ) clm3%g%l%c%p%pps%ftii(p,j) = rbuf2dp(j,p)
          if (flag == 'write') rbuf2dp(j,p) = clm3%g%l%c%p%pps%ftii(p,j)
       end do
    end do
    if (flag == 'write') call wrtout (nio, rbuf2dp, clmlevel=namep)
    deallocate(rbuf2dp)

    ! pft type energy state variable - t_veg
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pes%t_veg(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pes%t_veg(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

    ! pft type water state variable - h2ocan
    if (flag == 'read') call readin (nio, rbuf1dp, clmlevel=namep)
!dir$ concurrent
!cdir nodep
    do p = begp,endp
       if (flag == 'read' ) clm3%g%l%c%p%pws%h2ocan(p) = rbuf1dp(p)
       if (flag == 'write') rbuf1dp(p) = clm3%g%l%c%p%pws%h2ocan(p)
    end do
    if (flag == 'write') call wrtout (nio, rbuf1dp, clmlevel=namep)

#if (!defined DGVM)
    ! For read only:
    ! Determine average over all column pfts for h2ocan, needed by begwb
    ! computation in routine driver.F90) - this needs to be done after the
    ! weights are reset in the DGVM case
    ! The following should not be vectorized
    if (flag == 'read' ) then
       do c = begc,endc
          clm3%g%l%c%cws%pws_a%h2ocan(c) = 0.
       end do
       do p = begp,endp
          c = clm3%g%l%c%p%column(p)
          clm3%g%l%c%cws%pws_a%h2ocan(c) =  clm3%g%l%c%cws%pws_a%h2ocan(c) &
               + clm3%g%l%c%p%pws%h2ocan(p) * clm3%g%l%c%p%wtcol(p)
       end do
    end if
#endif

    ! For read only:
    ! Determine volumetric soil water
    if (flag == 'read' ) then
       do j = 1,nlevsoi
!dir$ concurrent
!cdir nodep
          do c = begc,endc
             clm3%g%l%c%cws%h2osoi_vol(c,j) = &
                  clm3%g%l%c%cws%h2osoi_liq(c,j)/(clm3%g%l%c%cps%dz(c,j)*denh2o) &
                + clm3%g%l%c%cws%h2osoi_ice(c,j)/(clm3%g%l%c%cps%dz(c,j)*denice)
          end do
       end do
    endif

    deallocate (ibuf1dg)
    deallocate (ibuf1dl)
    deallocate (ibuf1dc)
    deallocate (ibuf1dp)

    deallocate (rbuf1dg)
    deallocate (rbuf1dl)
    deallocate (rbuf1dc)
    deallocate (rbuf1dp)

  end subroutine restart_biogeophys

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: restart_wrapup
!
! !INTERFACE:
  subroutine restart_wrapup (nio, flag)
!
! !DESCRIPTION:
! Close and archive restart file and write restart pointer file if
! in write mode, otherwise just close restart file if in read mode
!
! !USES:
    use clm_varctl, only : mss_irt, mss_wpass, archive_dir
#if (defined COUP_CSM)
    use clm_csmMod, only : csmstop_next
#endif
    use fileutils, only : putfil, relavu, set_filename
    use time_manager, only : is_last_step
!
! !ARGUMENTS:
    implicit none
    integer, intent(inout) :: nio          !restart unit
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: i                   !index
    logical :: lremove             !true => remove file after archive
    character(len=256) :: rem_fn   !remote (archive) filename
    character(len=256) :: rem_dir  !remote (archive) directory
!-----------------------------------------------------------------------

   if (masterproc) then
      if (flag == 'write') then
         ! close and archive restart file and write restart pointer file
         lremove = .true.
#if (defined OFFLINE) || (defined COUP_CAM)
         if (is_last_step()) lremove = .false.
#elif (defined COUP_CSM)
         if (csmstop_next) lremove = .false.
#endif
         call relavu (nio)
         write(6,*) 'Successfully wrote local restart file ',trim(locfnr)
         if (mss_irt > 0) then
            rem_dir = trim(archive_dir) // '/rest/'
            rem_fn = set_filename(rem_dir, locfnr)
            call putfil (locfnr, rem_fn, mss_wpass, mss_irt, lremove)
         endif
         call write_rest_pfile( )
         write(6,'(72a1)') ("-",i=1,60)
         write(6,*)
      else
         ! close file
         call relavu (nio)
         write(6,'(72a1)') ("-",i=1,60)
         write(6,*) 'Successfully read restart data for restart run'
         write(6,*)
      endif
   endif

  end subroutine restart_wrapup

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: write_rest_pfile
!
! !INTERFACE:
  subroutine write_rest_pfile ( )
!
! !DESCRIPTION:
! Open restart pointer file. Write names of current restart and
! history files. If using mass store, these are the mass store
! names except if mss_irt=0 (no mass store files written). Close.
!
! !USES:
    use clm_varctl, only : rpntdir, mss_irt, archive_dir, rpntfil
    use fileutils, only : set_filename, relavu
    use fileutils, only : getavu, opnfil
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: m                           !index
    integer :: nio                         !restart pointer file
    character(len=256) :: filename         !local file name
    character(len=256) :: rem_dir          !remote directory
!-----------------------------------------------------------------------

    nio = getavu()
    filename= trim(rpntdir) //'/'// trim(rpntfil)
    call opnfil (filename, nio, 'f')

    ! write name of restart file to pointer file
    if (mss_irt == 0) then
       write(nio,'(a)') locfnr
    else
       rem_dir = trim(archive_dir) // '/rest/'
       write(nio,'(a)') set_filename(rem_dir, locfnr)
    endif

    ! add comments to pointer file of all files that are needed for restart
    write(nio,*)'The following lines list files needed for restart - do not edit'

    ! only write names of open history files
!    do m = 1,nhist
!       if (locfnh(m) /= ' ') then
!          if (mss_irt == 0) then
!             filename = locfnh(m)
!          else
!             rem_dir = trim(archive_dir) // '/hist/'
!             filename = set_filename(rem_dir, locfnh(m))
!          endif
!          write(nio, '(a)') filename
!       end if
!    end do

    call relavu (nio)
    write(6,*)'Successfully wrote local restart pointer file'

  end subroutine write_rest_pfile

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: set_restart_filename
!
! !INTERFACE:
  character(len=256) function set_restart_filename ()
!
! !DESCRIPTION:
!
! !USES:
    use clm_varctl  , only : caseid
    use time_manager, only : get_curr_date
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine restart in this module
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: cdate       !date char string
    integer :: day                    !day (1 -> 31)
    integer :: mon                    !month (1 -> 12)
    integer :: yr                     !year (0 -> ...)
    integer :: sec                    !seconds into current day
!-----------------------------------------------------------------------

    call get_curr_date (yr, mon, day, sec)
    write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    set_restart_filename = "./"//trim(caseid)//".clm2.r."//trim(cdate)

  end function set_restart_filename

end module restFileMod



