!===============================================================================
! CVS $Id: shr_msg_mod.F90,v 1.2.4.1 2004/01/02 18:50:55 mvr Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/shared/csm_share/shr/shr_msg_mod.F90,v $
! CVS $Name: cam3_0_brnchT_release01 $
!===============================================================================

MODULE shr_msg_mod

   use shr_kind_mod   ! defines real & integer kinds
   use shr_sys_mod    ! shared system call wrappers

   IMPLICIT none
!
! PUBLIC: Interfaces and data
!
   public :: shr_msg_stdio
   public :: shr_msg_chdir
   public :: shr_msg_dirio

!===============================================================================
CONTAINS
!===============================================================================

SUBROUTINE shr_msg_stdio(model)

   !--- arguments ---
   character(len=*),intent(in) :: model ! used to construct env varible name

   !--- formats ---
   character(len=*),parameter :: subName = "(shr_msg_stdio)"
   character(len=*),parameter :: F00 = "('(shr_msg_stdio) ',4a)"

!-------------------------------------------------------------------------------
! PURPOSE:
!   1) change the cwd (current working directory) and 
!   2) redirect stdin & stdout (units 5 & 6) to named files,
!   where the desired cwd & files are specified by namelist file.
!
!   Normally this is done to work around limitations in the execution syntax
!   of common MPI implementations.  For example, SGI's mpirun syntax is not 
!   flexible enough to allow MPMD models to select different execution
!   directories or to redirect stdin & stdout on the command line.  
!   Such functionality is highly desireable for CCSM purposes.  
!   ie. mpirun can't handle this:
!   unix> cd /usr/tmp/jdoe/csm/case01/atm ; atm < atm.parm > atm.log &
!   unix> cd /usr/tmp/jdoe/csm/case01/cpl ; cpl < cpl.parm > cpl.log &
!   etc.
!
! ASSUMPTIONS:
! o if the cwd, stdin, or stdout are to be changed, there must be a namelist
!   file in the cwd named <model>_stdio.nml  where <model> is provided via
!   subroutine dummy argument. 
!-------------------------------------------------------------------------------

   call shr_msg_chdir(model)   ! changes cwd
   call shr_msg_dirio(model)   ! open units 5 & 6 to named files
 
END SUBROUTINE shr_msg_stdio

!===============================================================================

SUBROUTINE shr_msg_chdir(model)

   !--- arguments ---
   character(len=*),intent(in) :: model ! used to construct env varible name

   !--- local ---
   character(256) :: filename ! namelist file to read
   logical        :: exists   ! true iff file exists
   character(256) :: dir      ! directory to cd to
   character(256) :: stdin    ! open unit 5 to this file
   character(256) :: stdout   ! open unit 6 to this file
   integer(SHR_KIND_IN) :: rcode    ! return code

   namelist / stdio / dir,stdin,stdout

   !--- formats ---
   character(len=*),parameter :: subName = "(shr_msg_chdir)"
   character(len=*),parameter :: F00 = "('(shr_msg_chdir) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_msg_chdir) ',2a,i6)"

!-------------------------------------------------------------------------------
! PURPOSE:
!   change the cwd (current working directory), see shr_msg_stdio for notes
!-------------------------------------------------------------------------------

   dir    = "nochange"
   stdin  = "nochange"
   stdout = "nochange"

   filename = trim(model)//"_stdio.nml"  ! eg. file="cpl_stdio.nml"
   inquire(file=filename,exist=exists)

   if (.not. exists) then
      write(6,F01)
      write(6,F00) "file ",trim(filename),& 
      & " doesn't exist, cwd has *not* been changed"
   else
      open (10,file=filename)
      read (10,nml=stdio,iostat=rcode)
      close(10)
      if (rcode > 0) then
         write(6,F01) 'ERROR: reading ',trim(filename),', iostat=',rcode
         call shr_sys_abort(subName//" ERROR reading "//trim(filename) )
      else if (dir /= "nochange") then
         call shr_sys_chdir(dir ,rcode)
         write(6,F00) "read ",trim(filename),", changed cwd to ",trim(dir)
      else
         write(6,F00) "read ",trim(filename),", cwd has *not* been changed"
      endif
   endif
 
   call shr_sys_flush(6)
 
END SUBROUTINE shr_msg_chdir

!===============================================================================

SUBROUTINE shr_msg_dirio(model)

   !--- arguments ---
   character(len=*),intent(in) :: model ! used to construct env varible name

   !--- local ---
   character(256) :: filename ! namelist file to read
   logical        :: exists   ! true iff file exists
   character(256) :: dir      ! directory to cd to
   character(256) :: stdin    ! open unit 5 to this file
   character(256) :: stdout   ! open unit 6 to this file
   integer(SHR_KIND_IN) :: rcode    ! return code

   namelist / stdio / dir,stdin,stdout

   !--- formats ---
   character(len=*),parameter :: subName = "(shr_msg_dirio)"
   character(len=*),parameter :: F00 = "('(shr_msg_dirio) ',4a)"
   character(len=*),parameter :: F01 = "('(shr_msg_dirio) ',2a,i6)"

!-------------------------------------------------------------------------------
! PURPOSE:
!   change the stdin & stdout (units 5 & 6), see shr_msg_stdio for notes
!-------------------------------------------------------------------------------

   dir    = "nochange"
   stdin  = "nochange"
   stdout = "nochange"

   filename = trim(model)//"_stdio.nml"  ! eg. file="cpl_stdio.nml"
   inquire(file=filename,exist=exists)

   if (.not. exists) then
      write(6,F00) "file ",trim(filename),& 
      & " doesn't exist, units 5 & 6 have *not* been changed"
   else
      open (10,file=filename)
      read (10,nml=stdio,iostat=rcode)
      close(10)
      if (rcode > 0) then
         write(6,F01) 'ERROR: reading ',trim(filename),', iostat=',rcode
         call shr_sys_abort(subName//" ERROR reading "//trim(filename) )
      else
         if (stdout /= "nochange") then
            open(unit=6,file=stdout,position='APPEND')
            write(6,F00) "read ",trim(filename),', unit 6 connected to ',trim(stdout)
         else
            write(6,F00) "read ",trim(filename),', unit 6 has *not* been redirected'
         endif
         if (stdin  /= "nochange") then
            open(unit=5,file=stdin ,status='UNKNOWN')
            write(6,F00) "read ",trim(filename),', unit 5 connected to ',trim(stdin)
         else
            write(6,F00) "read ",trim(filename),', unit 5 has *not* been redirected'
         endif
      endif
   endif
 
   call shr_sys_flush(6)
 
END SUBROUTINE shr_msg_dirio

!===============================================================================

END MODULE shr_msg_mod




