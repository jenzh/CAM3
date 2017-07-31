#include <misc.h>

module ramp_scon

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for initializing and computing 
!          solar constant.  It handles ramping of the solar constant.  
! 
! Author:
!   Original routines:  CMS
!   Module:             T. Henderson, January 2004
!
! $Id: ramp_scon.F90,v 1.3.2.4 2004/01/27 17:10:59 rosinski Exp $
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: masterproc
   use abortutils, only: endrun

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make default access private
   save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public rampnl_scon        ! Initialize runtime ramp options
   public ramp_sconst        ! Compute ramping of solar constant

!-----------------------------------------------------------------------
! Public data ----------------------------------------------------------
!-----------------------------------------------------------------------
   character(len=256), public :: bndtvscon    ! filename for ramped data

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------
   integer :: ntim = -1               ! number of yearly data values
   integer, allocatable :: yrdata(:)  ! yearly data values
   real(r8), allocatable :: sconst(:) ! input time-varying solar const (W/m2)
   logical fixYear_scon   ! true => Ramped gases fixed at specified year.
   integer rampYear_scon  ! ramped gases fixed at this year

!-----------------------------------------------------------------------
! Subroutines and functions --------------------------------------------
!-----------------------------------------------------------------------
contains


subroutine ramp_scon_read
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read ramped solar constant data
! 
! Author:  T. Henderson
! 
!-----------------------------------------------------------------------

   use ioFileMod, only: getfil
#if ( defined SPMD )
   use mpishorthand, only: mpicom, mpiint, mpir8
#endif

   include 'netcdf.inc'

!---------------------------Local variables-----------------------------
   integer :: ncid
   integer :: scon_id
   integer :: date_id
   integer :: time_id
   integer :: ierror
   character(len=256) :: locfn          ! netcdf local filename to open

   if (masterproc) then
     call getfil (bndtvscon, locfn, 0)
     call wrap_open (trim(locfn), NF_NOWRITE, ncid)
     write(6,*)'RAMP_SCON_READ:  reading ramped solar constant data from file ',trim(locfn)
     call wrap_inq_varid( ncid, 'date', date_id )
     call wrap_inq_varid( ncid, 'scon', scon_id )
     call wrap_inq_dimid( ncid, 'time', time_id )
     call wrap_inq_dimlen( ncid, time_id, ntim )
   endif
#if (defined SPMD )
   call mpibcast (ntim, 1, mpiint, 0, mpicom)
#endif
   ! these arrays are never deallocated
   allocate ( yrdata(ntim), sconst(ntim), stat=ierror )
   if (ierror /= 0) then
     call endrun ('RAMP_SCON_READ:  ERROR, allocate() failed')
   endif
   if (masterproc) then
     call wrap_get_var_int   (ncid, date_id, yrdata )
     yrdata = yrdata / 10000 
     call wrap_get_var_realx (ncid, scon_id, sconst )
     call wrap_close (ncid)
     write(6,*)'RAMP_SCON_READ:  successfully read ramped solar constant data from years ',yrdata(1),' through ',yrdata(ntim)
   endif
#if (defined SPMD )
   call mpibcast (sconst, ntim, mpir8, 0, mpicom)
   call mpibcast (yrdata, ntim, mpiint, 0, mpicom)
#endif

   return

end subroutine ramp_scon_read

!##############################################################################

subroutine rampnl_scon( year )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize the ramp options that are controlled by namelist input.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: <Who is primarily responsible for the code> 
! 
!-----------------------------------------------------------------------

   integer, intent(in) :: year ! Ramped gases fixed at this year

   ! read ramped solar constant data
   call ramp_scon_read
   rampYear_scon = year
   fixYear_scon = .false.
   if ( year > 0 ) then
      fixYear_scon = .true.
      if (masterproc) &
         write(6,*) 'RAMPNL_SCON:  Ramped gases being fixed at year ',rampYear_scon
   end if
   return
end subroutine rampnl_scon

!##############################################################################

subroutine ramp_sconst
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes ramping of solar constant
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: <Who is primarily responsible for the code> 
! 
!-----------------------------------------------------------------------
   use time_manager, only: get_curr_date, get_curr_calday
   use timeinterp,   only: validfactors

#include <comsol.h>
!---------------------------Local variables-----------------------------

   integer yrmodel           ! model year
   integer nyrm              ! year index
   integer nyrp              ! year index
   integer :: yr, mon, day   ! components of a date
   integer :: ncdate         ! current date in integer format [yyyymmdd]
   integer :: ncsec          ! current time of day [seconds]

   real(r8) :: calday            ! current calendar day
   real(r8) doymodel             ! model day of year
   real(r8) doydatam             ! day of year for input data yrdata(nyrm)
   real(r8) doydatap             ! day or year for input data yrdata(nyrp)
   real(r8) deltat               ! delta time
   real(r8) fact1, fact2         ! time interpolation factors
!
! ---------------------------------------------------------------------
!
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day
!
! determine index into input data
!
   if ( fixYear_scon ) then
      yrmodel  = rampYear_scon
   else
      yrmodel  = ncdate/10000
   end if

   nyrm       = yrmodel - yrdata(1) + 1
   nyrp       = nyrm + 1
!
! if current date is outside range of ramp values, then quit
!
   if ((nyrm < 1) .or. (nyrp > ntim)) then
      write(6,*)'nyrm = ',nyrm,' nyrp= ',nyrp, ' ncdate= ', ncdate
      call endrun ('RAMP_SCONST:  data time index is out of bounds')
   endif
!
! determine time interpolation factors, check sanity
! of interpolation factors to within 32-bit roundoff
! assume that day of year is 1 for all input data
!
   doymodel = yrmodel*365.    + calday
   doydatam = yrdata(nyrm)*365. + 1.
   doydatap = yrdata(nyrp)*365. + 1.
   deltat   = doydatap - doydatam
   fact1    = (doydatap - doymodel)/deltat
   fact2    = (doymodel - doydatam)/deltat

   if (.not. validfactors (fact1, fact2)) then
      write(6,*)'RAMP_SCONST: Bad fact1 and/or fact2=',fact1,fact2
      call endrun ()
   end if
!
! do time interpolation:
!
   scon = (sconst(nyrm)*fact1 + sconst(nyrp)*fact2)

   return
end subroutine ramp_sconst

end module ramp_scon

