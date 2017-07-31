#include <misc.h>

module timeinterp
   
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils, only: endrun

   private
   save

   public getfactors
   public validfactors

CONTAINS

subroutine getfactors (cycflag, np1, cdayminus, cdayplus, cday, &
                       fact1, fact2, str)
!---------------------------------------------------------------------------
!
! Purpose: Determine time interpolation factors (normally for a boundary dataset)
!          for linear interpolation.
!
! Method:  Assume 365 days per year.  Output variable fact1 will be the weight to
!          apply to data at calendar time "cdayminus", and fact2 the weight to apply
!          to data at time "cdayplus".  Combining these values will produce a result 
!          valid at time "cday".  Output arguments fact1 and fact2 will be between 
!          0 and 1, and fact1 + fact2 = 1 to roundoff.
!
! Author:  Jim Rosinski
!
!--------------------------------------------------------------------------- 
   implicit none
!
! Arguments
!
   logical, intent(in) :: cycflag             ! flag indicates whether dataset is being cycled yearly

   integer, intent(in) :: np1                 ! index points to forward time slice matching cdayplus

   real(r8), intent(in) :: cdayminus          ! calendar day of rearward time slice
   real(r8), intent(in) :: cdayplus           ! calendar day of forward time slice
   real(r8), intent(in) :: cday               ! calenar day to be interpolated to
   real(r8), intent(out) :: fact1             ! time interpolation factor to apply to rearward time slice
   real(r8), intent(out) :: fact2             ! time interpolation factor to apply to forward time slice

   character(len=*), intent(in) :: str        ! string to be added to print in case of error (normally the callers name)
!
! Local workspace
!
   real(r8) :: deltat                         ! time difference (days) between cdayminus and cdayplus
   real(r8), parameter :: daysperyear = 365.  ! number of days in a year
!
! Initial sanity checks
!
   if (np1 == 1 .and. .not. cycflag) then
      call endrun ('GETFACTORS:'//str//' cycflag false and forward month index = Jan. not allowed')
   end if

   if (np1 < 1) then
      call endrun ('GETFACTORS:'//str//' input arg np1 must be > 0')
   end if

   if (cycflag) then
      if ((cday < 1.) .or. (cday > (daysperyear+1.))) then
         write(6,*) 'GETFACTORS:', str, ' bad cday=',cday
         call endrun ()
      end if
   else
      if (cday < 1.) then
         write(6,*) 'GETFACTORS:', str, ' bad cday=',cday
         call endrun ()
      end if
   end if         
!
! Determine time interpolation factors.  Account for December-January 
! interpolation if dataset is being cycled yearly.
!
   if (cycflag .and. np1 == 1) then                     ! Dec-Jan interpolation
      deltat = cdayplus + daysperyear - cdayminus
      if (cday > cdayplus) then                         ! We are in December
         fact1 = (cdayplus + daysperyear - cday)/deltat
         fact2 = (cday - cdayminus)/deltat
      else                                              ! We are in January
         fact1 = (cdayplus - cday)/deltat
         fact2 = (cday + daysperyear - cdayminus)/deltat
      end if
   else
      deltat = cdayplus - cdayminus
      fact1 = (cdayplus - cday)/deltat
      fact2 = (cday - cdayminus)/deltat
   end if

   if (.not. validfactors (fact1, fact2)) then
      write(6,*) 'GETFACTORS: ', str, ' bad fact1 and/or fact2=', fact1, fact2
      call endrun ()
   end if

   return
end subroutine getfactors

logical function validfactors (fact1, fact2)
!---------------------------------------------------------------------------
!
! Purpose: check sanity of time interpolation factors to within 32-bit roundoff
!
!---------------------------------------------------------------------------
   implicit none

   real(r8), intent(in) :: fact1, fact2           ! time interpolation factors

   validfactors = .true.

   ! The fact1 .ne. fact1 and fact2 .ne. fact2 comparisons are to detect NaNs.
   if (abs(fact1+fact2-1.) > 1.e-6 .or. &
       fact1 > 1.000001 .or. fact1 < -1.e-6 .or. &
       fact2 > 1.000001 .or. fact2 < -1.e-6 .or. &
       fact1 .ne. fact1 .or. fact2 .ne. fact2) then

      validfactors = .false.
   end if
   
   return
end function validfactors
end module timeinterp
