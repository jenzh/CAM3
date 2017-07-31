#include <misc.h>
#include <params.h>
 
module carbonscales
!-----------------------------------------------------------------------
!
! Purposes:
!       Read in scaling dataset.
!       Provide (linearly interpolated) scaling when requested
!
! Public routines:
!       get_carbonscale
!       init_scale
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use time_manager, only: get_curr_calday, get_curr_date
   use pmgrid, only: masterproc
   use abortutils, only: endrun

   implicit none

   public init_scale      ! read file data
   public get_carbonscale ! return scale at current time step

   private 
   save

   real(r8),allocatable :: fscale(:)
   real(r8),allocatable :: days(:)  ! date converted to days
   integer :: indextopreviousdate
   integer :: numberoffiledata

contains

subroutine init_scale()
!-----------------------------------------------------------------------
!  read scale file into module variables
!-----------------------------------------------------------------------
   use ioFileMod, only: getfil
   use filenames, only: bndtvcarbonscale

#if ( defined SPMD )
   use mpishorthand
#endif

   character(len=256) :: locfn 
   integer :: scalenid = -1           ! netcdf id for scale file
   integer :: datedimid, scaleid, dateid
   integer,allocatable :: date(:)
   integer :: i
   integer :: dimids(1)

   integer :: yr, mon, day,sec
   real(r8) :: calday                        ! calendar day of current timestep
   real(r8) :: caldayloc                     ! calendar day of current timestep

   if (masterproc) then
!
! find and open file; abort if fail (getfil(,,0)).
!
      call getfil (bndtvcarbonscale, locfn, 0)
      call wrap_open (locfn, 0, scalenid)
      write(6,*)'init_scale: reading scale dataset'
!
! Get and check dimension info
!
      call wrap_inq_dimid( scalenid, 'date', datedimid )
      call wrap_inq_dimlen( scalenid, datedimid, numberoffiledata)
      allocate( date(numberoffiledata) )
      allocate( fscale(numberoffiledata) )

      call wrap_inq_varid( scalenid, 'carbonscale', scaleid )
      call wrap_inq_varid( scalenid, 'date', dateid )

      call wrap_inq_vardimid (scalenid, scaleid, dimids)
!
! read in Population and dates
!
      call wrap_get_var_realx (scalenid, scaleid, fscale)
      call wrap_get_var_int (scalenid, dateid, date)
 
#if ( defined SPMD )
      call mpibcast (numberoffiledata, 1, mpir8, 0, mpicom)
   else
      call mpibcast (numberoffiledata, 1, mpir8, 0, mpicom)
      allocate( date(numberoffiledata) )
      allocate( fscale(numberoffiledata) )
#endif
   endif

!
! broadcast scale and date to nodes
!
#if ( defined SPMD )
   call mpibcast (date, numberoffiledata, mpiint, 0, mpicom)
   call mpibcast (fscale, numberoffiledata, mpir8, 0, mpicom)
#endif

!
! convert date to days
!
   allocate( days(numberoffiledata) )
   do i = 1, numberoffiledata
     call bnddyi(date(i), 0, days(i))
     days(i) = days(i) + (date(i) / 10000) * 365
   enddo
   deallocate( date )

   calday = get_curr_calday ()
   call get_curr_date(yr, mon, day, sec)
   caldayloc = calday + yr*365
                                                                                
   if (caldayloc < days(1)) then
     write(6,*) 'init_scale: calday, caldayloc, yr', calday, caldayloc, yr
     write(6,*) 'init_scale: days(1)',days(1)
     call endrun ('init_scale: calendar manager returns date previous to scale data')
   endif

   if (caldayloc >=  days(numberoffiledata)) then
     write(6,*) 'init_scale: calday, caldayloc, yr', calday, caldayloc, yr
     write(6,*) 'init_scale: days(numberoffiledata)',days(numberoffiledata)
     call endrun ('init_scale: calendar manager returns date after all scale data')
   endif

!
! find data which bounds to the left
!
   do indextopreviousdate = 1, numberoffiledata
      if(caldayloc < days(indextopreviousdate)) exit
   enddo
   indextopreviousdate = indextopreviousdate - 1
  
   write(6,*) 'init_scale: caldayloc',caldayloc,'dateindex,date',indextopreviousdate,days

end subroutine init_scale

subroutine get_carbonscale(scaleattime)
!-----------------------------------------------------------------------
!  compute scale for time indicated by time_manager
!-----------------------------------------------------------------------

   real(r8), intent(out) :: scaleattime

   integer :: yr, mon, day,sec
   real(r8) :: calday                        ! calendar day of current timestep
   real(r8) :: caldayloc                     ! calendar day of current timestep
   real(r8) :: deltat, fact2, fact1          ! time interpolation factors
   real(r8) :: cdayp, cdaym

   calday = get_curr_calday ()
   call get_curr_date(yr, mon, day, sec)
   caldayloc = calday + yr*365

   if (caldayloc >=  days(indextopreviousdate+1)) then
      indextopreviousdate = indextopreviousdate + 1
      if (indextopreviousdate + 1 > numberoffiledata) then
          write(6,*) 'get_carbonscale: date exceeds scale data', caldayloc, days(numberoffiledata)
          call endrun
      endif
   endif

   if ( (caldayloc <  days(indextopreviousdate)) .or. caldayloc >= days(indextopreviousdate+1) ) then
      write(6,*) 'simulation day, bounding days',caldayloc, days(indextopreviousdate), days(indextopreviousdate+1)
      call endrun ('get_carbonscale: problem with data in scale file')
   endif

   deltat = days(indextopreviousdate+1) - days(indextopreviousdate)
   fact1 = (days(indextopreviousdate+1) - caldayloc) / deltat
   fact2 = (caldayloc - days(indextopreviousdate  )) / deltat
  
   scaleattime = fscale(indextopreviousdate  ) * fact1 + &
                 fscale(indextopreviousdate+1) * fact2

end subroutine get_carbonscale

end module carbonscales
