#include <misc.h>

module dmsbnd

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! This code does time interpolation for DMS boundary data in a netCDF
   ! file.  Assumptions on the data in the netCDF file are:
   ! 1. Coordinates are ordered (lon,lat,time)
   ! 2. The time coordinate is in days, and the data is assumed to be periodic
   !    annually.
   !
   ! It is assumed that the model calling this interface has been
   ! compiled so that 8 byte real data are being used.  On many
   ! machines this implies compiling with a "-r8" flag.
   ! 
   ! Author: B. Eaton
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid, only: pcols, begchunk, endchunk
   use pmgrid, only: masterproc, plon, plat
   use abortutils, only: endrun
!nf90   use netcdf

   implicit none
   save
   private
   public :: &
      dmsbndini,      &! initialize dmsbnd module
      dmsbndint,      &! interpolate dmsbnd data to requested date/time
      dmsbndget        ! return latitude slice data at current date/time

   ! public module data
   
   public :: dmsems_data ! full pathname for time-variant DMS emissions dataset
   character(len=256) :: dmsems_data

   ! private module data

#include <netcdf.inc>

   real(r8), allocatable, dimension(:) :: &
      time            ! time coordinate (calander days + frac)
   real(r8), allocatable :: &
      dmsin(:,:,:)    ! input data (pcols,begchunk:endchunk,2)
   real(r8), allocatable :: &
      dms(:,:)        ! interpolated data (pcols,begchunk:endchunk)

   integer :: &
      ncid,          &! ID for netCDF file
      nrec,          &! number of records (time samples)
      lotim,         &! time(lotim) .le. current time
      hitim,         &! current time .lt. time(hitim)
      loin,          &! index into input data array containing time(lotim) data
      hiin,          &! index into input data array containing time(hitim) data
      start(3),      &! start vector for netCDF hyperslabs
      count(3)        ! count vector for netCDF hyperslabs

!##############################################################################
contains
!##############################################################################

   subroutine dmsbndini( calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Open netCDF file containing DMS emissions data.  Initialize arrays
      ! with the data to be interpolated to the current time.
      !
      ! It is assumed that the time coordinate is increasing and represents
      ! calendar days (range = [1.,366.)).
      ! 
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: alloc_err, handle_ncerr, handle_err
      use phys_grid,      only: scatter_field_to_chunk
      use ioFileMod, only: getfil
      use filenames, only: bndtvdms

#if ( defined SPMD )
      use mpishorthand
#endif

!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         calday  ! current time in calendar days + fraction.

      ! Local variables:
      integer ::&
         did,   &
         istat, &
         recid, &! record ID
         nlon,  &
         vid
      real(r8), allocatable :: dmslatlon(:,:,:)  ! used for netCDF input
      !-----------------------------------------------------------------------

      start(1) = 1
      start(2) = 1
      start(3) = 1
      count(1) = plon
      count(2) = plat
      count(3) = 1

      ! Get file name.  
      call getfil(bndtvdms, dmsems_data, 0)

      ! Open file.
!nf90      call handle_ncerr( nf90_open( trim(dmsems_data), NF_NOWRITE, ncid ), &
!nf90         'dmsbndini: error opening file '//trim(dmsems_data) )

      if (masterproc) then

        call handle_ncerr( nf_open( trim(dmsems_data), NF_NOWRITE, ncid ), &
           'dmsbndini: error opening file '//trim(dmsems_data) )

      ! get the record id
!nf90      call handle_ncerr( nf90_inquire( ncid, unlimiteddimid=recid), &
!nf90         'dmsbndini: no record variables ' )

        call handle_ncerr( nf_inq_unlimdim( ncid, recid), &
           'dmsbndini: no record variables ' )

      ! Check that input data is a right resolution.
!nf90      call handle_ncerr( nf90_inq_dimid( ncid, 'lon', did ), 'dmsbndini: ' )
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, did, len=nlon ), 'dmsbndini: ' )

        call handle_ncerr( nf_inq_dimid( ncid, 'lon', did ), 'dmsbndini: ' )
        call handle_ncerr( nf_inq_dimlen( ncid, did, nlon ), 'dmsbndini: ' )
        if ( nlon .ne. plon ) then
           write(6,*)'dmsbndini: model plon = ',plon,', dataset nlon = ',nlon
           call endrun('dmsbndini:  plon != nlon')
        end if

      ! Get size of unlimited dimension.
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, recid, len=nrec ), 'dmsbndini: ' )
        call handle_ncerr( nf_inq_dimlen( ncid, recid, nrec ), 'dmsbndini: ' )

      endif  ! end of masterproc

#if ( defined SPMD )
      ! broadcast nrec to nodes
      call mpibcast (nrec, 1, mpiint, 0, mpicom)
#endif

      ! Allocate space for time coordinate data.
      allocate( time(nrec), stat=istat )
      call alloc_err( istat, 'dmsbndini', 'time', nrec )

      if (masterproc) then

      ! Get time coordinate.
!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'time', vid ), &
!nf90         'dmsbndini: cannot find time coordinate variable' )
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, time ), &
!nf90         'dmsbndini: error getting time coordinate data' )

        call handle_ncerr( nf_inq_varid( ncid, 'time', vid ), &
           'dmsbndini: cannot find time coordinate variable' )
        call handle_ncerr( nf_get_var_double( ncid, vid, time ), &
           'dmsbndini: error getting time coordinate data' )

      endif  ! end of masterproc

#if ( defined SPMD )
      ! broadcast time to nodes
      call mpibcast (time, nrec, mpir8, 0, mpicom)
#endif

      ! Make sure the time coordinate looks like calander day, and is
      ! increasing.
      call chktime( time, nrec )

      ! Find indices for time samples that bound the current time.
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      ! Read data.
      loin = 1
      hiin = 2

      allocate( dmsin(pcols,begchunk:endchunk,2), stat=istat )
      call alloc_err( istat, 'dmsbndini', 'dmsin', pcols*((endchunk-begchunk)+1)*2 )
      allocate( dms(pcols,begchunk:endchunk), stat=istat )
      call alloc_err( istat, 'dmsbndini', 'dms', pcols*((endchunk-begchunk)+1) )

      allocate( dmslatlon(plon,plat,2), stat=istat )
      call alloc_err( istat, 'dmsbndini', 'dmslatlon', &
                      plon*plat*2 )

      if (masterproc) then

!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'DMS', vid ), &
!nf90         'dmsbndini: cannot find variable '//'DMS' )
        call handle_ncerr( nf_inq_varid( ncid, 'DMS', vid ), &
           'dmsbndini: cannot find variable '//'DMS' )

        start(3) = lotim
!nf90      call handle_ncerr( nf90_get_var( ncid, vid,  dmslatlon(:,:,loin), start, count ), &
!nf90         'dmsbndini: cannot read data for '//'DMS' )
        call handle_ncerr( nf_get_vara_double( ncid, vid, start, count,  dmslatlon(:,:,loin) ), &
           'dmsbndini: cannot read data for '//'DMS' )

        start(3) = hitim
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, dmslatlon(:,:,hiin), start, count ), &
!nf90         'dmsbndini: cannot read data for '//'DMS' )
        call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, dmslatlon(:,:,hiin) ), &
           'dmsbndini: cannot read data for '//'DMS' )

      endif  ! end of masterproc

      ! scatter data to chunked data structures
      call scatter_field_to_chunk (1, 1, 2, plon, dmslatlon, &
                                   dmsin(1,begchunk,1))

      if (masterproc) then
        write(6,*)'dmsbndini: calendar day = ',calday, ' : read data for days ', &
           time(lotim), ' and ',time(hitim)
      endif  ! end of masterproc

      deallocate( dmslatlon, stat=istat )
      call handle_err( istat, &
         'ERROR deallocating memory for dmslatlon in routine dmsbndini')

   end subroutine dmsbndini

!#######################################################################

   subroutine dmsbndint( calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Interpolate DMS data to the current time.  Update the input data
      ! as necessary.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: alloc_err, handle_ncerr, handle_err
      use phys_grid,      only: scatter_field_to_chunk

!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      real(r8), intent(in) ::&
         calday  ! current time in calendar days + fraction.

      ! Local variables:
      integer ::&
         oldhitim,        &
         vid
      real(r8) ::&
         dt, dt1, tint
      real(r8), allocatable :: dmslatlon(:,:)  ! used for netCDF input
      integer istat
      !-----------------------------------------------------------------------

      ! Check to see if model time is still bounded by dataset times.
      oldhitim = hitim
      call findplb( time, nrec, calday, lotim )
      hitim = mod( lotim, nrec ) + 1

      if ( hitim .ne. oldhitim ) then
         ! Read in new hitim data.  Replace old lotim data.
         loin = hiin
         hiin = mod( loin, 2 ) + 1

         allocate( dmslatlon(plon,plat), stat=istat )
         call alloc_err( istat, 'dmsbndini', 'dmslatlon', &
                         plon*plat )

         if (masterproc) then

           start(3) = hitim
!nf90         call handle_ncerr( nf90_inq_varid( ncid, 'DMS', vid ), &
!nf90            'dmsbndint: cannot find variable '//'DMS' )
!nf90         call handle_ncerr( nf90_get_var( ncid, vid, dmslatlon, start, count ), &
!nf90            'dmsbndint: cannot read data for '//'DMS' )
           call handle_ncerr( nf_inq_varid( ncid, 'DMS', vid ), &
              'dmsbndint: cannot find variable '//'DMS' )
           call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, dmslatlon ), &
              'dmsbndint: cannot read data for '//'DMS' )
           write(6,*)'dmsbndint: read data for day ',time(hitim)

         endif  ! end of masterproc

         ! scatter data to chunked data structures
         call scatter_field_to_chunk (1, 1, 1, plon, dmslatlon, &
                                      dmsin(1,begchunk,hiin))
         if ( lotim .ne. oldhitim ) then
           if (masterproc) then
              ! Read in new lotim data.  Replace old hitim data.
              start(3) = lotim
!nf90            call handle_ncerr( nf90_inq_varid( ncid, 'DMS', vid ), &
!nf90               'dmsbndint: cannot find variable '//'DMS' )
!nf90            call handle_ncerr( nf90_get_var( ncid, vid, dmslatlon, start, count), &
!nf90               'dmsbndint: cannot read data for '//'DMS' )
              call handle_ncerr( nf_inq_varid( ncid, 'DMS', vid ), &
                 'dmsbndint: cannot find variable '//'DMS' )
              call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, dmslatlon), &
                 'dmsbndint: cannot read data for '//'DMS' )
              write(6,*)'dmsbndint: read data for day ',time(lotim)
           endif  ! end of masterproc

           ! scatter data to chunked data structures
           call scatter_field_to_chunk (1, 1, 1, plon, dmslatlon, &
                                        dmsin(1,begchunk,loin))

         end if

         deallocate( dmslatlon, stat=istat )
         call handle_err( istat, &
                 'ERROR deallocating memory for dmslatlon in routine dmsbndint')

      end if


      ! Linear interpolation...  Start by computing the number of days between
      !                          the lower and upper bounds, and days between
      !                          the model time and lower bound.

      if( time(hitim) .lt. time(lotim) )then
         dt = 365. - time(lotim) + time(hitim)
         if( calday .le. time(hitim) )then
            dt1 = 365. - time(lotim) + calday
         else
            dt1 = calday - time(lotim)
         end if
      else
         dt = time(hitim) - time(lotim)
         dt1 = calday - time(lotim)
      end if

      tint = dt1/dt
      call linintp( size(dms), 0._r8, 1._r8, tint, &
                    dmsin(1,begchunk,loin),        &
                    dmsin(1,begchunk,hiin),        &
                    dms )

   end subroutine dmsbndint

!#######################################################################

   subroutine dmsbndget( ncol, lchnk, x )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return DMS emission data for the requested latitude.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      integer, intent(in) :: ncol           ! number of columns used
      integer, intent(in) :: lchnk          ! chunk index

      real(r8), intent(out) :: x(pcols)     ! DMS emissions (kg DMS/m2/s)

      ! Local variables:
      integer ::         i
      !-----------------------------------------------------------------------

      do i = 1, ncol
         x(i) = dms(i,lchnk)
      end do

   end subroutine dmsbndget

!#######################################################################

end module dmsbnd
