#include <misc.h>

module soxbnd

   !----------------------------------------------------------------------- 
   ! Purpose: 
   ! Interpolate the SOX emissions data from the GEIA-SMITH datasets.  This
   ! dataset contains seasonal average data for various years at roughly 15
   ! year intervals.  The interpolation scheme does linear interpolations to
   ! the requested calendar day in the seasonal cycle for each of
   ! the two years in the dataset that bound the requested year.  Finally a 
   ! linear interpolation is done to the requested year.
   !
   ! It is assumed that the model calling this interface has been
   ! compiled so that 8 byte real data are being used.  On many
   ! machines this implies compiling with a "-r8" flag.
   ! 
   ! Author: B. Eaton
   !
   ! Modified 11 Jan 2000 PJR: work with 4 or 8 byte floats 
   ! Modified 15 Jan 2004 D Bundy: can cycle one year of emission data
   !          using namelist variable rampyear_prognostic_sulfur
   !
   !----------------------------------------------------------------------- 

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: masterproc, plon, plat
   use ppgrid, only: pcols, begchunk, endchunk
   use abortutils, only: endrun
   use infnan, only: bigint

!nf90   use netcdf

   implicit none
   save
   private
   public :: &
      soxbndini,      &! initialize soxbnd module
      soxbndint,      &! interpolate soxbnd data to requested date/time
      soxbndget        ! return latitude slice data at current date/time

   ! public module data
   
   !  Variables set by namelist
   character(len=256) :: soxems_data
!++drb
   integer, public :: rampyear_prognostic_sulfur = bigint  ! if = bigint => unset, else = YYYY to cycle
   character(len=16), public :: scenario_prognostic_sulfur = 'RAMPED' ! only option available

   logical :: doramp_sox   ! at the moment this is the only option
   logical :: fixyear_sox  ! if true then cycle on year from the emissions dataset
!--drb

#include <netcdf.inc>

   integer, parameter ::&
      ndlev=2            ! number of levels in emissions data

   integer, allocatable, dimension(:) :: &
      date,             &! date coordinate (yyyymmdd format)
      yr                 ! years of annual cycle data
   real(r8), allocatable, dimension(:) :: &
      cdays              ! mid-season calendar days (fractional)
   real(r8), allocatable, dimension(:,:,:,:) :: &
      soxyrlo,          &! SOX input data for lower bound year (pcols,begchunk:endchunk,ndlev,ntpy)
      soxyrhi            ! SOX input data for upper bound year
   real(r8), allocatable, dimension(:,:,:) :: &
      sox                ! SOX interpolated data (pcols,begchunk:endchunk,ndlev)

   integer ::&
      ntpy,             &! number of time samples per year in emissions dataset
      nyr,              &! number of years in emissions dataset
      ncid,             &! ID for netCDF file
      loyri,            &! index in yr array for lower bound year
      start(4),         &! start vector for netCDF hyperslabs
      count(4)           ! count vector for netCDF hyperslabs

!##############################################################################
contains
!##############################################################################


  subroutine soxbndini( inyear )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Open netCDF file containing SOX emissions data.  Initialize arrays
      ! with the data to be interpolated to the current time.
      !
      ! It is assumed that each year for which data is available contains
      ! the same number of time samples.
      !
      ! It is assumed that the input data contains a "date" coordinate that is
      ! increasing and represents dates in the format yyyymmdd.
      ! 
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: alloc_err, handle_ncerr, handle_err
      use phys_grid,      only: scatter_field_to_chunk
      use ioFileMod, only: getfil
      use filenames, only: bndtvsox

#if ( defined SPMD )
      use mpishorthand
#endif

#ifdef MATCH
      use calendar
#endif

!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      integer, intent(in) ::&
         inyear    ! current year

      ! Local variables:
      integer year  ! year that is used (may be inyear or rampyear_prognostic_sulfur if fixed)
      integer ::&
         i, istat, did, nlon, vid, recid, nrec
      real(r8), allocatable :: soxlatlonhi(:,:,:,:)  ! used for netCDF input
      real(r8), allocatable :: soxlatlonlo(:,:,:,:)  ! used for netCDF input
      character*10 errmsg                            ! error messages

      ! Externals
!      real(r8) ::&
!         caldayr
      !-----------------------------------------------------------------------
      ! Get file name.  
      call getfil(bndtvsox, soxems_data, 0)

      !++drb Set logical flags according to namelist variables
      !    Currently the only scenario available is RAMPED      
      if ( scenario_prognostic_sulfur == 'RAMPED' ) then
         doramp_sox = .true.
         if ( masterproc ) &
              write(6,*)'soxbndini: scenario_prognostic_sulfur = RAMPED'
      else if  ( scenario_prognostic_sulfur == 'FIXED' ) then
         call endrun('soxbndini: input namelist SCENARIO_SOX must be set to RAMPED')
      else
         call endrun('soxbndini: input namelist SCENARIO_SOX must be set to RAMPED')
      end if

      ! Determine if namelist variable to cycle year has been set
      if ( rampyear_prognostic_sulfur .ne. bigint ) then
         fixyear_sox = .true.
         if ( masterproc ) &
              write(6,*) 'soxbndini: cycling emissions from year ',rampyear_prognostic_sulfur
      else
         fixyear_sox = .false.
      end if

      ! Set year to get out of file to either year from model (inyear) or fixed (rampyear_prognostic_sulfur)
      if ( .not. fixyear_sox ) then
         year = inyear
      else
         if ( masterproc ) &
              write(6,*)'soxbndini: WARNING- overriding model year with rampyear_prognostic_sulfur',rampyear_prognostic_sulfur
         year = rampyear_prognostic_sulfur
      endif
      !--drb fixed year 

      if (masterproc) then

        write(6,*)'soxbndini: called for model year:',year

      ! Open file.
!nf90      call handle_ncerr( nf90_open( soxems_data, NF_NOWRITE, ncid ), &
!nf90         'soxbndini: error opening file '//trim(soxems_data) )
        call handle_ncerr( nf_open( soxems_data, NF_NOWRITE, ncid ), &
           'soxbndini: error opening file '//trim(soxems_data) )

      ! get the record id
!nf90      call handle_ncerr( nf90_inquire( ncid, unlimiteddimid=recid), &
!nf90         'soxbndini: no record variables ' )
        call handle_ncerr( nf_inq_unlimdim( ncid, recid), &
           'soxbndini: no record variables ' )

      ! Check that input data is a right resolution.
!nf90      call handle_ncerr( nf90_inq_dimid( ncid, 'lon', did ), 'soxbndini: ' )
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, did, len=nlon ), 'soxbndini: ' )
        call handle_ncerr( nf_inq_dimid( ncid, 'lon', did ), 'soxbndini: ' )
        call handle_ncerr( nf_inq_dimlen( ncid, did, nlon ), 'soxbndini: ' )
        if ( nlon .ne. plon ) then
           write(unit=errmsg,fmt='(i10)') nlon
           call endrun('soxbndini: number of longitudes ('//trim(errmsg)//') doesn''t match model resolution.')
        end if

      ! Get size of unlimited dimension.
!nf90      call handle_ncerr( nf90_inquire_dimension( ncid, recid, len=nrec ), 'soxbndini: ' )
        call handle_ncerr( nf_inq_dimlen( ncid, recid, nrec ), 'soxbndini: ' )

      endif  ! end of masterproc

#if ( defined SPMD )
      ! broadcast nrec to nodes
      call mpibcast (nrec, 1, mpiint, 0, mpicom)
#endif

      allocate( date(nrec), stat=istat )
      call alloc_err( istat, 'soxbndini', 'date', nrec )

      ! Get date coordinate.
      if (masterproc) then

!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'date', vid ), &
!nf90         'soxbndini: cannot find date coordinate variable' )
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, date ), &
!nf90         'soxbndini: error getting date coordinate data' )
        call handle_ncerr( nf_inq_varid( ncid, 'date', vid ), &
           'soxbndini: cannot find date coordinate variable' )
        call handle_ncerr( nf_get_var_int( ncid, vid, date ), &
           'soxbndini: error getting date coordinate data' )

      endif  ! end of masterproc

#if ( defined SPMD )
      ! broadcast date to nodes
      call mpibcast (date, nrec, mpiint, 0, mpicom)
#endif

      ! Determine number of time samples per year.
      ntpy = 1
      do i = 2, nrec
         if ( date(i)/10000 .eq. date(1)/10000 ) ntpy = ntpy + 1
      end do

      ! Construct the years array.
      nyr = nrec/ntpy
      allocate( yr(nyr), stat=istat )
      call alloc_err( istat, 'soxbndini', 'yr', nyr )
      do i = 1, nyr
         yr(i) = date((i-1)*ntpy+1)/10000
      end do
      if ( masterproc ) &
           write(6,*)'soxbndini: years in emission dataset:',(yr(i),i=1,nyr)

      ! Construct array of calendar days for the annual cycle.
      allocate( cdays(ntpy), stat=istat )
      call alloc_err( istat, 'soxbndini', 'cdays', ntpy )
      do i = 1, ntpy
#ifdef MATCH         
         cdays(i) = caldayr( date(i), 0 )      ! match version
#else
         call bnddyi( date(i), 0, cdays(i) ) ! ccm version
#endif
      end do
      if ( masterproc ) &
           write(6,*)'soxbndini: calendar days in annual cycle:', (cdays(i),i=1,ntpy)

      if ( masterproc ) &
           write (6,*) 'soxbndini: searching for emissions for year ', year

      ! Check that requested year is contained in the data.
      if ( year .lt. yr(1) .or. year .gt. yr(nyr) ) then
         call endrun ('SOXBNDINI: requested year outside data limits in '//trim(soxems_data))
      end if

      ! Find index for the data year that is the lower bound of the
      ! interval containing the input year.
      do i = 1, nyr
         if ( yr(i) .gt. year ) then
            loyri = i - 1
            exit
         end if
      end do

      allocate( soxyrlo(pcols,begchunk:endchunk,ndlev,ntpy), stat=istat )
      call alloc_err( istat, 'soxbndini', 'soxyrlo', pcols*((endchunk-begchunk)+1)*ndlev*ntpy )
      allocate( soxyrhi(pcols,begchunk:endchunk,ndlev,ntpy), stat=istat )
      call alloc_err( istat, 'soxbndini', 'soxyrhi', pcols*((endchunk-begchunk)+1)*ndlev*ntpy )
      allocate( sox(pcols,begchunk:endchunk,ndlev), stat=istat )
      call alloc_err( istat, 'soxbndini', 'sox', pcols*((endchunk-begchunk)+1)*ndlev )
      soxyrlo(:,:,:,:) = 0._r8
      soxyrhi(:,:,:,:) = 0._r8
      sox(:,:,:) = 0._r8

      allocate( soxlatlonlo(plon,plat,ndlev,ntpy), &
                soxlatlonhi(plon,plat,ndlev,ntpy), stat=istat )
      call alloc_err( istat, 'soxbndini', 'soxlatlonlo and soxlatlonhi', &
                      plon*plat*ndlev*ntpy*2 )

      if (masterproc) then

        ! Read SOx data for years surrounding initial year.
        start(1) = 1
        start(2) = 1
        start(3) = 1
        start(4) = 1
        count(1) = plon
        count(2) = plat
        count(3) = ndlev
        count(4) = ntpy

!nf90      call handle_ncerr( nf90_inq_varid( ncid, 'SOx', vid ), &
!nf90         'soxbndini: cannot find variable '//'SOx' )
        call handle_ncerr( nf_inq_varid( ncid, 'SOx', vid ), &
           'soxbndini: cannot find variable '//'SOx' )

        start(4) = (loyri-1)*ntpy + 1
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, soxlatlonlo, start, count ), &
!nf90         'soxbndini: cannot read data for '//'SOx' )
        call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, soxlatlonlo ), &
           'soxbndini: cannot read data for '//'SOx' )

        start(4) = start(4) + ntpy
!nf90      call handle_ncerr( nf90_get_var( ncid, vid, soxlatlonhi, start, count ), &
!nf90         'soxbndini: cannot read data for '//'SOx' )
        call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, soxlatlonhi ), &
           'soxbndini: cannot read data for '//'SOx' )

      endif  ! end of masterproc

      ! scatter data to chunked data structures
      call scatter_field_to_chunk (1, 1, ndlev*ntpy, plon, soxlatlonlo, &
                                   soxyrlo(1,begchunk,1,1))
      call scatter_field_to_chunk (1, 1, ndlev*ntpy, plon, soxlatlonhi, &
                                   soxyrhi(1,begchunk,1,1))

      if ( masterproc ) then
        write(6,*)'soxbndini: read data for years; ', yr(loyri),' and ',yr(loyri+1)
      endif  ! end of masterproc

      deallocate( soxlatlonlo, soxlatlonhi, stat=istat )
      call handle_err( istat, &
         'ERROR deallocating memory for soxlatlonlo and soxlatlonhi in routine soxbndini')

   end subroutine soxbndini

!#######################################################################

    subroutine soxbndint( inyear, calday )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Interpolate SOX data to the current time.  Update the input data
      ! as necessary.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      use error_messages, only: alloc_err, handle_ncerr, handle_err
      use phys_grid,      only: scatter_field_to_chunk

#if ( defined SPMD )
      use mpishorthand
#endif

!-----------------------------------------------------------------------
#include <netcdf.inc>
!-----------------------------------------------------------------------

      integer, intent(in) ::&
           inyear    ! current year
      real(r8), intent(in) ::&
         calday  ! calendar day (w/ fraction) in current year, range = [1.,366.)

      ! Local variables:
      integer year  ! year that is used (may be inyear or rampyear_prognostic_sulfur if fixed)
      integer ::&
         i, j, k, n,    &
         vid,           &
         lotim, hitim
      real(r8) :: dt, dt1, tint
      real(r8), allocatable :: anncyclo(:,:,:)   ! low year data interpolated to calday
      real(r8), allocatable :: anncychi(:,:,:)   ! high year data interpolated to calday
      ! TBH:  Get rid of all the replicated code for I/O and scattering!  
      real(r8), allocatable :: soxlatlonhi(:,:,:,:)  ! used for netCDF input
      integer :: istat
      !-----------------------------------------------------------------------

      allocate( anncyclo(pcols,begchunk:endchunk,ndlev), &
                anncychi(pcols,begchunk:endchunk,ndlev), stat=istat )
      call alloc_err( istat, 'soxbndint', 'anncyclo and anncychi', &
                      pcols*((endchunk-begchunk)+1)*ndlev*2 )
      anncyclo(:,:,:) = 0._r8
      anncychi(:,:,:) = 0._r8

      !++drb fixed year 
      if ( .not. fixyear_sox ) then
         year = inyear
      else
         year = rampyear_prognostic_sulfur
      endif
      !--drb fixed year 


      ! Check to see if model year is still bounded by dataset years.
      if ( year .gt. yr(nyr) ) then
         write(6,*)'soxbndint: requested year = ',year, ' last dataset year = ',yr(nyr)
         call endrun
      end if

      if ( year .gt. yr(loyri+1) ) then
         loyri = loyri + 1
         soxyrlo = soxyrhi
         allocate( soxlatlonhi(plon,plat,ndlev,ntpy), stat=istat )
         call alloc_err( istat, 'soxbndini', 'soxlatlonhi', &
                         plon*plat*ndlev*ntpy )
         ! Read in new soxyrhi data.  Replace old soxyrlo data.
         if (masterproc) then
           start(4) = start(4) + ntpy
!nf90         call handle_ncerr( nf90_inq_varid( ncid, 'SOx', vid ), &
!nf90            'soxbndint: cannot find variable '//'SOx' )
!nf90         call handle_ncerr( nf90_get_var( ncid, vid, soxlatlonhi, start, count ), &
!nf90            'soxbndint: cannot read data for '//'SOx' )
           call handle_ncerr( nf_inq_varid( ncid, 'SOx', vid ), &
              'soxbndint: cannot find variable '//'SOx' )
           call handle_ncerr( nf_get_vara_double( ncid, vid, start, count, soxlatlonhi ), &
              'soxbndint: cannot read data for '//'SOx' )
         endif  ! end of masterproc

         ! scatter data to chunked data structures
         call scatter_field_to_chunk (1, 1, ndlev*ntpy, plon, soxlatlonhi, &
                                      soxyrhi(1,begchunk,1,1))

         if ( masterproc ) then
           write(6,*)'soxbndint: read data for year; ',yr(loyri+1)
         endif  ! end of masterproc
         deallocate( soxlatlonhi, stat=istat )
         call handle_err( istat, &
             'ERROR deallocating memory for soxlatlonhi in routine soxbndint')

      end if

      ! Linear interpolation...  Start by computing the number of days between
      !                          the lower and upper bounds, and days between
      !                          the model time and lower bound.

      if ( ntpy .gt. 1 ) then

         call findplb( cdays, ntpy, calday, lotim )
         hitim = mod( lotim, ntpy ) + 1

         if( cdays(hitim) .lt. cdays(lotim) )then
            dt = 365. - cdays(lotim) + cdays(hitim)
            if( calday .le. cdays(hitim) )then
               dt1 = 365. - cdays(lotim) + calday
            else
               dt1 = calday - cdays(lotim)
            end if
         else
            dt = cdays(hitim) - cdays(lotim)
            dt1 = calday - cdays(lotim)
         end if
         tint = dt1/dt
         ! Annual cycle interpolations.
!TBH         call linintp( size(anncyclo), 0._r8, 1._r8, tint, soxyrlo(1,1,1,lotim), &
!TBH                       soxyrlo(1,1,1,hitim), anncyclo )
!TBH         call linintp( size(anncychi), 0._r8, 1._r8, tint, soxyrhi(1,1,1,lotim), &
!TBH                       soxyrhi(1,1,1,hitim), anncychi )
         call linintp( size(anncyclo), 0._r8, 1._r8, tint,  &
                       soxyrlo(1,begchunk,1,lotim), &
                       soxyrlo(1,begchunk,1,hitim), &
                       anncyclo )
         call linintp( size(anncychi), 0._r8, 1._r8, tint,  &
                       soxyrhi(1,begchunk,1,lotim), &
                       soxyrhi(1,begchunk,1,hitim), &
                       anncychi )
      else
         anncyclo(:,:,:) = soxyrlo(:,:,:,1)
         anncychi(:,:,:) = soxyrhi(:,:,:,1)
      end if

      ! Interpolate between years for which annual cycle data is present
      dt = yr(loyri+1) - yr(loyri)
      dt1 = year - yr(loyri)
      tint = dt1/dt
      call linintp( size(sox), 0._r8, 1._r8, tint, anncyclo, anncychi, sox )

      ! deallocate memory
      deallocate( anncyclo, anncychi, stat=istat )
      call handle_err( istat, &
        'ERROR deallocating memory for anncyclo and anncychi in routine soxbndint')

   end subroutine soxbndint

!#######################################################################

   subroutine soxbndget( ncol, lchnk, x )

      !----------------------------------------------------------------------- 
      ! Purpose: 
      ! Return SOX emission data for the requested latitude.
      !
      ! Author: B. Eaton
      !----------------------------------------------------------------------- 

      integer, intent(in) :: ncol          ! number of columns used
      integer, intent(in) :: lchnk         ! chunk index

      real(r8), intent(out) :: x(pcols,ndlev)  ! SOx emissions in Tg S/m2/s

      ! Local variables:
      integer :: i
      !-----------------------------------------------------------------------

      do i = 1, ncol
         x(i,1) = sox(i,lchnk,1)
         x(i,2) = sox(i,lchnk,2)
      end do

   end subroutine soxbndget

!#######################################################################

end module soxbnd
