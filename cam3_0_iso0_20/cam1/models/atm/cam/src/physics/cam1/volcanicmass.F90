#include <misc.h>
#include <params.h>

module volcanicmass
!----------------------------------------------------------------------- 
! 
! Purpose: read, store, interpolate, and retrieve mass fields describing
!   volcanic emissions
! 
! Author: A. Conley, but mostly copied from ozone code base.
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plat, plon, plond, masterproc, plev, plevp
  use commap,       only: latdeg
  use ppgrid,       only: begchunk, endchunk, pcols, pver, pverp
  use phys_grid,    only: scatter_field_to_chunk, get_ncols_p
  use filenames,    only: bndtvvolc
  use abortutils,   only: endrun
  use timeinterp,   only: validfactors

#if ( defined SPMD )
  use mpishorthand
#endif

  implicit none

  private
  save

  public read_volcanic_mass
  public get_volcanic_mass
  public volcanic_initialize


  real(r8) cdayprev  ! dataset calendar day previous month
  real(r8) cdaynext  ! dataset calendar day next month

  integer nprev     ! Array indices for previous month volcanic data
  integer nnext     ! Array indices for next month volcanic data
  integer ntmp      ! temp storage for swaping nprev, nnext
  integer massid    ! netcdf id for volcanic variable
  integer levsiz    ! size of level dimension on volcanic dataset
  integer latsiz    ! size of latitude dimension on volcanic dataset
  integer timesiz   ! size of time dimension on volcanic dataset
  integer np1       ! current forward time index of volcanic dataset
  integer ncid_volc ! netcdf file id for volcanic file

  type volcanic_pointers 
    real(r8), dimension(:,:,:), pointer :: mass  ! (pcols,levsiz,begchunk:endchunk)
  end type volcanic_pointers
  type (volcanic_pointers) :: volcanic(2)      ! pointers to volcanic masses in 2 bounding months
  type (volcanic_pointers) :: volcanic2(2)     ! pointers to volcanic masses in 2 bounding months on CAM grid
  real(r8), allocatable :: volcanic_mass(:,:,:)! (pcols,plev,begchunk:endchunk)
  real(r8), allocatable :: pin(:)              ! volcanic pressure level (levsiz)
  real(r8), allocatable :: volclat(:)          ! Latitudes of volcanic dataset (latsiz)

  integer, allocatable :: date_volc(:)         ! Date on volcanic dataset (YYYYMMDD)
                                               ! (timesiz)
  integer, allocatable :: sec_volc(:)          ! seconds of date (0-86399)
                                               ! (timesiz)

contains

subroutine volcanic_initialize
!----------------------------------------------------------------------- 
! 
! Purpose: Read appropriate portion of time-variant volcanic boundary 
!          dataset, containing volcanic mixing ratios as a 
!          function of latitude and pressure.  Read two
!          consecutive months between which the current date lies.  
!
!          Another routine
!          then evaluates the two path length integrals (with and without
!          pressure weighting) from zero to the interfaces between the input
!          levels.  It also stores the contribution to the integral from each
!          layer.
! 
! Method: Call appropriate netcdf wrapper routines and interpolate to model grid
! 
! Author: CCM Core Group
! Modified: P. Worley, August 2003, for chunking and performance optimization
! Modified: A. Conley, Sept 2003, for volcanic aerosols
! 
!-----------------------------------------------------------------------
  use commap, only: latdeg, londeg
  use rgrid, only: nlon
  use time_manager, only: get_curr_date, get_perp_date, get_curr_calday, &
                           is_perpetual
  use ioFileMod, only: getfil

!-----------------------------------------------------------------------
  include 'netcdf.inc'
!-----------------------------------------------------------------------
!
! Local workspace
!
  integer dateid                          ! netcdf id for date variable
  integer secid                           ! netcdf id for seconds variable
  integer latdimid                        ! netcdf id for latitude dimension
  integer levdimid                        ! netcdf id for level dimension
  integer bindimid                        ! netcdf id for bin dimension
  integer latid                           ! netcdf id for latitude variable
  integer levid                           ! netcdf id for level variable
  integer dimids(nf_max_var_dims)         ! variable shape
  integer dim_cnt(4)                      ! array of counts for each dimension
  integer dim_strt(4)                     ! array of starting indices
  integer ilon, ilev, ilat, time          ! longitude, level, latitude, time indices
  integer ichunk                          ! index into chunk
  integer  :: yr, mon, day                ! components of a date
  integer  :: ncdate                      ! current date in integer format [yyyymmdd]
  integer  :: ncsec                       ! current time of day [seconds]
  real(r8) :: calday                      ! current calendar day
  real(r8) caldayloc                      ! calendar day (includes yr)
  real(r8), allocatable :: volc2D(:,:)    ! temporary volcanic arrays
  real(r8), allocatable :: volc3D(:,:,:)
  real(r8), allocatable :: volcbdyprev(:,:,:,:) ! volcanic data previous time sample
  real(r8), allocatable :: volcbdynext(:,:,:,:) ! volcanic data next time sample

  character(len=256) :: locfn          ! netcdf local filename to open


!
! find and open file; abort if fail (getfil(,,0)).
!

  nprev = 1
  nnext = 2
!
! SPMD: Master does all the work.  Sends needed info to slaves
!
  if (masterproc) then
    calday = get_curr_calday()
    call get_curr_date(yr, mon, day, ncsec)
    ncdate = yr*10000 + mon*100 + day
    caldayloc = calday + yr*365.
!      call bnddyi(ncdate, 0, caldayloc) !!!!

!
! Get and check dimension info
!
    call getfil(bndtvvolc, locfn, 0)
    call wrap_open(locfn, 0, ncid_volc)

    CALL WRAP_INQ_DIMID( ncid_volc, 'lev', levdimid   )
    CALL WRAP_INQ_DIMID( ncid_volc, 'lat', latdimid   )
    CALL WRAP_INQ_DIMID( ncid_volc, 'bin', bindimid )
    CALL WRAP_INQ_DIMID( ncid_volc, 'date', dateid  )

    CALL WRAP_INQ_DIMLEN( ncid_volc, levdimid, levsiz   )
    CALL WRAP_INQ_DIMLEN( ncid_volc, latdimid, latsiz   )
    CALL WRAP_INQ_DIMLEN( ncid_volc, dateid, timesiz   )

!    CALL WRAP_INQ_VARID( ncid_volc, 'date', dateid   )
!    CALL WRAP_INQ_VARID( ncid_volc, 'datesec', secid   )
    CALL WRAP_INQ_VARID( ncid_volc, 'MVOLC', massid   )
    CALL WRAP_INQ_VARID( ncid_volc, 'lat', latid   )
    CALL WRAP_INQ_VARID( ncid_volc, 'lev', levid   )

    CALL WRAP_INQ_VARDIMID (ncid_volc, massid, dimids)
    if (dimids(1) /= levdimid .or. dimids(2) /= latdimid .or. dimids(3) /= bindimid .or. dimids(4) /= dateid ) then
      call endrun ('VOLCANIC_INITIALIZE: Data must be ordered lev, lat, time')
    endif
  endif

#if (defined SPMD )
  call mpibcast( latsiz, 1, mpiint, 0, mpicom )
  call mpibcast( levsiz, 1, mpiint, 0, mpicom )
  call mpibcast( timesiz, 1, mpiint, 0, mpicom )
#endif

!
! Dynamically allocated memory for module
!
  allocate (date_volc(timesiz))
  allocate (sec_volc(timesiz))
  allocate (pin(levsiz))
  allocate (volcanic(nprev)%mass(pcols,levsiz,begchunk:endchunk))
  allocate (volcanic(nnext)%mass(pcols,levsiz,begchunk:endchunk))
  allocate (volcanic_mass(pcols,levsiz,begchunk:endchunk))
!
! Locally dynamic that will be deallocated before "return"
!
  allocate (volc3D(plond,levsiz,plat))

  if (masterproc) then
!
! More dynamically allocated memory for module 
! (for masterproc only)
!
    allocate (volclat(latsiz))
!
! More locally dynamic that will be deallocated before return
!
    allocate (volcbdyprev(levsiz,latsiz,1,1))
    allocate (volcbdynext(levsiz,latsiz,1,1))
    allocate (volc2D(levsiz,plat))
!
! Retrieve longitude, latitude and level arrays for interpolation.
!
    call wrap_get_var_realx (ncid_volc, latid, volclat)
    call wrap_get_var_realx (ncid_volc, levid, pin)

!
! Retrieve entire date and sec variables.
!
    call wrap_get_var_int (ncid_volc,dateid,date_volc)
!    call wrap_get_var_int (ncid_volc,secid,sec_volc)
    sec_volc = 0


    dim_strt(1) = 1
    dim_strt(2) = 1
    dim_strt(3) = 1
    dim_cnt(1)  = levsiz
    dim_cnt(2)  = latsiz
    dim_cnt(3)  = 1  ! only 1 bin
    dim_cnt(4)  = 1
!
! Normal interpolation between consecutive time slices.
!
    do time=1,timesiz-1
      np1 = time + 1
      call bnddyi(date_volc(time), sec_volc(time), cdayprev)
      call bnddyi(date_volc(np1 ), sec_volc(np1 ), cdaynext)
      yr = date_volc(time)/10000
      cdayprev = cdayprev + yr*365.
      yr = date_volc(np1)/10000
      cdaynext = cdaynext + yr*365.
      if (caldayloc > cdayprev .and. caldayloc <= cdaynext) then
        dim_strt(4) = time
        call wrap_get_vara_realx (ncid_volc,massid,dim_strt,dim_cnt,volcbdyprev)
        dim_strt(4) = np1
        call wrap_get_vara_realx (ncid_volc,massid,dim_strt,dim_cnt,volcbdynext)
        goto 10
      endif
    end do
    write(6,*)'volcanic_initialize: Failed to find dates bracketing ncdate, ncsec=', ncdate, ncsec
    call endrun
10  continue
    write(6,*)'volcanic_initialize: Read volcanic data for dates ',date_volc(time), &
                sec_volc(time),' and ',date_volc(np1),sec_volc(np1)
!
! Latitude interpolation for "prev" data. Expand to 3D after interpolation.
!
    call lininterp (volcbdyprev(:,:,1,1),volclat   ,levsiz  ,latsiz  ,volc2D(:,:), &
                    latdeg  ,plat    )
    deallocate (volcbdyprev)

    do ilat=1,plat
      do ilev=1,levsiz
        do ilon=1,nlon(ilat)
          volc3D(ilon,ilev,ilat) = volc2D(ilev,ilat)
        end do
      end do
    end do
  
    call scatter_field_to_chunk(1,levsiz,1,plond,volc3D,volcanic(nprev)%mass)

#if (defined SPMD )
  else

    call scatter_field_to_chunk(1,levsiz,1,plond,volc3D,volcanic(nprev)%mass)

  endif
  if (masterproc) then
#endif

!
! Latitude interpolation for "next" data. Expand to 3D after interpolation.
!
    call lininterp (volcbdynext(:,:,1,1),volclat   ,levsiz  ,latsiz  ,volc2D(:,:), &
                    latdeg  ,plat    )
    deallocate (volcbdynext)

    do ilat=1,plat
      do ilev=1,levsiz
        do ilon=1,nlon(ilat)
          volc3D(ilon,ilev,ilat) = volc2D(ilev,ilat)
        end do
      end do
    end do
    deallocate (volc2D)

    call scatter_field_to_chunk(1,levsiz,1,plond,volc3D,volcanic(nnext)%mass)

#if (defined SPMD )
  else

    call scatter_field_to_chunk(1,levsiz,1,plond,volc3D,volcanic(nnext)%mass)
#endif
  endif

!
! Deallocate dynamic memory for local workspace.  NOT for pointers in common.
!
  deallocate (volc3D)

#if (defined SPMD )
  call mpibcast (np1, 1, mpiint, 0, mpicom)
  call mpibcast (pin, levsiz, mpir8, 0, mpicom)
  call mpibcast (date_volc, timesiz, mpiint, 0, mpicom )
  call mpibcast (sec_volc, timesiz, mpiint, 0, mpicom )
  call mpibcast (cdayprev, 1, mpir8, 0, mpicom)
  call mpibcast (cdaynext, 1, mpir8, 0, mpicom)
#endif

!
!  accumulate masses from surface before interpolating
!  allocate memory for boundary datasets on cam vertical grid
!  interpolate from file levels to cam levels
!  deallocate memory for boundary dataset (previous) on file levels
!

  do ichunk = begchunk, endchunk
    do ilev = levsiz-1,1,-1
      volcanic(nprev)%mass(:,ilev,ichunk)=volcanic(nprev)%mass(:,ilev,ichunk) &
          + volcanic(nprev)%mass(:,ilev+1,ichunk)
      volcanic(nnext)%mass(:,ilev,ichunk)=volcanic(nnext)%mass(:,ilev,ichunk) &
          + volcanic(nnext)%mass(:,ilev+1,ichunk)
    enddo
  enddo

  allocate (volcanic2(nprev)%mass(pcols,plevp,begchunk:endchunk))
  allocate (volcanic2(nnext)%mass(pcols,plevp,begchunk:endchunk))

  call vert_interpolate(nprev)
  call vert_interpolate(nnext)

!  deallocate (volcanic(nprev)%mass)

  return
end subroutine volcanic_initialize


subroutine read_volcanic_mass()
  
  use rgrid,        only: nlon
  use time_manager, only: get_curr_date, get_curr_calday

  integer :: dim_cnt(4)                ! array of counts for each dimension
  integer :: dim_strt(4)               ! array of starting indices
  integer :: ilon,ilev,ilat,icol,ncols ! indices
  integer :: lchnk                     ! chunk index
  integer :: yr,mon,day                ! components of date
  integer :: ncdate                    ! formated date [yyyymmdd]
  integer :: ncsec                     ! seconds past midnight of current ncdate
  real(r8):: calday                    ! current calendar day
  real(r8):: caldayloc                 ! calendar day (including yr)
  real(r8):: deltat                    ! time (days) between times of dataset

  real(r8), allocatable :: volc2D(:,:)        ! temporary volcanic arrays
  real(r8), allocatable :: volc3D(:,:,:)      ! temporary volcanic arrays
  real(r8), allocatable :: volcbdynext(:,:,:,:) ! volcanic data next time sample
  
  include 'netcdf.inc'


  call get_curr_date(yr, mon, day, ncsec)
  ncdate = yr*10000 + mon*100 + day

  calday = get_curr_calday()
  caldayloc = calday + yr*365.

  if (caldayloc > cdaynext) then
    write(6,*) 'caldayloc',caldayloc,'cdaynext',cdaynext
    np1 = np1 + 1
    if (np1 > timesiz) then
      call endrun ('READ_VOLCANIC_MASS: Attempt to read past end of Volcanic dataset')
    end if
    cdayprev = cdaynext
    call bnddyi(date_volc(np1), sec_volc(np1), cdaynext)
    yr = date_volc(np1)/10000
    cdaynext = cdaynext + yr*365.
    if(caldayloc <= cdaynext) then
      ntmp = nprev
      nprev = nnext
      nnext = ntmp

!
! Allocate memory for dynamic local workspace
!
      allocate (volc3D(plond,levsiz,plat))
      if (masterproc) then
        dim_strt(1) = 1
        dim_strt(2) = 1
        dim_strt(3) = 1
        dim_strt(4) = np1
        dim_cnt(1)  = levsiz
        dim_cnt(2)  = latsiz
        dim_cnt(3)  = 1  ! only 1 bin
        dim_cnt(4)  = 1
!
! Allocate memory for more dynamic local workspace
!
        allocate (volcbdynext(levsiz,latsiz,1,1))
        allocate (volc2D(levsiz,plat))

        call wrap_get_vara_realx(ncid_volc,massid,dim_strt,dim_cnt,volcbdynext)
        write(6,*)'get_volcanic_mass: read volcanic masses for date (yyyymmdd) ', date_volc(np1), ' sec ',sec_volc(np1)

!
! Spatial interpolation. 
!
        call lininterp (volcbdynext(:,:,1,1),volclat   ,levsiz  ,latsiz  ,volc2D(:,:), &
                        latdeg  ,plat    )

        do ilat=1,plat
          do ilev=1,levsiz
            do ilon=1,nlon(ilat)
              volc3D(ilon,ilev,ilat) = volc2D(ilev,ilat)
            end do
          end do
        end do
        call scatter_field_to_chunk(1,levsiz,1,plond,volc3D,volcanic(nnext)%mass)

        deallocate (volcbdynext)
        deallocate (volc2D)
        deallocate (volc3D)

      else  ! masterproc

        call scatter_field_to_chunk(1,levsiz,1,plond,volc3D,volcanic(nnext)%mass)
        deallocate (volc3D)
        
      endif ! masterproc

      do lchnk = begchunk,endchunk
        do ilev = levsiz-1,1,-1
          volcanic(nnext)%mass(:,ilev,lchnk)=volcanic(nnext)%mass(:,ilev,lchnk) &
              + volcanic(nnext)%mass(:,ilev+1,lchnk)
        enddo
      enddo
      call vert_interpolate(nnext)

    else  ! if(caldayloc <= cdaynext) 
      if (masterproc) then
        write(6,*)'READ_VOLCANIC_MASS: date from dataset ',date_volc(np1),' does not exceed model date ', ncdate
      endif
      call endrun
    endif ! if(caldayloc <= cdaynext) 
  endif ! (caldayloc > cdaynext) 

end subroutine read_volcanic_mass

subroutine get_volcanic_mass(lchnk, aerosol)
  
  use rgrid,        only: nlon
  use time_manager, only: get_curr_date, get_curr_calday

  integer, intent(in)   :: lchnk           ! chunk id for local aerosol array
  real(r8), intent(inout) :: aerosol(:,:)  ! aerosol array to be filled
  integer :: dim_cnt(4)                ! array of counts for each dimension
  integer :: dim_strt(4)               ! array of starting indices
  integer :: ilon,ilev,ilat,icol,ncols ! indices
  integer :: yr,mon,day                ! components of date
  integer :: ncdate                    ! formated date [yyyymmdd]
  integer :: ncsec                     ! seconds past midnight of current ncdate
  real(r8):: fact1, fact2              ! time interpolation factors
  real(r8):: calday                    ! current calendar day
  real(r8):: caldayloc                 ! calendar day (including yr)
  real(r8):: deltat                    ! time (days) between times of dataset

  real(r8), allocatable :: volc2D(:,:)        ! temporary volcanic arrays
  real(r8), allocatable :: volc3D(:,:,:)      ! temporary volcanic arrays
  real(r8), allocatable :: volcbdynext(:,:,:,:) ! volcanic data next time sample
  
!
! time interpolate
!
  call get_curr_date(yr, mon, day, ncsec)
  calday = get_curr_calday()
  caldayloc = calday + yr*365.

  deltat = cdaynext - cdayprev
  fact1 = (cdaynext - caldayloc)/deltat
  fact2 = (caldayloc - cdayprev)/deltat
!
! Check sanity of time interpolation calculation to within 32-bit roundoff
!
  if (.not. validfactors (fact1, fact2)) then
     write(6,*)'get_volcanic_mass: Bad fact1 and/or fact2=',fact1,fact2
     call endrun ()
  end if   
!
! Interpolate (in time)!
!
  ncols = get_ncols_p(lchnk)
  do ilev=1,plev
    do icol=1,ncols
      aerosol(icol,ilev) = volcanic2(nprev)%mass(icol,ilev,lchnk)*fact1 + volcanic2(nnext)%mass(icol,ilev,lchnk)*fact2
    enddo
  enddo
return

end subroutine get_volcanic_mass

subroutine vert_interpolate(month)
!-----------------------------------------------------------------------------
!
! volcanic2(month) is assigned vertically interpolated values of
! volcanic(month).  volcanic is on vertical grid defined by file.
! volcanic2 is on vertical grid of CAM simulation.
!
!-----------------------------------------------------------------------------

  use pmgrid,       only: plat, plon, plond, masterproc, plev, plevp
#include <comhyb.h>

  integer, intent(in):: month   ! which month boundary should be interpolated?
  real(r8) :: model_level(pverp) 
  real(r8) :: AER_diff(pverp) 
  integer :: ncols              ! number of columns in chunk
  integer :: lchnk              ! chunk index
  integer :: kk                 ! boundary data index into file data
  integer :: ilat, ilon, ilev, icol ! indices
  real(r8) :: dpl, dpu          ! delta p across lower and upper parts of interpolant

  model_level(:) = 1000._r8 *(hyai(:) + hybi(:))
  do lchnk = begchunk, endchunk
    ncols = get_ncols_p(lchnk)
    volcanic2(month)%mass(1:ncols,1,lchnk) = volcanic(month)%mass(1:ncols,1,lchnk)
!
! At every pressure level, interpolate onto that pressure level
!
    kk = 2 
    do ilev = 2, pver
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top pressure level for at least some
! of the longitude points.
!
!      if(masterproc) then
!        write(6,*) 'model_level:',model_level
!        write(6,*) 'pin:',pin
!      endif

      if (model_level(ilev) .lt. pin(1)) then
        volcanic2(month)%mass(1:ncols,ilev,lchnk) = volcanic(month)%mass(1:ncols,1,lchnk)

      else if (model_level(ilev) .ge. pin(levsiz)) then
        volcanic2(month)%mass(1:ncols,ilev,lchnk) = 0._r8

      else
        do while (model_level(ilev) .ge. pin(kk) ) 
           kk = kk + 1
        enddo
        dpu = model_level(ilev) - pin(kk-1)
        dpl =           pin(kk) - model_level(ilev)
        volcanic2(month)%mass(1:ncols,ilev,lchnk)    = &
              (volcanic(month)%mass(1:ncols,kk-1,lchnk)*dpl + &
               volcanic(month)%mass(1:ncols,kk  ,lchnk)*dpu)/(dpl + dpu)
      end if

    end do
!    do ilev = 1, pver
!      do icol = 1, ncols
!        if(volcanic2(month)%mass(icol,ilev,lchnk) .lt. 0) then
!          write(6,*)'vert interp error'
!          write(6,*)'lchnk,ilev,icol',lchnk,ilev,icol
!          write(6,*)'volcanic', volcanic(month)%mass(icol,:,lchnk)
!          write(6,*)'volcanic2',volcanic2(month)%mass(icol,:,lchnk)
!          call endrun
!        endif
!      enddo
!    enddo

!
! aerosol mass beneath lowest interface (pverp) must be 0
!
    volcanic2(month)%mass(1:ncols, pverp, lchnk)=0._r8

!
! Set mass in layer to zero whenever it is less than
!   1.e-40 kg/m^2 in the layer
!
!    do ilev = 1, pver
!      do icol = 1, ncols
!        if(volcanic2(month)%mass(icol,ilev,lchnk) < 1.e-40) volcanic2(month)%mass(icol,ilev,lchnk) = 0._r8
!      enddo
!    enddo

!!
!! Extract mass per layer from cumulative mass.
!!
    do ilev = 1, pver
      do icol = 1, ncols
        volcanic2(month)%mass(icol,ilev,lchnk) = volcanic2(month)%mass(icol,ilev,lchnk)-volcanic2(month)%mass(icol,ilev+1,lchnk)
        if(volcanic2(month)%mass(icol,ilev,lchnk) .lt. 0._r8) then
          volcanic2(month)%mass(icol,ilev,lchnk) = 0._r8
        endif
      enddo
    enddo

    do ilev = 1, pver
      do icol = 1, ncols
        if(volcanic2(month)%mass(icol,ilev,lchnk) .lt. 0) then
          write(6,*)'mass per layer error'
          write(6,*)'lchnk,ilev,icol',lchnk,ilev,icol
          write(6,*)'volcanic', volcanic(month)%mass(icol,:,lchnk)
          write(6,*)'volcanic2',volcanic2(month)%mass(icol,:,lchnk)
          call endrun ('VERT_INTERPOLATE')
        endif
      enddo
    enddo

  enddo ! lchnk

  return

end subroutine vert_interpolate


end module volcanicmass
