#include <misc.h>

module ghg_surfvals

!-----------------------------------------------------------------------------------
! Purpose: Provides greenhouse gas (ghg) values at the Earth's surface.
!          These values may be time dependent.
!
! Author: Brian Eaton (assembled module from existing scattered code pieces)
!-----------------------------------------------------------------------------------

   use shr_kind_mod, only: r8=>shr_kind_r8
   use pmgrid,       only: masterproc
   use abortutils,   only: endrun

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make default access private
   save

! Public methods
   public ::&
      ghg_surfvals_init,        &! initialize options that depend on namelist input
      ghg_surfvals_set,         &! set ghg surface values when ramp option is on
      ghg_surfvals_ramp,        &! returns true when ramp option is on
      ghg_surfvals_get_co2mmr    ! computes and returns the co2 mass mixing ratio

! Public data
   real(r8), public :: co2vmr = 3.550e-4               ! co2   volume mixing ratio 
   real(r8), public :: n2ovmr = 0.311e-6               ! n2o   volume mixing ratio 
   real(r8), public :: ch4vmr = 1.714e-6               ! ch4   volume mixing ratio 
   real(r8), public :: f11vmr = 0.280e-9               ! cfc11 volume mixing ratio 
   real(r8), public :: f12vmr = 0.503e-9               ! cfc12 volume mixing ratio 
   character(len=16), public :: scenario_ghg = 'FIXED' ! 'FIXED','RAMPED' or 'RAMP_CO2_ONLY'
   integer, public  :: rampYear_ghg = 0                ! ramped gases fixed at this year (if > 0)
   character(len=256), public  :: bndtvghg             ! filename for ramped data
   integer, public  :: ramp_co2_start_ymd = 0          ! start date for co2 ramping (yyyymmdd)
   real(r8), public :: ramp_co2_annual_rate = 1.0      ! % amount of co2 ramping per yr; default is 1% 
   real(r8), public :: ramp_co2_cap = -9999.0          ! co2 ramp cap if rate>0, floor otherwise 
                                                       ! as multiple or fraction of inital value
                                                       ! ex. 4.0 => cap at 4x initial co2 setting 

! Private methods
   private ::&
      ghg_surfvals_set_all,     &! set all ghg surface values when ramp option is on
      ghg_surfvals_set_co2       ! set just co2 surface values when ramp option is on

! Private module data
   logical :: doRamp_ghg    ! true => turn on ramping for ghg
   logical :: ramp_just_co2 ! true => ramping to be done just for co2 and not other ghg's
   integer :: fixYear_ghg   ! year at which Ramped gases are fixed
   integer :: co2_start     ! date at which co2 begins ramping
   real(r8) :: co2_daily_factor    ! daily multiplier to achieve annual rate of co2 ramp
   real(r8) :: co2_limit    ! value of co2vmr where ramping ends
   real(r8) :: co2_base     ! initial co2 volume mixing ratio, before any ramping
   integer :: ntim = -1               ! number of yearly data values
   integer,  allocatable :: yrdata(:) ! yearly data values
   real(r8), allocatable :: co2(:)    ! co2 mixing ratios in ppmv 
   real(r8), allocatable :: ch4(:)    ! ppbv
   real(r8), allocatable :: n2o(:)    ! ppbv
   real(r8), allocatable :: f11(:)    ! pptv
   real(r8), allocatable :: f12(:)    ! pptv
   real(r8), allocatable :: adj(:)    ! unitless adjustment factor for f11 & f12

!=========================================================================================
contains
!=========================================================================================

subroutine ghg_surfvals_init()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize the ramp options that are controlled by namelist input.
! Set surface values at initial time.
! N.B. This routine must be called after the time manager has been initialized
!      since ghg_surfvals_set calls time manager methods.
! 
! Author: B. Eaton - merged code from parse_namelist and rampnl_ghg.
! 
!-----------------------------------------------------------------------
   use infnan,       only: inf

!---------------------------Local variables-----------------------------


   if (scenario_ghg == 'FIXED') then
      doRamp_ghg = .false.
      ramp_just_co2 = .false.
      if (masterproc) &
         write(6,*)'ghg_surfvals_init: ghg surface values are fixed as follows'

   else if (scenario_ghg == 'RAMPED') then
      doRamp_ghg = .true.
      ramp_just_co2 = .false.
      call ghg_ramp_read

      fixYear_ghg = rampYear_ghg     ! set private member to namelist var
      if (masterproc) then
         if ( fixYear_ghg > 0 ) then
            write(6,*) '  FIXED values from year ',fixYear_ghg
         else
            write(6,*) '  RAMPED values initialized to'
         end if
      end if
      call ghg_surfvals_set()

   else if (scenario_ghg == 'RAMP_CO2_ONLY') then
      if(ramp_co2_start_ymd == 0) then
         call endrun ('ghg_surfvals_init: RAMP_CO2_START_YMD must be set for SCENARIO_GHG=RAMP_CO2_ONLY')
      else
         co2_start = ramp_co2_start_ymd
      end if

      if(ramp_co2_annual_rate <= -100.0) then
         write(6,*) 'RAMP_CO2:  invalid ramp_co2_annual_rate= ',ramp_co2_annual_rate
         call endrun ('ghg_surfvals_init: RAMP_CO2_ANNUAL_RATE must be greater than -100.0')
      end if

      doRamp_ghg = .true.
      ramp_just_co2 = .true.
      co2_base = co2vmr        ! save initial setting 
      if (masterproc) &
           write(6,*) '  RAMPED values initialized to'

      co2_daily_factor = (ramp_co2_annual_rate*0.01_r8+1.0)**(1.0_r8/365.0_r8)

      if(ramp_co2_cap > 0.0) then  
         co2_limit = ramp_co2_cap * co2_base
      else                                  ! if no cap/floor specified, provide default
         if(ramp_co2_annual_rate < 0.0) then
            co2_limit = 0.0
         else
            co2_limit = inf
         end if
      end if
      if((ramp_co2_annual_rate<0.0 .and. co2_limit>co2_base) .or. &
         (ramp_co2_annual_rate>0.0 .and. co2_limit<co2_base)) then
         write(6,*) 'RAMP_CO2: ramp_co2_cap is unreachable'
         write(6,*) 'RAMP_CO2: ramp_co2_annual_rate= ',ramp_co2_annual_rate,' ramp_co2_cap= ',ramp_co2_cap
         call endrun('GHG_SURFVALS_INIT:  ramp_co2_annual_rate and ramp_co2_cap incompatible')
      end if

      call ghg_surfvals_set()
   else
      call endrun ('ghg_surfvals_init: input namelist SCENARIO_GHG must be set to either FIXED, RAMPED or RAMP_CO2_ONLY')
   endif

   if (masterproc) then
      write(6,*) '  co2 volume mixing ratio = ',co2vmr
      write(6,*) '  ch4 volume mixing ratio = ',ch4vmr
      write(6,*) '  n2o volume mixing ratio = ',n2ovmr
      write(6,*) '  f11 volume mixing ratio = ',f11vmr
      write(6,*) '  f12 volume mixing ratio = ',f12vmr
   end if

end subroutine ghg_surfvals_init

!=========================================================================================

subroutine ghg_ramp_read()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Read ramped greenhouse gas surface data.  
! 
! Author: T. Henderson
! 
!-----------------------------------------------------------------------

   use ioFileMod, only: getfil
#if ( defined SPMD )
   use mpishorthand, only: mpicom, mpiint, mpir8
#endif

   include 'netcdf.inc'

!---------------------------Local variables-----------------------------
   integer :: ncid
   integer :: co2_id
   integer :: ch4_id
   integer :: n2o_id
   integer :: f11_id
   integer :: f12_id
   integer :: adj_id
   integer :: date_id
   integer :: time_id
   integer :: ierror
   character(len=256) :: locfn          ! netcdf local filename to open

   if (masterproc) then
     call getfil (bndtvghg, locfn, 0)
     call wrap_open (trim(locfn), NF_NOWRITE, ncid)
     write(6,*)'GHG_RAMP_READ:  reading ramped greenhouse gas surface data from file ',trim(locfn)
     call wrap_inq_varid( ncid, 'date', date_id )
     call wrap_inq_varid( ncid, 'CO2', co2_id )
     call wrap_inq_varid( ncid, 'CH4', ch4_id )
     call wrap_inq_varid( ncid, 'N2O', n2o_id )
     call wrap_inq_varid( ncid, 'f11', f11_id )
     call wrap_inq_varid( ncid, 'f12', f12_id )
     call wrap_inq_varid( ncid, 'adj', adj_id )
     call wrap_inq_dimid( ncid, 'time', time_id )
     call wrap_inq_dimlen( ncid, time_id, ntim )
   endif
#if (defined SPMD )
   call mpibcast (ntim, 1, mpiint, 0, mpicom)
#endif
   ! these arrays are never deallocated
   allocate ( yrdata(ntim), co2(ntim), ch4(ntim), n2o(ntim),    &
                 f11(ntim), f12(ntim), adj(ntim), stat=ierror )
   if (ierror /= 0) then
     write(6,*)'GHG_RAMP_READ:  ERROR, allocate() failed!'
     call endrun
   endif
   if (masterproc) then
     call wrap_get_var_int   (ncid, date_id, yrdata )
     yrdata = yrdata / 10000
     call wrap_get_var_realx (ncid, co2_id, co2 )
     call wrap_get_var_realx (ncid, ch4_id, ch4 )
     call wrap_get_var_realx (ncid, n2o_id, n2o )
     call wrap_get_var_realx (ncid, f11_id, f11 )
     call wrap_get_var_realx (ncid, f12_id, f12 )
     call wrap_get_var_realx (ncid, adj_id, adj )
     call wrap_close (ncid)
     write(6,*)'GHG_RAMP_READ:  successfully read ramped greenhouse gas surface data from years ',yrdata(1),' through ',yrdata(ntim)
   endif
#if (defined SPMD )
   call mpibcast (co2, ntim, mpir8, 0, mpicom)
   call mpibcast (ch4, ntim, mpir8, 0, mpicom)
   call mpibcast (n2o, ntim, mpir8, 0, mpicom)
   call mpibcast (f11, ntim, mpir8, 0, mpicom)
   call mpibcast (f12, ntim, mpir8, 0, mpicom)
   call mpibcast (adj, ntim, mpir8, 0, mpicom)
   call mpibcast (yrdata, ntim, mpiint, 0, mpicom)
#endif

   return

end subroutine ghg_ramp_read

!=========================================================================================

function ghg_surfvals_ramp()
   logical :: ghg_surfvals_ramp
   ghg_surfvals_ramp = doRamp_ghg
end function ghg_surfvals_ramp

!=========================================================================================

function ghg_surfvals_get_co2mmr()
  use physconst,    only: mwdry, mwco2

  real(r8), parameter :: rmwco2 = mwco2/mwdry    ! ratio of molecular weights of co2 to dry air
  real(r8) :: ghg_surfvals_get_co2mmr            ! co2 mass mixing ratio

  ghg_surfvals_get_co2mmr = rmwco2 * co2vmr  

end function ghg_surfvals_get_co2mmr

!=========================================================================================


subroutine ghg_surfvals_set()

   use time_manager, only: get_curr_date, is_end_curr_day

!---------------------------Local variables-----------------------------

   integer  :: yr, mon, day, ncsec ! components of a date
   integer  :: ncdate              ! current date in integer format [yyyymmdd]

   if(ramp_just_co2) then
      call ghg_surfvals_set_co2()
   else
      call ghg_surfvals_set_all()
   end if

   if (masterproc .and. is_end_curr_day()) then
      call get_curr_date(yr, mon, day, ncsec)
      ncdate = yr*10000 + mon*100 + day
      write(6,*) 'ghg_surfvals_set: ncdate= ',ncdate,' co2vmr=',co2vmr
   end if

   return
end subroutine ghg_surfvals_set

!=========================================================================================

subroutine ghg_surfvals_set_all()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes greenhouse gas volume mixing ratios via interpolation of
! yearly input data.
! 
! Author: B. Eaton - updated ramp_ghg for use in ghg_surfvals module
! 
!-----------------------------------------------------------------------
   use time_manager, only: get_curr_date, get_curr_calday
   use timeinterp,   only: validfactors

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
   real(r8) cfcscl               ! cfc scale factor for f11

!
! ---------------------------------------------------------------------
!
   calday = get_curr_calday()
   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day
!
! determine index into input data
!
   if ( fixYear_ghg > 0) then
      yrmodel  = fixYear_ghg
   else
      yrmodel  = yr   
   end if

   nyrm       = yrmodel - yrdata(1) + 1
   nyrp       = nyrm + 1
!
! if current date is before yrdata(1), quit
!
   if (nyrm < 1) then
      write(6,*)'ghg_surfvals_set_all: data time index is out of bounds'
      write(6,*)'nyrm = ',nyrm,' nyrp= ',nyrp, ' ncdate= ', ncdate
      call endrun
   endif
!
! if current date later than yrdata(ntim), call endrun.
! if want to use ntim values - uncomment the following lines
! below and comment the call to endrun and previous write
!
   if (nyrp > ntim) then
      call endrun ('ghg_surfvals_set_all: error - current date is past the end of valid data')
!         write(6,*)'ghg_surfvals_set_all: using ghg data for ',yrdata(ntim)
!         co2vmr = co2(ntim)*1.e-06
!         ch4vmr = ch4(ntim)*1.e-09
!         n2ovmr = n2o(ntim)*1.e-09
!         f11vmr = f11(ntim)*1.e-12*(1.+cfcscl)
!         f12vmr = f12(ntim)*1.e-12
!         co2mmr = rmwco2 * co2vmr
!         return
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
      write(6,*)'ghg_surfvals_set_all: Bad fact1 and/or fact2=',fact1,fact2
   end if
!
! do time interpolation:
!   co2     in ppmv
!   n2o,ch4 in ppbv
!   f11,f12 in pptv
!
   co2vmr = (co2(nyrm)*fact1 + co2(nyrp)*fact2)*1.e-06
   ch4vmr = (ch4(nyrm)*fact1 + ch4(nyrp)*fact2)*1.e-09
   n2ovmr = (n2o(nyrm)*fact1 + n2o(nyrp)*fact2)*1.e-09

   cfcscl = (adj(nyrm)*fact1 + adj(nyrp)*fact2)
   f11vmr = (f11(nyrm)*fact1 + f11(nyrp)*fact2)*1.e-12*(1.+cfcscl)
   f12vmr = (f12(nyrm)*fact1 + f12(nyrp)*fact2)*1.e-12

   return
end subroutine ghg_surfvals_set_all

!=========================================================================================

subroutine ghg_surfvals_set_co2()
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Computes co2 greenhouse gas volume mixing ratio via ramping info 
! provided in namelist var's
! 
! Author: B. Eaton - updated ramp_ghg for use in ghg_surfvals module
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use time_manager, only: get_curr_date, timemgr_datediff

!---------------------------Local variables-----------------------------

   real(r8) :: daydiff             ! number of days of co2 ramping
   integer  :: yr, mon, day, ncsec ! components of a date
   integer  :: ncdate              ! current date in integer format [yyyymmdd]
!-----------------------------------------------------------------------

   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day

   call timemgr_datediff(co2_start, 0, ncdate, ncsec, daydiff)

   if (daydiff > 0.0) then

      co2vmr = co2_base*(co2_daily_factor)**daydiff

      if(co2_daily_factor < 1.0) then
         co2vmr = max(co2vmr,co2_limit)
      else
         co2vmr = min(co2vmr,co2_limit)
      end if
   end if

   return
end subroutine ghg_surfvals_set_co2


!=========================================================================================

end module ghg_surfvals
