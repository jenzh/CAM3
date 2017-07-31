#include <misc.h>
#include <params.h>

module tracers_suite

!---------------------------------------------------------------
!
! Implements artificial suite of tracers
!    1) low tracer with unit mixing ratio at about 800 mb
!    2) med tracer with unit mixing ratio at about 500 mb
!    3) high tracer with unit mixing ratio at about 200 mb
!    4) reverse med tracer with unit mixing ratio everywhere except about 500 mb
!    5) unit tracer with unit mixing ratio everywhere
!
!  D Bundy June 2003
!  modified Feb 2004 to include TT_UN and smoothing
!
! This is a module that contains a suite of tracers that is interfaced
! by tracers.F90. The details of the suite are contained entirely
! in this file, so the public routines are all very generic. 
!
!  ------------  calling tree --------------
!  Initialize the tracer mixing ratio field
!  	-> tracers.F90: tracers_init_cnst 
!  		-> tracers_suite.F90:init_cnst_tr
!  			-> init_cnst_lw
!  			-> init_cnst_md
!  			-> init_cnst_hi
!  			-> init_cnst_un
!
!  Initialize data set, things that need to be done at the beginning of a 
!  run (whether restart or initial)
!  	-> tracers.F90: tracers_init
!  		-> tracers_suite.F90:init_tr
!
!  Timestepping:
!  	-> tracers_timestep_init
!  		-> tracers_suite.F90:timestep_init_tr
!
!  	-> tracers_timestep_tend
!  		-> tracers_suite.F90:flux_tr
!  			-> flux_lw  
!  			-> flux_md
!  			-> flux_hi
!  			-> flux_un
!
!  		-> tracers_suite.F90:tend_tr
!  			-> tend_lw
!  			-> tend_md
!  			-> tend_hi
!  			-> tend_un
!
!---------------------------------------------------------------


  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid !,       only: plon, plev, plat
  use ppgrid !,       only: pcols, pver
  use abortutils, only: endrun


  implicit none
  private
!  integer setpindxtr

  save

! Public interfaces
  public get_tracer_name  ! store names of tracers
  public init_cnst_tr      ! initialize tracer fields
  public init_tr      ! initialize data sets need for tracer
  public tend_tr      ! tracer tendency
  public flux_tr      ! surface flux of tracer
  public timestep_init_tr  ! interpolate tracer emissions data set

! Public module data
  integer, public, parameter :: trac_ncnst=5

! Private module data
  logical, parameter :: smooth = .false.
  
contains

!======================================================================
function get_tracer_name(n)

!------------------------------------------------------------------------
! Purpose:
!
! The tracer names are only defined in this module. This function is for
! outside programs to grab the name for each tracer number. 
!    If n > trac_ncst calls endrun.
!
!------------------------------------------------------------------------

! -----------------------------Arguments---------------------------------

  integer, intent(in) :: n
  character(len=8) :: get_tracer_name  ! return value

! ----------------------------- Local ---------------------------------
  character(len=8), dimension(trac_ncnst), parameter :: & ! constituent names
       tracer_names  =  (/ 'TT_LW', 'TT_MD', 'TT_HI', 'TTRMD' , 'TT_UN'/)
  
!-----------------------------------------------------------------------

  if ( n > trac_ncnst ) then
     write(6,*) 'tracers_suite:get_tracer_name()','requested tracer',n
     write(6,*) 'only ',trac_ncnst,' tracers available'
     call endrun
  else
     get_tracer_name = tracer_names(n)
  endif

  return

end function get_tracer_name


!======================================================================
!======================================================================
subroutine init_cnst_tr(m,q)

!----------------------------------------------------------------------- 
!
! Purpose:
! calls initialization routine for tracer m, returns mixing ratio in q
!
! This routine must be consistent with trac_ncnst.
!
!----------------------------------------------------------------------- 

   implicit none

   real(r8), intent(out) :: q(plon,plev,plat)    ! kg tracer/kg dry air
   integer ,intent(in) :: m                      ! index of tracer

   if ( m > trac_ncnst ) then
      write(6,*) 'tracers_suite:init_cnst_tr()'
      write(6,*) ' asked to initialize tracer number ',m
      write(6,*) ' but there are only trac_ncnst = ',trac_ncnst,' tracers'
      call endrun
   endif

   if ( m == 1 ) then
      call init_cnst_lw(q)
   else if ( m == 2 ) then
      call init_cnst_md(q)
   else if ( m == 3 ) then
      call init_cnst_hi(q)
   else if ( m == 4 ) then
      call init_cnst_md(q,rev_in=1)
   else if ( m == 5 ) then
      call init_cnst_un(q)
   else
      write(6,*) 'tracers_suite:init_cnst_tr()'
      write(6,*) 'no initialization routine specified for tracer',m
      call endrun
   endif
      
end subroutine init_cnst_tr



!======================================================================
subroutine init_cnst_lw(q)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test tracer 1: low
! 
!-----------------------------------------------------------------------
   implicit none

!Arguments
   real(r8), intent(out) :: q(plon,plev,plat)    ! kg tracer/kg dry air

! Local
  integer indx

!-----------------------------------------------------------------------
!
! Initialize low tracer to zero except at 800 level
!

   indx = setpindxtr(800.)

  if ( smooth ) then 
     call setsmoothtr(indx,q,.876)
  else 
     q = 0.0
     q(:,indx,:) = 1.0
  endif

!   print *,'suite 1 init_cnst_lw min/max q',minval(q),maxval(q)

end subroutine init_cnst_lw

!======================================================================
subroutine init_cnst_md(q,rev_in)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test tracer 2: med
! 
!-----------------------------------------------------------------------
   implicit none

!Arguments
   real(r8), intent(out) :: q(plon,plev,plat)    ! kg tracer/kg dry air
   integer,  intent(in), optional :: rev_in         ! reverse the mixing ratio

! Local
  integer indx
  integer rev

!-----------------------------------------------------------------------
!
! Initialize med tracer to zero except at 500 level
!

  rev = 0
  if (present(rev_in)) then
     if (rev_in == 1) then
        rev = 1
     endif
  endif

  indx = setpindxtr(500.)

  if ( smooth ) then 
     call setsmoothtr(indx,q,.876,rev_in=rev)
  else
     if (rev == 1 ) then
        q = 1.0
        q(:,indx,:) = 0.0
     else
        q = 0.0
        q(:,indx,:) = 1.0
     endif
  endif
   
!   print *,'suite 1 init_cnst_md min/max q',minval(q),maxval(q)

end subroutine init_cnst_md

!======================================================================
subroutine init_cnst_hi(q)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test tracer 3: high
! 
!-----------------------------------------------------------------------

   implicit none

!Arguments
   real(r8), intent(out) :: q(plon,plev,plat)    ! kg tracer/kg dry air

! Local
  integer indx

!-----------------------------------------------------------------------
!
! Initialize high tracer to zero except at 200 level
!

  indx = setpindxtr(200.)

  if ( smooth ) then 
     call setsmoothtr(indx,q,.3)
  else
     q = 0.0
     q(:,indx,:) = 1.0
  endif
  !   print *,'suite 1 init_cnst_hi min/max q',minval(q),maxval(q)


end subroutine init_cnst_hi

!======================================================================
subroutine init_cnst_un(q)

!----------------------------------------------------------------------- 
!
! Purpose:
! Initialize test unit tracer.
!    2) conserved unit tracer
!
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
! copied from tracers.F90:initesttr()  D.Bundy Oct 15 2002
!-----------------------------------------------------------------------
   implicit none

   real(r8), intent(out) :: q(plon,plev,plat)    ! kg tracer/kg dry air
!-----------------------------------------------------------------------
! Initialize conserved unit tracer.

   q = 1.0

end subroutine init_cnst_un


!======================================================================

subroutine init_tr

  !----------------------------------------------------------------------- 
  ! Purpose:
  !  Initialize any datasets that the constituents need. 
  !  Currently does nothing for this suite
  ! D Bundy, May 30, 2003
  !-----------------------------------------------------------------------


end subroutine init_tr

!======================================================================

subroutine timestep_init_tr
!----------------------------------------------------------------------- 
! 
! Purpose: call tracer timestep init processes
!----------------------------------------------------------------------- 


end subroutine timestep_init_tr

!======================================================================

subroutine setsmoothtr(indx,q,width,rev_in)

  implicit none

#include <comhyb.h>

  !Arguments
  integer, intent(in)     :: indx               ! k index of pressure level
  real(r8), intent(inout) :: q(plon,plev,plat)  ! kg tracer/kg dry air
  real(r8), intent(in)    :: width              ! eta difference from unit level where q = 0.1
  integer,  intent(in), optional :: rev_in      ! reverse the mixing ratio


  !Local variables
  integer k
  real(r8) alpha ! guassian width, determined by width, T
  real(r8) pdist ! pressure distance (eta.e4) from k=indx
  real(r8) T     ! desired m.r. in level specified by pdiff from k=indx
  integer rev  ! = 1 then reverse (q = 1, q(k=indx) = 0 )




  rev = 0
  if (present(rev_in)) then
     if (rev_in == 1) then
        rev = 1
     endif
  endif


  write(6,*)'DRB TR SMOOTH indx hypm(indx)',indx,hypm(indx)
  write(6,67)'DRB TR SMOOTH ','k','hypm(k)','pdist','-a*(pd^2)','e^-a*(pd^2)'

  T = 0.1
  alpha = -log(T)/(width*1.e4)**2  ! s.t. in level width from indx, mr = T

!  alpha = 3.e-8  ! m.r. ~ 0.1 in adjacent levels, where change eta ~ 0.08

  do k=1,pver
     pdist = hypm(k) - hypm(indx)

     if ( rev == 1 ) then
        q(:,k,:) = 1.0 - exp(-alpha*(pdist**2))
     else
        q(:,k,:) =       exp(-alpha*(pdist**2))
  endif
     write(6,66)'DRB TR SMOOTH ',k,hypm(k),pdist,-alpha*pdist**2,q(1,k,1)
  end do
  
66 format (a15,i4,3f15.2,g15.6)
67 format (a15,a4,4a15)
  
!  call endrun('tracers_suite.F90 L278 DRBDBG')

end subroutine setsmoothtr


!======================================================================

integer function setpindxtr(pmb)

  ! find the index of layer nearest pmb

  use pmgrid, only: plev, plevp

  implicit none
  
#include <comhyb.h>
  
  real(r8) pmb
  real(r8) pmin, pdist
  integer indx, k
  
  indx = 0
  pmin = 1.e36
  pdist = 1.e36
  do k=1,pver
     pdist = abs(hypm(k) - pmb*100.)
     if (pdist < pmin) then
        indx = k
        pmin = pdist
     end if
  end do
  write (6,*) ' index near', pmb, ' is ', indx
  setpindxtr = indx

end function setpindxtr

!======================================================================
!======================================================================

subroutine tend_tr(m,ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of radon tracer
! 
! Method:
!-------------------------Code History----------------------------------
!
! tracers.F90:rndecay()
! Original version:  B. Eaton, 1995
! Standardized:      T. Acker, Feb 1996
!
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: m                 ! tracer number
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    ! radon mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg tracer /(s kg moist air))
!--------------------------Local Variables------------------------------

  if ( m > trac_ncnst ) then
      write(6,*) 'tracers_suite:tend_tr()'
      write(6,*) ' asked to calculate tendency for tracer number ',m
      write(6,*) ' but there are only trac_ncnst = ',trac_ncnst,' tracers'
      call endrun('tracers_suite.F90:tend_tr() L457')
   endif

   if ( m == 1 ) then
      call tend_lw(ncol, q, deltat, snk)
   else if ( m == 2 ) then
      call tend_md(ncol, q, deltat, snk)
   else if ( m == 3 ) then
      call tend_hi(ncol, q, deltat, snk)
   else if ( m == 4 ) then
      call tend_md(ncol, q, deltat, snk)
   else if ( m == 5 ) then
      call tend_un(ncol, q, deltat, snk)


   else
      write(6,*) 'tracers_suite:tend_tr()'
      write(6,*) 'no tendency routine specified for tracer',m
      call endrun
   endif
      


end subroutine tend_tr

!======================================================================

subroutine tend_lw(ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of test tracer 1 (low) (null)
! 
! Method: null
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    ! low mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg q /(s kg moist air))
!--------------------------Local Variables------------------------------
!
   snk = 0.
!   print *,'suite 1 tend_lw min/max snk',minval(snk),maxval(snk)

end subroutine tend_lw

!======================================================================

subroutine tend_md(ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of test tracer 2 (med) = null
! 
! Method: null
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    ! low mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg q /(s kg moist air))
!--------------------------Local Variables------------------------------
   snk = 0.;
!   print *,'suite 1 tend_md min/max snk',minval(snk),maxval(snk)


end subroutine tend_md

!======================================================================

subroutine tend_hi(ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of test tracer 3 (high) =  null
! 
! Method:  null 
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    !  mixing ratio (kg/(kg moist air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg tracer /(s kg moist air))
!--------------------------Local Variables------------------------------

   snk = 0.

!   print *,'suite 1 tend_hi min/max snk',minval(snk),maxval(snk)


end subroutine tend_hi

!======================================================================

subroutine tend_un(ncol, q, deltat, snk)
!----------------------------------------------------------------------- 
!
! Purpose: Calculate tendency (decay) of test tracer 2 (unit) = null
! 
! Method: null
!
!-----------------------------------------------------------------------
   implicit none
!-------------------------Arguments--------------------------------------
!
   integer,  intent(in)  :: ncol              ! number of atmospheric columns
   real(r8), intent(in)  :: q(pcols,pver)    ! radon mixing ratio (kg/(kg dry air))
   real(r8), intent(in)  :: deltat            ! time step
   real(r8), intent(out) :: snk(pcols,pver) ! conversion rate
                                              ! (kg rn /(s kg dry air))
!--------------------------Local Variables------------------------------
   snk = 0.;
   !	print *,'suite 2 tend_un min/max snk',minval(snk),maxval(snk)


end subroutine tend_un

!======================================================================
!======================================================================

subroutine flux_tr(m,ncol, lchnk, landfrac, flux )
!----------------------------------------------------------------------- 
!

!-----------------------------------------------------------------------
   implicit none
!--------------------------Arguments-------------------------------------

   integer,  intent(in)  :: m               ! tracer number
   integer,  intent(in)  :: ncol            ! number of atmospheric columns in chunk
   integer , intent(in)  :: lchnk           ! current identifier
   real(r8), intent(in)  :: landfrac(pcols) ! landfraction
   real(r8), intent(out) :: flux(pcols)     ! specified radon flux in kg/m^2/s

!--------------------------Local Variables------------------------------


   if ( m > trac_ncnst ) then
      write(6,*) 'tracers_suite:flux_tr()'
      write(6,*) ' asked to calculate flux for tracer number ',m
      write(6,*) ' but there are only trac_ncnst = ',trac_ncnst,' tracers'
      call endrun
   endif
   
   ! flux is null for all tracers
   flux = 0.


end subroutine flux_tr


end module tracers_suite
