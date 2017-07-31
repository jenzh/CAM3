!======================================================================
! This is an interface for 3 test tracers: tr1,tr2,tr3
!
! This uses the tracers_suite module to initialize 
!   mixing ratios, fluxes and calculate tendencies.
!   Details of calling tree below. All of the detailed information about the
!  tracers should be store in the suite file, including the number & names of tracers. 
!
! Author B. Eaton
! History  D. Bundy, June 2003 modified to the format of physics interface
!        
!
!---------------------------------------------------------------
!
!  ------------  calling tree --------------
!  Register the tracers as advected fields, pass names to model
!  initindx.F90:			call tracers_register()
!
!  Initialize the tracer mixing ratio field
!  inidat.F90:read_inidat
!  	-> tracers.F90: tracers_init_cnst 
!  		-> tracers_suite.F90:init_cnst_tr
!
!  Initialize data set, things that need to be done at the beginning of a 
!  run (whether restart or initial)
!  inti.F90
!  	-> tracers.F90: tracers_init
!  		-> tracers_suite.F90:init_tr
!  		-> addfld/add default for surface flux (SF)
!
!  Timestepping:
!  advnce.F90
!  	-> tracers_timestep_init
!  		-> tracers_suite.F90:timestep_init_tr
!
!  tphysac.F90
!  	-> tracers_timestep_tend
!  		-> tracers_suite.F90:flux_tr
!  		-> tracers_suite.F90:tend_tr
!
!======================================================================

#include <misc.h>
#include <params.h>

module tracers

  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  private
  save

! Public interfaces
  public tracers_register                  ! register constituent
  public tracers_implements_cnst           ! true if named constituent is implemented by this package
  public tracers_init_cnst                 ! initialize constituent field
  public tracers_init                      ! initialize history fields, datasets
  public tracers_timestep_tend             ! calculate tendencies
  public tracers_timestep_init             ! interpolate dataset for constituent each timestep
  public set_state_pdry                          ! calculate dry air masses in state variable
  public set_wet_to_dry                         ! calculate dry air masses in state variable
  public set_dry_to_wet

! Data from namelist variables
  logical, public :: tracers_flag  = .false.     ! true => turn on test tracer code, namelist variable

! Private module data

  integer :: ixtrct=-999                             ! index of 1st constituent
  logical :: debug = .false.
  
contains
!======================================================================
subroutine tracers_register
!----------------------------------------------------------------------- 
!
! Purpose: register advected tracers. Called by initindx.F90
!  The registration lets the model know what the tracer names are
!  and returns the index number ixtrct for the constituent array
! 
! Author: D. Bundy
!-----------------------------------------------------------------------

   use physconst,    only: mwdry, cpair
   use constituents, only: cnst_add, advected
   use tracers_suite, only: get_tracer_name, trac_ncnst
   
   implicit none
!---------------------------Local workspace-----------------------------
   integer :: mm,m                                 ! dummy
   character(len=8) :: name   ! constituent name
   real(r8) minc

!-----------------------------------------------------------------------
   if ( tracers_flag ) then 
      minc = 0        ! min mixing ratio (normal setting)
      minc = -1.e36   ! min mixing ratio (disable qneg3)
      
      do m = 1,trac_ncnst 
         name = get_tracer_name(m)  ! get name from suite file
         
         ! add constituent name to list of advected, save index number ixtrct
         call cnst_add(name, advected, mwdry, cpair, minc, mm, &  
              readiv=.false.,mixtype='dry')
         if ( m .eq. 1 ) ixtrct = mm  ! save index number of first tracer
         
      end do
   end if

end subroutine tracers_register
!======================================================================

function tracers_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this package
! 
! Author: B. Eaton
! 
!-----------------------------------------------------------------------

  use tracers_suite, only: trac_ncnst, get_tracer_name
  
  implicit none
!-----------------------------Arguments---------------------------------
  
  character(len=*), intent(in) :: name   ! constituent name
  logical :: tracers_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
   integer :: m
!-----------------------------------------------------------------------

   tracers_implements_cnst = .false.
   if ( tracers_flag ) then 
      do m = 1, trac_ncnst
         if (name == get_tracer_name(m)) then
            tracers_implements_cnst = .true.
            return
         end if
      end do
   end if
end function tracers_implements_cnst

!===============================================================================
subroutine tracers_init_cnst(name, q)

!----------------------------------------------------------------------- 
!
! Purpose: initialize test tracers mixing ratio fields 
!  This subroutine is called at the beginning of an initial run ONLY
!
!-----------------------------------------------------------------------

  use pmgrid,     only: plon, plev, plat
  use tracers_suite,   only: init_cnst_tr, trac_ncnst, get_tracer_name

  implicit none

  character(len=*), intent(in) :: name
  real(r8), intent(out), dimension(plon,plev,plat) :: q    ! kg tracer/kg dry air

! Local
  integer m
  if ( tracers_flag ) then 
     do m = 1, trac_ncnst
        if (name ==  get_tracer_name(m))  then
           call init_cnst_tr(m,q)
        endif
     end do
  end if

end subroutine tracers_init_cnst

!===============================================================================
subroutine tracers_init

!----------------------------------------------------------------------- 
!
! Purpose: declare history variables, initialize data sets
!  This subroutine is called at the beginning of an initial or restart run
!
!-----------------------------------------------------------------------

  use tracers_suite,   only: init_tr, trac_ncnst
  use history,    only: addfld, add_default, phys_decomp
  use pmgrid,     only:  plev, plat

  implicit none

! Local
  integer m

  if ( tracers_flag ) then     
     

     ! this is where we would do addfld calls
     ! add surface flux and tendency fields. The tracer mixing ratios are added in history.F90
     
     ! initialize datasets, etc, needed for constituents.
     call init_tr  
  endif
     
end subroutine tracers_init

!======================================================================

subroutine tracers_timestep_init
!----------------------------------------------------------------------- 
!
! Purpose: At the beginning of a timestep, there are some things to do
! that just the masterproc should do. This currently just interpolates
! the emissions boundary data set to the current time step.
!
!-----------------------------------------------------------------------

  use tracers_suite, only: timestep_init_tr

  if ( tracers_flag ) then 
     
     call timestep_init_tr
     
     if (debug) write(6,*)'tracers_timestep_init done'
  endif

end subroutine tracers_timestep_init

!======================================================================

subroutine tracers_timestep_tend(state, ptend, cflx, landfrac, deltat)

!----------------------------------------------------------------------- 
!
! Purpose: During the timestep, compute test tracer mixing ratio 
! tendencies and surface fluxes.
! 
! Author: D. Bundy
!-----------------------------------------------------------------------

  use physics_types, only: physics_state, physics_ptend, physics_ptend_init
  use ppgrid,        only: pcols, pver
  use constituents,  only: ppcnst, sflxnam
  use tracers_suite, only: flux_tr,tend_tr,trac_ncnst
  use history,    only: outfld

  implicit none

  ! Arguments
   type(physics_state), intent(in)  :: state          ! state variables
   type(physics_ptend), intent(out) :: ptend          ! package tendencies
   real(r8),            intent(in)  :: deltat         ! timestep
   real(r8),            intent(in)  :: landfrac(pcols) ! Land fraction
   real(r8), intent(inout) :: cflx(pcols,ppcnst) ! Surface constituent flux (kg/m^2/s)

! Local variables
   integer  :: m               ! tracer number (internal)

!-----------------------------------------------------------------------

     ! Initialize output tendency structure
     ! Do this even if not going to be used, to clear it and 
     !    pass runtime debug checks 
     call physics_ptend_init(ptend)

  if ( tracers_flag ) then 
     ptend%name  = 'tracers'
     
     do  m = 1,trac_ncnst 
        if (debug) print *,'tracers.F90 calling for tracer ',m
        
        !calculate flux
        call flux_tr(m,state%ncol,state%lchnk, landfrac, cflx(:,ixtrct+m-1))
        
        !calculate tendency
        call tend_tr(m,state%ncol, state%q(:,:,ixtrct+m-1), deltat, ptend%q(:,:,ixtrct+m-1))
        ptend%lq(ixtrct+m-1) = .true.    ! logical flag to update physics tendency
        
        !outfld calls could go here
        call outfld (sflxnam(ixtrct+m-1),cflx(:,ixtrct+m-1),pcols,state%lchnk)
        
        
     end do
     
     if ( debug ) then 
        do  m = 1,trac_ncnst 
           print *,'tracers_timestep_tend ixtrct,m,ixtrct+m-1',ixtrct,m,ixtrct+m-1
           print *,'tracers_timestep_tend min max flux',minval(cflx(:,ixtrct+m-1)),maxval(cflx(:,ixtrct+m-1))
           print *,'tracers_timestep_tend min max tend',minval(ptend%q(:,:,ixtrct+m-1)),maxval(ptend%q(:,:,ixtrct+m-1))
        end do
        write(6,*)'tracers_timestep_tend end'
     endif
  endif


end subroutine tracers_timestep_tend

!======================================================================

subroutine set_state_pdry (state)

  use physics_types, only: physics_state
  use ppgrid,  only: pcols, pver
  use pmgrid

  implicit none

#include <comhyb.h>

  type(physics_state), intent(inout) :: state
  integer ncol
  integer i, k
  
  ncol = state%ncol

#ifdef METHOD1
  state%psdry(:ncol) = ps0 * hyai(1)
  do k = 1, pver
     state%psdry(:ncol) = state%psdry(:ncol) + state%pdel(:ncol,k)*(1.-state%q(:ncol,k,1))
  end do
  call plevs0(ncol, pcols, pver, state%psdry, &
       state%pintdry, state%pmiddry, state%pdeldry)
#else
  state%psdry(:ncol) = ps0 * hyai(1)
  state%pintdry(:ncol,1) = ps0 * hyai(1)
  do k = 1, pver
     state%pdeldry(:ncol,k) = state%pdel(:ncol,k)*(1.-state%q(:ncol,k,1))
     state%pintdry(:ncol,k+1) = state%pintdry(:ncol,k)+state%pdeldry(:ncol,k)
     state%pmiddry(:ncol,k) = (state%pintdry(:ncol,k+1)+state%pintdry(:ncol,k))/2.
     state%psdry(:ncol) = state%psdry(:ncol) + state%pdeldry(:ncol,k)
  end do
#endif

  state%rpdeldry(:ncol,:) = 1./state%pdeldry(:ncol,:)
  state%lnpmiddry(:ncol,:) = log(state%pmiddry(:ncol,:))
  state%lnpintdry(:ncol,:) = log(state%pintdry(:ncol,:))

end subroutine set_state_pdry 

!======================================================================

subroutine set_wet_to_dry (state)

  use physics_types, only: physics_state
  use ppgrid,  only: pcols, pver
  use pmgrid
  use constituents,   only: cnst_type
  use constituents,  only: ppcnst

  implicit none

#include <comhyb.h>

  type(physics_state), intent(inout) :: state
  integer ncol
  integer i, k, m
  
  ncol = state%ncol

  do m = 1,ppcnst
     if (cnst_type(m).eq.'dry') then
        state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdel(:ncol,:)/state%pdeldry(:ncol,:)
     endif
  end do


end subroutine set_wet_to_dry 

!======================================================================

subroutine set_dry_to_wet (state)

  use physics_types, only: physics_state
  use ppgrid,  only: pcols, pver
  use pmgrid
  use constituents,  only: cnst_type
  use constituents,  only: ppcnst

  implicit none

#include <comhyb.h>

  type(physics_state), intent(inout) :: state
  integer ncol
  integer i, k, m
  
  ncol = state%ncol

  do m = 1,ppcnst
     if (cnst_type(m).eq.'dry') then
        if ( debug ) write(6,*)'DRBDBG set_dry_to_wet cnst',m
        state%q(:ncol,:,m) = state%q(:ncol,:,m)*state%pdeldry(:ncol,:)/state%pdel(:ncol,:)
     endif
  end do

end subroutine set_dry_to_wet

!======================================================================


end module tracers
