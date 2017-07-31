#include <misc.h>
#include <params.h>

module carbon_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! written by PJR (extensively modified from chemistry module)
!---------------------------------------------------------------------------------

  use abortutils,  only: endrun
    

  implicit none

  private          ! Make default type private to the module

  save

!
! Public interfaces
!
  public carbon_register_cnst                        ! register consituents
  public carbon_implements_cnst                      ! returns true if consituent is implemented by this package
  public carbon_init_cnst                            ! initialize mixing ratios if not read from initial file
  public carbon_initialize                           ! initialize (history) variables
  public carbon_wet_intr                             ! interface to wet deposition
  public carbon_emis_intr                            ! interface to emission
  public carbon_drydep_intr                          ! interface to tendency computation
  public carbon_time_interp                          ! interpolate oxidants and fluxes to current time

contains

!===============================================================================
  subroutine carbon_register_cnst
!----------------------------------------------------------------------- 
! 
! Purpose: register advected constituents for all aerosols
! 
! Method: 
!-----------------------------------------------------------------------

    call endrun('CARBON_REGISTER_CNST:  code under development')

    return
  end subroutine carbon_register_cnst



!=======================================================================
  function carbon_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this 
!          package
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: carbon_implements_cnst      ! return value
!-----------------------------------------------------------------------

     call endrun('CARBON_IMPLEMENTS_CNST:  code under development')

  end function carbon_implements_cnst


!=======================================================================
  subroutine carbon_init_cnst(name, q)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Set initial mass mixing ratios.
!
!-----------------------------------------------------------------------
    use shr_kind_mod,only: r8 => shr_kind_r8

    implicit none
!-----------------------------Arguments---------------------------------
    
    character(len=*), intent(in) :: name         ! constituent name
    
    real(r8), intent(out) :: q(:,:,:)            !  mass mixing ratio
!-----------------------------------------------------------------------
    
    call endrun('CARBON_INIT_CNST:  code under development')

  end subroutine carbon_init_cnst


!===============================================================================
  subroutine carbon_initialize 
!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterization of carbon chemistry
!          (declare history variables)
! 
! Method: 
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------

    call endrun('CARBON_INITIALIZE:  code under development')

  end subroutine carbon_initialize


!===============================================================================
  subroutine carbon_wet_intr (state, ptend, cflx, nstep, dt, lat, clat, cme, prain, &
       evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to we processing of aerosols (source and sinks).
! 
! Method: 
!-----------------------------------------------------------------------
    use shr_kind_mod,only: r8 => shr_kind_r8
    use physics_types, only: physics_state, physics_ptend
    use ppgrid,      only: pcols, pver
    use constituents,only: ppcnst

    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    integer, intent(in) :: nstep
    integer, intent(in) :: lat(pcols)                  ! latitude index for S->N storage
    real(r8), intent(in) :: clat(pcols)                    ! latitude 
    real(r8), intent(in) :: cme(pcols,pver)            ! local condensation of cloud water
    real(r8), intent(in) :: prain(pcols,pver)            ! production of rain
    real(r8), intent(in) :: evapr(pcols,pver)            ! evaporation of rain
    real(r8), intent(in) :: cldn(pcols,pver)            ! cloud fraction
    real(r8), intent(in) :: cldc(pcols,pver)            ! convective cloud fraction
    real(r8), intent(in) :: cldv(pcols,pver)            ! cloudy volume undergoing scavenging

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
    real(r8), intent(inout)  :: cflx(pcols,ppcnst)       ! Surface constituent flux (kg/m^2/s)

    real(r8), intent(inout) :: fracis(pcols,pver,ppcnst)         ! fraction of transported species that are insoluble

    real(r8), intent(in) :: conicw(pcols, pver)
    real(r8), intent(in) :: cmfdqr(pcols, pver)
    real(r8), intent(in) :: rainmr(pcols, pver) ! rain mixing ratio
    real(r8) :: calday        ! current calendar day

!-----------------------------------------------------------------------

    call endrun('CARBON_WET_INTR:  code under development')

    return

  end subroutine carbon_wet_intr

  subroutine carbon_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
       fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, month, landfrac, &
       icefrac, ocnfrac)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to parameterized greenhouse gas chemisty (source/sink).
! 
! Method: 
!-----------------------------------------------------------------------
    use shr_kind_mod,only: r8 => shr_kind_r8
    use physics_types, only: physics_state, physics_ptend
    use ppgrid,      only: pcols

    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    integer, intent(in) :: lat(pcols)                  ! latitude index for S->N storage
    real(r8), intent(in) :: clat(pcols)                 ! latitude 
    real(r8), intent(in) :: fsds(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: obklen(pcols)                 ! longwave down at sfc
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel
    real(r8), intent(in) :: ts(pcols)                     ! sfc temp
    real(r8), intent(in) :: landfrac(pcols)               ! land fraction
    real(r8), intent(in) :: icefrac(pcols)                ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)                ! ocean fraction
    real(r8), intent(in) :: hflx(pcols)                  ! sensible heat flux
    real(r8), intent(in) :: prect(pcols)                     ! prect
    real(r8), intent(in) :: snowh(pcols)                     ! snow depth
    real(r8), intent(in) :: pblh(pcols)                     ! pbl height
    integer, intent(in)  :: month
    real(r8), intent(in) :: wvflx(pcols)       ! water vapor flux

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies

!-----------------------------------------------------------------------

    call endrun('CARBON_DRYDEP_INTR:  code under development')

    return
  end subroutine carbon_drydep_intr

  subroutine carbon_time_interp

    call endrun('CARBON_TIME_INTERP:  code under development')

  end subroutine carbon_time_interp

  subroutine carbon_emis_intr (state, ptend, dt)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to emission of all carbons
! 
! Method: 
!-----------------------------------------------------------------------
    use shr_kind_mod,only: r8 => shr_kind_r8
    use physics_types, only: physics_state, physics_ptend

    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies

    call endrun('CARBON_EMIS_INTR:  code under development')

    return
  end subroutine carbon_emis_intr 

end module carbon_intr
