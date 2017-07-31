#include <misc.h>
#include <params.h>

module seasalt_intr

!---------------------------------------------------------------------------------
! Module to interface the aerosol parameterizations with CAM
! written by PJR (extensively modified from chemistry module)
!---------------------------------------------------------------------------------

  use abortutils,  only: endrun
    
  implicit none

  private          ! Make default type private to the module

  save

! Public interfaces
!
  public seasalt_register_cnst                        ! register consituents
  public seasalt_implements_cnst                      ! returns true if consituent is implemented by this package
  public seasalt_init_cnst                            ! initialize mixing ratios if not read from initial file
  public seasalt_initialize                           ! initialize (history) variables
  public seasalt_srcsnk                               ! production and loss of seasalt

contains

!===============================================================================
  subroutine seasalt_register_cnst
!----------------------------------------------------------------------- 
! 
! Purpose: register seasalt
! 
! Method: 
!-----------------------------------------------------------------------

    call endrun('SEASALT_REGISTER_CNST:  code under development')

    return
  end subroutine seasalt_register_cnst


!=======================================================================
  function seasalt_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this 
!          package
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: seasalt_implements_cnst     ! return value
!-----------------------------------------------------------------------

     call endrun('SEASALT_IMPLEMENTS_CNST:  code under development')

  end function seasalt_implements_cnst


!=======================================================================
  subroutine seasalt_init_cnst(name, q)
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
    
    call endrun('SEASALT_INIT_CNST:  code under development')

  end subroutine seasalt_init_cnst


!===============================================================================
  subroutine seasalt_initialize 
!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterization of seasalt 
!          (declare history variables)
! 
! Method: 
!-----------------------------------------------------------------------

    call endrun('SEASALT_INITIALIZE:  code under development')

  end subroutine seasalt_initialize


!===============================================================================
  subroutine seasalt_srcsnk (state, dt, u10, ocnfrac, ptend)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to emission of all seasalts
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
    type(physics_state), intent(in ) :: state          ! Physics state variables
    real(r8),            intent(in)  :: dt             ! time step
    real(r8),            intent(in)  :: u10(pcols)     ! 10 meter wind (m/s)
    real(r8),            intent(in)  :: ocnfrac(pcols) ! ocean fraction
    type(physics_ptend), intent(inout) :: ptend        ! indivdual parameterization tendencies

    call endrun('SEASALT_SRCSNK:  code under development')

  end subroutine seasalt_srcsnk 

end module seasalt_intr
