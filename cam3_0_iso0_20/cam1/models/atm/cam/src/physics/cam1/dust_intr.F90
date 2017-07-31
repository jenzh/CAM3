#include <misc.h>
#include <params.h>

module dust_intr

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
  public dust_register_cnst                        ! register consituents
  public dust_implements_cnst                      ! returns true if consituent is implemented by this package
  public dust_init_cnst                            ! initialize mixing ratios if not read from initial file
  public dust_initialize                           ! initialize (history) variables
  public dust_wet_intr                             ! interface to wet deposition
  public dust_emis_intr                            ! interface to emission
  public dust_drydep_intr                          ! interface to tendency computation
  public dust_time_interp                          ! interpolate oxidants and fluxes to current time
  public dust_idx1                                 ! allow other parts of code to know where dust is

contains

!===============================================================================
  subroutine dust_register_cnst
!----------------------------------------------------------------------- 
! 
! Purpose: register advected constituents for all aerosols
! 
! Method: 
!-----------------------------------------------------------------------

    call endrun('DUST_REGISTER_CNST:  code under development')

    return
  end subroutine dust_register_cnst



!=======================================================================
  function dust_implements_cnst(name)
!----------------------------------------------------------------------- 
! 
! Purpose: return true if specified constituent is implemented by this 
!          package
! 
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------

     character(len=*), intent(in) :: name   ! constituent name
     logical :: dust_implements_cnst        ! return value
!-----------------------------------------------------------------------

     call endrun('DUST_IMPLEMENTS_CNST:  code under development')

  end function dust_implements_cnst


!=======================================================================
  subroutine dust_init_cnst(name, q)
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
    
    call endrun('DUST_INIT_CNST:  code under development')

  end subroutine dust_init_cnst



  function dust_idx1()
    implicit none
    integer dust_idx1
    call endrun('DUST_IDX1:  code under development')
  end function dust_idx1


!===============================================================================
  subroutine dust_initialize 
!----------------------------------------------------------------------- 
! 
! Purpose: initialize parameterization of dust chemistry
!          (declare history variables)
! 
! Method: 
!-----------------------------------------------------------------------

    call endrun('DUST_INITIALIZE:  code under development')

  end subroutine dust_initialize


!===============================================================================
  subroutine dust_wet_intr (state, ptend, cflx, nstep, dt, lat, clat, cme, prain, &
       evapr, cldv, cldc, cldn, fracis, calday, cmfdqr, conicw, rainmr)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to wet processing of aerosols (source and sinks).
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

    call endrun('DUST_WET_INTR:  code under development')

  end subroutine dust_wet_intr

  subroutine dust_drydep_intr (state, ptend, wvflx, dt, lat, clat, &
       fsds, obklen, ts, ustar, prect, snowh, pblh, hflx, month, landfrac, &
       icefrac, ocnfrac,fv,ram1)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to dry deposition and sedimentation of dust
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
    real(r8), intent(in) :: obklen(pcols)                 ! obklen
    real(r8), intent(in) :: ustar(pcols)                  ! sfc fric vel--used over oceans and sea ice.
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
    real(r8), intent(in) :: fv(pcols)        ! for dry dep velocities from land model for dust
    real(r8), intent(in) :: ram1(pcols)       ! for dry dep velocities from land model for dust

    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies

!-----------------------------------------------------------------------

    call endrun('DUST_DRYDEP_INTR:  code under development')

    return
  end subroutine dust_drydep_intr

  subroutine dust_time_interp

    call endrun('DUST_TIME_INTERP:  code under development')

  end subroutine dust_time_interp

  subroutine dust_emis_intr (state, ptend, dt,cflx)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interface to emission of all dusts.
! Notice that the mobilization is calculated in the land model (need #define BGC) and
! the soil erodibility factor is applied here.
! 
! Method: 
!-----------------------------------------------------------------------
    use shr_kind_mod,only: r8 => shr_kind_r8
    use physics_types, only: physics_state, physics_ptend
    use ppgrid,      only: pcols
    use constituents,only: ppcnst

    implicit none
!-----------------------------------------------------------------------
!
! Arguments:
!
    real(r8),            intent(in)  :: dt             ! time step
    type(physics_state), intent(in ) :: state          ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend          ! indivdual parameterization tendencies
    real(r8),            intent(inout) :: cflx(pcols,ppcnst)

    call endrun('DUST_EMIS_INTR:  code under development')

    return
  end subroutine dust_emis_intr 

end module dust_intr
