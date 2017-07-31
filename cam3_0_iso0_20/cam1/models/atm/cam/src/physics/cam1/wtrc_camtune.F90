#include <params.h>
#include <misc.h>

module wtrc_camtune
!-----------------------------------------------------------------------
!
! Module containing tunable parameters for water isotope tracers
! These are specific to the moist parameterizations in cam.
! Notice, other tunable parameter that are not specific to cam, but
! are still specific to the isotope scheme are found in water_isotopes.F90
!
! All of these control the isotope physics. Controls on the numerics
! remain in the individual routines (qtiny for instance).
!
! Author: David Noone <dcn@colorado.edu> - Mon Aug 16 16:20:52 MDT 2004
!
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8

  implicit none
  save
  public
!-----------------------------------------------------------------------
!
! Zhang-McFarlane's deep convection
!  (zm_convr.F90, wtrc_zm_convr.F90)
!
   real(r8), parameter :: feq_upd  = 1.0    ! equilibration of drops in updraft (dicm)
   real(r8), parameter :: feq_dnd  = 1.0    ! equilibration of rain with downdraft (dicm)
!   real(r8), parameter :: feq_zme  = 0.35   ! equilibration of ZM rain with environment
!   real(r8), parameter :: feq_zme  = 0.55   ! equilibration of ZM rain with environment
   real(r8), parameter :: feq_zme  = 0.65   ! equilibration of ZM rain with environment
!   real(r8), parameter :: feq_hke  = 0.55   ! equilibration of Hack rain with environment
!   real(r8), parameter :: feq_hke  = 0.85   ! equilibration of Hack rain with environment
   real(r8), parameter :: feq_hke  = 0.95   ! equilibration of Hack rain with environment
!   real(r8), parameter :: feq_hke  = 1.00   ! equilibration of Hack rain with environment

!   real(r8), parameter :: fsnk_zmc = 1.0    ! fraction detrain in precip sink (oethwise retained)
!   real(r8), parameter :: fsnk_zmc = 0.5    ! fraction detrain in precip sink (oethwise retained)
   real(r8), parameter :: fsnk_zmc = 0.0    ! fraction detrain in precip sink (oethwise retained)
!
   real(r8), parameter :: ftc_dnd  = 0.8    ! fraction of cloud temp for fractn.
   real(r8), parameter :: fcloud   = 1.0    ! fraction of cloudy sky into which drops evaporate
!
   logical , parameter :: lupdcorr = .true. ! apply updraft corrections
   logical , parameter :: ltndcorr = .true. ! apply tendency corrections
!
! Hack's convective mass flux formulation of moist convective adjustment
! (moistconvect.F90)
!
   real(r8), parameter :: feq_cmf  = 1.00   ! equilibration of drops in updraft (dicm)
!   real(r8), parameter :: fsnk_cmf = 1.00   ! fraction detrained water going to isotope sink 
!   real(r8), parameter :: fsnk_cmf = 0.50   ! fraction detrained water going to isotope sink 
   real(r8), parameter :: fsnk_cmf = 0.00   ! fraction detrained water going to isotope sink 
!
! Prognostic cloud water
! (cldcond.F90)
!
   real(r8), parameter :: feq_cld  = 1.0    ! fraction cloud water equilibration with vapour
   real(r8), parameter :: fsnk_cld = 1.0    ! fraction of precipitation to isotope sink
!
! Isotopic equilibration of cloud liquid (for sedimentation)
! (wiso_cldliq_equilibrate.F90)
!
   real(r8), parameter :: feq_liq  = 1.0    ! liquid fraction equilibrated
   real(r8), parameter :: feq_ice  = 0.01   ! ice fraction equilibrated
!
!-----------------------------------------------------------------------
end module wtrc_camtune

