#include <misc.h>
#include <params.h>

module iop
#if ( defined BFB_CAM_SCAM_IOP ) || ( defined SCAM )
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: iop
! 
! !DESCRIPTION: 
! iop specific routines
!
! !USES:
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst, cnst_name,ppcnst
  use pmgrid, only: plond,plev,plat,beglat,endlat
  use string_utils, only: to_lower
!
! !PUBLIC TYPES:
  implicit none
   real(r8), allocatable :: betasav(:)
   real(r8), allocatable :: fixmassav(:)
   real(r8), allocatable :: alphasav(:,:)
   real(r8), allocatable :: clat_plond(:)          ! d(ps)/dt
   real(r8), allocatable :: alpha_plond(:,:)   
   real(r8), allocatable :: fixmas_plond(:)         
   real(r8), allocatable :: beta_plond(:)          
   real(r8), allocatable :: dqfx3sav(:,:,:)       
   real(r8), allocatable :: dqfx3savm1(:,:,:,:)       
   real(r8), allocatable :: divq3dsav(:,:,:,:)
   real(r8), allocatable :: divt3dsav(:,:,:)       
   real(r8), allocatable :: t3sav(:,:,:)       
   real(r8), allocatable :: u3sav(:,:,:)       
   real(r8), allocatable :: v3sav(:,:,:)       
   real(r8), allocatable :: t2sav(:,:,:)         ! temp tendency
   real(r8), allocatable :: q3sav(:,:,:,:)
   real(r8), allocatable :: pssav(:,:)
   real(r8), allocatable :: tssav(:,:)
   character(len=8) alphanam(pcnst)       ! alpha fixer
   character(len=8) dqfxnam(pcnst)       ! dq fixer

!
! !PUBLIC MEMBER FUNCTIONS:
  public :: init_iop_fields
!
! !REVISION HISTORY:
! Created by John Truesdale
!
!EOP
!
! !PRIVATE MEMBER FUNCTIONS:
!----------------------------------------------------------------------- 

contains
   subroutine init_iop_fields(ps, t3, u3, v3, q3)
!------------------------------------------------------------------------------
! Coupler for converting dynamics output variables into physics input variables
! also writes dynamics variables (on physics grid) to history file
!------------------------------------------------------------------------------
   use ppgrid,       only: begchunk, endchunk,pcols
   use pmgrid,       only: beglat,endlat,beglatex, endlatex,plond,numlats,plev
   use constituents, only: ppcnst

   implicit none

!------------------------------Arguments--------------------------------
    real(r8), intent(in) :: ps  (plond, beglat:endlat)            ! surface pressure
    real(r8), intent(in) :: t3  (plond, plev, beglatex:beglatex+numlats-1)  ! temperature
    real(r8), intent(in) :: u3  (plond, plev, beglatex:beglatex+numlats-1)  ! u-wind component
    real(r8), intent(in) :: v3  (plond, plev, beglatex:beglatex+numlats-1)  ! v-wind component
    real(r8), intent(in) :: q3  (plond, plev, ppcnst, beglatex:beglatex+numlats-1) ! constituents

!-----------------------------------------------------------------------
        
   if(.not.allocated(betasav)) allocate (betasav(beglat:endlat))
   if(.not.allocated(fixmassav)) allocate (fixmassav(beglat:endlat))
   if(.not.allocated(alphasav)) allocate (alphasav(pcnst,beglat:endlat))
   if(.not.allocated(clat_plond)) allocate (clat_plond(plond))          ! d(ps)/dt
   if(.not.allocated(alpha_plond)) allocate (alpha_plond(plond,pcnst))
   if(.not.allocated(fixmas_plond)) allocate (fixmas_plond(plond))
   if(.not.allocated(beta_plond)) allocate (beta_plond(plond))
   if(.not.allocated(dqfx3sav)) allocate (dqfx3sav(plond,plev,pcnst))
   if(.not.allocated(dqfx3savm1)) allocate (dqfx3savm1(plond,plev,pcnst,beglat:endlat))
   if(.not.allocated(divq3dsav)) allocate (divq3dsav(plond,plev,ppcnst,beglat:endlat))
   if(.not.allocated(divt3dsav)) allocate (divt3dsav(plond,plev,beglat:endlat))
   if(.not.allocated(t3sav)) allocate (t3sav(plond,plev,beglat:endlat))
   if(.not.allocated(u3sav)) allocate (u3sav(plond,plev,beglat:endlat))
   if(.not.allocated(v3sav)) allocate (v3sav(plond,plev,beglat:endlat))
   if(.not.allocated(t2sav)) allocate (t2sav(plond,plev,beglat:endlat))  ! temp tendency
   if(.not.allocated(q3sav)) allocate (q3sav(plond,plev,ppcnst,beglat:endlat))
   if(.not.allocated(pssav)) allocate (pssav(plond,beglat:endlat))
   if(.not.allocated(tssav)) allocate (tssav(plond,beglat:endlat))
   t3sav=t3
   u3sav=u3
   v3sav=v3
   q3sav=q3
   pssav=ps
  end subroutine init_iop_fields

#endif

end module iop

