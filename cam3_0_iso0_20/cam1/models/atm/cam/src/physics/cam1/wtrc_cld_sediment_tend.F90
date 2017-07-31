#include <params.h>
#include <misc.h>

subroutine wtrc_cld_sediment_tend (ncol  ,mtrc  ,dtime  ,pmid   ,pdel   ,cloud ,  &
                                   qvap  ,qliq  ,qice   ,trvap  ,trliq  ,trice ,  &
                                   fxliq ,fxice ,liqtend,icetend,wvtend ,  &
                                   sfliq ,sfice  )

!-----------------------------------------------------------------------
!
! Purpose: Compute sedimentation fluxes for triplets of water species
!          (a single water tracer [isotope] speices... etc)
!          Apply Cloud Particle Gravitational Sedimentation to Condensate
!
!          Here, no fractionation is done. Afer all sedimentation is
!          done, and equilibration is performed. This is approximately
!          correct, excep that the precipitation flux will not be
!          exactly at equilibtrim. 
!
!          Notice, water tracers (isotopes) are passive so don't 
!          contribute to "htend".
!
! Code History: 
!    See also cld_sediment_tend (in pkg_cld_sediment.F90)
!    Based on code skeleton from Andrew G. 
!    David Noone <dcn@caltech.edu> - Thu Aug  7 19:01:25 PDT 2003
!    Sun Wong  <swong@atmos.umd.edu> - (2/17/2003) for relaxation to
!                                              equilibrium for liq.
! Author: 
!    David Noone <dcn@colorado.edu> - Fri Jul  2 18:22:03 MDT 2004
!         (simpler form for cam3)
!
!-----------------------------------------------------------------------
    use shr_kind_mod,    only: r8=>shr_kind_r8
    use ppgrid,          only: pcols, pver, pverp
    use physconst,       only: gravit
    use water_tracers,   only: iwspec, wtrc_ratio
!-----------------------------------------------------------------------
    implicit none
!------------------------- Input Arguments -----------------------------
    integer,  intent(in)  :: ncol                      ! number of colums to process
    integer,  intent(in)  :: mtrc                      ! water tracer index
    real(r8), intent(in)  :: dtime                     ! time step
    real(r8), intent(in)  :: pmid  (pcols,pver)        ! midpoint pressures (Pa)
    real(r8), intent(in)  :: pdel  (pcols,pver)        ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: cloud (pcols,pver)        ! cloud fraction (fraction)
    real(r8), intent(in)  :: qvap(pcols,pver)          ! prognostic vapor  
    real(r8), intent(in)  :: qliq(pcols,pver)          ! prognostic liquid
    real(r8), intent(in)  :: qice(pcols,pver)          ! prognostic ice
    real(r8), intent(in)  :: trvap(pcols,pver)         ! tracer vapor
    real(r8), intent(in)  :: trliq(pcols,pver)         ! tracer liquid
    real(r8), intent(in)  :: trice(pcols,pver)         ! tracerice     
    real(r8), intent(in)  :: fxliq(pcols,pverp)        ! model liq. fluxes at the interfaces, liquid (+down)
    real(r8), intent(in)  :: fxice(pcols,pverp)        ! model ice fluxes at the interfaces, ice    (+down)
!------------------------ Output Arguments -----------------------------
    real(r8), intent(out) :: liqtend(pcols,pver)       ! liquid condensate tend
    real(r8), intent(out) :: icetend(pcols,pver)       ! ice condensate tend
    real(r8), intent(out) :: wvtend (pcols,pver)       ! water vapor tend
    real(r8), intent(out) :: sfliq  (pcols)            ! surface flux of liquid (rain, kg/m/s)
    real(r8), intent(out) :: sfice  (pcols)            ! surface flux of ice    (snow, kg/m/s)
!------------------------- Local Variables- ----------------------------
!
    real(r8) :: fxliqm(pcols,pverp)                    ! (limited) tracer liq. flux
    real(r8) :: fxicem(pcols,pverp)                    ! (limited) tracer ice flux
    real(r8) :: cldab(pcols)                           ! cloud in layer above
    real(r8) :: evapliq                                ! evap of cloud liquid into environment
    real(r8) :: evapice                                ! evap of cloud ice into environment
    real(r8) :: cldovrl                                ! cloud overlap factor
    real(r8) :: rat				       ! tracer ratio

    integer :: i,k

!-----------------------------------------------------------------------
    real(r8) :: mxsedfac = 1.0_r8-1.e-13		       ! trivial
!-----------------------------------------------------------------------
! initialize variables
    fxliqm (:ncol,:) = 0. ! flux at interfaces (liquid)
    fxicem (:ncol,:) = 0. ! flux at interfaces (ice)
    liqtend(:ncol,:) = 0. ! condensate tend (liquid)
    icetend(:ncol,:) = 0. ! condensate tend (ice)
    wvtend(:ncol,:)  = 0. ! environmental moistening
    sfliq(:ncol)     = 0. ! condensate sedimentation flux out bot of column (liquid)
    sfice(:ncol)     = 0. ! condensate sedimentation flux out bot of column (ice)
!
! Compute tracer advective flux by simple up-stream method
! Here we differ from main code in we allow up to all but a tiny bit of
! the mass to be sedimented if needed. This is because the mxsedfac is
! already included in fxliq and fxice.
!
    do k = 1,pver
       do i = 1,ncol
         fxliqm(i,k+1) = fxliq(i,k+1)*wtrc_ratio(trliq(i,k), qliq(i,k))
         fxliqm(i,k+1) = min( fxliqm(i,k+1), mxsedfac * trliq(i,k) * pdel(i,k) )

         fxicem(i,k+1) = fxice(i,k+1)*wtrc_ratio(trice(i,k), qice(i,k))
         fxicem(i,k+1) = min( fxicem(i,k+1), mxsedfac * trice(i,k) * pdel(i,k) )
       end do
    end do

! Now calculate the tendencies assuming that condensate evaporates when
! it falls into environment, and does not when it falls into cloud.
! All flux out of the layer comes from the cloudy part.
! Assume maximum overlap for stratiform clouds
!  if cloud above < cloud,  all water falls into cloud below
!  if cloud above > cloud,  water split between cloud  and environment
!
! ISOTOPES: Do a direct non-fractionating advection. 
! This is quite sound, as we can assume the evaporation in complete.
! Later, after next physics update, equilibrate the cloud liquid with
! the vapour. This will be VERY close to the more accurate "implict" calulation.
!
    do k = 1,pver
       cldab(:ncol) = 0.
       do i = 1,ncol
! cloud overlap cloud factor
          cldovrl  = min( cloud(i,k) / (cldab(i)+.0001), 1. )
          cldab(i) = cloud(i,k)

! evaporation into environment cause moistening and cooling
          evapliq = fxliqm(i,k) * (1.-cldovrl) / (dtime * pdel(i,k))  !  into env (kg/kg/s)
          evapice = fxicem(i,k) * (1.-cldovrl) / (dtime * pdel(i,k))  !  into env (kg/kg/s)
          wvtend(i,k) = evapliq + evapice                          !  evaporation into environment (kg/kg/s)

! net flux into cloud changes cloud liquid/ice (all flux is out of cloud)
          liqtend(i,k)  = (fxliqm(i,k)*cldovrl - fxliqm(i,k+1)) / (dtime * pdel(i,k))
          icetend(i,k)  = (fxicem(i,k)*cldovrl - fxicem(i,k+1)) / (dtime * pdel(i,k))
       end do
    end do

! convert flux out the bottom to mass units Pa -> kg/m2/s
    sfliq(:ncol) = fxliqm(:ncol,pverp) / (dtime*gravit)
    sfice(:ncol) = fxicem(:ncol,pverp) / (dtime*gravit)
!
  return
end subroutine wtrc_cld_sediment_tend


