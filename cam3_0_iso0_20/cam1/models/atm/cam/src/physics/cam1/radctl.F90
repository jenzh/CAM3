#include <misc.h>
#include <params.h>

subroutine radctl(lchnk   ,ncol    ,                   &
                  lwup    ,emis    ,          &
                  pmid    ,pint    ,pmln    ,piln    ,t       , &
                  qm1     ,cld     ,cicewp  ,cliqwp  ,coszrs  , &
                  asdir   ,asdif   ,aldir   ,aldif   ,pmxrgn  , &
                  nmxrgn  ,fsns    ,fsnt    ,flns    ,flnt    , &
                  qrs     ,qrl     ,flwds   ,rel     ,rei     , &
                  sols    ,soll    ,solsd   ,solld   , &
                  landfrac,zm      ,state, fsds)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Driver for radiation computation.
! 
! Method: 
! Radiation uses cgs units, so conversions must be done from
! model fields to radiation fields.
!
! Author: CCM1,  CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid
   use pspect
   use commap
   use history, only: outfld
   use constituents, only: ppcnst, cnst_get_ind
   use prescribed_aerosols, only: get_aerosol, naer_all, aerosol_diagnostics, &
      aerosol_indirect, get_rf_scales, get_int_scales, radforce, idxVOLC
   use physics_types, only: physics_state
   use wv_saturation, only: aqsat
   use chemistry,    only: trace_gas
   use physconst, only: cpair, epsilo
   use aer_optics, only: idxVIS
   use aerosol_intr, only: set_aerosol_from_prognostics


   implicit none

#include <ptrrgrid.h>
#include <comctl.h>
#include <comsol.h>
!
! Input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   integer nspint            ! Num of spctrl intervals across solar spectrum
   integer naer_groups       ! Num of aerosol groups for optical diagnostics
   parameter ( nspint = 19 )
   parameter ( naer_groups = 7 )    ! current groupings are sul, sslt, all carbons, all dust, background, and all aerosols


   real(r8), intent(in) :: lwup(pcols)          ! Longwave up flux at surface
   real(r8), intent(in) :: emis(pcols,pver)     ! Cloud emissivity
   real(r8), intent(in) :: pmid(pcols,pver)     ! Model level pressures
   real(r8), intent(in) :: pint(pcols,pverp)    ! Model interface pressures
   real(r8), intent(in) :: pmln(pcols,pver)     ! Natural log of pmid
   real(r8), intent(in) :: rel(pcols,pver)      ! liquid effective drop size (microns)
   real(r8), intent(in) :: rei(pcols,pver)      ! ice effective drop size (microns)
   real(r8), intent(in) :: piln(pcols,pverp)    ! Natural log of pint
   real(r8), intent(in) :: t(pcols,pver)        ! Model level temperatures
   real(r8), intent(in) :: qm1(pcols,pver,ppcnst) ! Specific humidity and tracers
   real(r8), intent(in) :: cld(pcols,pver)      ! Fractional cloud cover
   real(r8), intent(in) :: cicewp(pcols,pver)   ! in-cloud cloud ice water path
   real(r8), intent(in) :: cliqwp(pcols,pver)   ! in-cloud cloud liquid water path
   real(r8), intent(in) :: coszrs(pcols)        ! Cosine solar zenith angle
   real(r8), intent(in) :: asdir(pcols)         ! albedo shortwave direct
   real(r8), intent(in) :: asdif(pcols)         ! albedo shortwave diffuse
   real(r8), intent(in) :: aldir(pcols)         ! albedo longwave direct
   real(r8), intent(in) :: aldif(pcols)         ! albedo longwave diffuse
   real(r8), intent(in) :: landfrac(pcols)      ! land fraction
   real(r8), intent(in) :: zm(pcols,pver)       ! Height of midpoints (above surface)
   type(physics_state), intent(in) :: state     
   real(r8), intent(inout) :: pmxrgn(pcols,pverp) ! Maximum values of pmid for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pmid for
!    1st region, pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
   integer, intent(inout) :: nmxrgn(pcols)     ! Number of maximally overlapped regions

    real(r8) :: pmxrgnrf(pcols,pverp)             ! temporary copy of pmxrgn
    integer  :: nmxrgnrf(pcols)     ! temporary copy of nmxrgn

!
! Output solar arguments
!
   real(r8), intent(out) :: fsns(pcols)          ! Surface absorbed solar flux
   real(r8), intent(out) :: fsnt(pcols)          ! Net column abs solar flux at model top
   real(r8), intent(out) :: flns(pcols)          ! Srf longwave cooling (up-down) flux
   real(r8), intent(out) :: flnt(pcols)          ! Net outgoing lw flux at model top
   real(r8), intent(out) :: sols(pcols)          ! Downward solar rad onto surface (sw direct)
   real(r8), intent(out) :: soll(pcols)          ! Downward solar rad onto surface (lw direct)
   real(r8), intent(out) :: solsd(pcols)         ! Downward solar rad onto surface (sw diffuse)
   real(r8), intent(out) :: solld(pcols)         ! Downward solar rad onto surface (lw diffuse)
   real(r8), intent(out) :: qrs(pcols,pver)      ! Solar heating rate
   real(r8), intent(out) :: fsds(pcols)          ! Flux Shortwave Downwelling Surface
!
! Output longwave arguments
!
   real(r8), intent(out) :: qrl(pcols,pver)      ! Longwave cooling rate
   real(r8), intent(out) :: flwds(pcols)         ! Surface down longwave flux


!
!---------------------------Local variables-----------------------------
!
   integer i, k              ! index

   integer :: in2o, ich4, if11, if12 ! indexes of gases in constituent array

   real(r8) solin(pcols)         ! Solar incident flux
!  real(r8) fsds(pcols)          ! Flux Shortwave Downwelling Surface
   real(r8) fsntoa(pcols)        ! Net solar flux at TOA
   real(r8) fsntoac(pcols)       ! Clear sky net solar flux at TOA
   real(r8) fsnirt(pcols)        ! Near-IR flux absorbed at toa
   real(r8) fsnrtc(pcols)        ! Clear sky near-IR flux absorbed at toa
   real(r8) fsnirtsq(pcols)      ! Near-IR flux absorbed at toa >= 0.7 microns
   real(r8) fsntc(pcols)         ! Clear sky total column abs solar flux
   real(r8) fsnsc(pcols)         ! Clear sky surface abs solar flux
   real(r8) fsdsc(pcols)         ! Clear sky surface downwelling solar flux
   real(r8) flut(pcols)          ! Upward flux at top of model
   real(r8) lwcf(pcols)          ! longwave cloud forcing
   real(r8) swcf(pcols)          ! shortwave cloud forcing
   real(r8) flutc(pcols)         ! Upward Clear Sky flux at top of model
   real(r8) flntc(pcols)         ! Clear sky lw flux at model top
   real(r8) flnsc(pcols)         ! Clear sky lw flux at srf (up-down)
   real(r8) ftem(pcols,pver)     ! temporary array for outfld
   real(r8) fln200(pcols)        ! net longwave flux interpolated to 200 mb
   real(r8) fln200c(pcols)       ! net clearsky longwave flux interpolated to 200 mb
   real(r8) fns(pcols,pverp)     ! net shortwave flux
   real(r8) fcns(pcols,pverp)    ! net clear-sky shortwave flux
   real(r8) fsn200(pcols)        ! fns interpolated to 200 mb
   real(r8) fsn200c(pcols)       ! fcns interpolated to 200 mb
   real(r8) fnl(pcols,pverp)     ! net longwave flux
   real(r8) fcnl(pcols,pverp)    ! net clear-sky longwave flux

   real(r8) pbr(pcols,pverr)     ! Model mid-level pressures (dynes/cm2)
   real(r8) pnm(pcols,pverrp)    ! Model interface pressures (dynes/cm2)
   real(r8) o3vmr(pcols,pverr)   ! Ozone volume mixing ratio
   real(r8) o3mmr(pcols,pverr)   ! Ozone mass mixing ratio
   real(r8) eccf                 ! Earth/sun distance factor
   real(r8) n2o(pcols,pver)      ! nitrous oxide mass mixing ratio
   real(r8) ch4(pcols,pver)      ! methane mass mixing ratio
   real(r8) cfc11(pcols,pver)    ! cfc11 mass mixing ratio
   real(r8) cfc12(pcols,pver)    ! cfc12 mass mixing ratio
   real(r8) rh(pcols,pverr)      ! level relative humidity (fraction)
   real(r8) lwupcgs(pcols)       ! Upward longwave flux in cgs units

   real(r8) esat(pcols,pverr)    ! saturation vapor pressure
   real(r8) qsat(pcols,pverr)    ! saturation specific humidity

   real(r8) :: frc_day(pcols) ! = 1 for daylight, =0 for night colums
   real(r8) :: aertau(pcols,nspint,naer_groups) ! Aerosol column optical depth
   real(r8) :: aerssa(pcols,nspint,naer_groups) ! Aerosol column averaged single scattering albedo
   real(r8) :: aerasm(pcols,nspint,naer_groups) ! Aerosol column averaged asymmetry parameter
   real(r8) :: aerfwd(pcols,nspint,naer_groups) ! Aerosol column averaged forward scattering

   real(r8) aerosol(pcols, pver, naer_all) ! aerosol mass mixing ratios
   real(r8) scales(naer_all)               ! scaling factors for aerosols

!
! Interpolate ozone volume mixing ratio to model levels
!
   call radozn(lchnk   ,ncol    ,pmid    ,o3vmr   )
   call outfld('O3VMR   ',o3vmr ,pcols, lchnk)

!
! Set chunk dependent radiation input
!
   call radinp(lchnk   ,ncol    ,                                &
               pmid    ,pint    ,o3vmr   , pbr     ,&
               pnm     ,eccf    ,o3mmr   )
!
! Solar radiation computation
!
   if (dosw) then

!
! calculate heating with aerosols
!
      call aqsat(state%t, state%pmid, esat, qsat, pcols, &
                 ncol, pver, 1, pver)

      ! calculate relative humidity
      rh(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / qsat(1:ncol,1:pver) * &
         ((1.0 - epsilo) * qsat(1:ncol,1:pver) + epsilo) / &
         ((1.0 - epsilo) * state%q(1:ncol,1:pver,1) + epsilo)

      if (radforce) then

         pmxrgnrf = pmxrgn
         nmxrgnrf = nmxrgn

         call get_rf_scales(scales)

         call get_aerosol(lchnk, pint, aerosol, scales)

         ! overwrite with prognostics aerosols
         call set_aerosol_from_prognostics (state, aerosol)

         call aerosol_indirect(ncol,lchnk,landfrac,pmid,t,qm1,cld,zm,rel)
   
         call t_startf('radcswmx_rf')

         call radcswmx(lchnk   ,ncol ,                            &
                    pnm     ,pbr     ,qm1     ,rh      ,o3mmr   , &
                    aerosol ,cld     ,cicewp  ,cliqwp  ,rel     , &
                    rei     ,eccf    ,coszrs  ,scon    ,solin   , &
                    asdir   ,asdif   ,aldir   ,aldif   ,nmxrgnrf, &
                    pmxrgnrf,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsnirt  ,fsnrtc  ,fsnirtsq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   ,frc_day ,                   &
                    aertau  ,aerssa  ,aerasm  ,aerfwd  ,fns     , &
                    fcns    )

         call t_stopf('radcswmx_rf')

!
! Convert units of shortwave fields needed by rest of model from CGS to MKS
!

            do i = 1, ncol
            solin(i) = solin(i)*1.e-3
            fsnt(i)  = fsnt(i) *1.e-3
            fsns(i)  = fsns(i) *1.e-3
            fsntc(i) = fsntc(i)*1.e-3
            fsnsc(i) = fsnsc(i)*1.e-3
         end do
         ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair

!
! Dump shortwave radiation information to history tape buffer (diagnostics)
!
         call outfld('QRS_RF  ',ftem  ,pcols,lchnk)
         call outfld('FSNT_RF ',fsnt  ,pcols,lchnk)
         call outfld('FSNS_RF ',fsns  ,pcols,lchnk)
         call outfld('FSNTC_RF',fsntc ,pcols,lchnk)
         call outfld('FSNSC_RF',fsnsc ,pcols,lchnk)
 
      endif ! if (radforce)

      call get_int_scales(scales)
 
      call get_aerosol(lchnk, pint, aerosol, scales)

      ! overwrite with prognostics aerosols
      call set_aerosol_from_prognostics (state, aerosol)

      call aerosol_indirect(ncol,lchnk,landfrac,pmid,t,qm1,cld,zm,rel)

      call t_startf('radcswmx')

      call radcswmx(lchnk   ,ncol    ,                            &
                    pnm     ,pbr     ,qm1     ,rh      ,o3mmr   , &
                    aerosol ,cld     ,cicewp  ,cliqwp  ,rel     , &
                    rei     ,eccf    ,coszrs  ,scon    ,solin   , &
                    asdir   ,asdif   ,aldir   ,aldif   ,nmxrgn  , &
                    pmxrgn  ,qrs     ,fsnt    ,fsntc   ,fsntoa  , &
                    fsntoac ,fsnirt  ,fsnrtc  ,fsnirtsq,fsns    , &
                    fsnsc   ,fsdsc   ,fsds    ,sols    ,soll    , &
                    solsd   ,solld   ,frc_day ,                   &
                    aertau  ,aerssa  ,aerasm  ,aerfwd  ,fns     , &
                    fcns)

      call t_stopf('radcswmx')


! -- tls ---------------------------------------------------------------2

!  Output net fluxes at 200 mb

      call vertinterp(ncol, pcols, pverp, pint, 20000._r8, fcns, fsn200c)
      call vertinterp(ncol, pcols, pverp, pint, 20000._r8, fns, fsn200)

!
! Convert units of shortwave fields needed by rest of model from CGS to MKS
!
      do i=1,ncol
         solin(i) = solin(i)*1.e-3
         fsds(i)  = fsds(i)*1.e-3
         fsnirt(i)= fsnirt(i)*1.e-3
         fsnrtc(i)= fsnrtc(i)*1.e-3
         fsnirtsq(i)= fsnirtsq(i)*1.e-3
         fsnt(i)  = fsnt(i) *1.e-3
         fsns(i)  = fsns(i) *1.e-3
         fsntc(i) = fsntc(i)*1.e-3
         fsnsc(i) = fsnsc(i)*1.e-3
         fsdsc(i) = fsdsc(i)*1.e-3
         fsntoa(i)=fsntoa(i)*1.e-3
         fsntoac(i)=fsntoac(i)*1.e-3
         fsn200(i)  = fsn200(i)*1.e-3
         fsn200c(i) = fsn200c(i)*1.e-3
      end do
      ftem(:ncol,:pver) = qrs(:ncol,:pver)/cpair
!
! Dump shortwave radiation information to history tape buffer (diagnostics)
!

      call outfld('frc_day ', frc_day, pcols, lchnk)
      call outfld('SULOD_v ', aertau(:,idxVIS,1) ,pcols,lchnk)
      call outfld('SSLTOD_v', aertau(:,idxVIS,2) ,pcols,lchnk)
      call outfld('CAROD_v ', aertau(:,idxVIS,3) ,pcols,lchnk)
      call outfld('DUSTOD_v', aertau(:,idxVIS,4) ,pcols,lchnk)
      call outfld('BGOD_v  ', aertau(:,idxVIS,5) ,pcols,lchnk)
      call outfld('VOLCOD_v', aertau(:,idxVIS,6) ,pcols,lchnk)
      call outfld('AEROD_v ', aertau(:,idxVIS,7) ,pcols,lchnk)
      call outfld('AERSSA_v', aerssa(:,idxVIS,7) ,pcols,lchnk)
      call outfld('AERASM_v', aerasm(:,idxVIS,7) ,pcols,lchnk)
      call outfld('AERFWD_v', aerfwd(:,idxVIS,7) ,pcols,lchnk)
      call aerosol_diagnostics (state, aerosol)

      call outfld('QRS     ',ftem  ,pcols,lchnk)
      call outfld('SOLIN   ',solin ,pcols,lchnk)
      call outfld('FSDS    ',fsds  ,pcols,lchnk)
      call outfld('FSNIRTOA',fsnirt,pcols,lchnk)
      call outfld('FSNRTOAC',fsnrtc,pcols,lchnk)
      call outfld('FSNRTOAS',fsnirtsq,pcols,lchnk)
      call outfld('FSNT    ',fsnt  ,pcols,lchnk)
      call outfld('FSNS    ',fsns  ,pcols,lchnk)
      call outfld('FSNTC   ',fsntc ,pcols,lchnk)
      call outfld('FSNSC   ',fsnsc ,pcols,lchnk)
      call outfld('FSDSC   ',fsdsc ,pcols,lchnk)
      call outfld('FSNTOA  ',fsntoa,pcols,lchnk)
      call outfld('FSNTOAC ',fsntoac,pcols,lchnk)
      call outfld('SOLS    ',sols  ,pcols,lchnk)
      call outfld('SOLL    ',soll  ,pcols,lchnk)
      call outfld('SOLSD   ',solsd ,pcols,lchnk)
      call outfld('SOLLD   ',solld ,pcols,lchnk)
      call outfld('FSN200  ',fsn200,pcols,lchnk)
      call outfld('FSN200C ',fsn200c,pcols,lchnk)

   end if
!
! Longwave radiation computation
!
   if (dolw) then
!
! Convert upward longwave flux units to CGS
!
      do i=1,ncol
         lwupcgs(i) = lwup(i)*1000.
      end do
!
! Do longwave computation. If not implementing greenhouse gas code then
! first specify trace gas mixing ratios. If greenhouse gas code then:
!  o ixtrcg   => indx of advected n2o tracer
!  o ixtrcg+1 => indx of advected ch4 tracer
!  o ixtrcg+2 => indx of advected cfc11 tracer
!  o ixtrcg+3 => indx of advected cfc12 tracer
!
      if (trace_gas) then
         call cnst_get_ind('N2O'  , in2o)
         call cnst_get_ind('CH4'  , ich4)
         call cnst_get_ind('CFC11', if11)
         call cnst_get_ind('CFC12', if12)
         call t_startf("radclwmx")
         call radclwmx(lchnk   ,ncol    ,                            &
                       lwupcgs ,t       ,qm1(1,1,1)       ,o3vmr ,   &
                       pbr     ,pnm     ,pmln    ,piln    ,          &
                       qm1(1,1,in2o)    ,qm1(1,1,ich4)    ,          &
                       qm1(1,1,if11)    ,qm1(1,1,if12)    ,          &
                       cld     ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                       flns    ,flnt    ,flnsc   ,flntc   ,flwds   , &
                       flut    ,flutc   ,                            &
                       aerosol(:,:,idxVOLC)      ,fnl     , fcnl   )

         call t_stopf("radclwmx")
      else
         call trcmix(lchnk   ,ncol    , &
                     pmid    ,n2o     ,ch4     ,                     &
                     cfc11   ,cfc12   )

         call t_startf("radclwmx")
         call radclwmx(lchnk     ,ncol    ,                            &
                       lwupcgs   ,t       ,qm1(1,1,1)       ,o3vmr ,   &
                       pbr       ,pnm     ,pmln    ,piln    ,          &
                       n2o       ,ch4     ,cfc11   ,cfc12   ,          &
                       cld       ,emis    ,pmxrgn  ,nmxrgn  ,qrl     , &
                       flns      ,flnt    ,flnsc   ,flntc   ,flwds   , &
                       flut      ,flutc   ,                            &
                       aerosol(:,:,idxVOLC)        ,fnl     ,fcnl    )
         call t_stopf("radclwmx")
      endif
!
!  Output fluxes at 200 mb
!
   call vertinterp(ncol, pcols, pverp, pint, 20000._r8, fnl, fln200)
   call vertinterp(ncol, pcols, pverp, pint, 20000._r8, fcnl, fln200c)
!
! Convert units of longwave fields needed by rest of model from CGS to MKS
!
      do i=1,ncol
         flnt(i)  = flnt(i)*1.e-3
         flut(i)  = flut(i)*1.e-3
         flutc(i) = flutc(i)*1.e-3
         flns(i)  = flns(i)*1.e-3
         flntc(i) = flntc(i)*1.e-3
         fln200(i)  = fln200(i)*1.e-3
         fln200c(i) = fln200c(i)*1.e-3
         flnsc(i) = flnsc(i)*1.e-3
         flwds(i) = flwds(i)*1.e-3
         lwcf(i)=flutc(i) - flut(i)
         swcf(i)=fsntoa(i) - fsntoac(i)
      end do
!
! Dump longwave radiation information to history tape buffer (diagnostics)
!
      call outfld('QRL     ',qrl(:ncol,:)/cpair,ncol,lchnk)
      call outfld('FLNT    ',flnt  ,pcols,lchnk)
      call outfld('FLUT    ',flut  ,pcols,lchnk)
      call outfld('FLUTC   ',flutc ,pcols,lchnk)
      call outfld('FLNTC   ',flntc ,pcols,lchnk)
      call outfld('FLNS    ',flns  ,pcols,lchnk)
      call outfld('FLNSC   ',flnsc ,pcols,lchnk)
      call outfld('LWCF    ',lwcf  ,pcols,lchnk)
      call outfld('SWCF    ',swcf  ,pcols,lchnk)
      call outfld('FLN200  ',fln200,pcols,lchnk)
      call outfld('FLN200C ',fln200c,pcols,lchnk)
!
   end if
!
   return
end subroutine radctl
