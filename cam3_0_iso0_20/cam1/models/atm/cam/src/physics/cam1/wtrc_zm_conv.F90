#include <misc.h>
#include <params.h>

module wtrc_zm_conv

!-----------------------------------------------------------------------
!
! Purpose: Impliments water [isotope] tracers for the ZM convection, 
!     including evaporation of convective prcipitation. 
!
!  (wtrc_zm_conv_evap should be pasted to this module when happy)
!
! Notes:
!
!   One big whammy here is that this calulations is done AFTER the mass 
!   fluxes, etc, are converted to pressure (rather then geometric height) 
!   conordinates. Thus you will notice some of the signs are the oposite 
!   between here and cldprp (there is a factor of dz/dp missing). 
!   This is not very satisfying whn trying to compare codes, however, 
!   really can't avaoid it "as is", as we need to have the cu, evp, scaled
!   by mb also which occurs AFTER cldprp in the zm_convr routine.
!    
!
! Author: 
!
!   David Noone <dcn@caltech.edu> - Tue Oct 28 11:18:54 MST 2003
!
!-----------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: pcols, pver, pverp
   use constituents,    only: pcnst, pnats, qmin
   use water_tracers,   only: trace_water, lwtrczmlin, iwspec, wtrc_is_vap, ixwti,ixwtx,ixh2oq, &
                              wtrc_qchk1, wtrc_qchk2
   use moistconvection, only: grav,rgrav,limcnv
   use history,         only: outfld, addfld, add_default, phys_decomp
   use abortutils,      only: endrun

   use wtrc_camtune,    only: feq_upd, feq_dnd, feq_zme, feq_hke, ftc_dnd, fsnk_zmc, &
                              fcloud, lupdcorr, ltndcorr


!-----------------------------------------------------------------------
   implicit none
   save 
   private
!-----------------------------------------------------------------------

!
! PUBLIC interfaces
!
   public wtrc_zm_convi          ! constants and history fields
   public wtrc_zm_convr          ! convective scheme
   public wtrc_zm_conv_evap      ! post convective evaporation
!
! Private data
!
   real(r8), parameter :: msmall = 1.e-20  ! small flux for denominator (cldprp)
!   real(r8), parameter :: msmall = 1.e-16  ! bigger
!   real(r8), parameter :: msmall = 1.e-10  ! even bigger
!   real(r8), parameter :: qtiny  = 1.e-20  ! tinq q for irrelevant denominator
   real(r8), parameter :: qtiny  = 1.e-18  ! tinq q for irrelevant denominator
                                           ! (also trigger bypass calc mu/md is small)
   real(r8), parameter :: ftiny = 1.e-6    ! maximum rain out factor
!
!=======================================================================
CONTAINS

!=======================================================================
subroutine wtrc_zm_convi()
!-----------------------------------------------------------------------
!
! Purpose: Initialize model module data and output fields
!
!-----------------------------------------------------------------------
  use constituents,    only: cnst_name
!------------------------- Local Variables -----------------------------
  character(len=2) :: atr             ! label for output variable name
  integer m                             ! tracer index
!-----------------------------------------------------------------------
  atr= '  '             ! this line just for compiler check
!
! History field definitions 
!
   do m = ixwti, ixwtx
     if (wtrc_is_vap(m)) then
       write(atr,'(i2.2)') m

! ZM_convr (note DLF and DIF here differ from that in tphys which is hack+zm)
       call addfld('ZMUPD'//atr,'kg/kg ',pver, 'A',  &
          trim(cnst_name(m))//' vapour in updraft',phys_decomp)
       call addfld('ZMLIQ'//atr,'kg/kg ',pver, 'A',  &
          trim(cnst_name(m))//' liquid in updraft',phys_decomp)
       call addfld('ZMICE'//atr,'kg/kg ',pver, 'A',  &
          trim(cnst_name(m))//' ice in updraft',phys_decomp)
       call addfld('ZMDND'//atr,'kg/kg ',pver, 'A',  &
          trim(cnst_name(m))//' vapour in downdraft',phys_decomp)
       call addfld('ZMDLF'//atr,'kg/kg/s ',pver, 'A',  &
          trim(cnst_name(m))//' liquid detrain from updraft',phys_decomp)
       call addfld('ZMDIF'//atr,'kg/kg/s ',pver, 'A',  &
          trim(cnst_name(m))//' ice detrain from updraft',phys_decomp)

       call add_default('ZMUPD'//atr, 1,' ')
       call add_default('ZMLIQ'//atr, 1,' ')
       call add_default('ZMICE'//atr, 1,' ')
       call add_default('ZMDND'//atr, 1,' ')
       call add_default('ZMDLF'//atr, 1,' ')
       call add_default('ZMDIF'//atr, 1,' ')
     
! ZM_evap
       call addfld('ZMVIEV'//atr,'kg/m2/s ',  1, 'A',  &
          trim(cnst_name(m))//' integrated evaporated precipitation from ZM convection',phys_decomp)

       call addfld('ZMFXPR'//atr,'kg/m2/s ',pverp, 'A',  &
          trim(cnst_name(m))//' flux of precipitation from ZM convection',phys_decomp)
       call addfld('ZMFXSN'//atr,'kg/m2/s ',pverp, 'A',  &
          trim(cnst_name(m))//' flux of snow from ZM convection',phys_decomp)
       call addfld('ZMNPPD'//atr,'kg/m2/s ',pver , 'A',  &
          trim(cnst_name(m))//' net precipitation production from ZM convection',phys_decomp)
       call addfld('ZMNSPD'//atr,'kg/m2/s ',pver , 'A',  &
          trim(cnst_name(m))//' net snow production from ZM convection',phys_decomp)

       call add_default('ZMVIEV'//atr, 1,' ')
!!       call add_default('ZMFXPR'//atr, 1,' ')
!!       call add_default('ZMFXSN'//atr, 1,' ')
!!       call add_default('ZMNPPD'//atr, 1,' ')
!!       call add_default('ZMNSPD'//atr, 1,' ')
       
! HK_evap
       call addfld('HKVIEV'//atr,'kg/m2/s ',  1, 'A',  &
          trim(cnst_name(m))//' integrated evaporated precipitation from ZM convection',phys_decomp)

       call addfld('HKFXPR'//atr,'kg/m2/s ',pverp, 'A',  &
          trim(cnst_name(m))//' flux of precipitation from HK convection',phys_decomp)
       call addfld('HKFXSN'//atr,'kg/m2/s ',pverp, 'A',  &
          trim(cnst_name(m))//' flux of snow from HK convection',phys_decomp)
       call addfld('HKNPPD'//atr,'kg/m2/s ',pver , 'A',  &
          trim(cnst_name(m))//' net precipitation production from HK convection',phys_decomp)
       call addfld('HKNSPD'//atr,'kg/m2/s ',pver , 'A',  &
          trim(cnst_name(m))//' net snow production from HK convection',phys_decomp)

       call add_default('HKVIEV'//atr, 1,' ')
!!       call add_default('HKFXPR'//atr, 1,' ')
!!       call add_default('HKFXSN'//atr, 1,' ')
!!       call add_default('HKNPPD'//atr, 1,' ')
!!       call add_default('HKNSPD'//atr, 1,' ')
     end if
   end do
!
  return
end subroutine wtrc_zm_convi
 

!=======================================================================
subroutine wtrc_zm_convr(lchnk   ,ncol    , &
                    t       ,qh      ,prec    , &
                    qtnd    ,tu      ,td      , &
                    dpp     , &
                    delt    ,cme     ,          &
                    dlf     ,dif     ,pflx    ,rprd    , &
                    mu      ,md      ,du      ,eu      ,ed      , &
                    dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                    lengath ,ql      ,qi      ,qu      ,qd      ,cu, evp , &
                    zi      ,geos    ,rr      ,jd      ,ficed   ,ficeup)

!-----------------------------------------------------------------------
!
! Purpose: impliments water [isotope] tracers for the ZM convection.
!
! Method:
!
!   For both the updrafts and downdrafts, solve the budget equations
!   as in the ZM routines cldprp (convtran) and q1q2, except here also allow
!   isotopic fractionation. Reevaporation is handled externally.
!
!
! Notes:
!
!   It is important to have mass conservation for the column.
!   It is important to get PFLX right for output precipitation.
!
!   This code is based on that originally put into CCM3, and essentially
!   solves a diffusion equation using a forward difference (trouble!).
!   Because of this, and physicsally, we the entrained environment may
!   have already interacted with falling precip, we should iterate.
!
!   Code modified to track ZM sceme from cam2.x
!   All cloud "dynamics" are input on gathered arrays.
!   All cloud "properties" are input on scattered arrays.
!
!
! Author:
!
!   David Noone <dcn@caltech.edu> - Wed Aug 08 10:25:17 PDT 2001 (ccm3 masses)
!   David Noone <dcn@caltech.edu> - Fri Sep 27 14:58:35 PDT 2002 (ccm3 ratios)
!   David Noone <dcn@caltech.edu> - Mon Oct 20 16:44:01 PDT 2003 (cam2.x)
!
!----------------------------------------------------------------------
   use water_tracers,  only: wtrc_ratio
!-----------------------------------------------------------------------
   implicit none
!---------------------------- Arguments --------------------------------
! G: gathers, S: scattered
!
   integer , intent(in)    :: lchnk                         ! chunk identifier
   integer , intent(in)    :: ncol                          ! number of atmospheric columns
   integer , intent(in)    :: lengath                       ! number of gathered points
   integer , intent(in)    :: ideep(pcols)                  ! S: indicies of gathered points 
   integer , intent(in)    :: maxg(pcols)                   ! G: index of convection base
   integer , intent(in)    :: jt(pcols)                     ! G: index of convection top
   integer , intent(in)    :: jd(pcols)                     ! G: index of downdraft top

   real(r8), intent(in)    :: delt                          ! time step (s)
   real(r8), intent(in)    :: t(pcols,pver)                 ! S: temperature
   real(r8), intent(in)    :: qh(pcols,pver,pcnst+pnats)    ! S: specific humidity.
   real(r8), intent(in)    :: dpp(pcols,pver)               ! S: level thichness (Pa)
!
   real(r8), intent(in)    :: tu(pcols,pver)                ! S: updraft temperature
   real(r8), intent(in)    :: td(pcols,pver)                ! S: downdraft temperature

   real(r8), intent(in)    :: mu(pcols,pver)                ! G: updraft mass flax
   real(r8), intent(in)    :: eu(pcols,pver)                ! G: updraft entrainment rate
   real(r8), intent(in)    :: du(pcols,pver)                ! G: updraft detrainment rate
   real(r8), intent(in)    :: md(pcols,pver)                ! G: downdraft masss flux
   real(r8), intent(in)    :: ed(pcols,pver)                ! G: downdraft entrainment rate
   real(r8), intent(in)    :: dp(pcols,pver)                ! G: level thichnesss (mb)
   real(r8), intent(in)    :: dsubcld(pcols)                ! G: sub-cloud thichness ( mb)
   real(r8), intent(in)    :: zi(pcols,pver+1)              ! height of interfaces
   real(r8), intent(in)    :: geos(pcols)                   ! surface geopotential

   real(r8), intent(in)    :: ficed(pcols,pver)             ! G: ice fraction detrained
   real(r8), intent(in)    :: ficeup(pcols,pver)            ! G: ice fraction in updraft
!
   real(r8), intent(in)    :: cu(pcols,pver)                ! S: updraft condensation
   real(r8), intent(in)    :: rr(pcols,pver)                ! S: updraft rain rate
   real(r8), intent(in)    :: evp(pcols,pver)               ! S: downdraft evaporation
   real(r8), intent(inout) :: qu(pcols,pver,pcnst+pnats)    ! S: updraft water vapour
   real(r8), intent(inout) :: qd(pcols,pver,pcnst+pnats)    ! S: downdraft water vapour
   real(r8), intent(inout) :: ql(pcols,pver,pcnst+pnats)    ! S: updraft liquid water
   real(r8), intent(inout) :: qi(pcols,pver,pcnst+pnats)    ! S: updraft ice water
   real(r8), intent(inout) :: cme(pcols,pver,pcnst+pnats)   ! S: dqdt from C-E
   real(r8), intent(inout) :: rprd(pcols,pver,pcnst+pnats)  ! S: precipitation production
   real(r8), intent(inout) :: pflx(pcols,pverp,pcnst+pnats) ! S: precip flux at each level
   real(r8), intent(inout) :: prec(pcols,pcnst+pnats)       ! S: surface precipitataion flux

   real(r8), intent(inout) :: dlf(pcols,pver,pcnst+pnats)   ! S: detraining cld liquid h2o tend
   real(r8), intent(inout) :: dif(pcols,pver,pcnst+pnats)   ! S: detraining cld ice h2o tend
   real(r8), intent(inout) :: qtnd(pcols,pver,pcnst+pnats)  ! S: specific humidity tendency (kg/kg/s)

   character(len=2) :: atr	 	   ! label for tracer history fields

!------------------------- Local Variabless ----------------------------
   integer  :: i,k,m                          ! array indicies
   integer  :: msg                            ! upper level bound
!
   real(r8) :: tg(pcols,pver)                 ! G: temperature
   real(r8) :: tug(pcols,pver)                ! G: temperature in updraft 
   real(r8) :: tdg(pcols,pver)                ! G: temperature in downdraft
!
   real(r8) :: zf(pcols,pver+1)               ! S: geopotential heigh of interfaces
   real(r8) :: zfg(pcols,pver+1)              ! G: gathered zf
   real(r8) :: dzdp(pcols,pver)               ! G: height pressure conversion
!
   real(r8) :: q(pcols,pver,pcnst+pnats)      ! S: local copy of shums
   real(r8) :: qg(pcols,pver,pcnst+pnats)     ! G: shum 'qbar'
   real(r8) :: qhat(pcols,pver,pcnst+pnats)   ! G: interface shum
   real(r8) :: qhup(pcols,pver,pcnst+pnats)   ! G: interface shum FLUX
   real(r8) :: qhdn(pcols,pver,pcnst+pnats)   ! G: interface shum FLUX
   real(r8) :: qupd(pcols,pver,pcnst+pnats)   ! G: updraft vapour FLUX
   real(r8) :: qliq(pcols,pver,pcnst+pnats)   ! G: updraft liquid  FLUX
   real(r8) :: qice(pcols,pver,pcnst+pnats)   ! G: updraft ice  FLUX
   real(r8) :: qdnd(pcols,pver,pcnst+pnats)   ! G: downdraft vapour FLUX
   real(r8) :: cug(pcols,pver,pcnst+pnats)    ! G: condensation 
   real(r8) :: rrg(pcols,pver,pcnst+pnats)    ! G: rain rate
   real(r8) :: evpg(pcols,pver,pcnst+pnats)   ! G: evaporation rate

   real(r8) :: ldetg(pcols,pver,pcnst+pnats)  ! G: detraining liquid isotope mass
   real(r8) :: idetg(pcols,pver,pcnst+pnats)  ! G: detraining ice isotope mass
   real(r8) :: pflxg(pcols,pverp,pcnst+pnats) ! G: precip flux at each interface
   real(r8) :: rprdg(pcols,pver,pcnst+pnats)  ! G: precip production at each level

   real(r8) :: dqemc(pcols,pver,pcnst+pnats)  ! G: evap minus condensation tend
   real(r8) :: dqmdt(pcols,pver,pcnst+pnats)  ! G: final convective tendancy
   real(r8) :: dlmdt(pcols,pver,pcnst+pnats)  ! G: detraining cld liquid h2o tend
   real(r8) :: dimdt(pcols,pver,pcnst+pnats)  ! G: detraining cld ice h2o tend

   real(r8) :: qdifr			      ! difference for interpolation

!-----------------------------------------------------------------------
!!   write(*,*) '(wtrc_zm_convr) start'
!
! Set internal variable "msg" (convection limit) to "limcnv-1"
!
   msg = limcnv - 1
!
! initialize necessary arrays. 
! These are all local variables - intent(in)
!
   zfg(:ncol,:)    = 0.
   dzdp(:ncol,:)   = 0.
!
   q    (:ncol,:,:) = 0.
   qg   (:ncol,:,:) = 0.
   qhat (:ncol,:,:) = 0.
   qhup (:ncol,:,:) = 0.
   qhdn (:ncol,:,:) = 0.
   qupd (:ncol,:,:) = 0.
   qliq (:ncol,:,:) = 0.
   qice (:ncol,:,:) = 0.
   qdnd (:ncol,:,:) = 0.
   cug  (:ncol,:,:) = 0.
   rrg  (:ncol,:,:) = 0.
   evpg (:ncol,:,:) = 0.
   rprdg(:ncol,:,:) = 0. 
   pflxg(:ncol,:,:) = 0. 
   ldetg(:ncol,:,:) = 0.
   idetg(:ncol,:,:) = 0.

   dqemc(:ncol,:,:) = 0. 
   dqmdt(:ncol,:,:) = 0. 
   dlmdt(:ncol,:,:) = 0.
   dimdt(:ncol,:,:) = 0.
!
! Zero only tracer arrays (all species), intent(inout), m=1 (in), m=tracer (out)
! (could also be dimensioned ixwti:ixwtx)
!
   qu  (:ncol,:,2:) = 0.
   qd  (:ncol,:,2:) = 0.
   ql  (:ncol,:,2:) = 0.
   cme (:ncol,:,2:) = 0.
   rprd(:ncol,:,2:) = 0.
   pflx(:ncol,:,2:) = 0.
   prec(:ncol  ,2:) = 0.

   dlf (:ncol,:,2:) = 0. 
   dif (:ncol,:,2:) = 0. 
   qtnd(:ncol,:,2:) = 0.
!
! Compute height of levels
!
   do k = 1, pver+1
     zf(:ncol,k) = zi(:ncol,k) + grav*geos(:ncol)
   end do

! Partition incomming liquid to ice and liquid - for output
   
  do k = 1, pver
    do i = 1, lengath
       qi(ideep(i),k,1)= (   ficeup(i,k))*ql(ideep(i),k,1)
       ql(ideep(i),k,1)= (1.-ficeup(i,k))*ql(ideep(i),k,1)
    end do
  end do
!
! Preprocess for all water vapour tracerso
! For the sake of completeness, also do totals
!
   do m = 1, ixwtx
     if (m==1 .or. wtrc_is_vap(m)) then
!
! store incoming specific humidity field for subsequent calculation
! of precipitation (through change in storage).
! define dry static energy (normalized by cp).
!
       do k = 1,pver
          do i = 1,ncol
             q(i,k,m)    = qh(i,k,m)
             qhat(i,k,m) = q(i,k,m)             ! ideep?
          end do
       end do
!
! obtain gathered arrays necessary for ensuing calculations.
! Also regather prognostic water inputs to gathered arrays.
!
       do k = 1,pver
          do i = 1,lengath
             qg(i,k,m) = q(ideep(i),k,m)
          end do
       end do
!
! Compute interface (grid average) specific humidity.
! Optionally use linear rather than log interpolation.
!
       do k = msg + 2,pver
         do i = 1,lengath
           if (lwtrczmlin) then		! linear interpolation for midpoints
             qhat(i,k,m) = 0.5* (qg(i,k,m)+qg(i,k-1,m))
           else				! log. interpolation for midpoints
              qdifr = 0.
              if (qg(i,k,m) > 0. .or. qg(i,k-1,m) > 0.) &
                 qdifr = abs((qg(i,k,m)-qg(i,k-1,m))/max(qg(i,k-1,m),qg(i,k,m)))
              if (qdifr > 1.E-6) then
                 qhat(i,k,m) = log(qg(i,k-1,m)/qg(i,k,m))*qg(i,k-1,m)*qg(i,k,m)/(qg(i,k-1,m)-qg(i,k,m))
              else
                 qhat(i,k,m) = 0.5* (qg(i,k,m)+qg(i,k-1,m))
              end if
           end if
         end do
       end do
!
    end if      ! m is water vapour tracer
  end do        ! m tracer (and m=1) loop
!
! Gather non-tracer quantities
!
  do k = 1,pver
    do i = 1,lengath
       tg   (i,k)   = t   (ideep(i),k)
!!       tug  (i,k)   = 0.5*(t(ideep(i),k) + t(ideep(i),min(k+1,pver)))
!!       tdg  (i,k)   = 0.5*(t(ideep(i),max(k-1,1)) + t(ideep(i),k))
       tug  (i,k)   = tu  (ideep(i),k)		! about 1K warmer
       tdg  (i,k)   = td  (ideep(i),k)		! about 2K cooler

       qliq (i,k,1) = ql  (ideep(i),k,1)
       qice (i,k,1) = qi  (ideep(i),k,1)
       qupd (i,k,1) = qu  (ideep(i),k,1)
       qdnd (i,k,1) = qd  (ideep(i),k,1)
!
       cug  (i,k,1) = cu  (ideep(i),k)
       rrg  (i,k,1) = rr  (ideep(i),k)
       evpg (i,k,1) = evp (ideep(i),k)
       rprdg(i,k,1) = rprd(ideep(i),k,1)

       dqmdt(i,k,1) = qtnd(ideep(i),k,1) 
       dqemc(i,k,1) = -cme(ideep(i),k,1)        ! minus for q1q2/cldprp compadability
       dlmdt(i,k,1) = dlf (ideep(i),k,1)
       dimdt(i,k,1) = dif (ideep(i),k,1)

    end do
  end do
!
  do k = 1, pver+1
    do i = 1, lengath
       zfg (i,k)   = zf  (ideep(i),k)
       pflxg(i,k,1) = pflx(ideep(i),k,1)
    end do
  end do
!
! Compure height to pressure conversion
!
  do k = 1, pver
    do i = 1, lengath
      dzdp(i,k) = (zfg(i,k)-zfg(i,k+1))/dp(i,k)
    end do
  end do
!
! Compute incloud (updraft and downdraft) propertures 
! [This routine does most of the real work - everything else is housekeeping]
!
     call wtrc_zm_convr_cldprp ( lchnk   ,1       ,lengath ,  &
               msg     ,jt      ,maxg    ,dp      ,dsubcld ,  &
               mu      ,du      ,eu      ,md      ,ed      ,  &
               tg      ,tug     ,qg      ,qhat    ,qhup    ,qhdn,     &
               qupd    ,  &
               qliq    ,qice    ,qdnd    ,cug     ,evpg    ,rprdg   ,  &
               dzdp    ,pflxg   ,tdg     ,jd      ,rrg     ,  &
               ldetg   ,idetg   ,ficed   ,ficeup )
!
! Compute tendancies
!
     call wtrc_zm_convr_q1q2 ( lchnk   ,1       ,lengath ,  &
             msg     ,jt      ,maxg    ,dp      ,dsubcld ,  &
             mu      ,du      ,md      ,qg      ,  &
             qhup    ,qhdn    ,qupd    ,qliq    ,qice   , qdnd    , &
             cug     ,evpg    ,dqmdt   ,dqemc   ,dlmdt   ,  &
             dimdt   ,ldetg   ,idetg )
!
! Finalize all tracer quantities for output
!
   do m = ixwti,ixwtx
     if (wtrc_is_vap(m)) then
!
! Gather back to assign to output.
!  q (local) is updated to compute net precip)
! (use ratios simply to avoid dividimg by potentially small mass fluxes)
!
       do k = msg + 1,pver
          do i = 1,lengath
             q   (ideep(i),k,m) = qh(ideep(i),k,m) + 2.*delt*dqmdt(i,k,m)
             qtnd(ideep(i),k,m) = dqmdt(i,k,m)
             cme (ideep(i),k,m) = -dqemc(i,k,m)		! check sign on output
             dlf (ideep(i),k,m) = dlmdt(i,k,m)
             dif (ideep(i),k,m) = dimdt(i,k,m)
             ql  (ideep(i),k,m) = ql(ideep(i),k,1)*wtrc_ratio(qliq (i,k,m), qliq (i,k,1))
             qi  (ideep(i),k,m) = qi(ideep(i),k,1)*wtrc_ratio(qice (i,k,m), qice (i,k,1))
             qu  (ideep(i),k,m) = qu(ideep(i),k,1)*wtrc_ratio(qupd (i,k,m), qupd (i,k,1))
             qd  (ideep(i),k,m) = qd(ideep(i),k,1)*wtrc_ratio(qdnd (i,k,m), qdnd (i,k,1))
             rprd(ideep(i),k,m) = rprdg(i,k,m)
             pflx(ideep(i),k,m) = pflxg(i,k,m)
          end do
       end do
!
       do i = 1,lengath
          pflx(ideep(i),pverp,m) = pflxg(i,pverp,m)
       end do
!
! Compute precip by integrating change in water vapor minus detrained cloud water
       do k = pver,msg + 1,-1
          do i = 1,ncol
             prec(i,m) = prec(i,m) - &
                 dpp(i,k)*(q(i,k,m)-qh(i,k,m)) - dpp(i,k)*(dlf(i,k,m)+dif(i,k,m))*2*delt
          end do
       end do

! obtain final precipitation rate in m/s.
       do i = 1,ncol
          prec(i,m) = rgrav*max(prec(i,m),0._r8) / (2.*delt)/1000.
       end do
!
! Write diagnostics to history file
!
       write(atr,'(i2.2)') m
       call outfld('ZMUPD'//atr, qu(:,:,m)  , pcols, lchnk)
       call outfld('ZMLIQ'//atr, ql(:,:,m)  , pcols, lchnk)
       call outfld('ZMICE'//atr, qi(:,:,m)  , pcols, lchnk)
       call outfld('ZMDND'//atr, qd(:,:,m)  , pcols, lchnk)
       call outfld('ZMDLF'//atr, dlf(:,:,m) , pcols, lchnk)
       call outfld('ZMDIF'//atr, dif(:,:,m) , pcols, lchnk)
!
     end if     ! m is water vapour tracer
   end do       ! m tracer loop
!
   return
end subroutine wtrc_zm_convr

!=======================================================================

subroutine wtrc_zm_conv_evap(state  , ptend  , prdprec, cldfrc , deltat , & 
                             flxpr  , flxsn  , fevprec, fevsnow, & 
                             prec   , snow   , hackflag )
!-----------------------------------------------------------------------
! Purpose:
!
!   Tendencies for all isotopic species due to evaporation of precip. from
!   the ZM scheme. Only do for water tracers.
!
! Method:
!
!   Computation of isotope tendencies in parallel with the total water
!   tendencies stated in zm_conv_evap.
!   Sublimation of snow falling from layer above occurs without
!   fractionation, as the isotopic exchange is limited by high
!   molecular diffusivity in the ice. The evaporation of rain (and 
!   melted snow) falling from above assumes equilibrium with vapour at
!   that layer, and does NOT include precipitation formed AT the layer.
!   This is because we assume that evaporation occurs in the cloud free
!   part of the grid box. The precipitation produced at the present
!   layer will be allowed to evaporate in the layer below. 
!
!   The ratios are determined (or traced) according to the flux of the snow, 
!
!      Revapsnow = Rsnow = Fs/Ft
!
!   And equlibration of rain follows,  as
!
!      rain_i/rain_t = alpha (vapor_i/vapor_t)
!      totalwater = rainfallingfromabove + initialvapour
!
!   The fractionation factor uses a mean temperature, and is modified to
!   incorperate kinetic effects due to sub-saturation (Stewart style).
!   Since we have also the solid phase, we could additionally
!   incorperate some microphysics for the mixed cloud system.
!   Finally, the fractionation factor is further modified to incorperate
!   the fact that if drops are sufficiently large they will not be able
!   to reach equilibrium during the time it takes to pass through the
!   layer. Presently this is set as a constant value (in
!   water_isotope module), but could easily be modified to better
!   parametrize a mean drop size given some assumed drop distribution, 
!   and thus fall rate, etc.
!
!
! Notes:
!
!   3. We've assumed things mix fast in the precip. for simplification
!   4. When snow is produced, we assume the ratio is the same as the total 
!      precip. The snow formation is a fraction whose T dependence seems not
!      matching that for fractionation. (Need further work, perhaps??)     
!   5. Need different isotope rules for ice and liquid. As we have
!      both, we can account for the supersaturation explictly.
!   6. Does not account for possibility of refreeze during fall.
!   7. Does not account corectly for mixed phase microphysics (bergeron perocess).
!
!
! Change history:
!   o Results for "total" water passed in rather than recalculated (dcn 23/10/03)
!   o Removed species dependance so can have more than one tracer of
!     the same isotopic species (dcn 23/10/03)
!   o Removed "tendancy.name" flag, as we only do one physics_update for
!     all tracers (dcn 23/11/03)
!   1 Input prdprec is computed in ZM scheme and extended to all isotopes (dcn 24/11/03)
!   2 Output prec and snow NOW extended to include all isotopes (dcn 24/11/03)
!
!
! Author: Sun Wong <swong@atmos.umd.edu> - (9/3/2003)
!         David Noone <dcn@caltech.edu> - Thu Oct 23 17:05:38 PDT 2003
!
!-----------------------------------------------------------------------
!
   use shr_kind_mod,          only: r8 => shr_kind_r8
   use ppgrid,                only: pcols, pver, pverp
   use physics_types,         only: physics_state, physics_ptend
   use constituents,          only: pcnst, pnats
   use physconst,             only: tmelt, gravit, cpair
   use wv_saturation,         only: aqsat
   use history,               only: outfld, addfld, add_default
   use cldwat,                only: cldwat_fice
   use zm_conv,               only: ke
   use water_tracers,         only: iwspec, ixwti, ixwtx, wtrc_ratio, &
                                    wtrc_qchk1, ixh2oq
   use water_isotopes,        only: wiso_alpl, wiso_akel, wiso_heff, wiso_liqvap_equil
!
!-----------------------------------------------------------------------
   implicit none

!---------------------------- Arguments --------------------------------

   type(physics_state), intent(in)    :: state    ! Physics state variabiles
   type(physics_ptend), intent(inout) :: ptend    ! parameterized tendencies

   real(r8), intent(in)    :: prdprec(pcols,pver,pcnst+pnats) ! precip production
   real(r8), intent(in)    :: cldfrc(pcols,pver)  ! cloud fraction
   real(r8), intent(in)    :: deltat              ! time step

   real(r8), intent(inout) :: flxpr(pcols,pverp) ! interface flux of precipitataion
   real(r8), intent(inout) :: flxsn(pcols,pverp) ! interface flux of snow sublimated
   real(r8), intent(in)	   :: fevprec(pcols,pver) ! fraction of precip evaporated
   real(r8), intent(in)	   :: fevsnow(pcols,pver) ! fraction of snow sublimated
   logical,  intent(in)    :: hackflag            ! input is from Hack instead 
                                                  ! of ZM if true (for output only)

   real(r8), intent(inout) :: prec(pcols,pcnst+pnats)  ! Conv.-scale precip. rate
   real(r8), intent(out)   :: snow(pcols,pcnst+pnats)  ! Conv.-scale snowfall rate
!
!------------------------- Local Variables -----------------------------
!
   integer   :: i, k, m                    ! longitude, level and tracer indices
   integer   :: ncol, lchnk                ! no. of columns and chunk index
!
   real(r8)  :: hum0(pcols,pver)	  ! initial humidity 

   real(r8)  :: flxprec(pcols,pverp,pcnst+pnats) ! interface flux of precipitataion
   real(r8)  :: flxsnow(pcols,pverp,pcnst+pnats) ! interface flux of snow sublimated
   real(r8)  :: flxrain(pcols,pverp,pcnst+pnats) ! Conv. flx of rain from above (kg/m2/s)
   real(r8)  :: ntprprd(pcols,pver,pcnst+pnats)  ! net precip production (kg/kg/s)
   real(r8)  :: ntsnprd(pcols,pver,pcnst+pnats)  ! net snow production (kg/kg/s)
!
   real(r8)  :: evpvint(pcols,pcnst+pnats) ! vertical integral of evap.
   real(r8)  :: evpprec(pcols,pcnst+pnats) ! evap. of precip.     (kg/kg/s)
   real(r8)  :: evpsnow(pcols,pcnst+pnats) ! evap. of snow fall   (kg/kg/s)
   real(r8)  :: evprain(pcols,pcnst+pnats) ! evap. of rain fall   (kg/kg/s)
   real(r8)  :: snowmlt(pcols,pver,pcnst+pnats) ! snow melt tendency   (kg/kg/s) 
!
   real(r8)  :: flxsntm(pcols,pverp,pcnst+pnats) ! flux of snow NOT Melted
!
   real(r8)  :: vaptot(pcols,pver)	   ! total vapour
   real(r8)  :: liqtot(pcols,pver)	   ! total vapour
   real(r8)  :: vapiso		           ! isotope vapour
   real(r8)  :: liqiso                     ! isotope liquid
   real(r8)  :: dliqiso			   ! change in liquid from vapour
   real(r8)  :: fsnow_conv(pcols,pver)     ! snow fraction in precip production
   real(r8)  :: fice(pcols,pver)           ! ice fraction in precip production
   real(r8)  :: qsat(pcols,pver)           ! saturation specific humidity
   real(r8)  :: est(pcols,pver)            ! Saturation vapor pressure
   real(r8)  :: evplimit                   ! temp. variable for evap. limit 
   real(r8)  :: work1, work2               ! temporary variables
   real(r8)  :: fsnow			   ! fraction of precip flux that is snow
   real(r8)  :: tbar			   ! mean temperature for fractionation
   real(r8)  :: heff			   ! effective humidity for kinetic effects
   real(r8)  :: alpliq			   ! isotopic fractionation factor
   real(r8)  :: rat		           ! an isotope ratio
   real(r8)  :: qtot			   ! vapour plus liquid in equilibrium
   real(r8)  :: vnew			   ! "new" vapour in closed system
   real(r8)  :: fclr			   ! fraction of grid box with clear sky
   real(r8)  :: feq			   ! fraction rain equilibrated
!
   character(len=2) :: atr	 	   ! label for tracer history fields
!
!------------------------------------------------------------------------ 
!   real(r8)  :: fcloud = 1.0		   ! 1 - fraction of cloudy sky to include
!   real(r8)  :: fcloud = 0.5
!   real(r8)  :: fcloud = 0.0		   ! all sky
!-----------------------------------------------------------------------
!!   write(*,*) '(wtrc_zm_conv_evap) Start. Is hack=',hackflag
!
   ncol = state%ncol
   lchnk = state%lchnk
!
! Initilize and assign inputs to local variables
!
   flxprec(:ncol,:pverp,:) = 0.
   flxsnow(:ncol,:pverp,:) = 0.
   flxrain(:ncol,:pverp,:) = 0.
   flxsntm(:ncol,:pverp,:) = 0.
!
   flxprec(:ncol,:pverp,1) = flxpr(:ncol,:pverp)
   flxsnow(:ncol,:pverp,1) = flxsn(:ncol,:pverp)
!
   if (hackflag) then
      feq = feq_hke
   else
      feq = feq_zme
   endif
!
! This is needed to partition the output precipitation as ice or vapour.
! This SHOULD BE CHECKED FOR CONSISTENCY WITH CONVECTIVE FRACTIONATION.
!
  call cldwat_fice(ncol, state%t, fice, fsnow_conv)
!
! Determine saturation vapor pressure, and compute initial humidity
!
  call aqsat(state%t, state%pmid, est, qsat, pcols, ncol, pver, 1, pver)
  hum0(:ncol,:) = state%q(:ncol,:,1) / qsat(:ncol,:)
!
! Convert unit from m/s to kg/m2/s
!
   prec(:ncol,:) = prec(:ncol,:)*1000.
   vaptot(:ncol,:) = 0.
   liqtot(:ncol,:) = 0.
!
! Precompute total water rain flux for diagnostic checks, 
! Notice use of cloud factor means inconsistent with physics_update
!
   do k = 1, pver
     do i = 1, ncol
       if (state%t(i,k) > tmelt) then       ! melt all snow
         flxsntm(i,k,1) = 0.
         snowmlt(i,k,1) = flxsnow(i,k,1) * gravit/state%pdel(i,k)
       else                                 ! retain all snow
         flxsntm(i,k,1) = flxsnow(i,k,1)
         snowmlt(i,k,1) = 0.
       end if
       evpprec(i,1) = fevprec(i,k)*flxprec(i,k,1)*gravit/state%pdel(i,k)
       evpsnow(i,1) = fevsnow(i,k)*flxsntm(i,k,1)*gravit/state%pdel(i,k)
       evprain(i,1) = evpprec(i,1) - evpsnow(i,1)
       flxrain(i,k,1) = flxprec(i,k,1) - flxsnow(i,k,1)
!
       fclr = 1. - fcloud*cldfrc(i,k)
       qtot  = fclr*state%q(i,k,1) + &
                deltat*(evpsnow(i,1)+ snowmlt(i,k,1) + flxrain(i,k,1)*gravit/state%pdel(i,k))
!
! Save total vapour and liquid
!
       vaptot(i,k) =  fclr*state%q(i,k,1) + deltat*ptend%q(i,k,1)
       liqtot(i,k) = qtot - vaptot(i,k)
!  correct for numerics, conserving all mass
       if (liqtot(i,k) < 0.) then
          vaptot(i,k) = vaptot(i,k) - liqtot(i,k)
          liqtot(i,k) = 0.
       end if

     end do
   end do
!
! Loop over all isotope tracers
!
   do m = ixwti, ixwtx
     if (wtrc_is_vap(m)) then
       ptend%lq(m) = .TRUE.
!
! Zero tracer interface fluxes and midpoint integrals
!
       evpvint(:ncol,       m) = 0.
!
       do k = 1,pver			! top down integral
         do i = 1,ncol
!
! Melt snow falling into layer, if necessary.
!
           if (state%t(i,k) > tmelt) then	! melt all snow
             flxsntm(i,k,m) = 0.
             snowmlt(i,k,m) = flxsnow(i,k,m) * gravit/state%pdel(i,k)
           else					! retain all snow
             flxsntm(i,k,m) = flxsnow(i,k,m)
             snowmlt(i,k,m) = 0.
           end if
!
! Compute the evaporation of non-melted snow wthout isotopic fractionation
!
           evpsnow(i,m) = fevsnow(i,k)*flxsntm(i,k,m)*gravit/state%pdel(i,k)
!
! Assume isotopic equilibrium for rain and grid total vapour
! Do not add in precipitation produced at this layer. It will fall from
! the cloudy part, and treated in the layer below on next k loop.
! The should (possibly) be a (1-cldfrc), in here as we evaporate only in
! the cloud free part of the grid. Thus only the cloud free vapour
! should be interacting with the liquid. However, this would mean the
! entire grid cell will be out of equilibrium. 
!
           fclr = 1. - fcloud*cldfrc(i,k)
           tbar  = state%t(i,k) + 0.5*deltat*ptend%s(i,k)/cpair
           alpliq = wiso_alpl(iwspec(m), tbar)
           alpliq = wiso_akel(iwspec(m), tbar, hum0(i,k),alpliq)
!
           qtot  = fclr*state%q(i,k,m) + &
                deltat*(evpsnow(i,m)+ snowmlt(i,k,m) + flxrain(i,k,m)*gravit/state%pdel(i,k))
!
! Compute initial state assume no fractionation
!
           vapiso = fclr*state%q(i,k,m) + deltat*(evpsnow(i,m))
           liqiso = deltat*(snowmlt(i,k,m) + flxrain(i,k,m)*gravit/state%pdel(i,k))
!
! Equilibrate the liquid and vapour
!
           call wiso_liqvap_equil(alpliq, feq, vaptot(i,k), liqtot(i,k), &
                            vapiso, liqiso, dliqiso)
!
           evpprec(i,m) = (vapiso - fclr*state%q(i,k,m))/deltat 
           evprain(i,m) = evpprec(i,m) - evpsnow(i,m)
!
! Apply mass conservation limiter for evaporation...
!   o do not evap more than precipitation flux throu upper interface
!   o do not evap more than column precip (this SHOULD be a gimme)
!
           evplimit = min(flxprec(i,k,m), prec(i,m)-evpvint(i,m)) * gravit/state%pdel(i,k)
!
!! Limit usually order 1.e-18 to 1.e-16, occasional 1.e-15 and now and again 1.e-9!!
!! And this is when evplimit = 0.
!!           if (evpprec(i,m) - evplimit > 1.e-18) then
!!           if (evpprec(i,m) - evplimit > 1.e-12) then
!!           if (evpprec(i,m) - evplimit > 1.e-11) then	!  10*qmin(1)
!! This is getting ridiculous! this should be a very robust routine.
!! It should be re coded to give better numerical precision.
!!           if (evpprec(i,m) - evplimit > 1.e-10) then	!  100*qmin(1)
           if (evpprec(i,m) - evplimit > 1.e-9) then	!  1000*qmin(1)
             write(*,*) 'WTRC_ZM_EVAP - evap exceeded maximum: m=',m
!!             write(*,*)  ' - evpprec:',evpprec(i,m)
!!             write(*,*)  flxprec(i,k,m)*gravit/state%pdel(i,k), &
!!                         (prec(i,m)-evpvint(i,m))*gravit/state%pdel(i,k)
             write(*,*)  ' Qerr=:',i,k,evplimit-evpprec(i,m)
           end if
           evpprec(i,m) = min(evplimit, evpprec(i,m))
!
! Vertically integrated evaporation
!
           evpvint(i,m) = evpvint(i,m) + evpprec(i,m)*state%pdel(i,k)/gravit
!
! Calculate net production of isotopes in precip. and snow
!
           ntprprd(i,k,m) = prdprec(i,k,m) - evpprec(i,m)
!
! Compute net snow production.
! Assume process is similar to that for H2O. Guess the fraction of
! precip which is snow - this SHOULD depend on and be tracked by cloud physics.
! The precip should by tracked by type from the convective scheme as
! the precipitation output as ice and snow will have different isotopic
! concentrations due to the cloud microphysics. This needs to be fixed.
! For now, we realize that this will only affect the result
! of the isotope amount in the snow (i.e., flxsnow), not the vapor
! phase tracer, and as such will be Ok for a first guess (dcn+swong)
!   
#ifdef PERGRO
           work1 = min(max(0.,flxsnow(i,k,m)/(flxprec(i,k,m)+8.64e-11_r8)),1.)
           write(*,*) 'PERGRO NOT CHECKED in wtrc_zm_conv_evap'
           call endrun
#else
           if (flxprec(i,k,1) > 0.) then
              work1 = min(max(0.,flxsnow(i,k,1)/flxprec(i,k,1)),1.) ! gues snow fraction
           else
              work1 = 0.
           end if
#endif
           work2 = max(fsnow_conv(i,k),work1)	!  guess again 
           if (snowmlt(i,k,1) > 0.) work2 = 0
           ntsnprd(i,k,m) = work2*prdprec(i,k,m) - evpsnow(i,m)  - snowmlt(i,k,m)

! precipitation fluxes for isotopes
           flxprec(i,k+1,m) = flxprec(i,k,m) + ntprprd(i,k,m)*state%pdel(i,k)/gravit
           flxsnow(i,k+1,m) = flxsnow(i,k,m) + ntsnprd(i,k,m)*state%pdel(i,k)/gravit
           flxrain(i,k+1,m) = flxprec(i,k+1,m) - flxsnow(i,k+1,m)

! protect against rounding error
           flxprec(i,k+1,m) = max(flxprec(i,k+1,m), 0.)
           flxsnow(i,k+1,m) = max(flxsnow(i,k+1,m), 0.)

! Record the tendencies
           ptend%q(i,k,m) = evpprec(i,m)  
!
         end do      ! i-loop 
       end do        ! k-loop top down
!
! Set output isotope abundance in precip. and snow
!
       prec(:ncol,m) = flxprec(:ncol,pver+1,m)/1000.
       snow(:ncol,m) = flxsnow(:ncol,pver+1,m)/1000.     
!
! record history variables (these need to be added master field list)
       write(atr,'(i2.2)') m
       if (hackflag) then
         call outfld('HKVIEV'//atr, evpvint(:,m)  , pcols, lchnk)

         call outfld('HKFXPR'//atr, flxprec(:,:,m), pcols, lchnk)
         call outfld('HKFXSN'//atr, flxsnow(:,:,m), pcols, lchnk)
         call outfld('HKNPPD'//atr, ntprprd(:,:,m), pcols, lchnk)
         call outfld('HKNSPD'//atr, ntsnprd(:,:,m), pcols, lchnk)
       else
         call outfld('ZMVIEV'//atr, evpvint(:,m)  , pcols, lchnk)

         call outfld('ZMFXPR'//atr, flxprec(:,:,m), pcols, lchnk)
         call outfld('ZMFXSN'//atr, flxsnow(:,:,m), pcols, lchnk)
         call outfld('ZMNPPD'//atr, ntprprd(:,:,m), pcols, lchnk)
         call outfld('ZMNSPD'//atr, ntsnprd(:,:,m), pcols, lchnk)
       end if

     end if     ! m is water vapour (all species)
   end do      	! m, water tracer loop 
!
   return
end subroutine wtrc_zm_conv_evap      
!

!=======================================================================
subroutine wtrc_zm_convr_cldprp ( lchnk   ,il1g    ,il2g    ,  &
                msg     ,jt      ,mx      ,dp      ,dsubcld ,  &
                mu      ,du      ,eu      ,md      ,ed      ,  &
                tg      ,tug     ,qg      ,qhat    ,qhup    ,qhdn   , &
                qupd    ,  &
                qliq    ,qice    ,qdnd    ,cug     ,evpg    ,rprdg   ,  &
                dzdp    ,pflxg   ,tdg     ,jd      ,rrg     ,  &
                ldetg   ,idetg   ,ficed   ,ficeup)

!-----------------------------------------------------------------------
!
! Purpose:
!
!   Compute updraft and dowdraft cloud propoerties given the mass fluxes
!   and large scale isotopic ratios.
!
!   This is the guts of the convective calulations
!
!   Work in interface fluxes, rather than mixing ratios to avoid
!   divide by zero problems with small (trivial) mass fluxes.
!
! Method:
!
!   The isotopic fractionation model assumes the partial removal of
!   liquid for formation of liquid condensate, with a rayleigh
!   distiallation for precipitation and total equilibrium of in cloud
!   water, with kinetic effects. A direct rayleigh fractionation is used 
!   for formation of ice. Isotopic equilibrium is assumes for water
!   detrained into the environment (with respect to the environment).
!   Isotopic equilibrium is assumed to be with respect to the
!   environment, with kinetic effects, following evaporation of ambient
!   liquid to maintain downdrafts. 
!
!   This code now works entirely with isotope ratios. This is because
!   I coundn't get the numerics sufficiently stable (precise) with
!   input mixing ration. This also by-passes the problem with using
!   mixing ratios rather than mole ratios, by actually using mole
!   ratios for the isotopic fractionation, and thus the modelq can be
!   considered some scaled value (by the ratio of molecular weights)
!   of the actual mixing ration.
!
!
! Notes:
!
!   One shortcomming of this code is that all (updraft) rain water is
!   available for evaporation in the downdraft. This is consitent with
!   the underlying formulation where the downdraft mass fluxes have
!   been scaled uniformly, in association with the total precip and
!   evap (the alpha factor). However, this not as physically sensible
!   as having fractional reductions of the downdraft in associated with
!   local (rather than column integral) evaporation. For the isotopes,
!   theprecipitation water is assumed to be available to the downdraft
!   from the top down, which will essentially give a downward
!   distillation for liquid. This is probably close to what is expected
!   with more depleted water comming from lower in the downdraft.
!   This present formulation, although simplistic, means no nevative
!   evaporations can occurs, and that the there is always enough
!   rainwater present to give saturation (to numerical precision).
!   It is possible to do this one layer at a time (because most of the
!   rain forms high up in the plume... in which case there isn't
!   actually much difference. Most (as far as I can tell), numerical
!   problems eventuate because the mass fluxes (especially downdrafts)
!   are very small, and appear the denominator. Thus, setting the
!   the mixing ratios to trivail values will have only small impact
!   as they will be agin multiplied by the mass fluxes (numerator) when
!   the teneancies are computed.
!
!   This version assumes "condensed phase" is either liquid, OR ice
!   (which is not quite right). We SHOULD be tracking cloud ice and liquid
!   separatly, even if the core ZM scheme does not, because the fractionation
!   will differ. Thus, detrainid ice, vesus, liquid questions could be troubled.
!
!
!
! Author:
!   David Noone <dcn@caltech.edu> - Tue Oct 28 13:56:56 MST 2003
!
!
!-----------------------------------------------------------------------
   use water_tracers, only: ixwti, ixwtx, wtrc_is_vap

   implicit none

!---------------------------- Arguments --------------------------------
   integer , intent(in)    :: lchnk             ! chunk identifier
   integer , intent(in)    :: il1g              ! fist index of column
   integer , intent(in)    :: il2g              ! last index of column
   integer , intent(in)    :: msg               ! lowest level in gather
   integer , intent(in)    :: mx(pcols)         ! level of convection base
   integer , intent(in)    :: jt(pcols)         ! level of convection top
   integer , intent(in)    :: jd(pcols)         ! level of downdraft top

   real(r8), intent(in)    :: mu(pcols,pver)    ! updraft mass flax
   real(r8), intent(in)    :: eu(pcols,pver)    ! updraft entrainment rate
   real(r8), intent(in)    :: du(pcols,pver)    ! updraft detrainment rate
   real(r8), intent(in)    :: md(pcols,pver)    ! downdraft masss flux
   real(r8), intent(in)    :: ed(pcols,pver)    ! downdraft entrainment rate
   real(r8), intent(in)    :: dp(pcols,pver)    ! level thichnesss (mb)
   real(r8), intent(in)    :: dsubcld(pcols)    ! sub-cloud thichness ( mb)
   real(r8), intent(in)    :: dzdp(pcols,pver)
!
   real(r8), intent(in)    :: tg(pcols,pver)    ! environment temperature
   real(r8), intent(in)    :: tug(pcols,pver)   ! updraft temperature
   real(r8), intent(in)    :: tdg(pcols,pver)   ! downdraft temperature
   real(r8), intent(in)    :: qg(pcols,pver,pcnst+pnats)    ! shum 'qbar'
   real(r8), intent(in)    :: qhat(pcols,pver,pcnst+pnats)  ! interface shum FLUX
   real(r8), intent(inout) :: qupd(pcols,pver,pcnst+pnats)  ! updraft vapour isotope mass FLUX
   real(r8), intent(inout) :: qliq(pcols,pver,pcnst+pnats)  ! updraft liquid isotope mass FLUX
   real(r8), intent(inout) :: qice(pcols,pver,pcnst+pnats)  ! updraft ice isotope mass FLUX
   real(r8), intent(inout) :: qdnd(pcols,pver,pcnst+pnats)  ! downdraft vapour isotope mass FLUX
   real(r8), intent(inout) :: cug(pcols,pver,pcnst+pnats)   ! condensation isotope mass
   real(r8), intent(inout) :: rrg(pcols,pver,pcnst+pnats)   ! rain rate isotope mass
   real(r8), intent(inout) :: evpg(pcols,pver,pcnst+pnats)  ! evaporation isotope mass
   real(r8), intent(inout) :: rprdg(pcols,pver,pcnst+pnats) ! precip production at each level
   real(r8), intent(inout) :: pflxg(pcols,pverp,pcnst+pnats)! precip flux at each interface

   real(r8), intent(out)   :: qhup(pcols,pver,pcnst+pnats)  ! updraft interface shum FLUX
   real(r8), intent(out)   :: qhdn(pcols,pver,pcnst+pnats)  ! downdraft interface shum FLUX
   real(r8), intent(out)   :: ldetg(pcols,pver,pcnst+pnats) ! detraining liquid isotope mass
   real(r8), intent(out)   :: idetg(pcols,pver,pcnst+pnats) ! detraining liquid isotope mass
   real(r8), intent(in)    :: ficed(pcols,pver)             ! ice fraction detrained
   real(r8), intent(in)    :: ficeup(pcols,pver)            ! ice fraction in updraft
!
!------------------------- Local Variables -----------------------------

   integer  :: k,m                            ! level index

!-----------------------------------------------------------------------
   if (il1g /= 1) call endrun('wtrc_cldprp assumes il1g = 1')
!
! Check inputs
!      
   call wtrc_qchk2('WTRC_ZMCcld','qg  ',il2g,qg(:,:,ixh2oq) ,qg(:,:,1))
   call wtrc_qchk2('WTRC_ZMCcld','qhat',il2g,qhat(:,:,ixh2oq) ,qhat(:,:,1))
!
! Convert inputs to interface fluxes in updrafts and downdrafts
! (these are the rerms appearing in the tendency equation)
!
   do k = msg+1,pver
     qliq(:il2g,k,1) = mu(:il2g,k)*qliq(:il2g,k,1)		! > 0
     qice(:il2g,k,1) = mu(:il2g,k)*qice(:il2g,k,1)		! > 0
     qupd(:il2g,k,1) = mu(:il2g,k)*qupd(:il2g,k,1)		! > 0
     qhup(:il2g,k,1) = mu(:il2g,k)*qhat(:il2g,k,1)		! > 0
     qdnd(:il2g,k,1) = md(:il2g,k)*qdnd(:il2g,k,1)		! < 0
     qhdn(:il2g,k,1) = md(:il2g,k)*qhat(:il2g,k,1)		! < 0
   end do
!
! Initialize tracer arrays
!
   do m = ixwti, ixwtx
     if (wtrc_is_vap(m)) then
        do k = msg+1,pver
          qliq (:il2g,k,m) = 0.
          qice (:il2g,k,m) = 0.
          qupd (:il2g,k,m) = mu(:il2g,k)*qg  (:il2g,k,m)	! > 0
          qdnd (:il2g,k,m) = md(:il2g,k)*qg  (:il2g,k,m)	! < 0
          qhup (:il2g,k,m) = mu(:il2g,k)*qhat(:il2g,k,m)	! > 0
          qhdn (:il2g,k,m) = md(:il2g,k)*qhat(:il2g,k,m)	! < 0
          cug  (:il2g,k,m) = 0.
          rrg  (:il2g,k,m) = 0.
          evpg (:il2g,k,m) = 0.
          rprdg(:il2g,k,m) = 0.
        end do
     end if          ! m is water
   end do            ! m, tracers
!
! Updraft properties: output qupd, qlid, cug, rprdg
!
   call wtrc_zm_updraft( lchnk  , il1g    , il2g    , &
                msg     , mx     , jt      , mu      , eu      , &
                du      , dp     , dzdp    , tg      , tug     , qg  ,&
                qupd    , qliq   , qice    , cug     , rrg     , rprdg   , &
                ldetg   , idetg  , ficed   ,ficeup )
!
! Downdraft properties: output qnd, evpg, rprdg, pflxg
!
   call wtrc_zm_dndraft( lchnk  , il1g     , il2g    , msg     , &
               mx      , jd     , md       , ed      , dp      , &
               tg      , tdg    , qg       , qdnd    , rrg     , evpg    , &
               rprdg   , pflxg)
!
   return
  end subroutine wtrc_zm_convr_cldprp

!=======================================================================
  subroutine wtrc_zm_updraft( lchnk  , il1g    , il2g    , &
                    msg     , mx     , jt      , mu      , eu      , &
                    du      , dp     , dzdp    , tg      , tug     , qg ,&
                    qupd    , qliq   , qice    , cug     , rrg     , rprdg   , &
                    ldetg   , idetg  , ficed   , ficeup)
!-----------------------------------------------------------------------
!
! UPDRAFT PROPERTIES: compute tracer qu, ql, cu, rr
!
! Recompute statitics from total water to account for all mass.
! Specifically, compute the part which is detrainment plus condensation.
! Which is just the initial (pre-condensation) vapour for the layer.
! To do the fractionation, we need both the liquid and the updraft
! water, as we integrate up the plume(s). All isotopic exchange is now
! done in the external routine dicm.
!
! David Noone <dcn@colorado.edu> - Mon Jul 19 18:26:52 MDT 2004
! David Noone <dcn@colorado.edu> - Mon Jul 26 16:53:51 MDT 2004 (simplified code)
!
!-----------------------------------------------------------------------
   use water_tracers,      only: wtrc_ratio
   use water_isotopes,     only: wiso_dicm, wiso_delta
   use zm_conv,            only: c0

   implicit none
!-----------------------------------------------------------------------

!---------------------------- Arguments --------------------------------
   integer , intent(in)    :: lchnk             ! chunk identifier
   integer , intent(in)    :: il1g              ! fist index of column
   integer , intent(in)    :: il2g              ! last index of column
   integer , intent(in)    :: msg               ! lowest level in gather
   integer , intent(in)    :: mx(pcols)         ! level of convection base
   integer , intent(in)    :: jt(pcols)         ! level of convection top

   real(r8), intent(in)    :: mu(pcols,pver)    ! updraft mass flax
   real(r8), intent(in)    :: eu(pcols,pver)    ! updraft entrainment rate
   real(r8), intent(in)    :: du(pcols,pver)    ! updraft detrainment rate
   real(r8), intent(in)    :: dp(pcols,pver)    ! level thichnesss (mb)
   real(r8), intent(in)    :: dzdp(pcols,pver)
!
   real(r8), intent(in)    :: tg(pcols,pver)    ! environment temperature
   real(r8), intent(in)    :: tug(pcols,pver)   ! updraft temperature
   real(r8), intent(in)    :: qg(pcols,pver,pcnst+pnats)    ! shum 'qbar'
!
   real(r8), intent(inout) :: qupd(pcols,pver,pcnst+pnats)  ! updraft vapour isotope mass FLUX
   real(r8), intent(inout) :: qliq(pcols,pver,pcnst+pnats)  ! updraft liquid isotope mass FLUX
   real(r8), intent(inout) :: qice(pcols,pver,pcnst+pnats)  ! updraft liquid isotope mass FLUX
   real(r8), intent(inout) :: cug(pcols,pver,pcnst+pnats)   ! condensation isotope mass
   real(r8), intent(inout) :: rrg(pcols,pver,pcnst+pnats)   ! rain rate isotope mass (snow+rain)
   real(r8), intent(inout) :: rprdg(pcols,pver,pcnst+pnats) ! precip production at each level
   real(r8), intent(out)   :: ldetg(pcols,pver,pcnst+pnats) ! detraining liquid isotope mass flux
   real(r8), intent(out)   :: idetg(pcols,pver,pcnst+pnats) ! detraining ice isotope mass flux

   real(r8), intent(in)    :: ficed(pcols,pver)             ! fracition ice detrained 
   real(r8), intent(in)    :: ficeup(pcols,pver)            ! fracition ice in updraft

!------------------------- Local Variables -----------------------------
   integer  :: i,k,m                            ! level index

   real(r8) :: qxs(pcols,pver,pcnst+pnats)	! excess needed to balance budgets
   real(r8) :: fsnow(pcols,pver)                ! initial snow fraction in precip production

   real(r8) :: duqst(pcols,pver)                ! detraining vapour
   real(r8) :: rng(pcols,pver,pcnst+pnats)      ! rain rate isotope mass (rr = rn+sn)
   real(r8) :: sng(pcols,pver,pcnst+pnats)      ! snow rate isotope mass (rr = rn+sn)

   real(r8) :: netdet				! net (liquid+ice) detrainment
   real(r8) :: ntlsnk,trntlsnk			! net liquid sink
   real(r8) :: ntisnk,trntisnk			! net ice sink
   real(r8) :: netliq,trnetliq			! net liquid retained
   real(r8) :: netice,trnetice			! net liquid retained
   real(r8) :: fntld, fntid			! detrained fraction of retained liquid/ice
   real(r8) :: flskd, fiskd			! detrained fraction of liquid/ice sink 

   real(r8) :: rat
   real(r8) :: netall,trnetall
!
! Variables passed in/out of fractionation sub-model
!
   integer, parameter :: piso = 2       ! solve for one species at a time
   integer  :: isp(piso)			! species indicies (1 is water, 2 is isotope)
   real(r8) :: vapold(piso)			! initial vapour
   real(r8) :: liqold(piso)			! initial liquid
   real(r8) :: iceold(piso)			! initial ice
   real(r8) :: vapent(piso)			! vapour entrainment
   real(r8) :: vapnew(piso)			! final vapour
   real(r8) :: liqnew(piso)			! final liquid
   real(r8) :: icenew(piso)			! final ice
   real(r8) :: vapdet(piso)			! vapour detrainment
   real(r8) :: rainpr(piso)			! rain production
   real(r8) :: snowpr(piso)			! snow production
   real(r8) :: dliqmt				! liquid change due to ice melt
   real(r8) :: dicefz				! ice change due to liquid freeze
   real(r8) :: told, tnew			! start, end fractionation temperature
!
! debugging variables
!
  character(len=12) cerr			! name of toubled variable
!!  integer  :: idbg = 6
!!  integer  :: idbg = 0		! no messages, traps only
  integer  :: idbg = -1		! all messages and traps

#undef DBGUPDFT
#undef UPDIAGS

!-----------------------------------------------------------------------
!
! Initialize isotope ratios as necessary
!
   qxs(:il2g,:,:)    = 0.0
   duqst(:il2g,:)    = 0.0
   ldetg(:il2g,:,:)  = 0.0			! needed in q1q2 as isotope is not as qliq(k+1)
   idetg(:il2g,:,:)  = 0.0			! needed in q1q2 as isotope is not as qice(k+1)
   rng(:il2g,:,:)    = 0.0
   sng(:il2g,:,:)    = 0.0
!
! Assign snow fraction to be the same as the ice fraction.
! This will be repartitioned later in evaporation calulations to the
! correct fsnow
!
    fsnow(:,:) = ficed(:,:)
!
! Remove any dodgy negatives - save these excesses, which we'll add as "environment" air
   do k = msg, pver
     do i = 1, il2g
       qxs(i,k,1) = 0.0			! do summation for all water fluxes
       if (qupd(i,k,1) < 0.) then
         write(*,*) 'WTRC_UPDRAFT qupd(1) < 0:',i,k,qupd(i,k,1)
         qxs(i,k,1)  = qxs(i,k,1) - min(qupd(i,k,1),0.)
         qupd(i,k,1) = max(qupd(i,k,1), 0.)
       end if
       if (qliq(i,k,1) < 0.) then
         write(*,*) 'WTRC_UPDRAFT qliq(1) < 0:',i,k,qliq(i,k,1)
         qxs(i,k,1)  = qxs(i,k,1) - min(qliq(i,k,1),0.)
         qliq(i,k,1) = max(qliq(i,k,1), 0.)
       end if
       if (qice(i,k,1) < 0.) then
         write(*,*) 'WTRC_UPDRAFT qice(1) < 0:',i,k,qice(i,k,1)
         qxs(i,k,1)  = qxs(i,k,1) - min(qice(i,k,1),0.)
         qice(i,k,1) = max(qice(i,k,1), 0.)
       end if
     end do
   end do
!
!!! count the number of problems (about 300 per step at T42, that's 1/10000 gridpoints)
!!   if (count(qxs(:,:) /= 0.) > 0) then		! what to do, what to do...
!!     write(*,*) 'wtrc_updraft: number of negatives from zmc:',count(qxs(:,:) /= 0.)
!!   end if
!
! Extract the ACTUAL detrainment rates from known cloud base properties.
! (thre is actually only two interestring lines of code here - duqst, ldetg)
!
      do k = pver,msg+2,-1      ! loop, bottom up
        do i = 1, il2g
!
! Compute total vapour available for fractionation 
!  (upflux plus entrainment ... detrain and condense later)
! Also, diagnose the detraining vapour and liquid required for mass balance.
!
          if (k < mx(i) .and. k>=jt(i) .and. mu(i,mx(i)) > 0.) then ! jt for ql and cu, qu(jt)=0)
!
            duqst(i,k) = (qupd(i,k+1,1)-qupd(i,k,1))/dp(i,k) - cug(i,k,1) + eu(i,k)*qg(i,k,1)
            ldetg(i,k,1) = du(i,k)*qliq(i,k+1,1)/(mu(i,k+1) + msmall)       ! l+, as guang
            idetg(i,k,1) = du(i,k)*qice(i,k+1,1)/(mu(i,k+1) + msmall)       ! i+, as guang
!
! Solve budgets for rain and snow partitions
! (this is not just fice*rrg, as we need to include the convergence of
!  fice at interfaces)
!
!!            rng(i,k,1) = (qliq(i,k+1,1)-qliq(i,k,1))/dp(i,k) - ldetg(i,k,1) + (1.-ficed(i,k))*cug(i,k,1)
!!            sng(i,k,1) = (qice(i,k+1,1)-qice(i,k,1))/dp(i,k) - idetg(i,k,1) + (   ficed(i,k))*cug(i,k,1)
 
            rng(i,k,1) = (1.-ficed(i,k))*rrg(i,k,1)
            sng(i,k,1) = (   ficed(i,k))*rrg(i,k,1)
!
! Check that we could recompute rrg accurately.
!
            if (abs(rng(i,k,1)+sng(i,k,1)-rrg(i,k,1)) > qtiny) then
               write(*,*) '(wtrc_updraft) Failed to recompute rrg.'
            end if

            if (ldetg(i,k,1) < 0.) then
               write(*,*) 'LDETG < 0:',ldetg(i,k,1)
            end if
            if (idetg(i,k,1) < 0.) then
               write(*,*) 'IDETG < 0:',idetg(i,k,1)
            end if
! Checks
            if (abs(du(i,k)) < 1.e-12) then
              if (abs(duqst(i,k)) > 1.0e-20 ) then
                 write(*,*) '(wtrc_zm_conv_UPDRFT) QDET > 0:',i,k,du(i,k),duqst(i,k),qupd(i,k,1)
                 write(*,*) '   :',qupd(i,k  ,1)/dp(i,k),qupd(i,k+1,1)/dp(i,k), &
                                   eu(i,k)*qg(i,k,1),cug(i,k,1)
              endif
            end if
! The following will occure if qu is < 0...
            if (du(i,k)<0. .and. abs(duqst(i,k)) > qtiny) then
               write(*,*) '(wtrc_zm_convr) DUQST /= with du=0'
               write(*,2) 'DUQST:',i,k,duqst(i,k),du(i,k)
!!              call endrun
            end if
          end if        ! k, within plume, and mass flux > 0
!
! Do all check for physical consistancy...
! This is the part that pisses me off. The GCM value qu, CAN become
! negative. Which essentially means the mass fluxes (or detrainment)
! is incorrect. This is now fixed in cldprp to be a max(qu,0), 
! which should fix most problems here to numerical precision.
!  (note, qu<0 begats cu<0 begats rr<0...)
!
         if (qupd(i,k,1) < 0.) then
           write(*,2)'qupd< 0:i,k',i,k,jt(i),qupd(i,k,1),qupd(i,k+1,1),mu(i,k)
           qupd(i,k,1) = max(qupd(i,k,1),0.)
!!         call endrun('Negative qupd')
         end if
!
         if (qliq(i,k,1) < 0.) then
           write(*,2)'QLIQ< 0:i,k',i,k,jt(i),qliq(i,k,1),mu(i,k)
           qliq(i,k,1) = max(qliq(i,k,1),0.)
!!              call endrun('Negative qliq.')
         end if
        end do          ! i, columns
      end do            ! k, levels
  2   format(a16,3i3,9e12.4)
!
! Loop over all WATER TRACERS (m)
!
   do m = ixwti, ixwtx
     if (wtrc_is_vap(m)) then
!!       write(*,*) 'WTRC_ZM: m=',m
!
! Cloud base, no condensation (no liquid), no fractionation, only transport.
! Note: mu = dp*eu at the base.
!
       do i = 1, il2g 
         if (mu(i,mx(i)) > 0.) then
           qupd(i,mx(i),m) = (dp(i,mx(i))*eu(i,mx(i))*qg(i,mx(i),m))

! Truncate because negative qu can exist in cldprp
!           if (qupd(i,mx(i),1) <= 1.1e-22) qupd(i,mx(i),m) = 0.
           if (qupd(i,mx(i),m) <= 0.) then
               write(*,*) 'WTRC_UPDRAFT: qupd(mx) <0:',qupd(i,mx(i),m)
               qupd(i,mx(i),m) = 0.
           end if
           if (qupd(i,mx(i),1) <= 0.) then
               write(*,*) 'WTRC_UPDRAFT: qupd(mx,1)=0 set qupd(m)=0'
               qupd(i,mx(i),m) = 0.	! this one is trivial
           end if
#ifdef DBGUPDFT
           if (m==4 .and. i==idbg) then
             write(*,*) 'DEBUG at base  i,k:',i,mx(i)
             write(*,11) 'Cu, Rr :',cug (i,mx(i),1),rrg (i,mx(i),1) ! == 0, defined
             write(*,11) 'QG   MX:',qg  (i,mx(i),m),qg  (i,mx(i),1), &
                                    qg  (i,mx(i),m)-qg  (i,mx(i),1)
             write(*,11) 'QUPD MX:',qupd(i,mx(i),m),qupd(i,mx(i),1), &
                                    qupd(i,mx(i),m)-qupd(i,mx(i),1)
  11               format(a8,3e16.6)
             write(*,*) '------------------------------------------------'
           end if
#endif
!
         end if		! mu > 0
       end do		! i, columns
!
! Cloud plume from groud up, with known cloud base properties.
!
       do k = pver,msg+2,-1      ! loop, bottom up
         do i = 1, il2g
           if (k < mx(i) .and. k>=jt(i) .and. mu(i,mx(i)) > 0.) then ! jt for ql and cu, qu(jt)=0)
!
! SOLVE ISOTOPIC FRACTIONATION IN UPDRAFTS
!
! Use the mixed cloud model.
!  Distillation to ice, liquid in equilibrium with ice, rain and snow
!  autoconverted from liquid and ice. 
!
! In contrast to the underlying code which assumes the liquid detraining is the 
! same as that at the lower interface (which would lead to te assumption that the 
! detraining liquid and vapour does not undergo fractionation), the
! isotope model assumes that it is associated with liquid and vapour
! that HAS undegone fractionation before bing dumped to the environment. 
! This assumtion, while also being more correct, means that the updraft
! liquid available for fractionation can not go negative, and indeed is
! simple the liquid transported accross the lower interface.
! Essentially what dicm does is to compute cu for isotopes, which must be iterative.
!
! Two possible assumtions for detraining liquid/ice
!   Both are not accurate: detrainment is not the precipitation sink, as precip
!   forms from large drops in the cloud, while detrainment is small
!   drops. However, doing the detrainment of the final water is also not
!   accurate as it does not account for the evolution of drops during
!   the timestep.... either probably fine.
!
!
                netliq = qliq(i,k,1) + (1.-fsnk_zmc)*ldetg(i,k,1)*dp(i,k)
                netice = qice(i,k,1) + (1.-fsnk_zmc)*idetg(i,k,1)*dp(i,k)
                ntlsnk = (rng(i,k,1) + (   fsnk_zmc)*ldetg(i,k,1))*dp(i,k)
                ntisnk = (sng(i,k,1) + (   fsnk_zmc)*idetg(i,k,1))*dp(i,k)

                fntld = wtrc_ratio((1.-fsnk_zmc)*ldetg(i,k,1)*dp(i,k), netliq)
                fntid = wtrc_ratio((1.-fsnk_zmc)*idetg(i,k,1)*dp(i,k), netice)
                flskd = wtrc_ratio((   fsnk_zmc)*ldetg(i,k,1)*dp(i,k), ntlsnk)
                fiskd = wtrc_ratio((   fsnk_zmc)*idetg(i,k,1)*dp(i,k), ntisnk)
!
! Determine melt as that required to give the expected 
! rain/snow ratio ? Perhaps just direct melt up updafraft material?
! What about some random number.. like 1/2 snow production?
!
                dliqmt = 0.
                dicefz = 0.
                if (ficeup(i,k) > ficeup(i,k+1)) then	! freeze in updraft
                   dliqmt = 0.
                   dicefz =  0.5*(ficeup(i,k) - ficeup(i,k+1)) * qliq(i,k+1,1) 
                else					! melt in updraft?
                   dliqmt = -0.5*(ficeup(i,k) - ficeup(i,k+1)) * qice(i,k+1,1) 
                   dicefz = 0.
                end if
!
                isp(1) = iwspec(1)
                vapold(1) = qupd(i,k+1,1)
                liqold(1) = qliq(i,k+1,1)
                iceold(1) = qice(i,k+1,1)
                vapent(1) = (dp(i,k)*eu(i,k)*qg(i,k,1) + max(-dp(i,k)*duqst(i,k),0.))
!
                isp(2) = iwspec(m)
                vapold(2) = qupd(i,k+1,m)
                liqold(2) = qliq(i,k+1,m)
                iceold(2) = qice(i,k+1,m)
                vapent(2) = (dp(i,k)*eu(i,k)*qg(i,k,m))
!
                vapnew(1) = qupd(i,k,1)
                liqnew(1) = netliq
                icenew(1) = netice
                rainpr(1) = ntlsnk
                snowpr(1) = ntisnk
                vapdet(1) = max(dp(i,k)*duqst(i,k), 0.)
!
! ADjust inputs for special case when vapold and vapnew = 0, and precip
! does not (need a way to get from liquid to ice). Here, 
! just initially repartition the liquid and ice to be the same as the rain/snow.
! dicm is not applicable in this case.
!
                if (abs(vapold(1)) < qtiny .and. abs(vapnew(1)) < qtiny) then
!!                    write(*,*) '(wtrc_updraft) no vapour: all change homogeneous.'
                    netall   = vapold(1) + liqold(1) + iceold(1) + vapent(1)
                    trnetall = vapold(2) + liqold(2) + iceold(2) + vapent(2)
                    rat = wtrc_ratio(trnetall, netall)

                    vapnew(2) = vapold(2) + rat*(vapnew(1) - vapold(1))
                    liqnew(2) = liqold(2) + rat*(liqnew(1) - liqold(1))
                    icenew(2) = iceold(2) + rat*(icenew(1) - iceold(1))

!                    netall   = netall   - vapnew(1) - liqnew(1) - icenew(1)
!                    netall   = max (trnetall, 0.)		! should be
!                    trnetall = trnetall - vapnew(2) - liqnew(2) - icenew(2)
!                    trnetall = max (trnetall, 0.)		! should be
!                    rat = wtrc_ratio(trnetall, netall)
                    rainpr(2) = rat*rainpr(1)
                    snowpr(2) = rat*snowpr(1)

!                    netall   = netall   - rainpr(1) - snowpr(1)
!                    netall   = max (trnetall, 0.)		! should be
!                    trnetall = trnetall - rainpr(2) - snowpr(2)
!                    trnetall = max (trnetall, 0.)		! should be
!                    rat = wtrc_ratio(trnetall, netall)
                    vapdet(2) = rat*vapdet(1)			! det as residule

                    go to 111

!!                  if (vapnew(1) > qtiny) then		! try and make it ... might work?
!!                    write(*,*) '(wtrc_updraft) Daringly faking dicm inputs for zero vapour case'
!!                    liqold(1) = qliq(i,k+1,1)
!!                    iceold(1) = qice(i,k+1,1)
!!                    liqold(2) = qliq(i,k+1,m)
!!                    iceold(2) = qice(i,k+1,m)
!!                    liqnew(1) = qliq(i,k  ,1)
!!                    icenew(1) = qice(i,k  ,1)
!!                  else
!!                    vapnew(2) = 0.0
!!                    liqnew(2) = liqold(2)*wtrc_ratio(liqnew(1),liqold(1))
!!                    icenew(2) = iceold(2)*wtrc_ratio(icenew(1),iceold(1))
!!                    vapdet(2) = vapold(2) - vapnew(2) 		! i.e. trivially zero
!!                    rainpr(2) = liqold(2) - liqnew(2)
!!                    snowpr(2) = iceold(2) - icenew(2)
!!                    write(*,*) '(wtrc_updraft) NO VAPOUR - non-fractionation repartitition to bypass dicm'
!!                    go to 111			! BY-PASS DICM
!!                  end if
                end if
!
! Assign temperature - control for numberical problems from zm output (T very small)
!
                told = tug(i,k+1)
                told = max(told, 160.)
                told = min(told, 330.)
                tnew = tug(i,k)
                tnew = max(tnew, 160.)
                tnew = min(tnew, 330.)
!
! Solve for the fractionation system.
!   dicm assumes vapour is present, if not, just move stuff around. No
!
!!                call t_startf('wiso_dicm_zm_updraft')
                call wiso_dicm(piso   ,  &		! piso = 2
                               isp    , feq_upd, told   , tnew   , &
                               vapold , liqold , iceold , vapent , &
                               vapnew , liqnew , icenew , vapdet , &
                               rainpr , snowpr , dliqmt , dicefz )
!!                call t_stopf('wiso_dicm_zm_updraft')
!
! Unpack dicm outputs to updraft variables 
!
111             continue
                qupd(i,k,m)  = vapnew(2) 
                trnetliq     = liqnew(2) 
                trnetice     = icenew(2)
                trntlsnk     = rainpr(2) 
                trntisnk     = snowpr(2)

                ldetg(i,k,m) = (flskd*trntlsnk + fntld*trnetliq)/dp(i,k)
                idetg(i,k,m) = (fiskd*trntisnk + fntid*trnetice)/dp(i,k)
                rng  (i,k,m) = (1.-flskd)*trntlsnk/dp(i,k)
                sng  (i,k,m) = (1.-fiskd)*trntisnk/dp(i,k)
                qliq (i,k,m) = (1.-fntld)*trnetliq 
                qice (i,k,m) = (1.-fntid)*trnetice 
!
! Get cu from budget equation, as a residule (this is condensed ice and liquid)
!
                cug(i,k,m) = (qupd(i,k+1,m)-qupd(i,k,m))/dp(i,k) + (vapent(2) - vapdet(2))/dp(i,k)
!
! Correct for integral imprecision at top of plume, add excesses to cu and rr
! (error can be 1.e-8, dumping them to rr and cu means they will be removed form the model)
!
                if ( k == jt(i)) then		! qu(jt) = 0., ql(jt) = 0.
                   rng(i,k,m) = rng(i,k,m) + qliq(i,k,m) / dp(i,k)
                   rng(i,k,m) = max(rng(i,k,m), 0.)		! qliq possible < 0
                   sng(i,k,m) = sng(i,k,m) + qice(i,k,m) / dp(i,k)
                   sng(i,k,m) = max(sng(i,k,m), 0.)		! qliq possible < 0
                   qliq(i,k,m) = 0.
                   qice(i,k,m) = 0.
!
                   cug(i,k,m) = cug(i,k,m) + qupd(i,k,m) / dp(i,k)
                   cug(i,k,m) = max(cug(i,k,m), 0.)		! qupd possibly < 0
                   qupd(i,k,m) = 0.
                end if
!
! Recombine rain and snow for output to downdrafts/precipitation
!
                rrg(i,k,m) = rng(i,k,m) + sng(i,k,m)
!
! Check for negatives - should be numerical (warn is TOO negative, when zm fails)
                if (qupd(i,k,m) < -10*qtiny) write(*,*) '(WTRC_UPDRAFT) qupd < 0: m=',m, qupd(i,k,m)
                qupd(i,k,m) = max(qupd(i,k,m), 0.)
                if (qliq(i,k,m) < -10*qtiny) write(*,*) '(WTRC_UPDRAFT) qliq < 0: m=',m, qliq(i,k,m)
                qliq(i,k,m) = max(qliq(i,k,m), 0.)
                if (cug (i,k,m) < -10*qtiny) write(*,*) '(WTRC_UPDRAFT) cug  < 0: m=',m, cug (i,k,m)
                cug (i,k,m) = max(cug (i,k,m), 0.)
                if (rng (i,k,m) < -10*qtiny) write(*,*) '(WTRC_UPDRAFT) rng  < 0: m=',m, rng (i,k,m)
                rng (i,k,m) = max(rng (i,k,m), 0.)
                if (sng (i,k,m) < -10*qtiny) write(*,*) '(WTRC_UPDRAFT) sng  < 0: m=',m, sng (i,k,m)
                sng (i,k,m) = max(sng (i,k,m), 0.)
                if (rrg (i,k,m) < -10*qtiny) write(*,*) '(WTRC_UPDRAFT) rrg  < 0: m=',m, rrg (i,k,m)
                rrg (i,k,m) = max(rrg (i,k,m), 0.)
!
! Also check for consistency with total
! (this will fix negative qu allowed in cldprp, which makes no sence for tracers)
!
                if (qupd(i,k,1) <= 0.) then	! no total? no tracer
                   if (abs(qupd(i,k,m)) > qtiny) &
                       write(*,*) '(wtrc_updraft) set qupd = 0, as qupd(1) <= 0:',&
                              m,i,k,jt(i),qupd(i,k,1),qupd(i,k,m)
                   qupd(i,k,m) = 0.
                end if
                if (qliq(i,k,1) <= 0.) then	! no total? no tracer
                   if (abs(qliq(i,k,m)) > qtiny) &
                       write(*,*) '(wtrc_updraft) set qliq = 0, as qliq(1) <= 0:',&
                              m,i,k,jt(i),qliq(i,k,1),qliq(i,k,m)
                   qliq(i,k,m) = 0.
                end if
                if (qice(i,k,1) <= 0.) then	! no total? no tracer
                   if (abs(qice(i,k,m)) > qtiny) &
                       write(*,*) '(wtrc_updraft) set qice = 0, as qice(1) <= 0:',&
                              m,i,k,jt(i),qice(i,k,1),qice(i,k,m)
                   qice(i,k,m) = 0.
                end if
                if (cug(i,k,1) <= 0.) then	! no total? no tracer
                   if (abs(cug(i,k,m)) > 1.e-09) &	! all error
                       write(*,*) '(wtrc_updraft) set  cug = 0, as  cug(1) <= 0:',&
                              m,i,k,jt(i),cug(i,k,1),cug(i,k,m)
                   cug(i,k,m) = 0.
                end if
                if (rrg(i,k,1) <= 0.) then	! no total? no tracer
                   if (abs(rrg(i,k,m)) > 1.e-09) &	! all error
                       write(*,*) '(wtrc_updraft) set  rrg = 0, as  rrg(1) <= 0:',&
                              m,i,k,jt(i),rrg(i,k,1),rrg(i,k,m)
                   rrg(i,k,m) = 0.
                end if
!
! Dump some diagnostics
!
#ifdef UPDIAGS
!                if (m == 4 .and. i == idbg) then
                if (m == 4) then
                 write(*,*) 'UPFRAC:',i,k
                 write(*,*) 'FICE :',fice(i,k),fsnow(i,k)
!                 write(*,*) 'ALPHA :',alpliq,alpice
                 write(*,3) 'WT_ZM-UP QTULCR:',1, &
                   qg   (i,k,1)*1000., &
                   qupd (i,k,1)*1000., &
                   qliq (i,k,1)*1000., &
                   qice (i,k,1)*1000., &
                   cug  (i,k,1)*dp(i,k)*1000., &
                   rng  (i,k,1)*dp(i,k)*1000., &
                   sng  (i,k,1)*dp(i,k)*1000., &
                   rrg  (i,k,1)*dp(i,k)*1000.
                 write(*,3) 'WT_ZM-U dQTULCR:',m, &
                   wiso_delta(iwspec(m),qg   (i,k,m), qg   (i,k,1)), &
                   wiso_delta(iwspec(m),qupd (i,k,m), qupd (i,k,1)), &
                   wiso_delta(iwspec(m),qliq (i,k,m), qliq (i,k,1)), &
                   wiso_delta(iwspec(m),qice (i,k,m), qice (i,k,1)), &
                   wiso_delta(iwspec(m),rng  (i,k,m), rng  (i,k,1)), &
                   wiso_delta(iwspec(m),sng  (i,k,m), sng  (i,k,1)), &
                   wiso_delta(iwspec(m),rrg  (i,k,m), rrg  (i,k,1))
                   wiso_delta(iwspec(m),cug  (i,k,m), cug  (i,k,1))
  3              format(a16,i3,6f12.5)
               end if
#endif
!
! Debug diagnostics - with stops
!
#ifdef DBGUPDFT
               if (m==4) then
                 if (i == idbg .or. idbg < 0) then
                   write(*,*)  'DEBUG in plume i,k:',i,k
                   write(*,12) 'TUPD   :',tug(i,k),tug(i,k+1)
                   write(*,12) 'FICEUP :',ficeup(i,k), ficeup(i,k+1)
                   write(*,12) 'FICED  :',ficed(i,k)
                   write(*,*)  'FZ,MT  :',dicefz, dliqmt
                   write(*,12) 'FNETD  :',fntld,fntid
                   write(*,12) 'FSNKD  :',flskd,fiskd
                   write(*,12) 'MU     :',mu(i,k),mu(i,k+1)
                   write(*,12) 'DU,EU  :',du(i,k),eu(i,k)
                   write(*,12) 'QXS    :',qxs(i,k,1) ,qxs(i,k,m) ,qxs(i,k,m)-qxs(i,k,1)
                   write(*,12) 'QG     :',qg  (i,k,m),qg  (i,k,1),qg  (i,k,m)-qg  (i,k,1)
                   write(*,12) 'vapold I',vapold(2)  ,vapold(1)  ,vapold(2)-vapold(1)
                   write(*,12) 'liqold I',liqold(2)  ,liqold(1)  ,liqold(2)-liqold(1)
                   write(*,12) 'iceold I',iceold(2)  ,iceold(1)  ,iceold(2)-iceold(1)
                   write(*,12) 'vapent I',vapent(2)  ,vapent(1)  ,vapent(2)-vapent(1)
                   write(*,12) 'vapnew I',vapnew(2)  ,vapnew(1)  ,vapnew(2)-vapnew(1)
                   write(*,12) 'liqnew I',liqnew(2)  ,liqnew(1)  ,liqnew(2)-liqnew(1)
                   write(*,12) 'icenew I',icenew(2)  ,icenew(1)  ,icenew(2)-icenew(1)
                   write(*,12) 'vapdet I',vapdet(2)  ,vapdet(1)  ,vapdet(2)-vapdet(1)
                   write(*,12) 'rainpr I',rainpr(2)  ,rainpr(1)  ,rainpr(2)-rainpr(1)
                   write(*,12) 'snowpr I',snowpr(2)  ,snowpr(1)  ,snowpr(2)-snowpr(1)
                   write(*,12) 'netliq I',trnetliq   ,netliq     ,trnetliq -netliq
                   write(*,12) 'ntlsnk I',trntlsnk   ,ntlsnk     ,trntlsnk -ntlsnk
                   write(*,12) 'netice I',trnetice   ,netice     ,trnetice -netice
                   write(*,12) 'ntisnk I',trntisnk   ,ntisnk     ,trntisnk -ntisnk
                   write(*,12) 'ldetg  :',ldetg(i,k,m),ldetg(i,k,1),ldetg(i,k,m)-ldetg(i,k,1)
                   write(*,12) 'idetg  :',idetg(i,k,m),idetg(i,k,1),idetg(i,k,m)-idetg(i,k,1)
                   write(*,12) 'QUPD   :',qupd(i,k,m),qupd(i,k,1), qupd(i,k,m)-qupd(i,k,1)
                   write(*,12) 'QLIQ   :',qliq(i,k,m),qliq(i,k,1), qliq(i,k,m)-qliq(i,k,1)
                   write(*,12) 'QICE   :',qice(i,k,m),qice(i,k,1), qice(i,k,m)-qice(i,k,1)
                   write(*,12) 'RNG    :',rng (i,k,m),rng (i,k,1), rng (i,k,m)-rng (i,k,1)
                   write(*,12) 'SNG    :',sng (i,k,m),sng (i,k,1), sng (i,k,m)-sng (i,k,1)
                   write(*,12) 'RRG    :',rrg (i,k,m),rrg (i,k,1), rrg (i,k,m)-rrg (i,k,1)
                   write(*,12) 'RRDEC  :',rng (i,k,1)+sng(i,k,1),rrg (i,k,1), rng (i,k,1)+sng(i,k,1)-rrg (i,k,1)
                   write(*,12) 'CUG    :',cug (i,k,m),cug (i,k,1), cug (i,k,m)-cug (i,k,1)
 12                format(a8,3e16.6)
                   write(*,*) '------------------------------------------------'
                 end if	! idbg

!
! Check for errors in variables that are used later
!
                 cerr = 'none'
		 if (abs (idetg(i,k,m)-idetg(i,k,1))> 1.e-10) cerr = 'idetg'
		 if (abs (ldetg(i,k,m)-ldetg(i,k,1))> 1.e-10) cerr = 'ldetg'
                 if (abs ( cug (i,k,m)-cug (i,k,1)) > 1.e-12) cerr = 'cug'
                 if (abs ( rrg (i,k,m)-rrg (i,k,1)) > 1.e-12) cerr = 'rrg'
                 if (abs ( sng (i,k,m)-sng (i,k,1)) > 1.e-12) cerr = 'sng'
                 if (abs ( rng (i,k,m)-rng (i,k,1)) > 1.e-12) cerr = 'rng'
                 if (abs ( qice(i,k,m)-qice(i,k,1)) > 1.e-12) cerr = 'qice'
                 if (abs ( qliq(i,k,m)-qliq(i,k,1)) > 1.e-12) cerr = 'qliq'
                 if (abs ( qupd(i,k,m)-qupd(i,k,1)) > 1.e-12) cerr = 'qupd'

                 if (cerr /= 'none') then
                    write(*,*) '(wtrc_zm_updraft) conservation error for: '//trim(cerr),i,k
                    call endrun('wtrc_zm_updraft: debugging error')
                 endif

               end if
#endif
!
           end if    ! k in plume, and mass flux > 0
         end do      ! i, columns
!
       end do        ! k, levels
!
     end if          ! m is water
   end do            ! m, tracers

!
! Do a fix for updraft fluxes as they can be unstable occasionally
! Assume deficit has ratio of environment. This is a "qneg" type thing.'
! One can imagine this is a bit more entrainment.
! This is really just to ensure the qchks pass, and must be done after "m" loop.
!
   if (lupdcorr) then
     call updraft_correction('QUPD',il2g,pver,1,pver,1,ixh2oq,qupd,qg) ! for dqdt
     call updraft_correction('QLIQ',il2g,pver,1,pver,1,ixh2oq,qliq,qg) ! for dlf
     call updraft_correction('QICE',il2g,pver,1,pver,1,ixh2oq,qice,qg) ! for dlf
     call updraft_correction('CUG ',il2g,pver,1,pver,1,ixh2oq,cug ,qg) ! for dqemc
     call updraft_correction('RRG ',il2g,pver,1,pver,1,ixh2oq,rrg ,qg) ! for pflx
   endif
!
! After possible correction to rrg, save the rain rate in "rain production" 
! as per cldprp, for pflx below.
!
   rprdg(:il2g,:,2:pcnst+pnats) = rrg(:il2g,:,2:pcnst+pnats)
!
! Check conservation of total: qupd (qu), qcnd (as cu) and qliq (ql)
!
   call wtrc_qchk2('WTRC_ZMCucld','qupd',il2g,qupd(:,:,ixh2oq) ,qupd(:,:,1))
   call wtrc_qchk2('WTRC_ZMCucld','qliq',il2g,qliq(:,:,ixh2oq) ,qliq(:,:,1))
   call wtrc_qchk2('WTRC_ZMCucld','qice',il2g,qice(:,:,ixh2oq) ,qice(:,:,1))
   call wtrc_qchk2('WTRC_ZMCucld','cug ',il2g,cug (:,:,ixh2oq) ,cug (:,:,1))
   call wtrc_qchk2('WTRC_ZMCucld','rrg ',il2g,rrg (:,:,ixh2oq) ,rrg (:,:,1))
   call wtrc_qchk2('WTRC_ZMCucld','rng ',il2g,rng (:,:,ixh2oq) ,rng (:,:,1))
   call wtrc_qchk2('WTRC_ZMCucld','sng ',il2g,sng (:,:,ixh2oq) ,sng (:,:,1))
   call wtrc_qchk2('WTRC_ZMCucld','ldetg',il2g,ldetg(:,:,ixh2oq) ,ldetg(:,:,1))
   call wtrc_qchk2('WTRC_ZMCucld','idetg',il2g,idetg(:,:,ixh2oq) ,idetg(:,:,1))

   return
  end subroutine wtrc_zm_updraft

!=======================================================================
  subroutine wtrc_zm_dndraft( lchnk  , il1g     , il2g    , msg     , &
                    mx      , jd     , md       , ed      , dp      , &
                    tg      , tdg    , qg       , qdnd    , rrg     , evpg    , &
                    rprdg   , pflxg)
!-----------------------------------------------------------------------
!
! DOWNDRAFT PROPERTIES: compute tracer qd, ev, pflx
!
!-----------------------------------------------------------------------
!
! Do the isotopic part, assuming total equilibration with the
! falling liquid (with kinetic fractionartion), and no fractionation 
! for ice. That is, compute pflx as we go ensuring that it is in equilibrium
! with downdraft vapour.
!
! Assuming downdraft is saturated (as it SHOULD BE according to Z&M), 
! one could claim that there is isotopic equilibrium between downdraft and
! precipitation that exists the level. That is Rpflx(k+1) = alpha Rdnd(k).
! Given this, we can back out the isotopic evaporation. 
!
! (One catch is that the downdrafts need to be saturated if pflx is constrained
! to be positive definite - see mod in cldprp. One SHOULD then use some kind of
! kinetic fractionation, but assume this problems is small.... for now).
!
! David Noone <dcn@colorado.edu> - Mon Jul 19 18:26:45 MDT 2004
!
!-----------------------------------------------------------------------
   use cldwat,             only: cldwat_fice
   use water_tracers,      only: wtrc_ratio
   use water_isotopes,     only: wiso_alpl, wiso_akel, wiso_delta, &
                                 wiso_liqvap_equil

   implicit none
!---------------------------- Arguments --------------------------------
   integer , intent(in)    :: lchnk             ! chunk identifier
   integer , intent(in)    :: il1g              ! fist index of column
   integer , intent(in)    :: il2g              ! last index of column
   integer , intent(in)    :: msg               ! lowest level in gather
   integer , intent(in)    :: mx(pcols)         ! level of convection base
   integer , intent(in)    :: jd(pcols)         ! level of downdraft top

   real(r8), intent(in)    :: md(pcols,pver)    ! downdraft masss flux
   real(r8), intent(in)    :: ed(pcols,pver)    ! downdraft entrainment rate
   real(r8), intent(in)    :: dp(pcols,pver)    ! level thichnesss (mb)
!
   real(r8), intent(in)    :: tg(pcols,pver)    ! environment temperature
   real(r8), intent(in)    :: tdg(pcols,pver)   ! downdraft temperature
   real(r8), intent(in)    :: qg(pcols,pver,pcnst+pnats)    ! shum 'qbar'
   real(r8), intent(inout) :: qdnd(pcols,pver,pcnst+pnats)  ! downdraft vapour isotope mass FLUX
   real(r8), intent(inout) :: rrg(pcols,pver,pcnst+pnats)   ! rain rate isotope mass
   real(r8), intent(inout) :: evpg(pcols,pver,pcnst+pnats)  ! evaporation isotope mass
   real(r8), intent(inout) :: rprdg(pcols,pver,pcnst+pnats) ! precip production at each level
   real(r8), intent(inout) :: pflxg(pcols,pverp,pcnst+pnats)! precip flux at each interface

!------------------------- Local Variables -----------------------------
   integer i,k,m

   real(r8) :: teff                             ! effective temperature of phase change
   real(r8) :: hum0                             ! initial humidity in downdraft region
   real(r8) :: alpliq                           ! fractionation factor for liquid/ice
   real(r8) :: vaptot,liqtot                    ! vapour and liquid for prognostic water
   real(r8) :: vapiso,liqiso                    ! vapour and liquid for isotopic water
   real(r8) :: dliqiso                          ! change in isotope liquid from equilibration

   real(r8) :: qtot(pcols,pver,pcnst+pnats)     ! interface non-condensed downdraft vapour
   real(r8) :: qprc(pcols,pver,pcnst+pnats)     ! interface precipitation flux
   real(r8) :: fice (pcols,pver)                ! ice fraction in precip production
   real(r8) :: fsnow(pcols,pver)                ! snow fraction in precip production
   real(r8) :: feq                              ! equilibration fraction
!
! debugging variables
!
  character(len=12) cerr                        ! name of toubled variable
!!  integer  :: idbg = 6
!!  integer  :: idbg = 0                	! no messages, traps only
  integer  :: idbg = -1         		! all messages and traps


#undef DNDIAGS
#undef DBGDNDFT
!-----------------------------------------------------------------------
!
! Determine phase of falling precipitation 
!
   call cldwat_fice(il2g, tdg, fice, fsnow)
!
   qtot(:,:,:) = 0.
!
   do m = ixwti, ixwtx
     if (wtrc_is_vap(m)) then
!!          write(*,*) 'Downdraft for water tracer m=',m
!
! Loop from top down.
!   o Compute the precipitation thoughfall (pflx) assuming no evaporation.
!   o Compute downdraft mixing ratio assuming no evaporation.
!   o If there IS evaporation, resistribute mass to give isotopic equilibration
!     at the lower interface (k+1). 
!   o Mass balance to deduce evaporation
!   o Thus recompute, tracer precipitation flux and downdraft mixing ratio
! (Remember, we are working on k+1, not k as in updrafts)
!
       pflxg(:il2g,1,m) = 0.
       do k = 2,msg+2     ! should be zero above plume... but I'm paranoid.
         do i = 1, il2g
           pflxg(i,k,m) = pflxg(i,k-1,m) + dp(i,k-1)*rrg(i,k-1,m)*100./grav
         end do
       end do
!
! Downdraft initiation
!
       do i = 1, il2g           ! downdraft top is entrained environment air (see cldprp)
         qdnd(i,jd(i),m) = md(i,jd(i))*qg(i,jd(i),m)
       end do
!
       do k = msg+2,pver
         do i = 1,il2g
!
! Process for downdraft, above cloud base (mx... debug, notice this is jb in cldprp)
!
           if (k >= jd(i) .and. k < mx(i) .and. abs(md(i,mx(i))) > msmall) then
             if (k >= mx(i)) then
                write(*,*) '(WTRC_DNDRAFT) K >= mx in downdraft (md < 0).'
             end if
!
! Use an interum guess at the downdraft (no evap), and precipitation flux to 
! compute the total water in system (in a mixing ratio sence).
!
!!debug check -THIS IS CORRECT  
!!debug chekk     qdnd(i,k+1,m) = rmdikp*(md(i,k)*qdnd(i,k,m) - dp(i,k)*(evpg(i,k,1)+ed(i,k)*qg(i,k,1)))
!!                qdnd(i,k+1,m) = rmdikp*(md(i,k)*qdnd(i,k,m) - dp(i,k)*ed(i,k)*qg(i,k,m))
!!                write(*,*) 'QDND^:',qdnd(i,k+1,m),qdnd(i,k+1,1)+rmdikp*dp(i,k)*evpg(i,k,1), &
!!                                    qdnd(i,k+1,m)-qdnd(i,k+1,1)-rmdikp*dp(i,k)*evpg(i,k,1)

! OK                qtot(i,k+1,1) = qdnd(i,k+1,1) - rmdikp*pflxg(i,k+1,1)*grav/100.  ! qdnd has evp mass
! OK                qtot(i,k+1,m) = qdnd(i,k+1,m) - rmdikp*pflxg(i,k+1,m)*grav/100.  ! pflx has evp mass

                qtot(i,k+1,1) = qdnd(i,k,1) - pflxg(i,k,1)*grav/100. - &
                                        dp(i,k)*(ed(i,k)*qg(i,k,1)+rrg(i,k,1))
                qtot(i,k+1,m) = qdnd(i,k,m) - pflxg(i,k,m)*grav/100. - &
                                        dp(i,k)*(ed(i,k)*qg(i,k,m)+rrg(i,k,m))
!
! Compute net precipitation available for evapration (indexed as upper
! interface)
!
                qprc(i,k,1) = -pflxg(i,k,1)*grav/100. - dp(i,k)*rrg(i,k,1)
                qprc(i,k,m) = -pflxg(i,k,m)*grav/100. - dp(i,k)*rrg(i,k,m)
!  
! Initial estimate with no fractionation
!  
                evpg(i,k  ,m) = evpg (i,k,1)*wtrc_ratio(qprc(i,k,m), qprc(i,k,1))
                qdnd(i,k+1,m) = qdnd (i,k,m) - dp(i,k)*(evpg(i,k,m) + ed(i,k)*qg(i,k,m))
!  
! Update the precipitation flux with evaporation estiment
! (thus precip available to equilibration)
!
                qprc(i,k,1) = qprc(i,k,1) + dp(i,k)*evpg(i,k,1)
                qprc(i,k,m) = qprc(i,k,m) + dp(i,k)*evpg(i,k,m)
!
!!                write(*,*)'QDND:',qdnd(i,k+1,m),qdnd(i,k+1,m)-(qtot(i,k+1,m)-qprc(i,k,m))
!!                write(*,*)'TOT1:',qtot(i,k+1,1),qtot(i,k+1,1)-qprc(i,k,1)-qdnd(i,k+1,1)
!!                write(*,*)'TOTM:',qtot(i,k+1,m),qtot(i,k+1,m)-qprc(i,k,m)-qdnd(i,k+1,m)

!
! SOLVE ISOTOPIC FRATIONATION IN DOWNDRAFTS
!
! Set up the temperaure dependent fractionation.
!   o fractionation temperature halfway between in draft and ambient
!   o kinetic effect for liquid. 
!   o sublimation of hail/snow wth no fractionation.
!   o effectve humidity weighted mean of case with no evap, and  saturation
!
! This is only done if the precipitation is part liquid.
!
!!                hum0 = (qdnd(i,k+1,1) + rmdikp*dp(i,k)*evpg(i,k,1)) / max(qdnd(i,k+1,1),qtiny)
!!                hum0 = 1. + mu(i,k+1)*dp(i,k)*evpg(i,k,1)/max(qdnd(i,k+1,1),qtiny)
!!                hum0 = qtot(i,k+1)/(qdnd(i,k+1) - msmall*qtiny)
                if (fsnow(i,k) < 1.0) then
                  hum0 = 1. + dp(i,k)*evpg(i,k,1)/min(qdnd(i,k+1,1),-msmall*qtiny)
!
                  teff = ftc_dnd*tdg(i,k+1) + (1.0-ftc_dnd)*tg(i,k)
                  if (teff > 330) then
                     write(*,*) 'DOWNDRAFT teff > 330',teff
                  end if
                  if (teff < 180) then
                     write(*,*) 'DOWNDRAFT teff < 180',teff ! not impos. trop.strat.
                  end if
!
                  alpliq = wiso_alpl(iwspec(m),teff)
                  alpliq = wiso_akel(iwspec(m),teff,hum0,alpliq)
!
! Compute a modified equilibration fraction to account for
! non-fractionating snow.
!
                  feq = feq_dnd * (1.-fsnow(i,k))
!
! Equlibrate downdraft vapour with liquid precipitate. 
! Snow plays no further role, just use liquid part of precipitation
! flux.
! Here compute vapour, back out the evap, then resolve for total
! precipitation.
! (note minus signs as qfluxes are downward)
!
                  vaptot = -qdnd(i,k+1,1)
                  liqtot = -qprc(i,k  ,1)
                  if (liqtot < 0.) then
                    vaptot = vaptot - liqtot
                    liqtot = 0.
                  end if
!
                  vapiso = -qdnd(i,k+1,m)
                  liqiso = -qprc(i,k  ,m)
                  if (liqiso < 0.) then
                    vapiso = vapiso - liqiso
                    liqiso = 0.
                  end if

                  call wiso_liqvap_equil(alpliq, feq, vaptot, liqtot, &
                                         vapiso, liqiso, dliqiso)
!
                  qdnd(i,k+1,m) = -vapiso
!
                end if          ! liquid precip present (fsnow < 1)
!
! Back out the isotopic evaporation by conservation using downdraft budget equation.
! Also update rprd for consistancy with m=1 (now just a diagnostics).
!
                evpg(i,k,m) = (qdnd(i,k,m)-qdnd(i,k+1,m))/dp(i,k) - ed(i,k)*qg(i,k,m)
                rprdg(i,k,m) = rprdg(i,k,m) - evpg(i,k,m)
!
! Dump some diagnostics
!
#ifdef DNDIAGS
               if (m == 4) then
                 write(*,*) 'DNFRAC:',qdnd(i,k+1,1)/min(qtot(i,k+1,1),-msmall*qtiny)
                 write(*,*) 'DNCNDS:',teff,hum0

                 write(*,4) 'WT_ZM-DN QTDEpP:',1, &
                   qg   (i,k,1)*1000., &
                   -qtot (i,k+1,1)*1000., &
                   -qdnd (i,k+1,1)*1000., &
                   -qprc (i,k+1,1)*1000., &
                   evpg (i,k,1)*dp(i,k)*1000., &
                    pflxg(i,k,1)                      *1000.*grav/100.,&
                   (pflxg(i,k,1)+dp(i,k)*rprdg(i,k,1))*1000.*grav/100
!                   (pflxg(i,k,1)+dp(i,k)*(rrg(i,k,1)-evpg(i,k,1)))*(-rmdikp)*1000.*grav/100
!                   rprdg(i,k,1)*dp(i,k)*(-rmdikp)*1000., &
                 write(*,4) 'WT_ZM-D dQTDEpP:',m, &
                   wiso_delta(iwspec(m),qg   (i,k,m), qg   (i,k,1)), &
                   wiso_delta(iwspec(m),-qtot (i,k+1,m), -qtot (i,k+1,1)), &
                   wiso_delta(iwspec(m),-qdnd (i,k+1,m), -qdnd (i,k+1,1)), &
                   wiso_delta(iwspec(m),-qprc (i,k+1,m), -qprc (i,k+1,1)), &
                   wiso_delta(iwspec(m),evpg (i,k,m), evpg (i,k,1)), &
                   wiso_delta(iwspec(m),pflxg(i,k,m), pflxg(i,k,1)), &
                   wiso_delta(iwspec(m),pflxg(i,k,m)+dp(i,k)*rprdg(i,k,m), &
                                        pflxg(i,k,1)+dp(i,k)*rprdg(i,k,1))
  4              format(a16,i3,7f12.5)
               end if
#endif

#ifdef DBGDNDFT
               if (m==4) then
                 if (i == idbg .or. idbg < 0) then
                  write(*,*) 'DNDFT i,k:',i,k
                  write(*,*) 'FSNOW :',fsnow(i,k)
                  write(*,*) 'MD k,+:',  md(i,k), md(i,k+1), ed(i,k)
                  write(*,1) 'QTOT +:', qtot(i,k+1,m),qtot(i,k+1,1),qtot(i,k+1,m)-qtot(i,k+1,1)
                  write(*,1) 'QDND  :',  qdnd(i,k  ,m),qdnd(i,k ,1),qdnd(i,k  ,m)-qdnd(i,k  ,1)
                  write(*,1) 'QPRC  :',  qprc(i,k  ,m),qprc(i,k ,1),qprc(i,k  ,m)-qprc(i,k  ,1)
                  write(*,1) 'QDND +:', qdnd(i,k+1,m),qdnd(i,k+1,1),qdnd(i,k+1,m)-qdnd(i,k+1,1)
                  write(*,1) 'QG    :',  qg  (i,k  ,m),qg  (i,k  ,1),qg (i,k  ,m)-qg  (i,k  ,1)
                  write(*,1) 'RRG   :',  rrg (i,k  ,m),rrg (i,k  ,1),rrg (i,k  ,m)-rrg (i,k  ,1)
                  write(*,1) 'EVPG  :',  evpg(i,k  ,m),evpg(i,k ,1),evpg(i,k  ,m)-evpg(i,k  ,1)
                  write(*,1) 'PFXG  :',  pflxg(i,k ,m),pflxg(i,k ,1),pflxg(i,k ,m)-pflxg(i,k ,1)
                 endif

                 cerr = 'none'
                 if (abs ( qdnd(i,k+1,m)-qdnd(i,k+1,1)) > 1.e-12) cerr = 'qdnd'
                 if (abs ( evpg(i,k  ,m)-evpg(i,k  ,1)) > 1.e-12) cerr = 'evpg'

                 if (cerr /= 'none') then
                    write(*,*) '(wtrc_zm_downdraft) conservation error for: '//trim(cerr),i,k
                    call endrun('wtrc_zm_downdraft: debugging error')
                 endif

  1               format(a8,3e16.6)
               end if
#endif
!
           else         ! just a debug check (delete these lines when working)
             if (evpg(i,k,1) /= 0.) then
                write(*,*) 'wtrc_zm_downdraft EVP /= 0  in downdraft.', evpg(i,k,1)
             end if
           end if    ! k in downdraft
!
! Update precipitation flux for next level.
!
           pflxg(i,k+1,m) = pflxg(i,k,m) + dp(i,k)*(rrg(i,k,m)-evpg(i,k,m))*100./grav
#ifdef DBGDNDFT
           if (m == 4) then
             if (abs(pflxg(i,k+1,1)-pflxg(i,k+1,m)) > 1.e-12) then
               write(*,*) 'PFLXG:',k,pflxg(i,k+1,1),pflxg(i,k+1,m), &
                                     pflxg(i,k+1,1)-pflxg(i,k+1,m)
               write(*,*) 'WTRC_DNDRAFT - inconsistent pflx.'
!              call endrun
             endif
           end if
!
#endif
         end do      ! i, columns
       end do        ! k, levels 
!
     end if          ! m is water
   end do            ! m, tracers
!
! Do a fix for downdraft fluxes as they can be imprecise
! Assume deficit has ratio of environment. This is a "qneg" type thing.'
! One can imagine this is a bit more entrainment.
! This is really just to ensure the qchks pass, and must be done after "m" loop.
!
!!   if (lupdcorr) then
!!     call updraft_correction('QDND',il2g,pver ,1,pver ,1,ixh2oq,qdnd,qg) ! for dqdt
!!!!     call updraft_correction('EVPG',il2g,pver,1,pver,1,ixh2oq,evpg,qg) ! for dqemc
!!     call updraft_correction('PFLX',il2g,pverp,1,pverp,1,ixh2oq,pflxg,qg) ! for pflx
!!   endif

! Check conservation of total: qdnd (qu), qevp (as ev) 
! (numerical imprecision accumulates in pflx as is based on both updraft and
! downdraft mass integrals, updraft isotope iteration for rrg, and downdraft
! isotope equilibration in evpg - thus, be a bit generous on precision
! requirement)
!
    call wtrc_qchk2('WTRC_ZMCcldd','qtot',il2g,qtot (:,:,ixh2oq) ,qtot (:,:,1), 1.e-14)	! eeerh!
    call wtrc_qchk2('WTRC_ZMCcldd','qdnd',il2g,qdnd (:,:,ixh2oq) ,qdnd (:,:,1), 1.e-15)
    call wtrc_qchk2('WTRC_ZMCcldd','evpg',il2g,evpg (:,:,ixh2oq) ,evpg (:,:,1))
    call wtrc_qchk2('WTRC_ZMCcldd','pflx',il2g,pflxg(:,2:pverp,ixh2oq),pflxg(:,2:pverp,1), 1.e-14)
!
   return
  end subroutine wtrc_zm_dndraft

!=======================================================================
  subroutine wtrc_zm_convr_q1q2( lchnk   ,il1g    ,il2g    ,  &
               msg     ,jt      ,mx      ,dp      ,dsubcld ,  &
               mu      ,du      ,md      ,qg      , &
               qhup    ,qhdn    ,qupd    ,qliq    ,qice    ,qdnd    ,  &
               cug     ,evpg    ,dqmdt   ,dqemc   ,dlmdt   ,  & 
               dimdt   ,ldetg   ,idetg)

!-----------------------------------------------------------------------
!
! Purpose:  Compute mass flux tendancies (like q1q2_pjr)
!
! Method:
!
!   As per q1q2 for total water, but for isotopes. This is easy
!   as we have dont the hard work above getting the ratios right.
!   This includes the mass transport and the condensation terms.
!   No checking for mass "over removal" is done. This is assumed 
!   to be a "numerical" problem associated with the time step (i.e.
!   unstable diffusive equation). Handled elsewhere explictly.
!
! Author: 
!
!    David Noone <dcn@caltech.edu> - Tue Oct 28 13:28:26 MST 2003
!
!---------------------------- Arguments --------------------------------
   integer , intent(in)   :: lchnk                ! chunk identifier
   integer , intent(in)   :: il1g                 ! fist index of column
   integer , intent(in)   :: il2g                 ! last index of column
   integer , intent(in)   :: msg                  ! lowest level in gather
   integer , intent(in)   :: jt(pcols)            ! level of convection top
   integer , intent(in)   :: mx(pcols)            ! level of convection base
!
   real(r8), intent(in)   :: dp(pcols,pver)       ! layer thickness
   real(r8), intent(in)   :: dsubcld(pcols)       ! sub cloud thichness
   real(r8), intent(in)   :: mu(pcols,pver)       ! mass flux of updraft
   real(r8), intent(in)   :: du(pcols,pver)       ! detrainment rate of updraft
   real(r8), intent(in)   :: md(pcols,pver)       ! mass flux of downdraft
!
   real(r8), intent(in)   :: qg(pcols,pver,pcnst+pnats)    ! environment q
   real(r8), intent(in)   :: qhup(pcols,pver,pcnst+pnats)  ! large scale interface hsum FLUX
   real(r8), intent(in)   :: qhdn(pcols,pver,pcnst+pnats)  ! large scale interface hsum FLUX
   real(r8), intent(in)   :: qupd(pcols,pver,pcnst+pnats)  ! updraft vapour isotope ratio FLUX
   real(r8), intent(in)   :: qliq(pcols,pver,pcnst+pnats)  ! updraft liquid isotope ratio FLUX
   real(r8), intent(in)   :: qice(pcols,pver,pcnst+pnats)  ! updraft liquid isotope ratio FLUX
   real(r8), intent(in)   :: qdnd(pcols,pver,pcnst+pnats)  ! downdraft vapour isotope ratio FLUX
   real(r8), intent(in)   :: cug(pcols,pver,pcnst+pnats)   ! rate of condensation in updraft
   real(r8), intent(in)   :: evpg(pcols,pver,pcnst+pnats)  ! rate of evaporation into downdraft
   real(r8), intent(in)   :: ldetg(pcols,pver,pcnst+pnats) ! detraining liquid isotope mass
   real(r8), intent(in)   :: idetg(pcols,pver,pcnst+pnats) ! detraining ice isotope mass

   real(r8),intent(inout) :: dqmdt(pcols,pver,pcnst+pnats) ! tracer vapour tendancies
   real(r8),intent(inout) :: dqemc(pcols,pver,pcnst+pnats) ! tendancy due to E-C
   real(r8),intent(inout) :: dlmdt(pcols,pver,pcnst+pnats) ! tracer vapour tendancies
   real(r8),intent(inout) :: dimdt(pcols,pver,pcnst+pnats) ! tracer vapour tendancies

!------------------------- Local Variables -----------------------------
   integer i,k,m                        ! array indicies
   integer kbm                          ! k index for convection bottom
   integer ktm                          ! k index for convection top
!-----------------------------------------------------------------------
!
! find the highest level top and bottom levels of convection
!
   ktm = pver
   kbm = pver
   do i = 1, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do
!
! Do computation with respect to the vapour so as to get the
! flux forms exactly equal to the total water
!
   do m = ixwti,ixwtx
     if (wtrc_is_vap(m)) then
!
       do k = msg + 1,pver
          do i = 1,il2g
             dqemc(i,k,m) = 0.
             dqmdt(i,k,m) = 0.
             dlmdt(i,k,m) = 0.
             dimdt(i,k,m) = 0.
          end do
       end do
!
! Tendancies for levels above the cloud base layers
!
       do k = ktm,pver-1
          do i = 1,il2g
!
             dqemc(i,k,m) = evpg(i,k,m) - cug (i,k,m)

             dqmdt(i,k,m) = dqemc(i,k,m) + &
                    ( (qupd(i,k+1,m)-qhup(i,k+1,m)) - (qupd(i,k  ,m)-qhup(i,k  ,m)) &
                     +(qdnd(i,k+1,m)-qhdn(i,k+1,m)) - (qdnd(i,k  ,m)-qhdn(i,k  ,m)) )/dp(i,k)

             dlmdt(i,k,m) = ldetg(i,k,m)	! fractionated in plume
             dimdt(i,k,m) = idetg(i,k,m)	! fractionated in plume

          end do
       end do
!
! Tendancies in the subcloud layer are one-sided
!
       do k = kbm,pver
          do i = 1,il2g
             if (k == mx(i)) then
                dqmdt(i,k,m) = (1./dsubcld(i))* &
                        (-(qupd(i,k,m)-qhup(i,k,m)) - (qdnd(i,k,m)-qhdn(i,k,m)))
             else if (k > mx(i)) then
                dqmdt(i,k,m) = dqmdt(i,k-1,m)
             end if
          end do
       end do
!
     end if          ! m is water vapour tracer
   end do            ! m tracers
!
! If these has been a qu < 0 output from cldprp, the tendencies can be
! different. To fix this apply a bute force correction, preserving as
! much mass as possible
!
   if (ltndcorr) then
     call tendency_correction('DQDT',il2g,1,ixh2oq,dqmdt,qg)
   end if
!
! Conservation check
!
   call wtrc_qchk2('WTRC_ZMCq1q2','dEmC',il2g,dqemc(:,:,ixh2oq) ,dqemc(:,:,1))
   call wtrc_qchk2('WTRC_ZMCq1q2','dql ',il2g,dlmdt(:,:,ixh2oq) ,dlmdt(:,:,1))
   call wtrc_qchk2('WTRC_ZMCq1q2','dqdt',il2g,dqmdt(:,:,ixh2oq) ,dqmdt(:,:,1))
!
   return
  end subroutine wtrc_zm_convr_q1q2

!=======================================================================
  subroutine updraft_correction(mesg,il2g,nlevd,klo,khi,mone,mbase,qfld,qenv)
!-----------------------------------------------------------------------
! Applies a corerction when tracer diverges from total water
! (which occuers when qu < 0 is output from zm_scheme)
!-----------------------------------------------------------------------
  use water_tracers, only: wtrc_ratio
!-----------------------------------------------------------------------
  implicit none
!---------------------------- Arguments --------------------------------
  character(len=*), intent(in) :: mesg	! a message
  integer , intent(in) :: il2g
  integer , intent(in) :: nlevd		! level dimension (pver, pverp or 1)
  integer , intent(in) :: klo, khi
  integer , intent(in) :: mone
  integer , intent(in) :: mbase
  real(r8), intent(in) :: qenv(pcols,nlevd,pcnst+pnats)
!
  real(r8), intent(inout) :: qfld(pcols,nlevd,pcnst+pnats)
!------------------------- Local Variables -----------------------------
  integer i,k,m
  real(r8) :: qerr(pcols,nlevd)
  real(r8) :: Rwrk
!-----------------------------------------------------------------------
   qerr(:,:) = 0.
!
! Compute the error term
!
   qerr(:il2g,klo:khi) = qfld(:il2g,klo:khi,mbase) - qfld (:il2g,klo:khi,mone)
!
! Apply error correction for all tracers
!
   do m = ixwti, ixwtx
     if (wtrc_is_vap(m)) then
        do k = klo,khi
          do i = 1, il2g
!
             if (abs(qerr(i,k)) > qtiny) then
               if (qerr(i,k) > 0.) then
                 Rwrk = wtrc_ratio(qfld(i,k,m),qfld(i,k,mbase))
                 if (abs(qfld(i,k,mbase)) < qtiny) then		!  nothing there, suck it from environment
!!                   write(*,*) '(wtrc_zm_updraft_correction) using environment.'
                   Rwrk = wtrc_ratio(qenv(i,k,m),qenv(i,k,mbase))
                 end if
                 if (abs(qerr(i,k)) > 1.e-14) then
                   write(*,*) '(wtrc_zm_updraft_correction) Correcting excess:',i,k
                   write(*,*) trim(mesg)//':',qerr(i,k),Rwrk
                 end if
               else
                 Rwrk = wtrc_ratio(qenv(i,k,m),qenv(i,k,mbase))
                 if (abs(qerr(i,k)) > 1.e-14) then
                   write(*,*) '(wtrc_zm_updraft_correction) Correcting deficit:',i,k
                   write(*,*) trim(mesg)//':',qerr(i,k),Rwrk
                 end if
               end if
               qfld(i,k,m) = qfld(i,k,m) - Rwrk*qerr(i,k)
             end if

          end do
        end do
     end if
   end do
!
  return
  end subroutine updraft_correction

!=======================================================================
  subroutine tendency_correction(mesg,il2g,mone,mbase,dqdt,qenv)
!-----------------------------------------------------------------------
! Corrects tracer tendencies: Assume computed ratio of tendency is
! correct.  As this will be assed prior to computation of of net
! precipitation, we have infact managed to conserve all mass!
!-----------------------------------------------------------------------
  use water_tracers, only: wtrc_ratio
!-----------------------------------------------------------------------
  implicit none
!---------------------------- Arguments --------------------------------
  character(len=*), intent(in) :: mesg	! a message
  integer , intent(in) :: il2g
  integer , intent(in) :: mone
  integer , intent(in) :: mbase
  real(r8), intent(in) :: qenv(pcols,pver,pcnst+pnats)
!
  real(r8), intent(inout) :: dqdt(pcols,pver,pcnst+pnats)
!------------------------- Local Variables -----------------------------
  integer i,k,m
  real(r8) :: qerr(pcols,pver)
  real(r8) :: Rwrk
!-----------------------------------------------------------------------

   qerr(:il2g,1:pver) = dqdt(:il2g,1:pver,mbase) - dqdt(:il2g,1:pver,mone)
!
! Apply the correction
!
   do m = ixwti, ixwtx
     if (wtrc_is_vap(m)) then
        do k = 1,pver
          do i = 1, il2g

             if (abs(qerr(i,k)) > qtiny) then
                 Rwrk = wtrc_ratio(dqdt(i,k,m),dqdt(i,k,mbase))
                 if (Rwrk < 0. .or. abs(dqdt(i,k,mone)) < qtiny) then
                   Rwrk = wtrc_ratio(qenv(i,k,m),qenv(i,k,mbase))
                   if (abs(qerr(i,k)) > 1000.*qtiny) then	! warn if big (1.e-15, is qmin(1)*dtime)
                     write(*,*) '(wtrc_zm_tendency_correction) Correcting excess WITH ENVIRONMENT:',i,k
                     write(*,*) trim(mesg)//':',qerr(i,k),Rwrk
                   endif
                 else
                   if (abs(qerr(i,k)) > 1000.*qtiny) then	! warn if big (1.e-15, is qmin(1)*dtime)
                     write(*,*) '(wtrc_zm_tendency_correction) Correcting excess:',i,k
                     write(*,*) trim(mesg)//':',qerr(i,k),Rwrk
                   endif
                 end if

                 dqdt(i,k,m) = dqdt(i,k,m) - Rwrk*qerr(i,k)
             end if

          end do
        end do
     end if
   end do
   
  return
  end subroutine tendency_correction

!=======================================================================
end module wtrc_zm_conv

