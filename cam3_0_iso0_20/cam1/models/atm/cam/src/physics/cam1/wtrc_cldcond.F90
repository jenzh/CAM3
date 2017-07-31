#include <misc.h>
#include <params.h>

module wtrc_cldcond

!-----------------------------------------------------------------------
!
!  Purpose:
!
!    Provides the CAM interface to the isotopic versions of prognostic 
!    cloud water and ice parameterization
!
!  Author:
!
!    David Noone <dcn@caltech.edu> - Tue Oct 28 14:48:45 MST 2003
!
!-----------------------------------------------------------------------

  use shr_kind_mod,  only: r8=>shr_kind_r8
  use abortutils,    only: endrun
  use ppgrid,        only: pcols, pver

!-----------------------------------------------------------------------
  implicit none
  save
  private
!-----------------------------------------------------------------------
!
! public interfaces 
!
  public :: wtrc_cldcond_init           ! initializes module and history fields
  public :: wtrc_cldcond_tend           ! computes large scale tendancies
  public :: wtrc_cldcond_zmconv_detrain ! handles detained convective liquid
  public :: wtrc_cldcond_sediment       ! computes sedimentation
!
!=======================================================================
CONTAINS

!=======================================================================
subroutine wtrc_cldcond_init
!-----------------------------------------------------------------------
! initialize module and history files access
!-----------------------------------------------------------------------
  use water_tracers, only: wtrc_is_vap, ixwti, ixwtx
  use constituents,  only: cnst_name
  use history,       only: addfld, add_default, phys_decomp

  implicit none
!-----------------------------------------------------------------------
  integer m
!-----------------------------------------------------------------------

  do m = ixwti, ixwtx
    if (wtrc_is_vap(m)) then

! Register history variables
      call addfld ('DLF'//trim(cnst_name(m)),'m/s     ',pver, 'A',trim(cnst_name(m))//' detrained could water',phys_decomp)
      call addfld ('DIF'//trim(cnst_name(m)),'m/s     ',pver, 'A',trim(cnst_name(m))//' detrained could water',phys_decomp)

! Set defaults for primary history file
      call add_default ('DLF'//trim(cnst_name(m)), 1, ' ')
      call add_default ('DIF'//trim(cnst_name(m)), 1, ' ')
!
    end if
  end do
!
  return
end subroutine wtrc_cldcond_init

!=======================================================================
subroutine wtrc_cldcond_sediment(state, ptend, dtime, &
                                 cloud, fxliq, fxice, prec, snow)
!-----------------------------------------------------------------------
!
! Purpose: Perform vertical gravitation transport of ice and liquid
!
! Method:
!   Does what "cldcond_sediment" does, for for an arbitrary number of
!   water tracers.This routine should be called immediately after 
!   cldcond_sediment. physics_update flags set for tracer quantities.
!
!   Notice that no fractionation is applied here - it is just gravitational.
!   After the sedimentation, cloid liquid water is equilibrated wit the vapour
!   during the large scale condensation, so all is good at the END of the step.
!
! ASSUME TRACER ORDERED AS TRIPPLETS: vapour(m), liquid(m+1), ice(m+2).
!
! Code history:
!   Based on "Sedimentation Version 2" from Andrew Gettelman,
!   "cldcond_sediment" from B.A. Boville and P.J. Rasch.
!   David Noone <dcn@caltech.edu> - Thu Aug  7 15:53:25 PDT 2003
!   Sun Wong  <swong@atmos.umd.edu> - (2/17/2003) for relaxation to
!                                              equilibrium for liq.        
!
! Author: (cam3)
!    David Noone <dcn@colorado.edu> - Fri Jul  2 18:57:10 MDT 2004
!
!-----------------------------------------------------------------------
!
  use shr_kind_mod,  only: r8=>shr_kind_r8
  use physics_types, only: physics_state, physics_ptend
  use ppgrid,        only: pcols, pver, pverp
  use constituents,  only: pcnst, pnats, cnst_get_ind
  use history,       only: outfld
  use water_tracers, only: trace_water, wtrc_is_vap, ixwti, ixwtx

!-----------------------------------------------------------------------
  implicit none
!---------------------------- Arguments --------------------------------
!
    type(physics_state), intent(in)    :: state   ! state variables (pre adjusted)
    type(physics_ptend), intent(inout) :: ptend   ! package tendencies

    real(r8), intent(in)   :: dtime                ! timestep
    real(r8), intent(in)   :: cloud(pcols,pver)    ! cloud fraction

    real(r8), intent(in)   :: fxliq(pcols,pverp)   ! liq sed flux at interfaces (+ down)
    real(r8), intent(in)   :: fxice(pcols,pverp)   ! ice sed dlux at interfaces (+down)

    real(r8), intent(inout):: prec(pcols,pcnst+pnats)  ! surface flux of total cloud water
    real(r8), intent(inout):: snow(pcols,pcnst+pnats)  ! surface flux of cloud ice
!
!------------------------- Local Variables -----------------------------
!
    integer  :: m                 	! loop indexes
    integer  :: ncol                    ! number columns in chunk
    integer  :: lchnk                   ! chunk index
    integer  :: ixcldliq,ixcldice	! cloud liquid and ice indicies
    integer  :: ixwtvap,ixwtliq,ixwtice ! index of water tracer triplets

    real(r8) :: rain(pcols,pcnst+pnats) ! surface flux of cloud liquid
!
!-----------------------------------------------------------------------
   ncol = state%ncol
   
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
!
! Compute the sedimentation tendancies.
! LOOP OVER TRIPLETS - assumes order vap, liq, ice (MUST BE REGESTERED)
!
   do m = ixwti,ixwtx
     if (wtrc_is_vap(m)) then
       ixwtvap = m
       ixwtliq = m+1
       ixwtice = m+2
!
       ptend%lq(ixwtvap) = .TRUE.
       ptend%lq(ixwtliq) = .TRUE.
       ptend%lq(ixwtice) = .TRUE.
!
       call wtrc_cld_sediment_tend (ncol  ,ixwtvap  ,dtime  , &
                 state%pmid   ,state%pdel ,cloud ,  &
                 state%q(:,:,      1), state%q(:,:,ixcldliq), state%q(:,:,ixcldice), &
                 state%q(:,:,ixwtvap), state%q(:,:,ixwtliq) , state%q(:,:,ixwtice) , &
                 fxliq ,fxice , &
                 ptend%q(:,:,ixwtliq), ptend%q(:,:,ixwtice), ptend%q(:,:,ixwtvap),  &
                 rain(:,ixwtvap) , snow(:,ixwtvap)  )
     end if
   end do
!
! convert rain and snow from kg/m2 to m/s, and sum for precip
!
   do m = ixwti,ixwtx
     if (wtrc_is_vap(m)) then
       snow(:ncol,m) = snow(:ncol,m)/1000.
       rain(:ncol,m) = rain(:ncol,m)/1000.
       prec(:ncol,m) = rain(:ncol,m) + snow(:ncol,m)
     end if
   end do
!
! Dump interesting things to history file
! (need to think of sensible 8 char naming convention....)
!
!!   lchnk = state%lchnk
!!   do m = ixwti,ixwtx
!!     if (wtrc_is_vap(m)) then
!!       ixwtvap = m
!!       ixwtliq = m+1
!!       ixwtice = m+2
!!       call outfld('WTRCDQSED'   ,ptend%q(:,:,ixwtvap) , pcols,lchnk)
!!       call outfld('WTRCDISED'   ,ptend%q(:,:,ixwtliq)  , pcols,lchnk)
!!       call outfld('WTRCDLSED'   ,ptend%q(:,:,ixwtice  , pcols,lchnk)
!!       call outfld('PSDWTRCPRECSED' ,prec(:,m)         , pcols,lchnk)
!!       call outfld('PSDWTRCSNOWSED' ,snow(:,m)         , pcols,lchnk)
!!       call outfld('PSDWTRCRAINSED' ,rain(:,m)         , pcols,lchnk)
!!     end if
!!   end do
!
  return
end subroutine wtrc_cldcond_sediment


!=======================================================================
subroutine wtrc_cldcond_tend(state, ptend, dtime, Rtrc, dqlfz, &
              evapprec, prodprec, evapsnow, prodsnow, trprec, trsnow) 

!-----------------------------------------------------------------------
!
! Purpose: impliments water [isotope] tracers for in-cloud, 
!          large scale condensation using a bulk method.
!
!   THIS ROUTINE DOES NOT INCLUDE ALL MICROPHYSICS INFORMATION 
!
!
! Method:
!
!   Based on output of cldcond_tend, and in particular, pcond.
!   After the adjustment, any cloud and rain liquid is in isotopic 
!   equilibrium with the vapour. 
!
!   Isotopic calulations essentially follow the Ciais and Jouzel
!   "Mixed Cloud Isotopic Model" (also by Arnold and Gedzelman, and
!   Federer at al.). Thus we have for some layer of interest some
!   total amount of vapour, liquid and ice. A liquid and ice source
!   exists at the upper boundary. Thus the mass conservation
!   (given the Vapour, Liquid, iCe, Rain and Snow) is:
!
!      v' + l' + c' ( + r' + s') = v + l c + r + s
!
!     OR             dv + dl + dc + dr + ds = 0
!     For isotopes   dvi + dli + dci + dri + dsi = 0
!
!
!   This treats the bulk phases assuming initially all snow mixes with
!   cloud ice, and all falling rain mixes with cloud liquid. This allows
!   solution for the change in the vapour due to distillation to ice
!   while isotopic equilibrium is maintained. Thus the final ice and
!   liquid can be repartitioned into the snow and rain components.
!
!   Because of these simplifications, there are some substantial science 
!   limitatations. Specifically we no longer have full corroberation
!   between the microphysical exchanges and the the large scale
!   (e.g., can not have snow evaporation and ice formation... etc)
!   Similarly we fail to address questions: does the "cloud" vapour/liquid 
!   interact with the  non-cloudy part? Does the precipitation form at 
!   layer k "mix" with that from layer k-1? Which ice is evaporated? 
!   Incloud? Falling? Hopefully some of these will play out in version 2
!   of the isotope scheme, but I'm guessing it will show the
!   microphysics is wrong, as is. 
!
!
! Code history:
!     David Noone <dcn@caltech.edu> - Mon Oct 20 16:44:01 PDT 2003
!     Sun Wong <swong@atmos.umd.edu> - Thu Nov 6 14:52 PDT 2003
!     David Noone <dcn@colorado.edu> - Sat Jul  3 12:36:38 MDT 2004 (cam3)
!
!----------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use physics_types,   only: physics_state, physics_ptend
  use ppgrid,          only: pcols, pver, pverp
  use physconst,       only: gravit, cpair, tmelt
  use constituents,    only: pcnst, pnats, cnst_get_ind
  use history,         only: outfld
  use water_tracers,   only: trace_water, wtrc_is_vap, iwspec, ixwti, ixwtx, wtrc_ratio
  use water_isotopes,  only: wisotope, wiso_dicm, wiso_delta
  use wtrc_camtune,    only: feq_cld, fsnk_cld
!
!-----------------------------------------------------------------------
   implicit none
!---------------------------- Arguments --------------------------------
   type(physics_state), intent(inout) :: state   ! state variables
   type(physics_ptend), intent(inout) :: ptend   ! package tendencies
   real(r8), intent(in)  :: dtime                ! timestep

   real(r8), intent(in)  :: Rtrc(pcols,pver,pcnst+pnats)  ! pre-repartitioned tracer ratio
   real(r8), intent(in)  :: dqlfz(pcols,pver)    ! cloud liquid change due to freeze (repartitioning)
   real(r8), intent(in)  :: prodprec(pcols,pver) ! rate of conversion to precip (1/s)
   real(r8), intent(in)  :: prodsnow(pcols,pver) ! rate of conversion to snow (1/s)
   real(r8), intent(in)  :: evapprec(pcols,pver) ! rate of evaporation of precip (1/s)
   real(r8), intent(in)  :: evapsnow(pcols,pver) ! rate of evaporation of snow (1/s)

   real(r8), intent(out) :: trprec(pcols,pcnst+pnats)  ! sfc flux of precip (m/s)
   real(r8), intent(out) :: trsnow(pcols,pcnst+pnats)  ! sfc flux of snow   (m/s)

!------------------------- Local Variabless ----------------------------
   integer i,k,m, kk			! array indicies
   integer ncol				! number of columns to process
   integer lchnk			! chunk identifier for history
   integer ixcldliq,ixcldice		! index to prognostic waters
   integer ixwtvap,ixwtliq,ixwtice	! index to tracer waters
!
   real(r8) :: qtttl(pcols,pver)	! total water (v+l+i+r+s)
   real(r8) :: qtvapold(pcols,pver)	! "old" vapour
   real(r8) :: qtliqold(pcols,pver)	! "old" net liquid (cloud + rain)
   real(r8) :: qticeold(pcols,pver)	! "old" net ice (cloud + snow)
   real(r8) :: qtvap(pcols,pver)	! "new" vapour
   real(r8) :: qtliq(pcols,pver)	! "new" net liquid (cloud + rain)
   real(r8) :: qtice(pcols,pver)	! "new" net ice (cloud + snow)
   real(r8) :: qtrnw(pcols,pver)	! rain water leaving layer
   real(r8) :: qtsnw(pcols,pver)	! snow water leaving layer

   real(r8) :: precipt(pcols) 		! cumulative precip from above
   real(r8) :: rainabt(pcols) 		! cumulative rain from above
   real(r8) :: snowabt(pcols)		! cumulative snow from above
!
   real(r8) :: trttl(pcols,pver)	! total tracer water (v+l+i+r+s)
   real(r8) :: trvapold(pcols,pver)	! "old" tracer vapour
   real(r8) :: trliqold(pcols,pver)	! "old" tracer net liquid (cloud + rain)
   real(r8) :: triceold(pcols,pver)	! "old" tracer net ice (cloud + snow)
   real(r8) :: trvap(pcols,pver)	! "new" tracer vapour
   real(r8) :: trliq(pcols,pver)	! "new" tracer net liquid (cloud + rain)
   real(r8) :: trice(pcols,pver)	! "new" tracer net ice (cloud + snow)
   real(r8) :: trrnw(pcols,pver)	! tracer rain water leaving layer
   real(r8) :: trsnw(pcols,pver)	! tracer snow water leaving layer

   real(r8) :: trprecip(pcols) 		! tracer cumulative precip from above
   real(r8) :: trrainab(pcols) 		! tracer cumulative rain from above
   real(r8) :: trsnowab(pcols)		! tracer cumulative snow from above

   real(r8) :: dqrepart(pcols,pver)	! repartition change in liquid and ice
!
   real(r8) :: fnlsk			! fraction of net liquid going to sink
   real(r8) :: fnisk			! fraction of net ice going to sink
!  
   real(r8) :: qtiny = 1.e-18           ! tiny q for conservation check
!
! Variables passed in/out of fractionation sub-model
!
   integer, parameter :: piso = 2       ! solve for one species at a time
   integer  :: isp(piso)                ! species indicies (1 is water, 2 is isotope)
   real(r8) :: told, tnew		! old and new temperature
   real(r8) :: vapold(piso)             ! initial vapour
   real(r8) :: liqold(piso)             ! initial liquid
   real(r8) :: iceold(piso)             ! initial ice
   real(r8) :: vapent(piso)             ! vapour entrainment
   real(r8) :: vapnew(piso)             ! final vapour
   real(r8) :: liqnew(piso)             ! final liquid
   real(r8) :: icenew(piso)             ! final ice
   real(r8) :: vapdet(piso)             ! vapour detrainment
   real(r8) :: rainpr(piso)             ! rain production
   real(r8) :: snowpr(piso)             ! snow production
   real(r8) :: dliqmt                   ! liquid change due to ice melt
   real(r8) :: dicefz                   ! ice change due to liquid freeze
!
!-----------------------------------------------------------------------
#ifdef PERGRO
   write(6,*) 'wtrc_cldond_tend: This does not allow PERGRO. See pcond.'
   call endrun
#endif

   ncol  = state%ncol
   lchnk = state%lchnk
!
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
!
! Repartition liquid and ice due to in-cloud freeze/melt.
!  [this is numerically robust as Rtrc = qtrc/(qtot + qtiny)]
!   Loop over tripplets of water variables (vapour, liquid, ice)
!
   do m = ixwti,ixwtx 
     if (wtrc_is_vap(m)) then
       ixwtvap = m              !
       ixwtliq = m+1            !=- assume triplets are ordered v, l, i
       ixwtice = m+2            !
!
       where (dqlfz(:ncol,:pver) < 0.)          ! freeze dR(liq) liquid->ice
         dqrepart(:ncol,:pver) = Rtrc(:ncol,:pver,ixwtliq)*dqlfz(:ncol,:pver)
       elsewhere                                ! melt dR(ice) ice->liquid
         dqrepart(:ncol,:pver) = Rtrc(:ncol,:pver,ixwtice)*dqlfz(:ncol,:pver)
       endwhere
       state%q(:ncol,:pver,ixwtliq) = state%q(:ncol,:pver,ixwtliq) + dqrepart(:ncol,:pver)
       state%q(:ncol,:pver,ixwtice) = state%q(:ncol,:pver,ixwtice) - dqrepart(:ncol,:pver)
     end if
   end do

!
! Compute the "new" values, as would be found from physics_update
!  notice dqvap+dqliq+dqice + qrain(k-1) + qsno(k-1) = qrain+qsnow
!  The total is only the water involved in the transitions, and could
!  be limited to the falling precip, or just the vapour in the cloudy fraction.
!
   qtrnw(:ncol,:) = 0.0
   qtsnw(:ncol,:) = 0.0
   precipt(:ncol) = 0.0
   rainabt(:ncol) = 0.0
   snowabt(:ncol) = 0.0
!
   do k = 1, pver
      do i = 1, ncol
        if (state%t(i,k) >= tmelt) snowabt(i) = 0.        
        rainabt(i) = precipt(i) - snowabt(i)
!
! "old" values: mix in situ cloud, and precip from above
!
        qtvapold(i,k) = state%q(i,k,       1)
        qtliqold(i,k) = state%q(i,k,ixcldliq) + dtime*rainabt(i)*gravit/state%pdel(i,k)
        qticeold(i,k) = state%q(i,k,ixcldice) + dtime*snowabt(i)*gravit/state%pdel(i,k)
!
! Check for rounding
!
        if (qtvapold(i,k) < -100.*qtiny) &
               write(*,*) '(wtrc_cldcond_trend) qtvapold < qtiny:',qtvapold(i,k)
        qtvapold(i,k) = max(qtvapold(i,k), 0.)
        if (qtliqold(i,k) < -100.*qtiny) &
               write(*,*) '(wtrc_cldcond_trend) qtliqold < qtiny:',qtliqold(i,k)
        qtliqold(i,k) = max(qtliqold(i,k), 0.)
        if (qticeold(i,k) < -100.*qtiny) &
               write(*,*) '(wtrc_cldcond_trend) qticeold < qtiny:',qticeold(i,k)
        qticeold(i,k) = max(qticeold(i,k), 0.)
!
! "new" values: to a "physics_update"
!
        qtvap(i,k) = state%q(i,k,       1) + dtime*ptend%q(i,k       ,1)
        qtliq(i,k) = state%q(i,k,ixcldliq) + dtime*ptend%q(i,k,ixcldliq)
        qtice(i,k) = state%q(i,k,ixcldice) + dtime*ptend%q(i,k,ixcldice)
!
! Check for rounding
!
        if (qtvap(i,k) < -100.*qtiny) &
               write(*,*) '(wtrc_cldcond_trend) qtvap < qtiny:',qtvap(i,k)
        qtvap(i,k) = max(qtvap(i,k), 0.)
        if (qtliq(i,k) < -100.*qtiny) &
               write(*,*) '(wtrc_cldcond_trend) qtliq < qtiny:',qtliq(i,k)
        qtliq(i,k) = max(qtliq(i,k), 0.)
        if (qtice(i,k) < -100.*qtiny) &
               write(*,*) '(wtrc_cldcond_trend) qtice < qtiny:',qtice(i,k)
        qtice(i,k) = max(qtice(i,k), 0.)

! precip falling to next layer
        precipt(i) = precipt(i) + (prodprec(i,k) - evapprec(i,k))*state%pdel(i,k)/gravit
        snowabt(i) = snowabt(i) + (prodsnow(i,k) - evapsnow(i,k))*state%pdel(i,k)/gravit
        rainabt(i) = precipt(i) - snowabt(i)

!!        precipt(i) = max(precipt(i), 0.)
!!        snowabt(i) = min(max(snowabt(i), 0.), precipt(i))
!!        rainabt(i) = min(max(rainabt(i), 0.), precipt(i))
!
! Save the rain and snow fluxes
!
        qtrnw(i,k) = dtime*rainabt(i)*gravit/state%pdel(i,k)
        qtsnw(i,k) = dtime*snowabt(i)*gravit/state%pdel(i,k)
!
! Check, zero
!
        if (qtrnw(i,k) < -100.*qtiny) &
               write(*,*) '(wtrc_cldcond_trend) qtrnw < qtiny:',qtrnw(i,k)
        qtrnw(i,k) = max(qtrnw(i,k), 0.)
        if (qtsnw(i,k) < -100.*qtiny) &
               write(*,*) '(wtrc_cldcond_trend) qtsnw < qtiny:',qtsnw(i,k)
        qtsnw(i,k) = max(qtsnw(i,k), 0.)

      end do	! i, columns
   end do	! k, top down integration for precip
!
! budget check for precipitation (DEBUGGING)
!
  do i = 1,ncol
   if (trprec(i,1) > 0. .and. trprec(i,1)+precipt(i)/1000._r8 > 1.e-16) then
    if (abs(precipt(i)/1000._r8-trprec(i,1)) > 1.0e-18) then
     write(6,*) 'Non-consistent precip.: ', i, precipt(i)/1000._r8, trprec(i,1)
    endif
   endif
   if (trsnow(i,1) > 0. .and. trsnow(i,1)+snowabt(i)/1000._r8 > 1.e-16) then
    if (abs(snowabt(i)/1000._r8-trsnow(i,1)) > 1.0e-18) then
     write(6,*) 'Non-consistent snow.: ', i, snowabt(i)/1000._r8, trsnow(i,1)
     write(6,*) 'Non-consistent snow.: ', i, precipt(i)/1000._r8, trprec(i,1)
    endif
   endif
  enddo
!
! Loop over tripplets of water variables (vapour, liquid, ice)
!
   do m = ixwti,ixwtx
     if (wtrc_is_vap(m)) then
!
       ixwtvap = m		!
       ixwtliq = m+1		!=- assume triplets are ordered v, l, i
       ixwtice = m+2		!
!
       ptend%lq(ixwtvap) = .true.
       ptend%lq(ixwtice) = .true.
       ptend%lq(ixwtliq) = .true.
!
! Compute the "macroscale" component of moisture tendancies
! Start to down calulations for forming and falling condensation
!
       trrnw(:ncol,:)  = 0.0
       trsnw(:ncol,:)  = 0.0
       trprecip(:ncol) = 0.0
       trrainab(:ncol) = 0.0
       trsnowab(:ncol) = 0.0
!
       do k = 1, pver
         do i = 1, ncol
!
! For tracers, we only know the pre-adjusted quantities
! (getting the "final" values is the point of this routine... duh!)
!
           if (state%t(i,k) >= tmelt) trsnowab(i) = 0.        
           trrainab(i) = trprecip(i) - trsnowab(i)
!
! "old" values: mix precip from above, with local cloud particles
!
           trvapold(i,k) = state%q(i,k,ixwtvap)
           trliqold(i,k) = state%q(i,k,ixwtliq) + dtime*trrainab(i)*gravit/state%pdel(i,k)
           triceold(i,k) = state%q(i,k,ixwtice) + dtime*trsnowab(i)*gravit/state%pdel(i,k)
!
! Check for rounding
!
           if (trvapold(i,k) < -100.*qtiny) &
                  write(*,*) '(wtrc_cldcond_trend) trvapold < qtiny:',m,trvapold(i,k)
           trvapold(i,k) = max(trvapold(i,k), 0.)
           if (trliqold(i,k) < -100.*qtiny) &
                  write(*,*) '(wtrc_cldcond_trend) trliqold < qtiny:',m,trliqold(i,k)
           trliqold(i,k) = max(trliqold(i,k), 0.)
           if (triceold(i,k) < -100.*qtiny) &
                  write(*,*) '(wtrc_cldcond_trend) triceold < qtiny:',m,triceold(i,k)
           triceold(i,k) = max(triceold(i,k), 0.)
!
! Solve for isotopes using the comprehensive multiphase model
! Here we really should keep account of rain and snow production explictly.
! But for simplity, dump it all into liquid and ice (at top of box)
! with newly "produced" stuff coming out the bottom of the box.
!
            isp(1) = 1
            vapold(1) = qtvapold(i,k)
            liqold(1) = qtliqold(i,k)
            iceold(1) = qticeold(i,k)
            vapent(1) = 0.
!
            isp(2) = iwspec(m)
            vapold(2) = trvapold(i,k)
            liqold(2) = trliqold(i,k)
            iceold(2) = triceold(i,k)
            vapent(2) = 0.

            vapnew(1) = qtvap(i,k)
            liqnew(1) = qtliq(i,k) + (1.-fsnk_cld)*qtrnw(i,k)
            icenew(1) = qtice(i,k) + (1.-fsnk_cld)*qtsnw(i,k)
            vapdet(1) = 0.
            rainpr(1) = fsnk_cld*qtrnw(i,k)
            snowpr(1) = fsnk_cld*qtsnw(i,k)
!
            fnlsk = wtrc_ratio((1.-fsnk_cld)*qtrnw(i,k), liqnew(1))
            fnisk = wtrc_ratio((1.-fsnk_cld)*qtsnw(i,k), icenew(1))
!
            dliqmt = 0.
            dicefz = -dliqmt

            told = state%t(i,k)
            tnew = state%t(i,k) + dtime*ptend%s(i,k)/cpair

!!            call t_startf('wiso_dicm_cldwat')
            call wiso_dicm(piso   , &		! piso = 2
                           isp    , feq_cld, told   , tnew   , &
                           vapold , liqold , iceold , vapent , &
                           vapnew , liqnew , icenew , vapdet , &
                           rainpr , snowpr , dliqmt , dicefz )
!!            call t_stopf('wiso_dicm_cldwat')

            trvap(i,k) = vapnew(2) 
            trliq(i,k) = (1.-fnlsk)*liqnew(2)
            trice(i,k) = (1.-fnisk)*icenew(2)
            trrnw(i,k) = rainpr(2) + fnlsk*liqnew(2)
            trsnw(i,k) = snowpr(2) + fnisk*icenew(2)
!
! Check numerics: "zero" can be order -1e-15, more typically -1.e-22
!
            if (trvap(i,k) < -qtiny) &
               write(*,*) '(wtrc_cldcond_trend) trvap < triny:',trvap(i,k)
            trvap(i,k) = max(trvap(i,k), 0.)
            if (trliq(i,k) < -qtiny) &
               write(*,*) '(wtrc_cldcond_trend) trliq < triny:',trliq(i,k)
            trliq(i,k) = max(trliq(i,k), 0.)
            if (trice(i,k) < -qtiny) &
               write(*,*) '(wtrc_cldcond_trend) trice < triny:',trice(i,k)
            trice(i,k) = max(trice(i,k), 0.)
!
! Convert snow and rain water back to mass units for fll to next layer
!
            trrainab(i) = trrnw(i,k)*state%pdel(i,k)/(gravit*dtime)
            trsnowab(i) = trsnw(i,k)*state%pdel(i,k)/(gravit*dtime)
            trprecip(i) = trrainab(i) + trsnowab(i)

!            trprecip(i) = max(trprecip(i), 0.)
!            trsnowab(i) = min(max(trsnowab(i), 0.), trprecip(i))
!            trrainab(i) = min(max(trrainab(i), 0.), trprecip(i))
!
! Back out final tendancies for use in physics_update
!
            ptend%q(i,k,ixwtvap) = (trvap(i,k) - state%q(i,k,ixwtvap))/dtime
            ptend%q(i,k,ixwtliq) = (trliq(i,k) - state%q(i,k,ixwtliq))/dtime
            ptend%q(i,k,ixwtice) = (trice(i,k) - state%q(i,k,ixwtice))/dtime
!
! Do some debugging checks
!
#ifdef DBGDIAGS
  3     format(a12,4e16.6)
  33    format(a12,4f16.6)
!           if (wiso_delta(iwspec(m), trttl(i,k), qtttl(i,k)) > 2000) then
            write(*,*) 'WTRC_CLDCOND_TEND:',m,i,k
            write(*,3) 'V,L,I old', &
                              qtvapold(i,k),  &
                              qtliqold(i,k),  &
                              qticeold(i,k)
            write(*,3) 'V,L,I new', &
                              qtvap(i,k),  &
                              qtliq(i,k),  &
                              qtice(i,k)
            write(*,33) '  V,L,I  in', &
                              wiso_delta(iwspec(m), state%q(i,k,ixwtvap), state%q(i,k,       1)), &
                              wiso_delta(iwspec(m), state%q(i,k,ixwtliq), state%q(i,k,ixcldliq)), &
                              wiso_delta(iwspec(m), state%q(i,k,ixwtice), state%q(i,k,ixcldice)) 
            write(*,33) 'V,L,I old', &
                              wiso_delta(iwspec(m), trvapold(i,k), qtvapold(i,k)),  &
                              wiso_delta(iwspec(m), trliqold(i,k), qtliqold(i,k)),  &
                              wiso_delta(iwspec(m), triceold(i,k), qticeold(i,k))
            write(*,33) '  V,L,I new', &
                              wiso_delta(iwspec(m), trvap   (i,k), qtvap   (i,k)),  &
                              wiso_delta(iwspec(m), trliq   (i,k), qtliq   (i,k)),  &
                              wiso_delta(iwspec(m), trice   (i,k), qtice   (i,k))
!!          end if
#endif
!
#undef DEBUGIT
#ifdef DEBUGIT
!!
!! This block works well if the "tracer" is set equal to the total at the top
           if (m == 4) then
             if (abs(trliq(i,k) - qtliq(i,k)) > 1.e-14 .or. &
                 abs(trice(i,k) - qtice(i,k)) > 1.e-14 .or. &
                 abs(trvap(i,k) - qtvap(i,k)) > 1.e-14 .or. &
                 abs(ptend%q(i,k,       1) - ptend%q(i,k,ixwtvap)) > 1.e-12 .or. &
                 abs(ptend%q(i,k,ixcldliq) - ptend%q(i,k,ixwtliq)) > 1.e-12 .or. &
                 abs(ptend%q(i,k,ixcldice) - ptend%q(i,k,ixwtice)) > 1.e-12 ) then
               write(*,*) 'wtrc_cldcond_tend BUDGET ERROR --------------  m:',m,i,k
               write(*,3) 'QTOT :',qtttl(i,k),trttl(i,k),qtttl(i,k)-trttl(i,k)
               write(*,3) 'VOLD :',qtvapold(i,k),trvapold(i,k),qtvapold(i,k)-trvapold(i,k)
               write(*,3) 'LOLD :',qtliqold(i,k),trliqold(i,k),qtliqold(i,k)-trliqold(i,k)
               write(*,3) 'IOLD :',qticeold(i,k),triceold(i,k),qticeold(i,k)-triceold(i,k)

!               write(*,3) 'RAIN :',rainabt (i,k),trrainab(i,k),rainabt(i,k)-trvrainab(i,k)
!               write(*,3) 'SNOW :',snowabt (i,k),trsnowab(i,k),snowabt(i,k)-trvsnowab(i,k)
!               write(*,3) 'PREC :',precabt (i,k),trprecab(i,k),precabt(i,k)-trvprecab(i,k)

               write(*,3) 'VNEW :',qtvap   (i,k),trvap   (i,k),qtvap   (i,k)-trvap   (i,k)
               write(*,3) 'LNEW :',qtliq   (i,k),trliq   (i,k),qtliq   (i,k)-trliq   (i,k)
               write(*,3) 'INEW :',qtice   (i,k),trice   (i,k),qtice   (i,k)-trice   (i,k)
               write(*,3) 'CLIQ :',qtrnw   (i,k),trrnw   (i,k),qtrnw   (i,k)-trrnw   (i,k)
               write(*,3) 'CICE :',qtsnw   (i,k),trsnw   (i,k),qtsnw   (i,k)-trsnw   (i,k)

               write(*,3) 'dVAP :',ptend%q(i,k,       1),ptend%q(i,k,ixwtvap),ptend%q(i,k,       1)-ptend%q(i,k,ixwtvap)
               write(*,3) 'dLIQ :',ptend%q(i,k,ixcldliq),ptend%q(i,k,ixwtliq),ptend%q(i,k,ixcldliq)-ptend%q(i,k,ixwtliq)
               write(*,3) 'dICE :',ptend%q(i,k,ixcldice),ptend%q(i,k,ixwtice),ptend%q(i,k,ixcldice)-ptend%q(i,k,ixwtice)

               call endrun('wtrc_cldcond_tend BUDGET DEBUG')
 3             format(a8,3e14.6)
             end if
           end if		!  == 4
#endif

         end do		! i, columns
       end do		! k, top down levels
!
! Fix for numerics - no tracer precip if no total precip
!     (THIS COUPLES THE TRACERS TO THE PROGNOSTIC)
!!         do i = 1, ncol
!!           if (precip(i) == 0.) trprecip(i) = 0.
!!           if (snowab(i) == 0.) trsnowab(i) = 0.
!!           if (rainab(i) == 0.) trrainab(i) = 0.
!!         end do
!
! Convert surface flux to precipitation amounts in preparation for output
!
       trprec(:ncol,m) =  trprecip(:ncol) / 1000.
       trsnow(:ncol,m) =  trsnowab(:ncol) / 1000.
!
! Save diagnostics to history file
!
!!       call outfld ()

     end if	! m is vapour (first in triplet)
   end do	! m, water tracers
!
   return
end subroutine wtrc_cldcond_tend


!=======================================================================
subroutine wtrc_cldcond_zmconv_detrain(trdlf, trdif, cld, state, ptend)

!-----------------------------------------------------------------------
! 
! Purpose: dump water detrained from convection into cloud water
!
!   See cldcond_zmconv_detrain in cldcond.F90, by B. Boville. 
!   This routine should be called immediately after the main routine,
!   such that only the tracer physics_update flags are set.
!
! Author: David Noone <dcn@caltech.edu> - Thu Aug  7 17:50:27 PDT 2003
!
!-----------------------------------------------------------------------
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use ppgrid,        only: pcols, pver
   use physics_types, only: physics_state, physics_ptend
   use constituents,  only: pcnst, pnats, cnst_name
   use history,       only: outfld
   use water_tracers, only: wtrc_is_vap, ixwti, ixwtx
!-----------------------------------------------------------------------
   implicit none
!---------------------------- Arguments --------------------------------

    type(physics_state), intent(in  )  :: state   ! state variables
    type(physics_ptend), intent(inout) :: ptend   ! package tendencies

    real(r8), intent(in) :: trdlf(pcols,pver,pcnst+pnats) ! detrained tracers liquid
    real(r8), intent(in) :: trdif(pcols,pver,pcnst+pnats) ! detrained tracers ice
    real(r8), intent(in) :: cld(pcols,pver)       ! cloud fraction


!------------------------- Local Variables -----------------------------
    integer :: ixwtvap,ixwtliq,ixwtice            ! constituent index
    integer :: i,k,m                              ! loop indexes
!-----------------------------------------------------------------------
!
   do m = ixwti, ixwtx				! only in range of tracers
     if (wtrc_is_vap(m)) then
       ixwtvap = m
       ixwtliq = m+1
       ixwtice = m+2

       ptend%lq(ixwtliq) = .TRUE.			! flag update
       ptend%lq(ixwtice) = .TRUE.			! flag update

       do k = 1, pver
         do i = 1, state%ncol
           ptend%q(i,k,ixwtliq) = trdlf(i,k,m)		! assign to cloud water
           ptend%q(i,k,ixwtice) = trdif(i,k,m)		! assign to cloud water
         end do
       end do
       call outfld('DLF'//trim(cnst_name(m)),trdlf(:,:,m), pcols,state%lchnk)
       call outfld('DIF'//trim(cnst_name(m)),trdif(:,:,m), pcols,state%lchnk)
     end if
   end do
!
   return
end subroutine wtrc_cldcond_zmconv_detrain

!=======================================================================
end module wtrc_cldcond
