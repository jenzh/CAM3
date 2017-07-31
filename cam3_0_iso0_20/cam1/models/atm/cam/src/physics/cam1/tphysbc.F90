#include <misc.h>
#include <params.h>

subroutine tphysbc (ztodt,   pblht,   tpert,   ts,      sst,     &
                    qpert,   precl,   precc,   precsl,  precsc,  &
                    asdir,   asdif,   aldir,   aldif,   snowh,   &
                    qrs,     qrl,     flwds,   fsns,    fsnt,    &
                    flns,    flnt,    lwup,    srfrad,  sols,    &
                    soll,    solsd,   solld,   state,   tend,    &
                    pbuf,    prcsnw,  fsds ,   landm,   landfrac,&
		    ocnfrac, icefrac,&
                    wtprecl, wtprecc, wtprecsl,wtprecsc,wtprcsnw  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Tendency physics BEFORE coupling to land, sea, and ice models.
! 
! Method: 
! Call physics subroutines and compute the following:
!     o cloud calculations (cloud fraction, emissivity, etc.)
!     o radiation calculations
! Pass surface fields for separate surface flux calculations
! Dump appropriate fields to history file.
! 
! Author: CCM1, CMS Contact: J. Truesdale
!
! Modified for water isotopes:
!    David Noone <dcn@colorado.edu> - Fri Jul  2 15:59:00 MDT 2004
!
! 
!-----------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p, get_lon_all_p, get_lat_all_p
   use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use cldcond,         only: cldcond_tend, cldcond_zmconv_detrain, cldcond_sediment
   use param_cldoptics, only: param_cldoptics_calc
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
   use diagnostics,     only: diag_dynvar
   use history,         only: outfld
   use physconst,       only: gravit, latvap, cpair, tmelt, cappa, zvir, rair, rga
   use radheat,         only: radheat_net
   use constituents,    only: pcnst, pnats, ppcnst, qmin, cnst_name
   use constituents,    only: dcconnam, cnst_get_ind
   use zm_conv,         only: zm_conv_evap, zm_convr
   use time_manager,    only: is_first_step, get_nstep, get_curr_calday
   use moistconvection, only: cmfmca
   use check_energy,    only: check_energy_chng, check_energy_fix
   use check_energy,    only: check_tracers_data, check_tracers_init, check_tracers_chng
   use dycore,          only: dycore_is
   use cloudsimulator,  only: doisccp, cloudsimulator_run
   use aerosol_intr,    only: aerosol_wet_intr

   use water_tracers,   only: trace_water, ixwti, ixwtx, &
                              wtrc_is_wtrc, wtrc_is_vap, wtrc_is_liq, wtrc_is_ice, &
                              wtrc_ratio_all, wtrc_chkdelta, wtrc_check, wtrc_qchk1, wtrc_qchk2, ixh2oq
   use water_isotopes,  only: wisotope
   use wtrc_cldcond,    only: wtrc_cldcond_tend, wtrc_cldcond_zmconv_detrain, wtrc_cldcond_sediment
   use wtrc_zm_conv,    only: wtrc_zm_conv_evap, wtrc_zm_convr
   use wiso_histogram,  only: wiso_hist_diag

   implicit none

#include <comctl.h>
!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: ts(pcols)                      ! surface temperature
   real(r8), intent(in) :: sst(pcols)                     ! sea surface temperature
   real(r8), intent(inout) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(inout) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(inout) :: qpert(pcols,ppcnst)         ! Thermal humidity & constituent excess
   real(r8), intent(in) :: asdir(pcols)                  ! Albedo: shortwave, direct
   real(r8), intent(in) :: asdif(pcols)                  ! Albedo: shortwave, diffuse
   real(r8), intent(in) :: aldir(pcols)                  ! Albedo: longwave, direct
   real(r8), intent(in) :: aldif(pcols)                  ! Albedo: longwave, diffuse
   real(r8), intent(in) :: snowh(pcols)                  ! Snow depth (liquid water equivalent)
   real(r8), intent(inout) :: qrs(pcols,pver)            ! Shortwave heating rate
   real(r8), intent(inout) :: qrl(pcols,pver)            ! Longwave  heating rate
   real(r8), intent(inout) :: flwds(pcols)               ! Surface longwave down flux
   real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(in) :: lwup(pcols)                    ! Surface longwave up flux
   real(r8), intent(out) :: srfrad(pcols)                 ! Net surface radiative flux (watts/m**2)
   real(r8), intent(inout) :: sols(pcols)                   ! Direct beam solar rad. onto srf (sw)
   real(r8), intent(inout) :: soll(pcols)                   ! Direct beam solar rad. onto srf (lw)
   real(r8), intent(inout) :: solsd(pcols)                  ! Diffuse solar radiation onto srf (sw)
   real(r8), intent(inout) :: solld(pcols)                  ! Diffuse solar radiation onto srf (lw)
   real(r8), intent(out) :: precl(pcols)                  ! Large-scale precipitation rate
   real(r8), intent(out) :: precc(pcols)                  ! Convective-scale preciptn rate
   real(r8), intent(out) :: precsl(pcols)                 ! L.S. snowfall rate
   real(r8), intent(out) :: precsc(pcols)                 ! C.S. snowfall rate
   real(r8), intent(out) :: prcsnw(pcols)                 ! snowfall rate (precsl + precsc)
   real(r8), intent(out) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
   real(r8), intent(in) :: landfrac(pcols)                ! land fraction
   real(r8), intent(in) :: ocnfrac(pcols)                ! land fraction
   real(r8), intent(in) :: icefrac(pcols)                ! land fraction

   real(r8), intent(out) :: wtprecl(pcols,pcnst+pnats)   ! tracer Large-scale precipitation rate
   real(r8), intent(out) :: wtprecc(pcols,pcnst+pnats)   ! tracer Convective-scale preciptn rate
   real(r8), intent(out) :: wtprecsl(pcols,pcnst+pnats)  ! tracer L.S.  snowfall rate
   real(r8), intent(out) :: wtprecsc(pcols,pcnst+pnats)  ! tracer C.S.  snowfall rate
   real(r8), intent(out) :: wtprcsnw(pcols,pcnst+pnats)  ! snowfall rate (precsl + precsc)

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf
!
!---------------------------Local workspace-----------------------------
!
   real(r8) :: rhdfda(pcols,pver)            ! dRh/dcloud, old 
   real(r8) :: rhu00 (pcols,pver)            ! Rh threshold for cloud, old

   type(physics_ptend)   :: ptend                  ! indivdual parameterization tendencies

   integer :: nstep                          ! current timestep number
   integer      lat(pcols)                   ! current latitudes(indices)
   integer      lon(pcols)                   ! current longtitudes(indices)

   real(r8) :: calday                        ! current calendar day
   real(r8) :: clat(pcols)                   ! current latitudes(radians)
   real(r8) :: clon(pcols)                   ! current longitudes(radians)

   real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
   real(r8) :: ftem(pcols,pver)              ! Temporary workspace for outfld variables
   real(r8) :: zmrprd(pcols,pver)            ! rain production in ZM convection
   real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c
   real(r8) :: cmfsl(pcols,pver)             ! Moist convection lw stat energy flux
   real(r8) :: cmflq(pcols,pver)             ! Moist convection total water flux
   real(r8) :: dtcond(pcols,pver)            ! dT/dt due to moist processes
   real(r8) :: dqcond(pcols,pver,ppcnst)     ! dq/dt due to moist processes

   real(r8) cldst(pcols,pver)
   real(r8) cltot(pcols)                      ! Diagnostic total cloud cover
   real(r8) cllow(pcols)                      !       "     low  cloud cover
   real(r8) clmed(pcols)                      !       "     mid  cloud cover
   real(r8) clhgh(pcols)                      !       "     hgh  cloud cover
   real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation
   real(r8) cmfdqr2(pcols,pver)               ! dq/dt due to moist convective rainout
   real(r8) cmfmc2(pcols,pverp)               ! Moist convection cloud mass flux
   real(r8) cmfsl2(pcols,pver)                ! Moist convection lw stat energy flux
   real(r8) cmflq2(pcols,pver)                ! Moist convection total water flux
   real(r8) cnt(pcols)                        ! Top level of convective activity
   real(r8) cnb(pcols)                        ! Lowest level of convective activity
   real(r8) cnt2(pcols)                       ! Top level of convective activity
   real(r8) cnb2(pcols)                       ! Bottom level of convective activity
   real(r8) concld(pcols,pver)             
   real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
   real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from convection
   real(r8) dif(pcols,pver)                   ! Detraining cld H20 from convection
   real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev
   real(r8) prect(pcols)                      ! total (conv+large scale) precip rate
   real(r8) dlf2(pcols,pver)                   ! dq/dt due to rainout terms
   real(r8) dif2(pcols,pver)                   ! dq/dt due to rainout terms
   real(r8) qpert2(pcols,ppcnst)              ! Perturbation q
   real(r8) rtdt                              ! 1./ztodt
   real(r8) tpert2(pcols)                     ! Perturbation T
   real(r8) pmxrgn(pcols,pverp)               ! Maximum values of pressure for each
!                                             !    maximally overlapped region.
!                                             !    0->pmxrgn(i,1) is range of pressure for
!                                             !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!                                             !    2nd region, etc
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns

   integer nmxrgn(pcols)                      ! Number of maximally overlapped regions
   integer  i,k,m                             ! Longitude, level, constituent indices
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
                                           
!  real(r8) engt                              ! Thermal   energy integral
!  real(r8) engk                              ! Kinetic   energy integral
!  real(r8) engp                              ! Potential energy integral
   real(r8) rel(pcols,pver)                   ! Liquid cloud particle effective radius
   real(r8) rei(pcols,pver)                   ! Ice effective drop size (microns)
   real(r8) emis(pcols,pver)                  ! Cloud longwave emissivity
   real(r8) clc(pcols)                        ! Total convective cloud (cloud scheme)
   real(r8) :: cicewp(pcols,pver)             ! in-cloud cloud ice water path
   real(r8) :: cliqwp(pcols,pver)             ! in-cloud cloud liquid water path
!
   real(r8) dellow(pcols)                     ! delta p for bottom three levels of model
   real(r8) tavg(pcols)                       ! mass weighted average temperature for 

! physics buffer fields to compute tendencies for cloud condensation package
   integer itim, ifld
   real(r8), pointer, dimension(:,:) :: qcwat, tcwat, lcwat, cld

! physics buffer fields for total energy and mass adjustment
   real(r8), pointer, dimension(:  ) :: teout
   real(r8), pointer, dimension(:,:) :: qini
   real(r8), pointer, dimension(:,:) :: tini
!                                          
! Used for OUTFLD only                     
!                                          
   real(r8) icwmr1(pcols,pver)                ! in cloud water mixing ration for zhang scheme
   real(r8) icwmr2(pcols,pver)                ! in cloud water mixing ration for hack scheme
   real(r8) fracis(pcols,pver,ppcnst)         ! fraction of transported species that are insoluble
   real(r8) timestep(pcols)
!
!     Variables for doing deep convective transport outside of zm_convr
!
   real(r8) mu2(pcols,pver)
   real(r8) eu2(pcols,pver)
   real(r8) du2(pcols,pver)
   real(r8) md2(pcols,pver)
   real(r8) ed2(pcols,pver)
   real(r8) dp(pcols,pver)
   real(r8) dpdry(pcols,pver)
   real(r8) dsubcld(pcols)
   real(r8) conicw(pcols,pver)
   real(r8) cmfdqrt(pcols,pver)               ! dq/dt due to moist convective rainout

   real(r8) iuwmr1(pcols,pver)             ! ZMC: cloud updraft mixing ratio from zm scheme
   real(r8) idwmr1(pcols,pver)             ! ZMC: cloud downdraft mixing ratio from zm scheme
   real(r8) cu(pcols,pver)                 ! ZMC: updraft condensation 
   real(r8) rr(pcols,pver)                 ! ZMC: updraft rain rate
   real(r8) evp(pcols,pver)                ! ZMC: downdraft evaporation 
   real(r8) tu(pcols,pver)		   ! ZMC: updraft temperature
   real(r8) td(pcols,pver)		   ! ZMC: downdraft temperature

! stratiform precipitation variables
   real(r8) :: prec_pcw(pcols)                ! total precip from prognostic cloud scheme
   real(r8) :: snow_pcw(pcols)                ! snow from prognostic cloud scheme
   real(r8) :: prec_sed(pcols)                ! total precip from cloud sedimentation
   real(r8) :: snow_sed(pcols)                ! snow from cloud ice sedimentation

! convective precipitation variables
   real(r8) :: prec_zmc(pcols)                ! total precipitation from ZM convection
   real(r8) :: snow_zmc(pcols)                ! snow from ZM convection
   real(r8) :: prec_cmf(pcols)                ! total precipitation from Hack convection
   real(r8) :: snow_cmf(pcols)                ! snow from Hack convection

! energy checking variables
   real(r8) :: zero(pcols)                    ! array of zeros
   real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
   real(r8) :: rice(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
   real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
   real(r8) :: rice2(pcols)                   ! vertical integral of liquid from shallow scheme
   real(r8) :: flx_cnd(pcols)
   real(r8) :: flx_ice(pcols)
   real(r8) :: flx_heat(pcols)
   logical  :: conserve_energy = .true.       ! flag to carry (QRS,QRL)*dp across time steps
   type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes

   integer jt(pcols)
   integer jd(pcols)
   integer maxg(pcols)
   integer ideep(pcols)
   integer lengath
   real(r8) cldc(pcols,pver)
   real(r8) nevapr(pcols,pver)		! evap precip (rain+snow)
   real(r8) nevaps(pcols,pver)		! evap snow
   real(r8) qme(pcols,pver)
   real(r8) prain(pcols,pver)		! production of precip (rain+snow)
   real(r8) psnow(pcols,pver)		! production of snow
   real(r8) cflx(pcols,ppcnst)
!
! Additional variables for water tracers
!
   real(r8) qbot_tmp(pcols,ppcnst)	   ! qbot for srfxfer
   real(r8) Rtrc(pcols,pver,ppcnst)        ! tracer ratio to prognostic water
   real(r8) wtprect(pcols,ppcnst)	   ! tracer total precipitation
!
   real(r8) wtcme(pcols,pver,ppcnst)       ! ZMC: cmf condensation - evaporation
   real(r8) wtdlf(pcols,pver,ppcnst)       ! ZMC: Detraining tracer cld H20 from convection
   real(r8) wtdif(pcols,pver,ppcnst)       ! ZMC: Detraining tracer cld H20 from convection
   real(r8) wtql(pcols,pver,ppcnst)        ! ZMC: updraft cloud liquid water ZM scheme
   real(r8) wtqi(pcols,pver,ppcnst)        ! ZMC: updraft cloud ice water ZM scheme
   real(r8) wtqu(pcols,pver,ppcnst)        ! ZMC: updraft water vapour ZM scheme
   real(r8) wtqd(pcols,pver,ppcnst)        ! ZMC: downdraft water vapour ZM scheme
   real(r8) wtzmrprd(pcols,pver,ppcnst)    ! ZMC: rain production in ZM convection
   real(r8) wtpflx(pcols,pverp,ppcnst)     ! ZMC: Conv rain flux thru out btm of lev
   real(r8) ficed(pcols,pver)	      	   ! ZMC: ice fraction detrained
   real(r8) ficeup(pcols,pver)	      	   ! ZMC: ice fraction in updraft

   real(r8) wtdlf2(pcols,pver,ppcnst)      ! CMF: Detraining tracer cld H20 from convection
   real(r8) wtdif2(pcols,pver,ppcnst)      ! CMF: Detraining tracer cld H20 from convection
   real(r8) wtdqr2(pcols,pver,ppcnst)      ! CMF: Tendency due to rain out
   real(r8) flxprec(pcols,pverp)           ! EVP: Convective-scale flux of precip at interfaces (kg/m2/s)
   real(r8) flxsnow(pcols,pverp)           ! EVP: Convective-scale flux of snow   at interfaces (kg/m2/s)
   real(r8) fevprec(pcols,pver)            ! EVP: fraction of precip flux that evaps
   real(r8) fevsnow(pcols,pver)            ! EVP: fraction flux that is snow
   real(r8) dqlfz(pcols,pver)              ! PCW: change in liquid due to freeze (repartitioning)
   real(r8) fxliq(pcols,pverp)             ! SED: fluxes at the interfaces, liquid (positive = down)
   real(r8) fxice(pcols,pverp)             ! SED: fluxes at the interfaces, ice    (positive = down)
!
   real(r8) wtprec_zmc(pcols,ppcnst)       ! tracer total precipitation from ZM convection
   real(r8) wtsnow_zmc(pcols,ppcnst)       ! tracer snow from ZM convection
   real(r8) wtprec_cmf(pcols,ppcnst)       ! tracer total precip from Hack convection
   real(r8) wtsnow_cmf(pcols,ppcnst)       ! tracer snow from Hack convection
   real(r8) wtprec_sed(pcols,ppcnst)       ! tracer total precip from cloud sedimentation
   real(r8) wtsnow_sed(pcols,ppcnst)       ! tracer snow from cloud ice sedimentation
   real(r8) wtprec_pcw(pcols,ppcnst)       ! tracer total precip from prog.cloud scheme
   real(r8) wtsnow_pcw(pcols,ppcnst)       ! tracer snow from prognostic cloud scheme
!
! Local flags to turn on and off specific tracer routines (does not effect climate)
!
   logical :: lwtrczmc = .true.	  	   ! ZM convection
   logical :: lwtrccmf = .true.		   ! Hack convection
   logical :: lwtrcevp = .true.		   ! evaporative falling convective precip
   logical :: lwtrcdet = .true.		   ! accumulated detrained water
   logical :: lwtrcpcw = .true.		   ! prognostic cloud water
   logical :: lwtrcsed = .true.		   ! cloud liquid/ice edimentation
   logical :: lwisoeql = .true.		   ! cloud liquid/vapour equilibration
!
!-----------------------------------------------------------------------
   zero = 0.
!
   lchnk = state%lchnk
   ncol  = state%ncol

   rtdt = 1./ztodt

   nstep = get_nstep()
   calday = get_curr_calday()
!
! Output NSTEP for debugging
!
   timestep(:ncol) = nstep
   call outfld ('NSTEP   ',timestep, pcols, lchnk)

!  dry surface pressure
   call outfld ('PSDRY',  state%psdry, pcols, lchnk)

! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('QCWAT')
   qcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('TCWAT')
   tcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('LCWAT')
   lcwat => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld = pbuf_get_fld_idx('TEOUT')
   teout => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,itim)
   ifld = pbuf_get_fld_idx('QINI')
   qini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld = pbuf_get_fld_idx('TINI')
   tini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
!
! Initialize tracer variables 
! (needed as individual routined may be turned off)
!
   wtprec_zmc(:,:) = 0.
   wtsnow_zmc(:,:) = 0.
   wtprec_cmf(:,:) = 0.
   wtsnow_cmf(:,:) = 0.
   wtprec_sed(:,:) = 0.
   wtsnow_sed(:,:) = 0.
   wtprec_pcw(:,:) = 0.
   wtsnow_pcw(:,:) = 0.
   wtdlf(:,:,:) = 0.
   wtdif(:,:,:) = 0.

   dif(:,:) = 0.		! debugging
   dif2(:,:) = 0.		! debugging


   call wtrc_qchk2('TPH_W_top ','state%q',ncol,state%q(:,:,ixh2oq),state%q(:,:,1))
   call wtrc_chkdelta('tphystop', ncol, state%q)
   call wtrc_check('tphysbc_start', ncol, state%q)

#ifdef FIXWATER		/* for debugging */
   write(*,*) 'TPHYSBC - water fixed to total on input'
   state%q(:,:,4) = state%q(:,:,1) 		! vap
   state%q(:,:,5) = state%q(:,:,2) 		! liq
   state%q(:,:,6) = state%q(:,:,3) 		! ice
#endif
   call wtrc_qchk1('TPHYS_IN','qpert',ncol,qpert(:,ixh2oq),qpert(:,1),1.e-12)

!
! Set physics tendencies to 0
   tend %dTdt(:ncol,:pver)  = 0.
   tend %dudt(:ncol,:pver)  = 0.
   tend %dvdt(:ncol,:pver)  = 0.

   call physics_ptend_init (ptend) ! Initialize parameterization tendency structure

!
! Make sure that input tracers are all positive (probably unnecessary)
!
   call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
              ppcnst,qmin  ,state%q )
!
! Setup q and t accumulation fields
!
   dqcond(:ncol,:,:) = state%q(:ncol,:,:)
   dtcond(:ncol,:)   = state%s(:ncol,:)

   fracis (:ncol,:,1:ppcnst) = 1.

! compute mass integrals of input tracers state
   call check_tracers_init(state, tracerint)

!===================================================
! Global mean total energy fixer
!===================================================
   !*** BAB's FV heating kludge *** save the initial temperature
   tini(:ncol,:pver) = state%t(:ncol,:pver)
   if (dycore_is('LR')) then
      call check_energy_fix(state, ptend, nstep, flx_heat)
      call physics_update(state, tend, ptend, ztodt)
      call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
   end if
   qini(:ncol,:pver) = state%q(:ncol,:pver,1)
   call wtrc_check('energy fix', ncol, state%q)

   call outfld('TEOUT', teout       , pcols, lchnk   )
   call outfld('TEINP', state%te_ini, pcols, lchnk   )
   call outfld('TEFIX', state%te_cur, pcols, lchnk   )
!
!===================================================
! Dry adjustment
!===================================================

! Copy state info for input to dadadj
! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy
! (Modified to mix all tracers - dcn)

   ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
   ptend%q(:ncol,:pver,:) = state%q(:ncol,:pver,:)

   call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
                ptend%s, ptend%q)
   ptend%name  = 'dadadj'
   ptend%ls    = .TRUE.
   ptend%lq(:) = .TRUE.
   ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
   ptend%q(:ncol,:,:) = (ptend%q(:ncol,:,:) - state%q(:ncol,:,:))/ztodt
   call physics_update (state, tend, ptend, ztodt)
   call wtrc_check('dadadj', ncol, state%q)
!
!===================================================
! Moist convection
!===================================================
!
! Since the PBL doesn't pass constituent perturbations, they
! are zeroed here for input to the moist convection routine
!
   if (.not. trace_water) then		! water tracers have it!
     qpert(:ncol,2:ppcnst) = 0.0
   end if
   call wtrc_qchk1('TPHYS_ZMC','qpert',ncol,qpert(:,ixh2oq),qpert(:,1),1.e-12)
!
! Begin with Zhang-McFarlane (1996) convection parameterization
!
   call t_startf ('zm_convr')
   call zm_convr( lchnk,    ncol, &
                  state%t, state%q,  prec_zmc,   cnt,     cnb,      &
                  pblht,   state%zm, state%phis,    state%zi,   ptend%q(:,:,1),     &
                  ptend%s, state%pmid,   state%pint,  state%pdel,       &
                  .5*ztodt,cmfmc,    cmfcme,             &
                  tpert,   dlf,      dif,     pflx,    zdu,     zmrprd,   &
                  mu2,     md2,      du2,     eu2,     ed2,      &
                  dp,      dsubcld,  jt,      maxg,    ideep,    &
                  lengath, icwmr1,   iuwmr1,  idwmr1,  cu,       &
                  evp,     rr,       tu,      td,      jd,       &
                  rliq,    rice,     ficed,   ficeup)

   ptend%name  = 'zm_convr'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
!
   if (trace_water .and. lwtrczmc) then
     call t_startf ('wtrc_zm_convr')
!
! initialize tracer variables
       wtdlf(:ncol,:,:) = 0.0
       wtdif(:ncol,:,:) = 0.0
       wtcme(:ncol,:,:) = 0.0
       wtql(:ncol,:,:) = 0.
       wtqi(:ncol,:,:) = 0.
       wtqu(:ncol,:,:) = 0.
       wtqd(:ncol,:,:) = 0.
       wtzmrprd(:ncol,:,:) = 0.
       wtpflx(:ncol,:,:) = 0.
       wtprec_zmc(:ncol,:) = 0.
!
! assign totals 
       wtcme(:ncol,:,1) = cmfcme(:ncol,:)
       wtdlf(:ncol,:,1) = dlf(:ncol,:)
       wtdif(:ncol,:,1) = dif(:ncol,:)
       wtql(:ncol,:,1) = icwmr1(:ncol,:)	! partitioned to liquid on output
       wtqi(:ncol,:,1) = 0.0			! partitioned to ice on output
       wtqu(:ncol,:,1) = iuwmr1(:ncol,:)
       wtqd(:ncol,:,1) = idwmr1(:ncol,:)
       wtzmrprd(:ncol,:,1) = zmrprd(:ncol,:)
       wtpflx(:ncol,:,1) = pflx(:ncol,:)
       wtprec_zmc(:ncol,1) = prec_zmc(:ncol)

   call wtrc_check('pre wtrc_zm_convr', ncol, state%q)
   call wtrc_qchk2('TPH_W_ZMin ','state%q',ncol,state%q(:,:,ixh2oq),state%q(:,:,1))

       call wtrc_zm_convr(lchnk   ,ncol    , &
                      state%t       ,state%q      ,wtprec_zmc    , &
                      ptend%q    , tu     ,td     ,&
                      state%pdel     , &
                      .5*ztodt,wtcme     ,          &
                      wtdlf   ,wtdif   ,wtpflx    ,wtzmrprd    , &
                      mu2     ,md2     ,du2     ,eu2     ,ed2     , &
                      dp      ,dsubcld ,jt      ,maxg    ,ideep   , &
                      lengath ,wtql    ,wtqi    ,wtqu    ,wtqd    ,cu      , &
                      evp     ,state%zi,state%phis ,rr   , jd     , &
                      ficed   ,ficeup )

   call wtrc_check('post wtrc_zm_convr', ncol, state%q)
   call wtrc_qchk2('TPH_W_ZMout','state%q',ncol,state%q(:,:,ixh2oq),state%q(:,:,1))
   call wtrc_qchk2('TPH_W_ZMout','ptend%q',ncol,ptend%q(:,:,ixh2oq),ptend%q(:,:,1))
   call wtrc_qchk2('TPH_W_ZMout','cme    ',ncol,wtcme(:,:,ixh2oq),wtcme(:,:,1))
   call wtrc_qchk2('TPH_W_ZMout','rprd   ',ncol,wtzmrprd(:,:,ixh2oq),wtzmrprd(:,:,1))
   call wtrc_qchk2('TPH_W_ZMout','dlf    ',ncol,wtdlf(:,:,ixh2oq),wtdlf(:,:,1))
   call wtrc_qchk2('TPH_W_ZMout','dif    ',ncol,wtdif(:,:,ixh2oq),wtdif(:,:,1))

       do m = ixwti,ixwtx
         if (wtrc_is_vap(m)) ptend%lq(m) = .true.
       end do
     call t_stopf ('wtrc_zm_convr')
   end if
!
   cmfsl (:ncol,:) = 0. ! This is not returned from zm, hence it is zeroed.
   cmflq (:ncol,:) = 0. ! This is not returned from zm, hence it is zeroed.

   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('ZMDT    ',ftem           ,pcols   ,lchnk   )
   call outfld('ZMDQ    ',ptend%q(1,1,1) ,pcols   ,lchnk   )
   call t_stopf('zm_convr')

   call physics_update(state, tend, ptend, ztodt)
   call wtrc_check('zm_convr', ncol, state%q)
   call wtrc_chkdelta('post zm_convr', ncol, state%q)
!
! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
   call zm_conv_evap(state, ptend, zmrprd, cld, ztodt, &
                     flxprec, flxsnow, fevprec, fevsnow, &
                     prec_zmc, snow_zmc, .false.)

! Water tracers: evaporat falling precip, and assign state
   if (trace_water .and. lwtrcevp) then
     call wtrc_zm_conv_evap(state  , ptend  , wtzmrprd, cld , ztodt , &
                            flxprec, flxsnow, fevprec, fevsnow, &
                            wtprec_zmc, wtsnow_zmc, .false. )
   end if
!
   call physics_update(state, tend, ptend, ztodt)
   call wtrc_check('zm_evap', ncol, state%q)
! Check energy integrals, including "reserved liquid"
   flx_cnd(:ncol) = prec_zmc(:ncol) + rliq(:ncol)
   call check_energy_chng(state, tend, "zm_evap", nstep, ztodt, zero, flx_cnd, snow_zmc, zero)

! Transport cloud water and ice only
!
   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   ptend%name = 'convtran1'
   ptend%lq(ixcldice) = .true.
   ptend%lq(ixcldliq) = .true.
   if (trace_water) then
     do m = ixwti, ixwtx
       if (wtrc_is_liq(m) .or. wtrc_is_ice(m)) &
          ptend%lq(m) = .true.
     end do
   end if
   call convtran (lchnk,                                        &
                  ptend%lq,state%q, ppcnst,  mu2,     md2,   &
                  du2,     eu2,     ed2,     dp,      dsubcld,  &
                  jt,      maxg,    ideep,   1,       lengath,  &
                  nstep,   fracis,  ptend%q, dpdry  )
   call physics_update (state, tend, ptend, ztodt)
   call wtrc_check('convtran1', ncol, state%q)
!
! Convert mass flux from reported mb/s to kg/m^2/s
!
   cmfmc(:ncol,:pver) = cmfmc(:ncol,:pver) * 100./gravit
!
! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
!
   call t_startf('cmfmca')
   tpert2(:ncol  ) =0.
   qpert2(:ncol,:) = qpert(:ncol,:)  ! BAB Why is this not zero, if tpert2=0???
   call wtrc_qchk1('TPHYS_CMF','qpert2',ncol,qpert2(:,ixh2oq),qpert2(:,1),1.e-12)

   call cmfmca (lchnk,   ncol, &
                nstep,   ztodt,   state%pmid,  state%pdel,   &
                state%rpdel,   state%zm,      tpert2,  qpert2,  state%phis,     &
                pblht,   state%t,   state%q,   ptend%s,   ptend%q,      &
                cmfmc2,  cmfdqr2, cmfsl2,  cmflq2,  prec_cmf,   &
                dlf2,     dif2,   cnt2,    cnb2,    icwmr2   , rliq2,  rice2, & 
                state%pmiddry, state%pdeldry, state%rpdeldry, &
                wtprec_cmf, wtdlf2, wtdif2,   wtdqr2, lwtrccmf)
   ptend%name  = 'cmfmca'
   ptend%ls    = .TRUE.
   ptend%lq(:) = .TRUE.

! Add shallow cloud water detrainment to cloud water detrained from ZM
   dlf(:ncol,:pver) = dlf(:ncol,:pver) + dlf2(:ncol,:pver)
   dif(:ncol,:pver) = dif(:ncol,:pver) + dif2(:ncol,:pver)
   wtdlf(:ncol,:pver,:) = wtdlf(:ncol,:pver,:) + wtdlf2(:ncol,:pver,:)
   wtdif(:ncol,:pver,:) = wtdif(:ncol,:pver,:) + wtdif2(:ncol,:pver,:)
   rliq(:ncol) = rliq(:ncol) + rliq2(:ncol)
   rice(:ncol) = rice(:ncol) + rice2(:ncol)
   
   ftem(:ncol,:pver) = ptend%s(:ncol,:pver)/cpair
   call outfld('CMFDT   ',ftem          ,pcols   ,lchnk   )
   call outfld('CMFDQ   ',ptend%q(1,1,1),pcols   ,lchnk   )
   call t_stopf('cmfmca')
   call physics_update (state, tend, ptend, ztodt)
   call wtrc_check('cmfmca', ncol, state%q)
   call wtrc_chkdelta('ex cmfmca', ncol, state%q)

   call wtrc_qchk2('TPH_W_CMFout','state%q',ncol,state%q(:,:,ixh2oq),state%q(:,:,1))
!!   call wtrc_qchk2('TPH_W_CMFout','state%l',ncol,state%q(:,:,ixh2oq+1),state%q(:,:,2))
!!   call wtrc_qchk2('TPH_W_CMFout','state%i',ncol,state%q(:,:,ixh2oq+2),state%q(:,:,3))
   call wtrc_qchk2('TPH_W_CMFout','ptend%q',ncol,ptend%q(:,:,ixh2oq),ptend%q(:,:,1))
   call wtrc_qchk2('TPH_W_CMFout','dqr   ',ncol,wtdqr2(:,:,ixh2oq),cmfdqr2(:,:))
   call wtrc_qchk2('TPH_W_CMFout','dlf    ',ncol,wtdlf2(:,:,ixh2oq),dlf2(:,:))
   call wtrc_qchk2('TPH_W_CMFout','dif    ',ncol,wtdif2(:,:,ixh2oq),dif2(:,:))
!
! Determine the phase of the precipitation produced and add latent heat of fusion
   call zm_conv_evap(state, ptend, cmfdqr2, cld, ztodt, &
                     flxprec, flxsnow, fevprec, fevsnow, &
                     prec_cmf, snow_cmf, .true.)

! Water tracers: evaporat falling precip, and assign state
   if (trace_water .and. lwtrcevp) then
     call wtrc_zm_conv_evap(state  , ptend  , wtdqr2, cld , ztodt , &
                            flxprec, flxsnow, fevprec, fevsnow, &
                            wtprec_cmf, wtsnow_cmf, .true. )
   end if
!
   call physics_update(state, tend, ptend, ztodt)
   call wtrc_check('hk_evap', ncol, state%q)
   flx_cnd(:ncol) = prec_cmf(:ncol) + rliq2(:ncol)
!!   flx_ice(:ncol) = snow_cmf(:ncol) + rice2(:ncol)
   flx_ice(:ncol) = snow_cmf(:ncol) 
   call check_energy_chng(state, tend, "hk_evap", nstep, ztodt, zero, flx_cnd, flx_ice, zero)
!
! Merge shallow/mid-level output with prior results from Zhang-McFarlane
!
   do i=1,ncol
      if (cnt2(i) < cnt(i)) cnt(i) = cnt2(i)
      if (cnb2(i) > cnb(i)) cnb(i) = cnb2(i)
   end do
!
   cmfmc (:ncol,:pver) = cmfmc (:ncol,:pver) + cmfmc2 (:ncol,:pver)
   cmfsl (:ncol,:pver) = cmfsl (:ncol,:pver) + cmfsl2 (:ncol,:pver)
   cmflq (:ncol,:pver) = cmflq (:ncol,:pver) + cmflq2 (:ncol,:pver)
   call outfld('CMFMC' , cmfmc  , pcols, lchnk)
!  output new partition of cloud condensate variables, as well as precipitation 
   call outfld('QC      ',dlf2+dif2      ,pcols   ,lchnk   )
   call outfld('PRECSH  ',prec_cmf      ,pcols   ,lchnk       )
   call outfld('CMFDQR', cmfdqr2, pcols, lchnk)
   call outfld('CMFSL' , cmfsl  , pcols, lchnk)
   call outfld('CMFLQ' , cmflq  , pcols, lchnk)
   call outfld('DQP'   , dlf2    , pcols, lchnk)
   !call outfld('DQP_name   , dif2    , pcols, lchnk)

! Water isotope tracers: equilibrate cloud liquid with gridscale water.
! Here so sedimentation precipitation is in equilibrium with vapour. (dcn)
! (This is an "isotope-only" routine.)
   if (trace_water .and. wisotope .and. lwisoeql) then
      call t_startf('wiso_cldliq_equilibrate1')
      call wiso_cldliq_equilibrate(state, ptend, ztodt)
      call physics_update(state, tend, ptend, ztodt)
      call wtrc_check('wiso_cldliq_eq1', ncol, state%q)
      call wtrc_chkdelta('post eql1', ncol, state%q)
      call t_stopf('wiso_cldliq_equilibrate1')
   end if

! Allow the cloud liquid drops and ice particles to sediment
! Occurs before adding convectively detrained cloud water, because the phase of the
! of the detrained water is unknown.
   call t_startf('cldwat_sediment')
   call cldcond_sediment(state, ptend, ztodt,cld, icefrac, landfrac, ocnfrac, prec_sed, &
                         snow_sed, landm, snowh, fxliq, fxice)

! Water tracers: allow cloud liquid and ice to sediment 
   if (trace_water .and. lwtrcsed) then
     call t_startf('wtrc_cldcond_sediment')
     call wtrc_cldcond_sediment(state, ptend, ztodt, cld, fxliq, fxice, wtprec_sed, wtsnow_sed)
     call t_stopf('wtrc_cldcond_sediment')
   end if

   call physics_update(state, tend, ptend, ztodt)
   call wtrc_check('cldwat_sediment', ncol, state%q)
   call wtrc_chkdelta('post cldwat_sediment', ncol, state%q)
   call t_stopf('cldwat_sediment')

! check energy integrals
   call check_energy_chng(state, tend, "cldwat_sediment", nstep, ztodt, zero, prec_sed, snow_sed, zero)

! Water isotope tracers: equilibrate cloud liquid with gridscale water.
! Here so state is in equilibrium after phase change during sedimentation. 
! (This is an "isotope-only" routine.)
   if (trace_water .and. wisotope .and. lwisoeql) then
      call t_startf('wiso_cldliq_equilibrate2')
      call wiso_cldliq_equilibrate(state, ptend, ztodt)
      call physics_update(state, tend, ptend, ztodt)
      call wtrc_check('wiso_cldliq_eq2', ncol, state%q)
      call wtrc_chkdelta('post eql2', ncol, state%q)
      call t_stopf('wiso_cldliq_equilibrate2')
   end if

! Put the detraining cloud water from convection into the cloud and environment. 
   call cldcond_zmconv_detrain(dlf, dif, cld, state, ptend)

! Water tracers: put detraining tracers into cloud and environment
   if (trace_water .and. lwtrcdet) then
     call wtrc_cldcond_zmconv_detrain(wtdlf, wtdif, cld, state, ptend)
   end if

   call physics_update(state, tend, ptend, ztodt)
   call wtrc_check('cldwat_detrain', ncol, state%q)

! check energy integrals, reserved liquid has now been used
   flx_cnd(:ncol) = -rliq(:ncol)
   call check_energy_chng(state, tend, "cldwat_detrain", nstep, ztodt, zero, flx_cnd, zero, zero)
!
! cloud fraction after transport and convection,
! derive the relationship between rh and cld from 
! the employed cloud scheme
!
   call cldnrh(lchnk,   ncol,                                &
               state%pmid,    state%t,   state%q(1,1,1),   state%omega, &
               cnt,     cnb,     cld,    clc,     state%pdel,   &
               cmfmc,   cmfmc2, landfrac,snowh,   concld,  cldst,    &
               ts,      sst, state%pint(1,pverp),       zdu,  ocnfrac, &
               rhdfda,   rhu00 , state%phis)
   call outfld('CONCLD  ',concld, pcols,lchnk)
   call outfld('CLDST   ',cldst,  pcols,lchnk)
   call outfld('CNVCLD  ',clc,    pcols,lchnk)

! Water tracers: precompute tracer ratios before repartitioning in cldcond_tend
! (needed as repoartitioning not done via ptend and physics_update)
    if (trace_water) then
       call t_startf('wtrc_ratio_all')
       call wtrc_ratio_all(ncol,state%q, Rtrc)
       call t_stopf('wtrc_ratio_all')
    end if

! cloud water and ice parameterizations
   call t_startf('cldwat_tend')
   call cldcond_tend(state, ptend, ztodt, &
       tcwat, qcwat, lcwat, prec_pcw, snow_pcw, icefrac, rhdfda, rhu00, cld, &
       nevapr, prain, nevaps, psnow, qme, snowh, dqlfz)

! Water tracers: tendencies due to cloud liquid and ice exchange
   if (trace_water .and. lwtrcpcw) then
     call t_startf('wtrc_cldwat_tend')
     wtprec_pcw(:ncol,:) = 0.
     wtsnow_pcw(:ncol,:) = 0.
     wtprec_pcw(:ncol,1) = prec_pcw(:ncol)
     wtsnow_pcw(:ncol,1) = snow_pcw(:ncol)
     call wtrc_cldcond_tend(state, ptend, ztodt, Rtrc, dqlfz, &
           nevapr, prain, nevaps, psnow, wtprec_pcw, wtsnow_pcw)
     call t_stopf('wtrc_cldwat_tend')
   end if

   call physics_update (state, tend, ptend, ztodt)
   call wtrc_check('cldwat_tend', ncol, state%q)
   call wtrc_chkdelta('post wtrc_cldwat_tend', ncol, state%q)
   call t_stopf('cldwat_tend')

! check energy integrals
   call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_pcw, snow_pcw, zero)

! Save off q and t after cloud water
   do k=1,pver
      qcwat(:ncol,k) = state%q(:ncol,k,1)
      tcwat(:ncol,k) = state%t(:ncol,k)
      lcwat(:ncol,k) = state%q(:ncol,k,ixcldice) + state%q(:ncol,k,ixcldliq)
   end do
!
!  aerosol wet chemistry determines scavenging fractions, and transformations
   call get_lat_all_p(lchnk, ncol, lat)
   call get_lon_all_p(lchnk, ncol, lon)
   call get_rlat_all_p(lchnk, ncol, clat)
   conicw(:ncol,:) = icwmr1(:ncol,:) + icwmr2(:ncol,:)
   cmfdqrt(:ncol,:) = zmrprd(:ncol,:) + cmfdqr2(:ncol,:)
   call aerosol_wet_intr (state, ptend, cflx, nstep, ztodt, lat, clat, lon,&
        qme, prain, &
       nevapr, cldc, cld, fracis, calday, cmfdqrt, conicw)
   call physics_update (state, tend, ptend, ztodt)

!
!     Convective transport of all trace species except water vapor and
!     cloud liquid and ice done here because we need to do the scavenging first
!     to determine the interstitial fraction.
!
   ptend%name  = 'convtran2'
   ptend%lq(:) = .true.
   ptend%lq(ixcldice) = .false.
   ptend%lq(ixcldliq) = .false.
   if (trace_water) then
     do m = ixwti, ixwtx
       if (wtrc_is_wtrc(m)) &
          ptend%lq(m) = .false.
     end do
   end if
   dpdry = 0
   do i = 1,lengath
      dpdry(i,:) = state%pdeldry(ideep(i),:)/100.
   end do
   call convtran (lchnk,                                           &
                  ptend%lq,state%q, ppcnst,     mu2,     md2,      &
                  du2,     eu2,     ed2,        dp,      dsubcld,  &
                  jt,      maxg,    ideep,      1,       lengath,  &
                  nstep,   fracis,  ptend%q,    dpdry)

   call physics_update (state, tend, ptend, ztodt)

! check tracer integrals
   call check_tracers_chng(state, tracerint, "cmfmca", nstep, ztodt, cflx)

!
! Compute rates of temperature and constituent change due to moist processes
!
   dtcond(:ncol,:) = (state%s(:ncol,:) - dtcond(:ncol,:))*rtdt / cpair
   dqcond(:ncol,:,:) = (state%q(:ncol,:,:) - dqcond(:ncol,:,:))*rtdt
   call outfld('DTCOND  ',dtcond, pcols   ,lchnk   )
   do m=1,ppcnst
      call outfld(dcconnam(m),dqcond(1,1,m),pcols   ,lchnk )
   end do

! Compute total convective and stratiform precipitation and snow rates
   do i=1,ncol
      precc (i) = prec_zmc(i) + prec_cmf(i)
      precl (i) = prec_sed(i) + prec_pcw(i)
      precsc(i) = snow_zmc(i) + snow_cmf(i)
      precsl(i) = snow_sed(i) + snow_pcw(i)
! jrm These checks should not be necessary if they exist in the parameterizations
      if(precc(i).lt.0.) precc(i)=0.
      if(precl(i).lt.0.) precl(i)=0.
      if(precsc(i).lt.0.) precsc(i)=0.
      if(precsl(i).lt.0.) precsl(i)=0.
      if(precsc(i).gt.precc(i)) precsc(i)=precc(i)
      if(precsl(i).gt.precl(i)) precsl(i)=precl(i)
! end jrm
   end do
   prcsnw(:ncol) = precsc(:ncol) + precsl(:ncol)   ! total snowfall rate: needed by slab ocean model
!
! Water tracers: do all these again for the tracer quantities
! (were someone bold, the regular "prec" bvariable could be dimensioned "ppcnst" everywhere)
!
   if (trace_water) then
     wtprecc (:ncol,1) = precc(:ncol)
     wtprecl (:ncol,1) = precl(:ncol)
     wtprecsc(:ncol,1) = precsc(:ncol)
     wtprecsl(:ncol,1) = precsl(:ncol)
     wtprcsnw(:ncol,1) = prcsnw(:ncol)

     do m = ixwti,ixwtx
       if (wtrc_is_vap(m)) then
         do i=1,ncol
            wtprecc (i,m) = wtprec_zmc(i,m) + wtprec_cmf(i,m)
            wtprecl (i,m) = wtprec_sed(i,m) + wtprec_pcw(i,m)
            wtprecsc(i,m) = wtsnow_zmc(i,m) + wtsnow_cmf(i,m)
            wtprecsl(i,m) = wtsnow_sed(i,m) + wtsnow_pcw(i,m)
            if(wtprecc(i,m).lt.0.) wtprecc(i,m)=0.
            if(wtprecl(i,m).lt.0.) wtprecl(i,m)=0.
            if(wtprecsc(i,m).lt.0.) wtprecsc(i,m)=0.
            if(wtprecsl(i,m).lt.0.) wtprecsl(i,m)=0.
            if(wtprecsc(i,m).gt.wtprecc(i,m)) wtprecsc(i,m)=wtprecc(i,m)
            if(wtprecsl(i,m).gt.wtprecl(i,m)) wtprecsl(i,m)=wtprecl(i,m)
         end do
         wtprcsnw(:ncol,m) = wtprecsc(:ncol,m) + wtprecsl(:ncol,m)
       end if
     end do
   end if
!
!===================================================
! Moist physical parameteriztions complete: 
! send dynamical variables, and derived variables to history file
!===================================================
!
   call diag_dynvar (lchnk, ncol, state)
!
! Output water tracer histogram output
!
   if (trace_water .and. wisotope) then
     call wiso_hist_diag(lchnk,ncol, state, cld)
   end if
!
!===================================================
! Radiation computations
!===================================================
!
! Cosine solar zenith angle for current time step
!
   call get_rlat_all_p(lchnk, ncol, clat)
   call get_rlon_all_p(lchnk, ncol, clon)
   call zenith (calday, clat, clon, coszrs, ncol)

   if (dosw .or. dolw) then

! Compute cloud water/ice paths and optical properties for input to radiation
      call t_startf('cldoptics')
      call param_cldoptics_calc(state, cld, landfrac, landm,icefrac, &
                                cicewp, cliqwp, emis, rel, rei, pmxrgn, nmxrgn, snowh)
      call t_stopf('cldoptics')
!
! Complete radiation calculations
!
      call t_startf ('radctl')
      call radctl (lchnk, ncol, lwup, emis, state%pmid,             &
                   state%pint, state%lnpmid, state%lnpint, state%t, state%q,   &
                   cld, cicewp, cliqwp, coszrs, asdir, asdif,               &
                   aldir, aldif, pmxrgn, nmxrgn, fsns, fsnt    ,flns    ,flnt    , &
                   qrs, qrl, flwds, rel, rei,                       &
                   sols, soll, solsd, solld,                  &
                   landfrac, state%zm, state, fsds)
      call t_stopf ('radctl')
!
! Cloud cover diagnostics
! radctl can change pmxrgn and nmxrgn so cldsav needs to follow 
! radctl.
!
      call cldsav (lchnk, ncol, cld, state%pmid, cltot, &
                   cllow, clmed, clhgh, nmxrgn, pmxrgn)
!
! Dump cloud field information to history tape buffer (diagnostics)
!
      call outfld('CLDTOT  ',cltot  ,pcols,lchnk)
      call outfld('CLDLOW  ',cllow  ,pcols,lchnk)
      call outfld('CLDMED  ',clmed  ,pcols,lchnk)
      call outfld('CLDHGH  ',clhgh  ,pcols,lchnk)
      call outfld('CLOUD   ',cld    ,pcols,lchnk)
      if (doisccp) then
         call cloudsimulator_run(state, ts, concld, cld, cliqwp, &
                                 cicewp, rel, rei, emis, coszrs  )
      end if
   else

! convert radiative heating rates from Q*dp to Q for energy conservation
      if (conserve_energy) then
         do k =1 , pver
            do i = 1, ncol
               qrs(i,k) = qrs(i,k)/state%pdel(i,k)
               qrl(i,k) = qrl(i,k)/state%pdel(i,k)
            end do
         end do
      end if
         
   end if

!
! Compute net flux
! Since fsns, fsnt, flns, and flnt are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
   do i=1,ncol
      tend%flx_net(i) = fsnt(i) - fsns(i) - flnt(i) + flns(i)
   end do
!
! Compute net radiative heating
!
   call radheat_net (state, ptend, qrl, qrs)
!
! Add radiation tendencies to cummulative model tendencies and update profiles
!
   call physics_update(state, tend, ptend, ztodt)

! check energy integrals
   call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, tend%flx_net)
!
! Compute net surface radiative flux for use by surface temperature code.
! Note that units have already been converted to mks in RADCTL.  Since
! fsns and flwds are in the buffer, array values will be carried across
! timesteps when the radiation code is not invoked.
!
   srfrad(:ncol) = fsns(:ncol) + flwds(:ncol)
   call outfld('SRFRAD  ',srfrad,pcols,lchnk)
!
! Save atmospheric fields to force surface models
!
   call wtrc_check('srfxfer', ncol, state%q)
   qbot_tmp(:ncol,:) = state%q(:ncol,pver,:)
!
   call srfxfer (lchnk, ncol, state%ps, state%u(1,pver), state%v(1,pver),    &
                 state%t(1,pver), qbot_tmp, state%exner(1,pver), state%zm(1,pver), &
                    state%pmid,      &
                 state%rpdel(1,pver))

!---------------------------------------------------------------------------------------
! Save history variables. These should move to the appropriate parameterization interface
!---------------------------------------------------------------------------------------

   call outfld('PRECL   ',precl   ,pcols   ,lchnk       )
   call outfld('PRECC   ',precc   ,pcols   ,lchnk       )
   call outfld('PRECSL  ',precsl  ,pcols   ,lchnk       )
   call outfld('PRECSC  ',precsc  ,pcols   ,lchnk       )
   
   prect(:ncol) = precc(:ncol) + precl(:ncol)
   call outfld('PRECT   ',prect   ,pcols   ,lchnk       )
   call outfld('PRECTMX ',prect   ,pcols   ,lchnk       )
!
! Output water tracer precipitation rates (dcn)
!
   if (trace_water) then
     do m = ixwti,ixwtx
       if (wtrc_is_vap(m)) then
         call outfld('PRL'//trim(cnst_name(m)),wtprecl(:,m) ,pcols ,lchnk       )
         call outfld('PRC'//trim(cnst_name(m)),wtprecc(:,m) ,pcols ,lchnk       )
         call outfld('PSL'//trim(cnst_name(m)),wtprecsl(:,m),pcols ,lchnk       )
         call outfld('PSC'//trim(cnst_name(m)),wtprecsc(:,m),pcols ,lchnk       )

         wtprect(:ncol,m) = wtprecc(:ncol,m) + wtprecl(:ncol,m)
         call outfld('PRT'//trim(cnst_name(m)),wtprect(:,m) ,pcols ,lchnk       )
         call outfld('PRX'//trim(cnst_name(m)),wtprect(:,m) ,pcols ,lchnk       )
       endif
     end do
   end if

#if ( defined COUP_CSM )
   call outfld('PRECLav ',precl   ,pcols   ,lchnk   )
   call outfld('PRECCav ',precc   ,pcols   ,lchnk   )
#endif

#if ( defined BFB_CAM_SCAM_IOP )
   call outfld('Prec   ',prect   ,pcols   ,lchnk       )
#endif
!     
! Compute heating rate for dtheta/dt 
!
   do k=1,pver
      do i=1,ncol
         ftem(i,k) = (qrs(i,k) + qrl(i,k))/cpair * (1.e5/state%pmid(i,k))**cappa
      end do
   end do
   call outfld('HR      ',ftem    ,pcols   ,lchnk   )

! convert radiative heating rates to Q*dp for energy conservation
   if (conserve_energy) then
      do k =1 , pver
         do i = 1, ncol
            qrs(i,k) = qrs(i,k)*state%pdel(i,k)
            qrl(i,k) = qrl(i,k)*state%pdel(i,k)
         end do
      end do
   end if
!
   call wtrc_check('tphysbc_end', ncol, state%q)

   return
end subroutine tphysbc
