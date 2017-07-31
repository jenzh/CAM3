#include <misc.h>
#include <params.h>

subroutine tphysac (ztodt,   pblh,    qpert,   tpert,  shf,  &
                    taux,    tauy,    cflx,    sgh,    lhf,  &
                    landfrac,snowh,   tref,    precc,  precl,&
                    precsc, precsl,   state,   tend,    pbuf,&
                    ocnfrac, fsds, icefrac, fv, ram1 )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Tendency physics after coupling to land, sea, and ice models.
! Computes the following:
!   o Radon surface flux and decay (optional)
!   o Vertical diffusion and planetary boundary layer
!   o Multiple gravity wave drag
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: CCM1, CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use ppgrid,             only: pcols, pver
   use chemistry,          only: trace_gas, chem_timestep_tend
   use gw_drag,            only: gw_intr
   use vertical_diffusion, only: vd_intr
   use physics_types,      only: physics_state, physics_tend, physics_ptend, physics_update,    &
                                 physics_ptend_init, physics_dme_adjust
   use phys_buffer,        only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use constituents,       only: ppcnst, qmin
   use tracers,       only: tracers_timestep_tend, set_state_pdry, set_dry_to_wet
   use physconst,          only: zvir, gravit, rhoh2o, latvap,latice, cpair, rair
   use aerosol_intr,       only: aerosol_emis_intr, aerosol_drydep_intr, aerosol_srcsnk_intr
   use check_energy,    only: check_energy_chng
   use check_energy,       only: check_tracers_data, check_tracers_init, check_tracers_chng
   use time_manager,    only: get_nstep
   use abortutils,      only: endrun
   use dycore,          only: dycore_is
   
   use water_tracers,   only: trace_water, wtrc_check, wtrc_qchk1, ixh2oq
   use water_isotopes,  only: wisotope
   use wiso_chemistry,  only: lwiso_ch4ox, lwiso_o3mif, wiso_chem_htodcy_tend, &
                              wiso_chem_ch4ox_tend, wiso_chem_o3mif_tend


   implicit none

#include <comctl.h>
!
! Arguments
!
   real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)
   real(r8), intent(in) :: landfrac(pcols)        ! Land fraction
   real(r8), intent(in) :: icefrac(pcols)         ! ice fraction
   real(r8), intent(in) :: ocnfrac(pcols)         ! Land fraction
   real(r8), intent(in) :: snowh(pcols)           ! snow depth (liquid water equivalent)
   real(r8), intent(in) :: fv(pcols)              ! for dry deposition velocities for dust from land model
   real(r8), intent(in) :: ram1(pcols)            ! for dry deposition velocities for dust from land model
   real(r8), intent(in) :: tref(pcols)            ! 2m air temperature
   real(r8), intent(in) :: precc(pcols)           ! convective precipitation
   real(r8), intent(in) :: precl(pcols)           ! large-scale precipitation
   real(r8), intent(in) :: fsds(pcols)            ! down solar flux
   real(r8), intent(in) :: precsc(pcols)           ! convective snow
   real(r8), intent(in) :: precsl(pcols)           ! large-scale snow
   real(r8), intent(out) :: pblh(pcols)           ! Planetary boundary layer height
   real(r8), intent(out) :: qpert(pcols,ppcnst)   ! Moisture/constit. perturbation (PBL)
   real(r8), intent(out) :: tpert(pcols)          ! Temperature perturbation (PBL)
   real(r8), intent(inout) :: shf(pcols)          ! Sensible heat flux (w/m^2)
   real(r8), intent(in) :: taux(pcols)            ! X surface stress (zonal)
   real(r8), intent(in) :: tauy(pcols)            ! Y surface stress (meridional)
   real(r8), intent(inout) :: cflx(pcols,ppcnst)  ! Surface constituent flux (kg/m^2/s)
   real(r8), intent(in) :: sgh(pcols)             ! Std. deviation of orography for gwd
   real(r8), intent(inout) :: lhf(pcols)          ! Latent heat flux (w/m^2)

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf

   type(check_tracers_data):: tracerint             ! tracer mass integrals and cummulative boundary fluxes

!
!---------------------------Local workspace-----------------------------
!
   type(physics_ptend)     :: ptend               ! indivdual parameterization tendencies

   integer  :: nstep                              ! current timestep number
   real(r8) :: zero(pcols)                        ! array of zeros

   integer :: lchnk                                ! chunk identifier
   integer :: ncol                                 ! number of atmospheric columns
   integer i,k,m                 ! Longitude, level indices
   integer :: yr, mon, day, tod       ! components of a date

   logical :: labort                            ! abort flag

   real(r8) tvm(pcols,pver)           ! virtual temperature
   real(r8) prect(pcols)              ! total precipitation
   real(r8) surfric(pcols)              ! surface friction velocity
   real(r8) obklen(pcols)             ! Obukhov length
   real(r8) :: fh2o(pcols)            ! h2o flux to balance source from methane chemistry

! physics buffer fields for total energy and mass adjustment
   integer itim, ifld
   real(r8), pointer, dimension(:  ) :: teout
   real(r8), pointer, dimension(:,:) :: qini
   real(r8), pointer, dimension(:,:) :: tini
   real(r8), pointer, dimension(:,:) :: cld
!
   real(r8) qbot_tmp(pcols,ppcnst)

!
!-----------------------------------------------------------------------
   lchnk = state%lchnk
   ncol  = state%ncol

! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('TEOUT')
   teout => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,itim)
   ifld = pbuf_get_fld_idx('QINI')
   qini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('TINI')
   tini  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   ifld = pbuf_get_fld_idx('CLD')
   cld => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)
!
! accumulate fluxes into net flux array
! jrm Include latent heat of fusion for snow
!
   do i=1,ncol
      tend%flx_net(i) = tend%flx_net(i) + shf(i) + (precc(i) + precl(i))*latvap*rhoh2o &
           + (precsc(i) + precsl(i))*latice*rhoh2o
   end do

! Check water tracers on input 
   if (trace_water) then
      call wtrc_check('tphysac_start', ncol, state%q)
   end if

! Initialize parameterization tendency structure

   call physics_ptend_init(ptend)

! emission of aerosols at surface
   call aerosol_emis_intr (state, ptend, cflx, ztodt)
   call physics_update (state, tend, ptend, ztodt)

! get nstep and zero array for energy checker
   zero = 0.
   nstep = get_nstep()
   call check_tracers_init(state, tracerint)

! Check if latent heat flux exceeds the total moisture content of the
! lowest model layer, thereby creating negative moisture.

   call wtrc_qchk1('Tphysac_qnegin','cflx',ncol,cflx(:,ixh2oq),cflx(:,1))
   qbot_tmp(:ncol,:) = state%q(:ncol,pver,:)
   call qneg4('TPHYSAC '       ,lchnk               ,ncol  ,ztodt ,          &
              qbot_tmp       ,state%rpdel(1,pver) ,shf ,lhf ,cflx)
   call wtrc_qchk1('Tphysac_qnegout','cflx',ncol,cflx(:,ixh2oq),cflx(:,1))
   if (trace_water) then
      call wtrc_check('tphysac_qneg4', ncol, state%q)
   end if

!===================================================
! Source/sink terms for advected tracers.
!===================================================
! Test tracers

      call tracers_timestep_tend(state, ptend, cflx, landfrac, ztodt)      
      call physics_update (state, tend, ptend, ztodt)
      call check_tracers_chng(state, tracerint, "tracers_timestep_tend", nstep, ztodt, cflx)

! Advected greenhouse trace gases:

   if (trace_gas) then
      call chem_timestep_tend(state, ptend, cflx, ztodt, fh2o)
! Water tracer source from methane oxidation (optional)
      if (trace_water .and. lwiso_ch4ox) then
        call wiso_chem_ch4ox_tend(state, ptend, ztodt)
      end if
!
      call physics_update (state, tend, ptend, ztodt)
!
! Check energy integrals
      call check_energy_chng(state, tend, "chem", nstep, ztodt, fh2o, zero, zero, zero)
! Check tracer integrals
      call check_tracers_chng(state, tracerint, "chem_timestep_tend", nstep, ztodt, cflx)
   end if
!
! Water isotope tracers: 
!
   if (trace_water .and. wisotope) then
!    Radioactive decay of HTO
     call wiso_chem_htodcy_tend(state, ptend, ztodt)
     call physics_update (state, tend, ptend, ztodt)

!    17O and O18 anomaly from O3 mass independent fractionation (optional, and miniscule)
     if (lwiso_o3mif) then
        call wiso_chem_o3mif_tend(state, ptend)
        call physics_update (state, tend, ptend, ztodt)
     end if
   end if

!===================================================
! Vertical diffusion/pbl calculation
! Call vertical diffusion code (pbl, free atmosphere and molecular)
!===================================================

   call wtrc_qchk1('Tphys_vdintr','cflx',ncol,cflx(:,ixh2oq),cflx(:,1))
   call wtrc_check('tphysac_prevd', ncol, state%q)

   call vd_intr (ztodt    ,state    ,taux     ,tauy     , shf    ,&
                 cflx     ,pblh     ,tpert    ,qpert    , surfric  ,&
                 obklen   ,ptend    ,cld      ,ocnfrac  , landfrac, sgh )

   call physics_update (state, tend, ptend, ztodt)
! Check energy integrals
   call check_energy_chng(state, tend, "vdiff", nstep, ztodt, cflx(:,1), zero, zero, shf)
   call check_tracers_chng(state, tracerint, "vdiff", nstep, ztodt, cflx)
   if (trace_water) then
      call wtrc_check('tphysac_vd_intr', ncol, state%q)
   end if

!  aerosol dry deposition processes
   call aerosol_drydep_intr (state, ptend, cflx(:,1), ztodt, &
       fsds, obklen, tref, surfric, prect, snowh, pblh, shf, landfrac, &
       icefrac, ocnfrac, fv, ram1)

   call physics_update (state, tend, ptend, ztodt)

!===================================================
! Gravity wave drag
!===================================================

   call gw_intr (state   ,sgh     ,pblh    ,ztodt   , ptend , landfrac)
   call physics_update (state, tend, ptend, ztodt)
! Check energy integrals
   call check_energy_chng(state, tend, "gwdrag", nstep, ztodt, zero, zero, zero, zero)

! conversion of aerosols from one category to another by non-wet processes
   call aerosol_srcsnk_intr (state, ptend, ztodt, ocnfrac)
   call physics_update (state, tend, ptend, ztodt)

!-------------- Energy budget checks vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
   teout(:ncol) = state%te_cur(:ncol)

!*** BAB's FV heating kludge *** apply the heating as temperature tendency.
!*** BAB's FV heating kludge *** modify the temperature in the state structure
   state%t(:ncol,:pver) = tini(:ncol,:pver) + ztodt*tend%dtdt(:ncol,:pver)

!
! FV: convert dry-type mixing ratios to moist here because physics_dme_adjust
!     assumes moist. This is done in p_d_coupling for other dynamics. Bundy, Feb 2004.
!
   if ( dycore_is('LR') ) call set_dry_to_wet(state)    ! Physics had dry, dynamics wants moist


! Scale dry mass and energy (does nothing if dycore is EUL or SLD)
   call physics_dme_adjust(state, tend, qini, ztodt)
!!!   REMOVE THIS CALL, SINCE ONLY Q IS BEING ADJUSTED. WON'T BALANCE ENERGY. TE IS SAVED BEFORE THIS
!!!   call check_energy_chng(state, tend, "drymass", nstep, ztodt, zero, zero, zero, zero)
!-------------- Energy budget checks ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   if (aqua_planet) then
      labort = .false.
      do i=1,ncol
         if (ocnfrac(i) /= 1.) labort = .true.
      end do
      if (labort) then
         call endrun ('TPHYSAC error:  grid contains non-ocean point')
      endif
   endif

!
! Check water tracers at end
!
   if (trace_water) then
      call wtrc_check('tphysac_end', ncol, state%q)
   end if

end subroutine tphysac
