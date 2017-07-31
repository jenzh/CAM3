#include <misc.h>
#include <preproc.h>

module clmtype

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clmtype
!
! !DESCRIPTION:
! Define derived type hierarchy. Includes declaration of
! the clm derived type and 1d mapping arrays.
!
!   1 => default
! landunits types can have values of (see clm_varcon.F90)
!          (note shallow lakes not currently implemented)
!   1  => (istsoil) soil (vegetated or bare soil landunit)
!   2  => (istice)  land ice
!   3  => (istdlak) deep lake
!   5  => (istwet)  wetland
!   6  => (isturb)  urban
! column types can have values of
!   1 => in  compete mode
!   pft type values (see below) => in non-compete mode
! pft types can have values of
!   0  => not vegetated
!   1  => needleleaf evergreen temperate tree
!   2  => needleleaf evergreen boreal tree
!   3  => needleleaf deciduous boreal tree
!   4  => broadleaf evergreen tropical tree
!   5  => broadleaf evergreen temperate tree
!   6  => broadleaf deciduous tropical tree
!   7  => broadleaf deciduous temperate tree
!   8  => broadleaf deciduous boreal tree
!   9  => broadleaf evergreen shrub
!   10 => broadleaf deciduous temperate shrub
!   11 => broadleaf deciduous boreal shrub
!   12 => c3 arctic grass
!   13 => c3 non-arctic grass
!   14 => c4 grass
!   15 => corn
!   16 => wheat
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar
!
! !PUBLIC TYPES:
  implicit none
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!*******************************************************************************
!----------------------------------------------------
! Begin definition of conservation check structures
!----------------------------------------------------
! energy balance structure
!----------------------------------------------------
type energy_balance_type
   real(r8), pointer :: errsoi(:)        !soil/lake energy conservation error (W/m**2)
   real(r8), pointer :: errseb(:)        !surface energy conservation error (W/m**2)
   real(r8), pointer :: errsol(:)        !solar radiation conservation error (W/m**2)
   real(r8), pointer :: errlon(:)        !longwave radiation conservation error (W/m**2)
end type energy_balance_type

!----------------------------------------------------
! water balance structure
!----------------------------------------------------
type water_balance_type
   real(r8), pointer :: begwb(:)         !water mass begining of the time step
   real(r8), pointer :: endwb(:)         !water mass end of the time step
   real(r8), pointer :: errh2o(:)        !water conservation error (mm H2O)
end type water_balance_type

!----------------------------------------------------
! carbon balance structure
!----------------------------------------------------
type carbon_balance_type
   real(r8), pointer :: dummy_entry(:)
end type carbon_balance_type

!----------------------------------------------------
! nitrogen balance structure
!----------------------------------------------------
type nitrogen_balance_type
   real(r8), pointer :: dummy_entry(:)
end type nitrogen_balance_type
!----------------------------------------------------
! End definition of conservation check structures
!----------------------------------------------------
!*******************************************************************************

!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the pft_type level
!----------------------------------------------------
! pft physical state variables structure
!----------------------------------------------------
type pft_pstate_type
   integer , pointer :: frac_veg_nosno(:)       !fraction of vegetation not covered by snow (0 OR 1) [-]
   integer , pointer :: frac_veg_nosno_alb(:)   !fraction of vegetation not covered by snow (0 OR 1) [-]
   real(r8), pointer :: rsw(:)                  !soil water content for root zone
   real(r8), pointer :: emv(:)                  !vegetation emissivity
   real(r8), pointer :: z0mv(:)                 !roughness length over vegetation, momentum [m]
   real(r8), pointer :: z0hv(:)                 !roughness length over vegetation, sensible heat [m]
   real(r8), pointer :: z0qv(:)                 !roughness length over vegetation, latent heat [m]
   real(r8), pointer :: rootfr(:,:)             !fraction of roots in each soil layer  (nlevsoi)
   real(r8), pointer :: rootr(:,:)              !effective fraction of roots in each soil layer  (nlevsoi)
   real(r8), pointer :: dewmx(:)                !Maximum allowed dew [mm]
   real(r8), pointer :: rssun(:)                !sunlit stomatal resistance (s/m)
   real(r8), pointer :: rssha(:)                !shaded stomatal resistance (s/m)
   real(r8), pointer :: laisun(:)               !sunlit leaf area
   real(r8), pointer :: laisha(:)               !shaded leaf area
   real(r8), pointer :: btran(:)                !transpiration wetness factor (0 to 1)
   real(r8), pointer :: fsun(:)                 !sunlit fraction of canopy
   real(r8), pointer :: tlai(:)                 !one-sided leaf area index, no burying by snow
   real(r8), pointer :: tsai(:)                 !one-sided stem area index, no burying by snow
   real(r8), pointer :: elai(:)                 !one-sided leaf area index with burying by snow
   real(r8), pointer :: esai(:)                 !one-sided stem area index with burying by snow
   real(r8), pointer :: igs(:)                  !growing season index (0=off, 1=on)
   real(r8), pointer :: stembio(:)              !stem biomass (kg /m**2)
   real(r8), pointer :: rootbio(:)              !root biomass (kg /m**2)
   real(r8), pointer :: fwet(:)                 !fraction of canopy that is wet (0 to 1)
   real(r8), pointer :: fdry(:)                 !fraction of foliage that is green and dry [-] (new)
   real(r8), pointer :: dt_veg(:)               !change in t_veg, last iteration (Kelvin)
   real(r8), pointer :: htop(:)                 !canopy top (m)
   real(r8), pointer :: hbot(:)                 !canopy bottom (m)
   real(r8), pointer :: z0m(:)                  !momentum roughness length (m)
   real(r8), pointer :: displa(:)               !displacement height (m)
   real(r8), pointer :: albd(:,:)               !surface albedo (direct)                       (numrad)
   real(r8), pointer :: albi(:,:)               !surface albedo (diffuse)                      (numrad)
   real(r8), pointer :: fabd(:,:)               !flux absorbed by veg per unit direct flux     (numrad)
   real(r8), pointer :: fabi(:,:)               !flux absorbed by veg per unit diffuse flux    (numrad)
   real(r8), pointer :: ftdd(:,:)               !down direct flux below veg per unit dir flx   (numrad)
   real(r8), pointer :: ftid(:,:)               !down diffuse flux below veg per unit dir flx  (numrad)
   real(r8), pointer :: ftii(:,:)               !down diffuse flux below veg per unit dif flx  (numrad)
   real(r8), pointer :: u10(:)                  !10-m wind (m/s) (for dust model)
   real(r8), pointer :: ram1(:)                 !aerodynamical resistance (s/m)
   real(r8), pointer :: fv(:)                   !friction velocity (m/s) (for dust model)
end type pft_pstate_type

!----------------------------------------------------
! pft ecophysiological constants structure
!----------------------------------------------------
type pft_epc_type
   integer , pointer :: ncorn(:)                !value for corn
   integer , pointer :: nwheat(:)               !value for wheat
   integer , pointer :: noveg(:)                !value for not vegetated
   integer , pointer :: ntree(:)                !value for last type of tree
   real(r8), pointer :: smpmax(:)               !Wilting point potential in mm
   real(r8), pointer :: foln(:)                 !foliage nitrogen (%)
   real(r8), pointer :: dleaf(:)                !characteristic leaf dimension (m)
   real(r8), pointer :: c3psn(:)                !photosynthetic pathway: 0. = c4, 1. = c3
   real(r8), pointer :: vcmx25(:)               !max rate of carboxylation at 25C (umol CO2/m**2/s)
   real(r8), pointer :: mp(:)                   !slope of conductance-to-photosynthesis relationship
   real(r8), pointer :: qe25(:)                 !quantum efficiency at 25C (umol CO2 / umol photon)
   real(r8), pointer :: xl(:)                   !leaf/stem orientation index
   real(r8), pointer :: rhol(:,:)               !leaf reflectance: 1=vis, 2=nir   (numrad)
   real(r8), pointer :: rhos(:,:)               !stem reflectance: 1=vis, 2=nir   (numrad)
   real(r8), pointer :: taul(:,:)               !leaf transmittance: 1=vis, 2=nir (numrad)
   real(r8), pointer :: taus(:,:)               !stem transmittance: 1=vis, 2=nir (numrad)
   real(r8), pointer :: z0mr(:)                 !ratio of momentum roughness length to canopy top height (-)
   real(r8), pointer :: displar(:)              !ratio of displacement height to canopy top height (-)
   real(r8), pointer :: roota_par(:)            !CLM rooting distribution parameter [1/m]
   real(r8), pointer :: rootb_par(:)            !CLM rooting distribution parameter [1/m]
   real(r8), pointer :: sla(:)                  !specific leaf area [m2 leaf g-1 carbon]
end type pft_epc_type

!----------------------------------------------------
! pft DGVM-specific ecophysiological constants structure
!----------------------------------------------------
type pft_dgvepc_type
   real(r8), pointer :: respcoeff(:)       !maintenance respiration coefficient [-]
   real(r8), pointer :: flam(:)            !flammability threshold [units?]
   real(r8), pointer :: resist(:)          !fire resistance index [units?]
   real(r8), pointer :: l_turn(:)          !leaf turnover period [years]
   real(r8), pointer :: l_long(:)          !leaf longevity [years]
   real(r8), pointer :: s_turn(:)          !sapwood turnover period [years]
   real(r8), pointer :: r_turn(:)          !root turnover period [years]
   real(r8), pointer :: l_cton(:)          !leaf C:N (mass ratio)
   real(r8), pointer :: s_cton(:)          !sapwood C:N (mass ratio)
   real(r8), pointer :: r_cton(:)          !root C:N (mass ratio)
   real(r8), pointer :: l_morph(:)         !leaf morphology: 1=broad, 2=needle, 3=grass
   real(r8), pointer :: l_phen(:)          !leaf phenology: 1=everg, 2=summerg, 3=raing, 4=any
   real(r8), pointer :: lmtorm(:)          !leaf:root ratio under non-water stressed conditions
   real(r8), pointer :: crownarea_max(:)   !tree maximum crown area [m2]
   real(r8), pointer :: init_lai(:)        !sapling (or initial grass) LAI [-]
   real(r8), pointer :: x(:)               !sapling: (heart+sapwood)/sapwood [-]
   real(r8), pointer :: tcmin(:)           !minimum coldest monthly mean temperature [units?]
   real(r8), pointer :: tcmax(:)           !maximum coldest monthly mean temperature [units?]
   real(r8), pointer :: gddmin(:)          !minimum growing degree days (at or above 5 C)
   real(r8), pointer :: twmax(:)           !upper limit of temperature of the warmest month [units?]
   real(r8), pointer :: lm_sapl(:)
   real(r8), pointer :: sm_sapl(:)
   real(r8), pointer :: hm_sapl(:)
   real(r8), pointer :: rm_sapl(:)
   logical , pointer :: tree(:)
   logical , pointer :: summergreen(:)
   logical , pointer :: raingreen(:)
   real(r8), pointer :: reinickerp(:)      !parameter in allometric equation
   real(r8), pointer :: wooddens(:)        !wood density (gC/m3)
   real(r8), pointer :: latosa(:)          !ratio of leaf area to sapwood cross-sectional area (Shinozaki et al 1964a,b)
   real(r8), pointer :: allom1(:)          !parameter in allometric
   real(r8), pointer :: allom2(:)          !parameter in allometric
   real(r8), pointer :: allom3(:)          !parameter in allometric
end type pft_dgvepc_type

!----------------------------------------------------
! pft energy state variables structure
!----------------------------------------------------
type pft_estate_type
   real(r8), pointer :: t_ref2m(:)          !2 m height surface air temperature (Kelvin)
   real(r8), pointer :: t_ref2m_min(:)      !daily minimum of average 2 m height surface air temperature (K)
   real(r8), pointer :: t_ref2m_max(:)      !daily maximum of average 2 m height surface air temperature (K)
   real(r8), pointer :: t_ref2m_min_inst(:) !instantaneous daily min of average 2 m height surface air temp (K)
   real(r8), pointer :: t_ref2m_max_inst(:) !instantaneous daily max of average 2 m height surface air temp (K)
   real(r8), pointer :: q_ref2m(:)          !2 m height surface specific humidity (kg/kg)
   real(r8), pointer :: t_veg(:)            !vegetation temperature (Kelvin)
end type pft_estate_type

!----------------------------------------------------
! pft water state variables structure
!----------------------------------------------------
type pft_wstate_type
   real(r8), pointer :: h2ocan(:)         !canopy water (mm H2O)
end type pft_wstate_type

!----------------------------------------------------
! pft carbon state variables structure
!----------------------------------------------------
type pft_cstate_type
   real(r8), pointer :: dummy_entry(:)
end type pft_cstate_type

!----------------------------------------------------
! pft nitrogen state variables structure
!----------------------------------------------------
type pft_nstate_type
   real(r8), pointer :: dummy_entry(:)
end type pft_nstate_type

!----------------------------------------------------
! pft VOC state variables structure
!----------------------------------------------------
type pft_vstate_type
   real(r8), pointer :: dummy_entry(:)
end type pft_vstate_type

!----------------------------------------------------
! pft DGVM state variables structure
!----------------------------------------------------
type pft_dgvstate_type
   real(r8), pointer :: agdd0(:)               !accumulated growing degree days above 0 deg C
   real(r8), pointer :: agdd5(:)               !accumulated growing degree days above -5
   real(r8), pointer :: agddtw(:)              !accumulated growing degree days above twmax
   real(r8), pointer :: agdd(:)                !accumulated growing degree days above 5
   real(r8), pointer :: t10(:)                 !10-day running mean of the 2 m temperature (K)
   real(r8), pointer :: t_mo(:)                !30-day average temperature (Kelvin)
   real(r8), pointer :: t_mo_min(:)            !annual min of t_mo (Kelvin)
   real(r8), pointer :: fnpsn10(:)             !10-day running mean net photosynthesis
   real(r8), pointer :: prec365(:)             !365-day running mean of tot. precipitation
   real(r8), pointer :: agdd20(:)              !20-yr running mean of agdd
   real(r8), pointer :: tmomin20(:)            !20-yr running mean of tmomin
   real(r8), pointer :: t10min(:)              !annual minimum of 10-day running mean (K)
   real(r8), pointer :: tsoi25(:)              !soil temperature to 0.25 m (Kelvin)
   real(r8), pointer :: annpsn(:)              !annual photosynthesis (umol CO2 /m**2)
   real(r8), pointer :: annpsnpot(:)           !annual potential photosynthesis (same units)
   logical , pointer :: present(:)             !whether PFT present in patch
   real(r8), pointer :: dphen(:)               !phenology [0 to 1]
   real(r8), pointer :: leafon(:)              !leafon days
   real(r8), pointer :: leafof(:)              !leafoff days
   real(r8), pointer :: nind(:)                !number of individuals (#/m**2)
   real(r8), pointer :: lm_ind(:)              !individual leaf mass
   real(r8), pointer :: sm_ind(:)              !individual sapwood mass
   real(r8), pointer :: hm_ind(:)              !individual heartwood mass
   real(r8), pointer :: rm_ind(:)              !individual root mass
   real(r8), pointer :: lai_ind(:)             !LAI per individual
   real(r8), pointer :: fpcinc(:)              !foliar projective cover increment (fraction)
   real(r8), pointer :: fpcgrid(:)             !foliar projective cover on gridcell (fraction)
   real(r8), pointer :: crownarea(:)           !area that each individual tree takes up (m^2)
   real(r8), pointer :: bm_inc(:)              !biomass increment
   real(r8), pointer :: afmicr(:)              !annual microbial respiration
   real(r8), pointer :: firelength(:)          !fire season in days
   real(r8), pointer :: litterag(:)            !above ground litter
   real(r8), pointer :: litterbg(:)            !below ground litter
   real(r8), pointer :: cpool_fast(:)          !fast carbon pool
   real(r8), pointer :: cpool_slow(:)          !slow carbon pool
   real(r8), pointer :: k_fast_ave(:)          !decomposition rate
   real(r8), pointer :: k_slow_ave(:)          !decomposition rate
   real(r8), pointer :: litter_decom_ave(:)    !decomposition rate
   real(r8), pointer :: turnover_ind(:)        !
end type pft_dgvstate_type

!----------------------------------------------------
! pft energy flux variables structure
!----------------------------------------------------
type pft_eflux_type
   real(r8), pointer :: sabg(:)              !solar radiation absorbed by ground (W/m**2)
   real(r8), pointer :: sabv(:)              !solar radiation absorbed by vegetation (W/m**2)
   real(r8), pointer :: fsa(:)               !solar radiation absorbed (total) (W/m**2)
   real(r8), pointer :: fsr(:)               !solar radiation reflected (W/m**2)
   real(r8), pointer :: parsun(:)            !average absorbed PAR for sunlit leaves (W/m**2)
   real(r8), pointer :: parsha(:)            !average absorbed PAR for shaded leaves (W/m**2)
   real(r8), pointer :: dlrad(:)             !downward longwave radiation below the canopy [W/m2]
   real(r8), pointer :: ulrad(:)             !upward longwave radiation above the canopy [W/m2]
   real(r8), pointer :: eflx_lh_tot(:)       !total latent heat flux (W/m8*2)  [+ to atm]
   real(r8), pointer :: eflx_lh_grnd(:)      !ground evaporation heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_soil_grnd(:)    !soil heat flux (W/m**2) [+ = into soil]
   real(r8), pointer :: eflx_sh_tot(:)       !total sensible heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_sh_grnd(:)      !sensible heat flux from ground (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_sh_veg(:)       !sensible heat flux from leaves (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_lh_vege(:)      !veg evaporation heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_lh_vegt(:)      !veg transpiration heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: cgrnd(:)             !deriv. of soil energy flux wrt to soil temp [w/m2/k]
   real(r8), pointer :: cgrndl(:)            !deriv. of soil latent heat flux wrt soil temp [w/m**2/k]
   real(r8), pointer :: cgrnds(:)            !deriv. of soil sensible heat flux wrt soil temp [w/m2/k]
   real(r8), pointer :: eflx_gnet(:)         !net heat flux into ground (W/m**2)
   real(r8), pointer :: dgnetdT(:)           !derivative of net ground heat flux wrt soil temp (W/m**2 K)
   real(r8), pointer :: eflx_lwrad_out(:)    !emitted infrared (longwave) radiation (W/m**2)
   real(r8), pointer :: eflx_lwrad_net(:)    !net infrared (longwave) rad (W/m**2) [+ = to atm]
   real(r8), pointer :: fsds_vis_d(:)        !incident direct beam vis solar radiation (W/m**2)
   real(r8), pointer :: fsds_nir_d(:)        !incident direct beam nir solar radiation (W/m**2)
   real(r8), pointer :: fsds_vis_i(:)        !incident diffuse vis solar radiation (W/m**2)
   real(r8), pointer :: fsds_nir_i(:)        !incident diffuse nir solar radiation (W/m**2)
   real(r8), pointer :: fsr_vis_d(:)         !reflected direct beam vis solar radiation (W/m**2)
   real(r8), pointer :: fsr_nir_d(:)         !reflected direct beam nir solar radiation (W/m**2)
   real(r8), pointer :: fsr_vis_i(:)         !reflected diffuse vis solar radiation (W/m**2)
   real(r8), pointer :: fsr_nir_i(:)         !reflected diffuse nir solar radiation (W/m**2)
   real(r8), pointer :: fsds_vis_d_ln(:)     !incident direct beam vis solar radiation at local noon (W/m**2)
   real(r8), pointer :: fsds_nir_d_ln(:)     !incident direct beam nir solar radiation at local noon (W/m**2)
   real(r8), pointer :: fsr_vis_d_ln(:)      !reflected direct beam vis solar radiation at local noon (W/m**2)
   real(r8), pointer :: fsr_nir_d_ln(:)      !reflected direct beam nir solar radiation at local noon (W/m**2)
end type pft_eflux_type

!----------------------------------------------------
! pft momentum flux variables structure
!----------------------------------------------------
type pft_mflux_type
   real(r8),pointer ::  taux(:)           !wind (shear) stress: e-w (kg/m/s**2)
   real(r8),pointer ::  tauy(:)           !wind (shear) stress: n-s (kg/m/s**2)
end type pft_mflux_type

!----------------------------------------------------
! pft water flux variables structure
!----------------------------------------------------
type pft_wflux_type
   real(r8), pointer :: qflx_prec_intr(:) !interception of precipitation [mm/s]
   real(r8), pointer :: qflx_prec_grnd(:) !water onto ground including canopy runoff [kg/(m2 s)]
   real(r8), pointer :: qflx_rain_grnd(:) !rain on ground after interception (mm H2O/s) [+]
   real(r8), pointer :: qflx_snow_grnd(:) !snow on ground after interception (mm H2O/s) [+]
   real(r8), pointer :: qflx_snowcap(:)   !excess precipitation due to snow capping (mm H2O /s) [+]
   real(r8), pointer :: qflx_evap_veg(:)  !vegetation evaporation (mm H2O/s) (+ = to atm)
   real(r8), pointer :: qflx_tran_veg(:)  !vegetation transpiration (mm H2O/s) (+ = to atm)
   real(r8), pointer :: qflx_evap_can(:)  !evaporation from leaves and stems
   real(r8), pointer :: qflx_evap_soi(:)  !soil evaporation (mm H2O/s) (+ = to atm)
   real(r8), pointer :: qflx_evap_tot(:)  !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
   real(r8), pointer :: qflx_evap_grnd(:) !ground surface evaporation rate (mm H2O/s) [+]
   real(r8), pointer :: qflx_dew_grnd(:)  !ground surface dew formation (mm H2O /s) [+]
   real(r8), pointer :: qflx_sub_snow(:)  !sublimation rate from snow pack (mm H2O /s) [+]
   real(r8), pointer :: qflx_dew_snow(:)  !surface dew added to snow pack (mm H2O /s) [+]
end type pft_wflux_type

!----------------------------------------------------
! pft carbon flux variables structure
!----------------------------------------------------
type pft_cflux_type
   real(r8), pointer :: psnsun(:)         !sunlit leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha(:)         !shaded leaf photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: fpsn(:)           !photosynthesis (umol CO2 /m**2 /s)
   real(r8), pointer :: frm(:)            !total maintenance respiration (umol CO2 /m**2/s)
   real(r8), pointer :: frmf(:)           !leaf maintenance respiration  (umol CO2 /m**2 /s)
   real(r8), pointer :: frms(:)           !stem maintenance respiration  (umol CO2 /m**2 /s)
   real(r8), pointer :: frmr(:)           !root maintenance respiration  (umol CO2 /m**2 /s)
   real(r8), pointer :: frg(:)            !growth respiration (umol CO2 /m**2 /s)
   real(r8), pointer :: dmi(:)            !total dry matter production (ug /m**2 /s)
   real(r8), pointer :: fco2(:)           !net CO2 flux (umol CO2 /m**2 /s) [+ = to atm]
   real(r8), pointer :: fmicr(:)          !microbial respiration (umol CO2 /m**2 /s)
end type pft_cflux_type

!----------------------------------------------------
! pft nitrogen flux variables structure
!----------------------------------------------------
type pft_nflux_type
   real(r8), pointer :: dummy_entry(:)
end type pft_nflux_type

!----------------------------------------------------
! pft VOC flux variables structure
!----------------------------------------------------
type pft_vflux_type
   real(r8), pointer :: vocflx_tot(:)     !total VOC flux into atmosphere [ug C m-2 h-1]
   real(r8), pointer :: vocflx(:,:)       !(nvoc) VOC flux [ug C m-2 h-1]
   real(r8), pointer :: vocflx_1(:)       !vocflx(1) (for history output) [ug C m-2 h-1]
   real(r8), pointer :: vocflx_2(:)       !vocflx(2) (for history output) [ug C m-2 h-1]
   real(r8), pointer :: vocflx_3(:)       !vocflx(3) (for history output) [ug C m-2 h-1]
   real(r8), pointer :: vocflx_4(:)       !vocflx(4) (for history output) [ug C m-2 h-1]
   real(r8), pointer :: vocflx_5(:)       !vocflx(5) (for history output) [ug C m-2 h-1]
end type pft_vflux_type

!----------------------------------------------------
! pft dust flux variables structure
!----------------------------------------------------
type pft_dflux_type
   real(r8), pointer :: flx_mss_vrt_dst(:,:)    !(ndst)  !surface dust emission (kg/m**2/s) [ + = to atm]
   real(r8), pointer :: flx_mss_vrt_dst_tot(:)  !total dust flux into atmosphere
   real(r8), pointer :: vlc_trb(:,:)            !(ndst) turbulent deposition velocity (m/s)
   real(r8), pointer :: vlc_trb_1(:)            !turbulent deposition velocity 1(m/s)
   real(r8), pointer :: vlc_trb_2(:)            !turbulent deposition velocity 2(m/s)
   real(r8), pointer :: vlc_trb_3(:)            !turbulent deposition velocity 3(m/s)
   real(r8), pointer :: vlc_trb_4(:)            !turbulent deposition velocity 4(m/s)
end type pft_dflux_type

!----------------------------------------------------
! End definition of structures defined at the pft_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the column_type level
!----------------------------------------------------
! column physical state variables structure
!----------------------------------------------------
type column_pstate_type
   type(pft_pstate_type) :: pps_a            !pft-level pstate variables averaged to the column
   integer , pointer :: snl(:)                !number of snow layers
   integer , pointer :: isoicol(:)            !soil color class
   real(r8), pointer :: bsw(:,:)              !Clapp and Hornberger "b" (nlevsoi)
   real(r8), pointer :: watsat(:,:)           !volumetric soil water at saturation (porosity) (nlevsoi)
   real(r8), pointer :: hksat(:,:)            !hydraulic conductivity at saturation (mm H2O /s) (nlevsoi)
   real(r8), pointer :: sucsat(:,:)           !minimum soil suction (mm) (nlevsoi)
   real(r8), pointer :: csol(:,:)             !heat capacity, soil solids (J/m**3/Kelvin) (nlevsoi)
   real(r8), pointer :: tkmg(:,:)             !thermal conductivity, soil minerals  [W/m-K] (new) (nlevsoi)
   real(r8), pointer :: tkdry(:,:)            !thermal conductivity, dry soil (W/m/Kelvin) (nlevsoi)
   real(r8), pointer :: tksatu(:,:)           !thermal conductivity, saturated soil [W/m-K] (new) (nlevsoi)
   real(r8), pointer :: smpmin(:)             !restriction for min of soil potential (mm) (new)
   real(r8), pointer :: gwc_thr(:)            !threshold soil moisture based on clay content
   real(r8), pointer :: mss_frc_cly_vld(:)    ![frc] Mass fraction clay limited to 0.20
   real(r8), pointer :: mbl_bsn_fct(:)        !??
   logical , pointer :: do_capsnow(:)         !true => do snow capping
   real(r8), pointer :: snowdp(:)             !snow height (m)
   real(r8), pointer :: snowage(:)            !non dimensional snow age [-] (new)
   real(r8), pointer :: frac_sno(:)           !fraction of ground covered by snow (0 to 1)
   real(r8), pointer :: zi(:,:)               !interface level below a "z" level (m) (-nlevsno+0:nlevsoi)
   real(r8), pointer :: dz(:,:)               !layer thickness (m)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: z(:,:)                !layer depth (m) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: frac_iceold(:,:)      !fraction of ice relative to the tot water (new) (-nlevsno+1:nlevsoi)
   integer , pointer :: imelt(:,:)            !flag for melting (=1), freezing (=2), Not=0 (new) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: eff_porosity(:,:)     !effective porosity = porosity - vol_ice (nlevsoi)
   real(r8), pointer :: sfact(:)              !term for implicit correction to evaporation
   real(r8), pointer :: sfactmax(:)           !maximim of "sfact"
   real(r8), pointer :: emg(:)                !ground emissivity
   real(r8), pointer :: z0mg(:)               !roughness length over ground, momentum [m]
   real(r8), pointer :: z0hg(:)               !roughness length over ground, sensible heat [m]
   real(r8), pointer :: z0qg(:)               !roughness length over ground, latent heat [m]
   real(r8), pointer :: htvp(:)               !latent heat of vapor of water (or sublimation) [j/kg]
   real(r8), pointer :: beta(:)               !coefficient of convective velocity [-]
   real(r8), pointer :: zii(:)                !convective boundary height [m]
   real(r8), pointer :: albgrd(:,:)           !ground albedo (direct) (numrad)
   real(r8), pointer :: albgri(:,:)           !ground albedo (diffuse) (numrad)
   real(r8), pointer :: rootr_column(:,:)     !effective fraction of roots in each soil layer (nlevsoi)
   real(r8), pointer :: wf(:)                 !soil water as frac. of whc for top 0.5 m
end type column_pstate_type

!----------------------------------------------------
! column energy state variables structure
!----------------------------------------------------
type column_estate_type
   type(pft_estate_type):: pes_a              !pft-level energy state variables averaged to the column
   real(r8), pointer :: t_grnd(:)             !ground temperature (Kelvin)
   real(r8), pointer :: dt_grnd(:)            !change in t_grnd, last iteration (Kelvin)
   real(r8), pointer :: t_soisno(:,:)         !soil temperature (Kelvin)  (-nlevsno+1:nlevsoi)
   real(r8), pointer :: t_lake(:,:)           !lake temperature (Kelvin)  (1:nlevlak)
   real(r8), pointer :: tssbef(:,:)           !soil/snow temperature before update (-nlevsno+1:nlevsoi)
   real(r8), pointer :: t_snow(:)             !vertically averaged snow temperature
   real(r8), pointer :: thv(:)                !virtual potential temperature (kelvin)
   real(r8), pointer :: thm(:)                !intermediate variable (forc_t+0.0098*forc_hgt_t)
end type column_estate_type

!----------------------------------------------------
! column water state variables structure
!----------------------------------------------------
type column_wstate_type
   type(pft_wstate_type):: pws_a             !pft-level water state variables averaged to the column
   real(r8), pointer :: h2osno(:)             !snow water (mm H2O)
   real(r8), pointer :: h2osoi_liq(:,:)       !liquid water (kg/m2) (new) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: h2osoi_ice(:,:)       !ice lens (kg/m2) (new) (-nlevsno+1:nlevsoi)
   real(r8), pointer :: h2osoi_vol(:,:)       !volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]  (nlevsoi)
   real(r8), pointer :: h2osno_old(:)         !snow mass for previous time step (kg/m2) (new)
   real(r8), pointer :: qg(:)                 !ground specific humidity [kg/kg]
   real(r8), pointer :: dqgdT(:)              !d(qg)/dT
   real(r8), pointer :: snowice(:)            !average snow ice lens
   real(r8), pointer :: snowliq(:)            !average snow liquid water
end type column_wstate_type

!----------------------------------------------------
! column carbon state variables structure
!----------------------------------------------------
type column_cstate_type
   type(pft_cstate_type):: pcs_a        !pft-level carbon state variables averaged to the column
   real(r8), pointer :: soilc(:)        !soil carbon (kg C /m**2)
end type column_cstate_type

!----------------------------------------------------
! column nitrogen state variables structure
!----------------------------------------------------
type column_nstate_type
   type(pft_nstate_type):: pns_a        !pft-level nitrogen state variables averaged to the column
end type column_nstate_type

!----------------------------------------------------
! column VOC state variables structure
!----------------------------------------------------
type column_vstate_type
   type(pft_vstate_type):: pvs_a        !pft-level VOC state variables averaged to the column
end type column_vstate_type

!----------------------------------------------------
! column DGVM state variables structure
!----------------------------------------------------
type column_dgvstate_type
   type(pft_dgvstate_type):: pdgvs_a
end type column_dgvstate_type

!----------------------------------------------------
! column dust state variables structure
!----------------------------------------------------
type column_dstate_type
   real(r8), pointer :: dummy_entry(:)
end type column_dstate_type

!----------------------------------------------------
! column energy flux variables structure
!----------------------------------------------------
type column_eflux_type
   type(pft_eflux_type):: pef_a !pft-level energy flux variables averaged to the column
   real(r8), pointer :: eflx_snomelt(:) !snow melt heat flux (W/m**2)
   real(r8), pointer :: eflx_impsoil(:) !implicit evaporation for soil temperature equation
end type column_eflux_type

!----------------------------------------------------
! column momentum flux variables structure
!----------------------------------------------------
type column_mflux_type
   type(pft_mflux_type)::  pmf_a    !pft-level momentum flux variables averaged to the column
end type column_mflux_type

!----------------------------------------------------
! column water flux variables structure
!----------------------------------------------------
type column_wflux_type
   type(pft_wflux_type):: pwf_a !pft-level water flux variables averaged to the column
   real(r8), pointer :: qflx_infl(:)    !infiltration (mm H2O /s)
   real(r8), pointer :: qflx_surf(:)    !surface runoff (mm H2O /s)
   real(r8), pointer :: qflx_drain(:)   !sub-surface runoff (mm H2O /s)
   real(r8), pointer :: qflx_top_soil(:)!net water input into soil from top (mm/s)
   real(r8), pointer :: qflx_snomelt(:) !snow melt (mm H2O /s)
   real(r8), pointer :: qflx_qrgwl(:)   !qflx_surf at glaciers, wetlands, lakes
   real(r8), pointer :: qmelt(:)        !snow melt [mm/s]
end type column_wflux_type

!----------------------------------------------------
! column carbon flux variables structure
!----------------------------------------------------
type column_cflux_type
   type(pft_cflux_type):: pcf_a         !pft-level carbon flux variables averaged to the column
end type column_cflux_type

!----------------------------------------------------
! column nitrogen flux variables structure
!----------------------------------------------------
type column_nflux_type
   type(pft_nflux_type):: pnf_a         !pft-level nitrogen flux variables averaged to the column
end type column_nflux_type

!----------------------------------------------------
! column VOC flux variables structure
!----------------------------------------------------
type column_vflux_type
   type(pft_vflux_type):: pvf_a         !pft-level VOC flux variables averaged to the column
end type column_vflux_type

!----------------------------------------------------
! column dust flux variables structure
!----------------------------------------------------
type column_dflux_type
   type(pft_dflux_type):: pdf_a         !pft-level dust flux variables averaged to the column
end type column_dflux_type
!----------------------------------------------------
! End definition of structures defined at the column_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the landunit_type level
!----------------------------------------------------
! landunit physical state variables structure
! note - landunit type can be vegetated (includes bare soil), deep lake,
! shallow lake, wetland, glacier or urban
!----------------------------------------------------
type landunit_pstate_type
   type(column_pstate_type):: cps_a             !column-level physical state variables averaged to landunit
end type landunit_pstate_type

!----------------------------------------------------
! landunit energy state variables structure
!----------------------------------------------------
type landunit_estate_type
   type(column_estate_type):: ces_a          !column-level energy state variables averaged to landunit
end type landunit_estate_type

!----------------------------------------------------
! landunit water state variables structure
!----------------------------------------------------
type landunit_wstate_type
   type(column_wstate_type):: cws_a           !column-level water state variables averaged to landunit
end type landunit_wstate_type

!----------------------------------------------------
! landunit carbon state variables structure
!----------------------------------------------------
type landunit_cstate_type
   type(column_cstate_type):: ccs_a           !column-level carbon state variables averaged to landunit
end type landunit_cstate_type

!----------------------------------------------------
! landunit nitrogen state variables structure
!----------------------------------------------------
type landunit_nstate_type
   type(column_nstate_type):: cns_a           !column-level nitrogen state variables averaged to landunit
end type landunit_nstate_type

!----------------------------------------------------
! landunit VOC state variables structure
!----------------------------------------------------
type landunit_vstate_type
   type(column_vstate_type):: cvs_a           !column-level VOC state variables averaged to landunit
end type landunit_vstate_type

!----------------------------------------------------
! landunit DGVM state variables structure
!----------------------------------------------------
type landunit_dgvstate_type
   real(r8):: dummy_entry
end type landunit_dgvstate_type

!----------------------------------------------------
! landunit dust state variables structure
!----------------------------------------------------
type landunit_dstate_type
   type(column_dstate_type):: cds_a             !column-level dust state variables averaged to landunit
end type landunit_dstate_type

!----------------------------------------------------
! landunit energy flux variables structure
!----------------------------------------------------
type landunit_eflux_type
   type(column_eflux_type):: cef_a              !column-level energy flux variables averaged to landunit
end type landunit_eflux_type

!----------------------------------------------------
! landunit momentum flux variables structure
!----------------------------------------------------
type landunit_mflux_type
   type(pft_mflux_type):: pmf_a                 !pft-level momentum flux variables averaged to landunit
end type landunit_mflux_type

!----------------------------------------------------
! landunit water flux variables structure
!----------------------------------------------------
type landunit_wflux_type
   type(column_wflux_type):: cwf_a              !column-level water flux variables averaged to landunit
end type landunit_wflux_type

!----------------------------------------------------
! landunit carbon flux variables structure
!----------------------------------------------------
type landunit_cflux_type
   type(column_cflux_type):: ccf_a               !column-level carbon flux variables averaged to landunit
end type landunit_cflux_type

!----------------------------------------------------
! landunit nitrogen flux variables structure
!----------------------------------------------------
type landunit_nflux_type
   type(column_nflux_type):: cnf_a              !column-level nitrogen flux variables averaged to landunit
end type landunit_nflux_type

!----------------------------------------------------
! landunit VOC flux variables structure
!----------------------------------------------------
type landunit_vflux_type
   type(pft_vflux_type):: pvf_a                  !pft-level VOC flux variables averaged to landunit
end type landunit_vflux_type

!----------------------------------------------------
! landunit dust flux variables structure
!----------------------------------------------------
type landunit_dflux_type
   type(pft_dflux_type):: pdf_a                  !pft-level dust flux variables averaged to landunit
end type landunit_dflux_type
!----------------------------------------------------
! End definition of structures defined at the landunit_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the gridcell_type level
!----------------------------------------------------
! gridcell physical state variables structure
!----------------------------------------------------
type gridcell_pstate_type
   type(column_pstate_type):: cps_a   !column-level physical state variables averaged to gridcell
   real(r8), pointer :: wtfact(:)      !Fraction of model area with high water table
end type gridcell_pstate_type

!----------------------------------------------------
! atmosphere -> land state variables structure
!----------------------------------------------------
type atm2lnd_state_type
#if (defined OFFLINE)
   real(r8), pointer :: flfall(:)       !fraction of liquid water within falling precipitation
#endif
   real(r8), pointer :: forc_t(:)       !atmospheric temperature (Kelvin)
   real(r8), pointer :: forc_u(:)       !atmospheric wind speed in east direction (m/s)
   real(r8), pointer :: forc_v(:)       !atmospheric wind speed in north direction (m/s)
   real(r8), pointer :: forc_wind(:)    !atmospheric wind speed
   real(r8), pointer :: forc_q(:)       !atmospheric specific humidity (kg/kg)
   real(r8), pointer :: forc_hgt(:)     !atmospheric reference height (m)
   real(r8), pointer :: forc_hgt_u(:)   !observational height of wind [m] (new)
   real(r8), pointer :: forc_hgt_t(:)   !observational height of temperature [m] (new)
   real(r8), pointer :: forc_hgt_q(:)   !observational height of humidity [m] (new)
   real(r8), pointer :: forc_pbot(:)    !atmospheric pressure (Pa)
   real(r8), pointer :: forc_th(:)      !atmospheric potential temperature (Kelvin)
   real(r8), pointer :: forc_vp(:)      !atmospheric vapor pressure (Pa)
   real(r8), pointer :: forc_rho(:)     !density (kg/m**3)
   real(r8), pointer :: forc_co2(:)     !atmospheric CO2 concentration (Pa)
   real(r8), pointer :: forc_o2(:)      !atmospheric O2 concentration (Pa)
   real(r8), pointer :: forc_psrf(:)    !surface pressure (Pa)
end type atm2lnd_state_type

!----------------------------------------------------
! land -> atmosphere state variables structure
!----------------------------------------------------
type lnd2atm_state_type
   real(r8), pointer :: t_rad(:)        !radiative temperature (Kelvin)
   real(r8), pointer :: t_ref2m(:)      !2 m height surface air temperature (Kelvin)
   real(r8), pointer :: q_ref2m(:)      !2 m height surface specific humidity (kg/kg)
   real(r8), pointer :: h2osno(:)       !snow water (mm H2O)
   real(r8), pointer :: albd(:,:)       !(numrad) surface albedo (direct)
   real(r8), pointer :: albi(:,:)       !(numrad) surface albedo (diffuse)
end type lnd2atm_state_type

!----------------------------------------------------
! gridcell energy state variables structure
!----------------------------------------------------
type gridcell_estate_type
   type(column_estate_type):: ces_a     !column-level energy state variables averaged to gridcell
end type gridcell_estate_type

!----------------------------------------------------
! gridcell water state variables structure
!----------------------------------------------------
type gridcell_wstate_type
   type(column_wstate_type):: cws_a     !column-level water state variables averaged to gridcell
end type gridcell_wstate_type

!----------------------------------------------------
! gridcell carbon state variables structure
!----------------------------------------------------
type gridcell_cstate_type
   type(column_cstate_type):: ccs_a     !column-level carbon state variables averaged to gridcell
end type gridcell_cstate_type

!----------------------------------------------------
! gridcell nitrogen state variables structure
!----------------------------------------------------
type gridcell_nstate_type
   type(column_nstate_type):: cns_a     !column-level nitrogen state variables averaged to gridcell
end type gridcell_nstate_type

!----------------------------------------------------
! gridcell VOC state variables structure
!----------------------------------------------------
type gridcell_vstate_type
   type(column_vstate_type):: cvs_a     !column-level VOC state variables averaged to gridcell
end type gridcell_vstate_type

!----------------------------------------------------
! gridcell dust state variables structure
!----------------------------------------------------
type gridcell_dstate_type
   type(column_dstate_type):: cds_a     !column-level dust state variables averaged to gridcell
end type gridcell_dstate_type

!----------------------------------------------------
! gridcell DGVM state variables structure
!----------------------------------------------------
type gridcell_dgvstate_type
   real(r8), pointer :: afirefrac(:)   ! fraction of gridcell affected by fire
   real(r8), pointer :: acfluxfire(:)  ! C flux to atmosphere from biomass burning
   real(r8), pointer :: bmfm(:,:)      ! biomass (NPP) for each naturally-vegetated pft
   real(r8), pointer :: afmicr(:,:)    ! microbial respiration (Rh) for each naturally-vegetated pft
   real(r8), pointer :: begwater(:)
   real(r8), pointer :: endwater(:)
   real(r8), pointer :: begenergy(:)
   real(r8), pointer :: endenergy(:)
end type gridcell_dgvstate_type

!----------------------------------------------------
! atmosphere -> land flux variables structure
!----------------------------------------------------
type atm2lnd_flux_type
   real(r8), pointer :: forc_lwrad(:)           !downward infrared (longwave) radiation (W/m**2)
   real(r8), pointer :: forc_solad(:,:)         !direct beam radiation (vis=forc_sols , nir=forc_soll ) (numrad)
   real(r8), pointer :: forc_solai(:,:)         !diffuse radiation     (vis=forc_solsd, nir=forc_solld) (numrad)
   real(r8), pointer :: forc_solar(:)           !incident solar radiation
   real(r8), pointer :: forc_rain(:)            !rain rate [mm/s]
   real(r8), pointer :: forc_snow(:)            !snow rate [mm/s]
end type atm2lnd_flux_type

!----------------------------------------------------
! land -> atmosphere flux variables structure
!----------------------------------------------------
type lnd2atm_flux_type
   real(r8), pointer :: taux(:)                 !wind stress: e-w (kg/m/s**2)
   real(r8), pointer :: tauy(:)                 !wind stress: n-s (kg/m/s**2)
   real(r8), pointer :: eflx_lh_tot(:)          !total latent heat flux (W/m8*2)  [+ to atm]
   real(r8), pointer :: eflx_sh_tot(:)          !total sensible heat flux (W/m**2) [+ to atm]
   real(r8), pointer :: eflx_lwrad_out(:)       !emitted infrared (longwave) radiation (W/m**2)
   real(r8), pointer :: qflx_evap_tot(:)        !qflx_evap_soi + qflx_evap_veg + qflx_tran_veg
   real(r8), pointer :: fsa(:)                  !solar radiation absorbed (total) (W/m**2)
end type lnd2atm_flux_type

!----------------------------------------------------
! gridcell energy flux variables structure
!----------------------------------------------------
type gridcell_eflux_type
   type(column_eflux_type):: cef_a              !column-level energy flux variables averaged to gridcell
end type gridcell_eflux_type

!----------------------------------------------------
! gridcell momentum flux variables structure
!----------------------------------------------------
type gridcell_mflux_type
   type(pft_mflux_type):: pmf_a                 !pft-level momentum flux variables averaged to gridcell
end type gridcell_mflux_type

!----------------------------------------------------
! gridcell water flux variables structure
!----------------------------------------------------
type gridcell_wflux_type
!FTO
   real(r8), pointer :: qchan2(:)               !history file RTM river (channel) flow (m**3 H2O /s)
   real(r8), pointer :: qchocn2(:)              !history file RTM river (channel) flow into ocean (m**3/s)
!FTO
   type(column_wflux_type):: cwf_a              !column-level water flux variables averaged to gridcell
end type gridcell_wflux_type

!----------------------------------------------------
! gridcell carbon flux variables structure
!----------------------------------------------------
type gridcell_cflux_type
   type(column_cflux_type):: ccf_a              !column-level carbon flux variables averaged to gridcell
end type gridcell_cflux_type

!----------------------------------------------------
! gridcell nitrogen flux variables structure
!----------------------------------------------------
type gridcell_nflux_type
   type(column_nflux_type):: cnf_a              !column-level nitrogen flux variables averaged to gridcell
end type gridcell_nflux_type

!----------------------------------------------------
! gridcell VOC flux variables structure
!----------------------------------------------------
type gridcell_vflux_type
   type(pft_vflux_type):: pvf_a                 !pft-level VOC flux variables averaged to gridcell
end type gridcell_vflux_type

!----------------------------------------------------
! gridcell dust flux variables structure
!----------------------------------------------------
type gridcell_dflux_type
   type(pft_dflux_type):: pdf_a                 !pft-level dust flux variables averaged to gridcell
end type gridcell_dflux_type
!----------------------------------------------------
! End definition of structures defined at the gridcell_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of structures defined at the CLM level
!----------------------------------------------------
! CLM physical state variables structure
!----------------------------------------------------
type model_pstate_type
   type(column_pstate_type) :: cps_a    !column-level physical state variables globally averaged
end type model_pstate_type

!----------------------------------------------------
! CLM energy state variables structure
!----------------------------------------------------
type model_estate_type
   type(column_estate_type):: ces_a             !column-level energy state variables globally averaged
end type model_estate_type

!----------------------------------------------------
! CLM water state variables structure
!----------------------------------------------------
type model_wstate_type
   type(column_wstate_type):: cws_a             !column-level water state variables globally averaged
end type model_wstate_type

!----------------------------------------------------
! CLM carbon state variables structure
!----------------------------------------------------
type model_cstate_type
   type(column_cstate_type):: ccs_a             !column-level carbon state variables globally averaged
end type model_cstate_type

!----------------------------------------------------
! CLM nitrogen state variables structure
!----------------------------------------------------
type model_nstate_type
   type(column_nstate_type):: cns_a             !column-level nitrogen state variables globally averaged
end type model_nstate_type

!----------------------------------------------------
! CLM VOC state variables structure
!----------------------------------------------------
type model_vstate_type
   type(column_vstate_type):: cvs_a             !column-level VOC state variables globally averaged
end type model_vstate_type

!----------------------------------------------------
! CLM dust state variables structure
!----------------------------------------------------
type model_dstate_type
   type(column_dstate_type):: cds_a             !column-level dust state variables globally averaged
end type model_dstate_type

!----------------------------------------------------
! CLM energy flux variables structure
!----------------------------------------------------
type model_eflux_type
   type(column_eflux_type):: cef_a              !column-level energy flux variables globally averaged
end type model_eflux_type

!----------------------------------------------------
! CLM momentum flux variables structure
!----------------------------------------------------
type model_mflux_type
   type(pft_mflux_type):: pmf_a                 !pft-level momentum flux variables globally averaged
end type model_mflux_type

!----------------------------------------------------
! CLM water flux variables structure
!----------------------------------------------------
type model_wflux_type
   type(column_wflux_type):: cwf_a              !column-level water flux variables globally averaged
end type model_wflux_type

!----------------------------------------------------
! CLM carbon flux variables structure
!----------------------------------------------------
type model_cflux_type
   type(column_cflux_type):: ccf_a              !column-level carbon flux variables globally averaged
end type model_cflux_type

!----------------------------------------------------
! CLM nitrogen flux variables structure
!----------------------------------------------------
type model_nflux_type
   type(column_nflux_type):: cnf_a              !column-level nitrogen flux variables globally averaged
end type model_nflux_type

!----------------------------------------------------
! CLM VOC flux variables structure
!----------------------------------------------------
type model_vflux_type
   type(pft_vflux_type):: pvf_a                 !pft-level VOC flux variables globally averaged
end type model_vflux_type

!----------------------------------------------------
! CLM dust flux variables structure
!----------------------------------------------------
type model_dflux_type
   type(pft_dflux_type):: pdf_a                 !pft-level dust flux variables globally averaged
end type model_dflux_type

!----------------------------------------------------
! End definition of structures defined at the model_type level
!----------------------------------------------------
!*******************************************************************************


!*******************************************************************************
!----------------------------------------------------
! Begin definition of spatial scaling hierarchy
!----------------------------------------------------

!----------------------------------------------------
! define the pft structure
!----------------------------------------------------

type pft_type
   ! indices into higher levels in hierarchy
   integer, pointer :: column(:)        !index into column level quantities
   integer, pointer :: landunit(:)      !index into landunit level quantities
   integer, pointer :: gridcell(:)      !index into gridcell level quantities

   ! topological mapping functionality
   integer , pointer :: itype(:)        !pft vegetation
   real(r8), pointer :: area(:)         !total land area for this pft (km^2)
   real(r8), pointer :: wtcol(:)        !weight (relative to column) for this pft (0-1)
   real(r8), pointer :: wtlunit(:)      !weight (relative to landunit) for this pft (0-1)
   real(r8), pointer :: wtgcell(:)      !weight (relative to gridcell) for this pft (0-1)
   integer , pointer :: ixy(:)          !xy lon index
   integer , pointer :: jxy(:)          !xy lat index
   integer , pointer :: mxy(:)          !m index for laixy(i,j,m),etc.
   integer , pointer :: snindex(:)      !corresponding pft index in s->n and east->west order
   real(r8), pointer :: latdeg(:)       !latitude (degrees)
   real(r8), pointer :: londeg(:)       !longitude (degrees)

   ! conservation check structures for the pft level
   type(energy_balance_type)   :: pebal !energy balance structure
   type(water_balance_type)    :: pwbal !water balance structure
   type(carbon_balance_type)   :: pcbal !carbon balance structure
   type(nitrogen_balance_type) :: pnbal !nitrogen balance structure

   ! DGVM state variables
   type(pft_dgvstate_type) :: pdgvs           !pft DGVM state variables

   ! state variables defined at the pft level
   type(pft_pstate_type) :: pps         !physical state variables
   type(pft_estate_type) :: pes         !pft energy state
   type(pft_wstate_type) :: pws         !pft water state
   type(pft_cstate_type) :: pcs         !pft carbon state
   type(pft_nstate_type) :: pns         !pft nitrogen state
   type(pft_vstate_type) :: pvs         !pft VOC state

   ! flux variables defined at the pft level
   type(pft_eflux_type)  :: pef         !pft energy flux
   type(pft_mflux_type)  :: pmf         !pft momentum flux
   type(pft_wflux_type)  :: pwf         !pft water flux
   type(pft_cflux_type)  :: pcf         !pft carbon flux
   type(pft_nflux_type)  :: pnf         !pft nitrogen flux
   type(pft_vflux_type)  :: pvf         !pft VOC flux
   type(pft_dflux_type)  :: pdf         !pft dust flux
end type pft_type

!----------------------------------------------------
! define the column structure
!----------------------------------------------------

type column_type
   ! lower levels in hierarchy
   type(pft_type)   :: p                !plant functional type (pft) data structure
   integer , pointer :: pfti(:)         !beginning pft index for each column
   integer , pointer :: pftf(:)         !ending pft index for each column
   integer , pointer :: npfts(:)        !number of pfts for each column

   ! higher level in hierarchy
   integer , pointer :: landunit(:)     !index into landunit level quantities
   integer , pointer :: gridcell(:)     !index into gridcell level quantities

   ! topological mapping functionality
   integer , pointer :: itype(:)        !column type
   real(r8), pointer :: area(:)         !total land area for this column (km^2)
   real(r8), pointer :: wtgcell(:)      !weight (relative to gridcell) for this column (0-1)
   real(r8), pointer :: wtlunit(:)      !weight (relative to landunit) for this column (0-1)
   integer , pointer :: ixy(:)          !xy lon index
   integer , pointer :: jxy(:)          !xy lat index
   integer , pointer :: snindex(:)      !corresponding column index in s->n and east->west order
   real(r8), pointer :: latdeg(:)       !latitude (degrees)
   real(r8), pointer :: londeg(:)       !longitude (degrees)

   ! conservation check structures for the column level
   type(energy_balance_type)   :: cebal !energy balance structure
   type(water_balance_type)    :: cwbal !water balance structure
   type(carbon_balance_type)   :: ccbal !carbon balance structure
   type(nitrogen_balance_type) :: cnbal !nitrogen balance structure

   ! state variables defined at the column level
   type(column_pstate_type) :: cps      !column physical state variables
   type(column_estate_type) :: ces      !column energy state
   type(column_wstate_type) :: cws      !column water state
   type(column_cstate_type) :: ccs      !column carbon state
   type(column_nstate_type) :: cns      !column nitrogen state
   type(column_vstate_type) :: cvs      !column VOC state
   type(column_dstate_type) :: cds      !column dust state

   ! flux variables defined at the column level
   type(column_eflux_type) :: cef       !column energy flux
   type(column_mflux_type) :: cmf       !column momentum flux
   type(column_wflux_type) :: cwf       !column water flux
   type(column_cflux_type) :: ccf       !column carbon flux
   type(column_nflux_type) :: cnf       !column nitrogen flux
   type(column_vflux_type) :: cvf       !column VOC flux
   type(column_dflux_type) :: cdf       !column dust flux

   ! dgvm variables defined at the column level
   type (column_dgvstate_type) :: cdgvs !column DGVM structure
end type column_type

!----------------------------------------------------
! define the geomorphological land unit structure
!----------------------------------------------------

type landunit_type
   ! lower levels in hierarchy
   type(column_type) :: c                !column data structure (soil/snow/canopy columns)
   integer, pointer :: coli(:)           !beginning column index for each landunit
   integer, pointer :: colf(:)           !ending column index for each landunit
   integer, pointer :: ncolumns(:)       !number of columns for each landunit
   integer, pointer :: pfti(:)           !beginning pft index for each landunit
   integer, pointer :: pftf(:)           !ending pft index for each landunit
   integer, pointer :: npfts(:)          !number of pfts for each landunit

   ! higher level in hierarchy
   integer, pointer :: gridcell(:)       !index into gridcell level quantities

   ! topological mapping functionality
   integer , pointer :: itype(:)         !landunit type
   real(r8), pointer :: area(:)          !total land area for this landunit (km^2)
   real(r8), pointer :: wtgcell(:)       !weight (relative to gridcell) for this landunit (0-1)
   integer , pointer :: ixy(:)           !xy lon index
   integer , pointer :: jxy(:)           !xy lat index
   integer , pointer :: snindex(:)       !corresponding landunit index in s->n and east->west order
   real(r8), pointer :: latdeg(:)        !latitude (degrees)
   real(r8), pointer :: londeg(:)        !longitude (degrees)
   logical , pointer :: ifspecial(:)     !BOOL: true=>landunit is not vegetated
   logical , pointer :: lakpoi(:)        !BOOL: true=>lake point

   ! conservation check structures for the landunit level
   type(energy_balance_type)   :: lebal  !energy balance structure
   type(water_balance_type)    :: lwbal  !water balance structure
   type(carbon_balance_type)   :: lcbal  !carbon balance structure
   type(nitrogen_balance_type) :: lnbal  !nitrogen balance structure

   ! state variables defined at the land unit level
   type(landunit_pstate_type) :: lps     !land unit physical state variables
   type(landunit_estate_type) :: les     !average of energy states over all columns
   type(landunit_wstate_type) :: lws     !average of water states over all columns
   type(landunit_cstate_type) :: lcs     !average of carbon states over all columns
   type(landunit_nstate_type) :: lns     !average of nitrogen states over all columns
   type(landunit_vstate_type) :: lvs     !average of VOC states over all columns
   type(landunit_dstate_type) :: lds     !average of dust states over all columns

   ! flux variables defined at the landunit level
   type(landunit_eflux_type) :: lef      !average of energy fluxes over all columns
   type(landunit_mflux_type) :: lmf      !average of momentum fluxes over all columns
   type(landunit_wflux_type) :: lwf      !average of water fluxes over all columns
   type(landunit_cflux_type) :: lcf      !average of carbon fluxes over all columns
   type(landunit_nflux_type) :: lnf      !average of nitrogen fluxes over all columns
   type(landunit_vflux_type) :: lvf      !average of VOC fluxes over all columns
   type(landunit_dflux_type) :: ldf      !average of dust fluxes over all columns
end type landunit_type

!----------------------------------------------------
! define the gridcell structure
!----------------------------------------------------

type gridcell_type
   ! lower level in hierarchy
   type(landunit_type) :: l               !geomorphological landunits
   integer, pointer :: luni(:)            !beginning landunit index for each gridcell
   integer, pointer :: lunf(:)            !ending landunit index for each gridcell
   integer, pointer :: nlandunits(:)      !number of landunit for each gridcell
   integer, pointer :: coli(:)            !beginning column index for each gridcell
   integer, pointer :: colf(:)            !ending column index for each gridcell
   integer, pointer :: ncolumns(:)        !number of columns for each gridcell
   integer, pointer :: pfti(:)            !beginning pft index for each gridcell
   integer, pointer :: pftf(:)            !ending pft index for each gridcell
   integer, pointer :: npfts(:)           !number of pfts for each gridcell

   ! topological mapping functionality
   integer , pointer :: itype(:)          !gridcell type
   real(r8), pointer :: area(:)           !total land area for this gridcell (km^2)
   real(r8), pointer :: wtglob(:)         !weight for this gridcell relative to global area (0-1)
   integer , pointer :: ixy(:)            !xy lon index
   integer , pointer :: jxy(:)            !xy lat index
   integer , pointer :: snindex(:)        !corresponding gridcell index in s->n and east->west order
   real(r8), pointer :: lat(:)            !latitude (radians)
   real(r8), pointer :: lon(:)            !longitude (radians)
   real(r8), pointer :: latdeg(:)         !latitude (degrees)
   real(r8), pointer :: londeg(:)         !longitude (degrees)
   real(r8), pointer :: landfrac(:)       !fractional land for this gridcell

   ! conservation check structures for the gridcell level
   type(energy_balance_type)   :: gebal  !energy balance structure
   type(water_balance_type)    :: gwbal  !water balance structure
   type(carbon_balance_type)   :: gcbal  !carbon balance structure
   type(nitrogen_balance_type) :: gnbal  !nitrogen balance structure

   ! dgvm variables defined at the gridcell level
   type(gridcell_dgvstate_type):: gdgvs !gridcell DGVM structure

   ! state variables defined at the gridcell level
   type(gridcell_pstate_type) :: gps     !gridcell physical state variables
   type(gridcell_estate_type) :: ges     !average of energy states over all landunits
   type(gridcell_wstate_type) :: gws     !average of water states over all landunits
   type(gridcell_cstate_type) :: gcs     !average of carbon states over all landunits
   type(gridcell_nstate_type) :: gns     !average of nitrogen states over all landunits
   type(gridcell_vstate_type) :: gvs     !average of VOC states over all landunits
   type(gridcell_dstate_type) :: gds     !average of dust states over all landunits
   type(atm2lnd_state_type)   :: a2ls    !atmospheric state variables required by the land
   type(lnd2atm_state_type)   :: l2as    !land state variables required by the atmosphere

   ! flux variables defined at the gridcell level
   type(gridcell_eflux_type) :: gef       !average of energy fluxes over all landunits
   type(gridcell_wflux_type) :: gwf       !average of water fluxes over all landunits
   type(gridcell_cflux_type) :: gcf       !average of carbon fluxes over all landunits
   type(gridcell_nflux_type) :: gnf       !average of nitrogen fluxes over all landunits
   type(gridcell_vflux_type) :: gvf       !average of VOC fluxes over all landunits
   type(gridcell_dflux_type) :: gdf       !average of dust fluxes over all landunits
   type(atm2lnd_flux_type)   :: a2lf      !atmospheric flux variables required by the land
   type(lnd2atm_flux_type)   :: l2af      !land flux variables required by the atmosphere
end type gridcell_type

!----------------------------------------------------
! define the top-level (model) structure
!----------------------------------------------------

type model_type
   ! lower level in hierarch
   type(gridcell_type) :: g               !gridicell data structure
   integer  :: ngridcells                 !number of gridcells allocated for this process
   real(r8) :: area                       !total land area for all gridcells (km^2)

   ! conservation check structures for the clm (global) level
   type(energy_balance_type)   :: mebal  !energy balance structure
   type(water_balance_type)    :: mwbal  !water balance structure
   type(carbon_balance_type)   :: mcbal  !carbon balnace structure
   type(nitrogen_balance_type) :: mnbal  !nitrogen balance structure

   ! globally average state variables
   type(model_pstate_type) ::  mps       !clm physical state variables
   type(model_estate_type) ::  mes       !average of energy states over all gridcells
   type(model_wstate_type) ::  mws       !average of water states over all gridcells
   type(model_cstate_type) ::  mcs       !average of carbon states over all gridcells
   type(model_nstate_type) ::  mns       !average of nitrogen states over all gridcells
   type(model_vstate_type) ::  mvs       !average of VOC states over all gridcells
   type(model_dstate_type) ::  mds       !average of dust states over all gridcells

   ! globally averaged flux variables
   type(model_eflux_type) ::   mef       !average of energy fluxes over all gridcells
   type(model_wflux_type) ::   mwf       !average of water fluxes over all gridcells
   type(model_cflux_type) ::   mcf       !average of carbon fluxes over all gridcells
   type(model_nflux_type) ::   mnf       !average of nitrogen fluxes over all gridcells
   type(model_vflux_type) ::   mvf       !average of VOC fluxes over all gridcells
   type(model_dflux_type) ::   mdf       !average of dust fluxes over all gridcells
end type model_type
!----------------------------------------------------
! End definition of spatial scaling hierarchy
!----------------------------------------------------
!*******************************************************************************

!*******************************************************************************
!----------------------------------------------------
! Declare single instance of clmtype
!----------------------------------------------------
type(model_type), target, save :: clm3

!----------------------------------------------------
! Declare single instance of array of ecophysiological constant types
!----------------------------------------------------
type(pft_epc_type), target, save :: pftcon

!----------------------------------------------------
! Declare single instance of array of dgvm ecophysiological constant types
!----------------------------------------------------
type(pft_dgvepc_type), target, save :: dgv_pftcon

character(len=8), parameter :: nameg  = 'gridcell'  ! name of gridcells
character(len=8), parameter :: namel  = 'landunit'  ! name of landunits
character(len=8), parameter :: namec  = 'column'    ! name of columns
character(len=8), parameter :: namep  = 'pft'       ! name of pfts
character(len=8), parameter :: ocnrof = 'ocnrof'    ! name of river routing ocean runoff
character(len=8), parameter :: lndrof = 'lndrof'    ! name of river routing land channel runoff
!
!EOP
!-----------------------------------------------------------------------
end module clmtype
