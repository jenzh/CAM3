#include <misc.h>
#include <preproc.h>

module histFldsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: histFldsMod
!
! !DESCRIPTION:
! Module containing initialization of clm history fields and files
! This is the module that the user must modify in order to add new
! history fields or modify defaults associated with existing history
! fields.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
  public initHistFlds ! Build master field list of all possible history
                      ! file fields
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 03/2003
!
!EOP
!------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: initHistFlds
!
! !INTERFACE:
  subroutine initHistFlds()
!
! !DESCRIPTION:
! Build master field list of all possible fields in a history file.
! Each field has associated with it a ``long\_name'' netcdf attribute that
! describes what the field is, and a ``units'' attribute. A subroutine is
! called to add each field to the masterlist.
!
! !USES:
    use clmtype
    use clm_varcon , only : spval
    use clm_varctl , only : nsrest
#if (defined RTM)
    use RunoffMod  , only : runoff
#endif
    use histFileMod, only : add_subscript, add_fld1d, add_fld2d, &
                            masterlist_printflds, htapes_build
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Mariana Vertenstein: Created 03/2003
! Mariana Vertenstein: Updated interface to create history fields 10/2003
!
!EOP
!-----------------------------------------------------------------------

    ! Determine what subscripts to add
    ! (uncomment the following call and modify it appropriately)

    ! call add_subscript(subname='subscript_name', subdim=subscript_dim)

    ! NOTE: make a field not appear on the primary history tape by default -
    ! add the keyword to default='inactive' to the call to addfld_1d or addfld_2d

    ! Snow properties
    ! These will be vertically averaged over the snow profile

    call add_fld1d (fname='SNOWDP',  units='m',  &
         avgflag='A', long_name='snow height', &
         ptr_col=clm3%g%l%c%cps%snowdp)

    call add_fld1d (fname='SNOWAGE',  units='unitless',  &
         avgflag='A', long_name='snow age', &
         ptr_col=clm3%g%l%c%cps%snowage)

    call add_fld1d (fname='FSNO',  units='unitless',  &
         avgflag='A', long_name='fraction of ground covered by snow', &
         ptr_col=clm3%g%l%c%cps%frac_sno)

    ! Temperatures

    call add_fld1d (fname='TSA', units='K',  &
         avgflag='A', long_name='2m air temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m)

    call add_fld1d (fname='TREFMNAV', units='K',  &
         avgflag='A', long_name='daily minimum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_min)

    call add_fld1d (fname='TREFMXAV', units='K',  &
         avgflag='A', long_name='daily maximum of average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_ref2m_max)

    call add_fld1d (fname='TV', units='K',  &
         avgflag='A', long_name='vegetation temperature', &
         ptr_pft=clm3%g%l%c%p%pes%t_veg)

    call add_fld1d (fname='TG',  units='K',  &
         avgflag='A', long_name='ground temperature', &
         ptr_col=clm3%g%l%c%ces%t_grnd)

    call add_fld1d (fname='TSNOW',  units='K',  &
         avgflag='A', long_name='snow temperature', &
         ptr_col=clm3%g%l%c%ces%t_snow, set_lake=spval)

    call add_fld2d (fname='TSOI',  units='K', type2d='levsoi', &
         avgflag='A', long_name='soil temperature', &
         ptr_col=clm3%g%l%c%ces%t_soisno)

    call add_fld2d (fname='TLAKE',  units='K', type2d='levsoi', &
         avgflag='A', long_name='lake temperature', &
         ptr_col=clm3%g%l%c%ces%t_lake)

    ! Specific humidity

    call add_fld1d (fname='Q2M', units='kg/kg',  &
         avgflag='A', long_name='2m specific humidity', &
         ptr_pft=clm3%g%l%c%p%pes%q_ref2m)

    ! Surface radiation

    call add_fld1d (fname='SABV', units='watt/m^2',  &
         avgflag='A', long_name='solar rad absorbed by veg', &
         ptr_pft=clm3%g%l%c%p%pef%sabv)

    call add_fld1d (fname='SABG', units='watt/m^2',  &
         avgflag='A', long_name='solar rad absorbed by ground', &
         ptr_pft=clm3%g%l%c%p%pef%sabg)

    call add_fld1d (fname='FSDSVD', units='watt/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_d)

    call add_fld1d (fname='FSDSND', units='watt/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_nir_d)

    call add_fld1d (fname='FSDSVI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse vis incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_i)

    call add_fld1d (fname='FSDSNI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse nir incident solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_nir_i)

    call add_fld1d (fname='FSRVD', units='watt/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_vis_d)

    call add_fld1d (fname='FSRND', units='watt/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_nir_d)

    call add_fld1d (fname='FSRVI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse vis reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_vis_i)

    call add_fld1d (fname='FSRNI', units='watt/m^2',  &
         avgflag='A', long_name='diffuse nir reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_nir_i)

    call add_fld1d (fname='FSDSVDLN', units='watt/m^2',  &
         avgflag='A', long_name='direct vis incident solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_vis_d_ln)

    call add_fld1d (fname='FSDSNDLN', units='watt/m^2',  &
         avgflag='A', long_name='direct nir incident solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsds_nir_d_ln)

    call add_fld1d (fname='FSRVDLN', units='watt/m^2',  &
         avgflag='A', long_name='direct vis reflected solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_vis_d_ln)

    call add_fld1d (fname='FSRNDLN', units='watt/m^2',  &
         avgflag='A', long_name='direct nir reflected solar radiation at local noon', &
         ptr_pft=clm3%g%l%c%p%pef%fsr_nir_d_ln)

    call add_fld1d (fname='FSA', units='watt/m^2',  &
         avgflag='A', long_name='absorbed solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsa)

    call add_fld1d (fname='FSR', units='watt/m^2',  &
         avgflag='A', long_name='reflected solar radiation', &
         ptr_pft=clm3%g%l%c%p%pef%fsr)

    call add_fld1d (fname='FIRA', units='watt/m^2',  &
         avgflag='A', long_name='net infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_net)

    call add_fld1d (fname='FIRE', units='watt/m^2',  &
         avgflag='A', long_name='emitted infrared (longwave) radiation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lwrad_out)

    ! Surface energy fluxes

    call add_fld1d (fname='FCTR', units='watt/m^2',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_vegt, set_lake=0._r8)

    call add_fld1d (fname='FCEV', units='watt/m^2',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_vege, set_lake=0._r8)

    call add_fld1d (fname='FGEV', units='watt/m^2',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_lh_grnd)

    call add_fld1d (fname='FSH', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_tot)

    call add_fld1d (fname='FSH_V', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat from veg', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_veg, set_lake=0._r8)

    call add_fld1d (fname='FSH_G', units='watt/m^2',  &
         avgflag='A', long_name='sensible heat from ground', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_sh_grnd )

    call add_fld1d (fname='FGR', units='watt/m^2',  &
         avgflag='A', long_name='heat flux into soil', &
         ptr_pft=clm3%g%l%c%p%pef%eflx_soil_grnd)

    call add_fld1d (fname='FSM',  units='watt/m^2',  &
         avgflag='A', long_name='snow melt heat flux', &
         ptr_col=clm3%g%l%c%cef%eflx_snomelt)

    call add_fld1d (fname='TAUX', units='kg/m/s^2',  &
         avgflag='A', long_name='zonal surface stress', &
         ptr_pft=clm3%g%l%c%p%pmf%taux)

    call add_fld1d (fname='TAUY', units='kg/m/s^2',  &
         avgflag='A', long_name='meridional surface stress', &
         ptr_pft=clm3%g%l%c%p%pmf%tauy)

    ! Vegetation phenology

    call add_fld1d (fname='ELAI', units='m^2/m^2', &
          avgflag='A', long_name='exposed one-sided leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%elai)

    call add_fld1d (fname='ESAI', units='m^2/m^2', &
          avgflag='A', long_name='exposed one-sided stem area index', &
         ptr_pft=clm3%g%l%c%p%pps%esai)

    ! Canopy physiology

    call add_fld1d (fname='RSSUN', units='s/m',  &
         avgflag='M', long_name='sunlit leaf stomatal resistance', &
         ptr_pft=clm3%g%l%c%p%pps%rssun, set_lake=spval)

    call add_fld1d (fname='RSSHA', units='s/m',  &
         avgflag='M', long_name='shaded leaf stomatal resistance', &
         ptr_pft=clm3%g%l%c%p%pps%rssha, set_lake=spval)

    call add_fld1d (fname='BTRAN', units='unitless',  &
         avgflag='A', long_name='transpiration beta factor', &
         ptr_pft=clm3%g%l%c%p%pps%btran, set_lake=spval)

    call add_fld1d (fname='FPSN', units='umol/m2s',  &
         avgflag='A', long_name='photosynthesis', &
         ptr_pft=clm3%g%l%c%p%pcf%fpsn, set_lake=0._r8)

#if (defined VOC)
    call add_fld1d (fname='VOCFLXT', units='uG/M2/H',  &
         avgflag='A', long_name='total VOC flux into atmosphere', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_tot, set_lake=0._r8)

    call add_fld1d (fname='ISOPRENE', units='uG/M2/H',  &
         avgflag='A', long_name='isoprene flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_1, set_lake=0._r8)

    call add_fld1d (fname='MONOTERP', units='uG/M2/H',  &
         avgflag='A', long_name='monoterpene flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_2, set_lake=0._r8)

    call add_fld1d (fname='OVOC', units='uG/M2/H',  &
         avgflag='A', long_name='other VOC flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_3, set_lake=0._r8)

    call add_fld1d (fname='ORVOC', units='uG/M2/H',  &
         avgflag='A', long_name='other reactive VOC flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_4, set_lake=0._r8)

    call add_fld1d (fname='BIOGENCO', units='uG/M2/H',  &
         avgflag='A', long_name='biogenic CO flux', &
         ptr_pft=clm3%g%l%c%p%pvf%vocflx_5, set_lake=0._r8)
#endif

    ! Hydrology

    call add_fld1d (fname='H2OSNO',  units='mm',  &
         avgflag='A', long_name='snow depth (liquid water)', &
         ptr_col=clm3%g%l%c%cws%h2osno)

    call add_fld1d (fname='H2OCAN', units='mm',  &
         avgflag='A', long_name='intercepted water', &
         ptr_pft=clm3%g%l%c%p%pws%h2ocan, set_lake=0._r8)

    call add_fld2d (fname='H2OSOI',  units='mm3/mm3', type2d='levsoi', &
         avgflag='A', long_name='volumetric soil water', &
         ptr_col=clm3%g%l%c%cws%h2osoi_vol)

    call add_fld2d (fname='SOILLIQ',  units='kg/m2', type2d='levsoi', &
         avgflag='A', long_name='soil liquid water', &
         ptr_col=clm3%g%l%c%cws%h2osoi_liq)

    call add_fld2d (fname='SOILICE',  units='kg/m2', type2d='levsoi', &
         avgflag='A', long_name='soil ice', &
         ptr_col=clm3%g%l%c%cws%h2osoi_ice)

    call add_fld1d (fname='SNOWLIQ',  units='kg/m2',  &
         avgflag='A', long_name='snow liquid water', &
         ptr_col=clm3%g%l%c%cws%snowliq)

    call add_fld1d (fname='SNOWICE',  units='kg/m2',  &
         avgflag='A', long_name='snow ice', &
         ptr_col=clm3%g%l%c%cws%snowice)

    call add_fld1d (fname='QINFL',  units='mm/s',  &
         avgflag='A', long_name='infiltration', &
         ptr_col=clm3%g%l%c%cwf%qflx_infl)

    call add_fld1d (fname='QOVER',  units='mm/s',  &
         avgflag='A', long_name='surface runoff', &
         ptr_col=clm3%g%l%c%cwf%qflx_surf)

    call add_fld1d (fname='QRGWL',  units='mm/s',  &
         avgflag='A', long_name='surface runoff at glaciers, wetlands, lakes', &
         ptr_col=clm3%g%l%c%cwf%qflx_qrgwl)

    call add_fld1d (fname='QDRAI',  units='mm/s',  &
         avgflag='A', long_name='sub-surface drainage', &
         ptr_col=clm3%g%l%c%cwf%qflx_drain)

    call add_fld1d (fname='QINTR', units='mm/s',  &
         avgflag='A', long_name='interception', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_prec_intr, set_lake=0._r8)

    call add_fld1d (fname='QDRIP', units='mm/s',  &
         avgflag='A', long_name='throughfall', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_prec_grnd)

    call add_fld1d (fname='QMELT',  units='mm/s',  &
         avgflag='A', long_name='snow melt', &
         ptr_col=clm3%g%l%c%cwf%qflx_snomelt)

    call add_fld1d (fname='QSOIL', units='mm/s',  &
         avgflag='A', long_name='ground evaporation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_soi)

    call add_fld1d (fname='QVEGE', units='mm/s',  &
         avgflag='A', long_name='canopy evaporation', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_evap_can, set_lake=0._r8)

    call add_fld1d (fname='QVEGT', units='mm/s',  &
         avgflag='A', long_name='canopy transpiration', &
         ptr_pft=clm3%g%l%c%p%pwf%qflx_tran_veg, set_lake=0._r8)

#if (defined RTM)
    ! RTM River Routing

    call add_fld1d (fname='QCHANR', units='m3/s',  typexy='rof', &
         avgflag='A', long_name='RTM river flow', &
         ptr_roflnd=runoff%lnd)

    call add_fld1d (fname='QCHOCNR', units='m3/s', typexy='rof', &
         avgflag='A', long_name='RTM river discharge into ocean', &
         ptr_rofocn=runoff%ocn)
#endif

    ! Water and energy balance checks

    call add_fld1d (fname='ERRSOI',  units='watt/m^2',  &
         avgflag='A', long_name='soil/lake energy conservation error', &
         ptr_col=clm3%g%l%c%cebal%errsoi)

    call add_fld1d (fname='ERRSEB',  units='watt/m^2',  &
         avgflag='A', long_name='surface energy conservation error', &
         ptr_pft=clm3%g%l%c%p%pebal%errseb)

    call add_fld1d (fname='ERRSOL',  units='watt/m^2',  &
         avgflag='A', long_name='solar radiation conservation error', &
         ptr_pft=clm3%g%l%c%p%pebal%errsol)

    call add_fld1d (fname='ERRH2O', units='mm',  &
         avgflag='A', long_name='total water conservation error', &
         ptr_col=clm3%g%l%c%cwbal%errh2o)

    ! Atmospheric forcing

    call add_fld1d (fname='RAIN', units='mm/s',  &
         avgflag='A', long_name='atmospheric rain', &
         ptr_gcell=clm3%g%a2lf%forc_rain)

    call add_fld1d (fname='SNOW', units='mm/s',  &
         avgflag='A', long_name='atmospheric snow', &
         ptr_gcell=clm3%g%a2lf%forc_snow)

    call add_fld1d (fname='TBOT', units='K',  &
         avgflag='A', long_name='atmospheric air temperature', &
         ptr_gcell=clm3%g%a2ls%forc_t)

    call add_fld1d (fname='THBOT', units='K',  &
         avgflag='A', long_name='atmospheric air potential temperature', &
         ptr_gcell=clm3%g%a2ls%forc_th)

    call add_fld1d (fname='WIND', units='m/s',  &
         avgflag='A', long_name='atmospheric wind velocity magnitude', &
         ptr_gcell=clm3%g%a2ls%forc_wind)

    call add_fld1d (fname='QBOT', units='kg/kg',  &
         avgflag='A', long_name='atmospheric specific humidity', &
         ptr_gcell=clm3%g%a2ls%forc_q)

    call add_fld1d (fname='ZBOT', units='m',  &
         avgflag='A', long_name='atmospheric reference height', &
         ptr_gcell=clm3%g%a2ls%forc_hgt)

    call add_fld1d (fname='FLDS', units='watt/m^2',  &
         avgflag='A', long_name='atmospheric longwave radiation', &
         ptr_gcell=clm3%g%a2lf%forc_lwrad)

    call add_fld1d (fname='FSDS', units='watt/m^2',  &
         avgflag='A', long_name='atmospheric incident solar radiation', &
         ptr_gcell=clm3%g%a2lf%forc_solar)

#if (defined DGVM)
    ! History output of accumulation variables

    call add_fld1d (fname='FMICR', units='umol/m2s',  &
         avgflag='A', long_name='microbial respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%fmicr, set_lake=0._r8)

    call add_fld1d (fname='FRMS', units='umol/m2s',  &
         avgflag='A', long_name='stem maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%frms, set_lake=0._r8)

    call add_fld1d (fname='FRMR', units='umol/m2s',  &
         avgflag='A', long_name='root maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%frmr, set_lake=0._r8)

    call add_fld1d (fname='FRMF', units='umol/m2s',  &
         avgflag='A', long_name='foliage maintenance respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%frmf, set_lake=0._r8)

    call add_fld1d (fname='FRG', units='umol/m2s',  &
         avgflag='A', long_name='growth respiration', &
         ptr_pft=clm3%g%l%c%p%pcf%frg, set_lake=0._r8)

    call add_fld1d (fname='FCO2', units='umol/m2s',  &
         avgflag='A', long_name='net CO2 flux', &
         ptr_pft=clm3%g%l%c%p%pcf%fco2, set_lake=0._r8)

    call add_fld1d (fname='DMI', units='umol/m2s',  &
         avgflag='A', long_name='net primary production', &
         ptr_pft=clm3%g%l%c%p%pcf%dmi, set_lake=0._r8)

    call add_fld1d (fname='HTOP', units='m',  &
         avgflag='A', long_name='height of top of canopy', &
         ptr_pft=clm3%g%l%c%p%pps%htop)

    call add_fld1d (fname='HBOT', units='m',  &
         avgflag='A', long_name='height of bottom of canopy', &
         ptr_pft=clm3%g%l%c%p%pps%hbot)

    call add_fld1d (fname='TLAI', units='m^2/m^2', &
          avgflag='A', long_name='total one-sided leaf area index', &
         ptr_pft=clm3%g%l%c%p%pps%tlai)

    call add_fld1d (fname='TSAI', units='m^2/m^2', &
          avgflag='A', long_name='total one-sided stem area index', &
         ptr_pft=clm3%g%l%c%p%pps%tsai)

    call add_fld1d (fname='TDA', units='K',  &
         avgflag='A', long_name='daily average 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pdgvs%t_mo)

    call add_fld1d (fname='T10', units='K',  &
         avgflag='A', long_name='10-day running mean of 2-m temperature', &
         ptr_pft=clm3%g%l%c%p%pdgvs%t10)

    call add_fld1d (fname='AGDD0', units='K',  &
         avgflag='A', long_name='growing degree-days base 0C', &
         ptr_pft=clm3%g%l%c%p%pdgvs%agdd0)

    call add_fld1d (fname='AGDD5', units='K',  &
         avgflag='A', long_name='growing degree-days base 5C', &
         ptr_pft=clm3%g%l%c%p%pdgvs%agdd5)
#endif

    ! Print masterlist of history fields

    call masterlist_printflds ()

    ! Initialize active history fields if not a restart run
    ! If a restart run, then this information is obtained from the restart history file

    if (nsrest == 0 .or. nsrest == 3) call htapes_build ()

  end subroutine initHistFlds

end module histFldsMod
