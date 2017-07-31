#include <misc.h>
subroutine initindx
!----------------------------------------------------------------------- 
! 
! Purpose: Register constituents and physics buffer fields.
! 
! Author:    CSM Contact: M. Vertenstein, Aug. 1997
!            B.A. Boville, Oct 2001
!
! Modified for water isotope tracers
!  David Noone <dcn@colorado.edu> - Thu Jul  1 18:59:09 MDT 2004
!
! 
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use constituents, only: pcnst, ppcnst, cnst_add, advected, cnst_chk_dim, cnst_name
  use phys_buffer,  only: pbuf_init
  use chemistry,    only: trace_gas, chem_register
  use cldcond,      only: cldcond_register
  use physconst,    only: mwdry, cpair, mwh2o, cph2o
  use tracers, only: tracers_register
!  use constituents, only: dcconnam, sflxnam, hadvnam, vadvnam, fixcnam, 
  use constituents, only: dcconnam, sflxnam, tendnam, tottnam
  use check_energy, only: check_energy_register
  use aerosol_intr, only: aerosol_register_cnst
  use water_tracers, only: trace_water, wtrc_register

#if ( defined BFB_CAM_SCAM_IOP )
  use iop
#endif
  implicit none
!-----------------------------------------------------------------------
#include <comctl.h>
!---------------------------Local variables-----------------------------
!
  integer m            ! loop index
  integer mm           ! constituent index 
!-----------------------------------------------------------------------

! Initialize physics buffer
  call pbuf_init()

! Register water vapor.
! ***** N.B. ***** This must be the first call to cnst_add so that
!                  water vapor is constituent 1.
! (for clarity, set the mixtype to dry)
  call cnst_add('Q', advected, mwh2o, cph2o, 1.E-12_r8, mm, &
                longname='Specific humidity', readiv=.true.)
!
! Register cloud water
  call cldcond_register()
!
! Register water tracers (immediately after prognostic water)
  if (trace_water) then
     call wtrc_register()
  endif
!
! Register chemical constituents
! CO2 CH4 CFC11 CFC12
  if (trace_gas) then
     call chem_register()
  endif
!
! register aerosols
! SO2 SO4 DMS H2O2
  call aerosol_register_cnst()

! Register advected test tracers and determine starting index
  call tracers_register()

!
! All tracers registered, check that the dimensions are correct
  call cnst_chk_dim()
!
! Set default names for non-water advected and non-advected tracers
! Set names of advected and non-advected tracer diagnostics
!
  do m=1,ppcnst
     dcconnam(m) = 'DC'//cnst_name(m)
     sflxnam(m)  = 'SF'//cnst_name(m)
  end do

  do m=1,pcnst
!     hadvnam(m)  = 'HA'//cnst_name(m)
!     vadvnam(m)  = 'VA'//cnst_name(m)
!     fixcnam(m)  = 'DF'//cnst_name(m)
     tendnam(m)  = 'TE'//cnst_name(m)
     tottnam(m)  = 'TA'//cnst_name(m)
  end do

#if ( defined BFB_CAM_SCAM_IOP )
  do m=1,pcnst
     alphanam(m) = 'AFIX'//cnst_name(m)
     alphanam(m)=to_lower(alphanam(m))
     dqfxnam(m) = 'DQFX'//cnst_name(m)
     dqfxnam(m) = to_lower(dqfxnam(m))
  end do
#endif
  call check_energy_register()

end subroutine initindx
