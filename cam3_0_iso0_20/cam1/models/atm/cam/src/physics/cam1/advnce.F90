#include <misc.h>
#include <params.h>

subroutine advnce
!-----------------------------------------------------------------------
!
! Purpose: 
! Advance time information
!
! Method: 
!
! Author: CCM1, CMS Contact: J. Truesdale
!
!-----------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use chemistry, only: trace_gas, chem_timestep_init
  use aerosol_intr, only: aerosol_time_interp
  use ghg_surfvals, only: ghg_surfvals_ramp, ghg_surfvals_set
  use tracers, only: tracers_flag, tracers_timestep_init
  use time_manager, only: get_nstep
  use volcanicmass, only: read_volcanic_mass
  use prescribed_aerosols, only: aerint, strat_volcanic
  use ramp_scon, only: ramp_sconst

  implicit none

#include <comctl.h>
#include <comlun.h>
!
!-----------------------------------------------------------------------
!
! Local workspace
!
  integer :: nstep             ! current timestep number
!-----------------------------------------------------------------------

  nstep = get_nstep()
!
! Determine whether it is time for a shortwave or longwave radiation 
! calculation
!
  dosw = nstep.eq.0 .or. iradsw.eq.1 .or. (mod(nstep-1,iradsw).eq.0 .and. nstep.ne.1)
  dolw = nstep.eq.0 .or. iradlw.eq.1 .or. (mod(nstep-1,iradlw).eq.0 .and. nstep.ne.1)
!
! Determine whether it is time for an absorptivity/emissivity calculation
!
  doabsems = nstep.eq.0 .or. iradae.eq.1 .or. (mod(nstep-1,iradae).eq.0 .and. nstep.ne.1)
  aeres = (mod(nstep,iradae).ne.0)
!
! Update ozone and aerosol data on shortwave or longwave time step. 
! Note that the ozone data is not needed on a longwave time step unless the
! absorptivities are being updated ("doabsems").
!
  if (dosw .or. dolw) then
     call oznint ()
     call aerint ()
     if(strat_volcanic) then
       call read_volcanic_mass
     end if
  end if
!
! Ramping ghg if appropriate
!
  if (ghg_surfvals_ramp()) call ghg_surfvals_set()
!
! Ramp solar constant if appropraite
!
  if (doRamp_scon) call ramp_sconst
!
! Time interpolate for chemistry, if appropriate
!
  if (trace_gas) call chem_timestep_init

! get the aerosol surface fluxes for this time step
  call aerosol_time_interp
!
! Time interpolate for tracers, if appropriate
!
  if (tracers_flag) call tracers_timestep_init


  return
end subroutine advnce
