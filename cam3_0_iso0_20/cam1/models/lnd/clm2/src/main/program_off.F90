#include <misc.h>
#include <preproc.h>

#if (defined OFFLINE)

!================================================================================
! CCSM was developed in cooperation with the National Science Foundation,
! the Department of Energy Los Alamos National Laboratory (LANL), the
! National Aeronautics and Space Administration Data Assimilation Office,
! and the University Corporation for Atmospheric Research's National
! Center for Atmospheric Research.*  Except for POP, SCRIP, and CICE,
! which are segregable components of the model, CCSM is public domain
! software.
! 
! POP, SCRIP, and CICE were developed at LANL and are protected by the
! following copyright notice.
! 
! Copyright �2002 The Regents of the University of California. All Rights
! Reserved.
! 
! Permission to use, copy, modify, and distribute this software and its
! documentation for educational, research and non-profit purposes, without
! fee, and without a written agreement is hereby granted, provided that
! the above copyright notice, this paragraph and the following three
! paragraphs appear in all copies.
! 
! Permission to incorporate this software into commercial products may be
! obtained by contacting the University of California, Charles Rzeszutko,
! Campus Liaison Officer Office of Technology Transfer, 1111 Franklin
! Street, 5th Floor, Oakland, CA 94607, (510) 587-6063,
! charles.rzeszutko@ucop.edu.
! 
! THIS SOFTWARE PROGRAM AND DOCUMENTATION ARE COPYRIGHTED BY
! THE REGENTS OF THE UNIVERSITY OF CALIFORNIA. THE SOFTWARE
! PROGRAM AND DOCUMENTATION ARE SUPPLIED "AS IS", WITHOUT ANY
! ACCOMPANYING SERVICES FROM THE REGENTS. THE REGENTS DOES NOT
! WARRANT THAT THE OPERATION OF THE PROGRAM WILL BE
! UNINTERRUPTED OR ERROR-FREE. THE END-USER UNDERSTANDS THAT
! THE PROGRAM WAS DEVELOPED FOR RESEARCH PURPOSES AND IS
! ADVISED NOT TO RELY EXCLUSIVELY ON THE PROGRAM FOR ANY
! REASON.
! 
! IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY
! PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
! CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
! THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE
! UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
! SUCH DAMAGE.  THE UNIVERSITY OF CALIFORNIA SPECIFICALLY
! DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN
! "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO
! OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
! ENHANCEMENTS, OR MODIFICATIONS.
! 
! *The National Center for Atmospheric Research is sponsored by the
! National Science Foundation.  Any opinions, findings and conclusions or
! recommendations expressed in this publication are those of the author(s)
! and do not necessarily reflect the views of the National Science
! Foundation.
!================================================================================

!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: program_off
!
! !INTERFACE:
PROGRAM program_off
!
! !DESCRIPTION:
! "off-line" code to mimic coupling to an atmospheric model.
! This program is an "off-line" driver for clm2.
! This code can be used to run the clm2 uncoupled from any atmospheric model.
! The appropriate atmospheric forcing is provided in module [atmdrvMod.F90]
! o If running as an offline driver, the land surface model may use
!   a different grid than the input atmospheric data. The atmospheric
!   data is then interpolated to the land model grid inside the
!   atmospheric driver module [atmdrvMod.F90].
! o If running as part of cam, the land surface model must use the
!   same grid as the cam.
! o If running through the flux coupler, the land surface model grid
!   is interpolated to the atmospheric grid inside the flux coupler
! o To map from the atmospheric grid to the land grid, the atmospheric
!   model must provide latitudes and longitudes (degrees) for each grid
!   point and the North, East, South, and West edges of atmospheric grid.
!   Comparable data for the land grid are provided by the land model.
!   When mapping from land to atm grid, an atm grid cell that is part
!   land and part ocean (as defined by the land surface grid) will have
!   fluxes only based on the land portion.
! o The zenith angle calculation is for the NEXT time step rather
!   than the current time step. Make sure the calendar day is for
!   the NEXT time step. Make sure the calendar day is for Greenwich
!   time (see next comment).
! o The land surface model calculates its own net solar radiation and
!   net longwave radiation at the surface. The net longwave radiation
!   at the surface will differ somewhat from that calculated in the
!   atmospheric model because the atm model will use the upward
!   longwave flux (or radiative temperature) from the previous time
!   step whereas the land surface model uses the flux for the current
!   time step. The net solar radiation should equal that calculated
!   in the atmospheric model. If not, there is a problem in how
!   the models are coupled.
!
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8
  use clm_varpar
  use clm_varctl   , only : irad, nsrest
  use initializeMod, only : initialize
  use atmdrvMod    , only : atmdrv
#if (defined SPMD)
  use spmdMod      , only : masterproc, iam, mpicom, spmd_init
#else
  use spmdMod      , only : masterproc, iam
#endif
  use shr_orb_mod          !orbital parameters and routines
  use time_manager, only : is_last_step, advance_timestep, get_nstep, get_step_size
!
! !ARGUMENTS:
    implicit none
#include "gpt.inc"
!
! !REVISION HISTORY:
! Author: Gordon Bonan and Mariana Vertenstein
! 11/30/01 Peter Thornton : Added use globals, removed use clm_varctl
!
!EOP
!
! !LOCAL VARIABLES:
  logical doalb     !true if surface albedo calculation time step
  integer nstep     !time step
  integer ier       !error code

! Earth's orbital characteristics

  integer iyear_AD  !Year (AD) to simulate above earth's orbital parameters for
  real(r8) eccen    !Earth's eccentricity factor (unitless) (typically 0 to 0.1)
  real(r8) obliq    !Earth's obliquity angle (degree's) (-90 to +90) (typically 22-26)
  real(r8) mvelp    !Earth's moving vernal equinox at perhelion (degree's) (0 to 360.0)

! Orbital information after call to routine shr_orbit_params

  real(r8) obliqr   !Earth's obliquity in radians
  real(r8) lambm0   !Mean longitude (radians) of perihelion at the vernal equinox
  real(r8) mvelpp   !Earth's moving vernal equinox longitude
  logical log_print !true=> print diagnostics
!-----------------------------------------------------------------------

  ! -----------------------------------------------------------------
  ! Initialize timing library and mpi communication
  ! -----------------------------------------------------------------

  ! Initialize timing library.  2nd arg 0 means disable, 1 means enable

  call t_setoptionf (usrsys, 0)
  call t_initializef ()

#if (defined SPMD)
  ! Initialize intra-MPI communication info

  call spmd_init()
#endif

  ! -----------------------------------------------------------------
  ! Orbital parameters.
  ! variables obliq, eccen and mvelp determined based on value of
  ! iyear_AD
  ! -----------------------------------------------------------------

   if (masterproc) then
     log_print = .true.
   else
     log_print = .false.
   end if
   iyear_AD = 1950
   obliq    = SHR_ORB_UNDEF_REAL
   eccen    = SHR_ORB_UNDEF_REAL
   mvelp    = SHR_ORB_UNDEF_REAL
   call shr_orb_params (iyear_AD, eccen, obliq, mvelp, obliqr, &
                        lambm0, mvelpp, log_print)

   ! -----------------------------------------------------------------
   ! Land surface model initialization.
   ! o Run control parameters (length of integration, starting date, etc)
   !   are set in routine initialize via the clmexp namelist.
   !   Land surface model dataset names are set in the clmexp namelist.
   ! o Initializes albedos [asdirxy], [asdifxy], [aldirxy], [aldifxy],
   !   surface temperature [tsxy], upward longwave radiation [lwupxy],
   !   and snow [snowxy] for land points only on the atmospheric grid.
   !   These are zero for non-land points and must be set by
   !   the appropriate surface model.
   ! -----------------------------------------------------------------

  call initialize (eccen, obliqr, lambm0, mvelpp)

  ! -----------------------------------------------------------------
  ! Time stepping loop
  ! -----------------------------------------------------------------

  call t_startf('total')

  ! begin time stepping loop
  do

     ! Current atmospheric state and fluxes for all [atmlon] x [atmlat] points.
     ! When coupling to an atmospheric model: solar radiation depends on
     ! surface albedos from the previous time step (based on current
     ! surface conditions and solar zenith angle for next time step).
     ! Longwave radiation depends on upward longwave flux from previous
     ! time step.

     nstep = get_nstep()
     call atmdrv(nstep)

     ! doalb is true when the next time step is a radiation time step
     ! this allows for the fact that an atmospheric model may not do
     ! the radiative calculations every time step. for example:
     !      nstep dorad doalb
     !        1     F     F
     !        2     F     T
     !        3     T     F
     ! The following expression for doalb is specific to CAM

     doalb = (irad==1 .or. (mod(nstep,irad)==0 .and. nstep+1/=1))

     ! Call land surface model driver
     ! Note that surface fields used by the atmospheric model are zero for
     ! non-land points and must be set by the appropriate surface model

     call driver (doalb, eccen, obliqr, lambm0, mvelpp)

     ! determine if time to stop

     if (is_last_step()) exit

     ! increment time step

     call advance_timestep()

  end do
  call t_stopf('total')

  ! -----------------------------------------------------------------
  ! Exit gracefully
  ! -----------------------------------------------------------------

  if (masterproc) then
     write(6,*)'SUCCESFULLY TERMINATING CLM MODEL at nstep= ',get_nstep()
  endif
  call t_prf(iam)
#if (defined SPMD)
  call mpi_barrier (mpicom, ier)
  call mpi_finalize(ier)
#endif

  stop
end program program_off

#else

!The following is only here since empty file won't compile
subroutine program_off_stub
  write(6,*) 'PROGRAM_OFF: this routine should not be called'
  return
end subroutine program_off_stub

#endif
