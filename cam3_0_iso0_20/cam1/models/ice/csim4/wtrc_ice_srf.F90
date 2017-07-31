#include "misc.h"
!-----------------------------------------------------------------------
!
! Purpose:
!  Compute surface fluxes between atmos/sea ice; update the snow pack, and ice
!  thickness.
!
!
! THIS IS REALLY PLACE HPOLDER CODE - make this more realistic!
!
! Notes:
!   Work on ice points (within chunck)
!
!   Does not include blowing snow or sea spray.
!
!
! Author: David Noone <dcn@colorado.edu> - Mon Mar 15 17:08:41 MST 2004
!   (based on CC's ice_srf.F90)
!
!-----------------------------------------------------------------------
subroutine wtrc_ice_srf(npts   , indx   , dtime,  &
               wtsubi , wtsubs , wtdhs  , wtdhsf , wtdhsm , wtsnowfall, &
               fhsadj , fhsred , qbot   , wtsnowh, wthi   , qflx   )


  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols
  use constituents,    only: pcnst, pnats
  use abortutils,     only : endrun
  
  use ice_dh,          only: prognostic_icesnow, fixice
  use ice_constants,   only: rhofresh, rhos, rhoi
  use water_tracers,   only: ixwti, ixwtx, wtrc_is_vap, iwspec, wtrc_ratio, &
                             wtrc_qchk1, ixh2oq
  use water_isotopes,  only: wiso_get_roce, wiso_get_rsic

  implicit none

!---------------------------- Arguments --------------------------------

  integer , intent(in)    :: npts              ! number of ice points
  integer , intent(in)    :: indx(pcols)       ! indicies to ice points
  real(r8), intent(in)    :: dtime             ! model timestep
  real(r8), intent(in)    :: wtsubi(pcols)     ! sublimation from bare ice
  real(r8), intent(in)    :: wtsubs(pcols)     ! sublimation from snow
  real(r8), intent(in)    :: wtdhs(pcols)      ! melting (?) from snow
  real(r8), intent(in)    :: wtdhsf(pcols)     ! flooding from snow
  real(r8), intent(in)    :: wtdhsm(pcols)     ! freeze/melt from snow
  real(r8), intent(in)    :: wtsnowfall(pcols,pcnst+pnats) ! net snow
  real(r8), intent(in)    :: fhsadj(pcols)     ! height over max adjustment
  real(r8), intent(in)    :: fhsred(pcols)     ! mass reduction due to ice melt
  real(r8), intent(in)    :: qbot(pcols,pcnst+pnats)    ! atmospheric vapour

  real(r8), intent(inout) :: wtsnowh(pcols,pcnst+pnats) ! tracer snow height
  real(r8), intent(inout) :: wthi(pcols,pcnst+pnats) ! tracer ice height
  real(r8), intent(inout) :: qflx(pcols,pcnst+pnats) ! tracer flux

!------------------------- Local Variables -----------------------------

  integer i,ii                  ! ice point index
  integer m                     ! tracer index

  real(r8) Ratm                 ! isotope ratio of atmospheric vapour
  real(r8) Rocn                 ! isotope ratio of ocean
  real(r8) Rsno                 ! isotope ratio of snow pack
  real(r8) Rice                 ! isotope ratio of ice pack
  real(r8) Rbice                ! isotope ratio of bare ice
  real(r8) Rsubi                ! isotope ratio of ice sublimation
  real(r8) Rsubs                ! isotope ratio of snow sublimation
  real(r8) Rdhs                 ! isotope ratio of change in snow from melt/sublim
  real(r8) Rdhsf                ! isotope ratio of change in snow from flooding

  real(r8) wths(pcols,pcnst+pnats) ! tracer snow height
  real(r8) dsnow		! total snow tendency
  real(r8) fsnow		! fraction of snow loss that is snow water

  real(r8) hxs			! excess snow needed to magically appear
  real(r8) was          	! saved

!
#undef DBGDIAGS
   integer, parameter :: idbg = 3

!-----------------------------------------------------------------------
!!   write(*,*) 'WTRC_ICE_SRF - starting...'
!!   write(*,*) '     in :',wtsnowh(3,4),wtsnowh(3,1)

!
! Check runtime flags
!
  if (.not. prognostic_icesnow) then
     write(*,*) 'WTRC_ICE_SFC - must gave prognostics snow and ice'
     call endrun
  end if
!
! Loop over all tracers, and all ice points
!
  do m = ixwti,ixwtx
    if (wtrc_is_vap(m)) then
      Rocn  = wiso_get_roce(iwspec(m))         ! ratio of ocen water
      Rbice = wiso_get_rsic(iwspec(m))          ! ratio of bare ice

!!      write(*,*) 'TRACER:',m,  Rocn,  Rbice 
!
      do ii = 1, npts
        i = indx(ii)
!
        Ratm = wtrc_ratio(qbot(i,m), qbot(i,1))
!
! Update the ice thickness
!
!   - for now set this is a known value for isotopes, although some place-holder
!     code is included as the "if" along the lines on the snow budget, which
!     would could evaluate if you REALLY need to (should make trivial difference
!     to precipitation simulations)
!
!     This type of calulation is only valid if there is some way for the
!     ice to not have the composition of fronzen sea water. This could
!     be if there was some temperature dependent fractionation during
!     freezing, sublimation, OR flooding of (highly depleted) snow.
!     In any case, assuming constant ratio is probably fine for all but
!     the most extreme calulations (and possibly fresh water flux)
!     estimates.
!

!!        if (.not. fixice) then	! SOM
!!           dhi   = Rsubi*wtsubi(i) + Rhit*dhit(i) + Rhib*wtdhib(i) + Rbice*wtdhif(i)
!!           hi(i) = max (0._r8, hi(i)+dhi)
!!        else				! DOM
           wthi(i,m) = Rbice*wthi(i,1)
           Rice = wtrc_ratio(wthi(i,m),wthi(i,1))
!!        end if
!
! Accumulate falling snow
!   - also part do total for computing ratios 
!     (do this for all tracers as it is easy.. not fast)
!
#ifdef DBGDIAGS
        if (m == 4 .and. i == idbg) then
              write(*,*) '   i  :',wtsnowh(i,m),wtsnowh(i,1)
              write(*,*) '     snow  :',wtsnowfall(i,m)*dtime*rhofresh/rhos, &
                                        wtsnowfall(i,1)*dtime*rhofresh/rhos
        end if
#endif
        wths(i,1) = (wtsnowh(i,1) + wtsnowfall(i,1)*dtime)*rhofresh/rhos
        wths(i,m) = (wtsnowh(i,m) + wtsnowfall(i,m)*dtime)*rhofresh/rhos
#ifdef DBGDIAGS
        if (m == 4 .and. i == idbg) then
              write(*,*) '   s  :',wths(i,m),wths(i,1)
        end if
#endif
        wths(i,1) = fhsadj(i)*wths(i,1)
        wths(i,m) = fhsadj(i)*wths(i,m)
#ifdef DBGDIAGS
        if (m == 4 .and. i == idbg) then
              write(*,*) '   a  :',wths(i,m),wths(i,1)
        end if
#endif
!
! Assume trivial bits of snow have merged with ice pack (saves numerics)
!
!       if (wths(i,1) > 1.e-12) then		! not quite precise enough
       if (wths(i,1) > 1.e-9) then
!        if (wths(i,1) > 1.e-6) then		! 1/1000 of mm
          Rsno = wtrc_ratio(wths(i,m),wths(i,1))
        else
          Rsno = Rbice
        end if
!
! Apply and adjustment accounting for height change with mas constant,
! but changeing area of ice due to freeze/melt (SOM only)
!
       wths(i,1) = wths(i,1) + wtdhsm(i)
       wths(i,m) = wths(i,m) + Rsno*wtdhsm(i)
!
! Set up uni-directional flux tracer ratios for snow
! Need to check that the ice does not become bare within step.
!  if so, compute a weighted isotope ratio of flux
!
        if (.not. fixice) then
!          dsnow = wtsubs(i) +  wtdhs(i) + wtdhsf(i)
          dsnow = min(wtsubs(i),0.) +  min(wtdhs(i),0.) + min(wtdhsf(i),0.)
        else
!	   dsnow = wtsubs(i) +  wtdhs(i)
          dsnow = min(wtsubs(i),0.) +  min(wtdhs(i),0.) 
        end if
!
        if (dsnow < 0.) then
          fsnow = -wths(i,1)/dsnow
          fsnow = max(min(fsnow, 1.0),0.)
        else
          fsnow = 0.
        end if
#ifdef DBGDIAGS
        if (m == 4 .and. i == idbg) then
          write (*,*) 'dsnow:',dsnow, wths(i,1)
          write (*,*) 'fsnow:',fsnow,max(min(fsnow, 1.0),0.)
        end if
#endif
!
        if (wtdhs(i) < 0) then            ! melt
          Rdhs = Rsno
          Rsubs = fsnow*Rsno + (1.-fsnow)*Rbice
#ifdef DBGDIAGS
          if (m.eq.4 .and. i == idbg) write(*,*) 'RDHS = snow:',Rsno
#endif
        else                              ! freeze... of ocean?
          Rdhs = Rbice 
#ifdef DBGDIAGS
          if (m.eq.4 .and. i == idbg) write(*,*) 'RDHS = bare ice:',Rbice
#endif
        end if

        if (wtdhsf(i) < 0) then            ! melt
          Rdhsf = Rsno
          Rsubs = fsnow*Rsno + (1.-fsnow)*Rbice
#ifdef DBGDIAGS
          if (m.eq.4 .and. i == idbg) write(*,*) 'RDHSF = snow:',Rsno
#endif
        else                              ! freeze... of ocean?
          Rdhsf = Rocn
#ifdef DBGDIAGS
          if (m.eq.4 .and. i == idbg) write(*,*) 'RDHSF = ocean:',Rocn
#endif
        end if

        if (wtsubs(i) < 0.) then          ! sublimation
           Rsubs = fsnow*Rsno + (1.-fsnow)*Rbice
        else                              ! rime
           Rsubs = Ratm
        end if
!
! If more than all ice sublimated modify to get some ocean too.
        if (wtsubi(i) < 0.) then          ! sublimation
           Rsubi = Rice	
        else                              ! rime
           Rsubi = Ratm
        end if
!
! Update the snow height
!
!!        if (wths(i,1) > 0.) then
!!!!        if (wths(i,1) > 1.e-12) then    ! 1 trillionth of a meter
          was = wths(i,m)
          
          if (.not. fixice) then
            wths(i,1) = wths(i,1) +       wtsubs(i) +      wtdhs(i) +       wtdhsf(i)
            wths(i,m) = wths(i,m) + Rsubs*wtsubs(i) + Rdhs*wtdhs(i) + Rdhsf*wtdhsf(i)
          else
            wths(i,1) = wths(i,1) +       wtsubs(i) +      wtdhs(i) 
            wths(i,m) = wths(i,m) + Rsubs*wtsubs(i) + Rdhs*wtdhs(i) 
          end if
!!        else
!!           wths(i,m) = 0.
!!        end if
#ifdef DBGDIAGS
        if (m == 4 .and. i == idbg) then
              write(*,*) '   p  :',wths(i,m),wths(i,1)
              write(*,*) '  RATS:',Rsubs,           Rdhs,          Rbice
              write(*,*) ' terms:',      wtsubs(i),      wtdhs(i),       wtdhsf(i)
              write(*,*) 'Rterms:',Rsubs*wtsubs(i), Rdhs*wtdhs(i), Rbice*wtdhsf(i)
        end if
#endif
!
! Account for all post budget reductions (due to lateral melt of ice,
! tiny ice fraction, etc). If this is an increment, assume new "snow" is freezing
! ocean water. 
!
        if (fhsred(i) > 1.0) then		! SOM only (THIS IS A HACK)
           hxs = wths(i,1)*(fhsred(i) - 1.)
           wths(i,1) = wths(i,1) + hxs
           wths(i,m) = wths(i,m) + hxs*Rice
!!           wths(i,1) = fhsred(i)*wths(i,1)
!!           wths(i,m) = fhsred(i)*wths(i,m)
        else		! snow falls into ocean when ice melts
           wths(i,1) = fhsred(i)*wths(i,1)
           wths(i,m) = fhsred(i)*wths(i,m)
        end if

#ifdef DBGDIAGS
        if (m == 4 .and. i == idbg) then
              write(*,*) '   q  :',wths(i,m),wths(i,1)
        end if
#endif
!
! Do a test to make sure I didn't screw up
!
        if (m == 4 .and. wths(i,m) < -1.e-14) then	! still damn small
          write(*,*) '----------------------'
          write(*,*) 'WTRC_ICE_SRF - negative snow mass!',i,m,fixice
          write(*,*) 'SNOW, tot:',wtsnowfall(i,m)*dtime, wths(i,1)
          write(*,*) 'HS (old,new):',was, wths(i,m)
          write(*,*) 'R:',Rsubs,  Rdhs,  Rbice
          write(*,*) wtsubs(i), wtdhs(i), wtdhsf(i)
          write(*,*) 'ERROR you need to know is:',wths(i,m)-wths(i,1)
!!          call endrun('wtrc_ice - negative violoation')
        end if
        wths(i,m) = max(wths(i,m), 0.0)
!
! Do one more check to keep us on the rails (this is like a qneg)
!
!        if (wths(i,1) <= 0.) then	! ZERO
        if (wths(i,1) <= 1.e-22) then	! damn near ZERO (units are meters w.e.)
!           wths(i,m) = Rbice*wths(i,1)
           wths(i,m) = 0.
        end if
!
! Debug check 
!
#ifdef DBGDIAGS
        if (m == 4 .and.i.eq.idbg ) then
          write(*,*) 'TRACER CHECK:',i
          write(*,*) 'new',wths(i,m),wths(i,1)
          write(*,*) 'old',was
          write(*,*) 'fhsadj:',fhsadj(i)
        end if
#endif
!
! Compute constituent fluxes (-evap)
!
        qflx(i,m) = -(Rsubi*rhoi*wtsubi(i) + Rsubs*rhos*wtsubs(i))/dtime
!
! Copy local snow height mass back to external snow height as in calling code
!
        wtsnowh(i,m) = wths(i,m)*rhos/rhofresh
#ifdef DBGDIAGS
        if (m == 4 .and. i == idbg) then
           write(*,*) 'WTSH_out   :',wths(i,m),wths(i,1)
           write(*,*) 'WTSNOWH_out:',wtsnowh(i,m), wtsnowh(i,1)*rhos/rhofresh
        end if
#endif
!
      end do    ! i, ice point
!
    end if      ! m, is water vapour
  end do        ! m, tracers
!
! Assign locally computed total tracer for testing outout
!
  do ii = 1, npts
    i = indx(ii)
    wtsnowh(i,1) = wths(i,1)*rhos/rhofresh
  end do
!
! Do a check - this can be 1.e-15
!
  call wtrc_qchk1('wtrc_ice','snowh',npts, wtsnowh(:,ixh2oq) , wtsnowh(:,1), 1.e-14)
!
!!   write(*,*) 'WTRC_ICE_SRF - done.'
  return
end subroutine wtrc_ice_srf
