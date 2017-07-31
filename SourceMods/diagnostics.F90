#include <misc.h>
#include <params.h>

!--------------------------------------------------------------------------------
! @222:
! #ifdef TAGGING
!   use water_tracers, only: trace_water,ixwti,ixwtx,wtrc_is_wtrc, wtrc_is_vap,tagmask
! #else
!    use water_tracers, only: trace_water,ixwti,ixwtx,wtrc_is_wtrc, wtrc_is_vap
! #endif
!
! @356:
!   #ifdef TAGGING
!   use water_tracers, only: trace_water,ixwti,ixwtx,wtrc_is_wtrc, wtrc_is_vap,tagmask
!   #endif
!
!
!--------------------------------------------------------------------------------

module diagnostics

    use shr_kind_mod, only: r8 => shr_kind_r8
    use ppgrid,  only: pcols, pver, pvermx
    use history, only: outfld
    use constituents, only: pcnst, pnats, cnst_name

    implicit none

!---------------------------------------------------------------------------------
! Module to compute a variety of diagnostics quantities for history files
!---------------------------------------------------------------------------------

contains

!===============================================================================
  subroutine diag_dynvar(lchnk, ncol, state)

!----------------------------------------------------------------------- 
! 
! Purpose: record dynamics variables on physics grid
!
!-----------------------------------------------------------------------
    use physics_types, only: physics_state
    use physconst,     only: gravit, rga, rair
    use wv_saturation, only: aqsat
#if ( defined COUP_CSM )
    use ccsm_msg, only: psl   ! Store sea-level pressure for CCSM
#endif
#if ( defined SCAM )
use pmgrid,    only: plev, plevp
#include <comfrc.h>
#endif
!-----------------------------------------------------------------------
!
! Arguments
!
    integer,  intent(in) :: lchnk            ! chunk identifier
    integer,  intent(in) :: ncol             ! longitude dimension

    type(physics_state), intent(inout) :: state
!
!---------------------------Local workspace-----------------------------
!
    real(r8) ftem(pcols,pver) ! temporary workspace
    real(r8) psl_tmp(pcols)   ! Sea Level Pressure
    real(r8) z3(pcols,pver)   ! geo-potential height
    real(r8) p_surf(pcols)    ! data interpolated to a pressure surface
    real(r8) tem2(pcols,pver) ! temporary workspace

    integer k, m              ! index
!
!-----------------------------------------------------------------------
!
!This field has the same name as the one that is needed for BFB SCAM
!IOP datasets so don't outfield it here (outfield in tfiltmassfix)
!
    call outfld('T       ',state%t , pcols   ,lchnk   )
    call outfld('PS      ',state%ps, pcols   ,lchnk   )
    call outfld('U       ',state%u , pcols   ,lchnk   )
    call outfld('V       ',state%v , pcols   ,lchnk   )
    do m=1,pcnst+pnats
       call outfld(cnst_name(m),state%q(1,1,m),pcols ,lchnk )
    end do
    call outfld('PDELDRY ',state%pdeldry, pcols,   lchnk     )
    call outfld('PHIS    ',state%phis,    pcols,   lchnk     )
#if (defined BFB_CAM_SCAM_IOP )
    call outfld('phis    ',state%phis,    pcols,   lchnk     )
#endif
!
! Add height of surface to midpoint height above surface 
!
    do k = 1, pver
       z3(:ncol,k) = state%zm(:ncol,k) + state%phis(:ncol)*rga
    end do
    call outfld('Z3      ',z3,pcols,lchnk)
!           
! Output Z3 on 500mb, 300, 50 and 700 mb surface
!
    call vertinterp(ncol, pcols, pver, state%pmid, 70000._r8, z3, p_surf)
    call outfld('Z700    ', p_surf, pcols, lchnk)
    call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, z3, p_surf)
    call outfld('Z500    ', p_surf, pcols, lchnk)
    call vertinterp(ncol, pcols, pver, state%pmid,  5000._r8, z3, p_surf)
    call outfld('Z050    ', p_surf, pcols, lchnk)
    call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, z3, p_surf)
    call outfld('Z300    ', p_surf, pcols, lchnk)
!
! Quadratic height fiels Z3*Z3
!
    ftem(:ncol,:) = z3(:ncol,:)*z3(:ncol,:)
    call outfld('ZZ      ',ftem,pcols,lchnk)

    ftem(:ncol,:) = z3(:ncol,:)*state%v(:ncol,:)*gravit
    call outfld('VZ      ',ftem,  pcols,lchnk)
!
! Meridional advection fields
!
    ftem(:ncol,:) = state%v(:ncol,:)*state%t(:ncol,:)
    call outfld ('VT      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:)*state%q(:ncol,:,1)
    call outfld ('VQ      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:)**2
    call outfld ('VV      ',ftem    ,pcols   ,lchnk     )

    ftem(:ncol,:) = state%v(:ncol,:) * state%u(:ncol,:)
    call outfld ('VU      ',ftem    ,pcols   ,lchnk     )

! zonal advection

    ftem(:ncol,:) = state%u(:ncol,:)**2
    call outfld ('UU      ',ftem    ,pcols   ,lchnk     )

! Wind speed
    ftem(:ncol,:) = sqrt( state%u(:ncol,:)**2 + state%v(:ncol,:)**2)
    call outfld ('WSPEED  ',ftem    ,pcols   ,lchnk     )

! Vertical velocity and advection

#if ( defined SCAM )
    call outfld('OMEGA   ',wfld,    pcols,   lchnk     )
#else
    call outfld('OMEGA   ',state%omega,    pcols,   lchnk     )
#endif

#if (defined BFB_CAM_SCAM_IOP )
    call outfld('omega   ',state%omega,    pcols,   lchnk     )
#endif

    ftem(:ncol,:) = state%omega(:ncol,:)*state%t(:ncol,:)
    call outfld('OMEGAT  ',ftem,    pcols,   lchnk     )
    ftem(:ncol,:) = state%omega(:ncol,:)*state%u(:ncol,:)
    call outfld('OMEGAU  ',ftem,    pcols,   lchnk     )
!
! Output omega at 850 and 500 mb pressure levels
!
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%omega, p_surf)
    call outfld('OMEGA850', p_surf, pcols, lchnk)
    call vertinterp(ncol, pcols, pver, state%pmid, 50000._r8, state%omega, p_surf)
    call outfld('OMEGA500', p_surf, pcols, lchnk)
!     
! Mass of q, by layer and vertically integrated
!
    ftem(:ncol,:) = state%q(:ncol,:,1) * state%pdel(:ncol,:) * rga
    call outfld ('MQ      ',ftem    ,pcols   ,lchnk     )

    do k=2,pver
       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    call outfld ('TMQ     ',ftem, pcols   ,lchnk     )
!
! Relative humidity
!
    call aqsat (state%t    ,state%pmid  ,tem2    ,ftem    ,pcols   , &
         ncol ,pver  ,1       ,pver    )
    ftem(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100.
    call outfld ('RELHUM  ',ftem    ,pcols   ,lchnk     )
!
! Sea level pressure
!
    call cpslec (ncol, state%pmid, state%phis, state%ps, state%t,psl_tmp, gravit, rair) 
    call outfld ('PSL     ',psl_tmp  ,pcols, lchnk     )
#if ( defined COUP_CSM )
    psl(:ncol,lchnk) = psl_tmp(:ncol)
#endif
!
! Output T,q,u,v fields on pressure surfaces
!
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%t, p_surf)
    call outfld('T850    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 30000._r8, state%t, p_surf)
    call outfld('T300    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%q(1,1,1), p_surf)
    call outfld('Q850    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%q(1,1,1), p_surf)
    call outfld('Q200    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%u, p_surf)
    call outfld('U850    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%u, p_surf)
    call outfld('U200    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 85000._r8, state%v, p_surf)
    call outfld('V850    ', p_surf, pcols, lchnk )
    call vertinterp(ncol, pcols, pver, state%pmid, 20000._r8, state%v, p_surf)
    call outfld('V200    ', p_surf, pcols, lchnk )

    ftem(:ncol,:) = state%t(:ncol,:)*state%t(:ncol,:)
    call outfld('TT      ',ftem    ,pcols   ,lchnk   )

#if ( defined COUP_CSM )
!
! Output U, V, T, Q, P and Z at bottom level
!
    call outfld ('UBOT    ', state%u(1,pver)  ,  pcols, lchnk)
    call outfld ('VBOT    ', state%v(1,pver)  ,  pcols, lchnk)
    call outfld ('QBOT    ', state%q(1,pver,1),  pcols, lchnk)
    call outfld ('ZBOT    ', state%zm(1,pver) , pcols, lchnk)
#endif

    return
  end subroutine diag_dynvar

!===============================================================================

subroutine diag_surf (srfflx, surface, icefrac, ocnfrac, landfrac, &
                      sicthk, snowhland, snowhice, tsice, trefmxav, &
                      trefmnav, wtsicthk, wtsnowhland, wtsnowhice )

!----------------------------------------------------------------------- 
! 
! Purpose: record surface diagnostics
!
!-----------------------------------------------------------------------

   use comsrf, only: srfflx_state, surface_state, tsnam
   use constituents,  only: ppcnst, sflxnam, cnst_name
#ifdef TAGGING
   use water_tracers, only: trace_water,ixwti,ixwtx,wtrc_is_wtrc, wtrc_is_vap,tagmask
#else
   use water_tracers, only: trace_water,ixwti,ixwtx,wtrc_is_wtrc, wtrc_is_vap
#endif

#if ( defined COUP_CSM )
    use time_manager, only: is_end_curr_day
#endif
!-----------------------------------------------------------------------
!
! Input arguments
!
    type(srfflx_state),  intent(in) :: srfflx
    type(surface_state), intent(in) :: surface

    real(r8), intent(in) :: icefrac(pcols)   ! ice fraction
    real(r8), intent(in) :: ocnfrac(pcols)   ! ocn fraction
    real(r8), intent(in) :: landfrac(pcols)  ! land fraction
    real(r8), intent(in) :: sicthk(pcols)    ! sea-ice thickness
    real(r8), intent(in) :: snowhland(pcols) ! equivalent liquid water snow depth
    real(r8), intent(in) :: snowhice(pcols)  ! equivalent liquid water snow depth
    real(r8), intent(in) :: tsice(pcols)     ! surface T over seaice

    real(r8), intent(in) :: wtsicthk(pcols,ppcnst)    ! tracer sea-ice thickness
    real(r8), intent(in) :: wtsnowhland(pcols,ppcnst) ! tracer equivalent liquid water snow depth
    real(r8), intent(in) :: wtsnowhice(pcols,ppcnst)  ! tracer equivalent liquid water snow depth

    real(r8), intent(inout) :: trefmnav(pcols) ! daily minimum tref  
    real(r8), intent(inout) :: trefmxav(pcols) ! daily maximum tref
!
!---------------------------Local workspace-----------------------------
!
    integer i,k,m           ! indexes
    integer :: lchnk        ! chunk identifier
    integer :: ncol         ! longitude dimension
!
    real(r8) :: ftmp(pcols)
    character(len=2) atr    ! tracer number as a string
!
!-----------------------------------------------------------------------
!
    lchnk = srfflx%lchnk
    ncol  = srfflx%ncol

    call outfld('SHFLX',    srfflx%shf,       pcols, lchnk)
    call outfld('LHFLX',    srfflx%lhf,       pcols, lchnk)
    call outfld('QFLX',     srfflx%cflx(1,1), pcols, lchnk)

! add unidirectional flux diagnostics (dcn)
    ftmp(:) = 0.
    where (srfflx%cflx(:,1) > 0.) ftmp(:) = srfflx%cflx(:,1)
    call outfld('QFLXUP',   ftmp,  pcols, lchnk)
    ftmp(:) = 0.
    where (srfflx%cflx(:,1) < 0.) ftmp(:) = srfflx%cflx(:,1)
    call outfld('QFLXDN',   ftmp,  pcols, lchnk)

! Output water tracer surface flux
    if (trace_water) then
      do m = ixwti, ixwtx
        if (wtrc_is_wtrc(m)) then
          call outfld(sflxnam(m),srfflx%cflx(1,m), pcols, lchnk)
        end if

        if (wtrc_is_vap(m)) then
          ftmp(:) = 0.
          where (srfflx%cflx(:,m) > 0.) ftmp(:) = srfflx%cflx(:,m)
          call outfld('SFU'//trim(cnst_name(m)),ftmp, pcols, lchnk)
          ftmp(:) = 0.
          where (srfflx%cflx(:,m) < 0.) ftmp(:) = srfflx%cflx(:,m)
          call outfld('SFD'//trim(cnst_name(m)),ftmp, pcols, lchnk)
        end if
      end do
    end if
!
    call outfld('TAUX',     srfflx%wsx,       pcols, lchnk)
    call outfld('TAUY',     srfflx%wsy,       pcols, lchnk)
    call outfld('TREFHT  ', srfflx%tref,      pcols, lchnk)
    call outfld('TREFHTMX', srfflx%tref,      pcols, lchnk)
    call outfld('TREFHTMN', srfflx%tref,      pcols, lchnk)
#if ( defined COUP_CSM )
    call outfld('QREFHT',   srfflx%qref,      pcols, lchnk)
#endif
#if (defined BFB_CAM_SCAM_IOP )
    call outfld('shflx   ',srfflx%shf,   pcols,   lchnk)
    call outfld('lhflx   ',srfflx%lhf,   pcols,   lchnk)
    call outfld('trefht  ',srfflx%tref,  pcols,   lchnk)
#endif
!
! Ouput ocn and ice fractions
!
    call outfld('LANDFRAC', landfrac,         pcols, lchnk)
    call outfld('ICEFRAC',  icefrac,          pcols, lchnk)
    call outfld('OCNFRAC',  ocnfrac,          pcols, lchnk)

#if ( defined COUP_CSM )
!
! Compute daily minimum and maximum of TREF
!
    do i = 1,ncol
       trefmxav(i) = max(srfflx%tref(i),trefmxav(i))
       trefmnav(i) = min(srfflx%tref(i),trefmnav(i))
    end do
    if (is_end_curr_day()) then
       call outfld('TREFMXAV', trefmxav,pcols,   lchnk     )
       call outfld('TREFMNAV', trefmnav,pcols,   lchnk     )
       trefmxav(:ncol) = -1.0e36
       trefmnav(:ncol) =  1.0e36
    endif

#else

    do k=1,pvermx
       call outfld(tsnam(k), surface%tssub(1,k), pcols, lchnk)
    end do
    call outfld('SICTHK',   sicthk,           pcols, lchnk)
    call outfld('TSICE',    tsice,            pcols, lchnk)
#endif

    call outfld('TS',       srfflx%ts,        pcols, lchnk)
    call outfld('TSMN',     srfflx%ts,        pcols, lchnk)
    call outfld('TSMX',     srfflx%ts,        pcols, lchnk)
    call outfld('SNOWHLND', snowhland,        pcols, lchnk)
    call outfld('SNOWHICE', snowhice ,        pcols, lchnk)
    call outfld('TBOT',     surface%tbot,     pcols, lchnk)

    if (trace_water) then
      atr = '  '
      do m = ixwti, ixwtx
        write(atr,'(i2.2)') m
        if (wtrc_is_vap(m)) then
#ifndef COUP_CSM
          call outfld('ITK'//trim(cnst_name(m)), wtsicthk(:,m)   ,        pcols, lchnk)
#endif
          call outfld('SHL'//trim(cnst_name(m)), wtsnowhland(:,m),        pcols, lchnk)
          call outfld('SHI'//trim(cnst_name(m)), wtsnowhice(:,m) ,        pcols, lchnk)
        end if
      end do
#ifdef TAGGING
        call outfld('TAGMASK', tagmask, pcols, lchnk)
#endif
    end if

    call outfld('ASDIR',    srfflx%asdir,     pcols, lchnk)
    call outfld('ASDIF',    srfflx%asdif,     pcols, lchnk)
    call outfld('ALDIR',    srfflx%aldir,     pcols, lchnk)
    call outfld('ALDIF',    srfflx%aldif,     pcols, lchnk)
    call outfld('SST',      srfflx%sst,       pcols, lchnk)

end subroutine diag_surf

end module diagnostics
