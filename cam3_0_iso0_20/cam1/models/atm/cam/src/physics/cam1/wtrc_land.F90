#include <misc.h>
#include <params.h>

module wtrc_land
!-----------------------------------------------------------------------
!
! Purpose: provides simple land scheme for water tracers
!
! Author: David Noone <dcn@colorado.edu> - Tue Jul 13 12:01:22 MDT 2004
!
!-----------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8

  implicit none
  private
  save
!-----------------------------------------------------------------------

  public :: wtrc_land_init
  public :: wtrc_land_drv
!
! Deardorffs land parameters
!
   real(r8) :: wgmax = 0.4              ! ground max volumetic content (m^3/m^3)
   real(r8) :: wbmax = 0.32             ! bulk max volumetic content (m^3/m^3)
   real(r8) :: wtiny = 0.001		! trivial amount to retain ratio
   real(r8) :: d1    = 0.05		! depth of upper layer (meters)
   real(r8) :: d2    = 0.15             ! depth of both layers (metres)
!   real(r8) :: d1    = 0.1		! depth of upper layer (meters)
!   real(r8) :: d2    = 0.5             ! depth of both layers (metres)
!   real(r8) :: d1    = 0.2		! depth of upper layer (meters)
!   real(r8) :: d2    = 1.0             ! depth of both layers (metres)
!   real(r8) :: d1    = 0.1		! depth of upper layer (meters)
!   real(r8) :: d2    = 0.8              ! depth of both layers (metres)
   real(r8) :: tau1  = 86400.           ! upper damping scale (seconds)
   real(r8) :: rhow  = 1.0e+03          ! density of water (kg/m^3)

!=======================================================================
contains

!=======================================================================
subroutine wtrc_land_init
!-----------------------------------------------------------------------
! Initialized module variables and input values
! Set pools initially to be half full of smow
!-----------------------------------------------------------------------
  use water_isotopes,  only: wiso_get_rstd
  use constituents,    only: cnst_get_ind
  use comsrf,          only: wtlwb, wtlwg
!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------
  integer m			! constituent index
  integer c			! chunk index
  integer mbase			! base variable index
  integer ixcldice, ixcldliq    ! constituent indices for cloud liquid and ice water
  real(r8) riso			! an isotope ratio
!-----------------------------------------------------------------------
!
#ifdef INITWTRCLAND	/* else, buckts will just spin up in a few months */
  write(*,*) '(wtrc_land_init) initialize land isotopes to SMOW.'
!
  do c = begchnk, endchnk
    wtlwb(:,1,c) = wbmax
    wtlwg(:,1,c) = wgmax

    do m = ixwti, ixwtx
      riso = wiso_get_rstd(iwspec(m))

         mbase = 0
         if (wtrc_is_vap(m)) then
           mbase = 1
         else if (wtrc_is_liq(m)) then
           mbase = ixcldliq 
         else if (wtrc_is_ice(m)) then
           mbase = ixcldice
         end if
  
         if (mbase /= 0) then		! this is some kind of water
           wtlwb(:,m,c) = riso*wtlwg(:,1,c)
           wblwb(:,m,c) = riso*wtlwb(:,1,c)
         end if
    end do
  end do		! chunks
#endif
!
  return
end subroutine wtrc_land_init

!=======================================================================
subroutine wtrc_land_drv(ncol, lchnk, ztodt, state, parm, landfrac, &
                         wtprc, wtlwg, wtlwb)

!-----------------------------------------------------------------------
! 
! Purpose: compute isotopic fluxes from land points (as we don't have clm)
!  We work on one chunk.
!
!  1) prescribe SMOW
!  2) use bucket scheme
!
! Author: David Noone <dcn@colorado.edu> - Tue Jul 13 11:45:24 MDT 2004
!-----------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid,          only: pcols
   use constituents,    only: pcnst, pnats, cnst_name, cnst_get_ind
   use history,         only: outfld
   use water_tracers,   only: trace_water, ixwti, ixwtx, iwspec, &
                              wtrc_ratio, wtrc_is_vap, wtrc_is_liq, wtrc_is_ice, &
                              wtrc_qchk1, ixh2oq
   use water_isotopes,  only: wiso_get_rstd, wiso_delta

   use comsrf,          only: srfflx_parm, srfflx_state

!-----------------------------------------------------------------------
  implicit none
!-----------------------------------------------------------------------

!---------------------------- Arguments --------------------------------
  integer, intent(in)               :: ncol
  integer, intent(in)		    :: lchnk
  real(r8), intent(in)              :: ztodt
  type(srfflx_state), intent(inout) :: state
  type(srfflx_parm), intent(inout)  :: parm
  real(r8), intent(in)		    :: landfrac(pcols)
  real(r8), intent(in)              :: wtprc(pcols,pcnst+pnats) ! net prec
  real(r8), intent(inout)	    :: wtlwg(pcols,pcnst+pnats)
  real(r8), intent(inout)	    :: wtlwb(pcols,pcnst+pnats)

!------------------------- Local Variables -----------------------------

  integer i,ii			! column index
  integer m			! tracer index
  integer mbase			! tracer index for base species
  integer ixcldice, ixcldliq    ! constituent indices for cloud liquid and ice water
!
  integer npts
  integer indx(pcols)
!
  real(r8) riso(pcols,pcnst+pnats)	! an isotope ratio
  real(r8) dwb,dwg			! tendency terms
  real(r8) prec(pcols,pcnst+pnats)	! precipitation (correct units)
  real(r8) wtrun(pcols,pcnst+pnats)	! ranoff at surface
  real(r8) wtdrn(pcols,pcnst+pnats)	! deep soil drainage
  real(r8) wbxs(pcols,pcnst+pnats)
  real(r8) wgxs(pcols,pcnst+pnats)
!
!-----------------------------------------------------------------------
     call cnst_get_ind('CLDLIQ', ixcldliq)
     call cnst_get_ind('CLDICE', ixcldice)
!
     wtrun(:,:) = 0.
     wtdrn(:,:) = 0.
     wgxs(:,:) = 0.
     wbxs(:,:) = 0.
!
! Assign "prognostic" on input (as there is no "prognostic" version of this)
!
     wtlwg(:,1) = wtlwg(:,ixh2oq)
     wtlwb(:,1) = wtlwb(:,ixh2oq)
!
     riso(:,1)        = 1.
     riso(:,ixcldliq) = 1.
     riso(:,ixcldice) = 1.
!
! Convert precipitation unitf from m/s to kg/m2/s
!
    prec(:,:) = wtprc(:,:)*rhow
!
! Set up column index
!
    npts = 0
    do i = 1, ncol
!!      if (landfrac(i) > 0.) then	! do everywhere!
         npts = npts + 1
         indx(npts) = i
!!      end if
    end do
    if (npts == 0) return		! nothing to do
!
! 1) prescribe ratio for fluxes based on surface scheme
! ( this is actually quite disatisfactory, as the "wetness" does not
!  affect cflx...so buckets can drain dry, and still have evap)
     do m = ixwti, ixwtx

       mbase = 0
       if (wtrc_is_vap(m)) then
         mbase = 1
       else if (wtrc_is_liq(m)) then
         mbase = ixcldliq 
       else if (wtrc_is_ice(m)) then
         mbase = ixcldice
       end if

       if (mbase /= 0) then		! this is some kind of water
         do ii=1,npts
            i = indx(ii)

            if (wtlwg(i,1) > wtiny) then		! to is wet, assign ratio
               riso(i,m) = wtrc_ratio(wtlwg(i,m),wtlwg(i,1))
            else if (wtlwb(i,1) > 0.99*wtiny) then	! top is dry, try bulk
               riso(i,m) = wtrc_ratio(wtlwb(i,m),wtlwb(i,1))
            else					! buckets are dry, assign smow
               riso(i,m) = wiso_get_rstd(iwspec(m))
            end if

!
! DEBUGGUING
!
!!            riso(i,m) = wiso_get_rstd(iwspec(m))

            parm%cflx(i,m)  = riso(i,m)*parm%cflx(i,mbase)
         end do
       end if
    end do

    call wtrc_qchk1('wtrc_land','s%cflx',ncol,state%cflx(:,ixh2oq),state%cflx(:,1))
    call wtrc_qchk1('wtrc_land','p%cflx',ncol, parm%cflx(:,ixh2oq), parm%cflx(:,1))

!
! 2) Update surface reservoirs:  dbucket = P - E - drain - Runoff
!
    do ii=1,npts
      i = indx(ii)
      dwb = (prec(i,1))/(d2*rhow)
      dwg = (prec(i,1))/(d1*rhow)
      dwg = dwg - (wtlwg(i,1)-wtlwb(i,1))/tau1	! drain/recharge
!
      wtlwg(i,1) = max(wtlwg(i,1) + ztodt*dwg, 0.)
      wgxs(i,1)  = max(wtlwg(i,1) - wgmax, 0.)
      wtlwg(i,1) = wtlwg(i,1) - wgxs(i,1)
!
      wtlwb(i,1) = max(wtlwb(i,1) + ztodt*dwb, 0.)
      wbxs(i,1)  = max(wtlwb(i,1) - wbmax, 0.)
      wtlwb(i,1) = wtlwb(i,1) - wbxs(i,1)
!
      wtrun(i,1) = d1*rhow*wgxs(i,1)/ztodt
      wtdrn(i,1) = d2*rhow*wbxs(i,1)/ztodt
!
! If ground bucket has dried up completely, add a tiny mass so that
! we can save the long term isotope ratio.
!
      wtlwb(i,1) = max(wtlwb(i,1), wtiny)

!!      call outfld('LWGQ', wtlwg(:,1)  ,pcols   ,lchnk     )
!!      call outfld('LWBQ', wtlwg(:,1)  ,pcols   ,lchnk     )
!!      call outfld('LWRQ', wtrun(:,1)  ,pcols   ,lchnk     )
    end do

    do m = ixwti, ixwtx
      if (wtrc_is_vap(m)) then
        do ii=1,npts
          i = indx(ii)
!
            dwb = (prec(i,m))/(d2*rhow)
            dwg = (prec(i,m))/(d1*rhow)
            dwg = dwg - (wtlwg(i,m)-wtlwb(i,m))/tau1	! drain/recharge
!
            wtlwg(i,m) = max(wtlwg(i,m) + ztodt*dwg, 0.)
!            wgxs(i,m)  = wgxs(i,1)*wtrc_ratio(prec(i,m),prec(i,1))
            wgxs(i,m)  = wgxs(i,1)*riso(i,m)
            wtlwg(i,m) = wtlwg(i,m) - wgxs(i,m)
!
            wtlwb(i,m) = max(wtlwb(i,m) + ztodt*dwb, 0.)
!            wbxs(i,m)  = wbxs(i,1)*wtrc_ratio(prec(i,m),prec(i,1))
            wbxs(i,m)  = wbxs(i,1)*riso(i,m)
            wtlwb(i,m) = wtlwb(i,m) - wbxs(i,m)
!
! If total has dried up apply the saved isotope ratio so that we have
! something sensible to evaporate next time
!
            if (wtlwb(i,1) <= wtiny) then
              wtlwb(i,m) = riso(i,m)*wtlwb(i,1)
            end if
!
            wtrun(i,m) = d1*rhow*wgxs(i,m)/ztodt
            wtdrn(i,m) = d2*rhow*wbxs(i,m)/ztodt
         end do

!
!!         if (m == 7) then
!!           write(*,*) 'I:',i,landfrac(i)
!!           write(*,*) '  P,  E:',prec(i,1),parm%cflx(i,1)
!!           write(*,*) 'DWB,DWG:',dwb*ztodt,dwg*ztodt
!!           write(*,*) ' WB, WG:',wtlwb(i,1), wtlwg(i,1)
!!           write(*,*) ' P:',prec(i,1) ,wiso_delta(iwspec(m),prec(i,m),prec(i,1))
!!           write(*,*) 'WG:',wtlwg(i,1),wiso_delta(iwspec(m),wtlwg(i,m),wtlwg(i,1))
!!         end if

         call outfld('LWG'//trim(cnst_name(m)), wtlwg(:,m)  ,pcols   ,lchnk )
         call outfld('LWB'//trim(cnst_name(m)), wtlwb(:,m)  ,pcols   ,lchnk )
         call outfld('LWR'//trim(cnst_name(m)), wtrun(:,m)  ,pcols   ,lchnk )
         call outfld('LWD'//trim(cnst_name(m)), wtdrn(:,m)  ,pcols   ,lchnk )
      end if
    end do
!
   return
end subroutine wtrc_land_drv


!=======================================================================

end module wtrc_land

