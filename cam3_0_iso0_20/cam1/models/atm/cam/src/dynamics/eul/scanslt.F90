#include <misc.h>
#include <params.h>

subroutine scanslt(ztodt   ,pmap    ,etadot  ,kdpmpf  ,kdpmph  , &
                   lam     ,phi     ,dphi    ,sinlam  ,coslam  , &
                   lbasdy  ,lbasdz  ,lbassd  ,lbasiy  ,detam   , &
                   detai   ,dlam    ,cwava   ,etamid  ,etaint  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Driving routine for semi-lagrangian transport.
! 
! Method: 
! The latitude loop in this routine is multitasked.
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch,  March 1996
!
!-----------------------------------------------------------------------
!
! $Id: scanslt.F90,v 1.9.2.4.10.1 2004/05/20 18:36:06 mvr Exp $
! $Author: mvr $
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,      only: plon, plond, plev, plevp, plat, platd, beglat, endlat, beglatex, &
                          endlatex, i1, j1
   use constituents,only: pcnst
   use comslt,      only: qfcst, lammp, phimp, sigmp, hw1lat, trcavg
   use prognostics, only: u3, v3, qminus, ps, n3m2, n3m1, q3, n3
   use rgrid,       only: nlon, nlonex
   use commap,      only: w
   use dynconst,    only: ra
   use time_manager, only: get_nstep

#if ( defined SCAM )
   use physconst, only: gravit
#endif

#if (defined SPMD)
# if ( defined TIMING_BARRIERS )
   use mpishorthand, only: mpicom
# endif
#endif
!-----------------------------------------------------------------------
   implicit none
!------------------------------Parameters-------------------------------
   integer itermx  ! number of iterations to be used in departure
!                     ! point calculation for nstep = 0 and 1
   integer itermn  ! number of iterations to be used in departure
!                     ! point calculation for nstep > 1
   parameter(itermx=4,itermn=1)
!-----------------------------------------------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: ztodt              ! twice the time step unless nstep = 0
!
   integer, intent(in) :: pmap                ! dimension of artificial array
!                                 !  used to locate vertical interval
!                                 !  in which departure point falls
   real(r8), intent(in) :: etadot(plon,plevp,beglat:endlat)! vertical motion (slt)
!
   integer, intent(in) :: kdpmpf(pmap)        ! mapping from artificial grid to
!                                 ! model levels
   integer, intent(in) :: kdpmph(pmap)        ! mapping from artificial grid to
!                                 ! interfaces
   real(r8), intent(in) :: lam(plond,platd)       ! longitude coordinates of model grid
   real(r8), intent(in) :: phi(platd)             ! latitude  coordinates of model grid
   real(r8), intent(in) :: dphi(platd)            ! latitudinal grid increments
   real(r8), intent(in) :: sinlam(plond,platd)    ! sine of longitude
   real(r8), intent(in) :: coslam(plond,platd)    ! cosine of longitude
   real(r8), intent(in) :: lbasdy(4,2,platd)      ! basis functions for lat deriv est.
   real(r8), intent(in) :: lbasdz(4,2,plev)       ! basis functions for vert deriv est.
!                                 !  (levels)
   real(r8), intent(in) :: lbassd(4,2,plevp)  ! basis functions for vert deriv est.
!                                 !  (levels)
   real(r8), intent(in) :: lbasiy(4,2,platd)  ! basis functions for Lagrange interp
   real(r8), intent(in) :: detam(plev)        ! delta eta at levels
   real(r8), intent(in) :: detai(plevp)       ! delta eta at interfaces
   real(r8), intent(in) :: dlam(platd)        ! longitudinal grid increment
   real(r8), intent(in) :: cwava(plat)        ! weight for global water vapor int.
   real(r8), intent(in) :: etamid(plev)       ! eta at levels
   real(r8), intent(in) :: etaint(plevp)      ! eta at interfaces
!
!---------------------------Local workspace-----------------------------
!
   integer iter                ! number of iterations for
!                                 ! departure point calculation
   integer m
   integer lat                 ! latitude index
   integer irow                ! N/S latitude pair index
   integer jcen                ! lat index (extended grid) of forecast
   integer :: nstep            ! current timestep number

   real(r8) pmid(plond,plev)            ! pressure at model levels
   real(r8) pint(plond,plevp)           ! pressure at interfaces
   real(r8) pdel(plond,plev)            ! pressure difference between
!
! Dynamic (SPMD) vs stack (shared memory)
!
   real(r8) uxl(plond,plev,beglatex:endlatex)     ! left  x-deriv of u/v
   real(r8) uxr(plond,plev,beglatex:endlatex)     ! left  x-deriv of u/v
   real(r8) vxl(plond,plev,beglatex:endlatex)     ! left  x-deriv of u/v
   real(r8) vxr(plond,plev,beglatex:endlatex)     ! left  x-deriv of u/v
   real(r8) qxl(plond,plev,pcnst,beglatex:endlatex) ! left  x-deriv of constituents
   real(r8) qxr(plond,plev,pcnst,beglatex:endlatex) ! right  x-deriv of constituents
#if ( defined SCAM )
   real(r8) gw(plat)                ! Gaussian weights
   integer k
#endif
!
!-----------------------------------------------------------------------
#if ( !defined SCAM )
#if ( defined SPMD )
#ifdef TIMING_BARRIERS
   call t_startf ('sync_bndexch')
   call mpibarrier (mpicom)
   call t_stopf ('sync_bndexch')
#endif
!
! Communicate boundary information 
!
   call t_startf ('bndexch')
   call bndexch
   call t_stopf ('bndexch')
#endif

   nstep = get_nstep()
!
! Initialize extended arrays
!
   call t_startf('sltini')
   call sltini (dlam,    sinlam,  coslam,  uxl,     uxr, &
                vxl,     vxr,     qxl,     qxr,     u3,  &
                v3,      qminus,  n3m1)
   call t_stopf('sltini')
#else
!************************************************************************
!     IF surface pressure changes with time we need to remap the vertical
!     coordinate for the slt advection process.  It has been empirically
!     determined that we can get away with 500 for pmap (instead of 20000)
!     This is necessary to make the procedure computationally feasible
!
!
!   call grdini(pmap    ,etamid    ,etaint    ,gravit  ,dlam    , &
!        lam     ,phi     ,dphi    ,gw      ,sinlam  , &
!        coslam  ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
!        detam   ,detai   ,kdpmpf  ,kdpmph  ,cwava   )
      
!   ps_prev = ps
!     
! Initial guess for trajectory midpoints in spherical coords.
! nstep = 0:  use arrival points as initial guess for trajectory midpoints.
! nstep > 0:  use calculated trajectory midpoints from previous time
! step as first guess.
! NOTE:  reduce number of iters necessary for convergence after nstep = 1.
!
   do k=1,plev
      sigmp(1,k,beglat) = etamid(k)
   end do

   nstep = get_nstep()
#endif
   if (nstep .le. 1) then
      iter = itermx
   else
      iter = itermn
   end if
!
! Initialize moisture mass integrals.
!
   do m=1,pcnst
      do lat=1,plat
         hw1lat(m,lat) = 0.0
      end do
   end do
!
! Loop through latitudes producing forecast
!
!$OMP PARALLEL DO PRIVATE (LAT, IROW, JCEN, PINT, PMID, PDEL )

   do lat=beglat,endlat
      if(lat.le.plat/2) then
         irow = lat
      else
         irow = plat + 1 - lat
      end if
      jcen = j1 - 1 + lat
!
! Only pdel is needed inside SLT.  pint and pmid are not.
!
      call plevs0 (nlon(lat),plond,plev,ps(1,lat,n3m2), pint, pmid, pdel)
!
! Calculate mass of moisture in field being advected by slt.
!

!  q3     is plond,plev,pcnst+pnats,beglattex:endlatex,ptimelevs
!  qminus is plond,plev,pcnst+pnats,beglattex:endlatex
!  qfcst  is plond,plev,pcnst,beglat:endlat
      call qmassa (cwava(lat),w(irow) ,qminus(i1,1,1,jcen),pdel    , &
                   hw1lat(1,lat),nlon(lat), q3(i1,1,1,jcen,n3m2), lat)
!
! Call slt interface routine.
!
      call sltb1 (pmap    ,jcen    ,lat     ,ztodt   ,ra      , &
                  iter    ,uxl     ,uxr     ,vxl     ,vxr     , &
                  etadot(1,1,lat)  ,qxl     ,qxr     ,lam     , &
                  phi     ,dphi    ,etamid  ,etaint  ,detam   , &
                  detai   ,lbasdy  ,lbasdz  ,lbassd  ,lbasiy  , &
                  kdpmpf  ,kdpmph  ,lammp(1,1,lat), phimp(1,1,lat), sigmp(1,1,lat), &
                  qfcst(1,1,1,lat), u3      ,v3     ,qminus, n3m1, &
                  nlon(lat), nlonex  )
   end do
   return
end subroutine scanslt
