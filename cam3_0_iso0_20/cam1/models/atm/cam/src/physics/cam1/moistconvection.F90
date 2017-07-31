#include <misc.h>
#include <params.h>
module moistconvection
!
! Moist convection. Primarily data used by both Zhang-McFarlane convection
! and Hack shallow convective schemes.
!
! $Id: moistconvection.F90,v 1.1.4.14 2004/05/10 21:39:47 jmccaa Exp $
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use abortutils, only: endrun
   implicit none

   private
   save
!
! Public interfaces
!
   public mfinti   !  Initialization of data for moist convection
   public cmfmca   !  Hack shallow convection
!
! Public Data for moist convection
!
   real(r8), public :: cp          ! specific heat of dry air
   real(r8), public :: grav        ! gravitational constant       
   real(r8), public :: rgrav       ! reciprocal of grav
   real(r8), public :: rgas        ! gas constant for dry air
   integer, public :: limcnv       ! top interface level limit for convection
!
! Private data used for Hack shallow convection
!
   real(r8) :: hlat        ! latent heat of vaporization
   real(r8) :: c0          ! rain water autoconversion coefficient
   real(r8) :: betamn      ! minimum overshoot parameter
   real(r8) :: rhlat       ! reciprocal of hlat
   real(r8) :: rcp         ! reciprocal of cp
   real(r8) :: cmftau      ! characteristic adjustment time scale
   real(r8) :: rhoh2o      ! density of liquid water (STP)
   real(r8) :: dzmin       ! minimum convective depth for precipitation
   real(r8) :: tiny        ! arbitrary small num used in transport estimates
   real(r8) :: eps         ! convergence criteria (machine dependent)
   real(r8) :: tpmax       ! maximum acceptable t perturbation (degrees C)
   real(r8) :: shpmax      ! maximum acceptable q perturbation (g/g)           

   integer :: iloc         ! longitude location for diagnostics
   integer :: jloc         ! latitude  location for diagnostics
   integer :: nsloc        ! nstep for which to produce diagnostics
!
   logical :: rlxclm       ! logical to relax column versus cloud triplet
   logical :: liceliq	   ! flag export ice as well as liquid

contains

subroutine mfinti (rair    ,cpair   ,gravit  ,latvap  ,rhowtr  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialize moist convective mass flux procedure common block, cmfmca
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J. Hack
! 
!-----------------------------------------------------------------------
   use pmgrid, only: plev, plevp, masterproc
   use dycore, only: dycore_is, get_resolution
#include <comhyb.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: rair              ! gas constant for dry air
   real(r8), intent(in) :: cpair             ! specific heat of dry air
   real(r8), intent(in) :: gravit            ! acceleration due to gravity
   real(r8), intent(in) :: latvap            ! latent heat of vaporization
   real(r8), intent(in) :: rhowtr            ! density of liquid water (STP)

   integer k              ! vertical level index
!
!-----------------------------------------------------------------------
!
! Turn on/off export of ice and snow, as well as liquid and rain
!
  liceliq = .true.
!
! Initialize physical constants for moist convective mass flux procedure
!
   cp     = cpair         ! specific heat of dry air
   hlat   = latvap        ! latent heat of vaporization
   grav   = gravit        ! gravitational constant
   rgas   = rair          ! gas constant for dry air
   rhoh2o = rhowtr        ! density of liquid water (STP)
!
! Initialize free parameters for moist convective mass flux procedure
!
   if (dycore_is('LR'))then
      if ( get_resolution() == '4x5' ) then
         cmftau = 1800.         ! characteristic adjustment time scale
         c0     = 2.0e-4        ! rain water autoconversion coeff (1/m)
      else
         cmftau = 1800.         ! characteristic adjustment time scale
         c0     = 1.0e-4        ! rain water autoconversion coeff (1/m)
      endif
   else
      if(get_resolution() == 'T85')then
         cmftau = 1800.         ! characteristic adjustment time scale
         c0     = 1.0e-4        ! rain water autoconversion coeff (1/m)
      elseif(get_resolution() == 'T31')then
         cmftau = 1800.         ! characteristic adjustment time scale
         c0     = 5.0e-4        ! rain water autoconversion coeff (1/m)
      else
         cmftau = 1800.         ! characteristic adjustment time scale
         c0     = 2.0e-4        ! rain water autoconversion coeff (1/m)
      endif
   endif
   dzmin  = 0.0           ! minimum cloud depth to precipitate (m)
   betamn = 0.10          ! minimum overshoot parameter
!
! Limit convection to regions below 40 mb
!
   if (hypi(1) >= 4.e3) then
      limcnv = 1
   else
      do k=1,plev
         if (hypi(k) < 4.e3 .and. hypi(k+1) >= 4.e3) then
            limcnv = k
            goto 10
         end if
      end do
      limcnv = plevp
   end if

10 continue

   if (masterproc) then
      write(6,*)'MFINTI: Convection will be capped at intfc ',limcnv, &
                ' which is ',hypi(limcnv),' pascals'
   end if

   tpmax  = 1.50          ! maximum acceptable t perturbation (deg C)
   shpmax = 1.50e-3       ! maximum acceptable q perturbation (g/g)
   rlxclm = .true.        ! logical variable to specify that relaxation
!                                time scale should applied to column as
!                                opposed to triplets individually
!
! Initialize miscellaneous (frequently used) constants
!
   rhlat  = 1.0/hlat      ! reciprocal latent heat of vaporization
   rcp    = 1.0/cp        ! reciprocal specific heat of dry air
   rgrav  = 1.0/grav      ! reciprocal gravitational constant
!
! Initialize diagnostic location information for moist convection scheme
!
   iloc   = 1             ! longitude point for diagnostic info
   jloc   = 1             ! latitude  point for diagnostic info
   nsloc  = 1             ! nstep value at which to begin diagnostics
!
! Initialize other miscellaneous parameters
!
   tiny   = 1.0e-36       ! arbitrary small number (scalar transport)
   eps    = 1.0e-13       ! convergence criteria (machine dependent)
!
   return
end subroutine mfinti

subroutine cmfmca(lchnk   ,ncol    , &
                  nstep   ,ztodt     ,pmid    ,pdel    , &
                  rpdel   ,zm      ,tpert   ,qpert   ,phis    , &
                  pblht   ,t       ,q       ,cmfdt   ,dq      , &
                  cmfmc   ,cmfdqr  ,cmfsl   ,cmflq   ,precc   , &
                  qcl     ,qci     ,cnt     ,cnb     ,icwmr   ,rliq    , rice , & 
                  pmiddry ,pdeldry ,rpdeldry,trprecc ,trqcl   ,trqci   , &
                  trdqr   ,lwtrccmf)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Moist convective mass flux procedure:
! 
! Method: 
! If stratification is unstable to nonentraining parcel ascent,
! complete an adjustment making successive use of a simple cloud model
! consisting of three layers (sometimes referred to as a triplet)
!
! Code generalized to allow specification of parcel ("updraft")
! properties, as well as convective transport of an arbitrary
! number of passive constituents (see q array).  The code
! is written so the water vapor field is passed independently
! in the calling list from the block of other transported
! constituents, even though as currently designed, it is the
! first component in the constituents field.
! 
! Author: J. Hack
!
! BAB: changed code to report tendencies in cmfdt and dq, instead of
! updating profiles. Cmfdq contains water only, made it a local variable
! made dq (all constituents) the argument.
!
! Modified for water tracers.
!  David Noone <dcn@colorado.edu> - Mon Jul  5 12:11:54 MDT 2004
! 
!-----------------------------------------------------------------------

!#######################################################################
!#                                                                     #
!# Debugging blocks are marked this way for easy identification        #
!#                                                                     #
!#######################################################################
   use constituents,  only: pcnst, pnats
   use constituents,    only: cnst_get_type_byind
   use ppgrid,    only: pcols, pver, pverp
   use phys_grid, only: get_lat_all_p, get_lon_all_p
   use wv_saturation, only: aqsatd, vqsatd
   use water_tracers, only: trace_water, wtrc_is_vap, iwspec, wtrc_ratio
   use water_isotopes, only: wiso_dicm
   use wtrc_camtune,   only: feq_cmf, fsnk_cmf
   use cldwat,         only: cldwat_fice
   use physconst,      only: latice

   real(r8) ssfac               ! supersaturation bound (detrained air)
   parameter (ssfac = 1.001)

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in) :: ncol                 ! number of atmospheric columns
   integer, intent(in) :: nstep                ! current time step index

   real(r8), intent(in) :: ztodt               ! 2 delta-t (seconds)
   real(r8), intent(in) :: pmid(pcols,pver)    ! pressure
   real(r8), intent(in) :: pdel(pcols,pver)    ! delta-p
   real(r8), intent(in) :: pmiddry(pcols,pver)    ! pressure
   real(r8), intent(in) :: pdeldry(pcols,pver)    ! delta-p
   real(r8), intent(in) :: rpdel(pcols,pver)   ! 1./pdel
   real(r8), intent(in) :: rpdeldry(pcols,pver)   ! 1./pdel
   real(r8), intent(in) :: zm(pcols,pver)      ! height abv sfc at midpoints
   real(r8), intent(in) :: tpert(pcols)        ! PBL perturbation theta
   real(r8), intent(in) :: qpert(pcols,pcnst+pnats)  ! PBL perturbation specific humidity
   real(r8), intent(in) :: phis(pcols)         ! surface geopotential
   real(r8), intent(in) :: pblht(pcols)        ! PBL height (provided by PBL routine)
   real(r8), intent(in) :: t(pcols,pver)       ! temperature (t bar)
   real(r8), intent(in) :: q(pcols,pver,pcnst+pnats) ! specific humidity (sh bar)
!
! Output arguments
!
   real(r8), intent(out) :: cmfdt(pcols,pver)   ! dt/dt due to moist convection
   real(r8), intent(out) :: cmfmc(pcols,pverp)  ! moist convection cloud mass flux
   real(r8), intent(out) :: cmfdqr(pcols,pver)  ! dq/dt due to convective rainout
   real(r8), intent(out) :: cmfsl(pcols,pver )  ! convective lw static energy flux
   real(r8), intent(out) :: cmflq(pcols,pver )  ! convective total water flux
   real(r8), intent(out) :: precc(pcols)        ! convective precipitation rate
! JJH mod to explicitly export cloud water
   real(r8), intent(out) :: qcl(pcols,pver)      ! dq/dt due to export of cloud water liquid
   real(r8), intent(out) :: qci(pcols,pver)      ! dq/dt due to export of cloud water ice
   real(r8), intent(out) :: cnt(pcols)          ! top level of convective activity
   real(r8), intent(out) :: cnb(pcols)          ! bottom level of convective activity
   real(r8), intent(out) :: dq(pcols,pver,pcnst+pnats) ! constituent tendencies
   real(r8), intent(out) :: icwmr(pcols,pver)
   real(r8), intent(out) :: rliq(pcols) 
   real(r8), intent(out) :: rice(pcols) 

   real(r8), intent(out) :: trprecc(pcols,pcnst+pnats)     ! tracer precipitation rate
   real(r8), intent(out) :: trqcl(pcols,pver,pcnst+pnats)   ! dtr/dq due to export (like qcl)
   real(r8), intent(out) :: trqci(pcols,pver,pcnst+pnats)   ! dtr/dq due to export (like qci)
   real(r8), intent(out) :: trdqr(pcols,pver,pcnst+pnats)  ! net tracer rainout tendency (like cmfdqr)

   logical , intent(in)  :: lwtrccmf		! flag to do water tracers
!
!---------------------------Local workspace-----------------------------
!
   real(r8) pm(pcols,pver)    ! pressure
   real(r8) pd(pcols,pver)    ! delta-p
   real(r8) rpd(pcols,pver)   ! 1./pdel

   real(r8) cmfdq(pcols,pver)   ! dq/dt due to moist convection
   real(r8) gam(pcols,pver)     ! 1/cp (d(qsat)/dT)
   real(r8) sb(pcols,pver)      ! dry static energy (s bar)
   real(r8) hb(pcols,pver)      ! moist static energy (h bar)
   real(r8) shbs(pcols,pver)    ! sat. specific humidity (sh bar star)
   real(r8) hbs(pcols,pver)     ! sat. moist static energy (h bar star)
   real(r8) shbh(pcols,pverp)   ! specific humidity on interfaces
   real(r8) sbh(pcols,pverp)    ! s bar on interfaces
   real(r8) hbh(pcols,pverp)    ! h bar on interfaces
   real(r8) cmrh(pcols,pverp)   ! interface constituent mixing ratio
   real(r8) prec(pcols)         ! instantaneous total precipitation
   real(r8) dzcld(pcols)        ! depth of convective layer (m)
   real(r8) beta(pcols)         ! overshoot parameter (fraction)
   real(r8) betamx(pcols)       ! local maximum on overshoot
   real(r8) eta(pcols)          ! convective mass flux (kg/m^2 s)
   real(r8) etagdt(pcols)       ! eta*grav*dt
   real(r8) cldwtr(pcols)       ! cloud water (mass)
   real(r8) rnwtr(pcols)        ! rain water  (mass)
   real(r8) snwtr(pcols)        ! snow water  (mass)
   real(r8) dtwtr(pcols)        ! detraining water  (mass)
!  JJH extension to facilitate export of cloud liquid water
   real(r8) totcond(pcols)	! total condensate; mix of precip and cloud water (mass)
   real(r8) sc  (pcols)         ! dry static energy   ("in-cloud")
   real(r8) shc (pcols)         ! specific humidity   ("in-cloud")
   real(r8) hc  (pcols)         ! moist static energy ("in-cloud")
   real(r8) cmrc(pcols)         ! constituent mix rat ("in-cloud")
   real(r8) dq1(pcols)          ! shb  convective change (lower lvl)
   real(r8) dq2(pcols)          ! shb  convective change (mid level)
   real(r8) dq3(pcols)          ! shb  convective change (upper lvl)
   real(r8) ds1(pcols)          ! sb   convective change (lower lvl)
   real(r8) ds2(pcols)          ! sb   convective change (mid level)
   real(r8) ds3(pcols)          ! sb   convective change (upper lvl)
   real(r8) dcmr1(pcols)        ! q convective change (lower lvl)
   real(r8) dcmr2(pcols)        ! q convective change (mid level)
   real(r8) dcmr3(pcols)        ! q convective change (upper lvl)
   real(r8) estemp(pcols,pver)  ! saturation vapor pressure (scratch)
   real(r8) vtemp1(2*pcols)     ! intermediate scratch vector
   real(r8) vtemp2(2*pcols)     ! intermediate scratch vector
   real(r8) vtemp3(2*pcols)     ! intermediate scratch vector
   real(r8) vtemp4(2*pcols)     ! intermediate scratch vector
   integer indx1(pcols)     ! longitude indices for condition true
   logical etagt0           ! true if eta > 0.0
   real(r8) sh1                 ! dummy arg in qhalf statement func.
   real(r8) sh2                 ! dummy arg in qhalf statement func.
   real(r8) shbs1               ! dummy arg in qhalf statement func.
   real(r8) shbs2               ! dummy arg in qhalf statement func.
   real(r8) cats                ! modified characteristic adj. time
   real(r8) rtdt                ! 1./ztodt
   real(r8) qprime              ! modified specific humidity pert.
   real(r8) tprime              ! modified thermal perturbation
   real(r8) pblhgt              ! bounded pbl height (max[pblh,1m])
   real(r8) fac1                ! intermediate scratch variable
   real(r8) shprme              ! intermediate specific humidity pert.
   real(r8) qsattp              ! sat mix rat for thermally pert PBL parcels
   real(r8) dz                  ! local layer depth
   real(r8) temp1               ! intermediate scratch variable
   real(r8) b1                  ! bouyancy measure in detrainment lvl
   real(r8) b2                  ! bouyancy measure in condensation lvl
   real(r8) temp2               ! intermediate scratch variable
   real(r8) temp3               ! intermediate scratch variable
   real(r8) g                   ! bounded vertical gradient of hb
   real(r8) tmass               ! total mass available for convective exch
   real(r8) denom               ! intermediate scratch variable
   real(r8) qtest1              ! used in negative q test (middle lvl)
   real(r8) qtest2              ! used in negative q test (lower lvl)
   real(r8) fslkp               ! flux lw static energy (bot interface)
   real(r8) fslkm               ! flux lw static energy (top interface)
   real(r8) fqlkp               ! flux total water (bottom interface)
   real(r8) fqlkm               ! flux total water (top interface)
   real(r8) botflx              ! bottom constituent mixing ratio flux
   real(r8) topflx              ! top constituent mixing ratio flux
   real(r8) efac1               ! ratio q to convectively induced chg (btm lvl)
   real(r8) efac2               ! ratio q to convectively induced chg (mid lvl)
   real(r8) efac3               ! ratio q to convectively induced chg (top lvl)
   real(r8) tb(pcols,pver)      ! working storage for temp (t bar)
   real(r8) shb(pcols,pver)     ! working storage for spec hum (sh bar)
   real(r8) adjfac              ! adjustment factor (relaxation related)
   real(r8) rktp
   real(r8) rk

! Additional variables for water tracers
! (some of these are really superfluous, except they clarify things for me)
   real(r8) trb(pcols,pver,pcnst+pnats)         ! working storage for tracer spec hum (sh bar)
   real(r8) trdq(pcols,pver,pcnst+pnats)  	! net tracer rainout (like cmfdq)
!
   real(r8) iuwmr(pcols,pver)                   ! in-cloud updraft water mixing ratio
   real(r8) rnwmr(pcols,pver)                   ! precipitation water mixing ratio
   real(r8) rtlmr(pcols,pver)                   ! retained water mixing ratio
   real(r8) dtlmr(pcols,pver)                   ! detrained water mixing ratio
   real(r8) cndmr(pcols,pver)                   ! condenation water mixing ratio
   real(r8) truwmr(pcols,pver,pcnst+pnats)      ! tracer updraft water mixing ratio
   real(r8) trcwmr(pcols,pver,pcnst+pnats)      ! tracer CONDENSED (liq+ice) water mixing ratio 
   real(r8) trrnwmr(pcols,pver,pcnst+pnats)     ! tracer rain water mixing ratio
   real(r8) trsnwmr(pcols,pver,pcnst+pnats)     ! tracer snow water mixing ratio
   real(r8) trrtlmr(pcols,pver,pcnst+pnats)     ! tracer retained liquid mixing ratio
   real(r8) trrtimr(pcols,pver,pcnst+pnats)     ! tracer retained ice mixing ratio
   real(r8) trdtlmr(pcols,pver,pcnst+pnats)     ! tracer detrained liquid mixing ratio
   real(r8) trdtimr(pcols,pver,pcnst+pnats)     ! tracer detrained ice mixing ratio
   real(r8) trcndmr(pcols,pver,pcnst+pnats)     ! tracer condenation water mixing ratio
   real(r8) fprt(pcols,pver)                    ! PBL vapour perturbation factor
   real(r8) trdtlq(pcols)                       ! tracer detraining liquid (like (1-fice)*dtwtr)
   real(r8) trdtic(pcols)                       ! tracer detraining ice (like fice*dtwtr)
   real(r8) trrain(pcols)                       ! tracer rain water (like (1-fsnow)*rnwtr)
   real(r8) trsnow(pcols)                       ! tracer snow water (like    fsnow *rnwtr)
   real(r8) trcond(pcols)                       ! tracer condensed (like totcond)
   real(r8) tcld(pcols,pver)
   real(r8) tbar                                ! mean fractionation temperature
   real(r8) alpliq,alpice                       ! fractionation coefficients
   real(r8) fice(pcols,pver)                    ! fracition ice vs liquid formed
   real(r8) fsnow(pcols,pver)                   ! fracition snow vs rain formed


   real(r8) trprec(pcols,pcnst+pnats) 		! tracer precipitation
   real(r8) netliq, trnetliq                    ! liquid available for fractionation
   real(r8) netice, trnetice                    ! ice available for fractionation
   real(r8) ntlsnk, trntlsnk                    ! total liquid sink (rain + liq detrain)
   real(r8) ntisnk, trntisnk                    ! total ice sink (snow + ice detrain)
   real(r8) fntld, fntid                        ! fraction of net liquid, ice detrained
   real(r8) flskd, fiskd                        ! fraction of liquid, ice sink detrained
! 
   real(r8) qliq, qice				! liquid and ice wmr
   real(r8) qrnw, qsnw				! rain and snow wmr
   real(r8) ldet, idet				! detraining liquid and ice wmr
!
! Variables passed in/out of fractionation sub-model
!
   integer, parameter :: piso = 2       ! solve for one species at a time
   integer  :: isp(piso)                        ! species indicies (1 is water, 2 is isotope)
   real(r8) :: told, tnew			! old and new temperature
   real(r8) :: vapold(piso)                     ! initial vapour
   real(r8) :: liqold(piso)                     ! initial liquid
   real(r8) :: iceold(piso)                     ! initial ice
   real(r8) :: vapent(piso)                     ! vapour entrainment
   real(r8) :: vapnew(piso)                     ! final vapour
   real(r8) :: liqnew(piso)                     ! final liquid
   real(r8) :: icenew(piso)                     ! final ice
   real(r8) :: vapdet(piso)                     ! vapour detrainment
   real(r8) :: rainpr(piso)                     ! rain production
   real(r8) :: snowpr(piso)                     ! snow production
   real(r8) :: dliqmt                           ! liquid change due to ice melt
   real(r8) :: dicefz                           ! ice change due to liquid freeze

#if ( defined DIAGNS )
!
!  Following 7 real variables are used in diagnostics calculations
!
   real(r8) rh                  ! relative humidity
   real(r8) es                  ! sat vapor pressure
   real(r8) hsum1               ! moist static energy integral
   real(r8) qsum1               ! total water integral
   real(r8) hsum2               ! final moist static energy integral
   real(r8) qsum2               ! final total water integral
   real(r8) fac                 ! intermediate scratch variable
#endif
   integer i,k              ! longitude, level indices
   integer ii               ! index on "gathered" vectors
   integer len1             ! vector length of "gathered" vectors
   integer m                ! constituent index
   integer ktp              ! tmp indx used to track top of convective layer
#if ( defined DIAGNS )
   integer n                ! vertical index     (diagnostics)
   integer kp               ! vertical index     (diagnostics)
   integer kpp              ! index offset, kp+1 (diagnostics)
   integer kpm1             ! index offset, kp-1 (diagnostics)
   integer lat(pcols)       ! latitude indices
   integer lon(pcols)       ! longitude indices
#endif
!
!---------------------------Statement functions-------------------------
!
   real(r8) qhalf, qhalf_wtrc
   qhalf(sh1,sh2,shbs1,shbs2) = min(max(sh1,sh2),(shbs2*sh1 + shbs1*sh2)/(shbs1+shbs2))
   qhalf_wtrc(sh1,sh2) = 0.5*(sh1 + sh2)		! nicer for water tracers
!
!-----------------------------------------------------------------------

!** BAB initialize output tendencies here
!       copy q to dq; use dq below for passive tracer transport
   cmfdt(:ncol,:)  = 0.
   cmfdq(:ncol,:)  = 0.
   dq(:ncol,:,2:)  = q(:ncol,:,2:)
   cmfmc(:ncol,:)  = 0.
   cmfdqr(:ncol,:) = 0.
   cmfsl(:ncol,:)  = 0.
   cmflq(:ncol,:)  = 0.
   qcl(:ncol,:)     = 0.
   qci(:ncol,:)     = 0.
   rliq(:ncol)     = 0.
   rice(:ncol)     = 0.
!
! initialize water tracer arrays
!
   fice(:ncol,:)   = 0.
   fsnow(:ncol,:)  = 0.
   fprt(:ncol,:)   = 0.
   trdq (:ncol,:,:)= 0.
   trdqr(:ncol,:,:)= 0.
   trqcl(:ncol,:,:) = 0.
   trqci(:ncol,:,:) = 0.
   trprec(:ncol,:) = 0.
!
#if ( defined DIAGNS )
! Determine chunk latitudes and longitudes
   call get_lat_all_p(lchnk, ncol, lat)
   call get_lon_all_p(lchnk, ncol, lon)
#endif
!
! Ensure that characteristic adjustment time scale (cmftau) assumed
! in estimate of eta isn't smaller than model time scale (ztodt)
! The time over which the convection is assumed to act (the adjustment
! time scale) can be applied with each application of the three-level
! cloud model, or applied to the column tendencies after a "hard"
! adjustment (i.e., on a 2-delta t time scale) is evaluated
!
   if (rlxclm) then
      cats   = ztodt             ! relaxation applied to column
      adjfac = ztodt/(max(ztodt,cmftau))
   else
      cats   = max(ztodt,cmftau) ! relaxation applied to triplet
      adjfac = 1.0
   endif
   rtdt = 1.0/ztodt
!
! Move temperature and moisture into working storage
!
   do k=limcnv,pver
      do i=1,ncol
         tb (i,k) = t(i,k)
         shb(i,k) = q(i,k,1)
         tcld (i,k) = t(i,k)
      end do
   end do
   do m = 1, pcnst+pnats
     do k=limcnv,pver
        do i=1,ncol
           trb(i,k,m) = q(i,k,m)
        end do
     end do
   end do
   do k=1,pver
      do i=1,ncol
         icwmr(i,k) = 0.
         iuwmr(i,k) = 0.
      end do
   end do
   do m = 1, pcnst+pnats
     do k=1,pver
        do i=1,ncol
           trcwmr(i,k,m) = 0.
           truwmr(i,k,m) = 0.
        end do
     end do
   end do
!
! Compute sb,hb,shbs,hbs
!
   call aqsatd(tb      ,pmid    ,estemp ,shbs    ,gam     , &
               pcols   ,ncol    ,pver   ,limcnv  ,pver    )
!
   do k=limcnv,pver
      do i=1,ncol
         sb (i,k) = cp*tb(i,k) + zm(i,k)*grav + phis(i)
         hb (i,k) = sb(i,k) + hlat*shb(i,k)
         hbs(i,k) = sb(i,k) + hlat*shbs(i,k)
      end do
   end do
!
! Compute sbh, shbh
!
   do k=limcnv+1,pver
      do i=1,ncol
         sbh (i,k) = 0.5*(sb(i,k-1) + sb(i,k))
         if (trace_water) then
           shbh(i,k) = qhalf_wtrc(shb(i,k-1),shb(i,k))
         else 
           shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
         end if
         hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
      end do
   end do
!
! Specify properties at top of model (not used, but filling anyway)
!
   do i=1,ncol
      sbh (i,limcnv) = sb(i,limcnv)
      shbh(i,limcnv) = shb(i,limcnv)
      hbh (i,limcnv) = hb(i,limcnv)
   end do
!
! Zero vertically independent control, tendency & diagnostic arrays
!
   do i=1,ncol
      prec(i)  = 0.0
      dzcld(i) = 0.0
      cnb(i)   = 0.0
      cnt(i)   = float(pver+1)
   end do
#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
!#    output initial thermodynamic profile if debug diagnostics        #
!#                                                                     #
   do i=1,ncol
     if ((lat(i).eq.jloc) .and. (lon(i).eq.iloc) &
         .and. (nstep.ge.nsloc)) then
!#                                                                     #
!#       approximate vertical integral of moist static energy          #
!#       and total preciptable water                                   #
!#                                                                     #
      hsum1 = 0.0
      qsum1 = 0.0
      do k=limcnv,pver
         hsum1 = hsum1 + pdel(i,k)*rgrav*hb(i,k)
         qsum1 = qsum1 + pdel(i,k)*rgrav*shb(i,k)
      end do
!#                                                                     #
      write (6,8010)
      fac = grav*864.
      do k=limcnv,pver
         rh = shb(i,k)/shbs(i,k)
         write(6,8020) shbh(i,k),sbh(i,k),hbh(i,k),fac*cmfmc(i,k),cmfsl(i,k), cmflq(i,k)
         write(6,8040) tb(i,k),shb(i,k),rh,sb(i,k),hb(i,k),hbs(i,k),ztodt*cmfdt(i,k), &
                       ztodt*cmfdq(i,k),ztodt*cmfdqr(i,k)
      end do
      write(6, 8000) prec(i)
     end if
   enddo
#endif
!#                                                                     #
!#                                                                     #
!#######################################################################
!
! Begin moist convective mass flux adjustment procedure.
! Formalism ensures that negative cloud liquid water can never occur
!
   do 70 k=pver-1,limcnv+1,-1
      do 10 i=1,ncol
         etagdt(i) = 0.0
         eta   (i) = 0.0
         beta  (i) = 0.0
         ds1   (i) = 0.0
         ds2   (i) = 0.0
         ds3   (i) = 0.0
         dq1   (i) = 0.0
         dq2   (i) = 0.0
         dq3   (i) = 0.0
!
! Specification of "cloud base" conditions
!
         qprime    = 0.0
         tprime    = 0.0
!
! Assign tprime within the PBL to be proportional to the quantity
! tpert (which will be bounded by tpmax), passed to this routine by
! the PBL routine.  Don't allow perturbation to produce a dry
! adiabatically unstable parcel.  Assign qprime within the PBL to be
! an appropriately modified value of the quantity qpert (which will be
! bounded by shpmax) passed to this routine by the PBL routine.  The
! quantity qprime should be less than the local saturation value
! (qsattp=qsat[t+tprime,p]).  In both cases, tpert and qpert are
! linearly reduced toward zero as the PBL top is approached.
!
         pblhgt = max(pblht(i),1.0_r8)
         if ( (zm(i,k+1) <= pblhgt) .and. dzcld(i) == 0.0 ) then
            fac1   = max(0.0_r8,1.0-zm(i,k+1)/pblhgt)
            tprime = min(tpert(i),tpmax)*fac1
            qsattp = shbs(i,k+1) + cp*rhlat*gam(i,k+1)*tprime
            shprme = min(min(qpert(i,1),shpmax)*fac1,max(qsattp-shb(i,k+1),0.0_r8))
            qprime = max(qprime,shprme)

            if (abs(qpert(i,1)) > 1.e-20) fprt(i,k) = qprime/qpert(i,1)

         else
            tprime = 0.0
            qprime = 0.0
         end if
!
! Specify "updraft" (in-cloud) thermodynamic properties
!
         sc (i)    = sb (i,k+1) + cp*tprime
         tcld(i,k) = tb (i,k+1) + tprime
         shc(i)    = shb(i,k+1) + qprime
         hc (i)    = sc (i    ) + hlat*shc(i)
         vtemp4(i) = hc(i) - hbs(i,k)
         dz        = pdel(i,k)*rgas*tb(i,k)*rgrav/pmid(i,k)
         if (vtemp4(i) > 0.0) then
            dzcld(i) = dzcld(i) + dz
         else
            dzcld(i) = 0.0
         end if
10       continue
#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
!#    output thermodynamic perturbation information                    #
!#                                                                     #
         do i=1,ncol
           if ((lat(i)==jloc).and.(lon(i)==iloc).and.(nstep>=nsloc)) then
            write (6,8090) k+1,sc(iloc),shc(iloc),hc(iloc)
           end if
         enddo
!#                                                                     #
!#######################################################################
#endif
!
! Compute the ice and snow fractions from the cloud base temperature
!
         if (liceliq) then
           call cldwat_fice(ncol, tcld, fice, fsnow)
           fsnow(:,k) = fice(:,k) 	! expell from cloud at fice, 
                                        ! evap repartitions to true fsnow
         else
           fice (:,k) = 0.
           fsnow(:,k) = 0.
         end if
!
! Check on moist convective instability
! Build index vector of points where instability exists
!
         len1 = 0
         do i=1,ncol
            if (vtemp4(i) > 0.0) then
               len1 = len1 + 1
               indx1(len1) = i
            end if
         end do
         if (len1 <= 0) go to 70
!
! Current level just below top level => no overshoot
!
         if (k <= limcnv+1) then
            do ii=1,len1
               i = indx1(ii)
               temp1     = vtemp4(i)/(1.0 + gam(i,k))
               cldwtr(i) = max(0.0_r8,(sb(i,k) - sc(i) + temp1))
               beta(i)   = 0.0
               vtemp3(i) = (1.0 + gam(i,k))*(sc(i) - sbh(i,k))
            end do
         else
!
! First guess at overshoot parameter using crude buoyancy closure
! 10% overshoot assumed as a minimum and 1-c0*dz maximum to start
! If pre-existing supersaturation in detrainment layer, beta=0
! cldwtr is temporarily equal to hlat*l (l=> liquid water)
!
!cdir nodep
!DIR$ CONCURRENT
            do ii=1,len1
               i = indx1(ii)
               temp1     = vtemp4(i)/(1.0 + gam(i,k))
               cldwtr(i) = max(0.0_r8,(sb(i,k)-sc(i)+temp1))
               betamx(i) = 1.0 - c0*max(0.0_r8,(dzcld(i)-dzmin))
               b1        = (hc(i) - hbs(i,k-1))*pdel(i,k-1)
               b2        = (hc(i) - hbs(i,k  ))*pdel(i,k  )
               beta(i)   = max(betamn,min(betamx(i), 1.0 + b1/b2))
               if (hbs(i,k-1) <= hb(i,k-1)) beta(i) = 0.0
!
! Bound maximum beta to ensure physically realistic solutions
!
! First check constrains beta so that eta remains positive
! (assuming that eta is already positive for beta equal zero)
!
               vtemp1(i) = -(hbh(i,k+1) - hc(i))*pdel(i,k)*rpdel(i,k+1)+ &
                           (1.0 + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i))
               vtemp2(i) = (1.0 + gam(i,k))*(sc(i) - sbh(i,k))
               vtemp3(i) = vtemp2(i)
               if ((beta(i)*vtemp2(i) - vtemp1(i)) > 0.) then
                  betamx(i) = 0.99*(vtemp1(i)/vtemp2(i))
                  beta(i)   = max(0.0_r8,min(betamx(i),beta(i)))
               end if
            end do
!
! Second check involves supersaturation of "detrainment layer"
! small amount of supersaturation acceptable (by ssfac factor)
!
!cdir nodep
!DIR$ CONCURRENT
            do ii=1,len1
               i = indx1(ii)
               if (hb(i,k-1) < hbs(i,k-1)) then
                  vtemp1(i) = vtemp1(i)*rpdel(i,k)
                  temp2 = gam(i,k-1)*(sbh(i,k) - sc(i) + cldwtr(i)) -  &
                          hbh(i,k) + hc(i) - sc(i) + sbh(i,k)
                  temp3 = vtemp3(i)*rpdel(i,k)
                  vtemp2(i) = (ztodt/cats)*(hc(i) - hbs(i,k))*temp2/ &
                              (pdel(i,k-1)*(hbs(i,k-1) - hb(i,k-1))) + temp3
                  if ((beta(i)*vtemp2(i) - vtemp1(i)) > 0.) then
                     betamx(i) = ssfac*(vtemp1(i)/vtemp2(i))
                     beta(i)   = max(0.0_r8,min(betamx(i),beta(i)))
                  end if
               else
                  beta(i) = 0.0
               end if
            end do
!
! Third check to avoid introducing 2 delta x thermodynamic
! noise in the vertical ... constrain adjusted h (or theta e)
! so that the adjustment doesn't contribute to "kinks" in h
!
!cdir nodep
!DIR$ CONCURRENT
            do ii=1,len1
               i = indx1(ii)
               g = min(0.0_r8,hb(i,k) - hb(i,k-1))
               temp1 = (hb(i,k) - hb(i,k-1) - g)*(cats/ztodt)/(hc(i) - hbs(i,k))
               vtemp1(i) = temp1*vtemp1(i) + (hc(i) - hbh(i,k+1))*rpdel(i,k)
               vtemp2(i) = temp1*vtemp3(i)*rpdel(i,k) + (hc(i) - hbh(i,k) - cldwtr(i))* &
                           (rpdel(i,k) + rpdel(i,k+1))
               if ((beta(i)*vtemp2(i) - vtemp1(i)) > 0.) then
                  if (vtemp2(i) /= 0.0) then
                     betamx(i) = vtemp1(i)/vtemp2(i)
                  else
                     betamx(i) = 0.0
                  end if
                  beta(i) = max(0.0_r8,min(betamx(i),beta(i)))
               end if
            end do
         end if
!
! Calculate mass flux required for stabilization.
!
! Ensure that the convective mass flux, eta, is positive by
! setting negative values of eta to zero..
! Ensure that estimated mass flux cannot move more than the
! minimum of total mass contained in either layer k or layer k+1.
! Also test for other pathological cases that result in non-
! physical states and adjust eta accordingly.
!
!cdir nodep
!DIR$ CONCURRENT
         do ii=1,len1
            i = indx1(ii)
            beta(i) = max(0.0_r8,beta(i))
            temp1 = hc(i) - hbs(i,k)
            temp2 = ((1.0 + gam(i,k))*(sc(i) - sbh(i,k+1) + cldwtr(i)) - &
                      beta(i)*vtemp3(i))*rpdel(i,k) - (hbh(i,k+1) - hc(i))*rpdel(i,k+1)
            eta(i) = temp1/(temp2*grav*cats)
            tmass = min(pdel(i,k),pdel(i,k+1))*rgrav
            if (eta(i) > tmass*rtdt .or. eta(i) <= 0.0) eta(i) = 0.0
!
! Check on negative q in top layer (bound beta)
!
            if (shc(i)-shbh(i,k) < 0.0 .and. beta(i)*eta(i) /= 0.0) then
               denom = eta(i)*grav*ztodt*(shc(i) - shbh(i,k))*rpdel(i,k-1)
               beta(i) = max(0.0_r8,min(-0.999*shb(i,k-1)/denom,beta(i)))
            end if
!
! Check on negative q in middle layer (zero eta)
!
            qtest1 = shb(i,k) + eta(i)*grav*ztodt*((shc(i) - shbh(i,k+1)) - &
                     (1.0 - beta(i))*cldwtr(i)*rhlat - beta(i)*(shc(i) - shbh(i,k)))* &
	             rpdel(i,k)
            if (qtest1 <= 0.0) eta(i) = 0.0
!
! Check on negative q in lower layer (bound eta)
!
            fac1 = -(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
            qtest2 = shb(i,k+1) - eta(i)*grav*ztodt*fac1
            if (qtest2 < 0.0) then
               eta(i) = 0.99*shb(i,k+1)/(grav*ztodt*fac1)
            end if
            etagdt(i) = eta(i)*grav*ztodt
         end do
!
#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
         do i=1,ncol
           if ((lat(i)==jloc).and.(lon(i)==iloc).and.(nstep >= nsloc)) then
            write(6,8080) beta(iloc), eta(iloc)
           end if
         enddo
!#                                                                     #
!#######################################################################
#endif
!
! Calculate cloud water, rain water, and thermodynamic changes
!
!cdir nodep
!DIR$ CONCURRENT
         do 30 ii=1,len1
            i = indx1(ii)
            icwmr(i,k) = cldwtr(i)*rhlat
            cldwtr(i) = etagdt(i)*cldwtr(i)*rhlat*rgrav
! JJH changes to facilitate export of cloud liquid water --------------------------------
            totcond(i) = (1.0 - beta(i))*cldwtr(i)
            rnwtr(i) = min(totcond(i),c0*(dzcld(i)-dzmin)*cldwtr(i)) ! rain + snow
            snwtr(i) = fsnow(i,k)*rnwtr(i)		
            dtwtr(i) = totcond(i)-rnwtr(i)		! detraining water
            ds1(i) = etagdt(i)*(sbh(i,k+1) - sc(i))*rpdel(i,k+1)
            dq1(i) = etagdt(i)*(shbh(i,k+1) - shc(i))*rpdel(i,k+1)
            ds2(i) = (etagdt(i)*(sc(i) - sbh(i,k+1)) +  &
                     hlat*grav*cldwtr(i) - beta(i)*etagdt(i)*(sc(i) - sbh(i,k)))*rpdel(i,k)
! JJH change for export of cloud liquid water; must use total condensate 
! since rainwater no longer represents total condensate
            dq2(i) = (etagdt(i)*(shc(i) - shbh(i,k+1)) - grav*totcond(i) - beta(i)* &
                     etagdt(i)*(shc(i) - shbh(i,k)))*rpdel(i,k)
            ds3(i) = beta(i)*(etagdt(i)*(sc(i) - sbh(i,k)) - hlat*grav*cldwtr(i))* &
                     rpdel(i,k-1)
            dq3(i) = beta(i)*etagdt(i)*(shc(i) - shbh(i,k))*rpdel(i,k-1)
!
! Check budget interms of mixing ratios, as needed for water tracers (dcn)
! shc = qu + ql = qu + beta*ql + cn = qu + beta*ql + rr + qd .... right?
!
            iuwmr(i,k) = shc(i)-icwmr(i,k)		! should be saturated
            rtlmr(i,k) =         beta(i)*icwmr(i,k)	! retained liquid
            cndmr(i,k) = (1.0 - beta(i))*icwmr(i,k)	! total condensation to be lost
            rnwmr(i,k) = min(cndmr(i,k),c0*(dzcld(i)-dzmin)*icwmr(i,k)) ! rain 
            dtlmr(i,k) = cndmr(i,k) - rnwmr(i,k)	! detraining water
!!            write(*,*) 'BUDGET at i,k' , i,k
!!            write(*,3) 'A :',shc(i), iuwmr(i,k)+icwmr(i,k), &
!!                             shc(i)-(iuwmr(i,k)+icwmr(i,k))
!!            write(*,3) 'B :',shc(i), iuwmr(i,k)+rtlmr(i,k)+cndmr(i,k), &
!!                             shc(i)-(iuwmr(i,k)+rtlmr(i,k)+cndmr(i,k))
!!            write(*,3) 'C :',shc(i), iuwmr(i,k)+rtlmr(i,k)+rnwmr(i,k)+dtlmr(i,k), &
!!                             shc(i)-(iuwmr(i,k)+rtlmr(i,k)+rnwmr(i,k)+dtlmr(i,k))

!
! Isolate convective fluxes for later diagnostics
!
            fslkp = eta(i)*(sc(i) - sbh(i,k+1))
            fslkm = beta(i)*(eta(i)*(sc(i) - sbh(i,k)) - hlat*cldwtr(i)*rtdt)
            fqlkp = eta(i)*(shc(i) - shbh(i,k+1))
            fqlkm = beta(i)*eta(i)*(shc(i) - shbh(i,k))
!
! Update thermodynamic profile (update sb, hb, & hbs later)
!
            tb (i,k+1) = tb(i,k+1)  + ds1(i)*rcp
            tb (i,k  ) = tb(i,k  )  + ds2(i)*rcp
            tb (i,k-1) = tb(i,k-1)  + ds3(i)*rcp
            shb(i,k+1) = shb(i,k+1) + dq1(i)
            shb(i,k  ) = shb(i,k  ) + dq2(i)
            shb(i,k-1) = shb(i,k-1) + dq3(i)
!
! ** Update diagnostic information for final budget **
! Tracking precipitation, temperature & specific humidity tendencies,
! rainout term, convective mass flux, convective liquid
! water static energy flux, and convective total water flux
! The variable afac makes the necessary adjustment to the
! diagnostic fluxes to account for adjustment time scale based on
! how relaxation time scale is to be applied (column vs. triplet)
!
            prec(i)    = prec(i) + (rnwtr(i)/rhoh2o)*adjfac
!
! The following variables have units of "units"/second
!
            cmfdt (i,k+1) = cmfdt (i,k+1) + ds1(i)*rtdt*adjfac
            cmfdt (i,k  ) = cmfdt (i,k  ) + ds2(i)*rtdt*adjfac
            cmfdt (i,k-1) = cmfdt (i,k-1) + ds3(i)*rtdt*adjfac
            cmfdq (i,k+1) = cmfdq (i,k+1) + dq1(i)*rtdt*adjfac
            cmfdq (i,k  ) = cmfdq (i,k  ) + dq2(i)*rtdt*adjfac
            cmfdq (i,k-1) = cmfdq (i,k-1) + dq3(i)*rtdt*adjfac
! JJH changes to export cloud liquid water --------------------------------
            qcl    (i,k  ) =  (1.-fice(i,k))*(grav*dtwtr(i)*rpdel(i,k))*rtdt*adjfac
            qci    (i,k  ) =  (   fice(i,k))*(grav*dtwtr(i)*rpdel(i,k))*rtdt*adjfac
            cmfdqr(i,k  ) = cmfdqr(i,k  ) + (grav*rnwtr(i)*rpdel(i,k))*rtdt*adjfac
            cmfmc (i,k+1) = cmfmc (i,k+1) + eta(i)*adjfac
            cmfmc (i,k  ) = cmfmc (i,k  ) + beta(i)*eta(i)*adjfac
!
! The following variables have units of w/m**2
!
            cmfsl (i,k+1) = cmfsl (i,k+1) + fslkp*adjfac
            cmfsl (i,k  ) = cmfsl (i,k  ) + fslkm*adjfac
            cmflq (i,k+1) = cmflq (i,k+1) + hlat*fqlkp*adjfac
            cmflq (i,k  ) = cmflq (i,k  ) + hlat*fqlkm*adjfac
30          continue
!
! Next, convectively modify passive constituents
! For now, when applying relaxation time scale to thermal fields after
! entire column has undergone convective overturning, constituents will
! be mixed using a "relaxed" value of the mass flux determined above
! Although this will be inconsistant with the treatment of the thermal
! fields, it's computationally much cheaper, no more-or-less justifiable,
! and consistent with how the history tape mass fluxes would be used in
! an off-line mode (i.e., using an off-line transport model)
!
            do 50 m=2,pcnst+pnats    ! note: indexing assumes water is first field
               if (cnst_get_type_byind(m).eq.'dry') then
                  pd(:ncol,:) = pdeldry(:ncol,:)
                  rpd(:ncol,:) = rpdeldry(:ncol,:)
                  pm(:ncol,:) = pmiddry(:ncol,:)
               else
                  pd(:ncol,:) = pdel(:ncol,:)
                  rpd(:ncol,:) = rpdel(:ncol,:)
                  pm(:ncol,:) = pmid(:ncol,:)
               endif
!
! Water tracers: Allow condensation during transport
! This differs from the constituent transport in that (as well as
! condensation and fractionation) we compute the tendencies and update
! the profile such that the adjustment is identical to that total water.
! This makes a difference only if adjfac /= 1.
! This is optionally by-passed.
!
              if (trace_water .and. wtrc_is_vap(m) .and. lwtrccmf) then
                pd(:ncol,:)  = pdel(:ncol,:)	!\
                rpd(:ncol,:) = rpdel(:ncol,:)	! }- "wet" as per water
                pm(:ncol,:)  = pmid(:ncol,:)	!/

                do 41 ii = 1, len1
                  i = indx1(ii)
!
! If any of the reported values of the constituent is negative in
! the three adjacent levels, nothing will be done to the profile
! (trb check means dq + checks for last level updates ... also see logic below)
!
!!                  if ((dq(i,k+1,m) < 0.0) .or. (dq(i,k,m) < 0.0) .or. (dq(i,k-1,m) < 0.0)) go to 41
                  if ((trb(i,k+1,m) < 0.0) .or. (trb(i,k,m) < 0.0) .or. (trb(i,k-1,m) < 0.0)) then
                    write(*,*) '(cmfmca) bypassing negative value found in water tracer triplet.'
                    go to 41
                  end if

! Estimate interface values (notice geometric form differs from total water)
                  cmrh(i,k  ) = qhalf_wtrc(trb(i,k-1,m), trb(i,k  ,m))
                  cmrh(i,k+1) = qhalf_wtrc(trb(i,k  ,m), trb(i,k+1,m))

! Set in-cloud specific properties 
!   (should make use of PBL perturbation as is done for total water)
!   (here scale qpert to also account to saturation condition as for m=1)
                  cmrc(i) = trb(i,k+1,m) + fprt(i,k)*qpert(i,m)
! 
! Assign problem as being EITHER wet or frozen 
! (otherwise we will need to define the whole budget for both ice and liquid)
!
                  told = tcld(i,k)
                  tnew = tcld(i,k) - icwmr(i,k)*hlat/cp	! minus liquid water potential
                  tbar = 0.5*(told + tnew)          ! tamb + tcld
!!                  write(*,*) 'CMFMCA: T=',told,tnew-told
!
! Solve isotopic budgets to get the rain production in the updraft
! Budget given flux mass cmrc from level k+1: 
!        vapour is at saturation and over saturated watr is liquid
!        level k   - (1-beta) of mass:
!                  - vapour is saturated, and all cloud liquid is
!                  converted to rain
!        level k-1 -  beta times mass of BOTH liquid and vapour detrains
!
! For fractionation, assume all cloud liquid and ice is produced during
! this one substep (consistent with cloud model, but not ideal for isotopes)
! Similarly, no initial precipitation. Solve the mixed-cloud model
! assuming distillation to ice, and equilibration for liquid. No
! fractionation during coelecence to rain and snow.
! As the model does not understand ice/liquid, rain/snow differences, 
! just merge them for output precipitation, and for output detraining
! condensate.... again, less than ideal (need to add cloud ice budget to cam!)
!
! Assume detraining water is all retained for fractionation, then dumped
! out, after condensation processes are all done.
!
! Budget is this:
!   shc = qu + ql = qu + beta*ql + cn = qu + beta*ql + rr + qd 
!   totcond = rain + detrainment = shc - (iuwmr +  beta*icwmr), 
!   (1-beta)*icwmr is retained above
!
                  dliqmt = 0.
                  dicefz = 0.

                  isp(1) = 1
                  vapold(1) = shc(i)
                  liqold(1) = 0.
                  iceold(1) = 0.
                  vapent(1) = 0.

                  isp(2) = iwspec(m)
                  vapold(2) = cmrc(i)
                  liqold(2) = 0.
                  iceold(2) = 0.
                  vapent(2) = 0.
!
! For "new" state, compute the net liquid, and fraction of net liquid which will detrain
! NOTE : adjfac should appear in "new" otherwise we get too much distillation.
! There are two possabilities for detraining liquid: at isotopic
! composition of cloud liq/ice or at rain/snow. Both a justifiable given
! we only have two size distinktion (small, and big).
!  i.e.,  1) all stuff is dumped out... big ones fall, little ones stay as cloud
!         2) drops stay incloud, and detrain as an after thought
!
! This is tuned via fsnk (some mixture of two end memebers)
!
                  qliq = (1.-fice (i,k))*rtlmr(i,k)
                  qice = (   fice (i,k))*rtlmr(i,k)
                  qrnw = (1.-fsnow(i,k))*rnwmr(i,k)
                  qsnw = (   fsnow(i,k))*rnwmr(i,k)
                  ldet = (1.-fice (i,k))*dtlmr(i,k)
                  idet = (   fice (i,k))*dtlmr(i,k)
!
                  netliq = qliq + (1.-fsnk_cmf)*ldet
                  netice = qice + (1.-fsnk_cmf)*idet
                  ntlsnk = qrnw + (   fsnk_cmf)*ldet
                  ntisnk = qsnw + (   fsnk_cmf)*idet
!
                  fntld = wtrc_ratio((1.-fsnk_cmf)*ldet, netliq)
                  fntid = wtrc_ratio((1.-fsnk_cmf)*idet, netice)
                  flskd = wtrc_ratio((   fsnk_cmf)*ldet, ntlsnk)
                  fiskd = wtrc_ratio((   fsnk_cmf)*idet, ntisnk)
!
! FINALLY, assign total water output
!
                  vapnew(1) = iuwmr(i,k)
                  liqnew(1) = netliq
                  icenew(1) = netice
                  rainpr(1) = ntlsnk
                  snowpr(1) = ntisnk
                  vapdet(1) = 0.
!
! Scale dicm input by the adjustment factor so as to not over do the fractionation
!
                  dliqmt    = dliqmt    * adjfac
                  dicefz    = dicefz    * adjfac

                  vapold(1) = vapold(1) * adjfac
                  liqold(1) = liqold(1) * adjfac
                  iceold(1) = iceold(1) * adjfac
                  vapent(1) = vapent(1) * adjfac
!
                  vapold(2) = vapold(2) * adjfac
                  liqold(2) = liqold(2) * adjfac
                  iceold(2) = iceold(2) * adjfac
                  vapent(2) = vapent(2) * adjfac
!
                  vapnew(1) = vapnew(1) * adjfac 
                  liqnew(1) = liqnew(1) * adjfac 
                  icenew(1) = icenew(1) * adjfac 
                  vapdet(1) = vapdet(1) * adjfac
                  rainpr(1) = rainpr(1) * adjfac 
                  snowpr(1) = snowpr(1) * adjfac 

!!                  call t_startf('wiso_dicm_cmfmca')
                  call wiso_dicm(piso   , &		! piso = 2
                                 isp    , feq_cmf, told   , tnew   , & 
                                 vapold , liqold , iceold , vapent , &
                                 vapnew , liqnew , icenew , vapdet , &
                                 rainpr , snowpr , dliqmt , dicefz )

!!                  call t_stopf('wiso_dicm_cmfmca')

!
! De-scale dicm output by the adjustment factor (giving no additional fractionation)
!
                  vapnew(2) = vapnew(2) / adjfac 
                  liqnew(2) = liqnew(2) / adjfac 
                  icenew(2) = icenew(2) / adjfac 
                  vapdet(2) = vapdet(2) / adjfac
                  rainpr(2) = rainpr(2) / adjfac 
                  snowpr(2) = snowpr(2) / adjfac 
!
! Assign dicm output to cmfmca variables
!
                  truwmr(i,k,m)  = vapnew(2)                    !  vapour, obviously
                  trnetliq       = liqnew(2)
                  trnetice       = icenew(2)
                  trntlsnk       = rainpr(2)
                  trntisnk       = snowpr(2)
!
! Decompose net values back to cmf type variables
!  Reconstruction as:
!    1) detrain some fraction of cloud liquid/ice and 
!         some fraction of liquid/ice sink
!    2) Retain rest of cloud liquid/ice in the cloud
!    3) Remove the rest of the sink as rain/snow
!

                  trdtlmr(i,k,m) = (   fntld)*trnetliq + (   flskd)*trntlsnk 
                  trdtimr(i,k,m) = (   fntid)*trnetice + (   fiskd)*trntisnk 
                  trrnwmr(i,k,m) =                       (1.-flskd)*trntlsnk
                  trsnwmr(i,k,m) =                       (1.-fiskd)*trntisnk
                  trrtlmr(i,k,m) = (1.-fntld)*trnetliq
                  trrtimr(i,k,m) = (1.-fntid)*trnetice
!
! Thus compute total condensation, and back out "cloud liquid"
!
                  trcndmr(i,k,m) = trrnwmr(i,k,m) + trdtlmr(i,k,m) &
                                 + trsnwmr(i,k,m) + trdtimr(i,k,m) !  cond=rain+snow+detliq+detice
                  trcwmr(i,k,m)  = trcndmr(i,k,m) / (1.-beta(i)) 
!
! Convert needed quantities from mixing ratios to mass units
!
                  trrain(i) = etagdt(i)*(trrnwmr(i,k,m))*rgrav
                  trsnow(i) = etagdt(i)*(trsnwmr(i,k,m))*rgrav
                  trdtlq(i) = etagdt(i)*(trdtlmr(i,k,m))*rgrav
                  trdtic(i) = etagdt(i)*(trdtimr(i,k,m))*rgrav
                  trcond(i) = etagdt(i)*(trrnwmr(i,k,m) + trsnwmr(i,k,m) &
                                       + trdtlmr(i,k,m) + trdtimr(i,k,m))*rgrav
!
! Both rain and exported cloud water has same tracer ratio of condensed
! water
!
                  trqcl (i,k,m) = (grav*trdtlq(i)*rpdel(i,k))*rtdt*adjfac
                  trqci (i,k,m) = (grav*trdtic(i)*rpdel(i,k))*rtdt*adjfac
                  trdqr(i,k,m)  = (grav*(trrain(i)+trsnow(i))*rpdel(i,k))*rtdt*adjfac
!
! Compute the total interface fluxes.
! Top interface flux accounts for enrichment of rain.
!
                  botflx =         etagdt(i)*(cmrc(i) -cmrh(i,k+1))
                  topflx = beta(i)*etagdt(i)*(truwmr(i,k,m)+trcwmr(i,k,m)-cmrh(i,k  ))

!
! Compute the tendancies 
! Apply logic for positive definite state - this is needed for HDO about once per step.
! This is teypically due to the "rain" component. scaling the rain plus
! incloud flux is a sensible choice - trying to scale just the rain part
! first would also be a good choice.
! Strictly, if we change the fluxes, the fractionation should change...
! but this inconsistency should be small.
!
!!                  dcmr1(i) = -(botflx         )*rpdel(i,k+1)
!!                  dcmr2(i) =  (botflx - topflx)*rpdel(i,k  ) - grav*trcond(i)*rpdel(i,k)
!!                  dcmr3(i) =  (         topflx)*rpdel(i,k-1)

                  dcmr1(i) = -botflx*rpdel(i,k+1)
                  efac1    = 1.0
                  efac2    = 1.0
                  efac3    = 1.0
!
                  if (trb(i,k+1,m)+dcmr1(i) < 0.0) then
                     efac1 = max(tiny,abs(trb(i,k+1,m)/dcmr1(i)) - eps)
                  end if
!
                  if (efac1 == tiny .or. efac1 > 1.0) efac1 = 0.0
                  dcmr1(i) = -efac1*botflx*rpdel(i,k+1)
                  dcmr2(i) = (efac1*botflx - topflx)*rpdel(i,k) - grav*trcond(i)*rpdel(i,k)
!
                  if (trb(i,k,m)+dcmr2(i) < 0.0) then
                     efac2 = max(tiny,abs(trb(i,k  ,m)/dcmr2(i)) - eps)
                  end if
!
                  if (efac2 == tiny .or. efac2 > 1.0) efac2 = 0.0
                  dcmr2(i) = (efac1*botflx - efac2*topflx)*rpdel(i,k) &
                                           - efac2*grav*trcond(i)*rpdel(i,k)
                  dcmr3(i) = efac2*topflx*rpdel(i,k-1)
!
                  if (trb(i,k-1,m)+dcmr3(i) < 0.0) then
                     efac3 = max(tiny,abs(trb(i,k-1,m)/dcmr3(i)) - eps)
                  end if
!
                  if (efac3 == tiny .or. efac3 > 1.0) efac3 = 0.0
                  efac3    = min(efac2,efac3)
                  dcmr2(i) = (efac1*botflx - efac3*topflx)*rpdel(i,k) &
                                           - efac3*grav*trcond(i)*rpdel(i,k)
                  dcmr3(i) = efac3*topflx*rpdel(i,k-1)
!!!
!!! Report mass conservation adjuatment....
!!                  if (efac1/=1. .or. efac2/=1. .or. efac3/=1.) then
!!                    write(*,*) 'cmfmca efac/=0: m,i,k=',m, i,k
!!                    write(*,*) efac1,efac2,efac3
!!                  end if

!
! Update the tracer profile (not adjusted)
!
                  trb(i,k+1,m) = trb(i,k+1,m) + dcmr1(i)
                  trb(i,k  ,m) = trb(i,k  ,m) + dcmr2(i)
                  trb(i,k-1,m) = trb(i,k-1,m) + dcmr3(i)
!
! Update tendencies to be used by physics_update
!
                  trdq(i,k+1,m) = trdq(i,k+1,m) + dcmr1(i)*rtdt*adjfac
                  trdq(i,k  ,m) = trdq(i,k  ,m) + dcmr2(i)*rtdt*adjfac
                  trdq(i,k-1,m) = trdq(i,k-1,m) + dcmr3(i)*rtdt*adjfac
!
                  trprec(i,m) = trprec(i,m) + ((trrain(i) + trsnow(i))/rhoh2o)*adjfac

!
! debugging diagnostics
!
#define TRCDIAGS
#ifdef TRCDIAGS
!                  if (m == 4 .and.  i==8.and.k==25) then
                  if (m == 4 .and.  &
                      (abs(dq1(i)-dcmr1(i)) > 1.e-12 .or. &
                       abs(dq2(i)-dcmr2(i)) > 1.e-12 .or. &
                       abs(dq3(i)-dcmr3(i)) > 1.e-12 .or. &
                        abs(prec(i)-trprec(i,m)) > 1.e-12) ) then

3                   format(a10,3e16.6)
                    write(*,*)'WTRC_CMFMCA k:',i,k
                    write(*,*)'fntld   :',fntld
                    write(*,*)'adjfac  :',adjfac
                    write(*,3)'CMRH k  :',shbh(i,k  ),cmrh(i,k  ),shbh(i,k  )-cmrh(i,k  )
                    write(*,3)'CMRH k+1:',shbh(i,k+1),cmrh(i,k+1),shbh(i,k+1)-cmrh(i,k+1)
                    write(*,3)'CMRC    :',shc (i)    ,cmrc(i    ),shc (i    )-cmrc(i    )
                    write(*,3)'QPERT   :',qpert(i,1) ,qpert(i,m),qpert(i,1)-qpert(i,m)
                    write(*,3)'ICWMR   :',icwmr(i,k),trcwmr(i,k,m),icwmr(i,k)-trcwmr(i,k,m)
                    write(*,3)'IUWMR   :',iuwmr(i,k),truwmr(i,k,m),iuwmr(i,k)-truwmr(i,k,m)
                    write(*,3)'RNWMR   :',rnwmr(i,k),trrnwmr(i,k,m)+trsnwmr(i,k,m),rnwmr(i,k)-trrnwmr(i,k,m)-trsnwmr(i,k,m)
                    write(*,3)'RTLMR   :',rtlmr(i,k),trrtlmr(i,k,m),rtlmr(i,k)-trrtlmr(i,k,m)
                    write(*,3)'DTLMR   :',dtlmr(i,k),trdtlmr(i,k,m)+trdtimr(i,k,m),dtlmr(i,k)-trdtlmr(i,k,m)-trdtimr(i,k,m)
                    write(*,3)'CNDMR   :',cndmr(i,k),trcndmr(i,k,m),cndmr(i,k)-trcndmr(i,k,m)
                    write(*,3)'qcond   :', shc(i)- iuwmr(i,k)  -beta(i)* icwmr(i,k), &
                                          cmrc(i)-truwmr(i,k,m)-beta(i)*trcwmr(i,k,m), &
                                         ( shc(i)- iuwmr(i,k)  -beta(i)* icwmr(i,k)    ) &
                                        -(cmrc(i)-truwmr(i,k,m)-beta(i)*trcwmr(i,k,m))
                    write(*,3)'DTLTR   :',(1.-fice(i,k))*dtwtr(i),trdtlq(i),(1.-fice(i,k))*dtwtr(i)-trdtlq(i)
                    write(*,3)'DTITR   :',    fice(i,k) *dtwtr(i),trdtic(i),    fice(i,k) *dtwtr(i)-trdtic(i)
                    write(*,3)'RNWTR   :',(1.-fsnow(i,k))*rnwtr (i),trrain(i),(1.-fsnow(i,k))*rnwtr(i)-trrain(i)
                    write(*,3)'SNWTR   :',snwtr(i),trsnow(i),snwtr(i)-trsnow(i)

                    write(*,3)'TOTCOND :',totcond(i),trcond(i),totcond(i)-trcond(i)
                    write(*,3) 'vapold I',vapold(2)  ,vapold(1) ,vapold(2)-vapold(1)
                    write(*,3) 'liqold I',liqold(2)  ,liqold(1) ,liqold(2)-liqold(1)
                    write(*,3) 'iceold I',iceold(2)  ,iceold(1) ,iceold(2)-iceold(1)
                    write(*,3) 'vapent I',vapent(2)  ,vapent(1) ,vapent(2)-vapent(1)
                    write(*,3) 'vapnew I',vapnew(2)  ,vapnew(1) ,vapnew(2)-vapnew(1)
                    write(*,3) 'liqnew I',liqnew(2)  ,liqnew(1) ,liqnew(2)-liqnew(1)
                    write(*,3) 'icenew I',icenew(2)  ,icenew(1) ,icenew(2)-icenew(1)
                    write(*,3) 'vapdet I',vapdet(2)  ,vapdet(1) ,vapdet(2)-vapdet(1)
                    write(*,3) 'rainpr I',rainpr(2)  ,rainpr(1) ,rainpr(2)-rainpr(1)
                    write(*,3) 'snowpr I',snowpr(2)  ,snowpr(1) ,snowpr(2)-snowpr(1)

                    write(*,3)'qcl     :',qcl(i,k),trqcl(i,k,m),qcl(i,k)-trqcl(i,k,m)
                    write(*,3)'qci     :',qci(i,k),trqci(i,k,m),qci(i,k)-trqci(i,k,m)
                    write(*,3)'DQR     :',cmfdqr(i,k),trdqr(i,k,m),cmfdqr(i,k)-trdqr(i,k,m)
                    write(*,3)'DQ1     :',dq1(i),dcmr1(i),dq1(i)-dcmr1(i)
                    write(*,3)'DQ2     :',dq2(i),dcmr2(i),dq2(i)-dcmr2(i)
                    write(*,3)'DQ3     :',dq3(i),dcmr3(i),dq3(i)-dcmr3(i)
	            write(*,3)'PREC    :',prec(i),trprec(i,m),prec(i)-trprec(i,m)
!                    call endrun('cmfmca debug')
                    write(*,*) 'WARNING WARNING WARNING.... cmfmca isotope mass lost'
                  end if
#endif		/* TRCDIAGS */

 41             continue

              else		! m, not water vapour tracer
!cdir nodep
!DIR$ CONCURRENT
               do 40 ii=1,len1
                  i = indx1(ii)
!
! If any of the reported values of the constituent is negative in
! the three adjacent levels, nothing will be done to the profile
!
                  if ((dq(i,k+1,m) < 0.0) .or. (dq(i,k,m) < 0.0) .or. (dq(i,k-1,m) < 0.0)) go to 40
!
! Specify constituent interface values (linear interpolation)
!
                  cmrh(i,k  ) = 0.5*(dq(i,k-1,m) + dq(i,k  ,m))
                  cmrh(i,k+1) = 0.5*(dq(i,k  ,m) + dq(i,k+1,m))
!
! Specify perturbation properties of constituents in PBL
!
                  pblhgt = max(pblht(i),1.0_r8)
                  if ( (zm(i,k+1) <= pblhgt) .and. dzcld(i) == 0.0 ) then
                     fac1 = max(0.0_r8,1.0-zm(i,k+1)/pblhgt)
                     cmrc(i) = dq(i,k+1,m) + qpert(i,m)*fac1
                  else
                     cmrc(i) = dq(i,k+1,m)
                  end if
!
! Determine fluxes, flux divergence => changes due to convection
! Logic must be included to avoid producing negative values. A bit
! messy since there are no a priori assumptions about profiles.
! Tendency is modified (reduced) when pending disaster detected.
!
                  botflx   = etagdt(i)*(cmrc(i) - cmrh(i,k+1))*adjfac
                  topflx   = beta(i)*etagdt(i)*(cmrc(i)-cmrh(i,k))*adjfac
                  dcmr1(i) = -botflx*rpd(i,k+1)
                  efac1    = 1.0
                  efac2    = 1.0
                  efac3    = 1.0
!
                  if (dq(i,k+1,m)+dcmr1(i) < 0.0) then
                     efac1 = max(tiny,abs(dq(i,k+1,m)/dcmr1(i)) - eps)
                  end if
!
                  if (efac1 == tiny .or. efac1 > 1.0) efac1 = 0.0
                  dcmr1(i) = -efac1*botflx*rpd(i,k+1)
                  dcmr2(i) = (efac1*botflx - topflx)*rpd(i,k)
!
                  if (dq(i,k,m)+dcmr2(i) < 0.0) then
                     efac2 = max(tiny,abs(dq(i,k  ,m)/dcmr2(i)) - eps)
                  end if
!
                  if (efac2 == tiny .or. efac2 > 1.0) efac2 = 0.0
                  dcmr2(i) = (efac1*botflx - efac2*topflx)*rpd(i,k)
!                   dcmr3(i) = efac2*topflx*rpdel(i,k-1)	! bug fix (dcn 22/07/04)
                  dcmr3(i) = efac2*topflx*rpd(i,k-1)
!
                  if (dq(i,k-1,m)+dcmr3(i) < 0.0) then
                     efac3 = max(tiny,abs(dq(i,k-1,m)/dcmr3(i)) - eps)
                  end if
!
                  if (efac3 == tiny .or. efac3 > 1.0) efac3 = 0.0
                  efac3    = min(efac2,efac3)
                  dcmr2(i) = (efac1*botflx - efac3*topflx)*rpd(i,k)
                  dcmr3(i) = efac3*topflx*rpd(i,k-1)
!
                  dq(i,k+1,m) = dq(i,k+1,m) + dcmr1(i)
                  dq(i,k  ,m) = dq(i,k  ,m) + dcmr2(i)
                  dq(i,k-1,m) = dq(i,k-1,m) + dcmr3(i)
40                continue
                 endif			! water tracers
50              continue                ! end of m=2,pcnst+pnats loop
!
! Constituent modifications complete
!
                  if (k == limcnv+1) go to 60
!
! Complete update of thermodynamic structure at integer levels
! gather/scatter points that need new values of shbs and gamma
!
                  do ii=1,len1
                     i = indx1(ii)
                     vtemp1(ii     ) = tb(i,k)
                     vtemp1(ii+len1) = tb(i,k-1)
                     vtemp2(ii     ) = pmid(i,k)
                     vtemp2(ii+len1) = pmid(i,k-1)
                  end do
                  call vqsatd (vtemp1  ,vtemp2  ,estemp  ,vtemp3  , vtemp4  , &
                               2*len1   )    ! using estemp as extra long vector
!cdir nodep
!DIR$ CONCURRENT
                  do ii=1,len1
                     i = indx1(ii)
                     shbs(i,k  ) = vtemp3(ii     )
                     shbs(i,k-1) = vtemp3(ii+len1)
                     gam(i,k  ) = vtemp4(ii     )
                     gam(i,k-1) = vtemp4(ii+len1)
                     sb (i,k  ) = sb(i,k  ) + ds2(i)
                     sb (i,k-1) = sb(i,k-1) + ds3(i)
                     hb (i,k  ) = sb(i,k  ) + hlat*shb(i,k  )
                     hb (i,k-1) = sb(i,k-1) + hlat*shb(i,k-1)
                     hbs(i,k  ) = sb(i,k  ) + hlat*shbs(i,k  )
                     hbs(i,k-1) = sb(i,k-1) + hlat*shbs(i,k-1)
                  end do
!
! Update thermodynamic information at half (i.e., interface) levels
!
!DIR$ CONCURRENT
                  do ii=1,len1
                     i = indx1(ii)
                     sbh (i,k) = 0.5*(sb(i,k) + sb(i,k-1))
                     if (trace_water) then
                       shbh(i,k) = qhalf_wtrc(shb(i,k-1),shb(i,k))
                     else
                       shbh(i,k) = qhalf(shb(i,k-1),shb(i,k),shbs(i,k-1),shbs(i,k))
                     end if
                     hbh (i,k) = sbh(i,k) + hlat*shbh(i,k)
                     sbh (i,k-1) = 0.5*(sb(i,k-1) + sb(i,k-2))
                     if (trace_water) then
                       shbh(i,k-1) = qhalf_wtrc(shb(i,k-2),shb(i,k-1))
                     else
                       shbh(i,k-1) = qhalf(shb(i,k-2),shb(i,k-1),shbs(i,k-2),shbs(i,k-1))
                     end if
                     hbh (i,k-1) = sbh(i,k-1) + hlat*shbh(i,k-1)
                  end do
!
#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
!#    this update necessary, only if debugging diagnostics requested   #
!#                                                                     #
                  do i=1,ncol
                     if (lat(i) == jloc .and. nstep >= nsloc) then
                        call vqsatd(tb(i,k+1),pmid(i,k+1),es,shbs(i,k+1),gam(i,k+1), &
                                    1)
                        sb (i,k+1) = sb(i,k+1) + ds1(i)
                        hb (i,k+1) = sb(i,k+1) + hlat*shb(i,k+1)
                        hbs(i,k+1) = sb(i,k+1) + hlat*shbs(i,k+1)
                        kpp = k + 2
                        if (k+1 == pver) kpp = k + 1
                        do kp=k+1,kpp
                           kpm1 = kp-1
                           sbh(i,kp)  = 0.5*(sb(i,kpm1) + sb(i,kp))
                           shbh(i,kp) = qhalf(shb(i,kpm1),shb(i,kp),shbs(i,kpm1),shbs(i,kp))
                           hbh(i,kp)  = sbh(i,kp) + hlat*shbh(i,kp)
                        end do
                     end if
                  end do
!#                                                                     #
!#          diagnostic output                                          #
!#                                                                     #
                  do i=1,ncol
                    if ((lat(i)== jloc).and.(lon(i)==iloc).and.(nstep>=nsloc)) then
                     write(6, 8060) k
                     fac = grav*864.
                     do n=limcnv,pver
                        rh  = shb(i,n)/shbs(i,n)
                        write(6,8020)shbh(i,n),sbh(i,n),hbh(i,n),fac*cmfmc(i,n), &
                                     cmfsl(i,n), cmflq(i,n)
!--------------write(6, 8050)
!--------------write(6, 8030) fac*cmfmc(i,n),cmfsl(i,n), cmflq(i,n)
                        write(6, 8040) tb(i,n),shb(i,n),rh,sb(i,n),hb(i,n), &
                                       hbs(i,n), ztodt*cmfdt(i,n),ztodt*cmfdq(i,n), &
	                               ztodt*cmfdqr(i,n)
                     end do
                     write(6, 8000) prec(i)
                    end if
                  end do
!#                                                                     #
!#                                                                     #
!#######################################################################
#endif
!
! Ensure that dzcld is reset if convective mass flux zero
! specify the current vertical extent of the convective activity
! top of convective layer determined by size of overshoot param.
!
60                do i=1,ncol
                     etagt0 = eta(i).gt.0.0
                     if ( .not. etagt0) dzcld(i) = 0.0
                     if (etagt0 .and. beta(i) > betamn) then
                        ktp = k-1
                     else
                        ktp = k
                     end if
                     if (etagt0) then
                        rk=k
                        rktp=ktp
                        cnt(i) = min(cnt(i),rktp)
                        cnb(i) = max(cnb(i),rk)
                     end if
                  end do
70                continue                  ! end of k loop
!
! ** apply final thermodynamic tendencies **
!
!**BAB don't update input profiles
!!$                  do k=limcnv,pver
!!$                     do i=1,ncol
!!$                        t (i,k) = t (i,k) + cmfdt(i,k)*ztodt
!!$                        q(i,k,1) = q(i,k,1) + cmfdq(i,k)*ztodt
!!$                     end do
!!$                  end do
! Set output q tendencies 
! With water tracers, use computed tendency rather than backing out the change
      dq(:ncol,:,1 ) = cmfdq(:ncol,:)
      if (trace_water .and. lwtrccmf) then
        do m = 2, pcnst+pnats
          if (wtrc_is_vap(m)) then
            dq(:ncol,:,m) = trdq(:ncol,:,m)
          else
            dq(:ncol,:,m) = (dq(:ncol,:,m) - q(:ncol,:,m))/ztodt
          end if
        end do
      else		! normal way
        dq(:ncol,:,2:) = (dq(:ncol,:,2:) - q(:ncol,:,2:))/ztodt
      end if
!
! Kludge to prevent cnb-cnt from being zero (in the event
! someone decides that they want to divide by this quantity)
!
                  do i=1,ncol
                     if (cnb(i) /= 0.0 .and. cnb(i) == cnt(i)) then
                        cnt(i) = cnt(i) - 1.0
                     end if
                  end do
!
                  do i=1,ncol
                     precc(i) = prec(i)*rtdt
                  end do
                  if (trace_water) then
                    do i = 1, ncol
                       trprecc(i,1) =prec(i)*rtdt
                       trprecc(i,2:pcnst+pnats) = trprec(i,2:pcnst+pnats)*rtdt
                    end do
                  end if
!
! Ass a kludge, add a tiny bit more heat to account for the fact that we
! did no account for the latent heat of fusion when condensing to ice
!
!!  do k = 1, pver
!!    do i = 1, ncol
!!      cmfdt(i,k) = cmfdt(i,k) + latice*qci(i,k)/ztodt
!!        write(*,*) ' extra:',latice*qci(i,k)*pdel(i,k)/grav
!!    end do
!!  end do
!
! Compute reserved liquid (not yet in cldliq) for energy integrals.
! Treat rliq, and rice, as flux out bottom, to be added back later.
! (rliq is like precip, and Rice is like snow)
   do k = 1, pver
      do i = 1, ncol
         rliq(i) = rliq(i) + (qcl(i,k) + qci(i,k))*pdel(i,k)/grav
         rice(i) = rice(i) + (qci(i,k)           )*pdel(i,k)/grav
      end do
   end do
   rliq(:ncol) = rliq(:ncol) /1000.
   rice(:ncol) = rice(:ncol) /1000.

#if ( defined DIAGNS )
!#######################################################################
!#                                                                     #
!#    we're done ... show final result if debug diagnostics requested  #
!#                                                                     #
                  do i=1,ncol
                    if ((lat(i)==jloc).and.(lon(i)==iloc).and.(nstep>=nsloc)) then
                     fac = grav*864.
                     write(6, 8010)
                     do k=limcnv,pver
                        rh = shb(i,k)/shbs(i,k)
                        write(6, 8020) shbh(i,k),sbh(i,k),hbh(i,k),fac*cmfmc(i,k), &
                                       cmfsl(i,k), cmflq(i,k)
                        write(6, 8040) tb(i,k),shb(i,k),rh   ,sb(i,k),hb(i,k), &
                                       hbs(i,k), ztodt*cmfdt(i,k),ztodt*cmfdq(i,k), &
                                       ztodt*cmfdqr(i,k)
                     end do
                     write(6, 8000) prec(i)
!#                                                                     #
!#       approximate vertical integral of moist static energy and      #
!#       total preciptable water after adjustment and output changes   #
!#                                                                     #
                     hsum2 = 0.0
                     qsum2 = 0.0
                     do k=limcnv,pver
                        hsum2 = hsum2 + pdel(i,k)*rgrav*hb(i,k)
                        qsum2 = qsum2 + pdel(i,k)*rgrav*shb(i,k)
                     end do
!#                                                                     #
                     write (6,8070) hsum1, hsum2, abs(hsum2-hsum1)/hsum2, &
                                    qsum1, qsum2, abs(qsum2-qsum1)/qsum2
                    end if
                  enddo
!#                                                                     #
!#                                                                     #
!#######################################################################
#endif
                  return                 ! we're all done ... return to calling procedure
#if ( defined DIAGNS )
!
! Formats
!
8000              format(///,10x,'PREC = ',3pf12.6,/)
8010              format('1**        TB      SHB      RH       SB', &
                        '       HB      HBS      CAH      CAM       PRECC ', &
                        '     ETA      FSL       FLQ     **', /)
8020              format(' ----- ',     9x,3p,f7.3,2x,2p,     9x,-3p,f7.3,2x, &
                        f7.3, 37x, 0p,2x,f8.2,  0p,2x,f8.2,2x,f8.2, ' ----- ')
8030              format(' ----- ',  0p,82x,f8.2,  0p,2x,f8.2,2x,f8.2, &
                         ' ----- ')
8040              format(' - - - ',f7.3,2x,3p,f7.3,2x,2p,f7.3,2x,-3p,f7.3,2x, &
                        f7.3, 2x,f8.3,2x,0p,f7.3,3p,2x,f7.3,2x,f7.3,30x, &
                         ' - - - ')
8050              format(' ----- ',110x,' ----- ')
8060              format('1 K =>',  i4,/, &
                           '           TB      SHB      RH       SB', &
                           '       HB      HBS      CAH      CAM       PREC ', &
                           '     ETA      FSL       FLQ', /)
8070              format(' VERTICALLY INTEGRATED MOIST STATIC ENERGY BEFORE, AFTER', &
                        ' AND PERCENTAGE DIFFERENCE => ',1p,2e15.7,2x,2p,f7.3,/, &
                        ' VERTICALLY INTEGRATED MOISTURE            BEFORE, AFTER' &,
                        ' AND PERCENTAGE DIFFERENCE => ',1p,2e15.7,2x,2p,f7.3,/)
8080              format(' BETA, ETA => ', 1p,2e12.3)
8090              format (' k+1, sc, shc, hc => ', 1x, i2, 1p, 3e12.4)
#endif
!
end subroutine cmfmca
end module moistconvection
