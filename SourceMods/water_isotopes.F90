#include <misc.h>
#include <params.h>

module water_isotopes
!-----------------------------------------------------------------------
!
! Provides generic tools and constants for water isotope fracionation.
!
! All interface routine are identified by wico_*, etc.
!
! This code works over species indices, rather than the constituent indices
! used in the water_tracers module. As such, MAKE SURE you call these
! routines with species indicies! The tracer variable names do not need to 
! match the species names, which are privided just for diagnostic output.
!
! * This module MUST be includable by CAM and CLM * (be careful with uses)
!
!
! This routine has a bunch of "qtiny" - which could be standardized.
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 16:28:20 MDT 2003
!
!-----------------------------------------------------------------------
! all fractionation factors = 1
#undef NOFRAC
! all kinetic effects off
#undef NOKIN
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils,   only: endrun


  implicit none

  private
  save

! Public interfaces
  public :: wiso_init                  ! iitialize module parameters

  public :: wiso_alpl                  ! look-up liquid/vapor equil. fractn.
  public :: wiso_alpi                  ! look-up ice/vapor equil. fractn.
  public :: wiso_alps                  ! look-up ice/liquid equil. fractn.

  public :: wiso_kmol                  ! kinetic effect for drag coeffecient
!!  public :: wiso_kmolv10               ! kmol (as above) from 10 meter wind (crappy)
  public :: wiso_akel                  ! kinetic fractionation at liquid evaporation
  public :: wiso_akci                  ! kinetic fractnation at ice condensation
  public :: wiso_arlx                  ! fractionation given finite drop size
  public :: wiso_heff                  ! effective humidity function
  public :: wiso_ssatf                 ! supersaturation function

  public :: wiso_delta                 ! compute the delta value
!!  public :: wiso_ratio                 ! compute isotope ratio (use wtrc_ratio)
  public :: wiso_diffs                 ! compute isotope diffusivity
  public :: wiso_efac		        ! equilibrium implicit factor

  public :: wiso_decay                 ! tendency from radioactive decay (HTO)
  public :: wiso_liqvap_equil	       ! equilibrates liquid and vapour
  public :: wiso_vap_distil	       ! distils vapour for some increment
  public :: wiso_dicm		       ! core fractionation routine <- GOOD STUFF

  public :: wiso_get_rnat              ! retrive internal Rnat
  public :: wiso_get_rstd              ! retrive internal Rstd
  public :: wiso_get_roce              ! retrive internal Roce
  public :: wiso_get_rsic              ! retrive internal Rsic
  public :: wiso_get_rao2              ! retrive internal Rao2
  public :: wiso_get_fisub             ! retrive internal fisub
  public :: wiso_get_mwisp             ! retrive internal molecular weight
  public :: wiso_get_epsmw             ! retrive internal molecular weight ratio

! Namelist variable
  logical, public :: wisotope = .false. ! activate water isotopes [off]

! Species indicies - public so thay can be seen by water_tracers
  integer, parameter, public :: ispundef = 0    ! Undefined
  integer, parameter, public :: isph2o   = 1    ! H2O
  integer, parameter, public :: isphdo   = 2    ! HDO
  integer, parameter, public :: isph218o = 3    ! H218O
  integer, parameter, public :: isph217o = 4    ! H217O
  integer, parameter, public :: isphto   = 5    ! HTO
!
  real(r8), parameter :: hlhto = 3.88e+8 ! halflife of HTO (12.43 years) [s]

! Module parameters
  integer , parameter :: pwtspec = 5    ! number of water species (h2o,hdo,h218o,h217o,hto)
!
! Tunable prameters for fractionation scheme
  real(r8), parameter :: dkfac = 0.58   ! diffusive evap. kinetic power law
  real(r8), parameter :: fsata = 1.000  ! supersaturation peramater s = a + bTdegC (Hoffman)
  real(r8), parameter :: fsatb =-0.003  ! supersaturation parameter s = a + bTdegC (Hoffman)
!  real(r8), parameter :: fsata = 1.02   ! supersaturation parameter s = a + bTdegC (Ciais)
!  real(r8), parameter :: fsatb =-0.0038 ! supersaturation parameter s = a + bTdegC (Ciais)
  real(r8), parameter :: ssatmx = 2.00  ! maximum supersaturation
  real(r8), parameter :: fkhum = 0.25   ! effective humidity factor
  real(r8), parameter :: tzero = 273.16 ! supercooled water in stratiform
  real(r8), parameter :: tkini = 273.16 ! min temp. for kinetic effects as ice appears 
!
  real(r8), parameter :: recrit = 1.0   ! critical raynolds number for kmol

  integer , parameter :: nisoitr = 20	! number of iterations in dicm (10 < nitr < 1000)

! Private isotopic constants
!
  character(len=8), dimension(pwtspec), parameter :: & ! species names
      spnam  = (/ 'H216O'   ,'HD16O'   ,'H218O'   ,'H217O'   ,'HTO' /)
!
! Physical constants for isotopic molecules
! TODO: 为什么HDO和HTO的值都有是2？
  real(r8), dimension(pwtspec), parameter :: &  ! isotopic subs.
       fisub = (/ 1._r8    ,2._r8     ,1._r8     ,1._r8     ,2._r8  /)

  real(r8), dimension(pwtspec), parameter :: &  ! molecular weights
       mwisp = (/ 18.      ,19.       ,20.       ,19.       ,20.    /)

  real(r8), dimension(pwtspec), parameter :: &  ! mol. weight ratio
       epsmw = (/ 1._r8    ,19./18.   ,20./18.   ,19./18.   ,20./18 /)

! TODO: 不太明白
  real(r8), dimension(pwtspec), parameter :: &  ! diffusivity ratio (note D/H, not HDO/H2O)
       difrm = (/ 1._r8    ,0.9836504 ,0.9686999 ,0.9836504 ,0.9686999  /)    ! kinetic theory
!       difrm = (/ 1._r8    ,0.9836504 ,0.9686999 ,0.9686999 ,0.9836504  /)    ! this with expk
!       difrm = (/ 1._r8    ,0.9755    ,0.9723    ,0.9873    ,0.96796    /)    ! Merlivat 1978
!       difrm = (/ 1._r8    ,0.9839    ,0.9691    ,0.9873    ,0.96796    /)    ! Cappa etal 2003

! Isotopic ratios in natural abundance (SMOW)
! TODO:又不懂了，如果是自然丰度，为啥H2O是100%呢？
  real(r8), dimension(pwtspec), parameter :: &  ! SMOW isotope ratios
       rnat  = (/ 1._r8, 155.76e-6, 2005.20e-6, 402.00e-6, 77.88e-06 /)

! Prescribed isotopic ratios (largely arbitrary and tunable)
  real(r8), dimension(pwtspec), parameter :: &  ! model standard isotope ratio
       rstd  = (/ 1._r8    ,1.0_r8    ,1.0_r8    ,1.0_r8    ,1.0_r8     /)  ! best numerics
!!       rstd  = (/ 1._r8    ,0.5_r8   ,0.25_r8     ,0.2_r8    , 0.1_r8 /)  ! test numerics
!!       rstd  = (/ 1._r8    ,155.76e-6,2005.20e-6  ,402.00e-6 , Rhto/)     ! natural abundance

! Isotopic ratio of seaice formed from freezing water
  real(r8), dimension(pwtspec), parameter :: &  ! bare sea ice enrichent
       bsic  = (/ 1._r8    ,1.020     ,1.003     ,1.0015    ,1.0105     /)
!!       bsic = (/ 1._r8    ,1.0_r8   ,1.0_r8     ,1.0_r8      ,1.0_r8/)

! Isotope enrichment at ocean surface (better to be computed or read from file)
  real(r8), dimension(pwtspec), parameter :: &  ! mean ocean surface enrichent 
       boce  = (/ 1._r8    ,1.004     ,1.0005    ,1.00025   ,1.0021     /)
!!       boce  = (/ 1._r8    , 1.0128   ,1.0016    ,1.0008      ,1.00671/)  !! LGM
!!       boce  = (/ 1._r8    , 1.000    ,1.000     ,1.000       ,1.000 /)

! Isotope enrighment of atmospheric oxygen (Bender 1999, triple-isotope)
  real(r8), dimension(pwtspec), parameter :: &  ! mean ocean surface enrichent 
       bao2  = (/ 1._r8    ,1.0       ,0.97704   ,0.988222  ,1.0        /)


! Ocean surface kinetic effects paramterization for 10 m wind
  real(r8), parameter, dimension(pwtspec) :: &  ! surface kinetic exchange (kmol) [OBSOLETE]
       aksmc = (/0.        ,0.00528   ,0.006     ,0.00314   ,0.01056    /), &
       akrfa = (/0.        ,0.2508e-3 ,0.285e-3  ,0.1495e-3 ,0.5016e-3  /), &
       akrfb = (/0.        ,0.7216e-3 ,0.82e-3   ,0.430e-3  ,1.4432e-3  /)

! Coefficients for fractionation
  real(r8), parameter, dimension(pwtspec) :: &  ! liquid/vapour
       alpal = (/0.        , 24.844e+3, 1.137e+3 , 1.137e+3 , 24.844e+3 /), &
       alpbl = (/0.        ,-76.248   ,-0.4156   ,-0.4156   ,-76.248    /), &
       alpcl = (/0.        , 52.612e-3,-2.0667e-3,-2.0667e-3, 52.612e-3 /)

  real(r8), parameter, dimension(pwtspec) :: &  ! ice/vapour
       alpai = (/0.        , 16288.   , 0.       , 0.       , 16288.    /), &
       alpbi = (/0.        , 0.       , 11.839   , 11.839   , 0.        /), &
       alpci = (/0.        ,-9.34e-2  ,-28.224e-3,-28.224e-3,-9.34e-2   /)

  real(r8), parameter, dimension(pwtspec) :: &  ! solid/liquid
       alpas = (/1.0       ,1.02      ,1.003     ,1.003     ,1.02       /)

! Modifier for non-standard species EQUILIBRIUM 
  real(r8), parameter, dimension(pwtspec) :: &  ! specices fractionation modifier
       expa  = (/ 1._r8    ,1.0       ,1.0      ,0.525      ,2.0        /)
!       expa  = (/ 1._r8    ,1.0       ,1.0      ,0.52441    ,2.0        /)

! Modifier for non-standard species KINETIC (from lnD/lnDi)
  real(r8), parameter, dimension(pwtspec) :: &  ! specices fractionation modifier
       expk  = (/ 1._r8    ,1.0       ,1.0      ,1.0        ,1.0      /)
!       expk  = (/ 1._r8    ,1.0       ,1.0      ,0.51838    ,1.9291      /)

contains

!=======================================================================
  subroutine wiso_init
!-----------------------------------------------------------------------
! Purpose: Initialize module internal data arrays
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:26 MDT 2003
!-----------------------------------------------------------------------
    write(6,*) 'WISO_INIT: Initializing water isotopes.'
    return
  end subroutine wiso_init

!=======================================================================
  subroutine wiso_kmol(ncols,npts,indx,isp,rbot,zbot,ustar,alpkn)
!-----------------------------------------------------------------------
!
! Purpose: compute kinetic modifier for drag coefficient (Merlivat & Jouzel)
!
! Method:
!   Code solves Brutsaert equations for theturbulent layer using GCM computed
!   quantities.  Operates on a vector of points.
!
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 14:05:38 MDT 2003
!
!-----------------------------------------------------------------------
    use shr_kind_mod, only: r8 => shr_kind_r8
    use physconst,    only: gravit, karman

    implicit none

    real(r8), parameter :: difair = 2.36e-5 ! molecular diffusivity of air
    real(r8), parameter :: muair  = 1.7e-5  ! dynamic viscosity of air
                                            ! about 17 degC, 1.73 at STP (Salby)
!---------------------------- Arguments --------------------------------
    integer , intent(in)  :: ncols        ! number of columns dimension
    integer , intent(in)  :: npts         ! number of point to compute
    integer , intent(in)  :: indx(ncols)  ! indicies to points
    integer , intent(in)  :: isp          ! species flag
    real(r8), intent(in)  :: rbot(ncols)  ! density of lowest layer (kg/m3)
    real(r8), intent(in)  :: zbot(ncols)  ! height of lowest level (m)
    real(r8), intent(in)  :: ustar(ncols) ! Friction velocity (m/s)
!
    real(r8), intent(out) :: alpkn(ncols) ! kinetic fractionation factor (1-kmol)

!------------------------- Local Variables -----------------------------
    integer ii, i
    real(r8) z0                 ! roughness length (constant in cam 9.5e-5)
    real(r8) reno               ! surface reynolds number
    real(r8) tmr                ! ratio of turbulen to molecular resistance
    real(r8) enn		! diffusive power
    real(r8) sc                 ! Schmidt number (Prandtl number)
    real(r8) vmu                ! kinematic viscocity of air
    real(r8) difn               ! ratio of difusivities to the power of n
    real(r8) difrmj		! isotopic diffusion with substitutions

    real(r8) kmol               ! Merlivals k_mol
!-----------------------------------------------------------------------
!
!!    difrmj = difrm(isp)/fisub(isp)
    difrmj = difrm(isp)
!
! Loop over all points
!
    do ii = 1, npts
      i = indx(ii)
!
!!      z0 = zzocen                     ! CAM prescribed (wrong?)
      z0 = ustar(i)**2 / (81.1*gravit)  ! Charnock's equation
      vmu = muair / rbot(i)             ! kinematic viscosity
      Sc  = vmu/difair
      reno = ustar(i)*z0 / vmu       ! reynolds number
!
      if (reno < recrit) then        ! Smooth (Re < 0.13)
         enn = 2./3.
         tmr  = ( (1./karman)*log(ustar(i)*zbot(i) / (30. * vmu)) ) / (13.6 * Sc**(2./3))
      else                           ! Rough  (Re > 2)
         enn = 1./2.
         tmr  = ( (1./karman)*log(zbot(i)/z0) - 5.) / (7.3 * reno**(1./4.) * Sc**(1./2))
      end if

      difn = (1./difrmj)**enn        ! use D/Di, not Di/D
      kmol = (difn - 1.) / (difn + tmr)

      alpkn(i) = 1. - kmol
!!      write(*,*) 'KMOL:',kmol,alpkn(i)

!
! Modify for non standard species
!
!!      alpkn(i) = alpkn(i)**expk(isp)

#ifdef NOKIN
      alpkn(i) = 1.
#endif
!
    end do
!
    return
  end subroutine wiso_kmol

!=======================================================================
  subroutine wiso_kmolv10(ncols,npts,indx,isp,zbot,ubot,vbot,alpkn)
!-----------------------------------------------------------------------
!
! Purpose: compute kinetic modifier for drag coefficient (Merlivat &
! Jouzel)
!
! Method:
!    Uses everyones favorite empirical relation to 10 meter windspeed
!
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 14:05:38 MDT 2003
!
!-----------------------------------------------------------------------
    use shr_kind_mod, only: r8 => shr_kind_r8
    use physconst,    only: gravit, karman

    implicit none

!---------------------------- Arguments --------------------------------
    integer , intent(in)  :: ncols        ! number of columns dimension
    integer , intent(in)  :: npts         ! number of point to compute
    integer , intent(in)  :: indx(ncols)  ! indicies to points
    integer , intent(in)  :: isp          ! species flag
    real(r8), intent(in)  :: zbot(ncols)  ! height of bottom level (m)
    real(r8), intent(in)  :: ubot(ncols)  ! U wind (m/s)
    real(r8), intent(in)  :: vbot(ncols)  ! V wind (m/s)
!
    real(r8), intent(out) :: alpkn(ncols) ! kinetic fractionation fatcor

!------------------------- Local Variables -----------------------------
    integer ii, i
    real(r8) vmag		! windspeed
    real(r8) v10		! 10 meter winds
    real(r8) kmol               ! Merlivat's K_mol
!-----------------------------------------------------------------------
!
! Trivial case
!
    if (isp == isph2o) then
      do ii = 1, npts
        i = indx(ii)
          alpkn(i) = 1.0
        end do
      return
    end if
!
! Loop over points, compute 18O value, then adjust.
!
    do ii = 1, npts
      i = indx(ii)
!
      vmag = max(1.,sqrt(ubot(i)**2 + vbot(i)**2))
      v10 = vmag*alog(zbot(i)/10.)/karman
!
! Compute the value for o18
!
      if (v10 < 7.0) then               ! smooth regime
         kmol = aksmc(isp)
      else                              ! rough regime
         kmol = akrfa(isp)*v10 + akrfb(isp)
      end if
!
      alpkn(i) = 1. - kmol
!
#ifdef NOKIN
      alpkn(i) = 1.0
#endif
!
    end do
!
    return
  end subroutine wiso_kmolv10


!=======================================================================
  function wiso_alpl(isp,tk)
!-----------------------------------------------------------------------
! Purpose: return liquid/vapour fractionation from loop-up tables
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 10:59:13 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp  ! species indes
    real(r8), intent(in)        :: tk    ! temperature (k)
    real(r8) :: wiso_alpl               ! return fractionation
!-----------------------------------------------------------------------
!
    if (isp == isph2o) then
      wiso_alpl = 1.
      return
    end if
!
    wiso_alpl = exp(alpal(isp)/tk**2 + alpbl(isp)/tk + alpcl(isp))
    wiso_alpl = wiso_alpl**expa(isp)

#ifdef NOFRAC
    wiso_alpl = 1.
#endif
!
    return
  end function wiso_alpl

!=======================================================================
  function wiso_alpi(isp,tk)
!-----------------------------------------------------------------------
! Purpose: return ice/vapour fractionation from loop-up tables
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp  ! species indes
    real(r8), intent(in)        :: tk   ! temperature (k)
    real(r8) :: wiso_alpi               ! return fractionation
!-----------------------------------------------------------------------
    if (isp == isph2o) then
      wiso_alpi = 1.
      return
    end if

    wiso_alpi = exp(alpai(isp)/tk**2 + alpbi(isp)/tk + alpci(isp))
    wiso_alpi = wiso_alpi**expa(isp)

#ifdef NOFRAC
    wiso_alpi = 1.
#endif
!
    return
  end function wiso_alpi

!=======================================================================
  function wiso_alps(isp,tk)
!-----------------------------------------------------------------------
! Purpose: return ice/liquid fractionation from loop-up tables
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp  ! species indes
    real(r8), intent(in)        :: tk   ! temperature (k)
    real(r8) :: wiso_alps               ! return fractionation
!-----------------------------------------------------------------------
    if (isp == isph2o) then
      wiso_alps = 1.
      return
    end if
! 
    wiso_alps = alpas(isp)
    wiso_alps = wiso_alps**expa(isp)
!
#ifdef NOFRAC
    wiso_alps = 1.
#endif
!
    return
  end function wiso_alps

!=======================================================================
  function wiso_akel(isp,tk,hum0,alpeq)
!-----------------------------------------------------------------------
! Purpose: return modified fractination for kinetic effects during 
!          liquid evaporation into unsaturated air.
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp   ! species indes
    real(r8), intent(in)        :: tk    ! Temperature (K)
    real(r8), intent(in)        :: hum0  ! initial humidity ()
    real(r8), intent(in)        :: alpeq ! equilibrium fractionation factor
    real(r8) :: wiso_akel                ! return effective fractionation
    real(r8) :: h0			 ! humidity
    real(r8) :: heff			 ! effective humidity
    real(r8) :: difrmj			 ! diffusivity for isotopically substituted molecule
    real(r8) :: dondi			 ! (D / Di)^fdif, (rather than Di/D)
!-----------------------------------------------------------------------
!!    if (tk > tkinl) then		! also do it for supercooled water
      h0 = min(1.0000,hum0)
!!      difrmj = difrm(isp)/fisub(isp)
      difrmj = difrm(isp)
      heff = wiso_heff(h0)
      dondi = (1/difrmj)**dkfac
      wiso_akel = alpeq*heff / (alpeq*dondi*(heff-1.) + 1.)
!!    else
!!      wiso_akel = alpeq
!!    end if
!
! Modify for non-standard isotope
!
!!    wiso_akel = wiso_akel**expk(isp)

#ifdef NOKIN
    wiso_akel = alpeq
#endif

    return
  end function wiso_akel

!=======================================================================
  function wiso_akci(isp,tk,alpeq)
!-----------------------------------------------------------------------
! Purpose: return modified fractination for kinetic effects during 
!          condensation to ice.
!          Make use of supersaturation function. 
! Author:  David Noone <dcn@caltech.edu> - Tue Jul  1 12:02:24 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)        :: isp   ! species indes
    real(r8), intent(in)        :: tk    ! temperature (k)
    real(r8), intent(in)        :: alpeq ! equilibrium fractionation factor
    real(r8) :: wiso_akci               ! return effective fractionation
    real(r8) :: sat1			! super sturation
    real(r8) :: difrmj			! isotopic diffusion for substituted molecule
    real(r8) :: dondi			! D / Di, (rather than Di/D)
!-----------------------------------------------------------------------
!
    if (tk < tkini) then		! anytime below freezing
      sat1 = max(1.0000, wiso_ssatf(tk))
!!      difrmj = difrm(isp)/fisub(isp)
      difrmj = difrm(isp)
      dondi = 1/difrmj
      wiso_akci = alpeq*sat1 / (alpeq*dondi*(sat1-1.) + 1.)
    else
      wiso_akci = alpeq
    end if
!
! Modify for non-standard isotope
!
!!    wiso_akci = wiso_akci**expk(isp)

#ifdef NOKIN
    wiso_akci = alpeq
#endif
!
    return
  end function wiso_akci

  
!=======================================================================
  function wiso_arlx(isp,dtime,pmid,tk,hum,rdrop,qvap,qliq,alpeq)
!-----------------------------------------------------------------------
! Purpose: return an effective fractionation to account for only partial
! equilibration during some time interval (model time step).
! Following jouzel 1986, small drops (< ~ 30um) will totally equilibtate
! in around 1 sec. bigger drops take longer - calulation assumes a
! "characteric" drop size for the drop population
! Author: David Noone <dcn@colorado.edu> - Sun Jul 25 20:27:41 MDT 2004
  !-----------------------------------------------------------------------
    integer , intent(in) :: isp		! species index
    real(r8), intent(in) :: dtime	! time step (seconds)
    real(r8), intent(in) :: pmid	! pressure (Pa)
    real(r8), intent(in) :: tk		! temperature (kelvin)
    real(r8), intent(in) :: rdrop	! drop radius (microns!)
    real(r8), intent(in) :: hum		! humidity 
    real(r8), intent(in) :: qvap	! vapour mixing ratio
    real(r8), intent(in) :: qliq	! liquid mixing ratio
    real(r8), intent(in) :: alpeq	! equilibrium fractionation factor

    real(r8) wiso_arlx			! output effective fractionation

    real(r8) :: taurlx			! relaxation time scale
    real(r8) :: avent			! ventillation factor
    real(r8) :: difair			! isotope diffusivity in air
    real(r8) :: eps			! descrimination
    real(r8) :: feq			! fraction of equilibration

!-----------------------------------------------------------------------
    real(r8) :: taumin = 1.		! 1 second minimum
!-----------------------------------------------------------------------
    eps = alpeq - 1.			
!
! Given drop size, compute adjustment e-folding time
!
     difair = wiso_diffs(isp,tk,pmid)
     avent = wiso_ratio(isph2o, qliq*hum, qvap)	 ! need molecular weights?
     taurlx = (1.0e-6*rdrop)**2 *alpeq*avent/(3._r8*difair)
     taurlx = max(taurlx,taumin)
!
! Given efolding time, how far to equilibrium will we get in the time step
!
     feq = exp(-dtime/taurlx)
     feq = min(feq, 1.)		! needed, as dtime finite for small tau
     feq = max(feq, 0.)		! just to be sure
!
! Form the effective fractionation, by reducing it toward unity.
!
     wiso_arlx = 1. + feq*eps
!
     return
  end function wiso_arlx

!=======================================================================
  function wiso_heff(h0)
!-----------------------------------------------------------------------
! Purpose: Compute effective humidity (Jouzel type thing)
! Author: David Noone <dcn@caltech.edu> - Fri Oct 24 12:06:55 PDT 2003
!-----------------------------------------------------------------------
    real(r8), intent(in)  :: h0       ! initial humidity
    real(r8) :: wiso_heff             ! return humidity (subsaturation)
!-----------------------------------------------------------------------
    wiso_heff = min(1.0, fkhum*h0 + 1.0-fkhum)
    return
  end function wiso_heff

!=======================================================================
  function wiso_ssatf(tk)
!-----------------------------------------------------------------------
! Purpose: Compute supersaturation based on temperature parameterization.
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:14 MDT 2003
!-----------------------------------------------------------------------
    real(r8), intent(in)  :: tk           ! temperature
    real(r8) :: wiso_ssatf            ! return supersaturation
!-----------------------------------------------------------------------
#ifdef OLDWAY
    wiso_ssatf = max(1.0, fsata + fsatb*(tk-tzero))
#else
    wiso_ssatf = fsata + fsatb*(tk-tzero)
!!    wiso_ssatf = max(wiso_ssatf, fsata)
    wiso_ssatf = max(wiso_ssatf, 1.0)
    wiso_ssatf = min(wiso_ssatf, ssatmx)
#endif
    return
  end function wiso_ssatf

!=======================================================================
  function wiso_get_rnat(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal rnat variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:14 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp          ! species index
    real(r8) :: wiso_get_rnat             ! return isotope ratio
!-----------------------------------------------------------------------
    wiso_get_rnat = rnat(isp)
    return
  end function wiso_get_rnat

!=======================================================================
  function wiso_get_rstd(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal Rstd variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:14 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp          ! species index
    real(r8) :: wiso_get_rstd             ! return isotope ratio
!-----------------------------------------------------------------------
    wiso_get_rstd = rstd(isp)
    return
  end function wiso_get_rstd

!=======================================================================
  function wiso_get_roce(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal Roce variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:04 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp          ! species index
    real(r8) :: wiso_get_roce             ! return isotope ratio
!-----------------------------------------------------------------------
    wiso_get_roce = boce(isp)*rstd(isp)
    return
  end function wiso_get_roce

!=======================================================================
  function wiso_get_rsic(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal Rsic variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:28:52 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp          ! species index
    real(r8) :: wiso_get_rsic             ! return isotope ratio
!-----------------------------------------------------------------------
!!    wiso_get_rsic = bsic(isp)*rstd(isp)
    wiso_get_rsic = bsic(isp)*boce(isp)*rstd(isp)
    return
  end function wiso_get_rsic

!=======================================================================
  function wiso_get_rao2(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal Rao2 variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:29:04 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp          ! species index
    real(r8) :: wiso_get_rao2             ! return isotope ratio
!-----------------------------------------------------------------------
    wiso_get_rao2 = bao2(isp)*rstd(isp)
    return
  end function wiso_get_rao2

!=======================================================================
  function wiso_get_fisub(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal fisub variable, based on species index
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 20:28:52 MDT 2003
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp         ! species index
    real(r8) :: wiso_get_fisub           ! return number of substitutions
!-----------------------------------------------------------------------
    wiso_get_fisub = fisub(isp)
    return
  end function wiso_get_fisub

!=======================================================================
  function wiso_get_mwisp(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal mwwsp variable, based on species index
! Author: David Noone <dcn@colorado.edu> - Tue Feb  3 14:49:08 MST 2004
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp         ! species index
    real(r8) :: wiso_get_mwisp           ! return molecular weight
!-----------------------------------------------------------------------
    wiso_get_mwisp = mwisp(isp)
    return
  end function wiso_get_mwisp

!=======================================================================
  function wiso_get_epsmw(isp)
!-----------------------------------------------------------------------
! Purpose: Retrieve internal epsmw variable, based on species index
! Author: David Noone <dcn@colorado.edu> - Tue Feb  3 14:49:08 MST 2004
!-----------------------------------------------------------------------
    integer , intent(in)  :: isp         ! species index
    real(r8) :: wiso_get_epsmw           ! return molecular weight
!-----------------------------------------------------------------------
    wiso_get_epsmw = epsmw(isp)
    return
  end function wiso_get_epsmw

!=======================================================================
  function wiso_delta(isp,qiso,qtot)
!-----------------------------------------------------------------------
! Purpose: Compute isotopic delta value from masses
! Author David Noone <dcn@caltech.edu> - Tue Jul  1 08:32:45 MDT 2003
!-----------------------------------------------------------------------
    integer, intent(in)  :: isp         ! species index
    real(r8),intent(in)  :: qiso        ! isotopic mass
    real(r8),intent(in)  :: qtot        ! isotopic mass
    real(r8) :: wiso_delta              ! return value
!-----------------------------------------------------------------------
    wiso_delta = 1000.*(wiso_ratio(isp,qiso,qtot)/Rstd(isp) - 1.)
    return
  end function wiso_delta

!=======================================================================
  function wiso_ratio(isp,qiso,qtot)
!-----------------------------------------------------------------------
! Purpose: Compute isotopic ratio from masses, with numerical checks
! Author David Noone <dcn@caltech.edu> - Tue Jul  1 08:32:45 MDT 2003
!-----------------------------------------------------------------------
    integer, intent(in)  :: isp         ! species index
    real(r8),intent(in)  :: qiso        ! isotopic mass
    real(r8),intent(in)  :: qtot        ! isotopic mass
    real(r8) :: wiso_ratio              ! return value
!-----------------------------------------------------------------------
    real(r8) :: qtiny = 1.e-16
!-----------------------------------------------------------------------
    if (qtot > 0) then
      wiso_ratio = qiso/(qtot+qtiny)
    else
      wiso_ratio = qiso/(qtot-qtiny)
    end if
!!    wiso_ratio = espmw(isp)*wiso_ratio/fisum(isp)      ! correct!
  end function wiso_ratio

!=======================================================================
  function wiso_diffs(isp,tk,pk)
!-----------------------------------------------------------------------
! Purpose: Compute isotope diffusivity based on the temperature independent
!          ratios difrm. Diffusivity of water vapor in air is from
!          Pruppacher and Klett (1997)
! Author: Sun Wong  <swong@atmos.umd.edu> - (02/19/2004)
!-----------------------------------------------------------------------
    integer,  intent(in)  :: isp        ! species index
    real(r8), intent(in)  :: tk         ! temperature
    real(r8), intent(in)  :: pk         ! pressure
    real(r8) :: dwatvap                 ! diffusivity of water vapor in air
    real(r8) :: wiso_diffs              ! return value
!-----------------------------------------------------------------------
    dwatvap = 2.5007e-5*(tk/298.15_r8)**1.94/pk
    wiso_diffs = difrm(isp)*dwatvap
    return
  end function wiso_diffs 

!=======================================================================
  subroutine wiso_decay(isp,dtime,q,dqdcy)
!-----------------------------------------------------------------------
! Impliments radioactive decay (for tritium, etc)
!-----------------------------------------------------------------------
  integer , intent(in)  :: isp		! species index (MUST BE isphto, for now)
  real(r8), intent(in)  :: q		! mass of stuff
  real(r8), intent(in)  :: dtime	! time interval of decay
  real(r8), intent(out) :: dqdcy	! change in mass due to decay
!-----------------------------------------------------------------------
    if (isp /= isphto) call endrun('(wiso_decay) isp /= isphto: TRITIUM ONLY')
    dqdcy = q * (0.5**(dtime/hlhto) - 1.) / dtime
    return
  end subroutine wiso_decay

!=======================================================================
subroutine wiso_liqvap_equil(alpha, feq0, vaptot, liqtot, vapiso, liqiso, dliqiso)
!-----------------------------------------------------------------------
! Equilibrate vapour and liquid, assuming the input mass is right
! but distribution is now. liquid vapour partitioning is given by input
! total water. 
!
! Given budgets are unchanged   : v + l = q, and vi + li = qi
! And equilibrium of final state: li/l = alpha vi/v
! Solve for final isotope vapour (thus liquid from budget):
!  vi = efac*qi
!    efac = 1/F; F = alpha*(l/v) + 1
!   
! Nothing fancy here.
!   
! David Noone <dcn@colorado.edu> - Wed Jul 21 19:08:54 MDT 2004
!
!-----------------------------------------------------------------------
  implicit none
!---------------------------- Arguments --------------------------------
  real(r8), intent(in)      :: alpha        ! fractionation fatcor 
  real(r8), intent(in)      :: feq0         ! fraction fractionated
  real(r8), intent(in)      :: vaptot       ! total vapour
  real(r8), intent(in)      :: liqtot       ! total liquid
  real(r8), intent(inout)   :: vapiso       ! isotopic vapour
  real(r8), intent(inout)   :: liqiso       ! isotopic liquid
  real(r8), intent(out)     :: dliqiso      ! change in liquid
!------------------------- Local Variables -----------------------------
  real(r8) dviso		            ! change in vapour
  real(r8) qtot, qiso                       ! total mass of total and isotope
  real(r8) efac                             ! equilibration factor
  real(r8) ratio
!-----------------------------------------------------------------------
!  real(r8) :: qtiny = 1.e-22
  real(r8) :: qtiny = 1.e-36
!-----------------------------------------------------------------------
!
  dliqiso = 0.
  qtot = vaptot + liqtot                ! not used
  qiso = vapiso + liqiso
!
! Check for trivial cases - while these fall out in the algebra, 
! the numericas can be a little picky... so treat them explictly.
!
  if (qtot < qtiny) then                ! makes no sence, do nothing
!!     write(*,*) '(wiso_liqvap_equil) no (trivial) water - doing nothing'
     return
  end if
  if (qiso < qtiny) then
!!     write(*,*) '(wiso_liqvap_equil) no (trivial) isotope - doing nothing'
     return
  end if
!
  if (liqtot < qtiny) then
!!     write(*,*) '(wiso_liqvap_equil) no liquid - dump all isotope to vapor'
     dliqiso = -liqiso
     vapiso = vapiso - dliqiso
     liqiso = 0.
     return
  end if
  if (vaptot < qtiny) then
!!     write(*,*) '(wiso_liqvap_equil) no vapour - dump all isotope to liquid'
     dliqiso = vapiso
     vapiso = 0.
     liqiso = liqiso + dliqiso
     return
  end if
!
! Compute mass exchange from Rc = a Rv and v + l + q, in the normal way
!
!!!!  write(*,*) 'TOT:',qtot, qiso
!!#define USEWISOEFAC
!!#undef USEWISOEFAC
!!#ifdef USEWISOEFAC
!!  efac = wiso_efac(qvaptot,qtot,alpha)
!!#else
!!  efac = 1./ (alpha*wiso_ratio(isph2o,liqtot,vaptot) + 1.)
!!  efac = max(efac, 0.)
!!  efac = min(efac, 1.)
!!#endif
!!!
!!  dliqiso = vapiso - efac*qiso
!!!!  write(*,*) 'dliqiso:',dliqiso
!!!
!!! Correct for potential overflow (can't move what's not there)
!!! [this is paranoia, fixes for efac above should do it]
!!!
!!  if (dliqiso > 0.) then
!!    dliqiso = min(dliqiso ,  vapiso )   ! don't move more than available vapour
!!  else
!!    dliqiso = max(dliqiso , -liqiso )   ! don't move more than available liquid
!!  end if

  dviso = wiso_dqequil(alpha,feq0,vaptot,liqtot,vapiso,liqiso)
  dliqiso = -dviso

!!  write(*,*) 'dliqiso fixed:',dliqiso
!
! Update for outpout
!
  liqiso = liqiso + dliqiso
  vapiso = vapiso - dliqiso
!
  return
end subroutine wiso_liqvap_equil

!=======================================================================
function wiso_dqequil(alpha,feq0,vtotnew,ltotnew,visoold,lisoold)
!-----------------------------------------------------------------------
! Computes changes in vapour isotope due to equilibration (optionally
! partial), in a two-phase system. Nothing fancy here:
! Solve the mass balance under constraints dv+dl = 0; Rl = aRv
! Author: David Noone <dcn@colorado.edu> - Wed Aug 11 13:32:00 MDT 2004
!-----------------------------------------------------------------------
  implicit none
!---------------------------- Arguments --------------------------------
  real(r8), intent(in)  :: alpha      ! fractionation factor
  real(r8), intent(in)  :: feq0       ! fraction equilibrated
  real(r8) ,intent(in)  :: vtotnew    ! new vapour
  real(r8) ,intent(in)  :: ltotnew    ! new liquid
  real(r8) ,intent(in)  :: visoold    ! old isotope vapour
  real(r8) ,intent(in)  :: lisoold    ! old isotope liquid
!
  real(r8) wiso_dqequil               ! return value
!-----------------------------------------------------------------------
  real(r8) qiso                       ! total isotope
  real(r8) vieql                      ! new isotope vapour with equilibration
  real(r8) vinof                      ! new isotope vapour no fractionation
  real(r8) visonew                    ! new isotope vapour
  real(r8) dviso                      ! change in isotope vapour
!-----------------------------------------------------------------------
!
  qiso = visoold + lisoold

! fractionating: Rc = alpha Rv
  vieql = qiso * wiso_efac(alpha, vtotnew, ltotnew)

! non-fractionating: Rc = Rv
  vinof = qiso * wiso_efac(1.0  , vtotnew, ltotnew)

! Merge for partial equilibration
  visonew = feq0*vieql + (1.-feq0)*vinof

! Compute tendency
  dviso = visonew - visoold

! Check for numerical overflow
  if (dviso < 0.) then
    dviso = max(dviso , -visoold )   ! don't move more than available vapour
  else
    dviso = min(dviso ,  lisoold )   ! don't move more than available liquid
  end if
! Assign to output
  wiso_dqequil = dviso
end function wiso_dqequil

!=======================================================================
function wiso_efac(alpha, vapnew, liqnew)
!-----------------------------------------------------------------------
! Computes isotopic equilibration factor - there are different ways to
! compute this which have slightly different numerical properties.
! Author: David Noone <dcn@colorado.edu> - Fri Aug 13 11:26:19 MDT 2004
!-----------------------------------------------------------------------
  real(r8)  , intent(in) :: alpha       ! fractionation factor
  real(r8)  , intent(in) :: vapnew      ! new total vapour
  real(r8)  , intent(in) :: liqnew      ! new total liquid
  real(r8) wiso_efac                    ! return value
!-----------------------------------------------------------------------
  real(r8) efac                         ! equilibration factor
  real(r8) alov                         ! alpha times l on v
  real(r8) qtot                         ! total water
!-----------------------------------------------------------------------

#define DIRECTWAY
#ifdef DIRECTWAY         /* This, most obvious way is less precise */
  alov = alpha*wiso_ratio(isph2o, liqnew, vapnew)
#else			/* this one can have more precise numerics */
  alov = wiso_ratio(isph2o, vapnew, vapnew+liqnew)
  alov = alpha*(1.0_r8/alov - 1.0_r8)   ! this is l/v
#endif 
  efac = 1.0_r8/(alov + 1.0_r8)
! 
  efac = max(efac, 0.)
  efac = min(efac, 1.)
! 
  wiso_efac = efac
  return
end function wiso_efac


!=======================================================================
subroutine wiso_vap_distil(alpha,vtotold,vtotnew,visoold,visonew,dvapiso)
!-----------------------------------------------------------------------
! Performs a rayleigh distillation on some vapour increment assuming
! the parcel is in isolation (this is in integral form). 
!
! Given, dvi/dv = alpha (vi/v), integration from vold to vnew gives:
!     vinew = viold*(vnew/vold)**alpha
!
! Nothing fancy here.
!
! David Noone <dcn@colorado.edu> - Thu Jul 22 11:18:29 MDT 2004
!
!-----------------------------------------------------------------------
  implicit none
!---------------------------- Arguments --------------------------------
  real(r8), intent(in)      :: alpha        ! fractionation factor
  real(r8), intent(in)      :: vtotold      ! initial total vapour
  real(r8), intent(in)      :: vtotnew      ! final total vapour
  real(r8), intent(in)      :: visoold      ! initial isotope vapour
  real(r8), intent(out)     :: visonew      ! initial isotope vapour
  real(r8), intent(out)     :: dvapiso      ! change in isotope vapour (diagnostic)
!-----------------------------------------------------------------------
  real(r8) :: qtiny = 1.e-22
!-----------------------------------------------------------------------
  dvapiso = 0.
!
! Check for trivial cases
!
  if (vtotold < vtotnew) then
    write(*,*) '(wtrc_vap_distill) vapour increase, while expecting dectrease'
!   call endrun('ABORT')
  endif
  if (vtotold < qtiny) then
    write(*,*) '(wtrc_vap_distil) no (trivial) vapour - doing nothing.'
    visonew = visoold
    return
  end if
!
! Solve the logarithm, and back out change
!
  visonew = visoold*wiso_ratio(isph2o,vtotnew,vtotold)**alpha
  dvapiso = visonew - visoold
  return
end subroutine wiso_vap_distil

!=======================================================================
subroutine wiso_dicm( niso   , &
                      isp    , feq0   , told   , tnew   , &
                      vapold , liqold , iceold , vapent , &
                      vapnew , liqnew , icenew , vapdet , &
                      rainpr , snowpr , dliqmt , dicefz )
!-----------------------------------------------------------------------
!
! Differential (David's) Isotope Cloud Model (Mirophysics).
!
! This follows the type of model proposed by Merlivat and jouzel, with
! extension of mixed phase by Gendzelman and Arnold, and Ciais and
! Jousel. Used in a plume model, the results should look very much like
! the Frederer et al model. 
!
! Assume changes to total happen linearly.
! This routine is quite generic given some thought about the inputs.
!
!       dvap(m) + dliq(m) + dice(m) + prain(m) + psnow(m) = pvent - pvdet
!
!       dvap(m) = -pliq_vap(m) - pice_vap(m)                             - pvdet(m) + pvent(m)
!       dliq(m) =  pliq_vap(m)               - pice_liq(m) + pliq_ice(m) - prain(m) 
!       dice(m) =                pice_vap(m) + pice_liq(m) - pliq_ice(m) - psnow(m)
!
! For isotopes, there is also and exchange between rain and vapour)
!
! All exchanges are without fractionation EXCEPT:
!  o deplosition of vapour onto ice as disillation (with kinetic
!  fractionatoin)
!  o vapour and liquid are held in equilibrim (equilibrium
!    fractionation,  unless evaporation then kinetic effects)
!  o vapour and rain tend toward equilibrium at some slow rate
!
! Assumes no rain or snow initially. If this is not the case, 
! you might judiciously dump all initial rain into cloud water, 
! and specify an effective precipitation rate (?).
! If this is REALLY a problem, email me, it can be fixed.
! [note to self: snow sublimates no frac, rain as distillation is
! easiest] Similarly, detrainment and entrainment is easily include...
! bot it all costs flops!
!
! In CAM, this module is suitable for shallow and deep convection, 
! and stratiform cloud condensation processes (although the latter
! needs some preprocessing to treat snow melt to rain, and
! evaporation of rain and snow).
!
! After 10 iterations the vapour and liquid delta values are good to the
! 5th significant figure of delta values (8 sig figs for mass). After
! 100 iterations the snow (which is the least accurate) is good to
! almost 3 sig figs delta. I recommend at least 10, and more if you can
! afford it. 20 or 30 seems to be a good choice.
! However, there is a better analytic shortcut that may have better
! convergence properties.
!
! REMBER: FEQ WILL BE APPLIED DIFFERENTIALLY, so ends up something stringe
! (infact it has an exponential form which can be composed to give the
! correct result...)
!
! Note, magnitude of  accumulated numerical imprecision increases with number of
! iterations. Mass imbalance of 1.e-13 after 1000000 iterations - not
! bad.
!
! David Noone <dcn@colorado.edu> - Wed Jul 21 18:07:18 MDT 2004
!
!-----------------------------------------------------------------------
  implicit none
!------------------------- Input Arguments -----------------------------
  integer , intent(in)    :: niso           ! number of isotope species
  integer , intent(in)    :: isp(niso)	    ! species index

  real(r8), intent(in)    :: feq0           ! fraction fractionated
  real(r8), intent(in)    :: told           ! initial temperature
  real(r8), intent(in)    :: tnew           ! final temperature

  real(r8), intent(inout) :: vapold(niso)   ! initial vapour mass
  real(r8), intent(inout) :: liqold(niso)   ! initial cloud liquid mass
  real(r8), intent(inout) :: iceold(niso)   ! initial cloud ice mass
  real(r8), intent(inout) :: vapent(niso)   ! vapour entrainment

  real(r8), intent(inout) :: vapnew(niso)   ! final vapour mass
  real(r8), intent(inout) :: liqnew(niso)   ! final cloud liquid mass
  real(r8), intent(inout) :: icenew(niso)   ! final cloud ice mass
  real(r8), intent(inout) :: vapdet(niso)   ! vapour detrainment
  real(r8), intent(inout) :: rainpr(niso)   ! rain production 
  real(r8), intent(inout) :: snowpr(niso)   ! snow production
!
  real(r8), intent(in)    :: dliqmt         ! change in cloud liquid due to ice melt
  real(r8), intent(in)    :: dicefz         ! change in cloud ice due to liquid freeze

!------------------------- Local Variables -----------------------------

  integer m                         ! isotope number count
  integer itr                       ! iteration count
  integer nitr                      ! iterations for finite integration

  logical ldistdif		    ! differential/integral distillation

  real(r8) pliq_vap(niso)           ! prod. cloud liquid from vapour condensation
  real(r8) pice_vap(niso)           ! prod. cloud ice from vapour condensation
  real(r8) pice_liq(niso)           ! prod. cloud ice from cloud liquid freeze
  real(r8) pliq_ice(niso)           ! prod. cloud liquid from cloud ice melt
  real(r8) prnw_vap(niso)           ! prod. rain water from vapour condensation exchange

  real(r8) pvent(niso)              ! prod. vapour mass by entrainment
  real(r8) pvdet(niso)              ! prod. detrained vapour mass
  real(r8) psnow(niso)              ! prod. snow from conversion of cloud ice
  real(r8) prain(niso)              ! prod. rain from conversion of cloud liq

  real(r8) dvap(niso)               ! change in vapour
  real(r8) dliq(niso)               ! change in cloud liquid
  real(r8) dice(niso)               ! change in cloud ice

  real(r8) qvap(niso)               ! vapour intergrand
  real(r8) qliq(niso)               ! cloud liquid integrand
  real(r8) qice(niso)               ! cloud ice integrand
  real(r8) vdet(niso)               ! vapour detrainment accumulator
  real(r8) rain(niso)               ! rain accumulator
  real(r8) snow(niso)               ! snow accumulator

  real(r8) alpice(niso)             ! fractionation factor for ice (kinetic)
  real(r8) alpliq(niso)             ! fractionation factor for liquid

  real(r8) totalold, totalnew       ! totals for conservation checks
  real(r8) vtotnew, visonew         ! updates from distillation to ice
  real(r8) dvapiso                  ! change in vapour from distillation
  real(r8) tk			    ! temperature (kelvin)
  real(r8) fint			    ! fraction through integration

!-----------------------------------------------------------------------
  real(r8)  :: feq_rn = 0.0	  ! fraction rain vapour equulibration
!  real(r8)  :: feq_rn = 0.1	   ! fraction rain vapour equulibration [VERY SENSITIVE]
!  real(r8)  :: fieq = 0.          ! fraction of initial equilibration (slow convergence)
!  real(r8)  :: fieq = 1.          ! fraction of initial equilibration (better)
  real(r8)  :: fieq = 0.5         ! fraction of initial equilibration (like extra iteration)
!-----------------------------------------------------------------------
  real(r8)  :: qtolerr = 1.e-12	  ! small q for error checks
  real(r8)  :: qtiny   = 1.e-22	  ! essentially zero q
!-----------------------------------------------------------------------
   alpliq(1) = 1.
   alpice(1) = 1.
!
! Check input budget
!
   totalold = vapold(1) + liqold(1) + iceold(1) + vapent(1)
   totalnew = vapnew(1) + liqnew(1) + icenew(1) + vapdet(1) + rainpr(1) + snowpr(1)
   if (abs(totalold - totalnew) > qtolerr) then
      write(*,*) '(wtrc_dicm) WARNING: total budget does not balance: old /= new'
      write(*,2) 'old',vapold(1), liqold(1), iceold(1), vapent(1)
      write(*,2) 'new',vapnew(1), liqnew(1), icenew(1), vapdet(1), rainpr(1), snowpr(1)
      write(*,2) 'total',totalold, totalnew, totalold - totalnew
!!      call endrun('(wiso_dicm) ABORTED.')
   end if
!
! Check things we assume (all quantities positive definite, and temperature reasonable)
!  
!!   if (vapnew(1) < qtolerr) then
!!     write(*,*) '(wiso_dicm) VAPNEW < qtolerr - will cause problems for delta values.',vapnew(1)
!!   end if
!  
   do m = 1, niso
     if (vapold(m) < 0. .or. liqold(m) < 0. .or. iceold(m) < 0.) then
       write(*,*) '(wtrc_dicm) old values < 0.: m=',m
       write(*,2) 'old:',vapold(m),liqold(m),iceold(m)
       write(*,2) 'new:',vapnew(m),liqnew(m),icenew(m)
       vapold(1) = max(vapold(1),0.)
       liqold(1) = max(liqold(1),0.)
       iceold(1) = max(iceold(1),0.)
       call endrun('(wiso_dicm) Abort')
     end if
     if (vapent(m) < 0.) then
        write(*,*) '(wtrc_dicm) entrainment < 0. : m=',m
        write(*,*) 'vapent:',vapent(m)
        vapent(m) = max(vapent(m),0.)
     end if
  end do
  if (vapnew(1) < 0. .or. liqnew(1) < 0. .or. icenew(1) < 0.) then
       write(*,*) '(wtrc_dicm) new values < 0.'
       write(*,2) 'old:',vapold(1),liqold(1),iceold(1)
       write(*,2) 'new:',vapnew(1),liqnew(1),icenew(1)
       vapnew(1) = max(vapnew(1),0.)
       liqnew(1) = max(liqnew(1),0.)
       icenew(1) = max(icenew(1),0.)
     call endrun('(wiso_dicm) Abort')
  end if
  if (vapdet(1) < 0.) then
     write(*,*) '(wtrc_dicm) detrainment < 0. : m=',m
     write(*,*) 'vapdet:',vapdet(1)
     vapdet(1) = max(vapdet(1),0.)
  end if
  if (rainpr(1) < 0. .or. snowpr(1) < 0.) then
     write(*,*) '(wtrc_dicm) precipitation production < 0.'
     write(*,*) 'rain, snow:',rainpr(1),snowpr(1)
     call endrun('(wiso_dicm) Abort...')
     rainpr(1) = max(rainpr(1), 0.)
     snowpr(1) = max(snowpr(1), 0.)
  end if

  if (told > 330 .or. tnew > 330.  .or. &
      told < 130 .or. tnew < 130. ) then
     write(*,*) '(wiso_dicm) input tenperature unreasonable:',told,tnew
!!     call endrun('(wiso_dicm) Abort')
  end if
 2 format(a6,10e16.6)
!
! Given input net changes to total waters, compose small finte differeces.
! These are held constant over all iterations.
!  
   dvap(1)  = (vapnew(1) - vapold(1))
   dliq(1)  = (liqnew(1) - liqold(1))
   dice(1)  = (icenew(1) - iceold(1))
!
! Work out what's going to happen, and decide if we need to iterate.
! If iterating use differential form for distillation, as it is slightly faster.
!
  nitr = nisoitr		! default full iterations
  ldistdif = .true.		! default full differential 
!
  if (abs(dvap(1))   < qtiny .and. &
      abs(dliq(1))   < qtiny .and. &
      abs(dice(1))   < qtiny .and. &
      abs(rainpr(1)) < qtiny .and. &
      abs(snowpr(1)) < qtiny .and. &
      abs(dliqmt)    < qtiny .and. & 
      abs(dicefz)    < qtiny) then 		! nothing to do, just one pass
    nitr = 1
    ldistdif = .false.

  else if (iceold(1) < qtiny .and. icenew(1) < qtiny .and.  &
           snowpr(1) < qtiny .and. rainpr(1) < qtiny .and.  &
           abs(dliqmt)<qtiny .and. abs(dicefz)<qtiny ) then	! equilibrate vap/liq
    nitr = 1

  else if (liqold(1) < qtiny .and. liqnew(1) < qtiny .and.  &
           snowpr(1) < qtiny .and. rainpr(1) < qtiny .and.  &
           abs(dliqmt)<qtiny .and. abs(dicefz)<qtiny ) then	! distil to ice
    nitr = 1
    ldistdif = .false.
  end if
!
! Convert changes to iterative increments
!
   pvent(:niso) = vapent(:niso) / real(nitr)
!
   dvap(1)  = dvap(1)   / real(nitr)
   dliq(1)  = dliq(1)   / real(nitr)
   dice(1)  = dice(1)   / real(nitr)

   pvdet(1) = vapdet(1) / real(nitr)
   prain(1) = rainpr(1) / real(nitr)
   psnow(1) = snowpr(1) / real(nitr)

   pliq_ice(1) = dliqmt / real(nitr)
   pice_liq(1) = dicefz / real(nitr)
!
! Given the amount of cloud liquid/ice melt freeze, 
! we can back out the other terms
!
   pliq_vap(1) = dliq(1) + prain(1) + pice_liq(1) - pliq_ice(1)
   pice_vap(1) = dice(1) + psnow(1) - pice_liq(1) + pliq_ice(1)
   prnw_vap(1) = 0.
!
! Copy inputs to local variables to iterate on
!  
   qvap(:niso) = vapold(:niso)
   qliq(:niso) = liqold(:niso)
   qice(:niso) = iceold(:niso)
   vdet(:niso) = 0.
   rain(:niso) = 0.
   snow(:niso) = 0.
!
! Integrate over small changes, for all isotopes
! Here order is important... if number of iterations are small.
! The iteration is done in two substeps to improve accuracy 
! (although I'm not convinced it's much better than doubling the number
!  of iterations!)
!

   do itr = 1, nitr
!
! Compute temperature, and assign fractionation factors
!
     fint = (real(itr)-0.5) / real(nitr)
     tk = fint*tnew + (1-fint)*told
!
     do m = 2, niso
       alpliq(m) = wiso_alpl(isp(m),tk)		! Kinetic?
       alpice(m) = wiso_alpi(isp(m),tk)
       alpice(m) = wiso_akci(isp(m),tk,alpice(m))
     end do
!
! Start by entraining vapour from the environment
!
     do m = 1, niso
       qvap(m) = qvap(m) + pvent(m)
     end do
!
! Equilibrating liquid and vapour, for a fraction of iteration vapour/liquid production
!
     if (fieq > 0.) then
       qliq(1) = qliq(1) + fieq*pliq_vap(1)
       qvap(1) = qvap(1) - fieq*pliq_vap(1)
       do m = 2, niso
         call wiso_liqvap_equil(alpliq(m),feq0,qvap(1),qliq(1), &
                      qvap(m),qliq(m),pliq_vap(m))
       end do
     end if
!
! Compute the ice melt and liquid freeze, and move it
!
     do m = 2, niso
        pice_liq(m) = pice_liq(1)*wiso_ratio(isph2o,qliq(m),qliq(1))
        pice_liq(m) = min(pice_liq(m), qliq(m))

        pliq_ice(m) = pliq_ice(1)*wiso_ratio(isph2o,qice(m),qice(1))
        pliq_ice(m) = min(pliq_ice(m), qice(m))
     end do
!
     do m = 1, niso
       qice(m) = qice(m) + pice_liq(m) - pliq_ice(m)
       qliq(m) = qliq(m) - pice_liq(m) + pliq_ice(m)
     end do
!
! Compute the ice/vapour deposition/sublimation
! Both integral and differential schemes have about the same convergence
! characteristics.
!
     if (ldistdif) then		!  use differential form for distillation 
       do m = 2, niso
         if (pice_vap(1) > 0.) then ! ice deposition by distillation, kinetic
            pice_vap(m) = alpice(m)*pice_vap(1)*wiso_ratio(isph2o,qvap(m), qvap(1))
            pice_vap(m) = min(pice_vap(m),  qvap(m))
         else                       ! sublimation, no fractionation
            pice_vap(m) = pice_vap(1)*wiso_ratio(isph2o,qice(m), qice(1))
            pice_vap(m) = max(pice_vap(m), -qice(m))
         end if
       end do
     else			! use integral form for distillation
       if (pice_vap(1) > 0.) then   ! ice deposition by distillation, kinetic
         vtotnew = qvap(1) - pice_vap(1)
         vtotnew = max(vtotnew, 0.)
         do m = 2, niso
           call wiso_vap_distil(alpice(m),qvap(1),vtotnew,qvap(m),visonew,dvapiso)
           pice_vap(m) = -dvapiso
           pice_vap(m) = min(pice_vap(m),  qvap(m))
         end do
       else                         ! sublimation, no fractionation
         do m = 2, niso
           pice_vap(m) = pice_vap(1)*wiso_ratio(isph2o,qice(m), qice(1))
           pice_vap(m) = max(pice_vap(m), -qice(m))
         end do
       end if
     end if

     do m = 1, niso
       qice(m) = qice(m) + pice_vap(m)
       qvap(m) = qvap(m) - pice_vap(m)
     end do
!
! equilibrate liquid/vapour equilibrium with second increment
!
     qliq(1) = qliq(1) + (1.-fieq)*pliq_vap(1)
     qvap(1) = qvap(1) - (1.-fieq)*pliq_vap(1)
     do m = 2, niso
       call wiso_liqvap_equil(alpliq(m),feq0,qvap(1),qliq(1), &
                    qvap(m),qliq(m),pliq_vap(m))
     end do
!
! Detrain vapour
!
     do m = 2, niso
       pvdet(m) = pvdet(1)*wiso_ratio(isph2o,qvap(m),qvap(1))
       pvdet(m) = min(pvdet(m), qvap(m))
     end do
!
     do m = 1, niso
       qvap(m) = qvap(m) - pvdet(m)
       vdet(m) = vdet(m) + pvdet(m)
     end do
!
! Complete the iteration by removing rain and snow production
!
     do m = 2, niso             ! must do m=1 first
       prain(m) = prain(1)*wiso_ratio(isph2o,qliq(m),qliq(1))
       prain(m) = min(prain(m), qliq(m))
       psnow(m) = psnow(1)*wiso_ratio(isph2o,qice(m),qice(1))
       psnow(m) = min(psnow(m), qice(m))
     end do
!
     do m = 1, niso
       qliq(m) = qliq(m) - prain(m)
       rain(m) = rain(m) + prain(m)
       qice(m) = qice(m) - psnow(m)
       snow(m) = snow(m) + psnow(m)
     end do
!
! Allows some fraction of the rain to equilibrate with cloud vapour
! (isotope only exchange effect)
!
     if (feq_rn > 0.0) then
       do m = 2, niso
         call wiso_liqvap_equil(alpliq(m),feq_rn,qvap(1),rain(1), &
                      qvap(m),rain(m),prnw_vap(m))
       end do
     end if
!
  end do                ! itr, iteration loop
!
! For final state, ensure liquid and vapour are in equilibrium
! (might bee needed in calling code)
!
  do m = 2, niso
    call wiso_liqvap_equil(alpliq(m),feq0,qvap(1),qliq(1), &
                 qvap(m),qliq(m),pliq_vap(m))
  end do
!
! Assin final values to output, and check budget was done correctly
!
  do m = 2, niso        ! dont assign to total
!
    vapnew(m) = qvap(m)
    liqnew(m) = qliq(m)
    icenew(m) = qice(m)
    vapdet(m) = vdet(m)
    rainpr(m) = rain(m)
    snowpr(m) = snow(m)
!
! Final check for budgets 
! If not balanced, try more iterations
!  Also, could do an adjustment to enforce mass conservation while
!  preserving ratios...
!
    totalold = vapold(m) + liqold(m) + iceold(m) + vapent(m)
    totalnew = vapnew(m) + liqnew(m) + icenew(m) + vapdet(m) + rainpr(m) + snowpr(m)
    if (abs(totalold - totalnew) > qtolerr) then
       write(*,*) '(wtrc_dicm) WARNING - isotope budget not balanced.'
       write(*,*) totalold, totalnew, totalold - totalnew
!!       call endrun('(wiso_dicm) ABORTED.')
    end if
!
  end do

  return
end subroutine wiso_dicm


!=======================================================================
end module water_isotopes

