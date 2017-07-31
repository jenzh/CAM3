#include <misc.h>
#include <params.h>

module water_tracers
!-----------------------------------------------------------------------
!
! Provide core functionality for water tracers.
!
! This module works in tandem with "waterisotope" which specifically
! teats the istopic fractionation.
!
! All interface routine are identified by wtrc_*, etc.
!
! Indexing ASSUMES (as everywhere else in CAM), normal water vapour is
! m=1. Cloud liquid and ice are probably m=2 and m=3... but not
! assumed, as they are correctly registered.
!
! DEFAULT CONFIGURATION is for 3 additional tracers with each of 3 phases, to
! parallel the base CAM water prognosis. This default is invoked by setting the
! namelist variable "wisotope". If not set, don't know what to do, so
! complain and crash.
!
! Note total vapout (Q) is registered as wet, even though it is treated
! as dry. This means the PD coupling is diferent from the vertical
! diffusion. Do get around this we make all water species wet, but have
! spacial cased for vapour, as in main cam code.
!
!
! Code here based on chemistry.F90 and cldcond.F90.
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 15:22:48 MDT 2003
! Modified for tagging C Bitz Apr 24, 2009, only for h2o hdo o18
! with help from David
!
! User must customize number of tagged regions. Search ntags and wtrc_names
! adjust them accordingly. The input mask file is one 2d field that should
! look roughly like this for 4 tagged regions (but probably much larger):
!      00000000000000000000
!      00000000011111000000
!      00000000222211000000
!      00000002222333003300
!      00000000222233300000
!      00000044444444000000
!      00000004444444000000
!      00000000000000000000
!  the tagging sourcemods below will make 4 2d-masks and use them
!  to mask the evaporative fluxes over ice, land, and ocean
!-----------------------------------------------------------------------

! added -DTAGGING to configure at build time instead of here
!/* for tagging - ncnst hard coded */
!#define TAGGING
! define this to turn on h217o and hto, set pcnst=18 */
#undef O17HTO
! DEBUG: define to bypass conservation checks */
#undef NOCHECK
! DEBUG: define to give "sucess" message for checks */
#undef QCHKMESS
! DEBUG: define to terminate when "qchk" failed */
#undef QCHKTERM

!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plat, plev, plevp, plon, masterproc
  use ppgrid,       only: pcols, pver, begchunk, endchunk
  use constituents, only: pcnst,pnats
  use abortutils,   only: endrun
!
  use water_isotopes, only: wisotope

  implicit none
  private
  save

!------------------------ Module Interfaces -----------------------------
!
! Public interfaces
!
  public :: wtrc_init                            ! initialize module parameters
  public :: wtrc_register                        ! register constituents
  public :: wtrc_implements_cnst                 ! checks tracer implementation
  public :: wtrc_is_wtrc			 ! logical function for m = water tracer
  public :: wtrc_is_vap				 ! logical function for m = vapour
  public :: wtrc_is_liq				 ! logical function for m = cloud liquid
  public :: wtrc_is_ice				 ! logical function for m = cloud ice
  public :: wtrc_init_cnst                       ! sets data if not on IC file
  public :: wtrc_init_qpert			 ! initialize boundary layer perturbation
  public :: wtrc_setup_diag			 ! write tracer configuration
  public :: wtrc_qchk1                           ! compare 1d tracer with prognostic
  public :: wtrc_qchk2                           ! compare 2d tracer with prognostic
  public :: wtrc_qchk3                           ! compare all 2d with prognostsic
  public :: wtrc_check			 	 ! checks tracer with prognostic
  public :: wtrc_chkdelta			 ! checks delta values
  public :: wtrc_ratio			  	 ! calulates ratio to precision
  public :: wtrc_ratio_all			 ! calulates ratio for all tracers
!  public :: wtrc_rescale			 ! scaler routine

!------------------- Module Variable Declarations -----------------------

! Namelist variable
  logical, public :: trace_water = .false.       ! set true to activate [off]

  logical, public :: lwtrcland   = .true.        ! set to use simple land model

! Tracer physics control flags
!
  logical, parameter, public :: lh2oadj    = .true.     ! adjust tracer H20 to Q
  logical, parameter, public :: lnomfix    = .true.     ! do not apply usual mass fixer (eul core)
  logical, parameter, public :: lwtrczmlin = .false.	! linear interp for zm midpoitns (else log)
!
! Set cloud liqiud and ice as advected constituents (match this with cldcond)
!
  logical, parameter :: cldw_adv = .true.	        ! true for advected, false for non-advected
!
! Choose if check should be terminal
!
  logical, parameter :: lcheck_warn_only = .true.       ! true for message only, no endrun

! Water tracer type identifiers
  integer, parameter :: iwtundef = 0     ! not water type
  integer, parameter :: iwtvap   = 1     ! water type is vapour
  integer, parameter :: iwtliq   = 2     ! water type is liquid
  integer, parameter :: iwtice   = 3     ! water type is ice

! Water tracer module data .
! (this is the default set up when wisotope is set true)
  integer, parameter :: nwatr = 3		! 3 types of water (vap, liq, ice)
#ifdef O17HTO
  integer, parameter :: nspec = 5		! 5 species of water (h2o, hdo, o18, o17, hto)
#else
#ifdef TAGGING
! tagged region user must be sure to adjust ntags and wtrc_names together
  integer, parameter :: ntags = 10 !4		! N *species* of water (user chooses)
  integer, parameter :: nspec = 3+3*ntags	! N *species* of water (user chooses)
                                                !    h2o, hdo, o18  globally
                                                !    h2o1, hdo1, o181  in mask region 1
                                                !    h2o2, hdo2, o182  in mask region 2
                                                !    etc
#else
  integer, parameter :: nspec = 3		! 3 species of water (h2o, hdo, o18)
#endif
#endif
  integer, parameter :: ncnst = nwatr*nspec     ! 3 types per species
!
! Tracer vapour names by me no more than 5 characters for history files
  character(len=8), dimension(ncnst), parameter :: & ! constituent names
#ifdef O17HTO
     wtrc_names = (/'H2O   ', 'H2OL  ', 'H2OI  ', &     ! NCNST=3
                    'HDO   ', 'HDOL  ', 'HDOI  ', &     ! NCNST=6
                    'H218O ', 'H218OL', 'H218OI', &     ! NCNST=9
                    'H217O ', 'H217OL', 'H217OI', &     ! NCNST=12
                    'HTO   ', 'HTOL  ', 'HTOI  ' /)     ! NCNST=15
#else
#ifdef TAGGING		/* SORRY, these have to be hard coded.... */
     wtrc_names = (/'H2O   ', 'H2OL  ', 'H2OI  ', &     ! NCNST=3
                    'HDO   ', 'HDOL  ', 'HDOI  ', &     ! NCNST=6
                    'H218O ', 'H218OL', 'H218OI', &     ! NCNST=9
                    'H2O01 ', 'H2OL01', 'H2OI01', &     ! NCNST=12
                    'HDO01 ', 'HDOL01', 'HDOI01', &     ! NCNST=15
                    '18O01 ', '18OL01', '18OI01', &     ! NCNST=18(stop here if ntags=1)
                    'H2O02 ', 'H2OL02', 'H2OI02', &     ! NCNST=21
                    'HDO02 ', 'HDOL02', 'HDOI02', &     ! NCNST=24
                    '18O02 ', '18OL02', '18OI02', &     ! NCNST=27(stop here if ntags=2)
		            'H2O03 ', 'H2OL03', 'H2OI03', &     ! NCNST=30
                    'HDO03 ', 'HDOL03', 'HDOI03', &     ! NCNST=33
		            '18O03 ', '18OL03', '18OI03', &     ! NCNST=36(stop here if ntags=3)
		            'H2O04 ', 'H2OL04', 'H2OI04', &     ! NCNST=39
                    'HDO04 ', 'HDOL04', 'HDOI04', &     ! NCNST=42
		            '18O04 ', '18OL04', '18OI04', &     ! NCNST=45(stop here if ntags=4)
		            'H2O05 ', 'H2OL05', 'H2OI05', &     ! NCNST=48
                    'HDO05 ', 'HDOL05', 'HDOI05', &     ! NCNST=51
		            '18O05 ', '18OL05', '18OI05', &     ! NCNST=54(stop here if ntags=5)
		            'H2O06 ', 'H2OL06', 'H2OI06', &     ! NCNST=57
                    'HDO06 ', 'HDOL06', 'HDOI06', &     ! NCNST=60
		            '18O06 ', '18OL06', '18OI06', &     ! NCNST=63(stop here if ntags=6)
		            'H2O07 ', 'H2OL07', 'H2OI07', &     ! NCNST=66
                    'HDO07 ', 'HDOL07', 'HDOI07', &     ! NCNST=69
		            '18O07 ', '18OL07', '18OI07', &     ! NCNST=72(stop here if ntags=7)
		            'H2O08 ', 'H2OL08', 'H2OI08', &     ! NCNST=75
                    'HDO08 ', 'HDOL08', 'HDOI08', &     ! NCNST=78
		            '18O08 ', '18OL08', '18OI08', &     ! NCNST=81(stop here if ntags=8)
		            'H2O09 ', 'H2OL09', 'H2OI09', &     ! NCNST=84
                    'HDO09 ', 'HDOL09', 'HDOI09', &     ! NCNST=87
		            '18O09 ', '18OL09', '18OI09', &     ! NCNST=90(stop here if ntags=5)
		            'H2O10 ', 'H2OL10', 'H2OI10', &     ! NCNST=93
                    'HDO10 ', 'HDOL10', 'HDOI10', &     ! NCNST=96
		            '18O10 ', '18OL10', '18OI10' /)     ! NCNST=99(stop here if ntags=10)
		    
#else
     wtrc_names = (/'H2O   ', 'H2OL  ', 'H2OI  ', &     ! NCNST=3
                    'HDO   ', 'HDOL  ', 'HDOI  ', &     ! NCNST=6
                    'H218O ', 'H218OL', 'H218OI' /)     ! NCNST=9
#endif
#endif


  integer, public :: &
       ixh2oq  , ixh2ol  , ixh2oi  , &! H2O   vap, liq, ice tracer indicies
       ixhdoq  , ixhdol  , ixhdoi  , &! HDO   vap, liq, ice tracer indicies
       ixh218oq, ixh218ol, ixh218oi, &! H218O vap, liq, ice tracer indicies
       ixh217oq, ixh217ol, ixh217oi, &! H217O vap, liq, ice tracer indicies
       ixhtoq  , ixhtol  , ixhtoi     ! HTO   vap, liq, ice tracer indicies
#ifdef TAGGING
   real(r8), public, allocatable :: tagmask3d(:,:,:) 
   real(r8), public, allocatable :: tagmask(:,:) 
   integer, public :: &
       ixtagq(nspec), ixtagl(nspec), ixtagi(nspec) ! (HDO) tracer vap, liq, ice indices
#endif
!
! Configuration pointers/indicies
  integer, public :: iwater(pcnst+pnats)     ! flag for water type:
                                        ! 1 vapour, 2 liquid, 3 ice
                                        ! 0 not water
  integer, public :: iwspec(pcnst+pnats)     ! flag for water (isotope) species
                                        ! see water_isotopes for definitions
!
! Index arrays for all specified water tracers
!
  integer, public :: iavap(nspec)  ! index arrays for vapour
  integer, public :: ialiq(nspec)  ! index arrays for cloud liquid
  integer, public :: iaice(nspec)  ! index arrays for cloud ice
!
  integer, public :: ixwti, ixwtx     ! lowest and highest index to search

!-----------------------------------------------------------------------
!
! Default minimum difference to trigger check failure
!
!  real(r8) :: qchkmin = 1.e-15		! loss of 3 s.f.
  real(r8) :: qchkmin = 1.e-16		! loss of 2 s.f.
!  real(r8) :: qchkmin = 1.e-17		! loss of 1 s.f.
!  real(r8) :: qchkmin = 1.e-18	  	! qmin (model sees anything less as zero
!  real(r8) :: qchkmin = 1.e-19
!  real(r8) :: qchkmin = 1.e-20
!  real(r8) :: qchkmin = 1.e-21		! very strict
!
!-----------------------------------------------------------------------
contains


!=======================================================================
  subroutine wtrc_init
!-----------------------------------------------------------------------
!
! Purpose: initialize water_tracer parameterizations and indexing
!          (declare additional history field)
!
! Method:
!   Set up water indexing scheme (which must be done AFTER tracers
!   have been registered). Also, set up indexing for the prognostic
!   waters, just for completeness, although are not used.
!   Calls initialization of water isotope module.
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 18:01:52 MDT 2003
!
!-----------------------------------------------------------------------
!!  use history, only: add_fld, add_default
  use water_isotopes, only: wiso_init
!-----------------------------------------------------------------------
!
! We must have the isotopes set up
!
    if (.not. wisotope) then
      write(6,*) 'WTRC_INIT: Water tracers require water isotopes.'
      call endrun
    end if
    write(6,*) 'WTRC_INIT: Initializing water tracers.'
!
! Initialize isotope module
!
     call wiso_init
!
! Add further diagnostics to history file
! (including snow pack, precipitation, bi-directional fluxes)
!
!!    call add_fld()
!!    call add_default
!
#ifdef TAGGING	
   call tag_mask_init
#endif
    return
  end subroutine wtrc_init

!=======================================================================
  subroutine wtrc_register
!-----------------------------------------------------------------------
!
! Purpose: resister advected water tracer constituents
!
! Method:
!  Calls CAM constituent registration routines based on
!  water tracer species and phase indexing.
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 15:31:56 MDT 2003
!
!-----------------------------------------------------------------------
    use physconst,      only: mwdry, cpair, mwh2o, cph2o
    use constituents,   only: advected, nonadvec, cnst_get_ind
    use water_isotopes, only: ispundef, isph2o, isphdo, isph218o, isph217o, isphto

    integer ntrip	 ! accumulator for number of water species (register check)
    integer ntag

    integer ixwprg       ! constituent index of prognostic
    integer flag

    character(len=2) tno ! tag number
!-----------------------------------------------------------------------
!
! Initialize all tracers as nonwater, with unknwon species
!
      ntrip = 0
      iwater(:) = iwtundef
      iwspec(:) = ispundef
!
! Set the species of the total water as H2O, but DONT set them
! as water tracers, the they are prognostic (not "tracers")
!
     ixwti = 0
     ixwtx = 0
!
     call cnst_get_ind('Q     ', ixwprg)
     iwspec(ixwprg) = isph2o
     ixwti = max(ixwtx,ixwprg)

     call cnst_get_ind('CLDLIQ', ixwprg)
     iwspec(ixwprg) = isph2o
     ixwti = max(ixwtx,ixwprg)

     call cnst_get_ind('CLDICE', ixwprg)
     iwspec(ixwprg) = isph2o
     ixwti = max(ixwtx,ixwprg)

     ixwti = ixwti + 1		! (SHOULD BE 4 FOR CAM2.x)

! Set names of variable tendencies and declare them as history variables
! Notice, following CAM code, these can be advected or non-advected.

    if (wisotope) then

      if (cldw_adv) then
         flag = advected
      else
         flag = nonadvec
     endif

! H2O
      call wtrc_cnst_add('H2O'   ,advected, iwtvap, isph2o, mwh2o, cph2o, 0._r8, &
            ixh2oq, longname='H2O tracer specific humidity')
      call wtrc_cnst_add('H2OL'  ,flag    , iwtliq, isph2o, mwdry, cpair, 0._r8, &
            ixh2ol, longname='H2O tracer grid box avg. liquid condensate amount')
      call wtrc_cnst_add('H2OI'  ,flag    , iwtice, isph2o, mwdry, cpair, 0._r8, &
            ixh2oi, longname='H2O tracer grid box avg. ice condensate amount')
!
      ntrip = ntrip + 1
      iavap(ntrip) = ixh2oq
      ialiq(ntrip) = ixh2ol
      iaice(ntrip) = ixh2oi
!HDO
      call wtrc_cnst_add('HDO'   ,advected, iwtvap, isphdo, mwh2o, cph2o, 0._r8, &
            ixhdoq, longname='HDO tracer specific humidity')
      call wtrc_cnst_add('HDOL'  ,flag    , iwtliq, isphdo, mwdry, cpair, 0._r8, &
            ixhdol, longname='HDO tracer grid box avg. liquid condensate amount')
      call wtrc_cnst_add('HDOI'  ,flag    , iwtice, isphdo, mwdry, cpair, 0._r8, &
            ixhdoi, longname='HDO tracer grid box avg. ice condensate amount')
!
      ntrip = ntrip + 1
      iavap(ntrip) = ixhdoq
      ialiq(ntrip) = ixhdol
      iaice(ntrip) = ixhdoi
!H218O
      call wtrc_cnst_add('H218O' ,advected, iwtvap, isph218o, mwh2o, cph2o, 0._r8, &
            ixh218oq, longname='H218O tracer specific humidity')
      call wtrc_cnst_add('H218OL',flag    , iwtliq, isph218o, mwdry, cpair, 0._r8, &
            ixh218ol, longname='H218O tracer grid box avg. liquid condensate amount')
      call wtrc_cnst_add('H218OI',flag    , iwtice, isph218o, mwdry, cpair, 0._r8, &
            ixh218oi, longname='H218O tracer grid box avg. ice condensate amount')
!
      ntrip = ntrip + 1
      iavap(ntrip) = ixh218oq
      ialiq(ntrip) = ixh218ol
      iaice(ntrip) = ixh218oi

#ifdef TAGGING
      ixtagq(:) = 0
      ixtagl(:) = 0
      ixtagi(:) = 0
      ntag = 1
      do while (ntrip < nspec) 
        write(tno,'(i2.2)') ntag
!H2O
        call wtrc_cnst_add('H2O'//tno   ,advected, iwtvap, isph2o, mwh2o, cph2o, 0._r8, &
              ixtagq(ntrip), longname='H2O'//tno//' tracer specific humidity')
        call wtrc_cnst_add('H2OL'//tno  ,flag    , iwtliq, isph2o, mwdry, cpair, 0._r8, &
              ixtagl(ntrip), longname='H2O'//tno//' tracer grid box avg. liquid condensate amount')
        call wtrc_cnst_add('H2OI'//tno  ,flag    , iwtice, isph2o, mwdry, cpair, 0._r8, &
              ixtagi(ntrip), longname='H2O'//tno//' tracer grid box avg. ice condensate amount')
!
        ntrip = ntrip + 1
        iavap(ntrip) = ixtagq(ntrip)
        ialiq(ntrip) = ixtagl(ntrip)
        iaice(ntrip) = ixtagi(ntrip)

!HDO
        call wtrc_cnst_add('HDO'//tno   ,advected, iwtvap, isphdo, mwh2o, cph2o, 0._r8, &
              ixtagq(ntrip), longname='HDO'//tno//' tracer specific humidity')
        call wtrc_cnst_add('HDOL'//tno  ,flag    , iwtliq, isphdo, mwdry, cpair, 0._r8, &
              ixtagl(ntrip), longname='HDO'//tno//' tracer grid box avg. liquid condensate amount')
        call wtrc_cnst_add('HDOI'//tno  ,flag    , iwtice, isphdo, mwdry, cpair, 0._r8, &
              ixtagi(ntrip), longname='HDO'//tno//' tracer grid box avg. ice condensate amount')
!
        ntrip = ntrip + 1
        iavap(ntrip) = ixtagq(ntrip)
        ialiq(ntrip) = ixtagl(ntrip)
        iaice(ntrip) = ixtagi(ntrip)

!H218O
        call wtrc_cnst_add('18O'//tno   ,advected, iwtvap, isph218o, mwh2o, cph2o, 0._r8, &
              ixtagq(ntrip), longname='H218O'//tno//' tracer specific humidity')
        call wtrc_cnst_add('18OL'//tno  ,flag    , iwtliq, isph218o, mwdry, cpair, 0._r8, &
              ixtagl(ntrip), longname='H218O'//tno//' tracer grid box avg. liquid condensate amount')
        call wtrc_cnst_add('18OI'//tno  ,flag    , iwtice, isph218o, mwdry, cpair, 0._r8, &
              ixtagi(ntrip), longname='H218O'//tno//' tracer grid box avg. ice condensate amount')
!
        ntrip = ntrip + 1
        iavap(ntrip) = ixtagq(ntrip)
        ialiq(ntrip) = ixtagl(ntrip)
        iaice(ntrip) = ixtagi(ntrip)

        ntag = ntag + 1
      end do

!    write(6,*) 'WTRC_REGISTER: ixtagq', ixtagq(:)
!    write(6,*) 'WTRC_REGISTER: ixtagl', ixtagl(:)
!    write(6,*) 'WTRC_REGISTER: ntrip,ntag', ntrip,ntag

#else
!H217O - (PCNST=15)
    if (ntrip < nspec) then
    call wtrc_cnst_add('H217O' ,advected, iwtvap, isph217o, mwh2o, cph2o, 0._r8, &
          ixh217oq, longname='H217O tracer specific humidity')
    call wtrc_cnst_add('H217OL',flag    , iwtliq, isph217o, mwdry, cpair, 0._r8, &
          ixh217ol, longname='H217O tracer grid box avg. liquid condensate amount')
    call wtrc_cnst_add('H217OI',flag    , iwtice, isph217o, mwdry, cpair, 0._r8, &
          ixh217oi, longname='H217O tracer grid box avg. ice condensate amount')
!
      ntrip = ntrip + 1
      iavap(ntrip) = ixh217oq
      ialiq(ntrip) = ixh217ol
      iaice(ntrip) = ixh217oi
    end if

!HTO - (PCNST=18)
    if (ntrip < nspec) then
    call wtrc_cnst_add('HTO' ,advected, iwtvap, isphto, mwh2o, cph2o, 0._r8, &
          ixhtoq, longname='HTO tracer specific humidity')
    call wtrc_cnst_add('HTOL',flag    , iwtliq, isphto, mwdry, cpair, 0._r8, &
          ixhtol, longname='HTO tracer grid box avg. liquid condensate amount')
    call wtrc_cnst_add('HTOI',flag    , iwtice, isphto, mwdry, cpair, 0._r8, &
          ixhtoi, longname='HTO tracer grid box avg. ice condensate amount')
!
      ntrip = ntrip + 1
      iavap(ntrip) = ixhtoq
      ialiq(ntrip) = ixhtol
      iaice(ntrip) = ixhtoi
    end if
!
!/* tagging hdo */
#endif
    else
      write(6,*) 'WTRC_REGISTER: unknown water tracer configuration.'
      call endrun
    endif
!
! Check registry and modyule dimensions
!
    if (ntrip /= nspec) then
      write(*,*) '(WTRC_REGISTER): Number of registered triplets differs from module dimensions.'
      write(*,*) 'NTRIP = ',ntrip,'   NWATR=',nwatr, 'NSPEC=',nspec
      call endrun
    end if
!
! Request space on physics buffer for variables that persist across time steps
!
!!    call pbuf_add()
!
!
! Once the registration is done, report what we actually have just to make sure
!
      call wtrc_setup_diag
      write(*,*) 'WATER TRACERS m=',ixwti,ixwtx ! ixwti = 4 ixwtx = 102
    write(6,*) 'WTRC_REGISTER: done.'
!
    return
  end subroutine wtrc_register

!=======================================================================
  function wtrc_is_wtrc(m)
!-----------------------------------------------------------------------
! Returns true if tracer is vapour
!-----------------------------------------------------------------------
  integer, intent(in) :: m              ! constituent index
  logical wtrc_is_wtrc
!-----------------------------------------------------------------------
!!    wtrc_is_wtrc = wtrc_is_vap(m) .or. wtrc_is_liq(m) .or.  wtrc_is_ice(m)
    wtrc_is_wtrc = .false.
    if (iwater(m) /= iwtundef) wtrc_is_wtrc = .true.
  return
  end function wtrc_is_wtrc

!=======================================================================
  function wtrc_is_vap(m)
!-----------------------------------------------------------------------
! Returns true if tracer is vapour
!-----------------------------------------------------------------------
  integer, intent(in) :: m		! constituent index
  logical wtrc_is_vap
!-----------------------------------------------------------------------
    wtrc_is_vap = .false.
    if (iwater(m) == iwtvap) wtrc_is_vap = .true.
  return
  end function wtrc_is_vap

!=======================================================================
  function wtrc_is_liq(m)
!-----------------------------------------------------------------------
! Returns true if tracer is cloud liquid
!-----------------------------------------------------------------------
  integer, intent(in) :: m		! constituent index
  logical wtrc_is_liq
!-----------------------------------------------------------------------
    wtrc_is_liq = .false.
    if (iwater(m) == iwtliq) wtrc_is_liq = .true.
  return
  end function wtrc_is_liq

!=======================================================================
  function wtrc_is_ice(m)
!-----------------------------------------------------------------------
! Returns true if tracer is cloud ice
!-----------------------------------------------------------------------
  integer, intent(in) :: m		! constituent index
  logical wtrc_is_ice
!-----------------------------------------------------------------------
    wtrc_is_ice = .false.
    if (iwater(m) == iwtice) wtrc_is_ice = .true.
  return
  end function wtrc_is_ice

!=======================================================================
  subroutine wtrc_cnst_add(name, type, iwt, isp, mwc, cpc, qminc, ind, &
                           longname, readiv, mixtype)
!-----------------------------------------------------------------------
! Purpose: provide a wrapper for cnst_add with added index condifuration
!          for more details registration of water tracers
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 21:02:25 MDT 2003
!-----------------------------------------------------------------------
  use constituents, only: cnst_add
!---------------------------- Arguments --------------------------------
    character(len=*), intent(in) :: &
       name      ! constituent name for variable name in history file(8 char max)
    character(len=*), intent(in), optional :: &
       longname  ! long_name attribute in netcdf output (128 char max) [name]
    logical,          intent(in), optional :: &
       readiv    ! true => read initial values from initial file (default: true)
    character(len=*),         intent(in), optional :: &
       mixtype    ! mixing ratio type (dry, wet)

    integer, intent(in)    :: type   ! flag indicating advected or nonadvected
    integer, intent(in)    :: iwt    ! water type indicator
    integer, intent(in)    :: isp    ! water species indicator
    real(r8),intent(in)    :: mwc    ! const. molecular weight (kg/kmol)
    real(r8),intent(in)    :: cpc    ! const. spcfic heat  const press (J/kg/K)
    real(r8),intent(in)    :: qminc  ! minimum  mass mixing ratio (kg/kg)
!                                        normally 0., except water 1.E-12, for
!                                        radiation.

    integer, intent(out)   :: ind    ! global constituent index (in q array)

!-----------------------------------------------------------------------
!
! Pass aruments on to normal code
!
    call cnst_add(name, type, mwc, cpc, qminc, ind, longname, readiv, mixtype)
!
! Knowing the tracer index assign water type and species
!
    iwater(ind) = iwt
    iwspec(ind) = isp
    ixwtx = max(ixwtx,ind)
!
    return
  end subroutine wtrc_cnst_add

!=======================================================================
  function wtrc_implements_cnst(name)
!-----------------------------------------------------------------------
!
! Purpose: return true if specified constituent is implemented by this package
! Notice wtrc_names should be the same as hard coded calls to cnst_resister.
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 16:10:29 MDT 2003
!
!-----------------------------------------------------------------------
     implicit none
!-----------------------------Arguments---------------------------------
     character(len=*), intent(in) :: name   ! constituent name
     logical :: wtrc_implements_cnst        ! return value
!---------------------------Local workspace-----------------------------
     integer :: m
!-----------------------------------------------------------------------
     wtrc_implements_cnst = .false.
     do m = 1, ncnst
        if (name == wtrc_names(m)) then
           wtrc_implements_cnst = .true.
           return
        end if
     end do
     return
  end function wtrc_implements_cnst

!=======================================================================
  subroutine wtrc_init_cnst(name, qwtrc_tmp, q)
!-----------------------------------------------------------------------
!
! Initializes water tracers if not read from initial conditions file.
! Assign as some standard mass fraction of the prognostic waters
! (which  can assumed are set correctly). If using water isotope
! Set the standard ration, else set zero
!
!
! Author: David Noone <dcn@caltech.edu> - Sun Jun 29 18:24:57 MDT 2003
!
!-----------------------------------------------------------------------
    use water_isotopes, only: wisotope, wiso_get_rstd
    use constituents, only: pcnst, pnats, cnst_get_ind
!---------------------------- Arguments --------------------------------
    character(len=*),intent(in)  :: name                ! tracer name
!    real(r8),        intent(in)  :: q3_tmp(plond,plev,pcnst+pnats,plat) ! cam2
    real(r8),        intent(in)  :: qwtrc_tmp(plon,plev,plat,pcnst+pnats) ! cam3
    real(r8),        intent(out) :: q(plon,plev,plat)   ! mass mixing ratio
!------------------------- Local Variables -----------------------------
    integer ixwtrc              ! index of water tracer
    integer ixwprg              ! intext of water prognostic
    real(r8) rat                ! an isotope ratio
!-----------------------------------------------------------------------
!
! Retrieve the tracer index, and work out index of equivilent prognostic
!
   call cnst_get_ind(name, ixwtrc)      ! this SHOULD be m in calling routine
!
   if (.not. wtrc_is_wtrc(ixwtrc)) then
      call endrun( 'WTRC_INIT_CNST: non water tracer detected.')
   else if (wtrc_is_vap(ixwtrc)) then    ! vapour
     call cnst_get_ind('Q     ', ixwprg)
   else if (wtrc_is_liq(ixwtrc)) then    ! liquid
     call cnst_get_ind('CLDLIQ', ixwprg)
   else if (wtrc_is_ice(ixwtrc)) then    ! ice
     call cnst_get_ind('CLDICE', ixwprg)
   else
      call endrun('WTRC_INIT_CNST: water tracer set as unknown water type.')
   end if
!
! Assign tracer to be total, scaled by some standard ratio
!
    if (wisotope) then
      rat = wiso_get_rstd(iwspec(ixwtrc))
    else
!      rat = 0.
      rat = 1.
    endif
!
    q(1:plon,1:plev,1:plat) = rat*qwtrc_tmp(1:plon,1:plev,1:plat,ixwprg)
!
    return
  end subroutine wtrc_init_cnst

!=======================================================================
   subroutine wtrc_init_qpert(qpert)
!-----------------------------------------------------------------------
! Initialize constituent perturbation to something (smow?)
!-----------------------------------------------------------------------
    use water_isotopes, only: wiso_get_rstd

    real(r8), intent(inout) :: qpert(pcols,pcnst+pnats)

    integer m
    real(r8) rat
!-----------------------------------------------------------------------
    do m = ixwti, ixwtx
       qpert(:,m) = 0.
      if (wtrc_is_vap(m)) then

        if (wisotope) then
          rat = wiso_get_rstd(iwspec(m))
        else
!          rat = 0.
          rat = 1.
        endif
!
        qpert(:,m) = rat*qpert(:,1)
 
      end if
     end do
!
     return
   end subroutine wtrc_init_qpert

!=======================================================================
  subroutine  wtrc_setup_diag
!-----------------------------------------------------------------------
! Purpose: Writes configuration of water tracer scheme to standard output.
!-----------------------------------------------------------------------
    use constituents, only: cnst_name
    use water_isotopes, only: wiso_get_fisub, wiso_get_rstd
!------------------------- Local Variables -----------------------------
    integer m
!-----------------------------------------------------------------------

    write(6,*) ' '
    write(6,*) '---- Water isotopes tracer configurtaion ----'
    write(6,*) 'name      W  S  f     Rstd'
    do m = 1, pcnst+pnats
      if (wtrc_is_wtrc(m)) then
      write(6,1) cnst_name(m),iwater(m), iwspec(m),  &
             int(wiso_get_fisub(iwspec(m))), wiso_get_rstd(iwspec(m))
      end if
 1    format(a8,' ', i3,i3,i3,e16.5)
    end do
    write(6,*) '---------------------------------------------'
    write(6,*) ' '
!
    return
  end subroutine wtrc_setup_diag

!=======================================================================
  subroutine wtrc_qchk3(subr, vname, ncol, q, qmag0)
!-----------------------------------------------------------------------
! Checks that all tracers areRstd*prognostic
! (used for debuggin with no fractionation)
!-----------------------------------------------------------------------
    use water_isotopes, only: wisotope, wiso_get_rstd
    use constituents, only: pcnst, pnats, cnst_get_ind

!---------------------------- Arguments --------------------------------
    character(len=*),intent(in) :: subr   ! name of calling subroutine
    character(len=*),intent(in) :: vname  ! name of variable
    integer , intent(in) :: ncol          ! number of columns to scan
    real(r8), intent(in) :: q(pcols,pver,pcnst+pnats)   ! tarcers
    real(r8), intent(in), optional :: qmag0      ! minimum magnitude of qprg
!------------------------- Local Variables -----------------------------
    real(r8) rstd
    integer ixvap,ixliq,ixice
    integer mvap, mliq, mice
    integer m
!-----------------------------------------------------------------------

    call cnst_get_ind('Q'     , ixvap)
    call cnst_get_ind('CLDLIQ', ixliq)
    call cnst_get_ind('CLDICE', ixice)
!
    do m = ixwti,ixwtx
      if (wtrc_is_vap(m)) then
        mvap = m
        mliq = m + 1
        mice = m + 2
        rstd = wiso_get_rstd(iwspec(m))
        write(*,'(a40,3i6,g16.6)') 'WTRC_QCHK3 ('//trim(subr)//') - tracers:',mvap,mliq,mice, rstd
        call wtrc_qchk2(subr,trim(vname)//'_v',ncol,q(:,:,mvap),rstd*q(:,:,ixvap),qmag0)
        call wtrc_qchk2(subr,trim(vname)//'_l',ncol,q(:,:,mliq),rstd*q(:,:,ixliq),qmag0)
        call wtrc_qchk2(subr,trim(vname)//'_i',ncol,q(:,:,mvap),rstd*q(:,:,ixvap),qmag0)
      end if
    end do
!
    return
  end subroutine wtrc_qchk3

!=======================================================================
  subroutine wtrc_qchk2(subr,vname,ncol,qtrc,qprg,qmag0)
!-----------------------------------------------------------------------
! Purpose: Check the tracer water mass equal the prognostic
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 19:00:15 MDT 2003
!-----------------------------------------------------------------------

!---------------------------- Arguments --------------------------------
    character(len=*),intent(in) :: subr   ! name of calling subroutine
    character(len=*),intent(in) :: vname  ! name of variable
    integer , intent(in) :: ncol          ! number of columns to scan
    real(r8), intent(in) :: qtrc(pcols,pver)   ! tracer water
    real(r8), intent(in) :: qprg(pcols,pver)   ! prognostic water
    real(r8), intent(in), optional :: qmag0      ! minimum magnitude of qprg
!------------------------- Local Variables -----------------------------
    real(r8) etest                        ! test variable
    real(r8) qmag
    real(r8) qdw, etw                   ! worst values
    integer nbad                        ! number of bad values found
    integer i,k
!-----------------------------------------------------------------------
    nbad = 0
!    qmag = 0.
    qmag = qchkmin
    qdw = 0.
    etw = 0.
    if (present(qmag0)) qmag = qmag0
!
    do k = 1, pver
      do i = 1, ncol
       if (wtrc_qchk_one(qtrc(i,k),qprg(i,k),etest,qmag) > 0) then
#ifdef QCHKTERM
          write(6,1) 'WTRC_QCHK2: '//'('//trim(subr)//'.'//trim(vname)//')'// &
          ' q(m,1):',i,k,qtrc(i,k),qprg(i,k),etest
1       format(a36,2i4,2e12.3,e12.4)
#endif
         etw = max(etw, abs(etest))
         qdw = max(qdw, abs(qtrc(i,k)-qprg(i,k)))
         nbad = nbad + 1
       end if
      end do
    end do
!
    if (nbad /= 0) then
        write(6,*) 'WTRC_QCHK2: '//'('//trim(subr)//'.'//trim(vname)//')',&
              ' *** WARNING - chunk tracers /= Q =',nbad,etw,qdw
#ifdef QCHKTERM 
        call endrun('QCHK2 failed.')
#endif
    else

#ifdef QCHKMESS		/* print a sucess message */
        write(6,*) 'WTRC_QCHK2: '//'('//trim(subr)//'.'//trim(vname)//')',&
              ' All OK.'
#endif
    end if
!
    return
  end subroutine wtrc_qchk2

!=======================================================================
  subroutine wtrc_qchk1(subr,vname,ncol,qtrc,qprg,qmag0)
!-----------------------------------------------------------------------
! Purpose: Check the tracer water mass equal the prognostic
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 19:00:15 MDT 2003
!-----------------------------------------------------------------------
!---------------------------- Arguments --------------------------------
    character(len=*),intent(in) :: subr   ! name of calling subroutine
    character(len=*),intent(in) :: vname  ! name of variable
    integer , intent(in) :: ncol	  ! number of columns to scan
    real(r8), intent(in) :: qtrc(pcols)   ! tracer water
    real(r8), intent(in) :: qprg(pcols)   ! prognostic water
    real(r8), intent(in),optional :: qmag0 ! minimum q for fail test
!------------------------- Local Variables -----------------------------
    real(r8) etest                        ! test variable
    real(r8) qmag
    real(r8) qdw, etw                   ! worst values
    integer nbad 		    	  ! number of bad values found
    integer i
!-----------------------------------------------------------------------
    nbad = 0
!    qmag = 0.
    qmag = qchkmin
    qdw = 0.
    etw = 0.
    if (present(qmag0)) qmag = qmag0
!
    do i = 1, ncol
       if (wtrc_qchk_one(qtrc(i),qprg(i),etest,qmag) > 0) then
#ifdef QCHKTERM
          write(6,1) 'WTRC_QCHK1: '//'('//trim(subr)//'.'//trim(vname)//')'// &
          ' q(m,1):',i,qtrc(i),qprg(i),etest
1       format(a40,i4,'    ',2e12.4,e12.4)
#endif
         etw = max(etw, abs(etest))
         qdw = max(qdw, abs(qtrc(i)-qprg(i)))
         nbad = nbad + 1
       end if
    end do
!
    if (nbad /= 0) then
        write(6,*) 'WTRC_QCHK1: '//'('//trim(subr)//'.'//trim(vname)//')',&
              ' *** WARNING - chunk tracers /= Q =',nbad,etw,qdw
#ifdef QCHKTERM               /* terminate */
        call endrun('QCHK1 failed.')
#endif
    else
#ifdef QCHKMESS		/* print a sucess message */
        write(6,*) 'WTRC_QCHK1: '//'('//trim(subr)//'.'//trim(vname)//')',&
              ' All OK.'
#endif
    end if
!
    return
  end subroutine wtrc_qchk1


!=======================================================================
  function wtrc_qchk_one(qtrc,qprg,etest,qmag)
!-----------------------------------------------------------------------
! Purpose: Check the one tracer water mass equal the prognostic
! Author: David Noone <dcn@caltech.edu> - Mon Jun 30 19:00:15 MDT 2003
!-----------------------------------------------------------------------
!    real(r8), parameter :: elimit = 1.0e-16 ! precision required
!    real(r8), parameter :: elimit = 1.0e-14 ! precision required (q1q2 fails)
!    real(r8), parameter :: elimit = 1.0e-12 ! precision required
    real(r8), parameter :: elimit = 1.0e-10 ! precision required
!    real(r8), parameter :: qmin = 1.0e-18 ! precision required
    real(r8), parameter :: qmin = 1.e-26 ! precision required
!---------------------------- Arguments --------------------------------
    real(r8), intent(in) :: qtrc          ! tracer water
    real(r8), intent(in) :: qprg          ! prognostic water
    real(r8), intent(in) :: qmag          ! difference limit
    real(r8), intent(out) :: etest        ! test variable
    real(r8)              qdiff,qmabs
    integer wtrc_qchk_one		  ! return value
!-----------------------------------------------------------------------
    wtrc_qchk_one = 0

#ifdef NOCHECK
!
! By-pass all checking if running with fractionation
!
#else
    qmabs = max(abs(qprg),abs(qtrc))
    qdiff = abs(qtrc - qprg)
    if (qmabs > qmin .and. qdiff > qmag) then
      etest = qdiff / qmabs
      if (etest > elimit) then
         wtrc_qchk_one = 1
#ifdef QCHKTERM		/* if going to fail, write all diagnostics */
         write(*,'(a36,4e16.9)') '(WTRC_QCHK) FAILED TEST:',qtrc,qprg,qdiff,etest
#endif
      end if
    end if
#endif
!
    return
  end function wtrc_qchk_one

!=======================================================================
  subroutine wtrc_chkdelta(subr, ncol, q)
!-----------------------------------------------------------------------
! Checks the delta values of a 2d array (lon,lev)
!-----------------------------------------------------------------------
    use constituents, only: cnst_get_ind
    use water_isotopes, only: wiso_delta
!---------------------------- Arguments --------------------------------
    character(len=*), intent(in) :: subr	! name of calling routine/message
    integer , intent(in) :: ncol	  	! number of columns to scan
    real(r8), intent(in) :: q(pcols,pver,pcnst+pnats)	! tracer quantity
    real(r8) del, delbad
    integer i,k,m
    integer ixvap,ixliq,ixice		! prognostic water species
    integer mbase			! prognostic base for tracer m
    integer nbad
    real(r8) qbad
!-----------------------------------------------------------------------
!
  call cnst_get_ind('Q'     , ixvap)
  call cnst_get_ind('CLDLIQ', ixliq)
  call cnst_get_ind('CLDICE', ixice)
!
! Apply appropriate scaling
!
    do m = ixwti, ixwtx
!
      if (wtrc_is_vap(m)) then
         mbase = ixvap
         mbase = ixh2oq		! relative to tracer water
      else if (wtrc_is_liq(m)) then
         mbase = ixliq
         mbase = ixh2ol		! relative to tracer water
      else if (wtrc_is_ice(m)) then
         mbase = ixice
         mbase = ixh2oi		! relative to tracer water
      else
         write(*,*) '(WTRC_CHKDELTA) unknown tracer.'
         call endrun
      end if
!
     delbad = 0.
     qbad = 0.
     nbad = 0.
      do k = 1, pver
         do i = 1, ncol
            del = wiso_delta(iwspec(m), q(i,k,m), q(i,k,mbase))
            if (abs(del) > 1001.) then
              nbad = nbad + 1
              if (abs(del) > abs(delbad)) then
                 qbad = q(i,k,mbase)
                 delbad = del
              endif
!!              call endrun('(wtrc_chkdelta) Stopped.')
            end if
         end do
      end do

!!      if (nbad > 0) then
!!        write(*,*) trim(subr)//' Bad delta values for m=',m
!!        write(*,*) 'nbad = ',nbad, '  worst=',delbad,qbad
!!      end if
    end do
    return
  end subroutine wtrc_chkdelta

!=======================================================================
  subroutine wtrc_check(subr, ncol, q)
!-----------------------------------------------------------------------
! Checks H2O tracer (ice, liquid and vapour) is the same as the prognostic
! (optioanllly adjust)
!-----------------------------------------------------------------------
    use constituents, only: cnst_get_ind, qmin
    use water_isotopes, only: wiso_delta
!---------------------------- Arguments --------------------------------
    character(len=*), intent(in) :: subr        ! name of calling routine/message
    integer , intent(in) :: ncol                ! number of columns to scan
    real(r8), intent(inout) :: q(pcols,pver,pcnst+pnats) ! tracer quantity (optionally scaled)
    real(r8) dvap(pcols,pver)
    real(r8) dliq(pcols,pver)
    real(r8) dice(pcols,pver)
    integer i,k,m
    integer iw(2)			! indices of worst values
    integer ixvap,ixliq,ixice           ! prognostic water species
    integer mbase                       ! prognostic base for tracer m
    logical lerrors
!-----------------------------------------------------------------------
    real(r8) :: qtol = 1.e-12		! qmin(Q) = 1.e-12, qmin(L,I) = 0.)
!    real(r8) :: qtol = 1.e-15		! qmin(Q) = 1.e-12, qmin(L,I) = 0.)
!    real(r8) :: qtol = 1.e-17		! Too strict for zm_evap (review?)
!-----------------------------------------------------------------------
!
    lerrors = .false.
    call cnst_get_ind('Q'     , ixvap)
    call cnst_get_ind('CLDLIQ', ixliq)
    call cnst_get_ind('CLDICE', ixice)
!
!  Get all differences
!
    do k = 1, pver
      do i = 1, ncol
        dvap(i,k) = q(i,k,ixvap) - q(i,k,ixh2oq)
        dliq(i,k) = q(i,k,ixliq) - q(i,k,ixh2ol)
        dice(i,k) = q(i,k,ixice) - q(i,k,ixh2oi)
      end do
    end do
!
! Send reports
!
    if (.not. lcheck_warn_only) then
      do k = 1, pver
        do i = 1, ncol
          if (abs(dvap(i,k)) > max(qtol,qmin(ixvap))) then
            write(*,1) i,k,q(i,k,ixvap), q(i,k,ixh2oq), dvap(i,k)
1           format (2i3,3e20.10)
          end if
        end do
      end do
    end if
!
    if (count(abs(dvap(:ncol,:))>max(qtol,qmin(ixvap))) > 0) then
       write(*,*) '(wtrc_check) vapour differences: '//trim(subr)
       iw = maxloc(abs(dvap(:ncol,:)))
       write(*,2) count(abs(dvap(:ncol,:))>qmin(ixvap)),iw(1),iw(2),dvap(iw(1),iw(2)), q(iw(1),iw(2),ixh2oq)
       if (.not. lcheck_warn_only) &
           call endrun('wtrc_check: vapour check failed.')
    endif
    if (count(abs(dliq(:ncol,:))>max(qtol,qmin(ixliq))) > 0) then
       write(*,*) '(wtrc_check) liquid differences: '//trim(subr)
       iw = maxloc(abs(dliq(:ncol,:)))
       write(*,2) count(abs(dliq(:ncol,:))>qtol),iw(1),iw(2),dliq(iw(1),iw(2)), q(iw(1),iw(2),ixh2ol)
       if (.not. lcheck_warn_only) &
           call endrun('wtrc_check: cloud liquid check failed.')
    endif
    if (count(abs(dice(:ncol,:))>max(qtol,qmin(ixice))) > 0) then
       write(*,*) '(wtrc_check) ice differences: '//trim(subr)
       iw = maxloc(abs(dice(:ncol,:)))
       write(*,2) count(abs(dice(:ncol,:))>qtol),iw(1),iw(2),dice(iw(1),iw(2)), q(iw(1),iw(2),ixh2oi)
       if (.not. lcheck_warn_only) &
           call endrun('wtrc_check: cloud ice check failed.')
    endif
 2  format(i4,' point(s), worst (i,k):',2i5,2e16.6)
!
! Do any adjustments
!
    if (lh2oadj) then
!!        write(*,*) 'Applying rescaling to tracers to prohibit drift.'
        call wtrc_rescale(q,ncol)
    end if
!
#ifdef QCHKMESS
   write(*,*) '(wtrc_check) all OK: '//trim(subr)
#endif
!
    return
  end subroutine wtrc_check

!=======================================================================
  subroutine wtrc_rescale(q,ncol)
!-----------------------------------------------------------------------
! Purpose: Ensures tracer water mass is exactly the same as prognostic
! Author: David Noone <dcn@caltech.edu> - Mon Mar  8 16:22:30 PST 2004
!-----------------------------------------------------------------------
    use constituents,   only: cnst_get_ind
!---------------------------- Arguments --------------------------------
  integer, intent(in)     :: ncol
  real(r8), intent(inout) :: q(pcols,pver,pcnst+pnats)
!-----------------------------------------------------------------------
  real(r8) qerr(pcols,pver)
  real(r8) rat
  integer i,k,m
  integer ixvap,ixliq,ixice	! prognostic water species
  integer mbase			! tracer base for tracer m
  integer mprog			! prognostic base for tracer m
!-----------------------------------------------------------------------
!
  call cnst_get_ind('Q'     , ixvap)
  call cnst_get_ind('CLDLIQ', ixliq)
  call cnst_get_ind('CLDICE', ixice)
!
! Apply appropriate scaling
!
    do m = ixwti, ixwtx
!
      if (wtrc_is_vap(m)) then
         mprog = ixvap
         mbase = ixh2oq		! relative to tracer water
      else if (wtrc_is_liq(m)) then
         mprog = ixliq
         mbase = ixh2ol		! relative to tracer water
      else if (wtrc_is_ice(m)) then
         mprog = ixice
         mbase = ixh2oi		! relative to tracer water
      else
         write(*,*) '(WTRC_RESCALE) unknown traaver.'
         call endrun
      end if
!
! Compute the error
!
      qerr(:,:) = q(:,:,mbase) - q(:,:,mprog)
!
! Compute tracer ratio, consistent with tracers, then
! apply to total (prognostic) mass
!
      do k = 1, pver
        do i = 1, ncol
           rat = wtrc_ratio(q(i,k,m), q(i,k,mbase))
!!           q(i,k,m) =  rat*q(i,k,mprog) 		! direct scale
           q(i,k,m) =  q(i,k,m) - rat*qerr(i,k) 	! scale error
        end do
      end do
!
    end do
!
    return
  end subroutine wtrc_rescale

!=======================================================================
  function wtrc_ratio(qtrc,qtot)
!-----------------------------------------------------------------------
! Purpose: Compute tracer ratio from masses, with numerical checks
! Author David Noone <dcn@colorado.edu> - Sat Jul  3 18:52:40 MDT 2004
!-----------------------------------------------------------------------
    real(r8),intent(in)  :: qtrc        ! tracer water mass
    real(r8),intent(in)  :: qtot        ! "base" water mass
    real(r8) :: wtrc_ratio              ! return value
!-----------------------------------------------------------------------
!    real(r8) :: qtiny = 1.e-16		! bigger makes scheme more stable
    real(r8) :: qtiny = 1.e-22		! smaller makes scheme more accurate
!-----------------------------------------------------------------------
    if (abs(qtot) < qtiny) then
       wtrc_ratio = 0.
       return
    end if
!
    if (qtot > 0.) then
      wtrc_ratio = qtrc/(qtot+qtiny)
    else
      wtrc_ratio = qtrc/(qtot-qtiny)
    end if
    return
  end function wtrc_ratio

!=======================================================================
  subroutine wtrc_ratio_all(ncol,q,rat)
!-----------------------------------------------------------------------
! Computes ratios for all water tracers
!-----------------------------------------------------------------------
    use constituents,   only: cnst_get_ind
!---------------------------- Arguments --------------------------------
   integer , intent(in) :: ncol
   real(r8), intent(in) :: q(pcols,pver,pcnst+pnats)
   real(r8), intent(out) :: rat(pcols,pver,pcnst+pnats)
!-----------------------------------------------------------------------
   integer i,k,m
   integer ixvap,ixliq,ixice		! prognostic water species
   integer mbase			! prognostic base for tracer m
!-----------------------------------------------------------------------
   rat(:,:,:) = 0.
!
   call cnst_get_ind('Q'     , ixvap)
   call cnst_get_ind('CLDLIQ', ixliq)
   call cnst_get_ind('CLDICE', ixice)
!
   rat(:,:,ixvap) = 1.
   rat(:,:,ixliq) = 1.
   rat(:,:,ixice) = 1.
!
!  Compute ratios based on "parent"
!
    do m = ixwti, ixwtx
!
      if (wtrc_is_vap(m)) then
         mbase = ixvap
      else if (wtrc_is_liq(m)) then
         mbase = ixliq
      else if (wtrc_is_ice(m)) then
         mbase = ixice
      else
         write(*,*) '(WTRC_RAT_ALL) unknown traaver.'
         call endrun
      end if

      do k = 1, pver
        do i = 1, ncol
           rat(i,k,m) = wtrc_ratio(q(i,k,m), q(i,k,mbase))
        end do
      end do

    end do
!
    return
  end subroutine wtrc_ratio_all

!=======================================================================
#ifdef TAGGING
  subroutine tag_mask_init
!-----------------------------------------------------------------------
! read in the mask for tagging
!-----------------------------------------------------------------------
   use filenames, only: bndtag
   use ioFileMod, only: getfil
   use phys_grid,    only: scatter_field_to_chunk

   implicit none
      
#include <netcdf.inc>

   integer :: ncid_tag   
   integer :: latid, lonid, maskid
   integer :: latsiz, lonsiz
   character(len=256) :: locfn            ! local file
   real(r8) :: tagmask_tmp3d(plon,ntags,plat)   ! dims known at compile time
   real(r8) :: tagmask_tmp(plon,plat)     ! dims known at compile time
   integer nt,i,j

!-----------------------------------------------------------------------
! some dims may not be known at compile time
    allocate (tagmask3d(pcols,ntags,begchunk:endchunk))  
    allocate (tagmask(pcols,begchunk:endchunk))  
    tagmask3d(:,:,:)    = 0.0     ! be sure to initialize in case ncols<pcols
    tagmask(:,:)    = 0.0     ! be sure to initialize in case ncols<pcols

    if ( masterproc ) then

    write (6, '(2x, a)') '_______________________________________________________'
    write (6, '(2x, a)') '__________________ Reading Tagging Mask _______________'
    write (6, '(2x, a)') '_______________________________________________________'

    call getfil(bndtag, locfn)
    call wrap_open(locfn, 0, ncid_tag)
    write(6,*)'INITEXT: NCOPN returns id ',ncid_tag,' for file ',trim(locfn)

    call wrap_inq_dimid(ncid_tag, 'lat', latid)
    call wrap_inq_dimid(ncid_tag, 'lon', lonid)

    call wrap_inq_dimlen(ncid_tag, latid, latsiz)
    call wrap_inq_dimlen(ncid_tag, lonid, lonsiz)
      if (lonsiz /= plon) then
         write(6,*)'TAG_MASK_INIT: lonsiz=',lonsiz,' must = plon=',plon
         call endrun ()
      end if
      if (latsiz  /=  plat) then
         write(6,*)'TAG_MASK_INIT: latsiz=',latsiz,' must = plat=',plat
         call endrun ()
      end if

    write (6, '(2x, a, x, i2)') 'tagging mask initialize : latsiz =', latsiz
    write (6, '(2x, a, x, i2)') 'tagging mask initialize : lonsiz =', lonsiz

    call wrap_inq_varid (ncid_tag, 'MASK',  maskid)
    call wrap_get_var_realx (ncid_tag, maskid,  tagmask_tmp)

!  deconstruct the input mask
    tagmask_tmp3d(:,:,:)=0.0
    do i=1,plon
      do nt=1,ntags
        do j=1,plat
          if ( (tagmask_tmp(i,j) .gt. nt-.99) .and. (tagmask_tmp(i,j) .lt. nt+.01) ) &
             tagmask_tmp3d(i,nt,j)=tagmask_tmp(i,j)-nt+1.0
        enddo
      enddo
    enddo

   end if ! end if ( masterproc )

   call scatter_field_to_chunk(1,ntags,1,plon,tagmask_tmp3d,tagmask3d(1,1,begchunk))
   call scatter_field_to_chunk(1,1,1,plon,tagmask_tmp,tagmask(1,begchunk))
   if (masterproc) then
      write(6,*) 'read tag mask and scattered'
   end if

  end subroutine tag_mask_init
#endif

!=======================================================================

end module water_tracers
