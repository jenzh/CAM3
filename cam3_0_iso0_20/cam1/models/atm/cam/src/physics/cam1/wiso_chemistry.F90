#include <misc.h>
#include <params.h>

module wiso_chemistry
!-----------------------------------------------------------------------
!
! Provides routines to handle computation of water isotope specific humidity
! tendancies due to external chemistry. 
!
!
! Author: David Noone <dcn@colorado.edu> - Tue Feb  3 12:50:10 MST 2004
!
!-----------------------------------------------------------------------
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver
!-----------------------------------------------------------------------
  implicit none
  private
  save
!-----------------------------------------------------------------------
!
! Public interfaces
  public :: wiso_chem_init                   ! Initialized module
  public :: wiso_chem_ch4ox_tend             ! methane oxidation (HDO)
  public :: wiso_chem_htodcy_tend            ! radioactive decay (HTO)
  public :: wiso_chem_o3mif_tend             ! mass-independent ozone (H218O, H217O)

!
! [Namelist?] variables 
! 
  logical, public :: lwiso_ch4ox = .false.    ! flag for methane oxidation
  logical, public :: lwiso_o3mif = .false.    ! flag for MI fractionation
!
contains

!=======================================================================
  subroutine wiso_chem_init
!-----------------------------------------------------------------------
! Initializes module: adds field to history tape
!-----------------------------------------------------------------------
    use constituents,   only: cnst_name
    use chemistry,      only: trace_gas
    use water_tracers,  only: ixwti, ixwtx, wtrc_is_vap, iwspec
    use water_isotopes, only: wisotope, isph2o, isphdo
    use history,        only: addfld, add_default, phys_decomp

!---------------------------Local workspace-----------------------------
    integer :: m                                   ! tracer index
!-----------------------------------------------------------------------
  if (.not. wisotope) then
     !没有同位素功能，化学模块没有任何意义
     write(*,*) '(wiso_chemistry_init) Makes not sence if wisotope = .false.'
  end if
!
! Turn on methane oxidation of using trace gas
!
  if (trace_gas) lwiso_ch4ox = .true.

! Add source temms to history file
  call addfld ('CH3D' ,'kg/kg   ',pver, 'A','Monodeuterated mathane',phys_decomp)
  call add_default ('CH3D' , 1, ' ')
!
  do m = ixwti, ixwtx
    if (wtrc_is_vap(m)) then
       call addfld (trim(cnst_name(m))//'SRC','kg/kg/s ',pver, 'A', &
                trim(cnst_name(m))//' source/sink (methane oxidation)',phys_decomp)
       call add_default (trim(cnst_name(m))//'SRC', 1, ' ')
    end if
  end do
!
  end subroutine wiso_chem_init

!=======================================================================
  subroutine wiso_chem_ch4ox_tend (state, ptend, dt)
!-----------------------------------------------------------------------
!
! Purpose:
! Interface to calculate isotopic (HDO) source from
!    parameterized greenhouse gas chemisty (source/sink).
!    This also modifies the water budget.
!
! Method:
!  Based on input tendencies (computed) and mixing ratio for methane
!  use analytic expression from McCarthy et al 2004 for [CH3D]
!      [CH3D] = a * [CH4] + b
!
! Assume CH3T concentration is zero.
! Assume oxygen composition of newly produced water has the same
! isotopic composition as atmospheric oxygen (taken from Bender at al., 1999)
!  (The O2 values SHOULD be revised)
!
! Author: David Noone <dcn@colorado.edu> - Tue Feb  3 14:26:51 MST 2004
!         David Noone <dcn@colorado.edu> - Tue Aug  3 15:04:00 MDT 2004
!         David Noone <dcn@colorado.edu> - Tue Aug 10 13:29:45 MDT 2004
!
!-----------------------------------------------------------------------
    use ppgrid,          only: pcols, pver
    use physics_types,   only: physics_state, physics_ptend
    use physconst,       only: mwdry
    use constituents,    only: cnst_name, cnst_get_ind
    use water_tracers,   only: wtrc_is_vap, ixwti, ixwtx, iwspec, wtrc_ratio
    use water_isotopes,  only: isph2o, isphdo, wiso_get_epsmw, wiso_get_fisub, &
                               wiso_get_rstd, wiso_get_rao2, wiso_get_rnat, wiso_delta
    use history,         only: outfld
!---------------------------- Arguments --------------------------------

    real(r8),            intent(in)    :: dt        ! time step
    type(physics_state), intent(in)    :: state     ! Physics state variables
    type(physics_ptend), intent(inout) :: ptend     ! parameterization tendencies
!
!------------------------- Local Variables -----------------------------

    integer :: i,k
    integer :: m                        ! constituent index
    integer :: ncol			! number of columns
    integer :: lchnk			! chunk index
    integer :: ixch4                    ! index of ich4 = 1

    real(r8):: ch4vmr, ch3dvmr		! volumetic mixing ratios
    real(r8):: ch4mmr, ch3dmmr		! mass mixing ratios
    real(r8):: qch3d(pcols,pver)	! q for ch3d
    real(r8):: qchd(pcols,pver)	        ! q for hd
    real(r8):: Rch3d(pcols,pver)	! isotope ratio CH3D/CH4
    real(r8):: Rvap(pcols,pver)		! isotope ratio of vapour
    real(r8):: Rao2			! isotope ratio of atmospheric oxygen

!--------------------------- Parameters --------------------------------
!!    real(r8), parameter :: fwteq = 1.973	   ! H2O_eq = 1.973 CH4 + H2O
    real(r8), parameter :: mwch4  = 16.            ! molecular weight of CH4
    real(r8), parameter :: mwch3d = 17.            ! molecular weight of CH3D
    real(r8), parameter :: mwh2   = 2.             ! molecular weight of H2
    real(r8), parameter :: mwhd   = 3.             ! molecular weight of HD
    real(r8), parameter :: a1ch3d = 5.16e-4	   ! linear fit CH3D = a CH4 + b
    real(r8), parameter :: b1ch3d = 0.0908e-9      ! linear fit CH3D = a CH4 + b (units ppb)
!!    real(r8), parameter :: a1hd   = -6.32e-5	   ! linear fit HD = a CH4 + b
!!    real(r8), parameter :: b1hd   = 0.297e-9      ! linear fit HD = a CH4 + b (units ppb)
!
! 
    real(r8), parameter :: rmwch4  = mwch4/mwdry   ! ratio ch4 weight to dry air
    real(r8), parameter :: rmwch3d = mwch3d/mwdry  ! ratio ch3d weight to dry air
    real(r8), parameter :: ppb = 1.e+9		   ! parts per billion
!-----------------------------------------------------------------------
    if (.not. lwiso_ch4ox) then
       write(*,*) '(WISO_CHEM_CH4OX_TEND) Not computing methane oxidation term for HDO.'
       write(*,*) '(WISO_CHEM_CH4OX_TEND) YOU HAVE ATTEMPTED TO USE UNTESTED PLACEHOLDER CODE.'
       return
    end if

! set index for methane in state arrays
    
    call cnst_get_ind('CH4',ixch4)
!
    ncol = state%ncol
    lchnk = state%lchnk
!
! Compute the isotope ratio D/H in CH4
!  (constants from McCarthy et al 2004, which used mole, rather than mixing, ratio)
! Notice, consistent with the model, we compute the isotopic ratio from
! mass mixing ratios. 
!
! Notice I use the same scaling as the McCarthy paper, but I'm not
! convinced this is correct.... might be :)
!

    do k = 1, pver
      do i = 1, ncol
!!        h2ommr  = state%q(i,k,1)
!!        h2ovmr  = h2ommr/rmwh2o
!!        hdovmr  = h2ovmr*wtrc_ratio(state%q(i,k,isphdo), state%q(i,k,isph2o))
!
        ch4mmr  = state%q(i,k,ixch4)
        ch4vmr  = ch4mmr/rmwch4
!
! Compute the approximate CH3D
!
        ch3dvmr = a1ch3d*ch4vmr + b1ch3d
        ch3dmmr = ch3dvmr*rmwch3d
        qch3d(i,k) = ch3dmmr
!!        if (i==1) write(*,*) 'RCH3D:',ch3dvmr*1.e+9,ch4vmr*1.e+9	! CORRECT
!!!
!!! Compute the approximate HD, and H2
!!!
!!        h2vmr = a1h2*ch4vmr + b1h2
!!        h2mmr = h2vmr*rmwh2
!!        qh2(i,k) = h2mmr
!!!
!!        hdvmr = a1hd*ch4vmr + b1hd
!!        hdmmr = hdvmr*rmwhd
!!        qhd(i,k) = hdmmr
!!!
!!! Compute the total water from H and D budget
!!!
!!          qhtot(i,k) = h2ovmr + ch4vmr  + h2vmr
!!          qdtot(i,k) = hdovmr + ch3dvmr + hdvmr
!
! Compute ratios for water budget
!
!!        Rch3d(i,k) = wtrc_ratio(ch3dmmr,ch4mmr)
        Rch3d(i,k) = wtrc_ratio(ch3dvmr*ppb,ch4vmr*ppb)		        ! ppb for numerics
!
! Account for fact that there are two posisble substitutions
!
        Rch3d(i,k) = Rch3d(i,k)/wiso_get_fisub(isphdo)
!
! ANeed to have the effect of ch3d on hdo, thus use "equivilent water" 
!  (value should match values in chemistry.F90)
!
        Rch3d(i,k) = Rch3d(i,k)/1.973
!
! Model carries qisotope = q*R, and R = ni/n, which is qi/q*(mm/mmi)
! So we also need to convert to mol fraction
!
!!        Rch3d(i,k) = Rch3d(i,k)*wiso_get_epsmw(isphdo)
!
! Apply appropriate units conversion for model water tracer quantities
! (not needed unless rstd /= rnat)
!
        Rch3d(i,k) = Rch3d(i,k)*wiso_get_rstd(isphdo)/wiso_get_rnat(isphdo)
!
! Convert ch3d mixing ratio to an effective water value
! (thus delta ch3d, ch4 of output works like water in postprocessing)
!
!!        qch3d(i,k) = qch3d(i,k)/wiso_get_fisub(isphdo)
!!        qch3d(i,k) = qch3d(i,k)/1.973
!!        qch3d(i,k) = qch3d(i,k)*wiso_get_rstd(isphdo)/wiso_get_rnat(isphdo)
!!        qch3d(i,k) = qch3d(i,k)*wiso_get_epsmw(isphdo)
!
! At this point, we have now screwed around wit the ratio enough!
!
      end do
!! Checke that after all this, we get McCarthy's numbers... yep!
!!      write(*,*) 'dCHD:',k,ch4vmr*1.e+9,wiso_delta(isphdo,Rch3d(1,k),1.0)
    end do
!
    call outfld('CH3D', qch3d(:,:), pcols, lchnk)
!
! Compute methane oxidation tendency given:
!  dH2O = -1.943 * dCH4 (computed in ghg chem)
!  CH3D = a CH4 + b
!
    do m = ixwti,ixwtx
      if (wtrc_is_vap(m)) then

        ptend%lq(m) = .true.
        Rao2 = wiso_get_Rao2(iwspec(m))		! o2 depletion
        do k = 1, pver
          do i = 1, ncol
             Rvap(i,k) = wtrc_ratio(state%q(i,k,m)*ppb,state%q(i,k,1)*ppb)       ! ppb for numerics
             if (ptend%q(i,k,1) > 0) then
               ptend%q(i,k,m) = Rao2*ptend%q(i,k,1) 		! just copy from prognostic, with o2 composition
             else
               ptend%q(i,k,m) = Rvap(i,k)*ptend%q(i,k,1) 	! water loss no fractionation
             end if
          end do
        end do
!
! Modify HDO source by CH3D ratio (Rao2 = 1 for HDO)
!
        if (iwspec(m) == isphdo) then		! HDO source
          do k = 1, pver
            do i = 1, ncol
               if (ptend%q(i,k,1) > 0) then
                 ptend%q(i,k,m) = Rch3d(i,k)*ptend%q(i,k,m)	! source with ratio of ch3d
               end if
            end do
          end do
        end if
!
        call outfld(trim(cnst_name(m))//'SRC', ptend%q(:,:,m), pcols, lchnk)

      end if
    end do
!
! Check the delta values: remember qtend is TINY!
!
!!    do k = 1, pver
!!      write(*,*) 'dCHD:',k,wiso_delta(isphdo,ptend%q(1,k,7)*ppb,ptend%q(1,k,1)*ppb), &
!!                           wiso_delta(isphdo,Rch3d(1,k),1.0)
!!    end do
!
    return
  end subroutine wiso_chem_ch4ox_tend

!=======================================================================
  subroutine wiso_chem_htodcy_tend(state,ptend,decaytime)
!-----------------------------------------------------------------------
!
! Purpose: 
!   Computes tendancy associated with radioactive decay of tritium 
!   during some time interval (need not be one timestep).
!
! Author: David Noone <dcn@caltech.edu> - Tue Feb  3 12:49:47 MST 2004
!
!-----------------------------------------------------------------------
    use ppgrid,          only: pcols, pver
    use physics_types,   only: physics_state, physics_ptend
    use water_tracers,   only: ixwti, ixwtx, iwspec
    use water_isotopes,  only: isphto, wiso_decay
!---------------------------- Arguments --------------------------------
    type(physics_state), intent(in)    :: state
    type(physics_ptend), intent(inout) :: ptend
    real(r8)           , intent(in)    :: decaytime    ! time for decay (sec)
!------------------------- Local Variables -----------------------------
    integer  :: i,k
    integer  :: m                 ! tracer index
    real(r8) :: dqdcy
!-----------------------------------------------------------------------
!
! Loop over all water tracers, and decay tritium in all phases
!
    do m = ixwti,ixwtx
      if (iwspec(m) == isphto) then 
        do k = 1, pver
          do i = 1, state%ncol
            call wiso_decay(isphto, decaytime, state%q(i,k,m), dqdcy)
            ptend%q(i,k,m) = dqdcy
          end do
        end do
      end if    ! water tracer is a HTO
    end do      ! m water tracers
!
    return
  end subroutine wiso_chem_htodcy_tend


!=======================================================================
  subroutine wiso_chem_o3mif_tend(state,ptend)
!-----------------------------------------------------------------------
!
! Purpose: provide mass independent scrubbing of isotop region in the
! presence of O(1d), in the region of ozone photolysis.
! This is an incredibly small effect O(-6)...
!
! Author: David Noone <dcn@caltech.edu> - Tue Feb  3 12:49:47 MST 2004
!
!-----------------------------------------------------------------------
    use physics_types,   only: physics_state, physics_ptend
    use water_tracers,   only: wtrc_is_vap, iwspec, ixwti, ixwtx
    use water_isotopes,  only: isph2o, isph218o, isph217o
!---------------------------- Arguments --------------------------------
    type(physics_state), intent(in)    :: state
    type(physics_ptend), intent(inout) :: ptend     ! parameterization tendancies
!------------------------- Local Variables -----------------------------
    integer  :: m                               ! tracer index
    real(r8) :: dqozn				! water tendency due to ozone photolysis
!--------------------------- Parameters --------------------------------
    real(r8), parameter :: bmifh218o = 1.       ! o18 in ozone?
    real(r8), parameter :: bmifh217o = 1.       ! o17 in ozone?
!-----------------------------------------------------------------------
    if (.not. lwiso_o3mif) then
       write(*,*) '(WISO_CHEM_O3MIF_TEND) Not computing MI source of H217O and H218O'
       write(*,*) '(WISO_CHEM_O3MIF_TEND) YOU HAVE ATTEMPTED TO USE UNTESTED PLACEHOLDER CODE.'
       return
    end if
!
! Compute water production due to ozone desctruction
! (convert to units of water equivilent)
!
    dqozn = 0.0		! taken from boundary data
    ptend%q(:,:,1) = dqozn
!
! Add in water source to all species of water vapour (with O substitutions)
!  (HOW ARE WE TO DO THIS???)
!
    do m = ixwti,ixwtx
      if (wtrc_is_vap(m)) then
        if      (iwspec(m) == isph2o  ) then
          ptend%q(:,:,m) = ptend%q(:,:,1)
        else if (iwspec(m) == isph218o) then
          ptend%q(:,:,m) = ptend%q(:,:,1)*bmifh218o
        else if (iwspec(m) == isph217o) then 
          ptend%q(:,:,m) = ptend%q(:,:,1)*bmifh217o
        end if    ! water tracer is a HTO
      end if
    end do      ! m water tracers
!
    return
  end subroutine wiso_chem_o3mif_tend


!=======================================================================
end module wiso_chemistry
