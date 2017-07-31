#include <misc.h>
#include <params.h>

module wiso_histogram
!-----------------------------------------------------------------------
!
! Computes "histogram" diagnostics for output.
! Units are "fraction", with sum of all bins = 1.0
!
! For now only one species (HDO)
!
! David Noone <dcn@colorado.edu> - Wed Oct 20 13:57:20 MDT 2004
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
! Interfaces
!
  public :: wiso_hist_init	! initialize module and history file
  public :: wiso_hist_diag	! compute diagnostics and output
!
! Module parameters and variables
!
  integer, parameter :: pbin     =  20	      ! number of bins
!
  character(len=8)   :: sname    = 'HDO     ' ! name of species
  logical            :: lhistdef = .true.     ! put on primary histoy file by default
  real(r8)           :: delmax   =  0.        ! maximum delta value
  real(r8)           :: delmin   = -1000.     ! minimum delta value
!
  real(r8) bmid(pbin)			      ! bin value at midpoint
  real(r8) bint(pbin+1)			      ! bin balue at interfaces
!
!-----------------------------------------------------------------------
contains

!=======================================================================
subroutine wiso_hist_init()
!-----------------------------------------------------------------------
!
! Initializses histogram statistics on output file
!
!-----------------------------------------------------------------------
  use constituents,    only: pcnst, pnats, qmin
  use history,         only: addfld, add_default, phys_decomp
!-----------------------------------------------------------------------
  character(len=16) brange
  character(len=2) atr
  integer n
!-----------------------------------------------------------------------
  if (.not. lhistdef) then
    write(*,*) 'NO wtare tracer histogram diahnostics.'
    return
  end if
!
! Set the bin boundaries
!
  do n = 1,pbin+1
     bint(n) = delmin + real(n-1)*(delmax - delmin) / real(pbin)
  end do
  bmid(1:pbin) = 0.5*(bint(1:pbin) + bint(2:pbin+1))
  write(*,*) 'WISO_HISTOGRM bin interfaces:',bint
!
! Add a field for each bin
!
  do n = 1, pbin
       write(atr,'(i2.2)') n
       write(brange,'(f8.2,a1,f7.2)') bint(n),':',bint(n+1)
       call addfld('DVAPB'//atr,'fraction ',pver, 'A',  &
          trim(sname)//' vapour frequency bin '//trim(brange),phys_decomp)
       call addfld('DTOTB'//atr,'fraction ',pver, 'A',  &
          trim(sname)//' cloud total water frequency bin '//trim(brange),phys_decomp)

! Optionally add to primary history file by default
       if (lhistdef) then
          call add_default('DVAPB'//atr, 1,' ')
          call add_default('DTOTB'//atr, 1,' ')
       end if
  end do
!
end subroutine wiso_hist_init

!=======================================================================
subroutine wiso_hist_diag(lchnk, ncol, state, cld)
!-----------------------------------------------------------------------
!
! Computes histogram diagnostics and outputs
!
!-----------------------------------------------------------------------
  use physics_types,  only: physics_state
  use constituents,   only: cnst_get_ind,ppcnst
  use history,        only: outfld
  use water_tracers,  only: iwspec
  use water_isotopes, only: wiso_delta
!-----------------------------------------------------------------------
  integer, intent(in) :: lchnk
  integer, intent(in) :: ncol
  type(physics_state), intent(in) :: state
  real(r8), intent(in) :: cld(pcols,pver)	! cloud fraction
!-----------------------------------------------------------------------
  character(len=2) atr
  real(r8) q(pcols,pver,ppcnst)
  real(r8) qtott(pcols,pver)	! total cloud total water
  real(r8) qtoti(pcols,pver)	! isotope cloud total water
  real(r8) delvap(pcols,pver)	! delta of vapout
  real(r8) deltot(pcols,pver)	! delta of total (cloud) water
  real(r8) ftmp(pcols,pver)	! frecuency
  integer icliq, icice
  integer ixvap, ixliq, ixice
  integer n
  integer i,k
!-----------------------------------------------------------------------
  if (.not. lhistdef) return	! who cares
!
  q(:ncol,:,:) = state%q(:ncol,:,:)
!
! Obtain the constituent indicies
!
  call cnst_get_ind('CLDLIQ',icliq)
  call cnst_get_ind('CLDICE',icice)
!
  call cnst_get_ind(trim(sname)     ,ixvap)
  call cnst_get_ind(trim(sname)//'L',ixliq)
  call cnst_get_ind(trim(sname)//'I',ixice)
!
! Compute the total water in the cloudy region
!
  qtott(:ncol,:) = cld(:ncol,:)*q(:ncol,:,    1) + q(:ncol,:,icliq) + q(:ncol,:,icice)
  qtoti(:ncol,:) = cld(:ncol,:)*q(:ncol,:,ixvap) + q(:ncol,:,ixliq) + q(:ncol,:,ixice)
!
! Compute the delta values of the environment vapour and total water
!
  do k = 1, pver
    do i = 1, ncol
      delvap(i,k) = wiso_delta(iwspec(ixvap), q    (i,k,ixvap), q    (i,k,1))
      deltot(i,k) = wiso_delta(iwspec(ixvap), qtoti(i,k      ), qtott(i,k  ))
    end do
  end do
!
! Assign masking to delta balue, just for numerics
!
  where (delvap(:ncol,:) < -1000.) delvap(:ncol,:) = -1000.
  where (delvap(:ncol,:) > +1000.) delvap(:ncol,:) = +1000.
  where (deltot(:ncol,:) < -1000.) deltot(:ncol,:) = -1000.
  where (deltot(:ncol,:) > +1000.) deltot(:ncol,:) = +1000.
!
! Assign to each bin, writing as we go
!
  do n = 1, pbin
    write(atr,'(i2.2)') n
    ftmp(:,:) = 0.
    where (delvap(:ncol,:)>bint(n) .and. delvap(:ncol,:)< bint(n+1))
       ftmp(:ncol,:) = 1.0
    endwhere
    call outfld('DVAPB'//atr, ftmp , pcols, lchnk)
!
    ftmp(:,:) = 0.
    where (deltot(:ncol,:)>bint(n) .and. deltot(:ncol,:)<bint(n+1))
       ftmp(:ncol,:) = 1.0
    endwhere
    call outfld('DTOTB'//atr, ftmp , pcols, lchnk)
!
  end do
!
  return
end subroutine wiso_hist_diag

!=======================================================================
end module wiso_histogram
