#include <misc.h>
#include <preproc.h>

module subgridAveMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: subgridAveMod
!
! !DESCRIPTION:
! Utilities to perfrom subgrid averaging
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varcon, only : spval
  use abortutils, only : endrun

! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: p2c   ! Perfrom an average from pfts to columns
  public :: p2l   ! Perfrom an average from pfts to landunits
  public :: p2g   ! Perfrom an average from pfts to gridcells
  public :: c2l   ! Perfrom an average from columns to landunits
  public :: c2g   ! Perfrom an average from columns to gridcells
  public :: l2g   ! Perfrom an average from landunits to gridcells

  interface p2c
     module procedure p2c_1d
     module procedure p2c_2d
     module procedure p2c_1d_filter
     module procedure p2c_2d_filter
  end interface
  interface p2l
     module procedure p2l_1d
     module procedure p2l_2d
  end interface
  interface p2g
     module procedure p2g_1d
     module procedure p2g_2d
  end interface
  interface c2l
     module procedure c2l_1d
     module procedure c2l_2d
  end interface
  interface c2g
     module procedure c2g_1d
     module procedure c2g_2d
  end interface
  interface l2g
     module procedure l2g_1d
     module procedure l2g_2d
  end interface
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_1d
!
! !INTERFACE:
  subroutine p2c_1d (lbp, ubp, lbc, ubc, parr, carr, p2c_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to columns.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column
    real(r8), intent(in)  :: parr(lbp:ubp)         ! pft array
    real(r8), intent(out) :: carr(lbc:ubc)         ! column array
    character(len=*), intent(in) :: p2c_scale_type ! scale type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: p,c,index              ! indices
    real(r8) :: scale_p2c(lbp:ubp)     ! scale factor for column->landunit mapping
    logical  :: found                  ! temporary for error check
    real(r8) :: sumwt(lbc:ubc)         ! sum of weights
    real(r8), pointer :: wtcol(:)      ! weight of pft relative to column
    integer , pointer :: pcolumn(:)    ! column index of corresponding pft
!------------------------------------------------------------------------

    wtcol    => clm3%g%l%c%p%wtcol
    pcolumn  => clm3%g%l%c%p%column

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(6,*)'p2c_1d error: scale type ',p2c_scale_type,' not supported'
       call endrun()
    end if

    carr(lbc:ubc) = spval
    sumwt(lbc:ubc) = 0.
    do p = lbp,ubp
       if (parr(p) /= spval .and. wtcol(p) /= 0.) then
          c = pcolumn(p)
          if (sumwt(c) == 0.) carr(c) = 0.
          carr(c) = carr(c) + parr(p) * scale_p2c(p) * wtcol(p)
          sumwt(c) = sumwt(c) + wtcol(p)
       end if
    end do
    found = .false.
    do c = lbc,ubc
       if (sumwt(c) > 1.0_r8 + 1.e-6) then
          found = .true.
          index = c
          exit
       else if (sumwt(c) /= 0._r8) then
          carr(c) = carr(c)/sumwt(c)
       end if
    end do
    if (found) then
       write(6,*)'p2c error: sumwt is greater than 1.0 at c= ',index
       call endrun()
    end if

  end subroutine p2c_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_2d
!
! !INTERFACE:
  subroutine p2c_2d (lbp, ubp, lbc, ubc, num2d, parr, carr, p2c_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from landunits to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: parr(lbp:ubp,num2d)   ! pft array
    real(r8), intent(out) :: carr(lbc:ubc,num2d)   ! column array
    character(len=*), intent(in) :: p2c_scale_type ! scale type
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: j,p,c,index            ! indices
    real(r8) :: scale_p2c(lbp:ubp)        ! scale factor for column->landunit mapping
    logical  :: found                  ! temporary for error check
    real(r8) :: sumwt(lbc:ubc)         ! sum of weights
    real(r8), pointer :: wtcol(:)      ! weight of pft relative to column
    integer , pointer :: pcolumn(:)    ! column index of corresponding pft
!------------------------------------------------------------------------

    wtcol    => clm3%g%l%c%p%wtcol
    pcolumn  => clm3%g%l%c%p%column

    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(6,*)'p2c_2d error: scale type ',p2c_scale_type,' not supported'
       call endrun()
    end if

    carr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0.
       do p = lbp,ubp
          if (parr(p,j) /= spval .and. wtcol(p) /= 0.) then
             c = pcolumn(p)
             if (sumwt(c) == 0.) carr(c,j) = 0.
             carr(c,j) = carr(c,j) + parr(p,j) * scale_p2c(p) * wtcol(p)
             sumwt(c) = sumwt(c) + wtcol(p)
          end if
       end do
       found = .false.
       do c = lbc,ubc
          if (sumwt(c) > 1.0_r8 + 1.e-6) then
             found = .true.
             index = c
             exit
          else if (sumwt(c) /= 0._r8) then
             carr(c,j) = carr(c,j)/sumwt(c)
          end if
          if (found) exit
       end do
       if (found) then
          write(6,*)'p2c_2d error: sumwt is greater than 1.0 at c= ',index,' lev= ',j
          call endrun()
       end if
    end do

  end subroutine p2c_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_1d_filter
!
! !INTERFACE:
  subroutine p2c_1d_filter (numfc, filterc, pftarr, colarr)
!
! !DESCRIPTION:
! perform pft to column averaging for single level pft arrays
!
! !USES:
    use clm_varpar, only : max_pft_per_col
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    real(r8), pointer     :: pftarr(:)
    real(r8), pointer     :: colarr(:)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: fc,c,pi,p           ! indices
    integer , pointer :: npfts(:)
    integer , pointer :: pfti(:)
    integer , pointer :: pftf(:)
    real(r8), pointer :: wtcol(:)
!-----------------------------------------------------------------------

    npfts => clm3%g%l%c%npfts
    pfti  => clm3%g%l%c%pfti
    pftf  => clm3%g%l%c%pftf
    wtcol => clm3%g%l%c%p%wtcol

#if (defined CPP_VECTOR)
    do fc = 1,numfc
       c = filterc(fc)
       colarr(c) = 0._r8
    end do
    do pi = 1,max_pft_per_col
!dir$ concurrent
!cdir nodep
       do fc = 1,numfc
          c = filterc(fc)
          if ( pi <=  npfts(c) ) then
             p = pfti(c) + pi - 1
             colarr(c) = colarr(c) + pftarr(p) * wtcol(p)
          end if
       end do
    end do
#else
    do fc = 1,numfc
       c = filterc(fc)
       colarr(c) = 0._r8
       do p = pfti(c), pftf(c)
          colarr(c) = colarr(c) + pftarr(p) * wtcol(p)
       end do
    end do
#endif

  end subroutine p2c_1d_filter

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2c_2d_filter
!
! !INTERFACE:
  subroutine p2c_2d_filter (lev, numfc, filterc, pftarr, colarr)
!
! !DESCRIPTION:
! perform pft to column averaging for multi level pft arrays
!
! !USES:
    use clm_varpar, only : max_pft_per_col

! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lev
    integer , intent(in)  :: numfc
    integer , intent(in)  :: filterc(numfc)
    real(r8), pointer     :: pftarr(:,:)
    real(r8), pointer     :: colarr(:,:)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: fc,c,pi,p,j    ! indices
    integer , pointer :: npfts(:)
    integer , pointer :: pfti(:)
    integer , pointer :: pftf(:)
    real(r8), pointer :: wtcol(:)
!-----------------------------------------------------------------------

    npfts => clm3%g%l%c%npfts
    pfti  => clm3%g%l%c%pfti
    pftf  => clm3%g%l%c%pftf
    wtcol => clm3%g%l%c%p%wtcol

#if (defined CPP_VECTOR)
    do j = 1,lev
       do fc = 1,numfc
          c = filterc(fc)
          colarr(c,j) = 0._r8
       end do
       do pi = 1,max_pft_per_col
!dir$ concurrent
!cdir nodep
          do fc = 1,numfc
             c = filterc(fc)
             if ( pi <=  npfts(c) ) then
                p = pfti(c) + pi - 1
                colarr(c,j) = colarr(c,j) + pftarr(p,j) * wtcol(p)
             end if
          end do
       end do
    end do
#else
    do j = 1,lev
       do fc = 1,numfc
          c = filterc(fc)
          colarr(c,j) = 0._r8
          do p = pfti(c), pftf(c)
             colarr(c,j) = colarr(c,j) + pftarr(p,j) * wtcol(p)
          end do
       end do
    end do
#endif

  end subroutine p2c_2d_filter

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2l_1d
!
! !INTERFACE:
  subroutine p2l_1d (lbp, ubp, lbc, ubc, lbl, ubl, parr, larr, &
       p2c_scale_type, c2l_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to landunits
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft indices
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
    real(r8), intent(in)  :: parr(lbp:ubp)         ! input column array
    real(r8), intent(out) :: larr(lbl:ubl)         ! output landunit array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: p,c,l,index            ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: sumwt(lbl:ubl)         ! sum of weights
    real(r8) :: scale_p2c(lbc:ubc)     ! scale factor for pft->column mapping
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor for column->landunit mapping
    real(r8), pointer :: wtlunit(:)    ! weight of pft relative to landunit
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    integer , pointer :: plandunit(:)  ! landunit of corresponding pft
!------------------------------------------------------------------------

    wtlunit   => clm3%g%l%c%p%wtlunit
    pcolumn   => clm3%g%l%c%p%column
    plandunit => clm3%g%l%c%p%landunit

    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else
       write(6,*)'p2l_1d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if
    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(6,*)'p2l_1d error: scale type ',p2c_scale_type,' not supported'
       call endrun()
    end if

    larr(:) = spval
    sumwt(:) = 0.
    do p = lbp,ubp
       if (parr(p) /= spval .and. wtlunit(p) /= 0.) then
          c = pcolumn(p)
          l = plandunit(p)
          if (sumwt(l) == 0.) larr(l) = 0.
          larr(l) = larr(l) + parr(p) * scale_p2c(p) * scale_c2l(c) * wtlunit(p)
          sumwt(l) = sumwt(l) + wtlunit(p)
       end if
    end do
    found = .false.
    do l = lbl,ubl
       if (sumwt(l) > 1.0_r8 + 1.e-6) then
          found = .true.
          index = l
          exit
       else if (sumwt(l) /= 0._r8) then
          larr(l) = larr(l)/sumwt(l)
       end if
    end do
    if (found) then
       write(6,*)'p2l_1d error: sumwt is greater than 1.0 at l= ',index
       call endrun()
    end if

  end subroutine p2l_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2l_2d
!
! !INTERFACE:
  subroutine p2l_2d(lbp, ubp, lbc, ubc, lbl, ubl, num2d, parr, larr, &
       p2c_scale_type, c2l_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to landunits
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft indices
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: parr(lbp:ubp,num2d)   ! input pft array
    real(r8), intent(out) :: larr(lbl:ubl,num2d)   ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: j,p,c,l,index          ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: sumwt(lbl:ubl)         ! sum of weights
    real(r8) :: scale_p2c(lbc:ubc)     ! scale factor for pft->column mapping
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor for column->landunit mapping
    real(r8), pointer :: wtlunit(:)    ! weight of pft relative to landunit
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    integer , pointer :: plandunit(:)  ! landunit of corresponding pft
!------------------------------------------------------------------------

    wtlunit   => clm3%g%l%c%p%wtlunit
    pcolumn   => clm3%g%l%c%p%column
    plandunit => clm3%g%l%c%p%landunit

    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else
       write(6,*)'p2l_2d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if
    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(6,*)'p2l_2d error: scale type ',p2c_scale_type,' not supported'
       call endrun()
    end if

    larr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0.
       do p = lbp,ubp
          if (parr(p,j) /= spval .and. wtlunit(p) /= 0.) then
             c = pcolumn(p)
             l = plandunit(p)
             if (sumwt(l) == 0.) larr(l,j) = 0.
             larr(l,j) = larr(l,j) + parr(p,j) * scale_p2c(p) * scale_c2l(c) * wtlunit(p)
             sumwt(l) = sumwt(l) + wtlunit(p)
          end if
       end do
       found = .false.
       do l = lbl,ubl
          if (sumwt(l) > 1.0_r8 + 1.e-6) then
             found = .true.
             index = l
             exit
          else if (sumwt(l) /= 0._r8) then
             larr(l,j) = larr(l,j)/sumwt(l)
          end if
       end do
       if (found) exit
    end do
    if (found) then
       write(6,*)'p2l_2d error: sumwt is greater than 1.0 at l= ',index,' j= ',j
       call endrun()
    end if

  end subroutine p2l_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2g_1d
!
! !INTERFACE:
  subroutine p2g_1d(lbp, ubp, lbc, ubc, lbl, ubl, lbg, ubg, parr, garr, &
       p2c_scale_type, c2l_scale_type, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp            ! beginning and ending pft indices
    integer , intent(in)  :: lbc, ubc            ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl            ! beginning and ending landunit indices
    integer , intent(in)  :: lbg, ubg            ! beginning and ending gridcell indices
    real(r8), intent(in)  :: parr(lbp:ubp)       ! input pft array
    real(r8), intent(out) :: garr(lbg:ubg)       ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: p,c,l,g,index          ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_p2c(lbp:ubp)     ! scale factor
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of pfts relative to gridcells
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    integer , pointer :: plandunit(:)  ! landunit of corresponding pft
    integer , pointer :: pgridcell(:)  ! gridcell of corresponding pft
!------------------------------------------------------------------------

    wtgcell   => clm3%g%l%c%p%wtgcell
    pcolumn   => clm3%g%l%c%p%column
    pgridcell => clm3%g%l%c%p%gridcell
    plandunit => clm3%g%l%c%p%landunit

    if (l2g_scale_type == 'unity') then
       do l = lbl,ubl
          scale_l2g(l) = 1.0_r8
       end do
    else
       write(6,*)'p2g_1d error: scale type ',l2g_scale_type,' not supported'
       call endrun()
    end if
    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else
       write(6,*)'p2g_1d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if
    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(6,*)'p2g_1d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    garr(:) = spval
    sumwt(:) = 0.
    do p = lbp,ubp
       if (parr(p) /= spval .and. wtgcell(p) /= 0.) then
          c = pcolumn(p)
          l = plandunit(p)
          g = pgridcell(p)
          if (sumwt(g) == 0.) garr(g) = 0.
          garr(g) = garr(g) + parr(p) * scale_p2c(p) * scale_c2l(c) * scale_l2g(l) * wtgcell(p)
          sumwt(g) = sumwt(g) + wtgcell(p)
       end if
    end do
    found = .false.
    do g = lbg, ubg
       if (sumwt(g) > 1.0_r8 + 1.e-6) then
          found = .true.
          index = g
          exit
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       write(6,*)'p2g_1d error: sumwt is greater than 1.0 at g= ',index
       call endrun()
    end if

  end subroutine p2g_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: p2g_2d
!
! !INTERFACE:
  subroutine p2g_2d(lbp, ubp, lbc, ubc, lbl, ubl, lbg, ubg, num2d, &
       parr, garr, p2c_scale_type, c2l_scale_type, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from pfts to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbp, ubp              ! beginning and ending pft indices
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
    integer , intent(in)  :: lbg, ubg              ! beginning and ending gridcell indices
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: parr(lbp:ubp,num2d)   ! input pft array
    real(r8), intent(out) :: garr(lbg:ubg,num2d)   ! output gridcell array
    character(len=*), intent(in) :: p2c_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: j,p,c,l,g,index        ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_p2c(lbp:ubp)     ! scale factor
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of pfts relative to gridcells
    integer , pointer :: pcolumn(:)    ! column of corresponding pft
    integer , pointer :: plandunit(:)  ! landunit of corresponding pft
    integer , pointer :: pgridcell(:)  ! gridcell of corresponding pft
!------------------------------------------------------------------------

    wtgcell   => clm3%g%l%c%p%wtgcell
    pcolumn   => clm3%g%l%c%p%column
    pgridcell => clm3%g%l%c%p%gridcell
    plandunit => clm3%g%l%c%p%landunit

    if (l2g_scale_type == 'unity') then
       do l = lbl,ubl
          scale_l2g(l) = 1.0_r8
       end do
    else
       write(6,*)'p2g_2d error: scale type ',l2g_scale_type,' not supported'
       call endrun()
    end if
    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else
       write(6,*)'p2g_2d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if
    if (p2c_scale_type == 'unity') then
       do p = lbp,ubp
          scale_p2c(p) = 1.0_r8
       end do
    else
       write(6,*)'p2g_2d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    garr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0.
       do p = lbp,ubp
          if (parr(p,j) /= spval .and. wtgcell(p) /= 0.) then
             c = pcolumn(p)
             l = plandunit(p)
             g = pgridcell(p)
             if (sumwt(g) == 0.) garr(g,j) = 0.
             garr(g,j) = garr(g,j) + parr(p,j) * scale_p2c(p) * scale_c2l(c) * scale_l2g(l) * wtgcell(p)
             sumwt(g) = sumwt(g) + wtgcell(p)
          end if
       end do
       found = .false.
       do g = lbg, ubg
          if (sumwt(g) > 1.0_r8 + 1.e-6) then
             found = .true.
             index = g
             exit
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) exit
    end do
    if (found) then
       write(6,*)'p2g_2d error: sumwt is greater than 1.0 at g= ',index,' j= ',j
       call endrun()
    end if

  end subroutine p2g_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2l_1d
!
! !INTERFACE:
  subroutine c2l_1d (lbc, ubc, lbl, ubl, carr, larr, c2l_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from columns to landunits
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc      ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl      ! beginning and ending landunit indices
    real(r8), intent(in)  :: carr(lbc:ubc) ! input column array
    real(r8), intent(out) :: larr(lbl:ubl) ! output landunit array
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: c,l,index              ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_c2l(lbc:ubc)        ! scale factor for column->landunit mapping
    real(r8) :: sumwt(lbl:ubl)         ! sum of weights
    real(r8), pointer :: wtlunit(:)    ! weight of landunits relative to gridcells
    integer , pointer :: clandunit(:)  ! gridcell of corresponding column
!------------------------------------------------------------------------

    wtlunit   => clm3%g%l%c%wtlunit
    clandunit => clm3%g%l%c%landunit

    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else
       write(6,*)'c2l_1d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    larr(l:) = spval
    sumwt(:) = 0.
    do c = lbc,ubc
       if (carr(c) /= spval .and. wtlunit(c) /= 0.) then
          l = clandunit(c)
          if (sumwt(l) == 0.) larr(l) = 0.
          larr(l) = larr(l) + carr(c) * scale_c2l(c) * wtlunit(c)
          sumwt(l) = sumwt(l) + wtlunit(c)
       end if
    end do
    found = .false.
    do l = lbl,ubl
       if (sumwt(l) > 1.0_r8 + 1.e-6) then
          found = .true.
          index = l
          exit
       else if (sumwt(l) /= 0._r8) then
          larr(l) = larr(l)/sumwt(l)
       end if
    end do
    if (found) then
       write(6,*)'c2l_1d error: sumwt is greater than 1.0 at l= ',index
       call endrun()
    end if

  end subroutine c2l_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2l_2d
!
! !INTERFACE:
  subroutine c2l_2d (lbc, ubc, lbl, ubl, num2d, carr, larr, c2l_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from columns to landunits
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc            ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl            ! beginning and ending landunit indices
    integer , intent(in)  :: num2d               ! size of second dimension
    real(r8), intent(in)  :: carr(lbc:ubc,num2d) ! input column array
    real(r8), intent(out) :: larr(lbl:ubl,num2d) ! output landunit array
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: j,l,c,index            ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_c2l(lbc:ubc)        ! scale factor for column->landunit mapping
    real(r8) :: sumwt(lbl:ubl)         ! sum of weights
    real(r8), pointer :: wtlunit(:)    ! weight of column relative to landunit
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
!------------------------------------------------------------------------

    wtlunit    => clm3%g%l%c%wtlunit
    clandunit  => clm3%g%l%c%landunit

    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else
       write(6,*)'c2l_2d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    larr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0.
       do c = lbc,ubc
          if (carr(c,j) /= spval .and. wtlunit(c) /= 0.) then
             l = clandunit(c)
             if (sumwt(l) == 0.) larr(l,j) = 0.
             larr(l,j) = larr(l,j) + carr(c,j) * scale_c2l(c) * wtlunit(c)
             sumwt(l) = sumwt(l) + wtlunit(c)
          end if
       end do
       found = .false.
       do l = lbl,ubl
          if (sumwt(l) > 1.0_r8 + 1.e-6) then
             found = .true.
             index = l
             exit
          else if (sumwt(l) /= 0._r8) then
             larr(l,j) = larr(l,j)/sumwt(l)
          end if
       end do
       if (found) exit
    end do
    if (found) then
       write(6,*)'c2l_2d error: sumwt is greater than 1.0 at l= ',index,' lev= ',j
       call endrun()
    end if

  end subroutine c2l_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2g_1d
!
! !INTERFACE:
  subroutine c2g_1d(lbc, ubc, lbl, ubl, lbg, ubg, carr, garr, &
       c2l_scale_type, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from columns to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
    integer , intent(in)  :: lbg, ubg              ! beginning and ending landunit indices
    real(r8), intent(in)  :: carr(lbc:ubc)         ! input column array
    real(r8), intent(out) :: garr(lbg:ubg)         ! output gridcell array
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: c,l,g,index            ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of columns relative to gridcells
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: cgridcell(:)  ! gridcell of corresponding column
!------------------------------------------------------------------------

    wtgcell   => clm3%g%l%c%wtgcell
    clandunit => clm3%g%l%c%landunit
    cgridcell => clm3%g%l%c%gridcell

    if (l2g_scale_type == 'unity') then
       do l = lbl,ubl
          scale_l2g(l) = 1.0_r8
       end do
    else
       write(6,*)'c2l_1d error: scale type ',l2g_scale_type,' not supported'
       call endrun()
    end if
    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else
       write(6,*)'c2l_1d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    garr(:) = spval
    sumwt(:) = 0.
    do c = lbc,ubc
       if (carr(c) /= spval .and. wtgcell(c) /= 0.) then
          l = clandunit(c)
          g = cgridcell(c)
          if (sumwt(g) == 0.) garr(g) = 0.
          garr(g) = garr(g) + carr(c) * scale_c2l(c) * scale_l2g(l) * wtgcell(c)
          sumwt(g) = sumwt(g) + wtgcell(c)
       end if
    end do
    found = .false.
    do g = lbg, ubg
       if (sumwt(g) > 1.0_r8 + 1.e-6) then
          found = .true.
          index = g
          exit
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       write(6,*)'c2g_1d error: sumwt is greater than 1.0 at g= ',index
       call endrun()
    end if

  end subroutine c2g_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: c2g_2d
!
! !INTERFACE:
  subroutine c2g_2d(lbc, ubc, lbl, ubl, lbg, ubg, num2d, carr, garr, &
       c2l_scale_type, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from columns to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbc, ubc              ! beginning and ending column indices
    integer , intent(in)  :: lbl, ubl              ! beginning and ending landunit indices
    integer , intent(in)  :: lbg, ubg              ! beginning and ending gridcell indices
    integer , intent(in)  :: num2d                 ! size of second dimension
    real(r8), intent(in)  :: carr(lbc:ubc,num2d)   ! input column array
    real(r8), intent(out) :: garr(lbg:ubg,num2d)   ! output gridcell array
    character(len=*), intent(in) :: c2l_scale_type ! scale factor type for averaging
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: j,c,g,l,index          ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_c2l(lbc:ubc)     ! scale factor
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of columns relative to gridcells
    integer , pointer :: clandunit(:)  ! landunit of corresponding column
    integer , pointer :: cgridcell(:)  ! gridcell of corresponding column
!------------------------------------------------------------------------

    wtgcell   => clm3%g%l%c%wtgcell
    clandunit => clm3%g%l%c%landunit
    cgridcell => clm3%g%l%c%gridcell

    if (l2g_scale_type == 'unity') then
       do l = lbl,ubl
          scale_l2g(l) = 1.0_r8
       end do
    else
       write(6,*)'c2g_2d error: scale type ',l2g_scale_type,' not supported'
       call endrun()
    end if
    if (c2l_scale_type == 'unity') then
       do c = lbc,ubc
          scale_c2l(c) = 1.0_r8
       end do
    else
       write(6,*)'c2g_2d error: scale type ',c2l_scale_type,' not supported'
       call endrun()
    end if

    garr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0.
       do c = lbc,ubc
          if (carr(c,j) /= spval .and. wtgcell(c) /= 0.) then
             l = clandunit(c)
             g = cgridcell(c)
             if (sumwt(g) == 0.) garr(g,j) = 0.
             garr(g,j) = garr(g,j) + carr(c,j) * scale_c2l(c) * scale_l2g(l) * wtgcell(c)
             sumwt(g) = sumwt(g) + wtgcell(c)
          end if
       end do
       found = .false.
       do g = lbg, ubg
          if (sumwt(g) > 1.0_r8 + 1.e-6) then
             found = .true.
             index = g
             exit
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) exit
    end do
    if (found) then
       write(6,*)'c2g_2d error: sumwt is greater than 1.0 at g= ',index
       call endrun()
    end if

  end subroutine c2g_2d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: l2g_1d
!
! !INTERFACE:
  subroutine l2g_1d(lbl, ubl, lbg, ubg, larr, garr, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from landunits to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbl, ubl       ! beginning and ending sub landunit indices
    integer , intent(in)  :: lbg, ubg       ! beginning and ending gridcell indices
    real(r8), intent(in)  :: larr(lbl:ubl)  ! input landunit array
    real(r8), intent(out) :: garr(lbg:ubg)  ! output gridcell array
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: l,g,index              ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of landunits relative to gridcells
    integer , pointer :: lgridcell(:)  ! gridcell of corresponding landunit
!------------------------------------------------------------------------

    wtgcell   => clm3%g%l%wtgcell
    lgridcell => clm3%g%l%gridcell

    if (l2g_scale_type == 'unity') then
       do l = lbl,ubl
          scale_l2g(l) = 1.0_r8
       end do
    else
       write(6,*)'l2g_1d error: scale type ',l2g_scale_type,' not supported'
       call endrun()
    end if

    garr(:) = spval
    sumwt(:) = 0.
    do l = lbl,ubl
       if (larr(l) /= spval .and. wtgcell(l) /= 0.) then
          g = lgridcell(l)
          if (sumwt(g) == 0.) garr(g) = 0.
          garr(g) = garr(g) + larr(l) * scale_l2g(l) * wtgcell(l)
          sumwt(g) = sumwt(g) + wtgcell(l)
       end if
    end do
    found = .false.
    do g = lbg, ubg
       if (sumwt(g) > 1.0_r8 + 1.e-6) then
          found = .true.
          index = g
          exit
       else if (sumwt(g) /= 0._r8) then
          garr(g) = garr(g)/sumwt(g)
       end if
    end do
    if (found) then
       write(6,*)'l2g_1d error: sumwt is greater than 1.0 at g= ',index
       call endrun()
    end if

  end subroutine l2g_1d

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: l2g_2d
!
! !INTERFACE:
  subroutine l2g_2d(lbl, ubl, lbg, ubg, num2d, larr, garr, l2g_scale_type)
!
! !DESCRIPTION:
! Perfrom subgrid-average from landunits to gridcells.
! Averaging is only done for points that are not equal to "spval".
!
! !ARGUMENTS:
    implicit none
    integer , intent(in)  :: lbl, ubl             ! beginning and ending column indices
    integer , intent(in)  :: lbg, ubg             ! beginning and ending gridcell indices
    integer , intent(in)  :: num2d                ! size of second dimension
    real(r8), intent(in)  :: larr(lbl:ubl,num2d)  ! input landunit array
    real(r8), intent(out) :: garr(lbg:ubg,num2d)  ! output gridcell array
    character(len=*), intent(in) :: l2g_scale_type ! scale factor type for averaging
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein 12/03
!
!EOP
!
! LOCAL VARIABLES:
    integer  :: j,g,l,index            ! indices
    logical  :: found                  ! temporary for error check
    real(r8) :: scale_l2g(lbl:ubl)     ! scale factor
    real(r8) :: sumwt(lbg:ubg)         ! sum of weights
    real(r8), pointer :: wtgcell(:)    ! weight of landunits relative to gridcells
    integer , pointer :: lgridcell(:)  ! gridcell of corresponding landunit
!------------------------------------------------------------------------

    wtgcell   => clm3%g%l%wtgcell
    lgridcell => clm3%g%l%gridcell

    if (l2g_scale_type == 'unity') then
       do l = lbl,ubl
          scale_l2g(l) = 1.0_r8
       end do
    else
       write(6,*)'l2g_2d error: scale type ',l2g_scale_type,' not supported'
       call endrun()
    end if

    garr(:,:) = spval
    do j = 1,num2d
       sumwt(:) = 0.
       do l = lbl,ubl
          if (larr(l,j) /= spval .and. wtgcell(l) /= 0.) then
             g = lgridcell(l)
             if (sumwt(g) == 0.) garr(g,j) = 0.
             garr(g,j) = garr(g,j) + larr(l,j) * scale_l2g(l) * wtgcell(l)
             sumwt(g) = sumwt(g) + wtgcell(l)
          end if
       end do
       found = .false.
       do g = lbg,ubg
          if (sumwt(g) > 1.0_r8 + 1.e-6) then
             found = .true.
             index= g
             exit
          else if (sumwt(g) /= 0._r8) then
             garr(g,j) = garr(g,j)/sumwt(g)
          end if
       end do
       if (found) exit
    end do
    if (found) then
       write(6,*)'l2g_2d error: sumwt is greater than 1.0 at g= ',index,' lev= ',j
       call endrun()
    end if

  end subroutine l2g_2d

end module subgridAveMod

