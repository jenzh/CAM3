#include <misc.h>
#include <preproc.h>

module initGridcellsMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: initGridcellsMod
!
! !DESCRIPTION:
! Initializes sub-grid mapping for each land grid cell
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clmtype
  use clm_varpar, only : lsmlon, lsmlat, maxpatch, maxpatch_pft
  use clm_varsur, only : numlon, area, latixy, longxy, landfrac

  use shr_sys_mod, only : shr_sys_flush
  use spmdMod    , only : iam
  use abortutils, only : endrun
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public initGridcells      ! Initialize sub-grid gridcell mapping
!
! !PIVATE MEMBER FUNCTIONS:
  private landunit_veg_compete
  private landunit_veg_noncompete
  private landunit_special
  private landunit_crop_noncompete
  private initGridcellsGlob  ! Initialize global part clmtype
                             ! (topological info)
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL MODULE VARIABLES:
  type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
  type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
  type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
  type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initGridcells
!
! !INTERFACE:
  subroutine initGridcells (vegxy, wtxy)
!
! !DESCRIPTION:
! Initialize sub-grid mapping and allocates space for derived type
! hierarchy.  For each land gridcell determine landunit, column and
! pft properties.  Note that ngcells, nlunits, ncols and npfts are
! per-processor totals here and are currently not used for anything other
! than placeholders.  Determine if there are any vegetated landunits and
! if so---the weight of the vegetated landunit relative to the gridcell
! The first landunit contains all the vegetated patches (if any) For now,
! the vegetated patches will all be gathered on a single landunit, with
! each vegetated type having its own column on that landunit.  The special
! patches (urban, lake, wetland, glacier) each get their own landunit
! having a single column and one non-vegetated pfts
!
! !USES:
    use decompMod    , only : get_proc_bounds, get_proc_global, &
                              get_gcell_info, get_gcell_xyind, &
                              get_sn_index
    use shr_const_mod, only : SHR_CONST_PI
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) ! PFT type
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid patch
                                                          ! weights
!
! !REVISION HISTORY:
! Created by Peter Thornton and Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: g,i,j,m,n,gi,li,ci,pi ! indices
    integer :: ngcells     ! temporary dummy
    integer :: nlunits     ! temporary dummy
    integer :: ncols       ! temporary dummy
    integer :: npfts       ! temporary dummy
    integer :: nveg        ! number of pfts in naturally vegetated landunit
    real(r8):: wtveg       ! weight (relative to gridcell) of naturally vegetated landunit
    integer :: ncrop       ! number of crop pfts in crop landunit
    real(r8):: wtcrop      ! weight (relative to gridcell) of crop landunit
    integer :: begp, endp  ! per-proc beginning and ending pft indices
    integer :: begc, endc  ! per-proc beginning and ending column indices
    integer :: begl, endl  ! per-proc beginning and ending landunit indices
    integer :: begg, endg  ! per-proc gridcell ending gridcell indices
    integer :: numg        ! total number of gridcells across all processors
    integer :: numl        ! total number of landunits across all processors
    integer :: numc        ! total number of columns across all processors
    integer :: nump        ! total number of pfts across all processors
    integer :: ier         ! error status
    integer :: ilunits, icols, ipfts  ! temporaries
!------------------------------------------------------------------------

    ! Set pointers into derived types for this module

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)
    call get_gcell_xyind(begg, endg, gptr%ixy, gptr%jxy)

    ! Determine number of land gridcells on this processor

    clm3%ngridcells = endg - begg + 1

    ! Determine gridcell properties.
    ! Set area, weight, and type information for this gridcell.
    ! For now there is only one type of gridcell, value = 1
    ! Still need to resolve the calculation of area for the gridcell

    ngcells = begg-1
    nlunits = begl-1
    ncols   = begc-1
    npfts   = begp-1

    do gi = begg, endg

       ! Get 2d grid indices

       i = gptr%ixy(gi)
       j = gptr%jxy(gi)

       gptr%area(gi)   = area(i,j)
       gptr%itype(gi)  = 1
      !gptr%wtglob(g)  = gptr%area(g)/clm3%area
       gptr%lat(gi)    = latixy(i,j) * SHR_CONST_PI/180.
       gptr%lon(gi)    = longxy(i,j) * SHR_CONST_PI/180.
       gptr%latdeg(gi) = latixy(i,j)
       gptr%londeg(gi) = longxy(i,j)
       gptr%landfrac(gi) = landfrac(i,j)

       gptr%luni(gi) = nlunits + 1
       gptr%coli(gi) = ncols   + 1
       gptr%pfti(gi) = npfts   + 1

       call get_gcell_info(i, j, wtxy, nlunits=ilunits, ncols=icols, npfts=ipfts)

       ngcells = ngcells + 1
       nlunits = nlunits + ilunits
       ncols   = ncols   + icols
       npfts   = npfts   + ipfts

       gptr%lunf(gi) = nlunits
       gptr%colf(gi) = ncols
       gptr%pftf(gi) = npfts

       gptr%nlandunits(gi) = gptr%lunf(gi) - gptr%luni(gi) + 1
       gptr%ncolumns(gi)   = gptr%colf(gi) - gptr%coli(gi) + 1
       gptr%npfts(gi)      = gptr%pftf(gi) - gptr%pfti(gi) + 1

    end do

    ! For each land gridcell determine landunit, column and pft properties.

    ngcells = 0
    nlunits = 0
    ncols   = 0
    npfts   = 0

    li = begl - 1
    ci = begc - 1
    pi = begp - 1

    do gi = begg,endg

       ! Determine 2d lat and lon indices

       i = gptr%ixy(gi)
       j = gptr%jxy(gi)

       ! Obtain gridcell properties

       call get_gcell_info(i, j, wtxy, nveg=nveg, wtveg=wtveg, ncrop=ncrop, wtcrop=wtcrop)

       ! Determine naturally vegetated landunit

#if (defined NOCOMPETE)
       if (nveg > 0) call landunit_veg_noncompete(nveg, wtveg, wtxy, vegxy, i, j, gi, li, ci, pi)
#else
       if (nveg > 0) call landunit_veg_compete(nveg, wtveg, wtxy, vegxy, i, j, gi, li, ci, pi)
#endif

       ! Determine crop landunit.

       if (ncrop > 0) call landunit_crop_noncompete(ncrop, wtcrop, wtxy, vegxy, i, j, gi, li, ci, pi)

       ! Determine special landunits (urban, lake, wetland, glacier).

       do m = npatch_urban, npatch_glacier
          if (wtxy(i,j,m) > 0.) call landunit_special(wtxy, i, j, m, gi, li, ci, pi)
       end do

    end do

    ! Determine global values

    call initGridcellsGlob()

    ! Determine south->north indices (the following are global quantities)

    call get_sn_index(gptr%snindex, type1d=nameg)
    call get_sn_index(lptr%snindex, type1d=namel)
    call get_sn_index(cptr%snindex, type1d=namec)
    call get_sn_index(pptr%snindex, type1d=namep)

  end subroutine initGridcells

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: landunit_veg_compete
!
! !INTERFACE:
  subroutine landunit_veg_compete (nveg, wtveg, wtxy, vegxy, i, j, &
                                   gi, li, ci, pi)
!
! !DESCRIPTION:
! Initialize vegetated landunit with competition
!
! !USES:
    use clm_varcon, only : istsoil
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nveg   ! number of vegetated patches in gridcell
    real(r8), intent(in) :: wtveg  ! weight relative to gridcell of veg
                                   ! landunit
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid patch 
                                                          ! weights
    integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) ! PFT type
    integer , intent(in) :: i      ! 2d longitude index
    integer , intent(in) :: j      ! 2d latitude index
    integer , intent(in) :: gi     ! gridcell index
    integer , intent(inout) :: li  ! landunit index
    integer , intent(inout) :: ci  ! column index
    integer , intent(inout) :: pi  ! pft index
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                          ! indices
!------------------------------------------------------------------------

    ! Set landunit properties
    ! Increment landunits and set indices into lower levels in hierarchy and higher levels
    ! in hierarchy and topological mapping functionality

    li = li + 1
    lptr%ncolumns(li) = 1
    lptr%coli(li) = ci + 1
    lptr%colf(li) = ci + 1
    lptr%npfts(li) = nveg
    lptr%pfti(li) = pi + 1
    lptr%pftf(li) = pi + nveg
    lptr%area(li) = gptr%area(gi) * wtveg
    lptr%gridcell(li) = gi
    lptr%wtgcell(li) = wtveg
    lptr%ixy(li) = i
    lptr%jxy(li) = j
    lptr%latdeg(li) = latixy(i,j)
    lptr%londeg(li) = longxy(i,j)
    lptr%ifspecial(li) = .false.
    lptr%lakpoi(li) = .false.
    lptr%itype(li) = istsoil

    ! Set column properties for this landunit
    ! Increment column  - set only one column on compete landunit -  and set indices into
    ! lower levels in hierarchy, higher levels in hierarchy and topological mapping
    ! functionality (currently all columns have type 1)

    ci = ci + 1
    cptr%npfts(ci) = nveg
    cptr%pfti(ci) = pi + 1
    cptr%pftf(ci) = pi + nveg
    cptr%area(ci) = lptr%area(li)
    cptr%landunit(ci) = li
    cptr%gridcell(ci) = gi
    cptr%wtlunit(ci) = 1.0
    cptr%wtgcell(ci) = wtveg
    cptr%ixy(ci) = i
    cptr%jxy(ci) = j
    cptr%latdeg(ci) = latixy(i,j)
    cptr%londeg(ci) = longxy(i,j)
    cptr%itype(ci) = 1

    ! Set pft properties for this landunit
    ! Topological mapping functionality

!dir$ concurrent
!cdir nodep
    do m = 1,maxpatch_pft
       if (wtxy(i,j,m) > 0.) then
          pi = pi+1
          pptr%column(pi) = ci
          pptr%landunit(pi) = li
          pptr%gridcell(pi) = gi
          pptr%wtcol(pi) = wtxy(i,j,m) / wtveg
          pptr%wtlunit(pi) = wtxy(i,j,m) / wtveg
          pptr%wtgcell(pi) = wtxy(i,j,m)
          pptr%area(pi) = cptr%area(ci) * pptr%wtcol(pi)
          pptr%ixy(pi) = i
          pptr%jxy(pi) = j
          pptr%mxy(pi) = m
          pptr%latdeg(pi) = latixy(i,j)
          pptr%londeg(pi) = longxy(i,j)
          pptr%itype(pi) = vegxy(i,j,m)
       end if ! non-zero weight for this pft
    end do ! loop through maxpatch_pft

  end subroutine landunit_veg_compete

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: landunit_veg_noncompete
!
! !INTERFACE:
  subroutine landunit_veg_noncompete (nveg, wtveg, wtxy, vegxy, i, j, &
                                      gi, li, ci, pi)
!
! !DESCRIPTION:
! Initialize vegetated landunit without competition
!
! !USES:
    use clm_varcon, only : istsoil
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: nveg       ! number of vegetated patches in gridcell
    real(r8), intent(in) :: wtveg      ! weight relative to gridcell of veg landunit
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid patch weights
    integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) ! PFT type
    integer , intent(in) :: i          ! 2d longitude index
    integer , intent(in) :: j          ! 2d latitude index
    integer , intent(in) :: gi         ! gridcell index
    integer , intent(inout) :: li      ! landunit index
    integer , intent(inout) :: ci      ! column index
    integer , intent(inout) :: pi      ! pft index
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                          ! indices
    real(r8) :: wtlunit                    ! weight relative to landunit
!------------------------------------------------------------------------

    ! Set landunit properties
    ! Increment landunits and set indices into lower levels in hierarchy and higher levels
    ! in hierarchy and topological mapping functionality

    li = li + 1
    lptr%ncolumns(li) = nveg
    lptr%coli(li) = ci + 1
    lptr%colf(li) = ci + nveg
    lptr%npfts(li) = nveg
    lptr%pfti(li) = pi + 1
    lptr%pftf(li) = pi + nveg
    lptr%area(li) = gptr%area(gi) * wtveg
    lptr%gridcell(li) = gi
    lptr%wtgcell(li) = wtveg
    lptr%ixy(li) = i
    lptr%jxy(li) = j
    lptr%latdeg(li) = latixy(i,j)
    lptr%londeg(li) = longxy(i,j)
    lptr%ifspecial(li) = .false.
    lptr%lakpoi(li) = .false.
    lptr%itype(li) = istsoil

    ! Set column properties for this landunit
    ! Increment column  - each column has its own pft -  and set indices into
    ! lower levels in hierarchy, higher levels in hierarchy and topological mapping
    ! functionality (currently all columns have type 1)
    ! Set column and pft properties
    ! Loop through regular (vegetated) patches, assign one column for each
    ! vegetated patch with non-zero weight. The weights for each column on
    ! the vegetated landunit must add to one when summed over the landunit,
    ! so the wtxy(i,j,m) values are taken relative to the total wtveg.

!dir$ concurrent
!cdir nodep
    do m = 1, maxpatch_pft
       if (wtxy(i,j,m) > 0.) then

          ! Determine weight relative to landunit of pft/column

          wtlunit = wtxy(i,j,m) / wtveg

          ! Increment number of columns on landunit

          ci = ci + 1
          cptr%npfts(ci) = 1
          cptr%pfti(ci) = ci
          cptr%pftf(ci) = ci
          cptr%area(ci) = lptr%area(li) * wtlunit
          cptr%landunit(ci) = li
          cptr%gridcell(ci) = gi
          cptr%wtlunit(ci) = wtlunit
          cptr%wtgcell(ci) = cptr%area(ci) / gptr%area(gi)
          cptr%ixy(ci) = i
          cptr%jxy(ci) = j
          cptr%latdeg(ci) = latixy(i,j)
          cptr%londeg(ci) = longxy(i,j)
          cptr%itype(ci) = 1

          ! Increment number of pfts on this landunit
          ! Set area, weight (relative to column) and type information for this pft
          ! For now, a single pft per column, so weight = 1
          ! pft type comes from the m dimension of wtxy()
          ! Set grid index, weight (relative to grid cell)
          ! and m index (needed for laixy, etc. reference)

          pi = pi + 1
          pptr%column(pi) = ci
          pptr%landunit(pi) = li
          pptr%gridcell(pi) = gi
          pptr%wtcol(pi) = 1.0
          pptr%wtlunit(pi) = cptr%wtlunit(ci)
          pptr%area(pi) = cptr%area(ci)
          pptr%wtgcell(pi) = pptr%area(pi) / gptr%area(gi)
          pptr%ixy(pi) = i
          pptr%jxy(pi) = j
          pptr%mxy(pi) = m
          pptr%latdeg(pi) = latixy(i,j)
          pptr%londeg(pi) = longxy(i,j)
          pptr%itype(pi) = vegxy(i,j,m)

       end if   ! end if non-zero weight
    end do   ! end loop through the possible vegetated patch indices

  end subroutine landunit_veg_noncompete

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: landunit_special
!
! !INTERFACE:
  subroutine landunit_special (wtxy, i, j, m, gi, li, ci, pi)
!
! !DESCRIPTION:
! Initialize special landunits (urban, lake, wetland, glacier)
!
! !USES:
    use pftvarcon, only : noveg
    use clm_varcon, only : istice, istwet, istdlak, isturb
    use clm_varpar, only : npatch_lake, npatch_wet, npatch_urban, &
                           npatch_glacier
!
! !ARGUMENTS:
    implicit none
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  !subgrid patch weights
    integer, intent(in) :: i            !2-dim longitude index
    integer, intent(in) :: j            !2-dim latitude index
    integer, intent(in) :: m            !2-dim PFT patch index
    integer, intent(in) :: gi           !gridcell index
    integer, intent(inout) :: li        !landunit index
    integer, intent(inout) :: ci        !column index
    integer, intent(inout) :: pi        !pft index
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: c             !column loop index
    integer  :: ncols         !number of columns
    integer  :: npfts         !number of pfts
    integer  :: ier           !error status
    real(r8) :: weight        !temporary weight
    integer  :: itype         !landunit type
!------------------------------------------------------------------------

    ! Define landunit type

    if (m == npatch_lake) then         !deep lake (from pctlak)
       itype = istdlak
    else if (m == npatch_wet) then     !wetland (from pctwet)
       itype = istwet
    else if (m == npatch_glacier) then !glacier (from pctgla)
       itype = istice
    else if (m == npatch_urban) then   !urban (from pcturb)
       itype = isturb
    else                               !error
       write(6,*)'special landunit are currently only:', &
            ' deep lake, wetland, glacier or urban)'
       call endrun()
    endif

    ! Determine landunit index and landunit properties

    li = li + 1
    lptr%ncolumns(li) = 1
    lptr%coli(li) = ci + 1
    lptr%colf(li) = ci + 1
    lptr%npfts(li) = 1
    lptr%pfti(li) = pi + 1
    lptr%pftf(li) = pi + 1
    lptr%area(li) = gptr%area(gi) * wtxy(i,j,m)
    lptr%gridcell(li) = gi
    lptr%wtgcell(li) = lptr%area(li) / gptr%area(gi)
    lptr%ixy(li) = i
    lptr%jxy(li) = j
    lptr%latdeg(li) = latixy(i,j)
    lptr%londeg(li) = longxy(i,j)
    lptr%ifspecial(li) = .true.
    if (itype == istdlak) then
       lptr%lakpoi(li) = .true.
    else
       lptr%lakpoi(li) = .false.
    end if
    lptr%itype(li) = itype

    ! For the special landunits there currently is only one column
    ! Later, the age classes will be implemented on different columns within
    ! the same landunit, so the column type will correspond to an age class

    ncols = 1

    ! Loop through columns for this landunit and set the column properties
    ! We know that there is only one column for the special landunit - but
    ! the loop is included for future consistency.

    do c = 1,ncols

       ! Determine column index and column properties
       ! For now all columns have the same type, value = 1

       weight = 1.0/ncols

       ci = ci + c
       cptr%npfts(ci) = 1
       cptr%pfti(ci) = pi + 1
       cptr%pftf(ci) = pi + 1
       cptr%area(ci) = lptr%area(li) * weight
       cptr%landunit(ci) = li
       cptr%gridcell(ci) = gi
       cptr%wtlunit(ci) = weight
       cptr%wtgcell(ci) = cptr%area(ci) / gptr%area(gi)
       cptr%ixy(ci) = i
       cptr%jxy(ci) = j
       cptr%latdeg(ci) = latixy(i,j)
       cptr%londeg(ci) = longxy(i,j)
       cptr%itype(ci) = 1

       ! Determine pft index and pft properties
       ! Each column has one non-vegetated pft
       ! Set area, weight (relative to column), and type information
       ! for this non-vegetated pft
       ! Set grid index, weight (relative to grid cell) and
       ! m index (needed for laixy, etc. reference)

       npfts = 1
       weight = 1.0/npfts

       pi = pi + 1
       pptr%column(pi) = ci
       pptr%landunit(pi) = li
       pptr%gridcell(pi) = gi
       pptr%area(pi) = lptr%area(li) * weight
       pptr%wtcol(pi) = weight
       pptr%wtlunit(pi) = cptr%wtlunit(ci)
       pptr%wtgcell(pi) = pptr%area(pi) / gptr%area(gi)
       pptr%ixy(pi) = i
       pptr%jxy(pi) = j
       pptr%mxy(pi) = m
       pptr%latdeg(pi) = latixy(i,j)
       pptr%londeg(pi) = longxy(i,j)
       pptr%itype(pi) = noveg

    end do   ! end loop through ncolumns

  end subroutine landunit_special

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: landunit_crop_noncompete
!
! !INTERFACE:
  subroutine landunit_crop_noncompete (ncrop, wtcrop, wtxy, vegxy, i, j, &
                                       gi, li, ci, pi)
!
! !DESCRIPTION:
! Initialize crop landunit without competition
!
! !USES:
    use clm_varcon, only : istsoil
    use clm_varpar, only : npatch_crop
!
! !ARGUMENTS:
    implicit none
    integer , intent(in) :: ncrop       ! number of vegetated patches in gridcell
    real(r8), intent(in) :: wtcrop      ! weight relative to gridcell of veg landunit
    real(r8), intent(in) :: wtxy(lsmlon,lsmlat,maxpatch)  ! subgrid patch weights
    integer , intent(in) :: vegxy(lsmlon,lsmlat,maxpatch) ! PFT type
    integer , intent(in) :: i          ! 2d longitude index
    integer , intent(in) :: j          ! 2d latitude index
    integer , intent(in) :: gi         ! gridcell index
    integer , intent(inout) :: li      ! landunit index
    integer , intent(inout) :: ci      ! column index
    integer , intent(inout) :: pi      ! pft index
!
! !REVISION HISTORY:
! Created by Sam Levis
!
!EOP
!
! !LOCAL VARIABLES:
    integer  :: m                          ! indices
    real(r8) :: wtlunit                    ! weight relative to landunit
!------------------------------------------------------------------------

    ! Set landunit properties
    ! Increment landunits and set indices into lower levels in hierarchy and higher levels
    ! in hierarchy and topological mapping functionality

    li = li + 1
    lptr%ncolumns(li) = ncrop
    lptr%coli(li) = ci + 1
    lptr%colf(li) = ci + ncrop
    lptr%npfts(li) = ncrop
    lptr%pfti(li) = pi + 1
    lptr%pftf(li) = pi + ncrop
    lptr%area(li) = gptr%area(gi) * wtcrop
    lptr%gridcell(li) = gi
    lptr%wtgcell(li) = wtcrop
    lptr%ixy(li) = i
    lptr%jxy(li) = j
    lptr%latdeg(li) = latixy(i,j)
    lptr%londeg(li) = longxy(i,j)
    lptr%ifspecial(li) = .false.
    lptr%lakpoi(li) = .false.
    lptr%itype(li) = istsoil

    ! Set column properties for this landunit
    ! Increment column  - each column has its own pft -  and set indices into
    ! lower levels in hierarchy, higher levels in hierarchy and topological mapping
    ! functionality (currently all columns have type 1)
    ! Set column and pft properties
    ! Loop through regular (vegetated) patches, assign one column for each
    ! vegetated patch with non-zero weight. The weights for each column on
    ! the vegetated landunit must add to one when summed over the landunit,
    ! so the wtxy(i,j,m) values are taken relative to the total wtcrop.

!dir$ concurrent
!cdir nodep
    do m = npatch_glacier+1, npatch_crop
       if (wtxy(i,j,m) > 0.) then

          ! Determine weight of crop pft/column relative to crop landunit

          wtlunit = wtxy(i,j,m) / wtcrop

          ! Increment number of columns on landunit

          ci = ci + 1
          cptr%npfts(ci) = 1
          cptr%area(ci) = lptr%area(li) * wtlunit
          cptr%landunit(ci) = li
          cptr%gridcell(ci) = gi
          cptr%wtlunit(ci) = wtlunit
          cptr%wtgcell(ci) = cptr%area(ci) / gptr%area(gi)
          cptr%ixy(ci) = i
          cptr%jxy(ci) = j
          cptr%latdeg(ci) = latixy(i,j)
          cptr%londeg(ci) = longxy(i,j)
          cptr%itype(ci) = 1

          ! Increment number of pfts on this landunit
          ! Set area, weight (relative to column) and type information for this pft
          ! For now, a single pft per column, so weight relative to column is 1
          ! pft type comes from the m dimension of wtxy()
          ! Set grid index, weight relative to grid cell and m index (needed for laixy, etc.)

          pi = pi + 1
          pptr%column(pi) = ci
          pptr%landunit(pi) = li
          pptr%gridcell(pi) = gi
          pptr%wtcol(pi) = 1.0
          pptr%wtlunit(pi) = cptr%wtlunit(ci)
          pptr%area(pi) = cptr%area(ci)
          pptr%wtgcell(pi) = pptr%area(pi) / gptr%area(gi)
          pptr%ixy(pi) = i
          pptr%jxy(pi) = j
          pptr%mxy(pi) = m
          pptr%latdeg(pi) = latixy(i,j)
          pptr%londeg(pi) = longxy(i,j)
          pptr%itype(pi) = vegxy(i,j,m)

          ! Set pft indices for column

          cptr%pfti(ci) = pi
          cptr%pftf(ci) = pi

       end if   ! end if non-zero weight
    end do   ! end loop through the possible vegetated patch indices

  end subroutine landunit_crop_noncompete

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: initGridcellsGlob
!
! !INTERFACE:
  subroutine initGridcellsGlob()
#if (defined SPMD)
!
! !DESCRIPTION:
! Set up 1d array of weights and indices for xy mapping.
! Note: if DGVM is defined, weights are updated in DGVM mode
!
! !USES:
    use decompMod      , only : get_proc_bounds, get_proc_global
    use spmdMod        , only : masterproc, MPI_INTEGER, mpicom
    use spmdGathScatMod, only : gather_data_to_master
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    integer :: numg         ! total number of gridcells across all processors
    integer :: numl         ! total number of landunits across all processors
    integer :: numc         ! total number of columns across all processors
    integer :: nump         ! total number of pfts across all processors
    integer :: ier          ! error status
    real(r8), pointer :: rloc(:)       !temporaries for mpi gather
    integer , pointer :: iloc(:)       !temporaries for mpi gather
    real(r8), pointer :: rglob(:)      !temporaries for mpi gather
    integer , pointer :: iglob(:)      !temporaries for mpi gather
!-------------------------------------------------------------------------

    ! Determine necessary indices

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)
    call get_proc_global(numg, numl, numc, nump)

    ! gridcell gather

    allocate(rloc(numg), rglob(numg), iloc(numg), iglob(numg), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initGridcellsGlob(): rloc,rglob,iloc,iglob (numg) allocation error'
       call endrun
    end if

    rloc(:) = gptr%wtglob(:)  ! this just has begg:endg filled in here
    call gather_data_to_master(rloc, rglob, clmlevel=nameg)
    if (masterproc) gptr%wtglob(:) = rglob(:)

    iloc(:) = gptr%ixy(:)
    call gather_data_to_master(iloc, iglob, clmlevel=nameg)
    call mpi_bcast(iglob, size(iglob), MPI_INTEGER, 0, mpicom, ier)
    gptr%ixy(:) = iglob(:)

    iloc(:) = gptr%jxy(:)
    call gather_data_to_master(iloc, iglob, clmlevel=nameg)
    call mpi_bcast(iglob, size(iglob), MPI_INTEGER, 0, mpicom, ier)
    gptr%jxy(:) = iglob(:)

    rloc(:) = gptr%lat(:)
    call gather_data_to_master(rloc, rglob, clmlevel=nameg)
    if (masterproc) gptr%lat(:) = rglob(:)

    rloc(:) = gptr%lon(:)
    call gather_data_to_master(rloc, rglob, clmlevel=nameg)
    if (masterproc) gptr%lon(:) = rglob(:)

    rloc(:) = gptr%latdeg(:)
    call gather_data_to_master(rloc, rglob, clmlevel=nameg)
    if (masterproc) gptr%latdeg(:) = rglob(:)

    rloc(:) = gptr%londeg(:)
    call gather_data_to_master(rloc, rglob, clmlevel=nameg)
    if (masterproc) gptr%londeg(:) = rglob(:)

    deallocate(rloc, rglob, iloc, iglob)

    ! landunit gather

    allocate(rloc(numl), rglob(numl), iloc(numl), iglob(numl), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initGridcellsGlob(): rloc,rglob,iloc,iglob (numl) allocation error'
       call endrun
    end if

    rloc(:) = lptr%wtgcell(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namel)
    if (masterproc) lptr%wtgcell(:) = rglob(:)

    iloc(:) = lptr%ixy(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namel)
    if (masterproc) lptr%ixy(:) = iglob(:)

    iloc(:) = lptr%jxy(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namel)
    if (masterproc) lptr%jxy(:) = iglob(:)

    iloc(:) = lptr%gridcell(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namel)
    if (masterproc) lptr%gridcell(:) = iglob(:)

    rloc(:) = lptr%latdeg(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namel)
    if (masterproc) lptr%latdeg(:) = rglob(:)

    rloc(:) = lptr%londeg(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namel)
    if (masterproc) lptr%londeg(:) = rglob(:)

    iloc(:) = lptr%itype(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namel)
    if (masterproc) lptr%itype(:) = iglob(:)

    deallocate(rloc, rglob, iloc, iglob)

    ! column gather

    allocate(rloc(numc), rglob(numc), iloc(numc), iglob(numc), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initGridcellsGlob(): rloc,rglob,iloc,iglob (numc) allocation error'
       call endrun
    end if

    rloc(:) = cptr%wtgcell(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namec)
    if (masterproc) cptr%wtgcell(:) = rglob(:)

    rloc(:) = cptr%wtlunit(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namec)
    if (masterproc) cptr%wtlunit(:) = rglob(:)

    iloc(:) = cptr%ixy(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namec)
    if (masterproc) cptr%ixy(:) = iglob(:)

    iloc(:) = cptr%jxy(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namec)
    if (masterproc) cptr%jxy(:) = iglob(:)

    iloc(:) = cptr%gridcell(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namec)
    if (masterproc) cptr%gridcell(:) = iglob(:)

    iloc(:) = cptr%landunit(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namec)
    if (masterproc) cptr%landunit(:) = iglob(:)

    rloc(:) = cptr%latdeg(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namec)
    if (masterproc) cptr%latdeg(:) = rglob(:)

    rloc(:) = cptr%londeg(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namec)
    if (masterproc) cptr%londeg(:) = rglob(:)

    iloc(:) = cptr%itype(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namec)
    if (masterproc) cptr%itype(:) = iglob(:)

    deallocate(rloc, rglob, iloc, iglob)

    ! pft gather

    allocate(rloc(nump), rglob(nump), iloc(nump), iglob(nump), stat=ier)
    if (ier /= 0) then
       write (6,*) 'initGridcellsGlob(): rloc,rglob,iloc,iglob (nump) allocation error'
       call endrun
    end if

    rloc(:) = pptr%wtgcell(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namep)
    if (masterproc) pptr%wtgcell(:) = rglob(:)

    rloc(:) = pptr%wtlunit(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namep)
    if (masterproc) pptr%wtlunit(:) = rglob(:)

    rloc(:) = pptr%wtcol(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namep)
    if (masterproc) pptr%wtcol(:) = rglob(:)

    iloc(:) = pptr%ixy(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namep)
    if (masterproc) pptr%ixy(:) = iglob(:)

    iloc(:) = pptr%jxy(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namep)
    if (masterproc) pptr%jxy(:) = iglob(:)

    iloc(:) = pptr%mxy(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namep)
    if (masterproc) pptr%mxy(:) = iglob(:)

    iloc(:) = pptr%gridcell(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namep)
    if (masterproc) pptr%gridcell(:) = iglob(:)

    iloc(:) = pptr%landunit(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namep)
    if (masterproc) pptr%landunit(:) = iglob(:)

    iloc(:) = pptr%column(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namep)
    if (masterproc) pptr%column(:) = iglob(:)

    rloc(:) = pptr%latdeg(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namep)
    if (masterproc) pptr%latdeg(:) = rglob(:)

    rloc(:) = pptr%londeg(:)
    call gather_data_to_master(rloc, rglob, clmlevel=namep)
    if (masterproc) pptr%londeg(:) = rglob(:)

    iloc(:) = pptr%itype(:)
    call gather_data_to_master(iloc, iglob, clmlevel=namep)
    if (masterproc) pptr%itype(:) = iglob(:)

    deallocate(rloc, rglob, iloc, iglob)

#endif

  end subroutine initGridcellsGlob

end module initGridcellsMod
