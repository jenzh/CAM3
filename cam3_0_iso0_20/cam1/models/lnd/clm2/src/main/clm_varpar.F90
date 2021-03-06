#include <misc.h>
#include <preproc.h>

module clm_varpar

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varpar
!
! !DESCRIPTION:
! Module containing CLM parameters
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Define land surface 2-d grid. This sets the model resolution according
! to cpp directives LSMLON and LSMLAT in preproc.h.
!
  integer, parameter :: lsmlon = LSMLON  !maximum number of longitude points on lsm grid
  integer, parameter :: lsmlat = LSMLAT  !number of latitude points on lsm grid
!
! Define indices used in surface file read
! maxpatch_pft  = max number of vegetated pfts in naturally vegetated landunit
! maxpatch_crop = max number of crop pfts in crop landunit
!
#if (defined DGVM)
  integer, parameter :: maxpatch_pft   = 10
#else
  integer, parameter :: maxpatch_pft   = 4
#endif
  integer, parameter :: maxpatch_cft   = 2
  integer, parameter :: npatch_urban   = maxpatch_pft + 1
  integer, parameter :: npatch_lake    = npatch_urban + 1
  integer, parameter :: npatch_wet     = npatch_lake  + 1
  integer, parameter :: npatch_glacier = npatch_wet   + 1
  integer, parameter :: npatch_crop    = npatch_glacier + maxpatch_cft
  integer, parameter :: maxpatch       = npatch_crop
!
! Determine maximums in subgrid hierarchy
!
  integer, parameter :: max_pft_per_col     = maxpatch_pft
#if (defined NOCOMPETE)
  integer, parameter :: max_col_per_lunit   = maxpatch_pft
#else
  integer, parameter :: max_col_per_lunit   = 1
#endif
  integer, parameter :: max_lunit_per_gcell = 5            !(soil,urban,lake,wetland,glacier)
!
! Define number of levels
!
  integer, parameter :: nlevsoi     =  10   !number of soil layers
  integer, parameter :: nlevlak     =  10   !number of lake layers
  integer, parameter :: nlevsno     =   5   !maximum number of snow layers
!
! Define miscellaneous parameters
!
  integer, parameter :: numwat      =   5   !number of water types (soil, ice, 2 lakes, wetland)
  integer, parameter :: numpft      =  16   !number of plant types
  integer, parameter :: npftpar     =  32   !number of pft parameters (in LPJ) (DGVM only)
  integer, parameter :: numcol      =   8   !number of soil color types
  integer, parameter :: numrad      =   2   !number of solar radiation bands: vis, nir
  integer, parameter :: ndst        =   4   !number of dust size classes (BGC only)
  integer, parameter :: dst_src_nbr =   3   !number of size distns in src soil (BGC only)
  integer, parameter :: sz_nbr      = 200   !number of sub-grid bins in large bin of dust size distribution (BGC only)
  integer, parameter :: nvoc        =   5   !number of voc categories (BGC only)
!
! Define parameters for RTM river routing model
!
  integer, parameter :: rtmlon = 720  !number of rtm longitudes
  integer, parameter :: rtmlat = 360  !number of rtm latitudes
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

end module clm_varpar
