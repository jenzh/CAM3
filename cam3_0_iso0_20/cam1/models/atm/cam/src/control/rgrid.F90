#include <misc.h>

module rgrid

  use pmgrid, only: plat, platd
  use pspect, only: pmmax, pmax
  use infnan, only: bigint

  implicit none

  integer :: nlon(plat)        = bigint ! num longitudes per latitude
  integer :: nlonex(platd)     = bigint ! num longitudes per lat (extended grid)
  integer :: beglatpair(pmmax) = bigint
#if ( defined SCAM )
  integer :: nmmax(plat)     = bigint
#else
  integer :: nmmax(plat/2)     = bigint
#endif
  integer :: wnummax(plat)     = bigint ! cutoff Fourier wavenumber

end module rgrid
