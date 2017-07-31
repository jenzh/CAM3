!
! $Id: comsta.h,v 1.6.42.1 2004/01/27 17:10:50 rosinski Exp $
! $Author: rosinski $
!
!
! Diagnostic statistics integrals
!
      common/comsta/rmsz(plat)  ,rmsd (plat)  ,rmst(plat)  ,stq(plat)
      common/comsta/             psurf(plat)
!
      real(r8) rmsz             ! lambda/p sum of w*dp/ps times square vorticity
      real(r8) rmsd             ! lambda/p sum of w*dp/ps times square divergence
      real(r8) rmst             ! lambda/p sum of w*dp/ps times square temperature
      real(r8) stq              ! lambda/p sum of w*dp/ps times square moisture
      real(r8) psurf            ! lambda/p sum of w*dp/ps times square surface press
