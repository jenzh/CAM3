module spmd_utils

!----------------------------------------------------------------------- 
! 
! Purpose: This module is responsible for miscellaneous SPMD utilities
!          and information that are shared between dynamics and 
!          physics packages.  
! 
! Author:
!   Original routines:  CMS
!   Module:             T. Henderson, December 2003
!
! $Id: spmd_utils.F90,v 1.1.2.2 2004/03/03 19:53:39 pworley Exp $
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make the default access private
   save


!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public pair      ! $$$here...  originally from eul|sld/spmd_dyn
   public ceil2     ! $$$here...  originally from eul|sld/spmd_dyn


!-----------------------------------------------------------------------
! Public data ----------------------------------------------------------
!-----------------------------------------------------------------------
! physics-motivated dynamics decomposition request
   logical, parameter :: def_mirror = .false.                 ! default
   logical, public :: phys_mirror_decomp_req = def_mirror 
                    ! flag indicating whether latitudes and their
                    ! reflections across the equator should assigned 
                    ! to consecutive processes

!-----------------------------------------------------------------------
! Subroutines and functions --------------------------------------------
!-----------------------------------------------------------------------
contains

!========================================================================

   integer function pair(np,p,k)

      integer np,p,k,q
      q = ieor(p,k)
      if(q.gt.np-1) then
         pair = -1
      else
         pair = q
      endif
      return

   end function pair

!========================================================================

  integer function ceil2(n)
     integer n,p
     p=1
     do while(p.lt.n)
        p=p*2
     enddo
     ceil2=p
     return
  end function ceil2

!========================================================================

end module spmd_utils

