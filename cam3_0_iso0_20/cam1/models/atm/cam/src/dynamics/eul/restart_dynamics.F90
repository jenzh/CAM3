#include <misc.h>
#include <params.h>

module restart_dynamics

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use prognostics
   use ppgrid, only: pcols, pver
   use comslt
   use binary_io
#if ( defined BFB_CAM_SCAM_IOP )
   use iop
#endif
   use abortutils, only: endrun

   implicit none

CONTAINS

   subroutine write_restart_dynamics (nrg)

#include <comqfl.h>

!
! Input arguments
!
      integer, intent(in) :: nrg     ! Unit number
!
! Local workspace
!
      integer :: begj    ! starting latitude
      integer :: ioerr   ! error status
!
      call wrtout_r8 (nrg,vort(1,1,beglat,n3m1), plndlv)
      call wrtout_r8 (nrg,vort(1,1,beglat,n3m2), plndlv)

      call wrtout_r8 (nrg,div(1,1,beglat,n3m1) , plndlv)
      call wrtout_r8 (nrg,div(1,1,beglat,n3m2) , plndlv)

      call wrtout_r8 (nrg,dpsl  ,plond )
      call wrtout_r8 (nrg,dpsm  ,plond )
      call wrtout_r8 (nrg,dps   ,plond )
      call wrtout_r8 (nrg,phis  ,plond )
      call wrtout_r8 (nrg,omga  ,plndlv)
!
! Write fields u3,v3,t3,q3,ps at time indices n3 and n3m1
!
      begj = beglatex + numbnd

      call wrtout_r8 (nrg,u3(1,1,begj,n3m1)  ,plndlv)
      call wrtout_r8 (nrg,v3(1,1,begj,n3m1)  ,plndlv)
      call wrtout_r8 (nrg,t3(1,1,begj,n3m1)  ,plndlv)
      call wrtout_r8 (nrg,ps(1,beglat,n3m1)  ,plond)

      call wrtout_r8 (nrg,u3(1,1,begj,n3m2)  ,plndlv)
      call wrtout_r8 (nrg,v3(1,1,begj,n3m2)  ,plndlv)
      call wrtout_r8 (nrg,t3(1,1,begj,n3m2)  ,plndlv)
      call wrtout_r8 (nrg,ps(1,beglat,n3m2)  ,plond)
      
      call wrtout_r8 (nrg,q3(1,1,1,begj,n3m1),plndlv*(pcnst+pnats))
      call wrtout_r8 (nrg,q3(1,1,1,begj,n3m2),plndlv*(pcnst+pnats))
!
! Write slt arrays (trajectory mid-point coordinates and 
! slt forcast of moisture and constituents
!
      call wrtout_r8 (nrg,lammp,plnlv)
      call wrtout_r8 (nrg,phimp,plnlv)
      call wrtout_r8 (nrg,sigmp,plnlv)
      call wrtout_r8 (nrg,qfcst,plndlv*pcnst)
!
! Write global integrals
!
      if (masterproc) then
         write(nrg, iostat=ioerr) tmass0, fixmas, hw1,    hw2,  &
                                  hw3, alpha
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun ('WRITE_RESTART_DYNAMICS')
         end if
      end if

#if ( defined BFB_CAM_SCAM_IOP )
!
! Write scam values
!
     call wrtout_r8 (nrg,alphasav(1,beglat),pcnst)
     call wrtout_r8 (nrg,dqfx3savm1(1,1,1,beglat),plndlv*pcnst)       
     call wrtout_r8 (nrg,divq3dsav(1,1,1,beglat),plndlv*ppcnst)
     call wrtout_r8 (nrg,divt3dsav(1,1,beglat),plndlv)       
     call wrtout_r8 (nrg,t3sav(1,1,beglat),plndlv)       
     call wrtout_r8 (nrg,u3sav(1,1,beglat),plndlv)
     call wrtout_r8 (nrg,v3sav(1,1,beglat),plndlv)
     call wrtout_r8 (nrg,t2sav(1,1,beglat),plndlv)
     call wrtout_r8 (nrg,q3sav(1,1,1,beglat),plndlv*ppcnst)
     call wrtout_r8 (nrg,pssav(1,beglat),plond)
     call wrtout_r8 (nrg,tssav(1,beglat),plond)
     call wrtout_r8 (nrg,fixmassav(beglat),1)
     call wrtout_r8 (nrg,betasav(beglat),1)
#endif
      return
   end subroutine write_restart_dynamics

!#######################################################################

   subroutine read_restart_dynamics (nrg)

#if ( defined SPMD )
      use mpishorthand
#endif

#include <comqfl.h>
!
! Input arguments
!
      integer, intent(in) :: nrg     ! Unit number
!
! Local workspace
!
      integer :: begj    ! starting latitude
      integer :: ioerr   ! error status
!
      call initialize_prognostics
      call readin_r8 (nrg,vort(1,1,beglat,n3m1), plndlv)
      call readin_r8 (nrg,vort(1,1,beglat,n3m2), plndlv)

      call readin_r8 (nrg,div(1,1,beglat,n3m1) , plndlv)
      call readin_r8 (nrg,div(1,1,beglat,n3m2) , plndlv)

      call readin_r8 (nrg,dpsl  ,plond )
      call readin_r8 (nrg,dpsm  ,plond )
      call readin_r8 (nrg,dps   ,plond )
      call readin_r8 (nrg,phis  ,plond )
      call readin_r8 (nrg,omga  ,plndlv)
!
! Write fields u3,v3,t3,q3,ps at time indices n3 and n3m1
!
      begj = beglatex + numbnd

      call readin_r8 (nrg,u3(1,1,begj,n3m1)  ,plndlv)
      call readin_r8 (nrg,v3(1,1,begj,n3m1)  ,plndlv)
      call readin_r8 (nrg,t3(1,1,begj,n3m1)  ,plndlv)
      call readin_r8 (nrg,ps(1,beglat,n3m1)  ,plond)

      call readin_r8 (nrg,u3(1,1,begj,n3m2)  ,plndlv)
      call readin_r8 (nrg,v3(1,1,begj,n3m2)  ,plndlv)
      call readin_r8 (nrg,t3(1,1,begj,n3m2)  ,plndlv)
      call readin_r8 (nrg,ps(1,beglat,n3m2)  ,plond)
      
      call readin_r8 (nrg,q3(1,1,1,begj,n3m1),plndlv*(pcnst+pnats))
      call readin_r8 (nrg,q3(1,1,1,begj,n3m2),plndlv*(pcnst+pnats))
!
! Write slt arrays (trajectory mid-point coordinates and 
! slt forcast of moisture and constituents
!
      call initialize_comslt
      call readin_r8 (nrg,lammp,plnlv)
      call readin_r8 (nrg,phimp,plnlv)
      call readin_r8 (nrg,sigmp,plnlv)
      call readin_r8 (nrg,qfcst,plndlv*pcnst)
!
! Read global integrals
!
      if (masterproc) then
         read (nrg, iostat=ioerr) tmass0, fixmas, hw1,    hw2,  &
                                  hw3, alpha
         if (ioerr /= 0 ) then
            write (6,*) 'WRITE ioerror ',ioerr,' on i/o unit = ',nrg
            call endrun ('READ_RESTART_DYNAMICS')
         end if
      end if

#if ( defined SPMD )
   call mpibcast (tmass0,1         ,mpir8  ,0,mpicom)      
   call mpibcast (fixmas,1         ,mpir8  ,0,mpicom)
   call mpibcast (hw1   ,pcnst     ,mpir8  ,0,mpicom)
   call mpibcast (hw2   ,pcnst     ,mpir8  ,0,mpicom)
   call mpibcast (hw3   ,pcnst     ,mpir8  ,0,mpicom)   
   call mpibcast (alpha ,pcnst     ,mpir8  ,0,mpicom)
#endif
#if ( defined BFB_CAM_SCAM_IOP )
!
! Read scam values
!
   allocate (betasav(beglat:endlat))
   allocate (fixmassav(beglat:endlat))
   allocate (alphasav(pcnst,beglat:endlat))
   allocate (clat_plond(plond))          ! d(ps)/dt
   allocate (alpha_plond(plond,pcnst))
   allocate (fixmas_plond(plond))
   allocate (beta_plond(plond))
   allocate (dqfx3sav(plond,plev,pcnst))
   allocate (dqfx3savm1(plond,plev,pcnst,beglat:endlat))
   allocate (divq3dsav(plond,plev,ppcnst,beglat:endlat))
   allocate (divt3dsav(plond,plev,beglat:endlat))
   allocate (t3sav(plond,plev,beglat:endlat))
   allocate (u3sav(plond,plev,beglat:endlat))
   allocate (v3sav(plond,plev,beglat:endlat))
   allocate (t2sav(plond,plev,beglat:endlat))  ! temp tendency
   allocate (q3sav(plond,plev,ppcnst,beglat:endlat))
   allocate (pssav(plond,beglat:endlat))
   allocate (tssav(plond,beglat:endlat))



     call readin_r8 (nrg,alphasav(1,beglat),pcnst)
     call readin_r8 (nrg,dqfx3savm1(1,1,1,beglat),plndlv*pcnst)       
     call readin_r8 (nrg,divq3dsav(1,1,1,beglat),plndlv*ppcnst)
     call readin_r8 (nrg,divt3dsav(1,1,beglat),plndlv)       
     call readin_r8 (nrg,t3sav(1,1,beglat),plndlv)       
     call readin_r8 (nrg,u3sav(1,1,beglat),plndlv)
     call readin_r8 (nrg,v3sav(1,1,beglat),plndlv)
     call readin_r8 (nrg,t2sav(1,1,beglat),plndlv)
     call readin_r8 (nrg,q3sav(1,1,1,beglat),plndlv*ppcnst)
     call readin_r8 (nrg,pssav(1,beglat),plond)
     call readin_r8 (nrg,tssav(1,beglat),plond)
     call readin_r8 (nrg,fixmassav(beglat),1)
     call readin_r8 (nrg,betasav(beglat),1)
#endif
      return

   end subroutine read_restart_dynamics

end module restart_dynamics
