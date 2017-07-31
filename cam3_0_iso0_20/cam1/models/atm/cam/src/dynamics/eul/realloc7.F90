#include <misc.h>
#include <params.h>

subroutine realloc7 (vmax2d, vmax2dt, vcour)

!----------------------------------------------------------------------- 
! 
! Purpose: Reallocation routine for energy and log stats
! 
! Method: MPI_Allgatherv (or point-to-point implementation)
! 
! Author: J. Rosinski
! Modified: P. Worley, September 2002, December 2003
! 
!-----------------------------------------------------------------------

#ifdef SPMD
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: plat, plev, iam, numlats, beglat, endlat
   use mpishorthand
   use spmd_dyn
   use swap_comm
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comsta.h>
!------------------------------Parameters-------------------------------
!
   integer, parameter :: msgtag  = 3000
!---------------------------Input arguments-----------------------------
!
   real(r8), intent(inout) :: vmax2d(plev,plat)   ! Max. wind at each lvl, lat
   real(r8), intent(inout) :: vmax2dt(plev,plat)  ! Max. truncated wind at each lvl, lat
   real(r8), intent(inout) :: vcour(plev,plat)    ! Max. Courant number at each lvl, lat
!
!---------------------------Local workspace-----------------------------
!
   integer procid
   integer bufpos
   integer procj
   integer step, j, k, jstrt
   integer beglat_p, endlat_p, numlats_p, jstrt_p
   integer offset_r               ! receive displacement + 1
   integer sndcnt
   integer sndids(npes-1)         ! nonblocking MPI send request ids
   integer rcvids(npes-1)         ! nonblocking MPI recv request ids
   integer rcvcnts(0:npes-1), rdispls(0:npes-1)
   logical delayed_recv           ! local copy of delayed_swap_recv flag
!-----------------------------------------------------------------------
!
! Compute send count
   sndcnt = (plev*3 + 5)*numlats
!
! Compute recv counts and displacements
   rcvcnts(:) = 0
   do step=1,npes-1
      procid = allgather_proc(step)
      rcvcnts(procid) = (plev*3 + 5)*nlat_p(procid)
   enddo
   rcvcnts(iam) = (plev*3 + 5)*numlats
!   
   rdispls(0) = 0
   do procid=1,npes-1
     rdispls(procid) = rdispls(procid-1) + rcvcnts(procid-1)
   enddo
!
   if (dyn_allgather .eq. 0) then
!
! Fill send buffer
      jstrt = beglat - 1
      bufpos = 0
! psurf
      do j=1,numlats
         buf1(bufpos+j) = psurf(jstrt+j)
      enddo
      bufpos = bufpos + numlats
! stq
      do j=1,numlats
         buf1(bufpos+j) = stq(jstrt+j)
      enddo
      bufpos = bufpos + numlats
! rmst
      do j=1,numlats
         buf1(bufpos+j) = rmst(jstrt+j)
      enddo
      bufpos = bufpos + numlats
! rmsd
      do j=1,numlats
         buf1(bufpos+j) = rmsd(jstrt+j)
      enddo
      bufpos = bufpos + numlats
! rmsz
      do j=1,numlats
         buf1(bufpos+j) = rmsz(jstrt+j)
      enddo
      bufpos = bufpos + numlats
!vmax2d
      do j=beglat,endlat
         do k=1,plev
            buf1(bufpos+k) = vmax2d(k,j)
         enddo
         bufpos = bufpos + plev
      enddo
! vmax2dt
      do j=beglat,endlat
         do k=1,plev
            buf1(bufpos+k) = vmax2dt(k,j)
         enddo
         bufpos = bufpos + plev
      enddo
! vcour
      do j=beglat,endlat
         do k=1,plev
            buf1(bufpos+k) = vcour(k,j)
         enddo
         bufpos = bufpos + plev
      enddo
!
! Gather the data
!
      call mpiallgatherv(buf1, sndcnt, mpir8, &
                         buf2, rcvcnts, rdispls, mpir8, &
                         mpicom)
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_P, ENDLAT_P, NUMLATS_P, BUFPOS, J, K)
!DIR$ CONCURRENT
      do step=1,npes-1
         procid = allgather_proc(step)
         beglat_p = cut(1,procid)
         endlat_p = cut(2,procid)
         numlats_p = nlat_p(procid)
         bufpos = rdispls(procid)
! psurf
         jstrt_p  = beglat_p - 1
         do j=1,numlats_p
            psurf(jstrt_p+j) = buf2(bufpos+j)
         enddo
         bufpos = bufpos + numlats_p
! stq
         do j=1,numlats_p
            stq(jstrt_p+j) = buf2(bufpos+j)
         enddo
         bufpos = bufpos + numlats_p
! rmst
         do j=1,numlats_p
            rmst(jstrt_p+j) = buf2(bufpos+j)
         enddo
         bufpos = bufpos + numlats_p
! rmsd
         do j=1,numlats_p
            rmsd(jstrt_p+j) = buf2(bufpos+j) 
         enddo
         bufpos = bufpos + numlats_p
! rmsz
         do j=1,numlats_p
            rmsz(jstrt_p+j) = buf2(bufpos+j) 
         enddo
         bufpos = bufpos + numlats_p
! vmax2d
         do j=beglat_p,endlat_p
            do k=1,plev
               vmax2d(k,j) = buf2(bufpos+k)
            enddo
            bufpos = bufpos + plev
         enddo
! vmax2dt
         do j=beglat_p,endlat_p
            do k=1,plev
               vmax2dt(k,j) = buf2(bufpos+k)
            enddo
            bufpos = bufpos + plev
         enddo
! vcour
         do j=beglat_p,endlat_p
            do k=1,plev
               vcour(k,j) = buf2(bufpos+k)
            enddo
            bufpos = bufpos + plev
         enddo
!
      enddo
!
   else
!
      delayed_recv = delayed_swap_recv()
!
! Post receive requests
      call swap1m(npes-1, msgtag, allgather_proc, rcvcnts, & 
                  rdispls, spmdbuf_siz, buf2, rcvids)
!
! Fill send buffer
      jstrt = beglat - 1
      bufpos = 0
! psurf
      do j=1,numlats
         buf1(bufpos+j) = psurf(jstrt+j)
      enddo
      bufpos = bufpos + numlats
! stq
      do j=1,numlats
         buf1(bufpos+j) = stq(jstrt+j)
      enddo
      bufpos = bufpos + numlats
! rmst
      do j=1,numlats
         buf1(bufpos+j) = rmst(jstrt+j)
      enddo
      bufpos = bufpos + numlats
! rmsd
      do j=1,numlats
         buf1(bufpos+j) = rmsd(jstrt+j)
      enddo
      bufpos = bufpos + numlats
! rmsz
      do j=1,numlats
         buf1(bufpos+j) = rmsz(jstrt+j)
      enddo
      bufpos = bufpos + numlats
!vmax2d
      do j=beglat,endlat
         do k=1,plev
            buf1(bufpos+k) = vmax2d(k,j)
         enddo
         bufpos = bufpos + plev
      enddo
! vmax2dt
      do j=beglat,endlat
         do k=1,plev
            buf1(bufpos+k) = vmax2dt(k,j)
         enddo
         bufpos = bufpos + plev
      enddo
! vcour
      do j=beglat,endlat
         do k=1,plev
            buf1(bufpos+k) = vcour(k,j)
         enddo
         bufpos = bufpos + plev
      enddo
!
! Send data, receive data (depending on communication protocol)
      do step=1,npes-1
         procid = allgather_proc(step)
!
         offset_r = rdispls(procid)+1
         call swap2(msgtag, procid, sndcnt, buf1, sndids(step), &
                    rcvcnts(procid), buf2(offset_r), rcvids(step))
!
         if (.not. delayed_recv) then
            beglat_p = cut(1,procid)
            endlat_p = cut(2,procid)
            numlats_p = nlat_p(procid)
            bufpos = rdispls(procid)
! psurf
            jstrt_p  = beglat_p - 1
            do j=1,numlats_p
               psurf(jstrt_p+j) = buf2(bufpos+j)
            enddo
            bufpos = bufpos + numlats_p
! stq
            do j=1,numlats_p
               stq(jstrt_p+j) = buf2(bufpos+j)
            enddo
            bufpos = bufpos + numlats_p
! rmst
            do j=1,numlats_p
               rmst(jstrt_p+j) = buf2(bufpos+j)
            enddo
            bufpos = bufpos + numlats_p
! rmsd
            do j=1,numlats_p
               rmsd(jstrt_p+j) = buf2(bufpos+j) 
            enddo
            bufpos = bufpos + numlats_p
! rmsz
            do j=1,numlats_p
               rmsz(jstrt_p+j) = buf2(bufpos+j) 
            enddo
            bufpos = bufpos + numlats_p
! vmax2d
            do j=beglat_p,endlat_p
               do k=1,plev
                  vmax2d(k,j) = buf2(bufpos+k)
               enddo
               bufpos = bufpos + plev
            enddo
! vmax2dt
            do j=beglat_p,endlat_p
               do k=1,plev
                  vmax2dt(k,j) = buf2(bufpos+k)
               enddo
               bufpos = bufpos + plev
            enddo
! vcour
            do j=beglat_p,endlat_p
               do k=1,plev
                  vcour(k,j) = buf2(bufpos+k)
               enddo
               bufpos = bufpos + plev
            enddo
!
         endif

      enddo
!
! Wait for any outstanding send or receive requests to complete.
      call swap3m(npes-1, msgtag, allgather_proc, sndids, rcvcnts, & 
                  rdispls, spmdbuf_siz, buf2, rcvids)
!
      if ( delayed_recv ) then
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_P, ENDLAT_P, NUMLATS_P, BUFPOS, J, K)
!DIR$ CONCURRENT
         do step=1,npes-1
            procid = allgather_proc(step)
            beglat_p = cut(1,procid)
            endlat_p = cut(2,procid)
            numlats_p = nlat_p(procid)
            bufpos = rdispls(procid)
! psurf
            jstrt_p  = beglat_p - 1
            do j=1,numlats_p
               psurf(jstrt_p+j) = buf2(bufpos+j)
            enddo
            bufpos = bufpos + numlats_p
! stq
            do j=1,numlats_p
               stq(jstrt_p+j) = buf2(bufpos+j)
            enddo
            bufpos = bufpos + numlats_p
! rmst
            do j=1,numlats_p
               rmst(jstrt_p+j) = buf2(bufpos+j)
            enddo
            bufpos = bufpos + numlats_p
! rmsd
            do j=1,numlats_p
               rmsd(jstrt_p+j) = buf2(bufpos+j) 
            enddo
            bufpos = bufpos + numlats_p
! rmsz
            do j=1,numlats_p
               rmsz(jstrt_p+j) = buf2(bufpos+j) 
            enddo
            bufpos = bufpos + numlats_p
! vmax2d
            do j=beglat_p,endlat_p
               do k=1,plev
                  vmax2d(k,j) = buf2(bufpos+k)
               enddo
               bufpos = bufpos + plev
            enddo
! vmax2dt
            do j=beglat_p,endlat_p
               do k=1,plev
                  vmax2dt(k,j) = buf2(bufpos+k)
               enddo
               bufpos = bufpos + plev
            enddo
! vcour
            do j=beglat_p,endlat_p
               do k=1,plev
                  vcour(k,j) = buf2(bufpos+k)
               enddo
               bufpos = bufpos + plev
            enddo
!
         enddo
!
      endif
!
   endif
#endif
   return
end subroutine realloc7

