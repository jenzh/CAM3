#include <misc.h> 
#include <params.h>

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reallocation routines for the Fourier coefficients
! 
! Method: 
!   1) After FFT preceding Legendre analysis, reallocate fftbuf
!      to decompose over wavenumber, recombining latitudes.
!   2) Before FFT following Legendre synthesis, reallocate fftbuf
!      to recombine wavenumbers, decomposing over latitude.
! 
!-----------------------------------------------------------------------
!
! $Id: realloc4.F90,v 1.4.22.4 2004/03/03 19:54:11 pworley Exp $
! $Author: pworley $
!
!-----------------------------------------------------------------------

subroutine realloc4a(fftbuf_in, fftbuf_out )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reallocation routines for the Fourier coefficients
! 
! Method: 
!   After FFT preceding Legendre analysis, reallocate fftbuf
!   to decompose over wavenumber, combining latitudes.
! 
! Author: 
! Original version:  J. Rosinski
! Modified:          P. Worley, October 2002, December 2003
! 
!-----------------------------------------------------------------------

#ifdef SPMD

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use spmd_dyn
   use mpishorthand
   use swap_comm
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comsta.h>
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtag  = 1000
!---------------------------Input arguments-----------------------------
!
   real(r8), intent(in)  :: fftbuf_in(plond,plev,5,beglat:endlat) 
                            ! buffer used for in-place FFTs
   real(r8), intent(out) :: fftbuf_out(2*maxm,plev,5,plat) 
                            ! buffer used for reordered Fourier coefficients
!
!---------------------------Local workspace-----------------------------
!
! xxx_l: local decomposition
! xxx_r: remote decomposition
   integer procid
   integer length_r, length_l, locrm(2*maxm)
   integer bpos
   integer step, if, k, i
   integer lat_l, lat_r, beglat_r, endlat_r
   integer offset_s               ! send displacement + 1
   integer offset_r               ! receive displacement + 1
   integer sndids(realloc4_steps) ! nonblocking MPI send request ids
   integer rcvids(realloc4_steps) ! nonblocking MPI recv request ids
   integer sndcnts(0:npes-1), sdispls(0:npes-1)
   integer rcvcnts(0:npes-1), rdispls(0:npes-1)
   logical delayed_recv           ! local copy of delayed_swap_recv flag
!-----------------------------------------------------------------------
! Compute send/recv counts and displacements
   sndcnts(:) = 0
   do step=1,realloc4_steps
      procid = realloc4_proc(step)
      length_r = 2*numm(procid)
      sndcnts(procid) = length_r*(plev*4 + 1)*numlats
   enddo
!   
   sdispls(0) = 0
   do procid=1,npes-1
     sdispls(procid) = sdispls(procid-1) + sndcnts(procid-1)
   enddo
!
   length_l = 2*numm(iam)
   rcvcnts(:) = 0
   do step=1,realloc4_steps
      procid = realloc4_proc(step)
      rcvcnts(procid) = length_l*(plev*4 + 1)*nlat_p(procid)
   enddo
!   
   rdispls(0) = 0
   do procid=1,npes-1
     rdispls(procid) = rdispls(procid-1) + rcvcnts(procid-1)
   enddo
!
   if (dyn_alltoall .eq. 0) then
!     
! Copy local data to new location
      length_l = 2*numm(iam)
      do i=1,numm(iam)
         locrm(2*i-1) = 2*locm(i,iam)-1
         locrm(2*i)   = 2*locm(i,iam)
      enddo
      do lat_l=beglat,endlat
         do if=1,4
            do k=1,plev
               do i=1,length_l
                  fftbuf_out(i,k,if,lat_l) = fftbuf_in(locrm(i),k,if,lat_l)
               enddo
            enddo
         enddo
         do i=1,length_l
            fftbuf_out(i,1,5,lat_l) = fftbuf_in(locrm(i),1,5,lat_l)
         enddo
      enddo
!
! Fill message buffer
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, I, LOCRM, BPOS, LAT_L, IF, K)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
         do i=1,numm(procid)
            locrm(2*i-1) = 2*locm(i,procid)-1
            locrm(2*i)   = 2*locm(i,procid)
         enddo
!
         bpos = sdispls(procid)
         do lat_l=beglat,endlat
            do if=1,4
               do k=1,plev
                  do i=1,length_r
                     buf1(bpos+i) = fftbuf_in(locrm(i),k,if,lat_l)
                  enddo
                  bpos = bpos+length_r
               enddo
            enddo
            do i=1,length_r
               buf1(bpos+i) = fftbuf_in(locrm(i),1,5,lat_l)
            enddo
            bpos = bpos+length_r
         enddo
      enddo
!
! Get remote data
!
      call mpialltoallv(buf1, sndcnts, sdispls, mpir8, &
                        buf2, rcvcnts, rdispls, mpir8, &
                        mpicom)
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, IF, K, I)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         beglat_r = cut(1,procid)
         endlat_r = cut(2,procid)
         bpos = rdispls(procid)
         do lat_r=beglat_r,endlat_r
            do if=1,4
               do k=1,plev
                  do i=1,length_l
                     fftbuf_out(i,k,if,lat_r) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_l
               enddo
            enddo
            do i=1,length_l
               fftbuf_out(i,1,5,lat_r) = buf2(bpos+i)
            enddo
            bpos = bpos+length_l
         enddo
!
      end do
   else
!
      delayed_recv = delayed_swap_recv()
!     
! Post receive requests
      call swap1m(realloc4_steps, msgtag, realloc4_proc, rcvcnts, & 
                  rdispls, spmdbuf_siz, buf2, rcvids)
!     
! Copy local data to new location
      length_l = 2*numm(iam)
      do i=1,numm(iam)
         locrm(2*i-1) = 2*locm(i,iam)-1
         locrm(2*i)   = 2*locm(i,iam)
      enddo
      do lat_l=beglat,endlat
         do if=1,4
            do k=1,plev
               do i=1,length_l
                  fftbuf_out(i,k,if,lat_l) = fftbuf_in(locrm(i),k,if,lat_l)
               enddo
            enddo
         enddo
         do i=1,length_l
            fftbuf_out(i,1,5,lat_l) = fftbuf_in(locrm(i),1,5,lat_l)
         enddo
      enddo
!
! Send data
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
!
         do i=1,numm(procid)
            locrm(2*i-1) = 2*locm(i,procid)-1
            locrm(2*i)   = 2*locm(i,procid)
         enddo
!
         bpos = sdispls(procid)
         do lat_l=beglat,endlat
            do if=1,4
               do k=1,plev
                  do i=1,length_r
                     buf1(bpos+i) = fftbuf_in(locrm(i),k,if,lat_l)
                  enddo
                  bpos = bpos+length_r
               enddo
            enddo
            do i=1,length_r
               buf1(bpos+i) = fftbuf_in(locrm(i),1,5,lat_l)
            enddo
            bpos = bpos+length_r
         enddo

         offset_s = sdispls(procid)+1
         offset_r = rdispls(procid)+1
         call swap2(msgtag, procid, sndcnts(procid), buf1(offset_s), sndids(step), &
                    rcvcnts(procid), buf2(offset_r), rcvids(step))
!
         if (.not. delayed_recv) then
            beglat_r = cut(1,procid)
            endlat_r = cut(2,procid)
            bpos = rdispls(procid)
            do lat_r=beglat_r,endlat_r
               do if=1,4
                  do k=1,plev
                     do i=1,length_l
                        fftbuf_out(i,k,if,lat_r) = buf2(bpos+i)
                     enddo
                     bpos = bpos+length_l
                  enddo
               enddo
               do i=1,length_l
                  fftbuf_out(i,1,5,lat_r) = buf2(bpos+i)
               enddo
               bpos = bpos+length_l
            enddo
         endif
!
      enddo
!
! Wait for send requests to complete.
      call swap3m(realloc4_steps, msgtag, realloc4_proc, sndids, rcvcnts, & 
                  rdispls, spmdbuf_siz, buf2, rcvids)
!
      if ( delayed_recv ) then
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, IF, K, I)
         do step=1,realloc4_steps
            procid = realloc4_proc(step)
            beglat_r = cut(1,procid)
            endlat_r = cut(2,procid)
            bpos = rdispls(procid)
            do lat_r=beglat_r,endlat_r
               do if=1,4
                  do k=1,plev
                     do i=1,length_l
                        fftbuf_out(i,k,if,lat_r) = buf2(bpos+i)
                     enddo
                     bpos = bpos+length_l
                  enddo
               enddo
               do i=1,length_l
                  fftbuf_out(i,1,5,lat_r) = buf2(bpos+i)
               enddo
               bpos = bpos+length_l
            enddo
!
         end do
      endif
   endif
#endif
   return
   end subroutine realloc4a
!
subroutine realloc4b(fftbuf_in, fftbuf_out )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Reallocation routines for the Fourier coefficients
! 
! Method: 
!   Before FFT following Legendre synthesis, reallocate fftbuf
!   to combine wavenumbers, decomposing over latitude.
! 
! Author:  P. Worley, September 2002
! Modified: P. Worley, December 2003
! 
!-----------------------------------------------------------------------

#ifdef SPMD

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use spmd_dyn
   use mpishorthand
   use swap_comm
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <comsta.h>
!------------------------------Parameters-------------------------------
!
  integer, parameter :: msgtag  = 2000
!---------------------------Input arguments--------------------------
!
   real(r8), intent(in)  :: fftbuf_in(2*maxm,11,plevp,plat) 
                            ! buffer of Fourier coefficients to be reordered
   real(r8), intent(out) :: fftbuf_out(plond,11,plevp,beglat:endlat) 
                            ! buffer used for in-place FFTs
!
!---------------------------Local workspace-----------------------------
!
! xxx_l: local decomposition
! xxx_r: remote decomposition
   integer procid
   integer length_r, length_l, locrm(2*maxm)
   integer bpos
   integer step, if, k, i
   integer lat_l, lat_r
   integer beglat_r, endlat_r
   integer offset_s               ! send displacement + 1
   integer offset_r               ! receive displacement + 1
   integer sndids(realloc4_steps) ! nonblocking MPI send request ids
   integer rcvids(realloc4_steps) ! nonblocking MPI recv request ids
   integer sndcnts(0:npes-1), sdispls(0:npes-1)
   integer rcvcnts(0:npes-1), rdispls(0:npes-1)
   logical delayed_recv           ! local copy of delayed_swap_recv flag
!-----------------------------------------------------------------------
! Compute send/recv counts and displacements
   length_l = 2*numm(iam)
   sndcnts(:) = 0
   do step=1,realloc4_steps
      procid = realloc4_proc(step)
      sndcnts(procid) = length_l*(11*plev + 4)*nlat_p(procid)
   enddo
!   
   sdispls(0) = 0
   do procid=1,npes-1
     sdispls(procid) = sdispls(procid-1) + sndcnts(procid-1)
   enddo
!
   rcvcnts(:) = 0
   do step=1,realloc4_steps
      procid = realloc4_proc(step)
      length_r = 2*numm(procid)
      rcvcnts(procid) = length_r*(11*plev + 4)*numlats
   enddo
!   
   rdispls(0) = 0
   do procid=1,npes-1
     rdispls(procid) = rdispls(procid-1) + rcvcnts(procid-1)
   enddo
!
   if (dyn_alltoall .eq. 0) then
!     
! Copy local data to new location
      length_l = 2*numm(iam)
      do i=1,numm(iam)
         locrm(2*i-1) = 2*locm(i,iam)-1
         locrm(2*i)   = 2*locm(i,iam)
      enddo
      do lat_l=beglat,endlat
         do k=1,plev
            do if=1,11
               do i=1,length_l
                  fftbuf_out(locrm(i),if,k,lat_l) = fftbuf_in(i,if,k,lat_l)
               enddo
            enddo
         enddo

         do if=1,4
            do i=1,length_l
               fftbuf_out(locrm(i),if,plevp,lat_l) = fftbuf_in(i,if,plevp,lat_l)
            enddo
         enddo
      enddo
!
! Fill message buffer
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, K, IF, I)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         beglat_r = cut(1,procid)
         endlat_r = cut(2,procid)
         bpos = sdispls(procid)
!
         do lat_r=beglat_r,endlat_r
            do k=1,plev
               do if=1,11
                  do i=1,length_l
                     buf1(bpos+i) = fftbuf_in(i,if,k,lat_r)
                  enddo
                  bpos = bpos+length_l
               enddo
            enddo
            do if=1,4
               do i=1,length_l
                  buf1(bpos+i) = fftbuf_in(i,if,plevp,lat_r)
               enddo
               bpos = bpos+length_l
            enddo
         enddo
      enddo
!
! Get remote data
!
      call mpialltoallv(buf1, sndcnts, sdispls, mpir8, &
                        buf2, rcvcnts, rdispls, mpir8, &
                        mpicom)
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, I, LOCRM, LAT_L, K, IF)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
         bpos = rdispls(procid)

         do i=1,numm(procid)
            locrm(2*i-1) = 2*locm(i,procid)-1
            locrm(2*i)   = 2*locm(i,procid)
         enddo

         do lat_l=beglat,endlat
            do k=1,plev
               do if=1,11
                  do i=1,length_r
                     fftbuf_out(locrm(i),if,k,lat_l) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_r
               enddo
            enddo

            do if=1,4
               do i=1,length_r
                  fftbuf_out(locrm(i),if,plevp,lat_l) = buf2(bpos+i)
               enddo
               bpos = bpos+length_r
            enddo

         enddo
!
      end do
!     
   else
!
      delayed_recv = delayed_swap_recv()
!     
! Post receive requests
      call swap1m(realloc4_steps, msgtag, realloc4_proc, rcvcnts, & 
                  rdispls, spmdbuf_siz, buf2, rcvids)
!     
! Copy local data to new location
      length_l = 2*numm(iam)
      do i=1,numm(iam)
         locrm(2*i-1) = 2*locm(i,iam)-1
         locrm(2*i)   = 2*locm(i,iam)
      enddo
      do lat_l=beglat,endlat
         do k=1,plev
            do if=1,11
               do i=1,length_l
                  fftbuf_out(locrm(i),if,k,lat_l) = fftbuf_in(i,if,k,lat_l)
               enddo
            enddo
         enddo

         do if=1,4
            do i=1,length_l
               fftbuf_out(locrm(i),if,plevp,lat_l) = fftbuf_in(i,if,plevp,lat_l)
            enddo
         enddo
      enddo
!
! Send data
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
!
         do i=1,numm(procid)
            locrm(2*i-1) = 2*locm(i,procid)-1
            locrm(2*i)   = 2*locm(i,procid)
         enddo
!
         beglat_r = cut(1,procid)
         endlat_r = cut(2,procid)
         bpos = sdispls(procid)
         do lat_r=beglat_r,endlat_r
            do k=1,plev
               do if=1,11
                  do i=1,length_l
                     buf1(bpos+i) = fftbuf_in(i,if,k,lat_r)
                  enddo
                  bpos = bpos+length_l
               enddo
            enddo

            do if=1,4
               do i=1,length_l
                  buf1(bpos+i) = fftbuf_in(i,if,plevp,lat_r)
               enddo
               bpos = bpos+length_l
            enddo

         enddo

         offset_s = sdispls(procid)+1
         offset_r = rdispls(procid)+1
         call swap2(msgtag, procid, sndcnts(procid), buf1(offset_s), sndids(step), &
                    rcvcnts(procid), buf2(offset_r), rcvids(step))
!
         if (.not. delayed_recv) then
            length_r = 2*numm(procid)
            bpos = rdispls(procid)

            do i=1,numm(procid)
               locrm(2*i-1) = 2*locm(i,procid)-1
               locrm(2*i)   = 2*locm(i,procid)
            enddo

            do lat_l=beglat,endlat
               do k=1,plev
                  do if=1,11
                     do i=1,length_r
                        fftbuf_out(locrm(i),if,k,lat_l) = buf2(bpos+i)
                     enddo
                     bpos = bpos+length_r
                  enddo
               enddo

               do if=1,4
                  do i=1,length_r
                     fftbuf_out(locrm(i),if,plevp,lat_l) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_r
               enddo

            enddo
!
         endif
!
      enddo
!
! Wait for send requests to complete.
      call swap3m(realloc4_steps, msgtag, realloc4_proc, sndids, rcvcnts, & 
                  rdispls, spmdbuf_siz, buf2, rcvids)
!
      if ( delayed_recv ) then
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, I, LOCRM, LAT_L, K, IF)
         do step=1,realloc4_steps
            procid = realloc4_proc(step)
            length_r = 2*numm(procid)
            bpos = rdispls(procid)

            do i=1,numm(procid)
               locrm(2*i-1) = 2*locm(i,procid)-1
               locrm(2*i)   = 2*locm(i,procid)
            enddo

            do lat_l=beglat,endlat
               do k=1,plev
                  do if=1,11
                     do i=1,length_r
                        fftbuf_out(locrm(i),if,k,lat_l) = buf2(bpos+i)
                     enddo
                     bpos = bpos+length_r
                  enddo
               enddo

               do if=1,4
                  do i=1,length_r
                     fftbuf_out(locrm(i),if,plevp,lat_l) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_r
               enddo
!
            enddo
         enddo
      endif
   endif
!
#endif
   return
   end subroutine realloc4b

