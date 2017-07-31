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
! $Id: realloc4.F90,v 1.3.22.6 2004/04/28 23:43:56 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------
subroutine realloc4a(nlon_fft_in, nlon_fft_out, fftbuf_in, fftbuf_out )

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
! Standardized:      J. Rosinski, Oct 1995
!                    J. Truesdale, Feb. 1996
! Modified:          P. Worley, September 2002, December 2003
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
   integer, intent(in) :: nlon_fft_in      ! first dimension of input array
   integer, intent(in) :: nlon_fft_out     ! first dimension of output array
   real(r8), intent(in)  :: fftbuf_in(nlon_fft_in,plev,9,beglat:endlat) 
                            ! buffer used for in-place FFTs
   real(r8), intent(out) :: fftbuf_out(nlon_fft_out,plev,9,plat) 
                            ! buffer used for reordered Fourier coefficients
!
!---------------------------Local workspace-----------------------------
!
! xxx_l: local decomposition
! xxx_r: remote decomposition
   integer procid
   integer length_r, length_l
   integer bpos
   integer step, ifld, k, i
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
      sndcnts(procid) = length_r*(plev*8 + 1)*numlats
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
      rcvcnts(procid) = length_l*(plev*8 + 1)*nlat_p(procid)
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
!$omp parallel do private(lat_l, ifld, k, i)
!DIR$ NEXTSCALAR, NOSTREAM
      do lat_l=beglat,endlat
!DIR$ STREAM
         do ifld=1,8
!DIR$ PREFERVECTOR, PREFERSTREAM
            do k=1,plev
               do i=1,length_l
                  fftbuf_out(i,k,ifld,lat_l) = fftbuf_in(locrm(i,iam),k,ifld,lat_l)
               enddo
            enddo
         enddo
         do i=1,length_l
            fftbuf_out(i,1,9,lat_l) = fftbuf_in(locrm(i,iam),1,9,lat_l)
         enddo
      enddo
!
! Fill message buffer
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, LAT_L, IFLD, K, I)
!CSD$ PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, LAT_L, IFLD, K, I)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
!
         bpos = sdispls(procid)
         do lat_l=beglat,endlat
!DIR$ CONCURRENT, PREFERVECTOR
            do ifld=1,8
!DIR$ CONCURRENT, PREFERVECTOR
               do k=1,plev
                  do i=1,length_r
                     buf1(bpos+i) = fftbuf_in(locrm(i,procid),k,ifld,lat_l)
                  enddo
                  bpos = bpos+length_r
               enddo
            enddo
            do i=1,length_r
               buf1(bpos+i) = fftbuf_in(locrm(i,procid),1,9,lat_l)
            enddo
            bpos = bpos+length_r
         enddo
      enddo
!CSD$ END PARALLEL DO
!
! Get remote data
!
      call mpialltoallv(buf1, sndcnts, sdispls, mpir8, &
                        buf2, rcvcnts, rdispls, mpir8, &
                        mpicom)
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, IFLD, K, I)
!CSD$ PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, IFLD, K, I)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         beglat_r = cut(1,procid)
         endlat_r = cut(2,procid)
         bpos = rdispls(procid)
         do lat_r=beglat_r,endlat_r
!DIR$ CONCURRENT, PREFERVECTOR
            do ifld=1,8
!DIR$ CONCURRENT, PREFERVECTOR
               do k=1,plev
                  do i=1,length_l
                     fftbuf_out(i,k,ifld,lat_r) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_l
               enddo
            enddo
            do i=1,length_l
               fftbuf_out(i,1,9,lat_r) = buf2(bpos+i)
            enddo
            bpos = bpos+length_l
         enddo
!
      end do
!CSD$ END PARALLEL DO
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
!DIR$ NEXTSCALAR, NOSTREAM
      do lat_l=beglat,endlat
!DIR$ STREAM
         do ifld=1,8
!DIR$ PREFERVECTOR, PREFERSTREAM
            do k=1,plev
               do i=1,length_l
                  fftbuf_out(i,k,ifld,lat_l) = fftbuf_in(locrm(i,iam),k,ifld,lat_l)
               enddo
            enddo
         enddo
         do i=1,length_l
            fftbuf_out(i,1,9,lat_l) = fftbuf_in(locrm(i,iam),1,9,lat_l)
         enddo
      enddo
!
! Send data, receive data (depending on communication protocol)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
!
         bpos = sdispls(procid)
         do lat_l=beglat,endlat
!DIR$ CONCURRENT, PREFERVECTOR
            do ifld=1,8
!DIR$ CONCURRENT, PREFERVECTOR
               do k=1,plev
                  do i=1,length_r
                     buf1(bpos+i) = fftbuf_in(locrm(i,procid),k,ifld,lat_l)
                  enddo
                  bpos = bpos+length_r
               enddo
            enddo
            do i=1,length_r
               buf1(bpos+i) = fftbuf_in(locrm(i,procid),1,9,lat_l)
            enddo
            bpos = bpos+length_r
         enddo

         offset_s = sdispls(procid)+1
         offset_r = rdispls(procid)+1
         call swap2(msgtag, procid, sndcnts(procid), buf1(offset_s), sndids(step), &
                    rcvcnts(procid), buf2(offset_r), rcvids(step))

         if (.not. delayed_recv) then
            beglat_r = cut(1,procid)
            endlat_r = cut(2,procid)
            bpos = rdispls(procid)
            do lat_r=beglat_r,endlat_r
!DIR$ CONCURRENT, PREFERVECTOR
               do ifld=1,8
!DIR$ CONCURRENT, PREFERVECTOR
                  do k=1,plev
                     do i=1,length_l
                        fftbuf_out(i,k,ifld,lat_r) = buf2(bpos+i)
                     enddo
                     bpos = bpos+length_l
                  enddo
               enddo
               do i=1,length_l
                  fftbuf_out(i,1,9,lat_r) = buf2(bpos+i)
               enddo
               bpos = bpos+length_l
            enddo
         endif
!
      enddo
!
! Wait for any outstanding send or receive requests to complete.
      call swap3m(realloc4_steps, msgtag, realloc4_proc, sndids, rcvcnts, & 
                  rdispls, spmdbuf_siz, buf2, rcvids)
!
      if ( delayed_recv ) then
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, IFLD, K, I)
!CSD$ PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, IFLD, K, I)
         do step=1,realloc4_steps
            procid = realloc4_proc(step)
            beglat_r = cut(1,procid)
            endlat_r = cut(2,procid)
            bpos = rdispls(procid)
            do lat_r=beglat_r,endlat_r
!DIR$ CONCURRENT, PREFERVECTOR
               do ifld=1,8
!DIR$ CONCURRENT, PREFERVECTOR
                  do k=1,plev
                     do i=1,length_l
                        fftbuf_out(i,k,ifld,lat_r) = buf2(bpos+i)
                     enddo
                     bpos = bpos+length_l
                  enddo
               enddo
               do i=1,length_l
                  fftbuf_out(i,1,9,lat_r) = buf2(bpos+i)
               enddo
               bpos = bpos+length_l
            enddo
!
         end do
!CSD$ END PARALLEL DO
      endif
   endif
#endif
   return
   end subroutine realloc4a

subroutine realloc4b(nlon_fft_in, nlon_fft_out, fftbuf_in, fftbuf_out )

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
   integer, intent(in) :: nlon_fft_in      ! first dimension of input array
   integer, intent(in) :: nlon_fft_out     ! first dimension of output array
   real(r8), intent(in)  :: fftbuf_in(nlon_fft_in,8,plevp,plat) 
                            ! buffer of Fourier coefficients to be reordered
   real(r8), intent(out) :: fftbuf_out(nlon_fft_out,8,plevp,beglat:endlat) 
                            ! buffer used for in-place FFTs
!
!---------------------------Local workspace-----------------------------
!
! xxx_l: local decomposition
! xxx_r: remote decomposition
   integer procid
   integer length_r, length_l
   integer bpos
   integer step, ifld, k, i
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
      sndcnts(procid) = length_l*(8*plev + 4)*nlat_p(procid)
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
      rcvcnts(procid) = length_r*(8*plev + 4)*numlats
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
!$omp parallel do private(lat_l, ifld, k, i)
!DIR$ NEXTSCALAR, NOSTREAM
      do lat_l=beglat,endlat
!DIR$ STREAM
         do k=1,plev
!DIR$ PREFERVECTOR, PREFERSTREAM
            do ifld=1,8
               do i=1,length_l
                  fftbuf_out(locrm(i,iam),ifld,k,lat_l) = fftbuf_in(i,ifld,k,lat_l)
               enddo
            enddo
         enddo
!
         do ifld=1,4
            do i=1,length_l
               fftbuf_out(locrm(i,iam),ifld,plevp,lat_l) = fftbuf_in(i,ifld,plevp,lat_l)
            enddo
         enddo
      enddo
!
! Fill message buffer
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, K, IFLD, I)
!CSD$ PARALLEL DO PRIVATE (STEP, PROCID, BEGLAT_R, ENDLAT_R, BPOS, LAT_R, K, IFLD, I)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         beglat_r = cut(1,procid)
         endlat_r = cut(2,procid)
         bpos = sdispls(procid)
!
         do lat_r=beglat_r,endlat_r
!DIR$ CONCURRENT, PREFERVECTOR
            do k=1,plev
!DIR$ CONCURRENT, PREFERVECTOR
               do ifld=1,8
                  do i=1,length_l
                     buf1(bpos+i) = fftbuf_in(i,ifld,k,lat_r)
                  enddo
                  bpos = bpos+length_l
               enddo
            enddo
            do ifld=1,4
               do i=1,length_l
                  buf1(bpos+i) = fftbuf_in(i,ifld,plevp,lat_r)
               enddo
               bpos = bpos+length_l
            enddo
         enddo
      enddo
!CSD$ END PARALLEL DO
!
! Get remote data
!
      call mpialltoallv(buf1, sndcnts, sdispls, mpir8, &
                        buf2, rcvcnts, rdispls, mpir8, &
                        mpicom)
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, LAT_L, K, IFLD, I)
!CSD$ PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, LAT_L, K, IFLD, I)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
         bpos = rdispls(procid)

         do lat_l=beglat,endlat
!DIR$ CONCURRENT, PREFERVECTOR
            do k=1,plev
!DIR$ CONCURRENT, PREFERVECTOR
               do ifld=1,8
                  do i=1,length_r
                     fftbuf_out(locrm(i,procid),ifld,k,lat_l) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_r
               enddo
            enddo

            do ifld=1,4
               do i=1,length_r
                  fftbuf_out(locrm(i,procid),ifld,plevp,lat_l) = buf2(bpos+i)
               enddo
               bpos = bpos+length_r
            enddo

         enddo
!
      end do
!CSD$ END PARALLEL DO
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
!DIR$ NEXTSCALAR, NOSTREAM
      do lat_l=beglat,endlat
!DIR$ STREAM
         do k=1,plev
!DIR$ PREFERVECTOR, PREFERSTREAM
            do ifld=1,8
               do i=1,length_l
                  fftbuf_out(locrm(i,iam),ifld,k,lat_l) = fftbuf_in(i,ifld,k,lat_l)
               enddo
            enddo
         enddo

         do ifld=1,4
            do i=1,length_l
               fftbuf_out(locrm(i,iam),ifld,plevp,lat_l) = fftbuf_in(i,ifld,plevp,lat_l)
            enddo
         enddo
      enddo
!
! Send data, receive data (depending on communication protocol)
      do step=1,realloc4_steps
         procid = realloc4_proc(step)
         length_r = 2*numm(procid)
!
         beglat_r = cut(1,procid)
         endlat_r = cut(2,procid)
         bpos = sdispls(procid)
         do lat_r=beglat_r,endlat_r
!DIR$ CONCURRENT, PREFERVECTOR
            do k=1,plev
!DIR$ CONCURRENT, PREFERVECTOR
               do ifld=1,8
                  do i=1,length_l
                     buf1(bpos+i) = fftbuf_in(i,ifld,k,lat_r)
                  enddo
                  bpos = bpos+length_l
               enddo
            enddo

            do ifld=1,4
               do i=1,length_l
                  buf1(bpos+i) = fftbuf_in(i,ifld,plevp,lat_r)
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

!DIR$ NEXTSCALAR, NOSTREAM
            do lat_l=beglat,endlat
!DIR$ CONCURRENT, PREFERVECTOR, STREAM
               do k=1,plev
!DIR$ CONCURRENT, PREFERVECTOR
                  do ifld=1,8
                     do i=1,length_r
                        fftbuf_out(locrm(i,procid),ifld,k,lat_l) = buf2(bpos+i)
                     enddo
                     bpos = bpos+length_r
                  enddo
               enddo

               do ifld=1,4
                  do i=1,length_r
                     fftbuf_out(locrm(i,procid),ifld,plevp,lat_l) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_r
               enddo

            enddo
!
         endif
!
      enddo
!
! Wait for any outstanding send or receive requests to complete.
      call swap3m(realloc4_steps, msgtag, realloc4_proc, sndids, rcvcnts, & 
                  rdispls, spmdbuf_siz, buf2, rcvids)
!
      if ( delayed_recv ) then
!
! Copy out of message buffers
!
!$OMP PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, LAT_L, K, IFLD, I)
!CSD$ PARALLEL DO PRIVATE (STEP, PROCID, LENGTH_R, BPOS, LAT_L, K, IFLD, I)
         do step=1,realloc4_steps
            procid = realloc4_proc(step)
            length_r = 2*numm(procid)
            bpos = rdispls(procid)

            do lat_l=beglat,endlat
!DIR$ CONCURRENT, PREFERVECTOR
               do k=1,plev
!DIR$ CONCURRENT, PREFERVECTOR
                  do ifld=1,8
                     do i=1,length_r
                        fftbuf_out(locrm(i,procid),ifld,k,lat_l) = buf2(bpos+i)
                     enddo
                     bpos = bpos+length_r
                  enddo
               enddo

               do ifld=1,4
                  do i=1,length_r
                     fftbuf_out(locrm(i,procid),ifld,plevp,lat_l) = buf2(bpos+i)
                  enddo
                  bpos = bpos+length_r
               enddo
!
            enddo
         end do
!CSD$ END PARALLEL DO
      endif
   endif
!
#endif
   return
   end subroutine realloc4b

