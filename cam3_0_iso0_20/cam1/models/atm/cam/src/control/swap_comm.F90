#include <misc.h>
#include <params.h>

module swap_comm

!----------------------------------------------------------------------- 
! 
! Purpose: swap communication routines used in performance portable 
!          distributed transpose algorithms.
! 
! Entry points:
!      swap_comm_init          Initialize swap module.
!
!      swap_comm_defaultopts   Get default runtime options.
!      swap_comm_setopts       Set runtime options.
!      
!      swap1                   First of three routines that implement swap 
!                              using MPI point-to-point. Depending on
!                              communication protocol, posts nonblocking
!                              receive requests and sends ready send
!                              handshaking messages.
!      swap2                   Second of three routines that implement swap 
!                              using MPI point-to-point. Sends data and,
!                              depending on communication protocol, receives
!                              data.
!      swap3                   Third of three routines that implement swap 
!                              using MPI point-to-point. Completes all
!                              outstanding send and receive requests.
!      swap1m                  Same as swap1, but for multiple messages
!      swap3m                  Same as swap3, but for multiple messages.
!
!      do_swap1                Logical function that indicates whether
!                              swap1 needs to be called for the current
!                              communication protocol.
!      do_swap3                Logical function that indicates whether
!                              swap3 needs to be called for the current
!                              communication protocol.
!      delayed_swap_recv       Logical function that indicates whether
!                              message receives are delayed until swap3
!                              for the current communication protocol.
!
! Author: P. Worley
!-----------------------------------------------------------------------

#if (defined SPMD)
!-----------------------------------------------------------------------
!- use statements ------------------------------------------------------
!-----------------------------------------------------------------------

   use abortutils, only: endrun
   use mpishorthand

!-----------------------------------------------------------------------
!- module boilerplate --------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
   private                   ! Make the default access private
   save

!-----------------------------------------------------------------------
! Public interfaces ----------------------------------------------------
!-----------------------------------------------------------------------
   public swap_comm_init
   public swap_comm_defaultopts 
   public swap_comm_setopts 
   public swap1
   public swap2
   public swap3
   public swap1m
   public swap3m
   public do_swap1
   public do_swap3
   public delayed_swap_recv

!-----------------------------------------------------------------------
! Private data ---------------------------------------------------------
!-----------------------------------------------------------------------
! Swap communication order option:
!  0: simple swap: send/recv
!  1: ordered swap: [send/recv]|[recv/send]
!  2: delayed-recv swap: send ... recv
   integer, private, parameter :: min_comm_order = 0
   integer, private, parameter :: max_comm_order = 2
   integer, private, parameter :: def_comm_order = 0              ! default
   integer, private :: swap_comm_order = def_comm_order

! Swap communication protocol option:
!  1, 3, 5, 7, 9:                  nonblocking send
!  2, 3, 4, 5, 8, 9:               nonblocking receive
!  4, 5:                           ready send
!  6 .and. swap_comm_order .eq. 0: sendrecv
!  6 .and. swap_comm_order .eq. 1: explicitly synchronous  
!  7, 8, 9, .or. 10:               synchronous send          
   integer, private, parameter :: min0_comm_protocol =  1
   integer, private, parameter :: max0_comm_protocol =  9
   integer, private, parameter :: min1_comm_protocol =  0
   integer, private, parameter :: max1_comm_protocol = 10
   integer, private, parameter :: def_comm_protocol  =  6        ! default
   integer, private :: swap_comm_protocol = def_comm_protocol

! Swap communicators
   integer, private :: swap_com = mpi_comm_null
                                      ! primary MPI communicator
   integer, private :: handshake_com  = mpi_comm_null
                                      ! MPI communicator for 
                                      !  handshaking messages

!-----------------------------------------------------------------------
! Subroutines and functions --------------------------------------------
!-----------------------------------------------------------------------
contains

!
!========================================================================
!
   subroutine swap_comm_init()

!----------------------------------------------------------------------- 
! 
! Purpose: Create communicators to be used in swap communication.
! 
! Method: 
!
! Author:  P. Worley
!-----------------------------------------------------------------------
   use mpishorthand, only: mpicom
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments-----------------------------
!
!---------------------------Local workspace-----------------------------
!
   integer ier               ! return error status    
!
!-----------------------------------------------------------------------
!
   call mpi_comm_dup(mpicom, swap_com, ier)
   if (ier /= mpi_success) then
      write(6,*)                                         &
         'SWAP_COMM_INIT:  ERROR:  mpi_comm_dup failed with IER=', ier
      call endrun
   endif
   call mpi_comm_dup(mpicom, handshake_com, ier)
   if (ier /= mpi_success) then
      write(6,*)                                         &
         'SWAP_COMM_INIT:  ERROR:  mpi_comm_dup failed with IER=', ier
      call endrun
   endif
!
   return
   end subroutine swap_comm_init
!
!========================================================================
!
   subroutine swap_comm_defaultopts(swap_comm_order_out, &
                                    swap_comm_protocol_out )
!----------------------------------------------------------------------- 
! Purpose: Return default runtime options
! Author: P. Worley (modelled after Tom Henderson's code)
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
   ! swap module communication order option
   integer, intent(out), optional :: swap_comm_order_out
   ! swap module communication protocol option
   integer, intent(out), optional :: swap_comm_protocol_out
!-----------------------------------------------------------------------
   if ( present(swap_comm_order_out) ) then
      swap_comm_order_out = def_comm_order
   endif
   if ( present(swap_comm_protocol_out) ) then
      swap_comm_protocol_out = def_comm_protocol
   endif
!
   return
   end subroutine swap_comm_defaultopts
!
!========================================================================
!
   subroutine swap_comm_setopts(swap_comm_order_in, &
                                swap_comm_protocol_in )
!----------------------------------------------------------------------- 
! Purpose: Set runtime options
! Author: P. Worley (modelled after Tom Henderson's code)
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
     ! swap module communication order option
     integer, intent(in), optional :: swap_comm_order_in
     ! swap module communication protocol option
     integer, intent(in), optional :: swap_comm_protocol_in
!-----------------------------------------------------------------------
     if ( present(swap_comm_order_in) ) then
        swap_comm_order = swap_comm_order_in
        if ((swap_comm_order < min_comm_order) .or. &
            (swap_comm_order > max_comm_order)) then
           write(6,*)                                         &
              'SWAP_COMM_SETOPTS:  ERROR:  swap_comm_order=', &
              swap_comm_order,                                &
              ' is out of range.  It must be between ',       &
              min_comm_order,' and ',max_comm_order
           call endrun
        endif
     endif
!
     if ( present(swap_comm_protocol_in) ) then
        swap_comm_protocol = swap_comm_protocol_in
        if ((swap_comm_order .eq. 0) .or. &
            (swap_comm_order .eq. 2)) then
           if ((swap_comm_protocol < min0_comm_protocol) .or. &
               (swap_comm_protocol > max0_comm_protocol)) then
              write(6,*)                                            &
                 'SWAP_COMM_SETOPTS:  ERROR:  swap_comm_protocol=', &
                 swap_comm_protocol,                                &
                 ' is out of range.  It must be between ',          &
                 min0_comm_protocol,' and ',max0_comm_protocol,     &
                 ' when swap_comm_order= ', swap_comm_order
              call endrun
           endif
        else
           if ((swap_comm_protocol < min1_comm_protocol) .or. &
               (swap_comm_protocol > max1_comm_protocol)) then
              write(6,*)                                            &
                 'SWAP_COMM_SETOPTS:  ERROR:  swap_comm_protocol=', &
                 swap_comm_protocol,                                &
                 ' is out of range.  It must be between ',          &
                 min1_comm_protocol,' and ',max1_comm_protocol,     &
                 ' when swap_comm_order= ', swap_comm_order
              call endrun
           endif
        endif
     endif
!
     return
   end subroutine swap_comm_setopts
!
!========================================================================
!
   subroutine swap1(mtag, swapnode, rcvlth, rcvmsg, rcvid)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! First of three routines that implement swap using MPI point-to-point
! routines.
! 
! Method: 
! This subroutine begins a swap operation that will be completed by
! swap2 and swap3. It posts a receive and sends handshaking messages
! when ready sends are used. 
! if (swap_comm_order .eq. 0) simple swap: send/recv
! if (swap_comm_order .eq. 1) ordered swap: [send/recv]|[recv/send]
! if (swap_comm_order .eq. 2) delayed-recv swap: send ... recv
! if (swap_comm_protocol .eq. 1, 3, 5, 7, .or. 9) nonblocking send
! if (swap_comm_protocol .eq. 2, 3, 4, 5, 8, .or. 9)  nonblocking receive
! if (swap_comm_protocol .eq. 4 .or. 5)  ready send
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 0) sendrecv
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 1) explicitly synchronous  
! if (swap_comm_protocol .eq. 7, 8, 9, .or. 10)  synchronous send          
!
! Author of original version:  P. Worley
! Ported to CAM: P. Worley, December 2003
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: iam
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
   integer, intent(in)   :: mtag           ! MPI message tag
   integer, intent(in)   :: swapnode       ! MPI process id of swap partner
   integer, intent(in)   :: rcvlth         ! length of incoming message buffer
   integer, intent(out)  :: rcvid          ! receive request id
   real(r8), intent(out) :: rcvmsg(rcvlth) ! incoming message buffer
!
!---------------------------Local workspace-----------------------------
!
   real(r8) signal      ! ready send signal
!
!-------------------------------------------------------------------------------------
!
   signal = 1.0
!
! simple swap: send/recv
   if ((swap_comm_order .eq. 0) .or. (swap_comm_order .eq. 2)) then
!
      if (swap_comm_protocol <= 5) then
!
         if (swap_comm_protocol <= 1) then
!
! this procotol does not use nonblocking receive.
!
         elseif (swap_comm_protocol <= 3) then
!
! post the receive before the send, increasing odds that the
! receive will be posted before the message arrives.
            call mpiirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                           swap_com, rcvid )
!
         else
!
! post the receive before send to allow use of ready send.
            call mpiirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                           swap_com, rcvid )
            call mpisend( signal, 1, mpir8, swapnode, mtag, handshake_com )
!
         endif
!
      elseif (swap_comm_protocol <= 9) then
!
         if (swap_comm_protocol <= 7) then
!
! these procotols do not use nonblocking receive.
!
         else
!
! post the receive before the synchronous send, increasing odds that the
! receive will be posted before the message arrives.
            call mpiirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com, rcvid )
!
         endif
!
      else
!
          write (0,901) swap_comm_order, swap_comm_protocol
  901     format(/,' fatal error in subroutine swap1:',   &
                 /,' unknown communication protocol specified',/, &
                   ' swap_comm_order = ',i6, ' swap_comm_protocol = ',i6)
          call endrun
!
      endif
!
   elseif (swap_comm_order .eq. 1) then
! ordered swap:
! if (iam <= swapnode) send/recv
! if (iam >= swapnode) recv/send
!
      if (swap_comm_protocol <= 5) then
!
         if (swap_comm_protocol <= 1) then
!
! this procotol does not use nonblocking receive.
!
         elseif (swap_comm_protocol <= 3) then
!
! post the receive before the initial send, increasing odds 
! that the receive will be posted before the message arrives.
            call mpiirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                           swap_com, rcvid )
!
         else
!
! post the receive before the send to allow use of forcetypes. 
            if (iam <= swapnode) then
               call mpiirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                              swap_com, rcvid )
            else
               call mpiirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                              swap_com, rcvid )
               call mpisend ( signal, 1, mpir8, swapnode, mtag, handshake_com )
            endif
!
         endif
!
      elseif (swap_comm_protocol <= 10) then
!
         if (swap_comm_protocol <= 7) then
!
! these protocols do not use nonblocking receive.
!
         elseif (swap_comm_protocol <= 9) then
!
! post the receive before the initial synchronous send, increasing odds 
! that the receive will be posted before the message arrives.
            call mpiirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com, rcvid )
!
         else
!
! this protocol does not use nonblocking receive.
!
         endif
!
      else
!
! protocol error
         write (0,901) swap_comm_order, swap_comm_protocol
         call endrun
!
      endif
!
   else
!***********************************************************************
!       undefined swap option
!***********************************************************************
!
       write (0,900) swap_comm_order
  900  format(/,' fatal error in subroutine swap1:', &
              /,' unknown communication option specified',/, &
                ' swap_comm_order = ',i6)                                 
       call endrun                                             
!
   endif
!
   return
   end subroutine swap1
!
!========================================================================
!
   subroutine swap2(mtag, swapnode, sndlth, sndmsg, sndid, &
                    rcvlth, rcvmsg, rcvid)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Second of three routines that implement swap using MPI point-to-point
! routines,
! 
! Method: 
! This subroutine continues the swap operation begun in swap1. It
! initiates the send and waits for the receive to complete. 
! if (swap_comm_order .eq. 0) simple swap: send/recv
! if (swap_comm_order .eq. 1) ordered swap: [send/recv]|[recv/send]
! if (swap_comm_order .eq. 2) delayed-recv swap: send ... recv
! if (swap_comm_protocol .eq. 1, 3, 5, 7, .or. 9) nonblocking send
! if (swap_comm_protocol .eq. 2, 3, 4, 5, 8, .or. 9)  nonblocking receive
! if (swap_comm_protocol .eq. 4 .or. 5)  ready send
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 0) sendrecv
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 1) explicitly synchronous  
! if (swap_comm_protocol .eq. 7, 8, 9, .or. 10)  synchronous send          
!
! Author of original version:  P. Worley
! Ported to CAM: P. Worley, December 2003
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: iam
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
   integer, intent(in)   :: mtag           ! MPI message tag
   integer, intent(in)   :: swapnode       ! MPI process id of swap partner
   integer, intent(in)   :: sndlth         ! length of outgoing message
   integer, intent(out)  :: sndid          ! send request id
   integer, intent(in)   :: rcvlth         ! length of incoming message buffer
   integer, intent(inout):: rcvid          ! receive request id
   real(r8), intent(in)  :: sndmsg(sndlth) ! outgoing message buffer
   real(r8), intent(out) :: rcvmsg(rcvlth) ! incoming message buffer
!
!---------------------------Local workspace-----------------------------
!
   real(r8) signal                         ! ready send signal
   integer  status(MPI_STATUS_SIZE)        ! MPI status integer
!
!-------------------------------------------------------------------------------------
!
   signal = 1.0
!
   if (swap_comm_order .eq. 0) then
!
      if (swap_comm_protocol <= 5) then
!
         if (swap_comm_protocol <= 1) then
!
! do not block for the send, enabling overlap of communication with computation
            call mpiisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com, sndid )
            call mpirecv ( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
!
         elseif (swap_comm_protocol <= 3) then
!
            if (swap_comm_protocol .eq. 2) then
! complete outstanding receive
               call mpisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
               call mpiwait( rcvid, status )
            else
! also do not block for the send, enabling overlap of communication with computation
               call mpiisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com, &
                              sndid )
               call mpiwait ( rcvid, status )
            endif
!
         else
!    
            if (swap_comm_protocol .eq. 4) then
! complete receive of ready send
               call mpirecv ( signal, 1, mpir8, swapnode, mtag, handshake_com )
               call mpirsend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
               call mpiwait ( rcvid, status )
            else
! also do not block for send, enabling overlap of communication with computation.
               call mpirecv  ( signal, 1, mpir8, swapnode, mtag, handshake_com )
               call mpiirsend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                               swap_com, sndid )
               call mpiwait  ( rcvid, status )
            endif
!
         endif
!
      elseif (swap_comm_protocol <= 9) then
!
         if (swap_comm_protocol <= 7) then
!
            if (swap_comm_protocol .eq. 6) then
! native sendrecv
               call mpisendrecv( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                 rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                                 swap_com )
            else
! do not block for the synchronous send, enabling overlap of 
! communication with computation.
               call mpiissend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com, &
                               sndid )
               call mpirecv  ( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
            endif
!
         else
!
            if (swap_comm_protocol .eq. 8) then
! complete outstanding receive,
               call mpissend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
               call mpiwait ( rcvid, status )
            else
! also do not block for the synchronous send, enabling overlap of
! communication with computation
               call mpiissend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com, &
                               sndid )
               call mpiwait( rcvid, status )
            endif
!
         endif
!
      else
!
         write (0,901) swap_comm_order, swap_comm_protocol
  901    format(/,' fatal error in subroutine swap2:', &
                /,' unknown communication protocol specified',/,  &
                  ' swap_comm_order = ',i6, ' swap_comm_protocol = ',i6)
         call endrun
!
      endif
!
   elseif (swap_comm_order .eq. 1) then
! ordered swap:
! if (iam <= swapnode) send/recv
! if (iam >= swapnode) recv/send
!
      if (swap_comm_protocol <= 5) then
!
         if (swap_comm_protocol <= 1) then
!
            if (swap_comm_protocol .eq. 0) then
!
               if (iam <= swapnode) then
                  call mpisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
                  call mpirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
               else
                  call mpirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
                  call mpisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
               endif
!
            else
!
! do not block for the send, enabling overlap of communication with computation.
               if (iam <= swapnode) then
                  call mpiisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com, &
                                 sndid )
                  call mpirecv ( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
               else
                  call mpirecv ( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
                  call mpiisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com, &
                                 sndid )
               endif
!
            endif
!
         elseif (swap_comm_protocol <= 3) then
!
            if (swap_comm_protocol .eq. 2) then
!
! complete outstanding receive.
               if (iam <= swapnode) then
                  call mpisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
                  call mpiwait( rcvid, status )
               else
                  call mpiwait( rcvid, status )
                  call mpisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
               endif
!
            else
!
! also do not block for the send, enabling overlap of communication with computation.
               if (iam <= swapnode) then
                  call mpiisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com, &
                                  sndid )
                  call mpiwait ( rcvid, status )
               else
                  call mpiwait ( rcvid, status )
                  call mpiisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com, &
                                 sndid )
               endif
 
            endif
!
         else
!
            if (swap_comm_protocol .eq. 4) then
!
! complete forcetype receive.
               if (iam <= swapnode) then
                  call mpirecv ( signal, 1, mpir8, swapnode, mtag, handshake_com )
                  call mpirsend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                 swap_com )
                  call mpiwait ( rcvid, status )
               else
                  call mpiwait ( rcvid, status )
                  call mpirsend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                 swap_com )
               endif
!
            else
!
! also do not block for the send, enabling overlap of communication with computation.
               if (iam <= swapnode) then
                  call mpirecv  ( signal, 1, mpir8, swapnode, mtag, handshake_com )
                  call mpiirsend( sndmsg, sndlth, mpir8, swapnode, mtag, & 
                                  swap_com, sndid )
                  call mpiwait  ( rcvid, status )
               else
                  call mpiwait  ( rcvid, status )
                  call mpiirsend( sndmsg, sndlth, mpir8, swapnode, mtag, & 
                                  swap_com, sndid )
               endif
!
            endif
!
         endif
!
      elseif (swap_comm_protocol <= 10) then
!
         if (swap_comm_protocol <= 7) then
!
            if (swap_comm_protocol .eq. 6) then
!
! synchronous ordered swap 
               if (iam <= swapnode) then
                  call mpirecv( signal, 1, mpir8, swapnode, mtag, &
                                handshake_com )
                  call mpisend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                swap_com )
                  call mpirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                                swap_com )
               else
                  call mpisend( signal, 1, mpir8, swapnode, mtag, &
                                handshake_com )
                  call mpirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                                swap_com )
                  call mpisend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                swap_com )
               endif
!
            else
!
! do not block for the synchronous send, enabling overlap of communication
! with computation.
               if (iam <= swapnode) then
                  call mpiissend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                  swap_com, sndid)
                  call mpirecv  ( rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                                  swap_com )
               else
                  call mpirecv  ( rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                                  swap_com )
                  call mpiissend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                  swap_com, sndid )
               endif
!
            endif
!
         elseif (swap_comm_protocol <= 9) then
!
            if (swap_comm_protocol .eq. 8) then
!
! complete outstanding receive.
               if (iam <= swapnode) then
                  call mpissend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                 swap_com )
                  call mpiwait( rcvid, status )
               else
                  call mpiwait ( rcvid, status )
                  call mpissend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                 swap_com )
               endif
!
            else
!
! also do not block for the synchronous send, enabling overlap of communication
! with computation.
               if (iam <= swapnode) then
                  call mpiissend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                  swap_com, sndid)
                  call mpiwait  ( rcvid, status )
               else
                  call mpiwait  ( rcvid, status )
                  call mpiissend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                  swap_com, sndid)
               endif
!
            endif
!
         else
! ordered swap using synchronous sends
            if (iam <= swapnode) then
               call mpissend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
               call mpirecv ( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
            else
               call mpirecv ( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
               call mpissend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
            endif
!
         endif
!
      else
!
! protocol error
         write (0,901) swap_comm_order, swap_comm_protocol
         call endrun
!
      endif
!
   elseif (swap_comm_order .eq. 2) then
! delayed swap: send ... recv
!
      if (swap_comm_protocol <= 5) then
!
         if (swap_comm_protocol <= 1) then
!
! do not block for send, enabling overlap of communication with computation.
            call mpiisend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                           swap_com, sndid )
!
         elseif (swap_comm_protocol <= 3) then
!
            if (swap_comm_protocol .eq. 2) then
! swap send
               call mpisend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
            else
! do not block for send, enabling overlap of communication with computation.
               call mpiisend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                              swap_com, sndid )
            endif
!
         else
!    
            if (swap_comm_protocol .eq. 4) then
! use ready send
               call mpirecv ( signal, 1, mpir8, swapnode, mtag, handshake_com )
               call mpirsend( sndmsg, sndlth, mpir8, swapnode, mtag, swap_com )
            else
! do not block for send, enabling overlap of communication with computation.
               call mpirecv  ( signal, 1, mpir8, swapnode, mtag, &
                               handshake_com )
               call mpiirsend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                               swap_com, sndid )
            endif
!
         endif
!
      elseif (swap_comm_protocol <= 9) then
!
         if (swap_comm_protocol <= 7) then
!
            if (swap_comm_protocol .eq. 6) then
!
! native swap
               call mpisendrecv( sndmsg, sndlth, mpir8, swapnode, mtag, &
                                 rcvmsg, rcvlth, mpir8, swapnode, mtag, &
                                 swap_com )
            else
! do not block for send, enabling overlap of communication with computation.
               call mpiisend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                              swap_com, sndid )
            endif
!
         else
!
            if (swap_comm_protocol .eq. 8) then
! swap send using synchronous send
               call mpissend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                              swap_com )
            else
! do not block for synchronous send, enabling overlap of communication 
! with computation.
               call mpiissend( sndmsg, sndlth, mpir8, swapnode, mtag, &
                               swap_com, sndid )
            endif
!
         endif
!
      else
!
! protocol error
         write (0,901) swap_comm_order, swap_comm_protocol
         stop                                                   
!
      endif
!
   else
! undefined swap option
!
       write (0,900) swap_comm_order
  900  format(/,' fatal error in subroutine swap2:', &
              /,' unknown communication option specified',/, &
                ' swap_comm_order = ',i6)                                 
       call endrun                                            
!
   endif
!
   return
   end subroutine swap2
!
!========================================================================
!
   subroutine swap3(mtag, swapnode, sndid, rcvlth, rcvmsg, rcvid)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Third of three routines that implement swap using MPI point-to-point
! routines,
! 
! Method: 
! This subroutine completes the swap operation begun in swap1 and swap2.
! It waits until the send and receive request made in swap2 have
! completed.
! if (swap_comm_order .eq. 0) simple swap: send/recv
! if (swap_comm_order .eq. 1) ordered swap: [send/recv]|[recv/send]
! if (swap_comm_order .eq. 2) delayed-recv swap: send ... recv
! if (swap_comm_protocol .eq. 1, 3, 5, 7, .or. 9) nonblocking send
! if (swap_comm_protocol .eq. 2, 3, 4, 5, 8, .or. 9)  nonblocking receive
! if (swap_comm_protocol .eq. 4 .or. 5)  ready send
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 0) sendrecv
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 1) explicitly synchronous  
! if (swap_comm_protocol .eq. 7, 8, 9, .or. 10)  synchronous send          
!
! Author of original version:  P. Worley
! Ported to CAM: P. Worley, December 2003
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: iam
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
   integer, intent(in)   :: mtag           ! MPI message tag
   integer, intent(in)   :: swapnode       ! MPI process id of swap partner
   integer, intent(inout):: sndid          ! send request id
   integer, intent(in)   :: rcvlth         ! length of incoming message buffer
   integer, intent(inout):: rcvid          ! receive request id
   real(r8), intent(out) :: rcvmsg(rcvlth) ! incoming message buffer
!
!---------------------------Local workspace-----------------------------
!
   real(r8) signal                         ! ready send signal
   integer  status(MPI_STATUS_SIZE)        ! MPI status integer
!
!-------------------------------------------------------------------------------------
!
   signal = 1.0
!
   if (swap_comm_order .eq. 0) then
! simple swap: send/recv
!
! complete send for nonblocking send protocols.
      if ((swap_comm_protocol .eq. 1) .or. (swap_comm_protocol .eq. 3) .or. &
          (swap_comm_protocol .eq. 5) .or. (swap_comm_protocol .eq. 7) .or. &
          (swap_comm_protocol .eq. 9)) then
         call mpiwait( sndid, status )
      elseif (swap_comm_protocol .gt. 9) then
         write (0,901) swap_comm_order, swap_comm_protocol
  901    format(/,' fatal error in subroutine swap3:', &
                /,' unknown communication protocol specified',/,  &
                  ' swap_comm_order = ',i6, ' swap_comm_protocol = ',i6)
         call endrun                                                   
      endif
!
   elseif (swap_comm_order .eq. 1) then
! ordered swap:
! if (iam <= swapnode) send/recv
! if (iam >= swapnode) recv/send
!
! complete send for nonblocking send protocols.
      if ((swap_comm_protocol .eq. 1) .or. (swap_comm_protocol .eq. 3) .or. &
          (swap_comm_protocol .eq. 5) .or. (swap_comm_protocol .eq. 7) .or. &
          (swap_comm_protocol .eq. 9)) then
         call mpiwait( sndid, status )
      elseif (swap_comm_protocol .gt. 10) then
         write (0,901) swap_comm_order, swap_comm_protocol
         call endrun
      endif
!
   elseif (swap_comm_order .eq. 2) then
! delayed-recv swap: recvbegin ... send ... recvend
!
      if (swap_comm_protocol <= 5) then
!
         if (swap_comm_protocol <= 1) then
!
            if (swap_comm_protocol .eq. 0) then
! receive message.
               call mpirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
            else
! also complete send.
               call mpirecv(rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
               call mpiwait( sndid, status )
            endif
!
         elseif (swap_comm_protocol <= 3) then
!
            if (swap_comm_protocol .eq. 2) then
! complete receive.
               call mpiwait( rcvid, status )
            else
! also complete send.
               call mpiwait( rcvid, status )
               call mpiwait( sndid, status )
            endif
!
         else
!
            if (swap_comm_protocol .eq. 4) then
! complete receive.
               call mpiwait( rcvid, status )
            else
! also complete send.
               call mpiwait( rcvid, status )
               call mpiwait( sndid, status )
            endif
!
         endif
!
      elseif (swap_comm_protocol <= 9) then
!
         if (swap_comm_protocol <= 7) then
!
            if (swap_comm_protocol .eq. 6) then
! receive already complete in "native" swap
            else
! also complete send.
               call mpirecv( rcvmsg, rcvlth, mpir8, swapnode, mtag, swap_com )
               call mpiwait( sndid, status )
            endif
!
         else
!
            if (swap_comm_protocol .eq. 8) then
! complete receive.
               call mpiwait( rcvid, status )
            else
! also complete send.
               call mpiwait( rcvid, status )
               call mpiwait( sndid, status )
            endif
!
         endif
!
      else
!
! protocol error
         write (0,901) swap_comm_order, swap_comm_protocol
         call endrun
!
      endif
!
   else
! undefined swap option
      write (0,900) swap_comm_order
  900 format(/,' fatal error in subroutine swap3:', &
             /,' unknown communication option specified',/, &
               ' swap_comm_order = ',i6)                                 
      call endrun
!
   endif
!
   return
   end subroutine swap3
!
!========================================================================
!
   subroutine swap1m(cnt, mtag, swapnodes, &
                     rcvlths, rdispls, bufsiz, rcvbuf, rcvids)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! First of three routines that implement swap using MPI point-to-point
! routines. Variant of swap1 for multiple messages.
! 
! Method: 
! This subroutine begins a swap operation that will be completed by
! swap2 and swap3. It posts a receive and sends handshaking messages
! when ready sends are used. 
! if (swap_comm_order .eq. 0) simple swap: send/recv
! if (swap_comm_order .eq. 1) ordered swap: [send/recv]|[recv/send]
! if (swap_comm_order .eq. 2) delayed-recv swap: send ... recv
! if (swap_comm_protocol .eq. 1, 3, 5, 7, .or. 9) nonblocking send
! if (swap_comm_protocol .eq. 2, 3, 4, 5, 8, .or. 9)  nonblocking receive
! if (swap_comm_protocol .eq. 4 .or. 5)  ready send
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 0) sendrecv
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 1) explicitly synchronous  
! if (swap_comm_protocol .eq. 7, 8, 9, .or. 10)  synchronous send          
!
! Author of original version:  P. Worley
! Ported to CAM: P. Worley, December 2003
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: iam
   use spmd_dyn, only: npes
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
   integer, intent(in)   :: cnt            ! number of swaps to initiate
   integer, intent(in)   :: mtag           ! MPI message tag
   integer, intent(in)   :: swapnodes(cnt) ! MPI process id of swap partners
   integer, intent(in)   :: rcvlths(0:npes-1)
                                           ! length of incoming messages
   integer, intent(in)   :: rdispls(0:npes-1) 
                                           ! offset from beginning of receive 
                                           !  buffer where incoming messages
                                           !  should be placed
   integer, intent(in)   :: bufsiz         ! message buffer size
   integer, intent(inout) :: rcvids(cnt)   ! receive request ids
   real(r8), intent(out) :: rcvbuf(bufsiz) ! buffer for incoming messages
!
!---------------------------Local workspace-----------------------------
!
   integer  i                              ! loop index
   integer  p                              ! process index
   integer  offset_r                       ! index of message beginning in 
                                           !  receive buffer
   real(r8) signal                         ! ready send signal
!
!-------------------------------------------------------------------------------------
!
   signal = 1.0
!
! simple swap: send/recv
   if ((swap_comm_order .eq. 0) .or. (swap_comm_order .eq. 2)) then
!
      if (swap_comm_protocol <= 5) then
!
         if (swap_comm_protocol <= 1) then
!
! this procotol does not use nonblocking receive.
!
         elseif (swap_comm_protocol <= 3) then
!
! post the receive before the send, increasing odds that the
! receive will be posted before the message arrives.
            do i=1,cnt
               p = swapnodes(i)
               offset_r = rdispls(p)+1
               call mpiirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                              swap_com, rcvids(i) )
            enddo
!
         else
!
! post the receive before send to allow use of ready send.
            do i=1,cnt
               p = swapnodes(i)
               offset_r = rdispls(p)+1
               call mpiirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                              swap_com, rcvids(i) )
               call mpisend( signal, 1, mpir8, p, mtag, handshake_com )
            enddo
!
         endif
!
      elseif (swap_comm_protocol <= 9) then
!
         if (swap_comm_protocol <= 7) then
!
! these procotols do not use nonblocking receive.
!
         else
!
! post the receive before the synchronous send, increasing odds that the
! receive will be posted before the message arrives.
            do i=1,cnt
               p = swapnodes(i)
               offset_r = rdispls(p)+1
               call mpiirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                              swap_com, rcvids(i) )
            enddo
!
         endif
!
      else
!
          write (0,901) swap_comm_order, swap_comm_protocol
  901     format(/,' fatal error in subroutine swap1m:',   &
                 /,' unknown communication protocol specified',/, &
                   ' swap_comm_order = ',i6, ' swap_comm_protocol = ',i6)
          call endrun
!
      endif
!
   elseif (swap_comm_order .eq. 1) then
! ordered swap:
! if (iam <= swapnode) send/recv
! if (iam >= swapnode) recv/send
!
      if (swap_comm_protocol <= 5) then
!
         if (swap_comm_protocol <= 1) then
!
! this procotol does not use nonblocking receive.
!
         elseif (swap_comm_protocol <= 3) then
!
! post the receive before the initial send, increasing odds 
! that the receive will be posted before the message arrives.
            do i=1,cnt
               p = swapnodes(i)
               offset_r = rdispls(p)+1
               call mpiirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                              swap_com, rcvids(i) )
            enddo
!
         else
!
! post the receive before the send to allow use of forcetypes. 
            do i=1,cnt
               p = swapnodes(i)
               offset_r = rdispls(p)+1
               if (iam <= swapnodes(p)) then
                  call mpiirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                 swap_com, rcvids(i) )
               else
                  call mpiirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                 swap_com, rcvids(i) )
                  call mpisend( signal, 1, mpir8, p, mtag, handshake_com )
               endif
            enddo
!
         endif
!
      elseif (swap_comm_protocol <= 10) then
!
         if (swap_comm_protocol <= 7) then
!
! these protocols do not use nonblocking receive.
!
         elseif (swap_comm_protocol <= 9) then
!
! post the receive before the initial synchronous send, increasing odds 
! that the receive will be posted before the message arrives.
            do i=1,cnt
               p = swapnodes(i)
               offset_r = rdispls(p)+1
               call mpiirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                              swap_com, rcvids(i) )
            enddo
!
         else
!
! this protocol does not use nonblocking receive.
!
         endif
!
      else
!
! protocol error
         write (0,901) swap_comm_order, swap_comm_protocol
         call endrun
!
      endif
!
   else
!***********************************************************************
!       undefined swap option
!***********************************************************************
!
       write (0,900) swap_comm_order
  900  format(/,' fatal error in subroutine swap1m:', &
              /,' unknown communication option specified',/, &
                ' swap_comm_order = ',i6)                                 
       call endrun                                             
!
   endif
!
   return
   end subroutine swap1m
!
!========================================================================
!
   subroutine swap3m(cnt, mtag, swapnodes, sndids, &
                     rcvlths, rdispls, bufsiz, rcvbuf, rcvids)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Third of three routines that implement swap using MPI point-to-point
! routines. Variant of swap3 for multiple messages.
! 
! Method: 
! This subroutine completes swap operations begun in swap1 and swap2.
! It waits until the send and receive request made in swap2 have
! completed.
! if (swap_comm_order .eq. 0) simple swap: send/recv
! if (swap_comm_order .eq. 1) ordered swap: [send/recv]|[recv/send]
! if (swap_comm_order .eq. 2) delayed-recv swap: send ... recv
! if (swap_comm_protocol .eq. 1, 3, 5, 7, .or. 9) nonblocking send
! if (swap_comm_protocol .eq. 2, 3, 4, 5, 8, .or. 9)  nonblocking receive
! if (swap_comm_protocol .eq. 4 .or. 5)  ready send
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 0) sendrecv
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 1) explicitly synchronous  
! if (swap_comm_protocol .eq. 7, 8, 9, .or. 10)  synchronous send          
!
! Author of original version:  P. Worley
! Ported to CAM: P. Worley, December 2003
! 
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid, only: iam
   use spmd_dyn, only: npes
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
   integer, intent(in) :: cnt             ! number of swaps to complete
   integer, intent(in) :: mtag            ! MPI message tag
   integer, intent(in) :: swapnodes(cnt)  ! MPI process id of swap partners
   integer, intent(inout) :: sndids(cnt)  ! send request ids
   integer, intent(in) :: rcvlths(0:npes-1)    
                                          ! length of incoming messages
   integer, intent(in) :: rdispls(0:npes-1) 
                                          ! offset from beginning of receive 
                                          !  buffer where incoming messages
                                          !  should be placed
   integer, intent(inout) :: rcvids(cnt)  ! receive request ids
   integer, intent(in) :: bufsiz          ! message buffer size
   real(r8), intent(out) :: rcvbuf(bufsiz) 
                                          ! buffer for incoming messages
!
!---------------------------Local workspace-----------------------------
!
   integer  i                              ! loop index
   integer  p                              ! process index
   integer  offset_r                       ! index of message beginning in 
                                           !  receive buffer
   integer  status(MPI_STATUS_SIZE,cnt)    ! MPI status integers
   real(r8) signal                         ! ready send signal
!
!-------------------------------------------------------------------------------------
!
   signal = 1.0
!
   if (swap_comm_order .eq. 0) then
! simple swap: send/recv
!
! complete send for nonblocking send protocols.
      if ((swap_comm_protocol .eq. 1) .or. (swap_comm_protocol .eq. 3) .or. &
          (swap_comm_protocol .eq. 5) .or. (swap_comm_protocol .eq. 7) .or. &
          (swap_comm_protocol .eq. 9)) then
         call mpiwaitall ( cnt, sndids, status )
      elseif (swap_comm_protocol .gt. 9) then
         write (0,901) swap_comm_order, swap_comm_protocol
  901    format(/,' fatal error in subroutine swap3m:', &
                /,' unknown communication protocol specified',/,  &
                  ' swap_comm_order = ',i6, ' swap_comm_protocol = ',i6)
         call endrun                                                   
      endif
!
   elseif (swap_comm_order .eq. 1) then
! ordered swap:
! if (iam <= swapnode) send/recv
! if (iam >= swapnode) recv/send
!
! complete send for nonblocking send protocols.
      if ((swap_comm_protocol .eq. 1) .or. (swap_comm_protocol .eq. 3) .or. &
          (swap_comm_protocol .eq. 5) .or. (swap_comm_protocol .eq. 7) .or. &
          (swap_comm_protocol .eq. 9)) then
         call mpiwaitall ( cnt, sndids, status )
      elseif (swap_comm_protocol .gt. 10) then
         write (0,901) swap_comm_order, swap_comm_protocol
         call endrun
      endif
!
   elseif (swap_comm_order .eq. 2) then
! delayed-recv swap: recvbegin ... send ... recvend
!
      if (swap_comm_protocol <= 5) then
!
         if (swap_comm_protocol <= 1) then
!
            if (swap_comm_protocol .eq. 0) then
! receive message.
               do i=1,cnt
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  call mpirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag,&
                                swap_com )
               enddo
            else
! also complete send.
               do i=1,cnt
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  call mpirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                swap_com )
               enddo
               call mpiwaitall( cnt, sndids, status )
            endif
!
         elseif (swap_comm_protocol <= 3) then
!
            if (swap_comm_protocol .eq. 2) then
! complete receive.
               call mpiwaitall( cnt, rcvids, status )
            else
! also complete send.
               call mpiwaitall( cnt, rcvids, status )
               call mpiwaitall( cnt, sndids, status )
            endif
!
         else
!
            if (swap_comm_protocol .eq. 4) then
! complete receive.
               call mpiwaitall( cnt, rcvids, status )
            else
! also complete send.
               call mpiwaitall( cnt, rcvids, status )
               call mpiwaitall( cnt, sndids, status )
            endif
!
         endif
!
      elseif (swap_comm_protocol <= 9) then
!
         if (swap_comm_protocol <= 7) then
!
            if (swap_comm_protocol .eq. 6) then
! receive already complete in "native" swap
            else
! also complete send.
               do i=1,cnt
                  p = swapnodes(i)
                  offset_r = rdispls(p)+1
                  call mpirecv( rcvbuf(offset_r), rcvlths(p), mpir8, p, mtag, &
                                swap_com )
               enddo
               call mpiwaitall( cnt, sndids, status )
            endif
!
         else
!
            if (swap_comm_protocol .eq. 8) then
! complete receive.
               call mpiwaitall( cnt, rcvids, status )
            else
! also complete send.
               call mpiwaitall( cnt, rcvids, status )
               call mpiwaitall( cnt, sndids, status )
            endif
!
         endif
!
      else
!
! protocol error
         write (0,901) swap_comm_order, swap_comm_protocol
         call endrun
!
      endif
!
   else
! undefined swap option
      write (0,900) swap_comm_order
  900 format(/,' fatal error in subroutine swap3m:', &
             /,' unknown communication option specified',/, &
               ' swap_comm_order = ',i6)                                 
      call endrun
!
   endif
!
   return
   end subroutine swap3m
!
!========================================================================
!
   logical function do_swap1()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Indicates whether swap1 does anything when called with the current
! communication option and communication protocol.
! if (swap_comm_order .eq. 0) simple swap: send/recv
! if (swap_comm_order .eq. 1) ordered swap: [send/recv]|[recv/send]
! if (swap_comm_order .eq. 2) delayed-recv swap: send ... recv
! if (swap_comm_protocol .eq. 1, 3, 5, 7, .or. 9) nonblocking send
! if (swap_comm_protocol .eq. 2, 3, 4, 5, 8, .or. 9)  nonblocking receive
! if (swap_comm_protocol .eq. 4 .or. 5)  ready send
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 0) sendrecv
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 1) explicitly synchronous  
! if (swap_comm_protocol .eq. 7, 8, 9, .or. 10)  synchronous send          
! 
! Method: 
!
! Author: P. Worley
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
!-------------------------------------------------------------------------------------
!
   do_swap1 = .false.
   if ((swap_comm_order >= 0) .and. (swap_comm_order <= 2)) then
      if (((swap_comm_protocol >= 2) .and. (swap_comm_protocol <= 5)) .or. &
         (((swap_comm_protocol >= 8) .and. (swap_comm_protocol <= 9)))) then
        do_swap1 = .true.
      endif
   endif
!
   return
   end function do_swap1
!
!========================================================================
!
   logical function do_swap3()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Indicates whether swap3 does anything when called with the current
! communication option and communication protocol.
! if (swap_comm_order .eq. 0) simple swap: send/recv
! if (swap_comm_order .eq. 1) ordered swap: [send/recv]|[recv/send]
! if (swap_comm_order .eq. 2) delayed-recv swap: send ... recv
! if (swap_comm_protocol .eq. 1, 3, 5, 7, .or. 9) nonblocking send
! if (swap_comm_protocol .eq. 2, 3, 4, 5, 8, .or. 9)  nonblocking receive
! if (swap_comm_protocol .eq. 4 .or. 5)  ready send
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 0) sendrecv
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 1) explicitly synchronous  
! if (swap_comm_protocol .eq. 7, 8, 9, .or. 10)  synchronous send          
! 
! Method: 
!
! Author: P. Worley
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
!-------------------------------------------------------------------------------------
!
   do_swap3 = .false.
   if ((swap_comm_order >= 0) .and. (swap_comm_order <= 1)) then
      if ((swap_comm_protocol >= 1) .and. (swap_comm_protocol <= 9)) then
         if (mod(swap_comm_protocol,2) .eq. 1) then
            do_swap3 = .true.
         endif
      endif
   elseif (swap_comm_order .eq. 2) then
      if (((swap_comm_protocol >= 0) .and. (swap_comm_protocol <= 5)) .or. &
         (((swap_comm_protocol >= 7) .and. (swap_comm_protocol <= 9)))) then
        do_swap3 = .true.
      endif
   endif
!
   return
   end function do_swap3

!
!========================================================================
!
   logical function delayed_swap_recv()

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Indicates whether message receives occur in swap3 with the current
! communication option and communication protocol.
! if (swap_comm_order .eq. 0) simple swap: send/recv
! if (swap_comm_order .eq. 1) ordered swap: [send/recv]|[recv/send]
! if (swap_comm_order .eq. 2) delayed-recv swap: send ... recv
! if (swap_comm_protocol .eq. 1, 3, 5, 7, .or. 9) nonblocking send
! if (swap_comm_protocol .eq. 2, 3, 4, 5, 8, .or. 9)  nonblocking receive
! if (swap_comm_protocol .eq. 4 .or. 5)  ready send
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 0) sendrecv
! if (swap_comm_protocol .eq. 6 .and. swap_comm_order .eq. 1) explicitly synchronous  
! if (swap_comm_protocol .eq. 7, 8, 9, .or. 10)  synchronous send          
! 
! Method: 
!
! Author: P. Worley
! 
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   implicit none
!---------------------------Input arguments--------------------------
!
!-------------------------------------------------------------------------------------
!
   delayed_swap_recv = .false.
   if (swap_comm_order .eq. 2) then
      if (((swap_comm_protocol >= 0) .and. (swap_comm_protocol <= 5)) .or. &
         (((swap_comm_protocol >= 7) .and. (swap_comm_protocol <= 9)))) then
        delayed_swap_recv = .true.
      endif
   endif
!
   return
   end function delayed_swap_recv

#endif

end module swap_comm
