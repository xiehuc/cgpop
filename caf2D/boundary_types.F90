!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! This module contains data types for updating ghost cell boundaries 
!<
module boundary_types
    ! !USES:
    use kinds_mod, only: i4, r4, r8, log_kind
    use simple_type, only: distrb
    use communicate, only: mpitag_bndy_2d, my_task, MPI_DBL
    use constants, only: field_loc_necorner, field_loc_center, &
        field_loc_eface, field_loc_nface, field_type_vector, &
        field_type_angle, field_type_scalar, p5, c0, ALG_MPI2S_2D, &
        ALG_MPI2S_1D, ALG_CAF_SINGLE_PUSH_2D, ALG_CAF_SINGLE_PULL_2D, &
        boundary_exchange_algorithm
    use simple_blocks, only: nx_block, ny_block, nghost, nblocks_x, nblocks_y
    use exit_mod, only: exit_POP, sigAbort

    implicit none
    private
    save

    ! !PUBLIC TYPES:
    type, public :: bndy
        integer (i4) :: &
            communicator       ,&! communicator to use for update messages
            nmsg_ew_snd        ,&! number of messages to send for e-w update
            nmsg_ns_snd        ,&! number of messages to send for n-s update
            nmsg_ew_rcv        ,&! number of messages to recv for e-w update
            nmsg_ns_rcv        ,&! number of messages to recv for n-s update
            maxblocks_ew_snd   ,&! max num blocks involved in east-west sends
            maxblocks_ew_rcv   ,&! max num blocks involved in east-west recvs
            maxblocks_ns_snd   ,&! max num blocks involved in north-south sends
            maxblocks_ns_rcv   ,&! max num blocks involved in north-south recvs
            nlocal_ew          ,&! num local copies for east-west bndy update
            nlocal_ns            ! num local copies for east-west bndy update

         integer (i4), dimension(:), pointer :: &
           nblocks_ew_snd     ,&! num blocks in each east-west send msg
           nblocks_ns_snd     ,&! num blocks in each north-south send msg
           nblocks_ew_rcv     ,&! num blocks in each east-west recv msg
           nblocks_ns_rcv     ,&! num blocks in each north-south recv msg
           ew_snd_proc        ,&! dest   proc for east-west send message
           ew_rcv_proc        ,&! source proc for east-west recv message
           ns_snd_proc        ,&! dest   proc for north-south send message
           ns_rcv_proc        ,&! source proc for north-south recv message
           local_ew_src_block ,&! source block for each local east-west copy
           local_ew_dst_block ,&! dest   block for each local east-west copy
           local_ns_src_block ,&! source block for each local north-south copy
           local_ns_dst_block   ! dest   block for each local north-south copy

         integer (i4), dimension(:,:), pointer :: &
           local_ew_src_add   ,&! starting source address for local e-w copies
           local_ew_dst_add   ,&! starting dest   address for local e-w copies
           local_ns_src_add   ,&! starting source address for local n-s copies
           local_ns_dst_add     ! starting dest   address for local n-s copies

         integer (i4), dimension(:,:), pointer :: &
           ew_src_block       ,&! source block for sending   e-w bndy msg
           ew_dst_block       ,&! dest   block for receiving e-w bndy msg
           ns_src_block       ,&! source block for sending   n-s bndy msg
           ns_dst_block         ! dest   block for sending   n-s bndy msg

         integer (i4), dimension(:,:,:), pointer :: &
           ew_src_add         ,&! starting source address for e-w msgs
           ew_dst_add         ,&! starting dest   address for e-w msgs
           ns_src_add         ,&! starting source address for n-s msgs
           ns_dst_add           ! starting dest   address for n-s msgs
    end type bndy

    type (bndy), public :: &!  ghost cell update info
        bndy_tropic          ! block distribution for barotropic part
end module boundary_types
