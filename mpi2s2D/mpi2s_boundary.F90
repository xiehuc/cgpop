!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! This module contains data types and routines for updating ghost cell
!! boundaries using MPI calls
!<
module mpi2s_boundary
    ! !USES:
    use kinds_mod, only : i4, r4, r8, log_kind
    use simple_type, only : distrb
    use communicate, only : mpitag_bndy_2d, my_task, MPI_DBL
    use constants, only : field_loc_necorner, field_loc_center, &
        field_loc_eface, field_loc_nface, field_type_vector, &
        field_type_angle, field_type_scalar, p5, c0
    use simple_blocks
    use exit_mod, only : exit_POP, sigAbort

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

    integer(i4), public :: timer_mpi2s_boundary_create, &
        timer_mpi2s_boundary_2d_dbl

    ! !PUBLIC MEMBER FUNCTIONS:
    public :: mpi2s_create_boundary,  &
        mpi2s_destroy_boundary, &
        mpi2s_boundary_2d

    interface mpi2s_boundary_2d  ! generic interface
        module procedure mpi2s_boundary_2d_dbl
        module procedure mpi2s_boundary_2d_int
    end interface

    type (bndy), public :: &!  ghost cell update info
        bndy_tropic          ! block distribution for barotropic part

    !maltrud debug
    integer (i4) :: index_check

  contains

    !***********************************************************************
    !>
    !! This routine creates a boundary type with info necessary for
    !! performing a boundary (ghost cell) update based on the input block
    !! distribution.
    !!
    !! @param newbndy       a new boundary type with info for updates.
    !! @param dist          A data structure which describes the distribution
    !!                      of blocks to MPI tasks.
    !! @param blocks        A data structure for the application specific
    !!                      metadata.
    !! @param ns_bndy_type  type of boundary to use in ns dir
    !! @param ew_bndy_type  type of boundary to use in ew dir
    !! @param nx_global     global X extent of domain
    !! @param ny_global     global Y extent of domain
    !<
    subroutine mpi2s_create_boundary(newbndy, dist, blocks, &
        ns_bndy_type, ew_bndy_type, nx_global, ny_global)

        ! !INPUT PARAMETERS:
        type (distrb), intent(in) :: dist       

        integer (i4) :: blocks(:)

        character (*), intent(in) :: &
            ns_bndy_type,             &
            ew_bndy_type

        integer (i4), intent(in) :: &
            nx_global, ny_global

        ! !OUTPUT PARAMETERS:
        type (bndy), intent(out) :: newbndy    

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer (i4) ::           &
            i,j,k,kk,n,                  &! dummy counters
            iblock_src  , jblock_src  ,  &! i,j index of source block
            iblock_dst  , jblock_dst  ,  &! i,j index of dest   block
            iblock_north, jblock_north,  &! i,j index of north neighbor block
            iblock_south, jblock_south,  &! i,j index of south neighbor block
            iblock_east , jblock_east ,  &! i,j index of east  neighbor block
            iblock_west , jblock_west ,  &! i,j index of west  neighbor block
            src_block_loc,               &! local block location of source
            dst_block_loc,               &! local block location of dest
            imsg_ew_snd, imsg_ew_rcv,    &! counters for ew send/recv
            imsg_ns_snd, imsg_ns_rcv,    &! counters for ns send/recv
            nprocs,                      &! num of processors involved
            nblocks,                     &! total number of blocks
            bid, pid,                    &! block and processor locators
            iblk, imsg,                  &!
            iloc_ew, iloc_ns,            &!
            src_proc, dst_proc            ! src,dst processor for message


        integer (i4), dimension(:), allocatable :: &
            ew_snd_count,    &! array for counting blocks in each message
            ew_rcv_count,    &! array for counting blocks in each message
            ns_snd_count,    &! array for counting blocks in each message
            ns_rcv_count,    &! array for counting blocks in each message
            msg_ew_snd  ,    &! msg counter for each active processor
            msg_ew_rcv  ,    &! msg counter for each active processor
            msg_ns_snd  ,    &! msg counter for each active processor
            msg_ns_rcv        ! msg counter for each active processor

        integer(i4), parameter :: east  = 1, &! index of east neighbor
            west  = 2, &! index of west neighbor
            north = 3, &! index of north neighbor
            south = 4   ! index of south neighbor

        integer(i4), parameter :: NumNeigh = 4  ! Number of neighbors
        integer(i4) :: sbie,sbje,dbie,dbje

        !-------------------------------------------------------------------
        !
        !  Initialize some useful variables and return if this task not
        !  in the current distribution.
        !
        !-------------------------------------------------------------------
        nprocs = dist%nprocs

        if (my_task >= nprocs) return

        nblocks = size(dist%proc(:))
        newbndy%communicator = dist%communicator

        !-------------------------------------------------------------------
        !
        !  Count the number of messages to send/recv from each processor
        !  and number of blocks in each message.  These quantities are
        !  necessary for allocating future arrays.
        !
        !-------------------------------------------------------------------

        allocate (ew_snd_count(nprocs), ew_rcv_count(nprocs), &
        ns_snd_count(nprocs), ns_rcv_count(nprocs))

        ew_snd_count = 0
        ew_rcv_count = 0
        ns_snd_count = 0
        ns_rcv_count = 0

        block_loop0: do n=1,nblocks
            src_proc  = dist%proc(n)

            call get_block_parameter(n,iblock=iblock_src,jblock=jblock_src)

            !*** compute cartesian i,j block indices for each neighbor
            !*** use zero if off the end of closed boundary
            !*** use jnorth=nblocks_y and inorth < 0 for tripole boundary
            !***   to make sure top boundary communicated to all top
            !***   boundary blocks

            call GetEWBlockIndex(ew_bndy_type,iblock_src,jblock_src, &
                iblock_east,jblock_east,iblock_west,jblock_west)

            call GetNSBlockIndex(ns_bndy_type,iblock_src, &
                jblock_src,iblock_north,jblock_north,iblock_south,jblock_south)
        end do block_loop0

        block_loop1: do n=1,nblocks
            src_proc  = dist%proc(n)

            call get_block_parameter(n,iblock=iblock_src,jblock=jblock_src)

            !*** compute cartesian i,j block indices for each neighbor
            !*** use zero if off the end of closed boundary
            !*** use jnorth=nblocks_y and inorth < 0 for tripole boundary
            !***   to make sure top boundary communicated to all top
            !***   boundary blocks
            call GetEWBlockIndex(ew_bndy_type,iblock_src,jblock_src, &
            iblock_east,jblock_east,iblock_west,jblock_west)

            call GetNSBlockIndex(ns_bndy_type,iblock_src, &
            jblock_src,iblock_north,jblock_north,iblock_south,jblock_south)

            !***
            !*** if any boundary is closed boundary, create a local
            !*** copy pseudo-message to fill ghost cells
            !***
            if (iblock_east == 0) &
                call increment_message_counter(ew_snd_count, ew_rcv_count, &
                0, src_proc)
            if (iblock_west == 0) &
                call increment_message_counter(ew_snd_count, ew_rcv_count, &
                0, src_proc)
            if (jblock_north == 0) &
                call increment_message_counter(ns_snd_count, ns_rcv_count, &
                0, src_proc)
            if (jblock_south == 0) &
                call increment_message_counter(ns_snd_count, ns_rcv_count, &
                0, src_proc)

            !***
            !*** now look through all the blocks for the neighbors
            !*** of the source block and check whether a message is
            !*** required for communicating with the neighbor
            !***
            do kk=1,NumNeigh
                k = Neigh(kk,n) 
                if(k>0) then 
                    call get_block_parameter( &
                        k,iblock=iblock_dst,jblock=jblock_dst)

                    dst_proc = dist%proc(k)  ! processor that holds dst block

                    !***
                    !*** if this block is an eastern neighbor
                    !*** increment message counter
                    !***
                    if (iblock_dst == iblock_east .and. &
                        jblock_dst == jblock_east) then

                        call increment_message_counter( &
                            ew_snd_count, ew_rcv_count, src_proc, dst_proc)
                    endif

                    !***
                    !*** if this block is an western neighbor
                    !*** increment message counter
                    !***

                    if (iblock_dst == iblock_west .and. &
                        jblock_dst == jblock_west) then

                        call increment_message_counter( &
                            ew_snd_count, ew_rcv_count, src_proc, dst_proc)
                    endif

                    !***
                    !*** if this block is an northern neighbor
                    !*** find out whether a message is required
                    !*** for tripole, must communicate with all
                    !*** north row blocks (triggered by iblock_dst <0)
                    !***

                    if ((iblock_dst == iblock_north .or. iblock_north < 0) &
                        .and. jblock_dst == jblock_north) then

                            call increment_message_counter(ns_snd_count, &
                                ns_rcv_count, src_proc, dst_proc)
                    endif

                    !***
                    !*** if this block is an southern neighbor
                    !*** find out whether a message is required
                    !***
                    if (iblock_dst == iblock_south .and. &
                        jblock_dst == jblock_south) then

                        call increment_message_counter(ns_snd_count, &
                            ns_rcv_count, src_proc, dst_proc)
                    endif
                endif ! k>0
            end do  ! loop over neighbors 
        end do block_loop1

        !*** if messages are received from the same processor
        !*** the message is actually a local copy - count them
        !*** and reset to zero
        newbndy%nlocal_ew = ew_rcv_count(my_task+1)
        newbndy%nlocal_ns = ns_rcv_count(my_task+1)
        ew_snd_count(my_task+1) = 0
        ew_rcv_count(my_task+1) = 0
        ns_snd_count(my_task+1) = 0
        ns_rcv_count(my_task+1) = 0

        !*** now count the number of actual messages to be
        !*** sent and received
        newbndy%nmsg_ew_snd = count(ew_snd_count /= 0)
        newbndy%nmsg_ns_snd = count(ns_snd_count /= 0)
        newbndy%nmsg_ew_rcv = count(ew_rcv_count /= 0)
        newbndy%nmsg_ns_rcv = count(ns_rcv_count /= 0)

        !*** find the maximum number of blocks sent in any one
        !*** message to use as an array size parameter
        newbndy%maxblocks_ew_snd = maxval(ew_snd_count)
        newbndy%maxblocks_ew_rcv = maxval(ew_rcv_count)
        newbndy%maxblocks_ns_snd = maxval(ns_snd_count)
        newbndy%maxblocks_ns_rcv = maxval(ns_rcv_count)

        !***
        !*** create buffers for tracking which message
        !*** is sent/received from which processor
        !***
        allocate(msg_ew_snd(nprocs), msg_ew_rcv(nprocs), &
            msg_ns_snd(nprocs), msg_ns_rcv(nprocs))

        msg_ew_snd = 0
        msg_ew_rcv = 0
        msg_ns_snd = 0
        msg_ns_rcv = 0

        !***
        !*** assign a location in buffer for each message to a
        !*** different processor. scramble the processor order
        !*** using current task id as offset to prevent all
        !*** processors sending to the same processor at the
        !*** same time
        !***
        imsg_ew_snd = 0
        imsg_ew_rcv = 0
        imsg_ns_snd = 0
        imsg_ns_rcv = 0

        do n=1,nprocs
            dst_proc = modulo(my_task+n,nprocs) + 1
            if (ew_snd_count(dst_proc) /= 0) then
                imsg_ew_snd = imsg_ew_snd + 1
                msg_ew_snd(dst_proc) = imsg_ew_snd
            endif
            if (ew_rcv_count(dst_proc) /= 0) then
                imsg_ew_rcv = imsg_ew_rcv + 1
                msg_ew_rcv(dst_proc) = imsg_ew_rcv
            endif
            if (ns_snd_count(dst_proc) /= 0) then
                imsg_ns_snd = imsg_ns_snd + 1
                msg_ns_snd(dst_proc) = imsg_ns_snd
            endif
            if (ns_rcv_count(dst_proc) /= 0) then
                imsg_ns_rcv = imsg_ns_rcv + 1
                msg_ns_rcv(dst_proc) = imsg_ns_rcv
            endif
        end do

        deallocate(ew_snd_count, ew_rcv_count, ns_snd_count, ns_rcv_count)

        !-----------------------------------------------------------------------
        !
        !  allocate buffers and arrays necessary for boundary comms
        !
        !-----------------------------------------------------------------------
        allocate (newbndy%nblocks_ew_snd(newbndy%nmsg_ew_snd), &
            newbndy%nblocks_ns_snd(newbndy%nmsg_ns_snd), &
            newbndy%nblocks_ew_rcv(newbndy%nmsg_ew_rcv), &
            newbndy%nblocks_ns_rcv(newbndy%nmsg_ns_rcv))

        allocate (newbndy%local_ew_src_block(newbndy%nlocal_ew), &
            newbndy%local_ew_dst_block(newbndy%nlocal_ew), &
            newbndy%local_ns_src_block(newbndy%nlocal_ns), &
            newbndy%local_ns_dst_block(newbndy%nlocal_ns), &
            newbndy%local_ew_src_add(2,newbndy%nlocal_ew), &
            newbndy%local_ew_dst_add(2,newbndy%nlocal_ew), &
            newbndy%local_ns_src_add(2,newbndy%nlocal_ns), &
            newbndy%local_ns_dst_add(2,newbndy%nlocal_ns))

        allocate ( &
        newbndy%ew_snd_proc (newbndy%nmsg_ew_snd), &
        newbndy%ew_rcv_proc (newbndy%nmsg_ew_rcv), &
        newbndy%ns_snd_proc (newbndy%nmsg_ns_snd), &
        newbndy%ns_rcv_proc (newbndy%nmsg_ns_rcv), &
        newbndy%ew_src_block(newbndy%maxblocks_ew_snd,newbndy%nmsg_ew_snd), &
        newbndy%ew_dst_block(newbndy%maxblocks_ew_rcv,newbndy%nmsg_ew_rcv), &
        newbndy%ns_src_block(newbndy%maxblocks_ns_snd,newbndy%nmsg_ns_snd), &
        newbndy%ns_dst_block(newbndy%maxblocks_ns_rcv,newbndy%nmsg_ns_rcv), &
        newbndy%ew_src_add(2,newbndy%maxblocks_ew_snd,newbndy%nmsg_ew_snd), &
        newbndy%ew_dst_add(2,newbndy%maxblocks_ew_rcv,newbndy%nmsg_ew_rcv), &
        newbndy%ns_src_add(2,newbndy%maxblocks_ns_snd,newbndy%nmsg_ns_snd), &
        newbndy%ns_dst_add(3,newbndy%maxblocks_ns_rcv,newbndy%nmsg_ns_rcv))

        newbndy%nblocks_ew_snd = 0
        newbndy%nblocks_ns_snd = 0
        newbndy%nblocks_ew_rcv = 0
        newbndy%nblocks_ns_rcv = 0

        newbndy%ew_snd_proc = 0
        newbndy%ew_rcv_proc = 0
        newbndy%ns_snd_proc = 0
        newbndy%ns_rcv_proc = 0

        newbndy%local_ew_src_block = 0
        newbndy%local_ew_dst_block = 0
        newbndy%local_ns_src_block = 0
        newbndy%local_ns_dst_block = 0
        newbndy%local_ew_src_add = 0
        newbndy%local_ew_dst_add = 0
        newbndy%local_ns_src_add = 0
        newbndy%local_ns_dst_add = 0

        newbndy%ew_src_block = 0
        newbndy%ew_dst_block = 0
        newbndy%ns_src_block = 0
        newbndy%ns_dst_block = 0
        newbndy%ew_src_add = 0
        newbndy%ew_dst_add = 0
        newbndy%ns_src_add = 0
        newbndy%ns_dst_add = 0

        !-----------------------------------------------------------------------
        !
        !  now set up indices into buffers and address arrays
        !
        !-----------------------------------------------------------------------

        allocate (ew_snd_count(newbndy%nmsg_ew_snd), &
            ew_rcv_count(newbndy%nmsg_ew_rcv), &
            ns_snd_count(newbndy%nmsg_ns_snd), &
            ns_rcv_count(newbndy%nmsg_ns_rcv))

        ew_snd_count = 0
        ew_rcv_count = 0
        ns_snd_count = 0
        ns_rcv_count = 0

        iloc_ew = 0
        iloc_ns = 0

        !-----------------------------------------------------------------------
        !
        !  repeat loop through blocks but this time, determine all the
        !  required message information for each message or local copy
        !
        !-----------------------------------------------------------------------

        !   print *,'before call to block_loop2:'
        block_loop2: do n=1,nblocks

            src_proc  = dist%proc(n)    ! processor location for this block

            call get_block_parameter( &
                n,iblock=iblock_src,jblock=jblock_src,ie=sbie, je=sbje)

            if (src_proc /= 0) then
                src_block_loc = dist%local_block(n)  ! local block location
            else
                src_block_loc = 0  ! block is a land block
            endif

            !*** compute cartesian i,j block indices for each neighbor
            !*** use zero if off the end of closed boundary
            !*** use jnorth=nblocks_y and inorth < 0 for tripole boundary
            !***   to make sure top boundary communicated to all top
            !***   boundary blocks
            call GetEWBlockIndex(ew_bndy_type,iblock_src,jblock_src, &
                iblock_east,jblock_east,iblock_west,jblock_west)

            call GetNSBlockIndex(ns_bndy_type,iblock_src, &
                jblock_src, iblock_north,jblock_north,iblock_south,jblock_south)

            !***
            !*** if blocks are at closed boundary, create local copy
            !*** pseudo-message to fill ghost cells
            !***
            if (src_block_loc /= 0 .and. src_proc == my_task+1) then
                if (iblock_east == 0) then
                    iloc_ew = iloc_ew + 1
                    newbndy%local_ew_src_block(iloc_ew) = 0
                    newbndy%local_ew_src_add(1,iloc_ew) = 0
                    newbndy%local_ew_src_add(2,iloc_ew) = 0
                    newbndy%local_ew_dst_block(iloc_ew) = src_block_loc
                    newbndy%local_ew_dst_add(1,iloc_ew) = sbie + 1
                    newbndy%local_ew_dst_add(2,iloc_ew) = 1
                else if (iblock_west == 0) then
                    iloc_ew = iloc_ew + 1
                    newbndy%local_ew_src_block(iloc_ew) = 0
                    newbndy%local_ew_src_add(1,iloc_ew) = 0
                    newbndy%local_ew_src_add(2,iloc_ew) = 0
                    newbndy%local_ew_dst_block(iloc_ew) = src_block_loc
                    newbndy%local_ew_dst_add(1,iloc_ew) = 1
                    newbndy%local_ew_dst_add(2,iloc_ew) = 1
                else if (jblock_north == 0) then
                    iloc_ns = iloc_ns + 1
                    newbndy%local_ns_src_block(iloc_ns) = 0
                    newbndy%local_ns_src_add(1,iloc_ns) = 0
                    newbndy%local_ns_src_add(2,iloc_ns) = 0
                    newbndy%local_ns_dst_block(iloc_ns) = src_block_loc
                    newbndy%local_ns_dst_add(1,iloc_ns) = 1
                    newbndy%local_ns_dst_add(2,iloc_ns) = sbje + 1
                else if (jblock_south == 0) then
                    iloc_ns = iloc_ns + 1
                    newbndy%local_ns_src_block(iloc_ns) = 0
                    newbndy%local_ns_src_add(1,iloc_ns) = 0
                    newbndy%local_ns_src_add(2,iloc_ns) = 0
                    newbndy%local_ns_dst_block(iloc_ns) = src_block_loc
                    newbndy%local_ns_dst_add(1,iloc_ns) = 1
                    newbndy%local_ns_dst_add(2,iloc_ns) = 1
                endif
            endif

            !***
            !*** now search through blocks looking for neighbors to
            !*** the source block
            !***
            do kk=1,NumNeigh
                k = Neigh(kk,n) 
                if(k>0) then 
                    dst_proc      = dist%proc(k)  ! processor holding dst block

                    !***
                    !*** compute the rest only if this block is not a land block
                    !***
                    if (dst_proc /= 0) then
                        call get_block_parameter( &
                            k,iblock=iblock_dst,jblock=jblock_dst,ie=dbie, &
                            je=dbje)

                        ! local block location
                        dst_block_loc = dist%local_block(k)  

                        !***
                        !*** if this block is an eastern neighbor
                        !*** determine send/receive addresses
                        !***
                        if (iblock_dst == iblock_east .and. &
                            jblock_dst == jblock_east) then

                            if (src_proc == my_task+1 .and. &
                                src_proc == dst_proc) then
                                !*** local copy from one block to another
                                iloc_ew = iloc_ew + 1
                                newbndy%local_ew_src_block(iloc_ew) = &
                                    src_block_loc
                                newbndy%local_ew_src_add(1,iloc_ew) = &
                                    sbie - nghost + 1
                                newbndy%local_ew_src_add(2,iloc_ew) = 1
                                newbndy%local_ew_dst_block(iloc_ew) = &
                                    dst_block_loc
                                newbndy%local_ew_dst_add(1,iloc_ew) = 1
                                newbndy%local_ew_dst_add(2,iloc_ew) = 1
                            else if (src_proc == 0 .and. &
                                dst_proc == my_task+1) then
                                !*** source block is all land so treat as local
                                !*** copy with source block zero to fill ghost
                                !*** cells with zeroes
                                iloc_ew = iloc_ew + 1
                                newbndy%local_ew_src_block(iloc_ew) = 0
                                newbndy%local_ew_src_add(1,iloc_ew) = 0
                                newbndy%local_ew_src_add(2,iloc_ew) = 0
                                newbndy%local_ew_dst_block(iloc_ew) = &
                                    dst_block_loc
                                newbndy%local_ew_dst_add(1,iloc_ew) = 1
                                newbndy%local_ew_dst_add(2,iloc_ew) = 1
                            else if (src_proc == my_task+1 .and. &
                                dst_proc /= my_task+1) then
                                !*** an actual message must be sent
                                imsg = msg_ew_snd(dst_proc)
                                ew_snd_count(imsg) = ew_snd_count(imsg) + 1
                                iblk = ew_snd_count(imsg)
                                newbndy%ew_snd_proc (     imsg) = dst_proc
                                newbndy%ew_src_block(iblk,imsg) = src_block_loc
                                newbndy%ew_src_add(1,iblk,imsg) = &
                                    sbie - nghost + 1
                                newbndy%ew_src_add(2,iblk,imsg) = 1
                                newbndy%nblocks_ew_snd(imsg) = &
                                newbndy%nblocks_ew_snd(imsg) + 1
                            else if (dst_proc == my_task+1 .and. &
                                src_proc /= my_task+1) then
                                !*** must receive a message
                                imsg = msg_ew_rcv(src_proc)
                                ew_rcv_count(imsg) = ew_rcv_count(imsg) + 1
                                iblk = ew_rcv_count(imsg)
                                newbndy%ew_rcv_proc (     imsg) = src_proc
                                newbndy%ew_dst_block(iblk,imsg) = dst_block_loc
                                newbndy%ew_dst_add(:,iblk,imsg) = 1
                                newbndy%nblocks_ew_rcv(imsg) = &
                                newbndy%nblocks_ew_rcv(imsg) + 1
                            endif
                        endif ! east neighbor

                        !***
                        !*** if this block is a western neighbor
                        !*** determine send/receive addresses
                        !***
                        if (iblock_dst == iblock_west .and. &
                            jblock_dst == jblock_west) then

                            if (src_proc == my_task+1 .and. &
                                src_proc == dst_proc) then
                                !*** perform a local copy
                                iloc_ew = iloc_ew + 1
                                newbndy%local_ew_src_block(iloc_ew) = &
                                    src_block_loc
                                newbndy%local_ew_src_add(1,iloc_ew) = nghost + 1
                                newbndy%local_ew_src_add(2,iloc_ew) = 1
                                newbndy%local_ew_dst_block(iloc_ew) = &
                                    dst_block_loc
                                newbndy%local_ew_dst_add(1,iloc_ew) = dbie + 1
                                newbndy%local_ew_dst_add(2,iloc_ew) = 1
                            else if (src_proc == 0 .and. &
                                dst_proc == my_task+1) then
                                !*** neighbor is a land block so zero ghost
                                !*** cells
                                iloc_ew = iloc_ew + 1
                                newbndy%local_ew_src_block(iloc_ew) = 0
                                newbndy%local_ew_src_add(1,iloc_ew) = 0
                                newbndy%local_ew_src_add(2,iloc_ew) = 0
                                newbndy%local_ew_dst_block(iloc_ew) = &
                                    dst_block_loc
                                newbndy%local_ew_dst_add(1,iloc_ew) = dbie + 1
                                newbndy%local_ew_dst_add(2,iloc_ew) = 1
                            else if (src_proc == my_task+1 .and. &
                                dst_proc /= my_task+1) then
                                !*** message must be sent
                                imsg = msg_ew_snd(dst_proc)
                                ew_snd_count(imsg) = ew_snd_count(imsg) + 1
                                iblk = ew_snd_count(imsg)
                                newbndy%ew_snd_proc (     imsg) = dst_proc
                                newbndy%ew_src_block(iblk,imsg) = src_block_loc
                                newbndy%ew_src_add(1,iblk,imsg) = nghost + 1
                                newbndy%ew_src_add(2,iblk,imsg) = 1
                                newbndy%nblocks_ew_snd(imsg) = &
                                newbndy%nblocks_ew_snd(imsg) + 1
                            else if (dst_proc == my_task+1 .and. &
                                src_proc /= my_task+1) then
                                !*** message must be received
                                imsg = msg_ew_rcv(src_proc)
                                ew_rcv_count(imsg) = ew_rcv_count(imsg) + 1
                                iblk = ew_rcv_count(imsg)
                                newbndy%ew_rcv_proc (     imsg) = src_proc
                                newbndy%ew_dst_block(iblk,imsg) = dst_block_loc
                                newbndy%ew_dst_add(1,iblk,imsg) = dbie + 1
                                newbndy%ew_dst_add(2,iblk,imsg) = 1
                                newbndy%nblocks_ew_rcv(imsg) = &
                                newbndy%nblocks_ew_rcv(imsg) + 1
                            endif
                        endif ! west neighbor

                        !***
                        !*** if this block is a northern neighbor
                        !***  compute send/recv addresses
                        !*** for tripole, must communicate with all
                        !*** north row blocks (triggered by iblock_dst <0)
                        !***
                        if ((iblock_dst == iblock_north .or. &
                            iblock_north < 0) .and. &
                            jblock_dst == jblock_north) then

                            if (src_proc == my_task+1 .and. &
                                src_proc == dst_proc) then
                                !*** local copy
                                iloc_ns = iloc_ns + 1
                                newbndy%local_ns_src_block(iloc_ns) = &
                                    src_block_loc
                                newbndy%local_ns_src_add(1,iloc_ns) = 1
                                newbndy%local_ns_src_add(2,iloc_ns) = &
                                    sbje - nghost + 1
                                newbndy%local_ns_dst_block(iloc_ns) = &
                                    dst_block_loc
                                newbndy%local_ns_dst_add(1,iloc_ns) = 1
                                newbndy%local_ns_dst_add(2,iloc_ns) = 1

                            else if (src_proc == 0 .and. &
                                dst_proc == my_task+1) then
                                !*** source is land block so zero ghost cells
                                iloc_ns = iloc_ns + 1
                                newbndy%local_ns_src_block(iloc_ns) = 0
                                newbndy%local_ns_src_add(1,iloc_ns) = 0
                                newbndy%local_ns_src_add(2,iloc_ns) = 0
                                newbndy%local_ns_dst_block(iloc_ns) = &
                                    dst_block_loc
                                newbndy%local_ns_dst_add(1,iloc_ns) = 1
                                newbndy%local_ns_dst_add(2,iloc_ns) = 1

                            else if (src_proc == my_task+1 .and. &
                                dst_proc /= my_task+1) then
                                !*** message must be sent
                                imsg = msg_ns_snd(dst_proc)
                                ns_snd_count(imsg) = ns_snd_count(imsg) + 1
                                iblk = ns_snd_count(imsg)
                                newbndy%ns_snd_proc (     imsg) = dst_proc
                                newbndy%ns_src_block(iblk,imsg) = src_block_loc
                                newbndy%ns_src_add(1,iblk,imsg) = 1
                                newbndy%ns_src_add(2,iblk,imsg) = &
                                    sbje - nghost + 1
                                newbndy%nblocks_ns_snd(imsg) = &
                                newbndy%nblocks_ns_snd(imsg) + 1

                            else if (dst_proc == my_task+1 .and. &
                                src_proc /= my_task+1) then
                                !*** message must be received
                                imsg = msg_ns_rcv(src_proc)
                                ns_rcv_count(imsg) = ns_rcv_count(imsg) + 1
                                iblk = ns_rcv_count(imsg)
                                newbndy%ns_rcv_proc (     imsg) = src_proc
                                newbndy%ns_dst_block(iblk,imsg) = dst_block_loc
                                newbndy%ns_dst_add(1,iblk,imsg) = 1
                                newbndy%ns_dst_add(2,iblk,imsg) = 1
                                newbndy%nblocks_ns_rcv(imsg) = &
                                newbndy%nblocks_ns_rcv(imsg) + 1
                            endif
                        endif ! north neighbor

                        !***
                        !*** if this block is a southern neighbor
                        !*** determine send/receive addresses
                        !***
                        if (iblock_dst == iblock_south .and. &
                            jblock_dst == jblock_south) then

                            if (src_proc == my_task+1 .and. &
                                src_proc == dst_proc) then
                                !*** local copy
                                iloc_ns = iloc_ns + 1
                                newbndy%local_ns_src_block(iloc_ns) = &
                                    src_block_loc
                                newbndy%local_ns_src_add(1,iloc_ns) = 1
                                newbndy%local_ns_src_add(2,iloc_ns) = nghost + 1
                                newbndy%local_ns_dst_block(iloc_ns) = &
                                    dst_block_loc
                                newbndy%local_ns_dst_add(1,iloc_ns) = 1
                                newbndy%local_ns_dst_add(2,iloc_ns) = dbje + 1
                            else if (src_proc == 0 .and. &
                                dst_proc == my_task+1) then
                                !*** neighbor is a land block so zero ghost
                                !*** cells
                                iloc_ns = iloc_ns + 1
                                newbndy%local_ns_src_block(iloc_ns) = 0
                                newbndy%local_ns_src_add(1,iloc_ns) = 0
                                newbndy%local_ns_src_add(2,iloc_ns) = 0
                                newbndy%local_ns_dst_block(iloc_ns) = &
                                    dst_block_loc
                                newbndy%local_ns_dst_add(1,iloc_ns) = 1
                                newbndy%local_ns_dst_add(2,iloc_ns) = dbje + 1
                            else if (src_proc == my_task+1 .and. &
                                dst_proc /= my_task+1) then
                                !*** message must be sent
                                imsg = msg_ns_snd(dst_proc)
                                ns_snd_count(imsg) = ns_snd_count(imsg) + 1
                                iblk = ns_snd_count(imsg)
                                newbndy%ns_snd_proc (     imsg) = dst_proc
                                newbndy%ns_src_block(iblk,imsg) = src_block_loc
                                newbndy%ns_src_add(1,iblk,imsg) = 1
                                newbndy%ns_src_add(2,iblk,imsg) = nghost+1
                                newbndy%nblocks_ns_snd(imsg) = &
                                    newbndy%nblocks_ns_snd(imsg) + 1
                            else if (dst_proc == my_task+1 .and. &
                                src_proc /= my_task+1) then
                                !*** message must be received
                                imsg = msg_ns_rcv(src_proc)
                                ns_rcv_count(imsg) = ns_rcv_count(imsg) + 1
                                iblk = ns_rcv_count(imsg)
                                newbndy%ns_rcv_proc (     imsg) = src_proc
                                newbndy%ns_dst_block(iblk,imsg) = dst_block_loc
                                newbndy%ns_dst_add(1,iblk,imsg) = 1
                                newbndy%ns_dst_add(2,iblk,imsg) = dbje + 1
                                newbndy%nblocks_ns_rcv(imsg) = &
                                    newbndy%nblocks_ns_rcv(imsg) + 1
                            endif
                        endif ! south neighbor
                    endif  ! not a land block
                endif ! k>0
            end do ! loop over Neighbors
        end do block_loop2

        deallocate(ew_snd_count, ew_rcv_count, ns_snd_count, ns_rcv_count)
        deallocate(msg_ew_snd, msg_ew_rcv, msg_ns_snd, msg_ns_rcv)
    end subroutine mpi2s_create_boundary

    !***********************************************************************

    !>
    !! This routine destroys a boundary by deallocating all memory
    !! associated with the boundary and nullifying pointers.
    !!
    !! @param in_bndy   boundary structure to be destroyed
    !<
    subroutine mpi2s_destroy_boundary(in_bndy)
        ! !INPUT/OUTPUT PARAMETERS:
        type (bndy), intent(inout) :: in_bndy

        !-------------------------------------------------------------------
        !
        !  reset all scalars
        !
        !-------------------------------------------------------------------
        in_bndy%communicator      = 0
        in_bndy%nmsg_ew_snd       = 0
        in_bndy%nmsg_ns_snd       = 0
        in_bndy%nmsg_ew_rcv       = 0
        in_bndy%nmsg_ns_rcv       = 0
        in_bndy%maxblocks_ew_snd  = 0
        in_bndy%maxblocks_ew_rcv  = 0
        in_bndy%maxblocks_ns_snd  = 0
        in_bndy%maxblocks_ns_rcv  = 0
        in_bndy%nlocal_ew         = 0
        in_bndy%nlocal_ns         = 0

        !-------------------------------------------------------------------
        !
        !  deallocate all pointers
        !
        !-------------------------------------------------------------------
        deallocate(in_bndy%nblocks_ew_snd,     &
            in_bndy%nblocks_ns_snd,     &
            in_bndy%nblocks_ew_rcv,     &
            in_bndy%nblocks_ns_rcv,     &
            in_bndy%ew_snd_proc,        &
            in_bndy%ew_rcv_proc,        &
            in_bndy%ns_snd_proc,        &
            in_bndy%ns_rcv_proc,        &
            in_bndy%local_ew_src_block, &
            in_bndy%local_ew_dst_block, &
            in_bndy%local_ns_src_block, &
            in_bndy%local_ns_dst_block, &
            in_bndy%local_ew_src_add,   &
            in_bndy%local_ew_dst_add,   &
            in_bndy%local_ns_src_add,   &
            in_bndy%local_ns_dst_add,   &
            in_bndy%ew_src_block,       &
            in_bndy%ew_dst_block,       &
            in_bndy%ns_src_block,       &
            in_bndy%ns_dst_block,       &
            in_bndy%ew_src_add,         &
            in_bndy%ew_dst_add,         &
            in_bndy%ns_src_add,         &
            in_bndy%ns_dst_add )
    end subroutine mpi2s_destroy_boundary

    !***********************************************************************

    !>
    !! This routine updates ghost cells for an input array and is a
    !! member of a group of routines under the generic interface
    !! update\_ghost\_cells.  This routine is the specific interface
    !! for 2d horizontal arrays of double precision.
    !!
    !! @param ARRAY         array containing horizontal slab to update
    !! @param in_bndy       boundary update structure for the array
    !! @param grid_loc      id for location on horizontal grid (center,
    !!                      NEcorner, Nface, Eface)
    !! @param field_type    id for type of field (scalar, vector, angle)
    !<
    subroutine mpi2s_boundary_2d_dbl(ARRAY, in_bndy, grid_loc, field_type)
        ! !USER:
        include 'mpif.h'   ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        type (bndy), intent(in) ::  in_bndy

        integer (i4), intent(in) :: &
            field_type,              &
            grid_loc

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), dimension(:,:,:), intent(inout) :: ARRAY

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer (i4) ::           &
            i,j,k,m,n,                   &! dummy loop indices
            ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
            ib_dst,ie_dst,jb_dst,je_dst, &!
            nx_global,                   &! global domain size in x
            src_block,                   &! local block number for source
            dst_block,                   &! local block number for destination
            bufsize,                     &! buffer size for send/recv buffers
            xoffset, yoffset,            &! address shifts for tripole
            isign,                       &! sign factor for tripole grids
            ierr                          ! MPI error flag

        integer (i4), dimension(:), allocatable :: &
            snd_request,              &! MPI request ids
            rcv_request                ! MPI request ids

        integer (i4), dimension(:,:), allocatable :: &
            snd_status,               &! MPI status flags
            rcv_status                 ! MPI status flags

        real (r8), dimension(:,:,:,:), allocatable :: &
            buf_ew_snd,       &! message buffer for east-west sends
            buf_ew_rcv,       &! message buffer for east-west recvs
            buf_ns_snd,       &! message buffer for north-south sends
            buf_ns_rcv         ! message buffer for north-south recvs

        real (r8) :: &
            xavg               ! scalar for enforcing symmetry at U pts

        !-------------------------------------------------------------------
        !
        !  allocate buffers for east-west sends and receives
        !
        !-------------------------------------------------------------------
        allocate(buf_ew_snd(nghost, ny_block, &
            in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block, &
            in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

        allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

        !-------------------------------------------------------------------
        !
        !  post receives
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nmsg_ew_rcv
            bufsize = ny_block*nghost*in_bndy%nblocks_ew_rcv(n)

            call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                in_bndy%ew_rcv_proc(n)-1,                &
                mpitag_bndy_2d + in_bndy%ew_rcv_proc(n), &
                in_bndy%communicator, rcv_request(n), ierr)
        end do

        !-------------------------------------------------------------------
        !
        !  fill send buffer and post sends
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nmsg_ew_snd
            bufsize = ny_block*nghost*in_bndy%nblocks_ew_snd(n)

            do i=1,in_bndy%nblocks_ew_snd(n)
                ib_src    = in_bndy%ew_src_add(1,i,n)
                ie_src    = ib_src + nghost - 1
                src_block = in_bndy%ew_src_block(i,n)
                buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
            end do

            call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_dbl, &
                in_bndy%ew_snd_proc(n)-1, &
                mpitag_bndy_2d + my_task + 1, &
                in_bndy%communicator, snd_request(n), ierr)
        end do

        !-------------------------------------------------------------------
        !
        !  do local copies while waiting for messages to complete
        !  also initialize ghost cells to zero
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nlocal_ew
            src_block = in_bndy%local_ew_src_block(n)
            dst_block = in_bndy%local_ew_dst_block(n)

            ib_src = in_bndy%local_ew_src_add(1,n)
            ie_src = ib_src + nghost - 1
            ib_dst = in_bndy%local_ew_dst_add(1,n)
            ie_dst = ib_dst + nghost - 1

            if (src_block /= 0) then
                ARRAY(ib_dst:ie_dst,:,dst_block) = &
                ARRAY(ib_src:ie_src,:,src_block)
            else
                ARRAY(ib_dst:ie_dst,:,dst_block) = c0
            endif
        end do

        !-------------------------------------------------------------------
        !
        !  wait for receives to finish and then unpack the recv buffer into
        !  ghost cells
        !
        !-------------------------------------------------------------------
        if(in_bndy%nmsg_ew_rcv > 0)  &
            call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr)

        do n=1,in_bndy%nmsg_ew_rcv
            do k=1,in_bndy%nblocks_ew_rcv(n)
                dst_block = in_bndy%ew_dst_block(k,n)

                ib_dst = in_bndy%ew_dst_add(1,k,n)
                ie_dst = ib_dst + nghost - 1

                ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
            end do
        end do

        !-------------------------------------------------------------------
        !
        !  wait for sends to complete and deallocate arrays
        !
        !-------------------------------------------------------------------
        if(in_bndy%nmsg_ew_snd > 0)  &
            call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)

        deallocate(buf_ew_snd, buf_ew_rcv)
        deallocate(snd_request, rcv_request, snd_status, rcv_status)

        !-------------------------------------------------------------------
        !
        !  now exchange north-south boundary info
        !
        !-------------------------------------------------------------------
        allocate(buf_ns_snd(nx_block, nghost+1, &
            in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block, nghost+1, &
            in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

        allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

        !-------------------------------------------------------------------
        !
        !  post receives
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nmsg_ns_rcv

            bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

            call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_dbl,   &
                in_bndy%ns_rcv_proc(n)-1,                &
                mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                in_bndy%communicator, rcv_request(n), ierr)
        end do

        !-------------------------------------------------------------------
        !
        !  fill send buffer and post sends
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nmsg_ns_snd
            bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_snd(n)

            do i=1,in_bndy%nblocks_ns_snd(n)
                jb_src    = in_bndy%ns_src_add(2,i,n)
                je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
                src_block = in_bndy%ns_src_block(i,n)
                buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
            end do

            call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_dbl, &
                in_bndy%ns_snd_proc(n)-1, &
                mpitag_bndy_2d + my_task + 1, &
                in_bndy%communicator, snd_request(n), ierr)
        end do

        !-------------------------------------------------------------------
        !
        !  do local copies while waiting for messages to complete
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nlocal_ns
            src_block = in_bndy%local_ns_src_block(n)
            dst_block = in_bndy%local_ns_dst_block(n)

            if (dst_block > 0) then ! straight local copy

                jb_src = in_bndy%local_ns_src_add(2,n)
                je_src = jb_src + nghost - 1
                jb_dst = in_bndy%local_ns_dst_add(2,n)
                je_dst = jb_dst + nghost - 1

                if (src_block /= 0) then
                    ARRAY(:,jb_dst:je_dst,dst_block) = &
                    ARRAY(:,jb_src:je_src,src_block)
                else
                    ARRAY(:,jb_dst:je_dst,dst_block) = c0
                endif

            endif
        end do

        !-------------------------------------------------------------------
        !
        !  wait for receives to finish and then unpack the recv buffer into
        !  ghost cells
        !
        !-------------------------------------------------------------------
        if(in_bndy%nmsg_ns_rcv > 0)  &
            call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)

        do n=1,in_bndy%nmsg_ns_rcv
            do k=1,in_bndy%nblocks_ns_rcv(n)
                dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

                if (dst_block > 0) then  ! normal receive
                    jb_dst = in_bndy%ns_dst_add(2,k,n)
                    je_dst = jb_dst + nghost - 1

                    ARRAY(:,jb_dst:je_dst,dst_block) = &
                        buf_ns_rcv(:,1:nghost,k,n)
                endif
            end do
        end do

        !-------------------------------------------------------------------
        !
        !  wait for sends to complete and deallocate arrays
        !
        !-------------------------------------------------------------------
        if(in_bndy%nmsg_ns_snd > 0)  &
            call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)

        deallocate(buf_ns_snd, buf_ns_rcv)
        deallocate(snd_request, rcv_request, snd_status, rcv_status)
    end subroutine mpi2s_boundary_2d_dbl

    !***********************************************************************
    !>
    !! This routine updates ghost cells for an input array and is a
    !! member of a group of routines under the generic interface
    !! update\_ghost\_cells.  This routine is the specific interface
    !! for 2d horizontal arrays of double precision.
    !!
    !! @param ARRAY         array containing horizontal slab to update
    !! @param in_bndy       boundary update structure for the array
    !! @param grid_loc      id for location on horizontal grid (center,
    !!                      NEcorner, Nface, Eface)
    !! @param field_type    id for type of field (scalar, vector, angle)
    !<
    subroutine mpi2s_boundary_2d_int(ARRAY, in_bndy, grid_loc, field_type)
        ! !USER:
        include 'mpif.h'   ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        type (bndy), intent(in) :: in_bndy

        integer (i4), intent(in) :: &
            field_type,               &
            grid_loc

        ! !INPUT/OUTPUT PARAMETERS:
        integer (i4), dimension(:,:,:), intent(inout) :: ARRAY

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer (i4) ::           &
            i,j,k,m,n,                   &! dummy loop indices
            ib_src,ie_src,jb_src,je_src, &! beg,end indices for bndy cells
            ib_dst,ie_dst,jb_dst,je_dst, &!
            nx_global,                   &! global domain size in x
            src_block,                   &! local block number for source
            dst_block,                   &! local block number for destination
            bufsize,                     &! buffer size for send/recv buffers
            xoffset, yoffset,            &! address shifts for tripole
            isign,                       &! sign factor for tripole grids
            ierr                          ! MPI error flag

        integer (i4), dimension(:), allocatable :: &
            snd_request,              &! MPI request ids
            rcv_request                ! MPI request ids

        integer (i4), dimension(:,:), allocatable :: &
            snd_status,               &! MPI status flags
            rcv_status                 ! MPI status flags

        integer (i4), dimension(:,:,:,:), allocatable :: &
            buf_ew_snd,       &! message buffer for east-west sends
            buf_ew_rcv,       &! message buffer for east-west recvs
            buf_ns_snd,       &! message buffer for north-south sends
            buf_ns_rcv         ! message buffer for north-south recvs

        integer (i4) :: &
            xavg               ! scalar for enforcing symmetry at U pts

        !-------------------------------------------------------------------
        !
        !  allocate buffers for east-west sends and receives
        !
        !-------------------------------------------------------------------
        allocate(buf_ew_snd(nghost, ny_block, &
            in_bndy%maxblocks_ew_snd, in_bndy%nmsg_ew_snd),&
            buf_ew_rcv(nghost, ny_block, &
            in_bndy%maxblocks_ew_rcv, in_bndy%nmsg_ew_rcv))

        allocate(snd_request(in_bndy%nmsg_ew_snd), &
            rcv_request(in_bndy%nmsg_ew_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ew_rcv))

        !-------------------------------------------------------------------
        !
        !  post receives
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nmsg_ew_rcv

            bufsize = ny_block*nghost*in_bndy%nblocks_ew_rcv(n)

            call MPI_IRECV(buf_ew_rcv(1,1,1,n), bufsize, mpi_integer, &
                in_bndy%ew_rcv_proc(n)-1,                  &
                mpitag_bndy_2d + in_bndy%ew_rcv_proc(n),   &
                in_bndy%communicator, rcv_request(n), ierr)
        end do

        !-------------------------------------------------------------------
        !
        !  fill send buffer and post sends
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nmsg_ew_snd
            bufsize = ny_block*nghost*in_bndy%nblocks_ew_snd(n)

            do i=1,in_bndy%nblocks_ew_snd(n)
                ib_src    = in_bndy%ew_src_add(1,i,n)
                ie_src    = ib_src + nghost - 1
                src_block = in_bndy%ew_src_block(i,n)
                buf_ew_snd(:,:,i,n) = ARRAY(ib_src:ie_src,:,src_block)
            end do

            call MPI_ISEND(buf_ew_snd(1,1,1,n), bufsize, mpi_integer, &
                in_bndy%ew_snd_proc(n)-1, &
                mpitag_bndy_2d + my_task + 1, &
                in_bndy%communicator, snd_request(n), ierr)
        end do

        !-------------------------------------------------------------------
        !
        !  do local copies while waiting for messages to complete
        !  also initialize ghost cells to zero
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nlocal_ew
            src_block = in_bndy%local_ew_src_block(n)
            dst_block = in_bndy%local_ew_dst_block(n)

            ib_src = in_bndy%local_ew_src_add(1,n)
            ie_src = ib_src + nghost - 1
            ib_dst = in_bndy%local_ew_dst_add(1,n)
            ie_dst = ib_dst + nghost - 1

            if (src_block /= 0) then
                ARRAY(ib_dst:ie_dst,:,dst_block) = &
                ARRAY(ib_src:ie_src,:,src_block)
            else
                ARRAY(ib_dst:ie_dst,:,dst_block) = 0
            endif
        end do

        !-------------------------------------------------------------------
        !
        !  wait for receives to finish and then unpack the recv buffer into
        !  ghost cells
        !
        !-------------------------------------------------------------------
        if(in_bndy%nmsg_ew_rcv > 0)  &
            call MPI_WAITALL(in_bndy%nmsg_ew_rcv, rcv_request, rcv_status, ierr)

        do n=1,in_bndy%nmsg_ew_rcv
            do k=1,in_bndy%nblocks_ew_rcv(n)
                dst_block = in_bndy%ew_dst_block(k,n)

                ib_dst = in_bndy%ew_dst_add(1,k,n)
                ie_dst = ib_dst + nghost - 1

                ARRAY(ib_dst:ie_dst,:,dst_block) = buf_ew_rcv(:,:,k,n)
            end do
        end do

        !-------------------------------------------------------------------
        !
        !  wait for sends to complete and deallocate arrays
        !
        !-------------------------------------------------------------------
        if(in_bndy%nmsg_ew_snd > 0)  &
            call MPI_WAITALL(in_bndy%nmsg_ew_snd, snd_request, snd_status, ierr)

        deallocate(buf_ew_snd, buf_ew_rcv)
        deallocate(snd_request, rcv_request, snd_status, rcv_status)

        !-------------------------------------------------------------------
        !
        !  now exchange north-south boundary info
        !
        !-------------------------------------------------------------------
        allocate(buf_ns_snd(nx_block, nghost+1, &
            in_bndy%maxblocks_ns_snd, in_bndy%nmsg_ns_snd),&
            buf_ns_rcv(nx_block, nghost+1, &
            in_bndy%maxblocks_ns_rcv, in_bndy%nmsg_ns_rcv))

        allocate(snd_request(in_bndy%nmsg_ns_snd), &
            rcv_request(in_bndy%nmsg_ns_rcv), &
            snd_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_snd), &
            rcv_status(MPI_STATUS_SIZE,in_bndy%nmsg_ns_rcv))

        !-------------------------------------------------------------------
        !
        !  post receives
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nmsg_ns_rcv
            bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_rcv(n)

            call MPI_IRECV(buf_ns_rcv(1,1,1,n), bufsize, mpi_integer,   &
                in_bndy%ns_rcv_proc(n)-1,                &
                mpitag_bndy_2d + in_bndy%ns_rcv_proc(n), &
                in_bndy%communicator, rcv_request(n), ierr)
        end do

        !-------------------------------------------------------------------
        !
        !  fill send buffer and post sends
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nmsg_ns_snd

            bufsize = nx_block*(nghost+1)*in_bndy%nblocks_ns_snd(n)

            do i=1,in_bndy%nblocks_ns_snd(n)
                jb_src    = in_bndy%ns_src_add(2,i,n)
                je_src    = jb_src + nghost  ! nghost+1 rows needed for tripole
                src_block = in_bndy%ns_src_block(i,n)
                buf_ns_snd(:,:,i,n) = ARRAY(:,jb_src:je_src,src_block)
            end do

            call MPI_ISEND(buf_ns_snd(1,1,1,n), bufsize, mpi_integer, &
                in_bndy%ns_snd_proc(n)-1, &
                mpitag_bndy_2d + my_task + 1, &
                in_bndy%communicator, snd_request(n), ierr)
        end do

        !-------------------------------------------------------------------
        !
        !  do local copies while waiting for messages to complete
        !
        !-------------------------------------------------------------------
        do n=1,in_bndy%nlocal_ns
            src_block = in_bndy%local_ns_src_block(n)
            dst_block = in_bndy%local_ns_dst_block(n)

            if (dst_block > 0) then ! straight local copy

                jb_src = in_bndy%local_ns_src_add(2,n)
                je_src = jb_src + nghost - 1
                jb_dst = in_bndy%local_ns_dst_add(2,n)
                je_dst = jb_dst + nghost - 1

                if (src_block /= 0) then
                    ARRAY(:,jb_dst:je_dst,dst_block) = &
                    ARRAY(:,jb_src:je_src,src_block)
                else
                    ARRAY(:,jb_dst:je_dst,dst_block) = 0
                endif

            endif
        end do

        !-------------------------------------------------------------------
        !
        !  wait for receives to finish and then unpack the recv buffer into
        !  ghost cells
        !
        !-------------------------------------------------------------------
        if(in_bndy%nmsg_ns_rcv > 0)  &
            call MPI_WAITALL(in_bndy%nmsg_ns_rcv, rcv_request, rcv_status, ierr)

        do n=1,in_bndy%nmsg_ns_rcv
            do k=1,in_bndy%nblocks_ns_rcv(n)
                dst_block = in_bndy%ns_dst_block(k,n)  ! dest block

                if (dst_block > 0) then  ! normal receive
                    jb_dst = in_bndy%ns_dst_add(2,k,n)
                    je_dst = jb_dst + nghost - 1

                    ARRAY(:,jb_dst:je_dst,dst_block) = &
                        buf_ns_rcv(:,1:nghost,k,n)
                endif
            end do
        end do

        !-------------------------------------------------------------------
        !
        !  wait for sends to complete and deallocate arrays
        !
        !-------------------------------------------------------------------
        if(in_bndy%nmsg_ns_snd > 0)  &
            call MPI_WAITALL(in_bndy%nmsg_ns_snd, snd_request, snd_status, ierr)

        deallocate(buf_ns_snd, buf_ns_rcv)
        deallocate(snd_request, rcv_request, snd_status, rcv_status)
    end subroutine mpi2s_boundary_2d_int

    !***********************************************************************
    !>
    !! This is a utility routine to increment the arrays for counting
    !! whether messages are required.  It is used only for creating
    !! boundary structures for updating ghost cells.
    !!
    !! @param snd_counter   array for counting messages to be sent
    !! @param rcv_counter   array for counting messages to be received
    !! @param src_proc      source processor for communication
    !! @param dst_proc      destination processor for communication
    !<
    subroutine increment_message_counter(snd_counter, rcv_counter, &
        src_proc, dst_proc)

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
            src_proc,                &
            dst_proc

        ! !INPUT/OUTPUT PARAMETERS:
        integer (i4), dimension(:), intent(inout) :: &
            snd_counter,       &
            rcv_counter

        !-------------------------------------------------------------------
        !
        !  if destination all land (proc = 0), then no send is necessary, 
        !  so do the rest only for proc /= 0
        !
        !-------------------------------------------------------------------
        if (dst_proc /= 0) then
            !*** if the current processor is the source, must send
            !*** data (local copy if dst_proc = src_proc)

            if (src_proc == my_task + 1) snd_counter(dst_proc) = &
                snd_counter(dst_proc) + 1

            !*** if the current processor is the destination, must
            !*** receive data (local copy if dst_proc = src_proc)
            if (dst_proc == my_task + 1) then
                if (src_proc /= 0) then  
                    !*** the source block has ocean points so
                    !*** increment the number of messages from
                    !*** the source process
                    rcv_counter(src_proc) = rcv_counter(src_proc) + 1
                else
                    !*** the source block has no ocean points so
                    !*** count this as a local copy in order to
                    !*** fill ghost cells with zeroes
                    rcv_counter(dst_proc) = rcv_counter(dst_proc) + 1
                endif
            endif
        endif
    end subroutine increment_message_counter
end module mpi2s_boundary
