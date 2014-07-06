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
module caf_pull_boundary
    ! !USES:
    use kinds_mod, only: i4, r4, r8, log_kind
    use simple_type, only : distrb
    use communicate, only: mpitag_bndy_2d, my_task, MPI_DBL, MPI_COMM_OCN
    use constants, only: field_loc_necorner, field_loc_center, &
        field_loc_eface, field_loc_nface, field_type_vector, &
        field_type_angle, field_type_scalar, p5, c0
    use simple_blocks
    use exit_mod, only: exit_POP, sigAbort
    use boundary_types, only: bndy
    use domain_size, only: block_size_x, block_size_y

    implicit none
    private
    save

    ! !MODULE VARS:
    integer(i4) :: nblocks
    integer(i4) :: nLocalBlocks
    integer(i4), allocatable :: local_blocks(:)
    type (distrb) :: mDist

    integer(i4), public :: timer_caf_boundary_create, &
        timer_caf_boundary_2d_dbl

    ! !PUBLIC MEMBER FUNCTIONS:
    public :: caf_pull_create_boundary,  &
        caf_pull_destroy_boundary, &
        caf_pull_boundary_2d

    interface caf_pull_boundary_2d
        module procedure caf_boundary_2d_dbl
        module procedure caf_boundary_2d_int
    end interface

  contains

    !***********************************************************************
    !>
    !! This routine creates a boundary type with info necessary for
    !! performing a boundary (ghost cell) update based on the input block
    !! distribution.
    !!
    !! @param newbndy       a new boundary type with info for updates
    !! @param dist          distribution of blocks across procs
    !! @param blocks        A data structure for the application specific
    !!                      metadata.
    !! @param ns_bndy_type  type of boundary to use in ns dir
    !! @param ew_bndy_type  type of boundary to use in ew dir
    !! @param nx_global     global X extent of domain
    !! @param ny_global     global Y extent of domain
    !<
    subroutine caf_pull_create_boundary(newbndy, dist, blocks, &
        ns_bndy_type, ew_bndy_type, nx_global, ny_global)

        ! !INPUT PARAMETERS:

        type (distrb), intent(in) :: dist

        integer(i4) :: blocks(:) 

        character (*), intent(in) :: &
            ns_bndy_type,             &
            ew_bndy_type

        integer (i4), intent(in) :: &
            nx_global, ny_global

        ! !OUTPUT PARAMETERS:
        type (bndy), intent(out) :: newbndy

        !  local variables
        integer (i4) ::           &
            n,                    & ! loop counter
            src_proc                ! proccess that owns block n

        !-------------------------------------------------------------------
        !
        !  Initialize some useful variables and return if this task not
        !  in the current distribution.
        !
        !-------------------------------------------------------------------
        mDist = dist
        nblocks = size(mDist%proc(:))
        allocate(local_blocks(nblocks))

        if (my_task >= dist%nprocs) return

        !-------------------------------------------------------------------
        !
        !  Iterate through all blocks, pick out the ones that are local and
        !  store them in the local_blocks array
        !
        !-------------------------------------------------------------------
        nLocalBlocks = 0
        block_loop: do n=1,nblocks
            src_proc  = mDist%proc(n)

            if(src_proc /= this_image()) cycle block_loop

            local_blocks(nLocalBlocks+1) = n
            nLocalBlocks = nLocalBlocks + 1
        end do block_loop
    end subroutine caf_pull_create_boundary

    !***********************************************************************
    !>
    !! This routine destroys a boundary by deallocating all memory
    !! associated with the boundary and nullifying pointers.
    !!
    !! @param in_bndy   boundary structure to be destroyed
    !<
    subroutine caf_pull_destroy_boundary(in_bndy)
        ! !INPUT/OUTPUT PARAMETERS:
        type (bndy), intent(inout) :: in_bndy

        deallocate(local_blocks)
    end subroutine caf_pull_destroy_boundary

    !***********************************************************************
    subroutine caf_boundary_2d_dbl(ARRAY, in_bndy, grid_loc, field_type)
        ! !USER:
        include 'mpif.h'   ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        type (bndy), intent(in) :: &
            in_bndy                 ! boundary update structure for the array

        integer (i4), intent(in) :: &
            field_type,             & ! id for type of field (scalar, vector,
                                    & ! angle)
            grid_loc                  ! id for location on horizontal grid
                                      !  (center, NEcorner, Nface, Eface)

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), intent(inout) :: &
            ARRAY(:,:,:)[*]    ! array containing horizontal slab to update

        ! ! Locals:
            integer(i4) :: i, glbBlk
            integer(i4) :: ib, jb, ie, je ! bounds of block (nonghostcell) data
            integer(i4) :: proc_west, proc_east, proc_north, proc_south
            integer(i4) :: block_west, block_east, block_north, block_south

        ! ***
        ! *** Set border elements to zeroes
        ! ***
        do i=1,nLocalBlocks
            ! get block bounds information for the current block
            glbBlk=local_blocks(i)
            call get_block_parameter(glbBlk, ib=ib, ie=ie, jb=jb, je=je)

            ! clear ghost cells in each of the four borders
            ARRAY(1:nghost,  jb:je, i) = 0              ! west
            ARRAY(ie+1:ie+nghost,  jb:je, i) = 0        ! east
            ARRAY(1:ie+nghost,  1:jb-1, i) = 0          ! south
            ARRAY(1:ie+nghost,  je+1:je+nghost, i) = 0  ! north
        end do
        sync all

        ! ***
        ! *** Iterate through local blocks perform ew boundary exchange
        ! ***
        do i=1,nLocalBlocks
            ! ***
            ! *** Get information about block
            ! ***
            glbBlk=local_blocks(i)
            call get_block_parameter(glbBlk, ib=ib, ie=ie, jb=jb, je=je)

            if(Neigh(iwest, glbBlk) == 0) then
                proc_west = 0
                block_west = 0
            else
                proc_west  = mDist%proc(Neigh(iwest, glbBlk))
                block_west = mDist%local_block(Neigh(iwest, glbBlk))
            end if
            if(Neigh(ieast, glbBlk) == 0) then
                proc_east = 0
                block_east = 0
            else
                proc_east  = mDist%proc(Neigh(ieast, glbBlk))
                block_east = mDist%local_block(Neigh(ieast, glbBlk))
            end if

            ! ***
            ! *** Grab east
            ! ***
            if(proc_east /= 0) then
                ARRAY(ie+1:nx_block, jb:je, i) = ARRAY( &
                    ib:ib+nghost-1, jb:je, block_east)[proc_east]
            end if

            ! ***
            ! *** Grab west
            ! ***
            if(proc_west /= 0) then
                ARRAY(1:nghost, jb:je, i) = ARRAY( &
                    ie-nghost+1:ie, jb:je, block_west)[proc_west]
            end if
        end do

        sync all

        ! ***
        ! *** Iterate through local blocks perform ns boundary exchange
        ! ***
        do i=1,nLocalBlocks
            ! ***
            ! *** Get information about block
            ! ***
            glbBlk=local_blocks(i)
            call get_block_parameter(glbBlk, ib=ib, ie=ie, jb=jb, je=je)

            if(Neigh(inorth, glbBlk) == 0) then
                proc_north = 0
                block_north = 0
            else
                proc_north  = mDist%proc(Neigh(inorth, glbBlk))
                block_north = mDist%local_block(Neigh(inorth, glbBlk))
            end if

            if(Neigh(isouth, glbBlk) == 0) then
                proc_south = 0
                block_south = 0
            else
                proc_south  = mDist%proc(Neigh(isouth, glbBlk))
                block_south = mDist%local_block(Neigh(isouth, glbBlk))
            end if

            ! ***
            ! *** Grab south
            ! ***
            if(proc_south /= 0) then
                ARRAY(1:nx_block, 1:nghost, i) = ARRAY( &
                    1:nx_block, je-nghost+1:je, block_south)[proc_south]
            end if

            ! ***
            ! *** Send north
            ! ***
            if(proc_north /= 0) then
                ARRAY(1:nx_block, je+1:ny_block, i) = ARRAY( &
                    1:nx_block, jb:jb+nghost-1, block_north)[proc_north]
            end if
        end do

        sync all
    end subroutine caf_boundary_2d_dbl

    !***********************************************************************
    subroutine caf_boundary_2d_int(ARRAY, in_bndy, grid_loc, field_type)
        include 'mpif.h'   ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        type (bndy), intent(in) :: &
            in_bndy                 ! boundary update structure for the array

        integer (i4), intent(in) :: &
            field_type,         &! id for type of field (scalar, vector, angle)
            grid_loc             ! id for location on horizontal grid
                                 ! (center, NEcorner, Nface, Eface)

        ! !INPUT/OUTPUT PARAMETERS:
        integer (i4), intent(inout) :: &
            ARRAY(:,:,:)[*]    ! array containing horizontal slab to update

        ! ! Locals:
        integer(i4) :: i, glbBlk
        integer(i4) :: ib, jb, ie, je    ! bounds of block (non-ghost-cell) data
        integer(i4) :: proc_west, proc_east, proc_north, proc_south
        integer(i4) :: block_west, block_east, block_north, block_south

        ! ***
        ! *** Set border elements to zeroes
        ! ***
        do i=1,nLocalBlocks
            ! get block bounds information for the current block
            glbBlk=local_blocks(i)
            call get_block_parameter(glbBlk, ib=ib, ie=ie, jb=jb, je=je)

            ! clear ghost cells in each of the four borders
            ARRAY(1:nghost,  jb:je, i) = 0              ! west
            ARRAY(ie+1:ie+nghost,  jb:je, i) = 0        ! east
            ARRAY(1:ie+nghost,  1:jb-1, i) = 0          ! south
            ARRAY(1:ie+nghost,  je+1:je+nghost, i) = 0  ! north
        end do
        sync all

        ! ***
        ! *** Iterate through local blocks perform ew boundary exchange
        ! ***
        do i=1,nLocalBlocks
            ! ***
            ! *** Get information about block
            ! ***
            glbBlk=local_blocks(i)
            call get_block_parameter(glbBlk, ib=ib, ie=ie, jb=jb, je=je)

            if(Neigh(iwest, glbBlk) == 0) then
                proc_west = 0
                block_west = 0
            else
                proc_west  = mDist%proc(Neigh(iwest, glbBlk))
                block_west = mDist%local_block(Neigh(iwest, glbBlk))
            end if
            if(Neigh(ieast, glbBlk) == 0) then
                proc_east = 0
                block_east = 0
            else
                proc_east  = mDist%proc(Neigh(ieast, glbBlk))
                block_east = mDist%local_block(Neigh(ieast, glbBlk))
            end if

            ! ***
            ! *** Grab east
            ! ***
            if(proc_east /= 0) then
                ARRAY(ie+1:nx_block, jb:je, i) = ARRAY( &
                    ib:ib+nghost-1, jb:je, block_east)[proc_east]
            end if

            ! ***
            ! *** Grab west
            ! ***
            if(proc_west /= 0) then
                ARRAY(1:nghost, jb:je, i) = ARRAY( &
                    ie-nghost+1:ie, jb:je, block_west)[proc_west]
            end if
        end do

        sync all

        ! ***
        ! *** Iterate through local blocks perform ns boundary exchange
        ! ***
        do i=1,nLocalBlocks
            ! ***
            ! *** Get information about block
            ! ***
            glbBlk=local_blocks(i)
            call get_block_parameter(glbBlk, ib=ib, ie=ie, jb=jb, je=je)

            if(Neigh(inorth, glbBlk) == 0) then
                proc_north = 0
                block_north = 0
            else
                proc_north  = mDist%proc(Neigh(inorth, glbBlk))
                block_north = mDist%local_block(Neigh(inorth, glbBlk))
            end if

            if(Neigh(isouth, glbBlk) == 0) then
                proc_south = 0
                block_south = 0
            else
                proc_south  = mDist%proc(Neigh(isouth, glbBlk))
                block_south = mDist%local_block(Neigh(isouth, glbBlk))
            end if

            ! ***
            ! *** Grab south
            ! ***
            if(proc_south /= 0) then
                ARRAY(1:nx_block, 1:nghost, i) = ARRAY( &
                    1:nx_block, je-nghost+1:je, block_south)[proc_south]
            end if

            ! ***
            ! *** Send north
            ! ***
            if(proc_north /= 0) then
                ARRAY(1:nx_block, je+1:ny_block, i) = ARRAY( &
                    1:nx_block, jb:jb+nghost-1, block_north)[proc_north]
            end if
        end do

        sync all
    end subroutine caf_boundary_2d_int
end module caf_pull_boundary
