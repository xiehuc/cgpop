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
!! boundaries using MPI calls when using the 2D data structure implementation.
!<
module boundary
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
    use mpi2s_boundary, only: mpi2s_create_boundary, mpi2s_destroy_boundary, &
        mpi2s_boundary_2d, bndy

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS:
    public :: create_boundary,  &
        destroy_boundary,       &
        update_ghost_cells

    !>
    !! Routine to perform the boundary exchange operation when local data
    !! is stored in a 2-dimensional array.
    !<
    interface update_ghost_cells  ! generic interface
        module procedure boundary_2d_dbl,  boundary_2d_int
    end interface

  contains

    !***********************************************************************
    !>
    !! This routine creates a boundary type with info necessary for
    !! performing a boundary (ghost cell) update based on the input block
    !! distribution.
    !<
    subroutine create_boundary(newbndy, dist, blocks, ns_bndy_type, &
        ew_bndy_type, nx_global, ny_global)

        ! !INPUT PARAMETERS:
        type (distrb), intent(in) :: &
            dist       ! distribution of blocks across procs

        integer(i4) :: blocks(:)

        character (*), intent(in) :: &
            ns_bndy_type,             &! type of boundary to use in ns dir
            ew_bndy_type               ! type of boundary to use in ew dir

        integer (i4), intent(in) :: &
            nx_global, ny_global       ! global extents of domain

        ! !OUTPUT PARAMETERS:
        type (bndy), intent(out) :: &
            newbndy    ! a new boundary type with info for updates


        select case (boundary_exchange_algorithm)
            case (ALG_MPI2S_2D)
                call mpi2s_create_boundary( &
                    newbndy, dist, blocks, ns_bndy_type, ew_bndy_type, &
                    nx_global, ny_global)
        end select
    end subroutine create_boundary

    !***********************************************************************
    !>
    !! This routine destroys a boundary by deallocating all memory
    !! associated with the boundary and nullifying pointers.
    !<
    subroutine destroy_boundary(in_bndy)
        ! !INPUT/OUTPUT PARAMETERS:
        type (bndy), intent(inout) :: &
            in_bndy          ! boundary structure to be destroyed


        select case (boundary_exchange_algorithm)
            case (ALG_MPI2S_2D)
                call mpi2s_destroy_boundary(in_bndy)
        end select
    end subroutine destroy_boundary

    !***********************************************************************
    !>
    !! This routine updates ghost cells for an input array and is a
    !! member of a group of routines under the generic interface
    !! update\_ghost\_cells.  This routine is the specific interface
    !! for 2d horizontal arrays of double precision.
    !<
    subroutine boundary_2d_dbl(ARRAY, in_bndy, grid_loc, field_type)
        ! !USER:
        include 'mpif.h'   ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        type (bndy), intent(in) :: &
            in_bndy                 ! boundary update structure for the array

        integer (i4), intent(in) :: &
            field_type,             & ! type of field (scalar, vector, angle)
            grid_loc                  ! id for location on horizontal grid
                                      !  (center, NEcorner, Nface, Eface)

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), intent(inout) :: &
            ARRAY(:,:,:)    ! array containing horizontal slab to update

        
        select case (boundary_exchange_algorithm)
            case (ALG_MPI2S_2D)
                call mpi2s_boundary_2d(ARRAY, in_bndy, grid_loc, field_type)
        end select
    end subroutine boundary_2d_dbl

    !***********************************************************************
    !>
    !! This routine updates ghost cells for an input array and is a
    !! member of a group of routines under the generic interface
    !! update\_ghost\_cells.  This routine is the specific interface
    !! for 2d horizontal arrays of double precision.
    !<
    subroutine boundary_2d_int(ARRAY, in_bndy, grid_loc, field_type)
        ! !USER:

        include 'mpif.h'   ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        type (bndy), intent(in) :: &
        in_bndy                 ! boundary update structure for the array

        integer (i4), intent(in) :: &
            field_type,        & ! id for type of field (scalar, vector, angle)
            grid_loc             ! id for location on horizontal grid
                                 !  (center, NEcorner, Nface, Eface)

        ! !INPUT/OUTPUT PARAMETERS:
        integer (i4), intent(inout) :: &
            ARRAY(:,:,:)    ! array containing horizontal slab to update


        select case (boundary_exchange_algorithm)
            case (ALG_MPI2S_2D)
                call mpi2s_boundary_2d(ARRAY, in_bndy, grid_loc, field_type)
        end select
    end subroutine boundary_2d_int
end module boundary
