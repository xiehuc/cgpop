!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

!>
!! Routines and variables that are necessary for communicating between
!! processors.
!<
module communicate
    use kinds_mod, only: i4

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS:
    public  :: init_communicate,          &
        exit_message_environment,  &
        abort_message_environment, &
        get_num_procs

    ! !PUBLIC DATA MEMBERS:
    integer (i4), public :: &
        MPI_COMM_OCN,             &! MPI communicator for ocn comms
        mpi_dbl,                  &! MPI type for dbl_kind
        my_task,                  &! MPI task number for this task
        master_task                ! task number of master task

    integer (i4), parameter, public :: &
        mpitag_bndy_2d        = 1,    &! MPI tags for various
        mpitag_bndy_3d	    = 2,    &
        mpitag_gs             = 1000   ! communication patterns

    !***********************************************************************

    contains

    !***********************************************************************

    !>
    !! This routine sets up MPI environment and defines ocean
    !! communicator.
    !<
    subroutine init_communicate
        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        include 'mpif.h'   ! MPI Fortran include file
        integer (i4) :: ierr  ! MPI error flag

        !-------------------------------------------------------------------
        !
        !  initiate mpi environment and create communicator for internal
        !  ocean communications
        !
        !-------------------------------------------------------------------
        call MPI_INIT(ierr)
        master_task = 0
        MPI_COMM_OCN = MPI_COMM_WORLD    
        call MPI_COMM_RANK  (MPI_COMM_OCN, my_task, ierr)

        !-------------------------------------------------------------------
        !
        !  On some 64-bit machines where real_kind and dbl_kind are
        !  identical, the MPI implementation uses MPI_REAL for both.
        !  In these cases, set MPI_DBL to MPI_REAL.
        !
        !-------------------------------------------------------------------
        MPI_DBL = MPI_DOUBLE_PRECISION
    end subroutine init_communicate

    !***********************************************************************
    !>
    !! This function returns the number of processor assigned to
    !! MPI_COMM_OCN
    !<
    function get_num_procs()
        integer (i4) :: get_num_procs

        !-----------------------------------------------------------------------
        !
        !  local variables
        !
        !-----------------------------------------------------------------------
        integer (i4) :: ierr

        call MPI_COMM_SIZE(MPI_COMM_OCN, get_num_procs, ierr)
    end function get_num_procs

    !***********************************************************************

    !>
    !! This routine exits the message environment properly when model
    !! stops.
    !!
    !! @param ierr   MPI error flag
    !<
    subroutine exit_message_environment(ierr)
        include 'mpif.h'   ! MPI Fortran include file

        ! !OUTPUT PARAMETERS:
        integer (i4), intent(out) :: ierr

        call MPI_FINALIZE(ierr)
    end subroutine exit_message_environment

    !***********************************************************************

    !>
    !! This routine aborts the message environment when model stops.
    !! It will attempt to abort the entire MPI COMM WORLD.
    !!
    !! @param ierr   MPI error flag
    !<
    subroutine abort_message_environment(ierr)
        include 'mpif.h'   ! MPI Fortran include file

        ! !OUTPUT PARAMETERS:
        integer (i4), intent(out) :: ierr

        call MPI_BARRIER(MPI_COMM_OCN, ierr)
        call MPI_ABORT(MPI_COMM_WORLD, ierr)
        call MPI_FINALIZE(ierr)
    end subroutine abort_message_environment

end module communicate
