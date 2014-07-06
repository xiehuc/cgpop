!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! This module contains all the routines for performing global
!! reductions like global sums, minvals, maxvals, etc.
!<
module reductions
    use kinds_mod,   only: i4, r4, r8
    use communicate, only: my_task, MPI_COMM_OCN, get_num_procs
    use constants,   only: c0

    implicit none

    private

    public :: global_minval, global_maxval, global_sum

    !>
    !! Computes the global sum of either a distributed 
    !! 1-d vector or a scalar
    !<
    interface global_sum
        module procedure global_sum_scalar_dbl, &
            global_sum_scalar_int,  &
            global_sum_scalar_real, &
            global_sum_vector_dbl,  &
            global_sum_vector_int
    end interface

    !>
    !! Computes the minimum value of a scalar over all tasks.
    !<
    interface global_minval
        module procedure global_minval_scalar_dbl,    &
            global_minval_scalar_real,   &
            global_minval_scalar_int
    end interface

    !>
    !! Computes the maximum value of a scalar over all tasks.
    !<
    interface global_maxval
        module procedure global_maxval_scalar_int, &
            global_maxval_scalar_dbl, &
            global_maxval_scalar_real
    end interface

    contains
    
    !***********************************************************************
    function global_maxval_scalar_int (local_scalar)
        !-------------------------------------------------------------------
        !
        !  this function returns the maximum scalar value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        integer (i4), intent(inout) :: &
            local_scalar                ! local scalar to be compared

        integer (i4) :: &
            global_maxval_scalar_int   ! resulting global max

        integer (i4) :: ierr ! MPI error flag

        call MPI_ALLREDUCE(local_scalar, global_maxval_scalar_int, 1, &
            MPI_INTEGER, MPI_MAX, MPI_COMM_OCN, ierr)
    end function global_maxval_scalar_int

    !***********************************************************************
    function global_maxval_scalar_dbl (local_scalar)
        !-------------------------------------------------------------------
        !
        ! this function returns the maximum scalar value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        real (r8), intent(inout) :: &
            local_scalar                ! local scalar to be compared

        real (r8) :: &
            global_maxval_scalar_dbl   ! resulting global max

        integer (i4) :: ierr ! MPI error flag

        call MPI_ALLREDUCE(local_scalar, global_maxval_scalar_dbl, 1, &
            MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_OCN, ierr)
    end function global_maxval_scalar_dbl

    !***********************************************************************
    function global_maxval_scalar_real (local_scalar)
        !-------------------------------------------------------------------
        !
        !  this function returns the maximum scalar value across processors
        !
        !-------------------------------------------------------------------

        include 'mpif.h'  ! MPI Fortran include file

        real (r4), intent(inout) :: &
            local_scalar                ! local scalar to be compared

        real (r4) :: &
            global_maxval_scalar_real   ! resulting global max

        integer (i4) :: ierr ! MPI error flag

        call MPI_ALLREDUCE(local_scalar, global_maxval_scalar_real, 1, &
            MPI_REAL, MPI_MAX, MPI_COMM_OCN, ierr)
    end function global_maxval_scalar_real

    !***********************************************************************
    function global_minval_scalar_dbl (local_scalar)
        !-----------------------------------------------------------------------
        !
        !  this function returns the minimum scalar value across processors
        !
        !-----------------------------------------------------------------------

        include 'mpif.h'  ! MPI Fortran include file

        real (r8), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        real (r8) :: &
        global_minval_scalar_dbl   ! resulting global min

        integer (i4) :: ierr ! MPI error flag


        call MPI_ALLREDUCE(local_scalar, global_minval_scalar_dbl, 1, &
            MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_OCN, ierr)
    end function global_minval_scalar_dbl

    !***********************************************************************
    function global_minval_scalar_real (local_scalar)
        !-------------------------------------------------------------------
        !
        !  this function returns the minimum scalar value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        real (r4), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        real (r4) :: &
        global_minval_scalar_real   ! resulting global min

        integer (i4) :: ierr ! MPI error flag

        call MPI_ALLREDUCE(local_scalar, global_minval_scalar_real, 1, &
            MPI_REAL, MPI_MIN, MPI_COMM_OCN, ierr)
    end function global_minval_scalar_real

    !***********************************************************************
    function global_minval_scalar_int (local_scalar)
        !-------------------------------------------------------------------
        !
        !  this function returns the minimum scalar value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        integer (i4), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        integer (i4) :: &
        global_minval_scalar_int   ! resulting global min

        integer (i4) :: ierr ! MPI error flag

        call MPI_ALLREDUCE(local_scalar, global_minval_scalar_int, 1, &
            MPI_INTEGER, MPI_MIN, MPI_COMM_OCN, ierr)
    end function global_minval_scalar_int

    !***********************************************************************
    function global_sum_scalar_dbl(local_scalar)
        !-------------------------------------------------------------------
        !
        !  this function returns the sum of scalar value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        real (r8), intent(inout) :: &
            local_scalar                ! local scalar to be compared

        real (r8) :: &
            global_sum_scalar_dbl   ! resulting global sum

        integer (i4) :: ierr ! MPI error flag

        if ( get_num_procs() > 1) then
            if (my_task < get_num_procs()) then
                call MPI_ALLREDUCE(local_scalar, global_sum_scalar_dbl, 1, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_OCN, ierr)
            else
                global_sum_scalar_dbl = c0
            endif
        else
            global_sum_scalar_dbl = local_scalar
        endif
    end function global_sum_scalar_dbl

    !***********************************************************************
    function global_sum_scalar_real(local_scalar)
        !-------------------------------------------------------------------
        !
        !  this function returns the sum of scalar value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        real (r4), intent(inout) :: &
        local_scalar                ! local scalar to be compared

        real (r4) :: &
            global_sum_scalar_real   ! resulting global sum

        integer (i4) :: ierr ! MPI error flag

        if (get_num_procs() > 1) then
            if (my_task < get_num_procs()) then
                call MPI_ALLREDUCE(local_scalar, global_sum_scalar_real, 1, &
                    MPI_REAL, MPI_SUM, MPI_COMM_OCN, ierr)
            else
                global_sum_scalar_real = c0
            endif
        else
            global_sum_scalar_real = local_scalar
        endif
    end function global_sum_scalar_real

    !***********************************************************************
    function global_sum_scalar_int(local_scalar)
        !-------------------------------------------------------------------
        !
        !  this function returns the sum of scalar value across processors
        !
        !-------------------------------------------------------------------

        include 'mpif.h'  ! MPI Fortran include file

        integer (i4), intent(inout) :: &
            local_scalar                ! local scalar to be compared

        integer (i4) :: &
            global_sum_scalar_int   ! resulting global sum

        integer (i4) :: ierr ! MPI error flag

        if (get_num_procs() > 1) then
            if (my_task < get_num_procs()) then
                call MPI_ALLREDUCE(local_scalar, global_sum_scalar_int, 1, &
                    MPI_INTEGER, MPI_SUM, MPI_COMM_OCN, ierr)
            else
                global_sum_scalar_int = 0
            endif
        else
            global_sum_scalar_int = local_scalar
        endif
    end function global_sum_scalar_int

    !***********************************************************************
    function global_sum_vector_dbl(local_vector)
        !-------------------------------------------------------------------
        !
        !  this function returns the sum of vector value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        real (r8), intent(inout) :: &
            local_vector(:)                ! local vector to be compared

        real (r8), dimension(size(local_vector,dim=1)) :: &
            global_sum_vector_dbl   ! resulting global sum

        integer (i4) :: ierr ! MPI error flag
        integer (i4) :: vlen ! length of Vector sum

        !-----------------------------------------------------------------------

        vlen = size(local_vector,dim=1)
        if (get_num_procs() > 1) then
            if (my_task < get_num_procs()) then
                call MPI_ALLREDUCE(local_vector, global_sum_vector_dbl, vlen, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_OCN, ierr)
            else
                global_sum_vector_dbl = c0
            endif
        else
            global_sum_vector_dbl = local_vector
        endif
    end function global_sum_vector_dbl

    !***********************************************************************
    function global_sum_vector_int(local_vector)
        !-------------------------------------------------------------------
        !
        !  this function returns the sum of vector value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        integer (i4), intent(inout) :: &
            local_vector(:)                ! local vector to be compared

        integer(i4), dimension(size(local_vector,dim=1)) :: &
            global_sum_vector_int   ! resulting global sum

        integer (i4) :: ierr ! MPI error flag
        integer (i4) :: vlen ! length of Vector sum

        !-----------------------------------------------------------------------

        vlen = size(local_vector,dim=1)
        if (get_num_procs() > 1) then
            if (my_task < get_num_procs()) then
                call MPI_ALLREDUCE(local_vector, global_sum_vector_int, vlen, &
                    mpi_integer, MPI_SUM, MPI_COMM_OCN, ierr)
            else
                global_sum_vector_int = 0
            endif
        else
            global_sum_vector_int = local_vector
        endif
    end function global_sum_vector_int
end module reductions
