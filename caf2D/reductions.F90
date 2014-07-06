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
    ! !USES:
    use kinds_mod, only: i4, r4, r8, log_kind
    use simple_type, only: distrb
    use communicate, only: my_task, MPI_DBL, MPI_COMM_OCN
    use constants, only: c0, field_loc_necorner, field_loc_nface
    use simple_blocks, only: nblocks_y, nblocks_tot, get_block_parameter
    use domain_size, only: nx_global

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS:

    public :: global_sum
    public :: global_maxval
    public :: global_minval

    !-----------------------------------------------------------------------
    !
    !  generic interfaces for module procedures
    !
    !-----------------------------------------------------------------------
    !>
    !! Computes the global sum of either a distributed 2-d array
    !! which corresponds to the _physical domain_, 1-d vectors or a scalar.
    !<
    interface global_sum
    module procedure global_sum_dbl
    module procedure global_sum_vec_dbl
    module procedure global_sum_scalar_int
    module procedure global_sum_scalar_dbl
    module procedure global_sum_scalar_real
    end interface 

    !>
    !! Computes the maximum value of a scalar over all tasks.
    !<
    interface global_maxval
    module procedure global_maxval_scalar_int
    module procedure global_maxval_scalar_dbl
    module procedure global_maxval_scalar_real
    end interface

    !>
    !! Computes the minimum value of a scalar over all tasks.
    !<
    interface global_minval
    module procedure global_minval_scalar_int
    module procedure global_minval_scalar_dbl
    module procedure global_minval_scalar_real
    end interface

    !-----------------------------------------------------------------------
    !
    !  module variables
    !
    !-----------------------------------------------------------------------

  contains

    !***********************************************************************
    function global_maxval_scalar_int (local_scalar)
        !-----------------------------------------------------------------------
        !
        !  this function returns the maximum scalar value across processors
        !
        !-----------------------------------------------------------------------
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
        !  this function returns the maximum scalar value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        real (r8), intent(inout) :: &
            local_scalar                ! local scalar to be compared

        real (r8) :: &
            global_maxval_scalar_dbl   ! resulting global max

        integer (i4) :: ierr ! MPI error flag

        call MPI_ALLREDUCE(local_scalar, global_maxval_scalar_dbl, 1, &
            mpi_dbl, MPI_MAX, MPI_COMM_OCN, ierr)
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
        !-------------------------------------------------------------------
        !
        !  this function returns the minimum scalar value across processors
        !
        !-------------------------------------------------------------------
        include 'mpif.h'  ! MPI Fortran include file

        real (r8), intent(inout) :: &
            local_scalar                ! local scalar to be compared

        real (r8) :: &
            global_minval_scalar_dbl   ! resulting global min

        integer (i4) :: ierr ! MPI error flag

        call MPI_ALLREDUCE(local_scalar, global_minval_scalar_dbl, 1, &
            mpi_dbl, MPI_MIN, MPI_COMM_OCN, ierr)
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

        if (num_images() > 1) then
            if (my_task < num_images()) then
                call MPI_ALLREDUCE(local_scalar, global_sum_scalar_dbl, 1, &
                    mpi_dbl, MPI_SUM, MPI_COMM_OCN, ierr)
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

        if (num_images() > 1) then
            if (my_task < num_images()) then
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


        if (num_images() > 1) then
            if (my_task < num_images()) then
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
    !>
    !! computes the global sum of the _physical domain_ of a 2-d
    !! array.
    !! 
    !! REMARKS:
    !! This is actually the specific interface for the generic global_sum
    !! function corresponding to double precision arrays.  The generic
    !! interface is identical but will handle real and integer 2-d slabs
    !! and real, integer, and double precision scalars.
    !<
    function global_sum_dbl(X, dist, field_loc, MASK)
        ! !USES:
        include 'mpif.h'  ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        real (r8), dimension(:,:,:), intent(in) :: &
            X                    ! array to be summed
        
        type (distrb), intent(in) :: &
            dist                 ! block distribution for array X

        integer (i4), intent(in) :: &
            field_loc            ! location of field on staggered grid

        real (r8), dimension(size(X,dim=1), &
            size(X,dim=2), &
            size(X,dim=3)), intent(in), optional :: &
                MASK                 ! real multiplicative mask

        ! !OUTPUT PARAMETERS:
        real (r8) :: global_sum_dbl       ! resulting global sum

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        real (r8) :: local_sum           ! sum of all local blocks

        integer (i4) :: &
            i,j,n,             &! local counters
            ib,ie,jb,je,       &! beg,end of physical domain
            bid,               &! block location
            ierr                ! MPI error flag

        !-------------------------------------------------------------------
        !
        !  use this code for sums that are not reproducible for best
        !  performance
        !
        !-------------------------------------------------------------------
        local_sum = c0

        do n=1,nblocks_tot
            if (dist%proc(n) == my_task+1) then
                bid = dist%local_block(n)
                call get_block_parameter(n,ib=ib,ie=ie,jb=jb,je=je)
                if (present(MASK)) then
                    do j=jb,je
                        do i=ib,ie
                            local_sum = local_sum + X(i,j,bid)*MASK(i,j,bid)
                        end do
                    end do
                else ! no mask
                    do j=jb,je
                        do i=ib,ie
                            local_sum = local_sum + X(i,j,bid)
                        end do
                    end do
                endif
            endif
        end do !block loop

        if (dist%nprocs > 1) then
            if (my_task < dist%nprocs) then
                call MPI_ALLREDUCE(local_sum, global_sum_dbl, 1, &
                mpi_dbl, MPI_SUM, dist%communicator, ierr)
            else
                global_sum_dbl = c0
            endif
        else
            global_sum_dbl = local_sum
        endif
    end function global_sum_dbl

    !***********************************************************************
    !>
    !! computes a vector global sum of the _physical domain_ of a 2-d
    !! array.
    !!
    !! REMARKS:
    !! This is actually the specific interface for the generic global_sum
    !! function corresponding to double precision arrays.  The generic
    !! interface is identical but will handle real and integer 2-d slabs
    !! and real, integer, and double precision scalars.
    !<
    function global_sum_vec_dbl(X, dist, field_loc,MASK)
        ! !USES:
        include 'mpif.h'  ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        real (r8), dimension(:,:,:,:), intent(in) :: &
            X                    ! array to be summed

        type (distrb), intent(in) :: &
        dist                 ! block distribution for array X

        integer (i4), intent(in) :: &
        field_loc            ! location of field on staggered grid

        real (r8), dimension(size(X,dim=1), &
        size(X,dim=2), &
        size(X,dim=4)), intent(in), optional :: &
        MASK                 ! real multiplicative mask

        ! !OUTPUT PARAMETERS:
        real (r8) :: &
        global_sum_vec_dbl(size(X,dim=3))       ! resulting global sum

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        real (r8) :: local_sum(size(X,dim=3))   ! sum of all local blocks

        integer (i4) :: &
            i,j,n,             &! local counters
            ib,ie,jb,je,       &! beg,end of physical domain
            iv,                &! vector index
            bid,               &! block location
            vlen,              &! length of vector reduction
            ierr                ! MPI error flag

        !-------------------------------------------------------------------
        !
        ! use this code for sums that are not reproducible for best performance
        !
        !-------------------------------------------------------------------
        vlen = SIZE(X,dim=3)
        local_sum = c0

        do n=1,nblocks_tot
            if (dist%proc(n) == my_task+1) then
                bid = dist%local_block(n)
                call get_block_parameter(n,ib=ib,ie=ie,jb=jb,je=je)
                if (present(MASK)) then
                    do iv=1,vlen
                        do j=jb,je
                            do i=ib,ie
                                local_sum(iv) = &
                                    local_sum(iv) + X(i,j,iv,bid)*MASK(i,j,bid)
                            end do
                        end do
                    enddo
                else ! no mask
                    do iv=1,vlen
                        do j=jb,je
                            do i=ib,ie
                                local_sum(iv) = local_sum(iv) + X(i,j,iv,bid)
                            end do
                        end do
                    enddo
                endif
            endif
        end do !block loop

        if (dist%nprocs > 1) then
            if (my_task < dist%nprocs) then
                call MPI_ALLREDUCE(local_sum, global_sum_vec_dbl, vlen, &
                    mpi_dbl, MPI_SUM, dist%communicator, ierr)
            else
                global_sum_vec_dbl = c0
            endif
        else
            global_sum_vec_dbl = local_sum
        endif
    end function global_sum_vec_dbl
end module reductions
