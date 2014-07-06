!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! Routines to initialize and perform the boundary exchange operation by
!! pulling data in an unbuffered fashion.
!<
module caf_multi_pull_halo
    use kinds_mod, only: i4, r8

    implicit none
    private
    save

    include 'mpif.h'

    ! ==================
    ! array information:
    ! ==================
    integer(i4), allocatable :: &
        halo2grab(:)    ! Identifies where remote data is located for each
                        ! point in the halo that will be pulled.

    integer(i4), allocatable :: &
        linear2Proc(:)  ! Identifies which proccess owns the data that is

    integer(i4) :: &
        linearMaximum, &    ! Identifies the size of the largest array that
                       &    ! could be passed to the update routine.  (This
                       &    ! value will be globally replicated across all
                       &    ! procs).
        nTotalCA,      &    ! Identifies the size of the array that will be
                       &    ! passed to the update routine.  (This value may
                       &    ! be unique on all procs).
        nActiveCA           ! Identifies how many points in the passed array
                            ! are active ocean points (non halo points).

    ! ==========================
    ! Public subroutines:
    ! ==========================
    public :: caf_multi_pull_update_1d_dbl
    public :: caf_multi_pull_update_1d_int
    public :: caf_multi_pull_init

  contains

    !***********************************************************************
    !>
    !! This subroutine performs the boundary exchange operation in order to
    !! update the halo elements in the passed array.
    !!
    !! @param handle
    !! @param array      vector to update
    !<
    subroutine caf_multi_pull_update_1d_dbl(handle,array)
        ! !INPUT PARAMETERS:
        integer(i4), intent(in) :: handle
        real(r8), intent(inout) :: array(:) 

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer(i4) :: i ! dummy counter
        real(r8), allocatable, save :: arrayCA(:)[:] ! co-array version of array

        allocate(arrayCA(linearMaximum)[*])
        arrayCA(1:ubound(array,1)) = array(:)
        sync all

        !-----------------------------------------------------------------------
        !
        !  iterate through halo elements, grabbing fresh values from remote
        !  images as needed.
        !
        !-----------------------------------------------------------------------
        do i=nActiveCA+1,nTotalCA
            array(i) = arrayCA(halo2grab(i))[linear2Proc(i)+1]
        enddo

        deallocate(arrayCA)
    end subroutine caf_multi_pull_update_1d_dbl


    !***********************************************************************
    !>
    !! This subroutine performs the boundary exchange operation in order to
    !! update the halo elements in the passed array.
    !!
    !! @param handle
    !! @param array      vector to update
    !<
    subroutine caf_multi_pull_update_1d_int(handle,array)
        ! !INPUT PARAMETERS:
        integer(i4), intent(in)    :: handle
        integer(i4), intent(inout) :: array(:)

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer(i4) :: i ! dummy counter
        integer(i4), allocatable, save :: arrayCA(:)[:]

        !-------------------------------------------------------------------
        !
        !  in order to perform the boundary exchange operation, the array
        !  parameter must be a co-array.  It would be ideal to refactor the
        !  code so this is the case coming in.  For now, however, I'll simply
        !  make a copy within the function.
        !
        !-------------------------------------------------------------------
        allocate(arrayCA(linearMaximum)[*])
        arrayCA(1:ubound(array,1)) = array(:)
        sync all

        !-------------------------------------------------------------------
        !
        !  iterate through halo elements, grabbing fresh values from remote
        !  images as needed.
        !
        !-------------------------------------------------------------------
        do i=nActiveCA+1,nTotalCA
            array(i) = arrayCA(halo2grab(i))[linear2Proc(i)+1]
        enddo

        deallocate(arrayCA)

        sync all
    end subroutine caf_multi_pull_update_1d_int


    !***********************************************************************
    !>
    !! This function initilizes scheduling metadata needed for the update
    !! routines (the update routines perform the boundary exchange operation.)
    !!
    !! @param COMM          MPI communicator
    !! @param maxlinear     Maximum size of vectors
    !! @param nTotal        Total elements in local vector
    !! @param nActive       Number of active (non-halo) elements
    !! @param LinearGdof    Maps vector indices to GDOFs
    !! @param LinearProc    Maps halo indices to processors
    !<
    function caf_multi_pull_init( &
        COMM,maxlinear, nTotal,nActive,LinearGdof,LinearProc) result(handle)

        ! !INPUT PARAMETERS:
        integer(i4), intent(in) :: &
            COMM,                  &
            maxlinear,             &
            nTotal,                &
            nActive,               &
            LinearGdof(maxlinear), &
            LinearProc(maxlinear)
        integer(i4) :: handle

        !-----------------------------------------------------------------------
        !
        !  local variables
        !
        !-----------------------------------------------------------------------
        integer(i4), save, allocatable :: &
            LinearGdofCA(:)[:]  ! Co-array copy of LinearGdof
        integer(i4), allocatable :: &
            copy(:)        ! For temp local copies of remote LinearGdfCA arrays
        integer(i4) :: i, j     ! Dummy loop-index variables
        integer(i4) :: idx      ! Vector index of a remote gdof element

        linearMaximum = maxlinear
        nTotalCA = nTotal
        nActiveCA = nActive

        allocate(copy(maxlinear))
        allocate(linear2Proc(maxlinear))

        allocate(LinearGdofCA(maxlinear)[*])
        LinearGdofCA(:) = LinearGdof(:)
        sync all

        linear2Proc(:) = LinearProc(:)

        allocate(halo2grab(nActive+1:nTotal))

        sync all

        !-----------------------------------------------------------------------
        !
        !  iterate through halo points, constructing halo2grab array (halo2grab
        !  maps elements in this processor's halo to indices in a remote
        !  processor's vector)
        !
        !-----------------------------------------------------------------------
        do_i: do i=nActive+1,nTotal
            idx = -1

            copy(:) = LinearGdofCA(:)[LinearProc(i)+1]

            do_j: do j=1,maxlinear
                if(copy(j) == LinearGdofCA(i)) then
                    idx = j
                    exit do_j
                end if
            end do do_j

            halo2grab(i) = idx
        enddo do_i

        sync all

        deallocate(LinearGdofCA)
    end function caf_multi_pull_init
end module caf_multi_pull_halo
