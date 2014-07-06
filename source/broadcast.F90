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
!! Broadcast routines.  This particular version contains MPI versions of these
!! routines.
!<
module broadcast
    use kinds_mod, only: i4, r4, r8, log_kind, char_len
    use communicate, only: MPI_COMM_OCN, MPI_DBL

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS:

    public  :: broadcast_scalar

    !-----------------------------------------------------------------------
    !
    !  generic interfaces for module procedures
    !
    !-----------------------------------------------------------------------

    !>
    !! Broadcasts a scalar variable from one processor (root_pe)
    !! to all other processors.
    !<
    interface broadcast_scalar
        module procedure broadcast_scalar_dbl,  &
            broadcast_scalar_real, &
            broadcast_scalar_int,  &
            broadcast_scalar_log,  &
            broadcast_scalar_char
    end interface

    !***********************************************************************

    contains

    !***********************************************************************

    !>
    !! Broadcasts a scalar dbl variable from one processor (root_pe)
    !! to all other processors. This is a specific instance of the generic
    !! broadcast\_scalar interface.
    !<
    subroutine broadcast_scalar_dbl(scalar, root_pe)
        include 'mpif.h'  ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
            root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), intent(inout) :: &
        scalar               ! scalar to be broadcast

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer (i4) :: ierr  ! local MPI error flag

        call MPI_BCAST(scalar, 1, MPI_DBL, root_pe, MPI_COMM_OCN, ierr)
    end subroutine broadcast_scalar_dbl

    !***********************************************************************
    !>
    !! Broadcasts a scalar real variable from one processor (root_pe)
    !! to all other processors. This is a specific instance of the generic
    !! broadcast\_scalar interface.
    !<
    subroutine broadcast_scalar_real(scalar, root_pe)
        include 'mpif.h'  ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        real (r4), intent(inout) :: &
        scalar               ! scalar to be broadcast

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer (i4) :: ierr  ! local MPI error flag

        call MPI_BCAST(scalar, 1, MPI_REAL, root_pe, MPI_COMM_OCN, ierr)
    end subroutine broadcast_scalar_real

    !***********************************************************************
    !>
    !! Broadcasts a scalar integer variable from one processor (root_pe)
    !! to all other processors. This is a specific instance of the generic
    !! broadcast\_scalar interface.
    !<
    subroutine broadcast_scalar_int(scalar, root_pe)
        include 'mpif.h'  ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        integer (i4), intent(inout) :: &
        scalar                ! scalar to be broadcast

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer (i4) :: ierr  ! local MPI error flag

        call MPI_BCAST(scalar, 1, MPI_INTEGER, root_pe, MPI_COMM_OCN,ierr)
    end subroutine broadcast_scalar_int

    !***********************************************************************
    !>
    !! Broadcasts a scalar logical variable from one processor (root_pe)
    !! to all other processors. This is a specific instance of the generic
    !! broadcast\_scalar interface.
    !<
    subroutine broadcast_scalar_log(scalar, root_pe)
        include 'mpif.h'  ! MPI Fortran include file
        
        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
            root_pe              ! processor number to broadcast from
        
        ! !INPUT/OUTPUT PARAMETERS:
        logical (log_kind), intent(inout) :: &
            scalar               ! scalar to be broadcast
        
        !-----------------------------------------------------------------------
        !
        !  local variables
        !
        !-----------------------------------------------------------------------
        integer (i4) :: &
            itmp,               &! local temporary
            ierr                 ! MPI error flag

        if (scalar) then
            itmp = 1
        else
            itmp = 0
        endif

        call MPI_BCAST(itmp, 1, MPI_INTEGER, root_pe, MPI_COMM_OCN, ierr)

        if (itmp == 1) then
            scalar = .true.
        else
            scalar = .false.
        endif
    end subroutine broadcast_scalar_log

    !***********************************************************************
    !>
    !! Broadcasts a scalar character variable from one processor (root_pe)
    !! to all other processors. This is a specific instance of the generic
    !! broadcast\_scalar interface.
    !<
    subroutine broadcast_scalar_char(scalar, root_pe)
        include 'mpif.h'  ! MPI Fortran include file

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
            root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        character (*), intent(inout) :: &
            scalar               ! scalar to be broadcast

        !--------------------------------------------------------------------
        !
        !  local variables
        !
        !--------------------------------------------------------------------
        integer (i4) :: &
            clength,            &! length of character
            ierr                 ! MPI error flag

        clength = len(scalar)
        call MPI_BCAST( &
            scalar, clength, MPI_CHARACTER, root_pe, MPI_COMM_OCN, ierr)
    end subroutine broadcast_scalar_char
end module broadcast

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
