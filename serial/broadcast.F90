!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

! !DESCRIPTION:
!  This module contains all the broadcast routines.  This
!  particular version contains serial versions of these routines
!  which typically perform no operations since there is no need
!  to broadcast what is already known.
module broadcast
    use kinds_mod, only : i4, r4, r8, log_kind
    use communicate, only :

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS:
    public  :: broadcast_scalar, broadcast_array

    !-----------------------------------------------------------------------
    !
    !  generic interfaces for module procedures
    !
    !-----------------------------------------------------------------------

    interface broadcast_scalar
        module procedure broadcast_scalar_dbl,  &
            broadcast_scalar_real, &
            broadcast_scalar_int,  &
            broadcast_scalar_log,  &
            broadcast_scalar_char
    end interface

    interface broadcast_array
        module procedure broadcast_array_dbl_1d,  &
            broadcast_array_real_1d, &
            broadcast_array_int_1d,  &
            broadcast_array_log_1d,  &
            broadcast_array_dbl_2d,  &
            broadcast_array_real_2d, &
            broadcast_array_int_2d,  &
            broadcast_array_log_2d,  &
            broadcast_array_dbl_3d,  &
            broadcast_array_real_3d, &
            broadcast_array_int_3d,  &
            broadcast_array_log_3d
    end interface

    contains

    subroutine broadcast_scalar_dbl(scalar, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a scalar dbl variable from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_scalar interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), intent(inout) :: &
        scalar               ! scalar to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------

    end subroutine broadcast_scalar_dbl

    subroutine broadcast_scalar_real(scalar, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a scalar real variable from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_scalar interface.
        !
        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        real (r4), intent(inout) :: &
        scalar               ! scalar to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------

    end subroutine broadcast_scalar_real

    subroutine broadcast_scalar_int(scalar, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a scalar integer variable from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_scalar interface.
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        integer (i4), intent(inout) :: &
        scalar                ! scalar to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------

    end subroutine broadcast_scalar_int


    subroutine broadcast_scalar_log(scalar, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a scalar logical variable from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_scalar interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        logical (log_kind), intent(inout) :: &
        scalar               ! scalar to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------

    end subroutine broadcast_scalar_log

    subroutine broadcast_scalar_char(scalar, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a scalar character variable from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_scalar interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        character (*), intent(inout) :: &
        scalar               ! scalar to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------
    end subroutine broadcast_scalar_char


    subroutine broadcast_array_dbl_1d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a vector dbl variable from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe           ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), dimension(:), intent(inout) :: &
        array             ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------
    end subroutine broadcast_array_dbl_1d


    subroutine broadcast_array_real_1d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a real vector from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        real (r4), dimension(:), intent(inout) :: &
        array                ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------
    end subroutine broadcast_array_real_1d


    subroutine broadcast_array_int_1d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts an integer vector from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        integer (i4), dimension(:), intent(inout) :: &
        array              ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------
    end subroutine broadcast_array_int_1d


    subroutine broadcast_array_log_1d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a logical vector from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        logical (log_kind), dimension(:), intent(inout) :: &
        array                ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------
    end subroutine broadcast_array_log_1d

    subroutine broadcast_array_dbl_2d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a dbl 2d array from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe           ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), dimension(:,:), intent(inout) :: &
        array             ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------
    end subroutine broadcast_array_dbl_2d


    subroutine broadcast_array_real_2d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a real 2d array from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        real (r4), dimension(:,:), intent(inout) :: &
        array                ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------
    end subroutine broadcast_array_real_2d


    subroutine broadcast_array_int_2d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a 2d integer array from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        integer (i4), dimension(:,:), intent(inout) :: &
        array              ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------

    end subroutine broadcast_array_int_2d

    subroutine broadcast_array_log_2d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a logical 2d array from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        logical (log_kind), dimension(:,:), intent(inout) :: &
        array                ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------

    end subroutine broadcast_array_log_2d


    subroutine broadcast_array_dbl_3d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a double 3d array from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe           ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), dimension(:,:,:), intent(inout) :: &
        array             ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------

    end subroutine broadcast_array_dbl_3d


    subroutine broadcast_array_real_3d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a real 3d array from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
            integer (i4), intent(in) :: &
            root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
            real (r4), dimension(:,:,:), intent(inout) :: &
            array                ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------

    end subroutine broadcast_array_real_3d


    subroutine broadcast_array_int_3d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts an integer 3d array from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
            root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        integer (i4), dimension(:,:,:), intent(inout) :: &
            array              ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------
    end subroutine broadcast_array_int_3d

    subroutine broadcast_array_log_3d(array, root_pe)
        ! !DESCRIPTION:
        !  Broadcasts a logical 3d array from one processor (root_pe)
        !  to all other processors. This is a specific instance of the generic
        !  broadcast\_array interface.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: &
        root_pe              ! processor number to broadcast from

        ! !INPUT/OUTPUT PARAMETERS:
        logical (log_kind), dimension(:,:,:), intent(inout) :: &
        array                ! array to be broadcast

        !-----------------------------------------------------------------------
        !
        !  for serial codes, nothing is required
        !
        !-----------------------------------------------------------------------

    end subroutine broadcast_array_log_3d
end module broadcast
