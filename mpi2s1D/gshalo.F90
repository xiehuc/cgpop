!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! Proxy module for one of the other boundary exchange modules.  The value
!! of boundary_exchange_algorithm will determine which boundary exchange
!! algorithm to use.
!<
module gshalo
    use mpi2s_gshalo
    use constants, only: ALG_MPI2S_1D,ALG_CAF_MULTI_PULL_1D, &
        ALG_CAF_SINGLE_PULL_1D,ALG_CAF_SINGLE_PUSH_1D, &
        boundary_exchange_algorithm

    implicit none
    private

    integer, parameter, private ::  &
       i4 = selected_int_kind(6),  &
       r8 = selected_real_kind(13)
   
    interface GSHalo_update
       module procedure GSHalo_update_2d_int
       module procedure GSHalo_update_2d_dbl 
    end interface

    public :: GSHalo_init, GSHalo_update

  contains

    !<
    !! Updates the a 1-diminsional double precision array
    !!
    !! @param handle  Previous was a pointer to memory. [REMOVE]
    !! @param array   Array on which to perform the halo update.
    !>
    subroutine GSHalo_update_2d_dbl(handle,array)
        integer(i4), intent(in)    :: handle
        real(r8), intent(inout) :: array(:)

        select case (boundary_exchange_algorithm)
            case (ALG_MPI2S_1D)
                call mpi2s_gshalo_update_2d_dbl(handle, array)
        end select
    end subroutine GSHalo_update_2d_dbl

    !<
    !! Updates the a 1-diminsional integer array
    !!
    !! @param handle  Previous was a pointer to memory. [REMOVE]
    !! @param array   Array on which to perform the halo update.
    !>
    subroutine GSHalo_update_2d_int(handle,array)
        integer(i4), intent(in)    :: handle
        integer(i4), intent(inout) :: array(:)

        select case (boundary_exchange_algorithm)
            case (ALG_MPI2S_1D)
                call mpi2s_gshalo_update_2d_int(handle, array)
        end select
    end subroutine GSHalo_update_2d_int
    
    !<
    !! Construct the message passing metadata necessary to perform the MPI based
    !! boundary exchanges necessary to update the halo region. 
    !!
    !! @param COMM        MPI communicator
    !! @param maxlinear   The maximum number of gridpoints in the 1D data
    !!                    structure.
    !! @param nTotal      Total number of gridpoints including active and halo
    !!                    points.
    !! @param nActive     The total number of active ocean gridpoints.
    !! @param LinearGdof  A 1D array which contains the global degrees of
    !!                    freedom for this MPI task.
    !! @param LinearProc  A 1D array which contains the MPI rank or task for
    !!                    each gridpoint.
    !>

    function GSHalo_init(COMM,maxlinear, nTotal,nActive,LinearGdof,LinearProc) &
        result(handle)

        integer(i4), intent(in) :: COMM
        integer(i4), intent(in) :: maxlinear
        integer(i4), intent(in) :: nTotal,nActive
        integer(i4), intent(in) :: LinearGdof(maxlinear)
        integer(i4), intent(in) :: LinearProc(maxlinear)

        integer(i4) :: handle

        select case (boundary_exchange_algorithm)
            case (ALG_MPI2S_1D)
                handle = mpi2s_gshalo_init( &
                    COMM,maxlinear, nTotal,nActive,LinearGdof,LinearProc)
        end select
    end function GSHalo_init
end module gshalo
