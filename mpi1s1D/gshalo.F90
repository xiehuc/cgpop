! This acts as a proxy for one of the boundary exchange modules.  The value
! of boundary_exchange_algorithm will determine which boundary  exchange
! algorithm to use.
module gshalo
    use mpi1s_gshalo
    !use mpi1s_single_pull_halo
    use mpi1s_single_push_halo
    use constants, only: ALG_MPI1S_MULTI_PULL_1D, ALG_MPI1S_SINGLE_PULL_1D, &
        ALG_MPI1S_SINGLE_PUSH_1D, boundary_exchange_algorithm

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

    subroutine GSHalo_update_2d_dbl(handle,array)
        integer(i4), intent(in)    :: handle
        real(r8), intent(inout) :: array(:)

        select case (boundary_exchange_algorithm)
            !case (ALG_MPI1S_MULTI_PULL_1D)
            !    call mpi1s_gshalo_update_1d_dbl(handle, array)
            !case (ALG_MPI1S_SINGLE_PULL_1D)
            !    call mpi1s_single_pull_update_1d_dbl(handle, array)
            case (ALG_MPI1S_SINGLE_PUSH_1D)
                call mpi1s_single_push_update_1d_dbl(handle, array)
        end select
    end subroutine GSHalo_update_2d_dbl

    subroutine GSHalo_update_2d_int(handle,array)
        integer(i4), intent(in)    :: handle
        integer(i4), intent(inout) :: array(:)

        select case (boundary_exchange_algorithm)
            !case (ALG_MPI1S_MULTI_PULL_1D)
            !    call mpi1s_gshalo_update_1d_int(handle, array)
            !case (ALG_MPI1S_SINGLE_PULL_1D)
            !    call mpi1s_single_pull_update_1d_int(handle, array)
            case (ALG_MPI1S_SINGLE_PUSH_1D)
                call mpi1s_single_push_update_1d_int(handle, array)
        end select
    end subroutine GSHalo_update_2d_int

    function GSHalo_init(COMM,maxlinear, nTotal,nActive,LinearGdof,LinearProc) &
        result(handle)

        integer(i4), intent(in) :: COMM
        integer(i4), intent(in) :: maxlinear
        integer(i4), intent(in) :: nTotal,nActive
        integer(i4), intent(in) :: LinearGdof(maxlinear)
        integer(i4), intent(in) :: LinearProc(maxlinear)

        integer(i4) :: handle

        select case (boundary_exchange_algorithm)
            !case (ALG_MPI1S_MULTI_PULL_1D)
            !    handle = mpi1s_gshalo_init( &
            !        COMM,maxlinear, nTotal,nActive,LinearGdof,LinearProc)
            !case (ALG_MPI1S_SINGLE_PULL_1D)
            !    handle = mpi1s_single_pull_init( &
            !        COMM,maxlinear, nTotal,nActive,LinearGdof,LinearProc)
            case (ALG_MPI1S_SINGLE_PUSH_1D)
                handle = mpi1s_single_push_init( &
                    COMM,maxlinear, nTotal,nActive,LinearGdof,LinearProc)
        end select
    end function GSHalo_init
end module gshalo
