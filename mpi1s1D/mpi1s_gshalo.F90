! Multi pull version using 1-sided MPI
module mpi1s_gshalo
    implicit none
    private

    include 'mpif.h'

    integer, parameter, private :: &
        i4 = selected_int_kind(6),  &
        r8 = selected_real_kind(13)

    integer(i4), parameter :: NUMexpectMSG  = 101, &
        expectMSG     = 102, &
        UpdateHaloMSG = 103

    integer(i4) :: my_task,nprocs
    integer(i4), private, allocatable :: Rrequest(:), Srequest(:)
    integer(i4), private, allocatable :: Rstatus(:,:), Sstatus(:,:)
    integer(i4) :: mpi_COMM

    integer(i4) :: nPrint

    ! ==========================
    ! indirect addressing arrays
    ! ==========================
    integer(i4), allocatable :: halo2send(:)
    integer(i4), allocatable :: recv2halo(:)

    ! ======================
    ! Message buffer arrays
    ! ======================
    integer(i4), allocatable :: intBuffer(:)
    real(r8), allocatable :: realBuffer(:)

    ! Communication metadata
    integer(i4), allocatable :: halo2grab(:)
    integer(i4), allocatable :: linear2Proc(:)
    integer(i4) :: linearMaximum, nTotalCA, nActiveCA
    integer(i4) :: intWinID, realWinID
    public :: mpi1s_gshalo_update_1d_dbl
    public :: mpi1s_gshalo_update_1d_int
    public :: mpi1s_gshalo_init

    integer(i4), public :: timer_mpi1s_gshalo_init, &
        timer_mpi1s_gshalo_update_dbl

  contains

    !***********************************************************************
    !>
    !! This subroutine performs the boundary exchange operation in order to
    !! update the halo elements in the passed array.
    !<
    subroutine mpi1s_gshalo_update_1d_dbl(handle,array)
        ! !INPUT PARAMETERS:
        integer(i4), intent(in) :: handle
        real(r8), intent(inout) :: array(:) ! vector to update

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer(i4) :: i    ! counter
        integer(i4) :: ierr ! MPI Status code
        integer(kind=MPI_ADDRESS_KIND) :: pos

        !-------------------------------------------------------------------
        !
        !  in order to perform the boundary exchange operation, the array data
        !  must be put into the buffer that we've placed a window over.
        !  There may be a way of updating this function to not do this copy;
        !  but, for now, I'm going to go ahead and make the copy.
        !
        !   TODO: Determine if this copy matters or if we can get away without
        !         having it.
        !
        !-------------------------------------------------------------------
        realBuffer(1:ubound(array,1)) = array(:)

        !-------------------------------------------------------------------
        !
        !  iterate through halo elements, grabbing fresh values from remote
        !  images as needed.
        !
        !-------------------------------------------------------------------
        call MPI_Win_fence(0, realWinID, ierr)
        do i=nActiveCA+1,nTotalCA
            pos = int(halo2grab(i)-1, MPI_ADDRESS_KIND)
            call MPI_Get(array(i), 1, MPI_REAL8, linear2Proc(i), pos, 1, &
                MPI_REAL8, realWinID, ierr)
        enddo
        call MPI_Win_fence(0, realWinID, ierr)
    end subroutine mpi1s_gshalo_update_1d_dbl

    !***********************************************************************
    !>
    !! This subroutine performs the boundary exchange operation in order to
    !! update the halo elements in the passed array.
    !<
    subroutine mpi1s_gshalo_update_1d_int(handle,array)
        ! !INPUT PARAMETERS:
        integer(i4), intent(in)    :: handle
        integer(i4), intent(inout) :: array(:) ! vector to update

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer(i4) :: i ! dummy counter
        integer(i4) :: ierr ! MPI Status code
        integer(kind=MPI_ADDRESS_KIND) :: pos

        integer(i4) :: tmp

        !-------------------------------------------------------------------
        !
        !  in order to perform the boundary exchange operation, the array data
        !  must be put into the buffer that we've placed a window over.
        !  There may be a way of updating this function to not do this copy;
        !  but, for now, I'm going to go ahead and make the copy.
        !
        !-------------------------------------------------------------------
        intBuffer(1:ubound(array,1)) = array(:)

        !-------------------------------------------------------------------
        !
        !  iterate through halo elements, grabbing fresh values from remote
        !  images as needed.
        !
        !-------------------------------------------------------------------
        call MPI_Win_fence(0, intWinID, ierr)
        do i=nActiveCA+1,nTotalCA
            pos = int(halo2grab(i)-1, MPI_ADDRESS_KIND)
            call MPI_Get(array(i), 1, MPI_INTEGER, linear2Proc(i), pos, 1, &
            MPI_INTEGER, intWinID, ierr)
        enddo
        call MPI_Win_fence(0, intWinID, ierr)
    end subroutine mpi1s_gshalo_update_1d_int

    !***********************************************************************
    !>
    !! This function initilizes scheduling metadata needed for the update
    !! routines (the update routines perform the boundary exchange
    !! operation.)
    !<
    function mpi1s_gshalo_init(COMM,maxlinear, nTotal,nActive,LinearGdof, &
        LinearProc) result(handle)
        ! !INPUT PARAMETERS:
        integer(i4), intent(in) :: &
            COMM,                  &    ! MPI communicator
            maxlinear,             &    ! Maximum size of vectors
            nTotal,                &    ! Total elements in local vector
            nActive,               &    ! Number of active (non-halo) elements
            LinearGdof(maxlinear), &    ! Maps vector indices to GDOFs
            LinearProc(maxlinear)       ! Maps halo indices to processors
        integer(i4) :: handle

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer(i4), save, allocatable :: &
            LinearGdofCA(:)  ! Co-array copy of LinearGdof
        integer(i4), allocatable :: &
            copy(:)        ! For temp local copies of remote LinearGdfCA arrays
        integer(i4) :: i, j     ! Dummy loop-index variables
        integer(i4) :: idx      ! Vector index of a remote gdof element
        integer(i4) :: ierr
        integer(kind=MPI_ADDRESS_KIND) :: sz
        integer :: info 

        mpi_COMM = COMM
        linearMaximum = maxlinear
        nTotalCA = nTotal
        nActiveCA = nActive

        allocate(copy(maxlinear))
        allocate(linear2Proc(maxlinear))

        allocate(LinearGdofCA(maxlinear))
        LinearGdofCA(:) = LinearGdof(:)

        linear2Proc(:) = LinearProc(:)

        allocate(halo2grab(nActive+1:nTotal))

        !-------------------------------------------------------------------
        !
        !  iterate through halo points, constructing halo2grab array (halo2grab
        !  maps elements in this processor's halo to indices in a remote
        !  processor's vector)
        !
        !-------------------------------------------------------------------
        do_i: do i=nActive+1,nTotal
            idx = -1

            copy(:) = LinearGdofCA(:)

            do_j: do j=1,maxlinear
                if(copy(j) == LinearGdofCA(i)) then
                    idx = j
                    exit do_j
                end if
            end do do_j

            halo2grab(i) = idx
        enddo do_i

        deallocate(LinearGdofCA)

        ! -----------------------------------------------------------------
        ! Setup windows for 1-sided MPI communication
        ! -----------------------------------------------------------------
        sz = maxlinear * 4
        allocate(intBuffer(maxlinear))
        print *,'mpi1s_gshalo_init: before win_create #1'
        call MPI_Win_create(intBuffer, sz, 4, MPI_INFO_NULL, MPI_COMM, &
            intWinID, ierr)

        sz = maxlinear * 8
        allocate(realBuffer(maxlinear))
        print *,'mpi1s_gshalo_init: before win_create #2'
        call MPI_Win_create(realBuffer, sz, 8, MPI_INFO_NULL, MPI_COMM, &
            realWinID, ierr)

        call MPI_Barrier(MPI_COMM, ierr)
    end function mpi1s_gshalo_init
end module mpi1s_gshalo
