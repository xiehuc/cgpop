program CGpop
    ! !DESCRIPTION:
    !  This is the main driver for the standalone Parallel Ocean Program (POP) 
    !  conjugate gradient solver.

    ! !USES:
    use kinds_mod, only: i4, r8
    use simple_blocks, only: nx_block,ny_block,nblocks_tot
    use domain_size, only: nx_global,ny_global, max_blocks_tropic
    use communicate, only: init_communicate, my_task, master_task, &
        MPI_COMM_OCN, get_num_procs, exit_message_environment
    use simple_domain, only: init_domain_blocks, init_domain_distribution
    use solvers, only: esolver, init_solvers
    use constants, only: init_constants, boundary_exchange_algorithm, &
        ALG_MPI2S_1D, ALG_MPI2S_2D, ALG_CAF_MULTI_PULL_1D, &
        ALG_CAF_SINGLE_PULL_1D, ALG_CAF_SINGLE_PUSH_1D, &
        ALG_CAF_SINGLE_PUSH_2D, ALG_CAF_SINGLE_PULL_2D,&
        ALG_MPI1S_MULTI_PULL_1D, ALG_MPI1S_SINGLE_PULL_1D, &
        ALG_MPI1S_SINGLE_PUSH_1D!,ntrials
    use timers, only: get_timer,timer_start,timer_stop,timer_print_all, &
        init_timers
    use check, only: CheckAnswers

    implicit none

    include 'mpif.h'

    integer :: ntrials = 1

    real(r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
        PRESSI, RHS, PRESSF, PRESS

    real(r8) :: sum_diff,gdiff

    integer(i4) :: timer_esolver
    integer(i4) :: timer_solver_init
    integer(i4) :: timer_mpi1s_multi_pull_solver_1D, &
        timer_mpi1s_single_push_solver_1D, &
        timer_mpi1s_single_pull_solver_1D

    integer(i4), pointer :: reorder(:)
    integer(i4) :: winID

    real(r8) :: arrayBuf(1:10000)
    pointer (p,arrayBuf)
    integer :: linearMaximum
    integer(kind=MPI_ADDRESS_KIND) :: sz

    !-----------------------------------------------------------------------
    !
    !  local variables
    !
    !-----------------------------------------------------------------------
    integer (i4) :: nscan,nstep
    integer (i4) :: iblock, ierr
    integer (i4) :: n,nprocs

    logical, parameter :: &
        DO_MPI1S_SINGLE_PUSH_1D = .TRUE.

    !-----------------------------------------------------------------------
    !  initialize message-passing or other communication protocol
    !-----------------------------------------------------------------------
    call init_communicate

    nprocs = get_num_procs()

    !-----------------------------------------------------------------------
    !  initialize constants and i/o stuff
    !-----------------------------------------------------------------------
    call init_constants

    !-----------------------------------------------------------------------
    !  initialize domain and grid
    !-----------------------------------------------------------------------
    call init_domain_blocks(reorder)
    call init_domain_distribution(reorder)

    !-----------------------------------------------------------------------
    !  initialize timers and additional communication routines
    !-----------------------------------------------------------------------
    call init_timers()
    nscan = 0

    !-----------------------------------------------------------------------
    ! Initialize the solver
    !-----------------------------------------------------------------------
    !-----------------------------------------
    ! MPI1S single push with 1D data structure
    !-----------------------------------------
    if(DO_MPI1S_SINGLE_PUSH_1D) then 
        call get_timer( &
            timer_mpi1s_single_push_solver_1D,'MPI1S_SINGLE_PUSH_1D',1,nprocs)
        boundary_exchange_algorithm = ALG_MPI1S_SINGLE_PUSH_1D
        call init_solvers(RHS,PRESSI,PRESSF)

        call MPI_BARRIER(MPI_COMM_OCN, ierr)
        call timer_start(timer_mpi1s_single_push_solver_1D)
        do n=1,ntrials
            PRESS=PRESSI
            call esolver(RHS,PRESS)
        enddo
        call timer_stop(timer_mpi1s_single_push_solver_1D)

        !----------------------------------
        ! check the accuracy of the solver
        !----------------------------------
        call CheckAnswers('MPI1S_SINGLE_PUSH_1D',PRESSF,PRESS)
    endif

    if(my_task == master_task) then 
        write(*,*) 'Number of trials: ',ntrials
    endif
    call timer_print_all()

    call exit_message_environment(ierr)
end program CGpop
