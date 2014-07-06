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
!! This is the main driver for the standalone Parallel Ocean Program (POP) 
!! conjugate gradient solver.
!<
program CGpop
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
        ALG_CAF_SINGLE_PUSH_2D, ALG_CAF_SINGLE_PULL_2D!, ntrials
        use timers, only: get_timer,timer_start,timer_stop,timer_print_all, &
            init_timers
    use check, only: CheckAnswers

    implicit none

    real(r8), dimension(nx_block,ny_block,max_blocks_tropic)[*] :: &
        PRESSI, RHS, PRESSF, PRESS

    real(r8) :: sum_diff,gdiff

    integer(i4) :: timer_esolver
    integer(i4) :: timer_solver_init

    integer(i4) :: timer_mpi2s_solver_1D, timer_mpi2s_solver_2D
    integer(i4) :: timer_caf_multi_pull_solver_1D, &
        timer_caf_multi_push_solver_1D
    integer(i4) :: timer_caf_single_pull_solver_1D, &
        timer_caf_single_push_solver_1D
    integer(i4) :: timer_caf_single_push_solver_2D, &
        timer_caf_single_pull_solver_2D

    integer(i4), pointer :: reorder(:)
    integer(i4), parameter :: ntrials=1

    !-----------------------------------------------------------------------
    !
    !  local variables
    !
    !-----------------------------------------------------------------------
    integer (i4) :: nscan,nstep
    integer (i4) :: iblock, ierr
    integer (i4) :: n,nprocs

    integer(i4) :: timer_esolver
    integer(i4) :: timer_solver_init

    logical, parameter :: &
        DO_CAF_SINGLE_PULL_1D  = .TRUE., &
        DO_CAF_SINGLE_PUSH_1D  = .TRUE.

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
    ! CAF single push with 1D data structure
    !-----------------------------------------
    if(DO_CAF_SINGLE_PUSH_1D) then 
        call get_timer(timer_caf_single_push_solver_1D,'CAF_SINGLE_PUSH_1D', &
            1,nprocs)
        boundary_exchange_algorithm = ALG_CAF_SINGLE_PUSH_1D
        call init_solvers(RHS,PRESSI,PRESSF)
        call MPI_BARRIER(MPI_COMM_OCN, ierr)
        call timer_start(timer_caf_single_push_solver_1D)
        do n=1,ntrials
            PRESS=PRESSI
            call esolver(RHS,PRESS)
        enddo
        call timer_stop(timer_caf_single_push_solver_1D)
        !----------------------------------
        ! check the accuracy of the solver
        !----------------------------------
        call CheckAnswers('CAF_SINGLE_PUSH_1D',PRESSF,PRESS)
    endif

    !-----------------------------------------
    ! CAF single pull with 1D data structure
    !-----------------------------------------
    if(DO_CAF_SINGLE_PULL_1D) then 
        call get_timer(timer_caf_single_pull_solver_1D,'CAF_SINGLE_PULL_1D', &
            1,nprocs)

        boundary_exchange_algorithm = ALG_CAF_SINGLE_PULL_1D
        call init_solvers(RHS,PRESSI,PRESSF)

        call MPI_BARRIER(MPI_COMM_OCN, ierr)
        call timer_start(timer_caf_single_pull_solver_1D)
        do n=1,ntrials
            PRESS=PRESSI
            call esolver(RHS,PRESS)
        enddo
        call timer_stop(timer_caf_single_pull_solver_1D)

        !----------------------------------
        ! check the accuracy of the solver
        !----------------------------------
        call CheckAnswers('CAF_SINGLE_PULL_1D',PRESSF,PRESS)
    endif

    call timer_print_all()

    call exit_message_environment(ierr)
end program CGpop
