!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! This is the main driver for the standalone Parallel Ocean Program (POP) 
!! conjugate gradient solver.
!<
program CGpop
    use kinds_mod, only: i4, r8
    use simple_blocks, only: nx_block,ny_block,nblocks_tot
    use domain_size, only: nx_global,ny_global, max_blocks_tropic
    use communicate, only: init_communicate, my_task, master_task, &
        MPI_COMM_OCN, get_num_procs, exit_message_environment
    use simple_domain, only: init_domain_blocks, init_domain_distribution
    use solvers, only: esolver, init_solvers
    use constants, only: init_constants!, ntrials
    use timers, only: get_timer, timer_start, timer_stop, timer_print_all, &
        init_timers
    use check, only: CheckAnswers

    implicit none

    real(r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
        PRESSI, & ! Initial guess for pressure (x in the Ax=b equation)
        RHS,    & ! Right-hand side of equation Ax=b
        PRESSF, & ! Expected resulting pressure (expected value of x)
        PRESS     ! Solved result (calculated value of x)

    integer(i4) :: timer_mpi2s_solver_1D  ! Timer for how long it takes to do
                                          ! ntrials of the CG algorithm
    integer(i4), pointer :: reorder(:)
    integer(i4) :: ntrials=1

!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------
    integer (i4) :: nstep
    integer (i4) :: ierr
    integer (i4) :: n, nprocs

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

!-----------------------------------------------------------------------
! Initialize the solver
!-----------------------------------------------------------------------

    !-----------------------------------
    ! 2 sided MPI with 1D data structure
    !-----------------------------------
    call get_timer(timer_mpi2s_solver_1D,'MPI2S_1D',1,nprocs)

    call init_solvers(RHS,PRESSI,PRESSF)
    call MPI_BARRIER(MPI_COMM_OCN, ierr)
    call timer_start(timer_mpi2s_solver_1d)
    do n=1,ntrials
        PRESS=PRESSI
        call esolver(RHS,PRESS)
    enddo
    call timer_stop(timer_mpi2s_solver_1d)
    !----------------------------------
    ! check the accuracy of the solver
    !----------------------------------
    call CheckAnswers('MPI2S_1D',PRESSF,PRESS)

    if(my_task == master_task) then 
       write(*,*) 'Number of trials: ',ntrials 
    endif
    call timer_print_all()

    call exit_message_environment(ierr)
end program CGpop
