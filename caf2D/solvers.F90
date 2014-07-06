!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! This module contains routines and operators for solving the elliptic
!! system for surface pressure in the barotropic mode.
!<
module solvers
    ! !USES:
    use kinds_mod, only: i4, r8, char_len, log_kind
    use simple_type, only : distrb
    use simple_blocks, only: nx_block,ny_block, get_block_parameter
    use communicate, only: my_task, master_task
    use simple_domain, only: distrb_tropic, nblocks_tropic, blocks_tropic
    use simple_domain, only: read_solverioT

    use domain_size, only: max_blocks_tropic, nx_global, ny_global
    use constants, only: c0, c1, field_type_scalar, field_loc_center, p25, &
        eps, blank_fmt, delim_fmt, boundary_exchange_algorithm, &
        ALG_MPI2S_1D, ALG_MPI2S_2D, solv_max_iters, ALG_CAF_MULTI_PULL_1D, &
        ALG_CAF_SINGLE_PULL_1D, ALG_CAF_SINGLE_PUSH_1D, &
        ALG_CAF_SINGLE_PUSH_2D, ALG_CAF_SINGLE_PULL_2D
    use boundary, only: update_ghost_cells, create_boundary
    use broadcast, only: broadcast_scalar
    use IOUnitsMod, only: stdout, nml_filename, nmlin
    use exit_mod, only: sigAbort, exit_POP
    use boundary_types, only: bndy_tropic
    use reductions, only: global_sum

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS:
    public :: init_solvers
    public :: esolver
    public :: pcg_chrongear

    ! !PUBLIC DATA MEMBERS:
    integer (i4), public :: &
        solv_sum_iters      ! accumulated no of iterations (diagnostic)

    real (r8), public ::  &
        rms_residual        ! residual (also a diagnostic)

    logical :: strip

    integer (i4), save :: &
        timer_tolinear,           &! timer for the conversion to linear form
        timer_to2D,               &! time for the conversion to 2D form
        timer_matvec,timer_gsum,timer_halo, &
        timer_initDOF,            &!  Time to setup the DOF based communicator 
        timer_solver               ! time to perform the CG solver

    !-----------------------------------------------------------------------
    !
    !  other operator and preconditioner weights for barotropic operator
    !
    !-----------------------------------------------------------------------
    real (r8), public, dimension &
        (nx_block,ny_block,max_blocks_tropic)[*] :: & 
            A0,AN,AE,ANE,         &! barotropic (9pt) operator coefficients
            RCALCT_B               ! land mask in barotropic distribution 

    real (r8), public, dimension(nx_block,ny_block,max_blocks_tropic) :: &
        A0_read,AN_read,AE_read,ANE_read,RCALCT_B_read,RHS_read, &
        PRESSI_read,PRESSF_read

    !-----------------------------------------------------------------------
    !
    !  scalar convergence-related variables
    !
    !-----------------------------------------------------------------------
    logical (log_kind) :: &
        lprecond            ! true if computed preconditioner to be used

    real (r8) ::          &
        sor,               &! for jacobi solver
        resid_norm          ! residual normalization

    integer (i4), parameter :: &
        solv_pcg = 1,      &! predefined solver types
        solv_cgr = 2,      &
        solv_jac = 3,      &
        solv_cg1 = 4,      &
        solv_lcg1 = 5,     &
        solv_lpcg = 6

    integer (i4) :: &
        solv_itype   ! integer solver method (1=pcg, 2=cgr, 3=jac)

    logical :: solv_linear        

  contains

    !***********************************************************************
    !>
    !! This subroutine calls the eliptic solver
    !!
    !! @param RHS   The right hand side or 'b' in the equation 'Ax=b'.
    !! @param PRESS The surface pressure or 'x' in the equation 'Ax=b'.
    !!
    !<
    subroutine esolver(RHS,PRESS)
        real(r8), dimension(nx_block,ny_block,max_blocks_tropic)[*] ::  &
            RHS,PRESS
        integer(i4) :: i

        call pcg_chrongear(PRESS,RHS) ! Chron-Gear pcg solver
    end subroutine esolver

    !***********************************************************************
    !>
    !! This routine initializes choice of solver, calculates the 
    !! coefficients of the 9-point stencils for the barotropic operator and
    !! reads in a preconditioner if requested.
    !!
    !! @param RHS     The right hand side or 'b' in the equation 'Ax=b'.
    !! @param PRESSI  The initial guess for the surface pressure 'x'.
    !! @param PRESSF  The final estimate for the surface pressure 'x'.
    !<
    subroutine init_solvers(RHS,PRESSI,PRESSF)
        !-------------------------------------------------------------------
        !
        !  local variables:
        !
        !       {X,Y}{NE,SE,NW,SW} = contribution to {ne,se,nw,sw} coefficients 
        !         from {x,y} components of divergence
        !       HU = depth at U points
        !
        !-------------------------------------------------------------------
        real (r8), dimension(nx_block,ny_block,max_blocks_tropic)[*] :: &
            RHS, PRESSI, PRESSF

        character (char_len) :: &
            solv_type,           &! user choice of solver method
            precond_file          ! file containing preconditioner

        character(len=80) :: tilefile
        namelist /input_nml/ tilefile

        integer (i4) :: &
            i,j,n,             &! dummy counter
            iblock,            &! block counter
            nu,                &! I/O unit number and status flag
            nml_error           ! namelist i/o error flag

        integer (i4),dimension(:), allocatable :: &
            icheck              ! check for PC/mask compatibility in block

        real(r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
            WORK0_read,RCALCT_read

        logical (log_kind) :: &
            mlandne, mlandnw, mlandse, mlandsw ! land mask at nbr points

        logical,parameter :: Debug = .FALSE.
        integer(i4) :: nstep

        real(r8) :: sum_work0,sum_rcalct
        integer(i4) :: ie,ib,je,jb
        real(r8) :: diff

        character (char_len) ::  &
            ew_boundary_type,   &
            ns_boundary_type

        !-------------------------------------------------------------------
        !
        !  read solver choice and solver constants from namelist input
        !  (namelist input file opened in initial.F)
        !
        !-------------------------------------------------------------------
        lprecond       = .false.
        solv_type      = 'cg1'
        precond_file   = 'empty'

        if (my_task == master_task) then
            open (nmlin, file=nml_filename, status='old',iostat=nml_error)
            if (nml_error /= 0) then
                nml_error = -1
            else
                nml_error =  1
            endif
            do while (nml_error > 0)
                read(nmlin, nml=input_nml,iostat=nml_error)
            end do
            if (nml_error == 0) close(nmlin)
        endif

        call broadcast_scalar(nml_error, master_task)
        if (nml_error /= 0) then
            call exit_POP(sigAbort,'ERROR reading solver_nml')
        endif

        if (my_task == master_task) then
            write(stdout,delim_fmt)
            write(stdout,blank_fmt)
            write(stdout,'(a35)') ' Solver options (barotropic solver)'
            write(stdout,blank_fmt)
            write(stdout,delim_fmt)

            write(stdout,'(a13)') ' Solver type:'
            select case(solv_type)
                case('cg1')
                    write(stdout,'(a38)') &
                        '  Chron-Gear preconditioned conjg grad'
                    solv_itype = solv_cg1
                case('lcg1')
                    write(stdout,'(a46)') &
                        ' Linear Chron-Gear preconditioner conjg grad'
                    solv_itype = solv_lcg1
                    solv_linear = .TRUE.
                case('pcg')
                    write(stdout,'(a35)') '  Preconditioned Conjugate Gradient'
                    solv_itype = solv_pcg
                case('lpcg')
                    write(stdout,'(a42)') &
                        ' Linear Preconditioned Conjugate Gradient'
                    solv_itype = solv_lpcg
                    solv_linear = .TRUE.
                case default
                    solv_itype = -1000
            end select

            write(stdout,'(a29,i6)') ' Solver maximum iterations = ', &
                solv_max_iters
        endif

        solv_linear = .true.
        call broadcast_scalar(tilefile,       master_task)

        if (solv_itype == -1000) then
            call exit_POP(sigAbort, &
                'unknown solver type: must be cg1, lcg1, pcg, lpcg')
        endif

        !-------------------------------------------------------------------
        !
        !  set sor for jacobi solver
        !
        !-------------------------------------------------------------------
        sor = p25     ! should be < 1/2

        !-------------------------------------------------------------------
        !
        !  compute nine point operator coefficients: compute on baroclinic
        !  decomposition first where grid info defined and redistribute
        !  to barotropic distribution
        !  leave A0,AC in baroclinic distribution to facilitate easy
        !  time-dependent changes in barotropic routine
        !
        !-------------------------------------------------------------------
        A0 = 0.0_r8
        AN = 0.0_r8
        AE = 0.0_r8
        ANE = 0.0_r8
        RHS = 0.0_r8
        PRESSI = 0.0_r8
        PRESSF = 0.0_r8
        RCALCT_B = 0.0_r8
        call read_solverioT(tilefile,nstep,A0,AN,AE,ANE,RHS, &
            PRESSI, PRESSF, RCALCT_B,WORK0_read)

        !-------------------------------------------------------------------
        !
        !  calculate normalization constant (darea,darea) for rms_residual
        !  in cgr routine.
        !
        !-------------------------------------------------------------------
        !   resid_norm = c1/global_sum(WORK0_read, distrb_tropic,
        !       field_loc_center, RCALCT_B)
        !-------------------------------------------------------------------
        !
        !  setup preconditioner if required
        !
        !-------------------------------------------------------------------
        if (lprecond) then
            call exit_POP(sigAbort,'This option not currently supported')
        else ! no preconditioner
            if (my_task == master_task) then
                write(stdout,'(a18)') ' No preconditioner'
            endif
        endif

        if(my_task == master_task) then 
            write(stdout,*) 'init_solvers: Using solver V10'
        endif

        !------------------------------------------------------------------
        !
        !  set up ghost cell updates for each distribution
        !  Boundary types are cyclic, closed, or tripole
        !
        !------------------------------------------------------------------
        ns_boundary_type = 'closed'
        ew_boundary_type = 'cyclic'
        call create_boundary(bndy_tropic, distrb_tropic, blocks_tropic,&
            trim(ns_boundary_type), trim(ew_boundary_type), nx_global, &
                ny_global)
    end subroutine init_solvers

    !***********************************************************************
    !>
    !!  This routine implements the Chronopoulos-Gear conjugate-gradient
    !!  solver with preconditioner for solving the linear system $Ax=b$.
    !!  It is a rearranged conjugate gradient solver that reduces the
    !!  number of inner products per iteration from two to one (not
    !!  counting convergence check). Both the operator $A$ and
    !!  preconditioner are nine-point stencils. If no preconditioner has
    !!  been supplied, a diagonal preconditioner is applied.
    !!
    !!
    !!  References:
    !!     Dongarra, J. and V. Eijkhout. LAPACK Working Note 159.
    !!        Finite-choice algorithm optimization in conjugate gradients.
    !!        Tech. Rep. ut-cs-03-502. Computer Science Department.
    !!        University of Tennessee, Knoxville. 2003.
    !!
    !!     D Azevedo, E.F., V.L. Eijkhout, and C.H. Romine. LAPACK Working
    !!        Note 56. Conjugate gradient algorithms with reduced
    !!        synchronization overhead on distributed memory multiprocessors.
    !!        Tech. Rep. CS-93-185.  Computer Science Department.
    !!        University of Tennessee, Knoxville. 1993.
    !!
    !! @param X     on input,  an initial guess for the solution;
    !!              on output, solution of the linear system
    !! @param B     right hand side of linear system
    !<
    subroutine pcg_chrongear(X,B)
        ! !INPUT PARAMETERS:
        real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
            intent(in) :: B

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
            intent(inout) :: X

        !-------------------------------------------------------------------
        character (char_len) :: &
            noconvrg           ! error message for no convergence

        integer (i4) :: &
            i,j,m,             &! local iteration counter
            iblock              ! local block     counter

        real (r8) :: & ! scalar results
            cg_alpha, cg_beta, cg_sigma, cg_delta, cg_rho_old, cg_rho, rr

        real (r8), save, dimension(nx_block,ny_block,max_blocks_tropic)[*] :: &
            R,                  &! residual (b-Ax)
            S,                  &! conjugate direction vector
            Q,Z,AZ,WORK0,       &! various cg intermediate results
            A0R                  ! diagonal preconditioner

        real (r8) :: rr0

        real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: WORKJMD

        real (r8), dimension(nx_block,ny_block,2,max_blocks_tropic) :: &
            WORKN              ! WORK array

        real (r8), dimension(2) :: &
            sumN               ! global sum results for multiple arrays

        integer (i4) :: gid

        !-------------------------------------------------------------------
        !
        !  initialize some scalars
        !
        !-------------------------------------------------------------------
        cg_rho = c1
        solv_sum_iters = solv_max_iters

        !-------------------------------------------------------------------
        !
        !  compute initial residual and initialize other arrays
        !
        !-------------------------------------------------------------------
        !$OMP PARALLEL DO PRIVATE(iblock,gid)
        do iblock=1,nblocks_tropic
            gid = blocks_tropic(iblock)

            R    (:,:,iblock) = c0
            S    (:,:,iblock) = c0
            Z    (:,:,iblock) = c0
            Q    (:,:,iblock) = c0
            AZ   (:,:,iblock) = c0
            WORK0(:,:,iblock) = c0
            WORKN(:,:,:,iblock) = c0

            !--- diagonal preconditioner if preconditioner not specified
            do j=1,ny_block
                do i=1,nx_block
                    if (A0(i,j,iblock) /= c0) then
                        A0R(i,j,iblock) = c1/A0(i,j,iblock)
                    else
                        A0R(i,j,iblock) = c0
                    endif
                end do
            end do

            ! use S as a temp here for Ax
            call btrop_operator(S,X,gid,iblock)
            R(:,:,iblock) = B(:,:,iblock) - S(:,:,iblock) ! b-Ax
            WORKJMD(:,:,iblock) = R(:,:,iblock)*R(:,:,iblock)
        end do ! block loop
        !$OMP END PARALLEL DO

        rr0 = &
            sqrt(global_sum(WORKJMD, distrb_tropic, field_loc_center, RCALCT_B))

        call update_ghost_cells(R, bndy_tropic, field_loc_center, &
            field_type_scalar)

        !-------------------------------------------------------------------
        !
        !    take one pass of standard algorithm
        !
        !-------------------------------------------------------------------
        !$OMP PARALLEL DO PRIVATE(iblock,gid)
        do iblock=1,nblocks_tropic
            gid = blocks_tropic(iblock)

            !---- calculate (PC)r store in Z
            Z(:,:,iblock) = R(:,:,iblock)*A0R(:,:,iblock)

            !---- Compute intermediate result for dot product
            WORKN(:,:,1,iblock) = R(:,:,iblock)*Z(:,:,iblock)


            !---- update conjugate direction vector S
            S(:,:,iblock) =  Z(:,:,iblock)

            !---- compute Q = A * S
            call btrop_operator(Q,S,gid,iblock)

            !---- compute intermediate result for dot product
            WORKN(:,:,2,iblock) = S(:,:,iblock)*Q(:,:,iblock)
        end do
        !$OMP END PARALLEL DO

        call update_ghost_cells(Q, bndy_tropic, field_loc_center,&
            field_type_scalar)

        !---- Form dot products
        sumN = global_sum(WORKN, distrb_tropic, field_loc_center, RCALCT_B)

        cg_rho_old = sumN(1) !(r,PCr)
        cg_sigma   = sumN(2) !(s,As)
        cg_alpha   = cg_rho_old/cg_sigma

        !---- compute first solution and residual

        !$OMP PARALLEL DO PRIVATE(iblock)
        do iblock=1,nblocks_tropic
            X(:,:,iblock) = X(:,:,iblock) + cg_alpha*S(:,:,iblock)
            R(:,:,iblock) = R(:,:,iblock) - cg_alpha*Q(:,:,iblock)
        end do
        !$OMP END PARALLEL DO

        !-------------------------------------------------------------------
        !
        !     iterate
        !
        !-------------------------------------------------------------------
        iter_loop: do m = 1, solv_max_iters

            !---------------------------------------------------------------
            !
            !     calculate (PC)r and A*(Pc)r
            !
            !---------------------------------------------------------------
            !$OMP PARALLEL DO PRIVATE(iblock,gid)
            do iblock=1,nblocks_tropic
                gid = blocks_tropic(iblock)

                Z(:,:,iblock) = R(:,:,iblock)*A0R(:,:,iblock)

                call btrop_operator(AZ,Z,gid,iblock)

                !--- intermediate results for inner products
                WORKN(:,:,1,iblock) =  R(:,:,iblock)*Z(:,:,iblock)
                WORKN(:,:,2,iblock) = AZ(:,:,iblock)*Z(:,:,iblock)
            end do
            !$OMP END PARALLEL DO

            call update_ghost_cells(AZ, bndy_tropic, field_loc_center, &
                field_type_scalar)

            sumN = global_sum(WORKN, distrb_tropic,field_loc_center, RCALCT_B)

            cg_rho     = sumN(1)   ! (r,(PC)r)
            cg_delta   = sumN(2)   ! (A (PC)r,(PC)r)
            cg_beta    = cg_rho/cg_rho_old
            cg_sigma   = cg_delta - (cg_beta**2)*cg_sigma
            cg_alpha   = cg_rho/cg_sigma
            cg_rho_old = cg_rho

            !---------------------------------------------------------------
            !
            !     compute S and Q
            !     compute next solution and residual
            !
            !---------------------------------------------------------------
            !$OMP PARALLEL DO PRIVATE(iblock, gid)
            do iblock=1,nblocks_tropic
                S(:,:,iblock) =  Z(:,:,iblock) + cg_beta *S(:,:,iblock)
                Q(:,:,iblock) = AZ(:,:,iblock) + cg_beta *Q(:,:,iblock)
                X(:,:,iblock) =  X(:,:,iblock) + cg_alpha*S(:,:,iblock)
                R(:,:,iblock) =  R(:,:,iblock) - cg_alpha*Q(:,:,iblock)
            end do
            !$OMP END PARALLEL DO
        end do iter_loop

        rms_residual = sqrt(rr*resid_norm)

        if (solv_sum_iters == solv_max_iters) then
        endif
    end subroutine pcg_chrongear

    !***********************************************************************
    !>
    !! This routine applies the nine-point stencil operator for the
    !! barotropic solver.  It takes advantage of some 9pt weights being 
    !! shifted versions of others.
    !!
    !! @param AX        nine point operator result (AX)
    !! @param X         array to be operated on 
    !! @param gid       global block id for ths block
    !! @param bid       local block address for this block
    !<
    subroutine btrop_operator(AX,X,gid,bid)
        ! !INPUT PARAMETERS:
        real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
            intent(in) :: X
        integer (i4), intent(in) :: bid
        integer (i4), intent(in) :: gid

        ! !OUTPUT PARAMETERS:
        real (r8), dimension(nx_block,ny_block,max_blocks_tropic), &
            intent(out) :: AX

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer (i4) :: i,j                ! dummy counters
        integer (i4) :: ib,ie,jb,je


        AX(:,:,bid) = c0

        call get_block_parameter(gid,ib=ib,ie=ie,jb=jb,je=je)
        do j=jb,je
            do i=ib,ie
                AX(i,j,bid) = A0 (i  ,j  ,bid)*X(i  ,j  ,bid) + &
                AN (i  ,j  ,bid)*X(i  ,j+1,bid) + &
                AN (i  ,j-1,bid)*X(i  ,j-1,bid) + &
                AE (i  ,j  ,bid)*X(i+1,j  ,bid) + &
                AE (i-1,j  ,bid)*X(i-1,j  ,bid) + &
                ANE(i  ,j  ,bid)*X(i+1,j+1,bid) + &
                ANE(i  ,j-1,bid)*X(i+1,j-1,bid) + &
                ANE(i-1,j  ,bid)*X(i-1,j+1,bid) + &
                ANE(i-1,j-1,bid)*X(i-1,j-1,bid)
            end do
        end do
    end subroutine btrop_operator
end module solvers
