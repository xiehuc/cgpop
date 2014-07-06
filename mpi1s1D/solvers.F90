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
    use matrix_mod, only : Matrix_t, matvec
    use reductions, only : global_sum

    use domain_size, only: max_blocks_tropic
    use constants, only: c0, c1, field_type_scalar, field_loc_center, p25, &
        eps, blank_fmt, delim_fmt, boundary_exchange_algorithm,solv_max_iters
    use broadcast, only: broadcast_scalar
    use IOUnitsMod, only: stdout, nml_filename, nmlin
    use exit_mod, only: sigAbort, exit_POP
    use linear, only: linearmask, ntotal, nactive, max_linear, &
    convert2dtolinear, update_halo, convertlinearto2d, initdof
    use matrix_mod, only: convertstencil, matvec, A

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS:
    public :: init_solvers
    public :: esolver

    ! !PUBLIC DATA MEMBERS:
    integer (i4), public :: &
        solv_sum_iters      ! accumulated no of iterations (diagnostic)

    real (r8), public ::  &
        rms_residual        ! residual (also a diagnostic)

    logical :: strip

    !-----------------------------------------------------------------------
    !
    !  other operator and preconditioner weights for barotropic operator
    !
    !-----------------------------------------------------------------------
    real (r8), public, dimension (nx_block,ny_block,max_blocks_tropic) &
        :: A0,AN,AE,ANE, & ! barotropic (9pt) operator coefficients
           RCALCT_B        ! land mask in barotropic distribution 

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
        solv_itype  ! integer solver method (1=pcg, 2=cgr, 3=jac)

    logical :: solv_linear        

  contains

    !***********************************************************************
    subroutine esolver(RHS,PRESS)
        real(r8), dimension(nx_block,ny_block,max_blocks_tropic) ::  RHS,PRESS

        real (r8), dimension(max_linear) :: &
            X_linear,    & ! The initial guess and solution in Linear form
            B_linear,    & ! The RHS in Linear form
            Mask_linear, & ! Mask for non-valid points
            Minv_linear    ! Diagonal preconditioner

        integer(i4) :: i

        !==============================
        ! The initial guess
        !==============================
        call Convert2DtoLinear(X_linear,PRESS)
        call update_halo(X_linear)

        !==============================
        ! the RHS
        !==============================
        call Convert2DtoLinear(B_linear,RHS)
        call update_halo(B_linear)

        !==============================
        ! The preconditioner
        !==============================
        call Convert2DtoLinear(Minv_linear,A0)
        do i=1,nActive
            Minv_linear(i) = 1.0D0/Minv_linear(i)
        enddo

        !===============================
        ! Change the form of update_halo
        !===============================
        call update_halo(Minv_linear)

        !==============================
        ! The A matrix
        !==============================
        call ConvertStencil(A,A0,AN,AE,ANE)

        if (my_task < distrb_tropic%nprocs) then
            call pcg_chrongear_linear(nTotal,nActive+1, &
                A,X_linear,B_linear,Minv_linear,LinearMask,solv_sum_iters, &
                    rms_residual)
        endif

        call ConvertLinearto2D(PRESS,X_linear)
    end subroutine esolver

    !***********************************************************************
    !>
    !! This routine initializes choice of solver, calculates the 
    !! coefficients of the 9-point stencils for the barotropic operator and
    !! reads in a preconditioner if requested.
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
        real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
            RHS, PRESSI, PRESSF

        character(len=80) :: fname

        character (char_len) :: &
            solv_type,           &! user choice of solver method
            precond_file          ! file containing preconditioner

        character (char_len) :: tilefile

        namelist /input_nml/ tilefile

        integer (i4) :: &
            i,j,n,      & ! dummy counter
            iblock,     & ! block counter
            ncheck,     & ! scalar for checking PC/mask compatibility
            nu,         & ! I/O unit number and status flag
            nml_error     ! namelist i/o error flag

        integer (i4),dimension(:), allocatable :: &
            icheck              ! check for PC/mask compatibility in block

        real(r8), dimension(nx_block,ny_block,max_blocks_tropic) :: WORK1

        logical (log_kind) :: &
            mlandne, mlandnw, mlandse, mlandsw ! land mask at nbr points

        logical,parameter :: Debug = .FALSE.
            integer(i4) :: nstep

        real(r8) :: sum_work0,sum_rcalct
        integer(i4) :: ie,ib,je,jb
        real(r8) :: diff

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
                    write(stdout,'(a35)') &
                        '  Preconditioned Conjugate Gradient'
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
        call broadcast_scalar(lprecond,       master_task)
        call broadcast_scalar(solv_itype,     master_task)
        call broadcast_scalar(precond_file,   master_task)
        call broadcast_scalar(solv_linear,    master_task)
        call broadcast_scalar(tilefile,       master_task)

        if (solv_itype == -1000) then
            call exit_POP(sigAbort, &
                'unknown solver type: must be cg1, lcg1, pcg, lpcg')
        endif
        if (my_task == master_task) then
            print *,'Input tilefile: ',TRIM(tilefile)
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
        PRESSI, PRESSF, RCALCT_B,WORK1)

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

        ! setup the degree of freedom for the linear data structure
        strip = .TRUE.  ! if true strips outs land points
        call initDOF(tilefile,strip,RCALCT_B)

        !====================================
        ! Create the full matrix
        !====================================
        A%maxNZ = 9*nTotal
        A%n = nActive
        call init_solvers1D(solv_max_iters,resid_norm)
    end subroutine init_solvers

    !***********************************************************************
    subroutine init_solvers1D(maxiters, norm)
        integer(i4) :: maxiters,ncheck
        real(r8) :: convrg,norm

        resid_norm = norm
    end subroutine init_solvers1D

    !***********************************************************************
    subroutine pcg_chrongear_linear(n,iptrHalo,A,X,B,Minv2,mask, &
        solv_sum_iters, rms_residual)
        ! !DESCRIPTION:
        !  This routine implements the Chronopoulos-Gear conjugate-gradient
        !  solver with preconditioner for solving the linear system $Ax=b$.
        !  It is a rearranged conjugate gradient solver that reduces the
        !  number of inner products per iteration from two to one (not
        !  counting convergence check). Both the operator $A$ and
        !  preconditioner are nine-point stencils. If no preconditioner has
        !  been supplied, a diagonal preconditioner is applied.  Convergence
        !  is checked every {\em ncheck} steps.
        !
        !
        !  References:
        !     Dongarra, J. and V. Eijkhout. LAPACK Working Note 159.
        !        Finite-choice algorithm optimization in conjugate gradients.
        !        Tech. Rep. ut-cs-03-502. Computer Science Department.
        !        University of Tennessee, Knoxville. 2003.
        !
        !     D Azevedo, E.F., V.L. Eijkhout, and C.H. Romine. LAPACK Working
        !        Note 56. Conjugate gradient algorithms with reduced
        !        synchronization overhead on distributed memory multiprocessors.
        !        Tech. Rep. CS-93-185.  Computer Science Department.
        !        University of Tennessee, Knoxville. 1993.

        ! !INPUT PARAMETERS:
        integer (i4), intent(in) :: n
        integer (i4), intent(in) :: iptrHalo


        type (Matrix_t) :: A
        real (r8), dimension(:), intent(in) :: &
            B                         ! right hand side of linear system

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), dimension(:), &
        intent(inout) :: X  ! on input,  an initial guess for the solution
                            ! on output, solution of the linear system

        real(r8), dimension(:), intent(in) :: mask
        real(r8), dimension(:), intent(in) :: Minv2

        integer(i4) :: solv_sum_iters
        real(r8) :: rms_residual

        !-----------------------------------------------------------------------
        character (char_len) :: &
            noconvrg           ! error message for no convergence

        integer (i4) :: &
            i,j,m,             &! local iteration counter
            iblock              ! local block     counter

        real (r8) :: & ! scalar results
            cg_alpha, cg_beta, cg_sigma, cg_delta, &
            cg_rho_old, cg_rho, rr, rr0, rr_local

        real (r8), dimension(max_linear) :: &
            R,                  &! residual (b-Ax)
            S,                  &! conjugate direction vector
            Q,Z,AZ,WORK0,       &! various cg intermediate results
            A0R                  ! diagonal preconditioner

        real (r8), dimension(max_linear,2) :: &
            WORKN              ! WORK array

        real (r8), dimension(3) :: &
            sumN               ! global sum results for multiple arrays

        real (r8), dimension(3) :: sumN_local

        real (r8) :: sumN1,sumN2,sumN3

        real (r8) :: stmp,qtmp

        integer  :: nActive
        integer  :: rem,unroll,n2

        !-----------------------------------------------------------------------
        !
        !  initialize some scalars
        !
        !-----------------------------------------------------------------------
        nActive = iptrHalo-1

        cg_rho = c1
        solv_sum_iters = solv_max_iters

        !-----------------------------------------------------------------------
        !
        !  compute initial residual and initialize other arrays
        !
        !-----------------------------------------------------------------------
        call matvec(n,A,S,X)

        rr_local = 0.0d0
        do i=1,n
            R(i) = B(i) - S(i)
            rr_local = rr_local + R(i)*R(i)
        enddo
        rr0 = sqrt(global_sum(rr_local))

        call update_halo(R)

        !-----------------------------------------------------------------------
        !
        !    take one pass of standard algorithm
        !
        !-----------------------------------------------------------------------
        sumN1=c0
        do i=1,nActive
            Z(i) = Minv2(i)*R(i)
            sumN1 = sumN1 + R(i)*Z(i)
            S(i) = Z(i)
        enddo

        ! ===========================================
        ! Apply the preconditioner to the halo region
        ! ===========================================
        do i=iptrHalo,n
            Z(i) = Minv2(i)*R(i)
            S(i) = Z(i)
        enddo

        sumN2=c0
        call matvec(n,A,Q,S)
        do i=1,nActive
            sumN2 = sumN2 + S(i)*Q(i)
        enddo

        sumN_local(1) = sumN1
        sumN_local(2) = sumN2
        sumN_local(3) = c0
        sumN = global_sum(sumN_local)

        call update_halo(Q)

        cg_rho_old = sumN(1) !(r,PCr)
        cg_sigma   = sumN(2) !(s,As)
        cg_alpha   = cg_rho_old/cg_sigma

        do i=1,n
            X(i) = X(i) + cg_alpha*S(i)
            R(i) = R(i) - cg_alpha*Q(i)
        enddo

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
            sumN1=c0
            sumN3=c0
            do i=1,nActive
                Z(i) = Minv2(i)*R(i)
                sumN1 = sumN1 + R(i)*Z(i)
                sumN3 = sumN3 + R(i)*R(i)
            enddo
            do i=iptrHalo,n
                Z(i) = Minv2(i)*R(i)
            enddo

            call matvec(n,A,AZ,Z)

            sumN2=c0
            do i=1,nActive
                sumN2 = sumN2 + AZ(i)*Z(i)
            enddo

            call update_halo(AZ)

            sumN_local(1)=sumN1
            sumN_local(2)=sumN2
            sumN_local(3)=sumN3
            sumN = global_sum(sumN_local)

            cg_rho     = sumN(1)   ! (r,(PC)r)
            cg_delta   = sumN(2)   ! (A (PC)r,(PC)r)
            rr         = sumN(3)
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
            do i=1,n
                stmp = Z(i) + cg_beta *S(i)
                qtmp= AZ(i) + cg_beta *Q(i)
                X(i) = X(i) + cg_alpha *stmp
                R(i) = R(i) - cg_alpha *qtmp
                S(i) = stmp
                Q(i) = qtmp
            enddo
        end do iter_loop
        rms_residual = sqrt(rr*resid_norm)
    end subroutine pcg_chrongear_linear
end module solvers
