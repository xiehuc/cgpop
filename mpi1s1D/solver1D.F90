module solver1D
    use kinds_mod, only: i4,r8, char_len
    use communicate, only: my_task, master_task
    use simple_domain, only: distrb_tropic
    use constants, only: c0, c1
    use matrix_mod, only: Matrix_t,matvec
    use reductions, only: global_sum
    use linear, only: max_linear,update_halo
    use exit_mod, only: sigAbort, exit_POP
    use timers, only: timer_print_all

    implicit none
    private

    integer(i4) :: solv_max_iters
    integer(i4) :: solv_ncheck
    real(r8) :: solv_convrg
    real(r8) :: resid_norm

    public :: pcg_chrongear_linear

    public :: init_solvers1D

  contains 

    subroutine init_solvers1D(maxiters, ncheck, convrg, norm)
        integer(i4) :: maxiters,ncheck
        real(r8) :: convrg,norm

        solv_max_iters=maxiters
        solv_ncheck = ncheck
        solv_convrg = convrg
        resid_norm = norm
    end subroutine init_solvers1D


    !***********************************************************************
    !>
    !! This routine implements the Chronopoulos-Gear conjugate-gradient
    !! solver with preconditioner for solving the linear system $Ax=b$.
    !! It is a rearranged conjugate gradient solver that reduces the
    !! number of inner products per iteration from two to one (not
    !! counting convergence check). Both the operator $A$ and
    !! preconditioner are nine-point stencils. If no preconditioner has
    !! been supplied, a diagonal preconditioner is applied.  Convergence
    !! is checked every {\em ncheck} steps.
    !<
    subroutine pcg_chrongear_linear(n,iptrHalo,A,X,B,Minv2,mask, &
        solv_sum_iters, rms_residual)
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
        real (r8), dimension(:), &
            intent(in) :: B     ! right hand side of linear system

        ! !INPUT/OUTPUT PARAMETERS:
        real (r8), dimension(:), &
        intent(inout) ::  X   ! on input,  an initial guess for the solution
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

        !-------------------------------------------------------------------
        !
        !  initialize some scalars
        !
        !-------------------------------------------------------------------
        nActive = iptrHalo-1

        cg_rho = c1
        solv_sum_iters = solv_max_iters

        !-------------------------------------------------------------------
        !
        !  compute initial residual and initialize other arrays
        !
        !-------------------------------------------------------------------

        call matvec(n,A,S,X)

        rr_local = 0.0d0
        do i=1,n
            R(i) = B(i) - S(i)
            rr_local = rr_local + R(i)*R(i)
        enddo
        rr0 = sqrt(global_sum(rr_local))
        !   print *,'initial residual: ',rr0

        call update_halo(R)

        !-------------------------------------------------------------------
        !
        !    take one pass of standard algorithm
        !
        !-------------------------------------------------------------------
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

            !---------------------------------------------------------------
            !
            !     test for convergence if it is time
            !
            !---------------------------------------------------------------
            if(mod(m,solv_ncheck) == 0) then
                if (sqrt(rr)/rr0 < solv_convrg) then 
                    solv_sum_iters = m
                    exit iter_loop
                endif
            endif

        end do iter_loop
        rms_residual = sqrt(rr*resid_norm)
    end subroutine pcg_chrongear_linear
end module solver1D
