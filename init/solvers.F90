!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 module solvers
!BOP
! !MODULE: solvers
!
! !DESCRIPTION:
!  This module contains routines and operators for solving the elliptic
!  system for surface pressure in the barotropic mode.
!
! !REVISION HISTORY:
!  CVS:$Id: solvers.F90,v 1.43 2006/05/30 16:02:00 dennis Exp $
!  CVS:$Name:  $

! !USES:

   use kinds_mod, only: i4, r8, char_len, log_kind
   use simple_type, only : distrb
   use blocks, only: nx_block,ny_block, nblocks_tot, get_block_parameter
   use communicate, only: my_task, master_task
   use domain, only: distrb_tropic, nblocks_tropic, blocks_tropic, bndy_tropic, ConvertG2T
   use domain_size, only: max_blocks_tropic, nx_global, ny_global
   use constants, only: c0, c1, field_type_scalar, field_loc_center, p25, eps, &
       blank_fmt, delim_fmt, boundary_exchange_algorithm, ALG_MPI2S_1D, ALG_MPI2S_2D, &
	ALG_CAF_MULTI_PULL_1D, ALG_CAF_SINGLE_PULL_1D,ALG_CAF_SINGLE_PUSH_1D, &
	ALG_CAF_SINGLE_PUSH_2D, ALG_CAF_SINGLE_PULL_2D

   use boundary, only: boundary_2d
!   use global_reductions, only: global_sum
   use reductions, only: global_sum
!   use solver1D, only: pcg_chrongear_linear, init_solvers1D
!   use solver2D, only: pcg_chrongear, init_solvers2D
   use broadcast, only: broadcast_scalar
   use timers, only: get_timer, timer_start, timer_stop
   use IOUnitsMod, only: stdout, nml_filename, nmlin
   use exit_mod, only: sigAbort, exit_POP
   use linear, only: linearmask, nActive, nTotal,  max_linear, convert2dtolinear, &
	update_halo, convertlinearto2d, initdof
   use matrix_mod, only: convertstencil, matvec, A, initMatrix
   use io_serial, only: read_garray, write_tiles   
!   use io, only: read_solverio

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

    integer (i4), save :: &
       timer_tolinear,           &! timer for the conversion to linear form
       timer_to2D,               &! time for the conversion to 2D form
       timer_matvec,timer_gsum,timer_halo, &
       timer_initDOF,            &!  Time to setup the DOF based communicator 
       timer_solver               ! time to perform the CG solver

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  other operator and preconditioner weights for barotropic operator
!
!-----------------------------------------------------------------------

   real (r8), public, dimension(:,:,:),allocatable :: & 
      A0,AN,AE,ANE,         &! barotropic (9pt) operator coefficients
      RCALCT_B               ! land mask in barotropic distribution 


!-----------------------------------------------------------------------
!
!  scalar convergence-related variables
!
!-----------------------------------------------------------------------

   logical (log_kind) :: &
      lprecond            ! true if computed preconditioner to be used

   real (r8) ::          &
      solv_convrg,       &! convergence error criterion
      sor,               &! for jacobi solver
      resid_norm          ! residual normalization

   integer (i4), parameter :: &
      solv_pcg = 1,      &! predefined solver types
      solv_cgr = 2,      &
      solv_jac = 3,      &
      solv_cg1 = 4,	 &
      solv_lcg1 = 5,     &
      solv_lpcg = 6

   integer (i4) :: &
      solv_itype,        &! integer solver method (1=pcg, 2=cgr, 3=jac)
      solv_max_iters,    &! max number of iterations for solver
      solv_ncheck         ! check convergence every ncheck iterations

    logical :: solv_linear        

!EOC
!***********************************************************************

 contains

  subroutine esolver(RHS,PRESS)

     real(r8), dimension(nx_block,ny_block,max_blocks_tropic) ::  RHS,PRESS
    
     real (r8), dimension(max_linear) :: &
       X_linear,        &! The initial guess and solution in Linear form
       B_linear,        &! The RHS in Linear form
      Mask_linear,      &! Mask for non-valid points
      Minv_linear        ! Diagonal preconditioner

     integer(i4) :: i

     select case(boundary_exchange_algorithm)
       case(ALG_MPI2S_1D)
	 solv_linear = .true.
       case(ALG_CAF_MULTI_PULL_1D)
	 solv_linear = .true.
       case(ALG_CAF_SINGLE_PULL_1D)
	 solv_linear = .true.
       case(ALG_CAF_SINGLE_PUSH_1D)
	 solv_linear = .true.
       case(ALG_MPI2S_2D)
	 solv_linear = .false.
       case(ALG_CAF_SINGLE_PUSH_2D)
	 solv_linear = .false.
       case(ALG_CAF_SINGLE_PULL_2D)
	 solv_linear = .false.
     end select
       
    if(solv_linear) then 
      if(my_task == master_task) print *,'esolver: using the linear datastructure version'
      !==============================
      ! The initial guess
      !==============================
!      call Convert2DtoLinear(X_linear,PRESS)
!      call update_halo(X_linear)

      !==============================
      ! the RHS
      !==============================
!      call Convert2DtoLinear(B_linear,RHS)
!      call update_halo(B_linear)

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
!      call update_halo(Minv_linear)

      !==============================
      ! The A matrix
      !==============================
      call ConvertStencil(A,A0,AN,AE,ANE)

      if (my_task < distrb_tropic%nprocs) then
!         call pcg_chrongear_linear(nTotal,nActive+1, &
!                A,X_linear,B_linear,Minv_linear,LinearMask,solv_sum_iters,rms_residual)
      endif

      call ConvertLinearto2D(PRESS,X_linear)

      ! ==================================
      ! Update the Halo with the solution
      ! ==================================
      call boundary_2d(PRESS,bndy_tropic, &
                field_loc_center,field_type_scalar)
    else 
!      call pcg_chrongear(PRESS,RHS) ! Chron-Gear pcg solver
    endif
!-----------------------------------------------------------------------


 end subroutine esolver
!***********************************************************************
!BOP
! !IROUTINE: init_solvers
! !INTERFACE:

 subroutine init_solvers(ifname,ofname,RHS,PRESSI,PRESSF)

! !DESCRIPTION:
!  This routine initializes choice of solver, calculates the 
!  coefficients of the 9-point stencils for the barotropic operator and
!  reads in a preconditioner if requested.
!
! !REVISION HISTORY:
!  same as module

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables:
!
!       {X,Y}{NE,SE,NW,SW} = contribution to {ne,se,nw,sw} coefficients 
!         from {x,y} components of divergence
!       HU = depth at U points
!
!-----------------------------------------------------------------------

!   real (r8), dimension(nx_block,ny_block,max_blocks_tropic) :: RHS, PRESSI, PRESSF
   real (r8),dimension(:,:,:) :: RHS,PRESSI,PRESSF
   character(len=80) :: ifname
   character(len=80) :: ofname

   character (char_len) :: &
      solv_type,           &! user choice of solver method
      precond_file          ! file containing preconditioner

   namelist /solver_nml/ solv_convrg, solv_max_iters, solv_ncheck, &
                         lprecond, solv_type, precond_file

   integer (i4) :: &
      i,j,n,             &! dummy counter
      iblock,            &! block counter
      ncheck,            &! scalar for checking PC/mask compatibility
      nu,                &! I/O unit number and status flag
      nml_error           ! namelist i/o error flag

   integer (i4),dimension(:), allocatable :: &
      icheck              ! check for PC/mask compatibility in block

  real(r8), dimension(nx_block,ny_block,max_blocks_tropic) :: &
     WORK0_read,RCALCT_read, RTMPT

  integer(i4), dimension(nx_global,ny_global) :: ITMPG
  integer(i4), dimension(nx_block,ny_block,max_blocks_tropic) :: ITMPT

  real(r8), allocatable :: RTMPT2(:,:,:)

   logical (log_kind) :: &
      mlandne, mlandnw, mlandse, mlandsw ! land mask at nbr points

   real(r8), allocatable :: RMASKG(:,:),RTMPG(:,:)

   logical,parameter :: Debug = .FALSE.
   integer(i4) :: nstep

   real(r8) :: sum_work0,sum_rcalct

   integer (i4), dimension(:),pointer :: iglob, jglob
   integer (i4) :: jb, je, ib, ie
   integer (i4) :: gbid


!JMD   call initMatrix()
   allocate(A0(nx_block,ny_block,max_blocks_tropic))
   allocate(AN(nx_block,ny_block,max_blocks_tropic))
   allocate(AE(nx_block,ny_block,max_blocks_tropic))
   allocate(ANE(nx_block,ny_block,max_blocks_tropic))
   allocate(RCALCT_B(nx_block,ny_block,max_blocks_tropic))

   allocate(RTMPG(nx_global,ny_global),RMASKG(nx_global,ny_global))
    
   allocate(RTMPT2(nx_block,ny_block,max_blocks_tropic))
!   allocate(PRESSI(nx_block,ny_block,max_blocks_tropic))
!   allocate(PRESSF(nx_block,ny_block,max_blocks_tropic))
!-----------------------------------------------------------------------
!
!  read solver choice and solver constants from namelist input
!  (namelist input file opened in initial.F)
!
!-----------------------------------------------------------------------

   solv_convrg    = eps
   solv_max_iters = 1000
   solv_ncheck    = 10
!JMD   solv_linear    = .FALSE.
   lprecond       = .false.
   solv_type      = 'pcg'
   precond_file   = 'empty'

   if (my_task == master_task) then
      open (nmlin, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nmlin, nml=solver_nml,iostat=nml_error)
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
         write(stdout,'(a38)') '  Chron-Gear preconditioned conjg grad'
         solv_itype = solv_cg1
      case('lcg1')
         write(stdout,'(a46)') ' Linear Chron-Gear preconditioner conjg grad'
         solv_itype = solv_lcg1
         solv_linear = .TRUE.
      case('pcg')
         write(stdout,'(a35)') '  Preconditioned Conjugate Gradient'
         solv_itype = solv_pcg
      case('lpcg')
         write(stdout,'(a42)') ' Linear Preconditioned Conjugate Gradient'
         solv_itype = solv_lpcg
         solv_linear = .TRUE.
      case default
         solv_itype = -1000
      end select

      write(stdout,'(a28,1pe12.5)') ' Solver converged for err < ', &
                                      solv_convrg
      write(stdout,'(a29,i6)') ' Solver maximum iterations = ', &
                                 solv_max_iters
      write(stdout,'(a35,i6,a11)') ' Solver convergence checked every ',&
                                     solv_ncheck, ' iterations'

   endif

   solv_linear = .true.
   call broadcast_scalar(solv_convrg,    master_task)
   call broadcast_scalar(solv_max_iters, master_task)
   call broadcast_scalar(solv_ncheck,    master_task)
   call broadcast_scalar(lprecond,       master_task)
   call broadcast_scalar(solv_itype,     master_task)
   call broadcast_scalar(precond_file,   master_task)
   call broadcast_scalar(solv_linear,    master_task)

   if (solv_itype == -1000) then
      call exit_POP(sigAbort, &
                 'unknown solver type: must be cg1, lcg1, pcg, lpcg')
   endif

!-----------------------------------------------------------------------
!
!  set sor for jacobi solver
!
!-----------------------------------------------------------------------

   sor = p25     ! should be < 1/2

!-----------------------------------------------------------------------
!
!  compute nine point operator coefficients: compute on baroclinic
!  decomposition first where grid info defined and redistribute
!  to barotropic distribution
!  leave A0,AC in baroclinic distribution to facilitate easy
!  time-dependent changes in barotropic routine
!
!-----------------------------------------------------------------------

   A0     = 0.0_r8
   AN     = 0.0_r8
   AE     = 0.0_r8
   ANE    = 0.0_r8
   RHS    = 0.0_r8
   PRESSI = 0.0_r8
   PRESSF = 0.0_r8
   RCALCT_read = 0.0_r8
   WORK0_read = 0.0_r8
!   fname = 'cgpop_v3.nc'
!   call read_solverio(fname,nstep,A0,AN,AE,ANE,RHS,PRESSI,PRESSF,RCALCT_read,WORK0_read)
    call read_garray(ifname,'RCALCT',RMASKG)
#if 1
    call ConvertG2T(RMASKG,RCALCT_read)
#else
    allocate(iglob(nx_block),jglob(ny_block))
    do n=1,nblocks_tot
       call get_block_parameter(n,jb=jb,je=je,ib=ib,ie=ie,i_glob=iglob,j_glob=jglob)
       do j=jb,je
          if(jglob(j) > 0) then 
            do i=ib,ie
	       if(iglob(i) >0) then 
		  RCALCT_read(i,j,n) = RMASKG(iglob(i),jglob(j))
	       endif
	    enddo
	  endif
       enddo
    enddo
    deallocate(iglob,jglob)
#endif

!-----------------------------------------------------------------------
!
!  calculate normalization constant (darea,darea) for rms_residual
!  in cgr routine.
!
!-----------------------------------------------------------------------

!JMD   resid_norm = c1/global_sum(WORK0_read, distrb_tropic, field_loc_center, RCALCT_read)

!JMD   solv_convrg = solv_convrg**2/resid_norm

   RCALCT_B = RCALCT_read

!-----------------------------------------------------------------------
!
!  setup preconditioner if required
!
!-----------------------------------------------------------------------

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

!JMD   if(solv_linear) then 

     strip = .TRUE.  ! if true strips outs land points
!JMD     call get_timer(timer_initDOF,'initDOF',1,distrb_tropic%nprocs)
!JMD     call timer_start(timer_initDOF)
     call initDOF(ofname,strip,RCALCT_B)

   
    !-------------------------------------------------
    ! add all the state information to the output flie
    !-------------------------------------------------
    RTMPT2(:,:,:)=0
    do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       RTMPT2(:,:,gbid) = RCALCT_B(:,:,n)
    enddo
    call write_tiles(ofname,'RCALCT',RTMPT2)

 
    call read_garray(ifname,'PRESS_I',RTMPG)
    call ConvertG2T(RTMPG,RTMPT)
    RTMPT2(:,:,:)=0
    do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       RTMPT2(:,:,gbid) = RTMPT(:,:,n)
    enddo
    call write_tiles(ofname,'PRESS_I',RTMPT2)


    !---------------------------------------------------------- 
    ! convert the RHS array from global to ordered tiled output 
    !---------------------------------------------------------- 
    call read_garray(ifname,'RHS',RTMPG)
    call ConvertG2T(RTMPG,RTMPT)
    RTMPT2(:,:,:)=0
    do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       RTMPT2(:,:,gbid) = RTMPT(:,:,n)
    enddo
    call write_tiles(ofname,'RHS',RTMPT2)


    call read_garray(ifname,'PRESS_F',RTMPG)
    call ConvertG2T(RTMPG,RTMPT)
    RTMPT2(:,:,:)=0
    do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       RTMPT2(:,:,gbid) = RTMPT(:,:,n)
    enddo
    call write_tiles(ofname,'PRESS_F',RTMPT2)

    call read_garray(ifname,'A0',RTMPG)
    call ConvertG2T(RTMPG,RTMPT)
    RTMPT2(:,:,:)=0
    do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       RTMPT2(:,:,gbid) = RTMPT(:,:,n)
    enddo
    call write_tiles(ofname,'A0',RTMPT2)

    call read_garray(ifname,'AN',RTMPG)
    call ConvertG2T(RTMPG,RTMPT)
    RTMPT2(:,:,:)=0
    do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       RTMPT2(:,:,gbid) = RTMPT(:,:,n)
    enddo
    call write_tiles(ofname,'AN',RTMPT2)

    call read_garray(ifname,'AE',RTMPG)
    call ConvertG2T(RTMPG,RTMPT)
    RTMPT2(:,:,:)=0
    do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       RTMPT2(:,:,gbid) = RTMPT(:,:,n)
    enddo
    call write_tiles(ofname,'AE',RTMPT2)

    call read_garray(ifname,'ANE',RTMPG)
    call ConvertG2T(RTMPG,RTMPT)
    RTMPT2(:,:,:)=0
    do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       RTMPT2(:,:,gbid) = RTMPT(:,:,n)
    enddo
    call write_tiles(ofname,'ANE',RTMPT2)

    call read_garray(ifname,'TAREASQ',RTMPG)
    call ConvertG2T(RTMPG,RTMPT)
    RTMPT2(:,:,:)=0
    do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       RTMPT2(:,:,gbid) = RTMPT(:,:,n)
    enddo
    call write_tiles(ofname,'TAREASQ',RTMPT2)

    
!JMD    call read_garray(ifname,'KMT',ITMPG)
!JMD    call ConvertG2T(ITMPG,ITMPT)
!JMD    call write_tiles(ofname,'KMT',ITMPT)

!JMD    stop 'after call to initDOF'
!JMD     call timer_stop(timer_initDOF)

     !====================================
     ! Create the full matrix
     !====================================
!JMD     A%maxNZ = 9*nTotal
!JMD     A%n = nActive
!JMD#ifdef _USEESSL
!JMD     A%nz_essl=9
!JMD#endif

     !================================================
     ! Create the timers to measure conversion to and
     ! from the linear in the solver
     !================================================
!JMD     call get_timer(timer_tolinear,'CONVERT TO LINEAR',1, &
!JMD                        distrb_tropic%nprocs)

!JMD     call get_timer(timer_to2D,'CONVERT TO 2D',1, &
!JMD                        distrb_tropic%nprocs)

!JMD     call init_solvers1D(solv_max_iters,solv_ncheck, solv_convrg, resid_norm)
!JMD   else 

!JMD     call init_solvers2D(solv_max_iters,solv_ncheck, solv_convrg, resid_norm, &
!JMD	A0,AN,AE,ANE,RCALCT_B)

!JMD   endif

!JMD     call get_timer(timer_solver,'SOLVER',1, &
!JMD                        distrb_tropic%nprocs)

!EOC
   deallocate(A0,AN,AE,ANE,RCALCT_B)
   deallocate(RTMPG,RMASKG)
   deallocate(RTMPT2)

 end subroutine init_solvers

!***********************************************************************

 end module solvers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
