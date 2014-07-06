!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!BOP
! !ROUTINE: CGinit
! !INTERFACE:
!>
!! @defgroup CGinit
!<
!!
!>
!! @ingroup CGinit
!<
 program CGinit

! !DESCRIPTION:
!  This program generates a series of input files for the CGpop miniapp 
!
! !REVISION HISTORY:
!  CVS:$Id: POP.F90,v 1.4 2007/01/25 23:34:14 dennis Exp $
!  CVS:$Name:  $
! !USES:

  use kinds_mod, only: i4, r8
  use blocks, only: nx_block,ny_block, all_blocks, AppMD_t, nblocks_tot, &
	destroy_blocks
  use boundary, only: destroy_boundary 
  use domain_size, only: block_size_x, block_size_y, nx_global,ny_global, max_blocks_tropic
  use communicate, only: init_communicate, my_task, master_task, MPI_COMM_OCN
  use domain, only: nblocks_tropic, blocks_tropic, distrb_tropic, &
     init_domain_blocks, init_domain_distribution, &
     destroy_domain_blocks, destroy_domain_distribution
  use reductions, only: global_sum
  use solvers, only: A0, AN, AE, ANE, RCALCT_B, solv_sum_iters,esolver, init_solvers
  use constants, only: init_constants, boundary_exchange_algorithm, &
	ALG_MPI2S_1D, ALG_MPI2S_2D, ALG_CAF_MULTI_PULL_1D, ALG_CAF_SINGLE_PULL_1D, &
        ALG_CAF_SINGLE_PUSH_1D, ALG_CAF_SINGLE_PUSH_2D, ALG_CAF_SINGLE_PULL_2D
  use io_serial, only: read_garray, write_AppMD, read_AppMD
  use timers, only: get_timer,timer_start,timer_stop,timer_print_all,init_timers

   implicit none

   real(r8), allocatable, dimension(:,:,:) :: RHS, PRESSI, PRESS, PRESSF


   real(r8) :: sum_diff,gdiff

   integer(i4),allocatable :: KMT(:,:)

   character(len=80) :: ifname,ofname

   integer(i4) :: timer_mpi2s_solver_1D, timer_mpi2s_solver_2D

   type (AppMD_t), pointer :: all_blocks_read(:)

   integer (i4) :: i

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: nscan,nstep
   integer (i4) :: iblock, ierr
   integer (i4) :: n

  character(len=10) tilesize

   integer(i4) :: timer_esolver
   integer(i4) :: timer_solver_init
 
   integer(i4), pointer :: reorder(:)
   integer(i4), pointer :: reorder_read(:)


   logical, parameter :: &
       DO_MPI2S_1D            = .TRUE., &
       DO_MPI2S_2D            = .TRUE., &
       DO_CAF_MULTI_PULL_1D   = .TRUE., &
       DO_CAF_SINGLE_PULL_1D  = .TRUE., &
       DO_CAF_SINGLE_PUSH_1D  = .TRUE., &
       DO_CAF_SINGLE_PUSH_2D  = .TRUE., &
       DO_CAF_SINGLE_PULL_2D  = .TRUE.

    integer :: is

    integer(i4), parameter :: nsize = 8
    integer, dimension(nsize) :: bsx,bsy,maxb
    data bsx  /180, 120,  90,  60,  48,    36,   24,   18/
    data bsy  /120,  80,  60,  40,  32,    24,   16,   12/
    data maxb /400, 900,1600,3600,5625, 10000,22500,40000/

!    integer(i4), parameter :: nsize = 1
!    integer, dimension(nsize) :: bsx,bsy,maxb
!    data bsx  /180/
!    data bsy  /120/
!    data maxb /400/

!    integer(i4), parameter :: nsize = 6
!    integer, dimension(nsize) :: bsx,bsy,maxb
!    data bsx  /  80,  40,  30    20,   15,   12/
!    data bsy  /  80,  40,  30,   20,   15,   12/
!    data maxb /1350,5400,9600,21600,38400,60000/

!    integer(i4), parameter :: nsize = 6
!    integer, dimension(nsize) :: bsx,bsy,maxb
!    data bsx  /  80,  40,  30    20,   15,   12/
!    data bsy  /  80,  40,  30,   20,   15,   12/
!    data maxb /1350,5400,9600,21600,38400,60000/
!-----------------------------------------------------------------------
!  initialize message-passing or other communication protocol
!-----------------------------------------------------------------------

   call init_communicate
   print *,'cginit: point #1'

!-----------------------------------------------------------------------
!  initialize constants and i/o stuff
!-----------------------------------------------------------------------

   call init_constants
   print *,'cginit: point #2'

do is=1,nsize
   block_size_x=bsx(is)
   block_size_y=bsy(is)
   max_blocks_tropic = maxb(is)
!-----------------------------------------------------------------------
!  initialize domain and grid
!-----------------------------------------------------------------------
   call init_domain_blocks
   print *,'cginit: after call to init_domain_blocks'
   print *,'cginit: point #2.1'



   ifname = '../data/cgpopState.nc'

   print *,'nx_global,ny_global: ',nx_global,ny_global 
   allocate(KMT(nx_global,ny_global))
   call read_garray(ifname,'KMT',KMT)
   print *,'cginit: point #3'

   allocate(reorder(nblocks_tot))

   call init_domain_distribution(KMT,reorder)
   deallocate(KMT)

   print *,'cginit: point #4'

  
  if(block_size_x > 99 .and. block_size_y > 99) then 
      write(tilesize,101) block_size_x,block_size_y
  elseif(block_size_x < 100 .and. block_size_y > 99) then 
      write(tilesize,102) block_size_x,block_size_y
  elseif(block_size_x > 99 .and. block_size_y < 100) then 
      write(tilesize,103) block_size_x,block_size_y
  elseif(block_size_x < 100 .and. block_size_y < 100) then 
      write(tilesize,104) block_size_x,block_size_y
  endif

 101 format(i3,'x',i3)
 102 format(i2,'x',i3)
 103 format(i3,'x',i2)
 104 format(i2,'x',i2)
   
   ofname = '../data/cgpoptile_'//TRIM(tilesize)//'.nc'

   call write_AppMD(ofname,all_blocks,reorder)
   print *,'cginit: point #5'
!   call  read_AppMD(ofname,all_blocks_read,reorder_read)
!   deallocate(all_blocks_read,reorder_read)
   print *,'cginit: point #6'


!-----------------------------------------------------------------------
!  initialize timers and additional communication routines
!-----------------------------------------------------------------------

   ifname = '../data/cgpopState.nc'
!-----------------------------------------------------------------------
! Initialize the solver
!-----------------------------------------------------------------------
   allocate(RHS(nx_block,ny_block,max_blocks_tropic))
   allocate(PRESSI(nx_block,ny_block,max_blocks_tropic))
   allocate(PRESSF(nx_block,ny_block,max_blocks_tropic))
   allocate(PRESS(nx_block,ny_block,max_blocks_tropic))

   call init_solvers(ifname,ofname,RHS,PRESSI,PRESSF)
   print *,'cginit: point #7'

   deallocate(RHS,PRESSI,PRESSF,PRESS)
   call destroy_domain_blocks()
   call destroy_domain_distribution()

enddo


 end program CGinit
