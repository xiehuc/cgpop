!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module domain

!BOP
! !MODULE: domain
!
! !DESCRIPTION:
!  This module contains the model domain and routines for initializing
!  the domain.  It also initializes the decompositions and
!  distributions across processors/threads by calling relevent
!  routines in the block, distribution modules.
!
! !REVISION HISTORY:
!  CVS:$Id: domain.F90,v 1.7 2007/01/29 17:04:15 dennis Exp $
!  CVS:$Name:  $

! !USES:

   use kinds_mod, only: i4, r8, char_len, log_kind
   use simple_type, only: distrb
   use constants, only: blank_fmt, delim_fmt
   use communicate, only: my_task, master_task, get_num_procs, create_communicator
   use blocks, only: nx_block,ny_block,nblocks_tot, destroy_blocks, &
	nghost, create_blocks, get_block_parameter, set_block_parameter
   use broadcast, only: broadcast_scalar
   use boundary, only: bndy, create_boundary, destroy_boundary
   use distribution, only: samedistribution, create_distribution, &
	create_local_block_ids, destroy_distribution
   use exit_mod, only: sigAbort, exit_pop
   use IOUnitsMod, only: stdout,nmlin,nml_filename
   use domain_size, only: nx_global,ny_global,max_blocks_tropic
   use reductions, only: global_maxval
   use io_serial, only: write_tiles, read_tiles

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS

   public  :: init_domain_blocks ,&
              init_domain_distribution
   public  :: destroy_domain_blocks, &
	      destroy_domain_distribution

   public :: ConvertG2T

   interface ConvertG2T
      module procedure ConvertG2T_int
      module procedure ConvertG2T_dbl
   end interface

   public :: write_ordered_tiles
   interface write_ordered_tiles
      module procedure write_ordered_tiles_int
   end interface

   public :: read_ordered_tiles
   interface read_ordered_tiles
      module procedure read_ordered_tiles_int
   end interface

! !PUBLIC DATA MEMBERS:


   integer (i4), public :: &
      nblocks_tropic     ! actuale number of blocks in tropic distribution

   integer (i4), dimension(:), pointer, public :: &
      blocks_tropic      ! block ids for local blocks in barotropic dist

   type (distrb), public :: & !  block distribution info
      distrb_tropic      ! block distribution for barotropic part

   type (bndy), public :: &!  ghost cell update info
      bndy_tropic          ! block distribution for barotropic part

!EOP
!BOC
!-----------------------------------------------------------------------
!
!   module private variables - for the most part these appear as
!   module variables to facilitate sharing info between init_domain1
!   and init_domain2.
!
!-----------------------------------------------------------------------

    character (char_len) ::      &
       tropic_distribution_type, &!    blocks in each case
       ew_boundary_type,         &! type of domain bndy in each logical
       ns_boundary_type           !    direction (ew is i, ns is j)

    integer (i4) :: &! decomposition info
       nprocs_tropic       ! num of processors in barotropic dist

    integer(i4),public :: timer_detail
!EOC
!***********************************************************************

 contains

  subroutine destroy_domain_blocks()

     call destroy_blocks()

  end subroutine destroy_domain_blocks


!***********************************************************************
!BOP
! !IROUTINE: init_domain_blocks
! !INTERFACE:

 subroutine init_domain_blocks

! !DESCRIPTION:
!  This routine reads in domain information and calls the routine
!  to set up the block decomposition.
!
! !REVISION HISTORY:
!  same as module
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (i4) :: &
      nml_error          ! namelist read error flag

!----------------------------------------------------------------------
!
!  input namelists
!
!----------------------------------------------------------------------

   namelist /domain_nml/ nprocs_tropic, &
                         tropic_distribution_type, &
                         ew_boundary_type,         &
                         ns_boundary_type,         &
                         timer_detail

!----------------------------------------------------------------------
!
!  read domain information from namelist input
!
!----------------------------------------------------------------------

   nprocs_tropic = -1
   timer_detail  = 0
   tropic_distribution_type = 'cartesian'
   ew_boundary_type = 'cyclic'
   ns_boundary_type = 'closed'

   if (my_task == master_task) then
      open (nmlin, file=nml_filename, status='old',iostat=nml_error)
      if (nml_error /= 0) then
         nml_error = -1
      else
         nml_error =  1
      endif
      do while (nml_error > 0)
         read(nmlin, nml=domain_nml,iostat=nml_error)
      end do
      if (nml_error == 0) close(nmlin)
   endif

   call broadcast_scalar(nml_error, master_task)
   if (nml_error /= 0) then
      call exit_POP(sigAbort,'ERROR reading domain_nml')
   endif

   call broadcast_scalar(nprocs_tropic,            master_task)
   call broadcast_scalar(tropic_distribution_type, master_task)
   call broadcast_scalar(ew_boundary_type,         master_task)
   call broadcast_scalar(ns_boundary_type,         master_task)
   call broadcast_scalar(timer_detail,             master_task)

!----------------------------------------------------------------------
!
!  perform some basic checks on domain
!
!----------------------------------------------------------------------

   if (nx_global < 1 .or. ny_global < 1) then
      !***
      !*** domain size zero or negative
      !***
      call exit_POP(sigAbort,'Invalid domain: size < 1') ! no domain
   else if (nprocs_tropic /= get_num_procs()) then
      !***
      !*** input nprocs does not match system (eg MPI) request
      !***
      print *,'input_domain_blocks: nprocs_tropic, get_num_procs: ', nprocs_tropic, get_num_procs() 
      call exit_POP(sigAbort,'Input nprocs not same as system request')
   else if (nghost < 2) then
      !***
      !*** must have at least 2 layers of ghost cells
      !***
      call exit_POP(sigAbort,'Not enough ghost cells allocated')
   endif

!----------------------------------------------------------------------
!
!  compute block decomposition and details
!
!----------------------------------------------------------------------

   call create_blocks(nx_global, ny_global, trim(ew_boundary_type), &
                                            trim(ns_boundary_type))

!----------------------------------------------------------------------
!
!  Now we need grid info before proceeding further
!  Print some domain information
!
!----------------------------------------------------------------------

   if (my_task == master_task) then
     write(stdout,delim_fmt)
     write(stdout,blank_fmt)
     write(stdout,'(a18)') 'Domain Information'
     write(stdout,blank_fmt)
     write(stdout,delim_fmt)
     write(stdout,'(a26,i6)') '  Horizontal domain: nx = ',nx_global
     write(stdout,'(a26,i6)') '                     ny = ',ny_global
     write(stdout,'(a26,i6)') '  Block size:  nx_block = ',nx_block
     write(stdout,'(a26,i6)') '               ny_block = ',ny_block
     write(stdout,'(a29,i6)') '  Processors for barotropic: ', &
                                 nprocs_tropic
     write(stdout,'(a31,a9)') '  Distribution for barotropic: ', &
                                 trim(tropic_distribution_type)
     write(stdout,'(a25,i2)') '  Number of ghost cells: ', nghost
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine init_domain_blocks

 subroutine destroy_domain_distribution 

     call destroy_distribution(distrb_tropic)

     deallocate(blocks_tropic)

     call destroy_boundary(bndy_tropic)

 end subroutine destroy_domain_distribution

!***********************************************************************
!BOP
! !IROUTINE: init_domain_distribution
! !INTERFACE:

 subroutine init_domain_distribution(KMTG, reorder)

! !DESCRIPTION:
!  This routine calls appropriate setup routines to distribute blocks
!  across processors and defines arrays with block ids for any local
!  blocks. Information about ghost cell update routines is also
!  initialized here through calls to the appropriate boundary routines.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (i4), dimension(nx_global,ny_global), intent(in) :: &
      KMTG             ! global KMT (topography) field

   integer (i4), pointer, intent(inout)  :: reorder(:)

!JMD   real (r8), dimension(nx_global, ny_global), intent(in)  :: rmaskg

   integer(i4), save :: timer_distrib
!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   character (char_len) :: outstring

   integer (i4), parameter :: &
      max_work_unit=10   ! quantize the work into values from 1,max

   integer (i4) :: &
      i,j,k,n              ,&! dummy loop indices
      count1, count2       ,&! dummy counters
      work_unit            ,&! size of quantized work unit
      nblocks_tmp          ,&! temporary value of nblocks
      nblocks_tmp_tropic   ,&! num blocks on proc for tropic
      nblocks_max_tropic     ! max blocks on proc for tropic

   integer (i4), dimension(:), allocatable :: &
      nocn               ,&! number of ocean points per block
      work_per_block       ! number of work units per block

   integer (i4), dimension(:), allocatable :: &
      noWork_per_block     ! =1 for all land blocks

   integer (i4), dimension(:), pointer :: iglob,jglob
   integer (i4) :: jb,je, ib, ie

!----------------------------------------------------------------------
!
!  estimate the amount of work per processor using the topography
!
!----------------------------------------------------------------------

   allocate(iglob(nx_block),jglob(ny_block))
   allocate(nocn(nblocks_tot))

   nocn = 0
   do n=1,nblocks_tot
      call get_block_parameter(n,jb=jb,je=je,ib=ib,ie=ie,i_glob=iglob,j_glob=jglob)
      do j=jb,je
         if (jglob(j) > 0) then
!            do i=1,nx_block
            do i=ib,ie
               if (iglob(i) > 0) then
                  if (KMTG(iglob(i),jglob(j)) > 0) nocn(n) = nocn(n) + 1
               endif
            end do
         endif
      end do

      !*** with array syntax, we actually do work on non-ocean
      !*** points, so where the block is not completely land,
      !*** reset nocn to be the full size of the block

!     if (nocn(n) > 0) nocn(n) = nx_block*ny_block
      call set_block_parameter(n,npoints=nocn(n))
   end do

   work_unit = maxval(nocn)/max_work_unit + 1

   !*** find number of work units per block

   allocate(work_per_block(nblocks_tot),noWork_per_block(nblocks_tot))

   where (nocn > 0)
     work_per_block = nocn/work_unit + 1
     noWork_per_block = 0
   elsewhere
     work_per_block = 0
     noWork_per_block = 1
   end where

   deallocate(nocn)
   deallocate(iglob,jglob)

!----------------------------------------------------------------------
!
!  determine the distribution of blocks across processors
!
!----------------------------------------------------------------------

   distrb_tropic = create_distribution(tropic_distribution_type, &
                                       nprocs_tropic, work_per_block,reorder)

   !-------------------------------------------------------
   ! Determie if the two distributions are the same or not
   !-------------------------------------------------------
   sameDistribution=.TRUE.
   
   deallocate(work_per_block,noWork_per_block)

!----------------------------------------------------------------------
!
!  allocate and determine block id for any local blocks in each
!  distribution.
!
!----------------------------------------------------------------------

   call create_local_block_ids(blocks_tropic, distrb_tropic)

   if (my_task < distrb_tropic%nprocs .and. &
       associated(blocks_tropic)) then
     nblocks_tropic = size(blocks_tropic)
   else
     nblocks_tropic = 0
   endif

   nblocks_max_tropic = global_maxval(nblocks_tropic)

   if (nblocks_max_tropic > max_blocks_tropic) then
     write(outstring,*) 'tropic blocks exceed max: increase max from ', &
	max_blocks_tropic, ' to ',nblocks_max_tropic
     call exit_POP(sigAbort,trim(outstring))
   else if (nblocks_max_tropic < max_blocks_tropic) then
     write(outstring,*) 'tropic blocks too large: decrease max to',&
                         nblocks_max_tropic
     if (my_task == master_task) write(stdout,*) trim(outstring)
     !call exit_POP(sigAbort,trim(outstring))
   endif

!----------------------------------------------------------------------
!
!  set up ghost cell updates for each distribution
!  Boundary types are cyclic, closed, or tripole
!
!----------------------------------------------------------------------

!   print *,'init_domain_distribution: before call to mpi2s_create_boundary'
   call create_boundary(bndy_tropic, distrb_tropic, blocks_tropic,&
                        trim(ns_boundary_type),     &
                        trim(ew_boundary_type),     &
                        nx_global, ny_global)


!----------------------------------------------------------------------
!EOC

 end subroutine init_domain_distribution

 subroutine read_ordered_tiles_int(fname,vname,var)

    character(len=80), intent(in) :: fname
    character(len=*), intent(in) :: vname
    integer(i4), dimension(:,:,:) :: var
    integer(i4), allocatable :: vartmp(:,:,:)
    integer(i4) :: n,n1,n2,n3,gid

     n1 = SIZE(vartmp,dim=1)
     n2 = SIZE(vartmp,dim=2)
     allocate(vartmp(n1,n2,nblocks_tot))
     var(:,:,:) = 0 
     call read_tiles(fname,vname,vartmp)
     do n=1,nblocks_tropic
	gid = blocks_tropic(n)
	var(:,:,n) = vartmp(:,:,gid)
     enddo
     deallocate(vartmp)

 end subroutine read_ordered_tiles_int

 subroutine write_ordered_tiles_int(fname,vname,var)

    character(len=80), intent(in) :: fname
    character(len=*), intent(in) :: vname
    integer(i4), dimension(:,:,:) :: var
    integer(i4), allocatable :: vartmp(:,:,:)
    integer(i4) :: n,n1,n2,n3,gid

     n1 = SIZE(vartmp,dim=1)
     n2 = SIZE(vartmp,dim=2)
     allocate(vartmp(n1,n2,nblocks_tot))
     vartmp(:,:,:)=0
     do n=1,nblocks_tropic
	gid = blocks_tropic(n)
	vartmp(:,:,gid) = var(:,:,n)        	
     enddo
     call write_tiles(fname,vname,vartmp)
     deallocate(vartmp)

 end subroutine write_ordered_tiles_int

 subroutine ConvertG2T_dbl(gvar,tvar)
   ! --------------------------------------------
   ! Convert between the global to tiled storage
   ! for a variable
   ! --------------------------------------------
   real(r8), intent(in)  :: gvar(:,:)
   real(r8), intent(out) :: tvar(:,:,:)

   !-----------------
   ! local variables
   !-----------------
   integer (i4), dimension(:), pointer :: iglob, jglob
   integer (i4) :: n,i,j,jb, je, ib, ie,gbid

   allocate(iglob(nx_block),jglob(ny_block))
   do n=1,nblocks_tropic
       gbid = blocks_tropic(n)
       call get_block_parameter(gbid,i_glob=iglob,j_glob=jglob)
       do j=1,ny_block
          if(jglob(j) > 0) then
            do i=1,nx_block
               if(iglob(i) >0) then
                  tvar(i,j,n) = gvar(iglob(i),jglob(j))
               endif
            enddo
          endif
       enddo
   enddo
   deallocate(iglob,jglob)


 end subroutine ConvertG2T_dbl

 subroutine ConvertG2T_int(gvar,tvar)
   ! --------------------------------------------
   ! Convert between the global to tiled storage
   ! for a variable
   ! --------------------------------------------
   integer(i4), intent(in)  :: gvar(:,:)
   integer(i4), intent(out) :: tvar(:,:,:)

   !-----------------
   ! local variables
   !-----------------
   integer (i4), dimension(:), pointer :: iglob, jglob
   integer (i4) :: n,i,j,jb, je, ib, ie, gbid

   allocate(iglob(nx_block),jglob(ny_block))
   do n=1,nblocks_tot
!JMD       call get_block_parameter(n,jb=jb,je=je,ib=ib,ie=ie,i_glob=iglob,j_glob=jglob)
       gbid = blocks_tropic(n)
       call get_block_parameter(gbid,i_glob=iglob,j_glob=jglob)
       do j=1,ny_block
          if(jglob(j) > 0) then
            do i=1,nx_block
               if(iglob(i) >0) then
                  tvar(i,j,n) = gvar(iglob(i),jglob(j))
               endif
            enddo
          endif
       enddo
   enddo
   deallocate(iglob,jglob)

 end subroutine ConvertG2T_int



 end module domain

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
