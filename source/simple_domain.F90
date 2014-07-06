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
!! Model domain and routines for initializing
!! the domain.  This module also initializes the decompositions and
!! distributions across processors/threads by calling relevent
!! routines in the block, distribution modules.
!<
module simple_domain
    use kinds_mod, only: i4, r8, char_len, log_kind
    use simple_type, only: distrb
    use constants, only: blank_fmt, delim_fmt
    use communicate, only: MPI_COMM_OCN, my_task, master_task, get_num_procs
    use simple_blocks, only: nx_block,ny_block,nblocks_tot,nblocks_x, &
        nblocks_y, nghost, get_block_parameter, set_block_parameter, &
        all_blocks, all_blocks_ij, Neigh, GetEWBlockIndex, GetNSBlockIndex, &
        get_block_id, NumNeigh, ieast, iwest, isouth, inorth, iseast, iswest, &
        ineast, inwest
    use broadcast, only: broadcast_scalar
    use exit_mod, only: sigAbort, exit_pop
    use IOUnitsMod, only: stdout,nmlin,nml_filename
    use domain_size, only: nx_global,ny_global,max_blocks_tropic, &
        block_size_x,block_size_y
    use reductions, only: global_maxval
    use io_serial, only: read_AppMD, read_tiles

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS
    public  :: init_domain_blocks ,&
        init_domain_distribution

    public :: read_solverioT
    public :: read_local_tiles 

    !>
    !! Reads blocks or tiles that are specific for each task. 
    !<
    interface read_local_tiles
        module procedure read_local_tiles_int
        module procedure read_local_tiles_dbl
    end interface

    ! !PUBLIC DATA MEMBERS:
    integer (i4), public :: &
        nblocks_tropic     ! actual number of blocks in tropic distribution

    integer (i4), dimension(:), pointer, public :: &
        blocks_tropic      ! block ids for local blocks in barotropic dist

    !------------------------------------------------------------
    type (distrb), public :: & !  block distribution info
        distrb_tropic      ! block distribution for barotropic part

    !-----------------------------------------------------------------------
    !
    !   module private variables - for the most part these appear as
    !   module variables to facilitate sharing info between init_domain1
    !   and init_domain2.
    !
    !-----------------------------------------------------------------------
    character (char_len) ::      &
        ew_boundary_type,        &! type of domain bndy in each logical
        ns_boundary_type           !    direction (ew is i, ns is j)

    integer (i4) :: &    ! decomposition info
        nprocs_tropic    ! num of processors in barotropic dist

    !***********************************************************************
    contains

    !>
    !!  This routine reads in domain information and calls the routine
    !!  to set up the block decomposition.
    !!
    !! @param reorder An integer array that reorders block from the cannonical 
    !!                 order in the input files to a space-filling curve order 
    !!		  which is used by the miniapp.
    !<
    subroutine init_domain_blocks(reorder)
        integer (i4), pointer :: reorder(:)

        !------------------------------------------------------------------
        !
        !  local variables
        !
        !------------------------------------------------------------------

        integer (i4) :: &
            nml_error          ! namelist read error flag

        integer (i4) :: i,j,n

        integer (i4) :: iblock_src,jblock_src
        integer (i4) :: iblock_east, jblock_east, &
            iblock_west, jblock_west, &
            iblock_south,jblock_south, &
            iblock_north,jblock_north, &
            iblock_swest,jblock_swest, &
            iblock_nwest,jblock_nwest, &
            iblock_seast,jblock_seast, &
            iblock_neast,jblock_neast

        character(len=80) :: tilefile

        !------------------------------------------------------------------
        !
        !  input namelists
        !
        !------------------------------------------------------------------
        namelist /input_nml/ tilefile

        !------------------------------------------------------------------
        !
        !  read domain information from namelist input
        !
        !------------------------------------------------------------------
        nprocs_tropic = get_num_procs()
        ew_boundary_type = 'cyclic'
        ns_boundary_type = 'closed'

        !------------------------------------------------------------------
        !
        !  perform some basic checks on domain
        !
        !------------------------------------------------------------------
        if (nx_global < 1 .or. ny_global < 1) then
            !***
            !*** domain size zero or negative
            !***
            call exit_POP(sigAbort,'Invalid domain: size < 1') ! no domain
        else if (nghost < 2) then
            !***
            !*** must have at least 2 layers of ghost cells
            !***
            call exit_POP(sigAbort,'Not enough ghost cells allocated')
        endif

        !------------------------------------------------------------------
        !
        !  compute block decomposition and details
        !
        !------------------------------------------------------------------

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
        call broadcast_scalar(tilefile,       master_task)

        call read_AppMD(tilefile,all_blocks,reorder)
        nblocks_x   = (nx_global-1)/block_size_x + 1
        nblocks_y   = (ny_global-1)/block_size_y + 1
        nblocks_tot = nblocks_x*nblocks_y

        allocate(all_blocks_ij(nblocks_x,nblocks_y))
        n=0
        do j=1,nblocks_y
        do i=1,nblocks_x
        n = n + 1
        all_blocks_ij(i,j) = n
        enddo
        enddo

        allocate(Neigh(NumNeigh,nblocks_tot))


        do n=1,nblocks_tot
        call get_block_parameter(n,iblock=iblock_src,jblock=jblock_src)

        !*** compute cartesian i,j block indices for each neighbor
        !*** use zero if off the end of closed boundary
        !*** use jnorth=nblocks_y and inorth < 0 for tripole boundary
        !***   to make sure top boundary communicated to all top
        !***   boundary blocks

        call GetEWBlockIndex(ew_boundary_type,iblock_src,jblock_src, &
        iblock_east,jblock_east,iblock_west,jblock_west)

        call GetNSBlockIndex(ns_boundary_type,iblock_src,jblock_src, &
        iblock_north,jblock_north,iblock_south,jblock_south)

        call GetEWBlockIndex(ew_boundary_type,iblock_south,jblock_south, &
        iblock_seast,jblock_seast,iblock_swest,jblock_swest)

        call GetEWBlockIndex(ew_boundary_type,iblock_north,jblock_north, &
        iblock_neast,jblock_neast,iblock_nwest,jblock_nwest)

        ! save Neighbor information
        Neigh(ieast,n)  = get_block_id(iblock_east,jblock_east)
        Neigh(iwest,n)  = get_block_id(iblock_west,jblock_west)
        Neigh(inorth,n) = get_block_id(iblock_north,jblock_north)
        Neigh(isouth,n) = get_block_id(iblock_south,jblock_south)
        Neigh(iseast,n) = get_block_id(iblock_seast,jblock_seast)
        Neigh(ineast,n) = get_block_id(iblock_neast,jblock_neast)
        Neigh(iswest,n) = get_block_id(iblock_swest,jblock_seast)
        Neigh(inwest,n) = get_block_id(iblock_nwest,jblock_nwest)
        enddo

        !------------------------------------------------------------------
        !
        !  Now we need grid info before proceeding further
        !  Print some domain information
        !
        !------------------------------------------------------------------

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
        write(stdout,'(a25,i2)') '  Number of ghost cells: ', nghost
        endif
    end subroutine init_domain_blocks

    !>
    !! This routine calls appropriate setup routines to distribute blocks
    !! across processors and defines arrays with block ids for any local
    !! blocks. Information about ghost cell update routines is also
    !! initialized here through calls to the appropriate boundary routines.
    !<
    subroutine init_domain_distribution(reorder)
        ! !INPUT PARAMETERS:
        integer (i4),pointer,  intent(inout) :: reorder(:)
        integer(i4), save :: timer_distrib

        !------------------------------------------------------------------
        !
        !  local variables
        !
        !------------------------------------------------------------------
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

        type (distrb) :: distrb_calc
        integer(i4), pointer :: Lreorder(:)

        !------------------------------------------------------------------
        !  estimate the amount of work per processor using the topography
        !------------------------------------------------------------------
        allocate(nocn(nblocks_tot))
        nocn = 0
        nocn(:) = all_blocks(:)%npoints

        !------------------------------------------------------------------
        !
        !  determine the distribution of blocks across processors
        !
        !------------------------------------------------------------------
        distrb_tropic = calc_distribution(nprocs_tropic,nocn,reorder)
        deallocate(nocn)

        !------------------------------------------------------------------
        !
        !  allocate and determine block id for any local blocks in each
        !  distribution.
        !
        !------------------------------------------------------------------
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
        endif
    end subroutine init_domain_distribution

    !>
    !!  This subroutine partitions the linear spacing-filling curve into 
    !!  a number of equal length segments.
    !!
    !! @param nprocs    The number of tasks over which to partition the problem.
    !! @param nocn      An integer array which indicates the number of active
    !!                  ocean points within a block.
    !! @param reorder   An integer array which reorders the blocks from
    !!                  cannonical to space-filling curve order.
    !<
    function  calc_distribution(nprocs,nocn,reorder) result(dist)
        integer (i4), intent(in) :: nprocs
        integer (i4), intent(in) :: nocn(:)
        integer (i4), intent(in) :: reorder(:)
        type (distrb) :: dist

        integer (i4) :: nb,nblocks,nActive,nblocksL,extra,s1
        integer (i4) :: n,il,nblocksT,pid
        integer (i4) :: ii,tmp1
        integer (i4), allocatable :: proc_tmp(:)

        nActive = COUNT(reorder(:) > 0) 
        nb = SIZE(nocn)

        dist%nprocs = nprocs
        dist%communicator = MPI_COMM_OCN 
        allocate(dist%proc(nb))
        allocate(dist%local_block(nb))

        nblocks=nActive  ! This is the total number of active blocks
        nblocksL = nblocks/nprocs

        ! every cpu gets nblocksL blocks, but the first 'extra' get nblocksL+1
        extra = mod(nblocks,nprocs)
        s1 = extra*(nblocksL+1)

        dist%proc(:)=0
        dist%local_block(:)=0

        !---------------------------------------
        ! use the cryptic SFC decomposition code 
        !---------------------------------------
        do n=1,nblocks_tot
        ii = reorder(n)
        if(ii>0) then 
        if(ii<=s1) then 
        ii=ii-1
        tmp1 = ii/(nblocksL+1)
        dist%proc(n) = tmp1+1
        else
        ii=ii-s1-1
        tmp1 = ii/nblocksL
        dist%proc(n) = extra + tmp1 + 1
        endif
        endif
        enddo

        allocate(proc_tmp(nprocs))
        proc_tmp = 0
        do n=1,nblocks_tot
        pid = dist%proc(n)
        if(pid>0) then
        proc_tmp(pid) = proc_tmp(pid) + 1
        dist%local_block(n) = proc_tmp(pid)
        endif
        enddo
        deallocate(proc_tmp)
    end function calc_distribution

    !>
    !! This routine determines which blocks in an input distribution are
    !! located on the local processor and creates an array of block ids
    !! for all local blocks.
    !<
    subroutine create_local_block_ids(block_ids, distribution)
        ! !INPUT PARAMETERS:
        type (distrb), intent(in) :: &
        distribution           ! input distribution for which local

        ! !OUTPUT PARAMETERS:
        integer (i4), dimension(:), pointer :: &
        block_ids              ! array of block ids for every block
                               ! that resides on the local processor

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !--------------------------------------------------------------------
        integer (i4) :: &
        n, bid, bcount        ! dummy counters

        !--------------------------------------------------------------------
        !
        !  first determine number of local blocks to allocate array
        !
        !--------------------------------------------------------------------
        bcount = 0
        do n=1,size(distribution%proc)
            if (distribution%proc(n) == my_task+1) bcount = bcount + 1
        end do

        if (bcount > 0) allocate(block_ids(bcount))

        !-------------------------------------------------------------------
        !
        !  now fill array with proper block ids
        !
        !-------------------------------------------------------------------
        if (bcount > 0) then
            do n=1,size(distribution%proc)
                if (distribution%proc(n) == my_task+1) then
                    block_ids(distribution%local_block(n)) = n
                endif
            end do
        endif
    end subroutine create_local_block_ids

    !>
    !! This subroutine reads in the blocks or tiles assigned to each task.
    !! Note that for simplicity this subroutine calls read_tiles which 
    !! reads in all blocks.  Blocks not owned by the task are discarded.
    !!
    !! @param fname   The name of the netCDF from which to read the blocks.
    !! @param vname   The name of the variable to read from the netCDF file.
    !! @param var     The variable into which data will be read.
    !<
    subroutine read_local_tiles_dbl(fname,vname,var)
        character(len=80), intent(in) :: fname
        character(len=*), intent(in) :: vname
        real(r8), dimension(:,:,:) :: var
        real(r8), allocatable :: tmp(:,:,:)
        integer(i4) :: n, gbid

        allocate(tmp(nx_block,ny_block,nblocks_tot))
        call read_tiles(fname,vname,tmp)

        ! just pull out tiles that
        do n=1,nblocks_tropic
            gbid = blocks_tropic(n)
            var(:,:,n) = tmp(:,:,gbid)
        enddo
        deallocate(tmp)
    end subroutine read_local_tiles_dbl

    !>
    !! This subroutine reads in the blocks or tiles assigned to each task.
    !! Note that for simplicity this subroutine calls read_tiles which 
    !! reads in all blocks.  Blocks not owned by the task are discarded.
    !!
    !! @param fname   The name of the netCDF from which to read the blocks.
    !! @param vname   The name of the variable to read from the netCDF file.
    !! @param var     The variable into which data will be read.
    !<
    subroutine read_local_tiles_int(fname,vname,var)
        character(len=80), intent(in) :: fname
        character(len=*), intent(in) :: vname
        integer(i4), dimension(:,:,:) :: var
        real(r8), allocatable :: tmp(:,:,:)
        integer(i4) :: n, gbid

        allocate(tmp(nx_block,ny_block,nblocks_tot))
        call read_tiles(fname,vname,tmp)

        ! just pull out tiles that
        do n=1,nblocks_tropic
            gbid = blocks_tropic(n)
            var(:,:,n) = tmp(:,:,gbid)
        enddo
        deallocate(tmp)
    end subroutine read_local_tiles_int

    !>
    !! This subroutine reads in all solver state necessary to initialize 
    !! the miniapp.
    !!
    !! @param fname    The filename of the netCDF tilefile which contains the
    !!                 solver state.
    !! @param nstep    The timestep in the simulation that the inputfile
    !!                 represents.
    !! @param A0       The diagonal coefficient of the stencil which forms the
    !!                 sparse matrix vector multiply.
    !! @param AN       The northern neighbor coefficient of the stencil.
    !! @param AE       The eastern neighbor coefficient of the stencil.
    !! @param ANE      The northeastern neighbor coefficient of the stencil.
    !! @param RHS      The right and side of the linear system. This
    !!                 corresponds to 'b' from the equation 'Ax=b'.
    !! @param PRESSI   The initial guess for the surface pressure.  This
    !!                 corresponds to 'x'.
    !! @param PRESSF   The final estimate for the surface pressure.
    !! @param RCALCT   The land mask.  Values greater than 0 indicate ocean
    !!                 points.
    !! @param TAREASQ  The area for each grid point.
    !! @param KMT      This integer array which indicates the number of
    !!                 vertical levels.
    !<
    subroutine read_solverioT(fname,nstep,A0,AN,AE,ANE, &
        RHS,PRESSI,PRESSF,RCALCT, TAREASQ,KMT)

        character(len=80), intent(in) :: fname
        integer(i4), intent(inout) :: nstep
        real(r8), intent(inout), &
            dimension(nx_block,ny_block,max_blocks_tropic) :: A0, AN, AE, ANE
        real(r8), intent(inout), &
            dimension(nx_block,ny_block,max_blocks_tropic) :: PRESSI,RHS,PRESSF
        real(r8), intent(inout), &
            dimension(nx_block,ny_block,max_blocks_tropic) :: RCALCT,TAREASQ
        integer(i4), intent(inout), optional, &
            dimension(nx_block,ny_block,max_blocks_tropic) :: KMT

        call read_local_tiles(fname,'A0',A0)
        call read_local_tiles(fname,'AN',AN)
        call read_local_tiles(fname,'AE',AE)
        call read_local_tiles(fname,'ANE',ANE)
        call read_local_tiles(fname,'PRESS_I',PRESSI)
        call read_local_tiles(fname,'RHS',RHS)
        call read_local_tiles(fname,'PRESS_F',PRESSF)
        call read_local_tiles(fname,'RCALCT',RCALCT)
        call read_local_tiles(fname,'TAREASQ',TAREASQ)

        if(present(KMT)) then
            call read_local_tiles(fname,'KMT',KMT)
        endif
    end subroutine read_solverioT
end module simple_domain

!||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
