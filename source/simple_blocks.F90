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
!! Data types and tools for decomposing a global
!! horizontal domain into a set of blocks.  This module contains a data type 
!! for describing each block and contains routines for creating and 
!! querying the block decomposition for a global domain.
!<
module simple_blocks
    use kinds_mod, only: i4, r8
    use simple_type, only: distrb
    use exit_mod, only: sigabort, exit_POP
    use domain_size, only: block_size_x, block_size_y

    implicit none
    private
    save

    ! !PUBLIC TYPES:
    type, public :: AppMD_t   ! block data type
        integer (i4) :: &
            block_id          ,&! global block number
            ib, ie, jb, je    ,&! begin,end indices for physical domain
            iblock, jblock      ! cartesian i,j position for bloc
        integer (i4) :: &
            npoints           ! number of actual ocean points in block
        integer (i4), dimension(:), pointer :: &
            i_glob, j_glob     ! global domain location for each point
    end type

    ! !PUBLIC MEMBER FUNCTIONS:
    public :: get_block_id        ,& 
        get_proc_id         ,&
        set_block_parameter ,&
        get_block_parameter
    public :: GetNSBlockIndex     ,&
        GetEWBlockIndex

    ! !DEFINED PARAMETERS:
    integer (i4), parameter, public :: &
        nghost = 2       ! number of ghost cells around each block
    integer (i4), parameter, public :: &! size of block domain in
        nx_block = block_size_x + 2*nghost,   &!  x,y dir including ghost
        ny_block = block_size_y + 2*nghost     !  cells 

    ! !PUBLIC DATA MEMBERS:
    integer (i4), public :: &
        nblocks_tot      ,&! total number of blocks in decomposition
        nblocks_x        ,&! tot num blocks in i direction
        nblocks_y          ! tot num blocks in j direction

    integer(i4), public, dimension(:),allocatable :: &
    ocn_per_block          ! number of ocean points per block 

    !--------------------
    !  The neighbor graph
    !--------------------
    integer(i4), public, dimension(:,:), allocatable :: &
        Neigh             ! array of block neighbors

    integer(i4), public, parameter :: &
        ieast  = 1, &! index of east neighbor
        iwest  = 2, &! index of west neighbor
        inorth = 3, &! index of north neighbor
        isouth = 4, &! index of south neighbor
        iseast = 5, &! index of southeast neighbor
        iswest = 6, &! index of southwest neighbor
        inwest = 7, &! index of northwest neighbor
        ineast = 8   ! index of northeast neighbor

    integer(i4), public, parameter :: NumNeigh = 8  ! Number of neighbors

    !-----------------------------------------------------------------------
    !
    !  module private data
    !
    !-----------------------------------------------------------------------
        type (AppMD_t), public, pointer, dimension(:) :: &
            all_blocks         ! block information for all blocks in domain

        integer (i4), public, pointer, dimension(:,:) :: &
            all_blocks_ij   ! The linear block index in 2D array

        integer (i4), dimension(:,:), allocatable, target :: &
            i_global,         &! global i index for each point in each block
            j_global           ! global j index for each point in each block

    contains

    !>
    !! This function sets the number of active ocean points in the all_block
    !! variable.
    !<
    subroutine set_block_parameter(bid,npoints)
        ! !INPUT PARAMETERS:

        integer(i4) :: bid                ! global block id
        integer(i4), intent(in), optional :: & ! number of valid ocean points 
        npoints              ! in a block

        if(present(npoints)) all_blocks(bid)%npoints=npoints
        !----------------------------------------------------------------------
    end subroutine set_block_parameter


    !>
    !! This function returns the linear index of a block given the global 
    !! x,y index of a block.
    !<
    function get_block_id(ix,iy)
        integer(i4), intent(in)  :: &
        ix,iy           ! x and y index of block

        ! !OUTPUT PARAMETERS:
        integer(i4) :: get_block_id  ! linear block index

        if(ix>0 .and. iy>0) then 
            get_block_id=all_blocks_ij(ix,iy)
        else
            get_block_id=0
        endif
    end function get_block_id

    !>
    !! This routine returns requested parts of the block data type
    !! for the block associated with the input block id
    !<
    subroutine get_block_parameter(block_id, local_id, ib, ie, jb, je, &
        iblock, jblock, npoints, i_glob, j_glob)

        ! !INPUT PARAMETERS:
            integer (i4), intent(in) :: &
                block_id  ! global block id for which parameters are requested

        ! !OUTPUT PARAMETERS:
            !(optional) parts of block data type to extract if requested
            integer (i4), intent(out), optional :: &
                local_id      ,&! local id assigned to block in current distrb
                ib, ie, jb, je,&! begin,end indices for physical domain
                iblock, jblock,&! cartesian i,j position for bloc
                npoints
            integer (i4), dimension(:), pointer, optional :: &
                i_glob, j_glob     ! global domain location for each point

        !----------------------------------------------------------------------
        !
        !  extract each component of data type if requested
        !
        !----------------------------------------------------------------------
        if (block_id < 1 .or. block_id > nblocks_tot) then
            call exit_POP(sigAbort,'get_block_parameter: invalid block_id')
        endif

        if (present(ib      )) ib       = all_blocks(block_id)%ib
        if (present(ie      )) ie       = all_blocks(block_id)%ie
        if (present(jb      )) jb       = all_blocks(block_id)%jb
        if (present(je      )) je       = all_blocks(block_id)%je
        if (present(npoints )) npoints  = all_blocks(block_id)%npoints 
        if (present(iblock  )) iblock   = all_blocks(block_id)%iblock
        if (present(jblock  )) jblock   = all_blocks(block_id)%jblock
        if (present(i_glob  )) i_glob   = all_blocks(block_id)%i_glob
        if (present(j_glob  )) j_glob   = all_blocks(block_id)%j_glob
    end subroutine get_block_parameter

    !>
    !! This subroutine deallocates the array with block information.
    !<
    subroutine destroy_blocks
        deallocate(all_blocks,all_blocks_ij)
    end subroutine destroy_blocks

    !>
    !!  This subroutine determins the i,j index of the east and west neighbors
    !!  of a source or current block
    !!
    !! @param bndy_type    Controls the type of grid boundary.  Either cyclic
    !!                     or closed. 
    !! @param iblock_src   The 'i' index of the source block.
    !! @param jblock_src   The 'j' index of the source block.
    !! @param iblock_east  The 'i' index of the block to the east of the source
    !!                     block.
    !! @param jblock_east  The 'j' index of the block to the east of the source
    !!                     block.
    !! @param iblock_west  The 'i' index of the block to the west of the source
    !!                     block.
    !! @param jblock_west  The 'j' index of the block to the west of the source
    !!                     block.
    !<
    subroutine GetEWBlockIndex(bndy_type,iblock_src,jblock_src, &
        iblock_east,jblock_east, &
        iblock_west,jblock_west)

        character(*), intent(in) :: bndy_type
        integer(i4), intent(in)  :: iblock_src,jblock_src
        integer(i4), intent(out) :: iblock_east,jblock_east, &
            iblock_west,jblock_west

        select case(bndy_type)
            case ('cyclic')
                iblock_east = mod(iblock_src,nblocks_x) + 1
                iblock_west = iblock_src - 1
                if (iblock_west == 0) iblock_west = nblocks_x
                jblock_east = jblock_src
                jblock_west = jblock_src
            case ('closed')
                iblock_east = iblock_src + 1
                iblock_west = iblock_src - 1
                if (iblock_east > nblocks_x) iblock_east = 0
                if (iblock_west < 1        ) iblock_west = 0
                jblock_east = jblock_src
                jblock_west = jblock_src
            case default
                call exit_POP(sigAbort, 'Unknown east-west boundary type')
        end select
    end subroutine GetEWBlockIndex

    !>
    !!  This subroutine determins the i,j index of the north and 
    !!  south neighbors of a source or current block
    !!
    !! @param bndy_type     Controls the type of grid boundary.  Either cyclic
    !!                      or closed. 
    !! @param iblock_src    The 'i' index of the source block.
    !! @param jblock_src    The 'j' index of the source block.
    !! @param iblock_north  The 'i' index of the block to the north of the
    !!                      source block.
    !! @param jblock_north  The 'j' index of the block to the north of the 
    !!                      source block.
    !! @param iblock_south  The 'i' index of the block to the south of the 
    !!                      source block.
    !! @param jblock_south  The 'j' index of the block to the south of the 
    !!                      source block.
    !<
    subroutine GetNSBlockIndex(bndy_type, &
        iblock_src,jblock_src, &
        iblock_north,jblock_north, &
        iblock_south,jblock_south)

        character(*), intent(in)       :: bndy_type
        integer(i4), intent(in)  :: iblock_src,jblock_src
        integer(i4), intent(out) :: iblock_north,jblock_north, &
        iblock_south,jblock_south

        select case(bndy_type)
            case ('cyclic')
                jblock_north = mod(jblock_src,nblocks_y) + 1
                jblock_south = jblock_src - 1
                if (jblock_south == 0) jblock_south = nblocks_y
                iblock_north = iblock_src
                iblock_south = iblock_src
            case ('closed')
                jblock_north = jblock_src + 1
                jblock_south = jblock_src - 1
                if (jblock_north > nblocks_y) jblock_north = 0
                if (jblock_south < 1        ) jblock_south = 0
                iblock_north = iblock_src
                iblock_south = iblock_src
            case default
                call exit_POP(sigAbort, 'Unknown north-south boundary type')
        end select
    end subroutine GetNSBlockIndex

    !>
    !! This determins the task or processor id for a neighboring block
    !!
    !! @param dist      A data structure which stores the distribution of
    !!                  blocks.  
    !! @param ineigh    An integer indicating the cardinal and ordinal
    !!                  directions.
    !! @param block_id  The current or source block identifier.
    !<
    function get_proc_id(dist,ineigh,block_id)
        type (distrb) :: dist
        integer (i4) :: ineigh
        integer (i4) :: block_id
        integer (i4) :: get_proc_id
        integer (i4) :: nbid

        nbid=Neigh(ineigh,block_id)
        if(nbid >0) then 
            if(all_blocks(nbid)%npoints>0) then
                get_proc_id = dist%proc(nbid)-1
            else
                get_proc_id = 0
            endif
        else 
            get_proc_id = 0
        endif
    end function get_proc_id

end module simple_blocks
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
