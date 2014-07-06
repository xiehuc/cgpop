!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module blocks

!BOP
! !MODULE: blocks
!
! !DESCRIPTION: 
!  This module contains data types and tools for decomposing a global
!  horizontal domain into a set of blocks.  It contains a data type 
!  for describing each block and contains routines for creating and 
!  querying the block decomposition for a global domain.
!
! !REVISION HISTORY:
!  CVS:$Id: blocks.F90,v 1.12 2006/08/23 02:24:06 dennis Exp $
!  CVS:$Name:  $
!
! !USES:

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
!JMD         local_id          ,&! local address of block in current distrib
         ib, ie, jb, je    ,&! begin,end indices for physical domain
         iblock, jblock      ! cartesian i,j position for bloc
         
      integer (i4) :: &
	  npoints	     ! number of actual ocean points in block
      integer (i4), dimension(:), pointer :: &
         i_glob, j_glob     ! global domain location for each point
   end type

! !PUBLIC MEMBER FUNCTIONS:

   public :: create_blocks       ,&
             destroy_blocks      ,&
	     get_block_id        ,& 
             get_proc_id         ,&
             set_block_parameter ,&
             get_block_parameter

   public :: GetNSBlockIndex     ,&
	     GetEWBlockIndex

#if 0
   public :: ConvertG2T

   interface ConvertG2T
      module procedure ConvertG2T_int
      module procedure ConvertG2T_dbl
   end interface
#endif
	        
  public ::  get_block
! !DEFINED PARAMETERS:
  

   integer (i4), parameter, public :: &
      nghost = 2       ! number of ghost cells around each block

   integer (i4), public ::  &! size of block domain in
      nx_block,ny_block

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
        ineast = 8

    integer(i4), public, parameter :: NumNeigh = 8  ! Number of neighbors

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module private data
!
!-----------------------------------------------------------------------

   type (AppMD_t), public, pointer, dimension(:) :: &
      all_blocks         ! block information for all blocks in domain

   integer (i4), public, pointer, dimension(:,:) :: &
	all_blocks_ij	! The linear block index in 2D array

   integer (i4), dimension(:,:), allocatable, target :: &
      i_global,         &! global i index for each point in each block
      j_global           ! global j index for each point in each block

!EOC
!***********************************************************************

contains

!***********************************************************************
!BOP
! !IROUTINE: SetBlock
! !INTERFACE:

  subroutine set_block_parameter(bid,npoints)

! !DESCRIPTION:
!  This function sets the number of active ocean points in the all_block
!   variable.
!
! !REVISION HISTORY:
!   J. Dennis (dennis@ucar.edu)  April, 2006
!
! !INPUT PARAMETERS:

    integer(i4) :: bid                ! global block id
    integer(i4), intent(in), optional :: & ! number of valid ocean points 
		npoints				 ! in a block
!EOP
!BOC

    if(present(npoints)) all_blocks(bid)%npoints=npoints

!EOC
!----------------------------------------------------------------------
  end subroutine set_block_parameter

!***********************************************************************
!BOP
! !IROUTINE: create_blocks
! !INTERFACE:

 subroutine create_blocks(nx_global, ny_global, ew_boundary_type, &
                                                ns_boundary_type)

! !DESCRIPTION:
!  This subroutine decomposes the global domain into blocks and
!  fills the data structures with all the necessary block information.
!
! !REVISION HISTORY: 
!  same as module
!
! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
      nx_global, ny_global           ! global domain size in x,y

   character (*), intent(in) :: &
      ew_boundary_type,  &! type of boundary in logical east-west dir
      ns_boundary_type    ! type of boundary in logical north-south dir

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (i4) :: &
      i, ip1, j, jp1, n    ,&! loop indices
      iblock, jblock       ,&! block loop indices
      is, ie, js, je         ! temp start, end indices

!----------------------------------------------------------------------
!
!  compute number of blocks and cartesian decomposition
!  if the requested block size does not divide the global domain
!  size evenly, add additional block space to accomodate padding
!
!----------------------------------------------------------------------

   nx_block = block_size_x + 2*nghost   !  x,y dir including ghost
   ny_block = block_size_y + 2*nghost   !  cells 
   nblocks_x   = (nx_global-1)/block_size_x + 1
   nblocks_y   = (ny_global-1)/block_size_y + 1
   nblocks_tot = nblocks_x*nblocks_y

!----------------------------------------------------------------------
!
!  allocate block arrays
!
!----------------------------------------------------------------------

   allocate(all_blocks(nblocks_tot))
   allocate(i_global(nx_block,nblocks_tot), &
            j_global(ny_block,nblocks_tot))

   allocate(ocn_per_block(nblocks_tot))
   ocn_per_block=0

   allocate(all_blocks_ij(nblocks_x,nblocks_y))

!----------------------------------------------------------------------
!
!  fill block data structures for all blocks in domain
!
!----------------------------------------------------------------------

   n = 0
   do jblock=1,nblocks_y
      js = (jblock-1)*block_size_y + 1
      je = js + block_size_y - 1
      if (js > ny_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: ny_block too large?')
      if (je > ny_global) je = ny_global ! pad array

      do iblock=1,nblocks_x
         n = n + 1  ! global block id

         is = (iblock-1)*block_size_x + 1
         ie = is + block_size_x - 1
         if (is > nx_global) call exit_POP(sigAbort, &
               'create_blocks: Bad block decomp: nx_block too large?')
         if (ie > nx_global) ie = nx_global

         all_blocks(n)%block_id = n
         all_blocks(n)%iblock   = iblock
         all_blocks(n)%jblock   = jblock
         all_blocks(n)%ib       = nghost + 1
         all_blocks(n)%jb       = nghost + 1
         all_blocks(n)%ie       = nx_block - nghost ! default value
         all_blocks(n)%je       = ny_block - nghost ! default value

         !JMD----------------------------------------------
         !JMD store the linear index in a rectangular array
         !JMD----------------------------------------------
         all_blocks_ij(iblock,jblock)=n

         do j=1,ny_block
            j_global(j,n) = js - nghost + j - 1


            !*** southern ghost cells

            if (j_global(j,n) < 1) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) + ny_global
               case ('closed')
                  j_global(j,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select
            endif

            !*** padding required

            if (j_global(j,n) > ny_global + nghost) then
               j_global(j,n) = 0   ! padding

            !*** northern ghost cells

            else if (j_global(j,n) > ny_global) then
               select case (ns_boundary_type)
               case ('cyclic')
                  j_global(j,n) = j_global(j,n) - ny_global
               case ('closed')
                  j_global(j,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown n-s bndy type')
               end select

            !*** set last physical point if padded domain

            else if (j_global(j,n) == ny_global .and. &
                     j > all_blocks(n)%jb) then
               all_blocks(n)%je = j   ! last physical point in padded domain
            endif
         end do

         all_blocks(n)%j_glob => j_global(:,n)

         do i=1,nx_block
            i_global(i,n) = is - nghost + i - 1

            !*** western ghost cells

            if (i_global(i,n) < 1) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) + nx_global
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select
            endif

            !*** padded domain - fill padded region with zero

            if (i_global(i,n) > nx_global + nghost) then
               i_global(i,n) = 0

            !*** eastern ghost cells

            else if (i_global(i,n) > nx_global) then
               select case (ew_boundary_type)
               case ('cyclic')
                  i_global(i,n) = i_global(i,n) - nx_global
               case ('closed')
                  i_global(i,n) = 0
               case default
                  call exit_POP(sigAbort, &
                                'create_blocks: unknown e-w bndy type')
               end select

            !*** last physical point in padded domain

            else if (i_global(i,n) == nx_global .and. &
                     i > all_blocks(n)%ib) then
               all_blocks(n)%ie = i
            endif
         end do

         all_blocks(n)%i_glob => i_global(:,n)

      end do
   end do

!EOC
!----------------------------------------------------------------------

end subroutine create_blocks

!***********************************************************************
!BOP
! !IROUTINE: get_block
! !INTERFACE:

 function get_block(block_id,local_id)

! !DESCRIPTION:
!  This function returns the block data structure for the block
!  associated with the input block id.
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
      block_id,   &! global block id for requested block info
      local_id     ! local  block id to assign to this block

! !OUTPUT PARAMETERS:

   type (AppMD_t) :: &
      get_block    ! block information returned for requested block

!EOP
!BOC
!----------------------------------------------------------------------
!
!  check for valid id.  if valid, return block info for requested block
!
!----------------------------------------------------------------------

!DBG   print *,'get_block: nblocks_tot:',nblocks_tot
   if (block_id < 1 .or. block_id > nblocks_tot) then
      call exit_POP(sigAbort,'get_block: invalid block_id')
   endif

   get_block = all_blocks(block_id)
!JMD   get_block%local_id = local_id

!----------------------------------------------------------------------
!EOC

 end function get_block

!***********************************************************************
!BOP
! !IROUTINE: get_block_id
! !INTERFACE:

 function get_block_id(ix,iy)

! !DESCRIPTION:
!  This function returns the linear index of a block given the global 
!  x,y index of a block.

   integer(i4), intent(in)  :: &
	ix,iy			! x and y index of block

! !OUTPUT PARAMETERS:

   integer(i4) :: get_block_id  ! linear block index

!EOP
!BOC
!----------------------------------------------------------------------
    
    if(ix>0 .and. iy>0) then 
      get_block_id=all_blocks_ij(ix,iy)
    else
      get_block_id=0
    endif

!----------------------------------------------------------------------
!EOC

 end function get_block_id

!**********************************************************************
!BOP
! !IROUTINE: get_block_parameter
! !INTERFACE:

 subroutine get_block_parameter(block_id, local_id, ib, ie, jb, je, &
                                iblock, jblock, npoints, i_glob, j_glob)

! !DESCRIPTION:
!  This routine returns requested parts of the block data type
!  for the block associated with the input block id
!
! !REVISION HISTORY:
!  same as module
!
! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
      block_id   ! global block id for which parameters are requested

! !OUTPUT PARAMETERS:

   !(optional) parts of block data type to extract if requested

   integer (i4), intent(out), optional :: &
      local_id         ,&! local id assigned to block in current distrb
      ib, ie, jb, je   ,&! begin,end indices for physical domain
      iblock, jblock   ,&! cartesian i,j position for bloc
      npoints

   integer (i4), dimension(:), pointer, optional :: &
      i_glob, j_glob     ! global domain location for each point

!EOP
!BOC
!----------------------------------------------------------------------
!
!  extract each component of data type if requested
!
!----------------------------------------------------------------------

   if (block_id < 1 .or. block_id > nblocks_tot) then
      call exit_POP(sigAbort,'get_block_parameter: invalid block_id')
   endif

!JMD   if (present(local_id)) local_id = all_blocks(block_id)%local_id
   if (present(ib      )) ib       = all_blocks(block_id)%ib
   if (present(ie      )) ie       = all_blocks(block_id)%ie
   if (present(jb      )) jb       = all_blocks(block_id)%jb
   if (present(je      )) je       = all_blocks(block_id)%je
   if (present(npoints )) npoints  = all_blocks(block_id)%npoints 
   if (present(iblock  )) iblock   = all_blocks(block_id)%iblock
   if (present(jblock  )) jblock   = all_blocks(block_id)%jblock
   if (present(i_glob  )) i_glob   = all_blocks(block_id)%i_glob
   if (present(j_glob  )) j_glob   = all_blocks(block_id)%j_glob

!----------------------------------------------------------------------
!EOC

 end subroutine get_block_parameter

!**********************************************************************
!BOP
! !IROUTINE: destroy_blocks
! !INTERFACE:

 subroutine destroy_blocks

! !DESCRIPTION:
!  This subroutine deallocates the array with block information.
!
! !REVISION HISTORY:
!  same as module
!EOP
!----------------------------------------------------------------------
!BOC

   deallocate(all_blocks,all_blocks_ij)
   deallocate(i_global,j_global)
   deallocate(ocn_per_block)

!EOC
!----------------------------------------------------------------------

 end subroutine destroy_blocks

 subroutine GetEWBlockIndex(bndy_type,iblock_src,jblock_src, &
                                iblock_east,jblock_east, &
                                iblock_west,jblock_west)

    character(*), intent(in)       :: bndy_type
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

!***********************************************************************

 function get_proc_id(dist,ineigh,block_id)

    type (distrb) :: dist
    integer (i4) :: ineigh
    integer (i4) :: block_id
    integer (i4) :: get_proc_id


    if(all_blocks(Neigh(ineigh,block_id))%npoints>0) then
      get_proc_id = dist%proc(Neigh(ineigh,block_id))-1
    else
      get_proc_id = 0
    endif

 end function get_proc_id

#if 0
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
#endif

 end module blocks

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
