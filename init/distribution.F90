!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 module distribution

!BOP
! !MODULE: distribution
!
! !DESCRIPTION:
!  This module provides data types and routines for distributing
!  blocks across processors.
!
! !REVISION HISTORY:
!  CVS:$Id: distribution.F90,v 1.11 2007/01/29 17:04:15 dennis Exp $
!  CVS:$Name:  $

! !USES:

   use kinds_mod, only: i4, r4
   use simple_type, only : distrb
   use communicate, only: my_task,create_communicator
   use blocks, only: nblocks_tot, nblocks_x, nblocks_y
   use spacecurve_mod, only: factor_t, isfactorable, genspacecurve,  &
	PrintCurve, factor, ProdFactor, MatchFactor, PrintFactor
   use exit_mod, only: sigabort, exit_POP

   implicit none
   private
   save

! !PUBLIC TYPES:

!JMD   type, public :: distrb  ! distribution data type
!JMD      integer (i4) :: &
!JMD        nprocs            ,&! number of processors in this dist
!JMD         communicator        ! communicator to use in this dist
!JMD
!JMD      integer (i4), dimension(:), pointer :: &
!JMD         proc              ,&! processor location for this block
!JMD         local_block         ! block position in local array on proc
!JMD   end type

   logical, public :: sameDistribution

! !PUBLIC MEMBER FUNCTIONS:

   public :: create_distribution, &
   	     destroy_distribution, &
             create_local_block_ids

!EOP
!BOC
!EOC
!***********************************************************************

 contains

   subroutine destroy_distribution(dist)

      type (distrb) :: dist

      deallocate(dist%proc)
      deallocate(dist%local_block)

   end subroutine destroy_distribution

!***********************************************************************
!BOP
! !IROUTINE: create_distribution
! !INTERFACE:

 function create_distribution(dist_type, nprocs, work_per_block, reorder)

! !DESCRIPTION:
!  This routine determines the distribution of blocks across processors
!  by call the appropriate subroutine based on distribution type
!  requested.  Currently only two distributions are supported:
!  2-d Cartesian distribution (cartesian) and a load-balanced
!  distribution (balanced) based on an input amount of work per
!  block.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      dist_type             ! method for distributing blocks
                            !  either cartesian or balanced

   integer (i4), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (i4), dimension(:), intent(inout) :: &
      work_per_block        ! amount of work per block

   integer (i4), dimension(:), intent(inout) :: reorder

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distribution   ! resulting structure describing
                            !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  select the appropriate distribution type
!
!----------------------------------------------------------------------

   select case (trim(dist_type))

   case('cartesian')

      create_distribution = create_distrb_cart(nprocs, work_per_block)

   case('spacecurve')

      create_distribution = create_distrb_spacecurve(nprocs, &
						   work_per_block, reorder)
   case default

      call exit_POP(sigAbort,'distribution: unknown distribution type')

   end select

!-----------------------------------------------------------------------
!EOC
!  print *,'create_distribution: dist%proc ',create_distribution%proc

 end function create_distribution

!***********************************************************************
!BOP
! !IROUTINE: create_local_block_ids
! !INTERFACE:

 subroutine create_local_block_ids(block_ids, distribution)

! !DESCRIPTION:
!  This routine determines which blocks in an input distribution are
!  located on the local processor and creates an array of block ids
!  for all local blocks.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   type (distrb), intent(in) :: &
      distribution           ! input distribution for which local
                             !  blocks required

! !OUTPUT PARAMETERS:

   integer (i4), dimension(:), pointer :: &
      block_ids              ! array of block ids for every block
                             ! that resides on the local processor

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (i4) :: &
      n, bid, bcount        ! dummy counters

!-----------------------------------------------------------------------
!
!  first determine number of local blocks to allocate array
!
!-----------------------------------------------------------------------

   bcount = 0
   do n=1,size(distribution%proc)
      if (distribution%proc(n) == my_task+1) bcount = bcount + 1
   end do

   if (bcount > 0) allocate(block_ids(bcount))

!-----------------------------------------------------------------------
!
!  now fill array with proper block ids
!
!-----------------------------------------------------------------------

   if (bcount > 0) then
      do n=1,size(distribution%proc)
         if (distribution%proc(n) == my_task+1) then
            block_ids(distribution%local_block(n)) = n
         endif
      end do
   endif

!EOC

 end subroutine create_local_block_ids

!***********************************************************************
!BOP
! !IROUTINE: create_distrb_cart
! !INTERFACE:

 function create_distrb_cart(nprocs, work_per_block)

! !DESCRIPTION:
!  This function creates a distribution of blocks across processors
!  using a 2-d Cartesian distribution.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
      nprocs     ! number of processors in this distribution

   integer (i4), dimension(:), intent(inout) :: &
      work_per_block        ! amount of work per block

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_cart  ! resulting structure describing Cartesian
                          !  distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (i4) :: &
      i, j, n               ,&! dummy loop indices
      iblock, jblock, nblck ,&!
      is, ie, js, je        ,&! start, end block indices for each proc
      local_block           ,&! block location on this processor
      nprocs_x              ,&! num of procs in x for global domain
      nprocs_y              ,&! num of procs in y for global domain
      nblocks_x_loc         ,&! num of blocks per processor in x
      nblocks_y_loc           ! num of blocks per processor in y

   type (distrb) :: dist  ! temp hold distribution

!----------------------------------------------------------------------
!
!  create communicator for this distribution
!
!----------------------------------------------------------------------

   call create_communicator(dist%communicator, nprocs)

!----------------------------------------------------------------------
!
!  try to find best processor arrangement
!
!----------------------------------------------------------------------

   dist%nprocs = nprocs

   call proc_decomposition(dist%nprocs, nprocs_x, nprocs_y)

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc       (nblocks_tot), &
             dist%local_block(nblocks_tot))

!----------------------------------------------------------------------
!
!  distribute blocks linearly across processors in each direction
!
!----------------------------------------------------------------------

   nblocks_x_loc = (nblocks_x-1)/nprocs_x + 1
   nblocks_y_loc = (nblocks_y-1)/nprocs_y + 1

   do j=1,nprocs_y
   do i=1,nprocs_x
      n = (j-1)*nprocs_x + i

      is = (i-1)*nblocks_x_loc + 1
      ie =  i   *nblocks_x_loc
      if (ie > nblocks_x) ie = nblocks_x
      js = (j-1)*nblocks_y_loc + 1
      je =  j   *nblocks_y_loc
      if (je > nblocks_y) je = nblocks_y

      local_block = 0
      do jblock = js,je
      do iblock = is,ie
         nblck = (jblock - 1)*nblocks_x + iblock
         if (work_per_block(nblck) /= 0) then
            local_block = local_block + 1
            dist%proc(nblck) = n
            dist%local_block(nblck) = local_block
         else
            dist%proc(nblck) = 0
            dist%local_block(nblck) = 0
         endif
      end do
      end do
   end do
   end do

!----------------------------------------------------------------------

   create_distrb_cart = dist  ! return the result

!----------------------------------------------------------------------
!EOC

 end function create_distrb_cart
!**********************************************************************
!BOP
! !IROUTINE: create_distrb_spacecurve
! !INTERFACE:

 function create_distrb_spacecurve(nprocs,work_per_block,reorder)

! !Description:
!  This function distributes blocks across processors in a 
!  load-balanced manner using space-filling curves
!
! !REVISION HISTORY:
!  added by J. Dennis 3/10/06 

! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
      nprocs                ! number of processors in this distribution

   integer (i4), dimension(:), intent(inout) :: &
      work_per_block        ! amount of work per block

   integer (i4), dimension(:), intent(inout) :: reorder

! !OUTPUT PARAMETERS:

   type (distrb) :: &
      create_distrb_spacecurve  ! resulting structure describing
                                ! load-balanced distribution of blocks

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (i4) :: &
      i,j,k,n              ,&! dummy loop indices
      pid                  ,&! dummy for processor id
      local_block          ,&! local block position on processor
      max_work             ,&! max amount of work in any block
      nprocs_x             ,&! num of procs in x for global domain
      nprocs_y               ! num of procs in y for global domain

   integer (i4), dimension(:),allocatable :: &
	idxT_i,idxT_j

   integer (i4), dimension(:,:),allocatable :: Mesh,Mesh2,Mesh3
   integer (i4) :: nblocksL,nblocks,ii,extra,i2,j2,tmp1,s1,ig

   integer (i4) :: ierr
   logical, parameter :: Debug = .TRUE.

   integer (i4), dimension(:), allocatable :: &
      priority           ,&! priority for moving blocks
      work_tmp           ,&! work per row or column for rake algrthm
      proc_tmp           ,&! temp processor id for rake algrthm
      block_count          ! counter to determine local block indx

   type (distrb) :: dist  ! temp hold distribution

   type (factor_t) :: xdim, ydim
   integer (i4) :: it,jj,itmp
   integer (i4) :: curveSize, sb_x, sb_y,numfac
   integer (i4) :: subNum, sfcNum
   logical :: foundX
   logical, parameter :: verbose = .false.
   
   integer (i4) :: icnt

   integer (i4), allocatable :: Lindx(:), Lindx2(:)
!----------------------------------------------------------------------
!
!  first set up as Cartesian distribution
!  retain the Cartesian distribution if nblocks_tot = nprocs
!  to avoid processors with no work
!
!----------------------------------------------------------------------
   !------------------------------------------------------
   ! Space filling curves only work if:
   ! 
   !    nblocks_x = nblocks_y	
   ! 	nblocks_x = 2^m 3^n 5^o where m,n,o are integers
   !------------------------------------------------------
   if((.not. IsFactorable(nblocks_y)) .or. (.not. IsFactorable(nblocks_x))) then 
     stop 'can not generate SFC partitioning'
     create_distrb_spacecurve = create_distrb_cart(nprocs, work_per_block)
     return
   endif


   !-----------------------------------------------
   ! Factor the numbers of blocks in each dimension
   !-----------------------------------------------
   xdim = Factor(nblocks_x)
!   if(verbose) print *,'IAM: ',my_task,'create_distrb_spacecurve: point #1.4'
   ydim = Factor(nblocks_y)
!   if(verbose) print *,'IAM: ',my_task,'create_distrb_spacecurve: point #1.5'
   numfac = xdim%numfact
    print *,'nblocks_x: ',nblocks_x
    print *,'nblocks_y: ',nblocks_y
!    call PrintFactor('xdim',xdim)
!    call PrintFactor('ydim',ydim)
!   if(verbose) print *,'IAM: ',my_task,'create_distrb_spacecurve: point #2'

   !---------------------------------------------
   ! Match the common factors to create SFC curve
   !---------------------------------------------
   curveSize=1
   do it=1,numfac
      call MatchFactor(xdim,ydim,itmp,foundX)
      curveSize = itmp*curveSize
   enddo
!   if(verbose) print *,'IAM: ',my_task,'create_distrb_spacecurve: point #3'
   !--------------------------------------
   ! determine the size of the sub-blocks
   ! within the space-filling curve
   !--------------------------------------
   sb_x = ProdFactor(xdim)
   sb_y = ProdFactor(ydim)

   print *,'sb_x: ',sb_x
   print *,'sb_y: ',sb_y
   print *,'curveSize: ',curveSize
   call create_communicator(dist%communicator, nprocs)

   dist%nprocs = nprocs

!----------------------------------------------------------------------
!
!  allocate space for decomposition
!
!----------------------------------------------------------------------

   allocate (dist%proc       (nblocks_tot), &
             dist%local_block(nblocks_tot))
   dist%proc=0
   dist%local_block=0


!----------------------------------------------------------------------
!  Create the array to hold the SFC
!----------------------------------------------------------------------
   allocate(Mesh(curveSize,curveSize))
   allocate(Mesh2(nblocks_x,nblocks_y))
   allocate(Mesh3(nblocks_x,nblocks_y))
   Mesh  = 0
   Mesh2 = 0
   Mesh3 = 0

   allocate(idxT_i(nblocks_tot),idxT_j(nblocks_tot))
   allocate(Lindx(nblocks_tot),Lindx2(nblocks_tot))


!----------------------------------------------------------------------
!  Generate the space-filling curve
!----------------------------------------------------------------------
   call GenSpaceCurve(Mesh)
   Mesh = Mesh + 1  ! make it 1-based indexing
!   if(Debug) then 
!     if(my_task ==0) call PrintCurve(Mesh)
!   endif 
   !-----------------------------------------------
   ! Reindex the SFC to address internal sub-blocks
   !-----------------------------------------------
   do j=1,curveSize
   do i=1,curveSize
      icnt = 0
      sfcNum = (Mesh(i,j) - 1)*(sb_x*sb_y) + 1
      do ii=1,sb_x
      do jj=1,sb_y
         subNum = (jj-1) + (ii-1)*sb_y
         i2 = (i-1)*sb_x + ii
         j2 = (j-1)*sb_y + jj
         Mesh2(i2,j2) = sfcNum + subNum
	 n = (j2-1)*nblocks_x + i2
         if(work_per_block(n) /=0 ) then  
	    icnt = icnt + 1
	 endif
      enddo
      enddo
      !--------------------------------------------------------------------------
      ! if there is any work in the set of sub-blocks then reset the entire group
      !--------------------------------------------------------------------------
      if(icnt > 0) then 
	 do ii=1,sb_x	
         do jj=1,sb_y
            i2 = (i-1)*sb_x + ii
            j2 = (j-1)*sb_y + jj
	    n = (j2-1)*nblocks_x + i2
	    work_per_block(n) = 1.0
	 enddo
	 enddo
      endif

   enddo
   enddo
!   call PrintCurve(mesh2)
!   stop 'after after creation of Mesh2 '
   !------------------------------------------------
   ! create a linear array of i,j coordinates of SFC
   !------------------------------------------------
   reorder=0;
   idxT_i=0;idxT_j=0;Lindx2=0;Lindx=0
   do j=1,nblocks_y
     do i=1,nblocks_x
	n = (j-1)*nblocks_x + i
	ig = Mesh2(i,j)
	if(work_per_block(n) /= 0) then 
	    idxT_i(ig)=i;idxT_j(ig)=j
	endif
        Lindx(n) = ig
     enddo
   enddo
   !-----------------------------
   ! Compress out the land blocks 
   !-----------------------------
   ii=0
   do i=1,nblocks_tot
      if(IdxT_i(i) .gt. 0) then 
	 ii=ii+1
	 Mesh3(idxT_i(i),idxT_j(i)) = ii
	 n = (idxT_j(i)-1)*nblocks_x + idxT_i(i)
	 Lindx2(n) = ii
      endif
   enddo
   nblocks=ii  
!   if(Debug) then 
!     if(my_task==0) call PrintCurve(Mesh2)
!   endif
   
   nblocksL = nblocks/nprocs
   ! every cpu gets nblocksL blocks, but the first 'extra' get nblocksL+1
   extra = mod(nblocks,nprocs)
   s1 = extra*(nblocksL+1)
   ! split curve into two curves:
   ! 1 ... s1  s2 ... nblocks
   !
   !  s1 = extra*(nblocksL+1)         (count be 0)
   !  s2 = s1+1
   !
   ! First region gets nblocksL+1 blocks per partition
   ! Second region gets nblocksL blocks per partition
   if(Debug) print *,'nprocs,extra,nblocks,nblocksL,s1: ', &
		nprocs,extra,nblocks,nblocksL,s1

   do j=1,nblocks_y
   do i=1,nblocks_x
      n = (j-1)*nblocks_x + i
!      i2 = idxT_i(n)
!      j2 = idxT_j(n)
      ii = Mesh3(i,j)
      if(ii>0) then 
        !DBG if(my_task ==0) print *,'i,j,ii:= ',i,j,ii
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
   enddo

!----------------------------------------------------------------------
!  Reset the dist data structure 
!----------------------------------------------------------------------

   allocate(proc_tmp(nprocs))
   proc_tmp = 0

   do n=1,nblocks_tot
      pid = dist%proc(n)
      if(pid>0) then 
        proc_tmp(pid) = proc_tmp(pid) + 1 	
        dist%local_block(n) = proc_tmp(pid)
      endif
   enddo

!   if(Debug) then 
!      if(my_task==0) print *,'dist%proc:= ',dist%proc
!      print *,'IAM: ',my_task,' SpaceCurve: Number of blocks {total,local} :=', &
!		nblocks_tot,nblocks,proc_tmp(my_task+1)
!   endif
   do n=1,nblocks_tot
	reorder(Lindx(n))=n
   enddo
   reorder(:) = Lindx2(:)

   deallocate(proc_tmp)
   deallocate(Lindx,Lindx2)
   ierr=1

   deallocate(Mesh,Mesh2,Mesh3)
   deallocate(idxT_i,idxT_j)
!----------------------------------------------------------------------
   create_distrb_spacecurve = dist  ! return the result

!----------------------------------------------------------------------
   print *,'Active ocean blocks: ',count(reorder>0)
!EOC

 end function create_distrb_spacecurve

!**********************************************************************
!BOP
! !IROUTINE: proc_decomposition
! !INTERFACE:

 subroutine proc_decomposition(nprocs, nprocs_x, nprocs_y)

! !DESCRIPTION:
!  This subroutine attempts to find an optimal (nearly square)
!  2d processor decomposition for a given number of processors.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (i4), intent(in) :: &
      nprocs                       ! total number or processors

! !OUTPUT PARAMETERS:

   integer (i4), intent(out) :: &
      nprocs_x, nprocs_y           ! number of procs in each dimension

!EOP
!BOC
!----------------------------------------------------------------------
!
!  local variables
!
!----------------------------------------------------------------------

   integer (i4) :: &
      iguess, jguess               ! guesses for nproc_x,y

   real (r4) :: &
      square                       ! square root of nprocs

!----------------------------------------------------------------------
!
!  start with an initial guess that is closest to square decomp
!
!----------------------------------------------------------------------

   square = sqrt(real(nprocs))
   nprocs_x = 0
   nprocs_y = 0

   iguess = nint(square)

!----------------------------------------------------------------------
!
!  try various decompositions to find the best
!
!----------------------------------------------------------------------

   proc_loop: do
      jguess = nprocs/iguess

      if (iguess*jguess == nprocs) then ! valid decomp

         !***
         !*** if the blocks can be evenly distributed, it is a
         !*** good decomposition
         !***

         if (mod(nblocks_x,iguess) == 0 .and. &
             mod(nblocks_y,jguess) == 0) then
            nprocs_x = iguess
            nprocs_y = jguess
            exit proc_loop

         !***
         !*** if the blocks can be evenly distributed in a
         !*** transposed direction, it is a good decomposition
         !***

         else if (mod(nblocks_x,jguess) == 0 .and. &
                mod(nblocks_y,iguess) == 0) then
            nprocs_x = jguess
            nprocs_y = iguess
            exit proc_loop

         !***
         !*** A valid decomposition, but keep searching for
         !***  a better one
         !***

         else
            if (nprocs_x == 0) then
               nprocs_x = iguess
               nprocs_y = jguess
            endif
            iguess = iguess - 1
            if (iguess == 0) then
               exit proc_loop
            else
               cycle proc_loop
            endif
         endif

      else ! invalid decomp - keep trying

         iguess = iguess - 1
         if (iguess == 0) then
            exit proc_loop
         else
            cycle proc_loop
         endif
      endif
   end do proc_loop

   if (nprocs_x == 0) then
      call exit_POP(sigAbort,'Unable to find 2d processor config')
   endif

!----------------------------------------------------------------------
!EOC

 end subroutine proc_decomposition

!***********************************************************************

end module distribution

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
