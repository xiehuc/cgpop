!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
module linear

    ! !DESCRIPTION:
    !  This module contains routines for converting to and from 
    !  the linear data structure in the solver. 

    ! !USES:
    use kinds_mod, only : i4, r4, r8, log_kind
    use simple_type, only: distrb
    use constants, only : field_type_scalar, field_loc_center, &
        field_loc_Nface, field_loc_NEcorner, c0
    use blocks, only : set_block_parameter, nx_block,ny_block, nblocks_y, &
        get_block_parameter, nblocks_tot
    use boundary, only : boundary_2d
    use domain, only : distrb_tropic, bndy_tropic, nblocks_tropic, &
        blocks_tropic, write_ordered_tiles, read_ordered_tiles
    use domain_size, only : block_size_x, block_size_y, max_blocks_tropic, &
        nx_global
    use io_serial, only : write_tiles, read_tiles

    implicit none
    private
    save

    ! !PUBLIC MEMBER FUNCTIONS:
    public :: Convert2DtoLinear
    public :: ConvertLinearto2D
    public :: initDOF
    public :: update_halo

    !-----------------------------------------------------------------------
    !
    !  overload module functions
    !
    !-----------------------------------------------------------------------
    interface Convert2DtoLinear
        module procedure Convert2DtoLinear_dbl
    end interface

    interface ConvertLinearto2D
        module procedure ConvertLinearto2D_dbl
    end interface

    interface update_halo
        module procedure update_halo_dbl
        module procedure update_halo2_dbl
    end interface

    !-----------------------------------------------------------------------
    !
    !  module variables
    !
    !-----------------------------------------------------------------------
    ! !PUBLIC DATA MEMBERS:
    integer (i4), dimension (:,:,:), allocatable, &
        public ::  &
        gdof,      & ! global degrees of freedom
        gdof_read, & ! global degrees of freedom (read in )
        ldof         ! local degrees of freedom

    real(r8), dimension(:,:,:), allocatable,  &
        public :: BarotropicMask

    integer (i4), public :: & ! number of linear points
        max_linear

    integer(i4), public :: num_linear
    integer(i4), public :: num_linear2

    integer(i4), public :: nTotal, nActive

    real(r8), dimension(:),allocatable, public :: LinearMask

    !***********************************************************************
    contains

    subroutine Convert2DtoLinear_dbl(dest,src)
        ! !DESCRIPTION:
        !  This subroutine converts a 2D data structure into 
        !  a linear array with the land points removed 
        !
        ! !REMARKS:
        !  This is the specific inteface for double precision arrays
        !  corresponding to the generic interface gather_global.  It is shown
        !  to provide information on the generic interface (the generic
        !  interface is identical, but chooses a specific inteface based
        !  on the data type of the input argument).

        ! !INPUT PARAMETERS:
        real(r8), intent(in), &
        dimension(nx_block,ny_block,max_blocks_tropic) &
            :: src  ! The 2D source array 

        ! !OUTPUT PARAMETERS:
        real(r8), intent(inout), & 
            dimension(max_linear) :: dest    ! The destination linear array

        !-----------------------------------------------------------------------
        !
        !  local variables
        !
        !-----------------------------------------------------------------------
        integer (i4) :: &
            i,j,ie,ib,je,jb,    &! dummy loop counters
            iblock,ig,npoints     ! loop tempories

        !-----------------------------------------------------------------------
        !
        !  copy local array into block decomposition
        !
        !-----------------------------------------------------------------------
        do iblock=1,nblocks_tropic
            call get_block_parameter( &
                blocks_tropic(iblock),npoints=npoints,jb=jb,je=je,ib=ib,ie=ie)
            if(npoints > 0) then
                do j=jb,je
                    do i=ib,ie
                        ig = ldof(i,j,iblock)
                        if( ig .gt. 0) dest(ig) = src(i,j,iblock)
                    enddo
                enddo
            endif
        enddo
    end subroutine Convert2DtoLinear_dbl


    subroutine ConvertLinearto2D_dbl(dest, src)
        real(r8), dimension(nx_block,ny_block,max_blocks_tropic) :: dest
        real(r8), intent(in), dimension(max_linear) :: src
        integer(i4) :: ie,ib,je,jb,i,j,iblock,ig, npoints

        do iblock=1,nblocks_tropic
            call get_block_parameter( &
                blocks_tropic(iblock),npoints=npoints,jb=jb,je=je,ib=ib,ie=ie)
            if(npoints > 0) then
                do j=jb,je
                    do i=ib,ie
                        ig = ldof(i,j,iblock)
                        if(ig .gt. 0) dest(i,j,iblock) = src(ig)
                    enddo
                enddo
            endif
        enddo
    end subroutine ConvertLinearto2D_dbl


    subroutine initDOF(ofname,strip,MASK)
        character(len=80) :: ofname
        logical, intent(in) :: strip
        real(r8), intent(in), &
        dimension(nx_block,ny_block,max_blocks_tropic) :: &
            MASK   ! land mask in arotropic distribution 

        integer(i4) :: ig,il,iblock,gbid
        integer(i4) :: i,j,ib,ie,jb,je, npoints
        integer(i4) :: bid
        logical, parameter :: Debug = .FALSE.
        character(len=80) :: fname

        integer(i4),allocatable :: vartmp(:,:,:)

        integer(i4) :: gid,n

        allocate(gdof(nx_block,ny_block,max_blocks_tropic))
        allocate(gdof_read(nx_block,ny_block,max_blocks_tropic))
        allocate(ldof(nx_block,ny_block,max_blocks_tropic))
        allocate(BarotropicMask(nx_block,ny_block,max_blocks_tropic))
        max_linear = block_size_x*block_size_y*max_blocks_tropic
        allocate(LinearMask(max_linear))

        ! Setup the BarotropicMask 
        do iblock = 1,nblocks_tropic
            gbid = blocks_tropic(iblock)

            !=====================
            ! Setup the local mask
            !=====================
            call local_mask(BarotropicMask(:,:,iblock), &
            distrb_tropic,field_loc_center,iblock, MASK(:,:,iblock))

            !======================================================
            ! Set the actual number of active ocean points in block
            !======================================================
            if(strip) then 
                npoints=COUNT(BarotropicMask(:,:,iblock) .gt. 0.0D0)
            else
                npoints=nx_block*ny_block
            endif

            !===================================
            ! set the npoints value in the block
            !===================================
            call set_block_parameter(gbid,npoints=npoints)
        enddo

        ig=1
        il=1
        do iblock=1,nblocks_tropic
            ! ================================================
            ! Zero out the global and local degrees of freedom
            ! ================================================
            gdof(:,:,iblock) = 0
            ldof(:,:,iblock) = 0

            ! =================
            ! Info on the block
            ! =================
            gbid = blocks_tropic(iblock)
            call get_block_parameter(gbid,ib=ib,ie=ie,jb=jb,je=je)

            if(strip) then 
                do j=jb,je
                    do i=ib,ie
                        if(BarotropicMask(i,j,iblock) .gt. 0.0D0) then
                            gdof(i,j,iblock) = ig
                            ldof(i,j,iblock) = il
                            ig=ig+1;il=il+1
                        endif
                    enddo
                enddo
            else
                do j=jb,je
                    do i=ib,ie
                        gdof(i,j,iblock) = ig
                        ldof(i,j,iblock) = il
                        ig=ig+1;il=il+1
                    enddo
                enddo
            endif
        enddo
        nActive = il-1
        nTotal = nActive
        num_linear=ig-1
        num_linear2 = num_linear

        call boundary_2d(gdof,bndy_tropic, field_loc_center, field_type_scalar)

        !JMD  stop 'after first call to mpi2s_boundary_2d'
        call boundary_2d(ldof,bndy_tropic, field_loc_center, field_type_scalar)

        bid=1
        call get_block_parameter(bid,ie=ie,ib=ib,je=je,jb=jb)

        bid=21
        call get_block_parameter(bid,ie=ie,ib=ib,je=je,jb=jb)

        if(Debug) then
            print *,'init_solvers: after update_ghost_cells'
            do iblock=1,nblocks_tropic
                gbid = blocks_tropic(iblock)
                call get_block_parameter(gbid,npoints=npoints)
                print *,'Block #:',gbid,' Ocean points: ', npoints
            enddo
        endif

        call Convert2DtoLinear(LinearMask,BarotropicMask)

        allocate(vartmp(nx_block,ny_block,nblocks_tot))
        vartmp(:,:,:)=0
        do n=1,nblocks_tropic
            gid = blocks_tropic(n)
            vartmp(:,:,gid) = gdof(:,:,n)
        enddo
        call write_tiles(ofname,'GDOF',vartmp)

        !-----------------------------------------------
        ! Lets read it back to check that it is the same
        !-----------------------------------------------
        call read_tiles(ofname,'GDOF',vartmp)
        do n=1,nblocks_tropic
            gid = blocks_tropic(n)
            gdof_read(:,:,n) = vartmp(:,:,gid)
        enddo

        deallocate(vartmp)
        print *,'SUM(gdof-gdof_read): ',SUM(gdof-gdof_read)
        deallocate(gdof,gdof_read,ldof)
        deallocate(BarotropicMask,LinearMask)
    end subroutine initDOF


    subroutine update_halo_dbl(n,array,mask)
        integer(i4), intent(in) :: n
        real(r8), dimension(:), intent(inout) :: array
        real(r8), dimension(:), intent(in) :: mask
    end subroutine update_halo_dbl

    subroutine update_halo2_dbl(n,iptrHalo,array)
        integer(i4), intent(in) :: n
        integer(i4), intent(in) :: iptrHalo
        real(r8), dimension(:), intent(inout) :: array
    end subroutine update_halo2_dbl

    subroutine local_mask(NewMask, dist, field_loc, n, MASK)
        ! !DESCRIPTION:
        !  computes the local sum of a single block
        !  array.
        !
        ! !REMARKS:
        !  This is actually the specific interface for the generic local_mask
        !  function corresponding to double precision arrays.  The generic
        !  interface is identical but will handle real and integer 2-d slabs
        !  and real, integer, and double precision scalars.

        ! !INPUT PARAMETERS:

        real (r8), dimension(:,:), intent(inout) :: &
            NewMask                   ! array to be summed

        type (distrb), intent(in) :: &
            dist                 ! block distribution for array X

        integer (i4), intent(in) :: &
            field_loc            ! location of field on staggered grid

        integer (i4), intent(in) :: &
            n                    ! local block index
            real (r8), dimension(size(NewMask,dim=1), &
            size(NewMask,dim=2)), intent(in) :: &
            MASK                 ! real multiplicative mask

        !-----------------------------------------------------------------------
        !
        !  local variables
        !
        !-----------------------------------------------------------------------
        real (r8) :: &
        local_block_sum      ! sum of local block domain

        integer (i4) :: &
        i,j,               &! local counters
        ib,ie,jb,je,       &! beg,end of physical domain
        bid                 ! block location

        NewMask = c0

        bid = dist%local_block(n)
        call get_block_parameter(n,ib=ib,ie=ie,jb=jb,je=je)
        local_block_sum = c0
        do j=jb,je
            do i=ib,ie
                NewMask(i,j) = MASK(i,j)
            end do
        end do
    end subroutine local_mask
end module linear
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
