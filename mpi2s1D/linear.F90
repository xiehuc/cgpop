!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! This module contains routines for converting to and from 
!! the linear data structure in the solver. 
!<
module linear
    use kinds_mod,     only: i4, r8
    use communicate,   only: MPI_COMM_OCN, my_task, get_num_procs, master_task
    use simple_blocks, only: nx_block,ny_block,get_proc_id, &
        set_block_parameter, get_block_parameter, &
        iwest,ieast,isouth,inorth,iswest,iseast,inwest,ineast
    use reductions,    only: global_sum, global_minval, global_maxval
    use simple_domain, only: nblocks_tropic, blocks_tropic, &
        distrb_tropic, read_local_tiles
    use domain_size,   only: max_blocks_tropic, block_size_x, block_size_y
    use mpi2s_gshalo,  only: mpi2s_gshalo_update_2d_dbl, &
        mpi2s_gshalo_update_2d_int, mpi2s_gshalo_init

    implicit none
    private
    save

    !-----------------------------------------------------------------------
    !
    !  overload module functions
    !
    !-----------------------------------------------------------------------
    ! !PUBLIC MEMBER FUNCTIONS:
    public :: Convert2DtoLinear
    public :: ConvertLinearto2D

    public :: CreateIndx2DtoLinear
    public :: initDOF
    public :: update_halo

    !-----------------------------------------------------------------------
    !
    !  module variables
    !
    !-----------------------------------------------------------------------
    ! !PUBLIC DATA MEMBERS:
    integer (i4), dimension (nx_block,ny_block,max_blocks_tropic), &
    public ::    &
        gdof,    &! global degrees of freedom
        ldof      ! local degrees of freedom

    integer (i4), parameter, public :: & ! number of linear points
        max_linear = nx_block*ny_block*max_blocks_tropic

    integer (i4),dimension(max_linear), public :: LinearGdof, LinearProc

    real(r8), dimension(max_linear), public :: LinearMask

    integer(i4),dimension(3,max_linear), public :: Indx2DtoLinear

    integer(i4), public :: nActive,nTotal

    integer(i4), public :: gsHalo_handle


  contains

    !***********************************************************************
    !>
    !! This subroutine converts a 2D data structure into  a linear array with
    !! the land points removed dest,src)
    !!
    !! @param dest  1D destination array 
    !! @param src   2D source array
    !<
    subroutine Convert2DtoLinear(dest,src)
        real(r8), intent(in), &
            dimension(nx_block,ny_block,max_blocks_tropic) :: src

        ! !OUTPUT PARAMETERS:
        real(r8), intent(inout), & 
            dimension(max_linear) :: dest    ! The destination linear array

        !-------------------------------------------------------------------
        !
        !  local variables
        !
        !-------------------------------------------------------------------
        integer (i4) :: i,j,ii,ib    ! dummy loop counters

        !-------------------------------------------------------------------
        !
        !  copy local array into block decomposition
        !
        !-------------------------------------------------------------------
        do ii=1,nActive
            i  = Indx2DtoLinear(1,ii)
            j  = Indx2DtoLinear(2,ii)
            ib = Indx2DtoLinear(3,ii)
            dest(ii)=src(i,j,ib)
        enddo
    end subroutine Convert2DtoLinear

    !***********************************************************************
    !>
    !! This subroutine creates the indirect addressing that converts the
    !! the 2-dimensional variable into 1-dimensional equivalent.
    !<
    subroutine CreateIndx2DtoLinear
        integer (i4) :: &
        i,j,ie,ib,je,jb,    &! dummy loop counters
        iblock,ig,npoints    ! loop tempories

        do iblock=1,nblocks_tropic
            call get_block_parameter( &
                blocks_tropic(iblock),npoints=npoints,jb=jb,je=je,ib=ib,ie=ie)

            if(npoints > 0) then
                do j=jb,je
                    do i=ib,ie
                        ig = ldof(i,j,iblock)
                        if( ig .gt. 0) then 
                            Indx2DtoLinear(1,ig) = i
                            Indx2DtoLinear(2,ig) = j
                            Indx2DtoLinear(3,ig) = iblock
                        endif
                    enddo
                enddo
            endif
        enddo
    end subroutine CreateIndx2DtoLinear

    !***********************************************************************
    !>
    !! This subroutine converts a linear array into a 2D data-structure
    !!
    !! @param dest  2D target array
    !! @param src   1D source vector
    !<
    subroutine ConvertLinearto2D(dest, src)
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
                        if(ig .gt. 0) &
                            dest(i,j,iblock) = src(ig)
                    enddo
                enddo
            endif
        enddo
    end subroutine ConvertLinearto2D

    !***********************************************************************
    !>
    !! This subroutine calculates the index arrays necessary to construct the
    !! halo update for 1-dimensional variables.
    !!
    !! @param fname   The name of the netCDF file which contains the metadata
    !!                necessary to construct the 1D halo update.
    !! @param strip   A logical variable which controls if land points are
    !!                removed from the 1D datastructure. [REMOVE]
    !! @param MASK    land mask in barotropic distribution 
    !>
    subroutine initDOF(fname,strip,MASK)
        character(len=80) :: fname
        logical, intent(in) :: strip
        real(r8), intent(in), &
            dimension(nx_block,ny_block,max_blocks_tropic) :: MASK

        integer (i4),allocatable :: ldof_linear(:)

        integer(i4) :: ig,il,iblock,gbid
        integer(i4) :: i,j,ib,ie,jb,je
        logical, parameter :: Debug = .FALSE.

        integer(i4) :: lsum,gsum

        integer(i4) :: nprocs

        integer :: ierr

        logical, parameter :: CheckGSHalo = .TRUE.

        integer(i4) :: nActive_min,nActive_max,nActive_avg,nActive_sum
        integer(i4) :: nTotal_sum, nTotal_tmp, nblocks_tmp

        integer(i4),dimension(max_linear)  :: LinearGdof2,LinearDiff

        integer(i4) :: npoints

        integer(i4) :: pid_South, pid_North, pid_East, pid_West, &
            pid_SouthEast, pid_SouthWest, pid_NorthEast, pid_NorthWest

        integer(i4) :: n

        nprocs = get_num_procs()

        LinearMask=0.d0
        do iblock = 1,nblocks_tropic
            gbid = blocks_tropic(iblock)
            call get_block_parameter(gbid,ib=ib,ie=ie,jb=jb,je=je)

            !======================================================
            ! Set the actual number of active ocean points in block
            !======================================================
            if(strip) then 
                npoints=COUNT(MASK(ib:ie,jb:je,iblock) .gt. 0.0D0)
            else
                npoints=nx_block*ny_block
            endif

            !===================================
            ! set the npoints value in the block
            !===================================
            call set_block_parameter(gbid,npoints=npoints)
        enddo

        LinearGdof=0
        LinearProc=my_task

        gdof(:,:,:) = 0
        !---------------------------------------------------
        ! read the global degree of freedom in from the file 
        !---------------------------------------------------
        call read_local_tiles(fname,'GDOF',gdof)

        il=1
        do iblock=1,nblocks_tropic
            ! ========================
            ! local degrees of freedom
            ! ========================
            ldof(:,:,iblock) = 0

            ! =================
            ! Info on the block
            ! =================
            gbid = blocks_tropic(iblock)
            call get_block_parameter(gbid,ib=ib,ie=ie,jb=jb,je=je)

            if(strip) then 
                do j=jb,je
                    do i=ib,ie
                        ig=gdof(i,j,iblock)
                        if(MASK(i,j,iblock) .gt. 0.0D0) then
                            ldof(i,j,iblock) = il
                            LinearMask(il)   = 1.0D0
                            ! Don't set the GDOF for interior points  
                            if((i .eq. ib) .or. (i .eq. ie) .or. &
                                (j .eq. jb) .or. (j .eq. je)) then  
                                LinearGdof(il)   = ig
                            endif
                            il=il+1
                        endif
                    enddo
                enddo
            else
                do j=jb,je
                    do i=ib,ie
                        ig=gdof(i,j,iblock)
                        ldof(i,j,iblock) = il
                        LinearMask(il)   = 1.0D0
                        LinearGdof(il)   = ig
                        il=il+1
                    enddo
                enddo
            endif
        enddo
        nActive=il-1

        do iblock=1,nblocks_tropic
            ! =================
            ! Info on the block
            ! =================
            gbid = blocks_tropic(iblock)
            call get_block_parameter(gbid,ib=ib,ie=ie,jb=jb,je=je)
            
            ! =================================
            ! Figure out the GDOFS for the halo
            ! =================================
            j=jb-1
            
            !--------------------------
            ! do the Southwest Neighbor
            !---------------------------
            pid_SouthWest=get_proc_id(distrb_tropic,iswest,gbid)
            call FormHaloDofCk( &
                ib-1,j,iblock,pid_SouthWest,gdof,ldof,LinearGdof,LinearProc,il)
            
            !----------------------
            ! do the South neighbor
            !----------------------
            pid_South=get_proc_id(distrb_tropic,isouth,gbid)
            do i = ib,ie
                call FormHaloDofCk( &
                    i,j,iblock,pid_South,gdof,ldof,LinearGdof,LinearProc,il)
            enddo
            
            !--------------------------
            ! do the Southeast Neighbor
            !---------------------------
            pid_SouthEast=get_proc_id(distrb_tropic,iseast,gbid)
            call FormHaloDofCk( &
                ie+1,j,iblock,pid_SouthEast,gdof,ldof,LinearGdof,LinearProc,il)
            j=je+1
            
            !--------------------------
            ! do the Northwest neighbor
            !--------------------------
            pid_NorthWest=get_proc_id(distrb_tropic,inwest,gbid)
            call FormHaloDofCk( &
                ib-1,j,iblock,pid_NorthWest,gdof,ldof,LinearGdof,LinearProc,il)
            
            !--------------------------
            ! do the North neighbor
            !--------------------------
            pid_North = get_proc_id(distrb_tropic,inorth,gbid)
            do i = ib,ie
                call FormHaloDofCk( &
                i,j,iblock,pid_North,gdof,ldof,LinearGdof,LinearProc,il)
            enddo
            
            !--------------------------
            ! do the Northeast neighbor
            !--------------------------
            i = ie+1
            pid_NorthEast = get_proc_id(distrb_tropic,ineast,gbid)
            call FormHaloDofCk( &
                ie+1,j,iblock,pid_NorthEast,gdof,ldof,LinearGdof,LinearProc,il)

            i=ib-1
            
            !--------------------------
            ! do the West neighbor
            !--------------------------
            pid_West = get_proc_id(distrb_tropic,iwest,gbid)
            do j = jb,je
                call FormHaloDofCk( &
                    i,j,iblock,pid_West,gdof,ldof,LinearGdof,LinearProc,il)
            enddo
            i=ie+1
            
            !--------------------------
            ! do the East neighbor
            !--------------------------
            pid_East = get_proc_id(distrb_tropic,ieast,gbid)
            do j = jb,je
                call FormHaloDofCk( &
                    i,j,iblock,pid_East,gdof,ldof,LinearGdof,LinearProc,il)
            enddo
        enddo
        nTotal = il-1

        nprocs = get_num_procs()

        !========================================
        ! Calculate stats for active ocean points 
        !========================================
        nActive_max = global_maxval(nActive)
        nActive_min = global_minval(nActive)
        nActive_sum = global_sum(nActive)
        nActive_avg = nActive_sum/nprocs

        ! ================================================
        ! Calculate stats for total points (active + halo)
        ! ================================================
        nTotal_tmp = block_size_x*block_size_y*nblocks_tropic
        nTotal_sum  = global_sum(nTotal_tmp)

        ! =================================
        ! Calculate stats for total blocks 
        ! =================================
        nblocks_tmp = global_sum(nblocks_tropic)

        if(my_task == master_task) then 
            write(6,*) 'initDOF: Active Ocean points {min,max,avg}: ', &
                nActive_min,nActive_max,nActive_avg
            write(6,*) 'initDOF: # of {total,active} Ocean points (global): ', &
                nTotal_sum, nActive_sum
            write(6,*) 'initDOF: # of total blocks: ', nblocks_tmp
        endif

        ! ================================
        ! form the communicator handle
        ! ================================
        gsHalo_handle = mpi2s_gshalo_init( &
            MPI_COMM_OCN,max_linear,nTotal,nActive,LinearGdof,LinearProc)

        if(CheckGSHalo) then 
            ! ============================
            ! Send through a test pattern 
            ! ============================
            LinearGdof2 = 0
            LinearDiff  = 0
            LinearGdof2(1:nActive) = LinearGdof(1:nActive)
            call mpi2s_gshalo_update_2d_int(gshalo_handle,LinearGdof2)

            LinearDiff(nActive+1:nTotal) = &
                LinearGdof(nActive+1:nTotal) - LinearGdof2(nActive+1:nTotal)

            lsum = SUM(LinearDiff)
            lsum = global_sum(lsum)
            if(master_task == my_task) then 
                print *,'initDOF: SUM(LinearDiff) := ',lsum 
            endif
        endif

        !-------------------------------------------------------------
        ! Create the index to convert from 2D to Linear data structure
        !-------------------------------------------------------------
        call CreateIndx2DtoLinear()
    end subroutine initDOF

    !***********************************************************************
    !>
    !! This calculates the global degree of freedom for constructing the 1D
    !! halo update.
    !!
    !! @param i         The 'i' index for the gridpoint.
    !! @param j         The 'j' index for the gridpoint.
    !! @param iblock    The local block index which owns the gridpoint. 
    !! @param pnum      The processor or task number which owns the gridpoint.
    !! @param gdof      The global degree of freedom in the 2D data-structure.
    !! @param ldof      The local 1D ordering for the gridpoint.
    !! @param linGdof   A 1D array which contains the global degree of freedom. 
    !! @param linProc   A 1D array which contains the task number which
    !!                  corresponds to the global degree of freedom.
    !! @param il        The pointer to the end last inserted GDOF in the
    !!                  linGdof array. 
    !>
    subroutine FormHaloDofCk(i,j,iblock,pnum, gdof,ldof,linGdof,linProc,il)
        integer(i4),intent(in) :: i,j,iblock
        integer(i4),intent(in) :: pnum
        integer(i4),intent(inout), &
            dimension(nx_block,ny_block,max_blocks_tropic)  :: gdof,ldof
        integer(i4), intent(inout) :: linGdof(max_linear),linProc(max_linear)
        integer(i4), intent(inout) :: il

        integer(i4) :: idx, ierr
        logical :: found

        if(pnum .ne. my_task) then
            if(gdof(i,j,iblock) .gt. 0) then
                call LinearFind( &
                    linGdof(nActive+1:il),gdof(i,j,iblock),found,idx)
                if(found) then 
                    ldof(i,j,iblock)=nActive+idx
                else
                    linGdof(il) = gdof(i,j,iblock)
                    linProc(il) = pnum
                    ldof(i,j,iblock) = il
                    il=il+1
                endif
            endif
        else 
            if(gdof(i,j,iblock) .gt. 0) then 
                call LinearFind(linGdof(1:nActive),gdof(i,j,iblock),found,idx)
                if(found) then 
                ldof(i,j,iblock) = idx
                endif
            endif
        endif
    end subroutine FormHaloDofCk

    !***********************************************************************
    !>
    !! This subroutine determines if a scalar value already exists in an array.
    !!
    !! @param array  A 1D integer array to be searched.
    !! @param value  The scalar value for which to search.
    !! @param found  A logical value which indicates the scalar value is found.
    !! @param indx   If found=.true. the location of the scalar in the array.
    !>
    subroutine LinearFind(array,value,found,indx)
        integer(i4), intent(in) :: array(:)
        integer(i4), intent(in) :: value
        logical, intent(out) :: found
        integer(i4), intent(out) :: indx

        integer(i4) :: i, n
        n = SIZE(array)
        found = .FALSE.
        do i=1,n
            if(array(i) == value) then
                found = .TRUE.
                indx = i
                return
            endif
        enddo
    end subroutine LinearFind

    !***********************************************************************
    !>
    !! Wrapper function which calls a subroutine to update the halo of a
    !! 1-dimensional variable.
    !!
    !! @param array  The variable whose halo should be updated. 
    !>
    subroutine update_halo(array)
        real(r8), dimension(:), intent(inout) :: array

        call mpi2s_gshalo_update_2d_dbl(gshalo_handle, array)
    end subroutine update_halo
end module linear

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
