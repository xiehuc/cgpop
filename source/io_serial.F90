!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!!  This module provides several methods to read and write 
!!  data needed by the miniapp
!<
module io_serial
    use kinds_mod, only: i4, r8
    use simple_blocks, only: AppMD_t
    use netcdf
    use exit_mod

    implicit none
    private

    public :: write_AppMD, read_AppMD  ! read and write application metadata

    public :: read_tiles

    !>
    !! Reads in variables decomposed into tiles
    !<
    interface read_tiles
        module procedure read_tiles_int
        module procedure read_tiles_dbl
    end interface

    public :: write_tiles

    !>
    !! Writes out variables decomposed into tiles
    !<
    interface write_tiles
        module procedure write_tiles_int
        module procedure write_tiles_dbl
    end interface

    contains 

    !>
    !! @details Reads in a double precision variable from netCDF tilefile. Note 
    !! that this subroutine reads in all tiles. 
    !! @param fname Filename of the tilefile from which to read the variable.
    !! @param vname Name of the variable to read.
    !! @param var   Double precision variable that is returned.
    !<
    subroutine read_tiles_dbl(fname,vname,var)
        character(len=80), intent(in) :: fname
        character(len=*), intent(in)  :: vname
        real(r8), dimension(:,:,:) :: var

        integer(i4) :: mode,fh,varid,ierr
        integer(i4) :: nx,ny
        integer(i4) :: nx_read,ny_read
        integer(i4) :: dimid_nx, dimid_ny
        character(len=80) :: error


        mode = NF90_nowrite
        nx = SIZE(var,dim=1)
        ny = SIZE(var,dim=2)
        ierr = sio_open(fname,mode,fh)

        ierr = nf90_inq_dimid(fh,'nx_block',dimid_nx)
        ierr = nf90_inquire_dimension(fh,dimid_nx,len=nx_read)
        ierr = nf90_inq_dimid(fh,'ny_block',dimid_ny)
        ierr = nf90_inquire_dimension(fh,dimid_ny,len=ny_read)
        if((nx_read .ne. nx) .or. (ny_read .ne. ny)) then 
        print *,'Error: mismatch between size of blocks'
        print *,'nx:= ',nx,'nx_read:= ',nx_read
        print *,'ny:= ',ny,'ny_read:= ',ny_read
        stop 'in read_tiles:dbl'
        endif


        ierr = nf90_inq_varid(fh,TRIM(vname),varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid: '//vname,__LINE__)
        ierr = nf90_get_var(fh,varid, var)
        if(ierr /= nf90_noerr) &
        call check_netcdf(fh,ierr,'nf90_get_var: '//vname,__LINE__)

        ierr = sio_close(fh) 
    end subroutine read_tiles_dbl

    !>
    !!
    !! Reads in a integer variable from netCDF tile file. Note 
    !! that this subroutine reads in all tiles. 
    !! 
    !! @param fname Filename of the tilefile from which to read the variable.
    !! @param vname Name of the variable to read 
    !! @param var   Integer variable that is returned.
    !<
    subroutine read_tiles_int(fname,vname,var)
        character(len=80), intent(in) :: fname
        character(len=*), intent(in)  :: vname
        integer(i4), dimension(:,:,:) :: var

        integer(i4) :: mode,fh,varid,ierr
        integer(i4) :: nx,ny
        integer(i4) :: nx_read,ny_read
        integer(i4) :: dimid_nx, dimid_ny

        mode = NF90_nowrite
        nx = SIZE(var,dim=1)
        ny = SIZE(var,dim=2)
        ierr = sio_open(fname,mode,fh)

        ierr = nf90_inq_dimid(fh,'nx_block',dimid_nx)
        ierr = nf90_inquire_dimension(fh,dimid_nx,len=nx_read)
        ierr = nf90_inq_dimid(fh,'ny_block',dimid_ny)
        ierr = nf90_inquire_dimension(fh,dimid_ny,len=ny_read)
        if((nx_read .ne. nx) .or. (ny_read .ne. ny)) then 
            print *,'Error: mismatch between size of blocks'
            print *,'nx:= ',nx,'nx_read:= ',nx_read
            print *,'ny:= ',ny,'ny_read:= ',ny_read
            stop 'in read_tiles:dbl'
        endif

        ierr = nf90_inq_varid(fh,TRIM(vname),varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid: '//vname,__LINE__)

        ierr = nf90_get_var(fh,varid, var)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var: '//vname,__LINE__)

        ierr = sio_close(fh) 
    end subroutine read_tiles_int

    !>
    !!
    !! Writes out an integer variable to a netCDF tile file. Note 
    !! that this subroutine writes out all tiles. 
    !! 
    !! @param fname Filename of the tilefile to which to write the variable.
    !! @param vname Name of the variable to write 
    !! @param var   Integer variable to write.
    !<
    subroutine write_tiles_int(fname,vname,var)
        character(len=80), intent(in) :: fname
        character(len=*), intent(in) :: vname
        integer(i4) :: var(:,:,:)

        integer(i4) :: dim_nx,dim_ny,dim_nb
        integer(i4) :: ncid, fh, ierr, mode
        integer(i4) :: nx_block,ny_block,nblocks

        mode = NF90_write
        ierr = sio_open(fname,mode,fh)

        nx_block = SIZE(var,dim=1)
        ny_block = SIZE(var,dim=2)
        nblocks  = SIZE(var,dim=3) 

        ierr = nf90_redef(fh)
        ierr = nf90_inq_dimid(fh,'nblocks',dimid=dim_nb)
        if(ierr /= nf90_noerr) then 
            ierr = nf90_def_dim(fh,name='nblocks',len=nblocks,dimid=dim_nb)
            if(ierr /= nf90_noerr) &
                call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
        endif

        ierr = nf90_inq_dimid(fh,'nx_block',dimid=dim_nx)
        if(ierr /= nf90_noerr) then 
            ierr = nf90_def_dim(fh,name='nx_block',len=nx_block,dimid=dim_nx)
            if(ierr /= nf90_noerr) &
                call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
        endif

        ierr = nf90_inq_dimid(fh,'ny_block',dimid=dim_ny)
        if(ierr /= nf90_noerr) then 
            ierr = nf90_def_dim(fh,name='ny_block',len=ny_block,dimid=dim_ny)
            if(ierr /= nf90_noerr) &
                call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
        endif

        ierr = nf90_def_var( &
            fh,name=TRIM(vname),xtype=NF90_INT,dimids= &
                (/dim_nx,dim_ny,dim_nb/),varid=ncid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        !------------------------
        ! end netCDF define mode
        !------------------------
        ierr = nf90_enddef(fh)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_enddef',__LINE__)

        !--------------------------------
        ! write the variables to the file
        !---------------------------------
        ierr = nf90_put_var(fh,ncid,var)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_close(fh)
    end subroutine write_tiles_int

    !>
    !!
    !! Writes out an double precision variable to a netCDF tile file. Note 
    !! that this subroutine writes out all tiles. 
    !! 
    !! @param fname Filename of the tilefile to which to write the variable.
    !! @param vname Name of the variable to write 
    !! @param var   Double precision variable to write.
    !<
    subroutine write_tiles_dbl(fname,vname,var)
        character(len=80), intent(in) :: fname
        character(len=*), intent(in) :: vname
        real(r8) :: var(:,:,:)

        integer(i4) :: dim_nx,dim_ny,dim_nb
        integer(i4) :: ncid,fh,ierr, mode
        integer(i4) :: nx_block,ny_block,nblocks

        mode = NF90_write
        ierr = sio_open(fname,mode,fh)

        nx_block = SIZE(var,dim=1)
        ny_block = SIZE(var,dim=2)
        nblocks  = SIZE(var,dim=3)

        ierr = nf90_redef(fh)
        ierr = nf90_inq_dimid(fh,'nblocks',dimid=dim_nb)
        if(ierr /= nf90_noerr) then
            ierr = nf90_def_dim(fh,name='nblocks',len=nblocks,dimid=dim_nb)
            if(ierr /= nf90_noerr) &
                call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
        endif

        ierr = nf90_inq_dimid(fh,'nx_block',dimid=dim_nx)
        if(ierr /= nf90_noerr) then
            ierr = nf90_def_dim(fh,name='nx_block',len=nx_block,dimid=dim_nx)
            if(ierr /= nf90_noerr) &
                call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
        endif

        ierr = nf90_inq_dimid(fh,'ny_block',dimid=dim_ny)
        if(ierr /= nf90_noerr) then
            ierr = nf90_def_dim(fh,name='ny_block',len=ny_block,dimid=dim_ny)
            if(ierr /= nf90_noerr) &
                call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
        endif

        ierr = nf90_def_var( &
            fh,name=TRIM(vname),xtype=NF90_DOUBLE, &
            dimids=(/dim_nx,dim_ny,dim_nb/),varid=ncid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        !------------------------
        ! end netCDF define mode
        !------------------------
        ierr = nf90_enddef(fh)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_enddef',__LINE__)

        !--------------------------------
        ! write the variables to the file
        !---------------------------------
        ierr = nf90_put_var(fh,ncid,var)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_close(fh)
    end subroutine write_tiles_dbl

    !>
    !! This subroutine writes application specific metadata to the tilefiles.
    !! Note that this metadata is common among all implementations of the
    !! miniapp
    !!
    !! @param fname   Filename which to write the application metadata.   
    !! @param blocks  The blocks data structure which contains the application
    !!                metadata
    !! @param reorder An interger array which reorder the cannonical block
    !!                layout  into a space-filling curve based layout.  
    !<
    subroutine write_AppMD(fname,blocks,reorder)
        character(len=80), intent(in) :: fname
        type (AppMD_t) :: blocks(:)
        integer (i4) :: reorder(:)

        integer(i4) :: mode, fh, ierr
        integer(i4) :: dim_nb,dim_nx,dim_ny
        integer(i4) :: ncid_block_id,ncid_ib,ncid_ie,ncid_jb,ncid_je, &
            ncid_iblock,ncid_jbloc,ncid_jblock,ncid_npoints,ncid_iglob, &
            ncid_jglob, ncid_reorder
        integer(i4), allocatable :: iglob(:,:), jglob(:,:)
        integer(i4) :: nblocks, iblock, jblock, i, nx_block, ny_block

        mode = NF90_CLOBBER
        ierr = sio_create(fname,mode,fh)

        nblocks = SIZE(blocks)
        nx_block = SIZE(blocks(1)%i_glob(:))
        ny_block = SIZE(blocks(1)%j_glob(:))

        allocate(iglob(nx_block,nblocks))
        allocate(jglob(ny_block,nblocks))

        !---------------------------------------------------------
        ! define the dimensions for the file
        !---------------------------------------------------------
        ierr = nf90_def_dim(fh,name='nblocks',len=nblocks,dimid=dim_nb)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
        ierr = nf90_def_dim(fh,name='nx_block',len=nx_block,dimid=dim_nx)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
        ierr = nf90_def_dim(fh,name='ny_block',len=ny_block,dimid=dim_ny)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)

        !---------------------------------------------------------
        ! define the variables which will be contained in the file
        !---------------------------------------------------------
        ierr = nf90_def_var( &
            fh,name='block_id', xtype = NF90_INT, &
            dimids=(/dim_nb/),varid=ncid_block_id)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='ib', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_ib)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='ie', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_ie)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='jb', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_jb)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='je', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_je)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='iblock', xtype = NF90_INT, dimids=(/dim_nb/), &
            varid=ncid_iblock)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='jblock', xtype = NF90_INT, dimids=(/dim_nb/), &
            varid=ncid_jblock)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='npoints', xtype = NF90_INT, dimids=(/dim_nb/), &
            varid=ncid_npoints)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='i_glob', xtype = NF90_INT, dimids=(/dim_nx,dim_nb/), &
            varid=ncid_iglob)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='j_glob', xtype = NF90_INT, dimids=(/dim_ny,dim_nb/), &
            varid=ncid_jglob)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        ierr = nf90_def_var( &
            fh,name='reorder', xtype = NF90_INT, dimids=(/dim_nb/), &
            varid=ncid_reorder)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

        !------------------------
        ! end netCDF define mode
        !------------------------
        ierr = nf90_enddef(fh)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_enddef',__LINE__)

        !--------------------------------
        ! write the variables to the file
        !---------------------------------
        ierr = nf90_put_var(fh,ncid_block_id,blocks(:)%block_id)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_put_var(fh,ncid_ib,blocks(:)%ib)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_put_var(fh,ncid_ie,blocks(:)%ie)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_put_var(fh,ncid_jb,blocks(:)%jb)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_put_var(fh,ncid_je,blocks(:)%je)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_put_var(fh,ncid_iblock,blocks(:)%iblock)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_put_var(fh,ncid_jblock,blocks(:)%jblock)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_put_var(fh,ncid_npoints,blocks(:)%npoints)
        if(ierr /= nf90_noerr)  &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)

        do i=1,nblocks
            iglob(:,i) = blocks(i)%i_glob(:)
            jglob(:,i) = blocks(i)%j_glob(:)
        enddo

        ierr = nf90_put_var(fh,ncid_iglob,iglob)        
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_put_var(fh,ncid_jglob,jglob)        
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_put_var(fh,ncid_reorder,reorder)        
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
        ierr = nf90_close(fh)
    end subroutine write_AppMD

    !>
    !! This subroutine reads application specific metadata from the tilefiles.
    !! Note that this metadata is common among all implementations of the
    !! miniapp
    !!
    !! @param fname   Filename from which to read the application metadata.   
    !! @param blocks  The blocks data structure which contains the application
    !! metadata
    !! @param reorder An interger array which reorders the cannonical block
    !! layout  into a space-filling curve based layout.  
    !<
    subroutine read_AppMD(fname,blocks,reorder)
        character(len=80), intent(in) :: fname
        type (AppMD_t), pointer, intent(inout)  :: blocks(:)
        integer (i4),pointer :: reorder(:)

        integer (i4) :: mode, fh, ierr, i, j
        integer (i4) :: nblocks, nx_block,ny_block
        integer (i4) :: dim_nb, dim_nx, dim_ny, varid

        integer (i4), allocatable :: iglob(:,:), jglob(:,:)

        mode = NF90_nowrite
        ierr = sio_open(fname,mode,fh)

        ierr = nf90_inq_dimid(fh,'nblocks',dimid=dim_nb)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_dimid',__LINE__)

        ierr = nf90_inq_dimid(fh,'nx_block',dimid=dim_nx)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_dimid',__LINE__)

        ierr = nf90_inq_dimid(fh,'ny_block',dimid=dim_ny)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_dimid',__LINE__)

        ierr = nf90_inquire_dimension(fh,dim_nb,len=nblocks)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inquire_dimension',__LINE__)

        ierr = nf90_inquire_dimension(fh,dim_nx,len=nx_block)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inquire_dimension',__LINE__)

        ierr = nf90_inquire_dimension(fh,dim_ny,len=ny_block)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inquire_dimension',__LINE__)

        allocate(blocks(nblocks))

        !--------------------------------------------------------------------
        ! Read and load the different componets of the AppMD_t data structure
        !--------------------------------------------------------------------
        ierr = nf90_inq_varid(fh,'block_id',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

        ierr = nf90_get_var(fh,varid, blocks(:)%block_id)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

        ierr = nf90_inq_varid(fh,'ib',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

        ierr = nf90_get_var(fh,varid, blocks(:)%ib)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

        ierr = nf90_inq_varid(fh,'ie',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

        ierr = nf90_get_var(fh,varid, blocks(:)%ie)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

        ierr = nf90_inq_varid(fh,'jb',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

        ierr = nf90_get_var(fh,varid, blocks(:)%jb)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

        ierr = nf90_inq_varid(fh,'je',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

        ierr = nf90_get_var(fh,varid, blocks(:)%je)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

        ierr = nf90_inq_varid(fh,'iblock',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

        ierr = nf90_get_var(fh,varid, blocks(:)%iblock)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

        ierr = nf90_inq_varid(fh,'jblock',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

        ierr = nf90_get_var(fh,varid, blocks(:)%jblock)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

        ierr = nf90_inq_varid(fh,'npoints',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

        ierr = nf90_get_var(fh,varid, blocks(:)%npoints)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

        allocate(iglob(nx_block,nblocks))
        allocate(jglob(ny_block,nblocks))

        ierr = nf90_inq_varid(fh,'i_glob',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
        ierr = nf90_get_var(fh,varid, iglob)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)
        ierr = nf90_inq_varid(fh,'j_glob',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
        ierr = nf90_get_var(fh,varid, jglob)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)
        
        do i=1,nblocks
            allocate(blocks(i)%i_glob(nx_block))
            allocate(blocks(i)%j_glob(ny_block))
            blocks(i)%i_glob(:) = iglob(:,i)
            blocks(i)%j_glob(:) = jglob(:,i)
        enddo
        
        deallocate(iglob,jglob)
        allocate(reorder(nblocks))
        
        ierr = nf90_inq_varid(fh,'reorder',varid)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
        
        ierr = nf90_get_var(fh,varid, reorder)
        if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)
    end subroutine read_AppMD

    !>
    !! Creates a netCDF file
    !!
    !! @param fname  Filename of the file to be created.
    !! @param mode   The netCDF file mode  
    !! @param fh     An integer which identifies the created netCDF file.
    !<
    function sio_create(fname,mode,fh)
        character(len=*), intent(in) :: fname
        integer(i4), intent(in) :: mode
        integer(i4), intent(out) :: fh
        integer(i4) :: sio_create
        integer(i4) :: fmode, ierr
        
        ierr = nf90_create(fname,mode,fh)
        ierr = nf90_set_fill(fh,NF90_NOFILL,fmode)
        
        sio_create = ierr
    end function sio_create

    !>
    !! Opens an existing netCDF file
    !!
    !! @param fname  The name of the file to ope.
    !! @param mode   The netCDF file mode  
    !! @param fh     An integer which identifies the netCDF file.
    !<
    function sio_open(fname,mode,fh)
        character(len=*), intent(in) :: fname
        integer(i4), intent(in) :: mode
        integer(i4), intent(out) :: fh
        integer(i4) :: sio_open
        integer(i4) :: ierr

        ierr = nf90_open(fname,mode,fh)
        if(ierr /= nf90_noerr) then
            call check_netcdf(fh,ierr,'nf90_open(): '//fname,__LINE__)
            call exit_POP(1, "Issue encountered while opening netCDF file")
        end if
        sio_open= ierr
    end function sio_open

    !>
    !! Closes an open netCDF file
    !!
    !! @param ncid The netCDF file identifier. 
    !<
    function sio_close(ncid)
        integer(i4), intent(inout) :: ncid
        integer(i4) :: sio_close
        integer(i4) :: ierr

        ierr = nf90_sync(ncid)
        ierr = nf90_close(ncid)

        sio_close = ierr
    end function sio_close

    !>
    !! This subroutine will print out a netCDF error message
    !!
    !! @param fh       The netCDF file identifier. 
    !! @param status   The error return code.
    !! @param filestr  A user supplied string to identify the location of the
    !!                 netCDF error
    !! @param line     The file line number which generated the error  
    !<
    subroutine check_netcdf(fh,status,filestr,line)
        integer(i4) :: fh
        integer(i4) :: status
        character(len=*), intent(in) :: filestr
        integer, intent(in) :: line

        print *,'ERROR: ',TRIM(filestr),' ',TRIM(nf90_strerror(status)), &
            ' on line: ',line 
    end subroutine check_netcdf
end module io_serial
