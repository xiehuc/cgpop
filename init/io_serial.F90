!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

module io_serial

     use kinds_mod, only: i4, r8
     use blocks, only: AppMD_t
     use netcdf

     implicit none
     private

!     ierr = nf90_create(fname,nmode,fh)
!     ierr = nf90_set_fill(fh,NF90_NOFILL,nmode)

!     ierr = nf90_sync(fh)
!     ierr = nf90_close(fh)
 
!     ierr = nf90_get_var(fh,varid,IOBUF,temp_start,temp_count)

     
      public :: write_AppMD, read_AppMD

      public :: read_garray
      interface read_garray
        module procedure read_garray_int
        module procedure read_garray_dbl
      end interface

      
      public :: read_tiles
      interface read_tiles
	module procedure read_tiles_int
	module procedure read_tiles_dbl
      end interface

      public :: write_tiles
      interface write_tiles
	module procedure write_tiles_int
	module procedure write_tiles_dbl
      end interface


contains 

    subroutine read_tiles_dbl(fname,vname,var)
      
       character(len=80), intent(in) :: fname
       character(len=*), intent(in)  :: vname
       real(r8), dimension(:,:,:) :: var

       integer(i4) :: mode
       integer(i4) :: fh
       integer(i4) :: varid
       integer(i4) :: ierr

    
       mode = NF90_nowrite
       ierr = sio_open(fname,mode,fh)

       ierr = nf90_inq_varid(fh,TRIM(vname),varid)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      
       ierr = nf90_get_var(fh,varid, var)
       if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)
      
       ierr = sio_close(fh) 
         
    end subroutine read_tiles_dbl

    subroutine read_tiles_int(fname,vname,var)
      
       character(len=80), intent(in) :: fname
       character(len=*), intent(in)  :: vname
       integer(i4), dimension(:,:,:) :: var

       integer(i4) :: mode
       integer(i4) :: fh
       integer(i4) :: varid
       integer(i4) :: ierr

    
       mode = NF90_nowrite
       ierr = sio_open(fname,mode,fh)

       ierr = nf90_inq_varid(fh,TRIM(vname),varid)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      
       ierr = nf90_get_var(fh,varid, var)
       if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)
      
       ierr = sio_close(fh) 
         
    end subroutine read_tiles_int


    subroutine write_tiles_int(fname,vname,var)

      character(len=80), intent(in) :: fname
      character(len=*), intent(in) :: vname
      integer(i4) :: var(:,:,:)

      integer(i4) :: ncid,fh
      integer(i4) :: dim_nx,dim_ny,dim_nb
      integer(i4) :: ierr, mode
      integer(i4) :: nx_block,ny_block,nblocks

      mode = NF90_write
      ierr = sio_open(fname,mode,fh)

      nx_block = SIZE(var,dim=1)
      ny_block = SIZE(var,dim=2)
      nblocks  = SIZE(var,dim=3) 
       
      print *,'write_tiles_int: top of subroutine'
      ierr = nf90_redef(fh)
      ierr = nf90_inq_dimid(fh,'nblocks',dimid=dim_nb)
      if(ierr /= nf90_noerr) then 
         ierr = nf90_def_dim(fh,name='nblocks',len=nblocks,dimid=dim_nb)
         if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
      endif

      ierr = nf90_inq_dimid(fh,'nx_block',dimid=dim_nx)
      if(ierr /= nf90_noerr) then 
          ierr = nf90_def_dim(fh,name='nx_block',len=nx_block,dimid=dim_nx)
          if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
      endif

      ierr = nf90_inq_dimid(fh,'ny_block',dimid=dim_ny)
      if(ierr /= nf90_noerr) then 
          ierr = nf90_def_dim(fh,name='ny_block',len=ny_block,dimid=dim_ny)
          if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
      endif

      ierr = nf90_def_var(fh,name=TRIM(vname),xtype=NF90_INT,dimids=(/dim_nx,dim_ny,dim_nb/),varid=ncid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

      !------------------------
      ! end netCDF define mode
      !------------------------
      ierr = nf90_enddef(fh)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_enddef',__LINE__)

      !--------------------------------
      ! write the variables to the file
      !---------------------------------
      ierr = nf90_put_var(fh,ncid,var)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)

      ierr = sio_close(fh)

    end subroutine write_tiles_int

    subroutine write_tiles_dbl(fname,vname,var)

      character(len=80), intent(in) :: fname
      character(len=*), intent(in) :: vname
      real(r8) :: var(:,:,:)

      integer(i4) :: ncid,fh
      integer(i4) :: dim_nx,dim_ny,dim_nb
      integer(i4) :: ierr, mode
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
         if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
      endif

      ierr = nf90_inq_dimid(fh,'nx_block',dimid=dim_nx)
      if(ierr /= nf90_noerr) then
          ierr = nf90_def_dim(fh,name='nx_block',len=nx_block,dimid=dim_nx)
          if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
      endif

      ierr = nf90_inq_dimid(fh,'ny_block',dimid=dim_ny)
      if(ierr /= nf90_noerr) then
          ierr = nf90_def_dim(fh,name='ny_block',len=ny_block,dimid=dim_ny)
          if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
      endif

      ierr = nf90_def_var(fh,name=TRIM(vname),xtype=NF90_DOUBLE,dimids=(/dim_nx,dim_ny,dim_nb/),varid=ncid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

      !------------------------
      ! end netCDF define mode
      !------------------------
      ierr = nf90_enddef(fh)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_enddef',__LINE__)

      !--------------------------------
      ! write the variables to the file
      !---------------------------------
      ierr = nf90_put_var(fh,ncid,var)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)

      ierr = nf90_close(fh)

    end subroutine write_tiles_dbl

    subroutine write_AppMD(fname,blocks,reorder)

      character(len=80), intent(in) :: fname
      type (AppMD_t) :: blocks(:)
      integer (i4) :: reorder(:)

      integer(i4) :: ierr
      integer(i4) :: mode, fh
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
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
       ierr = nf90_def_dim(fh,name='nx_block',len=nx_block,dimid=dim_nx)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)
       ierr = nf90_def_dim(fh,name='ny_block',len=ny_block,dimid=dim_ny)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_dim',__LINE__)

       !---------------------------------------------------------
       ! define the variables which will be contained in the file
       !---------------------------------------------------------
       ierr = nf90_def_var(fh,name='block_id', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_block_id)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)
       ierr = nf90_def_var(fh,name='ib', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_ib)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)
       ierr = nf90_def_var(fh,name='ie', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_ie)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)
       ierr = nf90_def_var(fh,name='jb', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_jb)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)
       ierr = nf90_def_var(fh,name='je', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_je)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

       ierr = nf90_def_var(fh,name='iblock', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_iblock)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)
       ierr = nf90_def_var(fh,name='jblock', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_jblock)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

       ierr = nf90_def_var(fh,name='npoints', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_npoints)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)
       ierr = nf90_def_var(fh,name='i_glob', xtype = NF90_INT, dimids=(/dim_nx,dim_nb/),varid=ncid_iglob)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)
       ierr = nf90_def_var(fh,name='j_glob', xtype = NF90_INT, dimids=(/dim_ny,dim_nb/),varid=ncid_jglob)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

       ierr = nf90_def_var(fh,name='reorder', xtype = NF90_INT, dimids=(/dim_nb/),varid=ncid_reorder)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_def_var',__LINE__)

       !------------------------
       ! end netCDF define mode
       !------------------------
       ierr = nf90_enddef(fh)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_enddef',__LINE__)
	
       !--------------------------------
       ! write the variables to the file
       !---------------------------------
       ierr = nf90_put_var(fh,ncid_block_id,blocks(:)%block_id)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
       ierr = nf90_put_var(fh,ncid_ib,blocks(:)%ib)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
       ierr = nf90_put_var(fh,ncid_ie,blocks(:)%ie)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
       ierr = nf90_put_var(fh,ncid_jb,blocks(:)%jb)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
       ierr = nf90_put_var(fh,ncid_je,blocks(:)%je)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
       ierr = nf90_put_var(fh,ncid_iblock,blocks(:)%iblock)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
       ierr = nf90_put_var(fh,ncid_jblock,blocks(:)%jblock)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
       ierr = nf90_put_var(fh,ncid_npoints,blocks(:)%npoints)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)

       do i=1,nblocks
	  iglob(:,i) = blocks(i)%i_glob(:)
	  jglob(:,i) = blocks(i)%j_glob(:)
       enddo

       ierr = nf90_put_var(fh,ncid_iglob,iglob)        
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)
       ierr = nf90_put_var(fh,ncid_jglob,jglob)        
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)

       ierr = nf90_put_var(fh,ncid_reorder,reorder)        
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_put_var',__LINE__)

       ierr = nf90_close(fh)

       deallocate(iglob,jglob)
    
    end subroutine write_AppMD

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
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_dimid',__LINE__)
      ierr = nf90_inq_dimid(fh,'nx_block',dimid=dim_nx)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_dimid',__LINE__)
      ierr = nf90_inq_dimid(fh,'ny_block',dimid=dim_ny)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_dimid',__LINE__)

      ierr = nf90_inquire_dimension(fh,dim_nb,len=nblocks)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inquire_dimension',__LINE__)
      ierr = nf90_inquire_dimension(fh,dim_nx,len=nx_block)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inquire_dimension',__LINE__)
      ierr = nf90_inquire_dimension(fh,dim_ny,len=ny_block)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inquire_dimension',__LINE__)

      allocate(blocks(nblocks))

      !--------------------------------------------------------------------
      ! Read and load the different componets of the AppMD_t data structure
      !--------------------------------------------------------------------
	
      ierr = nf90_inq_varid(fh,'block_id',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, blocks(:)%block_id)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      ierr = nf90_inq_varid(fh,'ib',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, blocks(:)%ib)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      ierr = nf90_inq_varid(fh,'ie',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, blocks(:)%ie)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      ierr = nf90_inq_varid(fh,'jb',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, blocks(:)%jb)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      ierr = nf90_inq_varid(fh,'je',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, blocks(:)%je)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      ierr = nf90_inq_varid(fh,'iblock',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, blocks(:)%iblock)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      ierr = nf90_inq_varid(fh,'jblock',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, blocks(:)%jblock)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)
     
      ierr = nf90_inq_varid(fh,'npoints',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, blocks(:)%npoints)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      allocate(iglob(nx_block,nblocks))
      allocate(jglob(ny_block,nblocks))

      ierr = nf90_inq_varid(fh,'i_glob',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, iglob)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      ierr = nf90_inq_varid(fh,'j_glob',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, jglob)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      do i=1,nblocks
	allocate(blocks(i)%i_glob(nx_block))
	allocate(blocks(i)%j_glob(ny_block))
	blocks(i)%i_glob(:) = iglob(:,i)
	blocks(i)%j_glob(:) = jglob(:,i)
      enddo

      deallocate(iglob,jglob)
      allocate(reorder(nblocks))

      ierr = nf90_inq_varid(fh,'reorder',varid)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)
      ierr = nf90_get_var(fh,varid, reorder)
      if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)
      
    end subroutine read_AppMD

    subroutine read_garray_int(fname,vname, var)
      
       character(len=80), intent(in) :: fname
       character(len=*), intent(in) :: vname
       integer(i4), dimension(:,:) :: var

       integer(i4) :: mode
       integer(i4) :: fh
       integer(i4) :: varid
       integer(i4) :: ierr

      print *,'read_garray_int: before sio_open'
       mode = NF90_nowrite
       ierr = sio_open(fname,mode,fh)

      print *,'SIZE(var): ',SIZE(var)
      print *,'read_garray_int: before nf90_inq_varid'
       ierr = nf90_inq_varid(fh,TRIM(vname),varid)
      print *,'read_garray_int: vname is: ',TRIM(vname)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

      print *,'read_garray_int: SIZE(var): ',SIZE(var)
      print *,'read_garray_int: before nf90_get_var'
       ierr = nf90_get_var(fh,varid,var) 
      print *,'read_garray_int: after nf90_get_var'

       if(ierr /= nf90_noerr) &
	    call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

      print *,'read_garray_int: before nf90_get_var'
       ierr = sio_close(fh) 

    end subroutine read_garray_int

    subroutine read_garray_dbl(fname,vname,var)

       character(len=80), intent(in) :: fname
       character(len=*), intent(in) :: vname
       real(r8), dimension(:,:) :: var
       integer(i4) :: mode
       integer(i4) :: fh
       integer(i4) :: varid
       integer(i4) :: ierr

       mode = NF90_nowrite
       ierr = sio_open(fname,mode,fh)

       ierr = nf90_inq_varid(fh,TRIM(vname),varid)
       if(ierr /= nf90_noerr) call check_netcdf(fh,ierr,'nf90_inq_varid',__LINE__)

       ierr = nf90_get_var(fh,varid, var)
       if(ierr /= nf90_noerr) &
            call check_netcdf(fh,ierr,'nf90_get_var',__LINE__)

       ierr = sio_close(fh)

    end subroutine read_garray_dbl

    function sio_create(fname,mode,fh)

       character(len=*), intent(in) :: fname
       integer(i4), intent(in) :: mode
       integer(i4), intent(out) :: fh
       integer(i4) :: sio_create
       integer(i4) :: ierr
       integer(i4) :: fmode
  
       
       ierr = nf90_create(fname,mode,fh)
       ierr = nf90_set_fill(fh,NF90_NOFILL,fmode)

       sio_create = ierr

    end function sio_create

    function sio_open(fname,mode,fh)

       character(len=*), intent(in) :: fname
       integer(i4), intent(in) :: mode
       integer(i4), intent(out) :: fh
       integer(i4) :: sio_open
       integer(i4) :: ierr
       
       ierr = nf90_open(fname,mode,fh)
       if(ierr /= nf90_noerr) then
	    call check_netcdf(fh,ierr,'nf90_open',__LINE__)
        stop
       end if
       sio_open= ierr

    end function sio_open

    function sio_close(ncid)
       integer(i4), intent(inout) :: ncid
       integer(i4) :: sio_close
       integer(i4) :: ierr

       ierr = nf90_sync(ncid)
       ierr = nf90_close(ncid)

       sio_close = ierr

    end function sio_close

    subroutine check_netcdf(fh,status,filestr,line)
         integer(i4) :: fh
	 integer(i4) :: status
	 character(len=*), intent(in) :: filestr
	integer, intent(in) :: line
         print *,'ERROR: ',TRIM(filestr),' ',TRIM(nf90_strerror(status)), ' on line: ',line 
    end subroutine check_netcdf


end module io_serial
