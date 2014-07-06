!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! Routines to initialize and perform the boundary exchange operation.
!<
module mpi2s_gshalo
    use kinds_mod, only: i4,r8
    implicit none
    private
    
    include 'mpif.h'
    
    integer(i4), parameter :: NUMexpectMSG  = 101, &
        expectMSG     = 102, &
        UpdateHaloMSG = 103
    
    integer(i4) :: my_task,nprocs
    integer(i4), private, allocatable :: Rrequest(:), Srequest(:)
    integer(i4), private, allocatable :: Rstatus(:,:), Sstatus(:,:)
    
    type, private :: Schedule_t
        integer(i4) :: COMM
        integer(i4) :: nRecv,nSend
        integer(i4) :: lenRecvBuffer,lenSendBuffer
        integer(i4), pointer :: sNeigh(:),ptrSend(:),SendCnt(:)
        integer(i4), pointer :: rNeigh(:),ptrRecv(:),RecvCnt(:)
    end type
    
    type (Schedule_t) :: Sched
    
    integer(i4) :: nPrint
    
    ! ==========================
    ! indirect addressing arrays
    ! ==========================
    integer(i4), allocatable :: halo2send(:)
    integer(i4), allocatable :: recv2halo(:)
    
    ! ======================
    ! Message buffer arrays
    ! ======================
    real(r8), allocatable :: BufferRecvDbl(:),BufferSendDbl(:)
    
    public :: mpi2s_gshalo_update_2d_dbl
    public :: mpi2s_gshalo_update_2d_int
    public :: mpi2s_gshalo_init
    
    integer(i4), public :: timer_mpi2s_halo_init, timer_mpi2s_update_dbl

  contains
    
    subroutine mpi2s_gshalo_update_2d_dbl(handle,array)
        integer(i4), intent(in)    :: handle
        real(r8), intent(inout) :: array(:)

        integer(i4) :: lenRecvBuffer,lenSendBuffer
        integer(i4) :: nSend,nRecv

        integer(i4) :: src,dest,len,iptr,tag


        integer(i4) :: ierr,i

        nRecv = Sched%nRecv; lenRecvBuffer=Sched%lenRecvBuffer
        nSend = sched%nSend; lenSendBuffer=Sched%lenSendBuffer

        tag = UpdateHaloMSG
        do i=1,nRecv
            len = Sched%RecvCnt(i)
            iptr = Sched%ptrRecv(i)
            src  = Sched%rNeigh(i)
            call MPI_Irecv(BufferRecvDbl(iptr),len,MPI_REAL8,src, &
                tag, Sched%COMM, Rrequest(i), ierr)
        enddo

        !=====================================================================
        ! Indirect address into the 1D data structure to create the SendBuffer
        !=====================================================================
        do i=1,lenSendBuffer
            BufferSendDbl(i) = array(halo2send(i))
        enddo

        do i=1,nSend
            len  = Sched%SendCnt(i)
            iptr = Sched%ptrSend(i)
            dest  = Sched%sNeigh(i)
            call MPI_Isend(BufferSendDbl(iptr),len,MPI_REAL8,dest, &
                tag, Sched%COMM,Srequest(i), ierr)
        enddo
        if(nSend>0) call MPI_Waitall(nSend,Srequest,Sstatus,ierr)
        if(nRecv>0) call MPI_Waitall(nRecv,Rrequest,Rstatus,ierr)

        !===============================================================
        ! Indirect address from the RecvBuffer to the 1D data structure
        !===============================================================
        do i=1,lenRecvBuffer
            array(recv2halo(i))=BufferRecvDbl(i)
        enddo
    end subroutine mpi2s_gshalo_update_2d_dbl

    subroutine mpi2s_gshalo_update_2d_int(handle,array)
        integer(i4), intent(in)    :: handle
        integer(i4), intent(inout) :: array(:)

        integer(i4) :: lenRecvBuffer,lenSendBuffer
        integer(i4) :: nRecv,nSend

        integer(i4) :: src,dest,len,iptr,tag

        integer(i4), allocatable :: BufferRecvInt(:), BufferSendInt(:)

        integer(i4) :: ierr,i

        logical, parameter :: Debug = .FALSE.


        nRecv = Sched%nRecv; lenRecvBuffer=Sched%lenRecvBuffer
        nSend = Sched%nSend; lenSendBuffer=Sched%lenSendBuffer

        allocate(BufferRecvInt(lenRecvBuffer))
        allocate(BufferSendInt(lenSendBuffer))

        if(Debug) &
            print *,'IAM: ', my_task, &
            ' GSHalo_update_2d_int: halo2send(1:lenSendBuffer) :=', &
            halo2send(1:lenSendBuffer)
        if(Debug) &
            print *,'IAM: ', my_task, &
            ' GSHalo_update_2d_int: recv2halo(1:lenRecvBuffer) :=', &
            recv2halo(1:lenRecvBuffer)

        !=====================================================================
        ! Indirect address into the 1D data structure to create the SendBuffer
        !=====================================================================
        do i=1,lenSendBuffer
            BufferSendInt(i) = array(halo2send(i))
        enddo
        do i=1,nSend
            iptr = Sched%ptrSend(i)
            len  = Sched%SendCnt(i)
            dest = Sched%sNeigh(i)
            if(Debug) &
                print *,'IAM: ',my_task, ' to: ',dest, &
                ' BufferSendInt(:):= ', BufferSendInt(iptr:iptr+len-1)
        enddo

        tag = UpdateHaloMSG
        do i=1,nSend
            len  = Sched%SendCnt(i)
            iptr = Sched%ptrSend(i)
            dest  = Sched%sNeigh(i)
            if(Debug) &
                print *,'IAM: ',my_task,'Posting Send of length ',len, &
                    ' to: ',dest
            call MPI_Isend(BufferSendInt(iptr),len,MPI_INTEGER,dest, &
            tag, Sched%COMM,Srequest(i), ierr)
        enddo

        do i=1,nRecv
            len = Sched%RecvCnt(i)
            iptr = Sched%ptrRecv(i)
            src  = Sched%rNeigh(i)
            if(Debug) &
                print *,'IAM: ',my_task,'Posting Recv of length ',len, &
                    ' from: ',src
            call MPI_IRecv(BufferRecvInt(iptr),len,MPI_INTEGER,src, &
                tag, Sched%COMM, Rrequest(i), ierr)
        enddo
        if(nSend>0) &
            call MPI_Waitall(nSend,Srequest,Sstatus,ierr)
        if(nRecv>0) &
            call MPI_Waitall(nRecv,Rrequest,Rstatus,ierr)

        do i=1,nRecv
            iptr = Sched%ptrRecv(i)
            len  = Sched%RecvCnt(i)
            src  = Sched%rNeigh(i)
            if(Debug) print *,'IAM: ',my_task, ' From: ',src, &
                ' BufferRecvInt(:):= ', BufferRecvInt(iptr:iptr+len-1)
        enddo

        !===============================================================
        ! Indirect address from the RecvBuffer to the 1D data structure
        !===============================================================
        do i=1,lenRecvBuffer
            array(recv2halo(i))=BufferRecvInt(i)
        enddo
        deallocate(BufferRecvInt,BufferSendInt)
    end subroutine mpi2s_gshalo_update_2d_int

    function mpi2s_gshalo_init( &
        COMM,maxlinear, nTotal,nActive,LinearGdof,LinearProc) result(handle)

        integer(i4), intent(in) :: COMM
        integer(i4), intent(in) :: maxlinear
        integer(i4), intent(in) :: nTotal,nActive
        integer(i4), intent(in) :: LinearGdof(maxlinear)
        integer(i4), intent(in) :: LinearProc(maxlinear)

        integer(i4) :: handle

        integer(i4) :: sWords,rWords
        integer(i4) :: sWords_min,sWords_max,sWords_avg,sWords_total

        integer(i4) :: nNeigh
        integer(i4) :: nNeigh_max,nNeigh_min,nNeigh_avg

        integer(i4) :: dest,src,cnt

        integer(i4) :: nSend,nRecv
        integer(i4) :: i,lenHalo,maxNeigh,ierr

        integer(i4), allocatable :: ptrRecv(:)
        integer(i4), allocatable :: RecvCnt(:)

        integer(i4), allocatable :: ptrSend(:)

        integer(i4), allocatable :: rNeigh(:), sNeigh(:)
        integer(i4), allocatable :: SendCnt(:)
        integer(i4), allocatable :: SendTMP(:)

        integer(i4), allocatable :: ptrCnt(:)
        integer(i4) :: lenRecvBuffer,lenSendBuffer
        integer(i4) :: tag
        integer(i4) :: iptr

        logical, parameter :: Debug = .FALSE.
        logical, parameter :: Info = .FALSE.

        integer(i4), allocatable :: tmpBuf(:),rtmpBuf(:)
        integer(i4), allocatable :: sCount(:),sCount2(:)
        integer(i4) :: idx,ig,ip,len
        logical :: found

        integer(i4) :: one
        integer(i4)     :: errorcode,errorlen
        character*(80)  :: errorstring

        ! ==============================
        ! Get some MPI information 
        ! ==============================
        call MPI_COMM_SIZE(COMM,nprocs,ierr)
        call MPI_COMM_RANK(COMM,my_task,ierr)

        lenHalo=nTotal-nActive
        maxNeigh=MIN(lenHalo,nprocs)

        lenSendBuffer=0
        lenRecvBuffer=0

        !=================================
        ! Allocate an array of size nprocs
        !=================================
        allocate(sCount(0:nprocs-1),sCount2(0:nprocs-1))

        allocate(sNeigh(maxNeigh))
        allocate(rNeigh(maxNeigh))
        allocate(RecvCnt(maxNeigh))
        allocate(ptrRecv(maxNeigh))
        rNeigh=-1
        sNeigh=-1
        RecvCnt=0

        nSend=0
        nRecv=0

        ! ==========================
        ! Get a list of my neighbors
        ! ==========================
        one=1
        do i=nActive+1,nTotal
        call InsertIntoArray(rNeigh,RecvCnt,LinearProc(i),one,ierr)
        enddo

        ! ==========================
        ! Count neighbors
        ! ==========================
        nRecv=COUNT(rNeigh .ge. 0)

        sCount=0
        do i=1,nRecv
        sCount(rNeigh(i)) = 1
        enddo
        call MPI_Allreduce(sCount,sCount2,nprocs,MPI_INTEGER,MPI_SUM,COMM,ierr)

        nSend=sCount2(my_task)

        deallocate(sCount,sCount2)
        if(Debug) &
            print *,'Processor: ',my_task,' has R-neighbors: ',rNeigh(1:nRecv)
        rWords = SUM(RecvCnt(1:nRecv))
        if(Info) print *,'IAM: ',my_task,'nNeigh: ',nRecv, &
            ' Bytes received: ', 8*rWords,' : ', 8*RecvCnt(1:nRecv)

        allocate(SendCnt(maxNeigh))
        allocate(SendTMP(maxNeigh))

        SendCnt=0

        if(nRecv>0) then 
            ptrRecv(1)=1
            do i=2,nRecv
                ptrRecv(i) = ptrRecv(i-1)+RecvCnt(i-1)
            enddo
            lenRecvBuffer=ptrRecv(nRecv)+RecvCnt(nRecv)-1
        endif

        if(Debug) &
            print *,'Processor: ',my_task,' Receive pointer: ',ptrRecv(1:nRecv)
        
        ! ===============================================
        ! Pack all the info into the Schedule_t structure
        ! ===============================================
        Sched%nRecv        =  nRecv
        Sched%nSend        =  nSend

        allocate(Sched%ptrRecv(nRecv))
        Sched%ptrRecv(1:nRecv) = ptrRecv(1:nRecv)

        allocate(Sched%RecvCnt(nRecv))
        Sched%RecvCnt(1:nRecv)=RecvCnt(1:nRecv)

        allocate(Sched%rNeigh(nRecv))
        Sched%rNeigh(1:nRecv) = rNeigh(1:nRecv)

        Sched%lenRecvBuffer =  lenRecvBuffer

        allocate(Srequest(nSend),Rrequest(nRecv))
        allocate(Sstatus(MPI_STATUS_SIZE,nSend), &
            Rstatus(MPI_STATUS_SIZE, nRecv))

        !=====================================================
        ! Figure out the number of points to Send to Neighbor
        !=====================================================
        tag = NUMexpectMSG

        do i=1,nSend
            call MPI_Irecv(SendTMP(i),1,MPI_INTEGER,MPI_ANY_SOURCE, &
                tag, COMM, Srequest(i), ierr)
            if(ierr .ne. MPI_SUCCESS) then 
                errorcode=ierr
                call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
                print *,'GSHalo_init: After call to MPI_Irecv: ',errorstring
            endif
        enddo
        
        do i=1,nRecv
            call MPI_Isend(RecvCnt(i),1,MPI_INTEGER,rNeigh(i), &
                tag, COMM,Rrequest(i), ierr)
            if(ierr .ne. MPI_SUCCESS) then 
                errorcode=ierr
                call MPI_Error_String(errorcode,errorstring,errorlen,ierr)
                print *,'GSHalo_init: After call to MPI_Isend: ',errorstring
            endif
        enddo
        
        if(nRecv>0) call MPI_Waitall(nRecv,Rrequest,Rstatus,ierr)
        if(nSend>0) call MPI_Waitall(nSend,Srequest,Sstatus,ierr)

        call MPI_Barrier(comm,ierr)

        do i=1,nSend
            cnt = SendTmp(i)
            dest = Sstatus(MPI_SOURCE,i)
            if(Debug) &
                print *,'IAM: ',my_task,' GSHalo_init: Message source is: ', &
                    dest, ' Length: ',cnt
            call InsertIntoArray(sNeigh,SendCnt,dest,cnt,ierr)
        enddo

        if(Debug) &
            print *,'Processor: ',my_task,' has S-neighbors: ',sNeigh(1:nSend)
        if(nSend>0) then 
            sWords = SUM(SendCnt(1:nSend))
        else
            sWords = 0
        endif
        if(Info) &
            print *,'IAM: ',my_task,'nNeigh: ',nSend,' Bytes sent: ', &
                8*sWords,' : ',8*SendCnt(1:nSend)
        nNeigh = nSend

        if(Info) then 
            call MPI_Allreduce( &
                sWords,sWords_total,1,MPI_INTEGER,MPI_SUM,COMM,ierr)
            call MPI_Allreduce( &
                sWords,sWords_min,1,  MPI_INTEGER,MPI_MIN,COMM,ierr)
            call MPI_Allreduce( &
                sWords,sWords_max,1,  MPI_INTEGER,MPI_MAX,COMM,ierr)
            call MPI_Allreduce( &
                nNeigh,nNeigh_avg,1,  MPI_INTEGER,MPI_SUM,COMM,ierr)
            call MPI_Allreduce( &
                nNeigh,nNeigh_min,1,  MPI_INTEGER,MPI_MIN,COMM,ierr)
            call MPI_Allreduce( &
                nNeigh,nNeigh_max,1,  MPI_INTEGER,MPI_MAX,COMM,ierr)

            if( my_task == 0) then 
                sWords_avg = sWords_total/nprocs
                print *,'gsHalo_init: Bytes {total,min,avg,max}: ', &
                    8*sWords_total,8*sWords_min, &
                    8*sWords_avg,  8*sWords_max 
                nNeigh_avg=nNeigh_avg/nprocs
                print *, 'gsHalo_init: Neighbors {min,avg,max}: ', &
                    nNeigh_min,nNeigh_avg,nNeigh_max
            endif
        endif

        ! ============================================
        ! allocate and fill out the Send pointer array 
        ! ============================================
        allocate(Sched%sNeigh(nSend))
        allocate(ptrSend(nSend))
        if(nSend>0) then
            Sched%sNeigh(1:nSend) = sNeigh(1:nSend)
            ptrSend(1)=1
            do i=2,nSend
                ptrSend(i) = ptrSend(i-1)+SendCnt(i-1)
            enddo
            lenSendBuffer=ptrSend(nSend)+SendCnt(nSend)-1
        endif

        allocate(ptrCnt(nRecv))
        allocate(tmpBuf(lenRecvBuffer))
        allocate(rtmpBuf(lenSendBuffer))

        ptrCnt=0
        do i=nActive+1,nTotal
            ! find index into Neigh array corresponding to processor
            ! LinearProc(i)
            call LinearOrderedFind(rNeigh(1:nRecv),LinearProc(i),found,ip)
            ! global dof
            ig = LinearGdof(i)
            ! calculate index into message buffer
            idx = ptrRecv(ip)+ptrCnt(ip)
            ! Store the dof in the right spot
            tmpBuf(idx)=ig
            ! Increment the the message buffer offset to processor
            ! LinearProc(i)
            ptrCnt(ip)=ptrCnt(ip)+1
        enddo

        nPrint=lenRecvBuffer
        allocate(halo2send(lenSendBuffer))
        allocate(recv2halo(lenRecvBuffer))

        do i=1,lenRecvBuffer
            ig = tmpBuf(i)
            call LinearFind(LinearGdof(nActive+1:nTotal),ig,found,idx)
            recv2halo(i)=nActive+idx
        enddo

        do i=1,nSend
            len = SendCnt(i)
            src = sNeigh(i)
            if(Debug) print *,'IAM: ',my_task,' Posting Recv from: ',src, &
                ' Length: ',len
        enddo
        call MPI_Barrier(COMM,ierr)

        !============================================
        ! Figure out which points to send to Neighbor 
        !============================================
        tag = expectMSG
        do i=1,nSend
            len = SendCnt(i)
            iptr = ptrSend(i)
            src  = sNeigh(i)
            call MPI_Irecv(rtmpBuf(iptr),len,MPI_INTEGER,src, &
                tag, COMM,Srequest(i),ierr)
        enddo

        do i=1,nRecv
            len  = RecvCnt(i)
            iptr = ptrRecv(i)
            dest = rNeigh(i)
            call MPI_Isend(tmpBuf(iptr),len,MPI_INTEGER,dest, &
                tag, COMM,Rrequest(i),ierr)
        enddo

        if(nSend>0) call MPI_Waitall(nSend,Srequest,Sstatus,ierr)
        if(nRecv>0) call MPI_Waitall(nRecv,Rrequest,Rstatus,ierr)

        do i=1,lenSendBuffer
            ig = rtmpBuf(i)
            call LinearOrderedFind(LinearGdof(1:nActive),ig,found,idx)
            halo2send(i) = idx
        enddo

        ! ========================
        ! Assign the send arrays
        ! ========================
        allocate(Sched%SendCnt(nSend))
        Sched%SendCnt(1:nSend)  = SendCnt(1:nSend)

        allocate(Sched%ptrSend(nSend))
        Sched%ptrSend(1:nSend) = ptrSend(1:nSend)

        Sched%COMM          = COMM
        Sched%lenSendBuffer = lenSendBuffer

        ! ===================================
        ! allocate the Real message buffers
        ! ===================================
        allocate(BufferRecvDbl(lenRecvBuffer))
        allocate(BufferSendDbl(lenSendBuffer))
        
        if(Debug) print *,'IAM: ',my_task, &
            ' GSHalo_init: finished with subroutine'
        call MPI_Barrier(COMM,ierr)
        handle=1
    end function mpi2s_gshalo_init

    !***********************************************************************
    !>
    !! This subroutine inserts a unique value into an array in sortered order. 
    !! It also maintains associated satellite data.
    !!
    !! @param array   The array into which a value will be inserted.
    !! @param cnt     Satellite data associated with the array.
    !! @param value   A scalar value for insertion into the array 
    !! @param cvalue  The satellite data for insertion.
    !! @param ierr    Error return code which is not set. [REMOVE]
    !<
    subroutine InsertIntoArray(array,cnt, value,cvalue,ierr)
        integer(i4), intent(inout) :: array(:), cnt(:)
        integer(i4), intent(in)    :: value, cvalue
        integer(i4), intent(out) :: ierr
        integer(i4), allocatable :: tmp(:), ctmp(:)
        logical  :: found
        integer(i4)  :: n,indx
        logical, parameter  :: Debug = .FALSE.

        call LinearOrderedFind(array,value,found,indx)

        if(.not. found) then
            n = SIZE(array)
            allocate(tmp(n),ctmp(n))
            tmp = array
            ctmp  = cnt

            if(indx== 1) then
                array(1)   = value
                cnt(1)     = cvalue
                array(2:n) = tmp(1:n-1)
                cnt(2:n)   = ctmp(1:n-1)
            else
                array(1:indx-1) = tmp(1:indx-1)
                cnt(1:indx-1)   = ctmp(1:indx-1)
                array(indx)     = value
                cnt(indx)       = cvalue
                array(indx+1:n) = tmp(indx:n-1)
                cnt(indx+1:n)   = ctmp(indx:n-1)
            endif

            deallocate(tmp,ctmp)
        else
            cnt(indx) = cnt(indx) + cvalue
        endif
    end subroutine InsertIntoArray

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
    !! This subroutine determines if a scalar value already exists in a sorted
    !! array.
    !!
    !! @param array  A 1D integer array to be searched.
    !! @param value  The scalar value for which to search.
    !! @param found  A logical value which indicates the scalar value is found.
    !! @param indx   If found=.true. the location of the scalar in the array.
    !>
    subroutine LinearOrderedFind(array,value,found,indx)
        integer(i4), intent(in)  :: array(:)
        integer(i4), intent(in)  :: value
        logical, intent(out) :: found
        integer(i4), intent(out) :: indx

        integer(i4) :: n,nz,i
        n = SIZE(array)

        ! =============================================
        ! Array of all zeros... insert at the begining
        ! =============================================
        found = .FALSE.
        if((array(1) < 0 ) .or. (value < array(1)) ) then
            found = .FALSE.
            indx = 1
            return
        endif
        nz = COUNT(array .ge. 0)
        do i=1,nz
            if(array(i) == value) then
                ! ===============
                ! Already in list
                ! ===============
                found = .TRUE.
                indx = i
                return
            else if((array(i) < value) .and. (value < array(i+1)) ) then
                ! =====================================
                ! Insert it into the middle of the list
                ! =====================================
                found = .FALSE.
                indx = i+1
                return
            else if((array(i) < value) .and. (array(i+1) == -1) )then
                ! =====================================
                ! Insert it into the end of the List
                ! =====================================
                found = .FALSE.
                indx = i+1
                return
            endif
        enddo
    end subroutine LinearOrderedFind
end module mpi2s_gshalo
