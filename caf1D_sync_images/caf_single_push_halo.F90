!==============================================================================
! Copyright (C) 2010, University Corporation for Atmospheric Research,
!                     Colorado State University,
!                     Los Alamos National Security, LLC,
!                     United States Department of Energy
!
! All rights reserved.  See ../COPYING for copyright details
!==============================================================================

!>
!! Routines to initialize and perform the boundary exchange operation by
!! buffering and pushing data.
!<
module caf_single_push_halo
    use kinds_mod, only: i4, r8
    use, intrinsic :: iso_fortran_env, only: lock_type

    implicit none
    private

    include 'mpif.h'

    ! ========================================
    ! MPI Communication Information
    ! ========================================

    ! Tags that are used in MPI communication (in initialization routine):
    integer(i4), parameter :: &
        NUMexpectMSG  = 101, &
        expectMSG     = 102, &
        UpdateHaloMSG = 103

    integer(i4) :: &
        my_task, &       ! Proccess ID for the local task
        nprocs           ! Total number of proccesses

    integer(i4), private, allocatable :: &
        Rrequest(:), &   ! Stores MPI request objects for non-blocking receive
                     &   ! calls
        Srequest(:)      ! Stores MPI request objects for non-blocking send
                         ! calls
    integer(i4), private, allocatable :: &
        Rstatus(:,:), &  ! Stores MPI status objects for non-blocking receive
                      &  ! calls
        Sstatus(:,:)     ! Stores MPI status objects for non-blocking send
                         ! calls

    integer(i4) :: COMM  ! MPI communicator

    ! ==========================
    ! schedule
    ! ==========================
    integer(i4) :: &
        nRecv[*], &         ! Number of messages to locally recieve
        nSend[*]            ! Number of messages to send from this proccess

    integer(i4) :: &
        lenRecvBuffer,  &   ! How large the local receive buffer needs to be
        maxnRecv,       &
        maxRecvBufSize, &
        recvAmt,        &
        sendAmt

    integer(i4) :: &
        lenSendBuffer[*]    ! How large the local send buffer needs to be

    integer(i4), allocatable :: &
        sNeigh(:)[:],  &    ! List of neighbors  to send data to (by MPI rank)
        ptrSend(:)[:], &    ! Array of offsets in send buffer where data to be
                       &    ! sent to each neighbor in sNeigh should be placed.
        SendCnt(:)[:]       ! Number of values that should be received from
                            ! each neighbor in rNeigh

    integer(i4), allocatable :: &
        rNeigh(:)[:],  &    ! List of neighbors  to receive data from (by MPI
                       &    ! rank)
        ptrRecv(:)[:], &    ! Array of offsets in receive buffer where data
                       &    ! received from each neighbor in rNeigh should be
                       &    ! placed.
        RecvCnt(:)[:]       ! Number of values that should be received from
                            ! each neighbor in rNeigh

    integer(i4), allocatable :: &
        place(:)[:]         ! For each push, identifies where to place data in
                            ! recieve buffer

    integer(i4), allocatable :: &
        halo2grab(:)        ! Identifies where remote data is located for each
                            ! point in the halo that will be pulled.

    integer(i4), allocatable :: &
        linear2Proc(:)      ! Identifies which proccess owns the data that is

    logical, allocatable :: &
        sent(:)

    integer(i4), allocatable :: &
        IDOnRecvSide(:)[:]

    integer(i4), allocatable :: &
        allNeigh(:)

    integer(i4) :: nPrint

    ! ==========================
    ! indirect addressing arrays
    ! ==========================
    integer(i4), allocatable :: &
        halo2send(:)        ! Maps points in halo to send buffer

    integer(i4), allocatable :: &
        recv2halo(:)        ! Maps points in receive buffer to halo

    ! ======================
    ! Message buffer arrays
    ! ======================
    real(r8), allocatable ::  &
        BufferRecvDbl(:,:)[:], &  ! Receive buffer (co-array since values will
                               &  ! be pushed onto this)
        BufferSendDbl(:)          ! Send buffer

    type(lock_type), private :: recvLocks(100)[*]

    ! ==================
    ! array information
    ! ==================
    integer(i4) :: &
        linearMaximum, &    ! Identifies the size of the largest array that
                       &    ! could be passed to the update routine.  (This
                       &    ! value will be globally replicated across all
                       &    ! procs).
        nTotalCA,      &    ! Identifies the size of the array that will be
                       &    ! passed to the update routine.  (This value may
                       &    ! be unique on all procs).
        nActiveCA           ! Identifies how many points in the passed array
                            ! are active ocean points (non halo points).

    ! ==========================
    ! Public subroutines:
    ! ==========================
    public :: caf_single_push_update_1d_dbl_start
    public :: caf_single_push_update_1d_dbl_finish
    public :: caf_single_push_update_1d_int
    public :: caf_single_push_init

  contains

    !***********************************************************************
    !>
    !! This subroutine performs the boundary exchange operation in order to
    !! update the halo elements in the passed array.
    !!
    !! @param handle     Communication handle (generated by initialization
    !!                   routine).
    !! @param array      vector to update
    !<
    subroutine caf_single_push_update_1d_dbl_start(handle,array)
        ! !INPUT PARAMETERS:
        integer(i4), intent(in) :: handle
        real(r8), intent(inout) :: array(:)

        integer(i4) :: src,dest,len,iptr,tag

        integer(i4) :: ierr,i,lockID,placeval,j, numSent
        logical :: gotL


        !=====================================================================
        ! Indirect address into the 1D data structure to create the SendBuffer
        !=====================================================================
        do i=1,lenSendBuffer
            BufferSendDbl(i) = array(halo2send(i))
        enddo

        sent = .FALSE.
        numSent = 0

        sync images(allNeigh)
        do i=1,sendAmt
            dest  = sNeigh(i) + 1
            iptr  = ptrSend(i)
            len   = SendCnt(i)
            placeval = place(i)

            BufferRecvDbl(1:len, IDOnRecvSide(i))[dest] = &
            BufferSendDbl(iptr:iptr+len-1)
        enddo
    end subroutine caf_single_push_update_1d_dbl_start

    !***********************************************************************
    !>
    !! This subroutine performs the boundary exchange operation in order to
    !! update the halo elements in the passed array.
    !!
    !! @param handle     Communication handle (generated by initialization
    !!                   routine).
    !! @param array      vector to update
    !<
    subroutine caf_single_push_update_1d_dbl_finish(handle,array)
        ! !INPUT PARAMETERS:
        integer(i4), intent(in) :: handle
        real(r8), intent(inout) :: array(:)

        integer(i4) :: i, j, ierr, dest
        logical :: gotL

        sync images(allNeigh)

        do i=1,recvAmt
            do j=1,RecvCnt(i)
                array(recv2halo(j + ptrRecv(i)-1)) = &
                BufferRecvDbl(j, i)
            enddo
        enddo
    end subroutine caf_single_push_update_1d_dbl_finish

    !***********************************************************************
    !>
    !! This subroutine performs the boundary exchange operation in order to
    !! update the halo elements in the passed array.
    !!
    !! @param handle     Communication handle (generated by initialization
    !!                   routine).
    !! @param array      vector to update
    !<
    subroutine caf_single_push_update_1d_int(handle,array)
        ! !INPUT PARAMETERS:
        integer(i4), intent(in)    :: handle
        integer(i4), intent(inout) :: array(:) ! vector to update

        integer(i4) :: src,dest,len,iptr,tag

        integer(i4), save, allocatable :: BufferSendInt(:)
        integer(i4), save, allocatable :: BufferRecvIntCAF(:)[:]

        integer(i4) :: ierr,i, placeval, j

        logical, parameter :: Debug = .FALSE.

        allocate(BufferRecvIntCAF(linearMaximum)[*])
        allocate(BufferSendInt(lenSendBuffer))

        sync all

        BufferSendInt = 0

        BufferRecvIntCAF = 0

        !=====================================================================
        ! Indirect address into the 1D data structure to create the SendBuffer
        !=====================================================================
        do i=1,lenSendBuffer
            BufferSendInt(i) = array(halo2send(i))
        enddo

        do i=1,nSend
            iptr  = ptrSend(i)
            len   = SendCnt(i)
            dest  = sNeigh(i) + 1
            placeval = place(i)

            BufferRecvIntCAF(placeval:placeval+len-1)[dest] = &
            BufferSendInt(iptr:iptr+len-1)
        enddo

        sync all

        !===============================================================
        ! Indirect address from the RecvBuffer to the 1D data structure
        !===============================================================
        do i=1,lenRecvBuffer
            array(recv2halo(i))=BufferRecvIntCAF(i)
        enddo

        deallocate(BufferSendInt)
        deallocate(BufferRecvIntCAF)
    end subroutine caf_single_push_update_1d_int


    !***********************************************************************
    !>
    !! This function initilizes scheduling metadata needed for the update
    !! routines
    !! (the update routines perform the boundary exchange operation.)
    !!
    !! @param COMM          MPI communicator
    !! @param maxlinear     Maximum size of vectors
    !! @param nTotal        Total elements in local vector
    !! @param nActive       Number of active (non-halo) elements
    !! @param LinearGdof    Maps vector indices to GDOFs
    !! @param LinearProc    Maps halo indices to processors
    !<
    function caf_single_push_init( &
        COMMin,maxlinear, nTotal,nActive,LinearGdof,LinearProc) result(handle)

        ! !INPUT PARAMETERS:
        integer(i4), intent(in) :: COMMin
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

        integer(i4) :: i,lenHalo,maxNeigh,ierr

        integer(i4), allocatable :: SendTMP(:)

        integer(i4), allocatable :: ptrCnt(:)
        integer(i4) :: tag
        integer(i4) :: iptr

        logical, parameter :: Debug = .FALSE.
        logical, parameter :: Info = .FALSE.

        integer(i4), allocatable :: tmpBuf(:),rtmpBuf(:)
        integer(i4), allocatable :: sCount(:),sCount2(:)
        integer(i4) :: idx,ig,ip,len
        logical :: found

        integer(i4)     :: one
        integer(i4)     :: errorcode,errorlen
        character*(80)  :: errorstring


        integer(i4) :: j, neigh, recieverID

        linearMaximum = maxlinear

        ! ==============================
        ! Get some MPI information 
        ! ==============================
        COMM = COMMin
        call MPI_COMM_SIZE(COMM,nprocs,ierr)
        call MPI_COMM_RANK(COMM,my_task,ierr)

        lenHalo=nTotal-nActive
        maxNeigh=nprocs

        lenSendBuffer=0
        lenRecvBuffer=0

        !=================================
        ! Allocate an array of size nprocs
        !=================================
        allocate(sCount(0:nprocs-1),sCount2(0:nprocs-1))

        allocate(sNeigh(maxNeigh)[*])
        allocate(rNeigh(maxNeigh)[*])
        allocate(RecvCnt(maxNeigh)[*])
        allocate(ptrRecv(maxNeigh)[*])
        allocate(IDOnRecvSide(maxNeigh)[*])

        rNeigh=-1
        sNeigh=-1
        RecvCnt=0
        IDOnRecvSide=-1
        nSend=0
        nRecv=0

        ! ==========================
        ! Get a list of my neighbors
        ! ==========================
        one=1
        do i=nActive+1,nTotal
            call InsertIntoArray(rNeigh,RecvCnt,LinearProc(i),one)
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
        allocate(SendCnt(maxNeigh)[*])
        allocate(SendTMP(maxNeigh))

        SendCnt=0

        if(nRecv>0) then 
            ptrRecv(1)=1
            do i=2,nRecv
                ptrRecv(i) = ptrRecv(i-1)+RecvCnt(i-1)
            enddo
            lenRecvBuffer=ptrRecv(nRecv)+RecvCnt(nRecv)-1
        endif

        ! ===============================================
        ! Pack all the info into the Schedule_t structure
        ! ===============================================
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
            if(Debug) print *,'IAM: ',my_task, &
            'GSHalo_init: Message source is: ',dest, ' Length: ',cnt
            call InsertIntoArray(sNeigh,SendCnt,dest,cnt)
        enddo

        if(Debug) &
            print *,'Processor: ',my_task,' has S-neighbors: ',sNeigh(1:nSend)
        if(nSend>0) then 
            sWords = SUM(SendCnt(1:nSend))
        else
            sWords = 0
        endif
        if(Info) print *,'IAM: ',my_task,'nNeigh: ',nSend,' Bytes sent: ', &
            8*sWords,' : ',8*SendCnt(1:nSend)
        nNeigh = nSend

        ! ============================================
        ! allocate and fill out the Send pointer array 
        ! ============================================
        allocate(ptrSend(maxLinear)[*])
        if(nSend>0) then
            sNeigh(1:nSend) = sNeigh(1:nSend)
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
            ! Increment the the message buffer offset to processor LinearProc(i)
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
        enddo
        call MPI_Barrier(COMM,ierr)

        !============================================
        ! Figure out which points to send to Neighbor 
        !============================================
        tag = expectMSG
        do i=1,nSend
            len  = SendCnt(i)
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
        COMM          = COMMin

        call MPI_AllReduce(lenSendBuffer, maxRecvBufSize, 1, MPI_INTEGER, &
            MPI_MAX, COMM)
        call MPI_AllReduce(nRecv, maxnRecv, 1, MPI_INTEGER, &
            MPI_MAX, COMM)

        ! ===================================
        ! allocate the Real message buffers
        ! ===================================
        allocate(BufferRecvDbl(maxRecvBufSize, maxnRecv)[*])
        allocate(BufferSendDbl(lenSendBuffer))

        if(Debug) &
            print *,'IAM: ',my_task,' GSHalo_init: finished with subroutine'
        call MPI_Barrier(COMM,ierr)
        handle=1

        ! ====================================
        ! Create the list of neighbors
        ! ====================================
        allocate(sent(nSend))
        allocate(allNeigh(nRecv + nSend + 1))
        do i=1,nRecv
            allNeigh(i) = rNeigh(i) + 1
        end do
        do i=1,nSend
            allNeigh(i + nRecv) = sNeigh(i) + 1
        end do
        allNeigh(nSend + nRecv + 1) = this_image()

        ! ============================
        ! Place
        ! ============================
        sync all

        allocate(place(maxLinear)[*])

        sync all

        do i=1,nSend
            neigh = sNeigh(i)+1

            recieverID = -1
            do_j: do j=1,nRecv[neigh]
                if(this_image() == rNeigh(j)[neigh]+1) then
                    recieverID = j
                    exit do_j
                end if
            end do do_j

            IDOnRecvSide(i) = recieverID

            place(i) = ptrRecv(recieverID)[neigh]

            if(place(i) /= ptrRecv(recieverID)[neigh]) then
                print *, "Error on image", this_image(), "neighbor", i
            end if
        end do

        sync all

        recvAmt=nRecv
        sendAmt=nSend
    end function caf_single_push_init

    !***********************************************************************
    !>
    !! This subroutine inserts a unique value into an array in sortered order.
    !! It also maintains associated satellite data.
    !!
    !! @param array   The array into which a value will be inserted.
    !! @param cnt     Satellite data associated with the array.
    !! @param value   A scalar value for insertion into the array
    !! @param cvalue  The satellite data for insertion.
    !<
    subroutine InsertIntoArray(array,cnt, value,cvalue)
        integer(i4), intent(inout) :: array(:)
        integer(i4), intent(inout) :: cnt(:)
        integer(i4), intent(in)    :: value
        integer(i4), intent(in)    :: cvalue
        integer(i4), allocatable :: tmp(:)
        integer, allocatable :: ctmp(:)
        logical  :: found
        integer(i4)  :: indx
        integer(i4)  :: n
        logical, parameter  :: Debug = .FALSE.

        call LinearOrderedFind(array,value,found,indx)

        if(Debug) print *,'IAM: ',my_task,' InsertIntoArray: array is :',array
        if(Debug) print *,'IAM: ',my_task,' InsertIntoArray: value is :',value
        if(Debug) print *,'IAM: ',my_task,' InsertIntoArray: found is :',found
        if(Debug) print *,'IAM: ',my_task,' InsertIntoArray: indx is :',indx

        if(.not. found) then
            n = SIZE(array)
            allocate(tmp(n))
            allocate(ctmp(n))
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

            deallocate(tmp)
            deallocate(ctmp)
        else
            cnt(indx) = cnt(indx) + cvalue
        endif
    end subroutine InsertIntoArray

    !***********************************************************************
    !>
    !! Find the index in an array where a value is stored.
    !!
    !! @param array   Array to search through.
    !! @param value   Value to search for.
    !! @param found   On output set to TRUE if value was found in array.
    !! @param indx    Index of first instance of value in the array.
    !<
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
    !! Assuming an array is sorted by numeric value, find the index where a
    !! value is stored.
    !!
    !! @param array   Array to search through.
    !! @param value   Value to search for.
    !! @param found   On output set to TRUE if value was found in array.
    !! @param indx    Index of first instance of value in the array.
    !<
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
end module caf_single_push_halo
