module mpi1s_single_push_halo
    implicit none
    private

    include 'mpif.h'

    integer, parameter, private :: &
        i4 = selected_int_kind(6),  &
        r8 = selected_real_kind(13)

    integer(i4), parameter :: NUMexpectMSG  = 101, &
        expectMSG     = 102, &
        UpdateHaloMSG = 103

    integer(i4) :: my_task,nprocs
    integer(i4), private, allocatable :: Rrequest(:), Srequest(:)
    integer(i4), private, allocatable :: Rstatus(:,:), Sstatus(:,:)

    ! ==========================
    ! schedule
    ! ==========================
    integer(i4) :: COMM
    integer(i4) :: nRecv, nSend
    integer(i4) :: lenRecvBuffer, lenSendBuffer
    integer(i4), allocatable :: sNeigh(:),ptrSend(:),SendCnt(:)
    integer(i4), allocatable :: rNeigh(:),ptrRecv(:),RecvCnt(:)

    integer(i4), allocatable :: place(:)
    integer(i4), allocatable :: rNeighCopies(:, :)
    integer(i4) :: win_nRecv
    integer(i4) :: win_rNeigh, win_ptrRecv

    integer(i4) :: nPrint
    ! ==========================
    ! indirect addressing arrays
    ! ==========================
    integer(i4), allocatable :: halo2send(:)
    integer(i4), allocatable :: recv2halo(:)

    ! ======================
    ! Message buffer arrays
    ! ======================
    real(r8), allocatable :: BufferSendDbl(:)
    real(r8), allocatable :: BufferRecvDblMPI(:)
    integer(i4), allocatable :: BufferRecvIntMPI(:)
    integer(i4) :: winDbl, winInt

    integer(i4), allocatable :: halo2grab(:)
    integer(i4), allocatable :: linear2Proc(:)
    integer(i4) :: linearMaximum, nTotalCA, nActiveCA
    integer(i4) :: maxBufSize, timer_communication

    public :: mpi1s_single_push_update_1d_dbl
    public :: mpi1s_single_push_update_1d_int
    public :: mpi1s_single_push_init

    integer(i4), public :: timer_single_push_halo_init, &
        timer_single_push_update_dbl

  contains

    !***********************************************************************
    !>
    !! This subroutine performs the boundary exchange operation in order to
    !! update the halo elements in the passed array.
    !<
    subroutine mpi1s_single_push_update_1d_dbl(handle,array)
        ! !INPUT PARAMETERS:
        integer(i4), intent(in) :: handle
        real(r8), intent(inout) :: array(:) ! vector to update

        integer(i4) :: src,dest,len,iptr,tag

        integer(i4) :: ierr,i, placeval,j
        integer(kind=MPI_ADDRESS_KIND) :: pos

        integer :: global, origins, targets, origin_array(nRecv)

        !=====================================================================
        ! Indirect address into the 1D data structure to create the SendBuffer
        !=====================================================================
        do i=1,lenSendBuffer
            BufferSendDbl(i) = array(halo2send(i))
        enddo

        origin_array(1:nRecv) = rNeigh(1:nRecv)

        ! Make origins and targets group
        call MPI_COMM_GROUP(MPI_COMM_WORLD, global, ierr)
        call MPI_GROUP_INCL(global, nRecv, rNeigh, origins, ierr)
        call MPI_GROUP_INCL(global, nSend, sNeigh, targets, ierr)

        call MPI_BARRIER(comm, ierr)

        call MPI_WIN_POST(origins, 0, winDbl, ierr)
        call MPI_WIN_START(targets, 0, winDbl, ierr)

        do i=1,nSend
            iptr  = ptrSend(i)
            len   = SendCnt(i)
            dest  = sNeigh(i) + 1
            placeval = place(i)

            pos = int(placeVal-1, MPI_ADDRESS_KIND)
            call MPI_Put(BufferSendDbl(iptr), len, MPI_REAL8, &
            dest-1, pos, len, MPI_REAL8, winDbl, ierr)
        enddo

        call MPI_WIN_COMPLETE(winDbl, ierr)
        call MPI_WIN_WAIT(winDbl, ierr)

        !===============================================================
        ! Indirect address from the RecvBuffer to the 1D data structure
        !===============================================================
        do i=1,lenRecvBuffer
            array(recv2halo(i)) = BufferRecvDblMPI(i)
        enddo
    end subroutine mpi1s_single_push_update_1d_dbl


    !***********************************************************************
    !>
    !! This subroutine performs the boundary exchange operation in order to
    !! update the halo elements in the passed array.
    !<
    subroutine mpi1s_single_push_update_1d_int(handle,array)

        ! !INPUT PARAMETERS:
        integer(i4), intent(in)    :: handle
        integer(i4), intent(inout) :: array(:) ! vector to update


        integer(i4) :: src,dest,len,iptr,tag

        integer(i4), save, allocatable :: BufferSendInt(:)

        integer(i4) :: ierr,i, placeval, j

        logical, parameter :: Debug = .FALSE.

        integer(kind=MPI_ADDRESS_KIND) :: pos

        allocate(BufferSendInt(lenSendBuffer))
        BufferSendInt = 0

        !=====================================================================
        ! Indirect address into the 1D data structure to create the SendBuffer
        !=====================================================================
        do i=1,lenSendBuffer
            BufferSendInt(i) = array(halo2send(i))
        enddo

        call MPI_Win_fence(0, winInt, ierr)
        do i=1,nSend
            iptr  = ptrSend(i)
            len   = SendCnt(i)
            dest  = sNeigh(i) + 1
            placeval = place(i)

            pos = int(placeVal-1, MPI_ADDRESS_KIND)
            call MPI_Put(BufferSendInt(iptr), len, MPI_INTEGER, &
                dest-1, pos, len, MPI_INTEGER, winInt, ierr)
        enddo
        call MPI_Win_fence(0, winInt, ierr)

        !===============================================================
        ! Indirect address from the RecvBuffer to the 1D data structure
        !===============================================================
        do i=1,lenRecvBuffer
            array(recv2halo(i))=BufferRecvIntMPI(i)
        enddo

        deallocate(BufferSendInt)
    end subroutine mpi1s_single_push_update_1d_int

    !***********************************************************************
    !>
    !! This function initilizes scheduling metadata needed for the update
    !! routines
    !! (the update routines perform the boundary exchange operation.)
    !<
    function mpi1s_single_push_init( &
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

        integer(i4) :: j, jmax, neigh, receiverID, neighID
        integer(i4), allocatable :: jmaxCopies(:)

        integer(kind=MPI_ADDRESS_KIND) :: sz, pos

        maxBufSize = maxlinear

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
            call InsertIntoArray(sNeigh,SendCnt,dest,cnt,ierr)
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
        allocate(ptrSend(maxLinear))
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

        ! ===================================
        ! allocate the message buffers
        ! ===================================
        allocate(BufferSendDbl(lenSendBuffer))
        allocate(BufferRecvDblMPI(maxLinear))
        allocate(BufferRecvIntMPI(maxLinear))

        if(Debug) &
            print *,'IAM: ',my_task,' GSHalo_init: finished with subroutine'
        call MPI_Barrier(COMM,ierr)
        handle=1

        ! ==============================================
        ! Initialize windows for one-sided communication
        ! ==============================================
        sz = maxlinear * 4
        call MPI_Win_create(BufferRecvIntMPI, sz, 4, MPI_INFO_NULL, COMM, &
            winInt, ierr)
        sz = maxLinear * 8
        call MPI_Win_create(BufferRecvDblMPI, sz, 8, MPI_INFO_NULL, COMM, &
            winDbl, ierr)

        ! ============================
        ! Place
        ! ============================
        allocate(place(maxLinear))
        call MPI_Barrier(comm,ierr)

        !================================================
        ! Allocate windows for 1-sided MPI communication
        !================================================
        allocate(jmaxCopies(maxNeigh))
        allocate(rNeighCopies(maxNeigh, maxLinear))

        !nRecv
        sz = 4
        call MPI_Win_create(nRecv, sz, 4, MPI_INFO_NULL, COMM, &
        win_nRecv, ierr)

        !rNeigh(:)
        sz = 4 * maxNeigh
        call MPI_Win_create(rNeigh, sz, 4, MPI_INFO_NULL, COMM, &
        win_rNeigh, ierr)

        !ptrRecv(:)
        sz = 4 * maxNeigh
        call MPI_Win_create(ptrRecv, sz, 4, MPI_INFO_NULL, COMM, &
            win_ptrRecv, ierr)

        ! grab jmax value from neighbors
        call MPI_Win_fence(0, win_nRecv, ierr)
        do i=1,nSend
            neigh = sNeigh(i)

            !jmaxCopies(i) = nRecv[neigh]
            pos = int(0, MPI_ADDRESS_KIND)
            call MPI_Get(jmaxCopies(i), 1, MPI_INTEGER, neigh, pos, 1, &
                MPI_INTEGER, win_nRecv, ierr)
        end do
        call MPI_Win_fence(0, win_nRecv, ierr)

        ! grab rNeigh value from neighbors
        call MPI_Win_fence(0, win_rNeigh, ierr)
        do i=1,nSend
            neigh = sNeigh(i)
            jmax = jmaxCopies(i)

            do j=1,jmax
                !rNeighCopies(i,j) = rNeigh(j)[neigh]
                pos = int(j-1, MPI_ADDRESS_KIND)
                call MPI_Get(rNeighCopies(i,j), 1, MPI_INTEGER, neigh, pos, 1, &
                    MPI_INTEGER, win_rNeigh, ierr)
            end do
        end do
        call MPI_Win_fence(0, win_rNeigh, ierr)

        ! calculate push values by iterating through procs to send to
        ! and figuring out where data for this proc is in those proc's buffers
        call MPI_Win_fence(0, win_ptrRecv, ierr)
        do i=1,nSend
            neigh = sNeigh(i)

            receiverID = -1

            jmax = jmaxCopies(i)
            do_j: do j=1,jmax
                neighID = rNeighCopies(i, j)

                if(my_task == neighID) then
                    receiverID = j
                    exit do_j
                end if
            end do do_j

            ! place(i) = ptrRecv(receiverID)[neigh]
            pos = int(receiverID-1, MPI_ADDRESS_KIND)
            call MPI_Get(place(i), 1, MPI_INTEGER, neigh, pos, 1, &
                MPI_INTEGER, win_ptrRecv, ierr)
        end do
        call MPI_Win_fence(0, win_ptrRecv, ierr)
    end function mpi1s_single_push_init

    !***********************************************************************
    subroutine InsertIntoArray(array,cnt, value,cvalue,ierr)
        integer(i4), intent(inout) :: array(:)
        integer(i4), intent(inout) :: cnt(:)
        integer(i4), intent(in)    :: value
        integer(i4), intent(in)    :: cvalue
        integer(i4), intent(out) :: ierr
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
end module mpi1s_single_push_halo
