!C
!C*** 
!C*** module solver_SR_3
!C***
!C
      module solver_SR_3
      contains
!C
!C*** SOLVER_SEND_RECV
!C
      subroutine  SOLVER_SEND_RECV_3                                    &
     &                ( N, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,&
     &                                        EXPORT_INDEX, EXPORT_ITEM,&
     &                  WS, WR, X, my_rank)

      implicit REAL*8 (A-H,O-Z)
      include  'mpif.h'
      include  'precision.inc'

      integer(kind=kint )                , intent(in)   ::  N
      integer(kind=kint )                , intent(in)   ::  NEIBPETOT
      integer(kind=kint ), pointer :: NEIBPE      (:)
      integer(kind=kint ), pointer :: IMPORT_INDEX(:)
      integer(kind=kint ), pointer :: IMPORT_ITEM  (:)
      integer(kind=kint ), pointer :: EXPORT_INDEX(:)
      integer(kind=kint ), pointer :: EXPORT_ITEM  (:)
      real   (kind=kreal), dimension(3*N), intent(inout):: WS
      real   (kind=kreal), dimension(3*N), intent(inout):: WR
      real   (kind=kreal), dimension(3*N), intent(inout):: X
      integer                            , intent(in)   :: my_rank

      integer(kind=kint ), dimension(:,:), save, allocatable :: sta1
      integer(kind=kint ), dimension(:,:), save, allocatable :: sta2
      integer(kind=kint ), dimension(:  ), save, allocatable :: req1
      integer(kind=kint ), dimension(:  ), save, allocatable :: req2  

      integer(kind=kint ), save :: NFLAG
      data NFLAG/0/

!C
!C-- INIT.
      if (NFLAG.eq.0) then
        allocate (sta1(MPI_STATUS_SIZE,NEIBPETOT))
        allocate (sta2(MPI_STATUS_SIZE,NEIBPETOT))
        allocate (req1(NEIBPETOT))
        allocate (req2(NEIBPETOT))
        NFLAG= 1
      endif
       
!C
!C-- SEND
      do neib= 1, NEIBPETOT
        istart= EXPORT_INDEX(neib-1)
        inum  = EXPORT_INDEX(neib  ) - istart
        do k= istart+1, istart+inum
               ii   = 3*EXPORT_ITEM(k)
           WS(3*k-2)= X(ii-2)
           WS(3*k-1)= X(ii-1)
           WS(3*k  )= X(ii  )
        enddo

        call MPI_ISEND (WS(3*istart+1), 3*inum,MPI_DOUBLE_PRECISION,    &
     &                  NEIBPE(neib), 0, MPI_COMM_WORLD, req1(neib),    &
     &                  ierr)
      enddo

!C
!C-- RECEIVE
      do neib= 1, NEIBPETOT
        istart= IMPORT_INDEX(neib-1)
        inum  = IMPORT_INDEX(neib  ) - istart
        call MPI_IRECV (WR(3*istart+1), 3*inum, MPI_DOUBLE_PRECISION,   &
     &                  NEIBPE(neib), 0, MPI_COMM_WORLD, req2(neib),    &
     &                  ierr)
      enddo

      call MPI_WAITALL (NEIBPETOT, req2, sta2, ierr)
   
      do neib= 1, NEIBPETOT
        istart= IMPORT_INDEX(neib-1)
        inum  = IMPORT_INDEX(neib  ) - istart
        do k= istart+1, istart+inum
            ii   = 3*IMPORT_ITEM(k)
          X(ii-2)= WR(3*k-2)
          X(ii-1)= WR(3*k-1)
          X(ii  )= WR(3*k  )
        enddo
      enddo

      call MPI_WAITALL (NEIBPETOT, req1, sta1, ierr)

      end subroutine solver_send_recv_3
      end module     solver_SR_3
