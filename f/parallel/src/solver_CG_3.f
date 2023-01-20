!C
!C*** 
!C*** module solver_CG_3
!C***
!C
      module solver_CG_3
      contains
!C
!C*** CG_3
!C
!C    CG_3 solves the linear system Ax = b with 3*3 block matrix 
!C    using the Conjugate Gradient iterative method with the following
!C    preconditioners for SMP nodes:
!C
      subroutine CG_3                                                   &
     &   (N, NP, NPL, NPU, D, AL, INL, IAL, AU, INU, IAU,               &
     &    B, X, ALU, RESID, ITER, ERROR, my_rank,                       &
     &    NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,                 &
     &                       EXPORT_INDEX, EXPORT_ITEM,                 &
     &                       PRECOND, iterPREmax)     

      use  solver_SR_3

      implicit REAL*8(A-H,O-Z)
      include  'precision.inc'
      include  'mpif.h'

      integer(kind=kint ), intent(in):: N, NP, NPU, NPL, my_rank
      integer(kind=kint ), intent(in):: NEIBPETOT
      integer(kind=kint ), intent(in):: PRECOND

      integer(kind=kint ), intent(inout):: ITER, ERROR
      real   (kind=kreal), intent(inout):: RESID

      real(kind=kreal), dimension(3*NP)   , intent(inout):: B, X
      real(kind=kreal), dimension(9*NPL), intent(inout):: AL
      real(kind=kreal), dimension(9*NPU), intent(inout):: AU
      real(kind=kreal), dimension(9*N  ), intent(inout):: ALU
      real(kind=kreal), dimension(9*NP ), intent(inout):: D

      integer(kind=kint ), dimension(0:NP) ,intent(in) :: INU, INL
      integer(kind=kint ), dimension(  NPL),intent(in) :: IAL
      integer(kind=kint ), dimension(  NPU),intent(in) :: IAU

      integer(kind=kint ), pointer :: NEIBPE(:)
      integer(kind=kint ), pointer :: IMPORT_INDEX(:), IMPORT_ITEM(:)
      integer(kind=kint ), pointer :: EXPORT_INDEX(:), EXPORT_ITEM(:)

      real(kind=kreal), dimension(:),    allocatable       :: WS, WR
      real(kind=kreal), dimension(:,:),  allocatable       :: WW

      integer(kind=kint), parameter ::  R= 1
      integer(kind=kint), parameter ::  Z= 2
      integer(kind=kint), parameter ::  Q= 2
      integer(kind=kint), parameter ::  P= 3
      integer(kind=kint), parameter :: ZP= 4

      integer(kind=kint ) :: MAXIT, iterPREmax

      real   (kind=kreal) :: TOL, W, SS

!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      ERROR= 0

      COMMtime= 0.d0
      COMPtime= 0.d0

      allocate (WW(3*NP,4))
      allocate (WS(3*NP), WR(3*NP))

      MAXIT  = ITER
       TOL   = RESID           

      X = 0.d0
      WW= 0.d0
      WS= 0.d0
      WR= 0.d0
!C===

!C
!C +-----------------------+
!C | {r0}= {b} - [A]{xini} |
!C +-----------------------+
!C===

!C
!C-- INTERFACE data EXCHANGE
      call SOLVER_SEND_RECV_3                                           &
     &   ( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,            &
     &     EXPORT_INDEX, EXPORT_ITEM, WS, WR, X , my_rank)

!C
!C-- BEGIN calculation

      do j= 1, N
           X1= X(3*j-2)
           X2= X(3*j-1)
           X3= X(3*j  )
        WVAL1= B(3*j-2) - D(9*j-8)*X1 - D(9*j-7)*X2 - D(9*j-6)*X3
        WVAL2= B(3*j-1) - D(9*j-5)*X1 - D(9*j-4)*X2 - D(9*j-3)*X3
        WVAL3= B(3*j  ) - D(9*j-2)*X1 - D(9*j-1)*X2 - D(9*j  )*X3

        do k= INL(j-1)+1, INL(j)
          i= IAL(k)
          X1= X(3*i-2)
          X2= X(3*i-1)
          X3= X(3*i  )
          WVAL1= WVAL1 - AL(9*k-8)*X1 - AL(9*k-7)*X2 - AL(9*k-6)*X3
          WVAL2= WVAL2 - AL(9*k-5)*X1 - AL(9*k-4)*X2 - AL(9*k-3)*X3
          WVAL3= WVAL3 - AL(9*k-2)*X1 - AL(9*k-1)*X2 - AL(9*k  )*X3
        enddo

        do k= INU(j-1)+1, INU(j)
          i= IAU(k)
          X1= X(3*i-2)
          X2= X(3*i-1)
          X3= X(3*i  )
          WVAL1= WVAL1 - AU(9*k-8)*X1 - AU(9*k-7)*X2 - AU(9*k-6)*X3
          WVAL2= WVAL2 - AU(9*k-5)*X1 - AU(9*k-4)*X2 - AU(9*k-3)*X3
          WVAL3= WVAL3 - AU(9*k-2)*X1 - AU(9*k-1)*X2 - AU(9*k  )*X3
        enddo

        WW(3*j-2,R)= WVAL1
        WW(3*j-1,R)= WVAL2
        WW(3*j  ,R)= WVAL3
      enddo


      BNRM20= 0.d0
      do i= 1, N
        BNRM20= BNRM20+B(3*i-2)**2+B(3*i-1)**2+B(3*i)**2
      enddo

      call MPI_allREDUCE (BNRM20, BNRM2, 1, MPI_DOUBLE_PRECISION,       &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)

      if (BNRM2.eq.0.d0) BNRM2= 1.d0
      ITER = 0
!C===

      S1_TIME= MPI_WTIME()
      do iter= 1, MAXIT
!      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!C
!C************************************************* Conjugate Gradient Iteration

!C
!C +----------------+
!C | {z}= [Minv]{r} |
!C +----------------+
!C===

      if (PRECOND.eq.0) then
!C
!C== Block SSOR

      do i= 1, N
        WW(3*i-2,ZP)= WW(3*i-2,R)
        WW(3*i-1,ZP)= WW(3*i-1,R)
        WW(3*i  ,ZP)= WW(3*i  ,R)
      enddo

      do i= 1, NP
        WW(3*i-2,Z )= 0.d0
        WW(3*i-1,Z )= 0.d0
        WW(3*i  ,Z )= 0.d0
      enddo

      do iterPRE= 1, iterPREmax
        do i= 1+N, NP
          WW(3*i-2,ZP)= 0.d0
          WW(3*i-1,ZP)= 0.d0
          WW(3*i  ,ZP)= 0.d0
        enddo

!C
!C-- FORWARD

        do i= 1, N
          SW1= WW(3*i-2,ZP)
          SW2= WW(3*i-1,ZP)
          SW3= WW(3*i  ,ZP)
          isL= INL(i-1)+1
          ieL= INL(i)
          do j= isL, ieL
              k= IAL(j)
             X1= WW(3*k-2,ZP)
             X2= WW(3*k-1,ZP)
             X3= WW(3*k  ,ZP)
            SW1= SW1 - AL(9*j-8)*X1 - AL(9*j-7)*X2 - AL(9*j-6)*X3
            SW2= SW2 - AL(9*j-5)*X1 - AL(9*j-4)*X2 - AL(9*j-3)*X3
            SW3= SW3 - AL(9*j-2)*X1 - AL(9*j-1)*X2 - AL(9*j  )*X3
          enddo

          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,ZP)= X1
          WW(3*i-1,ZP)= X2
          WW(3*i  ,ZP)= X3
        enddo

!C
!C-- BACKWARD

        do i= N, 1, -1
          isU= INU(i-1) + 1
          ieU= INU(i) 
          SW1= 0.d0
          SW2= 0.d0
          SW3= 0.d0
          do j= isU, ieU
              k= IAU(j)
             X1= WW(3*k-2,ZP)
             X2= WW(3*k-1,ZP)
             X3= WW(3*k  ,ZP)
            SW1= SW1 + AU(9*j-8)*X1 + AU(9*j-7)*X2 + AU(9*j-6)*X3
            SW2= SW2 + AU(9*j-5)*X1 + AU(9*j-4)*X2 + AU(9*j-3)*X3
            SW3= SW3 + AU(9*j-2)*X1 + AU(9*j-1)*X2 + AU(9*j  )*X3
          enddo
          X1= SW1
          X2= SW2
          X3= SW3
          X2= X2 - ALU(9*i-5)*X1
          X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
          X3= ALU(9*i  )*  X3
          X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
          X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
          WW(3*i-2,ZP)=  WW(3*i-2,ZP) - X1
          WW(3*i-1,ZP)=  WW(3*i-1,ZP) - X2
          WW(3*i  ,ZP)=  WW(3*i  ,ZP) - X3
        enddo

!C
!C-- additive Schwartz

        call SOLVER_SEND_RECV_3                                         &
     &     ( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,          &
     &       EXPORT_INDEX, EXPORT_ITEM, WS, WR, WW(1,ZP), my_rank)

        indexA= ZP
        indexB= Z

      do i= 1, NP
        WW(3*i-2,indexB)= WW(3*i-2,indexB) + WW(3*i-2,indexA)         
        WW(3*i-1,indexB)= WW(3*i-1,indexB) + WW(3*i-1,indexA)         
        WW(3*i  ,indexB)= WW(3*i  ,indexB) + WW(3*i  ,indexA)         
      enddo

      if (iterPRE.eq.iterPREmax) goto 750

      do j= 1, N
         X1= WW(3*j-2,indexB)
         X2= WW(3*j-1,indexB)
         X3= WW(3*j  ,indexB)
        WV1= WW(3*j-2,R) - D(9*j-8)*X1 - D(9*j-7)*X2 - D(9*j-6)*X3
        WV2= WW(3*j-1,R) - D(9*j-5)*X1 - D(9*j-4)*X2 - D(9*j-3)*X3
        WV3= WW(3*j  ,R) - D(9*j-2)*X1 - D(9*j-1)*X2 - D(9*j  )*X3
        do k= INL(j-1)+1, INL(j)
            i= IAL(k)
           X1= WW(3*i-2,indexB)
           X2= WW(3*i-1,indexB)
           X3= WW(3*i  ,indexB)
          WV1= WV1 - AL(9*k-8)*X1 - AL(9*k-7)*X2 - AL(9*k-6)*X3
          WV2= WV2 - AL(9*k-5)*X1 - AL(9*k-4)*X2 - AL(9*k-3)*X3
          WV3= WV3 - AL(9*k-2)*X1 - AL(9*k-1)*X2 - AL(9*k  )*X3
        enddo
        do k= INU(j-1)+1, INU(j)
            i= IAU(k)
           X1= WW(3*i-2,indexB)
           X2= WW(3*i-1,indexB)
           X3= WW(3*i  ,indexB)
          WV1= WV1 - AU(9*k-8)*X1 - AU(9*k-7)*X2 - AU(9*k-6)*X3
          WV2= WV2 - AU(9*k-5)*X1 - AU(9*k-4)*X2 - AU(9*k-3)*X3
          WV3= WV3 - AU(9*k-2)*X1 - AU(9*k-1)*X2 - AU(9*k  )*X3
        enddo
        WW(3*j-2,indexA)= WV1
        WW(3*j-1,indexA)= WV2
        WW(3*j  ,indexA)= WV3
      enddo
 
 750  continue


      enddo
      endif

      if (PRECOND.ne.0) then
!C
!C== Block SCALING

      do i= 1, N
        WW(3*i-2,Z)= WW(3*i-2,R)
        WW(3*i-1,Z)= WW(3*i-1,R)
        WW(3*i  ,Z)= WW(3*i  ,R)
      enddo

      do i= 1, N
        X1= WW(3*i-2,Z)
        X2= WW(3*i-1,Z)
        X3= WW(3*i  ,Z)
        X2= X2 - ALU(9*i-5)*X1
        X3= X3 - ALU(9*i-2)*X1 - ALU(9*i-1)*X2
        X3= ALU(9*i  )*  X3
        X2= ALU(9*i-4)*( X2 - ALU(9*i-3)*X3 )
        X1= ALU(9*i-8)*( X1 - ALU(9*i-6)*X3 - ALU(9*i-7)*X2)
        WW(3*i-2,Z)= X1
        WW(3*i-1,Z)= X2
        WW(3*i  ,Z)= X3
      enddo
      endif
!C===
      
!C
!C +---------------+
!C | {RHO}= {r}{z} |
!C +---------------+
!C===
      RHO0= 0.d0

      do i= 1, N
        RHO0= RHO0 + WW(3*i-2,R)*WW(3*i-2,Z) + WW(3*i-1,R)*WW(3*i-1,Z)  &
     &             + WW(3*i  ,R)*WW(3*i  ,Z)
      enddo

      START_TIME= MPI_WTIME()
      call MPI_allREDUCE (RHO0, RHO, 1, MPI_DOUBLE_PRECISION,           &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
      END_TIME= MPI_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME
!C===

!C
!C +-----------------------------+
!C | {p} = {z} if      ITER=1    |
!C | BETA= RHO / RHO1  otherwise |
!C +-----------------------------+
!C===
      if ( ITER.eq.1 ) then

        do i= 1, N
          WW(3*i-2,P)= WW(3*i-2,Z)
          WW(3*i-1,P)= WW(3*i-1,Z)
          WW(3*i  ,P)= WW(3*i  ,Z)
        enddo
       else
         BETA= RHO / RHO1
         do i= 1, N
           WW(3*i-2,P)= WW(3*i-2,Z) + BETA*WW(3*i-2,P)
           WW(3*i-1,P)= WW(3*i-1,Z) + BETA*WW(3*i-1,P)
           WW(3*i  ,P)= WW(3*i  ,Z) + BETA*WW(3*i  ,P)
         enddo
      endif
!C===

!C
!C +-------------+
!C | {q}= [A]{p} |
!C +-------------+
!C===        

!C
!C-- INTERFACE data EXCHANGE
      START_TIME= MPI_WTIME()
      call SOLVER_SEND_RECV_3                                           &
     &   ( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,            &
     &     EXPORT_INDEX, EXPORT_ITEM, WS, WR, WW(1,P), my_rank)
      END_TIME= MPI_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME

!C
!C-- BEGIN calculation

      do j= 1, N
           X1= WW(3*j-2,P)
           X2= WW(3*j-1,P)
           X3= WW(3*j  ,P)
        WVAL1= D(9*j-8)*X1 + D(9*j-7)*X2 + D(9*j-6)*X3
        WVAL2= D(9*j-5)*X1 + D(9*j-4)*X2 + D(9*j-3)*X3
        WVAL3= D(9*j-2)*X1 + D(9*j-1)*X2 + D(9*j  )*X3
        do k= INL(j-1)+1, INL(j)
           i= IAL(k)
          X1= WW(3*i-2,P)
          X2= WW(3*i-1,P)
          X3= WW(3*i  ,P)
          WVAL1= WVAL1 + AL(9*k-8)*X1 + AL(9*k-7)*X2 + AL(9*k-6)*X3
          WVAL2= WVAL2 + AL(9*k-5)*X1 + AL(9*k-4)*X2 + AL(9*k-3)*X3
          WVAL3= WVAL3 + AL(9*k-2)*X1 + AL(9*k-1)*X2 + AL(9*k  )*X3
        enddo
        do k= INU(j-1)+1, INU(j)
           i= IAU(k)
          X1= WW(3*i-2,P)
          X2= WW(3*i-1,P)
          X3= WW(3*i  ,P)
          WVAL1= WVAL1 + AU(9*k-8)*X1 + AU(9*k-7)*X2 + AU(9*k-6)*X3
          WVAL2= WVAL2 + AU(9*k-5)*X1 + AU(9*k-4)*X2 + AU(9*k-3)*X3
          WVAL3= WVAL3 + AU(9*k-2)*X1 + AU(9*k-1)*X2 + AU(9*k  )*X3
        enddo

        WW(3*j-2,Q)= WVAL1
        WW(3*j-1,Q)= WVAL2
        WW(3*j  ,Q)= WVAL3

      enddo
!C===

!C
!C +---------------------+
!C | ALPHA= RHO / {p}{q} |
!C +---------------------+
!C===
      C10= 0.d0
      do i= 1, N
        C10= C10 + WW(3*i-2,P)*WW(3*i-2,Q) + WW(3*i-1,P)*WW(3*i-1,Q)    &
     &           + WW(3*i  ,P)*WW(3*i  ,Q)
      enddo

      START_TIME= MPI_WTIME()
      call MPI_allREDUCE (C10, C1, 1, MPI_DOUBLE_PRECISION,             &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
      END_TIME= MPI_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME

      ALPHA= RHO / C1
!C===

!C
!C +----------------------+
!C | {x}= {x} + ALPHA*{p} |
!C | {r}= {r} - ALPHA*{q} |
!C +----------------------+
!C===

      do i= 1, N
         X(3*i-2)  = X (3*i-2)   + ALPHA * WW(3*i-2,P)
         X(3*i-1)  = X (3*i-1)   + ALPHA * WW(3*i-1,P)
         X(3*i  )  = X (3*i  )   + ALPHA * WW(3*i  ,P)
        WW(3*i-2,R)= WW(3*i-2,R) - ALPHA * WW(3*i-2,Q)
        WW(3*i-1,R)= WW(3*i-1,R) - ALPHA * WW(3*i-1,Q)
        WW(3*i  ,R)= WW(3*i  ,R) - ALPHA * WW(3*i  ,Q)
      enddo

      DNRM20= 0.d0
      do i= 1, N
        DNRM20= DNRM20 + WW(3*i-2,R)**2 + WW(3*i-1,R)**2                &
     &                                  + WW(3*i  ,R)**2
      enddo
      START_TIME= MPI_WTIME()
      call MPI_allREDUCE (DNRM20, DNRM2, 1, MPI_DOUBLE_PRECISION,       &
     &                    MPI_SUM, MPI_COMM_WORLD, ierr)
      END_TIME= MPI_WTIME()
      COMMtime = COMMtime + END_TIME - START_TIME
      RESID= dsqrt(DNRM2/BNRM2)

!C##### ITERATION HISTORY
        if (my_rank.eq.0 .and. mod(ITER,20).eq.1) then
          write (*, 1000) ITER, RESID
        endif
 1000   format (i5, 1pe16.6)
!C#####

        if ( RESID.le.TOL   ) exit
        if ( ITER .eq.MAXIT ) ERROR= -300

        RHO1 = RHO                                                             
      enddo
!C===

      if (my_rank.eq.0) then
        write (*, 1000) ITER, RESID
      endif
!C
!C-- INTERFACE data EXCHANGE
   30 continue
      E1_TIME= MPI_WTIME()
      COMPtime= E1_TIME - S1_TIME

      R1= 100.d0 * ( 1.d0 - COMMtime/COMPtime )
      if (my_rank.eq.0) then
        write (*,'(/,a,1pe16.6  )') '### elapsed    : ', COMPtime
        write (*,'(  a,1pe16.6  )') '### comm.      : ', COMMtime
        write (*,'(  a,1pe16.6,/)') '### work ratio : ', R1
      endif

      call SOLVER_SEND_RECV_3                                           &
     &   ( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,            &
     &     EXPORT_INDEX, EXPORT_ITEM, WS, WR, X, my_rank)

      deallocate (WW, WR, WS)

      end subroutine        CG_3
      end module     solver_CG_3
