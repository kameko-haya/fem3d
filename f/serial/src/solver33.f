      module SOLVER33

      contains
        subroutine SOLVE33

        use pfem_util
        use solver_CG_3

        implicit REAL*8 (A-H,O-Z)

        real(kind=kreal), dimension(3,3) :: ALU
        real(kind=kreal), dimension(3)   :: PW  

        integer :: ERROR, ICFLAG
        character(len=char_length) :: BUF

        data ICFLAG/0/

!C
!C +------------+
!C | PARAMETERs |
!C +------------+
!C===
        ITER      = pfemIarray(1)
        METHOD    = pfemIarray(2)
        PRECOND   = pfemIarray(3)
        NSET      = pfemIarray(4)
        iterPREmax= pfemIarray(5)
        NREST     = pfemIarray(6)

        RESID     = pfemRarray(1)
        SIGMA_DIAG= pfemRarray(2)

        if (iterPREmax.lt.1) iterPREmax= 1
        if (iterPREmax.gt.4) iterPREmax= 4
!C===

!C
!C +-----------+
!C | BLOCK LUs |
!C +-----------+
!C===
        if (ICFLAG.eq.0) then
          allocate (ALUG(9*N))
          ICFLAG= 1

          BUF(1:30)= '### LINEAR SOLVER:  3x3 Block' 

          if (METHOD.eq.1) BUF(31:40)= 'ssCG' 
          if (METHOD.eq.2) BUF(31:40)= 'ssBiCGSTAB' 

          if (PRECOND.eq.0) then
            if (iterPREmax.eq.1) BUF(41:64)= 'BILU(0)-no ASDD' 
            if (iterPREmax.eq.2) BUF(41:64)= 'BILU(0)-1x ASDD' 
            if (iterPREmax.eq.3) BUF(41:64)= 'BILU(0)-2x ASDD' 
            if (iterPREmax.eq.4) BUF(41:64)= 'BILU(0)-3x ASDD' 
          endif

          if (PRECOND.ne.0) BUF(41:64)= 'Block Scaling' 

          write (*,'(a64)') BUF
        endif

        if (NSET.eq.0) then
          ALUG  = 0.d0

          do ii= 1, N
            ALU(1,1)= D(9*ii-8)*SIGMA_DIAG
            ALU(1,2)= D(9*ii-7)
            ALU(1,3)= D(9*ii-6)
            ALU(2,1)= D(9*ii-5)
            ALU(2,2)= D(9*ii-4)*SIGMA_DIAG
            ALU(2,3)= D(9*ii-3)
            ALU(3,1)= D(9*ii-2)
            ALU(3,2)= D(9*ii-1)
            ALU(3,3)= D(9*ii  )*SIGMA_DIAG
            do k= 1, 3
               L = k
              ALO= dabs(ALU(L,k))
              do i= k+1, 3
                if (dabs(ALU(i,k)).gt.ALO) then
                   L = i
                  ALO= dabs(ALU(L,k))
                endif
              enddo

              ALU(k,k)= 1.d0/ALU(k,k)
              do i= k+1, 3
                ALU(i,k)= ALU(i,k) * ALU(k,k)
                do j= k+1, 3
                  PW(j)= ALU(i,j) - ALU(i,k)*ALU(k,j)
                enddo
                do j= k+1, 3
                  ALU(i,j)= PW(j)
                enddo
              enddo
            enddo
            ALUG(9*ii-8)= ALU(1,1)
            ALUG(9*ii-7)= ALU(1,2)
            ALUG(9*ii-6)= ALU(1,3)
            ALUG(9*ii-5)= ALU(2,1)
            ALUG(9*ii-4)= ALU(2,2)
            ALUG(9*ii-3)= ALU(2,3)
            ALUG(9*ii-2)= ALU(3,1)
            ALUG(9*ii-1)= ALU(3,2)
            ALUG(9*ii  )= ALU(3,3)
          enddo
        endif
!C===

!C
!C +------------------+
!C | ITERATIVE solver |
!C +------------------+
!C===
      if (METHOD.eq.1) then       
        call CG_3                                                       &
     &     ( N, NPL, NPU, D, AL, indexL, itemL, AU, indexU, itemU,      &
     &       B, X,  ALUG, RESID, ITER, ERROR,                           &
     &       PRECOND, iterPREmax) 
      endif

      ITERactual= ITER
!C===

      end subroutine SOLVE33
      end module SOLVER33
