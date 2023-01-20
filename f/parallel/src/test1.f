      program SOLVER33_TEST_SCALAR

      use solver33
      use pfem_util

      implicit REAL*8(A-H,O-Z)
!C
!C +-------+
!C | INIT. |
!C +-------+
!C=== 
      call PFEM_INIT
      call INPUT_CNTL
      call INPUT_GRID
!C===

!C
!C +---------------------+
!C | matrix connectivity |
!C +---------------------+
!C===
      START_TIME= MPI_WTIME()

      call MAT_CON0
      call MAT_CON1

      END_TIME= MPI_WTIME()
      if (my_rank.eq.0) then
        write (*, '("*** matrix conn. ", 1pe16.6, " sec.")')            &
     &         END_TIME-START_TIME
      endif
!C===
      
!C
!C +-----------------+
!C | MATRIX assemble |
!C +-----------------+
!C===
      START_TIME= MPI_WTIME()

      call MAT_ASS_MAIN
      call MAT_ASS_BC

      END_TIME= MPI_WTIME()
      if (my_rank.eq.0) then
        write (*, '("*** matrix ass.  ", 1pe16.6, " sec.",/)')          &
     &         END_TIME-START_TIME
      endif
!C===

!C
!C +--------+
!C | SOLVER |
!C +--------+
!C===
      START_TIME= MPI_WTIME()

      call SOLVE33

      END_TIME= MPI_WTIME()

      i= N
      write (*,'(i8,3e16.6)') i,X(3*i-2),X(3*i-1),X(3*i)


      if (my_rank.eq.0) then
        write (*, '("*** real  COMP.  ", 1pe16.6, " sec.",/)')          &
     &         END_TIME-START_TIME
      endif
!C===

!C
!C +--------+
!C | OUTPUT |
!C +--------+
!C===
      call RECOVER_STRESS
      call OUTPUT_UCD
!C===

      call PFEM_FINALIZE

      end program SOLVER33_TEST_SCALAR

!C
!C***
!C*** ERROR_EXIT
!C***
!C
      subroutine ERROR_EXIT (IFLAG, my_rank)

      if (my_rank.eq.0) then

      if (IFLAG.eq. 101) then
        write (*,'(a)') '#### PFEM-SOL-E0101: PEsmpTOT must be >0'
        call MPI_FINALIZE (errno)
        stop
      endif

      if (IFLAG.eq. 102) then
        write (*,'(a)') '#### PFEM-SOL-E0102: METHOD must be 1 or 2'
        call MPI_FINALIZE (errno)
        stop
      endif

      if (IFLAG.eq. 103) then
        write (*,'(a)') '#### PFEM-SOL-E0103: PRECOND must be 0 or 1'
        call MPI_FINALIZE (errno)
        stop
      endif

      if (IFLAG.eq. 104) then
        write (*,'(a)') '#### PFEM-SOL-E0104: ITER must be >0'
        call MPI_FINALIZE (errno)
        stop
      endif

      if (IFLAG.eq. 111) then
        write (*,'(a)') '#### PFEM-SOL-W0111: iterPREmax must be >=1'
        write (*,'(a)') '                     iterPREmax set to   =1'
        return
      endif

      if (IFLAG.eq. 112) then
        write (*,'(a)') '#### PFEM-SOL-W0112: iterPREmax must be =<4'
        write (*,'(a)') '                     iterPREmax set to   =4'
        return
      endif

      endif

      if (IFLAG.eq. 201) then
        write (*,'(a)') '#### PFEM-SOL-E0201: too many PEs specified'
        write (*,'(a)') '                       in MPIRUN.'
        call MPI_ABORT (errno)
        stop
      endif

      if (IFLAG.eq. 202) then
        write (*,'(a,i8)') '#### PFEM-SOL-E0202: invalid mesh data',    &
     &                                           my_rank
        call MPI_FINALIZE (errno)
        stop
      endif

      return
      end
      

