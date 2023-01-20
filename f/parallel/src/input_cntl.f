!C
!C***
!C*** INPUT_CNTL
!C***
!C
      subroutine INPUT_CNTL
      use pfem_util

      implicit REAL*8 (A-H,O-Z)

      if (my_rank.eq.0) then
        open (11,file= 'INPUT.DAT', status='unknown')
        read (11,'(a80)') HEADER
        read (11,*) METHOD, PRECOND
        read (11,*) iterPREmax
        read (11,*) ITER
        close (11)


        if (METHOD.ne.1.and.METHOD.ne.2.and. METHOD.ne.3) then
          call ERROR_EXIT (102, my_rank)
        endif
        if (PRECOND.ne.0.and.PRECOND.ne.1) then
          call ERROR_EXIT (103, my_rank)
        endif
        if (ITER.le.0) then
          call ERROR_EXIT (104, my_rank)
        endif 
        if (iterPREmax.lt.1) then
          iterPREmax= 1
          call ERROR_EXIT (111, my_rank)
        endif 
        if (iterPREmax.gt.4) then
          iterPREmax= 4
          call ERROR_EXIT (112, my_rank)
        endif 
      endif     

      call MPI_BCAST (HEADER, 80, MPI_CHARACTER, 0, MPI_COMM_WORLD,ierr)

      call MPI_BCAST (METHOD , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (PRECOND, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (iterPREmax, 1,                                    &
     &                            MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST (ITER   , 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


      SIGMA_DIAG= 1.d0
      SIGMA     = 0.d0
      RESID     = 1.d-8
      NSET      = 0

      pfemRarray(1)= RESID
      pfemRarray(2)= SIGMA_DIAG
      pfemRarray(3)= SIGMA

      pfemIarray(1)= ITER
      pfemIarray(2)= METHOD
      pfemIarray(3)= PRECOND
      pfemIarray(4)= NSET
      pfemIarray(5)= iterPREmax

      return
      end
