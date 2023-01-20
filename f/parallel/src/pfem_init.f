!C
!C***
!C*** PFEM_INIT
!C***
!C
!C    INIT. PFEM-FEM process's
!C
      subroutine PFEM_INIT
      use pfem_util
      implicit REAL*8 (A-H,O-Z)

      call MPI_INIT      (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, PETOT, ierr )
      call MPI_COMM_RANK (MPI_COMM_WORLD, my_rank, ierr )
      write(*,*) 'my_rank:', my_rank

      pfemRarray= 0.d0
      pfemIarray= 0

      return
      end
