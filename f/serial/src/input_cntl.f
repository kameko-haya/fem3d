!C
!C***
!C*** INPUT_CNTL
!C***
!C
      subroutine INPUT_CNTL
      use pfem_util

      implicit REAL*8 (A-H,O-Z)

        open (11,file= 'INPUT.DAT', status='unknown')
        read (11,'(a80)') fname
        read (11,*) METHOD, PRECOND
        read (11,*) iterPREmax
        read (11,*) ITER
        read (11,*) ELAST, POISSON
        close (11)

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
