      program SOLVER33_TEST_SCALAR

      use solver33
      use pfem_util

      implicit REAL*8(A-H,O-Z)
!C
!C +-------+
!C | INIT. |
!C +-------+
!C=== 
      call INPUT_CNTL
      call INPUT_GRID
!C===

!C
!C +---------------------+
!C | matrix connectivity |
!C +---------------------+
!C===
      call MAT_CON0
      call MAT_CON1
!C===
      
!C
!C +-----------------+
!C | MATRIX assemble |
!C +-----------------+
!C===
      call MAT_ASS_MAIN 
      call MAT_ASS_BC
!C===

!C
!C +--------+
!C | SOLVER |
!C +--------+
!C===
      call SOLVE33
!C===

!C
!C +--------+
!C | OUTPUT |
!C +--------+
!C===
      call RECOVER_STRESS
      call OUTPUT_UCD
!C===

      end program SOLVER33_TEST_SCALAR
      

