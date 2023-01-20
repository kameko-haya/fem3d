!C
!C***
!C*** MAT_CON1
!C***
!C
      subroutine MAT_CON1
      use pfem_util
      implicit REAL*8 (A-H,O-Z)

      allocate (indexL(0:N), indexU(0:N))
      indexL= 0
      indexU= 0

      do i= 1, N
        indexL(i)= indexL(i-1) + INL(i)
        indexU(i)= indexU(i-1) + INU(i)
      enddo

      NPL= indexL(N)
      NPU= indexU(N)

      allocate (itemL(NPL), itemU(NPU))

      do i= 1, N
        do k= 1, INL(i)
                kk = k + indexL(i-1)
          itemL(kk)=        IAL(i,k)
        enddo
        do k= 1, INU(i)
                kk = k + indexU(i-1)
          itemU(kk)=        IAU(i,k)
        enddo
      enddo

      deallocate (INL, INU, IAL, IAU)

      end subroutine MAT_CON1
