!C
!C***
!C*** MAT_ASS_BC
!C***
!C
      subroutine MAT_ASS_BC
      use pfem_util
      implicit REAL*8 (A-H,O-Z)

      allocate (IWKX(NP,2))
      IWKX= 0

!C
!C== Z=Zmax
      do in= 1, NP
        IWKX(in,1)= 0
      enddo

      ib0= -1
      do ib0= 1, NODGRPtot
        if (NODGRP_NAME(ib0).eq.'Zmax') exit
      enddo

      do ib= NODGRP_INDEX(ib0-1)+1, NODGRP_INDEX(ib0)
        in= NODGRP_ITEM(ib)
        IWKX(in,1)= 1
      enddo

      do in= 1, NP
        if (IWKX(in,1).eq.1) then
          B(3*in-2)= B(3*in-2) - D(9*in-6)*1.d0
          B(3*in-1)= B(3*in-1) - D(9*in-3)*1.d0
          D(9*in-6)= 0.d0
          D(9*in-3)= 0.d0
          D(9*in-2)= 0.d0
          D(9*in-1)= 0.d0
          D(9*in  )= 1.d0
          B(3*in  )= 1.d0

          iS= indexL(in-1) + 1
          iE= indexL(in  )
          do k= iS, iE
            AL(9*k-2)= 0.d0
            AL(9*k-1)= 0.d0
            AL(9*k  )= 0.d0
          enddo

          iS= indexU(in-1) + 1
          iE= indexU(in  )
          do k= iS, iE
            AU(9*k-2)= 0.d0
            AU(9*k-1)= 0.d0
            AU(9*k  )= 0.d0
          enddo
        endif
      enddo

      do in= 1, NP
        iS= indexL(in-1) + 1
        iE= indexL(in  )
        do k= iS, iE
          if (IWKX(itemL(k),1).eq.1) then
            B(3*in-2)=  B(3*in-2) - AL(9*k-6)*1.d0
            B(3*in-1)=  B(3*in-1) - AL(9*k-3)*1.d0
            B(3*in  )=  B(3*in  ) - AL(9*k  )*1.d0
            AL(9*k-6)= 0.d0
            AL(9*k-3)= 0.d0
            AL(9*k  )= 0.d0
          endif
        enddo

        iS= indexU(in-1) + 1
        iE= indexU(in  )
        do k= iS, iE
          if (IWKX(itemU(k),1).eq.1) then
            B(3*in-2)=  B(3*in-2) - AU(9*k-6)*1.d0
            B(3*in-1)=  B(3*in-1) - AU(9*k-3)*1.d0
            B(3*in  )=  B(3*in  ) - AU(9*k  )*1.d0
            AU(9*k-6)= 0.d0
            AU(9*k-3)= 0.d0
            AU(9*k  )= 0.d0
          endif
        enddo
      enddo

!C
!C== Z=Zmin
      do in= 1, NP
        IWKX(in,1)= 0
      enddo

      ib0= -1
      do ib0= 1, NODGRPtot
        if (NODGRP_NAME(ib0).eq.'Zmin') exit
      enddo

      do ib= NODGRP_INDEX(ib0-1)+1, NODGRP_INDEX(ib0)
        in= NODGRP_ITEM(ib)
        IWKX(in,1)= 1
      enddo

      do in= 1, N
        if (IWKX(in,1).eq.1) then
          D(9*in-6)= 0.d0
          D(9*in-3)= 0.d0
          D(9*in-2)= 0.d0
          D(9*in-1)= 0.d0
          D(9*in  )= 1.d0
          B(  3*in)= 0.d0
          
          iS= indexL(in-1) + 1
          iE= indexL(in  )
          do k= iS, iE
            AL(9*k-2)= 0.d0
            AL(9*k-1)= 0.d0
            AL(9*k  )= 0.d0
          enddo

          iS= indexU(in-1) + 1
          iE= indexU(in  )
          do k= iS, iE
            AU(9*k-2)= 0.d0
            AU(9*k-1)= 0.d0
            AU(9*k  )= 0.d0
          enddo
        endif
      enddo

      do in= 1, NP
        iS= indexL(in-1) + 1
        iE= indexL(in  )
        do k= iS, iE
          if (IWKX(itemL(k),1).eq.1) then
            AL(9*k-6)= 0.d0
            AL(9*k-3)= 0.d0
            AL(9*k  )= 0.d0
          endif
        enddo

        iS= indexU(in-1) + 1
        iE= indexU(in  )
        do k= iS, iE
          if (IWKX(itemU(k),1).eq.1) then
            AU(9*k-6)= 0.d0
            AU(9*k-3)= 0.d0
            AU(9*k  )= 0.d0
          endif
        enddo
      enddo

!C
!C== X= Xmin
      do in= 1, NP
        IWKX(in,1)= 0
      enddo

      ib0= -1
      do ib0= 1, NODGRPtot
        if (NODGRP_NAME(ib0).eq.'Xmin') exit
      enddo
      
      do ib= NODGRP_INDEX(ib0-1)+1, NODGRP_INDEX(ib0)
        in= NODGRP_ITEM(ib)
        IWKX(in,1)= 1
      enddo

      do in= 1, N
        if (IWKX(in,1).eq.1) then
          D(9*in-8)= 1.d0
          D(9*in-7)= 0.d0
          D(9*in-6)= 0.d0
          D(9*in-5)= 0.d0
          D(9*in-2)= 0.d0
          B(3*in-2)= 0.d0
          
          iS= indexL(in-1) + 1
          iE= indexL(in  )
          do k= iS, iE
            AL(9*k-8)= 0.d0
            AL(9*k-7)= 0.d0
            AL(9*k-6)= 0.d0
          enddo

          iS= indexU(in-1) + 1
          iE= indexU(in  )
          do k= iS, iE
            AU(9*k-8)= 0.d0
            AU(9*k-7)= 0.d0
            AU(9*k-6)= 0.d0
          enddo
        endif
      enddo

      do in= 1, NP
        iS= indexL(in-1) + 1
        iE= indexL(in  )
        do k= iS, iE
          if (IWKX(itemL(k),1).eq.1) then
            AL(9*k-8)= 0.d0
            AL(9*k-5)= 0.d0
            AL(9*k-2)= 0.d0
          endif
        enddo

        iS= indexU(in-1) + 1
        iE= indexU(in  )
        do k= iS, iE
          if (IWKX(itemU(k),1).eq.1) then
            AU(9*k-8)= 0.d0
            AU(9*k-5)= 0.d0
            AU(9*k-2)= 0.d0
          endif
        enddo
      enddo

!C
!C== Y= Ymin
      do in= 1, NP
        IWKX(in,1)= 0
      enddo
      ib0= -1

      do ib0= 1, NODGRPtot
        if (NODGRP_NAME(ib0).eq.'Ymin') exit
      enddo

      do ib= NODGRP_INDEX(ib0-1)+1, NODGRP_INDEX(ib0)
        in= NODGRP_ITEM(ib)
        IWKX(in,1)= 1
      enddo

      do in= 1, N
        if (IWKX(in,1).eq.1) then
          D(9*in-7)= 0.d0
          D(9*in-4)= 1.d0
          D(9*in-1)= 0.d0
          D(9*in-5)= 0.d0
          D(9*in-3)= 0.d0
          B(3*in-1)= 0.d0
          
          iS= indexL(in-1) + 1
          iE= indexL(in  )
          do k= iS, iE
            AL(9*k-5)= 0.d0
            AL(9*k-4)= 0.d0
            AL(9*k-3)= 0.d0
          enddo

          iS= indexU(in-1) + 1
          iE= indexU(in  )
          do k= iS, iE
            AU(9*k-5)= 0.d0
            AU(9*k-4)= 0.d0
            AU(9*k-3)= 0.d0
          enddo
        endif
      enddo

      do in= 1, NP
        iS= indexL(in-1) + 1
        iE= indexL(in  )
        do k= iS, iE
          if (IWKX(itemL(k),1).eq.1) then
            AL(9*k-7)= 0.d0
            AL(9*k-4)= 0.d0
            AL(9*k-1)= 0.d0
          endif
        enddo

        iS= indexU(in-1) + 1
        iE= indexU(in  )
        do k= iS, iE
          if (IWKX(itemU(k),1).eq.1) then
            AU(9*k-7)= 0.d0
            AU(9*k-4)= 0.d0
            AU(9*k-1)= 0.d0
          endif
        enddo
      enddo

      return
      end
