      subroutine RHOSET (ISTART)
      use partitioner
!C
!C-- calc. RHO

      write (*,*) '*** calc. WEIGHTING FACTOR'
      do ie= 1, IEDGTOT
        in1= IEDGNOD(ie,1)
        in2= IEDGNOD(ie,2)

        RHO(in1)= RHO(in1) + 1
        RHO(in2)= RHO(in2) + 1
      enddo

!C
!C-- find START POINT.

      write (*,*) '*** FIND START POINT'

      RHOMAX= RHO(1)
      RHOMIN= RHO(1)

      ISTART  = 1
      do i= 2, N
        IRANK= RHO(i)
        if (IRANK.lt.RHOMIN) then
          RHOMIN= RHO(i)
          ISTART= i       
        endif
        RHOMAX= max (RHOMAX,IRANK)
      enddo

      write (*,*) ISTART

      return
      end





