      subroutine ORDERING (ISTART)
      use partitioner
!C
!C-- ORDERING

      IMASK(ISTART)= 1

      iter= 1

      do
        iter= iter +1
        do   ie= 1, IEDGTOT
          in1= IEDGNOD(ie,1)
          in2= IEDGNOD(ie,2)
          if (IMASK(in1).eq.iter-1 .and. IMASK(in2).eq.0)               &
     &        IMASK(in2)= iter
          if (IMASK(in2).eq.iter-1 .and. IMASK(in1).eq.0)               &
     &        IMASK(in1)= iter
        enddo

        icou= 0
        do i= 1, N
          if (IMASK(i).eq.0) icou= icou + 1
        enddo

        write (*,'(i8," ITERs for ORDERING", i8)') iter, icou
        if (icou.eq.0) exit
      enddo

      ITERC= iter
      write (*,*) '*** ORDERING : MARKED ALL NODES'

      ICOND1(1)     = ISTART
      ICOND2(ISTART)= 1

      icou= 1
      do iter= 2, ITERC
      do i= 1, N
        if (IMASK(i).eq.iter) then
          icou= icou + 1
          ICOND1(icou)= i
          ICOND2(i   )= icou
        endif
      enddo
      enddo

      return
      end
 
