      subroutine PREP_NEXT_PARTITION (ISTART)
      use partitioner
!C
!C-- find START POINT for NEXT PARTITION

        do i= 1, N
          if (IMASK (i).ne.2) IMASK (i)=  0
          if (IGROUP(i).eq.0) IMASK (i)=  0
                              ISTACK(i)= -1
        enddo

        do ie= 1, IEDGTOT
          in1= IEDGNOD(ie,1)
          in2= IEDGNOD(ie,2)

          IM1= IMASK(in1)
          IM2= IMASK(in2)

          if (IM1.eq.2 .and. IM2.ne.2) IMASK(in2)= 1
          if (IM2.eq.2 .and. IM1.ne.2) IMASK(in1)= 1
        enddo

        do i= 1, N
          if (IMASK(i).eq.1 .and. IGROUP(i).eq.0) then
            ISTACK(ICOND2(i))= RHO(i)
          endif
        enddo

        OUTER:                                                          &
     &  do ir= 0, RHOMAX
        do ic= 1, N
            if (ISTACK(ic).eq.ir) then
              is = ICOND1(ic) 
              exit OUTER
            endif
        enddo
        enddo OUTER

        ISTART= is

      return
      end


