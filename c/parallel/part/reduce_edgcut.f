      subroutine REDUCE_EDGCUT
      use partitioner

          ISTART_PREV= 0
      IEDGCUTMAX_PREV= 0


      do 700 iter= 1, 100
        write (*,*) iter
        do 600 i= 1, N
          IMASK(i)= 0
 600    continue

!C
!C-- FIND NODE with the MAX. EDGCUT

        do 610 ie= 1, IEDGTOT
          in1= IEDGNOD(ie,1)
          in2= IEDGNOD(ie,2)
          if (IGROUP(in1).ne.IGROUP(in2)) then
            IMASK(in1)= IMASK(in1) + 1
            IMASK(in2)= IMASK(in2) + 1
          endif
 610    continue

        IEDGCUTMAX= IMASK(1)
          ISTART= 1

        ISTACK(1)= 0
        do 620 i= 2, N
          if (IMASK(i).gt.IEDGCUTMAX) then
            if (IDEAD(i).eq.0) then
              IEDGCUTMAX= IMASK(i)
              ISTART    = i
            endif
          endif
          ISTACK(i)= 0
 620    continue

!C
!C-- search MOST FREQUENT GROUP
        icou= 0
        do 630 ie= 1, IEDGTOT
          in1= IEDGNOD(ie,1)
          in2= IEDGNOD(ie,2)

          ig1= IGROUP(in1)
          ig2= IGROUP(in2)

          if (ig1.ne.ig2) then
            if (in1.eq.ISTART) then
                     icou = icou + 1
               IMASK(icou)= in2
              ISTACK(ig2 )= ISTACK(ig2) + 1
            endif       

            if (in2.eq.ISTART .and. ig1.ne.ig2) then
                     icou = icou + 1
               IMASK(icou)= in1
              ISTACK(ig1 )= ISTACK(ig1) + 1
            endif 
          endif 
 630    continue

        ICMAX = -100
        ICOLOR=    0
        do 640 ip= 1, NP
          if (ISTACK(ip).gt.ICMAX) then
             ICMAX= ISTACK(ip)
            ICOLOR= ip
          endif
 640    continue

        if (ISTART_PREV.eq.ISTART .and.                                 &
     &      IEDGCUTMAX_PREV.eq.IEDGCUTMAX) then
          IDEAD(ISTART)= 1
        endif

!C
!C-- RE-GROUP the NODE with MAX. EDGCUT
            ISTART_PREV= ISTART
        IEDGCUTMAX_PREV= IEDGCUTMAX

        IGROUP(ISTART)= ICOLOR

 700  continue

      return
      end
