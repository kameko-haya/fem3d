      subroutine PARTITION (ip, NRTOT, ISTART)
      use partitioner
      integer(kind=kint), dimension(:), allocatable :: IW1

!C
!C-- define. START POINT

        allocate (IW1(N))
        IW1= 0

        IMASK (ISTART)= 3
        IGROUP(ISTART)= ip
        NPN_CUR       = 1
        NPN_PREV      = 1
        IFLAG         = 0

        do icel= 1, IELMTOT
          IWORK(icel)= 0
        enddo

        do i= 1, N
          if (IMASK (i).le.2) IMASK (i)= 0
          if (IGROUP(i).eq.0) IMASK (i)= 0
        enddo

!C
!C-- EXTEND CONNECTED REGION
        

        OUTER:                                                          &
     &  do iter= 1, 2*NPN(ip)

!C
!C-- init.
          if (iter.ne.1 .and. NPN_CUR.eq.NPN_PREV) then
            icou= 0
            do ie= 1, IEDGTOT
              in1= IEDGNOD(ie,1)
              in2= IEDGNOD(ie,2)
              IG1= IGROUP(in1)
              IG2= IGROUP(in2)

              if (IG1.eq.ip .and. IG2.eq.0) then
                IMASK(in2)= 1
                icou= icou + 1
              endif
              if (IG2.eq.ip .and. IG1.eq.0) then
                IMASK(in1)= 1
                icou= icou + 1
              endif
            enddo
            if (icou.eq.0) then
              do ie= 1, IEDGTOT
                in1= IEDGNOD(ie,1)
                in2= IEDGNOD(ie,2)
                IG1= IGROUP(in1)
                IG2= IGROUP(in2)

                if (IG1.ne.0 .and. IG2.eq.0) then
                  IMASK(in2)= 1
                  icou= icou + 1
                endif
                if (IG2.ne.0 .and. IG1.eq.0) then
                  IMASK(in1)= 1
                  icou= icou + 1
                endif
              enddo
            endif
          endif

          NPN_PREV= NPN_CUR

        if (NPN_CUR.ge.NPN(ip)) exit OUTER
        if (NRTOT  .ge.  N    ) exit OUTER

!C
!C-- "MARK" nodes which connected to the INTERFACE nodes
          do icel= 1, IELMTOT
          do    k= 1, NODELM(icel)
              in= ICELNOD(icel,k)
              if (IMASK(in).ge.3) IWORK(icel)= 1
            enddo
          enddo

          do icel= 1, IELMTOT
            iflag= IWORK(icel)
            do  k= 1, NODELM(icel)
              in= ICELNOD(icel,k)
              IM= IMASK  (in)
              if (iflag.eq.1 .and. IM.lt.2) IMASK(in)= 1
            enddo
          enddo

!C
!C-- count. marked nodes
          icou= 0
          do i= 1, N
            if (IMASK(i).eq.1 .and. IGROUP(i).eq.0) then
              icou= icou + 1
              ISTACK(icou)=  i
              IW1   (icou)=  RHO   (i)
            endif
          enddo

          IACTTOT= icou

!C
!C-- add marked nodes to the partition
!C   priority is set to nodes 
!C
!C     1. with smaller RHO
!C     2. with smaller node-ID (RENUMBERED) for same RHO
!C

          do ir= 0, RHOMAX
          do ic= 1, IACTTOT
            if (IW1(ic).eq.ir) then
                  is = ISTACK(ic)
              NPN_CUR= NPN_CUR + 1
              NRTOT  = NRTOT   + 1
              IGROUP(is)= ip
              IMASK (is)=  3
              if (NPN_CUR.eq.NPN(ip)) exit OUTER
              if (NRTOT  .eq.N      ) exit OUTER
            endif
          enddo
          enddo

!C
!C-- ACTIVE EDGEs

          do iae= 1, IACTEDGTOT
             ie= IACTEDG(iae)
            in1= IEDGNOD(ie,1)
            in2= IEDGNOD(ie,2)
            IM1= IMASK(in1)
            IM2= IMASK(in2)
            if (IM1.ge.2) then
              RHO     (in2)= RHO(in2)-1
              IEDGFLAG(ie) = 1
            endif
            if (IM2.ge.2) then
              RHO     (in1)= RHO(in1)-1
              IEDGFLAG(ie) = 1
            endif
          enddo

          do i= 1, N
            if (RHO(i).le.0 .and. IGROUP(i).eq.0) then
              NPN_CUR= NPN_CUR + 1
              NRTOT  = NRTOT   + 1
              IGROUP(i)= ip
              IMASK (i)= 3
              if (NPN_CUR.eq.NPN(ip)) exit OUTER
            endif
            if (IMASK (i).le.2) IMASK (i)= 0
            if (IGROUP(i).eq.0) IMASK (i)= 0
          enddo

        if (NPN_CUR.ge.NPN(ip)) exit OUTER         
        if (NRTOT  .ge.  N    ) exit OUTER               

        deallocate (IW1)
        allocate   (IW1(IEDGTOT))
        icou= 0
        do iae= 1, IACTEDGTOT
          ie= IACTEDG(iae)
          if (IEDGFLAG(ie).eq.0) then
                icou = icou + 1
            IW1(icou)= ie
          endif
        enddo

        IACTEDGTOT= icou
        do iae= 1, IACTEDGTOT
          IACTEDG(iae)= IW1(iae)
        enddo

!        write (*,'(5i10)')  ip, iter, NPN_CUR, NPN_PREV, IACTEDGTOT

        enddo OUTER

!        write (*,'(5i10)')  ip, iter, NPN_CUR, NPN_PREV, IACTEDGTOT

        do i= 1, N
          if (IMASK(i).ge.3) IMASK(i)= 2
        enddo

        deallocate (IW1)

      return
      end
