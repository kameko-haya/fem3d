      subroutine NEIB_PE
      use partitioner

      do ip= 1, NP
        NEIBPETOT(ip)= 0
      do is= ISTACKC(ip-1)+1, ISTACKC(ip)
        icel= NPCID(is)
      do  k= 1, NODELM(icel)
        in= ICELNOD(icel,k)
        ig= IGROUP (in)

        if (ig.ne.ip) call FIND_NEIBPE (ig, ip)
      enddo
      enddo
      enddo
      
      write ( *,'(/," PE/NEIB-PE#    NEIB-PEs")')
      write (21,'(/," PE/NEIB-PE#    NEIB-PEs")')
      do ip= 1, NP
        write ( *,'(i3,i4,5x, 31i4)') ip-1, NEIBPETOT(ip),              &
     &           (NEIBPE(ip,k)-1,k=1,NEIBPETOT(ip))
        write (21,'(i3,i4,5x, 31i4)') ip-1, NEIBPETOT(ip),              &
     &           (NEIBPE(ip,k)-1,k=1,NEIBPETOT(ip))
      enddo

      return
      end




      subroutine FIND_NEIBPE (ig, ip)
      use partitioner

      do inei= 1, NEIBPETOT(ip)
        if (ig.eq.NEIBPE(ip,inei)) return
      enddo

                NEIBPETOT(ip) = NEIBPETOT(ip) + 1
      NEIBPE(ip,NEIBPETOT(ip))= ig

      return
      end
