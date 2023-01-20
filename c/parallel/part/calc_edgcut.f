      subroutine CALC_EDGCUT
      use partitioner
!C
!C-- calc. EDGECUT

      IEDGCUT= 0
      do ie= 1, IEDGTOT
        in1= IEDGNOD(ie,1)
        in2= IEDGNOD(ie,2)
        ig1= IGROUP(in1)
        ig2= IGROUP(in2)
        if (ig1.ne.ig2) IEDGCUT= IEDGCUT + 1
      enddo

      write ( *,'(/,"TOTAL EDGE     #   ", i8)') IEDGTOT
      write ( *,'(  "TOTAL EDGE CUT #   ", i8)') IEDGCUT
      write (21,'(/,"TOTAL EDGE     #   ", i8)') IEDGTOT
      write (21,'(  "TOTAL EDGE CUT #   ", i8)') IEDGCUT

      return
      end



