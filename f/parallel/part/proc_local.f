      subroutine PROC_LOCAL
      use partitioner

      open (21,file='partition.log',status='unknown')

      if (NTYP.eq.1) then
        write ( *,'(/,"RECURSIVE COORDINATE BISECTION")')
        write (21,'(/,"RECURSIVE COORDINATE BISECTION")')
      endif

      if (NTYP.eq.2) then
        write ( *,'(/,"GREEDY")')
        write (21,'(/,"GREEDY")')
      endif

      if (NTYP.eq.3) then
        write ( *,'(/,"MeTiS")')
        write (21,'(/,"MeTiS")')
      endif

      write  ( *,'(/,"*** GRID  file   ", a80)')  GRIDFIL
      write  (21,'(/,"*** GRID  file   ", a80)')  GRIDFIL

      if (NTYP.eq.3) then
        write  ( *,'("*** MeTiS file   ", a80)')  METISFIL
        write  (21,'("*** MeTiS file   ", a80)')  METISFIL
      endif

      write ( *,'(/,i5, " PEs")') NP
      write (21,'(/,i5, " PEs")') NP

!C
!C-- create LOCAL DATA
      call CALC_EDGCUT
      deallocate (IEDGNOD)
      call CRE_LOCAL_DATA

!C
!C-- OVERLAPPED ELEMENTs

      do icel= 1, IELMTOT
        ISTACK(icel)= 0
      enddo

      do icel= 1, IELMTOT
        do k1= 1, NODELM(icel)
        do k2= 1, NODELM(icel)
          ig1= IGROUP(ICELNOD(icel,k1))
          ig2= IGROUP(ICELNOD(icel,k2))
          if (ig1.ne.ig2) ISTACK(icel)= 1
        enddo
        enddo
      enddo

      icou= 0
      do icel= 1, IELMTOT
        if (ISTACK(icel).eq.1) icou= icou + 1
      enddo
      write ( *,'(/,"OVERLAPPED ELEMENTS", i8)')  icou
      write (21,'(/,"OVERLAPPED ELEMENTS", i8)')  icou

!C
!C-- NEIGHBORING PEs
      call NEIB_PE

!C
!C-- INTERFACE info.
      call INTERFACE_NODES
      close (21)

!C
!C-- distributed Local DATA

      call DOUBLE_NUMBERING
      call LOCAL_DATA

      return
      end



