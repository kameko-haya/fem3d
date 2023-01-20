!C
!C***
!C*** KL 
!C***
!C
!C    Graph-Partition Optimization Routine using
!C    Kernighan-Lin Heuristics
!C

      subroutine KL (ip0, ip1, NODPART)
      use partitioner
      integer(kind=kint), dimension (:), allocatable ::  ACTEDG

      type KL_EDGE_LIST
          integer :: GAIN
          integer :: PE
          integer :: GlobalID
          integer :: PREV
      type(KL_EDGE_LIST), pointer :: next
      end type KL_EDGE_LIST

      type(KL_EDGE_LIST), pointer :: root, ptr1, ptr2, ptr1C, ptr2C,    &
     &                                                 ptr1N, ptr2N

      integer :: eof
      integer :: GAIN , PE , GlobalID, PREV
      integer :: GAIN0, PE0, GlobalID0,PREV0
      integer :: GAIN1, PE1, GlobalID1,PREV1
      integer :: GAINT, PET, GlobalIDT,PREVT

      integer :: GAINA(100), PEA(100), GlobalIDA(100)
      integer                                   ACTEDGTOT

      character*2 LINE(80), CHAR(2)
      data CHAR/' 1',' 2'/

      allocate (root)
      allocate (root%next)

!C
!C-- ACTIVE EDGEs

      allocate (NEIBNODTOT(N))
      allocate (NEIBNOD   (N,40))
      allocate (ACTEDG(IEDGTOT))

      ACTEDGTOT= 0
      do ie= 1, IEDGTOT
        ig1= IGROUP(IEDGNOD(ie,1))
        ig2= IGROUP(IEDGNOD(ie,2))
        if ((ig1.eq.ip0 .or. ig1.eq.ip1) .and.                          &
     &      (ig2.eq.ip0 .or. ig2.eq.ip1) ) then
                 ACTEDGTOT = ACTEDGTOT + 1
          ACTEDG(ACTEDGTOT)= ie
        endif
      enddo

!C
!C-- Initial EDGCUT

      call CALC_EDGCUT_KL (IEDGCUT)

      IEDGCUTMIN= IEDGCUT
      IEDGCUTINI= IEDGCUT

      do i= 1, N
        IMASK(i)= IGROUP(i)
      enddo
     
!C
!C-- find 1st NODE

      do i= 1, N
        if (IGROUP(i).eq.ip0 .or. IGROUP(i).eq.ip1) then
          if (IGROUP(i).eq.ip0) PE= 0
          if (IGROUP(i).eq.ip1) PE= 1
          GlobalID= i
          call CALC_GAIN_NODE (i)
          exit
        endif
      enddo

!C
!C-- create LISTs

      root%next%GAIN     = GAIN
      root%next%PE       = PE
      root%next%GlobalID = GlobalID
      root%next%PREV     = 0

      NODEINI            = GlobalID
      icou= 1
      NULLIFY (root%next%next)          

      write (*,*) '*****', NODEINI
      do i= 1, N
        if (IGROUP(i).eq.ip0 .or. IGROUP(i).eq.ip1 .and.                &
     &      i.ne.NODEINI) then
          if (IGROUP(i).eq.ip0) PE= 0
          if (IGROUP(i).eq.ip1) PE= 1
          GlobalID= i
          PREV    = 0
          call CALC_GAIN_NODE (i)
          call INSERT_BOT
        endif
      enddo

!C
!C +-------------+
!C | HANDLE list |
!C +-------------+
!C===

      ITERMAX= NODPART

      do iter= 1, ITERMAX

!C
!C-- PE0 -> PE1      
        ptr1 => root
        ptr2 => ptr1%next
        icou =  0

        do
          icou= icou + 1

          if (.not. associated(ptr2)) exit
          if (ptr2%PE.eq.0 .and. ptr2%PREV.eq.0) then
            GAIN0    = ptr2%GAIN
            PE0      = ptr2%PE
            PREV0    = ptr2%PREV
            GlobalID0= ptr2%GlobalID
            exit
          endif
          ptr1 => ptr2
          ptr2 => ptr1%next
        enddo

        do i= 1, N
          IDEAD(i)= 0
        enddo

        icoutot= 0
        do neib= 1, NEIBNODTOT(GlobalID0)
          inb= NEIBNOD(GlobalID0,neib)
          if (IGROUP(inb).eq.ip0 .or. IGROUP(inb).eq.ip1) then
            icoutot   = icoutot + 1
            IDEAD(inb)= 1
          endif
        enddo

!C
!C-- PE1 -> PE0      
        ptr1 => root
        ptr2 => ptr1%next

        icou= 0
        do
          if (.not. associated(ptr2)) exit
          if (ptr2%PE.eq.1 .and. ptr2%PREV.eq.0) then
            icou= icou + 1
            if (icou.eq.1) then
              GAINT    = ptr2%GAIN
              PET      = ptr2%PE
              PREVT    = ptr2%PREV
              GlobalIDT= ptr2%GlobalID
            endif
            if (IDEAD(ptr2%GlobalID).eq.0) then    
              GAIN1    = ptr2%GAIN
              PE1      = ptr2%PE
              PREV1    = ptr2%PREV
              GlobalID1= ptr2%GlobalID
              exit
            endif
          endif
          ptr1 => ptr2
          ptr2 => ptr1%next
        enddo

        ISWAPFLAG= 0
        if (GAINT-2.gt.GAIN1 .and. PREVT.eq.0 .and. PET.eq.1) then
          GAIN1    = GAINT
          PE1      = PET
          PREV1    = PREVT
          GlobalID1= GlobalIDT
          ISWAPFLAG= 1
        endif

!C
!C-- clear PREV 
        icou= 0
        ptr1 => root
        ptr2 => ptr1%next
        do
          icou= icou + 1
          if (.not. associated(ptr2)) exit
          ptr2%PREV= 0
          if (icou.eq.NODPART) exit
          ptr1 => ptr2
          ptr2 => ptr1%next
        enddo

!C
!C-- NODEs connected to PE0-NODE
        do i= 1, N
          IDEAD(i)= 0
        enddo

        icoutot= 0
        do neib= 1, NEIBNODTOT(GlobalID0)
          inb= NEIBNOD(GlobalID0,neib)
          if (IGROUP(inb).eq.ip0 .or. IGROUP(inb).eq.ip1) then
            icoutot   = icoutot + 1
            IDEAD(inb)= 1
          endif
        enddo

        ptr1C => root
        ptr2C => ptr1C%next
        icou= 0
        do
          if (.not. associated(ptr2C)) exit
          if (IDEAD(ptr2C%GlobalID).eq.1) then
                      icou = icou + 1
            GAINA    (icou)= ptr2C%GAIN
            PEA      (icou)= ptr2C%PE
            GlobalIDA(icou)= ptr2C%GlobalID
          endif
          if (icou.eq.icoutot) exit
          ptr1C => ptr2C
          ptr2C => ptr1C%next
        enddo

        do neib= 1, icoutot
          GAIN    = GAINA    (neib)
          PE      = PEA      (neib)
          PREV    = 0
          GlobalID= GlobalIDA(neib)
          call DELETE
            if (PEA(neib).eq.PE0) GAIN= GAIN + 2
            if (PEA(neib).ne.PE0) GAIN= GAIN - 2
          call INSERT_TOP
        enddo

!C
!C-- NODEs connected to PE1-NODE
        do i= 1, N
          IDEAD(i)= 0
        enddo

        icoutot= 0
        do neib= 1, NEIBNODTOT(GlobalID1)
          inb= NEIBNOD(GlobalID1,neib)
          if (IGROUP(inb).eq.ip0 .or. IGROUP(inb).eq.ip1) then
            icoutot   = icoutot + 1
            IDEAD(inb)= 1
          endif
        enddo

        ptr1C => root
        ptr2C => ptr1C%next
        icou= 0
        do
          if (.not. associated(ptr2C)) exit
          if (IDEAD(ptr2C%GlobalID).eq.1) then
                      icou = icou + 1
            GAINA    (icou)= ptr2C%GAIN
            PEA      (icou)= ptr2C%PE
            GlobalIDA(icou)= ptr2C%GlobalID
          endif
          ptr1C => ptr2C
          ptr2C => ptr1C%next
        enddo

        do neib= 1, icoutot
          GAIN    = GAINA    (neib)
          PE      = PEA      (neib)
          PREV    = 0
          GlobalID= GlobalIDA(neib)

          call DELETE
            if (PEA(neib).eq.PE1) GAIN= GAIN + 2
            if (PEA(neib).ne.PE1) GAIN= GAIN - 2
          call INSERT_TOP
        enddo

!C
!C-- SWAPPED NODEs

        if (ISWAPFLAG.eq.0) then
          PE      =       PE0
          PREV    =     PREV0
          GlobalID= GlobalID0
          GAIN    =     GAIN0

          call DELETE

          PE      = 1
          PREV    = PREV0 + 1
          IGROUP(GlobalID)= ip1
          call CALC_GAIN_NODE (GlobalID)
          call INSERT_BOT
         else
          ptr1 => root
          ptr2 => ptr1%next
          icou =  0
          do
            icou= icou + 1

            if (.not. associated(ptr2)) exit
            if (ptr2%GlobalID.eq.GlobalID0) then
              GAIN0    = ptr2%GAIN
              PE0      = ptr2%PE
              PREV0    = ptr2%PREV
              GlobalID0= ptr2%GlobalID
              exit
            endif
            ptr1 => ptr2
            ptr2 => ptr1%next
          enddo

          PE      = 1
          PREV    = PREV0 + 1
          GlobalID= GlobalID0
          IGROUP(GlobalID)= ip1
          call CALC_GAIN_NODE (GlobalID)
          call INSERT_BOT
        endif

!C
        if (ISWAPFLAG.eq.0) then
          PE      =       PE1
          PREV    =     PREV1
          GlobalID= GlobalID1
          GAIN    =     GAIN1

          call DELETE

          PE      = 0
          PREV    = PREV + 1
          IGROUP(GlobalID)= ip0
          call CALC_GAIN_NODE (GlobalID)
          call INSERT_BOT
         else
          ptr1 => root
          ptr2 => ptr1%next
          icou =  0
          do
            icou= icou + 1
            if (.not. associated(ptr2)) exit
            if (ptr2%GlobalID.eq.GlobalID1) then
              GAIN1    = ptr2%GAIN
              PE1      = ptr2%PE
              PREV1    = ptr2%PREV
              GlobalID1= ptr2%GlobalID
              exit
            endif
            ptr1 => ptr2
            ptr2 => ptr1%next
          enddo

          PE      = 0
          PREV    = PREV1 + 1
          GlobalID= GlobalID0
          IGROUP(GlobalID)= ip0
          call CALC_GAIN_NODE (GlobalID)
          call INSERT_BOT
        endif

!C
!C-- min. EDGCUT

        call CALC_EDGCUT_KL (IEDGCUT)

        if (IEDGCUT.lt.IEDGCUTMIN) then
          do i= 1, N
            IMASK(i)= IGROUP(i)
          enddo
          IEDGCUTMIN= IEDGCUT
        endif

!        write (*,*) iter, IEDGCUTMIN,GlobalID0,GlobalID1
!        read  (*,*) iiii

!          NYP1= 11
!          NXP1=  5
!          do jj= NYP1,1,-1
!            do ii= 1, NXP1
!              LINE(ii)= CHAR(IGROUP((jj-1)*NXP1+ii))
!            enddo
!            write (*,'(5a2)')( LINE(k),k=1,NXP1)
!          enddo
!          read (*,*) iii
           call OUTPUT
!           read (*,*) iii
      enddo

!C===
      do i= 1, N
        IGROUP(i)= IMASK(i)
      enddo
      
      deallocate (root)
      deallocate (ACTEDG)

      return
      contains

!C
!C***
!C***  INSERT_BOT
!C***
!C     insert pointer at the BOTTOM of the stack

        subroutine INSERT_BOT
        ptr1 => root
        ptr2 => ptr1%next

        icou= 0
        do
          if (.not. associated(ptr2)) exit
          if (GAIN .gt. ptr2%GAIN) exit
          if (GAIN .eq. ptr2%GAIN .and.     PE .lt.     ptr2%PE  ) exit

          ptr1 => ptr2
          ptr2 => ptr1%next
        enddo

        allocate (ptr1%next)
        ptr1 => ptr1%next
                ptr1%GAIN     = GAIN
                ptr1%PE       = PE
                ptr1%GlobalID = GlobalID
                ptr1%PREV     = PREV

        ptr1%next => ptr2

        end subroutine INSERT_BOT

!C
!C***
!C***  INSERT_TOP
!C***
!C     insert pointer at the TOP of the stack
 
       subroutine INSERT_TOP
        ptr1 => root
        ptr2 => ptr1%next

        icou= 0
        do
          if (.not. associated(ptr2)) exit
          if (GAIN .gt. ptr2%GAIN) exit
          if (GAIN .eq. ptr2%GAIN .and.     PE .eq.     ptr2%PE  ) exit
          if (GAIN .eq. ptr2%GAIN .and.     PE .lt.     ptr2%PE  ) exit

          ptr1 => ptr2
          ptr2 => ptr1%next
        enddo

        allocate (ptr1%next)
        ptr1 => ptr1%next
                ptr1%GAIN     = GAIN
                ptr1%PE       = PE
                ptr1%GlobalID = GlobalID
                ptr1%PREV     = PREV

        ptr1%next => ptr2

        end subroutine INSERT_TOP
!C
!C***
!C***  DELETE
!C***
        subroutine DELETE
        ptr1 => root
        
        icou= 0
        do
          if (.not. associated(ptr1)) then
            write (*,*) 'NOT FOUND'
            exit
          else if (                                                     &
     &      ptr1%next%GAIN     == GAIN      .and.                       &
     &      ptr1%next%PE       == PE        .and.                       &
     &      ptr1%next%GlobalID == GlobalID ) then

            ptr1%next => ptr1%next%next
            exit
          endif
          ptr1 => ptr1%next
        enddo

        end subroutine DELETE
!C        
!C***
!C***  OUTPUT
!C***
        subroutine OUTPUT
        ptr1 => root%next

        icou0= 0
        icou1= 0
        do 
          if (.not. associated(ptr1)) exit
          if (ptr1%PE.eq.0) icou0= icou0 + 1
          if (ptr1%PE.eq.1) icou1= icou1 + 1
            if ( ((ptr1%PE.eq.0).and.icou0.le.5) .or.                   &
     &           ((ptr1%PE.eq.1).and.icou1.le.5) )                      &
     &       write (*,'(16i5)') ptr1%GAIN, ptr1%PE, ptr1%GlobalID,       &
     &       IGROUP(ptr1%GlobalID),                                     &
     &      (NEIBNOD(ptr1%GlobalID,neib),                               &
     &       neib= 1, NEIBNODTOT(ptr1%GlobalID))

            ptr1 => ptr1%next
            if (icou0.eq.5 .and. icou1.eq.5) exit
        enddo
        end subroutine OUTPUT
!C
!C***
!C***  CALC_EDGCUT_KL
!C***
        subroutine CALC_EDGCUT_KL (IEDGCUT)
        IEDGCUT= 0
        do ie= 1, ACTEDGTOT
          ig1= IGROUP(IEDGNOD(ie,1))
          ig2= IGROUP(IEDGNOD(ie,2))
          if (ig1.eq.ip0 .and. ig2.eq.ip1) IEDGCUT= IEDGCUT+1
          if (ig2.eq.ip0 .and. ig1.eq.ip1) IEDGCUT= IEDGCUT+1
        enddo
        end subroutine CALC_EDGCUT_KL
!C
!C***
!C***  CALC_GAIN_NODE
!C***
!C    calc. GAIN value for each node
 
        subroutine CALC_GAIN_NODE (i)

        GAIN= 0
        do neib= 1, NEIBNODTOT(i)
          in= NEIBNOD(i,neib)
          if (IGROUP(in).eq.ip0 .and. PE.eq.0) GAIN= GAIN - 1
          if (IGROUP(in).eq.ip0 .and. PE.eq.1) GAIN= GAIN + 1
          if (IGROUP(in).eq.ip1 .and. PE.eq.0) GAIN= GAIN + 1
          if (IGROUP(in).eq.ip1 .and. PE.eq.1) GAIN= GAIN - 1
        enddo
        end subroutine CALC_GAIN_NODE

      end subroutine KL



