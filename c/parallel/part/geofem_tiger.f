!------------------------------------------------------------
!     Copyright 1999 by the Research Organization for 
!     Information Science & Technology (RIST)
!------------------------------------------------------------

      program geofem_tiger

      use  geofem_util
      use  analyzer
      use  partitioner

      real(kind=kreal), dimension(:), allocatable :: VAL
      real(kind=kint ), dimension(:), allocatable :: IS1, IS2

      integer(kind=kint )               :: my_rank
      integer(kind=kint )               :: errno

!        call  geofem_init    (my_rank, errno)
        call  init_analyzer  (errno)

!C
!C     N            total node #
!C     NP           total PE (partition) #
!C
!C     NPN   (ip)   totoal node # in each PE
!C     NPC   (ip)   totoal cell # in each PE
!C     NPNID (ii)   1D compressed array for node-to-PE relation
!C     NPCID (ii)   1D compressed array for cell-to-PE relation
!C
!C     RHO   (i)    connected EDGE # to each node
!C     IGROUP(i)    PE ID for each node
!C
!C     IMASK (i)    flag for node partition
!C                  =2 : already belongs to certain PE
!C                  =1 : under operation
!C
!C     ISTACK   (i)    work array
!C
!C     RHOMAX       Max. RHO
!C     RHOMIN       Min. RHO
!C    
!C     ICOND1(NEW)= OLD
!C     ICOND2(OLD)= NEW
!C
!C     NEIBPETOT(ip)    neighboring PE #
!C     NEIBPE   (ip,k)  neighboring PE ID
!C

!C
!C-- PRECONDITIONING

      call PARASET
      if (NTYP.eq.2) then
        call RHOSET   (ISTART)
        call ORDERING (ISTART)
      endif


      if (NTYP.eq.1) then
!C
!C +-----+
!C | RCB |
!C +-----+
!C===
      allocate (VAL(N))
      allocate (IS1(N), IS2(-N:+N))

      do i= 1, N
        IGROUP(i)= 1
      enddo

      do iter= 1, NPOWER

        idir= 300
        do 
          write (*,'(/,"#####",i3,"-th BiSECTION #####")') iter
          write (*,*)                                   ' '
          write (*,'(" in which direction ? X:1, Y:2, Z:3")') 
          write (*,'(/,">>>")')
          read  (*,*) idir
          if (idir.ge.1 .and. idir.le.3) exit
        enddo

        if (idir.eq.1) write (*,'(" X-direction")') 
        if (idir.eq.2) write (*,'(" Y-direction")') 
        if (idir.eq.3) write (*,'(" Z-direction")') 

        do ip0= 1, 2**(iter-1)
          icou= 0
          do i= 1, N
            if (IGROUP(i).eq.ip0) then
                icou= icou + 1
              IS1(icou)= i
              VAL(icou)= XYZ(i,idir)
            endif
          enddo

          call SORT (VAL, IS1, IS2, N, icou)

          do ic= 1, icou/2
            in= IS1(ic)
            IGROUP(in)= ip0 + 2**(iter-1)
          enddo

!          ic0= 0
!          do i= 1, N
!            if (IGROUP(i).eq.ip0) then
!              ic0= ic0 + 1
!              IGROUP(i)= ip0 + 2**(iter-1)
!              if (ic0.eq.icou/2) exit
!            endif
!          enddo

!C
!C-- KL optimization
         ip1= ip0 + 2**(iter-1)

!         call KL (ip0,ip1,icou)

!          NXP1=  5
!          NYP1= 11
!          do jj= NYP1,1,-1
!            do ii= 1, NXP1
!              LINE(ii)= CHAR(IGROUP((jj-1)*NXP1+ii))
!            enddo
!            write (*,'(5a2)')( LINE(k),k=1,NXP1)
!          enddo


        enddo

      enddo

      deallocate (VAL, IS1, IS2)

      endif
!C===

      if (NTYP.eq.2) then
!C
!C +--------+
!C | GREEDY |
!C +--------+
!C===
      NRTOT= 1
      do ip= 1, NP-1
        write (*,*) '*** PROCESSING ', ip, '-th PARTITION' 
        call PARTITION (ip, NRTOT, ISTART)
        if (ip.ne.NP) call PREP_NEXT_PARTITION (ISTART)
      enddo

      do i= 1, N
        if (IGROUP(i).le.0) IGROUP(i)= NP
      enddo
!C===
      endif

!C
!C-- create LOCAL DATA
      call PROC_LOCAL
      call OUTPUT_UCD

      stop ' * normal termination'

 998  continue
        call ERROR_EXIT (21,0)
 999  continue
        call ERROR_EXIT (22,0)

      end program geofem_tiger

      subroutine SORT (STEM, INUM, ISTACK, NP, N)
      real   (kind=8), dimension(NP)      ::  STEM
      integer(kind=4), dimension(NP)      ::  INUM
      integer(kind=4), dimension(-NP:+NP) ::  ISTACK

      M     = 100
      NSTACK= NP

      jstack= 0
      l     = 1
      ir    = N

      ip= 0
 1    continue
      ip= ip + 1

      if (ir-l.lt.M) then
        do 12 j= l+1, ir
          ss= STEM(j)
          ii= INUM(j)

          do 11 i= j-1,1,-1
            if (STEM(i).le.ss) goto 2
            STEM(i+1)= STEM(i)
            INUM(i+1)= INUM(i)
 11       continue
          i= 0

 2        continue
            STEM(i+1)= ss
            INUM(i+1)= ii
 12     continue

        if (jstack.eq.0) return
        ir = ISTACK(jstack)
         l = ISTACK(jstack-1)
        jstack= jstack - 2
       else

        k= (l+ir) / 2
            temp = STEM(k)
        STEM(k)  = STEM(l+1)
        STEM(l+1)= temp

              it = INUM(k)
        INUM(k)  = INUM(l+1)     
        INUM(l+1)= it

        if (STEM(l+1).gt.STEM(ir)) then
              temp = STEM(l+1)
          STEM(l+1)= STEM(ir)
          STEM(ir )= temp
                it = INUM(l+1)
          INUM(l+1)= INUM(ir)
          INUM(ir )= it
        endif

        if (STEM(l).gt.STEM(ir)) then
             temp = STEM(l)
          STEM(l )= STEM(ir)
          STEM(ir)= temp
               it = INUM(l)
          INUM(l )= INUM(ir)
          INUM(ir)= it
        endif

        if (STEM(l+1).gt.STEM(l)) then
              temp = STEM(l+1)
          STEM(l+1)= STEM(l)
          STEM(l  )= temp
                it = INUM(l+1)
          INUM(l+1)= INUM(l)
          INUM(l  )= it
        endif

        i= l + 1
        j= ir

        ss= STEM(l)
        ii= INUM(l)

 3      continue
          i= i + 1
          if (STEM(i).lt.ss) goto 3

 4      continue
          j= j - 1
          if (STEM(j).gt.ss) goto 4     

        if (j.lt.i)        goto 5

        temp   = STEM(i)
        STEM(i)= STEM(j)
        STEM(j)= temp

        it     = INUM(i)
        INUM(i)= INUM(j)
        INUM(j)= it
 
        goto 3
      
 5      continue

        STEM(l)= STEM(j)
        STEM(j)= ss
        INUM(l)= INUM(j)
        INUM(j)= ii

        jstack= jstack + 2

        if (jstack.gt.NSTACK) then
          write (*,*) 'NSTACK overflow'
          stop
        endif

        if (ir-i+1.ge.j-1) then
          ISTACK(jstack  )= ir
          ISTACK(jstack-1)= i
          ir= j-1
         else
          ISTACK(jstack  )= j-1
          ISTACK(jstack-1)= l
          l= i
        endif 

      endif     

      goto 1

      end


      subroutine ERROR_EXIT (IFLAG, nn)
      use  partitioner

      write (*,'(/,a)')                                                 &
     &        "********** MESSAGE from GeoFEM Partitioner **********"
      if (IFLAG.ge. 1001 .and. IFLAG.lt.2000) then
        write (*,'(/,a)')                                               &
     &        " ### ABORT : unexpected ZERO/minus in the orginal file"
        if (IFLAG.eq.1001) write (*,'(  a,/)')                          &
     &        "     TOTAL NODE and/or ELEMENT NUMBER"
        if (IFLAG.eq.1002) write (*,'(  a,i8/)')                        &
     &        "     BOUNDARY GROUP NUMBER (1:node, 2:elem, 3:suf)", nn
        if (IFLAG.eq.1003) write (*,'(  a,i8/)')                        &
     &        "     BOUNDARY info ITEMs   (1:node, 2:elem, 3:suf)", nn
        if (IFLAG.eq.1004) write (*,'(  a,i8/)')                        &
     &        "     ELEMENT type", nn
        if (IFLAG.eq.1005) write (*,'(  a,i8/)')                        &
     &        "     ELEMENT connectivity in ", nn
        stop
      endif

      if (IFLAG.eq. 2001) then
        write (*,'(/,a,i8/)')                                           &
     &        " ### ABORT : local node ID > N appears in ELEMENT", nn
        stop
      endif

      if (IFLAG.eq.2002) then
        write (*,'(/,a  )')                                             &
     &        " ### ABORT : local node ID > N appears in BOUNDARY"
        write (*,'(  a,i8/  )')                                         &
     &        "     (1:node, 2:elem, 3:suf)", nn
        stop
      endif

      if (IFLAG.eq. 2) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : ERROR in ORIGINAL GRID FILE : Parallel Info"
        stop
      endif

      if (IFLAG.eq.11) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : ERROR in ORIGINAL GRID FILE"
        stop
      endif

      if (IFLAG.eq.12) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : UNEXPECTED EOF in ORIGINAL GRID FILE"
        stop
      endif

      if (IFLAG.eq.21) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : ERROR in MeTiS PARTIRION FILE"
        stop
      endif

      if (IFLAG.eq.12) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : UNEXPECTED EOF in MeTiS PARTITION FILE"
        stop
      endif

      if (IFLAG.eq.31) then
        write (*,'(/,a,/)')                                             &
     &        " ### ABORT : MeTiS file INCONSISTENCY"
        stop
      endif

      if (IFLAG.eq.32) then
        write (*,'(/,a,2i8/)')                                          &
     &        " ### ABORT : INVALID PE  and node #", NP, nn
        stop
      endif

      if (IFLAG.eq.33) then
        write (*,'(/,a,i8/)')                                           &
     &        " ### ABORT : INVALID element type", nn
        stop
      endif

      end
