!------------------------------------------------------------
!     Copyright 1999 by the Research Organization for 
!     Information Science & Technology (RIST)
!------------------------------------------------------------

      subroutine LOCAL_DATA
      use partitioner

      character*80 LINE
      integer(kind=kint), dimension(:)  , allocatable :: IWKM
      integer(kind=kint), dimension(:)  , allocatable :: IMASK1, IMASK2
      integer(kind=kint), dimension(:)  , allocatable :: IELMLOCAL


!C
!C-- init.
      allocate (IWKM  (NP))

      N2= max (IELMTOT, N)
      allocate (IMASK1(N2))
      allocate (IMASK2(N2))

      allocate (IELMLOCAL(IELMTOT))

      INODTOTG= N       
      IELMTOTG= IELMTOT

      do ip= 1, NP

!C
!C +--------------------------+
!C | read INITIAL LOCAL files |
!C +--------------------------+
!C===
      open (11,file=WORKFIL(ip),status='unknown',form='unformatted')
        rewind (11)
          read (11) ID
          read (11) INODTOT
          read (11) N2
          read (11) N3         
          read (11) (IWKM(i),i=1,N3)

          do i= 1, INODTOT
            read (11) iip, INODLOCAL(i)
          enddo

          NN1= STACK_IMPORT (N3,ip)
          NN2= STACK_EXPORT (N3,ip)
          if (NN1.gt.INODTOT .or. NN2.gt.INODTOT) then
            write (*,'(a,i5)') '### too many communication in PE ', ip
          endif

          STACK_IMPORT(0,ip)= 0
          read (11) (STACK_IMPORT(k,ip), k=1, N3)
          do is= 1, STACK_IMPORT(N3,ip)
            read (11) NOD_IMPORT(is)          
          enddo

          STACK_EXPORT(0,ip)= 0
          allocate (NOD_EXPORT(STACK_EXPORT(N3,ip)))
          read (11) (STACK_EXPORT(k,ip), k=1, N3)
          do is= 1, STACK_EXPORT(N3,ip)
            read (11) NOD_EXPORT(is)          
          enddo

          read (11) IELMTOT
          do i= 1, IELMTOT
            read (11) iip, IELMLOCAL(i), IWORK(i),                      &
     &           (ICELNOD(i,k), k=1,NODELM(IELMLOCAL(i)))
          enddo
        close (11)
!C===

!C
!C +-----------------+
!C | LOCAL NUMBERING |
!C +-----------------+
!C===
          do i= 1, INODTOTG
            IMASK1(i)= 0
          enddo

          do i= 1, IELMTOTG
            IMASK2(i)= 0
          enddo

          do i= 1, INODTOT
            in= INODLOCAL(i)
            IMASK1(in)= i 
          enddo

          do i= 1, IELMTOT
            in= IELMLOCAL(i)
            IMASK2(in)= i 
          enddo

          NODGRPSTACK(0)= 0
          do ig= 1, NODGRPTOT
            icou= 0
            do is= NODGRPSTACKG(ig-1)+1,NODGRPSTACKG(ig)
              i= NODGRPITEMG(is)
              if (IMASK1(i).ne.0) then
                icou= icou + 1
                  in= NODGRPSTACK(ig-1) + icou
                NODGRPITEM(in)= IMASK1(i)
              endif
            enddo
            NODGRPSTACK(ig)= NODGRPSTACK(ig-1) + icou
          enddo
!C===

!C
!C +-------------------------+
!C | write FINAL LOCAL files |
!C +-------------------------+
!C===
        open (12,file=FILNAME(ip), status='unknown')

        rewind (12)


!C
!C-- section 1.
          write (* ,'("PE:", i5, 4i10)') ip-1, INODTOT, N2,             &
     &           STACK_IMPORT(N3,ip),                                   &
     &           STACK_EXPORT(N3,ip)

!          write(12,'(10i10)')  CoarseGridLevels
!          write(12,'(10i10)')  HOWmanyADAPTATIONs

          write(12,'(10i10)')  ip-1
          write(12,'(10i10)')  N3
          write(12,'(10i10)') (IWKM(inei)-1,inei=1,N3)
!C
!C-- section 2.: element
          write(12,'(10i10)') INODTOT, N2

          do i= 1, INODTOT
            in= INODLOCAL(i)
            write (12,'(i10,i5, 3( 1pe16.6 ))')                         &
     &      HOME_NODE(in,2), HOME_NODE(in,1), (XYZ(in,k),k=1,3)
          enddo

          write(12,'(10i10)') NPC(ip), nELEM_internal(ip)

          write(12,'(10i10)') (IELMTYP(IELMLOCAL(i)), i=1,NPC(ip))
!          write(12,'(5(1pe16.6))') (CONDW(IELMLOCAL(i)), i=1,NPC(ip))

          ELEM_internal_LIST= 0
          icou= 0
          do i= 1, IELMTOT
            icel= IELMLOCAL(i)
            if (HOME_ELEM(icel,1)+1 .eq. ip) then
              icou= icou + 1
              ELEM_internal_LIST(icou)= i
            endif

!            if (IELMMAT(IELMLOCAL(i)).ne.0) write (*,*)                 &
!     &                 HOME_ELEM(icel,2), HOME_ELEM(icel,1)

            write (12,'(i10,2i5,8i10)')                                 &
     &      HOME_ELEM(icel,2), HOME_ELEM(icel,1),                       &
     &      IELMMAT(IELMLOCAL(i)),                                      &
     &     (ICELNOD(i,k),k=1,NODELM(IELMLOCAL(i)))
          enddo

          write (12,'(10i10)')                                          &
     &         (ELEM_internal_LIST(i), i=1,nELEM_internal(ip))

!C
!C-- section 3.
        if (NP.ne.1) then       
          STACK_IMPORT(0,ip)= 0
            write (12,'(10i10)') (STACK_IMPORT(k,ip), k=1, N3)
            if (N3.ne.0)                                                &
     &      write (12,'( 2i10)') (NOD_IMPORT(is),                       &
     &                 HOME_NODE(INODLOCAL(NOD_IMPORT(is)),1),          &
     &                            is= 1, STACK_IMPORT(N3,ip))

          if (N3.ne.0)                                                  &
     &    write(12,'(10i10)') (STACK_EXPORT(inei,ip),                   &
     &                                     inei= 1, N3)
          write(12,'(  i10)') (NOD_EXPORT(is),                          &
     &                         is= 1, STACK_EXPORT(N3,ip))
        endif
          deallocate (NOD_EXPORT)

!C
!C-- section 4.
          call DATA_COMPRESS (NODGRPTOT, IGTOT, NODGRPSTACK, IWORK)
!          if (IWORK(IGTOT).eq.0) IGTOT= 0
          IGTOT= NODGRPTOT

          write (12, '(  i10)')  IGTOT
          write (12, '(10i10)') (IWORK(ig), ig=1,IGTOT)

          if (IGTOT.ne.0) then
          do ig= 1, NODGRPTOT
            nn= NODGRPSTACK(ig) - NODGRPSTACK(ig-1)
            write (12, '(a64)')    NODGRPNAME(ig)
            if (nn.ne.0) then
     
            write (12, '(10i10)') (NODGRPITEM(is),                      &
     &                           is= NODGRPSTACK(ig-1)+1,               &
     &                               NODGRPSTACK(ig) )
            endif
          enddo
          endif

        close (12)
!C===
      enddo

      return
      end

      subroutine DATA_COMPRESS (I1, I2, I_INN, I_OUT)
      integer  I_INN(0:I1), I_OUT(0:I1)

      ITOT= I1
      do i= 0, I1
        I_OUT(i)= I_INN(i)
      enddo

!      do i= 1, I1
!        ip= I_OUT(i)
!        if (ip.eq.I_OUT(i-1)) then
!          ITOT= ITOT - 1
!          do j= i, ITOT
!            I_OUT(j)= I_OUT(j+1)
!          enddo   
!        endif      
!        if (i.eq.ITOT) exit
!      enddo   
!      I2= ITOT

      return
      end
