      implicit REAL*8 (A-H,O-Z)
      real(kind=8), dimension(:,:), allocatable :: X , Y      
      real(kind=8), dimension(:,:), allocatable :: X0, Y0      
      character(len=80) :: GRIDFILE, HHH
      integer , dimension(:,:), allocatable :: IW
!C
!C +-------+
!C | INIT. |
!C +-------+
!C===
      write (*,*) 'NX ? (must be > 0)'
      read  (*,*)  NX

      NY= NX 
      NZ= NX 

      NXP1= NX + 1
      NYP1= NY + 1
      NZP1= NZ + 1

      DX= 1.d0

      INODTOT=  NXP1*NYP1*NZP1
      ICELTOT=  NX  *NY  *NZ  
      IBNODTOT= NXP1*NYP1

      allocate (IW(IBNODTOT,4))
      IW= 0

      icou= 0
      ib  = 1
      do k= 1, NXP1      
      do j= 1, NYP1
        i= 1
        icou= icou + 1
        ii  = (k-1)*IBNODTOT + (j-1)*NXP1 + i
        IW(icou,ib)= ii
      enddo
      enddo

      icou= 0
      ib  = 2
      do k= 1, NXP1      
        j= 1
      do i= 1, NXP1
        icou= icou + 1
        ii  = (k-1)*IBNODTOT + (j-1)*NXP1 + i
        IW(icou,ib)= ii
      enddo
      enddo

      icou= 0
      ib  = 3
      k= 1
      do j= 1, NYP1
      do i= 1, NXP1
        icou= icou + 1
        ii  = (k-1)*IBNODTOT + (j-1)*NXP1 + i
        IW(icou,ib)= ii
      enddo
      enddo

      icou= 0
      ib  = 4
      k= NXP1      
      do j= 1, NYP1
      do i= 1, NXP1
        icou= icou + 1
        ii  = (k-1)*IBNODTOT + (j-1)*NXP1 + i
        IW(icou,ib)= ii
      enddo
      enddo
!C===

!C
!C +-------------+
!C | GeoFEM data |
!C +-------------+
!C===
      NN= 0
      write (*,*)      'GeoFEM gridfile name ?'
      GRIDFILE= 'cube.0'

      open  (12, file= GRIDFILE, status='unknown',form='formatted')
          write(12,'(10i10)') INODTOT
       
          icou= 0
          do k= 1, NZP1
          do j= 1, NYP1
          do i= 1, NXP1
            XX= dfloat(i-1)*DX
            YY= dfloat(j-1)*DX
            ZZ= dfloat(k-1)*DX

            icou= icou + 1
            write (12,'(i10,3(1pe16.6))') icou, XX, YY, ZZ
          enddo
          enddo
          enddo

          write(12,'(i10)') ICELTOT

          IELMTYPL= 361
          write(12,'(10i10)') (IELMTYPL, i=1,ICELTOT)

          icou= 0
          imat= 1
          do k= 1, NZ
          do j= 1, NY
          do i= 1, NX
            icou= icou + 1
            in1 = (k-1)*IBNODTOT + (j-1)*NXP1 + i
            in2 = in1 + 1
            in3 = in2 + NXP1
            in4 = in3 - 1
            in5 = in1 + IBNODTOT
            in6 = in2 + IBNODTOT
            in7 = in3 + IBNODTOT
            in8 = in4 + IBNODTOT
            write (12,'(10i10)') icou, imat, in1, in2, in3, in4,        &
     &                                       in5, in6, in7, in8
          enddo
          enddo
          enddo

          IGTOT= 4
          write (12,'(10i10)') IGTOT
          write (12,'(10i10)') IBNODTOT, 2*IBNODTOT, 3*IBNODTOT,        &
     &                                   4*IBNODTOT

          HHH= 'Xmin'
          write (12,'(a80)')  HHH
          write (12,'(10i10)') (IW(ii,1), ii=1,IBNODTOT)
          HHH= 'Ymin'
          write (12,'(a80)')  HHH
          write (12,'(10i10)') (IW(ii,2), ii=1,IBNODTOT)
          HHH= 'Zmin'
          write (12,'(a80)')  HHH
          write (12,'(10i10)') (IW(ii,3), ii=1,IBNODTOT)
          HHH= 'Zmax'
          write (12,'(a80)')  HHH
          write (12,'(10i10)') (IW(ii,4), ii=1,IBNODTOT)

          deallocate (IW)

          IGTOT= 1
          write (12,'(10i10)') IGTOT
          write (12,'(10i10)') NX*NY
          HHH= 'ZminE'
          write (12,'(a80)')  HHH
          write (12,'(10i10)') (ii, ii=1, NX*NY)

          IGTOT= 1
          write (12,'(10i10)') IGTOT
          write (12,'(10i10)') NX*NY
          HHH= 'ZmaxS'
          write (12,'(a80)')  HHH
          iS= ICELTOT - NX*NY + 1
          isuf= 6
          write (12,'(10i10)') (ii  , ii=iS, ICELTOT)
          write (12,'(10i10)') (isuf, ii=1, NX*NY)

      close (12)
!C===
      stop
      end

