      subroutine OUTPUT_UCD 
      use pfem_util

      implicit REAL*8 (A-H,O-Z)

      real(kind=kreal), dimension(:), allocatable :: XYZ_G
      real(kind=kreal), dimension(:), allocatable :: VAL_G, WSarray
      integer, dimension(:), allocatable :: NODEflag
      integer, dimension(:), allocatable :: ICELNOD_L, ICELNOD_G
      integer, dimension(:), allocatable :: PEnode, PEelem

      integer, dimension(:), allocatable :: rcountsN, displsN      
      integer, dimension(:), allocatable :: rcountsE, displsE      

      character(len=6) :: ETYPE

!C
!C +----------+
!C | AVS file |
!C +----------+
!C===
        open (21 ,file= 'test.inp', status='unknown')

        N0= 0
        N1= 1
        N3= 3
        N4= 4
        ZERO= 0.d0

        write (21,'(i8)')  N1
        write (21,'(a)' )  'data'
        write (21,'(a)' )  'step1'
        write (21,'(5i8)')  N, ICELTOT
        do i= 1, N
          XX= XYZ(i,1)
          YY= XYZ(i,2)
          ZZ= XYZ(i,3)
          write (21,'(i8,3(1pe16.6))') i, XX, YY, ZZ
        enddo
        do ie= 1, ICELTOT
          ETYPE= 'hex   '
          in1= ICELNOD(ie,1)
          in2= ICELNOD(ie,2)
          in3= ICELNOD(ie,3)
          in4= ICELNOD(ie,4)
          in5= ICELNOD(ie,5)
          in6= ICELNOD(ie,6)
          in7= ICELNOD(ie,7)
          in8= ICELNOD(ie,8)

          write (21,'(i8,i3,1x,a6,1x,8i8)')                             &
     &      ie, N1, ETYPE, in1, in2, in3, in4, in5, in6, in7, in8

        enddo

        write (21,'(10i3)')  N4, N0
        write (21,'(10i3)')  N4, N1, N1, N1, N1
        write (21,'(a  )') 'disp-x, '
        write (21,'(a  )') 'disp-y, '
        write (21,'(a  )') 'disp-z, '
        write (21,'(a  )') 'sigma-zz, '

        igmax    = 0
        VAL_G_max= 0.d0
        do i= 1, N
          write (21,'(i8, 4(1pe16.6))')  i, (X(3*i-3+k),k=1,3),         &
     &                                       SIGMA_N(3*i)
        enddo
        close (21)

        i= N
        write (*,'(/a))')  '### DISPLACEMENT at (Xmax,Ymax,Zmax))'
        write (*,'(i8, 3(1pe16.6),//)')  i, (X(3*i-3+k), k=1,3)
!C===
      end subroutine output_ucd

