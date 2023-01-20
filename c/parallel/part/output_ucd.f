      subroutine OUTPUT_UCD
      use  geofem_util
      use  analyzer
      use  partitioner

      implicit REAL*8 (A-H,O-Z)
      character(len=6) :: ETYPE

!C
!C +----------+
!C | AVS file |
!C +----------+
!C===
        open (21 ,file= 'part.inp', status='unknown')

        N0= 0
        N1= 1
        N3= 3
        N4= 4
        ZERO= 0.d0

        write (21,'(i8)')  N1
        write (21,'(a)' )  'data'
        write (21,'(a)' )  'step1'
        write (21,'(5i8)')  N, IELMTOTG
        do i= 1, N
          XX= XYZ(i,1)
          YY= XYZ(i,2)
          ZZ= XYZ(i,3)
          write (21,'(i8,3(1pe16.6))') i, XX, YY, ZZ
        enddo
        do ie= 1, IELMTOTG
          ETYPE= 'hex   '
          in1= ICELNODg(ie,1)
          in2= ICELNODg(ie,2)
          in3= ICELNODg(ie,3)
          in4= ICELNODg(ie,4)
          in5= ICELNODg(ie,5)
          in6= ICELNODg(ie,6)
          in7= ICELNODg(ie,7)
          in8= ICELNODg(ie,8)
       

          write (21,'(i8,i3,1x,a6,1x,8i8)')                             &
     &      ie, N1, ETYPE, in1, in2, in3, in4, in5, in6, in7, in8

        enddo

        write (21,'(10i3)')  N0, N1
        write (21,'(10i3)')  N1, N1
        write (21,'(a  )') 'color, '


        do ie= 1, IELMTOTG
          ip1= HOME_NODE(ICELNODg(ie,1),1)
          ip2= HOME_NODE(ICELNODg(ie,2),1)
          ip3= HOME_NODE(ICELNODg(ie,3),1)
          ip4= HOME_NODE(ICELNODg(ie,4),1)
          ip5= HOME_NODE(ICELNODg(ie,5),1)
          ip6= HOME_NODE(ICELNODg(ie,6),1)
          ip7= HOME_NODE(ICELNODg(ie,7),1)
          ip8= HOME_NODE(ICELNODg(ie,8),1)

          if (ip1.eq.ip2 .and. ip1.eq.ip3 .and.
     &        ip1.eq.ip4 .and. ip1.eq.ip5 .and.
     &        ip1.eq.ip6 .and. ip1.eq.ip7 .and.
     &        ip1.eq.ip8) then
            NCC= ip1
           else
            NCC= NP * 2
          endif

          write (21,'(i8, 4(1pe16.6))')  ie, dfloat(NCC)
         enddo
        close (21)

!C===
      end subroutine output_ucd
