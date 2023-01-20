!C
!C***
!C*** MAT_ASS_MAIN
!C***
!C
      subroutine MAT_ASS_MAIN
      use pfem_util
      implicit REAL*8 (A-H,O-Z)
      integer(kind=kint), dimension(  8) :: nodLOCAL

      allocate (AL(9*NPL), AU(9*NPU))
      allocate (B(3*NP), D(9*NP), X(3*NP))

      
      AL= 0.d0
      AU= 0.d0
       B= 0.d0
       X= 0.d0
       D= 0.d0

      WEI(1)= +1.0000000000D+00
      WEI(2)= +1.0000000000D+00

      POS(1)= -0.5773502692D+00
      POS(2)= +0.5773502692D+00

!C
!C-- INIT.
!C     PNQ   - 1st-order derivative of shape function by QSI
!C     PNE   - 1st-order derivative of shape function by ETA
!C     PNT   - 1st-order derivative of shape function by ZET
!C
      valA=            POISSON /      (1.d0-POISSON)
      valB= (1.d0-2.d0*POISSON)/(2.d0*(1.d0-POISSON))
      valX= ELAST*(1.d0-POISSON)/((1.d0+POISSON)*(1.d0-2.d0*POISSON))

      valA= valA * valX
      valB= valB * valX

      do kp= 1, 2
      do jp= 1, 2
      do ip= 1, 2
        QP1= 1.d0 + POS(ip)
        QM1= 1.d0 - POS(ip)
        EP1= 1.d0 + POS(jp)
        EM1= 1.d0 - POS(jp)
        TP1= 1.d0 + POS(kp)
        TM1= 1.d0 - POS(kp)
        SHAPE(ip,jp,kp,1)= O8th * QM1 * EM1 * TM1
        SHAPE(ip,jp,kp,2)= O8th * QP1 * EM1 * TM1
        SHAPE(ip,jp,kp,3)= O8th * QP1 * EP1 * TM1
        SHAPE(ip,jp,kp,4)= O8th * QM1 * EP1 * TM1
        SHAPE(ip,jp,kp,5)= O8th * QM1 * EM1 * TP1
        SHAPE(ip,jp,kp,6)= O8th * QP1 * EM1 * TP1
        SHAPE(ip,jp,kp,7)= O8th * QP1 * EP1 * TP1
        SHAPE(ip,jp,kp,8)= O8th * QM1 * EP1 * TP1
        PNQ(jp,kp,1)= - O8th * EM1 * TM1
        PNQ(jp,kp,2)= + O8th * EM1 * TM1
        PNQ(jp,kp,3)= + O8th * EP1 * TM1
        PNQ(jp,kp,4)= - O8th * EP1 * TM1
        PNQ(jp,kp,5)= - O8th * EM1 * TP1
        PNQ(jp,kp,6)= + O8th * EM1 * TP1
        PNQ(jp,kp,7)= + O8th * EP1 * TP1
        PNQ(jp,kp,8)= - O8th * EP1 * TP1
        PNE(ip,kp,1)= - O8th * QM1 * TM1
        PNE(ip,kp,2)= - O8th * QP1 * TM1
        PNE(ip,kp,3)= + O8th * QP1 * TM1
        PNE(ip,kp,4)= + O8th * QM1 * TM1
        PNE(ip,kp,5)= - O8th * QM1 * TP1
        PNE(ip,kp,6)= - O8th * QP1 * TP1
        PNE(ip,kp,7)= + O8th * QP1 * TP1
        PNE(ip,kp,8)= + O8th * QM1 * TP1
        PNT(ip,jp,1)= - O8th * QM1 * EM1
        PNT(ip,jp,2)= - O8th * QP1 * EM1
        PNT(ip,jp,3)= - O8th * QP1 * EP1
        PNT(ip,jp,4)= - O8th * QM1 * EP1
        PNT(ip,jp,5)= + O8th * QM1 * EM1
        PNT(ip,jp,6)= + O8th * QP1 * EM1
        PNT(ip,jp,7)= + O8th * QP1 * EP1
        PNT(ip,jp,8)= + O8th * QM1 * EP1
      enddo
      enddo
      enddo


      do icel= 1, ICELTOT
        in1= ICELNOD(icel,1)
        in2= ICELNOD(icel,2)
        in3= ICELNOD(icel,3)
        in4= ICELNOD(icel,4)
        in5= ICELNOD(icel,5)
        in6= ICELNOD(icel,6)
        in7= ICELNOD(icel,7)
        in8= ICELNOD(icel,8)

!C
!C
!C== JACOBIAN & INVERSE JACOBIAN
        nodLOCAL(1)= in1
        nodLOCAL(2)= in2
        nodLOCAL(3)= in3
        nodLOCAL(4)= in4
        nodLOCAL(5)= in5
        nodLOCAL(6)= in6
        nodLOCAL(7)= in7
        nodLOCAL(8)= in8

        X1= 0.d0
        Y1= 0.d0
        Z1= 0.d0

        X2= 1.d0
        Y2= 0.d0
        Z2= 0.d0

        X3= 1.d0
        Y3= 1.d0
        Z3= 0.d0

        X4= 0.d0
        Y4= 1.d0
        Z4= 0.d0

        X5= 0.d0
        Y5= 0.d0
        Z5= 1.d0

        X6= 1.d0
        Y6= 0.d0
        Z6= 1.d0

        X7= 1.d0
        Y7= 1.d0
        Z7= 1.d0

        X8= 0.d0
        Y8= 1.d0
        Z8= 1.d0


!C        X1= XYZ(in1,1)
!C        X2= XYZ(in2,1)
!C        X3= XYZ(in3,1)
!C        X4= XYZ(in4,1)
!C        X5= XYZ(in5,1)
!C        X6= XYZ(in6,1)
!C        X7= XYZ(in7,1)
!C        X8= XYZ(in8,1)

!C       Y1= XYZ(in1,2)
!C       Y2= XYZ(in2,2)
!C       Y3= XYZ(in3,2)
!C       Y4= XYZ(in4,2)
!C       Y5= XYZ(in5,2)
!C       Y6= XYZ(in6,2)
!C       Y7= XYZ(in7,2)
!C       Y8= XYZ(in8,2)

!C       Z1= XYZ(in1,3)
!C       Z2= XYZ(in2,3)
!C       Z3= XYZ(in3,3)
!C       Z4= XYZ(in4,3)
!C       Z5= XYZ(in5,3)
!C       Z6= XYZ(in6,3)
!C       Z7= XYZ(in7,3)
!C       Z8= XYZ(in8,3)

        call JACOBI (DETJ, PNQ, PNE, PNT, PNX, PNY, PNZ,                &
     &               X1, X2, X3, X4, X5, X6, X7, X8,                    &
     &               Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8,                    &
     &               Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8 )

!C
!C== CONSTRUCT the GLOBAL MATRIX
        do ie= 1, 8
          ip = nodLOCAL(ie)
        do je= 1, 8
          jp = nodLOCAL(je)

          kk= 0
          if (jp.gt.ip) then
            iiS= indexU(ip-1) + 1
            iiE= indexU(ip  )
            do k= iiS, iiE
              if ( itemU(k).eq.jp ) then
                kk= k
                exit
              endif
            enddo
          endif

          if (jp.lt.ip) then
            iiS= indexL(ip-1) + 1
            iiE= indexL(ip  )
            do k= iiS, iiE
              if ( itemL(k).eq.jp) then
                kk= k
                exit
              endif
            enddo
          endif

          PNXi= 0.d0
          PNYi= 0.d0
          PNZi= 0.d0
          PNXj= 0.d0
          PNYj= 0.d0
          PNZj= 0.d0

          VOL= 0.d0
          
          do kpn= 1, 2
          do jpn= 1, 2
          do ipn= 1, 2
            coef= dabs(DETJ(ipn,jpn,kpn))*WEI(ipn)*WEI(jpn)*WEI(kpn)

            VOL= VOL + coef
            PNXi= PNX(ipn,jpn,kpn,ie)
            PNYi= PNY(ipn,jpn,kpn,ie)
            PNZi= PNZ(ipn,jpn,kpn,ie)

            PNXj= PNX(ipn,jpn,kpn,je)
            PNYj= PNY(ipn,jpn,kpn,je)
            PNZj= PNZ(ipn,jpn,kpn,je)

            a11= (valX*PNXi*PNXj+valB*(PNYi*PNYj+PNZi*PNZj))*coef
            a22= (valX*PNYi*PNYj+valB*(PNZi*PNZj+PNXi*PNXj))*coef
            a33= (valX*PNZi*PNZj+valB*(PNXi*PNXj+PNYi*PNYj))*coef

            a12= (valA*PNXi*PNYj + valB*PNXj*PNYi)*coef
            a13= (valA*PNXi*PNZj + valB*PNXj*PNZi)*coef
            a21= (valA*PNYi*PNXj + valB*PNYj*PNXi)*coef
            a23= (valA*PNYi*PNZj + valB*PNYj*PNZi)*coef
            a31= (valA*PNZi*PNXj + valB*PNZj*PNXi)*coef
            a32= (valA*PNZi*PNYj + valB*PNZj*PNYi)*coef

            if (jp.gt.ip) then
              AU(9*kk-8)= AU(9*kk-8) + a11
              AU(9*kk-7)= AU(9*kk-7) + a12
              AU(9*kk-6)= AU(9*kk-6) + a13
              AU(9*kk-5)= AU(9*kk-5) + a21
              AU(9*kk-4)= AU(9*kk-4) + a22
              AU(9*kk-3)= AU(9*kk-3) + a23
              AU(9*kk-2)= AU(9*kk-2) + a31
              AU(9*kk-1)= AU(9*kk-1) + a32
              AU(9*kk  )= AU(9*kk  ) + a33
            endif

            if (jp.lt.ip) then
              AL(9*kk-8)= AL(9*kk-8) + a11
              AL(9*kk-7)= AL(9*kk-7) + a12
              AL(9*kk-6)= AL(9*kk-6) + a13
              AL(9*kk-5)= AL(9*kk-5) + a21
              AL(9*kk-4)= AL(9*kk-4) + a22
              AL(9*kk-3)= AL(9*kk-3) + a23
              AL(9*kk-2)= AL(9*kk-2) + a31
              AL(9*kk-1)= AL(9*kk-1) + a32
              AL(9*kk  )= AL(9*kk  ) + a33
            endif

            if (jp.eq.ip) then
              D(9*ip-8)= D(9*ip-8) + a11
              D(9*ip-7)= D(9*ip-7) + a12
              D(9*ip-6)= D(9*ip-6) + a13
              D(9*ip-5)= D(9*ip-5) + a21
              D(9*ip-4)= D(9*ip-4) + a22
              D(9*ip-3)= D(9*ip-3) + a23
              D(9*ip-2)= D(9*ip-2) + a31
              D(9*ip-1)= D(9*ip-1) + a32
              D(9*ip  )= D(9*ip  ) + a33
            endif
          enddo
          enddo
          enddo
        enddo
        enddo
      enddo

      open (50,file='debug', status='unknown')
      rewind(50)
      do i=1, 100
        do j=indexL(i-1)+1, indexL(i)
          write(50,*) AL(9*j-8), AL(9*j-7), AL(9*j-6)
          write(50,*) AL(9*j-5), AL(9*j-4), AL(9*j-3)
          write(50,*) AL(9*j-2), AL(9*j-1), AL(9*j)
        enddo
        write(50,*) D(9*i-8), D(9*i-7), D(9*i-6)
        write(50,*) D(9*i-5), D(9*i-4), D(9*i-3)
        write(50,*) D(9*i-2), D(9*i-1), D(9*i)
        do j=indexU(i-1)+1, indexU(i)
          write(50,*) AU(9*j-8), AU(9*j-7), AU(9*j-6)
          write(50,*) AU(9*j-5), AU(9*j-4), AU(9*j-3)
          write(50,*) AU(9*j-2), AU(9*j-1), AU(9*j)
        enddo
      enddo
      close(50)



      return
      end
