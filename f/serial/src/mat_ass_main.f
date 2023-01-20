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
      allocate (B(3*N), D(9*N), X(3*N))
      
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

        X1= XYZ(in1,1)
        X2= XYZ(in2,1)
        X3= XYZ(in3,1)
        X4= XYZ(in4,1)
        X5= XYZ(in5,1)
        X6= XYZ(in6,1)
        X7= XYZ(in7,1)
        X8= XYZ(in8,1)

        Y1= XYZ(in1,2)
        Y2= XYZ(in2,2)
        Y3= XYZ(in3,2)
        Y4= XYZ(in4,2)
        Y5= XYZ(in5,2)
        Y6= XYZ(in6,2)
        Y7= XYZ(in7,2)
        Y8= XYZ(in8,2)

        Z1= XYZ(in1,3)
        Z2= XYZ(in2,3)
        Z3= XYZ(in3,3)
        Z4= XYZ(in4,3)
        Z5= XYZ(in5,3)
        Z6= XYZ(in6,3)
        Z7= XYZ(in7,3)
        Z8= XYZ(in8,3)

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
            coef= -dabs(DETJ(ipn,jpn,kpn))*WEI(ipn)*WEI(jpn)*WEI(kpn)

            VOL= VOL + coef
            PNXi= PNX(ipn,jpn,kpn,ie)
            PNYi= PNY(ipn,jpn,kpn,ie)
            PNZi= PNZ(ipn,jpn,kpn,ie)

            PNXj= PNX(ipn,jpn,kpn,je)
            PNYj= PNY(ipn,jpn,kpn,je)
            PNZj= PNZ(ipn,jpn,kpn,je)

            a11= (valX*(PNXi*PNXj)+valB*(PNYi*PNYj+PNZi*PNZj))*coef
            a22= (valX*(PNYi*PNYj)+valB*(PNZi*PNZj+PNXi*PNXj))*coef
            a33= (valX*(PNZi*PNZj)+valB*(PNXi*PNXj+PNYi*PNYj))*coef

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

      return
      end
