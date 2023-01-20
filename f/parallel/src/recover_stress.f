!C
!C***
!C*** RECOVER_STRESS
!C***
!C
      subroutine RECOVER_STRESS
      use pfem_util
      use solver_SR_3

      implicit REAL*8 (A-H,O-Z)

      integer(kind=kint), dimension(  8) :: nodLOCAL
      real(kind=kreal)  , dimension(:),allocatable :: WS, WR

      allocate (SIGMA_N(3*NP), TAU_N(3*NP))

      SIGMA_N= 0.d0
      TAU_N  = 0.d0
      B      = 0.d0

      do icel= 1, ICELTOT
        E0  = ELAST
        POI0= POISSON

        valA=            POI0 /      (1.d0-POI0)
        valB= (1.d0-2.d0*POI0)/(2.d0*(1.d0-POI0))
        valX= E0*  (1.d0-POI0)/((1.d0+POI0)*(1.d0-2.d0*POI0))

        valA= valA * valX
        valB= valB * valX

        in1= ICELNOD(icel,1)
        in2= ICELNOD(icel,2)
        in3= ICELNOD(icel,3)
        in4= ICELNOD(icel,4)
        in5= ICELNOD(icel,5)
        in6= ICELNOD(icel,6)
        in7= ICELNOD(icel,7)
        in8= ICELNOD(icel,8)

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

!C
!C-- JACOBIAN & inv-JACOBIAN
        call JACOBI (DETJ, PNQ, PNE, PNT, PNX, PNY, PNZ,                &
     &                     X1, X2, X3, X4, X5, X6, X7, X8,              &
     &                     Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8,              &
     &                     Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8 )

!C
!C-- MATRIX
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

          UUj= X(3*jp-2)
          VVj= X(3*jp-1)
          WWj= X(3*jp  )
          UUi= X(3*ip-2)
          VVi= X(3*ip-1)
          WWi= X(3*ip  )

          EPS_xx= 0.d0
          EPS_yy= 0.d0
          EPS_zz= 0.d0
          GAM_xy= 0.d0
          GAM_xz= 0.d0
          GAM_yz= 0.d0
          GAM_xy= 0.d0
          GAM_xz= 0.d0
          GAM_yz= 0.d0
          
          VOL   = 0.d0
          do kpn= 1, 2
          do jpn= 1, 2
          do ipn= 1, 2
            coef= dabs(DETJ(ipn,jpn,kpn))*WEI(ipn)*WEI(jpn)*WEI(kpn)

            SHi= SHAPE(ipn,jpn,kpn,ie) * coef
            SHj= SHAPE(ipn,jpn,kpn,je)

            EPS_xx= EPS_xx + SHi*PNX(ipn,jpn,kpn,je)
            EPS_yy= EPS_yy + SHi*PNY(ipn,jpn,kpn,je)
            EPS_zz= EPS_zz + SHi*PNZ(ipn,jpn,kpn,je)
            GAM_xy= GAM_xy + SHi*PNX(ipn,jpn,kpn,je) * VVj              &
     &                     + SHi*PNY(ipn,jpn,kpn,je) * UUj  
            GAM_xz= GAM_xz + SHi*PNX(ipn,jpn,kpn,je) * WWj              &
     &                     + SHi*PNZ(ipn,jpn,kpn,je) * UUj  
            GAM_yz= GAM_yz + SHi*PNY(ipn,jpn,kpn,je) * WWj              &
     &                     + SHi*PNZ(ipn,jpn,kpn,je) * VVj  
            VOL   = VOL    +SHi
          enddo
          enddo
          enddo

          EPS_xx= EPS_xx * UUj
          EPS_yy= EPS_yy * VVj
          EPS_zz= EPS_zz * WWj

          SIGMA_N(3*ip-2)= SIGMA_N(3*ip-2) + valX*EPS_xx + valA*EPS_yy   &
     &                                                   + valA*EPS_zz
          SIGMA_N(3*ip-1)= SIGMA_N(3*ip-1) + valA*EPS_xx + valX*EPS_yy   &
     &                                                   + valA*EPS_zz
          SIGMA_N(3*ip  )= SIGMA_N(3*ip  ) + valA*EPS_xx + valA*EPS_yy   &
     &                                                   + valX*EPS_zz

          TAU_N(3*ip-2)= TAU_N(3*ip-2) + GAM_xy*valB
          TAU_N(3*ip-1)= TAU_N(3*ip-1) + GAM_xz*valB
          TAU_N(3*ip  )= TAU_N(3*ip  ) + GAM_yz*valB

          if (ip.eq.jp) B(ip)= B(ip) + VOL
        enddo
        enddo
      enddo

!C
!C-- NODAL VALUE
      do i= 1, N
        RB= 1.d0/B(i)
        SIGMA_N(3*i-2)= SIGMA_N(3*i-2)*RB
        SIGMA_N(3*i-1)= SIGMA_N(3*i-1)*RB
        SIGMA_N(3*i  )= SIGMA_N(3*i  )*RB

        TAU_N(3*i-2)  = TAU_N(3*i-2)*RB
        TAU_N(3*i-1)  = TAU_N(3*i-1)*RB
        TAU_N(3*i  )  = TAU_N(3*i  )*RB
      enddo

!C
!C-- UPDATE
      allocate (WS(3*NP), WR(3*NP))
      call SOLVER_SEND_RECV_3                                           &
     &   ( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,            &
     &     EXPORT_INDEX, EXPORT_ITEM, WS, WR, SIGMA_N, my_rank)
      call SOLVER_SEND_RECV_3                                           &
     &   ( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,            &
     &     EXPORT_INDEX, EXPORT_ITEM, WS, WR, TAU_N, my_rank)

      deallocate (WS, WR)

      end subroutine RECOVER_STRESS
