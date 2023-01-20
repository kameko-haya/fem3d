/***
 *** RECOVER_STRESS
 ***/
/** PARALLEL VERSION **/
#include <math.h>
#include "pfem_util.h"
#include "allocate.h"
extern void JACOBI();
extern void SOLVER_SEND_RECV_3();
void RECOVER_STRESS()
{
	int i,k,kk,icel;
	int ie,je,ip,jp;
	int ipn,jpn,kpn;
	int iiS,iiE;
	double RB;
	double UUi,VVi,WWi,UUj,VVj,WWj;
	double valX,valA,valB,E0,POI0,VOL,coef;
	int in1,in2,in3,in4,in5,in6,in7,in8;
	double X1,X2,X3,X4,X5,X6,X7,X8;
	double Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8;
	double Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8;
	double SHi,SHj;
	double EPS_xx,EPS_yy,EPS_zz,GAM_xy,GAM_xz,GAM_yz;

	KINT nodLOCAL[8];
	KREAL *WS,*WR;

	SIGMA_N=(KREAL*)allocate_vector(sizeof(KREAL),3*NP);
	TAU_N  =(KREAL*)allocate_vector(sizeof(KREAL),3*NP);

	for(i=0;i<3*NP;i++) SIGMA_N[i]=0.0;
	for(i=0;i<3*NP;i++) TAU_N[i]=0.0;
	for(i=0;i<3*NP;i++) B[i]=0.0;
	
	for( icel=0;icel< ICELTOT;icel++){
		E0  = ELAST;
		POI0= POISSON;

        valA=            POI0 /      (1.e0-POI0);
        valB= (1.e0-2.e0*POI0)/(2.e0*(1.e0-POI0));
        valX= E0*  (1.e0-POI0)/((1.e0+POI0)*(1.e0-2.e0*POI0));

        valA= valA * valX;
        valB= valB * valX;

        in1= ICELNOD[icel][0];
        in2= ICELNOD[icel][1];
        in3= ICELNOD[icel][2];
        in4= ICELNOD[icel][3];
        in5= ICELNOD[icel][4];
        in6= ICELNOD[icel][5];
        in7= ICELNOD[icel][6];
        in8= ICELNOD[icel][7];

        nodLOCAL[0]= in1;
        nodLOCAL[1]= in2;
        nodLOCAL[2]= in3;
        nodLOCAL[3]= in4;
        nodLOCAL[4]= in5;
        nodLOCAL[5]= in6;
        nodLOCAL[6]= in7;
        nodLOCAL[7]= in8;

        X1= XYZ[in1-1][0];
        X2= XYZ[in2-1][0];
        X3= XYZ[in3-1][0];
        X4= XYZ[in4-1][0];
        X5= XYZ[in5-1][0];
        X6= XYZ[in6-1][0];
        X7= XYZ[in7-1][0];
        X8= XYZ[in8-1][0];

        Y1= XYZ[in1-1][1];
        Y2= XYZ[in2-1][1];
        Y3= XYZ[in3-1][1];
        Y4= XYZ[in4-1][1];
        Y5= XYZ[in5-1][1];
        Y6= XYZ[in6-1][1];
        Y7= XYZ[in7-1][1];
        Y8= XYZ[in8-1][1];

        Z1= XYZ[in1-1][2];
        Z2= XYZ[in2-1][2];
        Z3= XYZ[in3-1][2];
        Z4= XYZ[in4-1][2];
        Z5= XYZ[in5-1][2];
        Z6= XYZ[in6-1][2];
        Z7= XYZ[in7-1][2];
        Z8= XYZ[in8-1][2];

/**
	JACOBIAN & inv-JACOBIAN
**/
        JACOBI (DETJ, PNQ, PNE, PNT, PNX, PNY, PNZ,               
				X1, X2, X3, X4, X5, X6, X7, X8,              
				Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8,             
				Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8 );

/**
	MATRIX
**/
		for(ie=0;ie<8;ie++){
			ip = nodLOCAL[ie];
			for(je=0;je<8;je++){
				jp = nodLOCAL[je];

				kk=0;
				if( jp > ip ){
					iiS=indexU[ip-1]+1;
					iiE=indexU[ip  ];
					for(k=iiS;k<=iiE;k++){
						if( itemU[k-1] == jp ) {
							kk=k;
							break;
						}
					}
				}
	
				if( jp < ip ){
					iiS=indexL[ip-1]+1;
					iiE=indexL[ip  ];
					for(k=iiS;k<=iiE;k++){
						if( itemL[k-1] == jp ) {
							kk=k;
							break;
						}
					}
				}
  
				UUj= X[3*jp-3];
				VVj= X[3*jp-2];
				WWj= X[3*jp-1];
				UUi= X[3*ip-3];
				VVi= X[3*ip-2];
				WWi= X[3*ip-1];

				EPS_xx= 0.e0;
				EPS_yy= 0.e0;
				EPS_zz= 0.e0;
				GAM_xy= 0.e0;
				GAM_xz= 0.e0;
				GAM_yz= 0.e0;
				GAM_xy= 0.e0;
				GAM_xz= 0.e0;
				GAM_yz= 0.e0;
          
				VOL   = 0.e0;

				for( ipn=0;ipn<2;ipn++){
					for( jpn=0;jpn<2;jpn++){
						for( kpn=0;kpn<2;kpn++){
							coef= fabs(DETJ[ipn][jpn][kpn])*WEI[ipn]*WEI[jpn]*WEI[kpn];
							SHi= SHAPE[ipn][jpn][kpn][ie] * coef;
							SHj= SHAPE[ipn][jpn][kpn][je];

							EPS_xx+= SHi*PNX[ipn][jpn][kpn][je];
							EPS_yy+= SHi*PNY[ipn][jpn][kpn][je];
							EPS_zz+= SHi*PNZ[ipn][jpn][kpn][je];
							GAM_xy+= SHi*PNX[ipn][jpn][kpn][je]* VVj
							+ SHi*PNY[ipn][jpn][kpn][je] * UUj;
							GAM_xz+= SHi*PNX[ipn][jpn][kpn][je] * WWj
							+ SHi*PNZ[ipn][jpn][kpn][je] * UUj ;
							GAM_yz+= GAM_yz + SHi*PNY[ipn][jpn][kpn][je] * WWj
							+ SHi*PNZ[ipn][jpn][kpn][je] * VVj;  
							VOL   = VOL    +SHi;
						}
					}
				}
				EPS_xx= EPS_xx * UUj;
				EPS_yy= EPS_yy * VVj;
				EPS_zz= EPS_zz * WWj;

				SIGMA_N[3*ip-3]+=valX*EPS_xx + valA*EPS_yy + valA*EPS_zz;
				SIGMA_N[3*ip-2]+=valA*EPS_xx + valX*EPS_yy + valA*EPS_zz;
				SIGMA_N[3*ip-1]+=valA*EPS_xx + valA*EPS_yy + valX*EPS_zz;
				TAU_N[3*ip-3]+=GAM_xy*valB;
				TAU_N[3*ip-2]+=GAM_xz*valB;
				TAU_N[3*ip-1]+=GAM_yz*valB;
				if (ie==je) B[ip-1]+=VOL;
			}
		}
	}
/****
	NODAL VALUE
***/
	for(i=0;i<N;i++){
		RB=1.0e0/B[i];

        SIGMA_N[3*i  ]*=RB;
        SIGMA_N[3*i+1]*=RB;
        SIGMA_N[3*i+2]*=RB;

		TAU_N[3*i  ]*=RB;
        TAU_N[3*i+1]*=RB;
        TAU_N[3*i+2]*=RB;

	}
/****
	UPDATE
***/
	WS=(KREAL*)allocate_vector(sizeof(KREAL),3*NP);
	WR=(KREAL*)allocate_vector(sizeof(KREAL),3*NP);
	SOLVER_SEND_RECV_3( NP,NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,
		                EXPORT_INDEX,EXPORT_ITEM, WS,WR,SIGMA_N, my_rank);
	SOLVER_SEND_RECV_3( NP,NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,
		                EXPORT_INDEX,EXPORT_ITEM, WS,WR,TAU_N, my_rank);

	deallocate_vector((KREAL*)WS);
	deallocate_vector((KREAL*)WR);

}
