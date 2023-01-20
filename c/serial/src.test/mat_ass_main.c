/***
 *** MAT_ASS_MAIN
 **/
#include <stdio.h>
#include <math.h>
#include "pfem_util.h"
#include "allocate.h"
extern FILE *fp_log;
extern void JACOBI();
void MAT_ASS_MAIN()
{
	int i,k,kk;
	int ip,jp,kp;
	int ipn,jpn,kpn;
	int icel;
	int ie,je;
	int iiS,iiE;
	int in1,in2,in3,in4,in5,in6,in7,in8;
	double valA,valB,valX;
	double QP1,QM1,EP1,EM1,TP1,TM1;
	double X1,X2,X3,X4,X5,X6,X7,X8;
	double Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8;
	double Z1,Z2,Z3,Z4,Z5,Z6,Z7,Z8;
	double PNXi,PNYi,PNZi,PNXj,PNYj,PNZj;
	double VOL;
	double coef;
	double a11,a12,a13,a21,a22,a23,a31,a32,a33;

	KINT nodLOCAL[8];

	AL=(KREAL*) allocate_vector(sizeof(KREAL),9*NPL);
	AU=(KREAL*) allocate_vector(sizeof(KREAL),9*NPU);
	B =(KREAL*) allocate_vector(sizeof(KREAL),3*N );
	D =(KREAL*) allocate_vector(sizeof(KREAL),9*N);
	X =(KREAL*) allocate_vector(sizeof(KREAL),3*N);

	for(i=0;i<9*NPL;i++) AL[i]=0.0;
	for(i=0;i<9*NPU;i++) AU[i]=0.0;
	for(i=0;i<3*N ;i++) B[i]=0.0;
	for(i=0;i<9*N ;i++) D[i]=0.0;
	for(i=0;i<3*N ;i++) X[i]=0.0;

	WEI[0]= 1.0000000000e0;
	WEI[1]= 1.0000000000e0;

	POS[0]= -0.5773502692e0;
	POS[1]=  0.5773502692e0;

/***
	INIT.
	PNQ   - 1st-order derivative of shape function by QSI
	PNE   - 1st-order derivative of shape function by ETA
	PNT   - 1st-order derivative of shape function by ZET
***/
	valA=            POISSON /      (1.e0-POISSON);
	valB= (1.e0-2.e0*POISSON)/(2.e0*(1.e0-POISSON));
	valX= ELAST*(1.e0-POISSON)/((1.e0+POISSON)*(1.e0-2.e0*POISSON));

	valA= valA * valX;
	valB= valB * valX;

	for(ip=0;ip<2;ip++){
		for(jp=0;jp<2;jp++){
			for(kp=0;kp<2;kp++){
				QP1= 1.e0 + POS[ip];
				QM1= 1.e0 - POS[ip];
				EP1= 1.e0 + POS[jp];
				EM1= 1.e0 - POS[jp];
				TP1= 1.e0 + POS[kp];
				TM1= 1.e0 - POS[kp];
				SHAPE[ip][jp][kp][0]= O8th * QM1 * EM1 * TM1;
				SHAPE[ip][jp][kp][1]= O8th * QP1 * EM1 * TM1;
				SHAPE[ip][jp][kp][2]= O8th * QP1 * EP1 * TM1;
				SHAPE[ip][jp][kp][3]= O8th * QM1 * EP1 * TM1;
				SHAPE[ip][jp][kp][4]= O8th * QM1 * EM1 * TP1;
				SHAPE[ip][jp][kp][5]= O8th * QP1 * EM1 * TP1;
				SHAPE[ip][jp][kp][6]= O8th * QP1 * EP1 * TP1;
				SHAPE[ip][jp][kp][7]= O8th * QM1 * EP1 * TP1;
				PNQ[jp][kp][0]= - O8th * EM1 * TM1;
				PNQ[jp][kp][1]= + O8th * EM1 * TM1;
				PNQ[jp][kp][2]= + O8th * EP1 * TM1;
				PNQ[jp][kp][3]= - O8th * EP1 * TM1;
				PNQ[jp][kp][4]= - O8th * EM1 * TP1;
				PNQ[jp][kp][5]= + O8th * EM1 * TP1;
				PNQ[jp][kp][6]= + O8th * EP1 * TP1;
				PNQ[jp][kp][7]= - O8th * EP1 * TP1;
				PNE[ip][kp][0]= - O8th * QM1 * TM1;
				PNE[ip][kp][1]= - O8th * QP1 * TM1;
				PNE[ip][kp][2]= + O8th * QP1 * TM1;
				PNE[ip][kp][3]= + O8th * QM1 * TM1;
				PNE[ip][kp][4]= - O8th * QM1 * TP1;
				PNE[ip][kp][5]= - O8th * QP1 * TP1;
				PNE[ip][kp][6]= + O8th * QP1 * TP1;
				PNE[ip][kp][7]= + O8th * QM1 * TP1;
				PNT[ip][jp][0]= - O8th * QM1 * EM1;
				PNT[ip][jp][1]= - O8th * QP1 * EM1;
				PNT[ip][jp][2]= - O8th * QP1 * EP1;
				PNT[ip][jp][3]= - O8th * QM1 * EP1;
				PNT[ip][jp][4]= + O8th * QM1 * EM1;
				PNT[ip][jp][5]= + O8th * QP1 * EM1;
				PNT[ip][jp][6]= + O8th * QP1 * EP1;
				PNT[ip][jp][7]= + O8th * QM1 * EP1;
			}
		}
	}

	for( icel=0;icel< ICELTOT;icel++){
		in1=ICELNOD[icel][0];
		in2=ICELNOD[icel][1];
		in3=ICELNOD[icel][2];
		in4=ICELNOD[icel][3];
		in5=ICELNOD[icel][4];
		in6=ICELNOD[icel][5];
		in7=ICELNOD[icel][6];
		in8=ICELNOD[icel][7];
/**
 **
 ** JACOBIAN & INVERSE JACOBIAN
**/
		nodLOCAL[0]= in1;
		nodLOCAL[1]= in2;
		nodLOCAL[2]= in3;
		nodLOCAL[3]= in4;
		nodLOCAL[4]= in5;
		nodLOCAL[5]= in6;
		nodLOCAL[6]= in7;
		nodLOCAL[7]= in8;

		X1=XYZ[in1-1][0];
		X2=XYZ[in2-1][0];
		X3=XYZ[in3-1][0];
		X4=XYZ[in4-1][0];
		X5=XYZ[in5-1][0];
		X6=XYZ[in6-1][0];
		X7=XYZ[in7-1][0];
		X8=XYZ[in8-1][0];

		Y1=XYZ[in1-1][1];
		Y2=XYZ[in2-1][1];
		Y3=XYZ[in3-1][1];
		Y4=XYZ[in4-1][1];
		Y5=XYZ[in5-1][1];
		Y6=XYZ[in6-1][1];
		Y7=XYZ[in7-1][1];
		Y8=XYZ[in8-1][1];

		Z1=XYZ[in1-1][2];
		Z2=XYZ[in2-1][2];
		Z3=XYZ[in3-1][2];
		Z4=XYZ[in4-1][2];
		Z5=XYZ[in5-1][2];
		Z6=XYZ[in6-1][2];
		Z7=XYZ[in7-1][2];
		Z8=XYZ[in8-1][2];

		JACOBI(DETJ, PNQ, PNE, PNT, PNX, PNY, PNZ,     
				X1, X2, X3, X4, X5, X6, X7, X8,
				Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8,
				Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8);

/**
	CONSTRUCT the GLOBAL MATRIX
**/
		for(ie=0;ie<8;ie++){
			ip=nodLOCAL[ie];
			for(je=0;je<8;je++){
				jp=nodLOCAL[je];

				kk=0;
				if( jp > ip ){
					iiS=indexU[ip-1]+1;
					iiE=indexU[ip  ];
					for( k=iiS;k<=iiE;k++){
						if( itemU[k-1] == jp ){
							kk=k;
							break;
						}
					}
				}

				if( jp < ip ){
					iiS=indexL[ip-1]+1;
					iiE=indexL[ip  ];
					for( k=iiS;k<=iiE;k++){
						if( itemL[k-1] == jp ){
							kk=k;
							break;
						}
					}
				}

				PNXi= 0.e0;
				PNYi= 0.e0;
				PNZi= 0.e0;
				PNXj= 0.e0;
				PNYj= 0.e0;
				PNZj= 0.e0;

				VOL= 0.e0;

				for(kpn=0;kpn<2;kpn++){
					for(jpn=0;jpn<2;jpn++){
						for(ipn=0;ipn<2;ipn++){
							coef= -fabs(DETJ[ipn][jpn][kpn])*WEI[ipn]*WEI[jpn]*WEI[kpn];

							VOL+=coef;
							PNXi= PNX[ipn][jpn][kpn][ie];
							PNYi= PNY[ipn][jpn][kpn][ie];
							PNZi= PNZ[ipn][jpn][kpn][ie];

							PNXj= PNX[ipn][jpn][kpn][je];
							PNYj= PNY[ipn][jpn][kpn][je];
							PNZj= PNZ[ipn][jpn][kpn][je];

							a11= (valX*PNXi*PNXj+valB*(PNYi*PNYj+PNZi*PNZj))*coef;
							a22= (valX*PNYi*PNYj+valB*(PNZi*PNZj+PNXi*PNXj))*coef;
							a33= (valX*PNZi*PNZj+valB*(PNXi*PNXj+PNYi*PNYj))*coef;

							a12= (valA*PNXi*PNYj + valB*PNXj*PNYi)*coef;
							a13= (valA*PNXi*PNZj + valB*PNXj*PNZi)*coef;
							a21= (valA*PNYi*PNXj + valB*PNYj*PNXi)*coef;
							a23= (valA*PNYi*PNZj + valB*PNYj*PNZi)*coef;
							a31= (valA*PNZi*PNXj + valB*PNZj*PNXi)*coef;
							a32= (valA*PNZi*PNYj + valB*PNZj*PNYi)*coef;

							if (jp > ip) {
								AU[9*kk-9]+=a11;
								AU[9*kk-8]+=a12;
								AU[9*kk-7]+=a13;
								AU[9*kk-6]+=a21;
								AU[9*kk-5]+=a22;
								AU[9*kk-4]+=a23;
								AU[9*kk-3]+=a31;
								AU[9*kk-2]+=a32;
								AU[9*kk-1]+=a33;
							}

							if (jp < ip) {
								AL[9*kk-9]+=a11;
								AL[9*kk-8]+=a12;
								AL[9*kk-7]+=a13;
								AL[9*kk-6]+=a21;
								AL[9*kk-5]+=a22;
								AL[9*kk-4]+=a23;
								AL[9*kk-3]+=a31;
								AL[9*kk-2]+=a32;
								AL[9*kk-1]+=a33;
							}

							if (jp == ip) {
								D[9*ip-9]+=a11;
								D[9*ip-8]+=a12;
								D[9*ip-7]+=a13;
								D[9*ip-6]+=a21;
								D[9*ip-5]+=a22;
								D[9*ip-4]+=a23;
								D[9*ip-3]+=a31;
								D[9*ip-2]+=a32;
								D[9*ip-1]+=a33;
							}
						}
					}
				}
			}
		}
	}
}

