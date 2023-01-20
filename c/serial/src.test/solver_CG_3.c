/***
 *** CG_3
 ***/
#include <stdio.h>
#include <math.h>
#include "precision.h"
#include "allocate.h"
extern FILE *fp_log;
/***
	CG_3 solves the linear system Ax = b with 3*3 block matrix 
	using the Conjugate Gradient iterative method with the following
	preconditioners for SMP nodes:
 ***/
void  CG_3(
		KINT N,KINT NP,KINT NPL, KINT NPU,KREAL D[],
		KREAL AL[],KINT INL[], KINT IAL[],
		KREAL AU[],KINT INU[], KINT IAU[],
		KREAL B[],KREAL X[],KREAL ALU[],
		KREAL RESID,KINT ITER, KINT *ERROR,
		KINT  PRECOND, KINT iterPREmax)
{
	int i,j,k;
	int ieL,isL,ieU,isU;
	double X1,X2,X3;
	double WVAL1,WVAL2,WVAL3;
	double SW1,SW2,SW3;
	double WV1,WV2,WV3;
	double BNRM20,BNRM2,DNRM20,DNRM2;
	double S1_TIME,E1_TIME;
	double ALPHA,BETA;
	double C1,C10,RHO,RHO0,RHO1;
	int    iterPRE;
	int    indexA,indexB;

	KREAL **WW;

	KINT R=0,Z=1,Q=1,P=2,ZP=3;
	KINT MAXIT;
	KREAL TOL;

	double COMPtime;

/** debug **/
/**
	for(i=0;i<9*NPL;i++) {
		fprintf(fp_log,"(S)AL[%d]=%e\n",i,AL[i]);
	}
	for(i=0;i<9*NPU;i++) {
		fprintf(fp_log,"(S)AU[%d]=%e\n",i,AU[i]);
	}
	for(i=0;i<N+1;i++) {
		fprintf(fp_log,"(S)INL[%d]=%d\n",i,INL[i]);
	}
	for(i=0;i<N+1;i++) {
		fprintf(fp_log,"(S)INU[%d]=%d\n",i,INU[i]);
	}
	for(i=0;i<NPL;i++) {
		fprintf(fp_log,"(S)IAL[%d]=%d\n",i,IAL[i]);
	}
	for(i=0;i<NPU;i++) {
		fprintf(fp_log,"(S)IAU[%d]=%d\n",i,IAU[i]);
	}
	for(i=0;i<9*N;i++) {
		fprintf(fp_log,"(S)ALU[%d]=%e\n",i,ALU[i]);
	}
	exit(0);
**/
/***
	+-------+
	| INIT. |
	+-------+
***/
	ERROR= 0;

	WW=(KREAL**) allocate_matrix(sizeof(KREAL),3*N,4);
  
	MAXIT  = ITER;
	TOL   = RESID;          

	for(i=0;i<N;i++){
		X[i]=0.0;	
	}
	for(i=0;i<3*N;i++) for(j=0;j<4;j++) WW[i][j]=0.0;
/**
	+-----------------------+
	| {r0}= {b} - [A]{xini} |
	+-----------------------+
***/
/**
 ** BEGIN calculation
 **/
	for(j=1;j<=N;j++){
		X1= X[3*j-3];
		X2= X[3*j-2];
		X3= X[3*j-1];
        WVAL1= B[3*j-3] - D[9*j-9]*X1 - D[9*j-8]*X2 - D[9*j-7]*X3;
        WVAL2= B[3*j-2] - D[9*j-6]*X1 - D[9*j-5]*X2 - D[9*j-4]*X3;
        WVAL3= B[3*j-1] - D[9*j-3]*X1 - D[9*j-2]*X2 - D[9*j-1]*X3;
		
		for( k=INL[j-1]+1;k<=INL[j];k++){
			i= IAL[k-1];
			X1= X[3*i-3];
			X2= X[3*i-2];
			X3= X[3*i-1];
			WVAL1+=  - AL[9*k-9]*X1 - AL[9*k-8]*X2 - AL[9*k-7]*X3;
			WVAL2+=  - AL[9*k-6]*X1 - AL[9*k-5]*X2 - AL[9*k-4]*X3;
			WVAL3+=  - AL[9*k-3]*X1 - AL[9*k-2]*X2 - AL[9*k-1]*X3;
		}

		for( k=INU[j-1]+1;k<=INU[j];k++){
			i=IAU[k-1];
			X1= X[3*i-3];
			X2= X[3*i-2];
			X3= X[3*i-1];
			WVAL1+=  - AU[9*k-9]*X1 - AU[9*k-8]*X2 - AU[9*k-7]*X3;
			WVAL2+=  - AU[9*k-6]*X1 - AU[9*k-5]*X2 - AU[9*k-4]*X3;
			WVAL3+=  - AU[9*k-3]*X1 - AU[9*k-2]*X2 - AU[9*k-1]*X3;
		}
		WW[3*j-3][R]=WVAL1;
		WW[3*j-2][R]=WVAL2;
		WW[3*j-1][R]=WVAL3;
	}

	BNRM20= 0.e0;
	for(i=0;i<N;i++){
		BNRM20+= B[3*i]*B[3*i] +B[3*i+1]*B[3*i+1] +B[3*i+2]*B[3*i+2];
	}

	BNRM2= BNRM20;

/**debug 
	fprintf(fp_log,"###bnrm2=%e\n",BNRM2);fflush(fp_log);
**/
	if (BNRM2 == 0.e0) BNRM2= 1.e0;

	ITER = 0;

	for( ITER=1;ITER<= MAXIT;ITER++){
/**
	************************************************* Conjugate Gradient Iteration
**/
/**
	+----------------+
	| {z}= [Minv]{r} |
	+----------------+
**/
		if( PRECOND == 0 ){
/**
	 Block SSOR
**/
			for( i=1;i<=N;i++){
				WW[3*i-3][ZP]=WW[3*i-3][R];
				WW[3*i-2][ZP]=WW[3*i-2][R];
				WW[3*i-1][ZP]=WW[3*i-1][R];
			}

			for( i=1;i<=N;i++){
				WW[3*i-3][Z]=0.e0;
				WW[3*i-2][Z]=0.e0;
				WW[3*i-1][Z]=0.e0;
			}

/**
	FORWARD
**/
				for( i=1;i<=N;i++){
					SW1= WW[3*i-3][ZP];
					SW2= WW[3*i-2][ZP];
					SW3= WW[3*i-1][ZP];

					isL=INL[i-1]+1;
					ieL=INL[i];

					for(j=isL;j<=ieL;j++){
						k=IAL[j-1];
						X1= WW[3*k-3][ZP];
						X2= WW[3*k-2][ZP];
						X3= WW[3*k-1][ZP];
						SW1+=  - AL[9*j-9]*X1 - AL[9*j-8]*X2 - AL[9*j-7]*X3;
						SW2+=  - AL[9*j-6]*X1 - AL[9*j-5]*X2 - AL[9*j-4]*X3;
						SW3+=  - AL[9*j-3]*X1 - AL[9*j-2]*X2 - AL[9*j-1]*X3;
					}

					X1= SW1;
					X2= SW2;
					X3= SW3;
					X2= X2 - ALU[9*i-6]*X1;
					X3= X3 - ALU[9*i-3]*X1 - ALU[9*i-2]*X2;
					X3= ALU[9*i-1]*  X3;
					X2= ALU[9*i-5]*( X2 - ALU[9*i-4]*X3 );
					X1= ALU[9*i-9]*( X1 - ALU[9*i-7]*X3 - ALU[9*i-8]*X2);

					WW[3*i-3][ZP]= X1;
					WW[3*i-2][ZP]= X2;
					WW[3*i-1][ZP]= X3;
				}
/**
	BACKWARD
**/
				for(i=N;i>=1;i--){
					isU= INU[i-1]+1;
					ieU= INU[i];
					SW1= 0.e0;
					SW2= 0.e0;
					SW3= 0.e0;
					for(j=isU;j<=ieU;j++){
						k=IAU[j-1];
						X1=WW[3*k-3][ZP];
						X2=WW[3*k-2][ZP];
						X3=WW[3*k-1][ZP];
						SW1+=  + AU[9*j-9]*X1 + AU[9*j-8]*X2 + AU[9*j-7]*X3;
						SW2+=  + AU[9*j-6]*X1 + AU[9*j-5]*X2 + AU[9*j-4]*X3;
						SW3+=  + AU[9*j-3]*X1 + AU[9*j-2]*X2 + AU[9*j-1]*X3;
					}
					X1= SW1;
					X2= SW2;
					X3= SW3;
					X2= X2 - ALU[9*i-6]*X1;
					X3= X3 - ALU[9*i-3]*X1 - ALU[9*i-2]*X2;
					X3= ALU[9*i-1]*  X3;
					X2= ALU[9*i-5]*( X2 - ALU[9*i-4]*X3 );
					X1= ALU[9*i-9]*( X1 - ALU[9*i-7]*X3 - ALU[9*i-8]*X2);
					WW[3*i-3][ZP]+= -X1;
					WW[3*i-2][ZP]+= -X2;
					WW[3*i-1][ZP]+= -X3;
				}

				for( i=1;i<=N;i++){
					WW[3*i-3][Z]=WW[3*i-3][ZP];
					WW[3*i-2][Z]=WW[3*i-2][ZP];
					WW[3*i-1][Z]=WW[3*i-1][ZP];
				}

		}

		if (PRECOND != 0 ) {
/**
	Block SCALING
**/
			for(i=0;i<N;i++){
				WW[3*i  ][Z]=WW[3*i  ][R];
				WW[3*i+1][Z]=WW[3*i+1][R];
				WW[3*i+2][Z]=WW[3*i+2][R];
			}

			for(i=0;i<N;i++){
				X1=WW[3*i  ][Z];
				X2=WW[3*i+1][Z];
				X3=WW[3*i+2][Z];
				X2= X2 - ALU[9*i+3]*X1;
				X3= X3 - ALU[9*i+6]*X1 - ALU[9*i+7]*X2;
				X3= ALU[9*i+8]*  X3;
				X2= ALU[9*i+4]*( X2 - ALU[9*i+5]*X3 );
				X1= ALU[9*i  ]*( X1 - ALU[9*i+2]*X3 - ALU[9*i+1]*X2);
				WW[3*i  ][Z]= X1;
				WW[3*i+1][Z]= X2;
				WW[3*i+2][Z]= X3;
			}
		}
/**** ****/  
/***
	+---------------+
	| {RHO}= {r}{z} |
	+---------------+
***/
		RHO0= 0.e0;
	
		for(i=0;i<N;i++){
			RHO0+=WW[3*i][R]*WW[3*i][Z] +WW[3*i+1][R]*WW[3*i+1][Z] +WW[3*i+2][R]*WW[3*i+2][Z];
		}

		RHO= RHO0;

/**debug  
	fprintf(fp_log,"###iter %d rho=%e\n",ITER,RHO);fflush(fp_log);
**/

/*** **/
/***
	+-----------------------------+
	| {p} = {z} if      ITER=1    |
	| BETA= RHO / RHO1  otherwise |
	+-----------------------------+
***/
		if( ITER == 1 ){
			for(i=0;i<N;i++){
				WW[3*i  ][P]=WW[3*i  ][Z];
				WW[3*i+1][P]=WW[3*i+1][Z];
				WW[3*i+2][P]=WW[3*i+2][Z];
			}
		}else{
			BETA= RHO / RHO1;
			for(i=0;i<N;i++){
				WW[3*i  ][P]=WW[3*i  ][Z] +BETA*WW[3*i  ][P];
				WW[3*i+1][P]=WW[3*i+1][Z] +BETA*WW[3*i+1][P];
				WW[3*i+2][P]=WW[3*i+2][Z] +BETA*WW[3*i+2][P];
			}
		}
/*** **/
/***
	+-------------+
	| {q}= [A]{p} |
	+-------------+
***/      
/**
	 BEGIN calculation
**/
		for( j=0;j<N;j++){
			X1=WW[3*j  ][P];
			X2=WW[3*j+1][P];
			X3=WW[3*j+2][P];
			WVAL1= D[9*j  ]*X1 + D[9*j+1]*X2 + D[9*j+2]*X3;
			WVAL2= D[9*j+3]*X1 + D[9*j+4]*X2 + D[9*j+5]*X3;
			WVAL3= D[9*j+6]*X1 + D[9*j+7]*X2 + D[9*j+8]*X3;
			for(k=INL[j]+1;k<=INL[j+1];k++){
				i=IAL[k-1];
				X1=WW[3*i-3][P];
				X2=WW[3*i-2][P];
				X3=WW[3*i-1][P];
				WVAL1+=   AL[9*k-9]*X1 + AL[9*k-8]*X2 + AL[9*k-7]*X3;
				WVAL2+=   AL[9*k-6]*X1 + AL[9*k-5]*X2 + AL[9*k-4]*X3;
				WVAL3+=   AL[9*k-3]*X1 + AL[9*k-2]*X2 + AL[9*k-1]*X3;	
			}
			for(k=INU[j]+1;k<=INU[j+1];k++){
				i=IAU[k-1];
				X1=WW[3*i-3][P];
				X2=WW[3*i-2][P];
				X3=WW[3*i-1][P];
				WVAL1+=   AU[9*k-9]*X1 + AU[9*k-8]*X2 + AU[9*k-7]*X3;
				WVAL2+=   AU[9*k-6]*X1 + AU[9*k-5]*X2 + AU[9*k-4]*X3;
				WVAL3+=   AU[9*k-3]*X1 + AU[9*k-2]*X2 + AU[9*k-1]*X3;	
			}
			WW[3*j  ][Q]=WVAL1;
			WW[3*j+1][Q]=WVAL2;
			WW[3*j+2][Q]=WVAL3;
		}
/** **/
/***
	+---------------------+
	| ALPHA= RHO / {p}{q} |
	+---------------------+
***/
		C10= 0.e0;
		for(i=0;i<N;i++){
			C10+=WW[3*i][P]*WW[3*i][Q]+WW[3*i+1][P]*WW[3*i+1][Q]+WW[3*i+2][P]*WW[3*i+2][Q];
		}
		C1=C10;

		ALPHA= RHO / C1;

/**
	debug fprintf(fp_log,"###iter %d alpha=%e\n",ITER,ALPHA);fflush(fp_log);
**/

/** **/
/***
	+----------------------+
	| {x}= {x} + ALPHA*{p} |
	| {r}= {r} - ALPHA*{q} |
	+----------------------+
***/
		for(i=0;i<N;i++){
			X[3*i  ]+=ALPHA *WW[3*i  ][P];
			X[3*i+1]+=ALPHA *WW[3*i+1][P];
			X[3*i+2]+=ALPHA *WW[3*i+2][P];
			WW[3*i  ][R]+= -ALPHA *WW[3*i  ][Q];
			WW[3*i+1][R]+= -ALPHA *WW[3*i+1][Q];
			WW[3*i+2][R]+= -ALPHA *WW[3*i+2][Q];
		}

		DNRM20= 0.e0;
		for(i=0;i<N;i++){
			DNRM20+=WW[3*i][R]*WW[3*i][R]+WW[3*i+1][R]*WW[3*i+1][R]+WW[3*i+2][R]*WW[3*i+2][R];
		}
		DNRM2= DNRM20;
		RESID= sqrt(DNRM2/BNRM2);

/** ##### ITERATION HISTORY ***/
		fprintf(stdout,"%d %e\n",ITER,RESID);
		fprintf(fp_log,"%d %e\n",ITER,RESID);
/** ***/
        if ( RESID <= TOL   ) break;
        if ( ITER  == MAXIT ) *ERROR= -300;

        RHO1 = RHO ;                                                           
	}
/** **/
/***
	INTERFACE data EXCHANGE
***/

        free(WW);
}
