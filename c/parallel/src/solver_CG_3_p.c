/***
 *** CG_3
 ***/
/** PARALLEL VERSION **/
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "precision.h"
#include "allocate.h"
extern FILE *fp_log;
/*** external functions **/
extern void SOLVER_SEND_RECV_3 (); 
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
		KREAL RESID,KINT ITER, KINT *ERROR,int my_rank,
		int NEIBPETOT,int NEIBPE[],
		int IMPORT_INDEX[], int IMPORT_ITEM[],
		int EXPORT_INDEX[], int EXPORT_ITEM[],
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

	KREAL *WS,*WR;
	KREAL **WW;

	KINT R=0,Z=1,Q=1,P=2,ZP=3;
	KINT MAXIT;
	KREAL TOL,W,SS;

	double COMPtime,COMMtime,R1;
        double START_TIME,END_TIME;

/** debug **/
/**
	for(i=0;i<9*NPL;i++) {
		fprintf(fp_log,"(S)AL[%d]=%e\n",i,AL[i]);
	}
	for(i=0;i<9*NPU;i++) {
		fprintf(fp_log,"(S)AU[%d]=%e\n",i,AU[i]);
	}
	for(i=0;i<NP+1;i++) {
		fprintf(fp_log,"(S)INL[%d]=%d\n",i,INL[i]);
	}
	for(i=0;i<NP+1;i++) {
		fprintf(fp_log,"(S)INU[%d]=%d\n",i,INU[i]);
	}
	for(i=0;i<NPL;i++) {
		fprintf(fp_log,"(S)IAL[%d]=%d\n",i,IAL[i]);
	}
	for(i=0;i<NPU;i++) {
		fprintf(fp_log,"(S)IAU[%d]=%d\n",i,IAU[i]);
	}
	for(i=0;i<9*NP;i++) {
		fprintf(fp_log,"(S)ALU[%d]=%e\n",i,ALU[i]);
	}
**/
/***
	+-------+
	| INIT. |
	+-------+
***/
	ERROR= 0;

	COMPtime=0.0;
	COMMtime=0.0;

	WW=(KREAL**) allocate_matrix(sizeof(KREAL),4,3*NP);
	WS=(KREAL* ) allocate_vector(sizeof(KREAL),3*NP);
	WR=(KREAL* ) allocate_vector(sizeof(KREAL),3*NP);
  
	MAXIT  = ITER;
	TOL   = RESID;          

	for(i=0;i<NP;i++){
		X[i]=0.0;	
	}
	for(i=0;i<4;i++) for(j=0;j<3*NP;j++) WW[i][j]=0.0;
	for(i=0;i<3*NP;i++) WS[i]=0.0;
	for(i=0;i<3*NP;i++) WR[i]=0.0;
/**
	+-----------------------+
	| {r0}= {b} - [A]{xini} |
	+-----------------------+
***/
/**
 ** INTERFACE data EXCHANGE
 **/

	SOLVER_SEND_RECV_3
	( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,
	  EXPORT_INDEX, EXPORT_ITEM, WS, WR, X , my_rank);
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
		WW[R][3*j-3]=WVAL1;
		WW[R][3*j-2]=WVAL2;
		WW[R][3*j-1]=WVAL3;
	}

	BNRM20= 0.e0;
	for(i=0;i<N;i++){
		BNRM20+= B[3*i]*B[3*i] +B[3*i+1]*B[3*i+1] +B[3*i+2]*B[3*i+2];
	}

	MPI_Allreduce (&BNRM20, &BNRM2, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

	if (BNRM2 == 0.e0) BNRM2= 1.e0;

	ITER = 0;

	S1_TIME= MPI_Wtime();

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
				WW[ZP][3*i-3]=WW[R][3*i-3];
				WW[ZP][3*i-2]=WW[R][3*i-2];
				WW[ZP][3*i-1]=WW[R][3*i-1];
			}

			for( i=1;i<=NP;i++){
				WW[Z][3*i-3]=0.e0;
				WW[Z][3*i-2]=0.e0;
				WW[Z][3*i-1]=0.e0;
			}

			for(iterPRE=1;iterPRE<=iterPREmax;iterPRE++){
				for( i=1+N;i<=NP;i++){
					WW[ZP][3*i-3]=0.0;
					WW[ZP][3*i-2]=0.0;
					WW[ZP][3*i-1]=0.0;
				}
/**
	FORWARD
**/
				for( i=1;i<=N;i++){
					SW1= WW[ZP][3*i-3];
					SW2= WW[ZP][3*i-2];
					SW3= WW[ZP][3*i-1];

					isL=INL[i-1]+1;
					ieL=INL[i];

					for(j=isL;j<=ieL;j++){
						k=IAL[j-1];
						X1= WW[ZP][3*k-3];
						X2= WW[ZP][3*k-2];
						X3= WW[ZP][3*k-1];
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

					WW[ZP][3*i-3]= X1;
					WW[ZP][3*i-2]= X2;
					WW[ZP][3*i-1]= X3;
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
						X1=WW[ZP][3*k-3];
						X2=WW[ZP][3*k-2];
						X3=WW[ZP][3*k-1];
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
					WW[ZP][3*i-3]+= -X1;
					WW[ZP][3*i-2]+= -X2;
					WW[ZP][3*i-1]+= -X3;
				}
/**
	additive Schwartz
**/
				SOLVER_SEND_RECV_3
				( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,
				EXPORT_INDEX, EXPORT_ITEM, WS, WR, WW[ZP], my_rank);
				
				indexA= ZP;
				indexB= Z;

				for( i=1;i<=NP;i++){
					WW[indexB][3*i-3]+=WW[indexA][3*i-3];
					WW[indexB][3*i-2]+=WW[indexA][3*i-2];
					WW[indexB][3*i-1]+=WW[indexA][3*i-1];
				}

				if (iterPRE == iterPREmax ) goto LABEL750;

				for( j=1;j<=N;j++){
					X1= WW[indexB][3*j-3];
					X2= WW[indexB][3*j-2];
					X3= WW[indexB][3*j-1];
					WV1= WW[R][3*j-3] - D[9*j-9]*X1 - D[9*j-8]*X2 - D[9*j-7]*X3;
					WV2= WW[R][3*j-2] - D[9*j-6]*X1 - D[9*j-5]*X2 - D[9*j-4]*X3;
					WV3= WW[R][3*j-1] - D[9*j-3]*X1 - D[9*j-2]*X2 - D[9*j-1]*X3;
					for(k=INL[j-1]+1;k<=INL[j];k++){
						i=IAL[k-1];
						X1= WW[indexB][3*i-3];
						X2= WW[indexB][3*i-2];
						X3= WW[indexB][3*i-1];
						WV1+=  - AL[9*k-9]*X1 - AL[9*k-8]*X2 - AL[9*k-7]*X3;
						WV2+=  - AL[9*k-6]*X1 - AL[9*k-5]*X2 - AL[9*k-4]*X3;
						WV3+=  - AL[9*k-3]*X1 - AL[9*k-2]*X2 - AL[9*k-1]*X3;	
					}
					for(k=INU[j-1]+1;k<=INU[j];k++){
						i=IAU[k-1];
						X1= WW[indexB][3*i-3];
						X2= WW[indexB][3*i-2];
						X3= WW[indexB][3*i-1];
						WV1+=  - AU[9*k-9]*X1 - AU[9*k-8]*X2 - AU[9*k-7]*X3;
						WV2+=  - AU[9*k-6]*X1 - AU[9*k-5]*X2 - AU[9*k-4]*X3;
						WV3+=  - AU[9*k-3]*X1 - AU[9*k-2]*X2 - AU[9*k-1]*X3;	
					}
					WW[indexA][3*j-3]=WV1;
					WW[indexA][3*j-2]=WV2;
					WW[indexA][3*j-1]=WV3;
				}

				LABEL750:
				continue;
			}

		}

		if (PRECOND != 0 ) {
/**
	Block SCALING
**/
			for(i=0;i<N;i++){
				WW[Z][3*i  ]=WW[R][3*i  ];
				WW[Z][3*i+1]=WW[R][3*i+1];
				WW[Z][3*i+2]=WW[R][3*i+2];
			}

			for(i=0;i<N;i++){
				X1=WW[Z][3*i  ];
				X2=WW[Z][3*i+1];
				X3=WW[Z][3*i+2];
				X2= X2 - ALU[9*i+3]*X1;
				X3= X3 - ALU[9*i+6]*X1 - ALU[9*i+7]*X2;
				X3= ALU[9*i+8]*  X3;
				X2= ALU[9*i+4]*( X2 - ALU[9*i+5]*X3 );
				X1= ALU[9*i  ]*( X1 - ALU[9*i+2]*X3 - ALU[9*i+1]*X2);
				WW[Z][3*i  ]= X1;
				WW[Z][3*i+1]= X2;
				WW[Z][3*i+2]= X3;
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
			RHO0+=WW[R][3*i]*WW[Z][3*i] +WW[R][3*i+1]*WW[Z][3*i+1] +WW[R][3*i+2]*WW[Z][3*i+2];
		}

		START_TIME=MPI_Wtime();

		MPI_Allreduce (&RHO0, &RHO, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

		END_TIME=MPI_Wtime();
		COMMtime+=END_TIME-START_TIME;

/*** **/
/***
	+-----------------------------+
	| {p} = {z} if      ITER=1    |
	| BETA= RHO / RHO1  otherwise |
	+-----------------------------+
***/
		if( ITER == 1 ){
			for(i=0;i<N;i++){
				WW[P][3*i  ]=WW[Z][3*i  ];
				WW[P][3*i+1]=WW[Z][3*i+1];
				WW[P][3*i+2]=WW[Z][3*i+2];
			}
		}else{
			BETA= RHO / RHO1;
			for(i=0;i<N;i++){
				WW[P][3*i  ]=WW[Z][3*i  ]+BETA*WW[P][3*i  ];
				WW[P][3*i+1]=WW[Z][3*i+1]+BETA*WW[P][3*i+1];
				WW[P][3*i+2]=WW[Z][3*i+2]+BETA*WW[P][3*i+2];
			}
		}
/*** **/
/***
	+-------------+
	| {q}= [A]{p} |
	+-------------+
***/
/**
   INTERFACE data EXCHANGE
**/
		START_TIME= MPI_Wtime();

		SOLVER_SEND_RECV_3
		( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,
		  EXPORT_INDEX, EXPORT_ITEM, WS, WR, WW[P], my_rank);

		END_TIME= MPI_Wtime();
		COMMtime+=END_TIME-START_TIME;
/**
	 BEGIN calculation
**/
		for( j=0;j<N;j++){
			X1=WW[P][3*j  ];
			X2=WW[P][3*j+1];
			X3=WW[P][3*j+2];
			WVAL1= D[9*j  ]*X1 + D[9*j+1]*X2 + D[9*j+2]*X3;
			WVAL2= D[9*j+3]*X1 + D[9*j+4]*X2 + D[9*j+5]*X3;
			WVAL3= D[9*j+6]*X1 + D[9*j+7]*X2 + D[9*j+8]*X3;

			for(k=INL[j]+1;k<=INL[j+1];k++){
				i=IAL[k-1];
				X1=WW[P][3*i-3];
				X2=WW[P][3*i-2];
				X3=WW[P][3*i-1];
				WVAL1+=   AL[9*k-9]*X1 + AL[9*k-8]*X2 + AL[9*k-7]*X3;
				WVAL2+=   AL[9*k-6]*X1 + AL[9*k-5]*X2 + AL[9*k-4]*X3;
				WVAL3+=   AL[9*k-3]*X1 + AL[9*k-2]*X2 + AL[9*k-1]*X3;	
			}

			for(k=INU[j]+1;k<=INU[j+1];k++){
				i=IAU[k-1];
				X1=WW[P][3*i-3];
				X2=WW[P][3*i-2];
				X3=WW[P][3*i-1];
				WVAL1+=   AU[9*k-9]*X1 + AU[9*k-8]*X2 + AU[9*k-7]*X3;
				WVAL2+=   AU[9*k-6]*X1 + AU[9*k-5]*X2 + AU[9*k-4]*X3;
				WVAL3+=   AU[9*k-3]*X1 + AU[9*k-2]*X2 + AU[9*k-1]*X3;	
			}
			WW[Q][3*j  ]=WVAL1;
			WW[Q][3*j+1]=WVAL2;
			WW[Q][3*j+2]=WVAL3;
		}
/** **/
/***
	+---------------------+
	| ALPHA= RHO / {p}{q} |
	+---------------------+
***/
		C10= 0.e0;
		for(i=0;i<N;i++){
			C10+=WW[P][3*i]*WW[Q][3*i]+WW[P][3*i+1]*WW[Q][3*i+1]+WW[P][3*i+2]*WW[Q][3*i+2];
		}
		
		START_TIME=MPI_Wtime();

		MPI_Allreduce (&C10,&C1, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

		END_TIME=MPI_Wtime();
		COMMtime+=END_TIME-START_TIME;

		ALPHA= RHO / C1;

/** **/
/***
	+----------------------+
	| {x}= {x} + ALPHA*{p} |
	| {r}= {r} - ALPHA*{q} |
	+----------------------+
***/
		for(i=0;i<N;i++){
			X[3*i  ]+=ALPHA *WW[P][3*i  ];
			X[3*i+1]+=ALPHA *WW[P][3*i+1];
			X[3*i+2]+=ALPHA *WW[P][3*i+2];
			WW[R][3*i  ]+= -ALPHA *WW[Q][3*i  ];
			WW[R][3*i+1]+= -ALPHA *WW[Q][3*i+1];
			WW[R][3*i+2]+= -ALPHA *WW[Q][3*i+2];
		}

		DNRM20= 0.e0;
		for(i=0;i<N;i++){
			DNRM20+=WW[R][3*i]*WW[R][3*i]+WW[R][3*i+1]*WW[R][3*i+1]+WW[R][3*i+2]*WW[R][3*i+2];
		}

		START_TIME=MPI_Wtime();

		MPI_Allreduce (&DNRM20,&DNRM2, 1, MPI_DOUBLE,MPI_SUM, MPI_COMM_WORLD);

		END_TIME=MPI_Wtime();
		COMMtime+=END_TIME-START_TIME;
		RESID= sqrt(DNRM2/BNRM2);

/** ##### ITERATION HISTORY ***/
		if( my_rank == 0 ) fprintf(stdout,"%d %e\n",ITER,RESID);
		if( my_rank == 0 ) fprintf(fp_log,"%d %e\n",ITER,RESID);
/** ***/
        	if ( RESID <= TOL   ) break;
        	if ( ITER  == MAXIT ) *ERROR= -300;

        	RHO1 = RHO ;                                                           
	}
/** **/
/***
	INTERFACE data EXCHANGE
***/
	E1_TIME= MPI_Wtime();
	COMPtime= E1_TIME - S1_TIME;

	R1= 100.e0 * ( 1.e0 - COMMtime/COMPtime );

	if (my_rank == 0) {
		fprintf(stdout,"### elapsed    :%e\n",COMPtime);
		fprintf(stdout,"### comm.      :%e\n",COMMtime);
		fprintf(stdout,"### work ratio :%e\n",R1);
		fprintf(fp_log,"### elapsed    :%e\n",COMPtime);
		fprintf(fp_log,"### comm.      :%e\n",COMMtime);
		fprintf(fp_log,"### work ratio :%e\n",R1);
	}

	SOLVER_SEND_RECV_3
	( NP, NEIBPETOT, NEIBPE, IMPORT_INDEX, IMPORT_ITEM,
	  EXPORT_INDEX, EXPORT_ITEM, WS, WR, X, my_rank);

	free ( (KREAL**)WW);
	deallocate_vector ( (KREAL**)WR);
	deallocate_vector( (KREAL**)WS);
}
