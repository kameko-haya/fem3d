#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pfem_util.h"
#include "allocate.h"
extern FILE *fp_log;
extern void CG_3();
void SOLVE33()
{
	int i,j,k,ii,L;
	KREAL ALU[3][3];
	KREAL PW[3];
	double AL0;

	int  ERROR, ICFLAG=0;
	CHAR_LENGTH BUF;

/**
	+------------+
	| PARAMETERs |
	+------------+
**/
	ITER      = pfemIarray[0];
	METHOD    = pfemIarray[1];
	PRECOND   = pfemIarray[2];
	NSET      = pfemIarray[3];
    iterPREmax= pfemIarray[4];
	NREST     = pfemIarray[5];

	RESID     = pfemRarray[0];
	SIGMA_DIAG= pfemRarray[1];

	if( iterPREmax  < 1   ) iterPREmax= 1;
	if (iterPREmax  > 4   ) iterPREmax= 4;

/**
	+-----------+
	| BLOCK LUs |
	+-----------+
**/
	if( ICFLAG == 0 ){
		ALUG  =(KREAL*)allocate_vector(sizeof(KREAL),9*N);
		ICFLAG= 1;
		strcpy(BUF.name,"### LINEAR SOLVER:  3x3 Block" );

		if (METHOD == 1) strcat(BUF.name,"ssCG"); 
		if (METHOD == 2) strcat(BUF.name,"ssBiCGSTAB");

		if (PRECOND == 0) {
                    strcat(BUF.name,"BILU(0)-no ASDD");
		}

		if (PRECOND != 0) strcat(BUF.name,"Block Scaling");

			fprintf(stdout,"%s\n",BUF.name);
			fprintf(fp_log,"%s\n",BUF.name);
	}

	if( NSET == 0 ){
		for(i=0;i<9*N;i++) ALUG[i]=0.0;

		for( ii=0;ii<N;ii++){
			ALU[0][0]= D[9*ii  ]*SIGMA_DIAG;
            ALU[0][1]= D[9*ii+1];
            ALU[0][2]= D[9*ii+2];
            ALU[1][0]= D[9*ii+3];
            ALU[1][1]= D[9*ii+4]*SIGMA_DIAG;
            ALU[1][2]= D[9*ii+5];
            ALU[2][0]= D[9*ii+6];
            ALU[2][1]= D[9*ii+7];
            ALU[2][2]= D[9*ii+8]*SIGMA_DIAG;

			for(k=1;k<=3;k++){
				L=k;
				AL0=fabs(ALU[L-1][k-1]);
				for( i=k+1;i<=3;i++){
					if( fabs(ALU[i-1][k-1]) > AL0 ){
						L=i;
						AL0=fabs(ALU[L-1][k-1]);
					}
				}

				ALU[k-1][k-1]= 1.e0/ALU[k-1][k-1];
				for(i=k+1;i<=3;i++){
					ALU[i-1][k-1]*=ALU[k-1][k-1];
					for(j=k+1;j<=3;j++){
						PW[j-1]=ALU[i-1][j-1] - ALU[i-1][k-1]*ALU[k-1][j-1];
					}
					for(j=k+1;j<=3;j++){
						ALU[i-1][j-1]=PW[j-1];
					}
				}
			}
			ALUG[9*ii  ]=ALU[0][0];
			ALUG[9*ii+1]=ALU[0][1];
			ALUG[9*ii+2]=ALU[0][2];
			ALUG[9*ii+3]=ALU[1][0];
			ALUG[9*ii+4]=ALU[1][1];
			ALUG[9*ii+5]=ALU[1][2];
			ALUG[9*ii+6]=ALU[2][0];
			ALUG[9*ii+7]=ALU[2][1];
			ALUG[9*ii+8]=ALU[2][2];
		}
	}
/***
	+------------------+
	| ITERATIVE solver |
	+------------------+
***/
	if (METHOD == 1 ) {
		CG_3( N, NP, NPL, NPU, D, AL, indexL, itemL, AU, indexU, itemU,
			  B, X,  ALUG, RESID, ITER, &ERROR,
			  PRECOND, iterPREmax);
	}

	ITERactual= ITER;
}

