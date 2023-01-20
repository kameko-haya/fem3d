/***
 *** OUTPUT_UCD
 ***/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "pfem_util.h"
#include "allocate.h"
extern FILE *fp_log;
extern void JACOBI();
void OUTPUT_UCD()
{
	FILE *fp;
	int i,ie;
	int N0,N1,N3,N4;
	int in1,in2,in3,in4,in5,in6,in7,in8;
	double ZERO;
	double XX,YY,ZZ;
	double igmax;
	double  VAL_G_max;
  
	char ETYPE[7];
 
/***
	+----------+
	| AVS file |
	+----------+
***/
		if( (fp=fopen("test.inp","w")) == NULL){
			fprintf(stdout,"output file cannot be opened!\n");
			  exit(1);}

        N0= 0;
        N1= 1;
        N3= 3;
        N4= 4;
        ZERO= 0.e0;

		fprintf(fp,"%8d\n",N1);
		fprintf(fp,"data\n");
		fprintf(fp,"step1\n");
		fprintf(fp,"%8d%8d\n",N,ICELTOT);

		for(i=0;i<N;i++){
			XX=XYZ[i][0];
			YY=XYZ[i][1];
			ZZ=XYZ[i][2];
			fprintf(fp,"%8d%16.6e%16.6e%16.6e\n",i+1,XX,YY,ZZ);
		}

		for(ie=0;ie<ICELTOT;ie++){
			strcpy(ETYPE," hex  ");
			in1= ICELNOD[ie][0];
			in2= ICELNOD[ie][1];
			in3= ICELNOD[ie][2];
			in4= ICELNOD[ie][3];
			in5= ICELNOD[ie][4];
			in6= ICELNOD[ie][5];
			in7= ICELNOD[ie][6];
			in8= ICELNOD[ie][7];
			fprintf(fp,"%8d%3d%-8s%8d%8d%8d%8d%8d%8d%8d%8d\n",
				ie+1,N1,ETYPE,in1, in2, in3, in4, in5, in6, in7, in8);
		}

		fprintf(fp,"%3d%3d\n",N4, N0);
		fprintf(fp,"%3d%3d%3d%3d%3d\n", N4, N1, N1, N1, N1);
		fprintf(fp,"disp-x,\n");
		fprintf(fp,"disp-y,\n");
		fprintf(fp,"disp-z,\n");
		fprintf(fp,"sigma-zz,\n");

        igmax    = 0;
        VAL_G_max= 0.e0;

		for(i=0;i<N;i++){
			fprintf(fp,"%8d%16.6e%16.6e%16.6e%16.6e\n",i+1,X[3*i],X[3*i+1],X[3*i+2],SIGMA_N[3*i+2]);
		}
		fclose(fp);

        i= N-1;
		fprintf(stdout,"### DISPLACEMENT at (Xmax,Ymax,Zmax))\n");
		fprintf(stdout,"%d %e %e %e\n",i+1,X[3*i],X[3*i+1],X[3*i+2]);
		fprintf(fp_log,"### DISPLACEMENT at (Xmax,Ymax,Zmax))\n");
		fprintf(fp_log,"%d %e %e %e\n",i+1,X[3*i],X[3*i+1],X[3*i+2]);
}

