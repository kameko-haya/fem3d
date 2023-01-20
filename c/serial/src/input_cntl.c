/***
 *** INPUT_CNTL
 ***/
#include <stdio.h>
#include <stdlib.h>
#include "pfem_util.h"
/** **/
void INPUT_CNTL()
{
	FILE *fp;
	
	if( (fp=fopen("INPUT.DAT","r")) == NULL){
		fprintf(stdout,"input file cannot be opened!\n");
		exit(1);
	}
	fscanf(fp,"%s",fname);
	fscanf(fp,"%d %d",&METHOD,&PRECOND);
	fscanf(fp,"%d",&iterPREmax);
	fscanf(fp,"%d",&ITER);
	fscanf(fp, "%lf %lf", &ELAST, &POISSON);
	fclose(fp);

	if( ( iterPREmax < 1 ) ){
		iterPREmax= 1;
	}
  	if( ( iterPREmax > 4 ) ){
		iterPREmax= 4;
	}

	SIGMA_DIAG= 1.0;
	SIGMA     = 0.0;
	RESID     = 1.e-8;
	NSET      = 0;

	pfemRarray[0]= RESID;
	pfemRarray[1]= SIGMA_DIAG;
	pfemRarray[2]= SIGMA;

	pfemIarray[0]= ITER;
	pfemIarray[1]= METHOD;
	pfemIarray[2]= PRECOND;
	pfemIarray[3]= NSET;
	pfemIarray[4]= iterPREmax;
}


