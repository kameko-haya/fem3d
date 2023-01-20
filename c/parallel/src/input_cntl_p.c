/***
 *** INPUT_CNTL
 ***/
/** PARALLEL VERSION **/
#include <stdio.h>
#include <stdlib.h>
#include "pfem_util.h"
/*** external functions **/
extern void ERROR_EXIT (int, int);
/** **/
void INPUT_CNTL()
{
	FILE *fp;
	
	if( my_rank == 0 ){
		if( (fp=fopen("INPUT.DAT","r")) == NULL){
			fprintf(stdout,"input file cannot be opened!\n");
			exit(1);
		}
		fscanf(fp,"%s",HEADER);
		fscanf(fp,"%d %d",&METHOD,&PRECOND);
		fscanf(fp,"%d",&iterPREmax);
		fscanf(fp,"%d",&ITER);
		fclose(fp);

		if( ( METHOD != 1 ) && ( METHOD != 2 ) &&( METHOD != 3) ){
			ERROR_EXIT(102,my_rank);
		}
		if( ( PRECOND != 0 ) && ( PRECOND != 1 ) ){
			ERROR_EXIT (103, my_rank);
		}
		if( ( ITER <= 0 ) ){
			ERROR_EXIT (104, my_rank);
		}
		if( ( iterPREmax < 1 ) ){
			iterPREmax= 1;
			ERROR_EXIT (111, my_rank);
		}
  		if( ( iterPREmax > 4 ) ){
			iterPREmax= 4;
			ERROR_EXIT (112, my_rank);
		}
	}

	MPI_Bcast(HEADER ,80,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&METHOD , 1,MPI_INTEGER,0,MPI_COMM_WORLD);
	MPI_Bcast(&PRECOND, 1,MPI_INTEGER,0,MPI_COMM_WORLD);
	MPI_Bcast(&iterPREmax, 1,MPI_INTEGER,0,MPI_COMM_WORLD);
	MPI_Bcast(&ITER   , 1,MPI_INTEGER,0,MPI_COMM_WORLD);

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


