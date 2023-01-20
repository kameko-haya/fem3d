/***
	program SOLVER33_TEST_SCALAR
***/
#include <stdio.h>
#include <stdlib.h>
FILE* fp_log;
#define GLOBAL_VALUE_DEFINE
#include "pfem_util.h"
//#include "solver33.h"
extern void INPUT_CNTL();
extern void INPUT_GRID();
extern void MAT_CON0();
extern void MAT_CON1();
extern void MAT_ASS_MAIN();
extern void MAT_ASS_BC();
extern void SOLVE33();
extern void RECOVER_STRESS();
extern void OUTPUT_UCD();
int main()
{

/** Logfile for debug **/
	if( (fp_log=fopen("log.log","w")) == NULL){
		fprintf(stdout,"input file cannot be opened!\n");
		exit(1);
	}

/**
	+-------+
	| INIT. |
	+-------+
**/ 
	INPUT_CNTL();
	INPUT_GRID();
/**
	+---------------------+
	| matrix connectivity |
	+---------------------+
**/

	MAT_CON0();
	MAT_CON1();
/**
	+-----------------+
	| MATRIX assemble |
	+-----------------+
**/
	MAT_ASS_MAIN();
	MAT_ASS_BC()  ;
/**
	+--------+
	| SOLVER |
	+--------+
**/
	SOLVE33();

/**
	+--------+
	| OUTPUT |
	+--------+
**/
	RECOVER_STRESS();
	OUTPUT_UCD()    ;
}

      

