/**
 ** PFEM_INIT
 **/
#include "pfem_util.h"
/**
	INIT. PFEM-FEM process's
**/
void PFEM_INIT(int argc, char* argv[])
{
	int i;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&PETOT);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

	for(i=0;i<100;i++)  pfemRarray[i]=0.0;
	for(i=0;i<100;i++)  pfemIarray[i]=0;
}
