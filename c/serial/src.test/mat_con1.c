/***
 *** MAT_CON1
 ***/
#include <stdio.h>
#include "pfem_util.h"
#include "allocate.h"
extern FILE* fp_log;
void MAT_CON1()
{
	int i,k,kk;

	indexL=(KINT*)allocate_vector(sizeof(KINT),N+1);
	indexU=(KINT*)allocate_vector(sizeof(KINT),N+1);
	for(i=0;i<N+1;i++) indexL[i]=0;
	for(i=0;i<N+1;i++) indexU[i]=0;

	for(i=0;i<N;i++){
		indexL[i+1]=indexL[i]+INL[i];
		indexU[i+1]=indexU[i]+INU[i];
	}

	NPL=indexL[N];
	NPU=indexU[N];

	itemL=(KINT*)allocate_vector(sizeof(KINT),NPL);
	itemU=(KINT*)allocate_vector(sizeof(KINT),NPU);

	for(i=0;i<N;i++){
		for(k=0;k<INL[i];k++){
			kk=k+indexL[i];
			itemL[kk]=IAL[i][k];
		}
		for(k=0;k<INU[i];k++){
			kk=k+indexU[i];
			itemU[kk]=IAU[i][k];
		}
	}
 

	deallocate_vector(INL);
	deallocate_vector(INU);
	deallocate_vector(IAL);
	deallocate_vector(IAU);
}
