/***
 *** MAT_ASS_BC
 ***/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pfem_util.h"
#include "allocate.h"
extern FILE *fp_log;
void MAT_ASS_BC()
{
	int i,j,k,in,ib,ib0,icel;
	int in1,in2,in3,in4,in5,in6,in7,in8;
	int iq1,iq2,iq3,iq4,iq5,iq6,iq7,iq8;
	int iS,iE;
	double STRESS,VAL;

	IWKX=(KINT**) allocate_matrix(sizeof(KINT),NP,2);
	for(i=0;i<NP;i++) for(j=0;j<2;j++) IWKX[i][j]=0;

/**
 Z=Zmax
**/
        for(in=0;in<NP;in++) IWKX[in][0]=0;

        ib0=-1;

        for( ib0=0;ib0<NODGRPtot;ib0++){
                if( strcmp(NODGRP_NAME[ib0].name,"Zmax") == 0 ) break;
        }

        for( ib=NODGRP_INDEX[ib0];ib<NODGRP_INDEX[ib0+1];ib++){
                in=NODGRP_ITEM[ib];
                IWKX[in-1][0]=1;
        }

        for(in=0;in<NP;in++){
                if( IWKX[in][0] == 1 ){
                        B[3*in  ]= B[3*in  ] - D[9*in+2]*1.e0;
                        B[3*in+1]= B[3*in+1] - D[9*in+5]*1.e0;
                        D[9*in+2]= 0.e0;
                        D[9*in+5]= 0.e0;
                        D[9*in+6]= 0.e0;
                        D[9*in+7]= 0.e0;
                        D[9*in+8]= 1.e0;
                        B[3*in+2]= 1.e0;
                        iS= indexL[in]+1;
                        iE= indexL[in+1];
                        for(k=iS;k<=iE;k++){
                                AL[9*k-3]= 0.e0;
                                AL[9*k-2]= 0.e0;
                                AL[9*k-1]= 0.e0;
                        }
                        iS= indexU[in]+1;
                        iE= indexU[in+1];
                        for(k=iS;k<=iE;k++){
                                AU[9*k-3]= 0.e0;
                                AU[9*k-2]= 0.e0;
                                AU[9*k-1]= 0.e0;
                        }
                }
        }
        for(in=0;in<NP;in++){
                iS= indexL[in]+1;
        iE= indexL[in+1];
                for(k=iS;k<=iE;k++){
                        if (IWKX[itemL[k-1]-1][0] == 1 ) {
                                B[3*in  ]= B[3*in  ] - AL[9*k-7]*1.e0;
                                B[3*in+1]= B[3*in+1] - AL[9*k-4]*1.e0;
                                B[3*in+2]= B[3*in+2] - AL[9*k-1]*1.e0;
                                AL[9*k-7]= 0.e0;
                                AL[9*k-4]= 0.e0;
                                AL[9*k-1]= 0.e0;
                        }
                }
                iS= indexU[in]+1;
        iE= indexU[in+1];

                for(k=iS;k<=iE;k++){
                        if (IWKX[itemU[k-1]-1][0] == 1 ) {
                                B[3*in  ]= B[3*in  ] - AU[9*k-7]*1.e0;
                                B[3*in+1]= B[3*in+1] - AU[9*k-4]*1.e0;
                                B[3*in+2]= B[3*in+2] - AU[9*k-1]*1.e0;
                                AU[9*k-7]= 0.e0;
                                AU[9*k-4]= 0.e0;
                                AU[9*k-1]= 0.e0;
                        }
                }
        }
/**
	Z=Zmin
**/
	for(in=0;in<NP;in++) IWKX[in][0]=0;

	ib0=-1;

	for( ib0=0;ib0<NODGRPtot;ib0++){
		if( strcmp(NODGRP_NAME[ib0].name,"Zmin") == 0 ) break;
	}
   
	for( ib=NODGRP_INDEX[ib0];ib<NODGRP_INDEX[ib0+1];ib++){
		in=NODGRP_ITEM[ib];
		IWKX[in-1][0]=1;
	}

	for(in=0;in<N;in++){
		if( IWKX[in][0] == 1 ){
			D[9*in+2]= 0.e0;
			D[9*in+5]= 0.e0;
			D[9*in+6]= 0.e0;
			D[9*in+7]= 0.e0;
			D[9*in+8]= 1.e0;
			B[3*in+2]= 0.e0;
			iS= indexL[in]+1;
			iE= indexL[in+1];
			for(k=iS;k<=iE;k++){
				AL[9*k-3]= 0.e0;
				AL[9*k-2]= 0.e0;
				AL[9*k-1]= 0.e0;
			}
			iS= indexU[in]+1;
			iE= indexU[in+1];
			for(k=iS;k<=iE;k++){
				AU[9*k-3]= 0.e0;
				AU[9*k-2]= 0.e0;
				AU[9*k-1]= 0.e0;
			}
		}
	}
	for(in=0;in<NP;in++){
		iS= indexL[in]+1;
        iE= indexL[in+1];
		for(k=iS;k<=iE;k++){
			if (IWKX[itemL[k-1]-1][0] == 1 ) {
				AL[9*k-7]= 0.e0;
				AL[9*k-4]= 0.e0;
				AL[9*k-1]= 0.e0;
			}
		}
		iS= indexU[in]+1;
        iE= indexU[in+1];

		for(k=iS;k<=iE;k++){
			if (IWKX[itemU[k-1]-1][0] == 1 ) {
				AU[9*k-7]= 0.e0;
				AU[9*k-4]= 0.e0;
				AU[9*k-1]= 0.e0;
			}
		}
	}
/**
	X=Xmin
**/
	for(in=0;in<NP;in++) IWKX[in][0]=0;

	ib0=-1;

	for( ib0=0;ib0<NODGRPtot;ib0++){
		if( strcmp(NODGRP_NAME[ib0].name,"Xmin") == 0 ) break;
	}
   
	for( ib=NODGRP_INDEX[ib0];ib<NODGRP_INDEX[ib0+1];ib++){
		in=NODGRP_ITEM[ib];
		IWKX[in-1][0]=1;
	}

	for(in=0;in<N;in++){
		if( IWKX[in][0] == 1 ){
			D[9*in  ]= 1.e0;
			D[9*in+1]= 0.e0;
			D[9*in+2]= 0.e0;
			D[9*in+3]= 0.e0;
			D[9*in+6]= 0.e0;
			B[3*in  ]= 0.e0;
		
			iS= indexL[in]+1;
			iE= indexL[in+1];
			for(k=iS;k<=iE;k++){
				AL[9*k-9]= 0.e0;
				AL[9*k-8]= 0.e0;
				AL[9*k-7]= 0.e0;
			}

			iS= indexU[in]+1;
			iE= indexU[in+1];
			for(k=iS;k<=iE;k++){
				AU[9*k-9]= 0.e0;
				AU[9*k-8]= 0.e0;
				AU[9*k-7]= 0.e0;
			}
		}
	}

	for(in=0;in<NP;in++){
		iS= indexL[in]+1;
        iE= indexL[in+1];
		for(k=iS;k<=iE;k++){
			if (IWKX[itemL[k-1]-1][0] == 1 ) {
				AL[9*k-9]= 0.e0;
				AL[9*k-6]= 0.e0;
				AL[9*k-3]= 0.e0;
			}
		}

		iS= indexU[in]+1;
        iE= indexU[in+1];
		for(k=iS;k<=iE;k++){
			if (IWKX[itemU[k-1]-1][0] == 1 ) {
				AU[9*k-9]= 0.e0;
				AU[9*k-6]= 0.e0;
				AU[9*k-3]= 0.e0;
			}
		}
	}
/**
	Y=Ymin
**/
	for(in=0;in<NP;in++) IWKX[in][0]=0;

	ib0=-1;

	for( ib0=0;ib0<NODGRPtot;ib0++){
		if( strcmp(NODGRP_NAME[ib0].name,"Ymin") == 0 ) break;
	}
   

	for( ib=NODGRP_INDEX[ib0];ib<NODGRP_INDEX[ib0+1];ib++){
		in=NODGRP_ITEM[ib];
		IWKX[in-1][0]=1;
	}

	for(in=0;in<N;in++){
		if( IWKX[in][0] == 1 ){
			D[9*in+1]= 0.e0;
			D[9*in+4]= 1.e0;
			D[9*in+7]= 0.e0;
			D[9*in+3]= 0.e0;
			D[9*in+5]= 0.e0;
			B[3*in+1]= 0.e0;

			iS= indexL[in]+1;
			iE= indexL[in+1];
			for(k=iS;k<=iE;k++){
				AL[9*k-6]= 0.e0;
				AL[9*k-5]= 0.e0;
				AL[9*k-4]= 0.e0;
			}

			iS= indexU[in]+1;
			iE= indexU[in+1];
			for(k=iS;k<=iE;k++){
				AU[9*k-6]= 0.e0;
				AU[9*k-5]= 0.e0;
				AU[9*k-4]= 0.e0;
			}
		}
	}

	for(in=0;in<NP;in++){
		iS= indexL[in]+1;
		iE= indexL[in+1];
		for(k=iS;k<=iE;k++){
			if (IWKX[itemL[k-1]-1][0] == 1 ) {
				AL[9*k-8]= 0.e0;
				AL[9*k-5]= 0.e0;
				AL[9*k-2]= 0.e0;
			}
		}
	
		iS= indexU[in]+1;
		iE= indexU[in+1];
		for(k=iS;k<=iE;k++){
			if (IWKX[itemU[k-1]-1][0] == 1 ) {
				AU[9*k-8]= 0.e0;
				AU[9*k-5]= 0.e0;
				AU[9*k-2]= 0.e0;
			}
		}
	}
/** debug **/
//	for(i=0;i<9*NPL;i++) {
//		fprintf(fp_log,"(A)AL[%d]=%e\n",i,AL[i]);
//	}
//	for(i=0;i<9*NPU;i++) {
//		fprintf(fp_log,"(A)AU[%d]=%e\n",i,AU[i]);
//	}
//	exit(0);

}
