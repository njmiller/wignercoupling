/* This is the wrapper routine for calling my Wigner 3j/6j/9j routines
 * */

#include <stdio.h>
#include <stdlib.h>
#include "idl_export.h"
#include "WignerCoupling.h"

//Since everything is a scalar or 1-d array, don't need to do lots of stuff needed 
//for wig3j6j. This is the Wigner3jvect calling function
IDL_VPTR wig3jv(int argc, IDL_VPTR argv[], char *argk) {
	IDL_VPTR result, src, tmp, lngargv[4], plainargv[4];
	IDL_LONG *ptjm, tjm[4];
	float *wigvals;
	double NormVal, *wigtemp, *wigvalsd;
	int i, j, conv, tjmax, tjmin, Nj;
	int n_dim;
	IDL_ARRAY_DIM dim;

	//keyword processing
	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD;
		IDL_LONG dodouble;
	} KW_RESULT;

	static IDL_KW_PAR kw_pars[] = {
		IDL_KW_FAST_SCAN, 
		{ "DOUBLE", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, IDL_KW_OFFSETOF(dodouble)}, 
		{ NULL }
	};
	KW_RESULT kw;
	IDL_KWProcessByOffset(argc,argv,argk,kw_pars,plainargv,1,&kw);

	//make sure we have 4 arguments
	if (argc != 4) {
		i = 1;
		//do something which says it needs 4 inputs
		//possible for 1 input with array of length 4
	}

	//convert to a long if needed
	//values are 2*j and 2*m
	if (argv[0]->type != IDL_TYP_LONG) {
		conv = 1;
		src = IDL_CvtLng(1,&argv[0]);
		lngargv[0] = src;
		src = IDL_CvtLng(1,&argv[1]);
		lngargv[1] = src;
		src = IDL_CvtLng(1,&argv[2]);
		lngargv[2] = src;
		src = IDL_CvtLng(1,&argv[3]);
		lngargv[3] = src;
	} else {
		conv = 0;
		lngargv[0] = argv[0];
		lngargv[1] = argv[1];
		lngargv[2] = argv[2];
		lngargv[3] = argv[3];
	}

	//put tj and tm values in tjm array
	//can be an array or a scalar
	//if an array, we only use the first value
	for (i=0;i<4;i++) {
		if (argv[i]->flags & IDL_V_ARR) {
			//then it is a 1-d array
			ptjm = (IDL_LONG *) lngargv[i]->value.arr->data;
			tjm[i] = *ptjm;
		} else {
			//scalar
			tjm[i] = lngargv[i]->value.l;
		}
	}
	
	//calculate size of array and make it (double or float depending on keywords
	tjmax = tjm[0]+tjm[1];
	tjmin = (tjm[0] > tjm[1]) ? tjm[0]-tjm[1] : tjm[1]-tjm[0];
	Nj = 1 + ((tjmax-tjmin) / 2);
	n_dim = 1;
	dim[0] = Nj;
	if (kw.dodouble) {
		wigvalsd = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,n_dim,dim,IDL_ARR_INI_NOP, &result);
	} else {
		wigvals = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT,n_dim,dim,IDL_ARR_INI_NOP, &result);
	}

	//wigvals = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT,n_dim,dim,IDL_ARR_INI_NOP, &result);

	//get results and normalize them
	wigtemp = Wigner3jVect(tjm[0],tjm[1],tjm[2],tjm[3],&NormVal);
	for (i=0;i<Nj;i++) {
		if (kw.dodouble) {
			wigvalsd[i] = (NormVal*wigtemp[i]);
		} else {
			wigvals[i] = (float) (NormVal*wigtemp[i]);
		}
	}

	//Delete temporary variables if needed
	if (conv == 1) {
		IDL_Deltmp(lngargv[0]);
		IDL_Deltmp(lngargv[1]);
		IDL_Deltmp(lngargv[2]);
		IDL_Deltmp(lngargv[3]);
	}

	free(wigtemp);
	return result;
}

// If no keyword double, then covert returned array to float
// Need to work on just scalar inputs, currently only works if it is an array
IDL_VPTR wig3j6j(int argc, IDL_VPTR argv[], char *argk, int type) {
	IDL_VPTR result, src, tmp, lngargv[6], plainargv[6];
	IDL_LONG *ptjm[6], jm[6], jm2[6]; //array of pointers to data elements, and holder for input for Wigner3j
	float *wigvals;
	double *wigvalsd;
	int i, j, nelts, temp, temp2, conv, scals[6], numscals, IDL_RET_TYPE;
	FILE *ttt;

	//keyword processing
	typedef struct {
		IDL_KW_RESULT_FIRST_FIELD;
		IDL_LONG dodouble;
	} KW_RESULT;

	static IDL_KW_PAR kw_pars[] = {
		IDL_KW_FAST_SCAN, 
		{ "DOUBLE", IDL_TYP_LONG, 1, IDL_KW_ZERO, 0, IDL_KW_OFFSETOF(dodouble)}, 
		{ NULL }
	};
	KW_RESULT kw;
	IDL_KWProcessByOffset(argc,argv,argk,kw_pars,plainargv,1,&kw);

	//convert inputs to type LONG 32-bit (4-byte) integer
	if (argv[0]->type != IDL_TYP_LONG) {
		conv = 1;
		src = IDL_CvtLng(1,&argv[0]);
		lngargv[0] = src;
		src = IDL_CvtLng(1,&argv[1]);
		lngargv[1] = src;
		src = IDL_CvtLng(1,&argv[2]);
		lngargv[2] = src;
		src = IDL_CvtLng(1,&argv[3]);
		lngargv[3] = src;
		src = IDL_CvtLng(1,&argv[4]);
		lngargv[4] = src;
		src = IDL_CvtLng(1,&argv[5]);
		lngargv[5] = src;
	} else {
		conv = 0;
		lngargv[0] = argv[0];
		lngargv[1] = argv[1];
		lngargv[2] = argv[2];
		lngargv[3] = argv[3];
		lngargv[4] = argv[4];
		lngargv[5] = argv[5];
	}

	//check if any are scalars and put the input arrays into correct arrays and 
	//scalars into correct scalars
	numscals = 0;
	for (i=0;i<6;i++) {
		if (lngargv[i]->flags & IDL_V_ARR) {
			//if it is an array
			ptjm[i] = (IDL_LONG *) lngargv[i]->value.arr->data;
			scals[i]=0;
		} else {
			//if it is a scalar
			jm2[i] = lngargv[i]->value.l;
			scals[i]=1;
			numscals++;
		}
	}

	if (numscals == 6) {
		//then there are no arrays and can output a scalar
		//Everything is a scalar
		jm[0]=jm2[0]; jm[1]=jm2[1]; jm[2]=jm2[2]; jm[3]=jm2[3]; jm[4]=jm2[4]; jm[5]=jm2[5];
		if (kw.dodouble) {
			wigvalsd = (double *) malloc(sizeof(double));
			switch (type) {
				case 0:
					wigvalsd[0] = (double) Wigner3j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 1:
					wigvalsd[0] = (double) ClebschGordon(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 2:
					wigvalsd[0] = (double) RacahV(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 3:
					wigvalsd[0] = (double) Wigner6j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 4:
					wigvalsd[0] = (double) RacahW(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				default:
					wigvalsd[0] = (double) Wigner3j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
			}
		} else {
			wigvals = (float *) malloc(sizeof(float));
			switch (type) {
				case 0:
					wigvals[0] = (float) Wigner3j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 1:
					wigvals[0] = (float) ClebschGordon(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 2:
					wigvals[0] = (float) RacahV(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 3:
					wigvals[0] = (float) Wigner6j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 4:
					wigvals[0] = (float) RacahW(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				default:
					wigvals[0] = (float) Wigner3j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
			}
		}
		result = IDL_GettmpLong(0); //create a temporary long, since there doesn't seem to be one for float
		if (kw.dodouble) {
			result->type = IDL_TYP_DOUBLE;
			result->value.d = wigvalsd[0];
		} else {
			result->type = IDL_TYP_FLOAT;
			result->value.f = wigvals[0];
		}
		if (conv == 1) {
			IDL_Deltmp(lngargv[0]);
			IDL_Deltmp(lngargv[1]);
			IDL_Deltmp(lngargv[2]);
			IDL_Deltmp(lngargv[3]);
			IDL_Deltmp(lngargv[4]);
			IDL_Deltmp(lngargv[5]);
		}
		return result;
	}
	
	//dimensions of output array are that of one with smallest number
	//of elements, ignore scalars or 1-element arrays
	nelts = -1;
	temp2 = -1;
	for (i=0;i<6;i++) {
		if (scals[i] == 0) {
			//only check for arrays
			temp = lngargv[i]->value.arr->n_elts;
			if ((temp < nelts) || (nelts < 0)) {
				//if this input has less elements, then use this as the output dimensions
				if (temp == 1) {
					//effectively a scalar, but output a 1-element array
					scals[i] = 2;
					jm2[i] = ptjm[i][0]; //stick value where it is needed
					temp2 = i; //so I can easily access this is everything is a 1-element array
				} else {
					nelts = temp;
					temp2 = i;
				}
			}
		}
	}

	//if nelts < 0 then everything is a 1-element array or a scalar (with at least 
	//one 1-element array. temp2 will point to one of these 1-element arrays
	if (nelts < 0)
		nelts = 1; //for case where everything is scalar

	if (kw.dodouble) {
		wigvalsd = (double *) IDL_MakeTempArray(IDL_TYP_DOUBLE,lngargv[temp2]->value.arr->n_dim,
				lngargv[temp2]->value.arr->dim,IDL_ARR_INI_NOP, &result);
	} else {
		wigvals = (float *) IDL_MakeTempArray(IDL_TYP_FLOAT,lngargv[temp2]->value.arr->n_dim,
				lngargv[temp2]->value.arr->dim,IDL_ARR_INI_NOP, &result);
	}

#pragma omp parallel default(shared) private(jm,i,j)
	{
#pragma omp for schedule(static)
		for (i=0;i<nelts;i++) {
			for (j=0;j<6;j++) {
				if (scals[j] == 0) jm[j] = ptjm[j][i]; else jm[j] = jm2[j]; //stick everything into private data
			}
			switch (type) {
				case 0:
					//wigner 3j
					if (kw.dodouble)
						wigvalsd[i] = (double) Wigner3j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					else
						wigvals[i] = (float) Wigner3j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 1:
					//Clebsch-gordon
					if (kw.dodouble)
						wigvalsd[i] = (double) ClebschGordon(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					else
						wigvals[i] = (float) ClebschGordon(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 2:
					//Racah-V
					if (kw.dodouble)
						wigvalsd[i] = (double) RacahV(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					else
						wigvals[i] = (float) RacahV(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 3:
					//Wigner 6j
					if (kw.dodouble)
						wigvalsd[i] = (double) Wigner6j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					else
						wigvals[i] = (float) Wigner6j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				case 4:
					//Racah-W
					if (kw.dodouble)
						wigvalsd[i] = (double) RacahW(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					else
						wigvals[i] = (float) RacahW(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					break;
				default:
					//default is Wigner3j
					if (kw.dodouble)
						wigvalsd[i] = (double) Wigner3j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
					else
						wigvals[i] = (float) Wigner3j(jm[0],jm[1],jm[2],jm[3],jm[4],jm[5]);
			}
		}
	}

	//if we had to convert to longs, then delete the temporary variables we had to create
	if (conv == 1) {
		IDL_Deltmp(lngargv[0]);
		IDL_Deltmp(lngargv[1]);
		IDL_Deltmp(lngargv[2]);
		IDL_Deltmp(lngargv[3]);
		IDL_Deltmp(lngargv[4]);
		IDL_Deltmp(lngargv[5]);
	}

	IDL_KW_FREE;
	return result;

}

IDL_VPTR wig3j(int argc, IDL_VPTR argv[], char *argk) {
	int type;

	type = 0; //tells wig3j6j to do correct thing
	return wig3j6j(argc, argv, argk, type);
}

IDL_VPTR cgfun(int argc, IDL_VPTR argv[], char *argk) {
	int type;

	type = 1;
	return wig3j6j(argc, argv, argk, type);
}

IDL_VPTR rvfun(int argc, IDL_VPTR argv[], char *argk) {
	int type;

	type = 2;
	return wig3j6j(argc, argv, argk, type);
}

IDL_VPTR wig6j(int argc, IDL_VPTR argv[], char *argk) {
	int type;

	type = 3; //tells wig3j6j to do 6j symbols
	return wig3j6j(argc, argv, argk, type);
}

IDL_VPTR rwfun(int argc, IDL_VPTR argv[], char *argk) {
	int type;

	type = 4; //tells wig3j6j to do correct thing
	return wig3j6j(argc, argv, argk, type);
}

int IDL_Load(void) {

	int retval;
	//will need to modify first 0 when routines start taking keywords
	static IDL_SYSFUN_DEF2 wig3j_addr[] = {
		{wig3j,"WIGNER3J",6,6,IDL_SYSFUN_DEF_F_KEYWORDS,0},};
	static IDL_SYSFUN_DEF2 wig6j_addr[] = {
		{wig6j,"WIGNER6J",6,6,IDL_SYSFUN_DEF_F_KEYWORDS,0},};
	/*static IDL_SYSFUN_DEF2 wig9j_addr[] = {
		{wig9j,"WIGNER9J",9,9,0,0},};*/
	static IDL_SYSFUN_DEF2 cg_addr[] = {
		{cgfun,"CLEBSCHGORDON",6,6,IDL_SYSFUN_DEF_F_KEYWORDS,0},};
	static IDL_SYSFUN_DEF2 rv_addr[] = {
		{rvfun,"RACAHV",6,6,IDL_SYSFUN_DEF_F_KEYWORDS,0},};
	static IDL_SYSFUN_DEF2 rw_addr[] = {
		{rwfun,"RACAHW",6,6,IDL_SYSFUN_DEF_F_KEYWORDS,0},};
	static IDL_SYSFUN_DEF2 wig3jv_addr[] = {
		{wig3jv,"WIGNER3JVECT",4,4,IDL_SYSFUN_DEF_F_KEYWORDS,0},};
	/*static IDL_SYSFUN_DEF2 wig6jv_addr[] = {
		{wig6jv,"WIGNER6JVECT",5,5,0,0},};*/

	//IDL_SysRnAdd returns True if it succeeds in adding the routine, False
	//in the event of an error
	retval = IDL_SysRtnAdd(wig3j_addr, TRUE, IDL_CARRAY_ELTS(wig3j_addr)) && 
		IDL_SysRtnAdd(cg_addr, TRUE, IDL_CARRAY_ELTS(cg_addr)) &&
		IDL_SysRtnAdd(rv_addr, TRUE, IDL_CARRAY_ELTS(rv_addr)) &&
		IDL_SysRtnAdd(wig6j_addr, TRUE, IDL_CARRAY_ELTS(wig6j_addr)) &&
		IDL_SysRtnAdd(rw_addr, TRUE, IDL_CARRAY_ELTS(rw_addr)) &&
		IDL_SysRtnAdd(wig3jv_addr, TRUE, IDL_CARRAY_ELTS(wig3jv_addr)); /* &&
		IDL_SysRtnAdd(wig9j_addr, TRUE, IDL_CARRAY_ELTS(wig9j_addr)) &&
		IDL_SysRtnAdd(wig6jv_addr, TRUE, IDL_CARRAY_ELTS(wig6jv_addr));*/
	
	return retval;
}
