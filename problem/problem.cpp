#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "algebra.h"
#include "problem.h"
#include "mystring.h"
#include "lsmop.h"
#include "shape.h"
#include "wfg.h"
#include "dtlz.h"
#include "zdt.h"
#include "mop.h"
#include "bk.h"
#include "dpam.h"
#include "dgo.h"
#include "fa.h"
#include "far.h"
#include "fes.h"
#include "ff.h"
#include "ikk.h"
#include "im.h"
#include "jos.h"
#include "kur.h"
#include "lrs.h"
#include "ltdz.h"
#include "le.h"
#include "mhhm.h"
#include "mlf.h"
#include "qv.h"
#include "sch.h"
#include "sp.h"
#include "ssfyy.h"
#include "sk.h"
#include "tkly.h"
#include "vu.h"
#include "vfm.h"
#include "zlt.h"
#include "parameter.h"
#include "uf.h"
#include "smop.h"
#include "bt.h"
#include "lmf.h"
#include "nntp.h"

static Problem_t *P;

Problem_t *Problem_new () { 
	Parameter_t*	parameter = NULL;
	char*		dup = NULL;
	char*		title = NULL; 
	int 		numObj; 
	int 		numVar; 
	int		i;

	// 1. get the parameter objective
	if (NULL == (parameter = Parameter_get ())) {
		fprintf (stderr, "ERROE:%s:%d: The parameter objective is NULL\n", __FILE__, __LINE__);
		exit (0);
	} 

	title  = parameter->problem; 
	numObj = parameter->numObj; 
 	numVar = parameter->numVar; 

	dup = strdup (title);
	i=strlen (dup) - 1;
	while (i >=0 && dup[i] >= '0' && dup[i] <='9') {
		dup[i] = '\0';
		i--;
	}

	if (!strcmp (dup, "LSMOP")) 	P = LSMOP_new (title, numObj, numVar);
	else if (!strcmp (dup, "WFG"))  P = WFG_new (title, numObj, numVar);
	else if (!strcmp (dup, "DTLZ")) P = DTLZ_new (title, numObj, numVar);
	else if (!strcmp (dup, "ZDT")) 	P = ZDT_new (title, numObj, numVar);
	else if (!strcmp (dup, "MOP")) 	P = MOP_new (title, numObj, numVar);
	else if (!strcmp (dup, "BK")) 	P = BK_new (title, numObj, numVar);
	else if (!strcmp (dup, "DPAM"))	P = DPAM_new (title, numObj, numVar);
	else if (!strcmp (dup, "DGO"))	P = DGO_new (title, numObj, numVar);
	else if (!strcmp (dup, "FA"))	P = FA_new (title, numObj, numVar);
	else if (!strcmp (dup, "FAR"))	P = FAR_new (title, numObj, numVar);
	else if (!strcmp (dup, "FES"))	P = FES_new (title, numObj, numVar);
	else if (!strcmp (dup, "FF"))	P = FF_new (title, numObj, numVar);
	else if (!strcmp (dup, "IKK"))	P = IKK_new (title, numObj, numVar);
	else if (!strcmp (dup, "IM"))	P = IM_new (title, numObj, numVar);
	else if (!strcmp (dup, "JOS"))	P = JOS_new (title, numObj, numVar);
	else if (!strcmp (dup, "KUR"))	P = KUR_new (title, numObj, numVar);
	else if (!strcmp (dup, "LRS"))	P = LRS_new (title, numObj, numVar);
	else if (!strcmp (dup, "LTDZ"))	P = LTDZ_new (title, numObj, numVar);
	else if (!strcmp (dup, "LE"))	P = LE_new (title, numObj, numVar);
	else if (!strcmp (dup, "MHHM"))	P = MHHM_new (title, numObj, numVar);
	else if (!strcmp (dup, "MLF"))	P = MLF_new (title, numObj, numVar);
	else if (!strcmp (dup, "QV"))	P = QV_new (title, numObj, numVar);
	else if (!strcmp (dup, "SCH"))	P = SCH_new (title, numObj, numVar);
	else if (!strcmp (dup, "SP"))	P = SP_new (title, numObj, numVar);
	else if (!strcmp (dup, "SSFYY"))P = SSFYY_new (title, numObj, numVar);
	else if (!strcmp (dup, "SK"))   P = SK_new (title, numObj, numVar);
	else if (!strcmp (dup, "TKLY")) P = TKLY_new (title, numObj, numVar);
	else if (!strcmp (dup, "VU"))   P = VU_new (title, numObj, numVar);
	else if (!strcmp (dup, "VFM"))  P = VFM_new (title, numObj, numVar);
	else if (!strcmp (dup, "ZLT"))  P = ZLT_new (title, numObj, numVar);
	else if (!strcmp (dup, "UF"))   P = UF_new (title, numObj, numVar);
	else if (!strcmp (dup, "SMOP")) P = SMOP_new (title, numObj, numVar);
	else if (!strcmp (dup, "BT"))   P = BT_new (title, numObj, numVar);
	else if (!strcmp (dup, "LMF"))  P = LMF_new (title, numObj, numVar);
	else if (!strcmp (dup, "NNTP")) P = NNTP_new (title, numObj, numVar);
	else {
		fprintf (stderr, "Problem:%s have not been implemented.", title);
		exit (0);
	}

	free (dup);

	P->fitness  = 0.0;			// record the number of fitness
	P->lifetime = parameter->lifetime;	// terminal condition
	P->ideal_point = (double *)malloc ((numObj+10)*sizeof (double));	// ideal point
	for (i=0; i<numObj; i++) { P->ideal_point[i] = 1.0e+100; }

	return P;
}


Problem_t *Problem_get () { 
	if (P == NULL) {
		fprintf (stderr, "The routine Problem_new should been invoked before Problem_get\n");
		exit (0);
	}
	return P;
}

Matrix_t *Problem_sample (char *title, int numObj) {
	char*	dup = NULL;
	int 	No = 0;	
	int 	i;

	dup = strdup (title);
	i=strlen (dup) - 1;
	while (i >=0 && dup[i] >='0' && dup[i] <='9') {
		dup[i] = '\0';
		i--;
	}
	sscanf (title+i+1, "%d", &No);

	if (!strcmp (dup, "LSMOP")) 	return LSMOP_sample (No, numObj);
	if (!strcmp (dup, "WFG"))	return WFG_sample (No, numObj);
	if (!strcmp (dup, "DTLZ")) 	return DTLZ_sample (No, numObj);
	if (!strcmp (dup, "ZDT")) 	return ZDT_sample (No, numObj);
	if (!strcmp (dup, "MOP")) 	return MOP_sample (No, numObj);
	if (!strcmp (dup, "BK")) 	return BK_sample (No, numObj);
	if (!strcmp (dup, "DPAM")) 	return DPAM_sample (No, numObj);
	if (!strcmp (dup, "DGO")) 	return DGO_sample (No, numObj);
	if (!strcmp (dup, "FA")) 	return FA_sample (No, numObj);
	if (!strcmp (dup, "FAR")) 	return FAR_sample (No, numObj);
	if (!strcmp (dup, "FES")) 	return FES_sample (No, numObj);
	if (!strcmp (dup, "FF")) 	return FF_sample (No, numObj);
	if (!strcmp (dup, "IKK")) 	return IKK_sample (No, numObj);
	if (!strcmp (dup, "IM")) 	return IM_sample (No, numObj);
	if (!strcmp (dup, "JOS")) 	return JOS_sample (No, numObj);
	if (!strcmp (dup, "KUR")) 	return KUR_sample (No, numObj);
	if (!strcmp (dup, "LRS")) 	return LRS_sample (No, numObj);
	if (!strcmp (dup, "LTDZ")) 	return LTDZ_sample (No, numObj);
	if (!strcmp (dup, "LE")) 	return LE_sample (No, numObj);
	if (!strcmp (dup, "MHHM")) 	return MHHM_sample (No, numObj);
	if (!strcmp (dup, "MLF")) 	return MLF_sample (No, numObj);
	if (!strcmp (dup, "QV")) 	return QV_sample (No, numObj);
	if (!strcmp (dup, "SCH")) 	return SCH_sample (No, numObj);
	if (!strcmp (dup, "SP")) 	return SP_sample (No, numObj);
	if (!strcmp (dup, "SSFYY")) 	return SSFYY_sample (No, numObj);
	if (!strcmp (dup, "SK")) 	return SK_sample (No, numObj);
	if (!strcmp (dup, "TKLY")) 	return TKLY_sample (No, numObj);
	if (!strcmp (dup, "VU")) 	return VU_sample (No, numObj);
	if (!strcmp (dup, "VFM")) 	return VFM_sample (No, numObj);
	if (!strcmp (dup, "ZLT")) 	return ZLT_sample (No, numObj);
	if (!strcmp (dup, "UF")) 	return UF_sample (No, numObj);
	if (!strcmp (dup, "SMOP")) 	return SMOP_sample (No, numObj);
	if (!strcmp (dup, "BT")) 	return BT_sample (No, numObj);
	if (!strcmp (dup, "LMF")) 	return LMF_sample (No, numObj);
	if (!strcmp (dup, "NNTP")) 	return NNTP_sample (No, numObj);

	fprintf (stderr, "The Sample of Problem:%s have not been completed\n", title);
	exit (0);
	return NULL;
}

void Problem_evaluate (double *var, int numVar, double *obj, int numObj) {
	int 	i;
	P->evaluate (var, numVar, obj, numObj);	
	P->fitness += 1;
	for (i=0; i<numObj; i++) if (obj[i] < P->ideal_point[i]){
		P->ideal_point[i] = obj[i]; 
	}
}

static double bar_history = -1;
void Problem_progressbar () {
	double 	t;
	int 	i;

	t = 100.0*P->fitness/P->lifetime;
	if (0.01 > t - bar_history) {
		return;
	} else {
		bar_history = t;
	}
	
	printf ("\r[PROGRESS BAR (%lld/%lld)]", P->fitness, P->lifetime);
	for (i=0; i<t; i+=2) {
		printf ("-");
	}
	printf ("-> %.2f%%", t);
	for (;i<100; i+=2) {
		printf ("|");
	}
	fflush (stdout);
}

static double percentage=0;	
int Problem_isTick () {
	if (P->fitness >= percentage*P->lifetime && percentage < 1.0) {
		percentage += 0.01;
		return 1;
	} else {
		return 0;
	}
}

int Problem_isEnd () {
	long long int d;

	d = P->fitness - P->lifetime;
	if (d >= 0 ) 	
		return 1;
	else
		return 0;
}

long long  Problem_getFitness () {
	return P->fitness;
}

long long  Problem_getLifetime () {
	return P->lifetime;
}

double* Problem_getLowerBound () {
	return P->lowBound;
}

double* Problem_getUpperBound () {
	return P->uppBound;
}

double* Problem_getIdealPoint () {
	return P->ideal_point;
}

void Problem_xToOne (Matrix_t* X) {
	int 	i, j;
	int 	rowDim = X->rowDim;
	int 	colDim = X->colDim;
	double 	value;
	
	if (X->colDim != P->numVar) {
		fprintf (stderr, "ERROE:%s:%d: The size of X does't match %d != %d\n", __FILE__, __LINE__, colDim, P->numVar);
		exit (0);
	}
	
	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			value = (X->elements[i*colDim+j] - P->lowBound[j]) / (P->uppBound[j] - P->lowBound[j]);
			X->elements[i*colDim+j] = value;
		}
	}
}

void Problem_xFromOne (Matrix_t* X) {
	int 	i, j;
	int 	rowDim = X->rowDim;
	int 	colDim = X->colDim;
	double 	value;
	
	if (X->colDim != P->numVar) {
		fprintf (stderr, "ERROE:%s:%d: The size of X does't match %d != %d\n", __FILE__, __LINE__, colDim, P->numVar);
		exit (0);
	}
	
	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			value = X->elements[i*colDim+j] * (P->uppBound[j] - P->lowBound[j]) + P->lowBound[j];
			X->elements[i*colDim+j] = value;
		}
	}
}

void Problem_xToOne (double* X, int numVar) {
	int 	i;
	double 	value;
	
	if (numVar != P->numVar) {
		fprintf (stderr, "ERROE:%s:%d: The size of X does't match %d != %d\n", __FILE__, __LINE__, numVar, P->numVar);
		exit (0);
	}
	
	for (i=0; i<numVar; i++) {
		value = (X[i] - P->lowBound[i]) / (P->uppBound[i] - P->lowBound[i]);
		X[i] = value;
	}
}

void Problem_xFromOne (double*  X, int numVar) {
	int 	i;
	double 	value;
	
	if (numVar != P->numVar) {
		fprintf (stderr, "ERROE:%s:%d: The size of X does't match %d != %d\n", __FILE__, __LINE__, numVar, P->numVar);
		exit (0);
	}
	
	for (i=0; i<numVar; i++) {
		value = X[i] * (P->uppBound[i] - P->lowBound[i]) + P->lowBound[i];
		X[i] = value;
	}
}
