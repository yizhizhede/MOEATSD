#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lmf.h"
#include "matrix.h"
#include "lsmop.h"
#include "landscape.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4
//
// static int 	type_fm = 2;//type for formulation model, 0: addition, 1: multiplication, 2: mixed model
static int 	type_lk = 2;//type of variable linkage, 0: linear linkage, 1: nonlinear linkage, 2: mixed/hybrid
// static int 	type_dg = 1;//type for deep grouping distance-related variables, 0 : even grouping; 1 : nonuniform grouping
// static int 	type_cv = 1;//type for contribution of variables, 0: balanced contribution; 1: unbalanced contribution

static int 	nip[100];	//Number of independent variables in each position-related variable group
static int 	nop[100];	//Number of overlapping variables in each position-related variable group
static int 	nsp;		//Number of shared variables in each position-related variable group
static int 	nid[100];	//Number of independent variables in each distance-related variable group
static int 	nod[100];	//Number of overlapping variables in each distance-related variable group
static int 	nsd;		// Number of shared variables in each distance-related variable group
static int 	K;		// Number of positon varialbes
static int 	L;		// Number of distance variables

//
static int 	gp[10][100];	// groups of position variables
static double 	gp_val[10][100];// groups of position variables
static int 	gp_len;		// Length 
static int 	gp_len_len[10];	// Length 

//
static int 	gd[10][1000000];	// groups of position variables
static int 	gd_len;			// groups of position variables
static int 	gd_len_len[10];		// Length of groups of position variables

//
static int 	dgd[10][1000][1000];	//
static double 	dgd_val[10][1000][1000];	//
static int	dgd_len;
static int	dgd_len_len[10];
static int 	dgd_len_len_len[10][1000];

//
static char 	hType[64]; 	
static char 	gType[64]; 	

//
static void 	getWeightsLogisticMap (double c1, double r, int num, double *w);
static double 	sum_avg_abs (double *x, int n);
static void 	overlapGrouping_gp ();
static void 	overlapGrouping_gd ();
static void 	deepGrouping ();

//
static void 	evaluate (double *x, int n, double *f, int m);
static void 	evalH(double *xI, double *h, int numObj);
static void 	evalG(double *g);

/**********************************************************************************************************************/
/***********  LMF_new () **********************************************************************************************/
/**********************************************************************************************************************/

Problem_t *LMF_new (char *title, int numObj, int numVar) {
	// common variable
	int 	i;
	size_t 	size;
	double* lowBound = NULL; 
	double* uppBound = NULL;
	int	numberOfObjectives_ = numObj;
	int 	numberOfVariables_  = numVar;
	int 	sumNIP, sumNID;
	int 	nsd_leas = 2;	//make sure the number of nsd is at least 2
	int 	proportion = 5;	//control the number of nod as a percentage of the number of nid
	double	w[100];		//

	// 0. allocating memory for a problem
        Problem_t *problem = (Problem_t *)malloc (sizeof (Problem_t));
        if (problem == NULL) {
        	fprintf (stderr, "Allocating memory failed\n");         
               	exit (-1);
        }
	strcpy (problem->title, title);
	problem->numObj = numObj; 
	problem->numVar = numVar;

	// 1. set nip, nop, nsp;
	sumNIP = 0;
	for(i=0; i<numberOfObjectives_-1; i++){
		nip[i] = 3;
		if(numberOfObjectives_ == 2){
			nop[i] = 0;
		}else{
			nop[i] = 1;
		}
		sumNIP += nip[i];
	}
	nsp = 1;
	K = nsp + sumNIP;
	L = numberOfVariables_ - K;

	// 2. set nid, nod, nsd;
	sumNID = 0;
	getWeightsLogisticMap(0.342, 3.8, numObj, w);
	for(i = 0; i < numberOfObjectives_; i++){
		nid[i] = (int)floor(w[i]*(L - nsd_leas));
		sumNID += nid[i];
	}
	for(i = 0; i < numberOfObjectives_; i++){
		if(i == 0){
			nod[i] = nid[numberOfObjectives_-1]/proportion;
		}else{
			nod[i] = nid[i-1]/proportion;
		}
	}
	nsd = L - sumNID;	

	// 3. overlapGrouping
	gp_len = numObj - 1;
	gd_len = numObj;
	overlapGrouping_gp();
	overlapGrouping_gd();

	// 4. deep Grouping
	deepGrouping(); 

	// print grouping information to screem
	printf ("-------------   Test Suite LMF Grouping Information ----------------------\n");
	for (i=0; i<gp_len; i++) {
		printf ("gp[%d]: ", i);
		for (int j=0; j<gp_len_len[i]; j++) {
			printf ("%d ", gp[i][j]);	
		}
		printf ("\n");
		printf ("--------------------------------------------------------------------------\n");
	}

	for (i=0; i<dgd_len; i++) {
		for (int j=0; j<dgd_len_len[i]; j++) {
			printf ("dgd[%d][%d]: ", i, j);
			for (int k=0; k<dgd_len_len_len[i][j]; k++) {
				printf ("%d ", dgd[i][j][k]);	
			}
			printf ("\n");
		}
		printf ("--------------------------------------------------------------------------\n");
	}

	// 5. setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<K; i++) {
		lowBound[i] = -1.0;
		uppBound[i] = 1.0;
	}
	for (i=K; i<numVar; i++) {
		lowBound[i] = 0.0;
		uppBound[i] = 10.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "LMF1"))  	  { strcpy(hType, "concave"); strcpy (gType, "LMF1"); } 
	else if (!strcmp (title, "LMF2")) { strcpy(hType, "convex");  strcpy (gType, "LMF2"); } 
	else if (!strcmp (title, "LMF3")) { strcpy(hType, "linear");  strcpy (gType, "LMF3"); } 
	else if (!strcmp (title, "LMF4")) { strcpy(hType, "concave"); strcpy (gType, "LMF4"); } 
	else if (!strcmp (title, "LMF5")) { strcpy(hType, "convex");  strcpy (gType, "LMF5"); }	
	else if (!strcmp (title, "LMF6")) { strcpy(hType, "linear");  strcpy (gType, "LMF6"); } 
	else if (!strcmp (title, "LMF7")) { strcpy(hType, "concave"); strcpy (gType, "LMF7"); } 
	else if (!strcmp (title, "LMF8")) { strcpy(hType, "linear");  strcpy (gType, "LMF8"); } 
	else if (!strcmp (title, "LMF9")) { strcpy(hType, "inverted_concave"); strcpy (gType, "LMF9"); } 
	else if (!strcmp (title, "LMF10")){ strcpy(hType, "inverted_linear");  strcpy (gType, "LMF10"); } 
	else if (!strcmp (title, "LMF11")){ strcpy(hType, "inverted_concave"); strcpy (gType, "LMF11"); } 
	else if (!strcmp (title, "LMF12")){ strcpy(hType, "inverted_linear");  strcpy (gType, "LMF12"); } 
	else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}
	problem->evaluate = evaluate; 

	return problem;
}

/**********************************************************************************************************************/
/***********  LMF_sample () **********************************************************************************************/
/**********************************************************************************************************************/

Matrix_t *LMF1_sample (int numObj);
Matrix_t *LMF2_sample (int numObj);
Matrix_t *LMF3_sample (int numObj);
Matrix_t *LMF9_sample (int numObj);
Matrix_t *LMF10_sample (int numObj);

Matrix_t *LMF_sample (int No, int numObj) {
	switch (No) {
		case 1:
		case 4:
		case 7:
			return LMF1_sample (numObj);
		case 2:
		case 5:
			return LMF2_sample (numObj);
		case 3:
		case 6:
		case 8:
			return LMF3_sample (numObj);


		case 9:
		case 11:
			return LMF9_sample (numObj);
		case 10:
		case 12:
			return LMF10_sample (numObj);
		default:
		       	fprintf (stderr, "LMF%d have been not on consideration\n", No);
               		exit (0);
	}
	return NULL;
}

Matrix_t *LMF1_sample (int numObj) {
	return LSMOP_sample (5, numObj);
}

Matrix_t *LMF2_sample (int numObj) {
	Matrix_t *sample = Matrix_new ();
	Matrix_t *Grid = NULL;
        int 	i, H, M = numObj;
        size_t 	size;
         
        sample->colDim = M;
	for (H=1; grid(M-1, H) < NUM_SAMPLE; H++){};
        sample->rowDim = grid(M-1, H);
         
        size = sample->rowDim * sample->colDim * sizeof (double);
        sample->elements = (double *)malloc (size);
	
	Grid = Grid_sample (M-1, H);

	for (i=0; i<sample->rowDim; i++) {
		evalH(Grid->elements+i*(M-1), sample->elements+i*M, M);
	}

	Matrix_free (&Grid);
	return sample;
}

Matrix_t *LMF3_sample (int numObj) {
	return LSMOP_sample (1, numObj);
}

Matrix_t *LMF9_sample (int numObj) {
	Matrix_t*	sample = NULL;
	int 		i, j;
	int 		rowDim, colDim;
	
	sample = LSMOP_sample (5, numObj);
	rowDim = sample->rowDim;
	colDim = sample->colDim;
	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			sample->elements[i*colDim+j] = 1.0 - sample->elements[i*colDim+j];
		}
	}
	return sample;
}
Matrix_t *LMF10_sample (int numObj) {
	Matrix_t*	sample = NULL;
	int 		i, j;
	int 		rowDim, colDim;
	
	sample= LSMOP_sample (1, numObj);
	rowDim = sample->rowDim;
	colDim = sample->colDim;
	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			sample->elements[i*colDim+j] = 1.0 - sample->elements[i*colDim+j];
		}
	}
	return sample;
}
/**********************************************************************************************************************/
/**********************************************************************************************************************/
/**********************************************************************************************************************/

static void getWeightsLogisticMap(double c1, double r, int num, double *w) {
	int 	i;
	double 	sum = c1;

	w[0] = c1;
	for(i = 0; i<num-1; i++) {
		w[i+1] = r*w[i]*(1-w[i]);
		sum += w[i+1];
	}
	
	for(i = 0; i<num; i++) {
		w[i] = w[i]/sum;
	}
}

static double sum_avg_abs(double* x, int n) {
	int 	i;
    	double 	avg = 0;

	for(i=0; i<n; i++){
		avg += x[i];
	}
	return fabs(avg/n);
} 

static void overlapGrouping_gp () {
	int 	i, j;
	int 	pointer = 0;
	int 	len, len_len;

	len = gp_len;
	for(i=0; i<len; i++) {		//assign the independent variables to the last part of each group
		gp_len_len[i] = nip[i] + nop[i] + nsp;
		len_len = gp_len_len[i];
		for(j=nop[i]+nsp; j<len_len; j++) {
			gp[i][j] = pointer;
			pointer = pointer+1;
		}
	}
	
	for(i=0; i<len; i++){		//assign the overlap variables to the first part of each group
		for(j=0; j<nop[i]; j++){
			if(i == 0) {
				gp[i][j] = gp[len-1][gp_len_len[len-1]-1-j];
			} else {
				gp[i][j] = gp[i-1][gp_len_len[i-1]-j-1];
			}
		}
	}
	
	for(j=0; j<nsp; j++){		//assign the shared variables to the middle part of each group
		for(i=0; i<len; i++){
			gp[i][j+nop[i]] = pointer;
		}
		pointer = pointer+1; 
	}
}

static void overlapGrouping_gd () {
	int 	i, j;
	int 	pointer = K;
	int 	len, len_len;

	len = gd_len;
	for(i=0; i<len; i++) {		//assign the independent variables to the last part of each group
		gd_len_len[i] = nid[i] + nod[i] + nsd;
		len_len = gd_len_len[i];
		for(j=nod[i]+nsd; j<len_len; j++) {
			gd[i][j] = pointer;
			pointer = pointer+1;
		}
	}
	
	for(i=0; i<len; i++){		//assign the overlap variables to the first part of each group
		for(j=0; j<nod[i]; j++){
			if(i == 0) {
				gd[i][j] = gd[len-1][gd_len_len[len-1]-1-j];
			} else {
				gd[i][j] = gd[i-1][gd_len_len[i-1]-j-1];
			}
		}
	}
	
	for(j=0; j<nsd; j++){		//assign the shared variables to the middle part of each group
		for(i=0; i<len; i++){
			gd[i][j+nod[i]] = pointer;
		}
		pointer = pointer+1; 
	}
}

//The first term of the arithmetic sequence is a, and the common difference is d
static void deepGrouping() {
	int 	i, j, k, a, d, len;
	int 	span, remain, t = 0;

	len = gd_len;
	dgd_len = len;
	for(i=0; i<len; i++) {
		a = 5;
		d = (i%2 == 0) ? 1 : 2;

		// 
	 	span = a;
    	 	remain = gd_len_len[i];
    	 	t = 0;
		k = 0;
		
		while(remain > span + a) {
			for(j=0; j<span; j++) {
				dgd[i][k][j] = gd[i][t];
				t = t + 1;
			}
			dgd_len_len_len[i][k] = span;
		    	remain = remain - span;
		    	span = span + d;
			k++;
		}
	    	if(remain > 0){
			for(j=0; j<remain; j++) {
				dgd[i][k][j] = gd[i][t];
				t = t + 1;
			}
			dgd_len_len_len[i][k] = remain;
			k++;
	    	}
		dgd_len_len[i] = k;
	}
}

//
static void evaluate (double *x, int n, double *f, int m) {
	int 	i, j, k, index;
	int 	len, len_len, len_len_len;
	double 	y[100], sb, h[100], g[100];

	// 1. update gp_val
	len = gp_len;
	for (i = 0; i < len; i++){
		len_len = gp_len_len[i];
		for(j=0; j<len_len; j++){
			gp_val[i][j] = x[gp[i][j]];
		}
	}

	// 2. update y
	len = gp_len;
	for (i = 0; i < len; i++){
		len_len = gp_len_len[i];
		y[i] = sum_avg_abs((double *)gp_val[i], len_len);
	}

	// 3. update dgd_val
	len = dgd_len;
	for (i = 0; i < len; i++){
		len_len = dgd_len_len[i];
		for (j=0; j<len_len; j++){
			len_len_len = dgd_len_len_len[i][j];
			for (k=0; k<len_len_len; k++){
				sb = x[dgd[i][j][k]];
				index = dgd[i][j][k];

				//
				if (type_lk == 0) {
				// 	from LMF code
				//	dgd_val[i][j][k] = (1.0 + (index+1.0)/(L))*sb - 10.0*y[0];

				// 	from LSMOP
					dgd_val[i][j][k] = (1.0 + (index+1.0)/(n))*sb - 10.0*y[0];
				} else if(type_lk == 1) {
				// 	from LMF code
				//	dgd_val[i][j][k] = (1.0 +cos(0.5*PI*((index+1.0)/(L))))*sb - 10.0*y[0];

				// 	from LSMOP
					dgd_val[i][j][k] = (1.0 +cos(0.5*PI*((index+1.0)/(n))))*sb - 10.0*y[0];
				} else if (type_lk == 2) {
					if(i % 2 == 0) {
					//	from LMF code
					//	dgd_val[i][j][k] = (1.0 + (index+1.0)/(L))*sb - 10.0*y[0];

					// 	from LSMOP
						dgd_val[i][j][k] = (1.0 + (index+1.0)/(n))*sb - 10.0*y[0];
					} else {
					//	from LMF code
					//	dgd_val[i][j][k] = (1.0 + cos(0.5*PI*((index+1.0)/(L))))*sb - 10.0*y[0];

					// 	from LSMOP
						dgd_val[i][j][k] = (1.0 + cos(0.5*PI*((index+1.0)/(n))))*sb - 10.0*y[0];
					}
				} else {
					// dgd_val[i][j][k] = sb;
					// dgd_val[i][j][k] = sb - 10 * 0.35;
					// dgd_val[i][j][k] = sb - 10 * (index+1.0)/(n);
					// dgd_val[i][j][k] = sb - 10 * y[0];
					dgd_val[i][j][k] = (1.0 + (index+1.0)/(n))*sb - 10 * y[0];
				}
			}
		}
	}	

	// 4. evalue H and G
	evalH(y, h, m);
	evalG(g);	
	
	// 5. update  f
	for (i = 0; i<m; i++) {
		if(i % 2 == 0) {
			f[i] = h[i]*(1 + g[i]);
		}else {
			f[i] = h[i] + g[i];
		}
	}
}


static void evalH(double *xI, double *h, int numObj) {
	int 	i, j;
	int 	numberOfObjectives_ = numObj;
	int	aux;

	if (strcmp (hType, (char *)"linear") == 0) {
		for (i = 0; i < numberOfObjectives_; i++) {
			h[i] = 1.0;
			for (j = 0; j < numberOfObjectives_ - (i + 1); j++)
				h[i] *= xI[j];
			if (i != 0) {
				aux = numberOfObjectives_ - (i + 1);
				h[i] *= (1 - xI[aux]);
			} // if
		} // for
	}else if(strcmp (hType, (char *)"concave") == 0){
		for (i = 0; i < numberOfObjectives_; i++) {
			h[i] = 1.0;
			for (j = 0; j < numberOfObjectives_ - (i + 1); j++)
				h[i] *= cos(xI[j] * 0.5 * PI);
			if (i != 0) {
				aux = numberOfObjectives_ - (i + 1);
				h[i] *= sin(xI[aux] * 0.5 * PI);
			} // if
		} // for
	}else if(strcmp (hType, (char *)"convex") == 0){
		for (i = 0; i < numberOfObjectives_; i++) {
			h[i] = 1.0;
			for (j = 0; j < numberOfObjectives_ - (i + 1); j++)
				h[i] *= cos(xI[j] * 0.5 * PI);
			if (i != 0) {
				int aux = numberOfObjectives_ - (i + 1);
				h[i] *= sin(xI[aux] * 0.5 * PI);
			} // if
			if(i != numberOfObjectives_-1)
				h[i] = pow(h[i], 4);
			else
				h[i] = pow(h[i], 2);
		} // for
	}if (strcmp (hType, (char *)"inverted_linear") == 0){
		for (i = 0; i < numberOfObjectives_; i++) {
			h[i] = 1.0;
			for (j = 0; j < numberOfObjectives_ - (i + 1); j++)
				h[i] *= xI[j];
			if (i != 0) {
				aux = numberOfObjectives_ - (i + 1);
				h[i] *= (1 - xI[aux]);
			} // if
			h[i] = 1.0 - h[i];
		} // for
	}else if(strcmp (hType, (char *)"inverted_concave") == 0){
		for (i = 0; i < numberOfObjectives_; i++) {
			h[i] = 1.0;
			for (j = 0; j < numberOfObjectives_ - (i + 1); j++)
				h[i] *= cos(xI[j] * 0.5 * PI);
			if (i != 0) {
				aux = numberOfObjectives_ - (i + 1);
				h[i] *= sin(xI[aux] * 0.5 * PI);
			} // if
			h[i] = 1.0 - h[i];
		} // for
	}
}

static void getGLMF1(double *g);
static void getGLMF2(double *g);
static void getGLMF3(double *g);
static void getGLMF4(double *g);
static void getGLMF5(double *g);
static void getGLMF6(double *g);
static void getGLMF7(double *g);
static void getGLMF8(double *g);
static void getGLMF9(double *g);
static void getGLMF10(double *g);
static void getGLMF11(double *g);
static void getGLMF12(double *g);

static void evalG(double *g) {
	if (strcmp (gType, (char *)"LMF1") == 0) 	getGLMF1(g);
	else if (strcmp (gType, (char *)"LMF2") == 0) 	getGLMF2(g);
	else if (strcmp (gType, (char *)"LMF3") == 0)	getGLMF3(g);
	else if (strcmp (gType, (char *)"LMF4") == 0)	getGLMF4(g);
	else if (strcmp (gType, (char *)"LMF5") == 0)	getGLMF5(g);
	else if (strcmp (gType, (char *)"LMF6") == 0)	getGLMF6(g);
	else if (strcmp (gType, (char *)"LMF7") == 0)	getGLMF7(g);
	else if (strcmp (gType, (char *)"LMF8") == 0)	getGLMF8(g);
	else if (strcmp (gType, (char *)"LMF9") == 0)	getGLMF9(g);
	else if (strcmp (gType, (char *)"LMF10") == 0)	getGLMF10(g);
	else if (strcmp (gType, (char *)"LMF11") == 0)	getGLMF11(g);
	else if (strcmp (gType, (char *)"LMF12") == 0) 	getGLMF12(g);
	else {
		printf ("Error: g function type %s invalid\n", gType);
		exit(0);
	}
}

//
static double getSphere(double *x, int len);
static double getSchwefel1(double *x, int len);
static double getSchwefel2(double *x, int len);
static double getAckley(double *x, int len);
static double getRastrigin(double *x, int len);
static double getRosenbrock(double *x, int len);


/*LMF1: GFunction used in LMF1*/
static void getGLMF1(double *g) {
	int 	i, j;
	double	w1[100], w2[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if (i%2 == 0) {
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for (j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if (j % 2 == 0) {
					g[i] += w1[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				} else {
					g[i] += w1[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}
			}
		} else {
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for (j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w2[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}
	
/*LMF2: GFunction used in LMF2*/
static void getGLMF2 (double *g) {
	int 	i, j;
	double	w1[100], w2[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if (i%2 == 0){
			getWeightsLogisticMap (0.23, 3.7, len_len, w1);
			for (j=0; j<len_len; j++) {
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w1[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		} else {
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w2[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}
	
/*LMF3: GFunction used in LMF3*/
static void getGLMF3(double *g) {
	int 	i, j;
	double	w1[100], w2[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if(i%2 == 0){
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w1[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		} else {
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w2[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getRastrigin((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}
	
/*LMF4: GFunction used in LMF4*/
static void getGLMF4(double *g) {
	int 	i, j;
	double	w1[100], w2[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if (i%2 == 0) {
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w1[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}
			}
		} else {
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w2[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getRosenbrock((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}
	
/*LMF5: GFunction used in LMF5*/
static void getGLMF5(double *g) {
	int 	i, j;
	double	w1[100], w2[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if (i%2 == 0) {
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 3 == 0) {
					g[i] += w1[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else if(j % 3 == 1) {
					g[i] += w1[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		}else{
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 3 == 0) {
					g[i] += w2[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else if(j % 3 == 1) {
					g[i] += w2[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getRastrigin((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}
	
/*LMF6: GFunction used in LMF6*/
static void getGLMF6(double* g) {
	int 	i, j;
	double	w1[100], w2[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if(i%2 == 0){
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len; j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 3 == 0) {
					g[i] += w1[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else if(j % 3 == 1) {
					g[i] += w1[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		}else{
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 3 == 0) {
					g[i] += w2[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else if(j % 3 == 1) {
					g[i] += w2[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getRastrigin((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}
	
/*LMF7: GFunction used in LMF7*/
static void getGLMF7(double *g) {
	int 	i, j;
	double	w1[100], w2[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if(i%2 == 0){
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 3 == 0) {
					g[i] += w1[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else if(j % 3 == 1) {
					g[i] += w1[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getRosenbrock((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		}else{
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 3 == 0) {
					g[i] += w2[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else if(j % 3 == 1) {
					g[i] += w2[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getRosenbrock((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}

	
/*LMF8: GFunction used in LMF8*/
static void getGLMF8(double *g) {
	int 	i, j;
	double	w1[100], w2[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if(i%2 == 0){
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 3 == 0) {
					g[i] += w1[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else if(j % 3 == 1) {
					g[i] += w1[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getRosenbrock((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		}else{
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 3 == 0) {
					g[i] += w2[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}else if(j % 3 == 1) {
					g[i] += w2[j] * getRastrigin((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getRosenbrock((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}
	
/*LMF9: GFunction used in LMF9*/
static void getGLMF9(double *g) {
	int 	i, j;
	double	w1[100], w2[100], w3[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if(i%3 == 0){
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w1[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		}else if(i%3 == 1){
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w2[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}
			}
		}else {
			getWeightsLogisticMap(0.23, 3.7, len_len, w3);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w3[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w3[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}
	
/*LMF10: GFunction used in LMF10*/
static void getGLMF10(double *g) {
	int 	i, j;
	double	w1[100], w2[100], w3[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if(i%3 == 0){
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w1[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		}else if(i%3 == 1){
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w2[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getRastrigin((double *)dgd_val[i][j], len_len_len);
				}
			}
		}else {
			getWeightsLogisticMap(0.23, 3.7, len_len, w3);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w3[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w3[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}

/*LMF11: GFunction used in LMF11*/
static void getGLMF11(double *g) {
	int 	i, j;
	double	w1[100], w2[100], w3[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if(i%3 == 0){
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w1[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getRastrigin((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		}else if(i%3 == 1){
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w2[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}
			}
		}else {
			getWeightsLogisticMap(0.23, 3.7, len_len, w3);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w3[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w3[j] * getRastrigin((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}

/*LMF12: GFunction used in LMF12*/
static void getGLMF12(double *g) {
	int 	i, j;
	double	w1[100], w2[100], w3[100];
	int 	len, len_len, len_len_len;

	len = dgd_len;
	for (i = 0; i < len; i++) {
		g[i] = 0.0;
		len_len = dgd_len_len[i];
		if(i%3 == 0){
			getWeightsLogisticMap(0.23, 3.7, len_len, w1);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w1[j] * getSphere((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w1[j] * getRosenbrock((double *)dgd_val[i][j], len_len_len);
				}
				
			}
		}else if(i%3 == 1){
			getWeightsLogisticMap(0.23, 3.75, len_len, w2);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w2[j] * getSchwefel2((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w2[j] * getRastrigin((double *)dgd_val[i][j], len_len_len);
				}
			}
		}else {
			getWeightsLogisticMap(0.23, 3.7, len_len, w3);
			for(j=0; j<len_len;j++){
				len_len_len = dgd_len_len_len[i][j];
				if(j % 2 == 0) {
					g[i] += w3[j] * getSchwefel1((double *)dgd_val[i][j], len_len_len);
				}else {
					g[i] += w3[j] * getAckley((double *)dgd_val[i][j], len_len_len);
				}
			}
		}
	}
}

/*bf1: Sphere Function*/
static double getSphere(double *x, int len) {
	double 	sum = 0.0;
	int 	i;

	for (i = 0; i < len; i++)
		sum += x[i]*x[i];		
	return sum/len;
}//Unimodal and Separable

/*bf2: Schwefel 1.2 Function*/
static double getSchwefel1(double *x, int len) {
	double 	sum = 0.0;
	double 	mid = 0.0;
	int 	i, j;

	for (i = 0; i < len; i++) {
		mid = 0.0;
		for(j = 0; j < i; j++) {
			mid += x[j];
		}
		sum += mid*mid;
	}			
	return sum/len;
}//Unimodal and Non-Separable

/*bf3: Schwefel 2.21 Function*/
static double getSchwefel2(double *x, int len) {
	double 	max = fabs(x[0]);
	int 	i;

	for (i = 1; i < len; i++) {
		if(max < fabs(x[i])) {
			max = fabs(x[i]);
		}
	}	
	return max/len;
}//Unimodal and Non-Separable

/*bf4: Ackley Function*/
static double getAckley(double *x, int len) {
	double 	sum1 = 0;
	double 	sum2 = 0;
	int 	i;

	for (i = 0; i < len; i++) {
		sum1 += ((x[i] * x[i]) / len);
		sum2 += (cos(2 * PI * x[i]) / len);
	}
	return (-20 * exp(-0.2 * sqrt(sum1)) - exp(sum2) + 20 + exp(1.0))/len;
}//Multimodal and Separable	

/*bf5: Rastrigin Function*/
static double getRastrigin(double *x, int len) {
	double 	result = 0.0;
	double 	a = 10.0;
	double 	w = 2 * PI;
	int 	i;

	for (i = 0; i < len; i++) {
		result += x[i] * x[i] - a * cos(w * x[i]);
	}
	result += a * len;
	return result/len;
}//Multimodal and Separable	

/*bf6: Rosenbrock Function*/
static double getRosenbrock(double *x, int len) {
	double 	sum = 0, t;
	int 	i;

	for (i = 0; i < len - 1; i++) {
		// original formula
		// t = 100 * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1]) + (1 - x[i]) * (1 - x[i]);

		// revised formula
		t = 100 * (x[i] * x[i] - x[i + 1]) * (x[i] * x[i] - x[i + 1]) + (x[i]) * (x[i]);
		sum += t;
	}
	return sum/len;
}//Multimodal and Non-Separable	
