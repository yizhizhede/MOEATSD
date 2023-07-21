#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lsmop.h"
#include "matrix.h"
#include "landscape.h"
#include "shape.h"
#include "linkage.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

/*
 R: 	group size of variables
 I: 	group index of variables
 ns: 	grouping paramter
 nf: 	grouping paramter
 nk: 	grouping paramter
 C: 	correlation matrix
 L: 	linkage funciton pointer
 gI:  	landscape funciton for odd index
 gII: 	landscape function for even index
 H:	shape funciton
 P:	test problem
*/

// Non-uniform grouping of decision variables
static int R[1000];
static int I[1000];
static int ns; 
static int nf;
static int nk;

// Correlation Matrix C
static int *C, *C1, *C2, *C3;

// Variable linkage function pointer
static void (*L) (double *x, int n, Problem_t *problem);

// landscape funciton pointer
static double (*gI)  (double *x, int n);
static double (*gII) (double *x, int n);

// shape function pointer
static void (*H) (double *x, double *g, double *f, int M);

// test problem
static Problem_t *P;

static void grouping (int ns, int nk, int numObj);
static void evaluate (double *x, int n, double *f, int m);

Problem_t *LSMOP_new (char *title, int numObj, int numVar) {
	// initial fixed parameters;
	nf = numObj - 1;
	ns = numVar -nf;	//	ns = 100 * numObj;
	nk = 5;

	// common variable
	int 	i, j, k;
	size_t 	size;
	double 	*lowBound = NULL, *uppBound = NULL;

	// allocating memory for a problem
        Problem_t *problem = (Problem_t *)malloc (sizeof (Problem_t));
        if (problem == NULL) {
        	fprintf (stderr, "Allocating memory failed\n");         
               	exit (-1);
        }
	strcpy (problem->title, title);
	problem->numObj = numObj; 
	problem->numVar = numVar;

	// non-uniform grouping of decision varibles: set I, R
	grouping (ns, nk, numObj);
	
	// setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<nf; i++) {
		lowBound[i] = 0.0;
		uppBound[i] = 1.0;
	}
	for (; i<numVar; i++) {
		lowBound[i] = 0.0;
		uppBound[i] = 10.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// Correlation Matrix 
	size = numObj * numObj * sizeof (int);
	C1 = (int *)malloc (size);
	C2 = (int *)malloc (size);
	C3 = (int *)malloc (size);
	for (i=0, k=0; i<numObj; i++) {
		for (j=0; j<numObj; j++, k++) {
			C1[k] = C2[k] = C3[k] = 0;
			if (i==j) {
				C1[k] = 1;
				C2[k] = 1;
			}
			if (i+1 == j) {
				C2[k] = 1;
			}
			C3[k] = 1;
		}
	}
	
	// landscape funcition
	double (*eta1)(double*,int) = Sphere;
  	double (*eta2)(double*,int) = Schwefel;
  	double (*eta3)(double*,int) = Rosenbrock;
  	double (*eta4)(double*,int) = Rastrigin;
	double (*eta5)(double*,int) = Griewank;
	double (*eta6)(double*,int) = Ackley;

	// 
	if (!strcmp (title, "LSMOP1")) {
		P = problem; L = L1; H = H1; gI = eta1; gII = eta1; C = C1;
	} else if (!strcmp (title, "LSMOP2")) {
		P = problem; L = L1; H = H1; gI = eta5; gII = eta2; C = C1;
	} else if (!strcmp (title, "LSMOP3")) {
		P = problem; L = L1; H = H1; gI = eta4; gII = eta3; C = C1;
	} else if (!strcmp (title, "LSMOP4")) {
		P = problem; L = L1; H = H1; gI = eta6; gII = eta5; C = C1;
	} else if (!strcmp (title, "LSMOP5")) {
		P = problem; L = L2; H = H2; gI = eta1; gII = eta1; C = C2;
	} else if (!strcmp (title, "LSMOP6")) {
		P = problem; L = L2; H = H2; gI = eta3; gII = eta2; C = C2;
	} else if (!strcmp (title, "LSMOP7")) {
		P = problem; L = L2; H = H2; gI = eta6; gII = eta3; C = C2;
	} else if (!strcmp (title, "LSMOP8")) {
		P = problem; L = L2; H = H2; gI = eta5; gII = eta1; C = C2;
	} else if (!strcmp (title, "LSMOP9")) {
		P = problem; L = L2; H = H3; gI = eta1; gII = eta6; C = C3;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	problem->evaluate = evaluate;
	return problem;
}

static void grouping (int ns, int nk, int numObj) {
	// initial fixed parameters
	double 	alpha = 3.8, c0 = 0.1;	

	// ordinate variable
	double 	sum = 0;
	double 	sequence[numObj+1];
	int 	i;

	sequence[0] = alpha * c0 * (1.0 - c0);
	sum += sequence[0];
	for (i=1; i<numObj; i++) {
		sequence[i] = alpha * sequence[i-1] * (1.0 - sequence[i-1]);
		sum += sequence[i];
	}

	for (i=0; i<numObj; i++) {
		sequence[i] = sequence[i] / sum;

		// floor() is used here base on PlatEMO, while the ceil() is used in paper.
		R[i] = floor(sequence[i]*ns/nk)*nk;	
	}

	I[0] = numObj -1;
	for (i=1; i<numObj; i++) {
		I[i] = I[i-1] + R[i-1];
	}
}

static void evaluate (double *x, int n, double *f, int m) {
	double 	g[m], g_bar[m*m];
	double 	h[m];
	int 	i, j, k, offset, len;
	double 	xx[n];

	// copy x to xx
	memcpy (xx, x, n*sizeof (double));

	// linkage
	L (xx, n, P);

	// landscape: g_bar
	for (i=0; i<m; i++) { 	 
		for (j=0; j<m; j++) {
			g_bar[i*m+j] = 0;
			offset = I[j];
			len = R[j] / nk;
			for (k=0; k<nk; k++, offset += len) {
				g_bar[i*m+j] += ((i&1) == 0) ? (gI(xx+offset, len)/R[j]) : (gII(xx+offset, len)/R[j]);
			}
		}
	} 

	// g
	for (i=0; i<m; i++) { 	 
		g[i] = 0;
		for (j=0; j<m; j++) {  
			// base on paper
			// g[i] += C[i*m+j] * g_bar[i*m+j];

			// base on PlatEMO
			g[i] += C[i*m+j] * g_bar[j*m+j];
		}
	} 

	// shape
	H (xx, g, h, m);	

	// fitness
	for (i=0; i<m; i++) {
		f[i] = h[i] * (1 + g[i]);
	}
}


Matrix_t *LSMOP_sample (int No, int numObj) {
	int H;
	switch (No) {
		case 1:
		case 2:
		case 3:
		case 4:
        		for (H=1; conb (H+numObj-1, numObj-1) < NUM_SAMPLE; H++){};
		        return H1_sample (numObj, H);	
		case 5:
		case 6:
		case 7:
		case 8:
        		for (H=1; conb (H+numObj-1, numObj-1) < NUM_SAMPLE; H++){};
		        return H2_sample (numObj, H);	
		case 9:
        		for (H=1; grid (numObj-1, H) < NUM_SAMPLE; H++){};
		        return H3_sample (numObj, H);	
		default:
		       	fprintf (stderr, "LSMOP%d have been not on consideration\n", No);
               		exit (0);
	}
	return NULL;
}
