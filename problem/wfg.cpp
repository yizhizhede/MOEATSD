#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "wfg.h"
#include "matrix.h"
#include "transition.h"
#include "shape.h"
#include "algebra.h"

#define NUM_SAMPLE 1.0e+4

#define PI 3.14159265358979323846264338327950288419716939937510
#define MAX(a,b) ((a>b)?(a):(b))
#define MIN(a,b) ((a<b)?(a):(b))

// constants
static double D;	// distance-related variable's scaling constant
static double S[100];	// scaling constants
static double A[100];	// degeneracy constants
static int    nk;	// the number of position-ralatived variables
static int    nl;	// the numer of distance-relatived variables

// test problem
static Problem_t *P;

// 
static void wfg1 (double *x, int n, double *f, int m);
static void wfg2 (double *x, int n, double *f, int m);
static void wfg3 (double *x, int n, double *f, int m);
static void wfg4 (double *x, int n, double *f, int m);
static void wfg5 (double *x, int n, double *f, int m);
static void wfg6 (double *x, int n, double *f, int m);
static void wfg7 (double *x, int n, double *f, int m);
static void wfg8 (double *x, int n, double *f, int m);
static void wfg9 (double *x, int n, double *f, int m);

Problem_t *WFG_new (char *title, int numObj, int numVar) {
	int 		i;
	size_t 		size;
	double* 	lowBound = NULL;
	double* 	uppBound = NULL;
	Problem_t*	problem = NULL;

	// initial constants: nk, nl, D, S, A.
	nk = 2*(numObj-1);		
	nl = numVar - nk;
	if (nl < 1) {
		nl = 1;
	}
	if (!strcmp ("WFG2", title) || !strcmp ("WFG3", title)) {
		nl = ((nl&1) == 0) ? (nl) : (nl+1);
	}

	D = 1;
	for (i=0; i<numObj; i++) {
		S[i] = 2.0*(i+1);
		A[i] = 1.0;
	}
	if (!strcmp (title, "WFG3")) {
		memset (A, 0, 100*sizeof (double));
		A[0] = 1.0;
	}

	// allocating memory for a problem
        if (NULL == (problem = (Problem_t *)calloc (1, sizeof (Problem_t)))) {
        	fprintf (stderr, "ERROR:%s:%d:Allocating memory failed\n", __FILE__, __LINE__);         
               	exit (0);
        }
	strcpy (problem->title, title);
	problem->numObj = numObj; 

	numVar = nk + nl;
	problem->numVar = numVar;
	
	// setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = 0.0;
		uppBound[i] = 2.0*(i+1);
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "WFG1")) {
		P = problem; problem->evaluate = wfg1;
	} else if (!strcmp (title, "WFG2")) {
		P = problem; problem->evaluate = wfg2;
	} else if (!strcmp (title, "WFG3")) {
		P = problem; problem->evaluate = wfg3;
	} else if (!strcmp (title, "WFG4")) {
		P = problem; problem->evaluate = wfg4;
	} else if (!strcmp (title, "WFG5")) {
		P = problem; problem->evaluate = wfg5;
	} else if (!strcmp (title, "WFG6")) {
		P = problem; problem->evaluate = wfg6;
	} else if (!strcmp (title, "WFG7")) {
		P = problem; problem->evaluate = wfg7;
	} else if (!strcmp (title, "WFG8")) {
		P = problem; problem->evaluate = wfg8;
	} else if (!strcmp (title, "WFG9")) {
		P = problem; problem->evaluate = wfg9;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void wfg1 (double *x, int n, double *f, int m) {
	int 	i;
	double 	z[n+10];
	double 	h[m+10];
	double 	w[n+10];

	// copy x to z 
	memcpy (z, x, n*sizeof (double));

	// t^0 
	for (i=0; i<n; i++) {
		z[i] /= P->uppBound[i];
		w[i] = 2.0*(i+1);
	}

	// t^1
	for (i=nk; i<n; i++) {
		s_linear (z+i, 0.35);
	}

	// t^2
	for (i=nk; i<n; i++) {
		b_flat (z+i, 0.8, 0.75, 0.85);
	}

	// t^3
	for (i=0; i<n; i++) {
		b_poly (z+i, 0.02);
	}


	// t^4
	for (i=0; i<m-1; i++) {
		z[i] = r_sum (z+i*nk/(m-1), w+i*nk/(m-1), nk/(m-1));
	}
	z[m-1] = r_sum (z+nk, w+nk, nl);

	// t^(p+1)
	for (i=0; i<m-1; i++) {
		z[i] = MAX(z[m-1],A[i])*(z[i]-0.5)+0.5;
	}

	// shape: convex & mixed
	convex (z, h, m);	
	mixedM (z, h, m);

	// fitness
	for (i=0; i<m; i++) {
		f[i] = D*z[m-1] + S[i]*h[i];
	}
}


static void wfg2 (double *x, int n, double *f, int m) {
	int 	i;
	double 	z[n+10];
	double 	h[m+10];
	double 	w[n+10];

	// copy x to z 
	memcpy (z, x, n*sizeof (double));

	// t^0 
	for (i=0; i<n; i++) {
		z[i] /= P->uppBound[i];
		w[i] = 1.0;
	}

	// t^1
	for (i=nk; i<n; i++) {
		s_linear (z+i, 0.35);
	}

	// t^2
	for (i=nk; i<nk+nl/2; i++) {
		z[i] = r_nonsep (z+nk+2*(i-nk), 2, 2);
	}

	// t^3
	for (i=0; i<m-1; i++) {
		z[i] = r_sum (z+i*nk/(m-1), w+i*nk/(m-1), nk/(m-1));
	}
	z[m-1] = r_sum (z+nk, w+nk, nl/2);

	// t^(p+1)
	for (i=0; i<m-1; i++) {
		z[i] = MAX(z[m-1],A[i])*(z[i]-0.5)+0.5;
	}

	// shape: convex & discM
	convex(z, h, m);	
	discM (z, h, m);

	// fitness
	for (i=0; i<m; i++) {
		f[i] = D*z[m-1] + S[i]*h[i];
	}
}

static void wfg3 (double *x, int n, double *f, int m) {
	int 	i;
	double 	z[n+10];
	double 	h[m+10];
	double 	w[n+10];

	// copy x to z 
	memcpy (z, x, n*sizeof (double));

	// t^0 
	for (i=0; i<n; i++) {
		z[i] /= P->uppBound[i];
		w[i] = 1.0;
	}

	// t^1
	for (i=nk; i<n; i++) {
		s_linear (z+i, 0.35);
	}

	// t^2
	for (i=nk; i<nk+nl/2; i++) {
		z[i] = r_nonsep (z+nk+2*(i-nk), 2, 2);
	}

	// t^3
	for (i=0; i<m-1; i++) {
		z[i] = r_sum (z+i*nk/(m-1), w+i*nk/(m-1), nk/(m-1));
	}
	z[m-1] = r_sum (z+nk, w+nk, nl/2);

	// t^(p+1)
	for (i=0; i<m-1; i++) {
		z[i] = MAX(z[m-1],A[i])*(z[i]-0.5)+0.5;
	}

	// shape
	linear (z, h, m);	

	// fitness
	for (i=0; i<m; i++) {
		f[i] = D*z[m-1] + S[i]*h[i];
	}
}

static void wfg4 (double *x, int n, double *f, int m) {
	int 	i;
	double 	z[n+10];
	double 	h[m+10];
	double 	w[n+10];

	// copy x to z 
	memcpy (z, x, n*sizeof (double));

	// t^0 
	for (i=0; i<n; i++) {
		z[i] /= P->uppBound[i];
		w[i] = 1.0;
	}

	// t^1
	for (i=0; i<n; i++) {
		s_multi (z+i, 30, 10, 0.35);
	}

	// t^2
	for (i=0; i<m-1; i++) {
		z[i] = r_sum(z+i*nk/(m-1), w+i*nk/(m-1), nk/(m-1));
	}
	z[m-1] = r_sum (z+nk, w+nk, nl);

	// t^(p+1)
	for (i=0; i<m-1; i++) {
		z[i] = MAX(z[m-1],A[i])*(z[i]-0.5)+0.5;
	}

	// shape
	concave (z, h, m);	

	// fitness
	for (i=0; i<m; i++) {
		f[i] = D*z[m-1] + S[i]*h[i];
	}
}

static void wfg5 (double *x, int n, double *f, int m) {
	int 	i;
	double 	z[n+10];
	double 	h[m+10];
	double 	w[n+10];

	// copy x to z 
	memcpy (z, x, n*sizeof (double));

	// t^0 
	for (i=0; i<n; i++) {
		z[i] /= P->uppBound[i];
		w[i] = 1.0;
	}

	// t^1
	for (i=0; i<n; i++) {
		s_decept (z+i, 0.35, 0.001, 0.05);
	}

	// t^2
	for (i=0; i<m-1; i++) {
		z[i] = r_sum(z+i*nk/(m-1), w+i*nk/(m-1), nk/(m-1));
	}
	z[m-1] = r_sum (z+nk, w+nk, nl);

	// t^(p+1)
	for (i=0; i<m-1; i++) {
		z[i] = MAX(z[m-1],A[i])*(z[i]-0.5)+0.5;
	}

	// shape: concave 
	concave (z, h, m);	

	// fitness
	for (i=0; i<m; i++) {
		f[i] = D*z[m-1] + S[i]*h[i];
	}
}

static void wfg6 (double *x, int n, double *f, int m) {
	int 	i;
	double 	z[n+10];
	double 	h[m+10];

	// copy x to z 
	memcpy (z, x, n*sizeof (double));

	// t^0 
	for (i=0; i<n; i++) {
		z[i] /= P->uppBound[i];
	}

	// t^1
	for (i=nk; i<n; i++) {
		s_linear (z+i, 0.35);
	}

	// t^2
	for (i=0; i<m-1; i++) {
		z[i] = r_nonsep(z+i*nk/(m-1), nk/(m-1), nk/(m-1));
	}
	z[m-1] = r_nonsep (z+nk, nl, nl);

	// t^(p+1)
	for (i=0; i<m-1; i++) {
		z[i] = MAX(z[m-1],A[i])*(z[i]-0.5)+0.5;
	}

	// shape
	concave (z, h, m);	

	// fitness
	for (i=0; i<m; i++) {
		f[i] = D*z[m-1] + S[i]*h[i];
	}
}

static void wfg7 (double *x, int n, double *f, int m) {
	int 	i;
	double 	z[n+10];
	double 	b[n+10];
	double 	h[m+10];
	double 	w[n+10];
	double 	u;

	// copy x to z 
	memcpy (z, x, n*sizeof (double));

	// t^0 
	for (i=0; i<n; i++) {
		z[i] /= P->uppBound[i];
		w[i] = 1.0;
	}

	// backup z
	memcpy (b, z, n*sizeof (double));

	// t^1
	for (i=nk+1, u=0; i<n; i++) {
		u += b[i];
	}
	for (i=nk-1; i>=0; i--) {
		u += b[i+1];					// reduce funciton
		b_param (z+i, u/(n-i-1), 0.98/49.98, 0.02, 50);	
	}

	// t^2
	for (i=nk; i<n; i++) {
		s_linear (z+i, 0.35);
	}

	// t^3
	for (i=0; i<m-1; i++) {
		z[i] = r_sum(z+i*nk/(m-1), w+i*nk/(m-1), nk/(m-1));
	}
	z[m-1] = r_sum (z+nk, w+nk, nl);

	// t^(p+1)
	for (i=0; i<m-1; i++) {
		z[i] = MAX(z[m-1],A[i])*(z[i]-0.5)+0.5;
	}

	// shape
	concave (z, h, m);	

	// fitness
	for (i=0; i<m; i++) {
		f[i] = D*z[m-1] + S[i]*h[i];
	}
}

static void wfg8 (double *x, int n, double *f, int m) {
	int 	i;
	double 	z[n+10];
	double 	b[n+10];	// backup
	double 	h[m+10];
	double 	w[n+10];
	double 	u;

	// copy x to z 
	memcpy (z, x, n*sizeof (double));

	// t^0 
	for (i=0; i<n; i++) {
		z[i] /= P->uppBound[i];
		w[i] = 1.0;
	}

	// backup z
	memcpy (b, z, n*sizeof (double));

	// t^1
	for (i=0, u=0; i<nk-1; i++) {
		u += b[i];	
	}
	for (i=nk; i<n; i++) {
		u += b[i-1];					// reduce function
		b_param (z+i, u/i, 0.98/49.98, 0.02, 50);	
	}

	// t^2
	for (i=nk; i<n; i++) {
		s_linear (z+i, 0.35);
	}

	// t^3
	for (i=0; i<m-1; i++) {
		z[i] = r_sum (z+i*nk/(m-1), w+i*nk/(m-1), nk/(m-1));
	}
	z[m-1] = r_sum (z+nk, w+nk, nl);

	// t^(p+1)
	for (i=0; i<m-1; i++) {
		z[i] = MAX(z[m-1],A[i])*(z[i]-0.5)+0.5;
	}

	// shape:
	concave (z, h, m);	

	// fitness
	for (i=0; i<m; i++) {
		f[i] = D*z[m-1] + S[i]*h[i];
	}
}

static void wfg9 (double *x, int n, double *f, int m) {
	int 	i;
	double 	z[n+10];
	double 	b[n+10];	// backup
	double 	h[m+10];
	double	u;

	// copy x to z 
	memcpy (z, x, n*sizeof (double));

	// t^0 
	for (i=0; i<n; i++) {
		z[i] /= P->uppBound[i];
	}

	// backup z
	memcpy (b, z, n*sizeof (double));

	// t^1
	for (i=n-2, u=0; i>=0; i--) {
		u += b[i+1];					// reduce function
		b_param (z+i, u/(n-i-1), 0.98/49.98, 0.02, 50);	
	}

	// t^2
	for (i=0; i<nk; i++) {
		s_decept (z+i, 0.35, 0.001, 0.05);
	}
	for (i=nk; i<n; i++) {
		s_multi (z+i, 30, 95, 0.35);
	}

	// t^3
	for (i=0; i<m-1; i++) {
		z[i] = r_nonsep(z+i*nk/(m-1), nk/(m-1), nk/(m-1));
	}
	z[m-1] = r_nonsep (z+nk, nl, nl);	

	// t^(p+1)
	for (i=0; i<m-1; i++) {
		z[i] = MAX(z[m-1],A[i])*(z[i]-0.5)+0.5;
	}

	// shape:
	concave (z, h, m);	

	// fitness
	for (i=0; i<m; i++) {
		f[i] = D*z[m-1] + S[i]*h[i];
	}
}

static Matrix_t *WFG1_sample (int M);
static Matrix_t *WFG2_sample (int M);
static Matrix_t *WFG3_sample (int M);
static Matrix_t *WFG4_sample (int M);

Matrix_t *WFG_sample (int No, int M) {
	Matrix_t * sample = NULL;
	
	switch (No) {
		case 1:
			return WFG1_sample (M);
		case 2:
			return WFG2_sample (M);
		case 3:
			return WFG3_sample (M);
		case 4:
		case 5:
		case 6:
		case 7:
		case 8:
		case 9:
			return WFG4_sample (M);
		default:
			fprintf (stderr, "WFG%d have not been implemented now\n", No);
			exit(0);
	}

	return sample;
}

static Matrix_t *WFG1_sample (int M) {
	Matrix_t *sample = Matrix_new ();
	Matrix_t *Grid = NULL;
	Matrix_t *front = NULL;
        int 	i, j, H;
        size_t 	size;
         
        sample->colDim = M;
	for (H=1; grid(M-1, H) < NUM_SAMPLE; H++){};
        sample->rowDim = grid(M-1, H);
         
        size = sample->rowDim * sample->colDim * sizeof (double);
        sample->elements = (double *)malloc (size);
	
	Grid = Grid_sample (M-1, H);

	for (i=0; i<sample->rowDim; i++) {
		convex (Grid->elements+i*(M-1), sample->elements+i*M, M);
		mixedM (Grid->elements+i*(M-1), sample->elements+i*M, M);
		for (j=0; j<M; j++) {
			sample->elements[i*M+j] *= 2*(j+1);
		}
	}

	front = Matrix_front (sample);
	Matrix_free (&Grid);
	Matrix_free (&sample);
	return front;
}

static Matrix_t *WFG2_sample (int M) {
	Matrix_t *sample = Matrix_new ();
	Matrix_t *Grid = NULL;
	Matrix_t *front = NULL;
        int 	i, j, H;
        size_t 	size;
         
        sample->colDim = M;
	for (H=1; grid(M-1, H) < NUM_SAMPLE; H++){};
        sample->rowDim = grid(M-1, H);
         
        size = sample->rowDim * sample->colDim * sizeof (double);
        sample->elements = (double *)malloc (size);
	
	Grid = Grid_sample (M-1, H);

	for (i=0; i<sample->rowDim; i++) {
		convex (Grid->elements+i*(M-1), sample->elements+i*M, M);
		discM (Grid->elements+i*(M-1), sample->elements+i*M, M);
		for (j=0; j<M; j++) {
			sample->elements[i*M+j] *= 2*(j+1);
		}
	}

	front = Matrix_front (sample);
	Matrix_free (&Grid);
	Matrix_free (&sample);
	return front;
}

static Matrix_t *WFG3_sample (int M) {
	Matrix_t *sample = Matrix_new ();
	Matrix_t *L = NULL;
        int 	i, j, H;
        size_t 	size;
         
        sample->colDim = M;
	for (H=1;conb(H+2-1, 2-1) < NUM_SAMPLE+1; H++){};
        sample->rowDim = conb(H+2-1, 2-1);
         
        size = sample->rowDim * sample->colDim * sizeof (double);
        sample->elements = (double *)malloc (size);
        memset (sample->elements, 0, size);

	L = H1_sample (2, H);

	for (i=0; i<sample->rowDim; i++) {
		for (j=M-2; j<M; j++) {
			sample->elements[i*M+j] = 2*(j+1)*L->elements[i*2+j-(M-2)];
		}
	}

	Matrix_free (&L);
	return sample;
}

static Matrix_t *WFG4_sample (int M) {
	Matrix_t *sample = NULL;
        int 	i, j, H;
         
	for (H=1;conb(H+M-1, M-1) < NUM_SAMPLE; H++){};
	sample = H2_sample (M, H);

	for (i=0; i<sample->rowDim; i++) {
		for (j=0; j<M; j++) {
			sample->elements[i*M+j] *= 2*(j+1);
		}
	}

	return sample;
}
