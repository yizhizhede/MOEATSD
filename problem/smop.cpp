#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "smop.h"
#include "matrix.h"
#include "shape.h"
#include "algebra.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  smop1 (double *x, int n, double *f, int m);
static void  smop2 (double *x, int n, double *f, int m);
static void  smop3 (double *x, int n, double *f, int m);
static void  smop4 (double *x, int n, double *f, int m);
static void  smop5 (double *x, int n, double *f, int m);
static void  smop6 (double *x, int n, double *f, int m);
static void  smop7 (double *x, int n, double *f, int m);
static void  smop8 (double *x, int n, double *f, int m);

// landscape function
static double g1 (double x, double t);
static double g2 (double x, double t);
static double g3 (double x, double t);
static double g4 (double x, double t);

// Common parameter
static double theta = 0.1;

Problem_t *SMOP_new (char *title, int numObj, int numVar) {
	// common variable
	int 	i;
	size_t 	size;
	double *lowBound = NULL, *uppBound = NULL;

	// allocating memory for a problem
        Problem_t *problem = (Problem_t *)malloc (sizeof (Problem_t));
        if (problem == NULL) {
        	fprintf (stderr, "Allocating memory failed\n");         
               	exit (-1);
        }
	strcpy (problem->title, title);
	problem->numObj = numObj; 
	problem->numVar = numVar;

	// setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numObj-1; i++) {
		lowBound[i] = 0.0;
		uppBound[i] = 1.0;
	}
	for (i=numObj-1; i<numVar; i++) {
		lowBound[i] = -1.0;
		uppBound[i] = 2.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "SMOP1")) {
		P = problem; problem->evaluate = smop1;
	} else if (!strcmp (title, "SMOP2")) {
		P = problem; problem->evaluate = smop2;
	} else if (!strcmp (title, "SMOP3")) {
		P = problem; problem->evaluate = smop3;
	} else if (!strcmp (title, "SMOP4")) {
		P = problem; problem->evaluate = smop4;
	} else if (!strcmp (title, "SMOP5")) {
		P = problem; problem->evaluate = smop5;
	} else if (!strcmp (title, "SMOP6")) {
		P = problem; problem->evaluate = smop6;
	} else if (!strcmp (title, "SMOP7")) {
		P = problem; problem->evaluate = smop7;
	} else if (!strcmp (title, "SMOP8")) {
		P = problem; problem->evaluate = smop8;
	} else {
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void smop1(double *x,const int n,double *f,const int m) {
	int 	k = ceil(theta * (n - m + 1));	 
	double 	g = 0.0;
	int 	i;

	// compute g
	for (i=m-1; i<m+k-1; i++) {
		g += g1(x[i], PI/3);		
	}
	for (i=m+k-1; i<n; i++) {
		g += g2(x[i], 0);
	}
	g = 1 + g/(n-m+1);

	// shape 
	linear (x, f, m);

	// formula
	for (i=0; i<m; i++) {
		f[i] *= g;
	}
}

static void smop2(double *x,const int n,double* f,const int m) {
	int 	k = ceil(theta * (n - m + 1));	 
	double 	g = 0.0;
	int 	i;

	// compute g
	for (i=m-1; i<m+k-1; i++) {
		g += g2(x[i], PI/3);		
	}
	for (i=m+k-1; i<n; i++) {
		g += g3(x[i], 0);
	}
	g = 1 + g/(n-m+1);

	// shape 
	linear (x, f, m);

	// formula
	for (i=0; i<m; i++) {
		f[i] *= g;
	}
}

static void smop3(double *x,const int n,double *f,const int m) {
	int 	k = ceil(theta * (n - m + 1));	 
	double 	g = 0.0, gns=0, gs = 0;
	int 	i, j;

	// compute g
	for (i=m-1; i<m+k-1; i++) {
		gns += g1(x[i], PI/3);		
	}
	for (i=m+k-1; i<n; i+=10) {
		for (j=i, g=0; j<i+10 && j < n; j++) {
			g += g1(x[j], 0);
		}
		g = (g > 0) ? (50 - g ): 0;
		gs += g;

	}
	g = 1 + (gns + gs)/(n-m+1);

	// shape 
	linear (x, f, m);

	// formula
	for (i=0; i<m; i++) {
		f[i] *= g;
	}
}

static int comp(const void*a,const void*b)
{
	return *(double*)a-*(double*)b;
}

static void  smop4(double *x,const int n,double *f,const int m) {
	int 	k = ceil(theta * (n - m + 1));	 
	double 	g = 0.0;
	int 	i;
	double 	gint[n+10];

	// compute g
	for (i=m-1; i<n; i++) {
		gint[i-m+1]= g3(x[i], PI/3);		
	}

	qsort (gint, n-m+1, sizeof(double), comp);

	for (i=0, g=0; i<n-m-k+1; i++) {
		g += gint[i];
	}

	g = 1 + g/(n-m+1);

	// shape 
	convex (x, f, m);

	// formula
	for (i=0; i<m; i++) {
		f[i] *= g;
	}

}

static void  smop5(double *x,const int n,double *f,const int m) {
	int 	k = ceil(theta * (n - m + 1));	 
	double 	g = 0.0, gns=0, gs=0;
	int 	i;

	// compute g
	for (i=m-1; i<n; i++) {
		gns += g1(x[i], PI/3)*g2(x[i], 0);		
	}

	for (i=m-1, g=0; i<n; i++) {
		gs += (x[i] < 1.0e-10 && x[i] > -1.0e-10) ? 0 : 1;
	}
	gs = fabs(k - gs);

	g = 1 + (gns + gs)/(n-m+1);

	// shape 
	convex (x, f, m);

	// formula
	for (i=0; i<m; i++) {
		f[i] *= g;
	}
}

static void  smop6(double *x,const int n,double *f,const int m) {
	int 	k = ceil(theta * (n - m + 1));	 
	double 	g = 0.0;
	int 	i;
	double 	gmul[n+10];

	// compute g
	for (i=m-1; i<n; i++) {
		gmul[i-m+1]= g4(x[i], (i-m+1.0)/(n-m));		
		if (x[i] <1.0e-10 && x[i] > -1.0e-10) {
			gmul[i-m+1] = 0;
		}
	}

	qsort (gmul, n-m+1, sizeof(double), comp);

	for (i=0, g=0; i<k; i++) {
		g += gmul[i];
	}

	g = 1 + g/(n-m+1);

	// shape 
	convex (x, f, m);

	// formula
	for (i=0; i<m; i++) {
		f[i] *= g;
	}
}


static void smop7(double *x,const int n,double *f,const int m) {
	int 	k = ceil(theta * (n - m + 1));	 
	double 	g = 0.0, gns=0, gs=0;
	int 	i;

	// compute g
	for (i=m-1; i<m+k-1; i++) {
		gns += g2(x[i], PI/3);		
	}
	for (i=m+k-1; i<n-1; i++) {
		gs += g2(x[i], 0.9*x[i+1]);
	}
	gs += g2(x[n-1], 0.9*x[m+k-1]);

	g = 1 + (gns+gs)/(n-m+1);

	// shape 
	H2(x, NULL, f, m);

	// formula
	for (i=0; i<m; i++) {
		f[i] *= g;
	}
}

static void  smop8 (double *x, int n, double *f, int m) {
	int 	k = ceil(theta * (n - m + 1));	 
	double 	g = 0.0, gns=0, gs=0, t;
	int 	i;

	// compute g
	for (i=m-1; i<m+k-1; i++) {
		t = (x[i+1]+PI) - floor((x[i+1]+PI)/2)*2;
		gns += g2(x[i], t);		
	}
	for (i=m+k-1; i<n-1; i++) {
		gs += g2(x[i], 0.9*x[i+1]);
	}
	gs += g2(x[n-1], 0.9*x[m+k-1]);

	g = 1 + (gns+gs)/(n-m+1);

	// shape 
	H2(x, NULL, f, m);

	// formula
	for (i=0; i<m; i++) {
		f[i] *= g;
	}
}

//**************************************************************************************************************
//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t *SMOP4_sample (int M);

Matrix_t *SMOP_sample (int No, int numObj) {
	int 	H;

        switch (No) {
         	case 1:
             	case 2:
                case 3:
                	for (H=2; conb (H+numObj-1, numObj-1) < NUM_SAMPLE; H++){};
	                return H1_sample (numObj, H);
                case 4:
                case 5:
                case 6:
			return SMOP4_sample (numObj);			
                case 7:
                case 8:
                      	for (H=1; grid (numObj-1, H) < NUM_SAMPLE; H++){};
                       	return H2_sample (numObj, H);
                default:
                       fprintf (stderr, "SMOP%d have been not on consideration\n", No);
                       exit (0);
       	}
        return NULL;
}

static double g1 (double x, double t) {
	return (x-t)*(x-t);
}

static double g2 (double x, double t) {
	return 2*(x-t)*(x-t) + sin(2*PI*(x-t))*sin(2*PI*(x-t));
}

static double g3 (double x, double t) {
	return 4 - (x-t) - 4.0/exp(100*(x-t)*(x-t));
}

static double g4 (double x, double t) {
	return (x-PI/3)*(x-PI/3) + t*sin(6*PI*(x-PI/3))*sin(6*PI*(x-PI/3));
}

static Matrix_t *SMOP4_sample (int M) {
	Matrix_t *sample = Matrix_new ();
	Matrix_t *Grid = NULL;
        int 	i, H;
        size_t 	size;
         
        sample->colDim = M;
	for (H=1; grid(M-1, H) < NUM_SAMPLE; H++){};
        sample->rowDim = grid(M-1, H);
         
        size = sample->rowDim * sample->colDim * sizeof (double);
        sample->elements = (double *)malloc (size);
	
	Grid = Grid_sample (M-1, H);

	for (i=0; i<sample->rowDim; i++) {
		convex (Grid->elements+i*(M-1), sample->elements+i*M, M);
	}

	Matrix_free (&Grid);
	return sample;
}
