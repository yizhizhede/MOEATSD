#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "dpam.h"
#include "matrix.h"
#include "algebra.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  dpam1 (double *x, int n, double *f, int m);

Problem_t *DPAM_new (char *title, int numObj, int numVar) {
	// common variable
	int 	i;
	size_t 	size;
	double*	lowBound = NULL; 
	double*	uppBound = NULL;

	// allocating memory for a problem
        Problem_t *problem = (Problem_t *)malloc (sizeof (Problem_t));
        if (problem == NULL) {
        	fprintf (stderr, "Allocating memory failed\n");         
               	exit (-1);
        }
	strcpy (problem->title, title);
	problem->numObj = 2; 
	problem->numVar = numVar;
	
	// setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -0.3;
		uppBound[i] = 0.3;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "DPAM1")) {
		P = problem; problem->evaluate = dpam1;
	} else {
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void dpam1(double *x,const int n,double *f,const int m) {
	double 		g = 0.0;
	int 		i, j;
	Matrix_t*	R = Matrix_new (n, n);
	double		y[n+10];

	// initiaize R
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			if (i==j || i+1 ==j) {
				R->elements[i*n+j] = 1;
			} else {
				R->elements[i*n+j] = 0;
			}
		}
	}

	// y
	for (i=0; i<n; i++) {
		for (j=0, y[i]=0; j<n; j++) {
			y[i] += R->elements[i*n+j]*x[j];
		}
	}
	
	// compute g
	for (i=1, g=0; i<n; i++) {
		g += y[i]*y[i]-10*cos(4*PI*y[i]);
	}
	g = 1 + 10*(n-1) + g;

	// f
	f[0] = y[0];
	f[1] = g*exp(-y[0]/g);
}


//**************************************************************************************************************
//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t *dpam1_sample ();

Matrix_t *DPAM_sample (int No, int numObj) {
        switch (No) {
         	case 1:
                	return dpam1_sample ();
                default:
                       fprintf (stderr, "DPAM%d have been not on consideration\n", No);
                       exit (0);
       	}
        return NULL;
}

static Matrix_t *dpam1_sample () {
	int 		i;
	double		d = 0.6/NUM_SAMPLE, y1;
	Matrix_t*	sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		y1 = -0.3 + i*d;
		sample->elements[i*2+0] = y1;
		sample->elements[i*2+1] = exp(-y1);
	}
	
	return sample;
}
