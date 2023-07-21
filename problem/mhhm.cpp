#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mhhm.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  mhhm1 (double *x, int n, double *f, int m);
static void  mhhm2 (double *x, int n, double *f, int m);

Problem_t* MHHM_new (char *title, int numObj, int numVar) {
	// common variable
	int 	i;
	size_t 	size;
	double*	lowBound = NULL; 
	double* uppBound = NULL;

	// allocating memory for a problem
        Problem_t *problem = (Problem_t *)malloc (sizeof (Problem_t));
        if (problem == NULL) {
        	fprintf (stderr, "Allocating memory failed\n");         
               	exit (-1);
        }
	strcpy (problem->title, title);
	if (!strcmp (title, "MHHM1")) {
		problem->numObj = 3; 
		problem->numVar = 1;
	} else if (!strcmp (title, "MHHM2")){
		problem->numObj = 3; 
		problem->numVar = 2;
	}
	
	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = 0.0;
		uppBound[i] = 1.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "MHHM1")) {
		P = problem; problem->evaluate = mhhm1;
	} else if (!strcmp (title, "MHHM2")) {
		P = problem; problem->evaluate = mhhm2;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  mhhm1 (double *x, int n, double *f, int m) {
	f[0] = (x[0]-0.8)*(x[0]-0.8);
	f[1] = (x[0]-0.85)*(x[0]-0.85);
	f[2] = (x[0]-0.9)*(x[0]-0.9);
}

static void  mhhm2 (double *x, int n, double *f, int m) {
	f[0] = (x[0]-0.8)*(x[0]-0.8) + (x[1]-0.6)*(x[1]-0.6);
	f[1] = (x[0]-0.85)*(x[0]-0.85) + (x[1]-0.7)*(x[1]-0.7);
	f[2] = (x[0]-0.9)*(x[0]-0.9) + (x[1]-0.6)*(x[1]-0.6);
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* MHHM1_sample ();
static Matrix_t* MHHM2_sample ();

Matrix_t* MHHM_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	MHHM1_sample ();
		case 2:
			return	MHHM2_sample ();
		default:
                       fprintf (stderr, "MHHM%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* MHHM1_sample () {
	int 	  i;
	double    d = 0.1/NUM_SAMPLE;
	double	  x[10];
	Matrix_t* sample = Matrix_new (NUM_SAMPLE+1, 3);

	for (i=0; i<=NUM_SAMPLE; i++) {
		x[0] = 0.8 + i * d;
		mhhm1 (x, 1, sample->elements+3*i, 3);
	}
	
	return sample;
}
static Matrix_t* MHHM2_sample () {
	int       i, j;
	int 	  n = sqrt (NUM_SAMPLE); 
	double 	  d = 0.1/n;
	double 	  x[10];
	Matrix_t* sample = Matrix_new ((n+1)*(n+1), 3);

	for (i=0; i<=n; i++) {
		for (j=0; j<=n; j++) {
			x[0] = 0.8 + i*d;
			x[1] = 0.6 + i*d;
			mhhm2 (x, 2, sample->elements+3*i, 3);
		}
	}
	
	return sample;
}
