#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ssfyy.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  ssfyy1 (double *x, int n, double *f, int m);
static void  ssfyy2 (double *x, int n, double *f, int m);

Problem_t* SSFYY_new (char *title, int numObj, int numVar) {
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
	problem->numObj = 2; 
	problem->numVar = 2;
	if (!strcmp (title, "SSFYY2")) {
		problem->numVar = 1;
	}

	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -100;
		uppBound[i] = 100;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "SSFYY1")) {
		P = problem; problem->evaluate = ssfyy1;
	} else if (!strcmp (title, "SSFYY2")) {
		P = problem; problem->evaluate = ssfyy2;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  ssfyy1 (double *x, int n, double *f, int m) {
	f[0] = x[0]*x[0]+x[1]*x[1];
	f[1] = (x[0]-1)*(x[0]-1)+(x[1]-2)*(x[1]-2);
}

static void  ssfyy2 (double *x, int n, double *f, int m) {
	f[0] = 10 + x[0]*x[0] - 10*cos(x[0]*PI/2.0);
	f[1] = (x[0] - 4.0)*(x[0] - 4.0);
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* SSFYY1_sample ();
static Matrix_t* SSFYY2_sample ();

Matrix_t* SSFYY_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	SSFYY1_sample ();
		case 2:
			return	SSFYY2_sample ();
		default:
                       fprintf (stderr, "SSFYY%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* SSFYY1_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}

static Matrix_t* SSFYY2_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}
