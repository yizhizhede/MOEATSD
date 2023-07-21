#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "zlt.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  zlt1 (double *x, int n, double *f, int m);

Problem_t* ZLT_new (char *title, int numObj, int numVar) {
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
	problem->numObj = numObj; 
	problem->numVar = numVar;
	if (numVar < numObj) {
		problem->numVar = numObj;
	}

	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -1000.0;
		uppBound[i] = 1000.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "ZLT1")) {
		P = problem; problem->evaluate = zlt1;
	} else { 
		fprintf (stderr, "%s:%d:error: %s is undefined\n", __FILE__, __LINE__, title);
		exit (0);
	}

	return problem;
}

static void  zlt1 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g;

	for (i=0, g=0; i<n; i++) {
		g += x[i]*x[i];
	}

	for (i=0; i<m; i++) {
		f[i] = g + 1.0 - 2.0*x[i];
	}
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* ZLT1_sample (int numObj);

Matrix_t* ZLT_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	ZLT1_sample (numObj);
		default:
                       fprintf (stderr, "ZLT%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* ZLT1_sample (int numObj) {
	int i;
	Matrix_t* sample = Matrix_new (1, numObj);

	for (i=0; i<numObj; i++) {
		sample->elements[i] = 0;
	}
	
	return sample;
}

