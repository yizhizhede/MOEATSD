#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "lrs.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  lrs1 (double *x, int n, double *f, int m);

Problem_t* LRS_new (char *title, int numObj, int numVar) {
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
	
	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -50;
		uppBound[i] = 50;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "LRS1")) {
		P = problem; problem->evaluate = lrs1;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  lrs1 (double *x, int n, double *f, int m) {
	f[0] = x[0]*x[0] + x[1]*x[1];
	f[0] = (x[0]+2)*(x[0]+2) + x[1]*x[1];
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* LRS1_sample ();

Matrix_t* LRS_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	LRS1_sample ();
		default:
                       fprintf (stderr, "LRS%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* LRS1_sample () {
	int 		i;
	double 		d = 2.0/NUM_SAMPLE;
	double		x[10];
	Matrix_t*	sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		x[0] = -2.0 + i * d;
		x[1] =0;
		lrs1 (x, 2, sample->elements+2*i, 2);
	}
	return sample;
}
