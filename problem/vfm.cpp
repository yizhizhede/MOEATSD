#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vfm.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  vfm1 (double *x, int n, double *f, int m);

Problem_t* VFM_new (char *title, int numObj, int numVar) {
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
	problem->numObj = 3; 
	problem->numVar = 2;

	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -2.0;
		uppBound[i] = 2.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "VFM1")) {
		P = problem; problem->evaluate = vfm1;
	} else { 
		fprintf (stderr, "%s:%d:error: %s is undefined\n", __FILE__, __LINE__, title);
		exit (0);
	}

	return problem;
}

static void  vfm1 (double *x, int n, double *f, int m) {
	f[0] = x[0]*x[0] + (x[1]-1)*(x[1]-1);
	f[1] = x[0]*x[0] + (x[1]+1)*(x[1]+1) + 1;
	f[2] = (x[0]-1)*(x[0]-1) + x[1]*x[1] + 2;
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* VFM1_sample ();

Matrix_t* VFM_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	VFM1_sample ();
		default:
                       fprintf (stderr, "VFM%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* VFM1_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}

