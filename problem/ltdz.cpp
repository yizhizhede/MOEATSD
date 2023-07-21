#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ltdz.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  ltdz1 (double *x, int n, double *f, int m);

Problem_t* LTDZ_new (char *title, int numObj, int numVar) {
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
	problem->numVar = 3;
	
	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = 0;
		uppBound[i] = 1;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "LTDZ1")) {
		P = problem; problem->evaluate = ltdz1;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  ltdz1 (double *x, int n, double *f, int m) {
	f[0] = (1+x[2])*cos(x[0]*PI/2.0)*cos(x[1]*PI/2.0);
	f[1] = (1+x[2])*cos(x[0]*PI/2.0)*sin(x[1]*PI/2.0);
	f[2] = (1+x[2])*sin(x[0]*PI/2.0);
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* LTDZ1_sample ();

Matrix_t* LTDZ_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	LTDZ1_sample ();
		default:
                       fprintf (stderr, "LTDZ%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* LTDZ1_sample () {
	Matrix_t* sample = Matrix_new (1, 3);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	sample->elements[2] = 0;

	return sample;
}
