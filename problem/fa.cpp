#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fa.h"
#include "matrix.h"

#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  fa1 (double *x, int n, double *f, int m);

Problem_t *FA_new (char *title, int numObj, int numVar) {
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
	problem->numObj = 3; 
	problem->numVar = 3;
	
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
	if (!strcmp (title, "FA1")) {
		P = problem; problem->evaluate = fa1;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  fa1 (double *x, int n, double *f, int m) {
	f[0] = (1-exp(-4.0*x[0]))/(1-exp(-4.0));
	f[1] = (x[1]+1.0)*(1.0-pow(f[0]/(x[1]+1.0),0.5));
	f[2] = (x[2]+1.0)*(1.0-pow(f[0]/(x[2]+1.0),0.1));
}


//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t *FA1_sample ();

Matrix_t *FA_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	FA1_sample ();
		default:
                       fprintf (stderr, "BK%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t *FA1_sample () {
	Matrix_t* sample = Matrix_new (1, 3);

	sample->elements[0]=0;
	sample->elements[1]=0;
	sample->elements[2]=0;
	
	return sample;
}
