#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "far.h"
#include "matrix.h"

#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  far1 (double *x, int n, double *f, int m);

Problem_t *FAR_new (char *title, int numObj, int numVar) {
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
	problem->numVar = 2;
	
	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -1.0;
		uppBound[i] = 1.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "FAR1")) {
		P = problem; problem->evaluate = far1;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  far1 (double *x, int n, double *f, int m) {
	f[0] = -2.0*exp(15*(-(x[0]-0.1)*(x[0]-0.1)-x[1]*x[1]))  
	       -exp(20*(-(x[0]-0.6)*(x[0]-0.6)-(x[1]-0.6)*(x[1]-0.6))) 
	       +exp(20*(-(x[0]+0.6)*(x[0]+0.6)-(x[1]-0.6)*(x[1]-0.6))) 
	       +exp(20*(-(x[0]-0.6)*(x[0]-0.6)-(x[1]+0.6)*(x[1]+0.6))) 
	       +exp(20*(-(x[0]+0.6)*(x[0]+0.6)-(x[1]+0.6)*(x[1]+0.6)));
	f[1] = 2.0*exp(20*(-x[0]*x[0]-x[1]*x[1]))  
	       +exp(20*(-(x[0]-0.4)*(x[0]-0.4)-(x[1]-0.6)*(x[1]-0.6)))
	       -exp(20*(-(x[0]+0.5)*(x[0]+0.5)-(x[1]-0.7)*(x[1]-0.7))) 
	       -exp(20*(-(x[0]-0.5)*(x[0]-0.5)-(x[1]+0.7)*(x[1]+0.7))) 
	       +exp(20*(-(x[0]+0.4)*(x[0]+0.4)-(x[1]+0.8)*(x[1]+0.8)));
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t *FAR1_sample ();

Matrix_t *FAR_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	FAR1_sample ();
		default:
                       fprintf (stderr, "BK%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t *FAR1_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0]=0;
	sample->elements[1]=0;
	
	return sample;
}
