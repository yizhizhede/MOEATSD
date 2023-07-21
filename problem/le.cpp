#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "le.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  le1 (double *x, int n, double *f, int m);

Problem_t* LE_new (char *title, int numObj, int numVar) {
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
		lowBound[i] = -5.0;
		uppBound[i] = 10.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "LE1")) {
		P = problem; problem->evaluate = le1;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  le1 (double *x, int n, double *f, int m) {
	f[0] = pow(x[0]*x[0]+x[1]*x[1], 0.125);
	f[0] = pow((x[0]-0.5)*(x[0]-0.5)+(x[1]-0.5)*(x[1]-0.5), 0.25);
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* LE1_sample ();

Matrix_t* LE_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	LE1_sample ();
		default:
                       fprintf (stderr, "LE%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* LE1_sample () {
	int 		i;
	double 		d = 0.5/NUM_SAMPLE;
	double		x;
	Matrix_t*	sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		x = i * d;
		sample->elements[2*i+0] = pow(2.0*x*x, 0.125);
		sample->elements[2*i+1] = pow(2.0*(x-0.5)*(x-0.5), 0.25);
	
	}
	return sample;
}
