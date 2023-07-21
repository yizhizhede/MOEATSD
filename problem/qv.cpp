#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "qv.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  qv1 (double *x, int n, double *f, int m);

Problem_t* QV_new (char *title, int numObj, int numVar) {
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
	problem->numVar = numVar;
	
	// setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -5.12;
		uppBound[i] = 5.12;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "QV1")) {
		P = problem; problem->evaluate = qv1;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  qv1 (double *x, int n, double *f, int m) {
	double  g; 
	int 	i;

	for (i=0, g=0; i<n; i++) {
		g += x[i]*x[i] - 10*cos(2.0*PI*x[i]) + 10;
	}
	f[0]= pow(g/n, 0.25);

	for (i=0, g=0; i<n; i++) {
		g += (x[i]-1.5)*(x[i]-1.5) - 10*cos(2.0*PI*(x[i]-1.5)) + 10;
	}
	f[1]= pow(g/n, 0.25);
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* QV1_sample ();

Matrix_t* QV_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	QV1_sample ();
		default:
                       fprintf (stderr, "QV%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* QV1_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}
