#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mlf.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  mlf1 (double *x, int n, double *f, int m);
static void  mlf2 (double *x, int n, double *f, int m);

Problem_t* MLF_new (char *title, int numObj, int numVar) {
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
	if (!strcmp (title, "MLF1")) {
		problem->numObj = 2; 
		problem->numVar = 1;
	} else if (!strcmp (title, "MLF2")) {
		problem->numObj = 2; 
		problem->numVar = 2;
	}
	
	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = 0.0;
		uppBound[i] = 20.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "MLF1")) {
		P = problem; problem->evaluate = mlf1;
	} else if (!strcmp (title, "MLF2")) {
		P = problem; problem->evaluate = mlf2;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  mlf1 (double *x, int n, double *f, int m) {
	f[0] = (1.0 + x[0]/20)*sin(x[0]);
	f[0] = (1.0 + x[0]/20)*cos(x[0]);
}

static void  mlf2 (double *x, int n, double *f, int m) {
	f[0] = -5 + ((x[0]*x[0]+x[1]-11)*(x[0]*x[0]+x[1]-11)+(x[0]+x[1]*x[1]-7)*(x[0]+x[1]*x[1]-7))/200.0;
	f[1] = -5 + ((4.0*x[0]*x[0]+2.0*x[1]-11)*(4.0*x[0]*x[0]+2.0*x[1]-11)+
		     (2.0*x[0]+4.0*x[1]*x[1]-7)*(2.0*x[0]+4.0*x[1]*x[1]-7))/200.0;
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* MLF1_sample ();
static Matrix_t* MLF2_sample ();

Matrix_t* MLF_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	MLF1_sample ();
		case 2:
			return	MLF2_sample ();
		default:
                       fprintf (stderr, "MLF%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* MLF1_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}
static Matrix_t* MLF2_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}
