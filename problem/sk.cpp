#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sk.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  sk1 (double *x, int n, double *f, int m);
static void  sk2 (double *x, int n, double *f, int m);

Problem_t* SK_new (char *title, int numObj, int numVar) {
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
	if (!strcmp (title, "SK1")) {
		problem->numObj = 2; 
		problem->numVar = 1;
	} else if (!strcmp (title, "SK2")) {
		problem->numObj = 2; 
		problem->numVar = 4;
	}

	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -10;
		uppBound[i] = 10;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "SK1")) {
		P = problem; problem->evaluate = sk1;
	} else if (!strcmp (title, "SK2")) {
		P = problem; problem->evaluate = sk2;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  sk1 (double *x, int n, double *f, int m) {
	f[0] = x[0]*x[0]*x[0]*x[0]+3.0*x[0]*x[0]*x[0]-10*x[0]*x[0]-10*x[0]-10;
	f[1] = 0.5*x[0]*x[0]*x[0]*x[0]-2.0*x[0]*x[0]*x[0]-10*x[0]*x[0]+10*x[0]-5;
}

static void  sk2 (double *x, int n, double *f, int m) {
	f[0] = (x[0]-2)*(x[0]-2) + (x[1]+3)*(x[1]+3) + (x[2]-5)*(x[2]-5) + (x[3]-4)*(x[3]-4) - 5;
	f[1] = -(sin(x[0]) + sin(x[1]) + sin(x[2]) + sin(x[3]))/(1+(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3])/100.0);
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* SK1_sample ();
static Matrix_t* SK2_sample ();

Matrix_t* SK_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	SK1_sample ();
		case 2:
			return	SK2_sample ();
		default:
                       fprintf (stderr, "SK%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* SK1_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}

static Matrix_t* SK2_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}
