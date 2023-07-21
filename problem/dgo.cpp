#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "dgo.h"
#include "matrix.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  dgo1 (double *x, int n, double *f, int m);
static void  dgo2 (double *x, int n, double *f, int m);

Problem_t *DGO_new (char *title, int numObj, int numVar) {
	// common variable
	size_t 	size;
	double* lowBound = NULL; 
	double* uppBound = NULL;

	// allocating memory for a problem
        Problem_t *problem = (Problem_t *)malloc (sizeof (Problem_t));
        if (problem == NULL) {
        	fprintf (stderr, "Allocating memory failed\n");         
               	exit (-1);
        }
	strcpy (problem->title, title);
	problem->numObj = 2; 
	problem->numVar = 1;
	
	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	if (!strcmp (title, "DGO1")) {
		lowBound[0] = -10.0;
		uppBound[0] = 13.0;
	} else if (!strcmp (title, "DGO2")) {
		lowBound[0] = -9.0;
		uppBound[0] = 9.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "DGO1")) {
		P = problem; problem->evaluate = dgo1;
	} else if (!strcmp (title, "DGO2")) {
		P = problem; problem->evaluate = dgo2;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}


	return problem;
}

static void  dgo1 (double *x, int n, double *f, int m) {
	f[0] = sin(x[0]);
	f[1] = sin(x[0]+0.7);
}

static void  dgo2 (double *x, int n, double *f, int m) {
	f[0] = x[0]*x[0];
	f[1] = 9 - sqrt (81-x[0]*x[0]);
}


//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t  *DGO1_sample ();
static Matrix_t  *DGO2_sample ();

Matrix_t* DGO_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	DGO1_sample ();
		case 2:
			return	DGO2_sample ();
		default:
                       fprintf (stderr, "DGO%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* DGO1_sample () {
	int 		i;
	double 		d = 0.7/NUM_SAMPLE; 
	double 		x[10];
	Matrix_t*	sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		x[0] = i * d;
		dgo1(x, 1, sample->elements+i*2, 2);
	}
	
	return sample;
}
static Matrix_t* DGO2_sample () {
	Matrix_t *sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}
