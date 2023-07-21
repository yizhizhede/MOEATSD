#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "sch.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  sch1 (double *x, int n, double *f, int m);

Problem_t* SCH_new (char *title, int numObj, int numVar) {
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
	problem->numVar = 1;
	
	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -100;
		uppBound[i] = 100;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "SCH1")) {
		P = problem; problem->evaluate = sch1;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  sch1 (double *x, int n, double *f, int m) {
	if (x[0] < 1)
		f[0] = -x[0];
	else if (x[0] < 3)
		f[0] = -2 + x[0];
	else if (x[0] < 4)
		f[0] = 4 - x[0];
	else 
		f[0] = -4 + x[0];
	f[1] = (x[0]-5)*(x[0]-5);
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* SCH1_sample ();

Matrix_t* SCH_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	SCH1_sample ();
		default:
                       fprintf (stderr, "SCH%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* SCH1_sample () {
	Matrix_t* sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}
