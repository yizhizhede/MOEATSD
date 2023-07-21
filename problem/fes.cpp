#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "fes.h"
#include "matrix.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  fes1 (double *x, int n, double *f, int m);
static void  fes2 (double *x, int n, double *f, int m);
static void  fes3 (double *x, int n, double *f, int m);

Problem_t* FES_new (char *title, int numObj, int numVar) {
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
	if (!strcmp (title, "FES1")) {
		problem->numObj = 2; 
	} else if (!strcmp (title, "FES2")) {
		problem->numObj = 3; 
	} else if (!strcmp (title, "FES3")) {
		problem->numObj = 4; 
	} else {
		fprintf (stdout, "%s:%d: It must be one of FES1,FES2 or FES3\n", __FILE__, __LINE__);
		exit (0);
	}
	problem->numVar = numVar;
	
	// setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = 0.0;
		uppBound[i] = 1.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "FES1")) {
		P = problem; problem->evaluate = fes1;
	} else if (!strcmp (title, "FES2")) {
		P = problem; problem->evaluate = fes2;
	} else if (!strcmp (title, "FES3")) {
		P = problem; problem->evaluate = fes3;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  fes1 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g, t;

	// f1 
	for (i=0, g=0; i<n; i++) {
		t = x[i] - exp((i+1.0)*(i+1.0)/(n*n))/3.0;
		if (t > DBL_EPSILON)
			g += pow (t, 0.5);
		else if (t > -DBL_EPSILON) {
			g+= 0;
		} else {
			g += pow(-t,0.5);
		}
	}
	f[0] = g;

	// f2 
	for (i=0, g=0; i<n; i++) {
		t = x[i] - 0.5*cos(10.0*PI*(i+1.0)/n)-0.5;
		g += t*t;
	}
	f[1] = g;
}

static void  fes2 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0;

	// f1
	for (i=0, g=0; i<n; i++) {
		t = x[i] - 0.5*cos(10.0*PI*(i+1.0)/n) - 0.5;
		g += t*t;
	}
	f[0] = g;

	// f2 
	for (i=0, g=0; i<n; i++) {
		t = x[i] - sin(i)*sin(i)*cos(i)*cos(i);
		if (t > DBL_EPSILON)
			g += pow (t, 0.5);
		else if (t > -DBL_EPSILON) {
			g+= 0;
		} else {
			g += pow(-t,0.5);
		}
	}
	f[1] = g;

	// f3
	for (i=0, g=0; i<n; i++) {
		t = x[i] - 0.25*cos(i)*cos(2.0*i) - 0.5;
		if (t > DBL_EPSILON) {
			g += pow (t, 0.5);
		} else if (t > -DBL_EPSILON) {
			g += 0;
		} else {
			g += pow(-t,0.5);
		}
	}
	f[2] = g;
}

static void  fes3 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0;

	// f1 
	for (i=0, g=0; i<n; i++) {
		t = x[i] - exp((i+1.0)*(i+1.0)/(n*n))/3.0;
		g += (t>0)?pow(t,0.5):pow(-t,0.5);
	}
	f[0] = g;

	// f2 
	for (i=0, g=0; i<n; i++) {
		t = x[i] - sin(i)*sin(i)*cos(i)*cos(i);
		g += (t>0)?pow(t,0.5):pow(-t,0.5);
	}
	f[1] = g;

	// f3
	for (i=0, g=0; i<n; i++) {
		t = x[i] - 0.25*cos(i)*cos(2.0*i) - 0.5;
		g += (t>0)?pow(t,0.5):pow(-t,0.5);
	}
	f[2] = g;

	// f4
	for (i=0, g=0; i<n; i++) {
		t = x[i] - 0.5*sin(1000*PI*(i+1)/n) - 0.5;
		g += t*t;
	}
	f[3] = g;
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* FES1_sample ();
static Matrix_t* FES2_sample ();
static Matrix_t* FES3_sample ();

Matrix_t* FES_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	FES1_sample ();
		case 2:
			return	FES2_sample ();
		case 3:
			return	FES3_sample ();
		default:
                       fprintf (stderr, "FES%d_sample have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* FES1_sample () {
	Matrix_t *sample = Matrix_new (1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	
	return sample;
}
static Matrix_t* FES2_sample () {
	Matrix_t *sample = Matrix_new (1, 3);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	sample->elements[2] = 0;
	
	return sample;
}
static Matrix_t* FES3_sample () {
	Matrix_t *sample = Matrix_new (1, 4);

	sample->elements[0] = 0;
	sample->elements[1] = 0;
	sample->elements[2] = 0;
	sample->elements[3] = 0;
	
	return sample;
}
