#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "zdt.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  zdt1 (double *x, int n, double *f, int m);
static void  zdt2 (double *x, int n, double *f, int m);
static void  zdt3 (double *x, int n, double *f, int m);
static void  zdt4 (double *x, int n, double *f, int m);
static void  zdt6 (double *x, int n, double *f, int m);

Problem_t* ZDT_new (char *title, int numObj, int numVar) {
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
		lowBound[i] = 0.0;
		uppBound[i] = 1.0;
	}
	if (!strcmp (title, "ZDT4")) {
		for (i=1; i<numVar; i++) {
			lowBound[i] = -0.5;
			uppBound[i] = 5.0;
		}
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "ZDT1")) {
		P = problem; problem->evaluate = zdt1;
	} else if (!strcmp (title, "ZDT2")) {
		P = problem; problem->evaluate = zdt2;
	} else if (!strcmp (title, "ZDT3")) {
		P = problem; problem->evaluate = zdt3;
	} else if (!strcmp (title, "ZDT4")) {
		P = problem; problem->evaluate = zdt4;
	} else if (!strcmp (title, "ZDT6")) {
		P = problem; problem->evaluate = zdt6;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}


	return problem;
}

static void  zdt1 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, h = 0;

	// f1
	f[0] = x[0];
	
	// g
	for (i=1; i<n; i++) {
		g += x[i];
	}
	g = 1 + 9.0*g/(n-1);
	
	// h
	h = 1 - sqrt (f[0]/g);

	// f2
	f[1] = g*h;
}

static void  zdt2 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, h = 0;

	// f1
	f[0] = x[0];
	
	// g
	for (i=1; i<n; i++) {
		g += x[i];
	}
	g = 1 + 9.0*g/(n-1);
	
	// h
	h = 1 - (f[0]/g)*(f[0]/g);

	// f2
	f[1] = g*h;
}

static void  zdt3 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, h = 0;

	// f1
	f[0] = x[0];
	
	// g
	for (i=1; i<n; i++) {
		g += x[i];
	}
	g = 1 + 9.0*g/(n-1);
	
	// h
	h = 1 - sqrt (f[0]/g)-(f[0]/g)*sin(10*PI*f[0]);

	// f2
	f[1] = g*h;
}

static void  zdt4 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, h = 0;

	// f1
	f[0] = x[0];
	
	// g
	for (i=1; i<n; i++) {
		g += (x[i]*x[i]-10*cos(4.0*PI*x[i]));
	}
	g = 1 + 10*(n-1)+ g;
	
	// h
	h = 1 - sqrt (f[0]/g);

	// f2
	f[1] = g*h;
}

static void  zdt6 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, h = 0;

	// f1
	f[0] = 1.0 - exp(-4.0*x[0])*pow(sin(6.0*PI*x[0]), 6.0);
	
	// g
	for (i=1; i<n; i++) {
		g += x[i];
	}
	g = 1 + 9.0*pow(g/(n-1), 0.25);
	
	// h
	h = 1 - (f[0]/g)*(f[0]/g);

	// f2
	f[1] = g*h;
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* ZDT1_sample ();
static Matrix_t* ZDT2_sample ();
static Matrix_t* ZDT3_sample ();
static Matrix_t* ZDT6_sample ();

Matrix_t* ZDT_sample (int No, int numObj) {
	switch (No) {
		case 1:
		case 4:
			return	ZDT1_sample ();
		case 2:
			return	ZDT2_sample ();
		case 3:
			return	ZDT3_sample ();
		case 6:
			return	ZDT6_sample ();
		default:
                       fprintf (stderr, "ZDT%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* ZDT1_sample () {
	int 	  i;
	double    d = 1.0/NUM_SAMPLE, f1, f2, y1;
	Matrix_t* sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		y1 = i * d;
		if (i % 2 == 0) {
			f1 = y1;
			f2 = 1.0 - sqrt (f1);
		} else {
			f2 = y1;
			f1 = (1-f2)*(1-f2);
		}
		sample->elements[2*i] = f1;
		sample->elements[2*i+1] = f2;
	}
	
	return sample;
}
static Matrix_t* ZDT2_sample () {
	int       i;
	double 	  d = 1.0/NUM_SAMPLE, f1, f2, y1;
	Matrix_t* sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		y1 = i * d;
		if (i%2==0) {
			f1 = y1;
			f2 = 1.0 - f1*f1;
		} else {
			f2 = y1;
			f1 = sqrt (1.0-f2);
		}
		sample->elements[2*i] = f1;
		sample->elements[2*i+1] = f2;
	}
	
	return sample;
}

static Matrix_t* ZDT3_sample () {
	int 	  i;
	double    d = 1.0/NUM_SAMPLE, f1, f2, y1;
	Matrix_t* sample = Matrix_new (NUM_SAMPLE+1, 2);
	Matrix_t* front = NULL;

	for (i=0; i<=NUM_SAMPLE; i++) {
		y1 = i * d;
		f1 = y1;
		f2 = 1.0 - sqrt (f1) - f1*sin(10*PI*f1);
		sample->elements[2*i] = f1;
		sample->elements[2*i+1] = f2;
	}
	
	front = Matrix_front (sample);
	Matrix_free (&sample);

	return front;
}
static Matrix_t* ZDT6_sample () {
	int 		i;
	double 		d = 1.0/NUM_SAMPLE, f1, f2, y1;
	Matrix_t*	sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		y1 = i * d;
		f1 = 1.0 - exp(-4.0*y1)*pow(sin(6*PI*y1), 6);
		f2 = 1.0 - f1*f1;
		sample->elements[2*i] = f1;
		sample->elements[2*i+1] = f2;
	}
	
	return sample;
}
