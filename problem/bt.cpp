#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bt.h"
#include "dtlz.h"
#include "lsmop.h"
#include "matrix.h"
#include "shape.h"
#include "uf.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  bt1 (double *x, int n, double *f, int m);
static void  bt2 (double *x, int n, double *f, int m);
static void  bt3 (double *x, int n, double *f, int m);
static void  bt4 (double *x, int n, double *f, int m);
static void  bt5 (double *x, int n, double *f, int m);
static void  bt6 (double *x, int n, double *f, int m);
static void  bt7 (double *x, int n, double *f, int m);
static void  bt8 (double *x, int n, double *f, int m);
static void  bt9 (double *x, int n, double *f, int m);

Problem_t* BT_new (char *title, int numObj, int numVar) {
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
	if (!strcmp (title, "BT9")) {
		problem->numObj = 3; 
	}
	
	// setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = 0.0;
		uppBound[i] = 1.0;
	}
	if (!strcmp (title, "BT7")) {
		for (i=1; i<numVar; i++) {
			lowBound[i] = -1.0;
			uppBound[i] = 1.0;
		}
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "BT1")) {
		P = problem; problem->evaluate = bt1;
	} else if (!strcmp (title, "BT2")) {
		P = problem; problem->evaluate = bt2;
	} else if (!strcmp (title, "BT3")) {
		P = problem; problem->evaluate = bt3;
	} else if (!strcmp (title, "BT4")) {
		P = problem; problem->evaluate = bt4;
	} else if (!strcmp (title, "BT5")) {
		P = problem; problem->evaluate = bt5;
	} else if (!strcmp (title, "BT6")) {
		P = problem; problem->evaluate = bt6;
	} else if (!strcmp (title, "BT7")) {
		P = problem; problem->evaluate = bt7;
	} else if (!strcmp (title, "BT8")) {
		P = problem; problem->evaluate = bt8;
	} else if (!strcmp (title, "BT9")) {
		P = problem; problem->evaluate = bt9;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static double D1(double g, double theta) {
	return g*g + (1 - exp(-g*g/theta)) / 5.0;
}

static double D2(double g, double theta) {
	return g*g + pow(fabs(g), theta) / 5.0;
}

static double S1(double x, double gama) {
	return pow(fabs(x), gama);
}

static double S2(double x, double gama) {
	if (x < 0.25)
		return (1 - pow(1-4.0*x, gama))/4;
	if (x < 0.50)
		return (1 + pow(4.0*x-1, gama))/4;
	if (x < 0.75)
		return (3 - pow(3.0-4.0*x, gama))/4;

	return (3 + pow(4.0*x-3.0, gama))/4;
}

static void  bt1 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0, y;

	// f1
	for (i=1, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D1 (y, 1.0e-10);
		g += t; 
	}
	f[0] = x[0] + g;
	
	// f2
	for (i=2, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D1 (y, 1.0e-10);
		g += t; 
	}
	f[1] = 1 - sqrt(x[0]) + g;
}

static void  bt2 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0, y;

	// f1
	for (i=1, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D2 (y, 0.2);
		g += t; 
	}
	f[0] = x[0] + g;
	
	// f2
	for (i=2, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D2 (y, 0.2);
		g += t; 
	}
	f[1] = 1 - sqrt(x[0]) + g;
}

static void  bt3 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0, y;

	// f1
	for (i=1, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D1 (y, 1.0e-8);
		g += t; 
	}
	f[0] = S1(x[0],0.02) + g;
	
	// f2
	for (i=2, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D1 (y, 1.0e-8);
		g += t; 
	}
	f[1] = 1 - sqrt(S1(x[0],0.02)) + g;
}

static void  bt4 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0, y;

	// f1
	for (i=1, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D1 (y, 1.0e-8);
		g += t; 
	}
	f[0] = S2(x[0],0.06) + g;
	
	// f2
	for (i=2, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D1 (y, 1.0e-8);
		g += t; 
	}
	f[1] = 1 - sqrt(S2(x[0],0.06)) + g;

}

static void  bt5 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0, y;

	// f1
	for (i=1, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D1 (y, 1.0e-10);
		g += t; 
	}
	f[0] = x[0] + g;
	
	// f2
	for (i=2, g=0; i < n; i += 2) {
		y = x[i] - sin ((i+1)*PI/(2*n)); 
		t = D1 (y, 1.0e-10);
		g += t; 
	}
	f[1] = (1 - x[0])*(1 - x[0]*sin(8.5*PI*x[0])) + g;
}

static void  bt6 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0, y;

	// f1
	for (i=1, g=0; i < n; i += 2) {
 		y = x[i] - pow(x[0], 0.5+1.5*(i+1-1.0)/(n-1));
		t = D1 (y, 1.0e-4);
		g += t; 
	}
	f[0] = x[0] + g;
	
	// f2
	for (i=2, g=0; i < n; i += 2) {
 		y = x[i] - pow(x[0], 0.5+1.5*(i+1-1.0)/(n-1));
		t = D1 (y, 1.0e-4);
		g += t; 
	}
	f[1] = 1 - sqrt(x[0]) + g;
}

static void  bt7 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0, y;

	// f1
	for (i=1, g=0; i < n; i += 2) {
		y = x[i] - sin(6*PI*x[0]);
		t = D1 (y, 1.0e-3);
		g += t; 
	}
	f[0] = x[0] + g;
	
	// f2
	for (i=2, g=0; i < n; i += 2) {
		y = x[i] - sin(6*PI*x[0]);
		t = D1 (y, 1.0e-3);
		g += t; 
	}
	f[1] = 1 - sqrt(x[0]) + g;
}

static void  bt8 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0, y;

	// f1
	for (i=1, g=0; i < n; i += 2) {
 		y = x[i] - pow(x[0], 0.5+1.5*(i+1-1.0)/(n-1));
		t = D1 (y, 1.0e-3);
		g += 4*t*t - cos(8*PI*t) + 1; 
	}
	f[0] = x[0] + g;
	
	// f2
	for (i=2, g=0; i < n; i += 2) {
 		y = x[i] - pow(x[0], 0.5+1.5*(i+1-1.0)/(n-1));
		t = D1 (y, 1.0e-3);
		g += 4*t*t - cos(8*PI*t) + 1; 
	}
	f[1] = 1 - sqrt(x[0]) + g;

}

static void  bt9 (double *x, int n, double *f, int m) {
	int 	i;
	double 	g = 0, t = 0, y;

	// f1
	for (i=2, g=0; i < n; i+=3) {
		y = x[i] - sin((i+1)*PI/(2*n));
		t = D1(y, 1.0e-9);
		g += t; 
	}
	f[0] = cos(0.5*x[0]*PI)*cos(0.5*x[1]*PI) + 10*g; 
	
	// f2
	for (i=3, g=0; i < n; i+=3) {
		y = x[i] - sin((i+1)*PI/(2*n));
		t = D1(y, 1.0e-9);
		g += t; 
	}
	f[1] = cos(0.5*x[0]*PI)*sin(0.5*x[1]*PI) + 10*g; 
	
	// f3
	for (i=4, g=0; i < n; i+=3) {
		y = x[i] - sin((i+1)*PI/(2*n));
		t = D1(y, 1.0e-9);
		g += t; 
	}
	f[2] = sin(0.5*x[0]*PI) + 10*g; 
}

//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* BT1_sample ();
static Matrix_t* BT5_sample ();
static Matrix_t* BT9_sample ();

Matrix_t* BT_sample (int No, int numObj) {
	switch (No) {
		case 1:
		case 2:
		case 3:
		case 4:
		case 6:
		case 7:
		case 8:
			return	BT1_sample ();
		case 5:
			return	BT5_sample ();
		case 9:
			return	BT9_sample ();
		default:
                       fprintf (stderr, "BT%d_sample have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* BT1_sample () {
	return 	UF_sample (1, 2);
}

static Matrix_t* BT5_sample () {
	int 	  i;
	double    d = 1.0/NUM_SAMPLE, f1, f2, y1;
	Matrix_t* sample = Matrix_new (NUM_SAMPLE+1, 2);
	Matrix_t* T = NULL;

	for (i=0; i<=NUM_SAMPLE; i++) {
		y1 = i * d;
		f1 = y1;
		f2 = (1.0 - f1) * (1.0 - f1*sin(8.5*PI*f1));
		sample->elements[2*i] = f1;
		sample->elements[2*i+1] = f2;
	}
	T  = sample;
	sample = Matrix_front (T);
	Matrix_free (&T);
	
	return sample;
}

static Matrix_t* BT9_sample () {
	return DTLZ_sample (2, 3);
}
