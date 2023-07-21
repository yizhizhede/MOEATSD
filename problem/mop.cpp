#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mop.h"
#include "matrix.h"
#include "shape.h"
#include "problem.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// 
static void  mop1 (double *x, int n, double *f, int m);
static void  mop2 (double *x, int n, double *f, int m);
static void  mop3 (double *x, int n, double *f, int m);
static void  mop4 (double *x, int n, double *f, int m);
static void  mop5 (double *x, int n, double *f, int m);
static void  mop6 (double *x, int n, double *f, int m);
static void  mop7 (double *x, int n, double *f, int m);

// genereate a problem
Problem_t *MOP_new (char *title, int numObj, int numVar) {
	// common variable
	int i;
	size_t size;
	double *lowBound = NULL, *uppBound = NULL;

	// allocating memory for a problem
        Problem_t *problem = (Problem_t *)malloc (sizeof (Problem_t));
        if (problem == NULL) {
        	fprintf (stderr, "Allocating memory failed\n");         
               	exit (-1);
        }

	// set title 
	strcpy (problem->title, title);

	// 
	if (!strcmp (title, "MOP1")) {
		// set Obj. & Var.
		problem->numObj = 2; 
		problem->numVar = 1;

		// setting the bound of varibles
		size = numVar * sizeof (double);
		lowBound = (double *)malloc (size);
		uppBound = (double *)malloc (size);
		for (i=0; i<numVar; i++) {
			lowBound[i] = -1.0e+5;
			uppBound[i] = 1.0e+5;
		}
		problem->lowBound = lowBound;
		problem->uppBound = uppBound;

		// set evaluate function
		problem->evaluate = mop1;
	} else if (!strcmp (title, "MOP2")) {
		// set Obj. & Var.
		problem->numObj = 2; 
		problem->numVar = numVar;

		// setting the bound of varibles
		size = numVar * sizeof (double);
		lowBound = (double *)malloc (size);
		uppBound = (double *)malloc (size);
		for (i=0; i<numVar; i++) {
			lowBound[i] = -4.0;
			uppBound[i] = 4.0;
		}
		problem->lowBound = lowBound;
		problem->uppBound = uppBound;

		// set evaluate function
		problem->evaluate = mop2;
	} else if (!strcmp (title, "MOP3")) {
		// set Obj. & Var.
		problem->numObj = 2; 
		problem->numVar = 2;

		// setting the bound of varibles
		size = numVar * sizeof (double);
		lowBound = (double *)malloc (size);
		uppBound = (double *)malloc (size);
		for (i=0; i<numVar; i++) {
			lowBound[i] = -PI;
			uppBound[i] = PI;
		}
		problem->lowBound = lowBound;
		problem->uppBound = uppBound;

		// set evaluate function
		problem->evaluate = mop3;
	} else if (!strcmp (title, "MOP4")) {
		// set Obj. & Var.
		problem->numObj = 2; 
		problem->numVar = 3;

		// setting the bound of varibles
		size = numVar * sizeof (double);
		lowBound = (double *)malloc (size);
		uppBound = (double *)malloc (size);
		for (i=0; i<numVar; i++) {
			lowBound[i] = -5.0;
			uppBound[i] = 5.0;
		}
		problem->lowBound = lowBound;
		problem->uppBound = uppBound;

		// set evaluate function
		problem->evaluate = mop4;
	} else if (!strcmp (title, "MOP5")) {
		// set Obj. & Var.
		problem->numObj = 3; 
		problem->numVar = 2;

		// setting the bound of varibles
		size = numVar * sizeof (double);
		lowBound = (double *)malloc (size);
		uppBound = (double *)malloc (size);
		for (i=0; i<numVar; i++) {
			lowBound[i] = -30.0;
			uppBound[i] = 30.0;
		}
		problem->lowBound = lowBound;
		problem->uppBound = uppBound;

		// set evaluate function
		problem->evaluate = mop5;
	} else if (!strcmp (title, "MOP6")) {
		// set Obj. & Var.
		problem->numObj = 2; 
		problem->numVar = 2;

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

		// set evaluate function
		problem->evaluate = mop6;
	} else if (!strcmp (title, "MOP7")) {
		// set Obj. & Var.
		problem->numObj = 3; 
		problem->numVar = 2;

		// setting the bound of varibles
		size = numVar * sizeof (double);
		lowBound = (double *)malloc (size);
		uppBound = (double *)malloc (size);
		for (i=0; i<numVar; i++) {
			lowBound[i] = -400.0;
			uppBound[i] = 400.;
		}
		problem->lowBound = lowBound;
		problem->uppBound = uppBound;

		// set evaluate function
		problem->evaluate = mop7;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void  mop1 (double *x, int n, double *f, int m) {
	f[0] = x[0] * x[0];
	f[1] = (x[0]-2.0) * (x[0]-2.0);
}

static void  mop2 (double *x, int n, double *f, int m) {
	int i;
	double s = 0;

	// f1
	for (i=0, s=0; i<n; i++) {
		s += (x[i] - 1.0/sqrt(n))*(x[i] - 1.0/sqrt(n));
	}
	f[0] = 1.0 - exp (-s);
	
	// f2
	for (i=0, s=0; i<n; i++) {
		s += (x[i] + 1.0/sqrt(n))*(x[i] + 1.0/sqrt(n));
	}
	f[1] = 1.0 - exp (-s);
}

static void  mop3 (double *x, int n, double *f, int m) {
	double A1, A2, B1, B2;

	A1 = 0.5*sin(1.0) - 2.0*cos(1.0) + sin(2.0) - 1.5*cos(2.0);
	A2 = 1.5*sin(1.0) - cos(1.0) + 2.0*sin(2.0) - 0.5*cos(2.0);

	B1 = 0.5*sin(x[0]) - 2.0*cos(x[0]) + sin(x[1]) - 1.5*cos(x[1]);
	B2 = 1.5*sin(x[0]) - cos(x[0]) + 2.0*sin(x[1]) - 0.5*cos(x[1]);

	f[0] = 1.0 + (A1-B1)*(A1-B1) + (A2-B2)*(A2-B2);
	f[1] = (x[0]+3)*(x[0]+3)+(x[1]+1)*(x[1]+1);
}

static void  mop4 (double *x, int n, double *f, int m) {
	int i;
	double s;

	for (i=0, s=0; i<2; i++) {
		s += -10.0*exp (-0.2*sqrt(x[i]*x[i]+x[i+1]*x[i+1]));
	}
	f[0] = s;

	for (i=0, s=0; i<3; i++) {
		s += pow(fabs(x[i]),0.8)+5.0*sin(pow(x[i],3));
	}
	f[1] = s;
}

static void  mop5 (double *x, int n, double *f, int m) {
	f[0] = 0.5*(x[0]*x[0]+x[1]*x[1]) + sin(x[0]*x[0]+x[1]*x[1]);
	f[1] = (3.0*x[0]-2.0*x[1]+4.0)*(3.0*x[0]-2.0*x[1]+4.0) / 8.0 + (x[0]-x[1]+1.0)*(x[0]-x[1]+1.0)/27.0 + 15.0;
	f[2] = 1.0 / (x[0]*x[0]+x[1]*x[1]+1.0) - 1.1*exp(-x[0]*x[0]-x[1]*x[1]);
}

static void  mop6 (double *x, int n, double *f, int m) {
	f[0] = x[0];
	f[1] = (1.0+10.0*x[1])*(1.0-(x[0]/(1.0+10*x[1]))*(x[0]/(1.0+10*x[1]))-x[0]/(1.0+10*x[1])*sin(8.0*PI*x[0]));
}

static void  mop7 (double *x, int n, double *f, int m) {
	f[0] = (x[0]-2)*(x[0]-2)/2 + (x[1]+1)*(x[1]+1)/13 + 3;
	f[1] = (x[0]+x[1]-3)*(x[0]+x[1]-3)/36 + (-x[0]+x[1]+2)*(-x[0]+x[1]+2)/8 - 17;
	f[2] = (x[0]+2*x[1]-1)*(x[0]+2*x[1]-1)/175 + (-x[0]+2*x[1])*(-x[0]+2*x[1])/17 - 13;
}

//**************************************************************************************************************
//**************************************************************************************************************
Matrix_t  *MOP1_sample ();
Matrix_t  *MOP2_sample ();
Matrix_t  *MOP3_sample ();
Matrix_t  *MOP4_sample ();
Matrix_t  *MOP5_sample ();
Matrix_t  *MOP6_sample ();
Matrix_t  *MOP7_sample ();

Matrix_t  *MOP_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return  MOP1_sample ();	
		case 2:
			return	MOP2_sample ();
		case 3:
			return	MOP3_sample ();
		case 4:
			return	MOP4_sample ();
		case 5:
			return	MOP5_sample ();
		case 6:
			return	MOP6_sample ();
		case 7:
			return	MOP7_sample ();
		default:
                       fprintf (stderr, "MOP%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

Matrix_t  *MOP1_sample () {
	Matrix_t* M = Matrix_new (NUM_SAMPLE+1, 2);
	int 	  i;
	double 	  d = 2.0/NUM_SAMPLE;
	double 	  f1, f2;

	for (i=0; i<=NUM_SAMPLE; i++) {
		f1 = (i*d)*(i*d);
		f2 = (i*d - 2.0)*(i*d - 2.0);
		M->elements[i*2+0] = f1;
		M->elements[i*2+1] = f2;
	}

	return M;
}

Matrix_t  *MOP2_sample () {
	return NULL;
}

Matrix_t  *MOP3_sample () {
	return NULL;
}

Matrix_t  *MOP4_sample () {
	return NULL;
}

Matrix_t  *MOP5_sample () {
	return NULL;
}

Matrix_t  *MOP6_sample () {
	return NULL;
}

Matrix_t  *MOP7_sample () {
	return NULL;
}
