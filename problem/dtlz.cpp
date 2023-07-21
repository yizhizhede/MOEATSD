#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "dtlz.h"
#include "matrix.h"
#include "shape.h"
#include "algebra.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  dtlz1 (double *x, int n, double *f, int m);
static void  dtlz2 (double *x, int n, double *f, int m);
static void  dtlz3 (double *x, int n, double *f, int m);
static void  dtlz4 (double *x, int n, double *f, int m);
static void  dtlz5 (double *x, int n, double *f, int m);
static void  dtlz6 (double *x, int n, double *f, int m);
static void  dtlz7 (double *x, int n, double *f, int m);

Problem_t *DTLZ_new (char *title, int numObj, int numVar) {
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
	strcpy (problem->title, title);
	problem->numObj = numObj; 

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
	if (!strcmp (title, "DTLZ1")) {
		P = problem; problem->evaluate = dtlz1;
	} else if (!strcmp (title, "DTLZ2")) {
		P = problem; problem->evaluate = dtlz2;
	} else if (!strcmp (title, "DTLZ3")) {
		P = problem; problem->evaluate = dtlz3;
	} else if (!strcmp (title, "DTLZ4")) {
		P = problem; problem->evaluate = dtlz4;
	} else if (!strcmp (title, "DTLZ5")) {
		P = problem; problem->evaluate = dtlz5;
	} else if (!strcmp (title, "DTLZ6")) {
		P = problem; problem->evaluate = dtlz6;
	} else if (!strcmp (title, "DTLZ7")) {
		P = problem; problem->evaluate = dtlz7;
	} else {
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}

static void dtlz1(double *x,const int n,double *f,const int m) {
	int 	k = n-m +1;	  // recommanded value k = 5
	double 	g = 0.0;
	int 	i, j;
	
	// compute g
	for (i=m-1; i<n; i++) {
		g += ((x[i]-0.5)*(x[i]-0.5)-cos(20*PI*(x[i]-0.5)));			
	}
	g = 100*(k + g);

	// f
	for (j=0; j<m; j++) {
		f[j] = 0.5*(1+g);
		for (i=0; i<m-j-1; i++) {
			f[j] *= x[i];		
		}
		if (j>0) {
			f[j] *= (1-x[i]);
		}
	}
}

static void dtlz2(double *x,const int n,double* f,const int m) {
	int 	i, j;
	double  g = 0.0;

	// compute g
	for (i=m-1; i<n; i++) {
		g += (x[i]-0.5)*(x[i]-0.5);		
	}

	// compute f
	for (j=0; j<m; j++) {
		f[j] = 1.0 + g;
		for (i=0; i<m-j-1; i++) {
			f[j] *= cos (0.5*PI*x[i]);
		}
		if (j>0) {
			f[j] *= sin(0.5*PI*x[i]);
		}
	}
}

static void dtlz3(double *x,const int n,double *f,const int m) {
	int 	i, j, k;
	double 	g=0.0;
	k = n-m+1;	// recommanded value k = 10
	
	// compute g
	for (i=m-1; i<n; i++) {
		g += ((x[i]-0.5)*(x[i]-0.5) - cos(20*PI*(x[i]-0.5)));			
	}
	g = 100*(k+g);

	// compute f
	for (j=0; j<m; j++) {
		f[j] = 1.0+g;
		for (i=0; i<m-j-1; i++) {
			f[j] *= cos(0.5*PI*x[i]);
		}
		if (j>0) {
			f[j] *= sin(0.5*PI*x[i]);
		}
	}
}

static void  dtlz4(double *x,const int n,double *f,const int m) {
	int 	i, j;
	double 	g = 0.0;
	double 	alpha = 100.0;

	// compute g
	for (i=m-1; i<n; i++) {
		g += (x[i]-0.5)*(x[i]-0.5);		
	}

	// compute f
	for (j=0; j<m; j++) {
		f[j] = 1.0 + g;
		for (i=0; i<m-j-1; i++) {
			f[j] *= cos (0.5*PI*pow(x[i],alpha));
		}
		if (j>0) {
			f[j] *= sin (0.5*PI*pow(x[i],alpha));
		}
	}
}

static void  dtlz5(double *x,const int n,double *f,const int m) {
	int 	i, j;
	double 	g = 0.0;

	// compute g
	for (i=m-1; i<n; i++) {
		g += (x[i]-0.5)*(x[i]-0.5);		
	}

	// compute f
	for (j=0; j<m; j++) {
		f[j] = 1.0 + g;
		for (i=0; i<m-j-1; i++) {
			if (0 == i) {
			 	f[j] *= cos(0.5*PI*x[i]);
			} else {
				f[j] *= cos((0.25*pi/(1+g))*(1+2*g*x[i]));
			}
		}
		if (j>0) {
			if (0 == i) {
			 	f[j] *= sin(0.5*PI*x[i]);
			} else {
				f[j] *= sin((0.25*pi/(1+g))*(1+2*g*x[i]));
			}
		}
	}
}

static void  dtlz6(double *x,const int n,double *f,const int m) {
	int 	i, j;
	double 	g = 0.0;

	// compute g
	for (i=m-1; i<n; i++) {
	 	g += pow(x[i], 0.1);				 	
	}

	// compute f
	for (j=0; j<m; j++) {
		f[j] = 1.0 + g;
		for (i=0; i<m-j-1; i++) {
			if (0 == i) {
			 	f[j] *= cos(0.5*PI*x[i]);
			} else {
				f[j] *= cos((0.25*pi/(1+g))*(1+2*g*x[i]));
			}
		}
		if (j>0) {
			if (0 == i) {
			 	f[j] *= sin(0.5*PI*x[i]);
			} else {
				f[j] *= sin((0.25*pi/(1+g))*(1+2*g*x[i]));
			}
		}
	}
}


static void dtlz7(double *x,const int n,double *f,const int m) {
	int 	i, j;
	double 	g=0.0, h=0.0;

	// compute f[0,1...m-2]
	for (j=0; j<m-1; j++) {
		f[j] = x[j];
	}

	// compute g
	for (i=m-1; i<n; i++) {
		g += x[i]; 								
	}
	g = 1+9.0*g/(n-m+1);

	// compute h
	for (i=0; i<m-1; i++) {
		h += (f[i]/(1+g))*(1+sin(3*pi*f[i]));
	}
	h = m - h;
	
	// compute f[m-1]
	f[m-1]=(1+g)*h;
}

//**************************************************************************************************************
//**************************************************************************************************************
//**************************************************************************************************************
Matrix_t *DTLZ_sample (int No, int numObj) {
	int H;
	Matrix_t *sample=NULL; 
	Matrix_t *half=NULL;

        switch (No) {
         	case 1:
                	for (H=2; conb (H+numObj-1, numObj-1) < NUM_SAMPLE; H++){};
	                sample = H1_sample (numObj, H);
	                half = scale_product (0.5, sample);
	                Matrix_free (&sample);
                	return half;
             	case 2:
                case 3:
                case 4:
                      	for (H=1; conb (H+numObj-1, numObj-1) < NUM_SAMPLE; H++){};
                    	return H2_sample (numObj, H);
                case 5:
                case 6:
			return H4_sample (numObj, NUM_SAMPLE);			
                case 7:
                      	for (H=1; grid (numObj-1, H) < NUM_SAMPLE; H++){};
                       	return H3_sample (numObj, H);
                default:
                       fprintf (stderr, "DTLZ%d have been not on consideration\n", No);
                       exit (0);
       	}
        return NULL;
}
