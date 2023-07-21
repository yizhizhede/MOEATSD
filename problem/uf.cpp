#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "uf.h"
#include "dtlz.h"
#include "lsmop.h"
#include "matrix.h"
#include "shape.h"

#define PI 3.14159265358979323846264338327950288419716939937510
#define pi 3.14159265358979323846264338327950288419716939937510
#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  uf1 (double *x, int n, double *f, int m);
static void  uf2 (double *x, int n, double *f, int m);
static void  uf3 (double *x, int n, double *f, int m);
static void  uf4 (double *x, int n, double *f, int m);
static void  uf5 (double *x, int n, double *f, int m);
static void  uf6 (double *x, int n, double *f, int m);
static void  uf7 (double *x, int n, double *f, int m);
static void  uf8 (double *x, int n, double *f, int m);
static void  uf9 (double *x, int n, double *f, int m);
static void  uf10 (double *x, int n, double *f, int m);

Problem_t* UF_new (char *title, int numObj, int numVar) {
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
	if (!strcmp (title, "UF8") ||
	    !strcmp (title, "UF9") ||
	    !strcmp (title, "UF10")) {
		problem->numObj = 3; 
	}
	
	// setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	lowBound[0] = 0.0;
	uppBound[0] = 1.0;
	for (i=1; i<numVar; i++) {
		lowBound[i] = -1.0;
		uppBound[i] = 1.0;
	}
	if (!strcmp (title, "UF3")) {
		for (i=1; i<numVar; i++) {
			lowBound[i] = 0;
			uppBound[i] = 1.0;
		}
	} else if (!strcmp (title, "UF4")) {
		for (i=1; i<numVar; i++) {
			lowBound[i] = -2.0;
			uppBound[i] = 2.0;
		}
	} else if (!strcmp (title, "UF8") || 
	 	   !strcmp (title, "UF9") ||
	 	   !strcmp (title, "UF10")) {
		lowBound[1] = 0;
		uppBound[1] = 1;
		for (i=2; i<numVar; i++) {
			lowBound[i] = -2.0;
			uppBound[i] = 2.0;
		}
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "UF1")) {
		P = problem; problem->evaluate = uf1;
	} else if (!strcmp (title, "UF2")) {
		P = problem; problem->evaluate = uf2;
	} else if (!strcmp (title, "UF3")) {
		P = problem; problem->evaluate = uf3;
	} else if (!strcmp (title, "UF4")) {
		P = problem; problem->evaluate = uf4;
	} else if (!strcmp (title, "UF5")) {
		P = problem; problem->evaluate = uf5;
	} else if (!strcmp (title, "UF6")) {
		P = problem; problem->evaluate = uf6;
	} else if (!strcmp (title, "UF7")) {
		P = problem; problem->evaluate = uf7;
	} else if (!strcmp (title, "UF8")) {
		P = problem; problem->evaluate = uf8;
	} else if (!strcmp (title, "UF9")) {
		P = problem; problem->evaluate = uf9;
	} else if (!strcmp (title, "UF10")) {
		P = problem; problem->evaluate = uf10;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}


	return problem;
}

static void  uf1 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0;
	double 	g = 0, t = 0;

	// f1
	for (i=2, g=0; i < n; i += 2) {
		J1++;
		t = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		g += t * t; 
	}
	f[0] = x[0] + 2 * g / J1;
	
	// f2
	for (i=1, g=0; i < n; i += 2) {
		J2++;
		t = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		g += t * t; 
	}
	f[1] = 1 - sqrt(x[0]) + 2 * g / J2;
	
}

static void  uf2 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0;
	double 	g = 0, t = 0;

	// f1
	for (i=2, g=0; i < n; i += 2) {
		J1++;
		t = x[i] -  (0.3*x[0]*x[0]*cos(24*PI*x[0]+4*(i+1)*PI/n) + 0.6*x[0])*cos(6.0*PI*x[0]+ (i+1)*PI/n);
		g += t * t; 
	}
	f[0] = x[0] + 2 * g / J1;
	
	// f2
	for (i=1, g=0; i < n; i += 2) {
		J2++;
		t = x[i] -  (0.3*x[0]*x[0]*cos(24*PI*x[0]+4*(i+1)*PI/n) + 0.6*x[0])*sin(6.0*PI*x[0]+ (i+1)*PI/n);
		g += t * t; 
	}
	f[1] = 1 - sqrt(x[0]) + 2 * g / J2;
}

static void  uf3 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0;
	double 	t = 0, yi;
	double 	sum, prod;

	// f1
	for (i=2, sum=0, prod=1; i < n; i += 2) {
		J1++;
		yi = x[i] - pow(x[0],0.5*(1.0+3.0*((i+1)-2.0)/(n-2.0)));	
		t  = cos(20.0*yi*PI/sqrt(i+1));
		sum += yi*yi;
		prod *= t;
	}
	f[0] = x[0] + 2 * (4*sum - 2*prod + 2) / J1;
	
	// f2
	for (i=1, sum=0, prod=1; i < n; i += 2) {
		J2++;
		yi = x[i] - pow(x[0],0.5*(1.0+3.0*((i+1)-2.0)/(n-2.0)));	
		t  = cos(20.0*yi*PI/sqrt(i+1));
		sum += yi*yi;
		prod *= t;
	}
	f[1] = 1 - sqrt(x[0]) + 2 * (4*sum - 2*prod + 2) / J2;
}

static void  uf4 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0;
	double 	g = 0, t = 0, y, h;

	// f1
	for (i=2, g=0; i < n; i += 2) {
		J1++;
		y = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		t = y > 0 ? y : -y;
		h = t / (1 + exp(2*t));
		g += h; 
	}
	f[0] = x[0] + 2 * g / J1;
	
	// f2
	for (i=1, g=0; i < n; i += 2) {
		J2++;
		y = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		t = y > 0 ? y : -y;
		h = t / (1 + exp(2*t));
		g += h; 
	}
	f[1] = 1 - x[0]*x[0] + 2 * g / J2;
}

static void  uf5 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0;
	double 	g = 0, t = 0, y, h;
	int	N = 10;
	double 	epsilon = 0.1;

	// f1
	for (i=2, g=0; i < n; i += 2) {
		J1++;
		y = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		t = y;
		h = 2*t*t - cos(4*PI*t) + 1;
		g += h; 
	}
	f[0] = x[0] + (0.5/N + epsilon)*fabs(sin(2.0*N*PI*x[0])) + 2 * g / J1;
	
	// f2
	for (i=1, g=0; i < n; i += 2) {
		J2++;
		y = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		t = y;
		h = 2*t*t - cos(4*PI*t) + 1;
		g += h; 
	}
	f[1] = 1 - x[0] + (0.5/N + epsilon)*fabs(sin(2.0*N*PI*x[0])) + 2 * g / J1 + 2 * g / J2;
	
}

static void  uf6 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0;
	double 	t = 0, y;
	double 	sum, prod;
	int	N = 2;
	double 	epsilon = 0.1;

	// f1
	for (i=2, sum=0, prod=1; i < n; i += 2) {
		J1++;
		y = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		t  = cos(20.0*y*PI/sqrt(i+1));
		sum += y*y;
		prod *= t;
	}
	t = 2.0*(0.5/N + epsilon)*sin(2.0*N*PI*x[0]);
	t = t > 0 ? t : 0;
	f[0] = x[0] + t + 2 * (4*sum - 2*prod + 2) / J1;
	
	// f2
	for (i=1, sum=0, prod=1; i < n; i += 2) {
		J2++;
		y = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		t  = cos(20.0*y*PI/sqrt(i+1));
		sum += y*y;
		prod *= t;
	}
	t = 2.0*(0.5/N + epsilon)*sin(2.0*N*PI*x[0]);
	t = t > 0 ? t : 0;
	f[1] = 1 - x[0] + t + 2 * (4*sum - 2*prod + 2) / J2;
}

static void  uf7 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0;
	double 	g = 0, y;

	// f1
	for (i=2, g=0; i < n; i += 2) {
		J1++;
		y = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		g += y*y; 
	}
	f[0] = pow(x[0], 0.2) + 2 * g / J1;
	
	// f2
	for (i=1, g=0; i < n; i += 2) {
		J2++;
		y = x[i] - sin(6*PI*x[0] + (i+1)*PI/n);
		g += y*y; 
	}
	f[1] = 1 - pow(x[0],0.2) + 2 * g / J2;
	
}

static void  uf8 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0, J3=0;
	double 	g = 0, t = 0;

	// f1
	for (i=3, g=0; i < n; i+=3) {
		J1++;
		t =  x[i] - 2.0*x[1]*sin(2.0*PI*x[0]+(i+1)*PI/n);
		g += t * t; 
	}
	f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2 * g / J1;
	
	// f2
	for (i=4, g=0; i < n; i+=3) {
		J2++;
		t =  x[i] - 2.0*x[1]*sin(2.0*PI*x[0]+(i+1)*PI/n);
		g += t * t; 
	}
	f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) +  2 * g / J2;
	
	// f3
	for (i=2, g=0; i < n; i+=3) {
		J3++;
		t =  x[i] - 2.0*x[1]*sin(2.0*PI*x[0]+(i+1)*PI/n);
		g += t * t; 
	}
	f[2] = sin(0.5*PI*x[0]) + 2 * g / J3;
}

static void  uf9 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0, J3=0;
	double 	g = 0, t = 0;
	double 	epsilon = 0.1;

	// f1
	for (i=3, g=0; i < n; i+=3) {
		J1++;
		t =  x[i] - 2.0*x[1]*sin(2.0*PI*x[0]+(i+1)*PI/n);
		g += t * t; 
	}
	t =  (1.0+epsilon)*(1.0-4.0*(2.0*x[0]-1.0)*(2.0*x[0]-1.0));
	t = t > 0 ? t : 0;
	f[0] = 0.5*(t + 2*x[0])*x[1] + 2 * g / J1;
	
	// f2
	for (i=4, g=0; i < n; i+=3) {
		J2++;
		t =  x[i] - 2.0*x[1]*sin(2.0*PI*x[0]+(i+1)*PI/n);
		g += t * t; 
	}
	t =  (1.0+epsilon)*(1.0-4.0*(2.0*x[0]-1.0)*(2.0*x[0]-1.0));
	t = t > 0 ? t : 0;
	f[1] = 0.5*(t - 2*x[0] + 2)*x[1] +  2 * g / J2;
	
	// f3
	for (i=2, g=0; i < n; i+=3) {
		J3++;
		t =  x[i] - 2.0*x[1]*sin(2.0*PI*x[0]+(i+1)*PI/n);
		g += t * t; 
	}
	f[2] = 1 - x[1] + 2 * g / J3;
}

static void  uf10 (double *x, int n, double *f, int m) {
	int 	i, J1=0, J2=0, J3=0, y;
	double 	g = 0, t = 0;

	// f1
	for (i=3, g=0; i < n; i+=3) {
		J1++;
		y = x[i] - 2.0*x[1]*sin(2.0*PI*x[0]+(i+1)*PI/n);
		t = 4.0*y*y - cos(8.0*PI*y) + 1.0;
		g += t; 
	}
	f[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2 * g / J1;
	
	// f2
	for (i=4, g=0; i < n; i+=3) {
		J2++;
		y = x[i] - 2.0*x[1]*sin(2.0*PI*x[0]+(i+1)*PI/n);
		t = 4.0*y*y - cos(8.0*PI*y) + 1.0;
		g += t; 
	}
	f[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) +  2 * g / J2;
	
	// f3
	for (i=2, g=0; i < n; i+=3) {
		J3++;
		y = x[i] - 2.0*x[1]*sin(2.0*PI*x[0]+(i+1)*PI/n);
		t = 4.0*y*y - cos(8.0*PI*y) + 1.0;
		g += t; 
	}
	f[2] = sin(0.5*PI*x[0]) + 2 * g / J3;
	
}


//**************************************************************************************************************
//**************************************************************************************************************
static Matrix_t* UF1_sample ();
static Matrix_t* UF4_sample ();
static Matrix_t* UF5_sample ();
static Matrix_t* UF6_sample ();
static Matrix_t* UF7_sample ();
static Matrix_t* UF8_sample ();
static Matrix_t* UF9_sample ();

Matrix_t* UF_sample (int No, int numObj) {
	switch (No) {
		case 1:
		case 2:
		case 3:
			return	UF1_sample ();
		case 4:
			return	UF4_sample ();
		case 5:
			return	UF5_sample ();
		case 6:
			return	UF6_sample ();
		case 7:
			return	UF7_sample ();
		case 8:
		case 10:
			return	UF8_sample ();
		case 9:
			return	UF9_sample ();
		default:
                       fprintf (stderr, "UF%d_sample have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

static Matrix_t* UF1_sample () {
	int 	  i;
	double    d = 1.0/NUM_SAMPLE, f1, f2, y1;
	Matrix_t* sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		y1 = i * d;
		if (0 == (i&1) ) {
			f1 = y1;
			f2 = 1.0 - sqrt (f1);
		} else {
			f2 = y1;
			f1 = (1.0-f2)*(1.0-f2);
		}
		sample->elements[2*i] = f1;
		sample->elements[2*i+1] = f2;
	}
	
	return sample;
}
static Matrix_t* UF4_sample () {
	int 	  i;
	double    d = 1.0/NUM_SAMPLE, f1, f2, y1;
	Matrix_t* sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		y1 = i * d;
		f1 = y1;
		f2 = 1.0 - f1 * f1;
		sample->elements[2*i] = f1;
		sample->elements[2*i+1] = f2;
	}
	
	return sample;
}

static Matrix_t* UF5_sample () {
	int 	  i, N = 10;
	double    d = 1.0/(2*N), f1, f2, y1;
	Matrix_t* sample = Matrix_new (2*N+1, 2);

	for (i=0; i<=2*N; i++) {
		y1 = i * d;
		f1 = y1;
		f2 = 1.0 - f1;
		sample->elements[2*i] = f1;
		sample->elements[2*i+1] = f2;
	}
	
	return sample;
}

static Matrix_t* UF6_sample () {
	int 	  i, j;
	double    d = 1.0/NUM_SAMPLE, f1, f2, y1;
	Matrix_t* sample = Matrix_new (NUM_SAMPLE+1, 2);

	sample->elements[0] = 0;
	sample->elements[1] = 1;
	for (i=1, j=1; i<=NUM_SAMPLE; i++) {
		y1 = i * d;
		if ((y1 > 0.25 && y1 < 0.5 ) || y1 > 0.75) {
			f1 = y1;
			f2 = 1.0 - f1;
			sample->elements[2*j] = f1;
			sample->elements[2*j+1] = f2;
			j++;
		}
	}
	sample->rowDim = j;

	return sample;
}

static Matrix_t* UF7_sample () {
	int 	  i;
	double    d = 1.0/NUM_SAMPLE, f1, f2, y1;
	Matrix_t* sample = Matrix_new (NUM_SAMPLE+1, 2);

	for (i=0; i<=NUM_SAMPLE; i++) {
		y1 = i * d;
		f1 = y1;
		f2 = 1.0 - y1;
		sample->elements[2*i] = f1;
		sample->elements[2*i+1] = f2;
	}
	
	return sample;
}

static Matrix_t* UF8_sample () {
	return DTLZ_sample (2, 3);
}

static Matrix_t* UF9_sample () {
	int i, j;
	Matrix_t* sample = NULL;
	Matrix_t* S = NULL;
	
	S = LSMOP_sample (1, 3);
	sample = Matrix_new (S->rowDim, S->colDim);

	for (i=0, j=0; i<S->rowDim; i++) {
		if (4.0*S->elements[i*3 + 0] < 1 - S->elements[i*3 + 2]) {
			sample->elements[j*3] = S->elements[i*3];  
			sample->elements[j*3 + 1] = S->elements[i*3+1];  
			sample->elements[j*3 + 2] = S->elements[i*3+2];  
			j++;
		} else if (4.0*S->elements[i*3 + 0] > 3 * (1 - S->elements[i*3 + 2])) {
			sample->elements[j*3] = S->elements[i*3];  
			sample->elements[j*3 + 1] = S->elements[i*3+1];  
			sample->elements[j*3 + 2] = S->elements[i*3+2];  
			j++;
		}
	}

	Matrix_free (&S);
	sample->rowDim = j;	
	return sample;
}
