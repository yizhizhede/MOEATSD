#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "bk.h"
#include "matrix.h"

#define NUM_SAMPLE 1.0e+4

// test problem
static Problem_t *P;

// 
static void  bk1 (double *x, int n, double *f, int m);

Problem_t *BK_new (char *title, int numObj, int numVar) {
	// common variable
	int 	i;
	size_t 	size;
	double*	lowBound = NULL; 
	double*	uppBound = NULL;

	// allocating memory for a problem
        Problem_t *problem = (Problem_t *)malloc (sizeof (Problem_t));
        if (problem == NULL) {
        	fprintf (stderr, "Allocating memory failed\n");         
               	exit (-1);
        }
	strcpy (problem->title, title);
	problem->numObj = 2; 
	problem->numVar = 2;
	
	// setting the bound of varibles
	size = problem->numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -5.0;
		uppBound[i] = 10.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "BK1")) {
		P = problem; problem->evaluate = bk1;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}


	return problem;
}

static void  bk1 (double *x, int n, double *f, int m) {
	f[0] = x[0]*x[0] + x[1]*x[1];
	f[1] = (x[0]-5)*(x[0]-5) + (x[1]-5)*(x[1]-5);
}


//**************************************************************************************************************
//**************************************************************************************************************
Matrix_t  *BK1_sample ();

Matrix_t  *BK_sample (int No, int numObj) {
	switch (No) {
		case 1:
			return	BK1_sample ();
		default:
                       fprintf (stderr, "BK%d have been not implemented now\n", No);
                       exit (0);
	}
	return NULL;
}

Matrix_t  *BK1_sample () {
	int 		i, j, n; 
	double		d;
	double		x[10];
	Matrix_t*	sample = NULL;

	n = (int)sqrt (NUM_SAMPLE);
	if (n*n < NUM_SAMPLE) {
		n++;
	}
	d=5.0/n;
	sample = Matrix_new ((n+1)*(n+1), 2);

	for (i=0; i<=n; i++) {
		for (j=0; j<=n; j++) {
			x[0] = i*d;
			x[1] = j*d;
			bk1(x, 2, sample->elements+(i*(n+1)+j)*2, 2);
		}
	}
	
	return sample;
}
