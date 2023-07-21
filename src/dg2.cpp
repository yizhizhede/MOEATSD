#include "dg2.h"
#include "algebra.h"
#include "problem.h"
#include "myrandom.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX(a,b) ((a)>(b)?(a):(b))

static double array_max(double *array,int len);
static double gamma_func(double d);

Matrix_t *dg2 () {
	int 		i, j;
	int 		numVar = (Problem_get ()) -> numVar;
	int 		numObj = (Problem_get ()) -> numObj;
	double	 	f_base;
	Matrix_t*	f_hat = NULL;
	Matrix_t*	F = NULL;
	Matrix_t*	Lambda = NULL;
	Matrix_t*	Theta = NULL;
	double		x1[numVar+10];
	double		x2[numVar+10];
	double		m[numVar+10];
	double*		lowBound = Problem_getLowerBound ();
	double*		uppBound = Problem_getUpperBound ();
	double		buff[100], delta1, delta2;
	double		eInf, eSup, eta0=0, eta1=0, eps;

	f_hat = Matrix_new (numVar, 1);
	F = Matrix_new (numVar, numVar);
	Lambda = Matrix_new (numVar, numVar); 
	Theta = Matrix_new (numVar, numVar); 


	for (i=0; i<numVar; i++) {
		x1[i] = lowBound[i] +randu()*(uppBound[i] - lowBound[i]);
		m[i]  = lowBound[i] +randu()*(uppBound[i] - lowBound[i]);
	}

	// compute f_base
        Problem_evaluate (x1, numVar, buff, numObj);
	f_base = sum(buff, numObj);

	// compute f_hat
	for (i=0; i<numVar; i++) {
		memcpy (x2, x1, numVar*sizeof (double));
		x2[i] = m[i];
        	Problem_evaluate (x2, numVar, buff, numObj);
		f_hat->elements[i] = sum(buff, numObj);
	}

	// compuate F
	for (i=0; i<numVar-1; i++) {
		for (j=i+1; j<numVar; j++) {
			memcpy (x2, x1, numVar*sizeof (double));
			x2[i] = m[i];
			x2[j] = m[j];
			Problem_evaluate (x2, numVar, buff, numObj);
			F->elements[i*numVar+j] = sum(buff, numObj);
			F->elements[j*numVar+i] = F->elements[i*numVar+j]; 
		}
	}

	// compute Lambda: function ISM ()
	for (i=0; i<numVar-1; i++) {
		for (j=i+1; j<numVar; j++) {
			delta1 = f_hat->elements[i] - f_base;
			delta2 = F->elements[i*numVar+j] - f_hat->elements[j];
			Lambda->elements[i*numVar+j] = fabs (delta1 - delta2);
			Lambda->elements[j*numVar+i] =Lambda->elements[i*numVar+j];
		}
	}

	// compute Theta: function DSM () 
	for (i=0; i<numVar*numVar; i++) {
		Theta->elements[i] = 100;
	}
	for (i=0; i<numVar-1; i++) {
		for (j=i+1; j<numVar; j++) {
			buff[0]	= f_base;
			buff[1] = F->elements[i*numVar+j];
			buff[2] = f_hat->elements[i];
			buff[3] = f_hat->elements[j];

			eInf = gamma_func(2.0)*MAX(buff[0]+buff[1],buff[2]+buff[3]);
			eSup = gamma_func(sqrt(numVar))*array_max (buff, 4);	

			if (Lambda->elements[i*numVar+j] < eInf) {
				Theta->elements[i*numVar+j] = 0;
				Theta->elements[j*numVar+i] = 0;
				eta0 += 1;
			} else if (Lambda->elements[i*numVar+j] > eSup) {
				Theta->elements[i*numVar+j] = 1.0;
				Theta->elements[j*numVar+i] = 1.0;
				eta1 += 1;
			}
		}
	}
	for (i=0; i<numVar-1; i++) {
		for (j=i+1; j<numVar; j++) if (Theta->elements[i*numVar+j] > 2) {
			buff[0]	= f_base;
			buff[1] = F->elements[i*numVar+j];
			buff[2] = f_hat->elements[i];
			buff[3] = f_hat->elements[j];

			eInf = gamma_func(2.0)*MAX(buff[0]+buff[1],buff[2]+buff[3]);
			eSup = gamma_func(sqrt(numVar))*array_max (buff, 4);	
			eps  = (eta0*eInf + eta1*eSup) / (eta0 + eta1);

			if (Lambda->elements[i*numVar+j] > eps) {
				Theta->elements[i*numVar+j] = 1.0;
				Theta->elements[j*numVar+i] = 1.0;
			} else {
				Theta->elements[i*numVar+j] = 0;
				Theta->elements[j*numVar+i] = 0;
			}
		}
	}

	//
	for (i=0; i<numVar; i++) {
		Theta->elements[i*numVar+i] = 0;
	}


	Matrix_free (&f_hat);
	Matrix_free (&F);
	Matrix_free (&Lambda); 
	return Theta;
}

static double array_max(double *array,int len) {
        double 	t=array[len-1];
        int 	i;

        for (i=len-2; i>=0; i--) {
                if (array[i] > t) {
                        t = array[i];
		}
        }
        return t;
}

static double gamma_func(double d) {
     double muM = FLT_EPSILON /2.0;
     return (d * muM)/(1 - (d * muM));
}
     
