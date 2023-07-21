#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "nntp.h"
#include "matrix.h"
#include "shape.h"

// test problem
static Problem_t *P;

// 
static Matrix_t 	*Data;
static Matrix_t 	*TrainData;
static int		nHidden = 20;

// 
static void nntp (double *x, int n, double *f, int m);

//
Problem_t* NNTP_new (char *title, int numObj, int numVar) {
	// common variable
	int 	i, j;
	size_t 	size;
	double*	lowBound = NULL; 
	double* uppBound = NULL;

	//
	double*	aveD = NULL;
	double* stdD = NULL; 	
	int	rowDim;
	int	colDim;

	// allocating memory for a problem
        Problem_t *problem = (Problem_t *)malloc (sizeof (Problem_t));
        if (problem == NULL) {
        	fprintf (stderr, "Allocating memory failed\n");         
               	exit (-1);
        }
	strcpy (problem->title, title);

	// load data
	if (!strcmp (title, "NNTP1")) {Data = Matrix_read ((char *)"./problem/dataset/DatasetStatlogAustralian");}
	else if (!strcmp (title, "NNTP2")) {Data = Matrix_read ((char *)"./problem/dataset/DatasetClimate");} 
	else if (!strcmp (title, "NNTP3")) {Data = Matrix_read ((char *)"./problem/dataset/DatasetStatlogGerman");} 
	else if (!strcmp (title, "NNTP4")) {Data = Matrix_read ((char *)"./problem/dataset/DatasetConnectionistBenchSonar");} 

	// get TrainData	
	aveD 	= Matrix_mean (Data);
	stdD 	= Matrix_std (Data);
	rowDim 	= (int)ceil(Data->rowDim * 0.8);
	colDim 	= Data->colDim;
	TrainData = Matrix_new (rowDim, colDim);	
	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim-1; j++) {
			TrainData->elements[i*colDim+j] = (Data->elements[i*colDim+j] - aveD[j]) / stdD[j];
		}
		TrainData->elements[i*colDim+j] = Data->elements[i*colDim+j];
	}
	free (aveD);	
	free (stdD);

	// numObj, numVar
	problem->numObj = 2; 
	problem->numVar = (colDim * nHidden) + (nHidden + 1) * 1;

	//
	numObj = 2;
	numVar = problem->numVar;

	// setting the bound of varibles
	size = numVar * sizeof (double);
	lowBound = (double *)malloc (size);
	uppBound = (double *)malloc (size);
	for (i=0; i<numVar; i++) {
		lowBound[i] = -1.0;
		uppBound[i] = 1.0;
	}
	problem->lowBound = lowBound;
	problem->uppBound = uppBound;

	// 
	if (!strcmp (title, "NNTP1")) {
		P = problem; problem->evaluate = nntp;
	} else if (!strcmp (title, "NNTP2")) {
		P = problem; problem->evaluate = nntp;
	} else if (!strcmp (title, "NNTP3")) {
		P = problem; problem->evaluate = nntp;
	} else if (!strcmp (title, "NNTP4")) {
		P = problem; problem->evaluate = nntp;
	} else { 
		fprintf (stderr, "error: %s is undefined\n", title);
		exit (0);
	}

	return problem;
}



static void Predict  (Matrix_t* W1, Matrix_t* W2, Matrix_t* X, Matrix_t* Y, Matrix_t* Z);
static void Train    (Matrix_t* W1, Matrix_t* W2, Matrix_t* X, Matrix_t* Y, Matrix_t* Z, int nEpoch);
static void toW1W2   (Matrix_t* W1, Matrix_t* W2, double *x);
static void fromW1W2 (Matrix_t* W1, Matrix_t* W2, double *x);

static void nntp (double *x, int n, double *f, int m) {
	int 		i, j, k;
	double 		g = 0;
	Matrix_t*	W1 = Matrix_new (TrainData->colDim, nHidden);	
	Matrix_t*	W2 = Matrix_new (nHidden+1, 1);	
	Matrix_t*	Y  = Matrix_new (TrainData->rowDim, nHidden);
	Matrix_t*	Z  = Matrix_new (TrainData->rowDim, 1);

	// copy x to W1 and W2
	toW1W2 (W1, W2, x);

	// Train 
	Train (W1, W2, TrainData, Y, Z, 1);

	// copy W1 and W2 to x
	fromW1W2 (W1, W2, x);

	// Predict
	Predict  (W1, W2, TrainData, Y, Z);
	
	// f1
	for (i=0, g = 0; i<n; i++) if (fabs(x[i]) > DBL_EPSILON) {
		g += 1;	
	}
	f[0] = g / n;

	//f2
	for (i=0, g=0; i<Z->rowDim; i++) { 
		j = (int)(Z->elements[i] + 0.5); 
		k = (int)(TrainData->elements[i*TrainData->colDim+(TrainData->colDim-1)]);

		if (j != k) {
			g += 1;
		}
	}
	f[1] = g / Z->rowDim;

	//
	Matrix_free (&Z);
	Matrix_free (&Y);
	Matrix_free (&W2);
	Matrix_free (&W1);
}

static void toW1W2   (Matrix_t* W1, Matrix_t* W2, double *x) {
	int 	i, j, k;

	// copy x to W1 and W2
	for (j=0, k=0; j<W1->colDim; j++) {
		for (i=0; i<W1->rowDim; i++) {
			W1->elements[i*W1->colDim+j] = x[k];
			k++;
		}
	}
	for (i=0; i<W2->rowDim; i++) {
		W2->elements[i] = x[k];
		k++;
	}
}
static void fromW1W2 (Matrix_t* W1, Matrix_t* W2, double *x) {
	int 	i, j, k;

	// copy W1 and W2 to x
	for (j=0, k=0; j<W1->colDim; j++) {
		for (i=0; i<W1->rowDim; i++) {
			x[k] = W1->elements[i*W1->colDim+j];
			k++;
		}
	}
	for (i=0; i<W2->rowDim; i++) {
		x[k] = W2->elements[i];
		k++;
	}
}

static void Predict (Matrix_t* W1, Matrix_t* W2, Matrix_t* X, Matrix_t* Y, Matrix_t* Z) {
	int 	i, j, k;
	double	t;

	// Y
	for (i=0; i<Y->rowDim; i++) {
		for (j=0; j<Y->colDim; j++) {
			t = W1->elements[0*W1->colDim+j];
			for (k=1; k<W1->rowDim; k++) {
				t += X->elements[i*X->colDim+(k-1)] * W1->elements[k*W1->colDim+j];
			}
			//
			Y->elements[i*Y->colDim+j] = 1 - 2 / (1 + exp(t));
		}
	}

	// Z
	for (i=0; i<Z->rowDim; i++) {
		t = W2->elements[0];
		for (k=1; k<W2->rowDim; k++) {
			t += Y->elements[i*Y->colDim+(k-1)] * W2->elements[k];
		}
		//
		Z->elements[i] = 1 / (1 + exp(-t));
	}
}

static void Train (Matrix_t* W1, Matrix_t* W2, Matrix_t* X, Matrix_t* Y, Matrix_t* Z, int nEpoch) {
	int 		i, j;
	int		epoch, row;
	Matrix_t*	iP = Matrix_new (X->rowDim, 1);
	Matrix_t*	iQ = Matrix_new (X->rowDim, nHidden);
	Matrix_t*	D1 = Matrix_new (W1->rowDim, W1->colDim);
	Matrix_t*	D2 = Matrix_new (W2->rowDim, W2->colDim);

	for (epoch=0; epoch<nEpoch; epoch++) {
		// Predict
		Predict (W1, W2, X, Y, Z);

		// iP
		for (i=0; i<iP->rowDim; i++) {
			iP->elements[i] = (Z->elements[i] - X->elements[i*X->colDim+(X->colDim-1)]) * 
			(Z->elements[i]) * (1 - Z->elements[i]);
		}

		// iQ
		for (i=0; i<iQ->rowDim; i++) {
			for (j=0; j<iQ->colDim; j++) {
				iQ->elements[i*iQ->colDim+j] = iP->elements[i] * W2->elements[j+1] * 
					(1 - Y->elements[i*Y->colDim+j]*Y->elements[i*Y->colDim+j]);
			}
		}

		// D1, D2 
		for (row=0; row < (X->rowDim); row++) {
			// D1
			for (j=0; j<D1->colDim; j++) {
				D1->elements[0*D1->colDim+j] += iQ->elements[row*iQ->colDim+j];
			}
			for (i=1; i<D1->rowDim; i++) {
				for (j=0; j<D1->colDim; j++) {
			D1->elements[i*D1->colDim+j] += X->elements[row*X->colDim+(i-1)] * iQ->elements[row*iQ->colDim+j];
				}
			}
			// D2
			D2->elements[0] += iP->elements[row];
			for (i=1; i<D2->rowDim; i++) {
				D2->elements[i] += Y->elements[row*Y->colDim+(i-1)] * iP->elements[row];
			}
		}

		// W1
		for (i=0; i<W1->rowDim; i++) {
			for (j=0; j<W1->colDim; j++) {
				W1->elements[i*W1->colDim+j] -= (D1->elements[i*D1->colDim+j] / X->rowDim);
				if (W1->elements[i*W1->colDim+j] > 1) {
					W1->elements[i*W1->colDim+j] = 1;
				}
				
				if (W1->elements[i*W1->colDim+j] < -1) {
					W1->elements[i*W1->colDim+j] = -1;
				}

				//
				if (!(W1->elements[i*W1->colDim+j] < 1.1  && W1->elements[i*W1->colDim+j] > -1.1)) {
					printf ("%d, %d\n", W1->rowDim, W1->colDim);
					printf ("W1[%d, %d] = %f\n", i, j, W1->elements[i*W1->colDim+j]);
					printf ("D1[%d, %d] = %f\n", i, j, D1->elements[i*D1->colDim+j]);
					exit (0);
				}
			}
		}

		// W2
		for (i=0; i<W2->rowDim; i++) {
			W2->elements[i] -= D2->elements[i] / X->rowDim;
			if (W2->elements[i] > 1) {
				W2->elements[i] = 1;
			}
			
			if (W2->elements[i] < -1) {
				W2->elements[i] = -1;
			}

			//
			if (!(W2->elements[i] < 1.1  && W2->elements[i] > -1.1)) {
				printf ("W2[%d] = %f\n", i, W2->elements[i]);
				exit (0);
			}
		}


	}
	
	//
	Matrix_free (&D2);
	Matrix_free (&D1);
	Matrix_free (&iP);
	Matrix_free (&iQ);
}

//**************************************************************************************************************
//**************************************************************************************************************
Matrix_t *NNTP_sample (int No, int numObj) {
	Matrix_t* sample = Matrix_new (1, numObj);
	return 	sample;
}
