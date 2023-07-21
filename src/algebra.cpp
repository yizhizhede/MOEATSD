#include "algebra.h"
#include "matrix.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>

double mean (double *A, int n) {
	double 	tmp;
	int 	i;	

	for (i=0, tmp=0; i<n; tmp += A[i], i++);
	return tmp / n;
}

double VAR  (double *A, int n) {
	double 	s, v;
	int 	i;	

	for (i=0, s=0; i<n; i++) { s += A[i]; }
	for (i=0, v=0, s /= n; i<n; i++) { v += (A[i]-s)*(A[i]-s); }
	return v / n;
}

double STD (double *A, int n) {
	return sqrt (VAR (A, n));
}

double norm (double *A, int n) {
	double tmp;
	int i;	

	for (i=0, tmp=0; i<n; tmp += A[i]*A[i], i++);
	return sqrt (tmp);
}

double sum  (double *A, int n) {
	double tmp;
	int i;	

	for (i=0, tmp=0; i<n; tmp += A[i], i++);
	return tmp;
}

int sum  (int *A, int n) {
	int tmp;
	int i;	

	for (i=0, tmp=0; i<n; tmp += A[i], i++);
	return tmp;
}

double inner_product (double *A, double *B, int n) {
	double tmp;
	int i;	

	for (i=0, tmp=0; i<n; tmp += A[i]*B[i], i++);
	return tmp;
}

double* scale_product (double s, double *A, int n) {
	size_t size = n * sizeof (double);
	double *B= (double *)malloc (size);
	int i;	

	for (i=0; i<n; B[i] = A[i]*s, i++);
	return B;
}

double *vector_substract (double *A, double *B, int n) {
	size_t size = n * sizeof (double);
	double * C= (double *)malloc (size);
	int i;	

	for (i=0; i<n; C[i] = A[i] - B[i], i++);
	return C;
}

double distance_p2l (double *p, double *line, int n) {
	double *a, *b;
	double t1, t2;

	t1 = inner_product (p, line, n);
	t2 = norm (line, n);
	a = scale_product (t1/(t2*t2), line, n);
	b = vector_substract (p, a, n);
	t1 = norm (b, n);
	free (a);
	free (b);
	return t1;
}

double distance_p2p (double *p1, double *p2, int n) {
	double tmp;
	int i;	

	for (i=0, tmp=0; i<n; tmp += (p1[i] - p2[i])*(p1[i] - p2[i]), i++);
	return sqrt (tmp);
}

double vector_angle (double *v1, double *v2, int n) {
	double t1,t2;

        t1 = norm(v1, n) * norm(v2, n);
        if (!t1) return 0.0;
         
        t2 = inner_product(v1, v2, n);
        t1 = t2 / t1;
        
        if (t1>1.0)
                t1 = 1.0;
        else if(t1<-1.0)
                t1 = -1.0;
        return 90.0 * acos(t1) / acos(0.0);
}
 
Matrix_t* scale_product (double s, Matrix_t *M) {
	int i;
	Matrix_t *dupM = Matrix_dup (M);

	for (i=M->rowDim*M->colDim - 1; i>=0; dupM->elements[i] *= s, i--);
	return dupM;
}

double distance_p2m (double *point, Matrix_t *M) {
	double tmp = distance_p2p (point, M->elements, M->colDim), t;
	int i;

	if (NULL == point) return 1.0e+300;
	if (NULL == M || NULL == M->elements || 1 > M->rowDim) return 1.0e+300;

	for (i=1; i<M->rowDim; i++) {
		t = distance_p2p (point, M->elements+i*M->colDim, M->colDim);
		if (t < tmp)
			tmp = t;
	}
	return tmp;
}

double distance_m2m (Matrix_t *M1, Matrix_t *M2) {
	double sum;
	int i;

	if (NULL == M1 || NULL == M1->elements || 1 > M1->rowDim) return 1.0e+300;
	if (NULL == M2 || NULL == M2->elements || 1 > M2->rowDim) return 1.0e+300;

	for (i=0, sum=0; i<M1->rowDim; sum += distance_p2m (M1->elements+i*M1->colDim, M2), i++);
	return sum / M1->rowDim;
}


double Lp_norm (double *A, int n, double p) {
	int i;
	double sum=0;

	if (p < 0) {
		return 0;	
	} else if (p==0) {
		for (i=0; i<n; i++)
			sum += (A[i]!=0)?1:0;
	} else if (p==1) {
		for (i=0; i<n; i++)
			sum += A[i] > 0 ? A[i] : -A[i];
	} else if (p==2) {
		sum = norm (A, n);
	} else {
		for (i=0; i<n; i++) {
			sum += A[i]>0?pow (A[i], p):pow (-A[i], p);
		}
		sum = pow (sum, 1.0/p);
	}
	
	return sum;
}

int isLinear (double *X, double *Y, int n) {
	int 	i;
	double	Sxx=0, Sxy=0, Syy=0;
	double	x_bar = mean (X, n);
	double 	y_bar = mean (Y, n);
	double 	r2, d;
		
	for (i=0; i<n; i++) {
		Sxx += (X[i] - x_bar)*(X[i] - x_bar);
		Sxy += (X[i] - x_bar)*(Y[i] - y_bar);
		Syy += (Y[i] - y_bar)*(Y[i] - y_bar);
	}

	if (Sxx < DBL_EPSILON) {
	 	if (Syy < DBL_EPSILON) 
			return 1;
		else 
			return 0;
	}
	if (Syy < DBL_EPSILON) {
		return 1;
	}

	r2 = (Sxy*Sxy)/(Sxx*Syy);
	d = fabs(r2 - 1.0);
	if (d < FLT_EPSILON) {
		return 1;
	} else {
		return 0;
	}
}

// Y = X * belta
void make_line (Matrix_t *X,  double *Y, double* belta, double& delta) {
	int 		i, j, k;
	int 		rowDim = X->rowDim;
	int 		colDim = X->colDim;
	Matrix_t* 	T = Matrix_new (colDim, rowDim);
	Matrix_t* 	TX= Matrix_new (colDim, colDim);
	Matrix_t* 	V =  NULL; //  (colDim, colDim);
	Matrix_t* 	VT= Matrix_new (colDim, rowDim);
	double		epsilon[rowDim], RSS, t;

	// T
	for (i=0; i<colDim; i++) {
		for (j=0; j<rowDim; j++) {
			T->elements[i*rowDim+j] = X->elements[j*colDim+i];
		}
	}

	// TX = T * X
	for (i=0; i<colDim; i++) {
		for (j=0; j<colDim; j++) {
			for (k=0, t=0; k<rowDim; k++) {
				t += T->elements[i*rowDim+k] * X->elements[k*colDim+j];
			}
			TX->elements[i*colDim+j] = t;
		}
	}
	
	// inverse: V
	V = Matrix_inverse (TX);

	// VT = V * T
	for (i=0; i<colDim; i++) {
		for (j=0; j<rowDim; j++) {
			for (k=0, t=0; k<colDim; k++) {
				t += V->elements[i*colDim+k] * T->elements[k*rowDim+j];
			}
			VT->elements[i*rowDim+j] = t;
		}
	}

	// belta = VT * Y
	for (i=0; i<colDim; i++) {
		for (k=0, t=0; k<rowDim; k++) {
			t += VT->elements[i*rowDim+k] * Y[k];
		}
		belta[i] = t;
	}

	// epsilon = Y - X*Belta
	for (i=0; i<rowDim; i++) {
		for (j=0, t=0; j<colDim; j++) {
			t += X->elements[i*colDim+j] * belta[j];
		}
		epsilon[i] = Y[i] - t;
	}

	// RSS, delta
	for (i=0, RSS=0; i<rowDim; i++) {
		RSS += (epsilon[i])*(epsilon[i]);
	}
	if (rowDim > 2) {
		RSS /= (rowDim - 2);
	}
	delta = sqrt (RSS);

	// free
	Matrix_free (&T); Matrix_free (&TX); Matrix_free (&V); Matrix_free (&VT);
}

void mean (double *X, int rowDim, int colDim, double* X_bar) {
	int 	i, j;
	
	memset (X_bar, 0, colDim*sizeof (double));
	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			X_bar[j] += X[i*colDim+j];
		}
	}
	for (j=0; j<colDim; j++) {
		X_bar[j] = X_bar[j] /rowDim;
	}
}
