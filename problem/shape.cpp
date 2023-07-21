#include "shape.h"
#include "dominate.h"
#include "link.h"
#include "matrix.h"
#include "normalize.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.14159265358979323846264338327950288419716939937510 

/* A linear hyperplane, where SUM (f_0:M) = 1*/
void H1 (double *x, double *g, double *f, int M) {
	int i, j;

	for (i=0; i<M; i++) {
		f[i] = 1.0;
		for (j=0; j<(M-1-i); j++) {
			f[i] *= x[j];
		}
		if (i > 0) {
			f[i] *= (1.0 - x[j]);
		}
	}
}

void H2 (double *x, double *g, double *f, int M) {
	int i, j;

	for (i=0; i<M; i++) {
		f[i] = 1.0;
		for (j=0; j<(M-1-i); j++) {
			f[i] *= cos(0.5*PI*x[j]);
		}
		if (i > 0) {
			f[i] *= sin(0.5*PI*x[j]);
		}
	}
}

void H3 (double *x, double *g, double *f, int M) {
	int 	i;
	double 	sum = 0;

	for (i=0; i<M-1; i++) {
		f[i] = x[i] / (1 + g[i]);
		sum += (x[i] * (1 + sin (3*PI*x[i])));
	}
	sum = sum / (2.0 + g[M-1]);	

	f[M-1] = (M - sum) * ( 2.0 + g[M-1]) / (1.0 + g[M-1]);
}

/*
 the simplex lattice design approach
 */
Matrix_t *H1_sample (int M, int H) {
	Matrix_t *sample = Matrix_new ();
	double 	tmp;
	int 	i, j;
	size_t 	size;
	int 	solution[M], cur, sum;
	
	sample->colDim = M;
        for (i=0, tmp=1.0; i<M-1; i++) {
		tmp *= (M - 1.0 + H - i) / (M - 1.0 - i);
	}

	sample->rowDim = (int)(tmp + 0.5);

	size = sample->rowDim * M * sizeof (double);
	sample->elements = (double *)malloc (size);

	i = 0;
	cur=1;
	solution[0] = 0;
	sum = 0;
	while (cur > 0) {
		if (cur < M && sum <= H) {
			solution[cur] = 0;
			cur++;
		} else if (cur == M && sum < H) {
			solution[cur-1]++;
			sum++;
		} else if (cur == M && sum == H) {
			for (j=0; j<M; j++) {
				sample->elements[i++] = solution[j] * 1.0 / H;
			}
			
			// 
			cur--;
			sum -= solution[cur];
			solution[cur-1]++;	
			sum++;
		} else 	if (sum > H) {
			cur--;
			sum -= solution[cur];
			if (cur > 0) {
				solution[cur-1]++;
				sum++;
			}
		}
	}

	return sample;
}

Matrix_t *H2_sample (int M, int H) {
	Matrix_t *sample = NULL;
	Matrix_t *normS = NULL;

	sample = H1_sample (M, H);
	normS = normalize_sphere (sample);
	Matrix_free (&sample);

	return normS;
}

Matrix_t *H3_sample (int M, int H) {
	Matrix_t *grid = Grid_sample (M-1, H);
	Matrix_t *sample = Matrix_new (grid->rowDim, M);
	double tmp, t;
	int i, j;

	for (i=0; i<sample->rowDim; i++) {
		for (j=0, tmp=0; j<M-1; j++) {
			t = grid->elements [i*(M-1)+j];
			sample->elements[i*M+j] = t;
			tmp += t * (1 + sin (3*PI*t));
		}
		sample->elements[i*M+j] = 2.0 * M - tmp;
	}
	
	// remove the dominated points
	Matrix_t *subM = Matrix_front (sample);

	Matrix_free (&grid);
	Matrix_free (&sample);

	return subM;
}

Matrix_t *H4_sample (int M, int H) {
	Matrix_t *sample = Matrix_new ();
	int i, j, k;
	size_t size;
	
	sample->colDim = M;
	sample->rowDim = H+1;

	size = sample->rowDim * M * sizeof (double);
	sample->elements = (double *)malloc (size);

	for (i=0; i<=H; i++) {
		for (j=0; j<M; j++) {
			sample->elements[i*M+j] = 1.0;
			for (k=0; k<(M-1-j); k++) {
				if (k==0)
					sample->elements[i*M+j] *= cos(0.5*PI*i/H);
				else
					sample->elements[i*M+j] *= cos(0.5*PI*0.5);
			}
			if (j > 0) {
				if (k==0)
					sample->elements[i*M+j] *= sin(0.5*PI*i/H);
				else
					sample->elements[i*M+j] *= sin(0.5*PI*0.5);
			}
		}
	}
	return sample;
}


void linear (double *x, double *f, int M) {
	H1 (x, NULL, f, M);
}

void convex (double *x, double *f, int M) {
	int i, j;

	for (i=0; i<M; i++) {
		f[i] = 1.0;
		for (j=0; j<(M-1-i); j++) {
			f[i] *= (1.0 - cos (x[j]*PI/2));
		}
		if (i > 0) {
			f[i] *= (1.0 - sin (x[j]*PI/2));
		}
	}
}

void concave (double *x, double *f, int M) {
	int i, j;

	for (i=0; i<M; i++) {
		f[i] = 1.0;
		for (j=0; j<(M-1-i); j++) {
			f[i] *= sin(0.5*PI*x[j]);
		}
		if (i > 0) {
			f[i] *= cos(0.5*PI*x[j]);
		}
	}
}

void mixedM (double *x, double *f, int M) {
	f[M-1] = 1.0 - x[0] - cos (10*PI*x[0] + 0.5*PI) /(10*PI);
}

void discM (double *x, double *f, int M) {
	f[M-1] = 1.0 - x[0]*cos(5*PI*x[0])*cos(5*PI*x[0]);
}

Matrix_t *Grid_sample (int Dim, int H) {
	Matrix_t *sample = Matrix_new ();
	double tmp, t;
	int i, j;
	size_t size;
	int solution[Dim+1];
	
	sample->colDim = Dim;
	for (i=0, tmp=1.0; i<Dim; tmp *=(H+1), i++);
	sample->rowDim = (int)(tmp + 0.5);

	size = sample->rowDim * Dim * sizeof (double);
	sample->elements = (double *)malloc (size);
	
	for (i=0; i<=Dim; solution[i] = 0, i++);
	i=0;
	while (solution[Dim] <= 0) {
		for (j=0; j<Dim; j++) {
			t = (1.0 * solution[j]) / H;
			sample->elements[i++] = t;
		}

		solution[0]++;
		for (j=0; j<Dim; j++) {
			if (solution[j] > H) {
				solution[j] -= (H+1);
				solution[j+1]++;
			} 
		}
	}

	return sample;
}


int conb (int n, int m) {        // conbination: to select m individuals from n $C_n^m$
	double tmp=1.0;
	int i;

	for (i=0; i<m; i++) {
		tmp *= 1.0*(n-i)/(m-i);
	}
	return (int)(tmp+0.5);
}

int grid (int dim, int H) {
	double tmp=1.0;
	int i;

	for (i=0; i<dim; i++) {
		tmp *= (H+1);
	}
	return (int)(tmp+0.5);
}
