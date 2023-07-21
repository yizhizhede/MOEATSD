#ifndef _MATRIX_H
#define _MATRIX_H

#include <stdio.h>
#include "link.h"

/* data type Matrix_t */
typedef struct Mat_tag {
	int rowDim;
       	int colDim;
      	double *elements;
} Matrix_t;

typedef struct Vector_tag {
	int dim;
	double *elements;
} Vector_t;

/* 
 operation on Matrix_t, its impletament file is 'src/matrix.cpp'
*/
Matrix_t* Matrix_read (char *filename);
Matrix_t* Matrix_new ();
Matrix_t* Matrix_new (int nrow);
Matrix_t* Matrix_new (int nrow, int ncol);
Matrix_t* Matrix_sub (Matrix_t *M, int *subset, int n);
Matrix_t* Matrix_sub (Matrix_t *M, Link_t *link);
Matrix_t* Matrix_min (Matrix_t* M);
Matrix_t* Matrix_max (Matrix_t* M);
Matrix_t* Matrix_bound (Matrix_t *M);
Matrix_t* Matrix_norm (Matrix_t *M);
Matrix_t* Matrix_norm (Matrix_t *M, Matrix_t *bound);
Matrix_t* Matrix_norm (Matrix_t *M, Matrix_t *maxmum, Matrix_t *minmum);
Matrix_t* Matrix_norm (Matrix_t *M, double *maxmum, double *minmum);
Matrix_t* Matrix_dup (Matrix_t *M);
Matrix_t* Matrix_front (Matrix_t *M);
Matrix_t* Matrix_front_1k (Matrix_t *M);
Matrix_t* Matrix_trim (Matrix_t *M);
Matrix_t* Matrix_extreme (Matrix_t *M);
Matrix_t* Matrix_inverse (Matrix_t *M);
Matrix_t* Matrix_trans (Matrix_t *M);
Matrix_t* Matrix_compress (Matrix_t *M);	// remove the repeated items;
Matrix_t* Matrix_limited (Matrix_t *M, double *p);	// to use a point to constraint M, return a PF. 

double*   Matrix_mean (Matrix_t *M);		// average
double*   Matrix_var (Matrix_t *M);		// varianece 
double*   Matrix_std (Matrix_t *M);		// standard variance 
double	  Matrix_angle (double *p, Matrix_t* M);// angle between the p and the M

void Matrix_print (Matrix_t *M);
void Matrix_print (Matrix_t *M, char *fn);
void Matrix_print (Matrix_t *M, FILE *fp);
void Matrix_latex (Matrix_t *M);
void Matrix_octave (Matrix_t *M);
void Matrix_free (Matrix_t **M);
void Matrix_cat (Matrix_t **M1, Matrix_t *M2);

Vector_t* Vector_new ();
Vector_t* Matrix_sub (Matrix_t *M, int nrow); // nrow = 0, 1, 2, ...
void      Vector_free (Vector_t **v);

/*
  the functions below are impleted in 'src/mysort.cpp'
*/
int* sort (Matrix_t *M);
int* sort (Matrix_t *M, int byCol); 	// byCol = 0, 1, 2, ...
int* sort (Matrix_t *M, const char* order);
int* sort (Matrix_t *M, int byCol, const char* order);

// rank
int* Matrix_rank_by_hypervolume (Matrix_t * M);
int* Matrix_rank_by_crowdingdistance (Matrix_t * M);
int* Matrix_rank_by_sum (Matrix_t * M);

#endif

