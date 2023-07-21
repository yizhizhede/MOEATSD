#ifndef _SINGULAR_VALUE_DECOMPOSITION_H
#define _SINGULAR_VALUE_DECOMPOSITION_H

#include <stdio.h>

/***********************************************
 * transport A, B is the result.	       *
 **********************************************/
void 
Tra(double *A,int nrows,int ncols,double *B);


/*****************************************************
 * decomposite A to U D V,namely A = UDV',           *
 * A is n*m, 					     *	
 * if n >= m then U is n * m, V is m * m, D is m * 1 *
 * if n <  m then U is n * n, V is m * n, D is n * 1 *
 ****************************************************/
int 
svd(double *A,int nrows,int ncols,double *U,double *D,double *V);

/**********************************************
 * Inverse A, B is the result.		      *	
 **********************************************/
void
Inv(double *A,int nrows,double *B);

/**********************************************
 * compute the coordinate x of B ,based on A  *
 **********************************************/
void 
Coo(double *A,int nrows,double *B,double *x);


#endif
