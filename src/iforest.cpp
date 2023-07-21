/*
   Title of paper: Isolation-Based Anomaly Detection.
   DOI: 10.1145/2133360.2133363 
   Net Address: http://doi.acm.org/10.1145/2133360.2133363
*/
#include "iforest.h"
#include "myrandom.h"
#include <math.h>
#include <float.h>

static int _t   = 100;		// Number of Trees 
static int _psi = 256;		// Subsampling Size 

static double c (int psi);	//


typedef struct iNode_tag {
	iNode_tag* 	lchild;
	iNode_tag* 	rchild;
	int	   	splitAtt;
	double		splitValue;		
	int		size;
} iNode_t;

static Matrix_t* sample  (Matrix_t *X, int psi);
static Matrix_t* filterL (Matrix_t* X, int q, double p);
static Matrix_t* filterR (Matrix_t* X, int q, double p);
static iNode_t** iForest (Matrix_t *X, int t, int psi);
static iNode_t*  iTree   (Matrix_t *X);
static void  	 freeTree   (iNode_t *T);
static void  	 freeForest (iNode_t **forest, int t);
static double	 pathLength (double* x, iNode_t* T, int hlim, int e);


/****************************************************************************************************************/
double* iForest(Matrix_t* X) {
	int 		t      = _t;
	int 		psi    = _psi;
	iNode_t**	forest = NULL;
	int		i, j;
	double 		h;
	double*		x = NULL;
	int		rowDim = X->rowDim;
	int		colDim = X->colDim;
	double*		s = (double*)malloc (rowDim*sizeof(double));

	//
	forest = iForest (X, t, psi);

	//
	for (i=0; i<rowDim; i++) {
		x = X->elements + i*colDim;
		
		for (j=0, h=0; j<t; j++) {
			h += pathLength (x, forest[j], psi-1, 0);		
		}
		h /= t;
		s[i] = pow (2.0, -h/c (psi));
	}

	freeForest (forest, t);
	return s;
}


/****************************************************************************************************************/
static double c (int psi) {
	if (psi <=1) 
		return 0;
	if (psi == 2) 
		return 1;
	return 	2*(log (psi-1.0) + 0.5772156649) - 2.0*(psi-1.0)/psi;
}

static Matrix_t* sample  (Matrix_t *X, int psi) {
	int 	subset[psi+10];	
	int 	i;
	int 	n = X->rowDim;

	for (i=0; i<psi; i++) {
		subset[i] = rand () % n;
	}
	return Matrix_sub (X, subset, psi);
}

static Matrix_t* filterL (Matrix_t* X, int q, double p) {
	int	i;
	int	rowDim = X->rowDim;
	int	colDim = X->colDim;
	int 	subset[rowDim+10], n = 0;

	for (i=0; i<rowDim; i++) if (X->elements[i*colDim+q] < p) {
		subset[n++] = i;
	}
	return Matrix_sub (X, subset, n);
}

static Matrix_t* filterR (Matrix_t* X, int q, double p) {
	int	i;
	int	rowDim = X->rowDim;
	int	colDim = X->colDim;
	int 	subset[rowDim+10], n = 0;

	for (i=0; i<rowDim; i++) if (X->elements[i*colDim+q] >= p) {
		subset[n++] = i;
	}
	return Matrix_sub (X, subset, n);
}

static iNode_t** iForest (Matrix_t *X, int t, int psi) {
	int		i;
	Matrix_t*	S      = NULL;	
	iNode_t** 	Forest = (iNode_t**)calloc (t+10, sizeof(iNode_t*));

	for (i=0; i<t; i++) {
		S = sample (X, psi);	
		Forest[i] = iTree (S);
		Matrix_free (&S);
	}
	return Forest;
}

static iNode_t* iTree (Matrix_t *X) {
	iNode_t*	pNode  = NULL;
	int		rowDim = X->rowDim;
	int		colDim = X->colDim;
	int		q;
	int		Q[colDim+10], n=0;
	double 		p;
	Matrix_t*	maxV = NULL;
	Matrix_t*	minV = NULL;
	Matrix_t*	S = NULL;
	int		i;
	
	if (X == NULL || X->rowDim < 1 || X->elements == NULL)
		return NULL;
	if (rowDim == 1) {
		pNode = (iNode_t*)malloc (sizeof (iNode_t));
		pNode->lchild = NULL;
		pNode->rchild = NULL;
		pNode->splitAtt = rand () % colDim;
	        pNode->splitValue = X->elements[pNode->splitAtt];		
		pNode->size = 1;
		return pNode;
	} 

	maxV = Matrix_max (X); 	
	minV = Matrix_min (X);
	for (i=0, n=0; i<colDim; i++) {
		if (maxV->elements[i] - minV->elements[i] > DBL_EPSILON) {
			Q[n++] = i;
		}
	}
	
	if (n == 0) {
		Matrix_free (&maxV);
		Matrix_free (&minV);
		pNode = (iNode_t*)malloc (sizeof (iNode_t));
		pNode->lchild = NULL;
		pNode->rchild = NULL;
		pNode->splitAtt = rand () % colDim;
	        pNode->splitValue = X->elements[pNode->splitAtt];		
		pNode->size = rowDim;
		return pNode;
	}

	q = Q[rand()%n];
	p = randu() * (maxV->elements[q] - minV->elements[q]) + minV->elements[q];
	Matrix_free (&maxV);
	Matrix_free (&minV);

	pNode = (iNode_t*)malloc (sizeof (iNode_t));

	S = filterL (X, q, p);
	pNode->lchild = iTree (S);
	Matrix_free (&S);

	S = filterR (X, q, p);
	pNode->rchild = iTree (S);
	Matrix_free (&S);

	pNode->splitAtt = q;
        pNode->splitValue = p;
	pNode->size = rowDim;
	return pNode;
}

static void freeTree (iNode_t *T) {
	if (T == NULL)
		return;
	if (T->lchild == NULL || T->rchild == NULL) {
		free (T);
		T = NULL;
		return;
	}
	freeTree (T->lchild);
	T->lchild = NULL;
	freeTree (T->rchild);
	T->rchild = NULL;
	free (T);
	T = NULL;
}

static void freeForest (iNode_t **forest, int t) {
	int 	i;
	
	for (i=0; i<t; i++) {
		freeTree (forest[i]);
	}
	free (forest);
}

static double pathLength (double* x, iNode_t* T, int hlim, int e) {
	int	a;
	if (T == NULL)
		return 0;
	if (T->lchild == NULL || T->rchild == NULL || e >= hlim)
		return e + c(T->size); 
	a = T->splitAtt;	
	if (x[a] < T->splitValue) {
		return pathLength (x, T->lchild, hlim, e+1);
	} else {
		return pathLength (x, T->rchild, hlim, e+1);
	}
}
