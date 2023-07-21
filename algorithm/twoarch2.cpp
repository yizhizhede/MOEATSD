#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "algorithm.h"
#include "mystring.h"
#include "dominate.h"
#include "link.h"
#include "algebra.h"
#include "myrandom.h"
#include "iepsilon.h"
#include "terminal.h"

Population_t* two_arch2_updateCA (Population_t* Q, int nca);	// update CA
Population_t* two_arch2_updateDA (Population_t* Q, int nda);	// update DA

Population_t* two_arch2 (Problem_t *problem) {
	int nca = 10;
	int nda = 100;

	nda = Population_getSize (problem);
	nca = nda / 10;

	Population_t *CA = Population_new (problem, nca);
	Population_t *DA = Population_new (problem, nda);
	Population_t *Q=NULL, *U=NULL; 

	while ( !isTerminal(DA)) {
		Population_cat (&U, CA);
		Population_cat (&U, DA);
		Q = Population_reproduce(U, 'R'); 	// reproduction of CA and DA
		Population_free (&U);
		
		U = Population_mutation(CA); 		// mutation of CA
		Population_cat (&Q, U);
		Population_free (&U);

		Population_cat (&Q, CA);
		Population_cat (&Q, DA);
		Population_free (&CA);
		Population_free (&DA);

		CA = two_arch2_updateCA (Q, nca);	// update CA
		DA = two_arch2_updateDA (Q, nda);	// update DA

		Population_free (&Q);

	}
	
	return DA;
}

/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
Population_t* two_arch2_updateCA (Population_t* Q, int nca) { 	// update CA
	double 		*fitness = Iepsilon_fitness (Q->obj);
	Link_t 		*link=NULL;
	Population_t 	*CA=NULL;
	int sel[Q->obj->rowDim], i, num=Q->obj->rowDim, r;
	double xValue;

	for (i=0; i<Q->obj->rowDim; i++) sel[i] = 1;
	while (num > nca) {
		for (i=0, xValue=1.0e+100, r=-1; i<Q->obj->rowDim; i++) if (sel[i] == 1) {
			xValue=fitness[i];
			r = i;
			break;
		}
		for (i=r+1; i<Q->obj->rowDim; i++) if (sel[i] == 1) {
			if (fitness[i] < xValue) {
				xValue = fitness[i];
				r = i; 
			}
		}

		fitness[r] = 0;
		sel[r] = 0;
		Ipesilon_update  (Q->obj, fitness, sel, r);
		num--;
	}

	
	for (i=0; i<Q->obj->rowDim; i++) if (sel[i] == 1) {
		Link_add (&link, i);	
	}

	CA = Population_sub (Q, link);	

	Link_free (&link);
	free (fitness);
	return CA;
}

Population_t* two_arch2_updateDA (Population_t* Q, int nda) {	// update DA
	Population_t 	*DA=NULL;
	Population_t 	*U=Population_Front1 (Q);
	Matrix_t	*normU = Matrix_norm (U->obj);
	Link_t 		*link=NULL;
	int sel[Q->obj->rowDim], i, j, num=0, r;
	double xValue, *x1=NULL, *x2=NULL, *x3=NULL;
	double similarity[Q->obj->rowDim];

	if (U->obj->rowDim <= nda)	return U;

	memset (sel, 0, U->obj->rowDim*sizeof (int));
	for (j=0; j<U->obj->colDim; j++) {
		xValue = U->obj->elements[j];
		r = 0;
		for (i=1, x1=U->obj->elements+i*U->obj->colDim; i<U->obj->rowDim; i++, x1+=U->obj->colDim) {
			if (x1[j]<xValue) {
				xValue = x1[j];
				r = i;
			}
		}
		sel[r] = 1;
	}
	for (i=0, num=0; i<U->obj->rowDim; i++) num+=sel[i];
	
	for (i=0; i<U->obj->rowDim; i++) if (!sel[i]) {
		x1 = normU->elements+i*U->obj->colDim;
		for (j=0; j<U->obj->rowDim; j++) if (sel[j]) {
			x2 = normU->elements+j*U->obj->colDim;
			x3 = vector_substract (x1, x2, U->obj->colDim);
			similarity[i] = Lp_norm (x3, U->obj->colDim, 1.0/U->obj->colDim);	
			free (x3);
			break;
		}
		for (j++; j<U->obj->rowDim; j++) if (sel[j]) {
			x2 = normU->elements+j*U->obj->colDim;
			x3 = vector_substract (x1, x2, U->obj->colDim);
			xValue = Lp_norm (x3, U->obj->colDim, 1.0/U->obj->colDim);	
			free (x3);
			if (xValue < similarity[i]) {
				similarity[i] = xValue;
			}
		}
	}

	while (num < nda) {
		for (i=0, xValue=-1.0e+100, r=-1; i<U->obj->rowDim; i++) if (!sel[i]) {
			xValue = similarity[i];
			r = i;
			break;
		}
		for (i++; i<U->obj->rowDim; i++) if (!sel[i]) {
			if (similarity[i]>xValue) {
				xValue = similarity[i];
				r = i;
			}
		}
		sel[r] = 1;

		num++;
		similarity[r] = 0.0;
		x2 = normU->elements+r*U->obj->colDim;
		for (i=0; i<U->obj->rowDim; i++) if (!sel[i]) {
			x1 = normU->elements+i*U->obj->colDim;
			x3 = vector_substract (x1, x2, U->obj->colDim);
			xValue = Lp_norm (x3, U->obj->colDim, 1.0/U->obj->colDim);	
			free (x3);
			if (xValue < similarity[i]) {
				similarity[i] = xValue;
			}
		}
	}
	
	for (i=0; i<U->obj->rowDim; i++) if (sel[i] == 1) {
		Link_add (&link, i);	
	}
	DA = Population_sub (U, link);	

	Link_free (&link);
	Population_free (&U);
	Matrix_free (&normU);

	return DA;
}
