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
#include "ibea.h"

Population_t* ibea (Problem_t *problem) {
	int nPop = Population_getSize (problem);

	Population_t *pop=NULL; 
	Population_t *Q = Population_new (problem, nPop);

	while (true) {
		pop = ibea_environmental_select (Q, nPop);
		Population_free (&Q);
		if (isTerminal(pop))
			break;

		Q = Population_reproduce(pop, 'T'); 	  
		Population_cat (&Q, pop);
		Population_free (&pop);
	}
	
	return pop;
}

/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
static double c = 0.0;
static double k = 0.05;
static double (*I)(double*, double*, int);
Matrix_t *ibea_fitness_assign (Matrix_t *Obj) {
	Matrix_t *normObj = Matrix_norm (Obj);
	Matrix_t *fitness = Matrix_new (Obj->rowDim, 1);
	I = Iepsilon_plus;
	int i, j;
	double *x1 = NULL, *x2=NULL, t;

	for (i=0, x1=normObj->elements, c=0.0; i<Obj->rowDim; i++, x1 += Obj->colDim) {
		for (j=0, x2=normObj->elements; j<Obj->rowDim; j++, x2 += Obj->colDim) {
			if (i==j) 
				continue;
			t = I(x1, x2, Obj->colDim);
			if (t < 0)
				t = -t;
			if (t > c)
				c = t;
		}
	}

	for (i=0, x1=normObj->elements; i<Obj->rowDim; i++, x1 += Obj->colDim) {
		for (j=0, x2=normObj->elements; j<Obj->rowDim; j++, x2 += Obj->colDim) {
			if (i==j) 
				continue;
			fitness->elements[i] += -exp (-I(x2, x1, Obj->colDim)/(c*k));
		}
	}

	Matrix_free (&normObj);
	return fitness;
}

void  ibea_fitness_update (Matrix_t *Obj, Matrix_t *fitness, int r) {
	Matrix_t 	*normObj = Matrix_norm (Obj);
	int 		i;
	double 		*x1 = NULL, *x2=NULL;

	x2=normObj->elements+r*Obj->colDim;
	for (i=0, x1=normObj->elements, c=0.0; i<Obj->rowDim; i++, x1 += Obj->colDim) {
		fitness->elements[i] += exp (-I(x2, x1, Obj->colDim)/(c*k));
	}

	Matrix_free (&normObj);
}

Population_t* ibea_environmental_select (Population_t* Q, int nPop) {
	Population_t *subQ=NULL;
	Matrix_t *fitness = ibea_fitness_assign (Q->obj);	
	Link_t	*link = NULL;
	int 	sel[Q->obj->rowDim], nSel=Q->obj->rowDim, i, r;
	double 	xValue;
	int 	*index=NULL;

	for (i=0; i<Q->obj->rowDim; i++) sel[i] = 1;

	while (nSel > nPop) {
		for (i=0, xValue=1.0e+100, r=-1; i<Q->obj->rowDim; i++) if (sel[i] == 1) {
			xValue = fitness->elements[i];
			r = i;
			break;
		}
		for (i++; i<Q->obj->rowDim; i++) if (sel[i] == 1) {
			if (fitness->elements[i] < xValue) {
				xValue = fitness->elements[i];
				r = i;
			}
		}

		sel[r] = 0;
		ibea_fitness_update (Q->obj, fitness, r);
		nSel--;
	}

	index = sort (fitness, "DES");	
	for (i=0; i<Q->obj->rowDim; i++) {
		r = index[i];
		if (sel[r] == 1) {
			Link_add (&link, r);
		}
	}

	subQ = Population_sub (Q, link);		

	free (index);
	Link_free (&link);
	Matrix_free (&fitness);
	
	return subQ;
}
