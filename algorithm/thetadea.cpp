#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "algorithm.h"
#include "mystring.h"
#include "dominate.h"
#include "link.h"
#include "algebra.h"
#include "terminal.h"
#include "myrandom.h"

static void  thetadea_normalize (Matrix_t *S, Matrix_t *ideal, Matrix_t *nadir);
static List_t* theta_non_dominated_sort (Matrix_t *S, Matrix_t *Z);	
static int*  thetadea_select (List_t *list, int nPop);

Population_t* thetadea (Problem_t *problem) {
	int 	*Sel=NULL, *arr=NULL, i;
	Link_t 	*L1 = NULL, *L2=NULL, *L3=NULL;
	List_t  *list = NULL;
	Matrix_t* S = NULL;
	Matrix_t* Z = Population_reference (problem);
	int 	nPop = Z->rowDim;

	Population_t *P = Population_new (problem, nPop);
	Population_t *Q=NULL; 

	Matrix_t *ideal = Matrix_min (P->obj);
	Matrix_t *nadir = Matrix_max (P->obj);

	while (!isTerminal (P)) {
		Q = Population_reproduce(P, 'R');	// reproduction

		Population_cat (&Q, P);			// combination
		Population_free (&P);			

		Population_layer (Q, &L1, &L2, &L3);	// layer population 
		S = Matrix_sub (Q->obj, L3);
	
		Matrix_free (&ideal);
		ideal = Matrix_min (S);

		thetadea_normalize (S, ideal, nadir);	// normalize

		list = theta_non_dominated_sort (S, Z);	// theta_non_domination_sort

		Sel=thetadea_select (list, nPop);	// select offspring

		arr = Link2Array (L3);
		for (i=0; i<nPop; i++) {
			Sel[i] = arr[Sel[i]+1];
		}
		P = Population_sub (Q, Sel, nPop);	

		free (Sel);
		free (arr);
		Link_free (&L1);
		Link_free (&L2);
		Link_free (&L3);
		List_free (&list);
		Matrix_free (&S);
		Population_free (&Q);

	}

	Matrix_free (&nadir);
	Matrix_free (&ideal);
	Matrix_free (&Z);

	return P;
}

/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
static double ASF (double *f, int m, int j, Matrix_t *ideal, Matrix_t *nadir) {
	int i;	
	double xValue=0, t;

	for (i=0; i<m; i++) {
		t = nadir->elements[i] - ideal->elements[i];
		if (t < 1.0e-6) 
			continue;
		t = (i==j?1:1.0e+6)*(f[i] - ideal->elements[i]) / t;
		if (t<0)
			t = -t;
		if (t>xValue)
			xValue = t;
	}
	return xValue;
}

static Matrix_t* thetadea_extreme (Matrix_t *S, Matrix_t *ideal, Matrix_t *nadir) {
	Matrix_t *M=NULL;
	Link_t *link = NULL;
	double *f=NULL, xValue, t;
	int i, j, k;
	for (j=0; j<S->colDim; j++) {
		for (i=0, f=S->elements, k=0; i<S->rowDim; i++, f+=S->colDim) {
			xValue = ASF (f, S->colDim, j, ideal, nadir);
			k = i;
			break;
		}
		for (i++, f+=S->colDim; i<S->rowDim; i++, f+=S->colDim) {
			t = ASF (f, S->colDim, j, ideal, nadir);
			if (t < xValue) {
				xValue = t;
				k = i;
			}
		}
		Link_add (&link, k);
	}
	
	M = Matrix_sub (S, link);

	Link_free (&link);
	return M;
}

static void thetadea_normalize (Matrix_t *S, Matrix_t *ideal, Matrix_t *nadir) {
	int i, j;
	Matrix_t *E = thetadea_extreme (S, ideal, nadir);
	Matrix_t *invE = NULL;
	Matrix_t *invEu = Matrix_new (E->rowDim, 1);
	Matrix_t *xV = Matrix_max (S);
	double t;
	

	for (i=0; i<E->rowDim; i++) {
		for (j=0; j<E->colDim; j++) {
			E->elements[i*E->colDim+j] -= ideal->elements[j];
		}
	}

	invE = Matrix_inverse (E);

	for (i=0; i<E->rowDim; i++) {
		for (j=0; j<E->colDim; j++) {
			invEu->elements[i] += invE->elements[i*E->colDim+j];
		}
		
		if (invEu->elements[i]>=1.0e-6 ) {
			nadir->elements[i] =  1.0/invEu->elements[i] + ideal->elements[i];
		} else {
			nadir->elements[i] = xV->elements[i];
		}

	}
	
	for (i=0; i<S->rowDim; i++) {
		for (j=0; j<S->colDim; j++) {
			S->elements[i*S->colDim+j] -= ideal->elements[j];

			t = nadir->elements[j] - ideal->elements[j];
			if (t > 1.0e-6)
				S->elements[i*S->colDim+j] /= (nadir->elements[j] - ideal->elements[j]);
		}
	}

	Matrix_free (&xV);
	Matrix_free (&invEu);
	Matrix_free (&invE);
	Matrix_free (&E);
}

static void thetadea_clustering (Matrix_t *S, Matrix_t *Z, int *C, double *Fj) {
	int i, j, k;
	double *p=NULL, *l=NULL;
	double xValue = 0, t, d2, d1;

	for (i=0, p=S->elements; i<S->rowDim; i++, p+=S->colDim) {
		for (j=0, l=Z->elements, k=0; j<Z->rowDim; j++, l+=Z->colDim) {
			xValue = distance_p2l (p, l, Z->colDim);
			k = j;
			break;
		}
		for (j++, l+=Z->colDim; j<Z->rowDim; j++, l+=Z->colDim) {
			t = distance_p2l (p, l, Z->colDim);
			if (t < xValue) {
				xValue = t;
				k = j;
			}
		}
		C[i] = k;
		d2 = xValue;
		d1 = inner_product (p, Z->elements+k*Z->colDim, Z->colDim) / norm (Z->elements+k*Z->colDim, Z->colDim); 
		for (j=0; j<Z->colDim && Z->elements[k*Z->colDim+j]<1.0; j++);
		Fj[i] = d1 + (j<Z->colDim?1.0e+6:5.0) * d2;
	}
}

static List_t* theta_non_dominated_sort (Matrix_t *S, Matrix_t *Z) {
	int *	C  = (int *)malloc (S->rowDim*sizeof (int));
	double*	Fj = (double *)malloc (S->rowDim*sizeof (double));
	List_t* list = NULL;
	Link_t* link = NULL;
	int 	cnt=0, i, j, k;
	double 	xValue;


	thetadea_clustering (S, Z, C, Fj); 

	while (cnt < S->rowDim) {
		for (i=0; i<Z->rowDim; i++) {
			xValue = -1.0;
			k = -1; 
			for (j=0; j<S->rowDim; j++) if (C[j] == i) {
				xValue = Fj[j];	
				k = j;
				break;
			}
			for (j++; j<S->rowDim; j++) if (C[j] == i) {
				if (Fj[j] < xValue) {
					xValue = Fj[j];	
					k = j;
				}
			}
			if (xValue < 0) 
				continue;
			Link_add (&link, k);
			C[k] = -1;
			cnt++;
		}
		List_add (&list, link);
		link = NULL;

	}


	free (C);
	free (Fj);
	return list;	
}

static int*  thetadea_select (List_t *list, int nPop) {
	Link_t *link=NULL;
	Node_t *p = NULL;
	int *Sel = (int *)malloc (nPop*sizeof (int)), nSel = 0;
	int *index = NULL, *arr=NULL;
	Matrix_t *M = NULL;
	int i;



	link = list->list_head;
	while (link && link->nNode + nSel <= nPop) {
		p = link->link_head;
		while (p) {
			Sel[nSel++] = p->data;
			p = p->next;
		}
		link = link->next;
	}
	
	if (nSel < nPop && link) {
		arr = Link2Array (link);	
		M = Matrix_new (arr[0], 1);
		for (i=0; i<M->rowDim; i++) {
			M->elements[i] = randu ();
		}
		index=sort (M);

		for (i=nPop-nSel; i>0; i--) {
			Sel[nSel++]=arr[index[i]+1];
		}
	
		Matrix_free (&M);
		free (index);
		free (arr);
	}
	return Sel;
}
