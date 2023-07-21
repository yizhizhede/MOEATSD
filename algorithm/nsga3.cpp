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
#include "nsga3.h"

Population_t* nsga3 (Problem_t *problem) {
	Matrix_t* Z = Population_reference (problem);
	int 	nPop = Z->rowDim;

	Population_t *P = Population_new (problem, nPop);
	Population_t *Q=NULL; 

	while (!isTerminal (P)) {
		Q = Population_reproduce(P, 'R'); 

		Population_cat (&Q, P);
		Population_free (&P);

		P = nsga3_eliminate (Q, Z);
		Population_free (&Q);
	}

	return P;
}

/********************************************************************************************************/
/********************************************************************************************************/
/********************************************************************************************************/
static void nsga3_normalize (Matrix_t *S) {
	Matrix_t *zmin = Matrix_min (S);	
	Matrix_t *xp   = Matrix_extreme (S);
	Matrix_t *inv  = Matrix_inverse (xp);
	int i, j;
	double a[S->colDim], *p;

	for (i=0; i<S->colDim; i++) {
		for (j=0, a[i]=0; j<S->colDim; j++) {
			a[i] += inv->elements[i*inv->colDim+j];
		}
	}

	for (i=0, p=S->elements; i<S->rowDim; i++) {
		for (j=0; j<S->colDim; j++, p++) {
			*p = (*p - zmin->elements[j]) * a[j];
		}
	}

	Matrix_free (&inv);
	Matrix_free (&xp);
	Matrix_free (&zmin);
}

static void nsga3_associate (Matrix_t *S, Matrix_t *Z, int *pi, double *d, int *arr) {
	int i, j, k;
	double *s=NULL, *w=NULL;
	double t;

	for (i=0, s=S->elements; i<S->rowDim; i++, s+=S->colDim) {
		k = arr[i];
		pi[k] = 0;
		d[k]  = distance_p2l (s, Z->elements, S->colDim);
		for (j=1, w=Z->elements+Z->colDim; j<Z->rowDim; j++, w+=Z->colDim) {
			t = distance_p2l (s, w, S->colDim);
			if (t < d[k]) {
				pi[k] = j;
				d[k]  = t;
			}
		}
	}
}

static void nsga3_niching (int K, int *rho, int n, int *pi, double *d, int m, Link_t **L1, Link_t *L2) {
	int k=0, i, j;
	int J[n], nJ=0, jbar;
	int I[n], nI=0;
	int vis[n];
	Node_t *p=NULL;
	double xValue=0;
	
	memset (vis, 0, n*sizeof (int));
	while (k<K) {
		nJ = 0;
		for (i=0; i<=m; i++) {
			for (j=0; j<n; j++) if (!vis[j]) {
				if (rho[j] == i) {
					J[nJ++] = j;
				}
			}
			if (nJ>0) 
				break;
		}
		jbar = J[(int)(nJ*randu ())];	

		nI = 0;
		p = L2->link_head;
		while (p) {
			if (pi[p->data] == jbar) {
				I[nI++] = p->data;
			}
			p=p->next;
		}

		if (nI > 0) {
			if (rho[jbar]==0) {
				for (i=0, xValue=1.0e+100, j=-1; i<nI; i++) {
					if (d[I[i]] < xValue) {
						j = I[i];
						xValue = d[I[i]];
					}
				}
				Link_add (L1, j);
				Link_del (L2, j);
			} else {
				j = I[(int)(randu()*nI)];
				Link_add (L1, j);
				Link_del (L2, j);
			}

			rho[jbar]++;
			k++;
		} else {
			vis[jbar] = 1;
		}
	}
}

Population_t *nsga3_eliminate (Population_t* Q, Matrix_t *Z) {
	Link_t 	*L1 = NULL, *L2=NULL, *L3=NULL;
	int 	nPop = Q->obj->rowDim / 2; 
	int 	K=nPop; 
	int 	*arr=NULL;
	Population_t 	*P = NULL;
	Matrix_t 	*S = NULL;
	int    	pi[Q->obj->rowDim];
	double 	d[Q->obj->rowDim];
	int 	rho[Z->rowDim];
	Node_t 	*p = NULL;

	Population_layer (Q, &L1, &L2, &L3);

	if (L3->nNode == nPop) {
		P = Population_sub (Q, L3);
	} else {
		K = (L1==NULL) ? nPop : (nPop - L1->nNode);	

		arr = Link2Array (L3);	
		S = Matrix_sub (Q->obj, arr+1, arr[0]);
		nsga3_normalize (S);
		nsga3_associate (S, Z, pi, d, arr+1);
		free (arr);
		arr = NULL;

		memset (rho, 0, Z->rowDim * sizeof (int));	

		p = (L1==NULL) ? NULL: (L1->link_head);
		while (p != NULL) {
			rho[pi[p->data]]++;
			p=p->next;
		}

		nsga3_niching (K, rho, Z->rowDim, pi, d, Q->obj->rowDim, &L1, L2);
		P = Population_sub (Q, L1);
	}

	Matrix_free (&S);
	Link_free (&L1);
	Link_free (&L2);
	Link_free (&L3);
	
	return P;
}



