#include "reduce.h"
#include "rank.h"
#include "link.h"
#include "dominate.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>


static int ndp (Link_t *link, Matrix_t *M); 	// number of dominating points
static int chv (Matrix_t *M); 			// contributing hypervolume 

Population_t *Reduce_one_by_hv (Population_t *Q) {
	Population_t *P = NULL;
	List_t *list = ndSort (Q->obj);
	int sub[Q->obj->rowDim], top=0, r, i;
	
	if (list->nLink > 1) {
		r = ndp (list->list_tail, Q->obj);
	} else {
		r = chv (Q->obj);	
	}
	
	for (i=0; i<Q->obj->rowDim; i++) {
		if (i==r) 
			continue;
		sub[top++]=i;
	}
	P = Population_sub (Q, sub, top);	

	List_free (&list);
	return P;
}

Population_t *Reduce_hal_by_cd (Population_t *Q) {       // reduce half of individuals by crowing distance
	return rank_by_crowding (Q, Q->obj->rowDim / 2);
}

static int ndp (Link_t *link, Matrix_t *M) { 	// number of dominating points
	int i, j, d, r=0, xValue=0;
	Node_t *p = NULL;

	p = link->link_head;
	while (p) {
		i = p->data;
		d = 0;
		for (j=0; j<M->rowDim; j++) {
			if (i==j) continue;
			if (isDominate (M->elements+j*M->colDim, M->elements+i*M->colDim, M->colDim)) {
				d++;
			}
		}
		if (d > xValue) {
			xValue = d;
			r = i;
		}

		p = p->next;
	}
	return r;
}

static double *exchv3d (Matrix_t *M) {
	double *a = (double *)malloc ((M->rowDim+2)*sizeof (double));
	double *b = (double *)malloc ((M->rowDim+2)*sizeof (double));
	double *chv = (double *)malloc ((M->rowDim+2)*sizeof (double));

	int *index = NULL;
	int i, j, k, stem;
	double top1, top2, *p;

	index = sort (M, 0);
	for (i=0; i<M->rowDim; i++) {
		a[i] = M->elements[3*index[i]];
	}
	free (index);

	index = sort (M, 1);
	for (i=0; i<M->rowDim; i++) {
		b[i] = M->elements[3*index[i]+1];
	}
	free (index);
	
	a[i] = b[i] = 1.1;
	
	memset (chv, 0, (M->rowDim+2)*sizeof (double));
	for (i=0; i<M->rowDim; i++) {
		for (j=0; j<M->rowDim; j++) {
			top1 = top2 = 1.1;
			stem = -1;
			for (k=0, p=M->elements; k<M->rowDim; k++, p +=3) {
				if (p[0] - a[i] <= 0 && p[1] - b[j] <= 0) {
					if (p[2] < top1) {
						top2 = top1;
						top1 = p[2];
						stem = k;
					} else if (p[2] - top2 <= 0) {
						top2 = p[2];
					}
				}
			}
			if (stem >=0) 
				chv[stem] += (a[i+1]-a[i])*(b[j+1]-b[j])*(top2-top1);
		}
	}

	free (a);
	free (b);
	return chv;
}

static int chv (Matrix_t *M) { 			// contributing hypervolume 
	Matrix_t *normM = Matrix_norm (M);
	int i, *index = NULL, r=0;
	double exchv = 0, xValue=1.0e+10;
	double *buff=NULL;

	if (M->colDim == 2) {
		index = sort (normM);
		for (i=1; i<M->rowDim-1; i++) {
			exchv = (normM->elements[index[i+1]*2]-normM->elements[index[i]*2])*
							(normM->elements[index[i-1]*2+1]-normM->elements[index[i]*2+1]);
			if (exchv < xValue) {
				xValue = exchv;
				r = index[i];
			}
		}
		free (index);
	} else if (M->colDim == 3){
		buff = exchv3d (normM);
		for (i=0; i<M->rowDim; i++) {
			if (buff[i] < xValue) {
				xValue = buff[i];
				r = i;
			}
		}
		free (buff);
	} else {
		r=0;
	}

	Matrix_free (&normM);	
	return r;
}
	

