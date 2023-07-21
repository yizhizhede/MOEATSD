#include "hv.h"
#include "dominate.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

Matrix_t *limitset (Matrix_t *M, Vector_t *p);

double hv_wfg (Matrix_t *M) {
	typedef struct hv_node_tag {
		int vis;
		double exclhv; 	// exclhv = inclhv - hv;
		double inclhv; 	// product {p[j] - refPoint[j]} 
		double hv; 	// sum exclhv (subPro)	

		Vector_t *p;
		Matrix_t *Ps;
	} hv_node_t;

	hv_node_t* stack[1000];
	int 	top = 0;
	hv_node_t *curNode = NULL, *parent = NULL, *child = NULL;
	Vector_t *p = NULL;
	Matrix_t *Ps = NULL;
	int subset[1000], len;
	int i;
	double hv;
	size_t size = sizeof (hv_node_t);

	curNode = (hv_node_t *)malloc (size);
	memset (curNode, 0, size);
	curNode->vis = 0;
	curNode->Ps = Matrix_dup (M);
	stack[0] = curNode;
	top = 1;

	while (top) {
		curNode = stack[top-1];
		if (curNode->vis < curNode->Ps->rowDim) {
			size = sizeof (hv_node_t);
			child = (hv_node_t *)malloc (size);
			memset (child, 0, size);
			child->vis = 0;
	
			p = Matrix_sub (curNode->Ps, curNode->vis);
			child->p = p;
			
			child->inclhv = 1.0;
			for (i=0; i<p->dim; i++) {
				child->inclhv *= (1.1 - p->elements[i]);
			}
	
			len=0;
			for (i=(curNode->vis+1); i<(curNode->Ps->rowDim); i++) {
				subset[len++] = i;
			}

			if (len == 0) {
				child->hv = 0;
				child->Ps = Matrix_new ();
			} else {
				Ps = Matrix_sub (curNode->Ps, subset, len);
				child->Ps = limitset(Ps, p);
				Matrix_free (&Ps);
			}

			stack[top++] = child;
			curNode->vis++;	
		} else {
			curNode->exclhv = curNode->inclhv - curNode->hv;
			if (top > 1) {
				parent = stack[top-2];
				parent->hv += curNode->exclhv;
			}

			hv=curNode->hv;
			if (curNode->p != NULL)
				Vector_free (&curNode->p);
			if (curNode->Ps != NULL)
				Matrix_free (&curNode->Ps);
			free (curNode);
			curNode = NULL;
			top--;
		}
	}
	return hv;
}

Matrix_t *limitset (Matrix_t *M, Vector_t *p) {
	int 	i, j, index;

	for (i=0, index=0; i<M->rowDim; i++) {
		for (j=0; j<M->colDim; j++) {
			if (M->elements[index] < p->elements[j]) {
				M->elements[index] = p->elements[j];
			}
			index++;
		}
	}
	
	return Matrix_front (M);
}
