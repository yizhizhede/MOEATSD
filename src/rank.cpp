#include "rank.h"
#include "crowding.h"
#include "dominate.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

Population_t *rank_by_crowding (Population_t *pop) {
	return rank_by_crowding (pop, pop->obj->rowDim);
}
 
Population_t *rank_by_crowding (Population_t *pop, int T) {
	int*		index = rank_by_crowding (pop->obj, T);
	Population_t*	subPop = Population_sub (pop, index, T);

	free (index);
	return subPop;
}

int *rank_by_crowding (Matrix_t *M) {
	return rank_by_crowding (M, M->rowDim);
}

int *rank_by_crowding (Matrix_t *M, int T) {
	List_t*		list = ndSort (M);
	Link_t*		link = NULL;
	Matrix_t*	subM = NULL;
	Matrix_t*	disM = NULL; 
	int*		index = NULL;
	int*		global = NULL;
	int*		queue = (int *)malloc (M->rowDim *sizeof (int));
	int 		tail=0, i;

	link = list->list_head;
	while (link) {
		subM = Matrix_sub (M, link);
		disM = Matrix_new (subM->rowDim, 1);
		free (disM->elements);
		disM->elements = crowding_distance (subM);		
		index = sort (disM);
		global = Link2Array (link);

		for (i=disM->rowDim-1; i>=0; i--) {
			queue[tail++] = global[index[i]+1];
		}

		free (global);
		free (index);
		Matrix_free (&disM);
		Matrix_free (&subM);

		if (tail >= T)
			break;
		link = link->next;
	}
	
	List_free (&list);
	return queue;
}

int *rank_by_density (Matrix_t *M) { 			// non-dominated sort + Shift-based density estimation (SDE)
	int		rowDim 	= M->rowDim;
	int 		colDim 	= M->colDim;
	Matrix_t*	D	= Matrix_new (rowDim, 3);	// NO. of front, -SDE, NO.
	int*		I	= NULL;
	List_t*		list 	= ndSort (M);
	Link_t*		front 	= NULL;
	int*		arr 	= NULL;
	double		t, d, SDE;
	int		i, j, k, a, b, NF;

	NF = 1;
	front = list->list_head;
	while (front) {
		arr = Link2Array (front);
		// compute SDE
		for (a=arr[0]; a>0; a--) {
			i = arr[a];
			SDE = 1.0e+100;
			for (b=arr[0]; b>0; b--) {
				j = arr[b];
				if (i==j) continue;
				for (k=0, d=0; k<colDim; k++) {
					t = M->elements[j*colDim+k] - M->elements[i*colDim+k];
					if (t > 0) {
						d += t*t;
					}
				}
				d = sqrt (d);
				if (d < SDE) {
					SDE = d;
				}
			}
			D->elements[i*3+0] = NF;	
			D->elements[i*3+1] = -SDE;	
			D->elements[i*3+2] = i;	
		}
		// keep one with minimun 
		for (k=0; k<colDim; k++) {
			d = 1.0e+100;
			b = -1;
			for (a=arr[0]; a>0; a--) {
				i = arr[a];
				if (M->elements[i*colDim+k] < d) {
					d = M->elements[i*colDim+k];
					b = i;
				}
			}
			if (-1 != b) {
				D->elements[b*3+1] = -1.0e+100;	
			}
		}
		// free arr
		free (arr);

		//
		NF++;
		front = front->next;
	}
	I = sort (D);

	// free
	List_free (&list);
	Matrix_free (&D);

	return I;
}
