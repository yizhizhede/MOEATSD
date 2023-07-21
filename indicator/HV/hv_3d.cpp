#include "hv.h"
#include "avltree.h"
#include <stdlib.h>
#include <string.h>

double hv_3d (Matrix_t *M) {
	if (M->colDim != 3) return -1.0;
	int* 		index = sort (M, 2);
	avl_tree_t 	T = NULL;
	avl_node_t	*node=NULL, *succ = NULL, *pred = NULL, *s = NULL;
	double 		volume = 0, z=0, A=0, height, width;
	int		i;

	typedef struct point_tag {
		double x, y, z;
	} point_t;

	size_t size = sizeof (point_t);
	point_t *p = (point_t *)malloc (size);
	p->x = M->elements[index[0]*3];
	p->y = M->elements[index[0]*3 + 1];
	p->z = M->elements[index[0]*3 + 2];
	node = avl_insert (&T, p->x, p);
	// avl_visualize (T);
	A = (1.1 - p->x)*(1.1 - p->y);
	z = p->z;

	for (i=1; i<M->rowDim; i++) {
		p = (point_t *)malloc (size);
		p->x = M->elements[index[i]*3];
		p->y = M->elements[index[i]*3 + 1];
		p->z = M->elements[index[i]*3 + 2];
		node = avl_insert (&T, p->x, p);
		// avl_visualize (T);
		pred = avl_pred (node);
		if (pred==NULL || ((point_t *)pred->value)->y > p->y) {
			volume += A * (p->z - z);
			z = p->z;
			succ = avl_succ (node);
			while (succ && !(((point_t *)succ->value)->y < p->y)) {
				s = avl_succ (succ);
				if (s)
					width = ((point_t *)s->value)->x - ((point_t *)succ->value)->x;
				else 
					width = 1.1 - ((point_t *)succ->value)->x;
				if (pred)
					height = ((point_t *)pred->value)->y - ((point_t *)succ->value)->y;
				else
					height = 1.1 - ((point_t *)succ->value)->y;
				A -= width * height;	

				s = succ;
				succ = avl_succ (succ);
				avl_delete (&T, s);
				// avl_visualize (T);
			}
			if (succ)
				width = ((point_t *)succ->value)->x - p->x;
			else 
				width = 1.1 - p->x;
			if (pred)
				height = ((point_t *)pred->value)->y - p->y;
			else
				height = 1.1 - p->y;
			A += width * height;	
		} else {
			avl_delete (&T, node);
			// avl_visualize (T);
		}
	}
	volume += A * (1.1 - z);	

	free (index);
	avl_tree_free (&T);

	return volume;
}
