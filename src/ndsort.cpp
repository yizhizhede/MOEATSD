#include "dominate.h"
#include "algebra.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#define NDSORT_FNS	0
#define NDSORT_TENS 	1
#define NDSORT_TYPE 	1


static List_t *fns  (Matrix_t *M, int minNum);
static List_t *tens (Matrix_t *M, int minNum);
static List_t *tens_1k (Matrix_t *M, int minNum);

List_t *ndSort (Matrix_t *M) {
	return ndSort (M, M->rowDim);
}

List_t *ndSort (Matrix_t *M, int minNum) {
	switch (NDSORT_TYPE) {
		case NDSORT_FNS: 
			return fns (M, minNum);
		case NDSORT_TENS:
			return tens (M, minNum);
		default: 
			return NULL;
	}
}


static List_t *fns  (Matrix_t *M, int minNum) {
    	int i, j, nrow = M->rowDim, ncol = M->colDim;
        Link_t*	S[nrow];   
        int 	n[nrow];    

        Node_t *p = NULL, *q = NULL;
	List_t *list = NULL;
	Link_t *curFront = NULL, *nextFront=NULL;
 
        for (i=0; i<nrow; i++) {
		S[i] = NULL;
                n[i] = 0;
                for (j=0; j<nrow; j++) {
			if (i==j) 
				continue;
                       	if (isDominate(M->elements+i*ncol, M->elements+j*ncol, ncol) == 1) {
				Link_add (&S[i], j);
                        } else if (isDominate(M->elements+j*ncol, M->elements+i*ncol, ncol) == 1) {
                                n[i]++;
                        }
                }
                if (n[i]==0) {
			Link_add (&curFront, i);
               	}
       	}
	List_add (&list, curFront); 

        while (curFront != NULL && curFront->nNode > 0) {
             	while ((p = Link_read (curFront)) != NULL) {
                      	while ((q = Link_read (S[p->data])) != NULL ) {
                            	n[q->data]--;
                                if (n[q->data]==0) {
                                	Link_add( &nextFront, q->data);
                             	} 
                  	}
			Link_free (&S[p->data]);
       		}

		Link_free (&curFront);
		curFront = nextFront;
		nextFront = NULL;
		List_add (&list, curFront);
      	}
        return list;
}

typedef struct TENS_Node_tag {
	double*	data;
	struct 	TENS_Node_tag* child[100];
	struct 	TENS_Node_tag* parent;
	int 	vis; 
	int 	span;
}TENS_Node_t;

TENS_Node_t* TENS_Node_new () {
	TENS_Node_t*	node = NULL;
	node = (TENS_Node_t *)calloc (1, sizeof (TENS_Node_t));
	return node;
}

static List_t *tens (Matrix_t *M, int minNum) {
	int 		*Q = NULL;
	int 		*vis = NULL;
	int 		cnt = 0, index=0, m;
	List_t* 	list = NULL;
	Link_t* 	link = NULL;
	TENS_Node_t 	*tree=NULL, *p=NULL, *q=NULL;
	TENS_Node_t 	*add_pos=NULL;
	int 		add_no; 
	int 		i;
	double 		d;

	// if the number of points is less then 1000, use the stack memory to speed up.
	if (M->rowDim <= 1000) {
		return tens_1k (M, minNum);
	}

	//
	Q = sort (M);
	vis = (int *)malloc (M->rowDim*sizeof (int));
	memset (vis, 0, M->rowDim*sizeof (int));

	while (cnt < minNum && cnt < M->rowDim) {	// {F1, F2, ...}
		for (index=0; index<M->rowDim && vis[index]; index++);
		if (index >= M->rowDim) continue;
		tree = TENS_Node_new ();	
		tree->data = M->elements+M->colDim*Q[index];
		tree->span = M->colDim-1;
		tree->vis = 0;
		Link_add (&link, Q[index]);
		vis[index] = 1;
		cnt++;
			
		while (index < M->rowDim) {
			for (index++; index<M->rowDim && vis[index]; index++);
			if (index >=M->rowDim) continue;
			p = TENS_Node_new ();	
			p->data = M->elements+M->colDim*Q[index];
			p->span = M->colDim-1;
			p->vis = 0;

			// find position
			q = tree;
			while (q) {
				for (m=1; m<M->colDim && q->data[m] <= p->data[m]; m++);
				if (m >= M->colDim) {
					// there are the same solution
					for (i=0; i<M->colDim; i++) {
						d = q->data[i] - p->data[i];
						if (d < 0)
							d = -d;
						if (d > DBL_EPSILON)
							break;
					}
					// if (distance_p2p (q->data, p->data, M->colDim) < DBL_EPSILON) {
					if (i >= M->colDim) {
						Link_add (&link, Q[index]);
						vis[index] = 1;
						cnt++;
					}

					free (p);
					p = NULL;
					break; 
				} else {
					add_pos = q;
					add_no = m;
					q->vis = 1;	
					q->span = m;
					q = q->child[m];
				}
			}

			if (p==NULL) continue;
			
			// check other solution
			q = add_pos;
			while (q) {
				if (!q->vis) {
					q->vis = 1;
					for (m=1; m<M->colDim && q->data[m] <= p->data[m]; m++);
					if (m >= M->colDim) {
						// there are the same solution
						for (i=0; i<M->colDim; i++) {
							d = q->data[i] - p->data[i];
							if (d < 0)
								d = -d;
							if (d > DBL_EPSILON)
								break;
						}
						// if (distance_p2p (q->data, p->data, M->colDim) < DBL_EPSILON) {
						if (i >= M->colDim) {
							Link_add (&link, Q[index]);
							vis[index] = 1;
							cnt++;
						}
						free (p);
						p = NULL;
						break;
					}
				} 

				while (q->span > 0 && q->child[q->span]==NULL) {
					q->span--;
				}
				if (q->span > 0 && q->child[q->span]) { 
					q = q->child[q->span];
				} else {
					q->vis = 0;
					q->span = M->colDim - 1;
					q = q->parent;
					if (q) 
						q->span--;
				}
			}

			if (p != NULL) {
				p->parent = add_pos;
				add_pos->child[add_no] = p;
				Link_add (&link, Q[index]);				
				vis[index] = 1;
				cnt++;
			}
		}

		List_add (&list, link);
		Link_free (&link);
		link = NULL;
		
		q = tree;
		while (q) {
			for (m=1; m<M->colDim; m++) if (q->child[m]) {
				q = q->child[m];
				break;
			}
			if (m >= M->colDim) {
				if (q->parent) {
					for (m=1; q->parent->child[m]!= q; m++);
					p = q->parent;
					p->child[m]= NULL;
					free (q);
					q = p; 
				} else {
					free (q);
					q = NULL;
				}
			}
		}
		tree = NULL;

	}

	free (Q);
	free (vis);
	return list;
}


static List_t *tens_1k (Matrix_t *M, int minNum) {
	// the flag of visited opoints
	int 	vis[1010];
	int*	Q=NULL;
	int	index = 0;

	//
	int 	cnt = 0, m;
	List_t* list = NULL;
	Link_t* link = NULL;
	int  	tree, p, q;
	int  	add_pos = -1;
	int 	add_no = 0; 
	int 	i, j;
	double 	d;
	int	greater, equal;
	int 	rowDim = M->rowDim;
	int 	colDim = M->colDim;

	// the tree
	double	data[1010][20];
	int	parent[1010];
	int 	child[1010][20];
	int	span[1010];
	int 	tree_vis[1010];
	
	Q = sort (M);
	memset (vis, 0, 1010*sizeof (int));

	// initize data of tree 
	for (i=0, index=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			data[i][j] = M->elements[index++];
		}
	}

		
	while (cnt < minNum && cnt < rowDim) {	// {F1, F2, ...}
		for (index=0; index<rowDim && vis[index]; index++){};
		/* if (index >= rowDim) { 
				continue;
		} */
		tree = Q[index]; 	
		parent[tree] = -1;
		for (i=0; i<colDim; i++) {
			child[tree][i] = -1;
		}
		span[tree] = colDim - 1;
		tree_vis[tree] = 0;
		Link_add (&link, tree);
		vis[index] = 1;
		cnt++;
			
		while (index < rowDim) {
			for (index++; index<rowDim && vis[index]; index++){};
			if (index >=rowDim) {
				continue;
			}
			p = Q[index];
			parent[p] = -1;
			for (i=0; i<colDim; i++) {
				child[p][i] = -1;
			}
			span[p] = colDim - 1;
			tree_vis[p] = 0;

			// find position
			q = tree;
			while (q != -1) {
				for (m=1, greater=0, equal=0; m<colDim; m++) {
					d = data[q][m] - data[p][m];
					if (d < -DBL_EPSILON ) {
						greater++;
					} else if (d < DBL_EPSILON) {
						equal++;
					} else {
						break;
					}
				}
				if (m >= colDim) { 
					d = data[q][0] - data[p][0];
					if (d < -DBL_EPSILON ) {
						greater++;
					} else if (d < DBL_EPSILON) {
						equal++;
					}
					if (equal == colDim) {	// there are the same solution
						Link_add (&link, p);
						vis[index] = 1;
						cnt++;
					}

					p = -1;
					break; 
				} else {
					add_pos = q;
					add_no = m;
					tree_vis[q] = 1;
					span[q] = m;
					q = child[q][m];
				}
			}

			if (-1 == p) 
				continue;
			
			// check other solution
			q = add_pos;
			while (q != -1) {
				if (0 == tree_vis[q]) {
					tree_vis[q] = 1;
					for (m=1, greater=0, equal=0; m<colDim; m++) {
						d = data[q][m] - data[p][m];
						if (d < -DBL_EPSILON ) {
							greater++;
						} else if (d < DBL_EPSILON) {
							equal++;
						} else {
							break;
						}
					}
					if (m >= colDim) { 
						d = data[q][0] - data[p][0];
						if (d < -DBL_EPSILON ) {
							greater++;
						} else if (d < DBL_EPSILON) {
							equal++;
						}
						if (equal == colDim) {	// there are the same solution
							Link_add (&link, p);
							vis[index] = 1;
							cnt++;
						}

						p = -1;
						break; 
					} 
				}

				while (span[q] > 0 && -1 == child[q][span[q]]) {
					span[q]--;
				}
				if (span[q] > 0 && child[q][span[q]] != -1) { 
					q = child[q][span[q]];
				} else {
					tree_vis[q] = 0;
					span[q] = colDim - 1;
					q = parent[q];
					if (q != -1) 
						span[q]--;
				}
			}

			if (p != -1) {
				parent[p] = add_pos;
				child[add_pos][add_no] = p;
				Link_add (&link, p);				
				vis[index] = 1;
				cnt++;
			}
		}

		List_add (&list, link);
		Link_free (&link);
		link = NULL;
	}

	free (Q);
	return list;

}
