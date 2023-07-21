#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>

#include "matrix.h"
#include "dominate.h"
#include "link.h"
#include "singular_value_decomposition.h"
#include "algebra.h"
#include "hv.h"
#include "crowding.h"

static void strstrip (char *str) {
	int i=0;
	while (str[i] != '\0') {
		if (str[i] == '\r' || str[i] == '\n') {
			str[i] = '\0';
			break;
		}
		i++;
	}
}

Matrix_t* Matrix_read (char *filename) {
	Matrix_t *M = (Matrix_t *)malloc (sizeof (Matrix_t));
	FILE 	*fp = NULL;
	char 	line[1024], *tok;
	int  	row;
	int	index; 

	fp = fopen(filename,"r");
	if (fp == NULL) {
		fprintf (stderr,"File %s does not open\n", filename);
		exit (EXIT_FAILURE);
	}

	M->rowDim = M->colDim = 0;
	M->elements = NULL;

	row = 0;
	index = 0;
	while (fgets(line, sizeof(line), fp) != NULL) {
		strstrip (line);
		tok = strtok (line, " ,\t\n");
		M->elements = (double *)realloc ((void *)M->elements, (index + 1)*sizeof (double));	
		M->elements[index++] = atof (tok);

		while ((tok = strtok (NULL, " ,\t\n")) != NULL) {
			M->elements = (double *)realloc ((void *)M->elements, (index + 1)*sizeof (double));	
			M->elements[index++] = atof (tok);
		}
		row++;	
	}

	M->rowDim = row;
	if (row > 0)	
		M->colDim = index / row;
	else 
		M->colDim = 0;
	return M;
}

Matrix_t *Matrix_new () {
	Matrix_t *M = (Matrix_t *)calloc (1, sizeof (Matrix_t));
	return M;
}

Matrix_t *Matrix_new (int nrow) {
	return Matrix_new (nrow, nrow);
}
Matrix_t *Matrix_new (int nrow, int ncol) {
	Matrix_t *M = (Matrix_t *)malloc (sizeof (Matrix_t));

	M->rowDim = nrow;
	M->colDim = ncol;
	M->elements = (double *)calloc (nrow * ncol, sizeof (double));
	return M;
}

void Matrix_print (Matrix_t *M) {
	if (M == NULL || M->elements == NULL) {
		printf ("The matrix is empty\n");
		return;
	}
	Matrix_print (M, stdout);
}

void Matrix_print (Matrix_t *M, char *fn) {
	FILE *fp = NULL;

	if (M == NULL || M->elements == NULL) {
		printf ("The matrix is empty\n");
		return;
	}

	fp = fopen (fn, "w");		
	if (fp == NULL) {
		fprintf (stderr, "file %s does not open\n", fn);
		return;
	}
	Matrix_print (M, fp);
	fclose (fp);
}

void Matrix_print (Matrix_t *M, FILE *fp) {
	int i, j, index;

	if (M == NULL || M->elements == NULL) {
		printf ("The matrix is empty\n");
		return;
	}

       	for (i=0, index=0; i < M->rowDim; i++) {
      		for (j=0; j<M->colDim; j++) {
                       		fprintf (fp, "%.16f ", M->elements[index++]);
	              }
               fprintf (fp, "\n");
       	}
}

void Matrix_latex (Matrix_t *M) {
	char buf[256];
	char filename[1024];
	FILE *fp = NULL;
	int i, j;
	
	sprintf (buf, "m%ld%p", time (NULL), (char *)filename);
	sprintf (filename, "latex/%s.tex", buf);
	if ((fp = fopen (filename, "w")) == NULL) {
		fprintf (stderr, "Cannot open file %s\n", filename);	
		exit (-1);
	}
		
	fprintf (fp, "\\documentclass{article}\n\n");
	fprintf (fp, "\\usepackage{booktabs}\n");
	fprintf (fp, "\\usepackage{diagbox}\n\n");
	fprintf (fp, "\\begin{document}\n");
	fprintf (fp, "\\begin{figure}\n");
	fprintf (fp, "\\centering\n");
	fprintf (fp, "\\begin{tabular}{c|*{%d}{c}}\n", M->colDim);
	fprintf (fp, "\\toprule\n");
	fprintf (fp, "\\diagbox{R}{C}");
	for (i=0; i<M->colDim; i++) {
		fprintf (fp, " & C%02d", i+1);
	}	
	fprintf (fp, "\\\\\n");
	fprintf (fp, "\\midrule\n");
	for (i=0; i<M->rowDim; i++) {
		fprintf (fp, "A%02d", i+1);
		for (j=0; j<M->colDim; j++) {
			fprintf (fp, " & %lf", M->elements[i*M->colDim + j]);
		}
		fprintf (fp, "\\\\\n");
	}
	fprintf (fp, "\\bottomrule");
	fprintf (fp, "\\end{tabular}\n");
	fprintf (fp, "\\end{figure}\n");
	fprintf (fp, "\\end{document}\n");
	fclose (fp);
	
	sprintf (filename, "cd latex && xelatex %s.tex > /dev/null && open -a safari %s.pdf", buf, buf);
	if (system (filename)) {
		fprintf (stderr, "The command [%s] cannot be excuted\n", filename);
		exit (0);
	}
}

void Matrix_octave (Matrix_t *M) {
	char buf[256];
	char filename[1024];
	FILE *fp = NULL;
	
	sprintf (buf, "m%ld%p", time (NULL), (char *)filename);
	sprintf (filename, "matlab/%s.m", buf);
	if ((fp = fopen (filename, "w")) == NULL) {
		fprintf (stderr, "Cannot open file %s\n", filename);	
		exit (-1);
	}
	fprintf (fp, "M=[\n");	
	Matrix_print (M, fp);
	fprintf (fp, "];\n");
	if (M->colDim == 2)
		fprintf (fp, "scatter (M(:,1), M(:,2));\n");
	else if (M->colDim == 3)
		fprintf (fp, "scatter3 (M(:,1), M(:,2), M(:,3));\n");
	else {
		fprintf (fp, "X=[1:%d]';\n", M->colDim);
		fprintf (fp, "X=repmat (X,[1, %d]);\n", M->rowDim);
		fprintf (fp, "plot (X, M');\n");
	}
	fprintf (fp, "print ('%s', '-depsc');\n", buf);

	fclose (fp);
	sprintf (filename, "cd matlab && octave %s.m > /dev/null && open -a safari %s.eps", buf, buf);
	if (system (filename)) {
		fprintf (stderr, "The command [%s] cannot be excuted\n", filename);
		exit (0);
	}
}

Matrix_t* Matrix_sub (Matrix_t *M, int *subset, int n) {
	int 		i;
	Matrix_t*	subM = (Matrix_t *)malloc (sizeof (Matrix_t));

	subM->rowDim = n;
	subM->colDim = M->colDim;
	size_t size = M->colDim * n * sizeof (double);
	subM->elements = (double *)malloc (size);

	for (i=0; i<n; i++) {
		size = M->colDim * sizeof (double);
		memcpy (subM->elements + i * M->colDim, M->elements + subset[i] * M->colDim, size);
	}	
	return subM;
}

Matrix_t* Matrix_sub (Matrix_t *M, Link_t *link) {
	int *arr = Link2Array (link);
	Matrix_t *subM = Matrix_sub (M, arr+1, arr[0]);

	free (arr);
	return subM;
}

void Matrix_free (Matrix_t **M) {
	if ((*M) != NULL && (*M)->elements != NULL) { 
		free ((*M)->elements);
	} 
	if ((*M) != NULL) {
		free (*M);
		*(M) = NULL;
	}
}


Matrix_t* Matrix_min (Matrix_t *M) {
	int i, j;
	Matrix_t *p = (Matrix_t *)malloc (sizeof (Matrix_t));
	p->colDim = M->colDim;
	p->rowDim = 1;
	size_t size = p->colDim * sizeof (double);
	p->elements = (double *) malloc (size);

	memcpy (p->elements, M->elements, size);
	for (i=1; i<M->rowDim; i++) {
		for (j=0; j<M->colDim; j++) {
			if (M->elements[i * M->colDim + j] < p->elements[j] )
				p->elements[j] = M->elements[i*M->colDim + j];
		}
	}

	return p;
}

Matrix_t* Matrix_max (Matrix_t *M) {
	int i, j;
	Matrix_t *p = (Matrix_t *)malloc (sizeof (Matrix_t));
	p->colDim = M->colDim;
	p->rowDim = 1;
	size_t size = p->colDim * sizeof (double);
	p->elements = (double *) malloc (size);

	memcpy (p->elements, M->elements, size);
	for (i=1; i < M->rowDim; i++) {
		for (j=0; j < M->colDim; j++) {
			if (M->elements[i * M->colDim + j] > p->elements[j] )
				p->elements[j] = M->elements[i*M->colDim + j];
		}
	}

	return p;
}

Matrix_t* Matrix_bound (Matrix_t *M) {
	int i, j;
	Matrix_t *bound = (Matrix_t *)malloc (sizeof (Matrix_t));
	bound->colDim = M->colDim;
	bound->rowDim = 2;
	size_t size = 2 * bound->colDim * sizeof (double);
	bound->elements = (double *) malloc (size);

	memcpy (bound->elements, M->elements, size);
	for (i=0; i < M->rowDim; i++) {
		for (j=0; j < M->colDim; j++) {
			if (M->elements[i * M->colDim + j] < bound->elements[j] )
				bound->elements[j] = M->elements[i*M->colDim + j];
			if (M->elements[i * M->colDim + j] > bound->elements[bound->colDim + j])
				bound->elements[bound->colDim + j] = M->elements[i*M->colDim + j];
		}
	}
	return bound;
}


void Matrix_cat (Matrix_t **M1, Matrix_t *M2) {
	if (M2==NULL || M2->rowDim == 0) return;
	if ((*M1) == NULL) *(M1) = Matrix_new();
	if ((*M1)->elements == M2->elements) return;

	int rowDim = (*M1)->rowDim + M2->rowDim;
	int colDim = (*M1)->colDim;

	size_t size = rowDim *colDim * sizeof (double);
	(*M1)->elements = (double *)realloc ((*M1)->elements, size);

	size = M2->rowDim * M2->colDim * sizeof (double);
	memcpy ((*M1)->elements + (*M1)->rowDim * (*M1)->colDim, M2->elements, size);
	
	(*M1)->rowDim = rowDim;
	(*M1)->colDim = colDim;
}

Matrix_t* Matrix_norm (Matrix_t *M){
	Matrix_t *bound = NULL;
	Matrix_t *normM = NULL;

	bound = Matrix_bound (M);

	normM = Matrix_norm (M, bound);
	Matrix_free (&bound);

	return normM;
}

Matrix_t *Matrix_norm (Matrix_t *M, Matrix_t *bound) {
	int 		i, j, index;
	double 		d;
	Matrix_t*  	normM = NULL;
	size_t 		size;

	if (M == NULL || M->elements == NULL || M->rowDim < 1 || M->colDim < 1) {
		return NULL;
	}

 	normM = (Matrix_t *)malloc (sizeof (Matrix_t));
	normM->rowDim = M->rowDim;	
	normM->colDim = M->colDim;

	size = M->rowDim * M->colDim * sizeof (double);
	normM->elements = (double *)malloc (size);

	for (i=0, index=0; i<M->rowDim; i++) {
		for (j=0; j<M->colDim; j++) {
			d = bound->elements[bound->colDim + j] - bound->elements[j];
			if (d < 1.0e-6) {
				normM->elements[index] = 0;
				index++;
			} else {
				normM->elements[index]	=  (M->elements[index] - bound->elements[j] ) / d;
				index++;
			}
		}
	}

	return normM;
}

Matrix_t *Matrix_norm (Matrix_t *M, Matrix_t *maxmum, Matrix_t *minmum) {
	int 		i, j, index;
	double 		d;
	Matrix_t*  	normM = NULL;  
	size_t 		size;

	if (M == NULL || M->elements == NULL || M->rowDim < 1 || M->colDim < 1) {
		return NULL;
	}

	normM = (Matrix_t *)malloc (sizeof (Matrix_t));
	normM->rowDim = M->rowDim;	
	normM->colDim = M->colDim;

	size = M->rowDim * M->colDim * sizeof (double);
	normM->elements = (double *)malloc (size);

	for (i=0, index=0; i<M->rowDim; i++) {
		for (j=0; j<M->colDim; j++) {
			d = maxmum->elements[j] - minmum->elements[j];
			if (d < 1.0e-6) {
				normM->elements[index] = 0;
				index++;
			} else {
				normM->elements[index]	=  (M->elements[index] - minmum->elements[j] ) / d;
				index++;
			}
		}
	}

	return normM;
}

Matrix_t *Matrix_norm (Matrix_t *M, double *maxmum, double *minmum) {
	int 		i, j, index;
	double 		d;
	Matrix_t*  	normM = NULL;  
	size_t 		size;

	if (M == NULL || M->elements == NULL || M->rowDim < 1 || M->colDim < 1) {
		return NULL;
	}

	normM = (Matrix_t *)malloc (sizeof (Matrix_t));
	normM->rowDim = M->rowDim;	
	normM->colDim = M->colDim;

	size = M->rowDim * M->colDim * sizeof (double);
	normM->elements = (double *)malloc (size);

	for (i=0, index=0; i<M->rowDim; i++) {
		for (j=0; j<M->colDim; j++) {
			d = maxmum[j] - minmum[j];
			if (d < 1.0e-6) {
				normM->elements[index] = 0;
				index++;
			} else {
				normM->elements[index]	=  (M->elements[index] - minmum[j]) / d;
				index++;
			}
		}
	}

	return normM;
}

Matrix_t *Matrix_dup (Matrix_t *M) {
	Matrix_t*	dup = NULL;
	size_t 		size;

	if (M == NULL || M->elements == NULL || M->rowDim < 1 || M->colDim < 1) {
		return NULL;
	}

	dup = (Matrix_t *)malloc (sizeof (Matrix_t));
	dup->rowDim = M->rowDim;
	dup->colDim = M->colDim;

	size = dup->rowDim * dup->colDim * sizeof (double);
	dup->elements = (double *)malloc (size);
	memcpy (dup->elements, M->elements, size);

	return dup;
}


Vector_t *Vector_new () {
	Vector_t *v = NULL;
	if ((v = (Vector_t *)malloc (sizeof (Vector_t))) == NULL) {
		fprintf (stderr, "%s:%d Cannot allocate memary\n", __FILE__, __LINE__);
		exit (-1);
	}
	memset (v, 0, sizeof (Vector_t));

	return v;
}

Vector_t* Matrix_sub (Matrix_t *M, int nrow) {
	Vector_t *v = NULL;
	if (!(nrow < M->rowDim)) return NULL;

	if ((v = (Vector_t *)malloc (sizeof (Vector_t))) == NULL) {
		fprintf (stderr, "%s:%d Cannot allocate memary\n", __FILE__, __LINE__);
		exit (-1);
	}

	v->dim = M->colDim;
	v->elements = (double *)malloc (v->dim * sizeof (double));
	memcpy (v->elements, M->elements + nrow * v->dim, v->dim * sizeof (double));

	return v;
}

void Vector_free (Vector_t **v) {
	if (v == NULL || *v == NULL) return;
	if ((*v)->elements != NULL)
		free ((*v)->elements);
	free (*v);
	*v = NULL;
}

Matrix_t* Matrix_front (Matrix_t *M) {
	Matrix_t*	front = NULL;
	Node_t*		p = NULL; 			
	List_t*		list = NULL;
	int		index[M->rowDim], len=0;	

	if (M->rowDim <= 1000) {
		return Matrix_front_1k (M);
	}

	list = ndSort (M);
	p = list->list_head->link_head;
	while (p != NULL) {
		index[len++] = p->data;
		p=p->next;
	}
	front = Matrix_sub (M, index, len);
	List_free (&list);
	
	return front;
}

Matrix_t* Matrix_front_1k (Matrix_t *M) {
	int*	Q=NULL;
	int	index = 0;

	//
	int  	tree, p, q;
	int  	add_pos = -1;
	int 	add_no = 0; 
	int 	i, j, m;
	double 	d;
	int 	rowDim = M->rowDim;
	int 	colDim = M->colDim;

	// the tree
	double	data[1010][20];
	int	parent[1010];
	int 	child[1010][20];
	int	span[1010];
	int 	tree_vis[1010];

	// queue
	int	queue[1010];
	int	tail = 0;
	
	// sort
	Q = sort (M);

	// initize data of tree 
	for (i=0, index=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			data[i][j] = M->elements[index++];
		}
	}

	// 
	index=0;
	tree = Q[index]; 	
	parent[tree] = -1;
	for (i=0; i<colDim; i++) {
		child[tree][i] = -1;
	}
	span[tree] = colDim - 1;
	tree_vis[tree] = 0;
	queue[tail++] = tree;
			
	for (index=1; index < rowDim; index++) {
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
			for (m=1; m<colDim; m++) {
				d = data[q][m] - data[p][m];
				if (d > DBL_EPSILON) {
					break;
				}
			}
			if (m >= colDim) { 
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
				for (m=1; m<colDim; m++) {
					d = data[q][m] - data[p][m];
					if (d > DBL_EPSILON ) {
						break;
					}
				}
				if (m >= colDim) { 
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
			queue[tail++] = p;
		}
	}


	free (Q);
	return Matrix_sub (M, queue, tail);
}


Matrix_t* Matrix_trim (Matrix_t *M) {
	Matrix_t *trim = Matrix_new ();
	int	i;
	size_t 	size;

	trim->rowDim = M->rowDim;
	trim->colDim = M->colDim-1;
	size = trim->rowDim * trim->colDim * sizeof (double);
	trim->elements = (double *)malloc (size);

	size = trim->colDim * sizeof (double);
	for (i=0; i<trim->rowDim; i++)
		memcpy (trim->elements + i * trim->colDim, M->elements + i * M->colDim + 1, size);

	return trim;
}

Matrix_t* Matrix_extreme (Matrix_t *M) {
	int arr[M->colDim];
	int i, j;
	double t;

	for (j=0; j<M->colDim; j++) {
		arr[j] = 0;
		t = M->elements[j];
		for (i=1; i<M->rowDim; i++) {
			if (M->elements[i*M->colDim+j] > t) {
				t = M->elements[i*M->colDim+j];
				arr[j] = i;
			}
		}
	}

	return Matrix_sub (M, arr, M->colDim);
}


Matrix_t* Matrix_inverse (Matrix_t *M) {
	Matrix_t *I = Matrix_dup (M);
	
	if (M->rowDim != M->colDim) {
		fprintf (stderr, "The numbers of columes(%d) and rows(%d) are not equavalent\n", M->colDim, M->rowDim);
		exit (-1);
	}

	Inv (M->elements, M->rowDim, I->elements);
	
	return I;
}

Matrix_t* Matrix_trans (Matrix_t *M) {
	Matrix_t *T = Matrix_new (M->colDim, M->rowDim);
	int i, j;

	for (i=0; i<T->rowDim; i++) {
		for (j=0; j<T->colDim; j++) {
			T->elements[i*T->colDim+j] = M->elements[j*M->colDim+i];
		}
	}

	return T;
}

Matrix_t* Matrix_compress (Matrix_t *M) {
	Matrix_t*	subM = NULL;
	Link_t*		link = NULL;
	int 		i, j, k;
	Node_t*		p=NULL;
	double 		d;	

	Link_add (&link, 0);
	for (i=1; i<M->rowDim; i++) {
		p = link->link_head;
		while (p) {
			j = p->data;

			// 	
			for (k=0; k<M->colDim; k++) {
				d = M->elements[i*M->colDim+k] - M->elements[j*M->colDim+k];
				if (d < 0)
					d = -d;
				if (d > 1.0e-6)
					break;
			}
			// if (distance_p2p(M->elements+i*M->colDim, M->elements+j*M->colDim, M->colDim)<1.0e-6){
			if (k >= M->colDim) {
				break;
			} else {
				p = p->next;
			}
		}
		if (p==NULL) {
			Link_add (&link, i);
		}
	}
	
	subM = Matrix_sub (M, link);
	Link_free (&link);

	return subM;
}

double* Matrix_mean (Matrix_t *M) {
	double	*A = NULL; 
	int 	i, j;

	if (M == NULL || M->elements == NULL || M->rowDim == 0)
		return NULL;

	A = (double *)calloc (M->colDim, sizeof(double));
	for (j=0; j<M->colDim; j++) {
		for (i=0; i<M->rowDim; i++) {
			A[j] += M->elements[i*M->colDim+j];
		}
		A[j] /= M->rowDim;
	}

	return A;
}

double* Matrix_var (Matrix_t *M) {
	double *A = NULL;
	double *D = NULL;	
	int i, j;

	if (M == NULL || M->elements == NULL || M->rowDim == 0)
		return NULL;
	A = Matrix_mean (M);

	D = (double *)calloc (M->colDim, sizeof(double));
	for (j=0; j<M->colDim; j++) {
		for (i=0; i<M->rowDim; i++) {
			D[j] += (M->elements[i*M->colDim+j] - A[j])*(M->elements[i*M->colDim+j] - A[j]);
		}
		D[j] /= M->rowDim;
	}

	free (A);
	return D;
}

double* Matrix_std (Matrix_t *M) {
	double*	S = NULL;
	int 	i;

	S = Matrix_var (M);
	for (i=0; i<M->colDim; i++) {
		S[i] = sqrt (S[i]);
	}
	return S;
}

double Matrix_angle (double *p, Matrix_t* M) {// angle between the p and the M
	double 	xValue, t;
	int 	i;

	if (M==NULL || M->elements == NULL || M->rowDim == 0)
		return 0;

	for (i=0, xValue=1.0e+100; i<M->rowDim; i++) {
		t = vector_angle (p, M->elements+i*M->colDim, M->colDim);
		if (t < xValue)
			xValue = t;
	}

	return xValue;
}

Matrix_t* Matrix_limited (Matrix_t *M, double *p) {     // to use a point to constraint M, return a PF. 
	int 		i, j;
	int 		rowDim, colDim;
	Matrix_t* 	S = NULL; 
	Matrix_t* 	F1 = NULL;

	if (M == NULL || M->elements == NULL || M->rowDim < 1 || M->colDim < 1) {
		return NULL;
	}

	rowDim = M->rowDim;
	colDim = M->colDim;
	S = Matrix_dup (M);

	for (i=0; i<rowDim; i++) {
		for (j=0; j<colDim; j++) {
			if (S->elements[i*colDim+j] < p[j]) {
				S->elements[i*colDim+j] = p[j];
			}
		}
	}

	F1 = Matrix_front (S);
	Matrix_free (&S);
	return F1;
}

int* Matrix_rank_by_hypervolume (Matrix_t * M) {
	int		rowDim = M->rowDim;
	int 		colDim = M->colDim;
	int*		index  = NULL;
	List_t*		list   = ndSort (M);
	Link_t*		link   = NULL;
	Matrix_t*	I = Matrix_new (rowDim, 2);	// [PF, HV (or angle)]
	Matrix_t*	N = NULL;
	Matrix_t*	T = NULL;
	Matrix_t*	F = NULL;
	int		i, j, k, a, b;
	double 		y[colDim+10], exchv = 1, angle;
	int*		array = NULL;

	// 1. set front number
	i=1;
	link = list->list_head;
	while (link) {
		array = Link2Array (link);
		for (j=array[0]; j>0; j--) {
			I->elements[array[j]*2+0] = i;
		}
		free (array);
		i++;
		link=link->next;
	}
	List_free (&list);

	// 2. find the minimum
	for (j=0; j<colDim; j++) {
		k = 0;
		memcpy (y, M->elements, colDim*sizeof (double));
		for (i=1; i<rowDim; i++) {
			for (a=0; a<colDim; a++) {
				b = (j+a) % colDim;
				if (M->elements[i*colDim+b] < y[b]) {
					k = i;	
					memcpy (y, M->elements+i*colDim, colDim*sizeof (double));
					break;
				} else if (fabs (M->elements[i*colDim+b] - y[b]) < DBL_EPSILON) {
					continue;
				} else {
					break;
				}
			}
		}
		I->elements[k*2+0] = -1;
	}

	// 3. compute exchv
	N = Matrix_norm (M);
	for (i=0; i<rowDim; i++) {
		memcpy (y, N->elements+i*colDim, colDim*sizeof(double));
		T = Matrix_dup (N);	
		for (j=0; j<colDim; j++) {
			T->elements[i*colDim+j] = T->elements[(rowDim-1)*colDim+j];
		}
		T->rowDim -= 1;
		for (k=0; k<T->rowDim; k++) {
			for (j=0; j<colDim; j++) {
				if (T->elements[k*colDim+j] < y[j]) {
					T->elements[k*colDim+j] = y[j];
				}
			}
		}
		F = Matrix_front (T);
		for (j=0, exchv=1; j<colDim; j++) {
			exchv *= (1.1 - y[j]);
		}
		exchv -= hv (F);
		I->elements[i*2+1] = -exchv;
		Matrix_free (&F);
		Matrix_free (&T);
	}

	// 4. compute angle
	for (i=0; i<rowDim; i++) if (I->elements[i*2+0] > 1){
		I->elements[i*2+1] = 1.0e+100;
		for (j=0; j<rowDim; j++) if (I->elements[j*2+0] <= I->elements[i*2+0] ){
			angle = vector_angle (N->elements+i*colDim, N->elements+j*colDim, colDim);
			if (angle < I->elements[i*2+1]) {
				I->elements[i*2+1] = -angle;
			}
		}
	}
	Matrix_free (&N);
	
	// 5. rank
	index = sort (I);
	Matrix_free (&I);

	return index;
}

int* Matrix_rank_by_crowdingdistance (Matrix_t * M) {
	int		rowDim = M->rowDim;
	int 		colDim = M->colDim;
	int*		index  = NULL;
	List_t*		list   = ndSort (M);
	Link_t*		link   = NULL;
	Matrix_t*	I = Matrix_new (rowDim, 2);	// [PF, crowding distance (or angle)]
	Matrix_t*	F = Matrix_sub (M, list->list_head);
	int		i, j, k, a, b;
	double 		y[colDim+10], angle;
	int*		array = NULL;
	double*		crowding = NULL;

	// 1. set front number
	i=1;
	link = list->list_head;
	while (link) {
		array = Link2Array (link);
		for (j=array[0]; j>0; j--) {
			I->elements[array[j]*2+0] = i;
		}
		free (array);
		i++;
		link=link->next;
	}

	// 2. find the minimum
	for (j=0; j<colDim; j++) {
		k = 0;
		memcpy (y, M->elements, colDim*sizeof (double));
		for (i=1; i<rowDim; i++) {
			for (a=0; a<colDim; a++) {
				b = (j+a) % colDim;
				if (M->elements[i*colDim+b] < y[b]) {
					k = i;	
					memcpy (y, M->elements+i*colDim, colDim*sizeof (double));
					break;
				} else if (fabs (M->elements[i*colDim+b] - y[b]) < DBL_EPSILON) {
					continue;
				} else {
					break;
				}
			}
		}
		I->elements[k*2+0] = -1;
	}

	// 3. compute crowding distance 
	crowding = crowding_distance (F);
	array = Link2Array (list->list_head);
	for (i=1; i<=array[0]; i++) {
		j = array[i];
		I->elements[j*2+1] = -crowding[i-1];
	}
	free (array);
	free (crowding);
	Matrix_free (&F);
	List_free (&list);

	// 4. compute angle
	for (i=0; i<rowDim; i++) if (I->elements[i*2+0] > 1){
		I->elements[i*2+1] = 1.0e+100;
		for (j=0; j<rowDim; j++) if (I->elements[j*2+0] <= I->elements[i*2+0] ){
			angle = vector_angle (M->elements+i*colDim, M->elements+j*colDim, colDim);
			if (angle < I->elements[i*2+1]) {
				I->elements[i*2+1] = -angle;
			}
		}
	}
	
	// 5. rank
	index = sort (I);
	Matrix_free (&I);

	return index;
}

int* Matrix_rank_by_sum (Matrix_t * M) {
	int		rowDim = M->rowDim;
	int 		colDim = M->colDim;
	int*		index  = NULL;
	Matrix_t*	I = Matrix_new (rowDim, 2);	// [sum, No. ]
	int		i;

	for (i=0; i<rowDim; i++) {
		I->elements[i*2+0] = sum (M->elements+i*colDim, colDim);
		I->elements[i*2+1] = i+0.1;
	}

	// rank
	index = sort (I);
	Matrix_free (&I);

	return index;
}
