#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "avltree.h"

static void retrace_for_insert (avl_node_t *p);
static void retrace_for_delete (avl_node_t *p);

static void rotate_left (avl_node_t *p);
static void rotate_right (avl_node_t *p);
static void rotate_left_right (avl_node_t *p);
static void rotate_right_left (avl_node_t *p);

avl_node_t* avl_node_new () {
	size_t size = sizeof (avl_node_t);
	avl_node_t * p = (avl_node_t *)malloc (size);
	memset (p, 0, size);
	return p;
}

avl_node_t* avl_node_new (double key, void *value) {
	size_t size = sizeof (avl_node_t);
	avl_node_t * p = (avl_node_t *)malloc (size);
	memset (p, 0, size);
	p->key = key;
	p->value = value;
	return p;
}

void avl_node_free (avl_node_t **node) {
	if (*node == NULL) return;
	if ((*node)->value != NULL)
		free ((*node)->value);
	if (*node != NULL)	
		free (*node);
	*node = NULL;
}

void avl_tree_free (avl_tree_t *tree) {
	if (*tree == NULL) return;
	avl_node_t *p = NULL, *q = NULL;

	p = *tree;
	while (p) {
		if (p->leftchild != NULL) {
			p = p->leftchild;
		} else if (p->rightchild != NULL){
			p = p->rightchild;	
		} else {
			q = p->parent;
			if (q != NULL && q->leftchild == p)
				q->leftchild = NULL;
			else if (q != NULL)
				q->rightchild = NULL;
			avl_node_free (&p);
			p = q;
		}
	}
	*tree = NULL;
}

avl_node_t* avl_insert (avl_tree_t *tree, double key, void *value) {
	avl_node_t* node = NULL; 

	if (*tree == NULL) {
		node = avl_node_new (key, value);
		*tree = node;
		return NULL;
	}
	
	avl_node_t *p = *tree, *q = NULL;
	while (p) {
		if (p->key > key) {
			q = p;
			p = p->leftchild;
		} else {
			q = p;
			p = p->rightchild;
		}
	}

	node = avl_node_new (key, value);
	node->parent = q;
	if (q->key > key) 
		q->leftchild = node;			
	else 
		q->rightchild = node;

	retrace_for_insert (node); // retracing the added node

	while ((*tree)->parent)
		(*tree)=(*tree)->parent;

	return node;
}

void avl_delete (avl_tree_t *tree, avl_node_t *p) {
	avl_node_t *q = NULL, *t = NULL;
	if (p == NULL || *tree == NULL) 
		return;
	if ((*tree)->leftchild == NULL && (*tree)->rightchild == NULL) {
		avl_node_free (tree);
		return;
	}
	if (p->value)
		free (p->value);
	if (p->leftchild) {
		t = p;
		q = avl_pred (t);
		t->key = q->key;
		t->value = q->value;
		while (q->leftchild) {
			t = q;
			q = avl_pred (t);
			t->key = q->key;
			t->value = q->value;
		}
	} else if (p->rightchild) {
		t = p;
		q = avl_succ (t);
		t->key = q->key;
		t->value = q->value;
		while (q->rightchild) {
			t = q;
			q = avl_succ (t);
			t->key = q->key;
			t->value = q->value;
		}
	} else {
		q = p;
	} 

	t = q->parent;			// Now, t is the parent of q to be free
	if (t && t->leftchild == q) {
		t->factor++;
		t->leftchild = NULL;
	}
	else if (t) {
		t->factor--;
		t->rightchild = NULL;
	}
	q->parent = NULL;
	free (q);

	retrace_for_delete (t); 	// Now, q is the predecessor without child
}

void avl_delete (avl_tree_t *tree, double key) {
	avl_tree_t p = NULL;

	p = avl_search (*tree, key);
	while (p) {
		avl_delete (tree, p);
		p = avl_search (*tree, key);
	//	break;
	}
}

avl_node_t* avl_search (avl_tree_t tree, double key) {
	avl_node_t *p = tree;

	while (p) {
		if (p->key > key) {
			p = p->leftchild;
		} else if (p->key < key) {
			p = p->rightchild;
		} else {
			return p;
		}
	}
	return NULL;
}

avl_node_t* avl_succ (avl_node_t *p) {
	avl_node_t *q = NULL;
	if (p==NULL) return NULL;
	
	if (p->rightchild != NULL) {
		q = p->rightchild;
		while (q->leftchild)
			q = q->leftchild;
	} else {
		q = p->parent;
		while (q != NULL && q->key < p->key)
			q = q->parent;
	} 
	return q;
}

avl_node_t* avl_pred (avl_node_t *p) {
	avl_node_t *q = NULL;
	if (p==NULL) return NULL;
	
	if (p->leftchild != NULL) {
		q = p->leftchild;
		while (q->rightchild)
			q = q->rightchild;
	} else {
		q = p->parent;
		while (q != NULL && q->key > p->key)
			q = q->parent;
	} 
	return q;
}

static void retrace_for_insert (avl_node_t *p) {
	avl_node_t *q = p->parent;
	while (q) {
		if (q->leftchild == p)
			q->factor--;
		else 
			q->factor++;
		
		if (q->factor == 0) {
			break;
		} else 	if (q->factor == -2) {
			if (q->leftchild->factor == -1) 	
				rotate_right (q);
			else
				rotate_left_right (q);
			break;
		} else if (q->factor == 2) {
			if (q->rightchild->factor == 1)
				rotate_left (q);
			else
				rotate_right_left (q);
			break;
		} else {
			p = q;
			q = p->parent;
		}
	} 
}

static void retrace_for_delete (avl_node_t *p) {
	avl_node_t *q = NULL;

	if (p==NULL) return;

	if (p->factor == 1 || p->factor == -1) {
		return;
	} else 	if (p->factor == -2) {
		if (p->leftchild->factor == -1) 	
			rotate_right (p);
		else
			rotate_left_right (p);
	} else if (p->factor == 2) {
		if (p->rightchild->factor == 1)
			rotate_left (p);
		else
			rotate_right_left (p);
	} else {
		q = p->parent;
		while (q) {
			if (q->leftchild == p)
				q->factor++;
			else 
				q->factor--;

			if (q->factor == 1 || q->factor == -1) {
				return;	
			} else if (q->factor == -2) {
				if (q->leftchild->factor == -1) 	
					rotate_right (q);
				else
					rotate_left_right (q);
			} else if (q->factor == 2) {
				if (q->rightchild->factor == 1)
					rotate_left (q);
				else
					rotate_right_left (q);
			} else {
				p = q;
				q = p->parent;
			} 
		}
	}
}

static void rotate_left (avl_node_t *p) {
	avl_tree_t R = p->parent;
	avl_tree_t q = p->rightchild;

	// #1
	if (R && R->leftchild == p)
		R->leftchild = q;
	else if (R)
		R->rightchild = q;
	q->parent = R;

	// #2
	p->rightchild = q->leftchild;
	if (q->leftchild)
		q->leftchild->parent = p;

	// #3
	q->leftchild = p;
	p->parent = q;

	q->factor = 0;
	p->factor = 0;
		
}

static void rotate_right (avl_node_t *p) {
	avl_tree_t R = p->parent;
	avl_tree_t q = p->leftchild;

	// #1
	if (R && R->rightchild == p)
		R->rightchild = q;
	else if (R)
		R->leftchild = q;
	q->parent = R;

	// #2
	p->leftchild = q->rightchild;
	if (q->rightchild)
		q->rightchild->parent = p;

	// #3
	q->rightchild = p;
	p->parent = q;

	q->factor = 0;
	p->factor = 0;
}

static void rotate_left_right (avl_node_t *p) {
	avl_node_t *R = p->parent;
	avl_node_t *q = p->leftchild;
	avl_node_t *r = q->rightchild;

	if (r->factor==0) {
		p->factor = 0;
		q->factor = 0;
	} else if (r->factor == 1) {	
		p->factor = 0;
		q->factor = -1;
		r->factor = 0;
	} else {
		p->factor = 1;
		q->factor = 0;
		r->factor = 0;
	}

	// #1
	q->rightchild = r->leftchild;
	if (r->leftchild)
		r->leftchild->parent = q;
	// #2
	r->leftchild = q;
	q->parent = r;

	// #3
	p->leftchild = r->rightchild;
	if (r->rightchild)
		r->rightchild->parent = p;
	
	// #4
	r->rightchild = p;
	p->parent = r;

	// #5
	if (R && R->leftchild == p)
		R->leftchild = r;
	else if (R)
		R->rightchild = r;
	r->parent = R;
}

static void rotate_right_left (avl_node_t *p) {
	avl_node_t *R = p->parent;
	avl_node_t *q = p->rightchild;
	avl_node_t *r = q->leftchild;

	if (r->factor==0) {
		p->factor = 0;
		q->factor = 0;
	} else if (r->factor == 1) {	
		p->factor = -1;
		q->factor = 0;
		r->factor = 0;
	} else {
		p->factor = 0;
		q->factor = 1;
		r->factor = 0;
	}
	// #1
	q->leftchild = r->rightchild;
	if (r->rightchild)
		r->rightchild->parent = q;

	// #2
	r->rightchild = q;
	q->parent = r;

	// #3
	p->rightchild = r->leftchild;
	if (r->leftchild)
		r->leftchild->parent = p;

	// #4
	r->leftchild = p;
	p->parent = r;

	// #5
	if (R && R->leftchild == p)
		R->leftchild = r;
	else if (R)
		R->rightchild = r;
	r->parent = R;
}

avl_tree_t avl_dup (avl_tree_t T) {
	if (T==NULL) return NULL;
	avl_tree_t R = NULL;
	avl_node_t *p = NULL, *q = NULL;
	avl_node_t *pp = NULL, *qq = NULL;
	avl_node_t *stackP[10000], *stackQ[10000];
	int len=0;

	R = avl_node_new (T->key, NULL);
	R->factor = T->factor;
	p = T;
	q = R;
	stackP[len] = p;
	stackQ[len] = q;
	len++;
	p->vis = 0;
	while (len) {
		p = stackP[len-1];
		q = stackQ[len-1];
		if (p->vis==0) { 	// visite left child of p
			p->vis++;
			pp = p->leftchild;
			if (pp != NULL)	{
				qq = avl_node_new (pp->key, NULL);
				qq->factor = pp->factor;
				stackP[len] = pp;
				stackQ[len] = qq;
				len++;
				pp->vis = 0;
				qq->parent = q;
				q->leftchild = qq;
			}
		} else if (p->vis==1) {
			p->vis++;
			pp = p->rightchild;
			if (pp != NULL) {
				qq = avl_node_new (pp->key, NULL);
				qq->factor = pp->factor;
				stackP[len] = pp;
				stackQ[len] = qq;
				len++;
				pp->vis = 0;
				qq->parent = q;
				q->rightchild = qq;
			}
		} else {
			len--;
		}
	}

	return R;
}


#define HEIGHT 3 
#define WIDTH  2 
#define STEP 1.0e-3 
#define MAXGEN 100000
#define QUEUEL 10000
typedef struct point_tag {
	double x, y;
	int NO;
	avl_tree_t last;
	avl_tree_t next;
} point_t;

static void visualize_init (avl_tree_t Q);
static void visualize_latex (avl_tree_t Q, char* file);
static void visualize_evolute (avl_tree_t Q);
static void visualize_grow (avl_tree_t Q);
static void visualize_recover (avl_tree_t Q);

void avl_visualize (avl_tree_t T) {
	char 	filename[256];
	char*	pointer = (char *)filename;
	avl_tree_t Q = avl_dup (T);		

	visualize_grow (Q);
	visualize_init (Q);


	visualize_evolute (Q);
	visualize_recover (Q);

	sprintf (filename, "avltree%ld%p", time (NULL), (char *)pointer);
	visualize_latex (Q, filename);
}

static void visualize_init (avl_tree_t Q)
{
	avl_node_t *queue[QUEUEL];
	int head=0, tail=0;
	double xAxis[100];
	memset (xAxis, 0, 100*sizeof (double));
	avl_node_t *p = NULL, *q=NULL;
	int NO=0;

	if (Q==NULL) return;
	point_t *t = (point_t *)malloc (sizeof (point_t));
	memset (t, 0, sizeof (point_t));
	t->x = 0;
	t->y = 0;
	t->NO = NO++;
	Q->value = (void *)t;
	Q->vis = 0;
	Q->depth = 0;
	queue[tail++] = Q;

	while (tail - head) {
		p = queue[head];
		if (p->vis == 0) {
			p->vis++;
			if (p->leftchild) {
				q = p->leftchild; 
				q->depth = p->depth + 1;
				t = (point_t *)malloc (sizeof (point_t));
				t->x = xAxis[q->depth];
				if (t->x < ((point_t *)p->value)->x)
					t->x = ((point_t *)p->value)->x;
				xAxis[q->depth] = t->x + WIDTH;
				t->y = 0 - q->depth * HEIGHT;
				t->NO = NO++;
				q->value = (void *)t;
				q->vis = 0;
				queue[tail] = q;
				tail = (tail + 1) % QUEUEL;
				if ((tail + 1) % QUEUEL == head) {
					fprintf (stderr, "the queue is overflow\n");
					exit (-1);
				}
			}
		} else if (p->vis == 1) {
			p->vis++;
			if (p->rightchild) {
				q = p->rightchild;
				q->depth = p->depth + 1;
				t = (point_t *)malloc (sizeof (point_t));
				t->x = xAxis[q->depth];
				if (t->x < ((point_t *)p->value)->x + WIDTH)
					t->x = ((point_t *)p->value)->x + WIDTH;
				xAxis[q->depth] = t->x + WIDTH;
				t->y = 0 - q->depth * HEIGHT;
				t->NO = NO++;
				q->value = (void *)t;
				q->vis = 0;
				queue[tail] = q;
				tail = (tail + 1) % QUEUEL;
				if ((tail + 1) % QUEUEL == head) {
					fprintf (stderr, "the queue is overflow\n");
					exit (-1);
				}
			}
		} else {
			head = (head + 1) % QUEUEL;
			if (head != tail) {
				((point_t *)p->value)->next = queue[head];
				((point_t *)queue[head]->value)->last = p;
			}
		}
	}
}


static void visualize_latex (avl_tree_t Q, char* file) {
	char buff[1024];
	FILE *fp = NULL;
	
	sprintf (buff, "latex/%s.tex", file);
	fp =fopen (buff, "w");
	if (fp == NULL) {
		fprintf (stderr, "file %s cannot be opened\n",file);
		exit (-1);

	}

	fprintf (fp, "\\documentclass{article}\n\n");

	fprintf (fp, "\\usepackage{pdflscape}\n");
	fprintf (fp, "\\usepackage{geometry}\n");
	fprintf (fp, "\\usepackage{tikz}\n");
	fprintf (fp, "\\usepackage{xcolor}\n\n");
	fprintf (fp, "\\geometry{a4paper, scale=1.0, centering}\n");
	fprintf (fp, "\\usetikzlibrary{arrows,shapes,chains}\n\n");

	fprintf (fp, "\\begin{document}\n");
	fprintf (fp, "\\begin{landscape}\n");
	fprintf (fp, "\\begin{figure}\n");
	fprintf (fp, "\\centering\n");
	fprintf (fp, "\\scalebox{1.0}{");
	fprintf (fp, "\\begin{tikzpicture}[scale=1.0]\n");
	fprintf (fp, "\\tikzstyle{point}=[circle, draw=red]\n");
	
	if (Q==NULL) return;
	avl_node_t *p = NULL, *pp = NULL;
	avl_node_t *stack[1000];
	int len=0;
	point_t *t=NULL;

	t = (point_t *)(Q->value);
	fprintf (fp, "\\node[point] (P%d) at (%lf, %lf) {$%.2lf(%d)$};\n", t->NO, t->x, t->y, Q->key, Q->factor);
	stack[len] = Q;
	len++;
	Q->vis = 0;
	while (len) {
		p = stack[len-1];
		if (p->vis==0) { 	// visite left child of p
			p->vis++;
			pp = p->leftchild;
			if (pp != NULL)	{
				t = (point_t *)pp->value;
				fprintf (fp, "\\node[point] (P%d) at (%lf, %lf) {$%.2lf(%d)$};\n",
					t->NO, t->x, t->y, pp->key, pp->factor);
				fprintf (fp, "\\draw[green] (P%d) -- node[left] {L} (P%d);\n", 
					((point_t *)(p->value))->NO, t->NO);
				stack[len] = pp;
				len++;
				pp->vis = 0;
			}
		} else if (p->vis==1) {
			p->vis++;
			pp = p->rightchild;
			if (pp != NULL) {
				t = (point_t *)pp->value;
				fprintf (fp, "\\node[point] (P%d) at (%lf, %lf) {$%.2lf(%d)$};\n",
					t->NO, t->x, t->y, pp->key, pp->factor);
				fprintf (fp, "\\draw[blue] (P%d) -- node[right] {R} (P%d);\n", 
					((point_t *)(p->value))->NO, t->NO);
				stack[len] = pp;
				len++;
				pp->vis = 0;
			}
		} else {
			len--;
		}
	}

	fprintf (fp, "\\end{tikzpicture}}\n");
	fprintf (fp, "\\end{figure}\n");
	fprintf (fp, "\\end{landscape}\n");
	fprintf (fp, "\\end{document}\n");

	fclose (fp);

	sprintf (buff, "cd latex && xelatex %s.tex > /dev/null && open -a safari %s.pdf" , file, file);
	if (system (buff)) {
		fprintf (stderr, "The command [%s] cannot be excuted\n", buff);
		exit (0);
	}
}

static void evolute_align (avl_tree_t p);
static void evolute_justify (avl_tree_t p);
static void visualize_evolute (avl_tree_t Q) {
	if (Q==NULL) return;
	avl_node_t *p = NULL, *pp = NULL;
	avl_node_t *stack[1000];
	int len=0;
	int i=0;

	for (i=MAXGEN; i>=0; i--) {
		len = 0;
		stack[len] = Q;
		len++;
		Q->vis = 0;
		while (len) {
			p = stack[len-1];

			evolute_align (p);
			evolute_justify (p);

			if (p->vis==0) { 	// visite left child of p
				p->vis++;
				pp = p->leftchild;
				if (pp != NULL)	{
					stack[len] = pp;
					len++;
					pp->vis = 0;
				}
			} else if (p->vis==1) {
				p->vis++;
				pp = p->rightchild;
				if (pp != NULL) {
					stack[len] = pp;
					len++;
					pp->vis = 0;
				}
			} else {
				len--;
			}
		}
	}
}


static void evolute_align (avl_tree_t p) {
	point_t *t=NULL;
	double mid = 0;

	if (p->leftchild && p->rightchild) {
		mid = (((point_t *)p->leftchild->value)->x + ((point_t *)p->rightchild->value)->x ) / 2;
		t = (point_t *)p->value;
		if (t->x < mid) {
			t->x += STEP;
		} else if (t->x >mid) {
			t->x -= STEP;	
		}
	} else if (p->leftchild) {
		mid = ((point_t *)p->value)->x - HEIGHT * 0.618 ;
		t = (point_t *)p->leftchild->value;
		if (t->x < mid) {
			t->x += STEP;
		} else if (t->x > mid) {
			t->x -= STEP;	
		}
		
			
/*		mid = ((point_t *)p->leftchild->value)->x + HEIGHT * 0.618;
		t = (point_t *)p->value;
		if (t->x < mid) {
			t->x += STEP;
		} else if (t->x >mid) {
			t->x -= STEP;	
		}
*/	
	} else if (p->rightchild) {
		mid = ((point_t *)p->value)->x + HEIGHT * 0.618 ;
		t = (point_t *)p->rightchild->value;
		if (t->x < mid) {
			t->x += STEP;
		} else if (t->x > mid) {
			t->x -= STEP;	
		}
	
			
/*		mid = ((point_t *)p->rightchild->value)->x - HEIGHT * 0.618;
		t = (point_t *)p->value;
		if (t->x < mid) {
			t->x += STEP;
		} else if (t->x >mid) {
			t->x -= STEP;	
		}
*/	
	}
}

static void evolute_justify (avl_tree_t p) {
	point_t *t=NULL;
	double mid = 0;
	avl_tree_t last = NULL, next = NULL;

	avl_tree_t head = NULL, tail = NULL;
	double max = 1.0e+100, min = -1.0e+100;
	int num;

	last = ((point_t *)p->value)->last;
	next = ((point_t *)p->value)->next;
	if (last && next && last->depth == next->depth) {
		mid = (((point_t *)last->value)->x + 
				((point_t *)next->value)->x ) / 2;
		t = (point_t *)p->value;
		if (t->x < mid) {
			t->x += STEP;
		} else if (t->x >mid) {
			t->x -= STEP;	
		}
	}

	num = 1;
	tail = ((point_t *)p->value)->next;
	while (tail && tail->depth == p->depth) {
		num++;
		max = ((point_t *)tail->value)->x;
		tail = ((point_t *)tail->value)->next;
	}
	
	head = ((point_t *)p->value)->last;
	while (head && head->depth == p->depth) {
		num++;
		min = ((point_t *)head->value)->x;
		head = ((point_t *)head->value)->last;
	}

	if (num > 1 && (max - min) < (num -1) * WIDTH) {
		num = 0;
		tail = ((point_t *)p->value)->next;
		while (tail && tail->depth == p->depth) {
			num++;
			((point_t *)tail->value)->x += num * STEP;
			tail = ((point_t *)tail->value)->next;
		}
	
		num = 0;	
		head = ((point_t *)p->value)->last;
		while (head && head->depth == p->depth) {
			num++;
			((point_t *)head->value)->x -= num *STEP;
			head = ((point_t *)head->value)->last;
		}
	}
}


static void visualize_grow (avl_tree_t Q) {
	avl_node_t *p = NULL, *pp = NULL;
	avl_node_t *stack[1000];
	int len=0;
	int maxdepth = 0, maxNO;

	stack[0] = Q;
	len=1;
	Q->vis = 0;
	Q->depth = 0;
	Q->NO = 1;
	while (len) {
		p = stack[len-1];
		if (p->depth > maxdepth)
			maxdepth = p->depth;

		if (p->vis==0) { 	// visite left child of p
			p->vis++;
			pp = p->leftchild;
			if (pp != NULL)	{
				stack[len] = pp;
				len++;
				pp->vis = 0;
				pp->depth = p->depth + 1;
				pp->NO = p->NO * 2;
			}
		} else if (p->vis==1) {
			p->vis++;
			pp = p->rightchild;
			if (pp != NULL) {
				stack[len] = pp;
				len++;
				pp->vis = 0;
				pp->depth = p->depth + 1;
				pp->NO = p->NO * 2 + 1;
			}
		} else {
			len--;
		}
	}

	maxNO = (1<<(maxdepth+1));
	stack[0] = Q;
	len=1;
	Q->vis = 0;
	Q->depth = 0;
	Q->NO = 1;
	while (len) {
		p = stack[len-1];
		if (p->vis==0) { 	// visite left child of p
			p->vis++;
			pp = p->leftchild;
			if (pp != NULL)	{
				stack[len] = pp;
				len++;
				pp->vis = 0;
			} else if (p->NO * 2 < maxNO) {
				pp = avl_node_new ();
				pp->invalid = 1;
				p->leftchild = pp;
				pp->parent = p;

				stack[len] = pp;
				len++;
				pp->vis = 0;
				pp->depth = p->depth + 1;
				pp->NO = p->NO*2;
			}
		} else if (p->vis==1) {
			p->vis++;
			pp = p->rightchild;
			if (pp != NULL) {
				stack[len] = pp;
				len++;
				pp->vis = 0;
			} else if (p->NO*2 + 1 < maxNO){
				pp = avl_node_new ();
				pp->invalid = 1;
				p->rightchild = pp;
				pp->parent = p;

				stack[len] = pp;
				len++;
				pp->vis = 0;
				pp->depth = p->depth + 1;
				pp->NO = p->NO*2 + 1;
			}
		} else {
			len--;
		}
	}
}

static void visualize_recover (avl_tree_t Q) {
	avl_node_t *p = NULL, *pp = NULL;
	avl_node_t *stack[1000];
	int len=0;

	stack[0] = Q;
	len=1;
	Q->vis = 0;
	while (len) {
		p = stack[len-1];
		if (p->invalid == 1) {
			pp = p->parent;

			if (pp && pp->leftchild == p)
				pp->leftchild = NULL;
			else if (pp)
				pp->rightchild = NULL;
			p->parent = NULL;
			avl_tree_free (&p);
			len--;
			continue;
		}

		if (p->vis==0) { 	// visite left child of p
			p->vis++;
			pp = p->leftchild;
			if (pp != NULL)	{
				stack[len] = pp;
				len++;
				pp->vis = 0;
			}
		} else if (p->vis==1) {
			p->vis++;
			pp = p->rightchild;
			if (pp != NULL) {
				stack[len] = pp;
				len++;
				pp->vis = 0;
			}
		} else {
			len--;
		}
	}

}
