#include "hv.h"
#include "link.h"

typedef struct hso_node_tag {
	Matrix_t	*M;
	int		*index;	// order of M sorted by first objective
	int 		cursor;	// the current cursor pointed to processing point
	double 		depth; 	// 

	struct 	hso_node_tag *next;
} hso_node_t;

typedef struct hso_stack_tag {
	hso_node_t *stack_top;
} hso_stack_t;

hso_node_t *hso_node_new ();

hso_node_t *hso_stack_top (hso_stack_t *stack);
void hso_stack_pop (hso_stack_t *stack); 
void hso_stack_push (hso_stack_t **stack, hso_node_t *node);

double hv_hso (Matrix_t *M)
{
	double 	hv = 0;
	hso_node_t 	*node = NULL, *nnode = NULL;
	hso_stack_t 	*stack = NULL;
	Matrix_t * subM = NULL, *trimM = NULL;

	//
	node = hso_node_new ();
	node->M = Matrix_dup (M);
	node->index = sort (node->M);
	node->cursor = 0;
	node->depth = 1.0;
	hso_stack_push (&stack, node);

	while ((node = hso_stack_top (stack)) != NULL) {
		if (node->M->colDim > 1 && node->cursor < node->M->rowDim) {  // push 
			nnode = hso_node_new ();			

			subM = Matrix_sub (node->M, node->index, node->cursor + 1);
			trimM = Matrix_trim (subM);
			Matrix_free (&subM);
			nnode->M = Matrix_front (trimM);
			Matrix_free (&trimM);

			nnode->index = sort (nnode->M);
			nnode->cursor = 0;
			if (node->cursor < node->M->rowDim - 1)
				nnode->depth = (node->depth) * 
					((node->M->elements[node->index[(node->cursor+1)]*node->M->colDim]) - 
							(node->M->elements[node->index[node->cursor]*node->M->colDim]));
			else
				nnode->depth = (node->depth) * 
					(1.1 - node->M->elements[node->index[node->cursor]*node->M->colDim]);
			hso_stack_push (&stack, nnode);
			node->cursor++;
		} else { // pop
			if (node->M->colDim == 1) {
				hv += (node->depth) * (1.1 - node->M->elements[node->index[0]]);
			}
			hso_stack_pop (stack);
		}
	}

	free (stack);
	return hv;
}

hso_node_t *hso_node_new () {
	size_t size = sizeof (hso_node_t);
	hso_node_t *node = (hso_node_t *)malloc (size);
	memset (node, 0, size);
	return node;
}

hso_node_t *hso_stack_top (hso_stack_t *stack) {
	if (stack == NULL) 
		return NULL;
	return stack->stack_top;
}

void hso_stack_pop (hso_stack_t *stack) {
	if (stack == NULL || stack->stack_top == NULL) return;

	hso_node_t *p = stack->stack_top;
	stack->stack_top = p->next;
	
	Matrix_free (&(p->M));
	free (p->index);
	free (p);
}

void hso_stack_push (hso_stack_t **stack, hso_node_t *node) {
	if (*stack == NULL) { 
		*stack = (hso_stack_t *)malloc (sizeof (hso_stack_t));
		(*stack)->stack_top = NULL;
	}
	node->next = (*stack)->stack_top;
	(*stack)->stack_top = node;
}
