#include "stack.h"
#include <stdlib.h>

void mystack_push (mystack_t **S, void *q) {
	if (*S == NULL) {
		*S = (mystack_t *)malloc (sizeof (mystack_t));
		(*S)->head = NULL;
		(*S)->len = 0;
	}
	mystack_ele_t *ele = (mystack_ele_t *) malloc (sizeof (mystack_ele_t));
	ele->value = q;
	ele->next = (*S)->head;
	(*S)->head = ele;
	(*S)->len++;
}

void mystack_free (mystack_t **S) {
	while (mystack_top (*S))	
		mystack_pop (*S);
	free (*S);
	*S = NULL;
}

void mystack_pop  (mystack_t *S) {
	if (S==NULL || S->len == 0 || S->head == NULL) return;
	mystack_ele_t *p = S->head;
	S->head = p->next;
	S->len--;
	if (p->value)
		free (p->value);
	free (p);
}

mystack_ele_t* mystack_top (mystack_t *S) {
	if (S==NULL)
		return NULL;
	return S->head;
}

