#include "queue.h"
#include <stdlib.h>

void queue_push (queue_t **Q, void *q) {
	if (*Q == NULL) {
		*Q = (queue_t *)malloc (sizeof (queue_t));
		(*Q)->head = NULL;
		(*Q)->tail = NULL;
		(*Q)->len = 0;
	}
	queue_ele_t *ele = (queue_ele_t *) malloc (sizeof (queue_ele_t));
	ele->next = NULL;
	ele->value = q;

	if ((*Q)->len == 0) {
		(*Q)->head = ele;
		(*Q)->tail = ele;
		(*Q)->len = 1;
	} else {
		(*Q)->tail->next = ele;
		(*Q)->tail = ele;
		(*Q)->len++;
	}
}

void queue_free (queue_t **Q) {
	while (queue_top (*Q))	
		queue_pop (*Q);
	free (*Q);
	*Q = NULL;
}

void queue_pop  (queue_t *Q) {
	if (Q==NULL || Q->len == 0 || Q->head == NULL) return;
	queue_ele_t *p = Q->head;
	Q->head = p->next;
	Q->len--;
	if (p->value)
		free (p->value);
	free (p);
}

queue_ele_t* queue_top (queue_t *Q) {
	if (Q==NULL)
		return NULL;
	return Q->head;
}

