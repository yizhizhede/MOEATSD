#ifndef _QUEUE_H
#define _QUEUE_H

typedef struct queue_ele_tag {
	struct 	queue_ele_tag *next;
	void 	*value;
} queue_ele_t;

typedef struct queue_tag {
	queue_ele_t *head;
	queue_ele_t *tail;
	int	len;
} queue_t;

void queue_push (queue_t **Q, void *q);
void queue_free (queue_t **Q);
void queue_pop  (queue_t *Q);

queue_ele_t* queue_top (queue_t *Q);

#endif
