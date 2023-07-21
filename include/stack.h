#ifndef _STACK_H
#define _STACK_H

typedef struct mystack_ele_tag {
	struct 	mystack_ele_tag *next;
	void 	*value;
} mystack_ele_t;

typedef struct mystack_tag {
	mystack_ele_t *head;
	int	len;
} mystack_t;

void mystack_push (mystack_t **S, void *q);
void mystack_free (mystack_t **S);
void mystack_pop  (mystack_t *S);

mystack_ele_t* mystack_top (mystack_t *S);

#endif
