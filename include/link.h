#ifndef _LINK_H
#define _LINK_H

typedef struct Node_tag {
	long long data;
	struct 	Node_tag *next;
} Node_t;

typedef struct Link_tag {
	Node_t *link_head;
	Node_t *link_tail;
	int	nNode;
	struct Link_tag *next;
} Link_t;

typedef struct List_tag {
	Link_t *list_head;
	Link_t *list_tail;
	int	nLink;
} List_t;

Node_t *newNode (long long data);
Link_t *newLink ();
List_t *newList ();

Node_t *Node_dup (Node_t *node);
Link_t *Link_dup (Link_t *link);
List_t *List_dup (List_t *list);

void Link_show (Link_t *link);
void Link_add  (Link_t **link, long long data);
void Link_del  (Link_t *link, Node_t *pNode);
void Link_del  (Link_t *link, long long data);
void Link_free (Link_t **link);
void Link_cat  (Link_t **link1, Link_t *link2);

void List_show (List_t *list);
void List_add  (List_t **list, Link_t *link);
void List_del  (List_t *list, Link_t *pLink);
void List_free (List_t **list);
void List_cat  (List_t **list1, List_t *list2);

int *Link2Array (Link_t *link);
Link_t *Array2Link (int *arr, int len);

Node_t *Link_read (Link_t *link);
Node_t *Link_search (Link_t *link, long long data);

Link_t *List_read (List_t *list);

#endif
