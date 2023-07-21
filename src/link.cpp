#include "link.h"

#include <stdlib.h>
#include <stdio.h>

Node_t *newNode (long long data) {
	Node_t *node = (Node_t *)malloc (sizeof (Node_t));
	if (node == NULL) {
		fprintf (stderr, "Can't allocate memory\n");
		exit (EXIT_FAILURE);
	}
	node->data = data;
	node->next = NULL;
	return node;
}

Link_t *newLink (){
	Link_t *link = (Link_t *)malloc (sizeof (Link_t));
	if (link == NULL) {
		fprintf (stderr, "Can't allocate memory\n");
		exit (EXIT_FAILURE);
	}
	link->link_head = NULL;
	link->link_tail = NULL;
	link->nNode = 0;
	link->next = NULL;
	return link;
}

List_t *newList () {
	List_t *list = (List_t *)malloc (sizeof (List_t));
	if (list == NULL) {
		fprintf (stderr, "Can't allocate memory\n");
		exit (EXIT_FAILURE);
	}
	list->list_head = NULL;
	list->list_tail = NULL;
	list->nLink = 0;
	return list;
}

Node_t *Node_dup (Node_t *node) {
	if (node == NULL) return NULL;
	Node_t *pnode = (Node_t *)malloc (sizeof (Node_t));
	pnode->data = node->data;
	pnode->next = NULL;
	return pnode;
}
Link_t *Link_dup (Link_t *link) {
	if (link==NULL) return NULL;
	Link_t *pLink = newLink();
	Node_t *p = link->link_head;
	Node_t *pNode = NULL;
	int i;

	i = link->nNode;
	if (link->nNode > 0) {
		pNode = Node_dup (p);
		pLink->link_head = pNode;
		pLink->link_tail = pNode;
		pLink->nNode = 1;
		p=p->next;

		while (p != NULL && i > 0) {
			pNode = Node_dup (p);
			pLink->link_tail->next = pNode;
			pLink->link_tail = pNode;
			pLink->nNode++;
			p=p->next;
			i--;
		}
	}
	return pLink;	
}
List_t *List_dup (List_t *list) {
	if (list==NULL) return NULL;
	List_t *pList = newList();
	Link_t *p = list->list_head;
	Link_t *pLink = NULL;
	int i;

	i = list->nLink;
	if (list->nLink > 0) {
		pLink = Link_dup (p);
		pList->list_head = pLink;
		pList->list_tail = pLink;
		pList->nLink = 1;
		p=p->next;

		while (p != NULL && i > 0) {
			pLink = Link_dup (p);
			pList->list_tail->next = pLink;
			pList->list_tail = pLink;
			pList->nLink++;
			p=p->next;
			i--;
		}
	}
	return pList;	
}

void Link_show (Link_t *link) {
	int i; 
	Node_t *pNode = NULL;

	if (link == NULL || link->nNode == 0) {
		printf ("link is empty\n");
		return;
	}

	i = link->nNode;
	pNode = link->link_head;
	while (pNode != NULL && i > 0) {
		printf ("%lld ", pNode->data);
		pNode = pNode->next;
		i--;
	}
	printf ("\n");
}


void Link_add (Link_t **link, long long data) {
	Node_t *node = newNode (data);

	if ((*link) == NULL) (*link) = newLink ();

	if ((*link)->nNode == 0) {
		(*link)->link_head = node;
		(*link)->link_tail = node;
		(*link)->nNode = 1;
	} else {
		(*link)->link_tail->next = node;
		(*link)->link_tail = node;
		(*link)->nNode = (*link)->nNode + 1;	
	}
}

void Link_del (Link_t *link, Node_t *pNode) {
	Node_t *p = NULL;
	int i;

	if (link == NULL || link->nNode < 1 || pNode == NULL) return;

	i = link->nNode;
	if (link->link_head == pNode) {
		link->link_head = pNode->next;
		if (pNode->next == NULL)
			link->link_tail = NULL;
		link->nNode--;
	} else {
		p = link->link_head;
		while (p != NULL && i > 0 && p->next != pNode) {
			p = p->next;
			i--;
		}
		if (p == NULL || i < 1) return;
		p->next = pNode->next;
		if (pNode == link->link_tail) 
			link->link_tail = p;
		link->nNode--;
	}
}

void Link_del  (Link_t *link, long long data) {
	Link_del (link, Link_search (link, data));
}

void Link_free (Link_t **link) {
	if (link == NULL || (*link) == NULL) return;

	int i = (*link)->nNode;
	Node_t *p = (*link)->link_head;
	while (p != NULL && i > 0) {
		(*link)->link_head = p->next;
		free (p);
		p = (*link)->link_head;
		i--;
	}
	free (*link);
	*link = NULL;
}

void Link_cat (Link_t **link1, Link_t *link2) {
	Link_t *pLink = NULL;

	if (link2 == NULL)     return;
	if ((*link1) == NULL) 	  *(link1) = newLink();

	pLink = Link_dup (link2);
	if ((*link1)->nNode == 0) {
		(*link1)->link_tail = pLink->link_tail;
		(*link1)->link_head = pLink->link_head;
		(*link1)->nNode = pLink->nNode;
	} else if ((*link1)->nNode > 0){
		(*link1)->link_tail->next = pLink->link_head;
		(*link1)->link_tail = pLink->link_tail;
		(*link1)->nNode += pLink->nNode;
	}
}

void List_show (List_t *list) {
	int i;
	Link_t *pLink = NULL;

	if (list==NULL || list->nLink == 0) {
		printf ("list is empty\n");
		return;
	}

	i = list->nLink;
	pLink = list->list_head;
	while (i > 0 && pLink != NULL) {
		Link_show (pLink);
		pLink = pLink->next;
		i--;
	}
}

void List_add (List_t **list, Link_t *link) {
	if (link == NULL || link->nNode == 0) return;

	if ((*list) == NULL) (*list) = newList ();
	Link_t *pLink = Link_dup (link);

	if ((*list)->nLink == 0) {
		(*list)->list_head = pLink;
		(*list)->list_tail = pLink;
		(*list)->nLink = 1;
	} else {
		(*list)->list_tail->next = pLink;
		(*list)->list_tail = pLink;
		(*list)->nLink = (*list)->nLink + 1;	
	}
}

void List_del (List_t *list, Link_t *pLink) {
	if (pLink == NULL || list == NULL) return;
	Link_t *p = NULL;
	int i = list->nLink;

	if (list->list_head == pLink) {
		list->list_head = pLink->next;
		if (pLink->next == NULL)
			list->list_tail = NULL;
		list->nLink--;
	} else {
		p = list->list_head;
		while (p != NULL && i > 0 && p->next != pLink ) {
			p = p->next;
			i--;
		}
		if (p == NULL || i < 1) 
			return;
		p->next = pLink->next;
		if (pLink == list->list_tail) 
			list->list_tail = p;
		list->nLink--;
	}
}

void List_free (List_t **list) {
	if (list==NULL || (*list) == NULL) return;
	int i = (*list)->nLink;

	Link_t *p = (*list)->list_head;
	while (p != NULL && i > 0) {
		(*list)->list_head = p->next;
		Link_free (&p);
		p = (*list)->list_head;
		i--;
	}
	free (*list);
	*list = NULL;
}

void List_cat (List_t **list1, List_t *list2) {
	List_t * pList = NULL;

	if (list2 == NULL)     return;
	if ((*list1) == NULL) 	  *(list1) = newList();

	pList = List_dup (list2);
	if ((*list1)->nLink == 0) {
		(*list1)->list_tail = pList->list_tail;
		(*list1)->list_head = pList->list_head;
		(*list1)->nLink = pList->nLink;
	} else if ((*list1)->nLink > 0){
		(*list1)->list_tail->next = pList->list_head;
		(*list1)->list_tail = pList->list_tail;
		(*list1)->nLink += pList->nLink;
	}
}

int *Link2Array (Link_t *link) {
	size_t 	size;
	int  	*buf, i;
	Node_t 	*p = NULL;

	if (link==NULL) return NULL;

	size = (link->nNode + 2)*sizeof (int);
	buf = (int *)malloc (size);

	buf[0] = link->nNode;
	i = 1;
	p = link->link_head;	
	while (p != NULL && i <= link->nNode) {
		buf[i] = p->data;
		p = p->next;
		i++;
	}

	return buf;
}

Link_t *Array2Link (int *arr, int len) {
	Link_t *link = NULL;
	int i;

	for (i=len-1; i>=0; i--)
		Link_add (&link, arr[i]);
	return link;
}


static Node_t *nextNode[100];
static Link_t *curLink[100];
Node_t *Link_read (Link_t *link) {
	int i;
	Node_t *p = NULL;
	if (link == NULL || link->link_head == NULL) return NULL;

	for (i=0; i<100; i++) {
		if (curLink[i] == link) {
			p = nextNode[i];
			if (p != NULL)
				nextNode[i] = p->next;
			else
				curLink[i] = NULL;
			return p;
		}
	}

	for (i=0; i<100; i++) {
		if (curLink[i] == NULL) {
			p = link->link_head;
			curLink[i] = link;
			nextNode[i] = p->next; 
			return p;
		}
	}
	
	fprintf (stderr, "the number of the link concurrently read is more than 100\n");
	exit (EXIT_FAILURE);
}

static Link_t *nextLink[100];
static List_t *curList[100];
Link_t *List_read (List_t *list) {
	Link_t *p = NULL;
	int i;

	if (list == NULL || list->list_head == NULL) return NULL;

	for (i=0; i<100; i++) {
		if (curList[i] == list) {
			p = nextLink[i];
			if (p != NULL)
				nextLink[i] = p->next;
			else
				curList[i] = NULL;
			return p;
		}
	}

	for (i=0; i<100; i++) {
		if (curList[i] == NULL) {
			p = list->list_head;
			curList[i] = list;
			nextLink[i] = p->next; 
			return p;
		}
	}
	
	fprintf (stderr, "the number of the list concurrently read is more than 100\n");
	exit (EXIT_FAILURE);
}

Node_t *Link_search (Link_t *link, long long data) {
	Node_t *p = NULL;

	p = link->link_head;
	while (p != NULL && p->data != data) {
		p=p->next;
	}

	return p;
}
