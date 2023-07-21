#ifndef _AVLTREE_H
#define _AVLTREE_H

typedef struct avl_node_tag {
	//
	int	factor;		// balance factor	
	int	vis;		// visiting flag
	int	depth;		// depth 
	int 	NO; 		// number in binary search tree
	int	invalid;	//

	// relation 
	struct avl_node_tag *parent;
	struct avl_node_tag *leftchild;
	struct avl_node_tag *rightchild;

	// content
	double 	key;
	void*	value;
} avl_node_t;

typedef avl_node_t* avl_tree_t;


avl_node_t* avl_node_new ();
avl_node_t* avl_node_new (double key, void *value);

void avl_node_free (avl_node_t **node);
void avl_tree_free (avl_tree_t *tree);

avl_node_t* avl_insert (avl_tree_t *tree, double key, void *value);

void avl_delete (avl_tree_t *tree, avl_node_t *p);
void avl_delete (avl_tree_t *tree, double key);

avl_node_t* avl_search (avl_tree_t tree, double key);
avl_node_t* avl_succ   (avl_node_t *p);
avl_node_t* avl_pred   (avl_node_t *p);

avl_tree_t avl_dup (avl_tree_t T);
void avl_visualize (avl_tree_t T);

#endif
