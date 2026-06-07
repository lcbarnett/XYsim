#ifndef VLIST_H
#define VLIST_H

// A minimal "vortex list" (FILO)

typedef struct vlist_node {
	int                wnum; // winding number
	int                vidx; // vortex lattice index
	struct vlist_node* next;
} vlist_t;

static inline vlist_t* vlist_push(vlist_t* const head, const int wnum, const int vidx)
{
	vlist_t* pnode = malloc(sizeof(vlist_t)); // pointer to new list node
	pnode->wnum = wnum;
	pnode->vidx = vidx;
	pnode->next = head;
	return pnode;
}

// Clean up memory to keep things professional
static inline void vlist_free(vlist_t* head)
{
	while (head) {
		vlist_t* const ptmp = head;
		head = head->next;
		free(ptmp);
	}
}

#endif // VLIST_H
