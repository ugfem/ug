// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	"../main/defs.h"
#include	"../main/structs.h"


void add_edges(graph, new_edges, old_edges, old_ewgts, last_new_edges, last_old_edges, last_old_ewgts)
struct vtx_data **graph;	/* graph data structure */
struct edgeslist *new_edges;	/* list of edges connecting graph */
struct ilists **old_edges;	/* edges data overwritten for connecting */
struct flists **old_ewgts;	/* weights of edges overwritten */
struct edgeslist **last_new_edges;    /* last of new edges list*/
struct ilists **last_old_edges;/* last of old edges list */
struct flists **last_old_ewgts;	/* last in weights of edges overwritten */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   struct ilists *save_list;	/* space to save old edge list */
   struct ilists **end_list;	/* end of save_list */
   struct flists *save_ewgts;	/* space to save old edge weights */
   struct flists **end_ewgts;	/* end edge weights list */
   struct edgeslist *edges;	/* loops through new edges */
   float *new_ewgts;		/* new edge weights */
   int *new_list;		/* new edge list */
   int nedges;			/* number of edges a vertex has */
   int using_ewgts;		/* are edge weights being used? */
   int vtx, vtx2;		/* two vertices in edge to be added */
   int i, j;			/* loop counter */
   struct ilists *prev_old_edges;	/* prev in old edge list */
   struct flists *prev_old_ewgts;	/* space to save old edge weights */
   struct edgeslist *prev_new_edges;	/* loops through new edges */

   using_ewgts = (graph[1]->ewgts != NULL);

   *old_edges = NULL;
   *old_ewgts = NULL;
   end_list = old_edges;
   prev_old_edges = NULL;
   prev_old_ewgts = NULL;
   prev_new_edges = NULL;
   end_ewgts = old_ewgts;
   edges = new_edges;
   while (edges != NULL) {
      for (j=0; j<2; j++) { 
	 if (j==0) { vtx = edges->vtx1; vtx2 = edges->vtx2; }
	 else { vtx = edges->vtx2; vtx2 = edges->vtx1; }

         /* Copy old edge list to new edge list. */
         nedges = graph[vtx]->nedges;
         new_list = (int *) (MEM_OK = smalloc((unsigned) (nedges+1)*sizeof(int));
         if (!MEM_OK) return;
         for (i=0; i<nedges; i++) new_list[i] = graph[vtx]->edges[i];
         new_list[nedges] = vtx2;

         /* Save old edges. */
         save_list = (struct ilists *) (MEM_OK = smalloc((unsigned) sizeof(struct ilists));
         if (!MEM_OK) return;
	 save_list->list = graph[vtx]->edges;



	 /* Add new list at end of linked list. */
	save_list->prev = prev_old_edges;
	save_list->next = NULL;
	prev_old_edges = save_list;
	 *end_list = save_list;
	 end_list = &(save_list->next);



         /* Now modify graph to have new edges list. */
         graph[vtx]->nedges++;
         graph[vtx]->edges = new_list;

         /* If using edge weights, I have to modify those too. */
         if (using_ewgts) {
            new_ewgts = (float *) (MEM_OK = smalloc((unsigned) (nedges+1)*sizeof(float));
            if (!MEM_OK) return;
            for (i=1; i<nedges; i++) new_ewgts[i] = graph[vtx]->ewgts[i];
            new_ewgts[nedges] = 1;
	    new_ewgts[0] = graph[vtx]->ewgts[0] - new_ewgts[nedges];

            /* Save old edge weights. */
            save_ewgts = (struct flists *) (MEM_OK = smalloc((unsigned) sizeof(struct flists));
            if (!MEM_OK) return;
	    save_ewgts->list = graph[vtx]->ewgts;



            save_ewgts->prev = prev_old_ewgts;
            save_ewgts->next = NULL;
	    prev_old_ewgts = save_ewgts;
	    *end_ewgts = save_ewgts;
	    end_ewgts = &(save_ewgts->next);




            /* Finally, modify graph to have new edge weights. */
            graph[vtx]->ewgts = new_ewgts;
	 }
      }
      *last_old_edges = save_list;
      *last_old_ewgts = save_ewgts;

      *last_new_edges = edges;
      edges->prev = prev_new_edges;
      prev_new_edges = edges;
      edges = edges->next;
      
   }
   *last_new_edges = prev_new_edges;
}
