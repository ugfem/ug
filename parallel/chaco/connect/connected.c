// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	"../main/defs.h"
#include	"../main/structs.h"


void make_connected(graph, nvtxs, nedges, mark, vtxlist, cdata)
/* Add edges to make graph connected. */
struct vtx_data **graph;	/* graph data structure */
int nvtxs;			/* number of vertices in graph */
int *nedges;			/* number of edges in graph */
short *mark;			/* space for nvtxs+1 shorts */
int *vtxlist;			/* space for nvtxs ints */
struct connect_data **cdata;	/* space for connectivity data */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   struct edgeslist *new_edges;	/* list of edges connecting graph */
   int nadded;			/* number of edges being added */
   int find_edges();
   void add_edges();

   /* First find edges needed to make graph connected. */
   nadded = find_edges(graph, nvtxs, mark, vtxlist, &new_edges);
   if (!MEM_OK) return;

   /* Now add these needed edges to graph data structure if needed. */
   if (nadded == 0) {
      *cdata = NULL;
   }
   else {
      *cdata = (struct connect_data *) (MEM_OK = smalloc(sizeof(struct connect_data));
      if (!MEM_OK) return;
      (*cdata)->new_edges = new_edges;
      (*cdata)->old_edges = NULL;
      (*cdata)->old_ewgts = NULL;
      (*cdata)->last_new_edges = NULL;
      (*cdata)->last_old_edges = NULL;
      (*cdata)->last_old_ewgts = NULL;
      add_edges(graph, (*cdata)->new_edges, &(*cdata)->old_edges,
	            &(*cdata)->old_ewgts,&(*cdata)->last_new_edges,
                &(*cdata)->last_old_edges,&(*cdata)->last_old_ewgts);
      if (!MEM_OK) return;
      *nedges += nadded;
   }
}


void make_unconnected(graph, nedges, cdata)
/* Restore graph to its pristine state and free space for connectivity. */
struct vtx_data **graph;	/* graph data structure */
int *nedges;			/* number of edges in graph */
struct connect_data **cdata;	/* space for connectivity data */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   struct ilists *old_edges = NULL;	/* edges overwritten for connecting */
   struct flists *old_ewgts = NULL;	/* weights of edges overwritten */
   struct edgeslist *new_edges;	/* list of edges connecting graph */
   struct ilists *tempi;	/* used for freeing space */
   struct flists *tempf;	/* used for freeing space */
   struct edgeslist *tempe;	/* used for freeing edgelist space */
   struct edgeslist *edges;	/* loops through new edges */
   int using_ewgts;		/* are edge weights being used? */
   int vtx;			/* vertex in an added edge */
   int j;			/* loop counters */
   

   if (*cdata == NULL) return;

   old_edges = (*cdata)->last_old_edges;
   old_ewgts = (*cdata)->last_old_ewgts;
   new_edges = (*cdata)->last_new_edges;
   sfree((char *) *cdata);
   *cdata = NULL;

	if (graph[1]->ewgts != NULL)
    {
      using_ewgts = (graph[1]->ewgts != NULL);
	}
	else using_ewgts = 0;

   edges = new_edges;
   while (edges != NULL) {
      /* Restore edges and weights to original status. */
      (*nedges)--;
      for (j=0; j<2; j++) { 

	 if (j==0) vtx = edges->vtx2;
	 else vtx = edges->vtx1;

	 sfree((char *) graph[vtx]->edges);
	 graph[vtx]->edges = old_edges->list;
	 graph[vtx]->nedges--;
	 tempi = old_edges;
	 old_edges = old_edges->prev;
	 sfree((char *) tempi);

	 if (using_ewgts) {
	    sfree((char *) graph[vtx]->ewgts);
	    graph[vtx]->ewgts = old_ewgts->list;
	    tempf = old_ewgts;
	    old_ewgts = old_ewgts->prev;
	    sfree((char *) tempf);
	 }
      }
      tempe = edges;
      edges = edges->prev;
      sfree((char *) tempe);
   }
}


/* Print out the added edges. */
void print_connected(cdata)
struct connect_data *cdata;	/* space for connectivity data */
{
   struct edgeslist *edges;	/* loops through new edges */

   if (cdata == NULL){char buf[150]; sprintf(buf,"No phantom edges\n");UserWrite(buf);}
   else {
     {char buf[150]; sprintf(buf,"Phantom edges: ");UserWrite(buf);}
      edges = cdata->new_edges;
      while (edges != NULL) {
	{char buf[150]; sprintf(buf,"(%d,%d) ", edges->vtx1, edges->vtx2);UserWrite(buf);}
	 edges = edges->next;
      }
     {char buf[150]; sprintf(buf,"\n");UserWrite(buf);}
   }
}
