// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	<string.h>
#include	"../main/defs.h"
#include	"../main/params.h"


void input_graph(fin, inname, start, adjacency, nvtxs, vweights, eweights)
FILE *fin;			/* input file */
char *inname;			/* name of input file */
int **start;			/* start of edge list for each vertex */
int **adjacency;		/* edge list data */
int *nvtxs;			/* number of vertices in graph */
int **vweights;			/* vertex weight list data */
int **eweights;			/* edge weight list data */
{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern int DEBUG_INPUT;	/* echo that input file read successful? */
   int *adjptr;			/* loops through adjacency data */
   int *ewptr;			/* loops through edge weight data */
   char line[LINE_LENGTH+1];	/* line of date from file */
   char *scanptr;		/* mark position within line */	
   char *oldptr;		/* mark previous position within line */
   int narcs;			/* number of edges in graph */
   int nedge;			/* loops through edges for each vertex */
   int flag;			/* condition indicator */
   int vtx;			/* vertex in graph */
   int sum_edges;		/* total number of edges read so far */
   int option;			/* input option */
   int edgewts;			/* are edge weights in input file? */
   int vtxwts;			/* are vertex weights in input file? */
   int vtxnums;			/* are vertex numbers in input file? */
   int vertex;			/* current vertex being read */
   int weight;			/* weight being read */
   int eweight;			/* edgeweight being read */
   int neighbor;		/* neighbor of current vertex */
   int j;			/* loop counters */
   
   void exit();
   long strtol();

   /* Read first line  of input (= nvtxs, narcs, option). */
   /* The (decimal) digits of the option variable mean:
	 1's digit not zero => input edge weights
	 10's digit not zero => input vertex weights
	 100's digit not zero => include vertex numbers
   */
   flag = FALSE;
   while (!flag && fgets(line, LINE_LENGTH, fin) != NULL) {
      if (line[0] != '%') {
	 j = sscanf(line, "%d%d%d", nvtxs, &narcs, &option);
	 if (j < 3) option = 0;
	 flag = TRUE;
	 if (j == 1) narcs = 0;
      }
   }

   edgewts = option - 10*(option/10);
   option /= 10;
   vtxwts = option - 10*(option/10);
   option /= 10;
   vtxnums = option - 10*(option/10);

   /* Allocate space for rows and columns. */
   *start = (int *) (MEM_OK = smalloc((unsigned) (*nvtxs+2)*sizeof(int));
   if (!MEM_OK) return;
   if (narcs != 0) {*adjacency = (int *) (MEM_OK = smalloc((unsigned) 2*narcs*sizeof(int));
                    if (!MEM_OK) return;
                   }
   else *adjacency = NULL;

   if (vtxwts) {*vweights = (int *) (MEM_OK = smalloc((unsigned) (*nvtxs+1)*sizeof(int));
                if (!MEM_OK) return;
               }
   else *vweights = NULL;

   if (edgewts) {*eweights = (int *) (MEM_OK = smalloc((unsigned) 2*narcs*sizeof(int));
                 if (!MEM_OK) return;
                }
   else *eweights = NULL;

   adjptr = *adjacency;
   ewptr = *eweights;

   sum_edges = 0;
   for (vtx=1; fgets(line, LINE_LENGTH, fin) != NULL; ) {
      nedge = -1;
      if (line[0] != '%') {
	 scanptr = line;

	 /* If vertices are entered unordered, read vertex number. */
	 if (vtxnums) {
	    vertex = (int) strtol(scanptr,&scanptr,10);
	    if (scanptr == line) {
	      {char buf[150]; sprintf(buf,"Input error in file %s:", inname);UserWrite(buf);}
	      {char buf[150]; sprintf(buf," no vertex number in line %d.\n", vtx);UserWrite(buf);}
	       exit(1);
	    }
	 }
	 else vertex = vtx;

	 if (vertex > *nvtxs || vertex <= 0) {
	   {char buf[150]; sprintf(buf,"Input error in file %s:", inname);UserWrite(buf);}
	   {char buf[150]; sprintf(buf," nvtxs=%d, but data for vertex %d was read.\n", *nvtxs, vertex);UserWrite(buf);}
	    exit(1);
	 }

	 /* If vertices are weighted, read vertex weight. */
	 if (vtxwts) {
	    oldptr = scanptr;
	    weight = (int) strtol(scanptr,&scanptr,10);
	    if (scanptr == oldptr) {
	      {char buf[150]; sprintf(buf,"Input error in file %s:", inname);UserWrite(buf);}
	      {char buf[150]; sprintf(buf," no weight for vertex %d.\n", vertex);UserWrite(buf);}
	       exit(1);
	    }
	    if (weight < 0) {
	      {char buf[150]; sprintf(buf,"Input error in file %s:", inname);UserWrite(buf);}
	      {char buf[150]; sprintf(buf," negative weight entered for vertex %d.\n", vertex);UserWrite(buf);}
	       exit(1);
	    }
	    (*vweights)[vertex] = weight;
	 }

	 nedge = 0;

	 oldptr = scanptr;
	 /* Read number of adjacent vertex. */
	 neighbor = (int) strtol(scanptr,&scanptr,10);

	 while (scanptr != oldptr) {

	    /* Check for edge only entered once. */
	    if (!vtxnums && neighbor < vertex && neighbor != 0) {
	       option = FALSE;
	       for (j=(*start)[neighbor]; j<(*start)[neighbor+1]; j++) {
		  if ((*adjacency)[j] == vertex) option = TRUE;
	       }
	       if (!option){char buf[150]; sprintf(buf," Edge (%d, %d) entered but not (%d, %d)\n", vertex, neighbor, neighbor, vertex);UserWrite(buf);}
	    }

	    if (neighbor > *nvtxs) {
	      {char buf[150]; sprintf(buf,"Input error in file %s:", inname);UserWrite(buf);}
	      {char buf[150]; sprintf(buf," nvtxs=%d, but edge (%d,%d) was input.\n", *nvtxs, vertex, neighbor);UserWrite(buf);}
	       exit(1);
	    }
	    if (neighbor < 0) {
	      {char buf[150]; sprintf(buf,"Input error in file %s:", inname);UserWrite(buf);}
	      {char buf[150]; sprintf(buf," negative vertex in edge (%d,%d) was input.\n", vertex, neighbor);UserWrite(buf);}
	       exit(1);
	    }

	    if (neighbor != 0) {	/* Allow Julie's input format. */
	       *adjptr++ = neighbor;
	       nedge++;
	    }

	    if (edgewts) {	/* Read edge weight if it's being input. */
	       oldptr = scanptr;
	       eweight = (int) strtol(scanptr,&scanptr,10);

	       if (scanptr == oldptr) {
	         {char buf[150]; sprintf(buf,"Input error in file %s:", inname);UserWrite(buf);}
		 {char buf[150]; sprintf(buf," no weight for edge (%d,%d).\n", vertex, neighbor);UserWrite(buf);}
	          exit(1);
	       }
	       if (eweight < 0) {
	         {char buf[150]; sprintf(buf,"Input error in file %s:", inname);UserWrite(buf);}
	         {char buf[150]; sprintf(buf," negative weight entered for edge (%d,%d).\n", vertex, neighbor);UserWrite(buf);}
	          exit(1);
	       }

	       *ewptr++ = eweight;
	    }
	    oldptr = scanptr;
	    /* Read number of adjacent vertex. */
	    neighbor = (int) strtol(scanptr,&scanptr,10);
	 }

	 (*start)[vertex] = sum_edges;
	 sum_edges += nedge;
if (!vtxnums) (*start)[vertex+1] = sum_edges;
	 vtx++;
      }
   }
   (*start)[(*nvtxs)+1] = 2*narcs;

   if (vtx != 1) {		/* Normal file was read. */
      nedge = adjptr-(*adjacency);
      if (nedge != 2*narcs) {
        {char buf[150]; sprintf(buf,"Input error in file %s:", inname);UserWrite(buf);}
        {char buf[150]; sprintf(buf," narcs=%d, but %d edges were read.\n", narcs, nedge/2);UserWrite(buf);}
         exit(1);
      }
   }

   else {
      /* Graph was empty => must be using inertial method. */
      sfree((char *) *start);
      if (*adjacency != NULL) sfree((char *) *adjacency);
      if (*vweights != NULL) sfree((char *) *vweights);
      if (*eweights != NULL) sfree((char *) *eweights);
      *start = NULL;
      *adjacency = NULL;
   }

   fclose(fin);

   if (DEBUG_INPUT > 0) {
     {char buf[150]; sprintf(buf,"Done reading graph file %s; %d vertices, %d edges.\n", inname, *nvtxs, narcs);UserWrite(buf);}
   }
}
