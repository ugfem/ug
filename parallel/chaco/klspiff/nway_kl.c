// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include	<stdio.h>
#include	"../main/structs.h"
#include	"../main/params.h"
#include	"../main/defs.h"

/* Idea:
   'buckets[i][j]' is a set of buckets to sort moves from i to j.
   listspace[i] is space for lists in buckets[i][j].
   Loop through all nonequal pairs [i][j], taking the first element
   in each list.  Compare them all to find the largest allowed move.
   Make that move, and save it in movelist.
*/


void nway_kl(graph, nvtxs, buckets, listspace, tops, dvals, sets,
     maxdval, nsets, goal, hops, vwgt_max)
struct vtx_data **graph;	/* data structure for graph */
int nvtxs;			/* number of vtxs in graph */
struct bilist ****buckets;	/* array of lists for bucket sort */
struct bilist **listspace;	/* list data structure for each vertex */
int **tops;			/* 2-D array of top of each set of buckets */
int **dvals;			/* d-values for each transition */
short *sets;			/* processor each vertex is assigned to */
int maxdval;			/* maximum d-value for a vertex */
int nsets;			/* number of sets divided into */
double *goal;			/* desired set sizes */
short (*hops)[MAXSETS];		/* cost of set transitions */
int vwgt_max;			/* largest vertex weight */

/* Suaris and Kedem algorithm for quadrisection, generalized to an */
/* arbitrary number of sets, with intra-set cost function specified by hops. */
/* Note: this is for a single divide step. */
/* Also, sets contains an intial (possibly crummy) partitioning. */

{
   extern Heap   *heap;     /* pointer to heap of multigrid */
   extern double *MEM_OK;   /* variable for memory overflow exeception */
   extern double DOUBLE_MAX;	/* large double precision value */
   extern int KL_BAD_MOVES;	/* # bad moves in a row to stop KL */
   extern int DEBUG_KL;		/* debug flag for KL */
   extern int KL_RANDOM;	/* use randomness in KL? */
   extern int KL_NTRIES_BAD;	/* number of unhelpful passes before quitting */
   extern int KL_UNDO_LIST;	/* should I back out of changes or start over? */
   struct bilist *movelist;	/* list of vtxs to be moved */
   struct bilist **endlist;	/* end of movelists */
   struct bilist *bestptr;	/* best vertex in linked list */
   struct bilist *bptr;		/* loops through bucket list */
   float *ewptr;		/* loops through edge weights */
   double *locked;		/* weight of vertices locked in a set */
   double *loose;		/* weight of vtxs that can move from a set */
   int *bspace;			/* list of active vertices for bucketsort */
   double *weightsum;		/* sum of vweights for each partition */
   double *startweight;		/* initial sum of vweights for each set */
   int *edges;			/* edge list for a vertex */
   double delta;		/* desire of sets to change size */
   double bestdelta;		/* strongest delta value */
   double deltaplus;		/* largest negative deviation from goal size */
   double deltaminus;		/* largest negative deviation from goal size */
   int list_length;		/* how long is list of vertices to bucketsort? */
   int balanced;		/* is partition balanced? */
   int temp_balanced;		/* is intermediate partition balanced? */
   int bestvtx;			/* best vertex to move */
   int bestval;			/* best change in value for a vtx move */
   int bestfrom, bestto;	/* sets best vertex moves between */
   int vweight;			/* weight of best vertex */
   int gtotal;			/* sum of changes from moving */
   int improved;		/* total improvement from KL */
   double bestg;		/* maximum gtotal found in KL loop */
   int beststep;		/* step where maximum value occurred */
   int neighbor;		/* neighbor of a vertex */
   int step_cutoff;		/* number of negative steps in a row allowed */
   int cost_cutoff;		/* amount of negative d-values allowed */
   int neg_steps;		/* number of negative steps in a row */
   int neg_cost;		/* decrease in sum of d-values */
   int vtx;			/* vertex number */
   int dval;			/* dval of a vertex */
   int group;			/* set that a vertex is assigned to */
   int diff;			/* change in a d-value */
   int stuck1st, stuck2nd;	/* how soon will moves become disallowed? */
   int beststuck1, beststuck2;	/* best stuck values used for tie-breaking */
   int eweight;			/* a particular edge weight */
   int worth_undoing;		/* is it worth undoing list? */
   float undo_frac;		/* fraction of vtxs indicating worth of undoing */
   int step;			/* loops through movements of vertices */
   int parity;			/* sort forwards or backwards? */
   int done;			/* has termination criteria been achieved? */
   int nbad;			/* number of unhelpful passes in a row*/
   int npass;			/* total number of passes */
   int nbadtries;		/* number of unhelpful passes before quitting */
   int using_ewgts;		/* are edge weights being used? */
   int enforce_balance;		/* force a balanced partition? */
   int enforce_balance_hard;	/* really force a balanced partition? */
   int i, j, k, l;		/* loop counters */
   
   double drandom();
   int make_kl_list();
   void bucketsorts(), pbuckets(), removebilist(), movebilist();

	bestdelta = 0.0;

   using_ewgts = (graph[1]->ewgts != NULL);

   nbadtries = KL_NTRIES_BAD;

   enforce_balance = FALSE;
   enforce_balance_hard = FALSE;

   undo_frac = .3;

   bspace = (int *) (MEM_OK = smalloc((unsigned) nvtxs*sizeof(int));
   if (!MEM_OK) return;
   weightsum = (double *) (MEM_OK = smalloc((unsigned) nsets*sizeof(double));
   if (!MEM_OK) return;
   startweight = (double *) (MEM_OK = smalloc((unsigned) nsets*sizeof(double));
   if (!MEM_OK) return;
   locked = (double *) (MEM_OK = smalloc((unsigned) nsets*sizeof(double));
   if (!MEM_OK) return;
   loose = (double *) (MEM_OK = smalloc((unsigned) nsets*sizeof(double));
   if (!MEM_OK) return;

   step_cutoff = KL_BAD_MOVES;
   cost_cutoff = maxdval*step_cutoff/7;
   if (cost_cutoff < step_cutoff) cost_cutoff = step_cutoff;

   /* Compute the sum of vertex weights for each set. */
   for (i=0; i<nsets; i++) startweight[i] = 0;
   for (i=1; i<=nvtxs; i++) startweight[sets[i]] += graph[i]->vwgt;
   deltaminus = deltaplus = 0;
   for (i=0; i<nsets; i++) {
      if (startweight[i] - goal[i] > deltaplus) {
	 deltaplus = startweight[i] - goal[i];
      }
      else if (goal[i] - startweight[i] > deltaminus) {
	 deltaminus = goal[i] - startweight[i];
      }
   }
   balanced = (deltaplus + deltaminus <= vwgt_max);

   list_length = nvtxs;
   parity = FALSE;
   eweight = 1;
   nbad = 0;
   npass = 0;
   improved = 0;
   done = FALSE;
   while (!done) {
      npass++;

      /* Initialize various quantities. */
      for (i=0; i<nsets; i++) {
	 for (j=0; j<nsets; j++) tops[i][j] = 2*maxdval;
	 weightsum[i] = startweight[i];
	 loose[i] = weightsum[i];
	 locked[i] = 0;
      }

      gtotal = 0;
      bestg = -DOUBLE_MAX;
      beststep = -1;

      movelist = NULL;
      endlist = &movelist;

      neg_steps = 0;

      /* Compute the initial d-values, and bucket-sort them. */
/*{char buf[150]; sprintf(buf," Calling bucketsorts, nvtxs=%d, list_length=%d\n", nvtxs, list_length);UserWrite(buf);}*/
      bucketsorts(graph, nvtxs, buckets, listspace, dvals, sets, maxdval,
		  nsets, parity, hops, bspace, list_length, npass);
      parity = !parity;

      if (DEBUG_KL > 1) {
         pbuckets(buckets, listspace, maxdval, nsets);
      }

      /* Now determine the set of K-L moves. */

      for (step=1; ; step++) {

         /* Find the highest d-value in each set. */
	 /* But only consider moves from large to small sets. */
	 /* Break ties in some nonarbitrary manner. */
	 bestval = -maxdval -1;
	 for (i=0; i<nsets; i++)	/* Only move from large sets ... */
	  if (weightsum[i] > goal[i] ||
	      (weightsum[i] == goal[i] && !enforce_balance_hard)) {
	    for (j=0; j<nsets; j++)		/* ... to small sets. */
	     if (i!=j && ((weightsum[j] <= goal[j] && !enforce_balance_hard) ||
		      (weightsum[i]-goal[i]-weightsum[j]+goal[j]> vwgt_max))) {
	       /* Find the best move from i to j. */
	       k = tops[i][j];
	       bptr = buckets[i][j][k];
	       while (bptr == NULL && k >= 0) {
	          bptr = buckets[i][j][--k];
	       }
	       tops[i][j] = k;

	       /* Is it the best move seen so far? */
	       if (k-maxdval > bestval) {
		  bestval = k-maxdval;
	          l = (j > i) ? j-1 : j;
	          bestvtx = ((int) buckets[i][j][k] - (int) listspace[l])/
			 sizeof(struct bilist);
		  bestto = j;
		  /* DO I NEED ALL THIS DATA?  Just used in tie breaking. */
		  bestdelta = (weightsum[i]-goal[i]) - (weightsum[j]-goal[j]);
		  beststuck1 = min(loose[i], goal[j]-locked[j]);
		  beststuck2 = max(loose[i], goal[j]-locked[j]);
	       }

	       else if (k-maxdval == bestval) {
	          /* Tied.  Is it better balanced than current best? */
	          l = (j > i) ? j-1 : j;
	          vtx = ((int) buckets[i][j][k] - (int) listspace[l])/
			 sizeof(struct bilist);
		  /* Is it best balanced so far? */
		  /* If tied, move among sets with most freedom. */
		  delta = (weightsum[i]-goal[i]) - (weightsum[j]-goal[j]);
		  stuck1st = min(loose[i], goal[j]-locked[j]);
		  stuck2nd = max(loose[i], goal[j]-locked[j]);

		  /* NOTE: Randomization in this check isn't good if more */
		  /* than two guys are tied. */
		  if (delta > bestdelta ||
			(delta == bestdelta && (stuck1st > beststuck1 ||
			  (stuck1st == beststuck1 &&
			   (stuck2nd > beststuck2 || 
			    (stuck2nd == beststuck2 && 
			     (KL_RANDOM && drandom() < .5))))))) {
		     bestval = k-maxdval;
		     bestvtx = vtx;
		     bestto = j;
		     bestdelta = delta;
		     beststuck1 = stuck1st;
		     beststuck2 = stuck2nd;
		  }
	       }
	    }
	 }

	 if (bestval == -maxdval - 1) {		/* No allowed moves */
	    if (DEBUG_KL > 0) {
	      {char buf[150]; sprintf(buf,"No KL moves at step %d.  bestg = %g at step %d.\n", step, bestg, beststep);UserWrite(buf);}
	    }
	    break;
	 }

	 bestptr = &(listspace[0][bestvtx]);
	 bestfrom = sets[bestvtx];

	 vweight =  graph[bestvtx]->vwgt;
	 weightsum[bestto] += vweight;
	 weightsum[bestfrom] -= vweight;
	 loose[bestfrom] -= vweight;
	 locked[bestto] += vweight;

	 if (enforce_balance) {	/* Check if this partition is balanced. */
	    deltaminus = deltaplus = 0;
	    for (i=0; i<nsets; i++) {
	       if (weightsum[i] - goal[i] > deltaplus) {
	 	  deltaplus = weightsum[i] - goal[i];
	       }
	       else if (goal[i] - weightsum[i] > deltaminus) {
	 	  deltaminus = goal[i] - weightsum[i];
	       }
	    }
	    temp_balanced = (deltaminus + deltaplus <= vwgt_max);
	 }

         gtotal += bestval;
         if (gtotal > bestg && (!enforce_balance || temp_balanced)) {
	    bestg = gtotal;
	    beststep = step;
         }

	 if (DEBUG_KL > 1) {
	   {char buf[150]; sprintf(buf,"At KL step %d, bestvtx=%d, bestval=%d (%d-> %d)\n", step, bestvtx, bestval, bestfrom, bestto);UserWrite(buf);}
	 }

	 /* Monitor the stopping criteria. */
	 if (bestval < 0) {
	    neg_steps++;
	    if (bestg != -DOUBLE_MAX) neg_cost = bestg - gtotal;
	    else neg_cost = -maxdval - 1;
	    if (neg_steps > step_cutoff || neg_cost > cost_cutoff) {
	       if (DEBUG_KL > 0) {
		  if (neg_steps > step_cutoff) {
		    {char buf[150]; sprintf(buf,"KL step cutoff at step %d.  bestg = %g at step %d.\n", step, bestg, beststep);UserWrite(buf);}
		  }
		  else if (neg_cost > cost_cutoff) {
		    {char buf[150]; sprintf(buf,"KL cost cutoff at step %d.  bestg = %g at step %d.\n", step, bestg, beststep);UserWrite(buf);}
		  }
	       }
	       break;
	    }
	 }
	 else if (bestval > 0) {
	    neg_steps = 0;
	 }

         /* Remove vertex from its buckets, and flag it as finished. */
	 l = 0;
	 for (k=0; k<nsets; k++) {
	    if (k != bestfrom) {
	       dval = dvals[bestvtx][l] + maxdval;
               removebilist(&listspace[l][bestvtx],
			    &buckets[bestfrom][k][dval]);
	       l++;
	    }
	 }

	 /* Is there a better way to do this? */
         sets[bestvtx] = -sets[bestvtx] - 1;

	 /* Set up the linked list of moved vertices. */
	 bestptr->next = NULL;
	 bestptr->prev = (struct bilist *) bestto;
	 *endlist = bestptr;
	 endlist = &(bestptr->next);

         /* Now update the d-values of all the neighbors */
	 edges = &(graph[bestvtx]->edges[1]);
	 if (using_ewgts) ewptr = &(graph[bestvtx]->ewgts[1]);
	 for (j=1; j<graph[bestvtx]->nedges; j++){
	    neighbor = *edges++;
	    if (using_ewgts) eweight = .5 + *ewptr++;

	    /* First make sure neighbor is alive. */
	    if (sets[neighbor] >= 0) {
	       group = sets[neighbor];

	       l = 0;
	       for (k=0; k<nsets; k++) {
		  if (k != group) {
		     diff = eweight*(
			hops[k][bestfrom] - hops[group][bestfrom] +
			hops[group][bestto] - hops[k][bestto]);
		     dval = dvals[neighbor][l] + maxdval;
		     movebilist(&listspace[l][neighbor], 
				&buckets[group][k][dval],
				&buckets[group][k][dval+diff]);
		     dvals[neighbor][l] += diff;
		     dval += diff;
		     if (dval > tops[group][k]) tops[group][k] = dval;
		     l++;
		  }
	       }
	    }
         }
         if (DEBUG_KL > 2) {
            pbuckets(buckets, listspace, maxdval, nsets);
         }
      }

      /* Done with a pass; should we actually perform any swaps? */
      bptr = movelist;
      if (bestg > 0 || (bestg != -DOUBLE_MAX && !balanced && enforce_balance)) {
	 improved += bestg;
         for (i=1; i<=beststep; i++) {
	    vtx = ((int) bptr - (int) listspace[0])/sizeof(struct bilist);
	    bestto = (int) bptr->prev;
	    startweight[bestto] += graph[vtx]->vwgt;
	    startweight[-sets[vtx]-1] -= graph[vtx]->vwgt;
	    sets[vtx] = bestto;
	    bptr = bptr->next;
         }

/* Alternately, just make a list of all the vertices that were affected. */
/* Then redo dvals for only these vertices. */
/* This will involve first making such a list,
   Then removing all of their entries in buckets,
   Then calling buckets with this shorter list to work with.
*/
	 deltaminus = deltaplus = 0;
	 for (i=0; i<nsets; i++) {
	    if (startweight[i] - goal[i] > deltaplus) {
	       deltaplus = startweight[i] - goal[i];
	    }
	    else if (goal[i] - startweight[i] > deltaminus) {
	       deltaminus = goal[i] - startweight[i];
	    }
	 }
         balanced = (deltaplus + deltaminus <= vwgt_max);
      }
      else {
	 enforce_balance = TRUE;
	 nbad++;
      }
      if (bestg == -DOUBLE_MAX) enforce_balance_hard = TRUE;

      worth_undoing = (step < undo_frac*nvtxs);
      done = (nbad >= nbadtries && balanced);
      if (!done) {		/* Prepare for next pass. */
	 if (KL_UNDO_LIST && worth_undoing) {
	    /* Make a list of modified vertices for next bucketsort. */
	    /* Also, ensure these vertices are removed from their buckets. */
	    list_length = make_kl_list(graph, movelist, buckets, listspace,
				       sets, nsets, bspace, dvals, maxdval);
	 }
      }
      if (done || !(KL_UNDO_LIST && worth_undoing)) {
         /* Restore set numbers of remaining, altered vertices. */
	 while (bptr != NULL) {
	   vtx = ((int) bptr - (int) listspace[0])/sizeof(struct bilist);
	   sets[vtx] = -sets[vtx] - 1;
	   bptr = bptr->next;
	 }
	 list_length = nvtxs;
      }
   }

   if (DEBUG_KL > 0) {
     {char buf[150]; sprintf(buf,"   KL required %d passes to improve by %d.\n", npass, improved);UserWrite(buf);}
   }

   sfree((char *) loose);
   sfree((char *) locked);
   sfree((char *) startweight);
   sfree((char *) weightsum);
   sfree((char *) bspace);
}
