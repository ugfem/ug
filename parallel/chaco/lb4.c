// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*                                                                          */
/* File:      lb4.c                                                         */
/*                                                                          */
/* Purpose:   clustered partitioning using Chaco's partitioning methods     */
/*                                                                          */
/*            This module combines the clustering technique with different  */
/*            simple and high level partitioning strategies provided by     */
/*            Chaco.                                                        */ 
/*            It divides the elements into clusters according to their      */
/*            position in the element tree and then uses one of Chaco's     */
/*            partitioning schemes to distribute the multigrid onto the     */
/*            processors cluster per cluster. Chaco's partioning schemes    */
/*            are one of 5 spectral methods, a high-end multilevel strategy,*/
/*            inertial partitioning or one of 3 simple method's. The        */
/*            partioning may be optimized by a local KL method.             */
/*            Currently  arbitrary nxm processor configurations are         */
/*            possible.                                                     */
/*                                                                          */
/*                                                                          */
/* Authors:   Stefan Lang                                                   */
/*            Institut fuer Mathematische Maschinen und                     */
/*            Datenverarbeitung III                                         */
/*            Universitaet Erlangen-Nuernberg                               */
/*            Martensstrasse 3                                              */
/*            91058 Erlangen                                                */
/*                                                                          */
/*            Peter Bastian                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                       */
/*            Im Neuenheimer Feld 368                                       */
/*            6900 Heidelberg                                               */
/*            internet: bastian@iwr1.iwr.uni-heidelberg.de                  */
/*                                                                          */
/* History:   3 Jan 94 begin                                                */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"
#include "misc.h"
/* #include "verbose.h" */
#include "lb4.h"
#include "parallel.h"
#include "ppif.h"
/* #include "ugxfer.h" */
/* #include "ugrefine.h" */
#include "debug.h"
#include "refine.h"
#include "rm.h"
#include "cmdint.h"
#include "./main/defs.h"
#include "./main/structs.h"
/* #include "Chaco/code/util/smalloc.h" */

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

#define MAXCLUSTERS         10000       /* max number of clusters           */
#define MAXSETS             16          /* max number of cluster sets       */

#define SMALL_COORD         1.0E-5      /* resolution when comparing COORDs */


#define DESCENDENTS(e)          ((e)->ge.ptmp2)
#define SET_DESCENDENTS(e,n)    (e)->ge.ptmp2 = n
#define MY_CLUSTER(e)           ((CLUSTER *)((e)->ge.ptmp1))
#define SET_MY_CLUSTER(e,p)     (e)->ge.ptmp1 = ((unsigned INT) (p))
#define HAS_CLUSTER(e)          (((CLUSTER *)(e)->ge.ptmp1)!=NULL)
#define MY_CLUSTER_ID(e)        (((CLUSTER *)(e)->ge.ptmp1)->edges[0])
#define SET_PARTITION(e,p)		(PARTITION(e) = p)

/****************************************************************************/
/*                                                                          */
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
/*                                                                          */
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$", UG_RCS_STRING);


/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

static CLUSTER *clusters;                 /* start of cluster array         */
static CLUSTER **sort_clusters;           /* pointers for sorting clusters  */
static INT cluster_set[MAXSETS];          /* starts of cluster sets         */
static INT cluster_cnt[MAXSETS];          /* #clusters in each set          */
static int set_cnt;                       /* number of cluster sets stored  */
static int total_cnt;                     /* total number of clusters stored*/
static int startid;						  /* startid for unique clusternumb */
static INT *load;                         /* total load on all levels&proc!!*/
static INT MarkKeyTop, MarkKeyBottom;     /* mark-keys for mem management   */
static INT quiet;

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* Function:  InitClustering                                                */
/*                                                                          */
/* Purpose:   allocates storage for clusters and inits sizes and sets       */
/*                                                                          */
/* Input:     MULTIGRID *mg: whole data structure                           */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int InitClustering (MULTIGRID *mg, ELEMENT ***elements)
{
	int i;
	INT ne;

	/* allocate storage for clusters in each processor */
	PRINTDEBUG(dddif,1,("%d: InitClustering() GetMem bytes=%d\n",me,MAXCLUSTERS*sizeof(CLUSTER)));
	clusters = GetMemUsingKey(MGHEAP(mg),MAXCLUSTERS*sizeof(CLUSTER),FROM_BOTTOM,MarkKeyBottom);
	if (clusters==NULL) return(1);
	for (i=0; i<MAXSETS; i++)
	{
		cluster_set[i] = 0;
		cluster_cnt[i] = 0;
	}
	total_cnt = 0;
	startid   = 1;

	/* allocate array of cluster pointers for sorting */
	PRINTDEBUG(dddif,1,("%d: InitClustering() GetMem bytes=%d\n",me,MAXCLUSTERS*sizeof(clusters)));
	sort_clusters = (CLUSTER **) GetMemUsingKey(MGHEAP(mg),MAXCLUSTERS*sizeof(clusters),FROM_BOTTOM,MarkKeyBottom);
	if (sort_clusters==NULL) return(1);

	/* allocate memory for array of pointers to elements used in load transfer */
	/* This memory is allocated from top since bottom is released before       */
	/* load transfer.                                                          */
	ne = 0;
	for (i=0; i<=TOPLEVEL(mg); i++)
		if (GRID_ON_LEVEL(mg,i)!=NULL)
			ne += NT(GRID_ON_LEVEL(mg,i));
	PRINTDEBUG(dddif,1,("%d: InitClustering() GetMem nitems=%d bytes=%d\n",me,ne,ne*sizeof(ELEMENT *)));
	*elements = (ELEMENT **) GetMemUsingKey(MGHEAP(mg),ne*sizeof(ELEMENT *),FROM_TOP,MarkKeyTop);
	if (*elements==NULL && ne>0) return(1);

	return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  InitMemLB4                                                    */
/*                                                                          */
/* Purpose:   allocates memory used by lb4 algorithm                        */
/*            i.e. load for each processor per level                        */
/*                                                                          */
/* Input:     MULTIGRID *mg: whole data structure                           */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int InitMemLB1 (MULTIGRID *mg)
{
	int i;

	/* allocate memory for storing load per level in each processor */
	load = GetMem(MGHEAP(mg),procs*MAXLEVEL*sizeof(INT),FROM_BOTTOM);
	if (load==NULL) return(1);
	for (i=0; i<procs*MAXLEVEL; i++) load[i] = 0;

	return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  ComputeDescendents                                            */
/*                                                                          */
/* Purpose:   computes number of descendents in the tree for each elemement */
/*            This is done only locally !                                   */
/*                                                                          */
/* Input:     MULTIGRID *mg: whole data structure                           */
/*                                                                          */
/* Output:    none                                                          */
/*                                                                          */
/****************************************************************************/

static void ComputeDescendents (MULTIGRID *mg)
{
	int k;
	ELEMENT *e,*f;
	GRID *g;

	/* set size word to 1          */
	/* set cluster pointer to NULL */
	/* set destination to me       */
	for (k=TOPLEVEL(mg); k>=0; k--)
	{
		g = GRID_ON_LEVEL(mg,k);
		if (g==NULL) continue;

		for (e=FIRSTELEMENT(g); e!=NULL; e=SUCCE(e))
		{
			SET_MY_CLUSTER(e,NULL);
			SET_PARTITION(e,me);
			if (DDD_InfoPriority(PARHDRE(e))==PrioMaster)
			{
				if (EstimateHere(e)) 
					SET_DESCENDENTS(e,1+
						NSONS_OF_RULE(MARK2RULEADR(e,MARK(e))));
				else
					SET_DESCENDENTS(e,1);
			}
			else
				SET_DESCENDENTS(e,0);
		}
	}

	/* add my size to father's size word */
	for (k=TOPLEVEL(mg); k>0; k--)
	{
		g = GRID_ON_LEVEL(mg,k);
		if (g==NULL) continue;

		for (e=FIRSTELEMENT(g); e!=NULL; e=SUCCE(e))
		{
			f = EFATHER(e);
			if (f==NULL) continue;

			SET_DESCENDENTS(f,DESCENDENTS(f)+DESCENDENTS(e));
		}
	}
}

/****************************************************************************/
/*                                                                          */
/* Function:  new_cluster                                                   */
/*                                                                          */
/* Purpose:   get a pointer to a new cluster in set 0 (my set)              */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    pointer to new cluster                                        */
/*            NULL if not enough memory                                     */
/*                                                                          */
/****************************************************************************/

static CLUSTER *new_cluster (void)
{
	CLUSTER *c;
	int i;

	if (total_cnt>=MAXCLUSTERS)
	{
		printf( PFMT "new_cluster(): increase MAXCLUSTERS in lb4.c\n",me);
		PrintErrorMessage('E',"new_cluster","increase MAXCLUSTERS in lb4.c\n");
		return(NULL);
	}
	cluster_cnt[0]++;
	c = clusters+(total_cnt);
	c->source = me;
	c->destination = 0;
	c->minlevel = c->depth = 0;
	for (i=0; i<MAXDEPTH; i++) c->level_size[i] = 0;
	c->sx = c->sy = 0.0;
	c->edges[0] = total_cnt++;
	c->root_element = NULL;
	c->size = 0;
	c->nedges = 0;
	for (i=1; i<=MAX_SIDES_OF_ELEM; i++) 
	{
		c->edges[i] = 0;
		/* c->ewgts[i] = 0.0; */
	}
	
	return(c);
}


/****************************************************************************/
/*                                                                          */
/* Function:  Clustering                                                    */
/*                                                                          */
/* Purpose:   assign each of my master elements to a cluster.               */
/*                                                                          */
/* Input:     MULTIGRID *mg: whole data structure                           */
/*            int minlevel          lowest level to balance                 */
/*            int cluster_depth     min cluster depth                       */
/*            int threshold         min size to start a new cluster         */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int check_condition_from_refinement (ELEMENT *e)
{
	int i;
	
	/* this function returns TRUE if e must be stored in the  */
	/* same processor as its father.                          */

	/* level 0 elements can be everywhere */
	if (LEVEL(e)==0) return(0);

	/* copies and irregular elements must be in the fathers proc */
	if ((ECLASS(e)==YELLOW_CLASS)||(ECLASS(e)==GREEN_CLASS)) return(1);

	/* now e is a regular element, if it has no sons it must     */
	/* be with its father.                                       */
	if (NSONS(e)==0) return(1); 

	/* if one of the sons is copy or irregular, it must be       */
	/* with its father.                                          */
    {
        ELEMENT *SonList[MAX_SONS];

        if (GetSons(e,SonList)!= GM_OK) assert(0);
        for (i=0; i<NSONS(e); i++)
        {
			if (SonList[i]!=NULL)
				if (ECLASS(SonList[i])!=RED_CLASS) return(1);
			else
				break;
        }
    }

	return(0);
} 


static int Clustering (MULTIGRID *mg, int minlevel, int cluster_depth, int threshold)
{
	int k,i,j,l;
	ELEMENT *e,*s;
	GRID *g;
	CLUSTER *c;

	/* compute number of descendents of each element */
	ComputeDescendents(mg);

	/* now cluster bottom up, using the element tree */
	for (k=minlevel; k<=TOPLEVEL(mg); k++)
	{
		g = GRID_ON_LEVEL(mg,k);
		if (g==NULL) continue;

		for (e=FIRSTELEMENT(g); e!=NULL; e=SUCCE(e))
		{
			if (!(DDD_InfoPriority(PARHDRE(e))==PrioMaster)) continue;

			/* if e has no cluster allocate a new one */
			if (MY_CLUSTER(e)==NULL)
			{
				/* get a new cluster for e      */
				/* e is root element of cluster */
				c = new_cluster();
				if (c==NULL) return(1);
				SET_MY_CLUSTER(e,c);
				c->source = me;
				c->minlevel = LEVEL(e);
				c->sx = c->sy = 0.0;
				#ifdef __THREEDIM__
				c->sz = 0.0;
				#endif
				for (i=0; i<CORNERS_OF_ELEM(e); i++)
				{
					c->sx += XC(MYVERTEX(CORNER(e,i)));
					c->sy += YC(MYVERTEX(CORNER(e,i)));
					#ifdef __THREEDIM__
					c->sz += ZC(MYVERTEX(CORNER(e,i)));
					#endif
				}
				c->sx /= ((COORD)(CORNERS_OF_ELEM(e)));
				c->sy /= ((COORD)(CORNERS_OF_ELEM(e)));
				#ifdef __THREEDIM__
				c->sz /= ((COORD)(CORNERS_OF_ELEM(e)));
				#endif
				c->root_element = e;
			}

			/* put e in its cluster */
			c = MY_CLUSTER(e);
			c->depth = MAX(c->depth,LEVEL(e)-c->minlevel);
			if (c->depth>=MAXDEPTH) return(2);
			(c->size)++;
			(c->level_size[LEVEL(e)-c->minlevel])++;

			/* put future sons of e AFTER REFINEMENT in the same cluster */
			if ((EstimateHere(e)) && NSONS_OF_RULE(MARK2RULEADR(e,MARK(e)))>0)
			{
				c->depth = MAX(c->depth,LEVEL(e)+1-c->minlevel);
				if (c->depth>=MAXDEPTH) return(2);
/* TODO: j must be estimated for new green closure in 3D */
#ifdef __THREEDIM__
if (TAG(e)!=TETRAHEDRON && MARKCLASS(e)==GREEN_CLASS) assert(0);
#endif
				j = NSONS_OF_RULE(MARK2RULEADR(e,MARK(e)));
				c->size += j;
				c->level_size[LEVEL(e)+1-c->minlevel] += j;
			}
			
			/* assign the sons of e to clusters */
			{
				ELEMENT *SonList[MAX_SONS];

				if (GetSons(e,SonList)!= GM_OK) assert(0);

				for (i=0; i<NSONS(e); i++)
				{
					s = SonList[i];
					if (s==NULL) break;
					if (!(DDD_InfoPriority(PARHDRE(s))==PrioMaster)) continue;

					/* if condition implied by refinement algorithm is true */
					/* put s in same cluster as e.                          */
					if (check_condition_from_refinement(s))
					{
						SET_MY_CLUSTER(s,c);
						continue;
					}

					/* if max depth of a cluster is not reached, */
					/* put s in same cluster as e.               */
					if (LEVEL(s)-(c->minlevel)<=cluster_depth)
					{
						SET_MY_CLUSTER(s,c);
						continue;
					}

					/* if son has not enough descendents to build its own cluster,*/
					/* put s in same cluster as e.                                */
					if (DESCENDENTS(s)<threshold)
					{
						SET_MY_CLUSTER(s,c);
						continue;
					}

					/* if none of the above is true, the cluster pointer in s  */
					/* remains NULL. Therefore a new cluster will be allocated */
					/* for s.                                                  */
				}
			}
		}
	}
	return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  InitConcentrateClusters                                       */
/*                                                                          */
/* Purpose:   initialize transfer of all clusters to the master processor   */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int InitConcentrateClusters (void)
{
	int l,n;

	/* set 0 is mine, start with set  #1 to receive */
	set_cnt = 1;

	/* get clusternumbers from down tree */
	for (l=degree-1; l>=0; l--)
	{
		/* get number of clusters to receive */
		GetConcentrate(l,&n,sizeof(int));

		/* allocate cluster memory */
		cluster_set[set_cnt] = total_cnt;
		total_cnt += n;
		cluster_cnt[set_cnt] = n;
		set_cnt++;
	}

	/* send clusternumbers uptree */
	Concentrate(&total_cnt,sizeof(int));

	return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  BroadCastStartIds                                             */
/*                                                                          */
/* Purpose:   send startids to processors for global unique cluster         */
/*            numbering                                                     */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int BroadCastStartIds (void)
{
    int l,n,start;
    CLUSTER *c;

    /* get startid from uptree */
	if (me!=master) RecvSync(uptree,&startid,sizeof(int));

    /* send startids down tree */
	start = startid+total_cnt;
    for (l=0; l<degree; l++)
    {
		start -= cluster_cnt[set_cnt-l-1];
        SendSync(downtree[l],&start,sizeof(int));
	}

    return(0);
}


/****************************************************************************/
/*                                                                          */
/* Function:  ExchangeNeighborInfo                                          */
/*                                                                          */
/* Purpose:   exchange neighbor graph numbers of interface lists            */
/*            between processors                                            */
/*                                                                          */
/* Input:     MULTIGRID *mg: grid level                                     */
/*            int itemSize: #bytes to buffer for each interface object      */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int GatherClusterNumber (DDD_OBJ theCoupling, void *data)
{
	ELEMENT *theElement = (ELEMENT *) theCoupling;
	CLUSTER *cptr;

	cptr = MY_CLUSTER(theElement);

	/* if Element is no master_copy, id is assumed to be zero */
	if (DDD_InfoPriority(PARHDRE(theElement))==PrioMaster && cptr!=NULL)
	{

		/* Is theElement rootelement? */
		if (theElement == cptr->root_element)

			(*(int *)data) = MY_CLUSTER_ID(theElement);

		/* cluster has no neighbors on this level */
		else
			(*(int *)data) = 0;
	}
	else 
		(*(int *)data) = 0;

	return(0);
}

static int ScatterClusterNumber (DDD_OBJ theCoupling, void *data)
{
	int i=1;
    ELEMENT *theElement = (ELEMENT *) theCoupling;
	CLUSTER *cptr;
	char buffer[100]; 

	cptr = (CLUSTER *) MY_CLUSTER(theElement);

	/* scatter id only if Element is not a master_copy */
	if (DDD_InfoPriority(PARHDRE(theElement))==PrioMaster && cptr!=NULL)
	{
		if (cptr->root_element!=theElement) return(0);
		while (cptr->edges[i]!=0 && i<=SIDES_OF_ELEM(theElement)) i++;
		if (i>SIDES_OF_ELEM(theElement)) 
		{
			sprintf(buffer,"ScatterClusterNumber: neighbors=%d, i=%d", TAG(theElement),i);  
			UserWrite(buffer);
			return(1);
		}
		if ((*(int *)data) != 0)
		{
			cptr->edges[i] = (*(int *)data);
			cptr->nedges++;
		}
	}

	return(0);
}

static int ExchangeNeighborInfo(MULTIGRID* mg)
{
	char buffer[80];
	int error=0;

	DDD_IFExchange(ElementSymmIF, sizeof(int), GatherClusterNumber,
										ScatterClusterNumber);

	if (error>0) 
	{
		sprintf(buffer,"communication in ExchangeNeighborInfo failed, error=%d\n",
				error);
		PrintErrorMessage('E',99,buffer);
		return(2);
	}


	return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  ComputeGraphInfo                                              */
/*                                                                          */
/* Purpose:   number each cluster with its unique id and compute            */
/*            neighbor idinformation as input for Chaco                     */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int ComputeGraphInfo (MULTIGRID *mg)
{
	int i,j,l,n,error;
	ELEMENT *e,*theElement;
	CLUSTER *cptr,*c;
	char buffer[100];

	/* give processors global unique ids for cluster numbering */
	error = BroadCastStartIds();

	/* number each locally created Cluster */ 
	for (i=0; i<cluster_cnt[0]; i++)
	{
		(clusters+i)->edges[0] += startid;
	}

	/* determine the neighbor ids of each cluster */
	for (i=0; i<cluster_cnt[0]; i++)
	{
		cptr = clusters+i;
		theElement = cptr->root_element;
		n = j = 1;

		/* set index to first free entry */
		while (cptr->edges[j] != 0) { j++; n++;}
		if (j>1)
		{
			sprintf(buffer,"ComputeGraphInfo: already %d neighbors in list\n", j-1);
			UserWrite(buffer);
			return(1);
		} 

		for (l=0; l<SIDES_OF_ELEM(theElement); l++)
		{
			/* more sides than possible */
			if (j>MAX_SIDES_OF_ELEM) return(1);

			/* look for root element of neighbor cluster */
			e = NBELEM(theElement,l);

			/* Is e neighbor? */
			if (e != NULL) 
			{
				/* if element has cluster, it has id too! */
				if (HAS_CLUSTER(e) && DDD_InfoPriority(PARHDRE(e))==PrioMaster)
				{
					/* cluster of e is only neighbor, if e */
					/* is root element of this cluster!    */
					c = MY_CLUSTER(e);
					if (c->root_element == e)
					{
						(cptr->edges)[j++] = MY_CLUSTER_ID(e);
						n++; 
					}
				}
			}
		}
		/* set number of neighbors */
		cptr->nedges = n;
	}

	/* exchange cluster ids of boundary clusters between processors */
	error = ExchangeNeighborInfo(mg);

}		

/****************************************************************************/
/*                                                                          */
/* Function:  ConcentrateClusters                                           */
/*                                                                          */
/* Purpose:   transfer all clusters to the master processor                 */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int ConcentrateClusters (void)
{
	int l,n;
	CLUSTER *c;

	/* set 0 is mine, start with set  #1 to receive */
	set_cnt = 1;
	c = clusters+cluster_cnt[0];

	/* get clusters from down tree */
	for (l=degree-1; l>=0; l--)
	{
		/* receive clusters */
		n = cluster_cnt[set_cnt];
		GetConcentrate(l,c,n*sizeof(CLUSTER));
		c += cluster_cnt[set_cnt++];
	}

	/* send clusters uptree */
	Concentrate(clusters,total_cnt*sizeof(CLUSTER));

	return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  balance_ccptm                                                 */
/*                                                                          */
/* Purpose:   balance all clusters                                          */
/*                                                                          */
/* Input:     int Const  min #elements per proc used for configuration size */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if an error occured                                        */
/*                                                                          */
/****************************************************************************/

static INT sum_load (int l, int x0, int y0, int x, int y)
{
	int i,j;
	INT sum;

	sum = 0;
	for (i=x0; i<x0+x; i++)
		for (j=y0; j<y0+y; j++)
			sum += load[(j*DimX+i)*MAXLEVEL+l];
	return(sum);
}


static void fit_configuration (int *p, int *x, int *y)
{
	int i,xx,yy;

	xx = yy = 1;
	*x = xx; *y = yy;

	for (i=0; i<DimX+DimY-2; i++)
	{
		if (((xx<=yy)||(yy==DimY))&&(xx<DimX))
		{
			xx++;
		}
		else
		{
			if (yy<DimY) 
			{
				yy++;
			}
		}
		if (xx*yy<=*p)
		{
			*x = xx;
			*y = yy;
		}
		else
		{ /* now: x*y <= p */ 
			*p = xx*yy;
			break;
		}
	}
}

static double compute_goal (int level, double average, double N, double **goal,
                            int p, int x0, int y0, int x, int y)
{
	extern Heap   *heap;     /* pointer to heap of multigrid */
	extern double *MEM_OK;   /* variable for memory overflow exeception */
	int i,j,pgoal;
	double *gptr,balance,avg,max,ngoal;
	char buffer[100];

	*goal = (MEM_OK = (double *) smalloc((unsigned) p*sizeof(double));
	if (!MEM_OK) return;
	gptr = *goal;
	avg = average;
	ngoal = 0;

	do
	{
		balance = N;
		pgoal = 0;
		max = 0;
		/* compute total number of elements of processor get */
		/* clusters for load balance                         */
		for (i=x0; i<x0+x; i++)
		{
			for (j=y0; j<y0+y; j++)
			{
				if (avg > load[(j*DimX+i)*MAXLEVEL+level])
				{
					balance += load[(j*DimX+i)*MAXLEVEL+level];
					if (max < load[(j*DimX+i)*MAXLEVEL+level])
						max = load[(j*DimX+i)*MAXLEVEL+level];
					pgoal++;
				}
			}
		}

		/* mean balance value which can be gained with balancing  */
		/* in current load situation. balance can be smaller than */
		/* average, if there are too few elements in remaining    */
		/* clusters on this level.                                */                             
		balance /= pgoal;
		avg = balance;
	}
	while (max > avg);
		
	
	if (balance < average-0.05)
	{
		sprintf(buffer, "warning: cannot gain load balance: average=%f, balance=%f, n=%d\n",
		average, balance, pgoal);
		UserWrite(buffer);
	}

	/* compute goal in current situation */
	for (i=x0; i<x0+x; i++)
	{
		for (j=y0; j<y0+y; j++)
		{
			if (balance > load[(j*DimX+i)*MAXLEVEL+level])
			{
				*(gptr) = MAX(balance-(double)load[(j*DimX+i)*MAXLEVEL+level],0);
				ngoal += *(gptr++);
			}
			else 
				*(gptr++) = 0;
		}
	}
	if (ngoal>N+0.05 || ngoal<N-0.05)
	{
		sprintf(buffer, "fatal: sum of goal array differs from sum of clusterelements goal=%f, N=%f\n",
		        ngoal, N);
		UserWrite(buffer);
	}
	return(balance); 
}

static void check_assign (CLUSTER **clusters, double *goal, int dimx, int dimy, 
                          short *assign, int ncluster)
{
	extern Heap   *heap;     /* pointer to heap of multigrid */
	extern double *MEM_OK;   /* variable for memory overflow exeception */
	double total_goal;
	double g;
	int total_assigned;
	int *assigned;
	int a,procs,x,y;
	int i,j;
	char buffer[100];

	procs = dimx*dimy;
	assigned = (int *) (MEM_OK = smalloc( (unsigned) procs *sizeof(int));
	if (!MEM_OK) return;
	total_goal = 0.0;
	total_assigned = 0;

	for (i=0; i<procs; i++) assigned[i] = 0;
	for (i=1; i<=ncluster; i++) 
	{
		if (assign[i]<0 || assign[i]>=procs)
		{
			sprintf(buffer,"assignment wrong: assign=%d, procs=%d\n",
			        assign[i], procs);
			UserWrite(buffer);
		}
		assigned[assign[i]] += clusters[i-1]->level_size[clusters[i-1]->depth];
		x = assign[i]%dimx;
		y = assign[i]/dimx;
		for (j=0; j<=clusters[i-1]->depth; j++)
			load[(y*DimX+x)*MAXLEVEL+clusters[i-1]->minlevel+j]+=clusters[i-1]->level_size[j];
	}

	for (i=0; i<procs; i++)
	{
		g = floor(goal[i]);	
		a = assigned[i];
if (!quiet)
		if (g+1<a || g>assigned[i])
		{
			sprintf(buffer,"goal for processor %d failed: goal=%f, assigned=%d\n",i,
			        goal[i], a);
			UserWrite(buffer);
		}
		total_goal += goal[i];
		total_assigned += a;
	}

	g = floor(total_goal);

	if (g+1<total_assigned || g>total_assigned)
	{
		sprintf(buffer,"total goal failed: total_goal=%f, total_assigned=%d\n",
		       total_goal, total_assigned);
		UserWrite(buffer);
	}

	sfree((char *) assigned);
}
		
static int compare_maxlevel (const void *e1, const void *e2)
{
	CLUSTER *c1,*c2;
	int m1,m2;
	
	/* the two clusters to compare */
	c1 = *((CLUSTER **) e1);
	c2 = *((CLUSTER **) e2);
	
	m1 = c1->minlevel;
	m2 = c2->minlevel;

	if (m1<m2) return( 1);
	if (m1>m2) return(-1);

	m1 += c1->depth;
	m2 += c2->depth;

	if (m1<m2) return( 1);
	if (m1>m2) return(-1);

	if (c1->source<c2->source) return( 1);
	if (c1->source>c2->source) return(-1);

	return(0);
}

static int ArrayToProcid (int aid, int dimx, int dimy)
{
	int x;

	x = aid/dimx*DimX;
	return(x+(aid-x));
}

static INT max_load (int l, int x0, int y0, int x, int y)
{
    int i,j;
    INT max;

    max = 0;
    for (i=x0; i<x0+x; i++)
        for (j=y0; j<y0+y; j++)
                max = MAX(load[(j*DimX+i)*MAXLEVEL+l],max);
    return(max);
}

static void reset_timers (void) 
{
	extern double input_time, partition_time, reformat_time;
	extern double check_input_time, count_time, print_assign_time;
	extern double coarsen_time, match_time, make_cgraph_time;
	extern double make_cewgts_time, adj_cewgts_time;
	extern double lanczos_time, splarax_time, orthog_time;
	extern double ql_time, tevec_time, ritz_time;
	extern double evec_time, check_time, blas_time;
	extern double init_time, scan_time, debug_time;
	extern double probe_time, pause_time;
	extern double rqi_symmlq_time, refine_time;
	extern double kl_total_time, kl_init_time, nway_kl_time;
	extern double inertial_time, rcb_time, inertial_axis_time, median_time;

	input_time = partition_time = reformat_time = 0;
	check_input_time = count_time = print_assign_time = 0;

	coarsen_time = match_time = make_cgraph_time = 0;
	make_cewgts_time = adj_cewgts_time = 0;

	lanczos_time = splarax_time = orthog_time = 0;
	ql_time = tevec_time = ritz_time = 0;
	evec_time = check_time = blas_time = 0;
	init_time = scan_time = debug_time = 0;
	probe_time = pause_time = 0;

	rqi_symmlq_time = refine_time = 0;

	kl_total_time = kl_init_time = nway_kl_time = 0;

	inertial_time = rcb_time = inertial_axis_time = median_time = 0;
}

static int balance_ccptm (MULTIGRID *mg, int Const, int strategy, int eigen, 
                          int loc, int dims, int weights, int coarse, int mode)
{
	extern Heap   *heap;     /* pointer to heap of multigrid */
	extern double *MEM_OK;   /* variable for memory overflow exeception */
	int i,j,maxlevel,minlevel,first,last,x,y,p;
	INT nc,N,M,ne;
	unsigned INT NC;
	unsigned long max_size;
	void *heap_memory;
	CLUSTER *c;
	double average,balance,maximum,cmaximum;
	char buffer[100];
	short *assign,*assignment;
	double *goal; 
	int *glob2loc;
	
	/* this is done only at the master processor */
	if (me!=master) return(0);
	
	/* fill cluster pointer array */
	for (i=0; i<total_cnt; i++) sort_clusters[i] = clusters+i;
	
	/* sort clusters by max level */
	qsort((void *)sort_clusters,total_cnt,sizeof(clusters),compare_maxlevel);

	/* create new heap of type GENERAL_HEAP to substitute Chaco's */
	/* smalloc/sfree functions by GetMem/DisposeMem of ugp        */
	max_size = HeapSize(MGHEAP(mg)) - HeapUsed(MGHEAP(mg)) - sizeof(BLOCK);
	max_size = max_size - ((ALIGNMENT-((max_size)&(ALIGNMENT-1)))&(ALIGNMENT-1));
	heap_memory = GetMem(MGHEAP(mg), max_size, FROM_BOTTOM);
	heap = NewHeap(GENERAL_HEAP, max_size, (void *)heap_memory); 

	/* allocate map for glob to loc numbering */
	glob2loc = (int *) (MEM_OK = smalloc((unsigned) (total_cnt+1)*sizeof(int));
	if (!MEM_OK) return;
	
	for (i=0; i<=total_cnt; i++) glob2loc[i] = 0;

	/* reset timing variables */
	reset_timers();

	/* process levels top to bottom */
	i = 0;
	while (i<total_cnt)
	{
		/* next maxlevel*/
		last = first = i++;
		c = sort_clusters[first];
		maxlevel = c->minlevel + c->depth;
		minlevel = c->minlevel;
		
		/* find last cluster with same maxlevel */
		while (i<total_cnt)
		{
			c = sort_clusters[i];
			if (c->minlevel+c->depth!=maxlevel||c->minlevel!=minlevel) break;
			last = i++;
		}
		
		/* compute number of clusters */
		nc = last-first+1;
		
		/* compute total number of elements on level maxlevel */
		M = sum_load(maxlevel,0,0,DimX,DimY);
		N = 0;
		for (j=first; j<=last; j++)
		{
			c = sort_clusters[j];
			N += c->level_size[maxlevel-c->minlevel];
		}
		
		/* compute number of processors to use */
		p = MAX(MIN((N+M)/Const,procs),1);

		/* fit array configuration to processor number */
		fit_configuration(&p,&x,&y);

		/* compute average load per processor desired */
		average = ((double)(N+M))/((double)p);

		/* compute number of elements that has to be assigned */
		/* to each processor                                  */
		balance = compute_goal(maxlevel,average,N,&goal,p,0,0,x,y);

		/* allocate space for assignments */
		assign = (short *) (MEM_OK = smalloc((unsigned) (nc+1)*sizeof(short));
		if (!MEM_OK) return;
		assignment = assign;

		for (j=0; j<=total_cnt; j++)
			if (glob2loc[j]!=0)
			{
				sprintf(buffer,"balance_ccptm: maxlevel %d, minlevel %d, glob2loc %d=%d\n",
				        maxlevel, minlevel, j, glob2loc[j]);
				UserWrite(buffer);
				glob2loc[j]=0;
			}

		sprintf(buffer,"Balancing: %d minlevel, %d maxlevel, %d vertices\n", minlevel, maxlevel, nc);
		UserWrite(buffer);

		/* balance clusters with desired partitioning strategy */
		interface(sort_clusters+first,nc,assign,goal,p,strategy,eigen,loc,
		          dims,weights,coarse,mode,glob2loc,x,y); 

		/* check for memory error in interface() */
		if (!MEM_OK) return;

		/* check whether partitioning succeeded */
		check_assign(sort_clusters+first, goal, x, y, assign, nc);
		

		/* convert tree id's in clusters to processor id's */
		/* this is done here, since the configuration size */
		/* is needed !                                     */
		for (j=first; j<=last; j++)
		{
			c = sort_clusters[j];
			assign++;
			c->destination = ArrayToProcid(*assign,x,y);
		}

		/* compute maximum element number of processors in any level */
		maximum = (double) max_load(maxlevel,0,0,DimX,DimY);

		/* print statistics */
		/* if (quiet_level<1) */
		if (1)
		{
			sprintf(buffer,"lb4: level=%2d p=%3d nc=%5d avg=%10.3lg max=%10.3lg imbal=%10.3lg\n",
					maxlevel,p,nc,average,maximum,ABS(maximum/average-1.0)*100.0);
			UserWrite(buffer);
		}

	/* free space allocated in this loop */
	sfree((char *) goal);
	sfree((char *) assignment);
	}

	sfree((char *) glob2loc);
	return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  OptimizeMapping                                               */
/*                                                                          */
/* Purpose:   optimizes the mapping of the clusters to the processors       */
/*            using nway-KL. Only all clusters may be changed between       */
/*            two distinct processors to minimize the array hop costs.      */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int OptimizeMapping(void)
{
	extern Heap   *heap;     /* pointer to heap of multigrid */
	extern double *MEM_OK;
	struct vtx_data **graph;
	struct vtx_data *links;
	CLUSTER c, *nbcl;
	double *goal;
	float *eweights;
	int **p2p;
	int *mtx, *p2pmtx;
	int *edges;
	int i,j,k;
	int max_deg,sum,nedges,nbid,nprocs;
	int **hop_mtx;
	short *assign;
	short *map;
	SHORT nbdest;
	
	/* this is only done at the master processor */
	if (me!=master) return;

	nedges = nprocs = 0;


	/* determine neighborship relations between processors */
	mtx = (int *) (MEM_OK = smalloc((unsigned)(DimX*DimY*DimX*DimY)*
	                                  sizeof(int));
	p2pmtx = mtx;
	if (!MEM_OK) return;
	p2p = (int **) (MEM_OK = smalloc((unsigned)(DimX*DimY)*
	                                  sizeof(int));
	if (!MEM_OK) return;
	for (i=0; i<DimX*DimY; i++)
	{
		p2p[i] = mtx;
		mtx += DimX*DimY;
	}
		
		

	for (i=0; i<DimX*DimY; i++)
		for (j=0; j<DimX*DimY; j++)
			p2p[i][j]=0;

	for (i=0; i<total_cnt; i++)
	{
		c = clusters[i];
		for (j=0; j<MAX_SIDES_OF_ELEM; j++)
		{
			nbid = c.edges[j];
			if (nbid!=0)
			{
				nbdest = clusters[nbid-1].destination;
				if (p2p[c.destination][nbdest]==0) nedges++;
				p2p[c.destination][nbdest]++;
			}
		}
	}

	for (i=0; i<DimX*DimY; i++)
		if (p2p[i][i]!=0) nprocs++;

	/* generate some space */
	assign = (short *) (MEM_OK = smalloc((unsigned)(nprocs)*sizeof(short));
	if (!MEM_OK) return;
	map = (short *) (MEM_OK = smalloc((unsigned)(DimX*DimY)*sizeof(short));
	if (!MEM_OK) return;

	/* build up graph with procssor neighborship relations */
	graph = (struct vtx_data **) (MEM_OK = smalloc((unsigned)(nprocs+1)
	                                       *sizeof(struct vtx_data *));
	if (!MEM_OK) return;
	links = (struct vtx_data *) (MEM_OK = smalloc((unsigned)(nprocs)
	                                      *sizeof(struct vtx_data));
	if (!MEM_OK) return;

	for (i=1; i<=nprocs; i++)
	{
		graph[i] = links++;
	}

	edges = (int *) (MEM_OK = smalloc((unsigned)nedges*sizeof(int));
	if (!MEM_OK) return;

	eweights = (float *) (MEM_OK = smalloc((unsigned)nedges*sizeof(float));
	if (!MEM_OK) return;

	/* Now fill in all the data fields. */
	for (i=0,k=1; i<DimX*DimY; i++) 
	{
		if (p2p[i][i]==0) continue;

		graph[k]->vwgt = 1;
		graph[k]->edges = edges;
		*edges = k;
		map[i] = k;
		assign[k] = i;
		edges++;
		nedges = 1;
		for (j=0; j<DimX*DimY; j++)
		{
			if ((p2p[i][j]!=0) && (i!=j))
			{
				*edges = j+1;
				edges++;
				nedges++;
			}
		}
		graph[k]->nedges = nedges;
		graph[k]->ewgts = eweights;
		eweights++;
		sum = 0;
		for (j=0; j<DimX*DimY; j++)
		{
			if ((p2p[i][j]!=0) && (i!=j))
			{
				sum += p2p[i][j];
				*eweights++ = p2p[i][j];
			}
		}
		graph[k++]->ewgts[0] = -sum;
	}

	/* generate goal array */
	goal = (double *) (MEM_OK = smalloc((unsigned) (nprocs+1)*sizeof(double));
	if (!MEM_OK) return;
	for (i=1; i<=nprocs; i++) goal[i] = 1;

	
	/* generate hop matrix */
	mtx = (int *) (MEM_OK = smalloc((unsigned) (DimX*DimY*DimX*DimY)*
	                                      sizeof(int));
	if (!MEM_OK) return;
	hop_mtx = (int **) (MEM_OK = smalloc((unsigned) (DimX*DimY)*
	                                      sizeof(int));
	if (!MEM_OK) return;

	for (i=0; i<DimX*DimY; i++)
	{
		hop_mtx[i] = mtx;
		mtx += DimX*DimY;
	}

	for (i=0; i<nprocs; i++) 
	{
		for (j=0; j<nprocs; j++) 
		{
			if (i==j)
				hop_mtx[i][j] = 0;
			else
				hop_mtx[i][j] = absval(assign[i]%DimX - assign[j]%DimX) + 
				                absval(assign[i]/DimX - assign[j]/DimX);
		}
	}

	/* optimize array hop costs using KL */
	max_deg = find_maxdeg(graph, nprocs);
	klspiff(graph, nprocs, assign, nprocs, hop_mtx, goal, 1, max_deg);

	/* if something has changed in assignment , set all cluster destinations */
	/* to new processor id.                                                  */ 
	for  (i=0,k=1; i<DimX*DimY; i++)
	{
		/* has processor clusters? */
		if (p2p[i][i]!=0) 
			if (assign[k]!=map[k])
			{
				for (j=0; j<total_cnt; j++)
				{
					c = clusters[i];
					if (c.destination == map[k])
						c.destination = assign[k];
				}
				k++;
			}
	}
	sfree(hop_mtx);
	sfree(mtx);
	sfree(goal);
	free_graph(graph);
	sfree(map);
	sfree(assign);
	sfree(p2p);
	sfree(p2pmtx);
}
		
		

/****************************************************************************/
/*                                                                          */
/* Function:  BroadcastDestinations                                         */
/*                                                                          */
/* Purpose:   broadcast new destinations to all processors                  */
/*                                                                          */
/* Input:     none                                                          */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if not enough memory                                       */
/*                                                                          */
/****************************************************************************/

static int BroadcastDestinations (void)
{
	int l,n;
	CLUSTER *c;

	/* get clusters from uptree */
	if (me!=master) RecvSync(uptree,clusters,total_cnt*sizeof(CLUSTER));

	/* send clusters down tree */
	for (l=0; l<degree; l++)
	{
		c = clusters+cluster_set[set_cnt-l-1];
		n = cluster_cnt[set_cnt-l-1];
		SendSync(downtree[l],c,n*sizeof(CLUSTER));
	}

	return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  SetDestinations                                               */
/*                                                                          */
/* Purpose:   writes destination in each element and fills array of pointers*/
/*            to elements to be transfered. The order in this array is      */
/*            cluster per cluster to minimize intermediate interfaces       */
/*            during iterative load transfer.                               */
/*                                                                          */
/* Input:     MULTIGRID *mg         whole data structure                    */
/*            int minlevel          lowest level to balance                 */
/*            ELEMENT **elements    pointer array to be filled              */
/*            INT *ne               #entries used in pointer array          */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if an error occured                                        */
/*                                                                          */
/****************************************************************************/

static void recursive_set_dest (ELEMENT *e, CLUSTER *c, ELEMENT **elements, 
                                INT *ne)
{
	int i;
	
	/* if e is not in cluster c then return */
	if (MY_CLUSTER(e)!=c) return;
	
	/* now e is in cluster c, if dest!=me put it in elements */
	if (c->destination!=me)	elements[(*ne)++] = e;

	/* set destination in e */
	SET_PARTITION(e,c->destination);
	
	/* delete cluster pointer to indicate that e has been processed */
	SET_MY_CLUSTER(e,NULL);
	
	/* process sons of e recursively */
	{
		ELEMENT *SonList[MAX_SONS];

		if (GetSons(e,SonList) != GM_OK) assert(0);

		for (i=0; i<NSONS(e); i++)
		{
			ELEMENT *s = SonList[i];
			
			if (s!=NULL)
				recursive_set_dest(s,c,elements,ne);
			else
				break;
		}
	}
}
	
static int SetDestinations (MULTIGRID *mg, int minlevel, ELEMENT **elements, 
                            INT *ne)
{
	int k;
	GRID *g;
	ELEMENT *e;
	CLUSTER *c;
	
	/* reset number of elements */
	*ne = 0;
	
	/* find destination for all elements */
	for (k=minlevel; k<=TOPLEVEL(mg); k++)
	{
		g = GRID_ON_LEVEL(mg,k);
		if (g==NULL) continue;

		for (e=FIRSTELEMENT(g); e!=NULL; e=SUCCE(e))
		{
			/* process only master elements */
			if (!(DDD_InfoPriority(PARHDRE(e))==PrioMaster)) continue;
			
			/* elements with NULL cluster pointer have been processed already */
			if (MY_CLUSTER(e)==NULL) continue;
			
			/* set destinations recursively */
			c = MY_CLUSTER(e);
			recursive_set_dest(e,c,elements,ne);
		}
	}

	return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  BalanceCCPTM                                                  */
/*                                                                          */
/* Purpose:   Load Balancing Call for partitioning with one of Chaco's      */
/*            partitioning method's - Multiplicative Version                */
/*                                                                          */
/* Input:     MULTIGRID *mg:        whole data structure                    */
/*            int minlevel          lowest level to balance                 */
/*            int cluster_depth     min cluster depth                       */
/*            int threshold         min size to start a new cluster         */
/*            int Const             min #elements per proc used for configur*/
/*            int element_limit     for load transfer                       */
/*            int channel_limit     for load transfer                       */
/*            int iter              for load transfer                       */
/*                                                                          */
/* Output:    0: if OK                                                      */
/*            >0 if an error occured                                        */
/*                                                                          */
/****************************************************************************/

int Balance_CCPTM (MULTIGRID *mg,                                 /* data    */
                  int minlevel, int cluster_depth, int threshold,/* clusters*/
				  int Const,                                     /* balance */
				  int element_limit, int channel_limit,          /* transfer*/
				  int strategy,                                  /* strategy*/
				  int eigen,                                     /* eigensol*/
				  int loc,                                       /* KL?     */
				  int dims,                                      /*dimension*/
				  int weights,                                   /* weights?*/
				  int coarse,                                    /* #coarse */
				  int mode,                                      /* mode?   */
				  int iter)                                      /* transfer*/
{
	extern double *MEM_OK;   /* variable for memory overflow exeception */
	int error,l;
	ELEMENT **elements;
	INT ne,nc,i;
	DOUBLE Begin,End;
	char buf[60];

	/* quiet */
	quiet = (mode==512);

	Begin = CURRENT_TIME;
	
	/* mark heap, clusters are allocated from bottom, */
	/* element array for load transfer is allocated   */
	/* from top.                                      */
	Mark(MGHEAP(mg),FROM_BOTTOM, &MarkKeyBottom);
	Mark(MGHEAP(mg),FROM_TOP, &MarkKeyTop);

	/* reset error flag and memory error flag*/
	error = 0;
	MEM_OK = (double *)0x1L;

	/* allocate all memory we might need now */
	error = InitClustering(mg,&elements);
	if (error>0) goto stage1;
	error = InitMemLB1(mg);
	if (error>0) goto stage1;

	/* compute clusters locally in each processor */
	error = Clustering(mg,minlevel,cluster_depth,threshold);
	
stage1: /* compute total number of clusters or error */
	if (error>0) nc=-1; else nc=total_cnt;
	for (l=degree-1; l>=0; l--)
	{
		GetConcentrate(l,&i,sizeof(INT));
		if ((i<0)||(nc<0)) 
			nc = -1;
		else
			nc += i;
	}
	Concentrate(&nc,sizeof(INT));
	Broadcast(&nc,sizeof(INT));
	if (nc<0)
	{
		Release(MGHEAP(mg),FROM_TOP, MarkKeyTop);
		Release(MGHEAP(mg),FROM_BOTTOM, MarkKeyBottom);
		UserWrite("error in stage 1\n");
		return(1);
	}
	if (nc>=MAXCLUSTERS)
	{
		Release(MGHEAP(mg),FROM_TOP, MarkKeyTop);
		Release(MGHEAP(mg),FROM_BOTTOM, MarkKeyBottom);
		sprintf(buf,"Not enough cluster memory: MAXCLUSTERS=%d\n",
                        MAXCLUSTERS);
		UserWrite(buf);
		return(2);
	}

	/* no error occured so far, memory is enough */
	
	/* prepare cluster transfer to master processor */ 
	error = InitConcentrateClusters();

	/* number cluster locally, but uniquely and compute neighbor ids */
	error = ComputeGraphInfo(mg); 

	/* transfer all clusters to master processor */
	error = ConcentrateClusters();

	/* now increasingly numbered clusters should be in master processor */

	/* balance load on master processor */
	error = balance_ccptm(mg,Const,strategy,eigen,loc,dims,weights,coarse,mode);

	/* optimize cluster to processor mapping */
/*	error = OptimizeMapping(); */

	/* check for memory error */
	Broadcast(&MEM_OK,sizeof(double *));
	if (!MEM_OK) 
	{
		UserWrite("Sorry, not enough memory for partitioning\n");
		goto mem_err; 
	}

	/* write back results to the corresponding sources */
	error = BroadcastDestinations();

	/* now assign each element to its destination */
	error = SetDestinations(mg,minlevel,elements,&ne);

	/* release memory for clusters */
	Release(MGHEAP(mg),FROM_BOTTOM, MarkKeyBottom);

	End = CURRENT_TIME;
	sprintf(buf,"BALAN: T=%.2f\n",End-Begin);
	UserWrite(buf);

	/* transfer load */
	/* error = LoadTransfer(mg,ne,elements,element_limit,channel_limit,iter);*/

	TransferGridFromLevel(mg,minlevel);
/*
	TransferGridFromCoarse(mg);
*/

mem_err: /* if memory error has occured leave everything unchanged */
	if (!MEM_OK) Release(MGHEAP(mg),FROM_BOTTOM, MarkKeyBottom);
	
	/* release memory for elements array */
	Release(MGHEAP(mg),FROM_TOP, MarkKeyTop);
	
	return(error);
}



