// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  refine.c														*/
/*																			*/
/* Purpose:   unstructured grid refinement using a general element concept	*/
/*			  (dimension independent for 2/3D)								*/
/*																			*/
/* Author:	  Stefan Lang                         							*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   08.08.95 begin serial algorithm, ug version 3.0				*/
/*																			*/
/* Remarks:   - level 0 grid consists of red elements only					*/
/*			  - the only restriction in the element hierarchie is that 		*/
/*				green or yellow elements might not have sons of class 		*/
/*				green or red												*/
/*			  - the rule set for refinement consists of regular (red)		*/
/*				and irregular rules; regular rules create red elements		*/
/*				while irregular rules result in green elements				*/
/*				( green elements are needed for the closure of the grid,	*/
/*				yellow elements, which are from copy rules, save			*/
/* 				the numerical properties of the solver and are handsome for */
/*				the discretisation											*/
/*			  - if the rule set for the red rules is not complete for build-*/
/*				up a consistent red refined region the FIFO might be used 	*/
/*				for some (hopefully not too much) iterations to find a 		*/
/*				consistent one												*/
/*			  - in 2D: exists a complete rule set for grids of triangles 	*/
/*					   and quadrilaterals exclusivly						*/
/*			  - in 3D: exists a complete rule set for tetrahedrons			*/
/*					   and we assume after some analysation a complete set  */
/*					   of rules described by an algorithm for hexahedrons	*/
/*			  - for mixed element types in arbitrary dimension no 		    */
/*				rule set for the closure exists								*/
/*			  - BEFORE refinement we assume a situation where the error		*/
/*				estimator has detected and marked the leaf elements for 	*/
/*				further refinement											*/
/*			  - AFTER refinement all elements are refined by a rule in way	*/
/*				that no hanging nodes remain (this is the default mode)		*/
/*				or with hanging nodes (in the hanging node mode)			*/
/*				if you use inconsistent red refinement, you need to tell 	*/
/* 				the algorithm explicitly to use the FIFO					*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files 									*/
/*																			*/
/****************************************************************************/

/* standard C library */
#include <assert.h>
#include <math.h>
#include <stdio.h>

/* low module */
#include "compiler.h"
#include "debug.h"
#include "heaps.h"
#include "misc.h"

/* dev module */
#include "devices.h"

/* gm module */
#include "algebra.h"
#include "evm.h"
#include "gm.h"
#include "GenerateRules.h"
#include "refine.h"
#include "rm.h"
#include "switch.h"
#include "ugm.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MINVNCLASS      2           /* determines copies, dep. on discr. !  */

/* defines for side matching of elements 8 bits:                            */
/* _ _ _ _ (4 bits for corner of one element) _ _ _ _ (4 bits for the other)*/ 

#define LINEPOINTS		51			/* 0011 0011 */
#define TRIPOINTS		119			/* 0111 0111 */
#define QUADPOINTS		255			/* 1111 1111 */

#define MAX_GREEN_SONS		32 		/* maximal number of sons for green refinement */

/* TODO: delete next line */
/* #define GET_PATTERN(e)				((SIDEPATTERN(e)<<EDGEPATTERN_LEN) | EDGEPATTERN(e)) */
#define EDGE_IN_PATTERN(p,i)		(((p)[(i)]) & 0x1)
#define SIDE_IN_PATTERN(p,i)		(((p)[EDGEPATTERN_LEN+(i)]) & 0x1)
#define EDGE_IN_PAT(p,i)			(((p)>>(i)) & 0x1)
#define SIDE_IN_PAT(p,i)			(((p)>>(i)) & 0x1)

#define MARK_BISECT_EDGE(r,i)		(r->pattern[i]==1)

#define SIDES_OF_ELEMDESC(t)		(element_descriptors[t]->sides_of_elem)
#define EDGES_OF_ELEMDESC(t)		(element_descriptors[t]->edges_of_elem)
#define CORNERS_OF_ELEMDESC(t)		(element_descriptors[t]->corners_of_elem)

#define REF_TYPE_CHANGES(e)			((REFINE(e)!=MARK(e)) || (REFINECLASS(e)!=MARKCLASS(e)))

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef NODE *ELEMENTCONTEXT[MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM];


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/



/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

static int rFlag=GM_REFINE_TRULY_LOCAL; /* type of refine                   */
static int hFlag=0;						/* refine with hanging nodes?		*/
static int fifoFlag=0;					/* use fifo? 0=no 1=yes				*/
static int first;						/* fifo loop counter 				*/
static ELEMENT *fifo_first=NULL;		/* first element in fifo work list	*/
static ELEMENT *fifo_last=NULL;			/* last element in fifo	work list	*/
static ELEMENT *fifo_insertfirst=NULL;	/* first element in fifo insertlist */
static ELEMENT *fifo_insertlast=NULL;	/* last element in fifo insertlist	*/

/* determine number of edge from reduced (i.e. restricted to one side) edgepattern */
/* if there are two edges marked for bisection, if not deliver -1. If the edge-    */
/* is not reduced (i.e. marked edges lying on more than one side) deliver -2       */
static INT TriSectionEdge[64][2] = {  {-1,-1},{-1,-1},{-1,-1},{ 1, 0},{-1,-1},{ 0, 2},{ 2, 1},{-1,-1},
                               {-1,-1},{ 3, 0},{-2,-2},{-2,-2},{ 2, 3},{-2,-2},{-2,-2},{-2,-2},
                               {-1,-1},{ 0, 4},{ 4, 1},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},
                               { 4, 3},{-1,-1},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},
                               {-1,-1},{-2,-2},{ 1, 5},{-2,-2},{ 5, 2},{-2,-2},{-2,-2},{-2,-2},
                               { 3, 5},{-2,-2},{-2,-2},{-2,-2},{-1,-1},{-2,-2},{-2,-2},{-2,-2},
                               { 5, 4},{-2,-2},{-1,-1},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},
                               {-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2},{-2,-2}  };

/* the indices of the edges of each side */
static INT  CondensedEdgeOfSide[4] = {0x07,0x32,0x2C,0x19};

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

#ifdef ModelPTest
/****************************************************************************/
/*																			*/
/* Function:  MakeRefMarkandMarkClassConsistent								*/
/*																			*/
/* Purpose:	  exchange the MARK and MARKCLASS flags between elements on 	*/
/*			  the vertical boundary of one level.							*/
/*																			*/
/* Input:	  int level: level for which to make flags consistent 	 		*/
/*																			*/
/* Output:	  void 															*/
/*																			*/
/****************************************************************************/

int GetMarkandMarkClass (OBJECT obj, void *data)
{
}

int PutMarkandMarkClass (OBJECT obj, void *data)
{
}

void MakeRefMarkandMarkClassConsistent (int level)
{
	INTERFACE id; 

	/* get id for vertical interface downwards */
	id = 1;
	/* exchange marks */
	DDD_IFExchange(id,INT,GetMarkandMarkClass, PutMarkandMarkClass);
}
#endif


/****************************************************************************/
/*																			*/
/* Function:  DropMarks 													*/
/*																			*/
/* Purpose:   drop marks from leafelements to first regular, and reset		*/
/*			  marks on all elements above (important for restrict marks) 	*/
/*																			*/
/* Param:	  MULTIGRID *theMG												*/
/*																			*/
/* return:	  INT 0: ok 													*/
/*			  INT 1: error													*/
/*																			*/
/****************************************************************************/

static INT DropMarks (MULTIGRID *theMG)
{
	INT k, Mark;
	GRID *theGrid;
	ELEMENT *theElement, *FatherElement;

	for (k=theMG->topLevel; k>0; k--)
	{
		theGrid = GRID_ON_LEVEL(theMG,k);
		for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
			#ifdef ModelP
			if ((MARKCLASS(theElement) == RED) && (ECLASS(theElement) != RED))
			#else
			if (LEAFELEM(theElement) && (MARKCLASS(theElement) == RED) && (ECLASS(theElement) != RED))
			#endif
			{
				Mark = MARK(theElement);
				/* TODO: marks must be changed if element type changes */
				if (TAG(theElement)!=HEXAHEDRON && TAG(EFATHER(theElement))==HEXAHEDRON)  Mark = HEXA_RED;
				FatherElement = theElement;

				#ifdef ModelP
				SETMARK(FatherElement,NO_REFINEMENT);
				SETMARKCLASS(FatherElement,0);
				FatherElement = EFATHER(FatherElement);
				#else
				while(ECLASS(FatherElement) != RED)
				{	
					SETMARK(FatherElement,NO_REFINEMENT);
					SETMARKCLASS(FatherElement,0);
					FatherElement = EFATHER(FatherElement);
				}
				#endif

				SETMARK(FatherElement,Mark);
				SETMARKCLASS(FatherElement,RED);

				#ifdef ModelPTest
				MakeRefMarkandMarkClassConsistent(k);
				#endif
			}
	}
	return (0);
}

#ifdef ModelPTest
/****************************************************************************/
/*																			*/
/* Function:  ExchangePatternOfMasterAndSlaves 								*/
/*																			*/
/* Purpose:	  exchange the PATTERN between elements on horizontal       	*/
/*			  boundary of one level.										*/
/*																			*/
/* Input:	  int level: level for which to make flags consistent 	 		*/
/*																			*/
/* Output:	  void 															*/
/*																			*/
/****************************************************************************/

int GetEdgePatternOfElement (OBJECT obj, void *data)
{
}

int PutEdgePatternOfElement (OBJECT obj, void *data)
{
}

/* TODO: perhaps it is better to exchange the PATTERN to check that they are */
/*       consistent then use the name ExchangePatternOfMasterToSlaves        */
void SendPatternFromMasterToSlaves(int level)
{
	int id;

	/* get interface id for horizontal interface */
	id = 1;
	DDD_IFExchange(id,INT,GetEdgePatternOfElement, PutEdgePatternOfElement);
}
#endif


/****************************************************************************/
/*																			*/
/* Function:  CloseGrid 													*/
/*																			*/
/* Purpose:   compute closure for next level. A closure can only be         */
/*			  determined if the rule set for the used elements is complete. */
/* 			  This means that for all side and edge patterns possible for   */
/*			  an element type exists a rule which closes the element.       */
/*	          In this case a FIFO for computing the closure is not needed   */
/*	          any more and the closure can be computed in one step.			*/
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  INT >0: elements will be refined								*/
/*			  INT 0: no elements will be refined							*/
/*			  INT -1: an error occured           		 					*/
/*																			*/
/****************************************************************************/

static int CloseGrid (GRID *theGrid)
{
	ELEMENT *theElement, *NbElement, *firstElement;
	EDGE *MyEdge, *NbEdge;
	INT i,j,k,cnt,n;
	INT Mark, MyEdgePattern, MySidePattern, MyEdgeNum, MyRule, MyPattern, NewPattern;
	INT NbEdgePattern, NbSidePattern, NbEdgeNum, NbSideMask;
	SHORT *myPattern;

	/* reset USED flag of elements and PATTERN and ADDPATTERN flag on the edges */
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
	{
		SETUSED(theElement,0);
		for (j=0; j<EDGES_OF_ELEM(theElement); j++)
		{
			MyEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
			assert (MyEdge != NULL);
			SETPATTERN(MyEdge,0);
			/* This is needed in RestrictMarks() */
			SETADDPATTERN(MyEdge,1);
		}
	}
	
	/* reset EDGE/SIDEPATTERN in elements, set SIDEPATTERN in elements for quadrilateral sides and set PATTERN on the edges */
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
	{
		if (MARKCLASS(theElement)==RED)
		{
			Mark = MARK(theElement);
			myPattern = MARK2PATTERN(theElement,Mark);

			for (i=0; i<EDGES_OF_ELEM(theElement); i++)
				/* TODO: delete this (old code) */
				/* if (myPattern[i]) */
				if (EDGE_IN_PATTERN(myPattern,i))
				{
					MyEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
					if (MyEdge != NULL)
						SETPATTERN(MyEdge,1);
					else 
						UserWriteF("CloseGrid(): ERROR edge i=%d of element e=%x not found!",i,theElement);
				}

			SETSIDEPATTERN(theElement,0);

			/* for 2D we are finished */
			if (DIM == 2) continue;

			for (i=0;i<SIDES_OF_ELEM(theElement); i++)
			{
				if (CORNERS_OF_SIDE(theElement,i)==4)
				{
					/* set SIDEPATTERN if side has node */
					if(SIDE_IN_PATTERN(myPattern,i))
						SETSIDEPATTERN(theElement,SIDEPATTERN(theElement) | 1<<i);
				}
			}					
		}
		else
		{
			SETSIDEPATTERN(theElement,0);
			SETMARKCLASS(theElement,0);
		}
		/* TODO: delete this */
		/* SETEDGEPATTERN(theElement,0);*/
	}
	
    firstElement = theGrid->elements;

	if (fifoFlag) {
		fifo_first=fifo_last=fifo_insertfirst=fifo_insertlast=NULL;
		first = 1;
		n = 0;
		UserWriteF("Using FIFO: loop 0\n");
	}

/* here starts the fifo loop */
FIFOSTART:

	/* set pattern (edge and side) on the elements */
	for (theElement=firstElement; theElement!=NULL; theElement=SUCCE(theElement))
	{
		/* make edgepattern consistent with pattern of edges */
		SETUSED(theElement,1);
		MyEdgePattern = 0;
		for (i=EDGES_OF_ELEM(theElement)-1; i>=0; i--)
		{
			MyEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
			MyEdgePattern = (MyEdgePattern<<1) | PATTERN(MyEdge);
		}
		
		/* for 2D we are finished */
		if (DIM == 2) continue;

		/* TODO: change this for red refinement of pyramids */
		if (DIM==3 && TAG(theElement)==PYRAMID) continue;

		/* make SIDEPATTERN consistent with neighbors	*/
		/* TODO: compute this in a separate function?	*/
		for (i=0; i<SIDES_OF_ELEM(theElement); i++)
		{
			NbElement = NBELEM(theElement,i);
			if (NbElement == NULL) continue;
			/* TODO: is this or only USED ok? */
			if (!USED(NbElement)) continue;
			/* now edgepattern from theelement and NbElement are in final state */

			/* search neighbors side */
			for (j=0; j<SIDES_OF_ELEM(NbElement); j++)
				if (NBELEM(NbElement,j) == theElement)
					break;
			assert(j<SIDES_OF_ELEM(NbElement));

			/* side is triangle or quadrilateral */
			switch (CORNERS_OF_SIDE(theElement,i))
			{
				case 3:
					/* because SIDEPATTERN is set to zero, I choose TriSectionEdgeÉÉ[0] */
					MyEdgeNum = TriSectionEdge[MyEdgePattern&CondensedEdgeOfSide[i]][0];
					if (MyEdgeNum == -2) {
						assert(0);
						return (-1);
					}
					if (MyEdgeNum == -1) continue;
					
					NbEdgePattern = 0;
					for (k=0; k<EDGES_OF_ELEM(NbElement); k++) {
						NbEdge=GetEdge(CORNER(NbElement,CORNER_OF_EDGE(NbElement,k,0)),CORNER(NbElement,CORNER_OF_EDGE(NbElement,k,1)));
						assert(NbEdge!=NULL);
						NbEdgePattern = NbEdgePattern | (PATTERN(NbEdge)<<k);
					}
						
					NbEdgeNum = TriSectionEdge[NbEdgePattern&CondensedEdgeOfSide[j]][0];
					if (NbEdgeNum == -2 || NbEdgeNum == -1) {
						assert(0);
						return (-1);
					}
					
					if (!( CORNER(theElement,CORNER_OF_EDGE(theElement,MyEdgeNum,0)) == CORNER(NbElement,CORNER_OF_EDGE(NbElement,NbEdgeNum,0)) &&
						   CORNER(theElement,CORNER_OF_EDGE(theElement,MyEdgeNum,1)) == CORNER(NbElement,CORNER_OF_EDGE(NbElement,NbEdgeNum,1))	) &&
						!( CORNER(theElement,CORNER_OF_EDGE(theElement,MyEdgeNum,0)) == CORNER(NbElement,CORNER_OF_EDGE(NbElement,NbEdgeNum,1)) &&
						   CORNER(theElement,CORNER_OF_EDGE(theElement,MyEdgeNum,1)) == CORNER(NbElement,CORNER_OF_EDGE(NbElement,NbEdgeNum,0))	)	 )
					{
						NbSidePattern = SIDEPATTERN(NbElement);
						NbSideMask = (1<<j);
						if ( NbSidePattern & NbSideMask )
							NbSidePattern &= ~NbSideMask;
						else
							NbSidePattern |= NbSideMask;
						SETSIDEPATTERN(NbElement,NbSidePattern);
					}
					break;
				case 4:
					/* if side of one of the neighboring elements has a sidenode, then both need a sidenode */
					NbSidePattern = SIDEPATTERN(NbElement);
					if (SIDE_IN_PAT(SIDEPATTERN(theElement),i))
					{
						SETSIDEPATTERN(NbElement,SIDEPATTERN(NbElement) | (1<<j));
					}
					else if (SIDE_IN_PAT(SIDEPATTERN(NbElement),j))
					{
						SETSIDEPATTERN(theElement,SIDEPATTERN(theElement) | (1<<i));
					}
					break;
				default:
					UserWriteF("CloseGrid(): ERROR: CORNER_OF_SIDE(e=%x,s=%d)=%d !\n",theElement,i,CORNERS_OF_SIDE(theElement,i));
					assert(0);
					return(-1);
			}
		}
	}

	#ifdef ModelPTest
	/* send the PATTERN flag from master to slave elements */
	SendPatternFromMasterToSlaves(GLEVEL(theGrid));
	#endif

	/* set refinement rules from edge- and sidepattern */
	cnt = 0;
	for (theElement=firstElement; theElement!=NULL; theElement=SUCCE(theElement))
	{
		MyEdgePattern = 0;
		for (i=EDGES_OF_ELEM(theElement)-1; i>=0; i--)
		{
			MyEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
			MyEdgePattern = (MyEdgePattern<<1) | PATTERN(MyEdge);
		}
		MySidePattern = SIDEPATTERN(theElement);
		MyPattern = MySidePattern<<EDGES_OF_ELEM(theElement) | MyEdgePattern;

		Mark = PATTERN2MARK(theElement,MyPattern);

		if (fifoFlag) {
			if (Mark == -1 && MARKCLASS(theElement)==RED) {
				/* there is no rule for this pattern, switch to red */
				Mark = RED;
			}
			else {
				assert(Mark != -1);
			}
		}			
		else if (hFlag==0 && MARKCLASS(theElement)!=RED) {
			Mark = 0;
		}
		else {
			assert(Mark != -1);
		}

		IFDEBUG(gm,1)
		UserWriteF("    ID=%d TAG=%d MyPattern=%d Mark=%d MyEdgePattern=%d MySidePattern=%d\n",ID(theElement),TAG(theElement),MyPattern,Mark,MyEdgePattern,MySidePattern);
		ENDDEBUG

		#ifdef __THREEDIM__
		if (TAG(theElement)==TETRAHEDRON && MARKCLASS(theElement) == RED)
			Mark = (*theFullRefRule)(theElement);
		#endif

		NewPattern = MARK2PAT(theElement,Mark);
		PRINTDEBUG(gm,1,("   MyPattern=%d NewPattern=%d Mark=%d\n",MyPattern,NewPattern,Mark));

		if (fifoFlag) {
          if (MARKCLASS(theElement)==RED && MyPattern != NewPattern) {
				#ifdef __TWODIM__
				for (j=0; j<EDGES_OF_ELEM(theElement); j++) {
                     if (EDGE_IN_PAT(MyPattern,j)==0 && EDGE_IN_PAT(NewPattern,j)) {
					     MyEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
					     if (MyEdge != NULL)
						    SETPATTERN(MyEdge,1);
					     else 
						    UserWriteF("CloseGrid(): ERROR edge i=%d of element e=%x not found!",i,theElement);
                      
						/* boundary case */
						if (SIDE(theElement,j) != NULL) continue;

						/* add the element sharing this edge to fifo_queue */
						NbElement = NBELEM(theElement,j);

                        if (NbElement==NULL) continue; 

						PRINTDEBUG(gm,1,("   ADDING to FIFO: NBID=%d\n",ID(NbElement)))

						/* unlink element from element list */
						if (PREDE(NbElement) != NULL)
							SUCCE(PREDE(NbElement)) = SUCCE(NbElement);
						if (SUCCE(NbElement) != NULL)
							PREDE(SUCCE(NbElement)) = PREDE(NbElement);
                        if (theGrid->elements == NbElement) theGrid->elements = SUCCE(NbElement);

						SUCCE(NbElement) = PREDE(NbElement) = NULL;
						/* insert into fifo */
						if (fifo_insertfirst == NULL) {
							fifo_insertfirst = fifo_insertlast = NbElement;
						}
						else {
							SUCCE(fifo_insertlast) = NbElement;
							PREDE(NbElement) = fifo_insertlast;
							fifo_insertlast = NbElement;
						}
					}
					if (EDGE_IN_PAT(MyPattern,j) && EDGE_IN_PAT(NewPattern,j)==0) {
						UserWriteF("CloseGrid(): ERROR EID=%d in fifo MyPattern=%d has edge=%d refined but NewPattern=%d NOT!\n",ID(theElement),MyPattern,j,NewPattern);
						return(-1);
					}
				}
				#endif
				#ifdef __THREEDIM__
				UserWriteF("CloseGrid(): ERROR fifo for 3D NOT implemented!\n");
				assert(0);
				return(-1);
				#endif
			}
		}

		/* TODO: why here swap from NOREFRULE to COPY_REFRULE */
		if (MARKCLASS(theElement)==RED && Mark==NO_REFINEMENT)
		{
			UserWriteF("CloseGrid(): MARKCLASS=RED && Mark==NO_REFINEMENT -> set Mark=COPY!\n");
			Mark = COPY;
		}
		/* TODO: delete or better  ... && MyPattern != 0 ?? */
		if (MARKCLASS(theElement)!=RED && Mark!=NO_REFINEMENT) {
			SETMARKCLASS(theElement,GREEN);
			IFDEBUG(gm,1)
			UserWriteF("   Switching MARKCLASS=%d for MARK=%d of EID=%d to GREEN\n",MARKCLASS(theElement),MARK(theElement),ID(theElement));
			ENDDEBUG
		}
		if (Mark)
			cnt++;
		SETMARK(theElement,Mark);
	}	

	if (fifoFlag) {
		/* insert fifo work list into elementlist */
		for (theElement=fifo_last; theElement!=NULL; theElement=PREDE(theElement)) {
			SUCCE(theElement) = theGrid->elements;
			PREDE(theGrid->elements) = theElement;
			theGrid->elements = theElement;
		}
		PREDE(theGrid->elements) = NULL;

		if (fifo_insertfirst!=NULL) {
			/* append fifo insert list to fifo work list */
			firstElement = fifo_first = fifo_insertfirst;
			fifo_last = fifo_insertlast;

			IFDEBUG(gm,2)
			UserWriteF(" FIFO Queue:");
			for (theElement=fifo_first; theElement!=NULL; theElement=SUCCE(theElement))
				UserWriteF(" %d\n", ID(theElement));
			ENDDEBUG

			fifo_insertfirst = fifo_insertlast = NULL;
			first = 0;
			n++;
			UserWriteF(" loop %d",n);
			goto FIFOSTART;
		}
	}


	#ifdef FIFO
	#ifdef ModelP
	do {
		/* exchange FIFO flag and PATTERN from slaves to master */
		IF_FIF0AndPat_S2M(GLEVEL(grid));

		/* add all master elements of horizontal interface to FIFO */
		for (theElement=firstElement; theElement!=NULL; theElement=SUCCE(theElement)) {
			if (IS_HOR_MASTER(theElement) && FIFO(theElement)) {
				AddToFIFIO(theElement);
				SETFIFO(theElement,0);
			}
		}

		/* check condition for termination of pattern adaptation */
		AllFIFOsEmpty = CheckGlobalFIFOStatus(fifo);
	}
	while (fifo==NULL && AllFIFOsEmpty==1)
	#endif /* ModelP */
	} /* end if (fifoFlag) */
	#endif

	
	/* set additional pattern on the edges */
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
	{
		if (MARKCLASS(theElement)!=RED) continue;
		for (j=0; j<EDGES_OF_ELEM(theElement); j++)
		{
			MyEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
			assert (MyEdge != NULL);
			/* ADDPATTERN is now set to 0 for all edges of red elements */
			SETADDPATTERN(MyEdge,0);
		}
	}
	
	/* build a green covering around the red elements */
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
	{
		if (MARKCLASS(theElement)==RED) continue;

		for (i=0; i<EDGES_OF_ELEM(theElement); i++) {
			MyEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
			assert (MyEdge != NULL);
			if (ADDPATTERN(MyEdge) == 0) {
				/* TODO: when does this case occur? */
				/* for HEXAHEDRA Patterns2Rules returns 0 for not red elements */
				/*       perhaps if the side of the red element has no edges markes and no side node */
				if (MARK(theElement) == NO_REFINEMENT) {
					SETMARK(theElement,COPY);
					IFDEBUG(gm,0)
					UserWriteF("   WARNING: EID=%d switching MARK from NO_REFINEMENT to COPY\n",ID(theElement));
					ENDDEBUG
				}
				SETMARKCLASS(theElement,GREEN);
				break;
			}
		}
	}	

	#ifdef ModelPTest
	/* send MARKCLASS from slave to masters */
	IF_S2M_MARKCLASS(GLEVEL(grid));
	#endif
	
	return(cnt);
}


/****************************************************************************/
/*																			*/
/* Function:  GetNeighborSons												*/
/*																			*/
/* Purpose:   fill SonList for theElement with a breadth first search		*/
/*																			*/
/* Param:	  ELEMENT *theElement:	father element							*/
/*			  ELEMENT *theSon: currently visited son						*/ 
/*			  ELEMENT *SonList[MAX_SONS]: the list of sons					*/
/*			  int 	count: son count										*/
/*			  int   nsons: number of sons									*/
/*																			*/
/* return:	  int: new son count											*/
/*																			*/
/****************************************************************************/

INT GetNeighborSons (ELEMENT *theElement, ELEMENT *theSon, ELEMENT *SonList[MAX_SONS], int count, int nsons)
{
	ELEMENT *NbElement;
	int i,j, startson, stopson;

	startson = count;

	for (i=0; i<SIDES_OF_ELEM(theSon); i++)
	{
		NbElement = NBELEM(theSon,i);
		if (NbElement == NULL) continue;
		if (EFATHER(NbElement) == theElement) 
		{
			/* is NbElement already in list */
			for (j=0; j<count; j++)
				if (SonList[j] == NbElement)
					break;
			if (j==count && count<nsons)
				SonList[count++] = NbElement;
		}
	}
	if (count == nsons) return(count);

	stopson = count;
	for (i=startson; i<stopson; i++)
	{
		if (count<nsons) count = GetNeighborSons(theElement,SonList[i],SonList,count,nsons);
		else return(count);
	}
	return(count);
}


/****************************************************************************/
/*																			*/
/* Function:  GetSons														*/
/*																			*/
/* Purpose:   fill SonList for theElement									*/
/*																			*/
/* Param:	  ELEMENT *theElement, ELEMENT *SonList[MAX_SONS]				*/
/*																			*/
/* return:	  0: ok 														*/
/*			  1: error														*/
/*																			*/
/****************************************************************************/

INT GetSons (ELEMENT *theElement, ELEMENT *SonList[MAX_SONS])
{
	REFRULE *theRule;
	ELEMENT *theSon;
	int SonID,PathPos,nsons;
	
	if (theElement==NULL) return (GM_ERROR);
	
	for (SonID=0; SonID<MAX_SONS; SonID++)
		SonList[SonID] = NULL;
	
	switch (TAG(theElement))
	{
		#ifdef __TWODIM__
		case (TRIANGLE):
			for (SonID=0;SonID<NSONS(theElement);SonID++)	
				SonList[SonID] = SON(theElement,SonID);
			break;
		case (QUADRILATERAL):
			for (SonID=0;SonID<NSONS(theElement);SonID++)	
				SonList[SonID] = SON(theElement,SonID);
			break;
		#endif
		#ifdef __THREEDIM__
		case (TETRAHEDRON):
			SonList[0] = SON(theElement,0);
			
			/* get other sons from path info in rules */
			theRule = MARK2RULEADR(theElement,REFINE(theElement));
			
			for (SonID=1; SonID<theRule->nsons; SonID++)
			{
				theSon = SonList[0];
				for (PathPos=0; PathPos<PATHDEPTH(theRule->sons[SonID].path); PathPos++)
					theSon = NBELEM(theSon,NEXTSIDE(theRule->sons[SonID].path,PathPos));
				
				if (theSon==NULL)
					return (GM_ERROR);
				
				SonList[SonID] = theSon;
			}
			break;
		case (PYRAMID):
			IFDEBUG(gm,0)
			if (REFINECLASS(theElement) != YELLOW)
				UserWriteF("GetSons(): ERROR PYRAMID has REFINECLASS=%d and MARK=%d\n",REFINECLASS(theElement),MARK(theElement));
			ENDDEBUG
			SonList[0] = SON(theElement,0);
			break;
		case (HEXAHEDRON):
			SonList[0] = SON(theElement,0);

			if (REFINECLASS(theElement) == GREEN)
			{
				if (NSONS(theElement)==0 || SonList[0]==NULL) return(GM_ERROR);
				nsons = 1;
				if (NSONS(theElement)>1)
					nsons = GetNeighborSons(theElement,SonList[0],SonList,1,NSONS(theElement));

				if (nsons != NSONS(theElement))
				{
					PRINTDEBUG(gm,2,("GetSons(): ERROR! Element ID=%d, NSONS=%d but nsons=%d\n",ID(theElement),NSONS(theElement),nsons))
					return(GM_ERROR);
				}
			}
			else
			{
				/* get other sons from path info in rules */
				theRule = MARK2RULEADR(theElement,REFINE(theElement));
				
				for (SonID=1; SonID<theRule->nsons; SonID++)
				{
					theSon = SonList[0];
					for (PathPos=0; PathPos<PATHDEPTH(theRule->sons[SonID].path); PathPos++)
						theSon = NBELEM(theSon,NEXTSIDEHEX(theRule->sons[SonID].path,PathPos));
					
					if (theSon==NULL)
						return (GM_ERROR);
					
					SonList[SonID] = theSon;
				}
			}
			break;
		#endif
		default:
			UserWriteF("GetSons(): ERROR TAG(e=%x)=%d !\n",theElement,TAG(theElement));
			return (GM_ERROR);
	}
	
	return (GM_OK);
}


/****************************************************************************/
/*																			*/
/* Function:  RestrictMarks 												*/
/*																			*/
/* Purpose:   restrict refinement marks when going down 					*/
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static void RestrictMarks (GRID *theGrid)
{
	ELEMENT *theElement,*SonList[MAX_SONS];
	EDGE *theEdge;
	int i,j,flag,CondensedPattern,Pattern;
	
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
	{
		if (GetSons(theElement,SonList)!=GM_OK) return;

		if (hFlag)
		{
			if (REFINE(theElement) == NO_REFINEMENT ||  /* if element is not refined anyway, then there are no restrictions to apply */
				ECLASS(theElement) == YELLOW ||   
				ECLASS(theElement) == GREEN ||      /* irregular elements are marked by estimator, because they are leaf elements */
				REFINECLASS(theElement) == YELLOW ) /* regular elements with YELLOW copies are marked by estimator, because the marks are dropped */
				continue;

			/* regular elements with GREEN refinement go to no refinement or red refinement */
			if (REFINECLASS(theElement)==GREEN)
			{
				for (i=0; i<NSONS(theElement); i++)
					/* Is the son marked for further refinement */ 
					if (MARK(SonList[i])>NO_REFINEMENT)
					{
						if (MARKCLASS(theElement)==RED)
						{
							/* TODO: this mark is from DropMarks()! */
							/* theElement is marked from outside */
							/* TODO: edit this for new element type or for different restrictions */
							switch (TAG(theElement)) {
								#ifdef __TWODIM__
									case TRIANGLE:
										SETMARK(theElement,T_RED);
										break;
									case QUADRILATERAL:
										SETMARK(theElement,Q_RED);
										break;
								#endif
								#ifdef __THREEDIM__
								case TETRAHEDRON:
									if (MARK(theElement)!=RED) 
										/* TODO: Is REFINE always as red rule available? */
										SETMARK(theElement,REFINE(theElement));
									break;
								case HEXAHEDRON:
									SETMARK(theElement,HEXA_RED);
									break;
								#endif
								default:
									UserWriteF("RestrictMarks(): for elementtype=%d mark restriction not implemented!\n",TAG(theElement));
									SETMARK(theElement,PATTERN2MARK(theElement,0));
									break;
							}
									
						}
						else
						{
							/* TODO: edit this for new element type or for different restrictions */
							switch (TAG(theElement)) {
								#ifdef __TWODIM__
									case TRIANGLE:
										SETMARK(theElement,T_RED);
										break;
									case QUADRILATERAL:
										SETMARK(theElement,Q_RED);
										break;
								#endif
								#ifdef __THREEDIM__
								case TETRAHEDRON:
									/* theElement is not marked from outside, so find a regular rule being consistent */
									/* with those neighbors of all sons of theElement which are marked for refine.	  */
									/* this choice will make sure these marks will not be distroyed.				  */
									Pattern = RULE2PAT(theElement,REFINE(theElement));
									for (j=0; j<EDGES_OF_ELEM(theElement); j++)
									{
										theEdge=GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,j,1)));
										assert (theEdge != NULL);
										/* TODO: What's on when MIDNODE exists?? */
										if (MIDNODE(theEdge)==NULL)
										{
											theEdge=GetEdge(SONNODE(CORNER(theElement,CORNER_OF_EDGE(theElement,j,0))),SONNODE(CORNER(theElement,CORNER_OF_EDGE(theElement,j,1))));
											assert(theEdge != NULL);
											/* TODO: Is ADDPATTERN needed for fitting with other green elements?? */ 
											if (ADDPATTERN(theEdge))
												Pattern |= (1<<j);
										}
									}
									SETMARK(theElement,PATTERN2MARK(theElement,Pattern));
									break;
								case HEXAHEDRON:
									SETMARK(theElement,HEXA_RED);
									break;
								#endif
								default:
									UserWriteF("RestrictMarks(): for elementtype=%d mark restriction not implemented!\n",TAG(theElement));
									SETMARK(theElement,PATTERN2MARK(theElement,0));
									break;
							}

							SETMARKCLASS(theElement,RED);
						}
						/* this must be done only once for each element "theElement" */
						break;
					}
				continue;
			}
					
			/* regular elements with regular refinement, are the only ones to coarsen */
			if (REFINECLASS(theElement) == RED)
			{
				SETMARK(theElement,REFINE(theElement));
				SETMARKCLASS(theElement,REFINECLASS(theElement));
			}
		}

		flag = 0;
		for (i=0; i<NSONS(theElement); i++)
			/* if not all sons are marked no unrefinement is possible */
			if (!COARSEN(SonList[i]))
			{
				flag = 1;
				break;
			}
			
		if (flag) continue;
	
		/* remove refinement */
		SETMARK(theElement,NO_REFINEMENT);
		SETMARKCLASS(theElement,0);
		SETCOARSEN(theElement,1);
	}
}


/****************************************************************************/
/*																			*/
/* Function:  ComputeCopies 												*/
/*																			*/
/* Purpose:   determine copy elements from node classes 					*/
/*																			*/
/* Param:	  GRID *theGrid: pointer to grid structure						*/
/*																			*/
/* return:	  0: ok 														*/
/*																			*/
/****************************************************************************/

static int ComputeCopies (GRID *theGrid)
{
	ELEMENT *theElement;
	int i,flag;
	
	/* set class of all dofs on next level to 0 */
	ClearNextVectorClasses(theGrid);
	
	/* seed dofs of regularly and irregularly refined elements to 3 */
	flag = 0;
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
		if (MARK(theElement)!=NO_REFINEMENT && (MARKCLASS(theElement)==RED || MARKCLASS(theElement)==GREEN))
		{
			SeedNextVectorClasses(theGrid,theElement);
			flag=1; /* there is at least one element to be refined */
		}

	/* copy all option or neighborhood */	
	if (rFlag==GM_COPY_ALL) 
	{
		if (flag)
			for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
				SeedNextVectorClasses(theGrid,theElement);
	}
	else
	{
		PropagateNextVectorClasses(theGrid);
	}
	
	/* an element is copied if it has a dof of class 2 and higher */
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement)) {
		if ((MARK(theElement)==NO_REFINEMENT)&&(MaxNextVectorClass(theGrid,theElement)>=MINVNCLASS))
		{
			SETMARK(theElement,COPY);
			SETMARKCLASS(theElement,YELLOW);
		}
	}

	return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  GetCurrentContext 											*/
/*																			*/
/* Purpose:   assemble references to objects which interact with the sons	*/
/*			  of the given element, as indicated by REFINE. 				*/
/*			  (i)	 corner nodes											*/
/*			  (ii)	 nodes at midpoints of edges							*/
/*			  (iii)	 nodes at center of element 							*/
/*																			*/
/* Param:	  ELEMENT *theElement: element to refine						*/
/*			  ELEMENTCONTEXT *theContext: context structure to fill 		*/
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static void GetCurrentContext (ELEMENT *theElement, NODE **theElementContext)
{
	ELEMENT *SonList[MAX_SONS];
	NODE *theNode,*theNode0,*theNode1;
	EDGE *theEdge;
	LINK *theLink0,*theLink1,*theLink;
	INT i,l,k,nodeindex;						
	REFRULE *rule;
	
	/* init Context to NULL */
	for(i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++)  
		theElementContext[i] = NULL;

	/* get nodes, can be NULL */
	for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
		theElementContext[i] = SONNODE(CORNER(theElement,i));

	if (DIM==3 && TAG(theElement)==HEXAHEDRON && REFINECLASS(theElement)==GREEN) {

		for (i=0; i<EDGES_OF_ELEM(theElement); i++)
		{
			theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
			if ((theNode = MIDNODE(theEdge)) != NULL)
				theElementContext[i+CORNERS_OF_ELEM(theElement)] = theNode;
			else
				theElementContext[i+CORNERS_OF_ELEM(theElement)] = NULL;
		}
		for (i=0; i<SIDES_OF_ELEM(theElement); i++)
		{
			theNode = NULL;
			theNode0 = theElementContext[EDGE_OF_SIDE(theElement,i,0)+CORNERS_OF_ELEM(theElement)];
			theNode1 = theElementContext[EDGE_OF_SIDE(theElement,i,2)+CORNERS_OF_ELEM(theElement)];
			if (theNode0 != NULL && theNode1 != NULL)
			  for (theLink0=START(theNode0); theLink0!=NULL; theLink0=NEXT(theLink0)) {
				for (theLink1=START(theNode1); theLink1!=NULL; theLink1=NEXT(theLink1)) 
					if (NBNODE(theLink0) == NBNODE(theLink1)) {
						IFDEBUG(gm,3)
						UserWriteF("    		possible node ID=%d\n",ID(NBNODE(theLink0)));
						ENDDEBUG

						l = 0;
						for (theLink=START(NBNODE(theLink0)); theLink!=NULL; theLink=NEXT(theLink)) 
							for (k=0; k<CORNERS_OF_SIDE(theElement,i); k++) 
								if (NBNODE(theLink) == SONNODE(CORNER(theElement,CORNER_OF_SIDE(theElement,i,k))))
									l = 1;
						if (l == 0) {
							theNode = NBNODE(theLink0);
							IFDEBUG(gm,3)
							UserWriteF("    		FOUND node ID=%d\n",ID(NBNODE(theLink)));
							continue;
							ENDDEBUG
							break;
						}
					}
				if (theNode != NULL) {
					IFDEBUG(gm,3)
					continue;
					ENDDEBUG
					break;
				}
			}
			theElementContext[i+CORNERS_OF_ELEM(theElement)+EDGES_OF_ELEM(theElement)] = theNode;
		}
		/* get center node */
		theNode = NULL;
		/* TODO: not independent of element corner numbering */
		theNode0 = SONNODE(CORNER(theElement,0));
		theNode1 = SONNODE(CORNER(theElement,6));
		for (theLink0=START(theNode0); theLink0!=NULL; theLink0=NEXT(theLink0)) {
			for (theLink1=START(theNode1); theLink1!=NULL; theLink1=NEXT(theLink1)) 
				if (NBNODE(theLink0) == NBNODE(theLink1)) {
					theNode = NBNODE(theLink0);
					break;
				}
			if (theNode != NULL) break;
		}
		assert(theNode != NULL);
		theElementContext[CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement)] = theNode;
	}
	else {
		if (GetSons(theElement,SonList)!=GM_OK) return;

		/* get midpoints */
		rule = MARK2RULEADR(theElement,REFINE(theElement));

		for (i=0; i<EDGES_OF_ELEM(theElement); i++)
		{								
			if (rule->pattern[i]==1)
				theElementContext[i+CORNERS_OF_ELEM(theElement)] = CORNER(SonList[rule->sonandnode[i][0]],rule->sonandnode[i][1]);
			else
				theElementContext[i+CORNERS_OF_ELEM(theElement)] = NULL;
		}

		/* get center node */
		nodeindex = CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement);
		if (rule->sonandnode[CENTER_NODE_INDEX(theElement)][0] == NO_CENTER_NODE)
			theElementContext[nodeindex] = NULL;
		else
			theElementContext[nodeindex] = CORNER(SonList[rule->sonandnode[CENTER_NODE_INDEX(theElement)][0]],rule->sonandnode[CENTER_NODE_INDEX(theElement)][1]);

		/* for 2D we are finished */
		if (DIM==2) return;

		/* get side nodes */
		for (i=EDGES_OF_ELEM(theElement); i<EDGES_OF_ELEM(theElement)+SIDES_OF_ELEM(theElement); i++)
		{								
			if (rule->pattern[i]==1)
				theElementContext[i+CORNERS_OF_ELEM(theElement)] = CORNER(SonList[rule->sonandnode[i][0]],rule->sonandnode[i][1]);
			else
				theElementContext[i+CORNERS_OF_ELEM(theElement)] = NULL;
		}
	}
}


/****************************************************************************/
/*																			*/
/* Function:  NodeNeededOnlyFor												*/
/*																			*/
/* Purpose:   check wether node is referenced only within theElementContext	*/
/*																			*/
/* Parameters:NODE *theNode													*/
/*            NODE **theElementContext										*/
/*																			*/
/* Return:    INT 0: node is needed also otherwise                          */
/*			  INT 1: node is needed only for theElementContext				*/
/*																			*/
/****************************************************************************/

static INT NodeNeededOnlyFor(NODE *theNode,NODE **theElementContext)
{
	LINK *theLink;
	INT i, found;
	
	for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
	{
		found=0;
		/* TODO: WRONG because NULL pointer might be equal ? */
		for (i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++)
			if (NBNODE(theLink)==theElementContext[i])
				found=1;
		if (!found)
			return (0);
	}
	
	return (1);
}


/****************************************************************************/
/*																			*/
/* Function:  UpdateContext 												*/
/*																			*/
/* Purpose:   assemble references to objects which interact with the sons	*/
/*			  of the given element, i.e.									*/
/*			  objects are allocated, kept or deleted as indicated by MARK	*/
/*			  (i)	 corner nodes											*/
/*			  (ii)	 nodes at midpoints of edges							*/
/*																			*/
/* Param:	  GRID *theGrid: grid level of the sons of theElement			*/
/*			  ELEMENT *theElement: element to refine						*/
/*			  ELEMENTCONTEXT *theContext: context structure to update		*/
/*																			*/
/* return:	  INT 0: ok 													*/
/*			  INT 1: fatal memory error 									*/
/*																			*/
/****************************************************************************/

static int UpdateContext (GRID *theGrid, ELEMENT *theElement, NODE **theElementContext) 
{
	NODE *theNode, **CenterNode; 			
	ELEMENT *theNeighbor,*theSon;			/* neighbor and a son of current elem.	*/
	ELEMENT *NeighborSonList[MAX_SONS];
	EDGE *theEdge,*fatherEdge;				/* temporary storage for an edge		*/
	INT i,j,r,Corner0, Corner1,candelete;	/* some integer variables				*/
	NODE **MidNodes;						/* nodes on refined edges				*/
	NODE **SideNodes;						/* nodes on refined sides				*/
	LINK *theLink;							/* scan through nodes neighbor list 	*/
	LINK *theLink0,*theLink1;
	NODE *Node0, *Node1;
	NODE *theNode0, *theNode1;
	VERTEX *CenterVertex;
	EDGEDATA *edata;
	REFRULE *rule, *nbrule;
	COORD *x, *xi, sx, sy;
	INT Mark,Refine,toBisect,toDelete,toCreate;
	INT k,l;
	char buffer[64];

	Mark = MARK(theElement);

	/* allocate corner nodes if necessary */
	if (Mark>0)
		for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
		{
			theNode = CORNER(theElement,i);
			if (SONNODE(theNode)==NULL)
			{
				SONNODE(theNode) = CreateNode(theGrid,NULL);
				if (SONNODE(theNode)==NULL) return(GM_ERROR);
				MYVERTEX(SONNODE(theNode)) = MYVERTEX(theNode);
				NFATHER(SONNODE(theNode)) = theNode;
				theElementContext[i] = SONNODE(theNode);
			}
		}

	/* delete corner nodes if possible */
	if (Mark==0)
		for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
			if (theElementContext[i]!=NULL)
			{
				if (NodeNeededOnlyFor(theElementContext[i],theElementContext))
				{
					DisposeNode(theGrid,theElementContext[i]);
					theElementContext[i] = NULL; 
					SONNODE(CORNER(theElement,i)) = NULL;
				}
			}
	
	/* allocate,keep, or delete midpoint nodes */
	/* allocate,keep, or delete corner/corner edges */
	MidNodes = theElementContext+CORNERS_OF_ELEM(theElement);
	for (i=0; i<EDGES_OF_ELEM(theElement); i++)
	{
		Corner0 = CORNER_OF_EDGE(theElement,i,0);
		Corner1 = CORNER_OF_EDGE(theElement,i,1);
		
		toBisect = 0;
		if (DIM==3 && TAG(theElement)==HEXAHEDRON && MARKCLASS(theElement)==GREEN) {
			theEdge = GetEdge(CORNER(theElement,Corner0),CORNER(theElement,Corner1));
			assert(theEdge != NULL);
			if (ADDPATTERN(theEdge) == 0)
				toBisect = 1;
		}
		else {
			rule = MARK2RULEADR(theElement,Mark);
			if (MARK_BISECT_EDGE(rule,i))
				toBisect = 1;
		}
		IFDEBUG(gm,2)
		if (MidNodes[i] == NULL)
		  UserWriteF("\n    MidNodes[%d]: toBisect=%d ID(Corner0)=%d ID(Corner1)=%d",i,toBisect,ID(CORNER(theElement,Corner0)),ID(CORNER(theElement,Corner1)));
		else
		  UserWriteF("\n    MidNodes[%d]: toBisect=%d ID(Corner0)=%d ID(Corner1)=%d ID(MidNode)=%d",i,toBisect,ID(CORNER(theElement,Corner0)),ID(CORNER(theElement,Corner1)),ID(MidNodes[i]));
		ENDDEBUG
		if (toBisect)
		{
			/* we need a midpoint node */
			if (MidNodes[i]!=NULL) continue;
			Node0 = CORNER(theElement,Corner0);
			Node1 = CORNER(theElement,Corner1);
			if ((theEdge = GetEdge(Node0,Node1))==NULL)
			 return (GM_ERROR);
			MidNodes[i] = MIDNODE(theEdge);
			if (MidNodes[i] == NULL)
			{
				MidNodes[i] = CreateMidNode(theGrid,theElement,i,theElementContext[Corner0]);
				if (MidNodes[i]==NULL) return(1);
				MIDNODE(theEdge) = MidNodes[i];
				if ((theEdge=CreateEdge(theGrid,theElementContext[Corner0],MidNodes[i]))==NULL) return(1);
				if ((theEdge=CreateEdge(theGrid,theElementContext[Corner1],MidNodes[i]))==NULL) return(1);
			} 
		}
		else
		{
			/* if we need a corner corner edge then allocate it */
			if (Mark>0)
			{
				if ((theEdge=CreateEdge(theGrid,theElementContext[Corner0],theElementContext[Corner1]))==NULL) return(1);
			}
		
			/* we don't need a midpoint node on that edge, lets see if it can be deleted */
			if (MidNodes[i]==NULL) continue;
			candelete = 1;

			/* This midnode can be deleted, if all of it's remaining links */
			/* are to the endpoints of the father edge. In this case all   */
			/* elements sharing that edge have been visited.			   */
			if (START(MidNodes[i]) == NULL)
			{
				DisposeVertex(theGrid,MYVERTEX(MidNodes[i]));
				DisposeNode(theGrid,MidNodes[i]);
				MidNodes[i] = NULL;
				theEdge = GetEdge(CORNER(theElement,Corner0),CORNER(theElement,Corner1));
				MIDNODE(theEdge) = NULL;
			}
		}
		IFDEBUG(gm,2)
		UserWriteF(" CHANGED ID(newMIDNODE)=%d\n",ID(MidNodes[i]));    
		ENDDEBUG
	}

	IFDEBUG(gm,2)
	UserWriteF("\n");    
	ENDDEBUG

	#ifdef __THREEDIM__
		SideNodes = theElementContext+CORNERS_OF_ELEM(theElement)+EDGES_OF_ELEM(theElement);
		for (i=0; i<SIDES_OF_ELEM(theElement); i++)
		{
			/* no side nodes for triangular sides yet */
			if (CORNERS_OF_SIDE(theElement,i) == 3) continue;

			toDelete = 0;
			/* if side node exists and is not needed, try to delete it */
			if (SideNodes[i]!=NULL) {
				if (DIM==3 && TAG(theElement)==HEXAHEDRON && MARKCLASS(theElement)==GREEN) {
					theNeighbor = NBELEM(theElement,i);

					if (theNeighbor!=NULL) {
						if (MARKCLASS(theNeighbor)==GREEN || MARKCLASS(theNeighbor)==YELLOW)
							toDelete = 1;
						else {
							for (j=0; j<SIDES_OF_ELEM(theNeighbor); j++) {
								if (NBELEM(theNeighbor,j) == theElement) break;			
							}
							assert(j<SIDES_OF_ELEM(theNeighbor));

							if (MARK2RULEADR(theNeighbor,MARK(theNeighbor))->sonandnode[EDGES_OF_ELEM(theNeighbor)+j][0]==-1) {
								toDelete = 1;
							}
						}
					}
					else
						assert(SIDE(theElement,i)!=NULL);
				}
				else if (MARK2RULEADR(theElement,Mark)->sonandnode[i][0]==-1) {
					toDelete = 1;
				}
			}

			toCreate = 0;
			/* if side node does not exist and is needed, allocate it */
			if (SideNodes[i]==NULL) {
				if (DIM==3 && TAG(theElement)==HEXAHEDRON && MARKCLASS(theElement)==GREEN) {
					theNeighbor = NBELEM(theElement,i);

					if (theNeighbor!=NULL) {
						if (MARKCLASS(theNeighbor)!=GREEN && MARKCLASS(theNeighbor)!=YELLOW) {
							for (j=0; j<SIDES_OF_ELEM(theNeighbor); j++) {
								if (NBELEM(theNeighbor,j) == theElement) break;			
							}
							assert(j<SIDES_OF_ELEM(theNeighbor));
							if (MARK2RULEADR(theNeighbor,MARK(theNeighbor))->sonandnode[EDGES_OF_ELEM(theNeighbor)+j][0]!=-1 &&
								(MARK2RULEADR(theNeighbor,REFINE(theNeighbor))->sonandnode[EDGES_OF_ELEM(theNeighbor)+j][0]==-1 ||
								 USED(theNeighbor)==0))
									toCreate = 1;
						}
					}
				}
				else if (MARK2RULEADR(theElement,Mark)->sonandnode[i][0]!=-1) {
					toCreate = 1;
				}
			}
			IFDEBUG(gm,2)
			if (SideNodes[i] == NULL)
			  UserWriteF("    SideNode[%d]: delete=%d create=%d old=%x",i,toDelete,toCreate,SideNodes[i]);
			else
			  UserWriteF("    SideNode[%d]: delete=%d create=%d old=%x oldID=%d",i,toDelete,toCreate,SideNodes[i],ID(SideNodes[i]));
			ENDDEBUG

			if (toDelete)
			{
				/* this node has no links any more delete it */
				if (START(SideNodes[i])==NULL)
				{
					DisposeVertex(theGrid,MYVERTEX(SideNodes[i]));
					DisposeNode(theGrid,SideNodes[i]);
					SideNodes[i] = NULL;
				}
			}

			if (toCreate)
			{
				theNeighbor = NBELEM(theElement,i);
				if (theNeighbor != NULL)
				  {
					IFDEBUG(gm,3)
					  UserWriteF("    ID(theNeighbor)=%d nbadr=%x:\n",ID(theNeighbor),theNeighbor);
					ENDDEBUG

					  /* check whether node exists already */
					  Refine = REFINE(theNeighbor);
				  }

				if (theNeighbor!=NULL && DIM==3 && TAG(theNeighbor)==HEXAHEDRON && 
					MARKCLASS(theNeighbor)==GREEN && USED(theNeighbor)==0) {

					IFDEBUG(gm,3)
					UserWriteF("    	Serching for side node allocated by green neighbor:\n");
					ENDDEBUG

					assert(SideNodes[i] == NULL);
					/* get the side node */
					theNode0 = theElementContext[EDGE_OF_SIDE(theElement,i,0)+CORNERS_OF_ELEM(theElement)];
					theNode1 = theElementContext[EDGE_OF_SIDE(theElement,i,2)+CORNERS_OF_ELEM(theElement)];
					for (theLink0=START(theNode0); theLink0!=NULL; theLink0=NEXT(theLink0)) {
						for (theLink1=START(theNode1); theLink1!=NULL; theLink1=NEXT(theLink1)) 
							if (NBNODE(theLink0) == NBNODE(theLink1)) {
								
								IFDEBUG(gm,3)
								UserWriteF("    		possible node ID=%d\n",ID(NBNODE(theLink0)));
								ENDDEBUG

								l = 0;
								for (theLink=START(NBNODE(theLink0)); theLink!=NULL; theLink=NEXT(theLink)) 
									for (k=0; k<CORNERS_OF_SIDE(theElement,i); k++) 
										if (NBNODE(theLink) == SONNODE(CORNER(theElement,CORNER_OF_SIDE(theElement,i,k))))
											l = 1;
								if (l == 0) {
									SideNodes[i] = NBNODE(theLink0);
									IFDEBUG(gm,3)
									UserWriteF("    		FOUND node ID=%d\n",ID(NBNODE(theLink)));
									continue;
									ENDDEBUG
									break;
								}
							}
						if (SideNodes[i] != NULL) {
							IFDEBUG(gm,3)
							continue;
							ENDDEBUG
							break;
						}
					}
					assert(SideNodes[i] != NULL);
				}
				else if (theNeighbor!=NULL &&  Refine>0 && !REF_TYPE_CHANGES(theNeighbor) &&
						 (DIM!=3 || TAG(theNeighbor)!=HEXAHEDRON || MARKCLASS(theNeighbor)!=GREEN))
				{
					/* get the side node */
					nbrule = MARK2RULEADR(theNeighbor,Refine);
					for (j=0; j<SIDES_OF_ELEM(theNeighbor); j++) {
						if (NBELEM(theNeighbor,j) == theElement) break;			
					}
					assert(j<SIDES_OF_ELEM(theNeighbor));
					GetSons(theNeighbor,NeighborSonList);
					theNeighbor = NeighborSonList[nbrule->sonandnode[EDGES_OF_ELEM(theNeighbor)+j][0]];
					SideNodes[i] = CORNER(theNeighbor,nbrule->sonandnode[EDGES_OF_ELEM(theNeighbor)+j][1]);
				}
				else
				{
					/* allocate the sidenode */
					/* TODO: this is dependent on element numbering */
					switch (i)
					{
						case (0):
						Node0 = theElementContext[8];   Node1 = theElementContext[10];
						break;
						case (1):
						Node0 = theElementContext[8];   Node1 = theElementContext[16];
						break;
						case (2):
						Node0 = theElementContext[9];   Node1 = theElementContext[17];
						break;
						case (3):
						Node0 = theElementContext[10];   Node1 = theElementContext[18];
						break;
						case (4):
						Node0 = theElementContext[11];   Node1 = theElementContext[19];
						break;
						case (5):
						Node0 = theElementContext[16];   Node1 = theElementContext[18];
						break;
					}
					assert ((Node0!=NULL)&&(Node1!=NULL));

					if ((SideNodes[i] = CreateSideNode(theGrid,theElement,Node0,Node1,NULL)) == NULL) return(GM_FATAL);;
				}
				assert (SideNodes[i]!=NULL);
				for (j=0; j<EDGES_OF_SIDE(theElement,i); j++)
				{
					fatherEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,EDGE_OF_SIDE(theElement,i,j),0)),CORNER(theElement,CORNER_OF_EDGE(theElement,EDGE_OF_SIDE(theElement,i,j),1)));
					Node0 = MIDNODE(fatherEdge);
					assert (Node0 != NULL);
					if ((theEdge=CreateEdge(theGrid,Node0,SideNodes[i]))==NULL) return(1);
				}
			}
			IFDEBUG(gm,2)
			if (SideNodes[i] != NULL) 
			  UserWriteF(" new=%x newID=%d\n",SideNodes[i],ID(SideNodes[i]));
			else
			  UserWriteF(" new=%x\n",SideNodes[i]);
			ENDDEBUG
		}
	#endif
	
	/* allocate/remove center node */
	CenterNode = theElementContext+CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement);

	toDelete = 0;
	if (CenterNode[0] != NULL) {
		if (DIM==3 && TAG(theElement)==HEXAHEDRON && MARKCLASS(theElement)==GREEN) {
			/* do nothing */
		}
		else if (MARK2RULEADR(theElement,Mark)->sonandnode[CENTER_NODE_INDEX(theElement)][0] == NO_CENTER_NODE) {
				toDelete = 1;
		}
	}
	toCreate = 0;
	if (CenterNode[0] == NULL) {
		if (DIM==3 && TAG(theElement)==HEXAHEDRON && MARKCLASS(theElement) == GREEN) {
			toCreate = 1;
		}
		else if (MARK2RULEADR(theElement,Mark)->sonandnode[CENTER_NODE_INDEX(theElement)][0] != NO_CENTER_NODE){
				toCreate = 1;
		}
	}

	if (toDelete)
	{
		/* existing center node has to be removed */
		DisposeVertex(theGrid,MYVERTEX(CenterNode[0]));
		DisposeNode(theGrid,CenterNode[0]);
		CenterNode[0] = NULL;
	}

	if (toCreate)
	{
		/* we need an interior node */
		if (DIM == 2)
		{
			CenterNode[0] = CreateNode(theGrid,NULL);
			if (CenterNode[0]==NULL) return(GM_ERROR);

			CenterVertex = CreateInnerVertex(theGrid,NULL);
			if (CenterVertex==NULL) return(GM_ERROR);

			SETUSED(theNode,1);
			theGrid->status |= 1;
			sx = sy = 0.0;
			for (j=0; j<CORNERS_OF_ELEM(theElement); j++)
			{
				sx += XC(MYVERTEX(CORNER(theElement,j)));
				sy += YC(MYVERTEX(CORNER(theElement,j)));
			}
			XC(CenterVertex) = 0.25*sx;
			YC(CenterVertex) = 0.25*sy;
			XI(CenterVertex) = 0.0;
			ETA(CenterVertex) = 0.0;
			VFATHER(CenterVertex) = theElement;
			NFATHER(CenterNode[0]) = NULL;
			MYVERTEX(CenterNode[0]) = CenterVertex;
			TOPNODE(CenterVertex) = CenterNode[0];
		}

		/* there is no center node, but we need one: so allocate it */
		if (DIM == 3)
		{
			CenterNode[0] = CreateNode(theGrid,NULL);
			if (CenterNode[0] == NULL) return (GM_ERROR);

			/* allocate center vertex and init local and global position */
			CenterVertex = CreateInnerVertex(theGrid,NULL); 
			if (CenterVertex == NULL) return (GM_ERROR);
			x = CVECT(CenterVertex);
			xi = LCVECT(CenterVertex);
			V3_CLEAR(x)
			V3_CLEAR(xi)
			for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
			{
				V3_LINCOMB(1.0, x, 1.0, CVECT(MYVERTEX(CORNER(theElement,i))), x);
				V3_LINCOMB(1.0, xi, 1.0, LOCAL_COORD_OF_ELEM(theElement,i), xi)
			}

			V3_SCALE(1.0/CORNERS_OF_ELEM(theElement), x)
			V3_SCALE(1.0/CORNERS_OF_ELEM(theElement), xi)
			/* TODO: delete this
			V3_SCALE(0.25, x)
			V3_SCALE(0.25, xi) */
			
			/* init ptrs */
			VFATHER(CenterVertex) = theElement;
			TOPNODE(CenterVertex) = CenterNode[0];
			MYVERTEX(CenterNode[0]) = CenterVertex;
		}
	}

	return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  UnrefineElement												*/
/*																			*/
/* Purpose:   remove previous refinement for an element	and all sonelements	*/
/*			  recursively, deletes:											*/
/*			  (i)	 all connections                                        */
/*			  (ii)	 all interior nodes and edges are deleted				*/
/*			  (iii)	 sons are deleted and references to sons reset to NULL	*/
/*																			*/
/* Param:	  GRID *theGrid: grid level of sons of theElement				*/
/*			  ELEMENT *theElement: element to refine						*/
/*			  ELEMENTCONTEXT *theContext: current context of element		*/
/*																			*/
/* return:	  none															*/
/*																			*/
/****************************************************************************/

static INT UnrefineElement (GRID *theGrid, ELEMENT *theElement, NODE **theElementContext)
{
	int i,j,s;
	EDGE *theEdge;
	NODE *CenterNode;
	REFRULE *rule;					/* current refinement rule of theElement*/
	ELEMENTCONTEXT sonContext;
	ELEMENT *theSon,*SonList[MAX_SONS];
	EDGEDATA *edata;

	/* something to do ? */
	if ((REFINE(theElement)==0)||(theGrid==NULL)) return(GM_OK);

	if (GetSons(theElement,SonList)!=GM_OK) return(GM_ERROR);
	
	for (s=0; s<NSONS(theElement); s++)
	{
		theSon = SonList[s];
		SETMARK(theSon,NO_REFINEMENT);
		if (REFINE(theSon)>0)
		{
				GetCurrentContext(theSon,sonContext);
				if (UnrefineElement(theGrid->finer,theSon,sonContext)) return (1);
				UpdateContext(theGrid->finer,theSon,sonContext);
		}
	}

	/* remove connections in neighborhood of sons */
	for (i=0; i<NSONS(theElement); i++)
		DisposeConnectionsInNeighborhood(theGrid,SonList[i]);

	/* remove son elements */
	IFDEBUG(gm,1)
	if (DIM!=3 || TAG(theElement)!=HEXAHEDRON || REFINECLASS(theElement)!=GREEN) {
		rule = MARK2RULEADR(theElement,REFINE(theElement));
		if (NSONS(theElement) != rule->nsons) 
			UserWriteF("ERROR: NSONS=%d but rule.sons=%d\n",NSONS(theElement),rule->nsons);
	}
	ENDDEBUG
	for (s=0; s<NSONS(theElement); s++)
	{
		/* dispose all edges not needed any more */
		if (DisposeEdgesFromElement(theGrid,SonList[s])) return (1);
		DisposeElement(theGrid,SonList[s]);
	}


	SETNSONS(theElement,0);
	SET_SON(theElement,0,NULL);

	return (0);
}


/****************************************************************************/
/*																			*/
/* Function:  RefineGreenElement											*/
/*																			*/
/* Purpose:   refine an element without context     	 					*/
/*			  (i)	 corner and midnodes are already allocated				*/
/*			  (ii)	 edges between corner and midnodes are ok				*/
/*			  (iii)  create interior nodes and edges						*/
/*			  (iv)	 create sons and set references to sons 				*/
/*																			*/
/* Param:	  GRID *theGrid: grid level of sons of theElement				*/
/*			  ELEMENT *theElement: element to refine						*/
/*																			*/
/* return:	  INT 0: ok 													*/
/*			  INT 1: fatal memory error 									*/
/*																			*/
/****************************************************************************/
		
static int RefineGreenElement (GRID *theGrid, ELEMENT *theElement, NODE **theContext)
{
	struct greensondata {
		short		tag;
		short		bdy;
		NODE		*corners[MAX_CORNERS_OF_ELEM];
		int			nb[MAX_SIDES_OF_ELEM];
		ELEMENT		*theSon;
	}; 
	typedef struct greensondata		GREENSONDATA;

	GREENSONDATA sons[MAX_GREEN_SONS];

	NODE *theNode, *theEdgeNode, *theNode0, *theNode1;
	NODE *theSideNodes[8];
	VERTEX *CenterVertex, *myvertex;
	VSEGMENT *vs;
	EDGE *theEdge;
	ELEMENT *theSon;
	ELEMENT *NbElement;
	ELEMENT *NbSonList[MAX_SONS];
	ELEMENT *NbSideSons[5];
	ELEMENTSIDE *oldSide, *newSide;
	NODE **CenterNode;
	REFRULE *NbRule;
	SONDATA *sdata2;
	COORD *x,*xi;
	int i,j,k,l,m,n,o,p,q,r,s,s2,t,found,points;
	int NbrSide,side,nbside,NbSonIndex,nelem,nedges,node,node0,nsi;
	int bdy,edge, sides[4], side0, side1;
	int tetNode0, tetNode1, tetNode2,
		tetSideNode0Node1, tetSideNode0Node2, tetSideNode1Node2,
		pyrSide, pyrNode0, pyrNode1, pyrNode2, pyrNode3,
		pyrSideNode0Node1, pyrSideNode1Node2, pyrSideNode2Node3, pyrSideNode0Node3;
	int elementsSide0[5], elementsSide1[5];

	IFDEBUG(gm,1)
	UserWriteF("RefineGreenElement(): ELEMENT ID=%d\n",ID(theElement)); 
	ENDDEBUG

	/* init son data array */
	for (i=0; i<MAX_GREEN_SONS; i++) {
		sons[i].tag = -1;
		sons[i].bdy = -1;
		for (j=0; j<MAX_CORNERS_OF_ELEM; j++)	sons[i].corners[j] = NULL;
		for (j=0; j<MAX_SIDES_OF_ELEM; j++)		sons[i].nb[j] = -1;
		sons[i].theSon = NULL;
	}

	IFDEBUG(gm,2)
	UserWriteF("         Element ID=%d actual CONTEXT is:\n",ID(theElement));
	for (i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++) UserWriteF(" %3d",i);
	UserWrite("\n");
	for (i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++) 
	  if (theContext[i] != NULL)
		UserWriteF(" %3d",ID(theContext[i]));
	UserWrite("\n");
	ENDDEBUG
	IFDEBUG(gm,3)
	for (i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++) {
		if (theContext[i] == NULL) continue;
		if (NDOBJ != OBJT(theContext[i])) UserWriteF(" ERROR NO NDOBJ(5) OBJT(i=%d)=%d ID=%d adr=%x\n",i,OBJT(theContext[i]),ID(theContext[i]),theContext[i]);
	}
	ENDDEBUG

	/* create inner edges */
	for (i=0; i<(MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM-1); i++) {
		if (theContext[i] == NULL) continue;
		if ((theEdge=CreateEdge(theGrid,theContext[i],theContext[CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement)]))==NULL) return(GM_FATAL);
	}

	/* create edges on edges */
	for (i=0; i<EDGES_OF_ELEM(theElement); i++) {
		if (theContext[i+CORNERS_OF_ELEM(theElement)] == NULL) {
			/* no midnode exists, create one corner corner edge */
			if ((theEdge=CreateEdge(theGrid,theContext[CORNER_OF_EDGE(theElement,i,0)],theContext[CORNER_OF_EDGE(theElement,i,1)]))==NULL) return(GM_FATAL);
		}
		else
		{
			/* midnode exists, create two corner midnode edges */
			if ((theEdge=CreateEdge(theGrid,theContext[i+CORNERS_OF_ELEM(theElement)],theContext[CORNER_OF_EDGE(theElement,i,0)]))==NULL) return(GM_FATAL);
			if ((theEdge=CreateEdge(theGrid,theContext[i+CORNERS_OF_ELEM(theElement)],theContext[CORNER_OF_EDGE(theElement,i,1)]))==NULL) return(GM_FATAL);
		}
	}

	/* init indices for son elements */
	/* outer side for tetrahedra is side 0 */
	tetNode0 = element_descriptors[TETRAHEDRON]->corner_of_side[0][0];
	tetNode1 = element_descriptors[TETRAHEDRON]->corner_of_side[0][1];
	tetNode2 = element_descriptors[TETRAHEDRON]->corner_of_side[0][2];

	for (i=0; i<element_descriptors[TETRAHEDRON]->edges_of_elem; i++)
		if(element_descriptors[TETRAHEDRON]->corner_of_edge[i][0]==tetNode0 && element_descriptors[TETRAHEDRON]->corner_of_edge[i][1]==tetNode1 ||
		   element_descriptors[TETRAHEDRON]->corner_of_edge[i][1]==tetNode0 && element_descriptors[TETRAHEDRON]->corner_of_edge[i][0]==tetNode1)
			break;
	for (j=0; j<element_descriptors[TETRAHEDRON]->sides_of_elem; j++)
		if (element_descriptors[TETRAHEDRON]->side_with_edge[i][j] != 0)
			break;
	tetSideNode0Node1 = element_descriptors[TETRAHEDRON]->side_with_edge[i][j];
	
	for (i=0; i<element_descriptors[TETRAHEDRON]->edges_of_elem; i++)
		if(element_descriptors[TETRAHEDRON]->corner_of_edge[i][0]==tetNode0 && element_descriptors[TETRAHEDRON]->corner_of_edge[i][1]==tetNode2 ||
		   element_descriptors[TETRAHEDRON]->corner_of_edge[i][1]==tetNode0 && element_descriptors[TETRAHEDRON]->corner_of_edge[i][0]==tetNode2)
			break;
	for (j=0; j<element_descriptors[TETRAHEDRON]->sides_of_elem; j++)
		if (element_descriptors[TETRAHEDRON]->side_with_edge[i][j] != 0)
			break;
	tetSideNode0Node2 = element_descriptors[TETRAHEDRON]->side_with_edge[i][j];
	
	for (i=0; i<element_descriptors[TETRAHEDRON]->edges_of_elem; i++)
		if(element_descriptors[TETRAHEDRON]->corner_of_edge[i][0]==tetNode1 && element_descriptors[TETRAHEDRON]->corner_of_edge[i][1]==tetNode2 ||
		   element_descriptors[TETRAHEDRON]->corner_of_edge[i][1]==tetNode1 && element_descriptors[TETRAHEDRON]->corner_of_edge[i][0]==tetNode2)
			break;
	for (j=0; j<element_descriptors[TETRAHEDRON]->sides_of_elem; j++)
		if (element_descriptors[TETRAHEDRON]->side_with_edge[i][j] != 0)
			break;
	tetSideNode1Node2 = element_descriptors[TETRAHEDRON]->side_with_edge[i][j];

	/* outer side for pyramid has 4 corners */
	for (i=0; i<element_descriptors[PYRAMID]->sides_of_elem; i++)
		if (element_descriptors[PYRAMID]->corners_of_side[i] == 4)
			break;
	pyrSide = i;
	pyrNode0 = element_descriptors[PYRAMID]->corner_of_side[i][0];
	pyrNode1 = element_descriptors[PYRAMID]->corner_of_side[i][1];
	pyrNode2 = element_descriptors[PYRAMID]->corner_of_side[i][2];
	pyrNode3 = element_descriptors[PYRAMID]->corner_of_side[i][3];

	for (k=0; k<element_descriptors[PYRAMID]->edges_of_elem; k++)
		if(element_descriptors[PYRAMID]->corner_of_edge[k][0]==pyrNode0 && element_descriptors[PYRAMID]->corner_of_edge[k][1]==pyrNode1 ||
		   element_descriptors[PYRAMID]->corner_of_edge[k][1]==pyrNode0 && element_descriptors[PYRAMID]->corner_of_edge[k][0]==pyrNode1)
			break;
	for (j=0; j<element_descriptors[PYRAMID]->sides_of_elem; j++)
		if (element_descriptors[PYRAMID]->side_with_edge[k][j] != i)
			break;
	pyrSideNode0Node1 = element_descriptors[PYRAMID]->side_with_edge[k][j];

	for (k=0; k<element_descriptors[PYRAMID]->edges_of_elem; k++)
		if(element_descriptors[PYRAMID]->corner_of_edge[k][0]==pyrNode1 && element_descriptors[PYRAMID]->corner_of_edge[k][1]==pyrNode2 ||
		   element_descriptors[PYRAMID]->corner_of_edge[k][1]==pyrNode1 && element_descriptors[PYRAMID]->corner_of_edge[k][0]==pyrNode2)
			break;
	for (j=0; j<element_descriptors[PYRAMID]->sides_of_elem; j++)
		if (element_descriptors[PYRAMID]->side_with_edge[k][j] != i)
			break;
	pyrSideNode1Node2 = element_descriptors[PYRAMID]->side_with_edge[k][j];

	for (k=0; k<element_descriptors[PYRAMID]->edges_of_elem; k++)
		if(element_descriptors[PYRAMID]->corner_of_edge[k][0]==pyrNode2 && element_descriptors[PYRAMID]->corner_of_edge[k][1]==pyrNode3 ||
		   element_descriptors[PYRAMID]->corner_of_edge[k][1]==pyrNode2 && element_descriptors[PYRAMID]->corner_of_edge[k][0]==pyrNode3)
			break;
	for (j=0; j<element_descriptors[PYRAMID]->sides_of_elem; j++)
		if (element_descriptors[PYRAMID]->side_with_edge[k][j] != i)
			break;
	pyrSideNode2Node3 = element_descriptors[PYRAMID]->side_with_edge[k][j];

	for (k=0; k<element_descriptors[PYRAMID]->edges_of_elem; k++)
		if(element_descriptors[PYRAMID]->corner_of_edge[k][0]==pyrNode3 && element_descriptors[PYRAMID]->corner_of_edge[k][1]==pyrNode0 ||
		   element_descriptors[PYRAMID]->corner_of_edge[k][1]==pyrNode3 && element_descriptors[PYRAMID]->corner_of_edge[k][0]==pyrNode0)
			break;
	for (j=0; j<element_descriptors[PYRAMID]->sides_of_elem; j++)
		if (element_descriptors[PYRAMID]->side_with_edge[k][j] != i)
			break;
	pyrSideNode0Node3 = element_descriptors[PYRAMID]->side_with_edge[k][j];

	/* create edges on inner of sides, create son elements and connect them */
	for (i=0; i<SIDES_OF_ELEM(theElement); i++) {
		theNode = theContext[CORNERS_OF_ELEM(theElement)+EDGES_OF_ELEM(theElement)+i];
		nedges = EDGES_OF_SIDE(theElement,i);

		bdy = 0;
		if (OBJT(theElement) == BEOBJ && SIDE(theElement,i)!= NULL)
			bdy = 1;
		nelem = 5*i;
		for (j=nelem; j<(nelem+5); j++)
			sons[j].bdy = bdy;

		k = 0;
		for (j=0; j<EDGES_OF_SIDE(theElement,i); j++) {
			edge = element_descriptors[HEXAHEDRON]->edge_of_side[i][j];
			for (l=0; l<MAX_SIDES_OF_ELEM; l++)
				if (element_descriptors[HEXAHEDRON]->side_with_edge[edge][l] != i) {
					sides[k++] = element_descriptors[HEXAHEDRON]->side_with_edge[edge][l]+MAX_GREEN_SONS;
					break;
				}
			assert(l<2);
		}
		
		k = 0;
		for (j=0; j<nedges; j++) {
			theSideNodes[2*j] = theContext[CORNER_OF_SIDE(theElement,i,j)];
			theSideNodes[2*j+1] = theContext[CORNERS_OF_ELEM(theElement)+EDGE_OF_SIDE(theElement,i,j)];
			if (theSideNodes[2*j+1] != NULL) k++;
		}

		IFDEBUG(gm,2)
		UserWriteF("    SIDE %d has %d nodes and sidenode=%x\n",i,k,theNode);   
		ENDDEBUG
		if (theNode == NULL) {
			switch (k) {
				case 0:
					sons[nelem].tag = PYRAMID;
					sons[nelem].corners[pyrNode0] = theSideNodes[0];
					sons[nelem].corners[pyrNode1] = theSideNodes[2];
					sons[nelem].corners[pyrNode2] = theSideNodes[4];
					sons[nelem].corners[pyrNode3] = theSideNodes[6]; 

					sons[nelem].nb[pyrSideNode0Node1] = sides[0];
					sons[nelem].nb[pyrSideNode1Node2] = sides[1];
					sons[nelem].nb[pyrSideNode2Node3] = sides[2];
					sons[nelem].nb[pyrSideNode0Node3] = sides[3];
					nelem++;
					
					break;
				case 1:
					for (j=0; j<nedges; j++) {
						node0 = 2*j+1; 
						if (theSideNodes[node0] != NULL) {
							if ((theEdge = CreateEdge(theGrid,theSideNodes[node0],theSideNodes[(node0+3)%(2*nedges)])) == NULL) return(GM_FATAL);
							if ((theEdge = CreateEdge(theGrid,theSideNodes[node0],theSideNodes[(node0+5)%(2*nedges)])) == NULL) return(GM_FATAL);

							/* define the son corners and inner side relations */
							sons[nelem].tag = TETRAHEDRON;
							sons[nelem].corners[tetNode0] = theSideNodes[node0];
							sons[nelem].corners[tetNode1] = theSideNodes[(node0+1)%(2*nedges)];
							sons[nelem].corners[tetNode2] = theSideNodes[(node0+3)%(2*nedges)]; 

							sons[nelem].nb[tetSideNode0Node1] = sides[j];
							sons[nelem].nb[tetSideNode1Node2] = sides[(j+1)%nedges];
							sons[nelem].nb[tetSideNode0Node2] = nelem+2;
							nelem++;

							sons[nelem].tag = TETRAHEDRON;
							sons[nelem].corners[tetNode0] = theSideNodes[node0];
							sons[nelem].corners[tetNode1] = theSideNodes[(node0+5)%(2*nedges)];
							sons[nelem].corners[tetNode2] = theSideNodes[(node0+7)%(2*nedges)]; 

							sons[nelem].nb[tetSideNode0Node1] = nelem+1;
							sons[nelem].nb[tetSideNode1Node2] = sides[(j+3)%nedges];
							sons[nelem].nb[tetSideNode0Node2] = sides[j];
							nelem++;

							sons[nelem].tag = TETRAHEDRON;
							sons[nelem].corners[tetNode0] = theSideNodes[node0];
							sons[nelem].corners[tetNode1] = theSideNodes[(node0+3)%(2*nedges)];
							sons[nelem].corners[tetNode2] = theSideNodes[(node0+5)%(2*nedges)]; 

							sons[nelem].nb[tetSideNode0Node1] = nelem-2;
							sons[nelem].nb[tetSideNode1Node2] = sides[(j+2)%nedges];
							sons[nelem].nb[tetSideNode0Node2] = nelem-1;
							nelem++;

							break;
						}
					}
					break;
				case 2:
					/* two cases: sidenodes are not on neighboring edges OR are on neighboring edges */
					for (j=0; j<nedges; j++) {
						node0 = 2*j+1; 
						if (theSideNodes[node0] != NULL) 
							break;
					}
					if (theSideNodes[(node0+6)%(2*nedges)] != NULL) {
						node0 = (node0+6)%(2*nedges);
						j = (j+3)%nedges;
					}
					if (theSideNodes[(node0+4)%(2*nedges)] == NULL) {
						if ((theEdge = CreateEdge(theGrid,theSideNodes[node0],theSideNodes[(node0+2)%(2*nedges)])) == NULL) return(GM_FATAL);
						if ((theEdge = CreateEdge(theGrid,theSideNodes[(node0+2)%(2*nedges)],theSideNodes[(node0+5)%(2*nedges)])) == NULL) return(GM_FATAL);
						if ((theEdge = CreateEdge(theGrid,theSideNodes[node0],theSideNodes[(node0+5)%(2*nedges)])) == NULL) return(GM_FATAL);

						sons[nelem].tag = TETRAHEDRON;
						sons[nelem].corners[tetNode0] = theSideNodes[node0];
						sons[nelem].corners[tetNode1] = theSideNodes[(node0+1)%(2*nedges)];
						sons[nelem].corners[tetNode2] = theSideNodes[(node0+2)%(2*nedges)]; 

						sons[nelem].nb[tetSideNode0Node1] = sides[(j)%nedges];
						sons[nelem].nb[tetSideNode1Node2] = sides[(j+1)%nedges];
						sons[nelem].nb[tetSideNode0Node2] = nelem+3;
						nelem++;

						sons[nelem].tag = TETRAHEDRON;
						sons[nelem].corners[tetNode0] = theSideNodes[node0];
						sons[nelem].corners[tetNode1] = theSideNodes[(node0+5)%(2*nedges)];
						sons[nelem].corners[tetNode2] = theSideNodes[(node0+7)%(2*nedges)]; 

						sons[nelem].nb[tetSideNode0Node1] = nelem+2;
						sons[nelem].nb[tetSideNode1Node2] = sides[(j+3)%nedges];
						sons[nelem].nb[tetSideNode0Node2] = sides[(j)%nedges];
						nelem++;

						sons[nelem].tag = TETRAHEDRON;
						sons[nelem].corners[tetNode0] = theSideNodes[(node0+2)%(2*nedges)];
						sons[nelem].corners[tetNode1] = theSideNodes[(node0+3)%(2*nedges)];
						sons[nelem].corners[tetNode2] = theSideNodes[(node0+5)%(2*nedges)]; 

						sons[nelem].nb[tetSideNode0Node1] = sides[(j+1)%nedges];
						sons[nelem].nb[tetSideNode1Node2] = sides[(j+2)%nedges];
						sons[nelem].nb[tetSideNode0Node2] = nelem+1;
						nelem++;

						sons[nelem].tag = TETRAHEDRON;
						sons[nelem].corners[tetNode0] = theSideNodes[node0];
						sons[nelem].corners[tetNode1] = theSideNodes[(node0+2)%(2*nedges)];
						sons[nelem].corners[tetNode2] = theSideNodes[(node0+5)%(2*nedges)]; 

						sons[nelem].nb[tetSideNode0Node1] = nelem-3;
						sons[nelem].nb[tetSideNode1Node2] = nelem-1;
						sons[nelem].nb[tetSideNode0Node2] = nelem-2;
						nelem++;
					}
					else {
						if ((theEdge = CreateEdge(theGrid,theSideNodes[node0],theSideNodes[(node0+4)%(2*nedges)])) == NULL) return(GM_FATAL);

						sons[nelem].tag = PYRAMID;
						sons[nelem].corners[pyrNode0] = theSideNodes[node0];
						sons[nelem].corners[pyrNode1] = theSideNodes[(node0+1)%(2*nedges)];
						sons[nelem].corners[pyrNode2] = theSideNodes[(node0+3)%(2*nedges)];
						sons[nelem].corners[pyrNode3] = theSideNodes[(node0+4)%(2*nedges)]; 

						sons[nelem].nb[pyrSideNode0Node1] = sides[(j)%nedges];
						sons[nelem].nb[pyrSideNode1Node2] = sides[(j+1)%nedges];
						sons[nelem].nb[pyrSideNode2Node3] = sides[(j+2)%nedges];
						sons[nelem].nb[pyrSideNode0Node3] = nelem+1;
						nelem++;

						sons[nelem].tag = PYRAMID;
						sons[nelem].corners[pyrNode0] = theSideNodes[(node0+4)%(2*nedges)];
						sons[nelem].corners[pyrNode1] = theSideNodes[(node0+5)%(2*nedges)];
						sons[nelem].corners[pyrNode2] = theSideNodes[(node0+7)%(2*nedges)];
						sons[nelem].corners[pyrNode3] = theSideNodes[(node0+8)%(2*nedges)]; 

						sons[nelem].nb[pyrSideNode0Node1] = sides[(j+2)%nedges];
						sons[nelem].nb[pyrSideNode1Node2] = sides[(j+3)%nedges];
						sons[nelem].nb[pyrSideNode2Node3] = sides[(j)%nedges];
						sons[nelem].nb[pyrSideNode0Node3] = nelem-1;
						nelem++;
					}
					break;
				case 3:
					for (j=0; j<nedges; j++) {
						node0 = 2*j+1; 
						if (theSideNodes[node0] == NULL) 
							break;
					}
					if ((theEdge = CreateEdge(theGrid,theSideNodes[(node0+2)%(2*nedges)],theSideNodes[(node0+6)%(2*nedges)])) == NULL) return(GM_FATAL);
					if ((theEdge = CreateEdge(theGrid,theSideNodes[(node0+2)%(2*nedges)],theSideNodes[(node0+4)%(2*nedges)])) == NULL) return(GM_FATAL);
					if ((theEdge = CreateEdge(theGrid,theSideNodes[(node0+4)%(2*nedges)],theSideNodes[(node0+6)%(2*nedges)])) == NULL) return(GM_FATAL);
					sons[nelem].tag = PYRAMID;
					sons[nelem].corners[pyrNode0] = theSideNodes[(node0+1)%(2*nedges)];
					sons[nelem].corners[pyrNode1] = theSideNodes[(node0+2)%(2*nedges)];
					sons[nelem].corners[pyrNode2] = theSideNodes[(node0+6)%(2*nedges)];
					sons[nelem].corners[pyrNode3] = theSideNodes[(node0+7)%(2*nedges)]; 

					sons[nelem].nb[pyrSideNode0Node1] = sides[(j+1)%nedges];
					sons[nelem].nb[pyrSideNode1Node2] = nelem+3;
					sons[nelem].nb[pyrSideNode2Node3] = sides[(j+3)%nedges];
					sons[nelem].nb[pyrSideNode0Node3] = sides[(j)%nedges];
					nelem++;

					sons[nelem].tag = TETRAHEDRON;
					sons[nelem].corners[tetNode0] = theSideNodes[(node0+2)%(2*nedges)];
					sons[nelem].corners[tetNode1] = theSideNodes[(node0+3)%(2*nedges)];
					sons[nelem].corners[tetNode2] = theSideNodes[(node0+4)%(2*nedges)]; 

					sons[nelem].nb[tetSideNode0Node1] = sides[(j+1)%nedges];
					sons[nelem].nb[tetSideNode1Node2] = sides[(j+2)%nedges];
					sons[nelem].nb[tetSideNode0Node2] = nelem+2;
					nelem++;

					sons[nelem].tag = TETRAHEDRON;
					sons[nelem].corners[tetNode0] = theSideNodes[(node0+4)%(2*nedges)];
					sons[nelem].corners[tetNode1] = theSideNodes[(node0+5)%(2*nedges)];
					sons[nelem].corners[tetNode2] = theSideNodes[(node0+6)%(2*nedges)]; 

					sons[nelem].nb[tetSideNode0Node1] = sides[(j+2)%nedges];
					sons[nelem].nb[tetSideNode1Node2] = sides[(j+3)%nedges];
					sons[nelem].nb[tetSideNode0Node2] = nelem+1;
					nelem++;

					sons[nelem].tag = TETRAHEDRON;
					sons[nelem].corners[tetNode0] = theSideNodes[(node0+2)%(2*nedges)];
					sons[nelem].corners[tetNode1] = theSideNodes[(node0+4)%(2*nedges)];
					sons[nelem].corners[tetNode2] = theSideNodes[(node0+6)%(2*nedges)]; 

					sons[nelem].nb[tetSideNode0Node1] = nelem-2;
					sons[nelem].nb[tetSideNode1Node2] = nelem-1;
					sons[nelem].nb[tetSideNode0Node2] = nelem-3;
					nelem++;

					break;
				case 4:
					for (j=0; j<nedges; j++) {
						node0 = 2*j+1; 
						if ((theEdge = CreateEdge(theGrid,theSideNodes[node0],theSideNodes[(node0+2)%(2*nedges)])) == NULL) return(GM_FATAL);

						sons[nelem].tag = TETRAHEDRON;
						sons[nelem].corners[tetNode0] = theSideNodes[node0];
						sons[nelem].corners[tetNode1] = theSideNodes[(node0+1)%(2*nedges)];
						sons[nelem].corners[tetNode2] = theSideNodes[(node0+2)%(2*nedges)];

						sons[nelem].nb[tetSideNode0Node1] = sides[(j)%nedges];
						sons[nelem].nb[tetSideNode1Node2] = sides[(j+1)%nedges];
						sons[nelem].nb[tetSideNode0Node2] = nelem+(nedges-j);
						nelem++;
					}

					sons[nelem].tag = PYRAMID;
					sons[nelem].corners[pyrNode0] = theSideNodes[1];
					sons[nelem].corners[pyrNode1] = theSideNodes[3];
					sons[nelem].corners[pyrNode2] = theSideNodes[5];
					sons[nelem].corners[pyrNode3] = theSideNodes[7];

					sons[nelem].nb[pyrSideNode0Node1] = nelem-4;
					sons[nelem].nb[pyrSideNode1Node2] = nelem-3;
					sons[nelem].nb[pyrSideNode2Node3] = nelem-2;
					sons[nelem].nb[pyrSideNode0Node3] = nelem-1;
					nelem++;
					break;

				default:
					return(GM_FATAL);
			}				
		}
		else {
			/* create the four side edges */
			for (j=0; j<nedges; j++) {
				node0 = 2*j+1; 
				if (theSideNodes[node0] == NULL) break;
				if ((theEdge = CreateEdge(theGrid,theNode,theSideNodes[node0])) == NULL) return(GM_FATAL);

				sons[nelem].tag = PYRAMID;
				sons[nelem].corners[pyrNode0] = theSideNodes[node0%(2*nedges)];
				sons[nelem].corners[pyrNode1] = theSideNodes[(node0+1)%(2*nedges)];
				sons[nelem].corners[pyrNode2] = theSideNodes[(node0+2)%(2*nedges)];
				sons[nelem].corners[pyrNode3] = theNode;

				sons[nelem].nb[pyrSideNode0Node1] = sides[(j)%nedges];
				sons[nelem].nb[pyrSideNode1Node2] = sides[(j+1)%nedges];
				if (j == 3)
					sons[nelem].nb[pyrSideNode2Node3] = nelem-3;
				else
					sons[nelem].nb[pyrSideNode2Node3] = nelem+1;
				if (j == 0)
					sons[nelem].nb[pyrSideNode0Node3] = nelem+3;
				else
					sons[nelem].nb[pyrSideNode0Node3] = nelem-1;
				nelem++;

			}
			assert(j==4);
		}
	}

	/* connect elements over edges */
	for (i=0; i<EDGES_OF_ELEM(theElement); i++) {
		side0 = SIDE_WITH_EDGE(theElement,i,0);
		side1 = SIDE_WITH_EDGE(theElement,i,1);

		if (theContext[i+CORNERS_OF_ELEM(theElement)] == NULL) {
			/* two elements share this edge */

			/* get son elements for this edge */
			found = 0;
			for (j=side0*5; j<(side0*5+5); j++) {
				for (k=0; k<MAX_SIDES_OF_ELEM; k++)
					if ((sons[j].nb[k]-MAX_GREEN_SONS)==side1) { 
						found = 1;
						break;
					}
				if (found) break;
			}
			assert(j<side0*5+5);

			found = 0;
			for (l=side1*5; l<side1*5+5; l++) {
				for (m=0; m<MAX_SIDES_OF_ELEM; m++)
					if ((sons[l].nb[m]-MAX_GREEN_SONS)==side0) {
						found = 1;
						break;
					}
				if (found) break;
			}
			assert(j<side1*5+5);

			sons[j].nb[k] = l;
			sons[l].nb[m] = j;
		} 
		else {
			/* four elements share this edge */

			/* get son elements for this edge */
			l = 0;
			for (j=side0*5; j<(side0*5+5); j++) {
				for (k=0; k<MAX_SIDES_OF_ELEM; k++)
					if ((sons[j].nb[k]-MAX_GREEN_SONS)==side1)  
						elementsSide0[l++] = j;
			}
			assert(l==2);

			l = 0;
			for (j=side1*5; j<(side1*5+5); j++) {
				for (m=0; m<MAX_SIDES_OF_ELEM; m++)
					if ((sons[j].nb[m]-MAX_GREEN_SONS)==side0)
						elementsSide1[l++] = j;
			}
			assert(l==2);

			/* determine neighboring elements */
			theNode0 = theContext[CORNERS_OF_ELEM(theElement)+i];
			for (j=0; j<CORNERS_OF_EDGE; j++) {
				theNode1 = theContext[CORNER_OF_EDGE(theElement,i,j)];
				found = 0;
				for (l=0; l<2; l++) {
					for (k=0; k<MAX_CORNERS_OF_ELEM; k++) {
						if (theNode1 == sons[elementsSide0[l]].corners[k]) {
							found = 1;
							break;
						}
					}
					if (found) break;
				}
				assert(k<MAX_CORNERS_OF_ELEM);
				assert(l<2);

				found = 0;
				for (m=0; m<2; m++) {
					for (k=0; k<MAX_CORNERS_OF_ELEM; k++) {
						if (theNode1 == sons[elementsSide1[m]].corners[k]) {
							found = 1;
							break;
						}
					}
					if (found) break;
				}
				assert(k<MAX_CORNERS_OF_ELEM);
				assert(m<2);

				/* init neighbor field */
				for (k=0; k<MAX_SIDES_OF_ELEM; k++)
					if ((sons[elementsSide0[l]].nb[k]-MAX_GREEN_SONS)==side1) 
						break;
				assert(k<MAX_SIDES_OF_ELEM);
				sons[elementsSide0[l]].nb[k] = elementsSide1[m];

				for (k=0; k<MAX_SIDES_OF_ELEM; k++)
					if ((sons[elementsSide1[m]].nb[k]-MAX_GREEN_SONS)==side0)
						break;
				assert(k<MAX_SIDES_OF_ELEM);
				sons[elementsSide1[m]].nb[k] = elementsSide0[l];
			}
		}
	}

	/* create son elements */
	IFDEBUG(gm,1)
	UserWriteF("    Creating SON elements for element ID=%d:\n",ID(theElement));
	ENDDEBUG
	n = 0;
	for (i=0; i<MAX_GREEN_SONS; i++) {
		if (sons[i].tag >= 0) {
			IFDEBUG(gm,2)
			if (i%5 == 0)
				UserWriteF("     SIDE %d:\n",i/5);
			ENDDEBUG
			if (sons[i].bdy == 1) 
				sons[i].theSon = CreateBoundaryElement(theGrid,NULL,sons[i].tag);
			else
				sons[i].theSon = CreateInnerElement(theGrid,NULL,sons[i].tag);
			if (sons[i].theSon==NULL) return(GM_FATAL);

			k = l = 0;
			IFDEBUG(gm,0)
			for (j=0; j<CORNERS_OF_ELEM(sons[i].theSon); j++) 
			  for (m=0; m<CORNERS_OF_ELEM(sons[i].theSon); m++) 
				if (sons[i].corners[j] == NULL || sons[i].corners[m] == NULL)
				  {
					if ((m!=j) && (sons[i].corners[j] == sons[i].corners[m]))
					  UserWriteF("     ERROR: son %d has equivalent corners %d=%d adr=%x adr=%x\n",n,j,m,sons[i].corners[j],sons[i].corners[m]); 
				  }
				else
				  if ((m!=j) && (sons[i].corners[j] == sons[i].corners[m] || 
					  (ID(sons[i].corners[j]) == ID(sons[i].corners[m])))) 
					UserWriteF("     ERROR: son %d has equivalent corners %d=%d  ID=%d ID=%d adr=%x adr=%x\n",n,j,m,ID(sons[i].corners[j]),ID(sons[i].corners[m]),sons[i].corners[j],sons[i].corners[m]); 
			ENDDEBUG

			IFDEBUG(gm,2)
			UserWriteF("      SONS[i=%d] ID=%d: CORNERS ",i,ID(sons[i].theSon)); 
			ENDDEBUG
			for (j=0; j<CORNERS_OF_ELEM(sons[i].theSon); j++) {
				if (sons[i].corners[j] != NULL) 
					SET_CORNER(sons[i].theSon,k++,sons[i].corners[j]);	
				else {
					SET_CORNER(sons[i].theSon,k++,theContext[CORNERS_OF_ELEM(theElement)+CENTER_NODE_INDEX(theElement)]);	
					l ++;
				}
				IFDEBUG(gm,2)
				if (sons[i].corners[j] != NULL)
				  UserWriteF(" %d",ID(sons[i].corners[j])); 
				ENDDEBUG
			}
			IFDEBUG(gm,2)
			UserWriteF("\n"); 
			ENDDEBUG

			assert(k == CORNERS_OF_ELEM(sons[i].theSon));
			assert(l == 1);

			SET_EFATHER(sons[i].theSon,theElement);
			SETECLASS(sons[i].theSon,GREEN);
			SETNSONS(theElement,NSONS(theElement)+1);
			if (i == 0) SET_SON(theElement,0,sons[i].theSon);
			for (s=0; s<SIDES_OF_ELEM(sons[i].theSon); s++) {
				SET_NBELEM(sons[i].theSon,s,NULL);
				if (sons[i].bdy == 1) SET_SIDE(sons[i].theSon,s,NULL);
			}

			for (j=0; j<EDGES_OF_ELEM(sons[i].theSon); j++)
			{
				theEdge = CreateEdge(theGrid,CORNER(sons[i].theSon,CORNER_OF_EDGE(sons[i].theSon,j,0)),CORNER(sons[i].theSon,CORNER_OF_EDGE(sons[i].theSon,j,1)));
				assert(theEdge!=NULL);
				if (NO_OF_ELEM(theEdge)<NO_OF_ELEM_MAX-1)
					INC_NO_OF_ELEM(theEdge);
				else
					return (GM_ERROR);
			}
			n++;
		}
	}
	IFDEBUG(gm,1)
	UserWriteF("    n=%d sons created NSONS=%d\n",n,NSONS(theElement)); 
	ENDDEBUG

	/* translate neighbor information */
	for (i=0; i<MAX_GREEN_SONS; i++) {
		if (sons[i].tag >= 0) {
			k = l = 0;
			IFDEBUG(gm,0)
			for (j=0; j<SIDES_OF_ELEM(sons[i].theSon); j++) {
				for (m=0; m<SIDES_OF_ELEM(sons[i].theSon); m++) {
					if (sons[i].nb[j] == sons[i].nb[m] && (m!=j))
						 UserWriteF("     ERROR: son %d has equivalent neighbors %d=%d  NB=%d\n",n,j,m,sons[i].nb[m]); 
				}
			}
			ENDDEBUG
			for (j=0; j<SIDES_OF_ELEM(sons[i].theSon); j++) {
				if (sons[i].nb[j] != -1)
					SET_NBELEM(sons[i].theSon,k++,sons[sons[i].nb[j]].theSon);
				else {
					l++;
					k++;
				}
			}
			assert(k == SIDES_OF_ELEM(sons[i].theSon));
			assert(l == 1);
		}
	}

	/* init outer side relations of son elements */
	for (i=0; i<SIDES_OF_ELEM(theElement); i++) {
		for (j=0; j<5; j++) {
			if (sons[i*5+j].tag < 0) continue; 
			theSon = sons[i*5+j].theSon;
			side = 0;
			if (sons[i*5+j].tag == PYRAMID) {
				for (k=0; k<element_descriptors[PYRAMID]->sides_of_elem; k++)
					if (element_descriptors[PYRAMID]->corners_of_side[k] == 4)
					break;
				side = k;
			}

			/* connect outer sides of son elements */
			if (sons[i*5+j].bdy == 1) {

				for (k=0; k<SIDES_OF_ELEM(sons[i*5+j].theSon); k++)
					SET_SIDE(sons[i*5+j].theSon,k,NULL);

				/* search boundary side */
				oldSide = SIDE(theElement,i);
				assert(oldSide != NULL);

				newSide = CreateElementSide(theGrid);
				if (newSide == NULL) return(GM_FATAL);
				SET_SIDE(sons[i*5+j].theSon,0,newSide);
				ES_PATCH(newSide) = ES_PATCH(oldSide);

				for (k=0; k<CORNERS_OF_SIDE(sons[i*5+j].theSon,side); k++) {
					node = CORNER_OF_SIDE(sons[i*5+j].theSon,side,k);
					for (l=0; l<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; l++)
						if (theContext[l] == CORNER(sons[i*5+j].theSon,node)) break;
					assert(l<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM);

					if (l<CORNERS_OF_ELEM(theElement)) { /* a corner of theElement */
						nsi = CORNER_OF_SIDE_INV(theElement,i,l); /* position of corner in side numeration */
						assert(nsi != -1);
						PARAM(newSide,k,0) = PARAM(oldSide,nsi,0);
						PARAM(newSide,k,1) = PARAM(oldSide,nsi,1);
					}
					else { /* a midpoint of an edge or side of theElement */
						myvertex = MYVERTEX(theContext[l]);
						assert(myvertex!=NULL);
						assert (OBJT(myvertex) == BVOBJ);
						assert (VSEG(myvertex) != NULL);

						/* find common boundary segment */
						for( vs=VSEG(myvertex); vs!=NULL; vs = NEXTSEG(vs) ) {
								if (VS_PATCH(vs) == ES_PATCH(oldSide)) break;
						}
						assert(vs!=NULL);
						PARAM(newSide,k,0) =  LAMBDA(vs,0);
						PARAM(newSide,k,1) =  LAMBDA(vs,1);
					}
				}
			}
			else {
				/* search neighboring element */
				NbElement = NBELEM(theElement,i);
				assert(NbElement != NULL);

				for (l=0; l<SIDES_OF_ELEM(NbElement); l++)
					if (NBELEM(NbElement,l) == theElement)
						break;
				assert(l<SIDES_OF_ELEM(NbElement));
				IFDEBUG(gm,2)
				UserWriteF("     SIDE=%d NBSIDE=%d ID(theElement)=%d ID(NbElement)=%d\n",i,l,ID(theElement),ID(NbElement));
				ENDDEBUG

				switch (MARKCLASS(NbElement)) {
					case YELLOW:
						if (REF_TYPE_CHANGES(NbElement)) continue;
						SET_NBELEM(SON(NbElement,0),l,sons[i*5+j].theSon);
						SET_NBELEM(sons[i*5+j].theSon,side,SON(NbElement,0));
						assert(sons[i*5+j+1].theSon == NULL);
						assert(j==0);
						break;
					case GREEN:
						/* green neighbor not refined yet */
						if (USED(NbElement) != 0) continue;
						/* TODO: delete this */
						/* if (NSONS(NbElement) == 0) continue; */
						if (j>0) continue;

						/* determine the sons of neighbor */
						if (GetSons(NbElement,NbSonList) != 0)
							return(1);

						/* determine son elements on current side */
						q = 0;
						for (m=0; m<MAX_SONS; m++) {
							if (NbSonList[m] == NULL) continue;
							p = 0;
							nbside = 0;
							if (TAG(NbSonList[m]) == PYRAMID) nbside =  pyrSide;
							for (n=0; n<CORNERS_OF_SIDE(NbSonList[m],nbside); n++) {
								for (o=0; o<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; o++) {
									if (theContext[o] == NULL) continue;
									if (CORNER(NbSonList[m],CORNER_OF_SIDE(NbSonList[m],nbside,n)) == theContext[o]) {
										IFDEBUG(gm,3)
										UserWriteF("                 NbSonList[m=%2d] matching CORNER theContext[%2d] ID=%d NBID=%d adr=%x nbadr=%x\n",m,o,ID(theContext[o]),ID(CORNER(NbSonList[m],CORNER_OF_SIDE(NbSonList[m],nbside,n))),theContext[o],CORNER(NbSonList[m],CORNER_OF_SIDE(NbSonList[m],nbside,n)));
										ENDDEBUG
										p++;
									}
								}
							}
							/* this son element is to connect */ 
							if (p == CORNERS_OF_SIDE(NbSonList[m],nbside)) {
								NbSideSons[q++] = NbSonList[m];
								IFDEBUG(gm,2)
								UserWriteF("             NbSonList[m=%d] matching SIDE p=%d for side=%d\n",m,p,nbside);
								ENDDEBUG
							}
						}
						assert(q>0 && q<=6);

						/* search corresponding elements */
						k = 0;
						for (m=0; m<q; m++) {
							found = 0;
							nbside = 0;
							if (TAG(NbSideSons[m]) == PYRAMID) nbside = pyrSide;
							for (n=i*5+j; n<i*5+q; n++) {
								r = 0;
								if (CORNERS_OF_SIDE(NbSideSons[m],nbside) != CORNERS_OF_SIDE(sons[n].theSon,nbside)) continue;
								for (o=0; o<CORNERS_OF_SIDE(NbSideSons[m],nbside); o++)
									for (p=0; p<CORNERS_OF_SIDE(sons[n].theSon,nbside); p++)
										if ((CORNER(NbSideSons[m],CORNER_OF_SIDE(NbSideSons[m],nbside,o))) == (CORNER(sons[n].theSon,CORNER_OF_SIDE(sons[n].theSon,nbside,p))))

											r++;
								if (r == CORNERS_OF_SIDE(sons[n].theSon,nbside)) {
									IFDEBUG(gm,2)
									UserWriteF("RefineGreenElement(): Matching sides for element ID=%d and neighbor ID=%d\n",ID(sons[n].theSon),ID(NbSideSons[m]));
									UserWriteF("                      m=%d n=%d i=%d j=%d q=%d r=%d\n",m,n,i,j,q,r);
									UserWrite("                      element node IDs:");
									for (o=0; o<CORNERS_OF_SIDE(NbSideSons[m],nbside); o++) UserWriteF(" %d",ID(CORNER(NbSideSons[m],CORNER_OF_SIDE(NbSideSons[m],nbside,o))));
									UserWrite("\n");
									UserWrite("                      neighbr node IDs:");
									for (p=0; p<CORNERS_OF_SIDE(sons[n].theSon,nbside); p++) UserWriteF(" %d",ID(CORNER(sons[n].theSon,CORNER_OF_SIDE(sons[n].theSon,nbside,p))));
									UserWrite("\n");
									ENDDEBUG
									found = 1;
									break;
								}
							}
							assert(found==1);
							assert(sons[n].theSon != NULL);
							assert(NbSideSons[m] != NULL);
							if (found) {
								SET_NBELEM(NbSideSons[m],nbside,sons[n].theSon);
								SET_NBELEM(sons[n].theSon,nbside,NbSideSons[m]);
								k++;
							}
						}
						assert(k == q);
						
						break;
					case RED:
						if (REF_TYPE_CHANGES(NbElement)) continue;
						NbRule = MARK2RULEADR(NbElement,REFINE(NbElement));
						found = 0;
						if (GetSons(NbElement,NbSonList)!=0)
							return (1);
						for (s2=0; s2<NbRule->nsons; s2++) {
							sdata2 = &(NbRule->sons[s2]);
							for (k=0; k<SIDES_OF_ELEM(NbSonList[s2]); k++) {
								if (sdata2->nb[k] != FATHER_SIDE_OFFSET+l) continue;

								IFDEBUG(gm,2)
									UserWriteF("elid=%3d: side:",ID(theSon));
									for (p=0; p<CORNERS_OF_SIDE(theSon,side); p++)
										UserWriteF(" %2d",ID(CORNER(theSon,CORNER_OF_SIDE(theSon,side,p))));
									UserWriteF(" OUTSIDE of father");
									UserWriteF("\nnbid=%3d: side:",ID(NbSonList[s2]));
									for (p=0; p<CORNERS_OF_SIDE(NbSonList[s2],k); p++)
										UserWriteF(" %2d",ID(CORNER(NbSonList[s2],CORNER_OF_SIDE(NbSonList[s2],k,p))));
									UserWriteF("\n\n");
								ENDDEBUG

								points=0;
								for (p=0; p<CORNERS_OF_SIDE(theSon,side); p++)
									for (q=0; q<CORNERS_OF_SIDE(NbSonList[s2],k); q++)
										if (CORNER(theSon,CORNER_OF_SIDE(theSon,side,p)) == CORNER(NbSonList[s2],CORNER_OF_SIDE(NbSonList[s2],k,q))) {
											/* look whether each corner of current neighbor side */
											/* has matching corner of current element side       */ 
											points |= ((1<<p) | (16<<q));
											break;
										}
								/* TODO: decrypt this expression and generalize for all side types 				*/ 
								/* (edges 2D, triangle and quad 3D) 							   				*/
								/* for tetra: 63 (0111 0111) and means all corners have matching corners        */
								/* for hexa: 127 (1111 1111) 													*/
								/* for quad and tri: 15 (0011 0011) 											*/
								switch (points) {
									#ifdef __TWODIM__
									case (LINEPOINTS):
									#endif
									#ifdef __THREEDIM__
									case (TRIPOINTS):
									case (QUADPOINTS): /* neighbor found */
									#endif
										/* no match for quadside with only three points */
										if (points==TRIPOINTS && CORNERS_OF_SIDE(theSon,side)==4) {
											PrintErrorMessage('E',"RefineGreenElement","quad side with 3 equal nodes");
											return(GM_FATAL);
										}

										IFDEBUG(gm,3)
											UserWriteF("Matching Sides:\n",ID(theSon));
											UserWriteF("elid=%3d: side:",ID(theSon));
											for (p=0; p<CORNERS_OF_SIDE(theSon,side); p++)
												UserWriteF(" %2d",ID(CORNER(theSon,CORNER_OF_SIDE(theSon,side,p))));
											UserWriteF(" OUTSIDE of father");
											UserWriteF("\nnbid=%3d: side:",ID(NbSonList[s2]));
											for (p=0; p<CORNERS_OF_SIDE(NbSonList[s2],k); p++)
												UserWriteF(" %2d",ID(CORNER(NbSonList[s2],CORNER_OF_SIDE(NbSonList[s2],k,p))));
											UserWriteF("\n\n");
										ENDDEBUG
										/* adjust pointers */
										SET_NBELEM(theSon,side,NbSonList[s2]);
										SET_NBELEM(NbSonList[s2],k,theSon);
								
										/* dispose doubled side vectors if */
                                        #ifdef __THREEDIM__
										if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR))										
										  if (DisposeDoubledSideVector(theGrid,theSon,side,NbSonList[s2],k))
											return (1);
                                        #endif
										found=1;
										break;
									default:
										break;
								}
								if (found) break;
							}
							if (found) break;
						}
						assert (found==1);

						break;
					default:
						break;
				}
			}
		}
	}
	return(GM_OK);
}


/****************************************************************************/
/*																			*/
/* Function:  RefineElement 												*/
/*																			*/
/* Purpose:   refine an element in the given context						*/
/*			  (i)	 corner and midnodes are already allocated				*/
/*			  (ii)	 edges between corner and midnodes are ok				*/
/*			  (iii)  create interior nodes and edges						*/
/*			  (iv)	 create sons and set references to sons 				*/
/*																			*/
/* Param:	  GRID *theGrid: grid level of sons of theElement				*/
/*			  ELEMENT *theElement: element to refine						*/
/*			  ELEMENTCONTEXT *theContext: current context of element		*/
/*																			*/
/* return:	  INT 0: ok 													*/
/*			  INT 1: fatal memory error 									*/
/*																			*/
/****************************************************************************/
		
static int RefineElement (GRID *theGrid, ELEMENT *theElement, NODE **theElementContext)
{
	INT i,j,l,s,s2,ni,pi,nsi,p,q,side,pyrSide,nbside;
	INT points,ni0,ni1,m,n,o,k,r;
	ELEMENT *theSon,*theNeighbor;
	EDGE *theEdge;
	ELEMENT *SonList[MAX_SONS],*SonList2[MAX_SONS];
	ELEMENT *NbSideSons[5];
	ELEMENTSIDE *oldSide,*newSide;
	VERTEX *myvertex;
	VSEGMENT *vs;
	NODE *mynode;
	COORD *v0,*v1,*v2, vd[3];
	COORD *x;
	INT boundaryelement, found;
	REFRULE *rule, *rule2;
	EDGEDATA *edata;
	SONDATA *sdata, *sdata2;
	
	/* is something to do ? */
	if (MARK(theElement)==0) return(0);

	if (DIM==3 && TAG(theElement)==HEXAHEDRON && MARKCLASS(theElement)==GREEN) {
		if (RefineGreenElement(theGrid,theElement,theElementContext) != GM_OK) return(GM_FATAL);
		return(GM_OK);
	}

	rule = MARK2RULEADR(theElement,MARK(theElement));

	/* create interior edges */
	for (s=0; s<MAX_NEW_EDGES_DIM; s++)
	{
		edata = &(rule->edges[s]);
		if (edata->type != INNER_EDGE) continue;
		if ((theEdge=CreateEdge(theGrid,theElementContext[edata->from], theElementContext[edata->to]))==NULL)
			return(1);
	}

	/* create elements */
	for (s=0; s<rule->nsons; s++)
	{
		boundaryelement = 0;
		/* TODO: how can boundary detection be generalized */
		if (OBJT(theElement) == BEOBJ)
			for (i=0; i<SIDES_OF_ELEMDESC(rule->sons[s].tag); i++)
					if ( (side = rule->sons[s].nb[i]) >= FATHER_SIDE_OFFSET )  /* exterior side */
						if (SIDE(theElement,side-FATHER_SIDE_OFFSET)!=NULL) 	/* at the boundary */
						{
							boundaryelement = 1;
							break;
						}

		if (boundaryelement)
				theSon = CreateBoundaryElement(theGrid,NULL,rule->sons[s].tag);
		else
				theSon = CreateInnerElement(theGrid,NULL,rule->sons[s].tag);
		if (theSon==NULL) return(GM_ERROR);

		/* fill in son data */
		SonList[s] = theSon;
		SETECLASS(theSon,MARKCLASS(theElement));
		SET_EFATHER(theSon,theElement);
		for (i=0; i<CORNERS_OF_ELEM(theSon); i++)
			SET_CORNER(theSon,i,theElementContext[rule->sons[s].corners[i]]);
		/* TODO: how can this be generalized */
		for (j=0; j<EDGES_OF_ELEM(theSon); j++)
		{
			theEdge = CreateEdge(theGrid,CORNER(theSon,CORNER_OF_EDGE(theSon,j,0)),CORNER(theSon,CORNER_OF_EDGE(theSon,j,1)));
			assert(theEdge!=NULL);
			if (NO_OF_ELEM(theEdge)<NO_OF_ELEM_MAX-1)
				INC_NO_OF_ELEM(theEdge);
			else
				return (GM_ERROR);
		}
	}
	SETNSONS(theElement,rule->nsons);
	#ifdef __TWODIM__
	for (i=0;i<NSONS(theElement); i++) SET_SON(theElement,i,SonList[i]);
	#endif
	#ifdef __THREEDIM__
	SET_SON(theElement,0,SonList[0]);
	#endif
	
	/* create element sides at the boundary */
	if (OBJT(theElement)==BEOBJ)
		for(s=0; s<rule->nsons; s++)
		{
				if (OBJT(SonList[s]) != BEOBJ) continue;
				for (j=0; j<SIDES_OF_ELEM(SonList[s]); j++)
				{
					SET_SIDE(SonList[s],j,NULL);
					if ((side = rule->sons[s].nb[j]) < FATHER_SIDE_OFFSET) continue;
					side -= FATHER_SIDE_OFFSET;
					if ((oldSide = SIDE(theElement,side)) == NULL) continue;

					newSide = CreateElementSide(theGrid);
					if (newSide==NULL) return(GM_FATAL);
					SET_SIDE(SonList[s],j,newSide);
					ES_PATCH(newSide) = ES_PATCH(oldSide);
	
					for (i=0; i<CORNERS_OF_SIDE(SonList[s],j); i++)
					{
							ni = CORNER_OF_SIDE(SonList[s],j,i);		  /* node index of son	 */
							pi = rule->sons[s].corners[ni];   /* point in theElement */

							if (pi<CORNERS_OF_ELEM(theElement)) /* a corner of theElement */
							{
								nsi = CORNER_OF_SIDE_INV(theElement,side,pi); /* position of corner in side numeration */
								assert(nsi != -1);
								PARAM(newSide,i,0) = PARAM(oldSide,nsi,0);
								PARAM(newSide,i,1) = PARAM(oldSide,nsi,1);
							}
							else /* a midpoint of an edge of theElement */
							{
								myvertex = MYVERTEX(theElementContext[pi]);
								assert (OBJT(myvertex) == BVOBJ);
								assert (VSEG(myvertex) != NULL);
					
								/* find common boundary segment */
								for( vs=VSEG(myvertex); vs!=NULL; vs = NEXTSEG(vs) )
								{
										if (VS_PATCH(vs) == ES_PATCH(oldSide)) break;
								}
								assert(vs!=NULL);
								PARAM(newSide,i,0) =  LAMBDA(vs,0);
								PARAM(newSide,i,1) =  LAMBDA(vs,1);
							}
					}
				}
		}
	
	/* connect elements */
	for (s=0; s<rule->nsons; s++)
	{
		sdata = &(rule->sons[s]);
		for (i=0; i<SIDES_OF_ELEM(SonList[s]); i++)
		{
			SET_NBELEM(SonList[s],i,NULL);

			/* an interior triangle face */
			if ( (side = sdata->nb[i]) < FATHER_SIDE_OFFSET )
			{
				SET_NBELEM(SonList[s],i,SonList[side]);

				IFDEBUG(gm,3)
					UserWriteF("elid=%3d: side:",ID(SonList[s]));
					for (p=0; p<CORNERS_OF_SIDE(SonList[s],i); p++)
						UserWriteF(" %2d",ID(CORNER(SonList[s],CORNER_OF_SIDE(SonList[s],i,p))));
					UserWriteF(" INSIDE of father");
					UserWriteF("\nnbid=%3d: side:",ID(SonList[side]));
					{
					int ss,qq,pp,f,pts;
					f=0;
					for (ss=0; ss<SIDES_OF_ELEM(SonList[side]); ss++)
					{
						pts = 0;
						for (pp=0; pp<CORNERS_OF_SIDE(SonList[s],i); pp++)
							for (qq=0; qq<CORNERS_OF_SIDE(SonList[side],ss); qq++)
								if (CORNER(SonList[s],CORNER_OF_SIDE(SonList[s],i,pp)) == CORNER(SonList[side],CORNER_OF_SIDE(SonList[side],ss,qq)))
								{
									pts |= ((1<<pp) | (16<<qq));
									break;
								}
						switch (pts)
						{
							#ifdef __TWODIM__
							case (LINEPOINTS):
							#endif
							#ifdef __THREEDIM__
							case (TRIPOINTS):
							case (QUADPOINTS): 
							#endif
								if (pts==TRIPOINTS && CORNERS_OF_SIDE(SonList[s],i)==4)
								{
									PrintErrorMessage('E',"RefineElement","quad side with 3 equal nodes");
									return(GM_FATAL);
								}
								f=1;
								break;
							default:
								break;
						}
						if (f) break;
					}
					assert (f=1);
					for (pp=0; pp<CORNERS_OF_SIDE(SonList[side],ss); pp++)
						UserWriteF(" %2d",ID(CORNER(SonList[side],CORNER_OF_SIDE(SonList[side],ss,pp))));
					}
					UserWriteF("\n\n");
				ENDDEBUG

				assert(SonList[side]!=NULL);
		
				/* dispose doubled side vectors if */
				#ifdef __THREEDIM__
				if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR))
				  {
					for (l=0; l<SIDES_OF_ELEM(theElement); l++)
					  if (rule->sons[side].nb[l]==s)
						break;
					if (DisposeDoubledSideVector(theGrid,SonList[s],i,SonList[side],l))
					return (1);
				  }
                #endif
				continue;
			}

			/* the boundary case */
			if ((OBJT(SonList[s]) == BEOBJ) && (SIDE(SonList[s],i) != NULL)) continue;

			/* check, if neighbor has been refined */
			side -= FATHER_SIDE_OFFSET;
			theNeighbor = NBELEM(theElement,side);

			if(theNeighbor==NULL)
			{
				IFDEBUG(gm,0)
				/* NULL only if theElement is an copy element */
				if (ECLASS(theElement)!=YELLOW)
				{
					PrintErrorMessage('E',"RefineElement","element has no neighbor, but is not on boundary and no yellow element!");
					return(GM_FATAL);
				}
				ENDDEBUG
				continue;
			}
	
			if (REF_TYPE_CHANGES(theNeighbor) || (REFINE(theNeighbor) == 0 ))
				continue;

			/* TODO: generalize this */
			#ifdef __TWODIM__
			for (l=0; l<EDGES_OF_ELEM(theNeighbor); l++) 
				if (NBELEM(theNeighbor,l) == theElement)
					break;
			#endif
			#ifdef __THREEDIM__
			for (l=0; l<SIDES_OF_ELEM(theNeighbor); l++) 
				if (NBELEM(theNeighbor,l) == theElement)
					break;
			#endif

			assert(l<SIDES_OF_ELEM(theNeighbor));			
	
			/* connect green hexahedral neighbor */
			if (DIM==3 && TAG(theNeighbor)==HEXAHEDRON && MARKCLASS(theNeighbor)==GREEN) {

				/* outer side for pyramid has 4 corners */
				for (n=0; n<element_descriptors[PYRAMID]->sides_of_elem; n++)
					if (element_descriptors[PYRAMID]->corners_of_side[n] == 4)
						break;
				pyrSide = n;

				/* green neighbor not refined yet */
				if (USED(theNeighbor) != 0) continue;

				/* determine the sons of neighbor */
				if (GetSons(theNeighbor,SonList2) != 0)
					return(1);

				/* determine son elements on current side */
				q = 0;
				for (m=0; m<MAX_SONS; m++) {
					if (SonList2[m] == NULL) continue;
					p = 0;
					nbside = 0;
					if (TAG(SonList2[m]) == PYRAMID) nbside = pyrSide;

					/* compare corners of neighbor elements with elementcontext */
					for (n=0; n<CORNERS_OF_SIDE(SonList2[m],nbside); n++) {
						for (o=0; o<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; o++) {
							if (theElementContext[o] == NULL) continue;
							if (CORNER(SonList2[m],CORNER_OF_SIDE(SonList2[m],nbside,n)) == theElementContext[o]) {
								IFDEBUG(gm,3)
								UserWriteF("                 SonList2[m=%2d] matching CORNER theElementContext[%2d] ID=%d NBID=%d adr=%x nbadr=%x\n",m,o,ID(theElementContext[o]),ID(CORNER(SonList2[m],CORNER_OF_SIDE(SonList2[m],nbside,n))),theElementContext[o],CORNER(SonList2[m],CORNER_OF_SIDE(SonList2[m],nbside,n)));
								ENDDEBUG
								p++;
							}
						}
					}
					/* this son element is to connect */ 
					if (p == CORNERS_OF_SIDE(SonList2[m],nbside)) {
						NbSideSons[q++] = SonList2[m];
						IFDEBUG(gm,2)
						UserWriteF("             SonList2[m=%d] matching SIDE p=%d for side=%d\n",m,p,nbside);
						ENDDEBUG
					}
				}
				assert(q>0 && q<=6);

				/* search corresponding elements */
				k = 0;
				for (m=0; m<q; m++) {
					found = 0;
					nbside = 0;
					if (TAG(NbSideSons[m]) == PYRAMID) nbside = pyrSide;

					/* compare corners of each side */
					r = 0;
					if (CORNERS_OF_SIDE(NbSideSons[m],nbside) != CORNERS_OF_SIDE(SonList[s],i)) continue;
					for (o=0; o<CORNERS_OF_SIDE(NbSideSons[m],nbside); o++)
						for (p=0; p<CORNERS_OF_SIDE(SonList[s],i); p++)
							if ((CORNER(NbSideSons[m],CORNER_OF_SIDE(NbSideSons[m],nbside,o))) == (CORNER(SonList[s],CORNER_OF_SIDE(SonList[s],i,p))))

								r++;
					
					/* matching sides */
					if (r == CORNERS_OF_SIDE(SonList[s],i)) {
						IFDEBUG(gm,2)
						UserWriteF("RefineElement(): Matching sides for element ID=%d and neighbor ID=%d\n",ID(SonList[s]),ID(NbSideSons[m]));
						UserWriteF("                      m=%d n=%d i=%d j=%d q=%d r=%d\n",m,n,i,j,q,r);
						UserWrite("                      neighbr node IDs:");
						for (o=0; o<CORNERS_OF_SIDE(NbSideSons[m],nbside); o++) UserWriteF(" %d",ID(CORNER(NbSideSons[m],CORNER_OF_SIDE(NbSideSons[m],nbside,o))));
						UserWrite("\n");
						UserWrite("                      element node IDs:");
						for (p=0; p<CORNERS_OF_SIDE(SonList[s],nbside); p++) UserWriteF(" %d",ID(CORNER(SonList[s],CORNER_OF_SIDE(SonList[s],i,p))));
						UserWrite("\n");
						ENDDEBUG
						found = 1;
					}
					assert(SonList[s] != NULL);
					assert(NbSideSons[m] != NULL);

					/* connect elements over these sides */
					if (found) {
						SET_NBELEM(NbSideSons[m],nbside,SonList[s]);
						SET_NBELEM(SonList[s],i,NbSideSons[m]);
						k++;
						break;
					}
				}
				assert(k == 1);
				
				continue;
			}

			rule2 = MARK2RULEADR(theNeighbor,REFINE(theNeighbor));
			found = 0;
			if (GetSons(theNeighbor,SonList2)!=0)
				return (1);
			for (s2=0; s2<rule2->nsons; s2++)
			{
				sdata2 = &(rule2->sons[s2]);
				for (j=0; j<SIDES_OF_ELEM(SonList2[s2]); j++)
				{
					if (sdata2->nb[j] != FATHER_SIDE_OFFSET+l) continue;

					IFDEBUG(gm,2)
						UserWriteF("elid=%3d: side:",ID(SonList[s]));
						for (p=0; p<CORNERS_OF_SIDE(SonList[s],i); p++)
							UserWriteF(" %2d",ID(CORNER(SonList[s],CORNER_OF_SIDE(SonList[s],i,p))));
						UserWriteF(" OUTSIDE of father");
						UserWriteF("\nnbid=%3d: side:",ID(SonList2[s2]));
						for (p=0; p<CORNERS_OF_SIDE(SonList2[s2],j); p++)
							UserWriteF(" %2d",ID(CORNER(SonList2[s2],CORNER_OF_SIDE(SonList2[s2],j,p))));
						UserWriteF("\n\n");
					ENDDEBUG

					points=0;
					for (p=0; p<CORNERS_OF_SIDE(SonList[s],i); p++)
						for (q=0; q<CORNERS_OF_SIDE(SonList2[s2],j); q++)
							if (CORNER(SonList[s],CORNER_OF_SIDE(SonList[s],i,p)) == CORNER(SonList2[s2],CORNER_OF_SIDE(SonList2[s2],j,q)))
							{
								/* look whether each corner of current neighbor side */
								/* has matching corner of current element side       */ 
								points |= ((1<<p) | (16<<q));
								break;
							}
					/* TODO: decrypt this expression and generalize for all side types 				*/ 
					/* (edges 2D, triangle and quad 3D) 							   				*/
					/* for tetra: 63 (0111 0111) and means all corners have matching corners        */
					/* for hexa: 127 (1111 1111) 													*/
					/* for quad and tri: 15 (0011 0011) 											*/
					switch (points)
					{
						#ifdef __TWODIM__
						case (LINEPOINTS):
						#endif
						#ifdef __THREEDIM__
						case (TRIPOINTS):
						case (QUADPOINTS): /* neighbor found */
						#endif
							/* no match for quadside with only three points */
							if (points==TRIPOINTS && CORNERS_OF_SIDE(SonList[s],i)==4)
							{
								PrintErrorMessage('E',"RefineElement","quad side with 3 equal nodes");
								return(GM_FATAL);
							}

							IFDEBUG(gm,3)
								UserWriteF("Matching Sides:\n",ID(SonList[s]));
								UserWriteF("elid=%3d: side:",ID(SonList[s]));
								for (p=0; p<CORNERS_OF_SIDE(SonList[s],i); p++)
									UserWriteF(" %2d",ID(CORNER(SonList[s],CORNER_OF_SIDE(SonList[s],i,p))));
								UserWriteF(" OUTSIDE of father");
								UserWriteF("\nnbid=%3d: side:",ID(SonList2[s2]));
								for (p=0; p<CORNERS_OF_SIDE(SonList2[s2],j); p++)
									UserWriteF(" %2d",ID(CORNER(SonList2[s2],CORNER_OF_SIDE(SonList2[s2],j,p))));
								UserWriteF("\n\n");
							ENDDEBUG
							/* adjust pointers */
							SET_NBELEM(SonList[s],i,SonList2[s2]);
							SET_NBELEM(SonList2[s2],j,SonList[s]);
					
							/* dispose doubled side vectors if */
                            #ifdef __THREEDIM__
							if (TYPE_DEF_IN_GRID(theGrid,SIDEVECTOR))
							  if (DisposeDoubledSideVector(theGrid,SonList[s],i,SonList2[s2],j))
								return (1);
                            #endif
							found=1;
							break;
						default:
							break;
					}
					if (found) break;
				}
				if (found) break;
			}
			assert (found==1);
		}
	}
	
	return(0);
}


/****************************************************************************/
/*																			*/
/* Function:  RefineGrid													*/
/*																			*/
/* Purpose:   refine one level of the grid									*/
/*																			*/
/* Param:	  GRID *theGrid: grid level to refine							*/
/*																			*/
/* return:	  INT 0: ok 													*/
/*			  INT 1: fatal memory error 									*/
/*																			*/
/****************************************************************************/

static int RefineGrid (GRID *theGrid)
{
	int i;
	ELEMENT *theElement;
	ELEMENTCONTEXT theContext;
	GRID *fineGrid;
	NODE *theNode;
	
	fineGrid = theGrid->finer;
	if (fineGrid==NULL) return(1);

	IFDEBUG(gm,1)
	UserWrite("RefineGrid():\n");
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
		UserWriteF("EID=%d TAG=%d ECLASS=%d RefineClass=%d MarkClass=%d Refine=%d Mark=%d Coarse=%d\n",ID(theElement),TAG(theElement),ECLASS(theElement),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement),COARSEN(theElement));
	ENDDEBUG
	
	/* refine elements */
	RESETGSTATUS(fineGrid,GRID_CHANGED);
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
	{
		if (REF_TYPE_CHANGES(theElement)||
			(DIM==3 && TAG(theElement)==HEXAHEDRON && MARKCLASS(theElement)==GREEN))
		{
			if (hFlag == 0 && MARKCLASS(theElement)!=RED) {
				SETMARK(theElement,NO_REF);
				SETMARKCLASS(theElement,0);
				continue; 
			}
			IFDEBUG(gm,1)
			UserWriteF("REFINING element ID=%d TAG=%d REFINECLASS=%d MARKCLASS=%d REFINE=%d MARK=%d\n",ID(theElement),TAG(theElement),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement));
			ENDDEBUG

			GetCurrentContext(theElement,theContext);

			IFDEBUG(gm,2)
			UserWrite("  CurrentContext is :\n");
			for(i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++)  
				UserWriteF("%3d",i);
			UserWrite("\n");
			for(i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++)  
			  if (theContext[i] != NULL)
				UserWriteF("%3d",ID(theContext[i]));
			UserWrite("\n");
			ENDDEBUG

			if (UnrefineElement(fineGrid,theElement,theContext))  return(1);
			if (UpdateContext(fineGrid,theElement,theContext)!=0) return(1);

			IFDEBUG(gm,2)
			UserWrite("  UpdateContext is :\n");
			for(i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++)  
				UserWriteF("%3d",i);
			UserWrite("\n");
			for(i=0; i<MAX_CORNERS_OF_ELEM+MAX_NEW_CORNERS_DIM; i++)  
			  if (theContext[i] != NULL)
				UserWriteF("%3d",ID(theContext[i]));
			UserWrite("\n");
			ENDDEBUG

			if (RefineElement(fineGrid,theElement,theContext)!=0) return(1);
			SETREFINE(theElement,MARK(theElement));
			SETREFINECLASS(theElement,MARKCLASS(theElement));
			SETGSTATUS(fineGrid,GRID_CHANGED);
			SETUSED(theElement,0);
		}
	}

	IFDEBUG(gm,1)
	UserWrite("RefineGrid():\n");
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
		UserWriteF("EID=%d TAG=%d ECLASS=%d RefineClass=%d MarkClass=%d Refine=%d Mark=%d Coarse=%d\n",ID(theElement),TAG(theElement),ECLASS(theElement),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement),COARSEN(theElement));
	ENDDEBUG
	
	/* reset coarse flags */
	for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement)) {
		SETCOARSEN(theElement,0);
	}

	/* set node class on next level */
	for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode))
		if (SONNODE(theNode)!=NULL)
		{
			SETCLASS(SONNODE(theNode),NCLASS(theNode));
			if (NCLASS(theNode)>=2) TOPNODE(MYVERTEX(theNode)) = SONNODE(theNode);
		}

	return(0);
}

			
   SYNOPSIS:
/*																			*/
/* Function:  RefineMultiGrid												*/
/*																			*/
/* Purpose:   refine whole multigrid structure								*/
/*																			*/
/* Param:	  MULTIGRID *theMG: multigrid to refine 						*/
/*			  INT 		flag:   flag for switching between different yellow */
/*								closures									*/
/*																			*/
/* return:	  INT 0: ok 													*/
/*			  INT 1: out of memory, but data structure as before			*/
/*			  INT 2: fatal memory error, data structure corrupted			*/
/*																			*/
	int level,toplevel,nrefined;
	int newlevel;
INT RefineMultiGrid (MULTIGRID *theMG, INT flag)
	GRID *theGrid, *FinerGrid;
	int j,k,r;

	DEBUG_TIME(0);

/*
	rFlag=flag & 0x03; 	 	/* copy local or all */
	/* set different flag */
	refine_seq = seq;

	No_Green_Update=0;
	
	/* prepare algebra (set internal flags correctly) */
	PrepareAlgebraModification(theMG);
		if (DropMarks(theMG))
			return(GM_ERROR);
	toplevel = TOPLEVEL(theMG);

	REFINE_MULTIGRID_LIST(1,theMG,"RefineMultiGrid()","","")
	
	j = theMG->topLevel;
	for (level=toplevel; level>0; level--)
	IFDEBUG(gm,1)
	UserWrite("RefineMultiGrid():\n");
	for (k=j; k>=0; k--)
	{
		theGrid = GRID_ON_LEVEL(theMG,k);
		for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
			UserWriteF("EID=%d TAG=%d ECLASS=%d RefineClass=%d MarkClass=%d Refine=%d Mark=%d Coarse=%d\n",ID(theElement),TAG(theElement),ECLASS(theElement),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement),COARSEN(theElement));
	}
	ENDDEBUG
		theGrid = GRID_ON_LEVEL(theMG,level);
		
	for (k=j; k>0; k--)
		{
		theGrid = GRID_ON_LEVEL(theMG,k);

		if (hFlag) {
			IFDEBUG(gm,1)
			UserWriteF("Begin CloseGrid(%d):\n",k);
			ENDDEBUG
			if ((r = CloseGrid(theMG->grids[k]))<0) {
				PrintErrorMessage('E',"RefineMultiGrid","error in CloseGrid");
				return(GM_ERROR);

			IFDEBUG(gm,1)
			UserWriteF("End CloseGrid(%d):\n",k);
			for (theElement=theMG->grids[k]->elements; theElement!=NULL; theElement=SUCCE(theElement))
				UserWriteF("EID=%d TAG=%d ECLASS=%d EFATHERID=%d RefineClass=%d MarkClass=%d Refine=%d Mark=%d Coarse=%d\n",ID(theElement),TAG(theElement),ECLASS(theElement),ID(EFATHER(theElement)),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement),COARSEN(theElement));
			ENDDEBUG

		RestrictMarks(theMG->grids[k-1]);
		IFDEBUG(gm,1)
		UserWriteF("End RestrictMarks(%d):\n",k-1);
		for (theElement=theMG->grids[k-1]->elements; theElement!=NULL; theElement=SUCCE(theElement))
		  if (k-1 == 0)
			UserWriteF("EID=%d TAG=%d ECLASS=%d RefineClass=%d MarkClass=%d Refine=%d Mark=%d Coarse=%d\n",ID(theElement),TAG(theElement),ECLASS(theElement),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement),COARSEN(theElement));
		  else
			UserWriteF("EID=%d TAG=%d ECLASS=%d EFATHERID=%d RefineClass=%d MarkClass=%d Refine=%d Mark=%d Coarse=%d\n",ID(theElement),TAG(theElement),ECLASS(theElement),ID(EFATHER(theElement)),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement),COARSEN(theElement));
		ENDDEBUG
	#endif


	for (k=0; k<=j; k++)
		if (level<toplevel) FinerGrid = GRID_ON_LEVEL(theMG,level+1); else FinerGrid = NULL;
		theGrid = GRID_ON_LEVEL(theMG,k);
		if (k<j) FinerGrid = GRID_ON_LEVEL(theMG,k+1); else FinerGrid = NULL;

		if (hFlag)
		{
		for (theNode=theGrid->firstNode; theNode!=NULL; theNode=SUCCN(theNode)) SETMODIFIED(theNode,0);
			for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
			{
				if ((ECLASS(theElement)==RED_CLASS) && MARKCLASS(theElement)==RED_CLASS) continue;
				SETMARK(theElement,NO_REFINEMENT);
			for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))

				if ((ECLASS(theElement)==RED) && MARKCLASS(theElement)==RED) continue;

			/* determine regular and irregular elements on next level */
			if ((nrefined = GridClosure(theGrid))<0)
				RETURN(GM_ERROR);
			IFDEBUG(gm,1)
			UserWriteF("Begin 2. CloseGrid(%d):\n",k);
			ENDDEBUG
			if ((r = CloseGrid(theGrid))<0) {
				PrintErrorMessage('E',"RefineMultiGrid","error in 2. CloseGrid");
				return(GM_ERROR);

			IFDEBUG(gm,1)
			UserWriteF("End 2. CloseGrid(%d):\n",k);
			for (theElement=theMG->grids[k]->elements; theElement!=NULL; theElement=SUCCE(theElement))
			  {
				if (k>0)
				  UserWriteF("EID=%d TAG=%d ECLASS=%d EFATHERID=%d RefineClass=%d MarkClass=%d Refine=%d Mark=%d Coarse=%d\n",ID(theElement),TAG(theElement),ECLASS(theElement),ID(EFATHER(theElement)),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement),COARSEN(theElement));
				else
				  UserWriteF("EID=%d TAG=%d ECLASS=%d RefineClass=%d MarkClass=%d Refine=%d Mark=%d Coarse=%d\n",ID(theElement),TAG(theElement),ECLASS(theElement),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement),COARSEN(theElement));
			  }
			ENDDEBUG
			ComputeCopies(theGrid);
			/* by the neighborhood of elements were MARK != REFINE. 					 */
			{
				for (theElement=FIRSTELEMENT(FinerGrid); theElement!=NULL; theElement=SUCCE(theElement))
				{
			if (k<j)
					if (REFINE(EFATHER(theElement))!=MARK(EFATHER(theElement))) 
						if (DisposeConnectionsInNeighborhood(FinerGrid,theElement)!=GM_OK)
							RETURN(GM_FATAL);
					assert(EFATHER(theElement) != NULL);
			}
		}
							return(GM_FATAL);
		/* TODO: bug fix to force new level creation */
		if (!hFlag)
		{
			/* set this variable>0 */
#endif
		if ( (r>0) && (k==j) )
		
			newlevel = 1;
			if (CreateNewLevel(theMG)==NULL)
				return(GM_FATAL);
			FinerGrid = GRID_ON_LEVEL(theMG,j+1);


#ifdef ModelP
		if ( k<j || newlevel )
			if (RefineGrid(theGrid)!=GM_OK) 
				return(GM_FATAL);
			
		if ((k<j)||(newlevel))
			
			/* and compute the vector classes on the new (or changed) level */
			ClearVectorClasses(FinerGrid);
			if (GridCreateConnection(FinerGrid)) return (GM_FATAL);
				if (ECLASS(theElement)>=GREEN_CLASS || (rFlag==GM_COPY_ALL)) 
				  SeedVectorClasses(FinerGrid,theElement);
			PropagateVectorClasses(FinerGrid);
				if (ECLASS(theElement)>=IRREGULAR_CLASS) 
		/* TODO: delete special debug */ PRINTELEMID(-1)
	DEBUG_TIME(0);

			
	#endif

	if (CreateAlgebra(theMG) != GM_OK)
	if (theMG->topLevel > 0) DisposeTopLevel(theMG);
	theMG->currentLevel = theMG->topLevel;

	/* set grid status of grid 0 */
	RESETGSTATUS(theMG->grids[0],GRID_CHANGED);

	IFDEBUG(gm,1)
	UserWrite("END RefineMultiGrid():\n");
	j = theMG->topLevel;
	for (k=j; k>=0; k--)
	{
		theGrid = GRID_ON_LEVEL(theMG,k);
		for (theElement=theGrid->elements; theElement!=NULL; theElement=SUCCE(theElement))
			UserWriteF("EID=%d TAG=%d ECLASS=%d RefineClass=%d MarkClass=%d Refine=%d Mark=%d Coarse=%d\n",ID(theElement),TAG(theElement),ECLASS(theElement),REFINECLASS(theElement),MARKCLASS(theElement),REFINE(theElement),MARK(theElement),COARSEN(theElement));
	}
	ENDDEBUG
}
