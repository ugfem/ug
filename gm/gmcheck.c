// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  gmcheck.c     												*/
/*																			*/
/* Purpose:   checks of the dtat structure   								*/
/*																			*/
/* Author:	  Stefan Lang                         							*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   Juli 1, 1997 moved from ugm.c     							*/
/*																			*/
/* Remarks:       															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/*		defines to exclude functions										*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files 									*/
/*																			*/
/****************************************************************************/

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <errno.h>

#include "compiler.h"
#include "heaps.h"
#include "ugenv.h"
#include "debug.h"
#include "general.h"
#include "fifo.h"

#include "devices.h"

#include "evm.h"
#include "gm.h"
#include "rm.h"
#include "misc.h"
#include "dlmgr.h"
#include "algebra.h"
#include "ugm.h"
#include "elements.h"
#include "shapes.h"
#include "refine.h"
#include "domain.h"

#ifdef ModelP
#include "parallel.h"
#endif

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define RESOLUTION       20     /* resolution for creating boundary midnode */
#define SMALL1 0.001

#define ORDERRES		1e-3	/* resolution for OrderNodesInGrid			*/
#define LINKTABLESIZE	32		/* max number of inks per node for ordering	*/

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

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

static char buffer[4*256];			/* general purpose text buffer			*/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

REP_ERR_FILE;

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

static INT CheckVertex (ELEMENT *theElement, NODE* theNode, VERTEX *theVertex)
{
	ELEMENT *theFather = VFATHER(theVertex);
	EDGE	*theEdge;
	INT		i,nerrors,n;
	DOUBLE  *global,*local,diff;
	DOUBLE_VECTOR global1;
	DOUBLE *x[MAX_CORNERS_OF_ELEM];

	nerrors = 0;
	if (theFather==NULL && MASTER(theNode) && LEVEL(theVertex)>0)
	{
        #ifdef ModelP
		if (!CORNERTYPE(theNode))
		{
		    nerrors = 0;
			IFDEBUG(gm,1)
		    nerrors = 1;
			ENDDEBUG
		}
        #endif
		if (nerrors == 0) return(nerrors);
		UserWriteF(PFMT "elem=" EID_FMTX " node=" ID_FMTX " vertex=" VID_FMTX
			" VFATHER=NULL vertex needs VFATHER\n",me,EID_PRTX(theElement),ID_PRTX(theNode),
			VID_PRTX(theVertex));
		return(nerrors++);
	}

	if (theFather!=NULL && HEAPCHECK(theFather))
	{
		UserWriteF(PFMT "elem=" EID_FMTX " node=" ID_FMTX " vertex=" VID_FMTX
			" VFATHER=%x is pointer to ZOMBIE\n",me,EID_PRTX(theElement),ID_PRTX(theNode),
			VID_PRTX(theVertex),theFather);
		return(nerrors++);
	}

	
	if (theFather!=NULL && MASTER(theNode) && EPRIO(theFather)==PrioHGhost)
	{
        #ifdef ModelP
		if (!CORNERTYPE(theNode))
		{
		    nerrors = 0;
			IFDEBUG(gm,1)
		    nerrors = 1;
			ENDDEBUG
		}
        #endif
		if (nerrors == 0) return(nerrors);
		UserWriteF(PFMT "elem=" EID_FMTX " node=" ID_FMTX " vertex=" VID_FMTX
			" VFATHER=" EID_FMTX " vertex needs VFATHER with prio master or vghost\n",
			me,EID_PRTX(theElement),ID_PRTX(theNode),VID_PRTX(theVertex),EID_PRTX(theFather));
		return(nerrors++);
	}

	if (theFather != NULL) {
	    CORNER_COORDINATES(theFather,n,x);
		global = CVECT(theVertex);
		local = LCVECT(theVertex);
		LOCAL_TO_GLOBAL(n,x,local,global1);
		V_DIM_EUKLIDNORM_OF_DIFF(global1,global,diff);
		if (diff > MAX_PAR_DIST) {
			nerrors++;
			#ifdef ModelP
			if (CORNERTYPE(theNode))
			{
				nerrors = 0;
				IFDEBUG(gm,1)
				nerrors = 1;
				ENDDEBUG
			}
			#endif
			if (nerrors >= 1)
			{
				UserWriteF(PFMT "elem=" EID_FMTX " node=" ID_FMTX "/%d vertex=" VID_FMTX
					" WARNING VFATHER=%x WARNING diff %f local and global coordinates don't match\n",
					me,EID_PRTX(theElement),ID_PRTX(theNode),
					NTYPE(theNode),VID_PRTX(theVertex),theFather,diff);
			}
		}
	}

	switch (NTYPE(theNode))
	{
		case (CORNER_NODE):
			if (LEVEL(theVertex)==0 && theFather != NULL)
			{
				UserWriteF(PFMT "EID=" EID_FMTX " NID=" ID_FMTX
					" VID=" VID_FMTX " CORNER_NODE has VFATHER\n",
					me,EID_PRTX(theElement),ID_PRTX(theNode),VID_PRTX(theVertex));
			}

			#ifdef ModelP
			IFDEBUG(gm,0)
			/* break for ghost nodes if debugging off */
			if (GHOST(theNode)) break;
			ENDDEBUG
			#endif

			if (LEVEL(theVertex)>0 && theFather == NULL)
			{
				UserWriteF(PFMT "EID=" EID_FMTX " NID=" ID_FMTX
					" VID=" VID_FMTX " CORNER_NODE has no VFATHER\n",
					me,EID_PRTX(theElement),ID_PRTX(theNode),VID_PRTX(theVertex));
			}
			break;

		case (MID_NODE):
			/* check ONEDGE and VFATHER */
			if (theFather == NULL)
			{
				#ifdef ModelP
				IFDEBUG(gm,0)
				/* break for ghost nodes if debugging off */
				if (GHOST(theNode)) break;
				ENDDEBUG
				#endif

				UserWriteF(PFMT "EID=" EID_FMTX " NID=" ID_FMTX 
					" VID=" VID_FMTX " MID_NODE VFATHER=NULL\n",
					me,EID_PRTX(theElement),ID_PRTX(theNode),VID_PRTX(theVertex));
				nerrors++;	
				break;
			}
			i = ONEDGE(theVertex);
			theEdge = GetEdge(CORNER(theFather,CORNER_OF_EDGE(theFather,i,0)),
							  CORNER(theFather,CORNER_OF_EDGE(theFather,i,1)));

			if (theEdge==NULL || theNode!=MIDNODE(theEdge))
			{
				nerrors++;	
                #ifdef ModelP
				if (EGHOST(theElement)) {
				    nerrors = 0;
					IFDEBUG(gm,1)
					    nerrors = 1;
					ENDDEBUG
				}
                #endif
				if (nerrors == 0) break;
				UserWriteF(PFMT "EID=" EID_FMTX " NID=" ID_FMTX " VID=" VID_FMTX 
					" ONEDGE and VFATHER incompatible edgeptr=%08x\n",
					me,EID_PRTX(theElement),ID_PRTX(theNode),
					VID_PRTX(theVertex),theEdge);
			}
			break;

		#ifdef __THREEDIM__
		case (SIDE_NODE):
			if (theFather == NULL)
			{
				nerrors++;	
                #ifdef ModelP
				if (EPRIO(theElement)==PrioHGhost) {
				    nerrors = 0;
					IFDEBUG(gm,1)
					    nerrors = 1;
					ENDDEBUG
				}
                #endif
				if (nerrors == 0) break;
				UserWriteF(PFMT "EID=" EID_FMTX " NID=" ID_FMTX 
					" VID=" VID_FMTX " SIDE_NODE VFATHER=NULL\n",
					me,EID_PRTX(theElement),ID_PRTX(theNode),VID_PRTX(theVertex));
				break;
			}
			else {
				if (GetSideNode(theFather,ONSIDE(theVertex)) != theNode) {
					nerrors = 1;
					UserWriteF(PFMT "EID=" EID_FMTX " NID=" ID_FMTX 
							   " VID=" VID_FMTX " inconsistent ONSIDE entry\n",
							   me,EID_PRTX(theElement),ID_PRTX(theNode),
							   VID_PRTX(theVertex));
				}
			}
			break;
        #endif

		case (CENTER_NODE):
			if (theFather == NULL)
			{
				nerrors++;	
                #ifdef ModelP
				if (EPRIO(theElement)==PrioHGhost) {
				    nerrors = 0;
					IFDEBUG(gm,1)
					    nerrors = 1;
					ENDDEBUG
				}
                #endif
				if (nerrors == 0) break;
				UserWriteF(PFMT "EID=" EID_FMTX " NID=" ID_FMTX 
					" VID=" VID_FMTX " CENTER_NODE VFATHER=NULL\n",
					me,EID_PRTX(theElement),ID_PRTX(theNode),VID_PRTX(theVertex));
				break;
			}
			break;
	}

	return(nerrors);
}

static INT CheckNode (ELEMENT *theElement, NODE* theNode, INT i) 
{
	VERTEX	*theVertex	= MYVERTEX(theNode);
	NODE	*FatherNode;
	EDGE	*FatherEdge;
	INT		nerrors		= 0;		

	SETUSED(theNode,1);

	if (OBJT(theNode) != NDOBJ)
	{
		UserWriteF(PFMT " node=" ID_FMTX " has wrong OBJ=%d\n",
			me,ID_PRTX(theNode),OBJT(theNode));
		return(nerrors++);
	}

	switch (NTYPE(theNode))
	{
	    case (LEVEL_0_NODE):
		    if (LEVEL(theNode) > 0) {
			    UserWriteF(PFMT " node=" ID_FMTX " has NTYPE=LEVEL_0_NODE"
						   " but is on level=%d\n",
						   me,ID_PRTX(theNode),LEVEL(theNode));
				return(nerrors++);
			}
		    break;
		case (CORNER_NODE):
			{
				FatherNode = (NODE *)NFATHER(theNode);
				if (FatherNode == NULL)
				{
					#ifdef ModelP
					if (MASTER(theNode))
					{
					#endif
						UserWriteF(PFMT " ERROR cornernode=" ID_FMTX " has no father level=%d\n",
							me,ID_PRTX(theNode),LEVEL(theNode));
						nerrors++;
					#ifdef ModelP
					}
					else
					{
						INT print = 0;
						IFDEBUG(gm,1)
						print = 1;
						ENDDEBUG
						if (print)
							UserWriteF(PFMT " WARN cornernode=" ID_FMTX " has no father level=%d\n",
								me,ID_PRTX(theNode),LEVEL(theNode));
					}
					#endif
				}
				if (FatherNode != NULL)
				{
					if (HEAPCHECK(FatherNode))
					{
						UserWriteF(PFMT "elem=" EID_FMTX " cornernode=%d NID=" ID_FMTX 
							" has father pointer to ZOMBIE\n",me,EID_PRTX(theElement),i,ID_PRTX(theNode));
						nerrors++;
						break;
					}

					if (OBJT(FatherNode) != NDOBJ)
					{
						UserWriteF(PFMT " cornernode=" ID_FMTX 
							" has father of wrong type=%d\n", 
							me,ID_PRTX(theNode),OBJT(FatherNode));
						nerrors++;
					}
					else
					{
						if (SONNODE(FatherNode) != theNode)
						{
							UserWriteF(PFMT " cornernode=" ID_FMTX 
								" has father with wrong backptr=%x\n", 
								me,ID_PRTX(theNode),SONNODE(FatherNode));
							nerrors++;
						}
					}
				}
			}
			break;

		case (MID_NODE):
			if (LEVEL(theNode)>0)
			{
				FatherEdge = (EDGE *)NFATHER(theNode);
				if (FatherEdge == NULL)
				{
					#ifdef ModelP
					if (MASTER(theNode))
					{
					#endif
						UserWriteF(PFMT " ERROR midnode=" ID_FMTX " has no father level=%d\n",
							me,ID_PRTX(theNode),LEVEL(theNode));
						nerrors++;
					#ifdef ModelP
					}
					else
					{
                        IFDEBUG(gm,1)
						UserWriteF(PFMT " WARN midnode=" ID_FMTX " has no father level=%d\n",
							me,ID_PRTX(theNode),LEVEL(theNode));
						ENDDEBUG 
					}
					#endif
				}
				if (FatherEdge != NULL)
				{
					if (HEAPCHECK(FatherEdge))
					{
						UserWriteF(PFMT "elem=" EID_FMTX " edge=%d/%x midnode NID=" ID_FMTX 
							" fatherpointer to edge=%d/%x is ZOMBIE\n",me,EID_PRTX(theElement),
							ID_PRTX(theNode),i,FatherEdge);
						nerrors++;
						break;
					}

						if (OBJT(FatherEdge) != EDOBJ)
					{
						UserWriteF(PFMT " midnode=" ID_FMTX 
							" has father of wrong type=%d obj=\n", 
							me,ID_PRTX(theNode),OBJT(FatherEdge));
						nerrors++;
					}
					else
					{

						if (MIDNODE(FatherEdge) != theNode)
						{
							UserWriteF(PFMT " midnode=" ID_FMTX 
								" has father with wrong backptr=%x\n", 
								me,ID_PRTX(theNode),MIDNODE(FatherEdge));
							nerrors++;
						}
					}
				}
			}
			else
			{
				UserWriteF(PFMT " node=" ID_FMTX " is midnode BUT on level=%d\n",
					me,ID_PRTX(theNode),LEVEL(theNode));
				nerrors++;
			}
			break;

		case (SIDE_NODE):
			break;

		case (CENTER_NODE):
			break;

		default:
			UserWriteF(PFMT " node=" ID_FMTX " has unrecognized NTYPE=%d\n",
				me,ID_PRTX(theNode),NTYPE(theNode));
			break;
	}

	if (theVertex != NULL)
	{
		CheckVertex(theElement,theNode,theVertex);
	}
	else
	{
		UserWriteF(PFMT "elem=" EID_FMTX " node[%d]=" ID_FMTX " vertex=NULL\n",
			me,EID_PRTX(theElement),i,ID_PRTX(theNode));
		nerrors++;	
	}

	return(nerrors);
}

static INT CheckEdge (ELEMENT *theElement, EDGE* theEdge, INT i)
{
	INT		nerrors = 0;
	NODE	*theNode;
	VERTEX	*theVertex;

	SETUSED(theEdge,1);

	theNode = MIDNODE(theEdge);
	if (theNode == NULL)
	{
#ifdef TET_RULESET
		if (((REFINE(theElement) == RED) && (TAG(theElement) != TETRAHEDRON))
			|| ((TAG(theElement) == TETRAHEDRON) && (NSONS(theElement) == 8)))
#else
		if (REFINE(theElement) == RED)
#endif
		{

        #ifdef ModelP
		IFDEBUG(gm,1)
	    #endif
				UserWriteF(PFMT "elem=" EID_FMTX " edge%d=" EDID_FMTX " midnode NID=NULL" 
				" BUT REFINE(elem)=RED\n",me,EID_PRTX(theElement),i,EDID_PRTX(theEdge));
		#ifdef ModelP
		ENDDEBUG
	    #endif
			return(nerrors++);
		}
		else
			return(nerrors);
	}

	if (HEAPCHECK(theNode))
	{
		UserWriteF(PFMT "elem=" EID_FMTX " edge=%d/%x midnode NID=" ID_FMTX 
			" is pointer to ZOMBIE\n",me,EID_PRTX(theElement),i,theEdge,ID_PRTX(theNode));
		return(nerrors++);
	}

	theVertex = MYVERTEX(theNode);
	if (theVertex == NULL)
	{
		UserWriteF(PFMT "elem=" EID_FMTX " edge=%d/%x midnode NID=" ID_FMTX " vertex=NULL\n",
			me,EID_PRTX(theElement),i,theEdge,ID_PRTX(theNode));
		return(nerrors++);	
	}

	if (VFATHER(theVertex) != theElement)
		return(nerrors);

	if (i != ONEDGE(theVertex))
	{
		if (EGHOST(theElement))
		{
			IFDEBUG(gm,1)
			UserWriteF(PFMT "EID=" EID_FMTX " VID=" VID_FMTX 
				" WARNING edgenumber of vertex wrong\n",
				me,EID_PRTX(theElement),VID_PRTX(theVertex));
			ENDDEBUG
		}
		else
		{
			UserWriteF(PFMT "EID=" EID_FMTX " VID=" VID_FMTX 
				" ERROR edgenumber of vertex wrong\n",
				me,EID_PRTX(theElement),VID_PRTX(theVertex));
			nerrors++;
		}
		return(nerrors);	
	}

	return(nerrors);
}

static INT CheckElement (GRID *theGrid, ELEMENT *theElement, INT *SideError, INT *EdgeError,
						 INT *NodeError, INT *ESonError, INT *NSonError)
{
	INT		i,j,k,l,n,nsons,bserror,nerrors;
	NODE	*theNode,*theNode1;
	EDGE	*theEdge;
	ELEMENT *NbElement,*theFather;
	ELEMENT *SonList[MAX_SONS];
	VERTEX	*theVertex,*Vertices[MAX_CORNERS_OF_ELEM];
PAR(
	DOUBLE  *x[MAX_CORNERS_OF_ELEM];
	DOUBLE_VECTOR center;
)ENDPAR
	
	*SideError = 0;
	*NodeError = 0;
	*EdgeError = 0; 
	*ESonError = 0;
	*NSonError = 0;
	nerrors    = 0;
	
	bserror = 0;

	/* check level */
	if (GLEVEL(theGrid) != LEVEL(theElement))
		UserWriteF(PFMT "elem=" EID_FMTX " ERROR level=%2d but gridlevel=%2d\n",
				me,EID_PRTX(theElement),LEVEL(theElement),LEVEL(theGrid));

	/* check side information */
	for (i=0; i<SIDES_OF_ELEM(theElement); i++)
	{
		NbElement = NBELEM(theElement,i);
		if (NbElement != NULL)
		{
			/* lets see if NbElement has the neighbor theElement */
			for (j=0; j<SIDES_OF_ELEM(NbElement); j++)
				if (NBELEM(NbElement,j) == theElement)
					break;
			if (j == SIDES_OF_ELEM(NbElement))
				*SideError |= (1<<i);
			else
			{
				/* if this is a boundary side it has to be an inner boundary 
				   and the neighbour side is also a boundary side */
				/* TODO: check boundary side for NbElement==NULL */
				if (OBJT(theElement) == BEOBJ)
					if (SIDE_ON_BND(theElement,i))
					{
						INT id,nbid,id_nb,nbid_nb,part;
						
						if (BNDS_BndSDesc(ELEM_BNDS(theElement,i),&id,&nbid,&part))
							bserror |= (1<<i);
						else
						{
							if ((id==0) || (nbid==0))
								/* no interior boundary */
								bserror |= (1<<i);
							if (id==nbid)
								/* should be avoided */
								bserror |= (1<<i);
							
							/* check neighbour */
							if (!SIDE_ON_BND(NbElement,j))
								bserror |= (1<<i);
							else
							{
								if (BNDS_BndSDesc(ELEM_BNDS(NbElement,j),&id_nb,&nbid_nb,&part))
									bserror |= (1<<i);
								else
								{
									if (id!=nbid_nb)
										bserror |= (1<<i);
									if (nbid!=id_nb)
										bserror |= (1<<i);
								}
							}
						}
					}
			}
			if (ECLASS(theElement)!=YELLOW_CLASS) 
			{
				n = CORNERS_OF_SIDE(theElement,i);
				for (k=0; k<n; k++)
					if (CORNER(theElement,CORNER_OF_SIDE(theElement,i,k))
						== CORNER(NbElement,CORNER_OF_SIDE(NbElement,j,0)))
						break;
				if (k == n)
					*SideError |= (1<<i);
				if (TAG(theElement)!=TETRAHEDRON 
				#ifdef Debug
				|| Debuggm>=1
				#endif
				)
				for (l=1; l<n; l++)
					if (CORNER(theElement,
						CORNER_OF_SIDE(theElement,i,(n+k-l)%n))
						!= CORNER(NbElement,CORNER_OF_SIDE(NbElement,j,l)))
					{
						*SideError |= (1<<i);
					}
			}
		}
		else
		{
			if (ECLASS(theElement)!=YELLOW_CLASS)
				if (OBJT(theElement) == IEOBJ)
				#ifdef ModelP
				if (EMASTER(theElement))
				#endif
					*SideError |= (1<<(i+MAX_SIDES_OF_ELEM));

			if (OBJT(theElement) == BEOBJ)
			{
				if (SIDE_ON_BND(theElement,i))
				{
					if (INNER_BOUNDARY(theElement,i)) {
						*SideError |= (1<<(i+2*MAX_SIDES_OF_ELEM));
						UserWriteF(PFMT "no nb Element for inner boundary, el =  " EID_FMTX "\n",
								   me,EID_PRTX(theElement));
					}
					for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
					{
						theVertex = MYVERTEX(CORNER(theElement,(k=CORNER_OF_SIDE(theElement,i,j))));
						if (OBJT(theVertex) == IVOBJ)
							*NodeError |= (1<<(k+MAX_CORNERS_OF_ELEM));
					}
				}
				else if (ECLASS(theElement)!=YELLOW_CLASS)
					#ifdef ModelP
					if (EMASTER(theElement))
					#endif
						*SideError |= (1<<(i+2*MAX_SIDES_OF_ELEM));
			}
		}
	}

	/* check node information */
	for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
	{
		theNode = CORNER(theElement,i);

		if (theNode != NULL)
			nerrors += CheckNode(theElement,theNode,i);
		else
			UserWriteF(PFMT "elem=" EID_FMTX " corner=%d nodeptr=NULL\n",
				me,EID_PRTX(theElement),i);
	}

	/* check edge information */
	for (i=0; i<EDGES_OF_ELEM(theElement); i++)
	{
		theNode	 = CORNER(theElement,CORNER_OF_EDGE(theElement,i,0));	
		theNode1 = CORNER(theElement,CORNER_OF_EDGE(theElement,i,1));

		if (theNode == NULL || theNode1 == NULL)
		{
			UserWriteF(PFMT "elem=" EID_FMTX " edge=%d n0ptr=NULL or n1ptr=NULL\n",
				me,EID_PRTX(theElement),i,theNode,theNode1);
			continue;
		}

		theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
						  CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));

		if (theEdge != NULL)
			nerrors += CheckEdge(theElement,theEdge,i);
		else
			UserWriteF(PFMT "elem=" EID_FMTX " edge=%d n0=" ID_FMTX " n1=" 
				ID_FMTX " edgeptr=NULL\n",
				me,EID_PRTX(theElement),i,ID_PRTX(theNode),ID_PRTX(theNode1));
	}
	
	/* check orientation */
	for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
		Vertices[i] = MYVERTEX(CORNER(theElement,i));
if (0)
	if (!CheckOrientation(CORNERS_OF_ELEM(theElement),Vertices))
	{
			UserWriteF(PFMT "elem=" EID_FMTX " wrong orientation",me,EID_PRTX(theElement));
			nerrors++;
	}

	/* check father information */
	theFather = EFATHER(theElement);
	if (theFather != NULL)
	{
		/* check MIDNODE information of father */
		for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
		{
			theNode = CORNER(theElement,i);
			if (NTYPE(theNode) == MID_NODE)
			{
				for (j=0; j<EDGES_OF_ELEM(theFather); j++)
				{
					theEdge = GetEdge(CORNER(theFather,
									  CORNER_OF_EDGE(theFather,j,0)),
									  CORNER(theFather,
									  CORNER_OF_EDGE(theFather,j,1)));
					if (MIDNODE(theEdge) == theNode) break;
				}
				if (j >= EDGES_OF_ELEM(theFather)) 
				{
					#ifdef ModelP
					if (EMASTER(theFather))
					#endif
						UserWriteF(PFMT "ELEM(" EID_FMTX ") ERROR MIDNODE=NULL"
							" for mid node[%d]=" ID_FMTX "\n",
							me,EID_PRTX(theFather),i,ID_PRTX(theNode));
					#ifdef ModelP
					else
                        IFDEBUG(gm,1) 
						UserWriteF(PFMT "ELEM(" EID_FMTX ") WARNING MIDNODE=NULL"
							" for mid node[%d]" ID_FMTX "\n",
							me,EID_PRTX(theFather),i,ID_PRTX(theNode));
					    ENDDEBUG 
					#endif
				}
			}
		}

		/* check son information of father     */
		if (GetAllSons(theFather,SonList))
		{
			UserWrite("cannot get sons\n");
			return (1);
		}
		for (i=0; i<NSONS(theFather); i++)
		{
			if (SonList[i] == theElement) break;
		}
		if (i == NSONS(theFather))
		{
			UserWriteF(PFMT "ELEM(" EID_FMTX ") FATHER(" EID_FMTX 
				")element is not in SonList NSONS=%d\n",
				me,EID_PRTX(theElement),EID_PRTX(theFather),
				NSONS(theFather));
		}
	}
	#ifdef ModelP
	else
	{
		if (LEVEL(theElement) > 0)
		{
			if (EMASTER(theElement))
				UserWriteF(PFMT "ELEM(" EID_FMTX ") ERROR father=NULL\n",
					me,EID_PRTX(theElement));
			else
			{
				CORNER_COORDINATES(theElement,n,x);
				V_DIM_CLEAR(center);
				for (i=0; i<n; i++)
					V_DIM_ADD(center,x[i],center)
				V_DIM_SCALE(1.0/(DOUBLE)n,center)

				/* search in element list of grid level below */
				{
				GRID		*downGrid;
				MULTIGRID	*theMG;
				INT			downlevel,downlevel1;

				theMG		= MYMG(theGrid);
				downlevel	= GLEVEL(theGrid)-1;
				downGrid	= GRID_ON_LEVEL(theMG,downlevel);

				theFather = FindElementFromPosition(downGrid,center);
				if (theFather != NULL)
					UserWriteF(PFMT "ELEM(" EID_FMTX ") has no fatherpointer but father="
						EID_FMTX "\n",me,EID_PRTX(theElement),EID_PRTX(theFather));

			if (0)
				for (theFather = PFIRSTELEMENT(downGrid);
					 theFather != NULL; theFather = SUCCE(theFather))
				{
					k = PointInElement(center,theFather);
					switch (k)
					{
						case 0: 
							UserWriteF(PFMT "ELEM(" EID_FMTX ") PointInElement() returned"
								" error",me,EID_PRTX(theElement));
							break;
						case 1:
						case 2:
						case 3:
						case 4: 
							/* point (nearly) in father */
							UserWriteF(PFMT "ELEM(" EID_FMTX ") has no fatherpointer but father="
								EID_FMTX "\n",me,EID_PRTX(theElement),EID_PRTX(theFather));
							break;
						case 5:
							/* point not in father */
							break;
						default:
							UserWriteF(PFMT "ELEM(" EID_FMTX ") PointInElement() unexpected"
								" Returncode=%d",me,EID_PRTX(theElement),k);
							break;
					}
				}
				}
			}
		}
	}
	#endif
	
	/* check son information */
	if (NSONS(theElement)!=0)
	{
		nsons = NSONS(theElement);

		if (GetAllSons(theElement,SonList))
		{
			UserWrite("cannot get sons\n");
			return (1);
		}
		for (i=0; SonList[i]!=NULL || i<nsons; i++)
		{
		    IFDEBUG(gm,1)
			if (REFINE(theElement)==0)
			{
				UserWriteF(PFMT "ELEM(" EID_FMTX "): element is not refined "
					"but has NSONS=%d\n",me,EID_PRTX(theElement),nsons);
			}
			ENDDEBUG

			if (i >= nsons)
			{
				UserWriteF(PFMT "ELEM(" EID_FMTX "): element has nsons=%d but "
					" son[%d]=" EID_FMTX " exists\n", me,EID_PRTX(theElement),
					NSONS(theElement),i,EID_PRTX(SonList[i]));
			}

			if (SonList[i] == NULL)
			{
				UserWriteF(PFMT "ELEM(" EID_FMTX "): element has nsons=%d but "
					" son[%d]=NULL\n", me,EID_PRTX(theElement),nsons,i);
				*ESonError |= (1<<i);
				continue;
			}
			if (EFATHER(SonList[i])!=theElement)
			{
				UserWriteF(PFMT "i=%d theElement=" EID_FMTX 
					" SonList[i]=" EID_FMTX "\n",
					me,i,EID_PRTX(theElement),EID_PRTX(SonList[i]));
				*ESonError |= (1<<i);
			}
		}
	}
	
	if (bserror)
		nerrors++;
	if (nerrors > 0)
		UserWriteF("ELEM(" EID_FMTX "): element has %d errors\n",
			EID_PRTX(theElement),nerrors);

	if (*SideError || *EdgeError || *NodeError || *ESonError || *NSonError)
		return (1);
		
	return (0);
}		

static INT CheckGeometry (GRID *theGrid)
{
	NODE *theNode;
	ELEMENT *theElement;
	EDGE *theEdge;
	LINK *theLink;
	int i,j;
	INT SideError, EdgeError, NodeError, ESonError, NSonError, count;
	INT errors = 0;

	/* reset used flags */
	for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
	{
		SETUSED(theNode,0);
		for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
			SETUSED(MYEDGE(theLink),0);
	}

	/* check elements */
	for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
		 theElement=SUCCE(theElement))
	{
		if (CheckElement(theGrid,theElement,&SideError,&EdgeError,
				&NodeError,&ESonError,&NSonError)==0) continue;

		UserWriteF("ELEM=" EID_FMTX "\n",EID_PRTX(theElement));

		/* evaluate side information */
		if (SideError)
			for (i=0; i<SIDES_OF_ELEM(theElement); i++)
			{
				/* back pointer failure */
				if (SideError & 1<<i)
				{
					errors++;

					UserWriteF("   SIDE %d=(",i);
					for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
					{
						UserWriteF(ID_FMTX,ID_PRTX(CORNER(theElement,
							CORNER_OF_SIDE(theElement,i,j))));

						if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
					}
					UserWriteF(") has neighbour=" EID_FMTX " but a backPtr does not exist\n",
						EID_PRTX(NBELEM(theElement,i)));
				}

				/* neighbor pointer failure */
				if (SideError & 1<<(i+MAX_SIDES_OF_ELEM))
				{
					errors++;

					UserWrite("   SIDE(");
					for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
					{
						UserWriteF(ID_FMTX,ID_PRTX(CORNER(theElement,
							CORNER_OF_SIDE(theElement,i,j))));

						if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
					}
					UserWrite(") has no neighbour but element is IEOBJ\n");
				}

				/* boundary failure */
				if (SideError & 1<<(i+2*MAX_SIDES_OF_ELEM))
				{
					errors++;

					UserWrite("   SIDE(");
					for (j=0; j<CORNERS_OF_SIDE(theElement,i); j++)
					{
						UserWriteF(ID_FMTX,ID_PRTX(CORNER(theElement,
							CORNER_OF_SIDE(theElement,i,j))));

						if (j<CORNERS_OF_SIDE(theElement,i)-1) UserWrite(",");
					}
					UserWrite(") has no neighbour, element is BEOBJ "
						"but there is no SIDE\n");
				}
			}

		/* evaluate edge information */
		if (EdgeError)
			for (i=0; i<EDGES_OF_ELEM(theElement); i++)
			{
				if (!(EdgeError & 1<<i)) continue;

				errors++;
				UserWriteF("   EDGE(" ID_FMTX " , " ID_FMTX ") is missing\n",
					ID_PRTX(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0))),
					ID_PRTX(CORNER(theElement,CORNER_OF_EDGE(theElement,i,1))));
			}

		/* evaluate node information */
		if (NodeError)
			for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
			{
				if (NodeError & (1<<i))
				{
					errors++;
					UserWriteF("   CORNER=" ID_FMTX " is BVOBJ," 
						" ids from elementside "
						"and vertexsegment are not consistent\n",
						ID_PRTX(CORNER(theElement,i)));
				}
				if (NodeError & (1<<(i+MAX_CORNERS_OF_ELEM)))
				{
					errors++;
					UserWriteF("   CORNER " ID_FMTX " is IVOBJ, but lies on "
						"elementside\n",ID_PRTX(CORNER(theElement,i)));
				}
			}

		/* evaluate son information */
		if (ESonError)
		{
			for (i=0; i<NSONS(theElement); i++)
			{
				if ((ESonError & 1<<i))
				{
					errors++;
					UserWriteF("   ESON(%d) has wrong EFATHER "
						"pointer\n",i);
				}
			}
		}

		if (NSonError)
		{
			for (i=0; i<MAX_CORNERS_OF_ELEM; i++)
			{
				if (NSonError & (1<<i))
				{
					errors++;
					UserWriteF("   SONNODE(CORNER %d) != CORNER(ESON)\n",i);
				}
				if (NSonError & (1<<(i+MAX_CORNERS_OF_ELEM)))
				{
					errors++;
					UserWriteF("   CORNER %d != EFATHER(CORNER(ESON))\n",i);
				}
			}
			
			for (i=0; i<MAX_EDGES_OF_ELEM; i++)
			{
				
				if (NSonError & (1<<(i+MAX_CORNERS_OF_ELEM)))
				{
					errors++;
					UserWriteF("   MIDNODE(edge %d) != CORNER(ESON)\n",i);
				}
			}

			if (NSonError & (1<<(MAX_EDGES_OF_ELEM+2*MAX_CORNERS_OF_ELEM)))
			{
				errors++;
				UserWriteF("   NFATHER(CENTERNODE(ESON)) != NULL\n");
			}
		}
	}

	/* look for dead edges */
	for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
	{
		for (theLink=START(theNode); theLink!=NULL; theLink=NEXT(theLink))
		{
			theEdge = MYEDGE(theLink);
			if (!USED(theEdge))
			{
				errors++;
				UserWriteF("edge between " ID_FMTX " and " ID_FMTX 
					" has no element, NO_OF_ELEM=%d \n",
					ID_PRTX(theNode),ID_PRTX(NBNODE(theLink)),
					NO_OF_ELEM(theEdge));

				#ifdef Debug
				{
					NODE *nb;
					LINK *theLink1;

					nb = NBNODE(theLink);
					UserWriteF("linklist of nbnode %d:",ID(nb));

					for (theLink1=START(nb); theLink1!=NULL; 
						 theLink1=NEXT(theLink1))
						UserWriteF(" %d-%d",ID(NBNODE(theLink1)),
							ID(NBNODE(REVERSE(theLink1))));
					UserWrite("\n");
				}
				#endif
			}
		}
	}

	/* look for dead nodes */
	for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
	{
		if (!USED(theNode))
		{
			errors++;
			UserWriteF("node=" ID_FMTX " is dead\n",ID_PRTX(theNode));
		}
		else
			SETUSED(theNode,0);
	}

	/* check number of elem and their pointers */
	count = 0;
	for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; 
		 theElement=SUCCE(theElement))
	{
		if (SUCCE(theElement)!=NULL)
		{
			if (OBJT(SUCCE(theElement))!=IEOBJ && 
				OBJT(SUCCE(theElement))!=BEOBJ)
			{
				errors++;
				UserWriteF("pointer of ELEM(" EID_FMTX ") (number %ld) "
					"to next element is no pointer to an element\n",
					EID_PRTX(theElement),(long)count);
				break;
			}
			if (PREDE(SUCCE(theElement))!=NULL)
			{
				if (PREDE(SUCCE(theElement))!=theElement)
				{
					errors++;
					UserWriteF("pointer of ELEM(" EID_FMTX ") (number %ld) "
						"to previous element is not the previous element\n",
						EID_PRTX(SUCCE(theElement)),(long)(count+1));
				}
			}
			#ifndef ModelP
			else
			{
				errors++;
				UserWriteF("pointer of ELEM(" EID_FMTX ") (number %ld) "
					"to previous element is NULL\n",
					EID_PRTX(SUCCE(theElement)),(long)(count+1));
			}
			#endif
		}
		count++;
	}

	if (FIRSTELEMENT(theGrid) != NULL)
		if (PREDE(FIRSTELEMENT(theGrid)) != NULL)
		{
			errors++;
			UserWriteF("first element of the grid has a previous 'element'\n");
		}

	if (LASTELEMENT(theGrid) != NULL)
		if (SUCCE(LASTELEMENT(theGrid)) != NULL)
		{
			errors++;
			UserWriteF("last element of the grid has a following 'element'\n");
		}

	if (count != NT(theGrid))
	{
		errors++;
		UserWriteF("there are %ld elements but %ld expected\n",(long)(count),
			(long)NT(theGrid));
	}
	
	return(errors);
}

static INT CheckElementList (GRID *theGrid)
{
	ELEMENT *theElement;

	if (LEVEL(theGrid) <= 0) return(0);

	for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL; 
		 theElement=SUCCE(theElement))
	{
		ELEMENT *Father	= EFATHER(theElement);
PAR(	INT 	prio 	= EPRIO(theElement);               )ENDPAR

		if (Father == NULL)  
		{
			UserWriteF(PFMT "ERROR: element=" EID_FMTX " has no father\n",
				me,EID_PRTX(theElement));
			continue;
		}
		if (theElement == SON(Father,PRIO2INDEX(prio)))
		{
			if (PREDE(theElement) != NULL)
				if (EFATHER(PREDE(theElement))==Father
PAR(				&& EPRIO(theElement)==EPRIO(PREDE(theElement)) )ENDPAR )
				{
					UserWriteF(PFMT " ERROR element=" EID_FMTX " is not first"
						"son in list pred elem=" EID_FMTX " father=" EID_FMTX
						"\n",me,EID_PRTX(theElement),EID_PRTX(PREDE(theElement)),
						EID_PRTX(Father));
				}
		}
		else
		{
			if (PREDE(theElement)==NULL || 
				EFATHER(PREDE(theElement))!=Father)
			{
				UserWriteF(PFMT " ERROR element=" EID_FMTX " has no"
					"PREDE with same father=" EID_FMTX
					"\n",me,EID_PRTX(theElement),EID_PRTX(Father));
			}
		}
    }
	return (0);
}

INT CheckLists (GRID *theGrid)
{
	int objs = 0;

	/* perform gm dependent check */
	CheckElementList(theGrid);

	/* perform standard list check */
	GRID_CHECK_ELEMENT_LIST(theGrid);
	GRID_CHECK_NODE_LIST(theGrid);
	GRID_CHECK_VERTEX_LIST(theGrid);
	GRID_CHECK_VECTOR_LIST(theGrid);

	return(GM_OK);
}

/****************************************************************************/
/*D
   CheckGrid - Check consistency of data structure

   SYNOPSIS:
   INT CheckGrid (GRID *theGrid, INT checkgeom, INT checkalgebra, INT checklists, INT checkif)

   PARAMETERS:
.  theGrid - grid to check
.  checkgeom - check geomtry
.  checkalgebra - check algebra
.  checklists - checklists
.  checkif - check the processor interfaces

   DESCRIPTION:
   This function checks the consistency of data structure.

   RETURN VALUE:
   INT
.n   GM_OK if ok
.n   GM_ERROR if an error occured.
D*/
/****************************************************************************/

#ifndef ModelP
INT CheckGrid (GRID *theGrid, INT checkgeom, INT checkalgebra, INT checklists)
#else
INT CheckGrid (GRID *theGrid, INT checkgeom, INT checkalgebra, INT checklists,
			   INT checkif)
#endif
{
	INT error		= 0;
	INT errors		= 0; 
	INT totalerrors	= 0; 

	/* check geometrical data structures */
	if (checkgeom)
	{
		UserWrite(" geometry:");
		fflush(stdout);

		if ((errors = CheckGeometry(theGrid)) != GM_OK)
		{
			totalerrors += errors;
			error++;
			UserWriteF(" geometry BAD: %d errors",errors);
			fflush(stdout);
		}
		else
			UserWrite(" ok");
	}

	/* check algebraic data structures */
	if (checkalgebra)
	{
		UserWrite(", algebra:");
		fflush(stdout);

		if ((errors = CheckAlgebra(theGrid)) != GM_OK)
		{
			totalerrors += errors;
			error++;
			UserWriteF(" algebra BAD: %d errors",errors);
			fflush(stdout);
		}
		else
			UserWrite(" ok");
	}

	/* check lists and counters */
	if (checklists)
	{
		UserWrite(", lists:");
		fflush(stdout);

		if ((errors = CheckLists(theGrid)) != GM_OK)
		{
			totalerrors += errors;
			error++;
			UserWriteF(" lists BAD: %d errors",errors);
			fflush(stdout);
		}
		else
			UserWrite(" ok");
	}

	#ifdef ModelP
	/* check interfaces to other procs */
	if (checkif)
	{
		UserWrite(", interface:");
		fflush(stdout);

		if (errors = CheckInterfaces(theGrid) != GM_OK)
		{
			totalerrors += errors;
			error++;
			UserWriteF(" interfaces BAD: %d errors",errors);
			fflush(stdout);
		}
		else
			UserWrite(" ok");
	}
	#endif

	if (totalerrors) 
		UserWriteF(", grid BAD: %d check(s) with %d totalerror(s)",
			error,totalerrors);
	else 
		UserWrite(", grid ok");

	return(error);
}
