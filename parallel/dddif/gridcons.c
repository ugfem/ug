// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  gridcons.c													*/
/*																			*/
/* Purpose:   basic functions for managing consistency of distributed grids */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken										*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: birken@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   960906 kb  begin 												*/
/*																			*/
/* Remarks: 																*/
/*																			*/
/****************************************************************************/

#ifdef ModelP

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files 									*/
/*																			*/
/****************************************************************************/

#include <stdlib.h>

#include "debug.h"
#include "parallel.h"
#include "general.h"
#include "gm.h"
#include "ugm.h"
#include "evm.h"
#include "shapes.h"
#include "devices.h"

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/* macros for merge new priority with objects existing one */
/* valid only for all types of ghost priorities            */
#define PRIO_CALC(e) ((USED(e) && THEFLAG(e)) ? PrioVHGhost :                \
						(THEFLAG(e)) ? PrioVGhost : (USED(e)) ?              \
						PrioHGhost : (assert(0),0))

/* macros for setting object priorities with related objects */
/* macros for setting object priorities with related objects */
#define NODE_PRIORITY_SET(g,n,prio)                                          \
		{                                                                    \
			/* set priorities of node */                                     \
			SETPRIOX(n,prio);                                                \
                                                                             \
			if (VEC_DEF_IN_OBJ_OF_GRID(g,NODEVEC))                           \
			    if (NVECTOR(n) != NULL)                                      \
				    SETPRIOX(NVECTOR(n),prio);                               \
		}

#ifdef __TWODIM__
#define PRIO_SET_EDGE(e,prio)
#endif
#ifdef __THREEDIM__
#define PRIO_SET_EDGE(e,prio)  SETPRIOX(e,prio);
#endif

#define EDGE_PRIORITY_SET(g,e,prio)                                          \
		{                                                                    \
			/* set priorities of node for 3D */                              \
			PRIO_SET_EDGE(e,prio)                                            \
                                                                             \
			/* set priority of edge vector */                                \
			if (VEC_DEF_IN_OBJ_OF_GRID(g,EDGEVEC))                           \
			    if (EDVECTOR(e) != NULL)                                     \
				    SETPRIOX(EDVECTOR(e),prio);                              \
		}

#define CHECK_OBJECT_PRIO(o,prio,master,ghost,id,s,_nerr_)                   \
	if (USED(o)==1 && ! master (o))                                          \
	{                                                                        \
			UserWriteF("MASTER %s=" id ## _FMTX " has WRONG prio=%d\n",      \
				s, id ## _PRTX(o),prio(o));                                  \
			_nerr_++;                                                        \
	}                                                                        \
	if (USED( o )==0 && ! ghost ( o ))                                       \
	{                                                                        \
			UserWriteF("GHOST %s=" id ## _FMTX " has WRONG prio=%d\n",       \
				s, id ## _PRTX( o ),prio(o));                                \
			_nerr_++;                                                        \
	}


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

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/*
	for all PrioMaster-nodes with remote copies, set exactly one
	to PrioMaster, the other copies to PrioBorder in order to establish
	the BorderNodeIF. this is done for one grid.
*/


static int ComputeNodeBorderPrios (DDD_OBJ obj)
{
	NODE    *node  = (NODE *)obj;
	int     *plist = DDD_InfoProcList(PARHDR(node));
	int      i, min_proc = procs;

	/*
		minimum processor number will get Master-node,
		all others get Border-nodes
	*/
	for(i=0; plist[i]>=0; i+=2)
	{
		if (plist[i+1]==PrioMaster && plist[i]<min_proc)
			min_proc = plist[i];
	}

	if (min_proc==procs)
		return(0);

	if (me!=min_proc)
		SETPRIO(node, PrioBorder);
}

static int ComputeVectorBorderPrios (DDD_OBJ obj)
{
	VECTOR  *vector  = (VECTOR *)obj;
	int     *plist = DDD_InfoProcList(PARHDR(vector));
	int      i, min_proc = procs;

	/*
		minimum processor number will get Master-node,
		all others get Border-nodes
	*/
	for(i=0; plist[i]>=0; i+=2)
	{
		if (plist[i+1]==PrioMaster && plist[i]<min_proc)
			min_proc = plist[i];
	}

	if (min_proc==procs)
		return(0);

	if (me!=min_proc)
		SETPRIO(vector, PrioBorder);
}

#ifdef __THREEDIM__
static int ComputeEdgeBorderPrios (DDD_OBJ obj)
{
	EDGE	*edge  =	(EDGE *)obj;
	int		*plist =	DDD_InfoProcList(PARHDR(edge));
	int		i, min_proc	= procs;

	/*
		minimum processor number will get Master-node,
		all others get Border-nodes
	*/
	for(i=0; plist[i]>=0; i+=2)
	{
		if (plist[i+1]==PrioMaster && plist[i]<min_proc)
			min_proc = plist[i];
	}

	if (min_proc==procs)
		return(0);

	if (me!=min_proc)
		SETPRIO(edge, PrioBorder);
}
#endif

void SetGhostObjectPriorities (GRID *theGrid)
{
	ELEMENT *theElement,*theNeighbor,*SonList[MAX_SONS];
	NODE	*theNode;
	EDGE	*theEdge;
	VECTOR	*theVector;
	INT 	i,prio,*proclist,hghost,vghost;

	/* reset USED flag for objects of ghostelements */ 
	for (theElement=PFIRSTELEMENT(theGrid);
		 theElement!=NULL;
		 theElement=SUCCE(theElement))
	{
		SETUSED(theElement,0); SETTHEFLAG(theElement,0);
		for (i=0; i<EDGES_OF_ELEM(theElement); i++)
		{
			theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
							  CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
			ASSERT(theEdge != NULL);
			SETUSED(theEdge,0); SETTHEFLAG(theEdge,0);
		}
		if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
			for (i=0; i<SIDES_OF_ELEM(theElement); i++)
			{
				theVector = SVECTOR(theElement,i);
				SETUSED(theVector,0); SETTHEFLAG(theVector,0);
			}
	}
	/* to reset also nodes which are at corners of the boundary */
	/* reset nodes through the node list                        */
	for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
	{
		SETUSED(theNode,0); SETTHEFLAG(theNode,0);
		SETMODIFIED(theNode,0);
	}

	/* set FLAG for objects of vertical overlap */
	for (theElement=PFIRSTELEMENT(theGrid);
		 theElement!=NULL;
		 theElement=SUCCE(theElement))
	{
		if (PARTITION(theElement) == me) continue;

		/* check for horizontal ghost */
		hghost = 0;
		for (i=0; i<SIDES_OF_ELEM(theElement); i++)
		{
			theNeighbor = NBELEM(theElement,i);
			if (theNeighbor == NULL) continue;

			if (PARTITION(theNeighbor) == me)
			{
				hghost = 1;
				break;
			}
		}

		/* check for vertical ghost */
		vghost = 0;
		GetAllSons(theElement,SonList);
		for (i=0; SonList[i]!=NULL; i++)
		{
			if (PARTITION(SonList[i]) == me)
			{
				vghost = 1;
				break;
			}
		}

		/* one or both of vghost and hghost should be true here   */
		/* except for elements which will be disposed during Xfer */

		if (vghost) SETTHEFLAG(theElement,1);
		if (hghost) SETUSED(theElement,1);
		for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
		{
			theNode = CORNER(theElement,i);
			if (vghost) SETTHEFLAG(theNode,1);
			if (hghost) SETUSED(theNode,1);
		}
		for (i=0; i<EDGES_OF_ELEM(theElement); i++)
		{
			theEdge = GetEdge(CORNER_OF_EDGE_PTR(theElement,i,0),
							  CORNER_OF_EDGE_PTR(theElement,i,1));
			ASSERT(theEdge != NULL);
			if (vghost) SETTHEFLAG(theEdge,1);
			if (hghost) SETUSED(theEdge,1);
		}
		if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
			for (i=0; i<SIDES_OF_ELEM(theElement); i++)
			{
				theVector = SVECTOR(theElement,i);
				if (vghost) SETTHEFLAG(theVector,1);
				if (hghost) SETUSED(theVector,1);
			}
	}

	/* set USED flag for objects of master elements */
	/* reset FLAG for objects of master elements  */
	for (theElement=PFIRSTELEMENT(theGrid);
		 theElement!=NULL;
		 theElement=SUCCE(theElement))
	{
		if (PARTITION(theElement) != me) continue;

		SETUSED(theElement,0); SETTHEFLAG(theElement,0);
		for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
		{
			theNode = CORNER(theElement,i);
			SETUSED(theNode,0); SETTHEFLAG(theNode,0);
			SETMODIFIED(theNode,1);
		}
		for (i=0; i<EDGES_OF_ELEM(theElement); i++)
		{
			theEdge = GetEdge(CORNER_OF_EDGE_PTR(theElement,i,0),
							  CORNER_OF_EDGE_PTR(theElement,i,1));
			ASSERT(theEdge != NULL);
			SETUSED(theEdge,0); SETTHEFLAG(theEdge,0);
		}
		if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
			for (i=0; i<SIDES_OF_ELEM(theElement); i++)
			{
				theVector = SVECTOR(theElement,i);
				SETUSED(theVector,0); SETTHEFLAG(theVector,0);
			}
	}

	/* set object priorities for ghostelements */
	for (theElement=PFIRSTELEMENT(theGrid);
		 theElement!=NULL;
		 theElement=SUCCE(theElement))
	{
		if (PARTITION(theElement) == me) continue;

		if (USED(theElement) || THEFLAG(theElement))
		{
			prio = PRIO_CALC(theElement);
			PRINTDEBUG(gm,1,("SetGhostObjectPriorities(): e=" EID_FMTX " new prio=%d\n",
				EID_PRTX(theElement),prio))
			SETEPRIOX(theElement,prio);

			if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,ELEMVEC))
			{
				theVector = EVECTOR(theElement);
				SETPRIOX(theVector,prio);
			}
		}

		if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,EDGEVEC) || DIM==3)
		{
			/* set edge priorities */
			for (i=0; i<EDGES_OF_ELEM(theElement); i++)
			{

				theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
								  CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
				ASSERT(theEdge != NULL);

				if (USED(theEdge) || THEFLAG(theEdge))
				{
					PRINTDEBUG(dddif,3,(PFMT " dddif_SetGhostObjectPriorities():"
						" downgrade edge=" EDID_FMTX " from=%d to PrioHGhost\n",
						me,EDID_PRTX(theEdge),prio)); 

					EDGE_PRIORITY_SET(theGrid,theEdge,PRIO_CALC(theEdge));
				}
			}

			#ifdef __THREEDIM__
			/* if one(all) of the side nodes is (are) a hghost (vghost) node   */
			/* then its a hghost (vghost) side vector                          */
			if (VEC_DEF_IN_OBJ_OF_GRID(theGrid,SIDEVEC))
				for (i=0; i<SIDES_OF_ELEM(theElement); i++)
				{
					if (USED(theVector) || THEFLAG(theVector))
						SETPRIOX(theVector,PRIO_CALC(theVector));
				}
			#endif
		}
	}
	/* to set also nodes which are at corners of the boundary   */
	/* set them through the node list                           */
	for (theNode=PFIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
	{
		/* check if its a master node */
		if (USED(theNode) || THEFLAG(theNode))
		{
			PRINTDEBUG(dddif,3,(PFMT " dddif_SetGhostObjectPriorities():"
				" downgrade node=" ID_FMTX " from=%d to PrioHGhost\n",
				me,ID_PRTX(theNode),prio)); 

			/* set node priorities of node to ghost */
			NODE_PRIORITY_SET(theGrid,theNode,PRIO_CALC(theNode))
		}
		else if (MODIFIED(theNode) == 0)
		{
			/* this is a node of the boundary without connection to master elements */
			NODE_PRIORITY_SET(theGrid,theNode,PrioHGhost)
		}
	}

}


INT SetBorderPriorities (GRID *theGrid)
{
	DDD_IFAExecLocal(BorderNodeSymmIF,GRID_ATTR(theGrid),
		ComputeNodeBorderPrios);

	DDD_IFAExecLocal(BorderVectorSymmIF,GRID_ATTR(theGrid),
		ComputeVectorBorderPrios);

#ifdef __THREEDIM__
	DDD_IFAExecLocal(BorderEdgeSymmIF,GRID_ATTR(theGrid),
		ComputeEdgeBorderPrios);
#endif

	return(GM_OK);
}

INT SetGridBorderPriorities (GRID *theGrid)
{
	/* set border priorities on next higher level */
	if (SetBorderPriorities(UPGRID(theGrid)) != GM_OK) return(GM_FATAL);

	return(GM_OK);
}

void ConstructConsistentGrid (GRID *theGrid)
{
	INT		i,j,k,l,m,o;
	DOUBLE  fac,*local;
	ELEMENT *theElement,*theFather,*theNb;
	NODE	*theNode;
	EDGE	*theEdge;
	VERTEX	*theVertex;

	/* the setting of the priorities has to be done in two waves after */
	/* completion of the grid transfer, since                          */
	/* - decisions about vghost prio can only be done if all sons are  */
	/*   available in SetGhostObjectPriorities()                       */
	/* - setting of the border priorities can only be done if all      */
	/*   ghost objects have their proper priority                      */
	DDD_XferBegin();
	SetGhostObjectPriorities(theGrid);
	DDD_XferEnd();
	DDD_XferBegin();
	SetBorderPriorities(theGrid);
	DDD_XferEnd();

    #ifdef __TWODIM__
	for (theVertex = PFIRSTVERTEX(theGrid); theVertex != NULL;
		 theVertex = SUCCV(theVertex)) {
	    if (VXGHOST(theVertex)) 
		    VFATHER(theVertex) = NULL;
		else if (OBJT(theVertex) == BVOBJ) 
		    if (MOVED(theVertex))
			{
			    INT n;
				DOUBLE *x[MAX_CORNERS_OF_ELEM];
				
				theElement = VFATHER(theVertex);
				if (theElement == NULL) continue;
				HEAPFAULT(theElement);
				CORNER_COORDINATES(theElement,n,x);			
				UG_GlobalToLocal(n,(const DOUBLE **)x,
								 CVECT(theVertex),LCVECT(theVertex));
			}
	}
	#endif

	/* reconstruct VFATHER pointers */
	for (theElement = FIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
	{
		for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
		{
			theNode = CORNER(theElement,i);
			if (CORNERTYPE(theNode)) continue;

			theVertex = MYVERTEX(theNode);
			theFather = EFATHER(theElement);

/* this is too few for arbitrary load balancing, since 
	VFATHER pointer may have changed (970828 s.l.)
   			if (VFATHER(theVertex)==NULL || EPRIO(VFATHER(theVertex))==PrioHGhost)
*/
   			if (VFATHER(theVertex)==NULL || EPRIO(theFather)!=PrioHGhost)
			{
				switch (NTYPE(theNode))
				{
					case (MID_NODE):
					{
						INT             co0,co1;

						for (j=0; j<EDGES_OF_ELEM(theFather); j++)
						{
							theEdge = GetEdge(CORNER(theFather,CORNER_OF_EDGE(theFather,j,0)),
							                  CORNER(theFather,CORNER_OF_EDGE(theFather,j,1)));
							if (MIDNODE(theEdge) == theNode) break;
						}
/* here should be an assertion, but not in each situation the 
   midnode pointer is set (970829 s.l.)
						ASSERT(j<EDGES_OF_ELEM(theFather));
*/
						if (j>=EDGES_OF_ELEM(theFather))
						{
						  for (j=0; j<EDGES_OF_ELEM(theFather); j++)
							{
							  theEdge = GetEdge(CORNER(theFather,CORNER_OF_EDGE(theFather,j,0)),
												CORNER(theFather,CORNER_OF_EDGE(theFather,j,1)));
							  if (theEdge->midnode != NULL)
							  PRINTDEBUG(dddif,0,
										 (PFMT " ConstructConsistentGrid(): elem=" EID_FMTX 
										  " i=%d n1=" ID_FMTX " n2=" ID_FMTX " midnode= " ID_FMTX  "\n",
										  me,theFather,EID_PRTX(theFather),j,
										  ID_PRTX(NBNODE(LINK0(theEdge))),
										  ID_PRTX(NBNODE(LINK1(theEdge))),
										  ID_PRTX(theEdge->midnode)))
							}



							PRINTDEBUG(dddif,0,
									   ("ConstructConsistentGrid(): WARN "
										" theNode= " ID_FMTX 
										" vertex= " VID_FMTX 
										" recalculation of VFATHER impossible\n",
										ID_PRTX(NBNODE(LINK0(theEdge))),
										VID_PRTX(theVertex)));
						  break;
						}

							
						/* reconstruct local coordinates of vertex */
						co0 = CORNER_OF_EDGE(theFather,j,0);
						co1 = CORNER_OF_EDGE(theFather,j,1);

						/* local coordinates have to be local towards pe */


						V_DIM_LINCOMB(0.5, LOCAL_COORD_OF_ELEM(theFather,co0),
									  0.5, LOCAL_COORD_OF_ELEM(theFather,co1),
									  LCVECT(theVertex));
						SETONEDGE(theVertex,j);

						break;
					}

					#ifdef __THREEDIM__
					case (SIDE_NODE):
						/* always compute new coords for this case! */
						k =  GetSideIDFromScratch(theElement,theNode);
						ASSERT(k < SIDES_OF_ELEM(theFather));

						SETONSIDE(theVertex,k);

						m = CORNERS_OF_SIDE(theFather,k);
						local = LCVECT(theVertex);
						fac = 1.0 / m;
						V_DIM_CLEAR(local);
						for (o=0; o<m; o++)
						{
							l = CORNER_OF_SIDE(theFather,k,o);
							V_DIM_LINCOMB(1.0,local,1.0,
										  LOCAL_COORD_OF_ELEM(theFather,l),local);
						}
						V_DIM_SCALE(fac,local);

						theNb = NBELEM(theFather,k);
						if (theNb != NULL)
						{
							for (j=0; j<SIDES_OF_ELEM(theNb); j++)
							{
								if (NBELEM(theNb,j) == theFather) break;
							}
							ASSERT(j < SIDES_OF_ELEM(theNb));
							SETONNBSIDE(theVertex,j);
						}
						else SETONNBSIDE(theVertex,MAX_SIDES_OF_ELEM);
						break;

					#endif
					case (CENTER_NODE):
					case (LEVEL_0_NODE):
						/* nothing to do */
						break;

					case (CORNER_NODE):
					default:
						assert(0);
						break;
				}
				VFATHER(theVertex) = theFather;

				if (OBJT(theVertex) == BVOBJ) 
				    if (MOVED(theVertex)) {
					    INT n;
						DOUBLE *x[MAX_CORNERS_OF_ELEM];
				
						CORNER_COORDINATES(theFather,n,x);			
						UG_GlobalToLocal(n,(const DOUBLE **)x,
										 CVECT(theVertex),LCVECT(theVertex));
					}
			}
		}
	}
}

INT CheckProcListCons (int *proclist, int uniqueTag)
{
	int nunique = 0;

	/* check uniqueness */
	while (*proclist != -1)
	{
		if (*(proclist+1) == uniqueTag) nunique++;
		proclist += 2;
	}

	/* nunique must be 1 for master elements   */
	/* nunique can  be 0/1 for (inner) nodes   */
	/*   with PrioBorder/PrioMaster            */
	return (nunique);
}

INT ListProcList (int *proclist, int uniqueTag)
{
	while (*proclist != -1)
	{
		if (*(proclist+1) == uniqueTag)
			UserWriteF(" proc=%d",*proclist);
		proclist += 2;
	}
	return(0);
}

INT CheckVectorPrio (ELEMENT *theElement, VECTOR *theVector)
{
	INT nmaster;
	INT nerrors = 0;

	/* check vector prio */
	CHECK_OBJECT_PRIO(theVector,PRIO,MASTER,GHOST,VINDEX,"Vector",nerrors)

	/* master copy has to be unique */
	if ((nmaster = CheckProcListCons(PROCLIST(theVector),PrioMaster)) > 1)	
	{
		UserWriteF("NODE=" ID_FMTX " ERROR: master copy not unique, nmaster=%d:",
			ID_PRTX(theVector),nmaster);
		ListProcList(PROCLIST(theVector),PrioMaster);
		UserWriteF("\n");
		nerrors++;
	}

	return(nerrors);
}

INT CheckNodePrio (ELEMENT *theElement, NODE *theNode)
{
	INT nmaster;
	INT nerrors = 0;

	/* check node prio */
	CHECK_OBJECT_PRIO(theNode,PRIO,MASTER,GHOST,ID,"NODE",nerrors)

	/* master copy has to be unique */
	if ((nmaster = CheckProcListCons(PROCLIST(theNode),PrioMaster)) > 1)	
	{
		UserWriteF("NODE=" ID_FMTX " ERROR: master copy not unique, nmaster=%d:",
			ID_PRTX(theNode),nmaster);
		ListProcList(PROCLIST(theNode),PrioMaster);
		UserWriteF("\n");
		nerrors++;
	}

	if (dddctrl.nodeData)
		nerrors += CheckVectorPrio(theElement,NVECTOR(theNode));

	return(nerrors);
}


INT CheckEdgePrio (ELEMENT *theElement, EDGE *theEdge)
{
	INT nmaster;
	INT nerrors = 0;

	#ifdef __THREEDIM__
	/* check edge prio */
	CHECK_OBJECT_PRIO(theEdge,PRIO,MASTER,GHOST,ID,"EDGE",nerrors)

	/* master copy has to be unique */
	if ((nmaster = CheckProcListCons(PROCLIST(theEdge),PrioMaster)) > 1)	
	{
		UserWriteF("EDGE=" EDID_FMTX " ERROR: master copy not unique, nmaster=%d:",
			EDID_PRTX(theEdge),nmaster);
		ListProcList(PROCLIST(theEdge),PrioMaster);
		UserWriteF("\n");
		nerrors++;
	}
	#endif

	if (dddctrl.edgeData)
	    if (EDVECTOR(theEdge) != NULL)
		    nerrors += CheckVectorPrio(theElement,EDVECTOR(theEdge));

	return(nerrors);
}

INT CheckElementPrio (ELEMENT *theElement)
{
	INT		i,nmaster,prio,valid_copy;
	INT		nerrors = 0;
	NODE	*theNode;
	EDGE	*theEdge;
	ELEMENT *SonList[MAX_SONS];

	if (PARTITION(theElement)==me && !EMASTER(theElement))
	{
		UserWriteF(PFMT "#FATAL# MASTER ELEM=" EID_FMTX " has WRONG part=%d prio=%d\n",
			me,EID_PRTX(theElement),PARTITION(theElement),EPRIO(theElement));
		nerrors++;
	}
	if (PARTITION(theElement)!=me && !EGHOST(theElement))
	{
		UserWriteF(PFMT "#FATAL# GHOST ELEM=" EID_FMTX " has WRONG part=%d prio=%d\n",
			me,EID_PRTX(theElement),PARTITION(theElement),EPRIO(theElement));
		nerrors++;

		/* test ghost prio */
		prio = 0;
		for (i=0; i<SIDES_OF_ELEM(theElement); i++)
		{
			if (EMASTER(NBELEM(theElement,i))) prio = PrioHGhost; 
		}
		if (GetSons(theElement,SonList) != 0) RETURN(1);
		if (SonList[0] != NULL) prio += PrioVGhost;

		if (EPRIO(theElement) != prio)
		{
			UserWriteF(PFMT "ERROR GHOST ELEM=" EID_FMTX 
				" has WRONG prio=%d should be prio=%d\n",
				me,EID_PRTX(theElement),EPRIO(theElement),prio);
			nerrors++;
		}
	}

	/* check element prio */
	CHECK_OBJECT_PRIO(theElement,EPRIO,EMASTER,EGHOST,EID,"ELEM",nerrors)

	/* master copy has to be unique */
	if ((nmaster = CheckProcListCons(EPROCLIST(theElement),PrioMaster)) != 1)	
	{
		UserWriteF("ELEM=" EID_FMTX " ERROR: master copy not unique, nmaster=%d:",
			EID_PRTX(theElement),nmaster);
		ListProcList(EPROCLIST(theElement),PrioMaster);
		UserWriteF("\n");
		nerrors++;
	}

	/* hghost copy needs to a master neighbor */
	if (EHGHOST(theElement))
	{
		valid_copy = 0;	
		for (i=0; i<SIDES_OF_ELEM(theElement); i++)
		{
			if (NBELEM(theElement,i)!=NULL && EMASTER(NBELEM(theElement,i)))
				valid_copy = 1;
		}
		if (!valid_copy) 
		{
			UserWriteF("ELEM=" EID_FMTX " ERROR: hghost copy with no master neighbor!\n",
				EID_PRTX(theElement));
			nerrors++;
		}
	}

	/* vghost copy needs to a master Son */
	if (EVGHOST(theElement))
	{
		if (GetSons(theElement,SonList) != 0) RETURN(1);
		if (SonList[0] == NULL) valid_copy = 0;
		else 					valid_copy = 1;
		if (!valid_copy) 
		{
			UserWriteF("ELEM=" EID_FMTX " ERROR: vghost copy with no master son!\n",
				EID_PRTX(theElement));
			nerrors++;
		}
	}

	if (dddctrl.elemData)
		nerrors += CheckVectorPrio(theElement,EVECTOR(theElement));

	if (dddctrl.sideData)
	{
		for (i=0; i<SIDES_OF_ELEM(theElement); i++)
			nerrors += CheckVectorPrio(theElement,SVECTOR(theElement,i));
	}

	for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
	{
		theNode = CORNER(theElement,i);
		nerrors += CheckNodePrio(theElement,theNode);
	}

	for (i=0; i<EDGES_OF_ELEM(theElement); i++)
	{
		theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
						  CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
		ASSERT(theEdge != NULL);
		nerrors += CheckEdgePrio(theElement,theEdge);
	}

	return (nerrors);
}

INT CheckInterfaces(GRID *theGrid)
{
	INT		i,j;
	ELEMENT	*theElement;
	NODE	*theNode;
	EDGE	*theEdge;
	VECTOR	*theVector;
	int		nerrors = 0;

	/* reset USED flag of all grid objects  */
	/* set USED flag of master grid objects */
	for (j=0; j<2; j++)
	{
		for (theElement =(j==0 ? PFIRSTELEMENT(theGrid) : FIRSTELEMENT(theGrid)); 
			 theElement!=NULL;
			 theElement=SUCCE(theElement))
		{
			SETUSED(theElement,j);
			if (dddctrl.elemData)
				SETUSED(EVECTOR(theElement),j);
			if (dddctrl.sideData)
			{
				for (i=0; i<SIDES_OF_ELEM(theElement); i++)
					SETUSED(SVECTOR(theElement,i),j);
			}

			for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
			{
				theNode = CORNER(theElement,i);
				SETUSED(theNode,j);
				if (dddctrl.nodeData)
					SETUSED(NVECTOR(theNode),j);
				SETUSED(MYVERTEX(theNode),j);
			}

			for (i=0; i<EDGES_OF_ELEM(theElement); i++)
			{
				theEdge = GetEdge(CORNER(theElement,CORNER_OF_EDGE(theElement,i,0)),
								  CORNER(theElement,CORNER_OF_EDGE(theElement,i,1)));
				ASSERT(theEdge != NULL);
				SETUSED(theEdge,j);
				if (dddctrl.edgeData)
				    if (EDVECTOR(theEdge) != NULL)
					    SETUSED(EDVECTOR(theEdge),j);
			}
		}
	}

	/* check validity of priorities */
	for (theElement=PFIRSTELEMENT(theGrid); theElement!=NULL;
		 theElement=SUCCE(theElement))
	{
		nerrors += CheckElementPrio(theElement);
	}

	/* check ddd interface consistency */
	DDD_SetOption(OPT_QUIET_CONSCHECK, OPT_ON);
	nerrors += DDD_ConsCheck();
	DDD_SetOption(OPT_QUIET_CONSCHECK, OPT_OFF);

	return(nerrors);
}

/****************************************************************************/

#endif

