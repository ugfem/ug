// NOTE: The current revision of this file was left untouched when the DUNE source files were reindented!
// NOTE: It contained invalid syntax that could not be processed by uncrustify.

/****************************************************************************/
/*																			*/
/* File:	  gridcons.c													*/
/*																			*/
/* Purpose:   functions for managing consistency of distributed grids       */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken										*/
/*			  Institut fuer Computeranwendungen III 						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
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
#include "refine.h"
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

REP_ERR_FILE;

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*
   ConstructConsistentGrid - 

   SYNOPSIS:
   void ConstructConsistentGrid (GRID *theGrid);

   PARAMETERS:
.  theGrid

   DESCRIPTION:
   Provide a consistent grid. Do not call for a single grid level, but 
   always for the whole multigrid from bottom to top, since otherwise
   consistency is not ensured!

   RETURN VALUE:
   void
*/
/****************************************************************************/

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

	/* this is the simplest fix for VFATHER zombies  */
	/* just reset all VFATHER pointers and set them  */
	/* only by master nodes of this or upper levels. */
	/* A more complicated fix would be to set the    */
	/* priorities of the vertices correctly.         */
	/* (980126 s.l.)                                 */
	for (theVertex = PFIRSTVERTEX(theGrid); theVertex != NULL;
		 theVertex = SUCCV(theVertex)) {
		VFATHER(theVertex) = NULL;
/*
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
*/
	}

	/* reconstruct VFATHER pointers and                */
	/* make ghost neighborships symmetric (only for 3d)*/
	for (theElement = PFIRSTELEMENT(theGrid); theElement!=NULL; theElement=SUCCE(theElement))
	{
/* TODO: delete now done in ElementObjMkCons()
		#ifdef __THREEDIM__
		if (EVGHOST(theElement))
		{
			ELEMENT *NbElement; 

			for (i=0; i<SIDES_OF_ELEM(theElement); i++)
			{
				NbElement = NBELEM(theElement,i);
				for (j=0; j<SIDES_OF_ELEM(NbElement); j++)
				{
					if (NBELEM(NbElement,j) == theElement) break;
				}
				if (j>=SIDES_OF_ELEM(NbElement))
					SET_NBELEM(theElement,i,NULL);
			}
		}
		#endif 
*/

		theFather = EFATHER(theElement);

		/* no reconstruction of VFATHER possible */
		if (theFather == NULL) continue;

		for (i=0; i<CORNERS_OF_ELEM(theElement); i++)
		{
			theNode = CORNER(theElement,i);
			if (CORNERTYPE(theNode)) continue;

			theVertex = MYVERTEX(theNode);

/* this is too few for arbitrary load balancing, since 
	VFATHER pointer may have changed (970828 s.l.)
   			if (VFATHER(theVertex)==NULL || EPRIO(VFATHER(theVertex))==PrioHGhost)
*/
/* this is too few for arbitrary load balancing, since
	VFATHER pointer may be already a zombie pointer (980126 s.l.)
   			if (VFATHER(theVertex)==NULL || EPRIO(theFather)!=PrioHGhost)
*/
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
							  PRINTDEBUG(dddif,1,
										 (PFMT " ConstructConsistentGrid(): elem=" EID_FMTX 
										  " i=%d n1=" ID_FMTX " n2=" ID_FMTX " midnode= " ID_FMTX  "\n",
										  me,theFather,EID_PRTX(theFather),j,
										  ID_PRTX(NBNODE(LINK0(theEdge))),
										  ID_PRTX(NBNODE(LINK1(theEdge))),
										  ID_PRTX(theEdge->midnode)))
							}



							PRINTDEBUG(dddif,1,
									   ("ConstructConsistentGrid(): WARN "
										" theNode= " ID_FMTX 
										" vertex= " VID_FMTX 
										" recalculation of VFATHER impossible\n",
										ID_PRTX(NBNODE(LINK0(theEdge))),
										VID_PRTX(theVertex)));
							/* if it couldn't  be recalculated reset it */
							VFATHER(theVertex) = NULL;
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
				#ifdef Debug
				if (theFather != NULL) HEAPFAULT(theFather);
				#endif
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

#endif

