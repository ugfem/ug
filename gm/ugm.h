// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  ugm.h                                                                                                                 */
/*																			*/
/* Purpose:   unstructured grid manager header file                                             */
/*			  internal interface in grid manager module                                     */
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de			                        */
/*																			*/
/* History:   09.03.92 begin, ug version 2.0								*/
/*			  13.12.94 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __UGM__
#define __UGM__

#ifndef __COMPILER__
#include "compiler.h"
#endif

#ifndef __SWITCH__
#include "switch.h"
#endif

#ifndef __GM__
#include "gm.h"
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

#define MAX_PAR_DIST    1.0E-5          /* max.dist between different parameter */

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* init */
INT              InitUGManager                  (void);

/* object handling */
INT              GetFreeOBJT                    (void);
INT              ReleaseOBJT                    (INT type);
void            *GetMemoryForObject             (MULTIGRID *theMG, INT size, INT type);
INT              PutFreeObject                  (MULTIGRID *theMG, void *object, INT size, INT type);
#ifdef ModelP
void            *GetMemoryLocal                 (MULTIGRID *theMG, INT size, INT type);
INT              PutFreeObjectLocal             (MULTIGRID *theMG, void *object, INT size, INT type);
#endif

/* create basic objects */
VERTEX          *CreateBoundaryVertex   (GRID *theGrid, VERTEX *after);
VERTEX          *CreateInnerVertex              (GRID *theGrid, VERTEX *after);
VSEGMENT        *CreateVertexSegment    (GRID *theGrid, VERTEX *vertex);
NODE            *CreateNode                     (GRID *theGrid, NODE *after);
EDGE            *CreateEdge                     (GRID *theGrid, NODE *from, NODE *to);
ELEMENT         *CreateBoundaryElement  (GRID *theGrid, ELEMENT *after, INT tag);
ELEMENT         *CreateInnerElement     (GRID *theGrid, ELEMENT *after, INT tag);
ELEMENTSIDE *CreateElementSide          (GRID *theGrid);
INT         CreateSonElementSide    (GRID *theGrid, ELEMENT *theElement,
                                     INT side, ELEMENT *theSon, INT son_side);

GRID            *CreateNewLevel                 (MULTIGRID *theMG);

/* dispose basic objects */
INT              DisposeEdge                    (GRID *theGrid, EDGE *theEdge);
INT              DisposeNode                    (GRID *theGrid, NODE *theNode);
INT              DisposeVertex                  (GRID *theGrid, VERTEX *theVertex);
INT              DisposeElementSide             (GRID *theGrid, ELEMENTSIDE *theElementSide);
INT              DisposeElement                 (GRID *theGrid, ELEMENT *theElement);
INT              DisposeTopLevel                (MULTIGRID *theMG);
INT                      DisposeEdgesFromElement(GRID *theGrid, ELEMENT *theElement);

/* miscellaneous */
INT              FindNeighborElement    (const ELEMENT *theElement, INT Side, ELEMENT **theNeighbor, INT *NeighborSide);
INT              PointInElement                 (const COORD*, const ELEMENT *theElement);
VIRT_HEAP_MGMT *GetGenMGUDM             (void);

NODE            *CreateMidNode                  (GRID *theGrid,ELEMENT *theElement,INT side,NODE *after);
NODE        *CreateCenterNode       (GRID *theGrid, ELEMENT *theElement);

#ifdef __TWODIM__
EDGE            *CreateAuxEdge                  (GRID *theGrid, NODE *from, NODE *to);
INT              DisposeAuxEdges                (GRID *theGrid);
#endif

#ifdef __THREEDIM__
INT          DisposeEdgesFromElement(GRID *theGrid, ELEMENT *theElement);
NODE            *CreateSideNode                 (GRID *theGrid, ELEMENT *theElement, INT side);
#endif

#endif
