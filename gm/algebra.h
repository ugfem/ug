// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  algebra.h                                                                                                     */
/*																			*/
/* Purpose:   header for algebraic structures								*/
/*			  internal interface for grid manager module					*/
/*																			*/
/* Author:	  Klaus Johannsen												*/
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 294										*/
/*			  6900 Heidelberg												*/
/*			  internet: johannse@iwr1.iwr.uni-heidelberg.de                                 */
/*																			*/
/* History:   1.12.93 begin, ug 3d											*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifdef __version3__

#ifndef __ALGEBRA__
#define __ALGEBRA__


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

/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

extern INT MatrixType[MAXVECTORS][MAXVECTORS];

/****************************************************************************/
/*																			*/
/* control word definitions                                                                                             */
/*																			*/
/****************************************************************************/

/* vector */
#define VBUILDCON_CE                            1
#define VBUILDCON_SHIFT                         3
#define VBUILDCON_LEN                           1
#define VBUILDCON(p)                            CW_READ(p,VBUILDCON_CE)
#define SETVBUILDCON(p,n)                       CW_WRITE(p,VBUILDCON_CE,n)

/* matrix */

/* element */
#define EBUILDCON_CE                            52
#define EBUILDCON_SHIFT                         11
#define EBUILDCON_LEN                           1
#define EBUILDCON(p)                            CW_READ(p,EBUILDCON_CE)
#define SETEBUILDCON(p,n)                       CW_WRITE(p,EBUILDCON_CE,n)


/****************************************************************************/
/*																			*/
/* function declarations													*/
/*																			*/
/****************************************************************************/

/* basic create and dispose functions */
INT             CreateVector                                    (GRID *theGrid, VECTOR *After, INT VectorType, VECTOR **VectorHandle);
CONNECTION      *CreateConnection                               (GRID *theGrid, VECTOR *from, VECTOR *to);
INT             DisposeVector                                   (GRID *theGrid, VECTOR *theVector);
INT             DisposeConnection                               (GRID *theGrid, CONNECTION *theConnection);

/* more create and dispose */
INT                     MGCreateConnection                              (MULTIGRID *theMG);
INT             CreateConnectionsInNeighborhood (GRID *theGrid, ELEMENT *theElement);
INT             InsertedElementCreateConnection (GRID *theGrid, ELEMENT *theElement);
INT             DisposeDoubledSideVector                (GRID *theGrid, ELEMENT *Elem0, INT Side0, ELEMENT *Elem1, INT Side1);
INT             DisposeConnectionFromVector     (GRID *theGrid, VECTOR *theVector);
INT             DisposeConnectionFromElement    (GRID *theGrid, ELEMENT *theElement);
INT             DisposeConnectionsInNeighborhood(GRID *theGrid, ELEMENT *theElement);

/* query functions */
INT             GetVectorsOfElement                     (const ELEMENT *theElement, INT *cnt, VECTOR **vList);
INT             GetVectorsOfSides                               (const ELEMENT *theElement, INT *cnt, VECTOR **vList);
INT             GetVectorsOfEdges                               (const ELEMENT *theElement, INT *cnt, VECTOR **vList);
INT             GetVectorsOfNodes                               (const ELEMENT *theElement, INT *cnt, VECTOR **vList);
INT             GetElementInfoFromSideVector    (const VECTOR *theVector, ELEMENT **Elements, INT *Sides);

/* gridwise functions */
INT             GridCreateConnection                    (GRID *theGrid);

/* check algebra */
ELEMENT         *ElementCheckConnection                 (GRID *theGrid, ELEMENT *theElement);

/* determination of vector classes */
INT             ClearVectorClasses                              (GRID *theGrid);
INT             SeedVectorClasses                               (ELEMENT *theElement);
INT             PropagateVectorClasses                  (GRID *theGrid);
INT             ClearNextVectorClasses                  (GRID *theGrid);
INT             SeedNextVectorClasses                   (ELEMENT *theElement);
INT             PropagateNextVectorClasses              (GRID *theGrid);
INT             MaxNextVectorClass                              (ELEMENT *theElement);

/* miscellaneous routines */
INT             PrepareAlgebraModification              (MULTIGRID *theMG);
INT             MoveVector                                              (GRID *theGrid, VECTOR *moveVector, VECTOR *destVector, INT after);

/* Initialization */
INT             InitAlgebra                                     (void);

#endif
#endif
