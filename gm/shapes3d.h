// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  shapes.h														*/
/*																			*/
/* Purpose:   header file for shape functions								*/
/*																			*/
/* Author:	  Peter Bastian                                                                                                 */
/*			  Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen	*/
/*			  Universitaet Heidelberg										*/
/*			  Im Neuenheimer Feld 368										*/
/*			  6900 Heidelberg												*/
/*			  internet: ug@ica3.uni-stuttgart.de					*/
/*																			*/
/* History:   29.01.92 begin, ug version 2.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* auto include mechanism and other include files							*/
/*																			*/
/****************************************************************************/

#ifndef __SHAPES3D__
#define __SHAPES3D__

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

#define LOCAL_TO_GLOBAL(A,B,C)          {C[0] = (1.0-B[0]-B[1]-B[2])*A[0][0]+B[0]*A[1][0]+B[1]*A[2][0]+B[2]*A[3][0];\
                                         C[1] = (1.0-B[0]-B[1]-B[2])*A[0][1]+B[0]*A[1][1]+B[1]*A[2][1]+B[2]*A[3][1];\
                                         C[2] = (1.0-B[0]-B[1]-B[2])*A[0][2]+B[0]*A[1][2]+B[1]*A[2][2]+B[2]*A[3][2];}


/****************************************************************************/
/*																			*/
/* data structures exported by the corresponding source file				*/
/*																			*/
/****************************************************************************/

typedef INT (*UpwindProcPtr)(COORD **, COORD_VECTOR *,DOUBLE_VECTOR *, COORD_VECTOR *);

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

DOUBLE  N                                       (const INT i, const COORD_VECTOR LocalCoord);
INT     GlobalToLocal3d         (const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);
INT     TetraDerivative         (const COORD **theCorners, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM]);
INT     TetraVolume             (const COORD **theCorners, COORD *volume);
INT     FV_TetInfo                      (const COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM]);
INT     FV_TetInfo_for_conv (const COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM], COORD_VECTOR conv);
INT     FV_AliTetInfo           (const COORD **CornerPoints, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv, COORD_VECTOR GIP[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     Side_TetInfo            (const COORD **theCorners, INT side, COORD_VECTOR Area, COORD_VECTOR GIP[3]);
INT     GetSkewedUIP            (const COORD_VECTOR *theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], const DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     GFUIP                           (const COORD **theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     GCUIP                           (const COORD **theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     COPYIP                          (const COORD **theCorners, const COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     TetraSideNormals        (const COORD **theCorners, COORD_VECTOR theNormals[MAX_SIDES_OF_ELEM]);
INT     TetMaxSideAngle         (const COORD **theCorners, COORD *MaxAngle);
INT     TetAngleAndLength       (const COORD **theCorners, COORD *Angle, COORD *Length);

#endif
