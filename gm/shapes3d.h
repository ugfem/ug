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

DOUBLE  N                                       (INT i, COORD_VECTOR LocalCoord);
INT     GlobalToLocal3d         (COORD **Corners,COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);
INT     TetraDerivative         (const COORD **theCorners, COORD_VECTOR theGradient[MAX_CORNERS_OF_ELEM]);
INT     TetraVolume             (COORD **theCorners, COORD *volume);
INT     FV_TetInfo                      (COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM]);
INT     FV_TetInfo_for_conv (COORD **theCorners, COORD_VECTOR Area[MAX_EDGES_OF_ELEM], COORD_VECTOR GIP[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM], COORD_VECTOR conv);
INT     Side_TetInfo            (COORD **theCorners, INT side, COORD_VECTOR Area, COORD_VECTOR GIP[3]);
INT     GSUIP                           (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     GFUIP                           (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     GCUIP                           (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     COPYIP                          (COORD **theCorners, COORD_VECTOR LIP[MAX_EDGES_OF_ELEM], DOUBLE_VECTOR conv[MAX_EDGES_OF_ELEM], COORD_VECTOR LUIP[MAX_EDGES_OF_ELEM]);
INT     TetraSideNormals        (COORD **theCorners, COORD_VECTOR theNormals[MAX_SIDES_OF_ELEM]);
INT     TetMaxSideAngle         (COORD **theCorners, COORD *MaxAngle);
INT     TetAngleAndLength       (COORD **theCorners, COORD *Angle, COORD *Length);

#endif
