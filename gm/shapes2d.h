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

#ifndef __SHAPES2D__
#define __SHAPES2D__

/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define DNDS(n,i,s,t,r)                 if (n==3)\
  {\
    switch (i)\
    {\
    case 0 : r=-1.0; break;\
    case 1 : r=1.0; break;\
    case 2 : r=0.0; break;\
    }\
  }\
  else if (n==4)\
  {\
    switch (i)\
    {\
    case 0 : r=-0.25*(1-t); break;\
    case 1 : r=0.25*(1-t); break;\
    case 2 : r=0.25*(1+t); break;\
    case 3 : r=-0.25*(1+t); break;\
    }\
  }

#define DNDT(n,i,s,t,r)             if (n==3)\
  {\
    switch (i)\
    {\
    case 0 : r=-1; break;\
    case 1 : r=0; break;\
    case 2 : r=1; break;\
    }\
  }\
  else if (n==4)\
  {\
    switch (i)\
    {\
    case 0 : r=-0.25*(1-s); break;\
    case 1 : r=-0.25*(1+s); break;\
    case 2 : r=0.25*(1+s); break;\
    case 3 : r=0.25*(1-s); break;\
    }\
  }



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

DOUBLE  N                               (int n, int i, DOUBLE s, DOUBLE t);
DOUBLE  dNds                    (int n, int i, DOUBLE s, DOUBLE t);
DOUBLE  dNdt                    (int n, int i, DOUBLE s, DOUBLE t);

int     Derivatives     (int n, const DOUBLE *px, const DOUBLE *py, DOUBLE ips, DOUBLE ipt, DOUBLE *dNdx, DOUBLE *dNdy, DOUBLE *detJ);
INT             Gradients               (INT n, const COORD **theCorners, DOUBLE ips, DOUBLE ipt, DOUBLE_VECTOR Gradient[MAX_CORNERS_OF_ELEM], DOUBLE *DetJ);

INT             L2GDerivative2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD *Derivative);
INT             LocalToGlobal2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR GlobalCoord);
INT     GlobalToLocal2d (INT n,const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);
INT             specialGlobalToLocal2d (INT n, const COORD **Corners, const COORD_VECTOR EvalPoint, COORD_VECTOR LocalCoord);

#endif
