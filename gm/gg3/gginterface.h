// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      gginterface.h                                                 */
/*                                                                          */
/* Purpose:   interface header file for netgen                                                  */
/*                                                                          */
/* Author:    Christian Wieners                                             */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart, Germany										*/
/*			  email: ug@ica3.uni-stuttgart.de		                                        */
/*																			*/
/* History:   18 March 96 begin, ug version 3.2                             */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __GGINTERFACE__
#define __GGINTERFACE__

/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

int AddInnerNode (double x, double y, double z);
int AddInnerNode2ug (double x, double y, double z);
int AddTetrahedron (int node0, int node1, int node2, int node3);
int AddSurfaceNode (int nodeid, double x, double y, double z);
int AddSurfaceTriangle (int node0, int node1, int node2);
int AddSurfaceTriangle2ug (int node0, int node1, int node2);
int InitNetgen (char *rulefilename);
int StartNetgen (double h,int smooth,int display);

int AddGeomPoint (int nodeid, double x, double y, double z);
int AddGeomElement (int node0, int node1, int node2);
int AddLinePoint (int id, double x, double y, double z);
int AddLineSegment (int i1,int i2);
int InitSurfaceNetgen (char *rulefilename);
int StartSurfaceNetgen (double h,int smooth,int display);
#endif
