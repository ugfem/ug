// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      heap.h														*/
/*																			*/
/* Purpose:   cmg heap class												*/
/*																			*/
/* Author:    Christian Wagner												*/
/*			  Institut fuer Computeranwendungen  III						*/
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  internet: chris@ica3.uni-stuttgart.de							*/
/*																			*/
/*																			*/
/* History:   November 97 begin, Stuttgart									*/
/*			  August 98 integration into ug (Christian Wrobel)				*/
/*																			*/
/* Remarks:																	*/
/*																			*/
/****************************************************************************/

#ifndef __CMG_HEAP__
#define __CMG_HEAP__

/* RCS_ID
   $Header$
 */

#define CMGMAXSTACK 8

class CMGHeap
{
public:
  void *GetMem(unsigned long, int);
  int Mark(int);
  int Release(int);
  CMGHeap(unsigned long );
  ~CMGHeap();
private:
  void *buffer;
  unsigned long top, bottom;
  int ntop, nbottom;
  unsigned long topstack[CMGMAXSTACK], bottomstack[CMGMAXSTACK];
};

#define CMGALIGNMENT 8
#define CMGCEIL(n) ((n)+((CMGALIGNMENT-((n)& (CMGALIGNMENT-1)))& (CMGALIGNMENT-1)))

#define CMG_FROM_TOP  1
#define CMG_FROM_BOTTOM  2


/****************************************************************************/
/*                                                                          */
/* Functions                                                                */
/*                                                                          */
/****************************************************************************/
void *CMGGetMem(unsigned long size, int mode);
void CMGSetHeap(CMGHeap *ptr);
CMGHeap *CMGGetHeap();
int CMGMarkHeap(int mode);
int CMGReleaseHeap(int mode);

#endif
