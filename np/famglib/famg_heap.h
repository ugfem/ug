// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:      famg_heap.h													*/
/*																			*/
/* Purpose:   famg heap class												*/
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

#ifndef __FAMG_HEAP__
#define __FAMG_HEAP__

/* RCS_ID
   $Header$
 */

#define FAMGMAXSTACK 8

class FAMGHeap
{
public:
  FAMGHeap(unsigned long );
  ~FAMGHeap();

  void *GetMem(unsigned long, int);
  int Mark(int);
  int Release(int);
  void PrintInfo();
private:
  void *buffer;
  unsigned long top, bottom, info_max_bottom, info_max_top, info_min_free, info_size;
  int ntop, nbottom;
  unsigned long topstack[FAMGMAXSTACK], bottomstack[FAMGMAXSTACK];
#ifdef USE_UG_DS
  int FAMGHeapMarkKey;
#endif
};

#define FAMGALIGNMENT 8
#define FAMGCEIL(n) ((n)+((FAMGALIGNMENT-((n)& (FAMGALIGNMENT-1)))& (FAMGALIGNMENT-1)))

#define FAMG_FROM_TOP  1
#define FAMG_FROM_BOTTOM  2


/****************************************************************************/
/*                                                                          */
/* Functions                                                                */
/*                                                                          */
/****************************************************************************/
void *FAMGGetMem(unsigned long size, int mode);
void FAMGSetHeap(FAMGHeap *ptr);
void FAMGFreeHeap();
FAMGHeap *FAMGGetHeap();
int FAMGMarkHeap(int mode);
int FAMGReleaseHeap(int mode);

#endif
