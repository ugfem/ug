// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      fifo.h														*/
/*                                                                          */
/* Purpose:   cmg fifo class												*/
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

#ifndef __CMG_FIFO__
#define __CMG_FIFO__

/* RCS_ID
   $Header$
 */

// integer fifo. for pointers: cast pointers to integer.

class CMGFifo
{
public:
  CMGFifo(void **buffer, int el_size);
  void Clear();
  int Empty();
  int Full();
  int In(void *el);
  void *Out();

private:
  int start;
  int end;
  int size;
  int used;
  void **elements;
};


inline CMGFifo::CMGFifo(void **buffer, int total_size)
{
  size = total_size / sizeof(int);
  elements =  buffer;
  start = 0;
  end = 0;
  used = 0;
}

inline void CMGFifo::Clear()
{
  start = 0;
  end = 0;
  used = 0;
}

inline int CMGFifo::Empty()
{
  return (used == 0);
}

inline int CMGFifo::Full()
{
  return (used == size);
}

inline int CMGFifo::In(void *el)
{
  if(used < size)
  {
    elements[end] = el;
    end = (end+1)%size;
    used++;
    return 0;
  }
  else
  {
    return 1;
  }
}

inline void* CMGFifo::Out()
{
  if (used == 0) return NULL;
  int i = start;
  start = (start+1)%size;
  used--;
  return elements[i];
}


#endif
