// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      sll.h                                                         */
/*                                                                          */
/* Purpose:   header file for sll templates                                 */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            internet: birken@ica3.uni-stuttgart.de                        */
/*                                                                          */
/* History:   960826 kb  begin                                              */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __SLL_H__
#define __SLL_H__


/****************************************************************************/

/*
        which structure components are needed for each item?
 */
#define SLL_INFO(T)   \
  int sll_n;                   /* unique index number */  \
  struct T *sll_next           /* linked list         */



/****************************************************************************/

/*
        macros for template support: variables
 */
#define _list(T) list ## T
#define list(T) _list(T)

#define _free(T) free ## T
#define free(T) _free(T)

#define _n(T) n ## T
#define n(T) _n(T)


/*
        macros for template support: functions
 */
#define _New(T) New ## T
#define New(T) _New(T)

#define _Free(T) Free ## T
#define Free(T) _Free(T)

#define _SortedArray(T) SortedArray ## T
#define SortedArray(T) _SortedArray(T)

#define _sort_OrigOrder(T) sort_OrigOrder ## T
#define sort_OrigOrder(T) _sort_OrigOrder(T)

#define _OrigOrder(T) OrigOrder ## T
#define OrigOrder(T) _OrigOrder(T)

#define _Unify(T) Unify ## T
#define Unify(T) _Unify(T)

#define _Init(T) Init ## T
#define Init(T) _Init(T)

#define _FreeAll(T) FreeAll ## T
#define FreeAll(T) _FreeAll(T)



/****************************************************************************/

#endif
