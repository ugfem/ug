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

#define SEGM_SIZE  256


/****************************************************************************/

/*
        which structure components are needed for each item?
 */
#define SLL_INFO(T)   \
  int sll_n;                   /* unique index number */  \
  struct T *sll_next           /* linked list         */



/****************************************************************************/

/*
        macros for template support: datatypes
 */
#define _aSegm(T) aSegm ## T
#define aSegm(T) _aSegm(T)

#define _Segm(T) Segm ## T
#define Segm(T) _Segm(T)


/*
        macros for template support: variables
 */
#define _segms(T) segms ## T
#define segms(T) _segms(T)

#define _list(T) list ## T
#define list(T) _list(T)

#define _n(T) n ## T
#define n(T) _n(T)


/*
        macros for template support: functions
 */
#define _NewSegm(T) NewSegm ## T
#define NewSegm(T) _NewSegm(T)

#define _FreeSegms(T) FreeSegms ## T
#define FreeSegms(T) _FreeSegms(T)

#define _New(T) New ## T
#define New(T) _New(T)

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
