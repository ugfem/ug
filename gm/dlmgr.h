// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      dlmgr.h                                                       */
/*                                                                          */
/* Purpose:   defines for dynamic linked list management                    */
/*                                                                          */
/* Author:    Stefan Lang                                                   */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70550 Stuttgart                                               */
/*            email: stefan@ica3.uni-stuttgart.de                           */
/*            phone: 0049-(0)711-685-7003                                   */
/*            fax  : 0049-(0)711-685-7000                                   */
/*                                                                          */
/* History:   960915 sl  start of dynamic list management					*/
/*                                                                          */
/* Remarks:                                                                 */
/*            Management of dynamic linked lists, which consist of          */
/*            several parts. An object can be mapped to a part of the list  */
/*            by its priority using PRIO2LISTPART().                        */
/*            The formulation on object basis allows for management of      */
/*            elements, nodes, vectors and vertices.                        */
/*            A list has the form:                                          */
/*                p0first-p0last->p1first-p1last->...->pnfirst..pnlast      */
/*            where each part p0,p1,..,pn is limited by two pointers,       */
/*            pxfirst and pxlast, x=0..n. The part numbers are ordered      */
/*            increasingly and connected in a manner that one can run       */
/*            through the whole list in increasing order (SUCC) but only    */
/*            through one listpart in decreasing order (PRED).              */
/*            Linking/Unlinking of objects in a list part is done in a      */
/*            way that preserves these conventions.                         */
/*                                                                          */
/****************************************************************************/

#ifndef __DLMGR_H__

#include "gm.h"

#ifdef ModelP
#include "parallel.h"
#endif

#ifdef ModelP

#define FIRSTPART_OF_LIST               0
#define LASTPART_OF_LIST(OTYPE) ((OTYPE ## PRIOS) -1)

#define GRID_UNLINK_OBJECT(Grid,Object,Prio,OTYPE,PRED,SUCC)\
  {\
    INT listpart = PRIO2LISTPART(OTYPE ## _LIST,Prio);\
    IFDEBUG(gm,1) \
    printf("%d: GRID_UNLINK_" # OTYPE "():" # Object \
           " has listpart=%d for prio=%d\n",me,listpart,Prio);\
    fflush(stdout);\
    ENDDEBUG \
\
    if (listpart<0 || listpart>LASTPART_OF_LIST(OTYPE)) {\
      printf("%d: GRID_UNLINK_" # OTYPE "(): ERROR " # Object \
             " has no valid listpart=%d for prio=%d\n",me,listpart,Prio);\
      fflush(stdout);\
      ASSERT(0);\
    }\
\
    switch (listpart) {\
\
    case FIRSTPART_OF_LIST :\
      if (PRED(Object)!=NULL)\
        SUCC(PRED(Object)) = SUCC(Object);\
      else\
        LISTPART_FIRST ## OTYPE(Grid,FIRSTPART_OF_LIST) = SUCC(Object);\
      if (SUCC(Object)!=LISTPART_FIRST ## OTYPE(Grid,FIRSTPART_OF_LIST+1))\
        PRED(SUCC(Object)) = PRED(Object);\
      else\
        LISTPART_LAST ## OTYPE(Grid,FIRSTPART_OF_LIST) = PRED(Object);\
      break;\
        \
    case LASTPART_OF_LIST(OTYPE) :\
      if (PRED(Object)!=NULL) \
        SUCC(PRED(Object)) = SUCC(Object);\
      else {\
        LISTPART_FIRST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)) = SUCC(Object);\
        if (LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)-1)!=NULL)\
          SUCC(LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)-1)) = SUCC(Object);\
      }\
      if (SUCC(Object)!=NULL) \
        PRED(SUCC(Object)) = PRED(Object);\
      else \
        LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)) = PRED(Object);\
      break;\
        \
    default :\
      /* unlink in middle of list */\
      if (PRED(Object)!=NULL) \
        SUCC(PRED(Object)) = SUCC(Object);\
      else {\
        INT listpartprev=listpart;\
\
        PRED(SUCC(Object)) = NULL;\
\
        do \
          listpartprev--;\
        while (listpartprev>FIRSTPART_OF_LIST && \
               LISTPART_LAST ## OTYPE(Grid,listpartprev)==NULL);\
\
        if (LISTPART_LAST ## OTYPE(Grid,listpartprev)!=NULL)\
          SUCC(LISTPART_LAST ## OTYPE(Grid,listpartprev)) = SUCC(Object);\
      }\
      if (LISTPART_LAST ## OTYPE(Grid,listpart) != Object) {\
        LISTPART_FIRST ## OTYPE(Grid,listpart) = SUCC(Object);\
        PRED(SUCC(Object)) = PRED(Object);\
      }\
      else {\
        LISTPART_FIRST ## OTYPE(Grid,listpart) = NULL;\
        LISTPART_LAST ## OTYPE(Grid,listpart) = NULL;\
      }\
\
      if (Object == LISTPART_LAST ## OTYPE(Grid,listpart))\
        LISTPART_LAST ## OTYPE(Grid,listpart) = PRED(Object);\
      break;\
    }\
  }

#else

#define GRID_UNLINK_OBJECT(Grid,Object,Prio,OTYPE,PRED,SUCC)\
  {\
    if (PRED(Object)!=NULL) \
      SUCC(PRED(Object)) = SUCC(Object);\
    else {\
      FIRST ## OTYPE(Grid) = SUCC(Object);\
      PRED(SUCC(Object)) = NULL;\
    }\
    if (SUCC(Object)!=NULL) \
      PRED(SUCC(Object)) = PRED(Object);\
    else {\
      LAST ## OTYPE(Grid) = PRED(Object);\
      SUCC(PRED(Object)) = NULL;\
    }\
  }
#endif


#ifdef ModelP
#define GRID_LINK_OBJECT(Grid,Object,Prio,OTYPE,PRED,SUCC)\
  {\
    INT listpart = PRIO2LISTPART(OTYPE ## _LIST,Prio);\
\
    ASSERT(Grid != NULL);\
    ASSERT(Object != NULL);\
    ASSERT(Prio >= 0);\
\
    IFDEBUG(gm,1) \
    printf("%d: GRID_LINK_" # OTYPE "():" # Object \
           " has listpart=%d for prio=%d\n",me,listpart,Prio);\
    fflush(stdout);\
    ENDDEBUG \
\
    if (listpart<0 || listpart>LASTPART_OF_LIST(OTYPE)) {\
      printf("%d: GRID_LINK_" # OTYPE "(): ERROR " # Object \
             " has no valid listpart=%d for prio=%d\n",me,listpart,Prio);\
      fflush(stdout);\
      ASSERT(0);\
    }\
\
    switch  (listpart) {\
      OTYPE *before;\
      OTYPE *after;\
\
    case FIRSTPART_OF_LIST :\
      before = LISTPART_FIRST ## OTYPE(Grid,FIRSTPART_OF_LIST);\
      LISTPART_FIRST ## OTYPE(Grid,FIRSTPART_OF_LIST) = Object;\
      if (LISTPART_LAST ## OTYPE(Grid,FIRSTPART_OF_LIST) == NULL)\
        LISTPART_LAST ## OTYPE(Grid,FIRSTPART_OF_LIST) = Object;\
      PRED(Object) = NULL;\
      if (before==NULL) {\
        SUCC(Object) = LISTPART_FIRST ## OTYPE(Grid,FIRSTPART_OF_LIST+1);\
      }\
      else {\
        SUCC(Object) = before;\
        PRED(before) = Object;\
      }\
      break;\
\
    case LASTPART_OF_LIST(OTYPE) : \
      after = LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE));\
      SUCC(Object) = NULL;\
      if (after==NULL)\
      {\
        PRED(Object) = NULL;\
        LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)) = Object;\
        LISTPART_FIRST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)) = Object;\
        if (LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)-1)!=NULL)\
          SUCC(LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)-1)) = Object;\
      }\
      else\
      {\
        PRED(Object) = after;\
        LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)) = Object;\
        SUCC(after) = Object;\
      }\
      break;\
\
    default : {\
      /* link in middle of list */\
      OTYPE *Object1 = LISTPART_FIRST ## OTYPE(Grid,listpart);\
      INT listpartprev = listpart;\
      INT listpartnext = listpart;\
\
      LISTPART_FIRST ## OTYPE(Grid,listpart) = Object;\
      SUCC(Object) = Object1;\
\
      /* empty list? */\
      if (Object1 == NULL) {\
        LISTPART_LAST ## OTYPE(Grid,listpart) = Object;\
        do\
          listpartnext++;\
        while (listpartnext<LASTPART_OF_LIST(OTYPE) &&\
               (Object1=LISTPART_FIRST ## OTYPE(Grid,listpartnext))==NULL);\
        if (Object1 != NULL)\
          SUCC(Object) = Object1;\
      }\
      else \
        PRED(Object1) == Object;\
\
      do\
        listpartprev--;\
      while (listpartprev>FIRSTPART_OF_LIST &&\
             (Object1=LISTPART_LAST ## OTYPE(Grid,listpartprev))==NULL);\
\
      if (Object1 != NULL)\
        SUCC(Object1) = Object;\
      break;\
    }\
    }\
  }
#else
#define GRID_LINK_OBJECT(Grid,Object,Prio,OTYPE,PRED,SUCC)   \
  {                                                        \
    OTYPE *after;                                        \
                                                             \
    after = LAST ## OTYPE(Grid);                         \
    SUCC(Object) = NULL;                                 \
    if (after==NULL) {                                   \
      PRED(Object) = NULL;                             \
      LAST ## OTYPE(Grid) = Object;                    \
      FIRST ## OTYPE(Grid) = Object;                   \
    }                                                    \
    else {                                               \
      PRED(Object) = after;                            \
      LAST ## OTYPE(Grid) = Object;                    \
      SUCC(after) = Object;                            \
    }                                                    \
  }
#endif

/********************************************************/
/*	macro definition for object managed in lists            */
/*  these definitions are based on templates from above	*/
/********************************************************/

/* element */
#define GRID_UNLINK_ELEMENT(Grid,Element,Prio) \
  GRID_UNLINK_OBJECT(Grid,Element,Prio,ELEMENT,PREDE,SUCCE)

#define GRID_LINK_ELEMENT(Grid,Element,Prio) \
  GRID_LINK_OBJECT(Grid,Element,Prio,ELEMENT,PREDE,SUCCE)

/* node */
#define GRID_UNLINK_NODE(Grid,Node,Prio) \
  GRID_UNLINK_OBJECT(Grid,Node,Prio,NODE,PREDN,SUCCN)

#define GRID_LINK_NODE(Grid,Node,Prio) \
  GRID_LINK_OBJECT(Grid,Node,Prio,NODE,PREDN,SUCCN)

/* vector */
#define GRID_UNLINK_VECTOR(Grid,Vector,Prio) \
  GRID_UNLINK_OBJECT(Grid,Vector,Prio,VECTOR,PREDVC,SUCCVC)

#define GRID_LINK_VECTOR(Grid,Vector,Prio) \
  GRID_LINK_OBJECT(Grid,Vector,Prio,VECTOR,PREDVC,SUCCVC)

#endif /* __DLMGR_H__ */
