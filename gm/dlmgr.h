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

#define HDRELEMENT      PARHDRE
#define HDRNODE         PARHDR
#define HDRVECTOR       PARHDR
#define HDR(OTYPE)      HDR ## OTYPE

#define GRID_UNLINK_OBJECT(Grid,Object,OTYPE,PRED,SUCC)\
  {\
    INT Prio = DDD_InfoPriority(HDR(OTYPE) ## (Object));\
    INT listpart = PRIO2LISTPART(OTYPE ## _LIST,Prio);\
    INT listpart1 = listpart;\
    OTYPE *Object1 = NULL;\
\
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
\
      if (Object != LISTPART_LAST ## OTYPE(Grid,FIRSTPART_OF_LIST))\
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
\
        do {\
          listpart1--;\
          Object1 = LISTPART_LAST ## OTYPE(Grid,listpart1);\
        }\
        while (listpart1>FIRSTPART_OF_LIST && Object1==NULL);\
\
        if (Object1!=NULL)\
          SUCC(Object1) = SUCC(Object);\
      }\
      if (SUCC(Object)!=NULL) \
        PRED(SUCC(Object)) = PRED(Object);\
      else {\
        LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)) = PRED(Object);\
        if (PRED(Object) != NULL) \
          SUCC(PRED(Object)) = NULL;\
      }\
      break;\
        \
    default :\
      /* unlink in middle of list */\
      if (PRED(Object)!=NULL) \
        SUCC(PRED(Object)) = SUCC(Object);\
      else {\
\
        PRED(SUCC(Object)) = NULL;\
\
        do {\
          listpart1--;\
          Object1 = LISTPART_LAST ## OTYPE(Grid,listpart1);\
        }\
        while (listpart1>FIRSTPART_OF_LIST && Object1==NULL);\
\
        if (Object1!=NULL)\
          SUCC(Object1) = SUCC(Object);\
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
    SUCC(Object) = PRED(Object) = NULL;\
  }

#else

#define GRID_UNLINK_OBJECT(Grid,Object,OTYPE,PRED,SUCC)\
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
    INT listpartprev = listpart;\
    INT listpartnext = listpart;\
    OTYPE *Object1 = NULL;\
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
\
    case FIRSTPART_OF_LIST :\
      Object1 = LISTPART_FIRST ## OTYPE(Grid,FIRSTPART_OF_LIST);\
      PRED(Object) = NULL;\
      LISTPART_FIRST ## OTYPE(Grid,FIRSTPART_OF_LIST) = Object;\
      if (Object1==NULL) {\
        LISTPART_LAST ## OTYPE(Grid,FIRSTPART_OF_LIST) = Object;\
        do {\
          listpartnext++;\
          Object1=LISTPART_FIRST ## OTYPE(Grid,listpartnext);\
        }\
        while (listpartnext<LASTPART_OF_LIST(OTYPE) && Object1==NULL);\
        SUCC(Object) = Object1;\
      }\
      else {\
        SUCC(Object) = Object1;\
        PRED(Object1) = Object;\
      }\
      break;\
\
    case LASTPART_OF_LIST(OTYPE) : \
      Object1 = LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE));\
      SUCC(Object) = NULL;\
      LISTPART_LAST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)) = Object;\
      if (Object1==NULL)\
      {\
        PRED(Object) = NULL;\
        LISTPART_FIRST ## OTYPE(Grid,LASTPART_OF_LIST(OTYPE)) = Object;\
\
        do {\
          listpartprev--;\
          Object1=LISTPART_LAST ## OTYPE(Grid,listpartprev);\
        }\
        while (listpartprev>FIRSTPART_OF_LIST && Object1==NULL);\
\
        if (Object1!=NULL)\
          SUCC(Object1) = Object;\
      }\
      else\
      {\
        PRED(Object) = Object1;\
        SUCC(Object1) = Object;\
      }\
      break;\
\
    default : \
      /* link in middle of list */\
      Object1 = LISTPART_FIRST ## OTYPE(Grid,listpart);\
\
      LISTPART_FIRST ## OTYPE(Grid,listpart) = Object;\
      SUCC(Object) = Object1;\
      PRED(Object) = NULL;\
\
      /* empty list? */\
      if (Object1 == NULL) {\
        LISTPART_LAST ## OTYPE(Grid,listpart) = Object;\
        do {\
          listpartnext++;\
          Object1=LISTPART_FIRST ## OTYPE(Grid,listpartnext);\
        }\
        while (listpartnext<LASTPART_OF_LIST(OTYPE) && Object1==NULL);\
        SUCC(Object) = Object1;\
      }\
      else \
        PRED(Object1) == Object;\
\
      do {\
        listpartprev--;\
        Object1=LISTPART_LAST ## OTYPE(Grid,listpartprev);\
      }\
      while (listpartprev>FIRSTPART_OF_LIST && Object1==NULL);\
\
      if (Object1 != NULL)\
        SUCC(Object1) = Object;\
      break;\
    }\
  }
#else
#define GRID_LINK_OBJECT(Grid,Object,Prio,OTYPE,PRED,SUCC)\
  {\
    OTYPE *after;\
\
    after = LAST ## OTYPE(Grid);\
    SUCC(Object) = NULL;\
    if (after==NULL) {\
      PRED(Object) = NULL;\
      LAST ## OTYPE(Grid) = Object;\
      FIRST ## OTYPE(Grid) = Object;\
    }\
    else {\
      PRED(Object) = after;\
      LAST ## OTYPE(Grid) = Object;\
      SUCC(after) = Object;\
    }\
  }
#endif

#ifdef ModelP
#define GRID_INIT_OBJECT_LIST(Grid,OTYPE) \
  {\
    INT i;\
    for (i=0; i<=LASTPART_OF_LIST(OTYPE); i++) {\
      LISTPART_FIRST  ## OTYPE(Grid,i) = NULL;\
      LISTPART_LAST ## OTYPE(Grid,i) = NULL;\
    }\
  }
#else
#define GRID_INIT_OBJECT_LIST(Grid,OTYPE) \
  FIRST ## OTYPE(Grid) = LAST ## OTYPE(Grid) = NULL;
#endif

#ifdef ModelP
/* TODO: define this */
#define GRID_CHECK_OBJECT_LIST(Grid,OTYPE)
#endif

/********************************************************/
/*	macro definition for object managed in lists            */
/*  these definitions are based on templates from above	*/
/********************************************************/

/* element */
#define GRID_UNLINK_ELEMENT(Grid,Element) \
  GRID_UNLINK_OBJECT(Grid,Element,ELEMENT,PREDE,SUCCE)

#define GRID_LINK_ELEMENT(Grid,Element,Prio) \
  GRID_LINK_OBJECT(Grid,Element,Prio,ELEMENT,PREDE,SUCCE)

#define GRID_INIT_ELEMENT_LIST(Grid)\
  GRID_INIT_OBJECT_LIST(Grid,ELEMENT)

/* node */
#define GRID_UNLINK_NODE(Grid,Node) \
  GRID_UNLINK_OBJECT(Grid,Node,NODE,PREDN,SUCCN)

#define GRID_LINK_NODE(Grid,Node,Prio) \
  GRID_LINK_OBJECT(Grid,Node,Prio,NODE,PREDN,SUCCN)

#define GRID_INIT_NODE_LIST(Grid)\
  GRID_INIT_OBJECT_LIST(Grid,NODE)

/* vector */
#define GRID_UNLINK_VECTOR(Grid,Vector) \
  GRID_UNLINK_OBJECT(Grid,Vector,VECTOR,PREDVC,SUCCVC)

#define GRID_LINK_VECTOR(Grid,Vector,Prio) \
  GRID_LINK_OBJECT(Grid,Vector,Prio,VECTOR,PREDVC,SUCCVC)

#define GRID_INIT_VECTOR_LIST(Grid)\
  GRID_INIT_OBJECT_LIST(Grid,VECTOR)

#endif /* __DLMGR_H__ */
