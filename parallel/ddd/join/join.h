// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      join.h                                                        */
/*                                                                          */
/* Purpose:   header file for join module                                   */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Institut fuer Computeranwendungen III                         */
/*            Universitaet Stuttgart                                        */
/*            Pfaffenwaldring 27                                            */
/*            70569 Stuttgart                                               */
/*            email: birken@ica3.uni-stuttgart.de                           */
/*            phone: 0049-(0)711-685-7007                                   */
/*            fax  : 0049-(0)711-685-7000                                   */
/*                                                                          */
/* History:   980126 kb  begin                                              */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/* RCS_ID
   $Header$
 */


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __JOIN_H__
#define __JOIN_H__


#define DebugJoin     10   /* 0 is all, 10 is off */


#include "basic/lowcomm.h"
#include "basic/oopp.h"    /* for object-orientated style via preprocessor */



/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/* some macros for customizing oopp */
#define _NEWPARAMS
#define _NEWPARAMS_OR_VOID    void

#define __INDENT(n)   { int i; for(i=0; i<n; i++) fputs("   ",fp);}
#define _PRINTPARAMS  , int indent, FILE *fp
#define _PRINTPARAMS_DEFAULT  ,0,stdout
#define _INDENT       __INDENT(indent)
#define _INDENT1      __INDENT(indent+1)
#define _INDENT2      __INDENT(indent+2)
#define _PRINTNEXT    , indent+1, fp
#define _PRINTSAME    , indent, fp

/* map memory allocation calls */
#define OO_Allocate  join_AllocTmp
#define OO_Free      join_FreeTmp


/* extra prefix for all join-related data structures and/or typedefs */
/*#define ClassPrefix*/


/* overall mode of transfer */
enum JoinMode {
  JMODE_IDLE = 0,                /* waiting for next DDD_JoinBegin() */
  JMODE_CMDS,                    /* after DDD_JoinBegin(), before DDD_JoinEnd() */
  JMODE_BUSY                     /* during DDD_JoinEnd() */
};



/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/* JIJoin: represents JoinObj command from application                      */
/****************************************************************************/

#define ClassName JIJoin
Class_Data_Begin
DDD_HDR hdr;                    /* local obj for which the join is requested    */
DDD_PROC dest;                  /* proc for joining                             */
DDD_GID new_gid;                /* gid of object on dest which should be joined */
Class_Data_End
void Method(Print)   (DefThis _PRINTPARAMS);
int  Method(Compare) (ClassPtr, ClassPtr);

#undef ClassName


/* define container class */
#define SetOf          JIJoin
#define Set_SegmSize   256
#define Set_BTreeOrder 32
#include "basic/ooppcc.h"



/****************************************************************************/
/* JIPartner: description of phase1 message on sender side                  */
/****************************************************************************/

typedef struct _JIPartner
{
  DDD_HDR hdr;                  /* local obj which has been contacted by join   */
  DDD_PROC proc;                /* proc which initiated join                    */
} JIPartner;



/****************************************************************************/
/* JIAddCpl: remote command to add cpl for join-obj during phase 2          */
/****************************************************************************/

typedef struct _TEAddCpl
{
  DDD_GID gid;                  /* gid of object to add coupling                */
  DDD_PROC proc;                /* proc of new coupling                         */
  DDD_PRIO prio;                /* prio of new coupling                         */

} TEAddCpl;


#define ClassName JIAddCpl
Class_Data_Begin
DDD_PROC dest;                  /* receiver of this item                        */
TEAddCpl te;                    /* table entry (for message)                    */
Class_Data_End
void Method(Print)   (DefThis _PRINTPARAMS);
int  Method(Compare) (ClassPtr, ClassPtr);

#undef ClassName


/* define container class */
#define SetOf          JIAddCpl
#define Set_SegmSize   256
#define Set_BTreeOrder 32
#include "basic/ooppcc.h"



/****************************************************************************/
/* JOINMSG1: description of phase1 message on sender side                   */
/****************************************************************************/

typedef struct _TEJoin
{
  DDD_GID gid;                  /* gid of distributed obj to join with          */
  DDD_PRIO prio;                /* prio of new local object which is joined     */

  DDD_HDR hdr;                  /* hdr of local DDD object (only for receiver)  */
} TEJoin;


typedef struct _JOINMSG1
{
  DDD_PROC dest;                 /* receiver of message                         */

  struct _JOINMSG1 *next;

  JIJoinPtr *arrayJoin;
  int nJoins;

  /* lowcomm message handle */
  LC_MSGHANDLE msg_h;

} JOINMSG1;



/****************************************************************************/
/* JOINMSG2: description of phase2 message on sender side                   */
/****************************************************************************/

typedef struct _JOINMSG2
{
  DDD_PROC dest;                 /* receiver of message                         */

  struct _JOINMSG2 *next;

  JIAddCplPtr *arrayAddCpl;
  int nAddCpls;

  /* lowcomm message handle */
  LC_MSGHANDLE msg_h;

} JOINMSG2;



/****************************************************************************/
/* JOINMSG3: description of phase3 message on sender side                   */
/****************************************************************************/

typedef struct _JOINMSG3
{
  DDD_PROC dest;                 /* receiver of message                         */

  struct _JOINMSG3 *next;

  JIAddCplPtr *arrayAddCpl;
  int nAddCpls;

  /* lowcomm message handle */
  LC_MSGHANDLE msg_h;

} JOINMSG3;




/****************************************************************************/
/* JOIN_GLOBALS: global data for join module                                */
/****************************************************************************/

typedef struct _JOIN_GLOBALS
{
  /* mode of join module */
  int joinMode;

  /* description for phase1 message */
  LC_MSGTYPE phase1msg_t;
  LC_MSGCOMP jointab_id;

  /* description for phase2 message */
  LC_MSGTYPE phase2msg_t;
  LC_MSGCOMP addtab_id;

  /* description for phase3 message */
  LC_MSGTYPE phase3msg_t;
  LC_MSGCOMP cpltab_id;


  /* entry points for global sets */
  JIJoinSet   *setJIJoin;
  JIAddCplSet *setJIAddCpl2;
  JIAddCplSet *setJIAddCpl3;

} JOIN_GLOBALS;


/* one instance of JOIN_GLOBALS */
extern JOIN_GLOBALS joinGlobals;


/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/


/* join.c, used only by cmds.c */
/*
   XICopyObj **CplClosureEstimate(XICopyObjPtrArray *, int *);
   int  PrepareObjMsgs(XICopyObjPtrArray *, XINewCpl **, int,
                XIOldCpl **, int, XFERMSG **, size_t *);
   void ExecLocalXIDelCmd(XIDelCmd  **, int);
   void ExecLocalXISetPrio(XISetPrioPtrArray *, XIDelObj  **,int, XICopyObj **,int);
   void ExecLocalXIDelObj(XIDelObj  **, int, XICopyObj **,int);
   void PropagateCplInfos(XISetPrio **, int, XIDelObj  **, int,
                TENewCpl *, int);
 */
void join_SetTmpMem (int);
int JoinMode (void);
int JoinStepMode(int);


/* pack.c,   used only by cmds.c */
/*
   void XferPackMsgs (XFERMSG *);
 */


/* unpack.c, used only by cmds.c */
/*
   void XferUnpack (LC_MSGHANDLE *, int, DDD_HDR *, int,
                XISetPrioPtrArray *, XIDelObj  **, int,
                XICopyObjPtrArray *, XICopyObj **, int);
 */


/* ctrl.c */
/*
   void XferDisplayMsg (char *comment, LC_MSGHANDLE);
 */



#endif
