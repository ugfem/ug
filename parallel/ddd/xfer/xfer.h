// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      xfer.h                                                        */
/*                                                                          */
/* Purpose:   header file for xfer module                                   */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/*                                                                          */
/* History:   931122 kb  begin                                              */
/*            950404 kb  copied from dddi.h                                 */
/*            960718 kb  introduced lowcomm-layer (sets of messages)        */
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

#ifndef __XFER_H__
#define __XFER_H__


#define DebugXfer     10   /* 0 is all, 10 is off */
#define DebugPack     6    /* off: 6 */
#define DebugUnpack   5    /* off: 5 */
#define DebugCmdMsg   10   /* off: 10 */
#define DebugCplMsg   10   /* off: 10 */




#include "basic/lowcomm.h"
#include "basic/oopp.h"    /* for object-orientated style via preprocessor */
#include "sll.h"           /* TODO: remove this in future versions */



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
#define OO_Allocate  AllocTmp
#define OO_Free      FreeTmp


/* extra prefix for all xfer-related data structures and/or typedefs */
/*#define ClassPrefix*/


/* types of 'newness' of incoming object (for xferunpack.c) */
enum XferNewType {
  NOTNEW     = 0x00,             /* object is not new                           */
  PARTNEW    = 0x01,             /* object has been updated partially           */
  PRUNEDNEW  = 0x02,             /* object is new due to PruneDel (see cmdmsg.c)*/
  TOTALNEW   = 0x04,             /* object has been updated completely          */
  OTHERMSG   = 0x08,             /* object is taken from another message        */
  THISMSG    = 0x10              /* object is taken from this msg, temp setting */
};


/* overall mode of transfer */
enum XferMode {
  XMODE_IDLE = 0,                /* waiting for next DDD_XferBegin() */
  XMODE_CMDS,                    /* after DDD_XferBegin(), before DDD_XferEnd() */
  XMODE_BUSY                     /* during DDD_XferEnd() */
};



/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/* XFERADDDATA: description of additional data on sender side               */
/****************************************************************************/

typedef struct xferaddinfo {
  int addCnt;
  DDD_TYPE addTyp;
  int addLen;                    /* length of additional data                   */
  int addNPointers;
  int      *sizes;

  struct   xferaddinfo *next;
} XFERADDDATA;


/****************************************************************************/
/* XICopyObj:                                                               */
/****************************************************************************/

#define ClassName XICopyObj
Class_Data_Begin
DDD_HDR hdr;                    /* local obj for which a copy should be created */
DDD_GID gid;                    /* gid of local object                          */
DDD_PROC dest;                  /* proc involved with operation, me -> proc     */
DDD_PRIO prio;                  /* priority for new object copy                 */
size_t size;                    /* V1.2: needed for variable-sized objects      */

int addLen;
XFERADDDATA *add;               /* additional data items                   */

int flags;
Class_Data_End
void Method(Print)   (DefThis _PRINTPARAMS);
int  Method(Compare) (ClassPtr, ClassPtr);

#undef ClassName


/* define container class */
#define SetOf          XICopyObj
#define Set_SegmSize   256
#define Set_BTreeOrder 32
#include "basic/ooppcc.h"


/* usage of flags in XICopyObj */
/* usage of 0x01 while PruneXIDelCpl, temporarily */
#define MASK_CO_SELF 0x00000001
#define CO_SELF(c) (((int)((c)->flags))&MASK_CO_SELF)
#define SET_CO_SELF(c,d)     ((c)->flags) =   \
  ((((c)->flags)&(~MASK_CO_SELF)) | ((d)&MASK_CO_SELF))

/* usage of 0x02 for new_owner:
        did destination proc know the object's gid before? */
#define MASK_CO_NEWOWNER 0x00000002
#define CO_NEWOWNER(c) (((int)((c)->flags))&MASK_CO_NEWOWNER)
#define SET_CO_NEWOWNER(c)     ((c)->flags) = (((c)->flags) | MASK_CO_NEWOWNER)
#define CLEAR_CO_NEWOWNER(c)   ((c)->flags) = (((c)->flags)&(~MASK_CO_NEWOWNER))



/****************************************************************************/
/* XIDelObj:                                                                */
/****************************************************************************/

/* XIDelCmd represents a XferDeleteObj-cmd by the application program */
typedef struct _XIDelCmd
{
  SLL_INFO(_XIDelCmd);
  DDD_HDR hdr;

        #ifdef CPP_FRONTEND
  DDD_Object *obj;
        #endif

} XIDelCmd;

/* include template */
#define T XIDelCmd
#include "sll.ht"
#undef T


/* XIDelObj represents an object-delete-action (cf. XferRegisterDelete()) */
typedef struct _XIDelObj
{
  SLL_INFO(_XIDelObj);
  DDD_GID gid;                  /* gid of local object                          */

  /* hdr is explicitly not stored here, because object may be deleted
     when this XIDelObj-item is evaluated. */

  struct _XIDelCpl *delcpls;        /* couplings of deleted object              */

} XIDelObj;

/* include template */
#define T XIDelObj
#include "sll.ht"
#undef T



/****************************************************************************/
/* XISetPrio:                                                               */
/****************************************************************************/

#define ClassName XISetPrio
Class_Data_Begin
DDD_HDR hdr;                    /* local obj for which prio should be set       */
DDD_GID gid;                    /* gid of local object                          */
DDD_PRIO prio;                  /* new priority                                 */

int is_valid;                   /* invalid iff there's a DelObj for same gid    */
Class_Data_End
void Method(Print)   (DefThis _PRINTPARAMS);
int  Method(Compare) (ClassPtr, ClassPtr);

#undef ClassName


/* define container class */
#define SetOf          XISetPrio
#define Set_SegmSize   256
#define Set_BTreeOrder 32
#include "basic/ooppcc.h"


/****************************************************************************/
/* XINewCpl:                                                                */
/****************************************************************************/

typedef struct _TENewCpl
{
  DDD_GID gid;                  /* obj-gid for which new copy will be created   */
  DDD_PROC dest;                /* destination of new object copy               */
  DDD_PRIO prio;                /* priority of new object copy                  */
  DDD_TYPE type;                /* ddd-type of gid for PriorityMerge on receiver*/

} TENewCpl;


typedef struct _XINewCpl
{
  SLL_INFO(_XINewCpl);
  DDD_PROC to;                  /* receiver of this item                        */
  TENewCpl te;                  /* table entry (for message)                    */

} XINewCpl;

/* include template */
#define T XINewCpl
#include "sll.ht"
#undef T



/****************************************************************************/
/* XIOldCpl:                                                                */
/****************************************************************************/

typedef struct _TEOldCpl
{
  DDD_GID gid;                  /* obj-gid of local copy                        */
  DDD_PROC proc;                /* owner of that local object                   */
  DDD_PRIO prio;                /* priority of that local object                */

} TEOldCpl;


typedef struct _XIOldCpl
{
  SLL_INFO(_XIOldCpl);
  DDD_PROC to;                  /* receiver of this item                        */
  TEOldCpl te;                  /* table entry (for message)                    */

} XIOldCpl;

/* include template */
#define T XIOldCpl
#include "sll.ht"
#undef T


/****************************************************************************/
/* XIAddCpl:                                                                */
/****************************************************************************/

typedef struct _TEAddCpl
{
  DDD_GID gid;                  /* obj-gid of new local object                  */
  DDD_PROC proc;                /* owner of new object copy                     */
  DDD_PRIO prio;                /* priority of new local object                 */

} TEAddCpl;


typedef struct _XIAddCpl
{
  SLL_INFO(_XIAddCpl);
  DDD_PROC to;                  /* receiver of this item                        */
  TEAddCpl te;                  /* table entry (for message)                    */

} XIAddCpl;

/* include template */
#define T XIAddCpl
#include "sll.ht"
#undef T



/****************************************************************************/
/* XIDelCpl:                                                                */
/****************************************************************************/

typedef struct _TEDelCpl
{
  DDD_GID gid;                  /* obj-gid of deleted local object              */

} TEDelCpl;


typedef struct _XIDelCpl
{
  SLL_INFO(_XIDelCpl);
  DDD_PROC to;                  /* receiver of this item                        */
  TEDelCpl te;                  /* table entry (for message)                    */

  DDD_PRIO prio;                /* prio of deleted coupling                     */
  struct _XIDelCpl *next;       /* linked list of XIDelCpls                     */

} XIDelCpl;

/* include template */
#define T XIDelCpl
#include "sll.ht"
#undef T




/****************************************************************************/
/* XIModCpl:                                                                */
/****************************************************************************/

typedef struct _TEModCpl
{
  DDD_GID gid;                  /* obj-gid of corresponding object              */
  DDD_PRIO prio;                /* new priority for this obj on sending proc    */

} TEModCpl;


typedef struct _XIModCpl
{
  SLL_INFO(_XIModCpl);
  DDD_PROC to;                  /* receiver of this item                        */
  TEModCpl te;                  /* table entry (for message)                    */

  DDD_TYPE typ;                 /* type of corresponding object                 */

} XIModCpl;

/* include template */
#define T XIModCpl
#include "sll.ht"
#undef T


/****************************************************************************/
/* XFERMSG: complete description about message on sender side               */
/****************************************************************************/

typedef struct _XFERMSG
{
  DDD_PROC proc;                 /* receiver of message                         */
  size_t size;                   /* size of message data                        */

  struct _XFERMSG *next;


  XICopyObjPtr *xferObjArray;
  int nObjItems;

  XINewCpl     **xferNewCpl;
  int nNewCpl;

  XIOldCpl     **xferOldCpl;
  int nOldCpl;

  int nPointers;
  int nObjects;


  /* lowcomm message handle */
  LC_MSGHANDLE msg_h;

} XFERMSG;




/****************************************************************************/
/* SYMTAB_ENTRY: single entry of symbol table inside message                */
/****************************************************************************/


typedef struct
{
  DDD_GID gid;

  union {
    DDD_HDR hdr;                      /* used on receiver side only */
    DDD_OBJ    *ref;                  /* used on sender side only   */
  } adr;                          /* undefined during transfer  */
} SYMTAB_ENTRY;



/****************************************************************************/
/* OBJTAB_ENTRY: single entry of object table inside message                */
/****************************************************************************/


/* TODO: this should only store the offset of the DDD_HDR inside the
   message. gid, typ, prio(?), attr can be derived from that and are
   stored twice at present. */

typedef struct
{
  DDD_GID gid;
  int offset;                     /* offset from beginObjMem */
  DDD_TYPE typ;                   /* type of object */
  DDD_PRIO prio;                  /* new priority of object */
  DDD_ATTR attr;                  /* attr of object */
  int addLen;                     /* length of additional data */
  size_t size;                    /* size of object, ==desc->len for
                                      fixed-sized objects */

  DDD_HDR hdr;                /* TODO this is probably not used on sender side */

  /* TODO: the following data is only used on receiver side */
  int is_new;
  DDD_PRIO oldprio;
} OBJTAB_ENTRY;



/****************************************************************************/
/* DELTAB_ENTRY: single entry of deletion table inside message              */
/****************************************************************************/

typedef struct {
  DDD_GID gid;
  DDD_PROC proc;
} DELTAB_ENTRY;



/****************************************************************************/
/* XFER_GLOBALS: global data for xfer module                                */
/****************************************************************************/

typedef struct _XFER_GLOBALS
{
  /* mode of xfer module */
  int xferMode;

  /* description for object message */
  LC_MSGTYPE objmsg_t;
  LC_MSGCOMP symtab_id, objtab_id;
  LC_MSGCOMP newcpl_id, oldcpl_id;
  LC_MSGCOMP objmem_id;


  /* entry points for global sets */
  XICopyObjSet *setXICopyObj;
  XISetPrioSet *setXISetPrio;

} XFER_GLOBALS;


/* one instance of XFER_GLOBALS */
extern XFER_GLOBALS xferGlobals;


/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/


/* supp.c */
XFERADDDATA *NewXIAddData (void);
void FreeAllXIAddData (void);
int *AddDataAllocSizes(int);
/* and others, via template mechanism */



/* cplmsg.c */
void CommunicateCplMsgs (XIDelCpl **, int,
                         XIModCpl **, int, XIAddCpl **, int, DDD_HDR *, int);
void CplMsgInit (void);
void CplMsgExit (void);


/* cmdmsg.c */
int  PruneXIDelCmd (XIDelCmd **, int, XICopyObjPtrArray *);
void CmdMsgInit (void);
void CmdMsgExit (void);


/* xfer.c, used only by cmds.c */
XICopyObj **CplClosureEstimate(XICopyObjPtrArray *, int *);
int  PrepareObjMsgs(XICopyObjPtrArray *, XINewCpl **, int,
                    XIOldCpl **, int, XFERMSG **, size_t *);
void ExecLocalXIDelCmd(XIDelCmd  **, int);
void ExecLocalXISetPrio(XISetPrioPtrArray *, XIDelObj  **,int, XICopyObj **,int);
void ExecLocalXIDelObj(XIDelObj  **, int, XICopyObj **,int);
void PropagateCplInfos(XISetPrio **, int, XIDelObj  **, int,
                       TENewCpl *, int);
int XferMode (void);
int XferStepMode(int);


/* pack.c,   used only by cmds.c */
void XferPackMsgs (XFERMSG *);


/* unpack.c, used only by cmds.c */
void XferUnpack (LC_MSGHANDLE *, int, DDD_HDR *, int,
                 XISetPrioPtrArray *, XIDelObj  **, int,
                 XICopyObjPtrArray *, XICopyObj **, int);


/* ctrl.c */
void XferDisplayMsg (char *comment, LC_MSGHANDLE);



#endif
