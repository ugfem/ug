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

/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __XFER_H__
#define __XFER_H__


#define DebugXfer  10  /* 0 is all, 10 is off */

#include "basic/lowcomm.h"
#include "sll.h"



/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/



/* types of 'newness' of incoming object (for xferunpack.c) */
enum XferNewType {
  NOTNEW     = 0x00,             /* object is not new                           */
  PARTNEW    = 0x01,             /* object has been updated partially           */
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

typedef struct _XICopyObj
{
  SLL_INFO(_XICopyObj);
  DDD_HDR hdr;                  /* local obj for which a copy should be created */
  DDD_GID gid;                  /* gid of local object                          */
  DDD_PROC dest;                /* proc involved with operation, me -> proc     */
  DDD_PRIO prio;                /* priority for new object copy                 */
  size_t size;                  /* V1.2: needed for variable-sized objects      */

  int new_owner;                /* did destination proc know this gid before?   */

  int addLen;
  struct   xferaddinfo *add;         /* additional data items                   */

} XICopyObj;

/* include template */
#define T XICopyObj
#include "sll.ht"



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



/****************************************************************************/
/* XISetPrio:                                                               */
/****************************************************************************/

typedef struct _XISetPrio
{
  SLL_INFO(_XISetPrio);
  DDD_HDR hdr;                  /* local obj for which prio should be set       */
  DDD_GID gid;                  /* gid of local object                          */
  DDD_PRIO prio;                /* new priority                                 */

  int is_valid;                 /* invalid iff there's a DelObj for same gid    */

} XISetPrio;


/* include template */
#define T XISetPrio
#include "sll.ht"




/****************************************************************************/
/* XINewCpl:                                                                */
/****************************************************************************/

typedef struct _TENewCpl
{
  DDD_GID gid;                  /* obj-gid for which new copy will be created   */
  DDD_PROC dest;                /* destination of new object copy               */
  DDD_PRIO prio;                /* priority of new object copy                  */

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


/****************************************************************************/
/* XFERMSG: complete description about message on sender side               */
/****************************************************************************/

typedef struct _XFERMSG
{
  DDD_PROC proc;                 /* receiver of message                         */
  size_t size;                   /* size of message data                        */

  struct _XFERMSG *next;


  XICopyObj **xferObjArray;
  int nObjItems;

  XINewCpl  **xferNewCpl;
  int nNewCpl;

  XIOldCpl  **xferOldCpl;
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


/* TODO was soll diese datenstruktur?? */

typedef struct _DELINFO
{
  DDD_GID gid;
  DDD_PROC proc;
  struct _DELINFO *next;
} DELINFO;




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
  /*
          LC_MSGCOMP  priotab_id, deltab_id;
   */

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


/* xfer.c, used only by cmds.c */
XICopyObj **CplClosureEstimate(XICopyObj **, int, int *);
int  PrepareObjMsgs(XICopyObj **, int, XINewCpl **, int,
                    XIOldCpl **, int, XFERMSG **, size_t *);
void ExecLocalXIDelCmd(XIDelCmd  **, int);
void ExecLocalXISetPrio(XISetPrio **,int, XIDelObj  **,int, XICopyObj **,int);
void ExecLocalXIDelObj(XIDelObj  **, int, XICopyObj **,int);
void PropagateCplInfos(XISetPrio **, int, XIDelObj  **, int,
                       TENewCpl *, int);
int XferMode (void);
int XferStepMode(int);


/* pack.c,   used only by cmds.c */
void XferPackMsgs (XFERMSG *);


/* unpack.c, used only by cmds.c */
void XferUnpack (LC_MSGHANDLE *, int, DDD_HDR *, int,
                 XISetPrio **, int, XIDelObj  **, int,
                 XICopyObj **, int, XICopyObj **, int);


/* ctrl.c */
void XferDisplayMsg (char *comment, LC_MSGHANDLE);



#endif
