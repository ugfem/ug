// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/*            DDD V1.7.3                                                    */
/*                                                                          */
/* File:      ddd.h                                                         */
/*                                                                          */
/* Purpose:   header file for ddd module                                    */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/*                                                                          */
/* History:   93/11/22 kb  begin                                            */
/*            94/06/27 kb  revision for usage in C++ context                */
/*            94/08/22 kb  major revision, corresp. to Ref.Man. V1.0        */
/*            94/09/13 kb  added DDD_Status                                 */
/*            95/01/13 kb  added range functionality                        */
/*            95/01/16 kb  added Statistics features                        */
/*            95/01/18 kb  moved Statistics to dddaddon.h                   */
/*            95/02/06 kb  added extended Ifs                               */
/*            95/03/08 kb  added UserData features                          */
/*            95/03/21 kb  added variable-sized Objects                     */
/*            95/05/22 kb  added var-sized AddData features                 */
/*            95/11/04 kb  complete redesign of ObjMgr and Registering      */
/*            95/11/06 kb  changed parameters of DDD_Init()                 */
/*            95/11/15 kb  ddd_hdr with arbitrary offset (started)          */
/*            96/01/08 kb  renamed range to attr                            */
/*            96/05/12 kb  xfer-unpack rewritten                            */
/*            96/06/05 kb  changed handler management                       */
/*            96/09/06 kb  xfer-module completely rewritten                 */
/*            96/11/28 kb  merged F_FRONTEND functionality from code branch */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __DDD__
#define __DDD__
#ifdef __cplusplus
extern "C" {
#endif


/****************************************************************************/
/*                                                                          */
/* settings for switching C_FRONTEND/F_FRONTEND and for                     */
/* creating portability across different F77-platforms.                     */
/*                                                                          */
/****************************************************************************/

/* default frontend is C_FRONTEND (at least as long as F_FRONTEND is being
   developed */
#ifndef C_FRONTEND
 #ifndef F_FRONTEND
  #define C_FRONTEND
 #endif
#endif

#define DDD_VERSION    "1.7.3"


/* F77SYM(lsym,usym) macro is defined in compiler.h. 961127 KB */


/* for size_t */
#include <stddef.h>



/****************************************************************************/
/*                                                                          */
/* compile time constants defining static data size (i.e. arrays)           */
/* other constants                                                          */
/*                                                                          */
/****************************************************************************/


/* types of elements for StructRegister */
/* (use negative values for combination with positive DDD_TYPEs) */
enum ElemType {
  EL_DDDHDR   =  0,                     /* element type: DDD header             */
  EL_GDATA    = -1,                     /* element type: global data            */
  EL_LDATA    = -2,                     /* element type: local data             */
  EL_DATAPTR  = -3,                     /* element type: data pointer           */
  EL_OBJPTR   = -4,                     /* element type: object pointer         */
  EL_CONTINUE = -5,                     /* continued element definition list    */
  EL_END      = -6                      /* end of element definition list       */
};



/* types of handlers for HandlerRegister */
enum HandlerType {
  HANDLER_LDATACONSTRUCTOR = 0,         /* handler: construct ldata             */
  HANDLER_UPDATE,                       /* handler: update object internals     */
  HANDLER_DESTRUCTOR,                   /* handler: remove obj from appl's data */
  HANDLER_OBJMKCONS,                    /* handler: make obj consistent         */
  HANDLER_XFERCOPY,                     /* handler: copy cmd during xfer        */
  HANDLER_XFERDELETE,                   /* handler: delete cmd during xfer      */
  HANDLER_XFERGATHER,                   /* handler: send additional data        */
  HANDLER_XFERSCATTER,                  /* handler: recv additional data        */
  HANDLER_XFERGATHERX,                  /* handler: send add. data, var-sized   */
  HANDLER_XFERSCATTERX,                 /* handler: recv add. data, var-sized   */
  HANDLER_COPYMANIP,                    /* handler: manipulate incoming copy    */
#ifdef F_FRONTEND
  HANDLER_ALLOCOBJ,                                     /* handler: get an OBJ_ID (Array index) */
  HANDLER_FREEOBJ,                                      /* handler: free an OBJ_ID				*/
#endif

  HANDLER_DELETE,                       /* handler: delete object and free mem  */
  HANDLER_SETPRIORITY,                  /* handler: change priority of object   */
  HANDLER_MAX,                          /* current number of handlers           */
  HANDLER_END = 999                     /* end of handler list                  */
};



/* options for DDD_SetOption */
enum OptionType {
  OPT_WARNING_VARSIZE_OBJ=0,       /* warning on differing obj sizes            */
  OPT_WARNING_SMALLSIZE,           /* warning on obj sizes smaller than declared*/
  OPT_WARNING_PRIOCHANGE,          /* warning on inconsistency in prio-change   */
  OPT_WARNING_DESTRUCT_HDR,        /* warning on inconsistency in HdrDestructor */
  OPT_DEBUG_XFERMESGS,             /* print debug info for xfer messages        */
  OPT_QUIET_CONSCHECK,             /* do ConsCheck in a quiet manner            */
  OPT_IDENTIFY_MODE,               /* one of the IDMODE_xxx constants           */
  OPT_WARNING_REF_COLLISION,       /* warning on collision in reference-localize*/
  OPT_END
};


enum SwitchType {
  OPT_OFF = 0,
  OPT_ON,

  IDMODE_LISTS,             /* ordering of each identify-tupel is relevant      */
  IDMODE_SETS               /* ordering of each identify-tupel is not sensitive */
};




/* direction of interface communication (DDD_IFOneway) */
enum IFDirection {
  IF_FORWARD  = 1,                     /* communicate from A to B               */
  IF_BACKWARD = 2                      /* communicate from B to A               */
};


/* ID of (predefined) standard interface */
enum IFConstants {
  STD_INTERFACE = 0
};


/* DDD_TYPE DDD_USER_DATA: send stream of bytes with XferAddData */
enum XferConstants {
  /* DDD_TYPE DDD_USER_DATA: send stream of bytes with XferAddData */
  DDD_USER_DATA = 0x8000,

  /* second parameter for MKCONS handler */
  XFER_UPGRADE = 0x9000,          /* object has been upgraded due to RULE C3 */
  XFER_NEW                        /* object is totally new */
};





/****************************************************************************/
/*                                                                          */
/* data structures and new types                                            */
/*                                                                          */
/****************************************************************************/


/*
        DDD object header, include this into all parallel object structures

        Some remarks:

           - don't touch the member elements of DDD_HEADER in the
             application program, they will be changed in further
             DDD versions!

           - use DDD functional interface for accessing the header fields;
             elements which are not accessible via the DDD functional interface
             should not be accessed by the application program anyway,
 */
typedef struct _DDD_HEADER
{
  /* control word elements */
  unsigned char typ;
  unsigned char prio;
  unsigned char attr;
  unsigned char flags;

  unsigned int myIndex;         /* global object array index */
  unsigned int gid;             /* global id */

  /*	char          empty[4];  *//* 4 unused bytes in current impl. */
} DDD_HEADER;



/*
        new DDD types, used during access of DDD functional interface
 */
typedef unsigned int DDD_TYPE;
typedef unsigned int DDD_IF;
typedef unsigned int DDD_GID;
typedef unsigned int DDD_PROC;
typedef unsigned int DDD_PRIO;
typedef unsigned int DDD_ATTR;
#ifdef C_FRONTEND
typedef char           * DDD_OBJ;
#else
typedef int DDD_OBJ;
#endif
typedef DDD_HEADER     * DDD_HDR;
typedef unsigned int DDD_OPTION;
typedef void (*HandlerPtr)();

#ifdef C_FRONTEND
typedef int (*ExecProcPtr)(DDD_OBJ);
typedef int (*ExecProcXPtr)(DDD_OBJ, DDD_PROC, DDD_PRIO);
typedef int (*ComProcPtr)(DDD_OBJ, void *);
typedef int (*ComProcXPtr)(DDD_OBJ, void *, DDD_PROC, DDD_PRIO);
#else
typedef int (*ExecProcPtr)(DDD_OBJ *);
typedef int (*ExecProcXPtr)(DDD_OBJ *, DDD_PROC *, DDD_PRIO *);
typedef int (*ComProcPtr)(DDD_OBJ *, void *);
typedef int (*ComProcXPtr)(DDD_OBJ *, void *, DDD_PROC *, DDD_PRIO *);
#endif



/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/


/*
        temporary macro in order to access header

        the application program should implement a
        DDD_HEADER access function itself.

        this will be removed in future versions.
 */
#define DDD_OBJ(o)     (&((o)->ddd))



/*
        external access of elements in DDD_HEADER

        C++ users should implement these access functions
        as inline functions.
 */
#ifndef __cplusplus
#ifdef C_FRONTEND
#define DDD_InfoPriority(ddd_hdr)    ((ddd_hdr)->prio)
#endif
#define DDD_InfoGlobalId(ddd_hdr)    ((ddd_hdr)->gid)
#define DDD_InfoAttr(ddd_hdr)       ((ddd_hdr)->attr)
#define DDD_InfoType(ddd_hdr)       ((ddd_hdr)->typ)
#endif


/****************************************************************************/
/*                                                                          */
/* declaration of DDD functional interface                                  */
/*                                                                          */
/****************************************************************************/


/*
        General DDD Module
 */
void     DDD_Init (int *argcp, char ***argvp);
void     DDD_Exit (void);
void     DDD_Status (void);
void     DDD_SetOption (DDD_OPTION, int);

#ifdef C_FRONTEND
DDD_PROC DDD_InfoMe (void);
DDD_PROC DDD_InfoMaster (void);
DDD_PROC DDD_InfoProcs (void);
#endif


/*
        Redirect line-oriented output, new in V1.2
 */
void     DDD_LineOutRegister (void (*func)(char *s));


/*
        Type Manager Module
 */
#ifdef C_FRONTEND
DDD_TYPE DDD_TypeDeclare (char *name);
void     DDD_TypeDefine (DDD_TYPE, ...);
void     DDD_TypeDisplay (DDD_TYPE);
void     DDD_HandlerRegister (DDD_TYPE, ...);
int      DDD_InfoTypes (void);
int      DDD_InfoHdrOffset (DDD_TYPE);
#else
#define  DDD_TypeDeclare     F77SYM(ddd_typedeclare,DDD_TYPEDECLARE)
#define  DDD_TypeDefine      F77SYM(ddd_typedefine,DDD_TYPEDEFINE)
#define  DDD_TypeDisplay     F77SYM(ddd_typedisplay,DDD_TYPEDISPLAY)
#define  DDD_HandlerRegister F77SYM(ddd_handlerregister,DDD_HANDLERREGISTER)
#define  DDD_InfoTypes       F77SYM(ddd_infotypes,DDD_INFOTYPES)

void     DDD_TypeDeclare (char *name, int *size, DDD_TYPE *type);
void     DDD_TypeDefine (DDD_TYPE *, ...);
void     DDD_TypeDisplay (DDD_TYPE *);
void     DDD_HandlerRegister (DDD_TYPE *, ...);
int      DDD_InfoTypes (void);
#endif



/*
        Object Properties
 */
void     DDD_PrioritySet (DDD_HDR, DDD_PRIO);
void     DDD_AttrSet (DDD_HDR, DDD_ATTR); /* this shouldn't be allowed */
int  *   DDD_InfoProcList (DDD_HDR);
DDD_PROC DDD_InfoProcPrio (DDD_HDR, DDD_PRIO);
int      DDD_InfoIsLocal (DDD_HDR);
int      DDD_InfoNCopies (DDD_HDR);

#ifdef F_FRONTEND
#define  DDD_InfoPriority    F77SYM(ddd_infopriority,DDD_INFOPRIORITY)
#endif


/*
        Identification Module
 */
#ifdef C_FRONTEND
void     DDD_IdentifyBegin (void);
void     DDD_IdentifyEnd (void);
void     DDD_IdentifyNumber (DDD_HDR, DDD_PROC, int);
void     DDD_IdentifyString (DDD_HDR, DDD_PROC, char *);
void     DDD_IdentifyObject (DDD_HDR, DDD_PROC, DDD_HDR);
#else
#define  DDD_IdentifyBegin  F77SYM(ddd_identifybegin,DDD_IDENTIFYBEGIN)
#define  DDD_IdentifyEnd    F77SYM(ddd_identifyend,DDD_IDENTIFYEND)
#define  DDD_IdentifyNumber F77SYM(ddd_identifynumber,DDD_IDENTIFYNUMBER)
#define  DDD_IdentifyString F77SYM(ddd_identifystring,DDD_IDENTIFYSTRING)
#define  DDD_IdentifyObject F77SYM(ddd_identifyobject,DDD_IDENTIFYOBJECT)

void     DDD_IdentifyBegin (void);
void     DDD_IdentifyEnd (void);
void     DDD_IdentifyNumber (DDD_TYPE *, DDD_OBJ *, DDD_PROC *, int *);
void     DDD_IdentifyString (DDD_TYPE *, DDD_OBJ *, DDD_PROC *, char *);
void     DDD_IdentifyObject (DDD_TYPE *, DDD_OBJ *, DDD_PROC *, DDD_TYPE *, DDD_OBJ *);
#endif


/*
        Interface Module
 */
#ifdef C_FRONTEND
DDD_IF   DDD_IFDefine (int, DDD_TYPE O[], int, DDD_PRIO A[], int, DDD_PRIO B[]);
void     DDD_IFDisplay (void);

void     DDD_IFExchange   (DDD_IF,             size_t, ComProcPtr,ComProcPtr);
void     DDD_IFOneway     (DDD_IF,         int,size_t, ComProcPtr,ComProcPtr);
void     DDD_IFExecLocal  (DDD_IF,                     ExecProcPtr);
void     DDD_IFAExchange  (DDD_IF,DDD_ATTR,    size_t, ComProcPtr,ComProcPtr);
void     DDD_IFAOneway    (DDD_IF,DDD_ATTR,int,size_t, ComProcPtr,ComProcPtr);
void     DDD_IFAExecLocal (DDD_IF,DDD_ATTR,            ExecProcPtr);
void     DDD_IFExchangeX  (DDD_IF,             size_t, ComProcXPtr,ComProcXPtr);
void     DDD_IFOnewayX    (DDD_IF,         int,size_t, ComProcXPtr,ComProcXPtr);
void     DDD_IFExecLocalX (DDD_IF,                     ExecProcXPtr);
void     DDD_IFAExchangeX (DDD_IF,DDD_ATTR,    size_t, ComProcXPtr,ComProcXPtr);
void     DDD_IFAOnewayX   (DDD_IF,DDD_ATTR,int,size_t, ComProcXPtr,ComProcXPtr);
void     DDD_IFAExecLocalX(DDD_IF,DDD_ATTR,            ExecProcXPtr);


#else

void     DDD_IFDefine (int *, DDD_TYPE *, int *, DDD_PRIO *, int *, DDD_PRIO *, DDD_IF *);
void     DDD_IFDisplay (void);

void     DDD_IFExchange   (DDD_IF *,                 size_t *, ComProcPtr,ComProcPtr);
void     DDD_IFOneway     (DDD_IF *,           int *,size_t *, ComProcPtr,ComProcPtr);
void     DDD_IFExecLocal  (DDD_IF *,                           ExecProcPtr);
void     DDD_IFAExchange  (DDD_IF *,DDD_ATTR *,      size_t *, ComProcPtr,ComProcPtr);
void     DDD_IFAOneway    (DDD_IF *,DDD_ATTR *,int *,size_t *, ComProcPtr,ComProcPtr);
void     DDD_IFAExecLocal (DDD_IF *,DDD_ATTR *,                ExecProcPtr);
void     DDD_IFExchangeX  (DDD_IF *,                 size_t *, ComProcXPtr,ComProcXPtr);
void     DDD_IFOnewayX    (DDD_IF *,           int *,size_t *, ComProcXPtr,ComProcXPtr);
void     DDD_IFExecLocalX (DDD_IF *,                           ExecProcXPtr);
void     DDD_IFAExchangeX (DDD_IF *,DDD_ATTR *,      size_t *, ComProcXPtr,ComProcXPtr);
void     DDD_IFAOnewayX   (DDD_IF *,DDD_ATTR *,int *,size_t *, ComProcXPtr,ComProcXPtr);
void     DDD_IFAExecLocalX(DDD_IF *,DDD_ATTR *,                ExecProcXPtr);

#define DDD_IFDefine     F77SYM(ddd_ifdefine,DDD_IFDEFINE)
#define DDD_IFDisplay    F77SYM(ddd_ifdisplay,DDD_IFDISPLAY)
#define DDD_IFExchange   F77SYM(ddd_ifexchange,DDD_IFEXCHANGE)
#define DDD_IFOneway     F77SYM(ddd_ifoneway,DDD_IFONEWAY)
#define DDD_IFExecLocal  F77SYM(ddd_ifexeclocal,DDD_IFEXECLOCAL)
#define DDD_IFAExchange  F77SYM(ddd_ifaexchange,DDD_IFAEXCHANGE)
#define DDD_IFAOneway    F77SYM(ddd_ifaoneway,DDD_IFAONEWAY)
#define DDD_IFAExecLocal F77SYM(ddd_ifaexeclocal,DDD_IFAEXECLOCAL)
#define DDD_IFExchangeX   F77SYM(ddd_ifexchangex,DDD_IFEXCHANGEX)
#define DDD_IFOnewayX     F77SYM(ddd_ifonewayx,DDD_IFONEWAYX)
#define DDD_IFExecLocalX  F77SYM(ddd_ifexeclocalx,DDD_IFEXECLOCALX)
#define DDD_IFAExchangeX  F77SYM(ddd_ifaexchangex,DDD_IFAEXCHANGEX)
#define DDD_IFAOnewayX    F77SYM(ddd_ifaonewayx,DDD_IFAONEWAYX)
#define DDD_IFAExecLocalX F77SYM(ddd_ifaexeclocalx,DDD_IFAEXECLOCALX)
#endif



/*
        Transfer Module
 */
#ifdef C_FRONTEND
void     DDD_XferBegin (void);
void     DDD_XferEnd (void);
void     DDD_XferCopyObj (DDD_HDR, DDD_PROC, DDD_PRIO);
void     DDD_XferCopyObjX (DDD_HDR, DDD_PROC, DDD_PRIO, size_t);
void     DDD_XferAddData (int cnt, DDD_TYPE);
void     DDD_XferAddDataX (int cnt, DDD_TYPE, size_t sizes[]);
void     DDD_XferDeleteObj (DDD_HDR);
#else
#define DDD_XferBegin     F77SYM(ddd_xferbegin,DDD_XFERBEGIN)
#define DDD_XferEnd       F77SYM(ddd_xferend,DDD_XFEREND)
#define DDD_XferCopyObj   F77SYM(ddd_xfercopyobj,DDD_XFERCOPYOBJ)
#define DDD_XferDeleteObj F77SYM(ddd_xferdeleteobj,DDD_XFERDELETEOBJ)

void     DDD_XferBegin (void);
void     DDD_XferEnd (void);
void     DDD_XferCopyObj (DDD_TYPE *, DDD_OBJ *, DDD_PROC *, DDD_PRIO *);
/* void     DDD_XferCopyObjX (DDD_HDR, DDD_PROC, DDD_PRIO, size_t); */
/* void     DDD_XferAddData (int cnt, DDD_TYPE); */
/* void     DDD_XferAddDataX (int cnt, DDD_TYPE, size_t sizes[]); */
void     DDD_XferDeleteObj (DDD_TYPE *, DDD_OBJ *);
#endif

/*
        Object Manager
 */
DDD_OBJ  DDD_ObjNew (size_t, DDD_TYPE, DDD_PRIO, DDD_ATTR);
void     DDD_ObjDelete (DDD_OBJ, size_t, DDD_TYPE);
void     DDD_HdrConstructor (DDD_HDR, DDD_TYPE, DDD_PRIO, DDD_ATTR);
void     DDD_HdrDestructor (DDD_HDR);
void     DDD_HdrConstructorMove (DDD_HDR, DDD_HDR);
#ifdef C_FRONTEND
DDD_OBJ  DDD_ObjGet (size_t, DDD_TYPE, DDD_PRIO, DDD_ATTR);
void     DDD_ObjUnGet (DDD_HDR, size_t);
#else
#define DDD_ObjGet   F77SYM(ddd_objget,DDD_OBJGET)
#define DDD_ObjUnGet F77SYM(ddd_objunget,DDD_OBJUNGET)

void     DDD_ObjGet (DDD_TYPE *, DDD_PRIO *, DDD_ATTR *, DDD_OBJ *);
void     DDD_ObjUnGet (DDD_OBJ *, DDD_TYPE *);
#endif


/*
        Special Fortran functions
 */
#ifdef F_FRONTEND
#define DDD_SetConfig F77SYM(ddd_setconfig,DDD_SETCONFIG)

void     DDD_SetConfig (int *, int *, int *);
#endif


/*
        Maintainance & Debugging
 */
#ifdef F_FRONTEND
#define DDD_ConsCheck        F77SYM(ddd_conscheck,DDD_CONSCHECK)
#define DDD_ListLocalObjects F77SYM(ddd_listlocalobjects,DDD_LISTLOCALOBJECTS)
#endif
int      DDD_ConsCheck (void);  /* returns total #errors since V1.6.6 */
void     DDD_ListLocalObjects (void);
DDD_HDR  DDD_SearchHdr (DDD_GID);


/****************************************************************************/


#ifdef __cplusplus
}
#endif
#endif
