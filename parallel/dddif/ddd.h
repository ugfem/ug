// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/*            DDD V1.6                                                      */
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
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __DDD__
#define __DDD__
#ifdef __cplusplus
extern "C" {
#endif


/* temporary setting, set default = C_FRONTEND */
/* (as long as F_FRONTEND is on its way)       */
#ifndef C_FRONTEND
#define C_FRONTEND
#endif



#define DDD_VERSION    "1.6"


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
  HANDLER_DESTRUCTOR,                   /* handler: destruct object             */
  HANDLER_OBJMKCONS,                    /* handler: make obj consistent         */
  HANDLER_XFERCOPY,                     /* handler: copy cmd during xfer        */
  HANDLER_XFERDELETE,                   /* handler: delete cmd during xfer      */
  HANDLER_XFERGATHER,                   /* handler: send additional data        */
  HANDLER_XFERSCATTER,                  /* handler: recv additional data        */
  HANDLER_XFERGATHERX,                  /* handler: send add. data, var-sized   */
  HANDLER_XFERSCATTERX,                 /* handler: recv add. data, var-sized   */
  HANDLER_COPYMANIP,                    /* handler: manipulate incoming copy    */
  HANDLER_END                           /* end of handler list                  */
};



/* options for DDD_SetOption */
enum OptionType {
  OPT_WARNING_VARSIZE_OBJ = 0,        /* option: warning on differing obj sizes */
  OPT_END
};


enum SwitchType {
  OPT_OFF = 0,
  OPT_ON  = 1
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


/* DDD_TYPE USER_DATA: send stream of bytes with XferAddData */
enum XferConstants {
  DDD_USER_DATA = -1
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

  char empty[4];                /* 4 unused bytes in current impl. */
} DDD_HEADER;



/*
        new DDD types, used during access of DDD functional interface
 */
typedef int DDD_TYPE;
typedef int DDD_IF;
typedef unsigned int DDD_GID;
typedef int DDD_PROC;
typedef int DDD_PRIO;
typedef int DDD_ATTR;
#ifdef C_FRONTEND
typedef char           * DDD_OBJ;
#else
typedef int DDD_OBJ;
#endif
typedef DDD_HEADER     * DDD_HDR;
typedef int DDD_OPTION;
typedef void (*HandlerPtr)();
typedef int (*ComProcPtr)(DDD_OBJ, void *);
typedef int (*ComProcXPtr)(DDD_OBJ, void *, DDD_PROC, DDD_PRIO);



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
#define DDD_InfoPriority(ddd_hdr)    ((ddd_hdr)->prio)
#define DDD_InfoGlobalId(ddd_hdr)    ((ddd_hdr)->gid)
#define DDD_InfoAttr(ddd_hdr)       ((ddd_hdr)->attr)
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
#ifdef __sgi
#define  DDD_TypeDeclare ddd_typedeclare_
#define  DDD_TypeDefine ddd_typedefine_
#define  DDD_TypeDisplay ddd_typedisplay_
#define  DDD_HandlerRegister ddd_handlerregister_
#define  DDD_InfoTypes ddd_infotypes_
#endif

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
/* void     DDD_AttrSet (DDD_HDR, DDD_ATTR); not allowed */
int  *   DDD_InfoProcList (DDD_HDR);



/*
        Identification Module
 */
void     DDD_IdentifyBegin (void);
void     DDD_IdentifyEnd (void);
void     DDD_IdentifyNumber (DDD_HDR, DDD_PROC, int);
void     DDD_IdentifyString (DDD_HDR, DDD_PROC, char *);
void     DDD_IdentifyObject (DDD_HDR, DDD_PROC, DDD_HDR);



/*
        Interface Module
 */
DDD_IF   DDD_IFDefine (int nO, int O[], int nA, int A[], int nB, int B[]);
void     DDD_IFDisplay (void);
void     DDD_IFExchange (DDD_IF, size_t, ComProcPtr GatherP, ComProcPtr ScatterP);
void     DDD_IFOneway (DDD_IF, int dir,
                       size_t, ComProcPtr GatherP, ComProcPtr ScatterP);

/*
        Attr-oriented communication, new in V1.1
 */
void     DDD_IFRExchange (DDD_IF, DDD_ATTR,
                          size_t, ComProcPtr GatherP, ComProcPtr ScatterP);
void     DDD_IFROneway (DDD_IF, DDD_ATTR, int dir,
                        size_t, ComProcPtr GatherP, ComProcPtr ScatterP);

/*
        Extended communication, new in V1.1
 */
void     DDD_IFExchangeX (DDD_IF,
                          size_t, ComProcXPtr GatherP, ComProcXPtr ScatterP);
void     DDD_IFOnewayX (DDD_IF, int dir,
                        size_t, ComProcXPtr GatherP, ComProcXPtr ScatterP);



/*
        Transfer Module
 */
void     DDD_XferBegin (void);
void     DDD_XferEnd (void);
void     DDD_XferCopyObj (DDD_HDR, DDD_PROC, DDD_PRIO);
void     DDD_XferCopyObjX (DDD_HDR, DDD_PROC, DDD_PRIO, size_t);
void     DDD_XferAddData (int cnt, DDD_TYPE);
void     DDD_XferAddDataX (int cnt, DDD_TYPE, size_t sizes[]);
void     DDD_XferDeleteObj (DDD_HDR);


/*
        Object Manager
 */
DDD_OBJ  DDD_ObjNew (size_t, DDD_TYPE, DDD_PRIO, DDD_ATTR);
void     DDD_ObjDelete (void *, size_t, DDD_TYPE);
void     DDD_HdrConstructor (DDD_HDR, DDD_TYPE, DDD_PRIO, DDD_ATTR);
void     DDD_HdrDestructor (DDD_HDR);
void     DDD_HdrConstructorMove (DDD_HDR, DDD_HDR);
DDD_OBJ  DDD_ObjGet (size_t, DDD_TYPE, DDD_PRIO, DDD_ATTR);
void     DDD_ObjUnGet (DDD_HDR, size_t);



/*
        Maintainance & Debugging
 */
void     DDD_ConsCheck (void);
void     DDD_ListLocalObjects (void);



#ifdef __cplusplus
}
#endif
#endif
