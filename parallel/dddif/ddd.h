// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/*            DDD V1.3                                                      */
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
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __DDD__
#define __DDD__
#ifdef __cplusplus
extern "C" {
#endif


#define DDD_VERSION    "1.3beta"


/****************************************************************************/
/*                                                                          */
/* compile time constants defining static data size (i.e. arrays)           */
/* other constants                                                          */
/*                                                                          */
/****************************************************************************/


/* types of elements for StructRegister */
enum ElemType {
  EL_GDATA    = 1,                      /* element type: global data            */
  EL_LDATA    = 2,                      /* element type: local data             */
  EL_DATAPTR  = 3,                      /* element type: data pointer           */
  EL_OBJPTR   = 4,                      /* element type: object pointer         */
  EL_END      = 0                       /* end of element definition list       */
};



/* types of handlers for HandlerRegister */
enum HandlerType {
  HANDLER_OBJINIT = 0,                  /* handler: init object                 */
  HANDLER_OBJUPDATE,                    /* handler: update objects internals    */
  HANDLER_OBJMKCONS,                    /* handler: make obj consistent         */
  HANDLER_OBJDELETE,                    /* handler: delete object               */
  HANDLER_XFERCOPY,                     /* handler: copy cmd during xfer        */
  HANDLER_XFERDELETE,                   /* handler: delete cmd during xfer      */
  HANDLER_XFERGATHER,                   /* handler: send additional data        */
  HANDLER_XFERSCATTER,                  /* handler: recv additional data        */
  HANDLER_XFERGATHERX,                  /* handler: send add. data, var-sized   */
  HANDLER_XFERSCATTERX,                 /* handler: recv add. data, var-sized   */
  HANDLER_COPYMANIP,                    /* handler: manipulate incoming copy    */
  HANDLER_LDATARESTORE,                 /* handler: restore LDATA after copying */
  HANDLER_END                           /* end of handler list                  */
};



/* direction of interface communication (DDD_IFOneway) */
enum IFDirection {
  IF_FORWARD  = 1,                     /* communicate from A to B               */
  IF_BACKWARD = 2                      /* communicate from B to A               */
};



/* no of (predefined) standard interface */
#define STD_INTERFACE     0


/* DDD_TYPE USER_DATA: send stream of bytes with XferAddData */
#define USER_DATA        -1



/****************************************************************************/
/*                                                                          */
/* data structures and new types                                            */
/*                                                                          */
/****************************************************************************/


/*
        DDD object header, include this into all parallel object structures

        Some remarks:

           - include macro DDD_HEADER as first line in all structure
             definitions in the application program

           - don't touch the member elements of OBJ_HEADER in the
             application program, they will be changed in further
             DDD versions!

           - use DDD functional interface for accessing the header fields;
             elements which are not accessible via the DDD functional interface
             should not be accessed by the application program anyway,
 */
typedef struct obj_header {
  /* control word elements */
  unsigned char typ;
  unsigned char prio;
  unsigned char range;
  unsigned char flags;

  unsigned int myIndex;         /* global object array index */
  unsigned int gid;             /* global id */

  char empty[4];                /* 4 unused bytes in current impl. */
} OBJ_HEADER;



/*
        new DDD types, used during access of DDD functional interface
 */
typedef int DDD_TYPE;
typedef int INTERFACE;
typedef int PROC_ID;
typedef OBJ_HEADER       *OBJECT;
typedef OBJ_HEADER       *HEADER;
typedef void (*HandlerPtr)();
typedef int (*ComProcPtr)(OBJECT theObject, void *data);
typedef int (*ComProcXPtr)(OBJECT theObject, void *data, int proc, int prio);



/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/


/* referencing header from object */
#define DDD_OBJ(o)    ((OBJ_HEADER *) o)


/* DDD_HEADER for declaration of shared data structs */
#define DDD_HEADER    OBJ_HEADER ddd_hdr


/* external usage of fields in ddd_hdr */
#define DDD_InfoPriority(o)    ((o)->ddd_hdr.prio)
#define DDD_InfoGlobalId(o)    ((o)->ddd_hdr.gid)
#define DDD_InfoRange(o)       ((o)->ddd_hdr.range)



/****************************************************************************/
/*                                                                          */
/* declaration of DDD functional interface                                  */
/*                                                                          */
/****************************************************************************/


/*
        General DDD Module
 */
void      DDD_Init (int argc, char **argv);
void      DDD_Exit (void);
void      DDD_Status (void);

/*
        Redirect line-oriented output, new in V1.2
 */
void      DDD_LineOutRegister (void (*func)(char *s));


/*
        Definition Module
 */
DDD_TYPE  DDD_StructRegister (char *name, OBJ_HEADER *hdr, ...);
void      DDD_HandlerRegister (DDD_TYPE id, ...);
void      DDD_StructDisplay (DDD_TYPE id);



/*
        Object Properties
 */
void      DDD_PrioritySet (OBJECT obj, int prio);
void      DDD_RangeSet (OBJ_HEADER *obj, int range);
int  *    DDD_InfoProcList (OBJECT obj);



/*
        Identification Module
 */
void      DDD_IdentifyBegin (void);
void      DDD_IdentifyEnd (void);
void      DDD_IdentifyNumber (OBJECT obj, PROC_ID proc, int ident);
void      DDD_IdentifyString (OBJECT obj, PROC_ID proc, char *ident);
void      DDD_IdentifyObject (OBJECT obj, PROC_ID proc, OBJECT ident);



/*
        Interface Module
 */
INTERFACE DDD_IFDefine (int nO, int O[], int nA, int A[], int nB, int B[]);
void      DDD_DisplayIF (void);
void      DDD_IFExchange (INTERFACE ifId,
                          int itemSize, ComProcPtr GatherProc, ComProcPtr ScatterProc);
void      DDD_IFOneway (INTERFACE ifId, int dir,
                        int itemSize, ComProcPtr GatherProc, ComProcPtr ScatterProc);

/*
        Range-oriented communication, new in V1.1
 */
void      DDD_IFRExchange (INTERFACE ifId, int range,
                           int itemSize, ComProcPtr GatherProc, ComProcPtr ScatterProc);
void      DDD_IFROneway (INTERFACE ifId, int range, int dir,
                         int itemSize, ComProcPtr GatherProc, ComProcPtr ScatterProc);

/*
        Extended communication, new in V1.1
 */
void      DDD_IFExchangeX (INTERFACE ifId,
                           int itemSize, ComProcXPtr GatherProc, ComProcXPtr ScatterProc);
void      DDD_IFOnewayX (INTERFACE ifId, int dir,
                         int itemSize, ComProcXPtr GatherProc, ComProcXPtr ScatterProc);



/*
        Transfer Module
 */
void      DDD_XferBegin (void);
void      DDD_XferEnd (void);
void      DDD_XferCopyObj (OBJECT obj, PROC_ID proc, int prio);
void      DDD_XferCopyObjX (OBJECT obj, PROC_ID proc, int prio, int size);
void      DDD_XferAddData (int cnt, DDD_TYPE typ);
void      DDD_XferAddDataX (int cnt, DDD_TYPE typ, int sizes[]);
void      DDD_XferDeleteObj (OBJECT obj);


/*
        Object Manager
 */
void *    DDD_ObjGet (DDD_TYPE id);
void *    DDD_ObjGetX (DDD_TYPE id, int size);       /* new in V1.2 */
void      DDD_ObjDelete (OBJECT obj);
void      DDD_ObjMoveTo (OBJECT obj, char *mem);
void      DDD_ObjUnlink (OBJECT obj);


/*
        Maintainance & Debugging
 */
void      DDD_ConsCheck (void);
void      DDD_ListLocalObjects (void);



#ifdef __cplusplus
}
#endif
#endif
