// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
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
/*            97/02/10 kb  started with CPP_FRONTEND implementation         */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#ifndef __DDD__
#define __DDD__


#define DDD_VERSION    "1.8.7"


/****************************************************************************/
/*                                                                          */
/* settings for switching C_FRONTEND/CPP_FRONTEND/F_FRONTEND                */
/*                                                                          */
/****************************************************************************/

#ifdef DDD_FRONTEND_C
#define C_FRONTEND
#endif

#ifdef DDD_FRONTEND_CPP
#define CPP_FRONTEND
#endif

#ifdef DDD_FRONTEND_F
#define F_FRONTEND
#endif

/* check FRONTEND-setting for plausibility */
#if defined(C_FRONTEND) && defined(CPP_FRONTEND)
#error DDD Configuration Error: C_FRONTEND and CPP_FRONTEND are set.
#endif

#if defined(C_FRONTEND) && defined(F_FRONTEND)
#error DDD Configuration Error: C_FRONTEND and F_FRONTEND are set.
#endif

#if defined(CPP_FRONTEND) && defined(F_FRONTEND)
#error DDD Configuration Error: CPP_FRONTEND and F_FRONTEND are set.
#endif



/* default frontend is C_FRONTEND */
#ifndef C_FRONTEND
 #ifndef CPP_FRONTEND
  #ifndef F_FRONTEND
   #define C_FRONTEND
  #endif
 #endif
#endif


/* helpful macros for FRONTEND switching, will be #undef'd when ddd.h ends */
#ifdef F_FRONTEND
#define _FPTR     *
#define _OBJREF   DDD_TYPE *, DDD_OBJ *
#else
#define _FPTR
#define _OBJREF   DDD_HDR
#endif


/* F77SYM(lsym,usym) macro is defined in compiler.h. 961127 KB */


#ifdef __cplusplus
 #ifndef CPP_FRONTEND
extern "C" {
 #endif
#endif



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
/* NOTE: changes must be also done in fddd.f */
enum ElemType {
  EL_DDDHDR   =  0,                     /* element type: DDD header             */
  EL_GDATA    = -1,                     /* element type: global data            */
  EL_LDATA    = -2,                     /* element type: local data             */
  EL_GBITS    = -3,                     /* element type: bitwise, 1=global      */
  EL_DATAPTR  = -4,                     /* element type: data pointer           */
  EL_OBJPTR   = -5,                     /* element type: object pointer         */
  EL_CONTINUE = -6,                     /* continued element definition list    */
  EL_END      = -7                      /* end of element definition list       */
};




/* options for DDD_SetOption */
/* NOTE: changes must be also done in fddd.f */
enum OptionType {
  OPT_IDENTIFY_MODE=0,             /* one of the IDMODE_xxx constants           */

  OPT_WARNING_VARSIZE_OBJ=8,       /* warning on differing obj sizes            */
  OPT_WARNING_SMALLSIZE,           /* warning on obj sizes smaller than declared*/
  OPT_WARNING_PRIOCHANGE,          /* warning on inconsistency in prio-change   */
  OPT_WARNING_DESTRUCT_HDR,        /* warning on inconsistency in HdrDestructor */
  OPT_WARNING_REF_COLLISION,       /* warning on collision in reference-localize*/
  OPT_WARNING_OLDSTYLE,            /* warning on usage of old-style ddd-funcs   */

  OPT_QUIET_CONSCHECK=16,          /* do ConsCheck in a quiet manner            */
  OPT_DEBUG_XFERMESGS,             /* print debug info for xfer messages        */
  OPT_INFO_XFER,                   /* display some statistical info during xfer */
  OPT_INFO_IF_WITH_ATTR,           /* display interfaces detailed (with attrs)  */

  OPT_XFER_PRUNE_DELETE,           /* prune del-cmd in del/xfercopy-combination */

  OPT_IF_REUSE_BUFFERS,            /* reuse interface buffs as long as possible */
  OPT_IF_CREATE_EXPLICIT,          /* dont (re-)create interfaces automatically */

  OPT_END
};


/* NOTE: changes must be also done in fddd.f */
enum SwitchType {
  OPT_OFF = 0,
  OPT_ON,

  IDMODE_LISTS = 10,        /* ordering of each identify-tupel is relevant      */
  IDMODE_SETS,              /* ordering of each identify-tupel is not sensitive */

  XFER_SHOW_NONE     = 0x0000,        /* show no statistical infos              */
  XFER_SHOW_OBSOLETE = 0x0100,        /* show #obsolete xfer-commands           */
  XFER_SHOW_MEMUSAGE = 0x0200         /* show sizes of message buffers          */
};




/* direction of interface communication (DDD_IFOneway) */
/* NOTE: changes must be also done in fddd.f */
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
  /* application may add small integers in order to get more
     stream-of-byte channels, up to DDD_USER_DATA_MAX */
  DDD_USER_DATA     = 0x4000,
  DDD_USER_DATA_MAX = 0x4fff,


  /* additional parameter for user defined handlers in xfer */

  /* object has been rejected due to RULE C3 */
  XFER_REJECT  = 0x9000,

  /* object has been upgraded due to RULE C3 */
  XFER_UPGRADE,

  /* object has been downgraded due to PruneDel */
  XFER_DOWNGRADE,

  /* object is totally new */
  XFER_NEW
};


/* several default modes for priority handling */
enum PrioMatrixDefaults {
  PRIOMERGE_MAXIMUM = 0,
  PRIOMERGE_MINIMUM
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

        #ifdef C_FRONTEND
  char empty[4];                 /* 4 unused bytes in current impl. */
        #endif
} DDD_HEADER;



/*
        new DDD types, used during access of DDD functional interface
 */
typedef unsigned int DDD_GID;
/*
   typedef unsigned short   DDD_TYPE;
   typedef unsigned short   DDD_IF;
   typedef unsigned short   DDD_PROC;
   typedef unsigned short   DDD_PRIO;
   typedef unsigned short   DDD_ATTR;
 */
typedef unsigned int DDD_TYPE;
typedef unsigned int DDD_IF;
typedef unsigned int DDD_PROC;
typedef unsigned int DDD_PRIO;
typedef unsigned int DDD_ATTR;

#ifdef CPP_FRONTEND
typedef unsigned int DDD_INDEX;
#endif
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
typedef char           * DDD_OBJ;
#endif
#ifdef F_FRONTEND
typedef int DDD_OBJ;
#endif

typedef DDD_HEADER     * DDD_HDR;
typedef unsigned int DDD_OPTION;


/* NULL values for DDD types */
#define DDD_TYPE_NULL  0
#define DDD_PROC_NULL  0
#define DDD_PRIO_NULL  0
#define DDD_ATTR_NULL  0


/* types of handlers for HandlerRegister */
/* this is supported for downward compatibility only,
   will be removed in further versions. kb 970212, ddd-1.7.8 */

enum Handlers {
  HANDLER_LDATACONSTRUCTOR = 0,
  HANDLER_DESTRUCTOR,
  HANDLER_DELETE,
  HANDLER_UPDATE,
  HANDLER_OBJMKCONS,
  HANDLER_SETPRIORITY,
  HANDLER_XFERCOPY,
  HANDLER_XFERDELETE,
  HANDLER_XFERGATHER,
  HANDLER_XFERSCATTER,
  HANDLER_XFERGATHERX,
  HANDLER_XFERSCATTERX,
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  HANDLER_XFERCOPYMANIP,
#endif
#ifdef F_FRONTEND
  HANDLER_ALLOCOBJ,        /* handler: get an OBJ_ID (Array index) */
  HANDLER_FREEOBJ,         /* handler: free an OBJ_ID              */
#endif
  HANDLER_END=999
};


/* handler prototypes */

/* handlers related to certain DDD_TYPE (i.e., member functions) */
typedef void (*HandlerLDATACONSTRUCTOR)(DDD_OBJ _FPTR);
typedef void (*HandlerDESTRUCTOR)(DDD_OBJ _FPTR);
typedef void (*HandlerDELETE)(DDD_OBJ _FPTR);
typedef void (*HandlerUPDATE)(DDD_OBJ _FPTR);
typedef void (*HandlerOBJMKCONS)(DDD_OBJ _FPTR, int _FPTR);
typedef void (*HandlerSETPRIORITY)(DDD_OBJ _FPTR, DDD_PRIO _FPTR);
typedef void (*HandlerXFERCOPY)(DDD_OBJ _FPTR, DDD_PROC _FPTR, DDD_PRIO _FPTR);
typedef void (*HandlerXFERDELETE)(DDD_OBJ _FPTR);
typedef void (*HandlerXFERGATHER)(DDD_OBJ _FPTR, int _FPTR, DDD_TYPE _FPTR, void *);
typedef void (*HandlerXFERSCATTER)(DDD_OBJ _FPTR, int _FPTR, DDD_TYPE _FPTR, void *, int _FPTR);
typedef void (*HandlerXFERGATHERX)(DDD_OBJ _FPTR, int _FPTR, DDD_TYPE _FPTR, char **);
typedef void (*HandlerXFERSCATTERX)(DDD_OBJ _FPTR, int _FPTR, DDD_TYPE _FPTR, char **, int _FPTR);
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
typedef void (*HandlerXFERCOPYMANIP)(DDD_OBJ _FPTR);
#endif
#ifdef F_FRONTEND
typedef void (*HandlerALLOCOBJ)(DDD_OBJ _FPTR);
typedef void (*HandlerFREEOBJ)(DDD_OBJ _FPTR);
#endif

/* handlers not related to DDD_TYPE (i.e., global functions) */
typedef DDD_TYPE (*HandlerGetRefType)(DDD_OBJ _FPTR, DDD_OBJ _FPTR);



#if defined(C_FRONTEND) || defined(F_FRONTEND)
typedef int (*ExecProcPtr)(DDD_OBJ _FPTR);
typedef int (*ExecProcXPtr)(DDD_OBJ _FPTR, DDD_PROC _FPTR, DDD_PRIO _FPTR);
typedef int (*ComProcPtr)(DDD_OBJ _FPTR, void *);
typedef int (*ComProcXPtr)(DDD_OBJ _FPTR, void *, DDD_PROC _FPTR, DDD_PRIO _FPTR);
#endif



/* special feature: hybrid reftype at TypeDefine-time */
#define DDD_TYPE_BY_HANDLER   127   /* must be > MAX_TYPEDESC */


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
/*#define DDD_OBJ(o)     (&((o)->ddd))*/



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


#ifdef CPP_FRONTEND

class DDD_Library
{
public:
  DDD_Library (int *, char ***);
  ~DDD_Library ();

  // DDD_Library is a Singleton
  static DDD_Library* Instance();

  void Status (void);
  void SetOption (DDD_OPTION, int);

  DDD_PROC InfoMe (void);
  DDD_PROC InfoMaster (void);
  DDD_PROC InfoProcs (void);
  int IsMaster (void)  {
    return InfoMe()==InfoMaster();
  }

  void LineOutRegister (void (*func)(char *s));

  // from TypeManager
  /* TODO, class DDD_Type with derived classes DDD_TypeStruct/Index?? */
  // for struct-like objects
  DDD_TYPE TypeDeclareStruct (char* n="");
  DDD_TYPE TypeDeclare (char* n="")
  {
    return TypeDeclareStruct(n);
  }

  // for array-like objects
  DDD_TYPE TypeDeclareIndex (int s=0, char* n="");
  DDD_TYPE TypeDeclare (int s, char* n="")
  {
    return TypeDeclareIndex(s,n);
  }

  // for struct- or array-like objects
  void TypeDefine (DDD_TYPE, ...);
  void TypeChangeName (DDD_TYPE, char*);
  void TypeDisplay (DDD_TYPE);
  int InfoTypes (void);
  int InfoHdrOffset (DDD_TYPE);

  // handler registration
  void SetHandlerLDATACONSTRUCTOR (DDD_TYPE, HandlerLDATACONSTRUCTOR);
  void SetHandlerDESTRUCTOR       (DDD_TYPE, HandlerDESTRUCTOR);
  void SetHandlerDELETE           (DDD_TYPE, HandlerDELETE);
  void SetHandlerUPDATE           (DDD_TYPE, HandlerUPDATE);
  void SetHandlerOBJMKCONS        (DDD_TYPE, HandlerOBJMKCONS);
  void SetHandlerSETPRIORITY      (DDD_TYPE, HandlerSETPRIORITY);
  void SetHandlerXFERCOPY         (DDD_TYPE, HandlerXFERCOPY);
  void SetHandlerXFERDELETE       (DDD_TYPE, HandlerXFERDELETE);
  void SetHandlerXFERGATHER       (DDD_TYPE, HandlerXFERGATHER);
  void SetHandlerXFERSCATTER      (DDD_TYPE, HandlerXFERSCATTER);
  void SetHandlerXFERGATHERX      (DDD_TYPE, HandlerXFERGATHERX);
  void SetHandlerXFERSCATTERX     (DDD_TYPE, HandlerXFERSCATTERX);
  void SetHandlerXFERCOPYMANIP    (DDD_TYPE, HandlerXFERCOPYMANIP);

  // from PrioManager
  void PrioMergeDefault (DDD_TYPE, int);
  void PrioMergeDefine (DDD_TYPE, DDD_PRIO, DDD_PRIO, DDD_PRIO);
  DDD_PRIO PrioMerge (DDD_TYPE, DDD_PRIO, DDD_PRIO);
  void PrioMergeDisplay (DDD_TYPE);

  // Identification
  void IdentifyBegin (void);
  void IdentifyEnd (void);

  // Transfer
  void XferBegin (void);
  void XferEnd (void);

  /* TODO some are missing here */

  // miscellaneous
  int ConsCheck (void);
  void ListLocalObjects (void);
  DDD_HDR SearchHdr (DDD_GID);


private:
  void ddd_TypeMgrInit (void);

private:
  static DDD_Library* _instance;
};


class DDD_Object
{
public:
  DDD_Object (DDD_TYPE, DDD_PRIO, DDD_ATTR a=0);
  void Init (DDD_TYPE, DDD_PRIO, DDD_ATTR a=0);
  DDD_Object ();
  virtual ~DDD_Object ();


  // object properties
  void PrioritySet (DDD_PRIO);
  void AttrSet (DDD_ATTR);               /* this shouldn't be allowed */
  int* InfoProcList (void);
  DDD_PROC InfoProcPrio (DDD_PRIO);
  int InfoIsLocal (void);
  int InfoNCopies (void);

  DDD_GID  InfoGlobalId (void) {
    return (DDD_GID) _hdr.gid;
  }
  DDD_TYPE InfoType (void)     {
    return (DDD_TYPE)_hdr.typ;
  }
  DDD_PRIO InfoPriority (void) {
    return (DDD_PRIO)_hdr.prio;
  }
  DDD_ATTR InfoAttr (void)     {
    return (DDD_ATTR)_hdr.attr;
  }

  /* Transfer */
  void XferCopyObj (DDD_PROC, DDD_PRIO);
  void XferDeleteObj (void);

  /* Identification */
  void IdentifyNumber (DDD_PROC, int);
  void IdentifyString (DDD_PROC, char *);
  void IdentifyObject (DDD_PROC, DDD_Object*);

  friend void DDD_Library::ddd_TypeMgrInit (void);


  // DDD Handlers as virtual functions
  // TODO: is this general enough?
  virtual void HandlerLDATACONSTRUCTOR (void) { }
  //virtual void HandlerDESTRUCTOR       (void) { }
  //virtual void HandlerDELETE           (void) { }
  virtual void HandlerUPDATE           (void) { }
  //virtual void HandlerOBJMKCONS        (int) { }
  //virtual void HandlerSETPRIORITY      (DDD_PRIO) { }
  virtual void HandlerXFERCOPY         (DDD_PROC, DDD_PRIO) { }
  //virtual void HandlerXFERDELETE       (void) { }
  //virtual void HandlerXFERGATHER       (int, DDD_TYPE, void *) { }
  //virtual void HandlerXFERSCATTER      (int, DDD_TYPE, void *, int) { }
  //virtual void HandlerXFERGATHERX      (int, DDD_TYPE, char **) { }
  //virtual void HandlerXFERSCATTERX     (int, DDD_TYPE, char **, int) { }

private:
  DDD_HEADER _hdr;
};



template<class T>
class DDD_ObjectOf : public DDD_Object
{
public:
  DDD_ObjectOf<T> (DDD_PRIO p, DDD_ATTR a=0)
  : DDD_Object(GetType(),p,a)
  { }

  inline static DDD_TYPE GetType (void)
  {
    if (_dddtype==0)
    {
      //T* p=0;
      _dddtype = DDD_Library::Instance()->TypeDeclareStruct();
    }
    return _dddtype;
  }

private:
  static DDD_TYPE _dddtype;
};

template<class T>
DDD_TYPE DDD_ObjectOf<T>::_dddtype = 0;



class DDD_IndexObject : public DDD_Object
{
public:
  DDD_IndexObject (DDD_TYPE, DDD_INDEX, DDD_PRIO, DDD_ATTR a=0);
  virtual ~DDD_IndexObject() {}

  DDD_INDEX Index()    {
    return _index;
  }
  operator int()       {
    return _index;
  }

private:
  DDD_INDEX _index;
};


template<class T>
class DDD_IndexObjectOf : public DDD_IndexObject
{
public:
  DDD_IndexObjectOf<T> (DDD_INDEX i, DDD_PRIO p, DDD_ATTR a=0)
  : DDD_IndexObject(GetType(),i,p,a)
  { }

  static DDD_TYPE GetType (void)
  {
    if (_dddtype==0)
    {
      //T* p=0;
      _dddtype = DDD_Library::Instance()->TypeDeclareIndex();
    }
    return _dddtype;
  }

private:
  static DDD_TYPE _dddtype;
};

template<class T>
DDD_TYPE DDD_IndexObjectOf<T>::_dddtype = 0;



class DDD_GatherScatter
{
public:
  virtual int Gather  (DDD_Object*, void*) = 0;
  virtual int Scatter (DDD_Object*, void*) = 0;
};

class DDD_GatherScatterX
{
public:
  virtual int Gather  (DDD_Object*, void*, DDD_PROC, DDD_PRIO) = 0;
  virtual int Scatter (DDD_Object*, void*, DDD_PROC, DDD_PRIO) = 0;
};

class DDD_Exec
{
public:
  virtual int Exec (DDD_Object*) = 0;
};

class DDD_ExecX
{
public:
  virtual int Exec (DDD_Object*, DDD_PROC, DDD_PRIO) = 0;
};


class DDD_Interface
{
public:
  DDD_Interface (int, DDD_TYPE O[],
                 int, DDD_PRIO A[], int, DDD_PRIO B[], char* n="");
  DDD_Interface (DDD_TYPE, DDD_PRIO, DDD_PRIO, char* n="");
  void SetName (char *);

  static void DisplayAll (void);
  void Display (void);

  static size_t InfoMemoryAll (void);
  size_t InfoMemory (void);


  void Exchange  (             size_t, DDD_GatherScatter*);
  void Exchange  (             size_t, DDD_GatherScatter&);
  void Oneway    (         int,size_t, DDD_GatherScatter*);
  void Oneway    (         int,size_t, DDD_GatherScatter&);
  /* TODO: NIY
                  void ExecLocal (                     DDD_Exec*);
                  void AExchange (DDD_ATTR,    size_t, DDD_GatherScatter*);
                  void AOneway   (DDD_ATTR,int,size_t, DDD_GatherScatter*);
                  void AExecLocal(DDD_ATTR,            DDD_Exec*);

                  void Exchange  (             size_t, DDD_GatherScatterX*);
                  void Oneway    (         int,size_t, DDD_GatherScatterX*);
                  void ExecLocal (                     DDD_ExecX*);
                  void AExchange (DDD_ATTR,    size_t, DDD_GatherScatterX*);
                  void AOneway   (DDD_ATTR,int,size_t, DDD_GatherScatterX*);
                  void AExecLocal(DDD_ATTR,            DDD_ExecX*);
   */

private:
  void Init (int,DDD_TYPE*, int,DDD_PRIO*, int,DDD_PRIO*, char*);

private:
  DDD_IF _id;
};

#endif



/*
        General DDD Module
 */
#if defined(C_FRONTEND) || defined(F_FRONTEND)
void     DDD_Init (int *argcp, char ***argvp);
void     DDD_Exit (void);
void     DDD_Status (void);
void     DDD_SetOption (DDD_OPTION _FPTR, int _FPTR);
#endif

#if defined(C_FRONTEND)
DDD_PROC DDD_InfoMe (void);
DDD_PROC DDD_InfoMaster (void);
DDD_PROC DDD_InfoProcs (void);
#endif


/*
        Redirect line-oriented output, new in V1.2
 */
#if defined(C_FRONTEND)
void     DDD_LineOutRegister (void (*func)(char *s));
#endif


/*
        Type Manager Module
 */
#ifdef F_FRONTEND
#define  DDD_TypeDeclare     F77SYM(ddd_typedeclare,DDD_TYPEDECLARE)
#define  DDD_TypeDefine      F77SYM(ddd_typedefine,DDD_TYPEDEFINE)
#define  DDD_TypeDisplay     F77SYM(ddd_typedisplay,DDD_TYPEDISPLAY)
#define  DDD_InfoTypes       F77SYM(ddd_infotypes,DDD_INFOTYPES)
#define  DDD_HandlerRegister F77SYM(ddd_handlerregister,DDD_HANDLERREGISTER)
/* TODO hier fehlen DDD_SetHandler-umsetzungen fuer fortran */

void     DDD_TypeDeclare (char *name, int *size, DDD_TYPE *type);
#endif
#ifdef C_FRONTEND
DDD_TYPE DDD_TypeDeclare (char *name);
int      DDD_InfoHdrOffset (DDD_TYPE);
#endif
#if defined(C_FRONTEND) || defined(F_FRONTEND)
void     DDD_TypeDefine (DDD_TYPE _FPTR, ...);
void     DDD_TypeDisplay (DDD_TYPE _FPTR);

/* oldstyle setting of DDD-handlers, will be removed in later versions */
void     DDD_HandlerRegister (DDD_TYPE _FPTR, ...);
#endif
int      DDD_InfoTypes (void);


#if defined(C_FRONTEND) || defined(F_FRONTEND)
/* newstyle, type-secure setting of handlers */
void     DDD_SetHandlerLDATACONSTRUCTOR(DDD_TYPE _FPTR, HandlerLDATACONSTRUCTOR);
void     DDD_SetHandlerDESTRUCTOR      (DDD_TYPE _FPTR, HandlerDESTRUCTOR);
void     DDD_SetHandlerDELETE          (DDD_TYPE _FPTR, HandlerDELETE);
void     DDD_SetHandlerUPDATE          (DDD_TYPE _FPTR, HandlerUPDATE);
void     DDD_SetHandlerOBJMKCONS       (DDD_TYPE _FPTR, HandlerOBJMKCONS);
void     DDD_SetHandlerSETPRIORITY     (DDD_TYPE _FPTR, HandlerSETPRIORITY);
void     DDD_SetHandlerXFERCOPY        (DDD_TYPE _FPTR, HandlerXFERCOPY);
void     DDD_SetHandlerXFERDELETE      (DDD_TYPE _FPTR, HandlerXFERDELETE);
void     DDD_SetHandlerXFERGATHER      (DDD_TYPE _FPTR, HandlerXFERGATHER);
void     DDD_SetHandlerXFERSCATTER     (DDD_TYPE _FPTR, HandlerXFERSCATTER);
void     DDD_SetHandlerXFERGATHERX     (DDD_TYPE _FPTR, HandlerXFERGATHERX);
void     DDD_SetHandlerXFERSCATTERX    (DDD_TYPE _FPTR, HandlerXFERSCATTERX);
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
void     DDD_SetHandlerXFERCOPYMANIP   (DDD_TYPE _FPTR, HandlerXFERCOPYMANIP);
#endif
#ifdef F_FRONTEND
void     DDD_SetHandlerALLOCOBJ        (DDD_TYPE _FPTR, HandlerALLOCOBJ);
void     DDD_SetHandlerFREEOBJ         (DDD_TYPE _FPTR, HandlerFREEOBJ);
#endif
#endif


#ifdef C_FRONTEND
void     DDD_PrioMergeDefault (DDD_TYPE, int);
void     DDD_PrioMergeDefine (DDD_TYPE, DDD_PRIO, DDD_PRIO, DDD_PRIO);
DDD_PRIO DDD_PrioMerge (DDD_TYPE, DDD_PRIO, DDD_PRIO);
void     DDD_PrioMergeDisplay (DDD_TYPE);
#endif



/*
        Object Properties
 */
#ifndef CPP_FRONTEND
void     DDD_PrioritySet (DDD_HDR, DDD_PRIO);
void     DDD_AttrSet (DDD_HDR, DDD_ATTR); /* this shouldn't be allowed */
int  *   DDD_InfoProcList (DDD_HDR);
DDD_PROC DDD_InfoProcPrio (DDD_HDR, DDD_PRIO);
int      DDD_InfoIsLocal (DDD_HDR);
int      DDD_InfoNCopies (DDD_HDR);
size_t   DDD_InfoCplMemory (void);
#endif

#ifdef F_FRONTEND
#define  DDD_InfoPriority    F77SYM(ddd_infopriority,DDD_INFOPRIORITY)
#endif


/*
        Identification Module
 */
#ifdef F_FRONTEND
#define  DDD_IdentifyBegin  F77SYM(ddd_identifybegin,DDD_IDENTIFYBEGIN)
#define  DDD_IdentifyEnd    F77SYM(ddd_identifyend,DDD_IDENTIFYEND)
#define  DDD_IdentifyNumber F77SYM(ddd_identifynumber,DDD_IDENTIFYNUMBER)
#define  DDD_IdentifyString F77SYM(ddd_identifystring,DDD_IDENTIFYSTRING)
#define  DDD_IdentifyObject F77SYM(ddd_identifyobject,DDD_IDENTIFYOBJECT)
#endif
#if defined(C_FRONTEND) || defined(F_FRONTEND)
void     DDD_IdentifyBegin (void);
void     DDD_IdentifyEnd (void);
void     DDD_IdentifyNumber (_OBJREF, DDD_PROC _FPTR, int _FPTR);
void     DDD_IdentifyString (_OBJREF, DDD_PROC _FPTR, char *);
void     DDD_IdentifyObject (_OBJREF, DDD_PROC _FPTR, _OBJREF);
#endif


/*
        Interface Module
 */
#ifdef F_FRONTEND
#define DDD_IFDefine        F77SYM(ddd_ifdefine,DDD_IFDEFINE)
#define DDD_IFDisplayAll    F77SYM(ddd_ifdisplayall,DDD_IFDISPLAYALL)
#define DDD_IFDisplay       F77SYM(ddd_ifdisplay,DDD_IFDISPLAY)
#define DDD_IFInfoMemoryAll F77SYM(ddd_ifinfomemoryall,DDD_IFINFOMEMORYALL)
#define DDD_IFInfoMemory    F77SYM(ddd_ifinfomemory,DDD_IFINFOMEMORY)
#define DDD_IFRefreshAll    F77SYM(ddd_ifrefreshall,DDD_IFREFRESHALL)
#define DDD_IFExchange      F77SYM(ddd_ifexchange,DDD_IFEXCHANGE)
#define DDD_IFOneway        F77SYM(ddd_ifoneway,DDD_IFONEWAY)
#define DDD_IFExecLocal     F77SYM(ddd_ifexeclocal,DDD_IFEXECLOCAL)
#define DDD_IFAExchange     F77SYM(ddd_ifaexchange,DDD_IFAEXCHANGE)
#define DDD_IFAOneway       F77SYM(ddd_ifaoneway,DDD_IFAONEWAY)
#define DDD_IFAExecLocal    F77SYM(ddd_ifaexeclocal,DDD_IFAEXECLOCAL)
#define DDD_IFExchangeX     F77SYM(ddd_ifexchangex,DDD_IFEXCHANGEX)
#define DDD_IFOnewayX       F77SYM(ddd_ifonewayx,DDD_IFONEWAYX)
#define DDD_IFExecLocalX    F77SYM(ddd_ifexeclocalx,DDD_IFEXECLOCALX)
#define DDD_IFAExchangeX    F77SYM(ddd_ifaexchangex,DDD_IFAEXCHANGEX)
#define DDD_IFAOnewayX      F77SYM(ddd_ifaonewayx,DDD_IFAONEWAYX)
#define DDD_IFAExecLocalX   F77SYM(ddd_ifaexeclocalx,DDD_IFAEXECLOCALX)

void     DDD_IFDefine (int *, DDD_TYPE *, int *, DDD_PRIO *, int *, DDD_PRIO *, DDD_IF *);
#endif

#ifdef C_FRONTEND
DDD_IF   DDD_IFDefine (int, DDD_TYPE O[], int, DDD_PRIO A[], int, DDD_PRIO B[]);
void     DDD_IFSetName (DDD_IF, char *);
#endif

#if defined(C_FRONTEND) || defined(F_FRONTEND)
void     DDD_IFDisplayAll (void);
void     DDD_IFDisplay (DDD_IF _FPTR);
size_t   DDD_IFInfoMemoryAll (void);
size_t   DDD_IFInfoMemory (DDD_IF _FPTR);
void     DDD_IFRefreshAll (void);

void     DDD_IFExchange   (DDD_IF _FPTR,                         size_t _FPTR, ComProcPtr,ComProcPtr);
void     DDD_IFOneway     (DDD_IF _FPTR,               int _FPTR,size_t _FPTR, ComProcPtr,ComProcPtr);
void     DDD_IFExecLocal  (DDD_IF _FPTR,                                       ExecProcPtr);
void     DDD_IFAExchange  (DDD_IF _FPTR,DDD_ATTR _FPTR,          size_t _FPTR, ComProcPtr,ComProcPtr);
void     DDD_IFAOneway    (DDD_IF _FPTR,DDD_ATTR _FPTR,int _FPTR,size_t _FPTR, ComProcPtr,ComProcPtr);
void     DDD_IFAExecLocal (DDD_IF _FPTR,DDD_ATTR _FPTR,                        ExecProcPtr);
void     DDD_IFExchangeX  (DDD_IF _FPTR,                         size_t _FPTR, ComProcXPtr,ComProcXPtr);
void     DDD_IFOnewayX    (DDD_IF _FPTR,               int _FPTR,size_t _FPTR, ComProcXPtr,ComProcXPtr);
void     DDD_IFExecLocalX (DDD_IF _FPTR,                                       ExecProcXPtr);
void     DDD_IFAExchangeX (DDD_IF _FPTR,DDD_ATTR _FPTR,          size_t _FPTR, ComProcXPtr,ComProcXPtr);
void     DDD_IFAOnewayX   (DDD_IF _FPTR,DDD_ATTR _FPTR,int _FPTR,size_t _FPTR, ComProcXPtr,ComProcXPtr);
void     DDD_IFAExecLocalX(DDD_IF _FPTR,DDD_ATTR _FPTR,                        ExecProcXPtr);
#endif


/*
        Transfer Module
 */
#ifdef F_FRONTEND
#define DDD_XferBegin     F77SYM(ddd_xferbegin,DDD_XFERBEGIN)
#define DDD_XferEnd       F77SYM(ddd_xferend,DDD_XFEREND)
#define DDD_XferCopyObj   F77SYM(ddd_xfercopyobj,DDD_XFERCOPYOBJ)
#define DDD_XferDeleteObj F77SYM(ddd_xferdeleteobj,DDD_XFERDELETEOBJ)
#endif
#ifdef C_FRONTEND
void     DDD_XferAddData (int _FPTR, DDD_TYPE _FPTR);
void     DDD_XferAddDataX (int _FPTR, DDD_TYPE _FPTR, size_t sizes[]);
#endif
#if defined(C_FRONTEND) || defined(F_FRONTEND)
void     DDD_XferBegin (void);
void     DDD_XferEnd (void);
void     DDD_XferCopyObj (_OBJREF, DDD_PROC _FPTR, DDD_PRIO _FPTR);
void     DDD_XferCopyObjX (_OBJREF, DDD_PROC _FPTR, DDD_PRIO _FPTR, size_t _FPTR);
void     DDD_XferDeleteObj (_OBJREF);
#endif

/*
        Object Manager
 */
#ifdef F_FRONTEND
#define DDD_ObjGet   F77SYM(ddd_objget,DDD_OBJGET)
#define DDD_ObjUnGet F77SYM(ddd_objunget,DDD_OBJUNGET)

void     DDD_ObjGet (DDD_TYPE *, DDD_PRIO *, DDD_ATTR *, DDD_OBJ *);
void     DDD_ObjUnGet (DDD_OBJ *, DDD_TYPE *);
#endif
#ifdef C_FRONTEND
DDD_OBJ  DDD_ObjNew (size_t, DDD_TYPE, DDD_PRIO, DDD_ATTR);
void     DDD_ObjDelete (DDD_OBJ, size_t, DDD_TYPE);
void     DDD_HdrConstructor (DDD_HDR, DDD_TYPE, DDD_PRIO, DDD_ATTR);
void     DDD_HdrConstructorMove (DDD_HDR, DDD_HDR);
void     DDD_HdrDestructor (DDD_HDR);
DDD_OBJ  DDD_ObjGet (size_t, DDD_TYPE, DDD_PRIO, DDD_ATTR);
void     DDD_ObjUnGet (DDD_HDR, size_t);
#endif


/*
        Special functions for F_FRONTEND
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
#if defined(C_FRONTEND) || defined(F_FRONTEND)
int      DDD_ConsCheck (void);  /* returns total #errors since V1.6.6 */
void     DDD_ListLocalObjects (void);
DDD_HDR  DDD_SearchHdr (DDD_GID _FPTR);
#endif


/****************************************************************************/

#undef _FPTR
#undef _OBJREF


#ifdef __cplusplus
 #ifndef CPP_FRONTEND
}
 #endif
#endif

/****************************************************************************/

#endif
