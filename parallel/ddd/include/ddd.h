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
/*            98/01/28 kb  new ddd-environment: Join.                       */
/*            98/05/14 kb  redesigned memory handling.                      */
/*            98/07/20 kb  new ddd-environment: Prio.                       */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/* RCS_ID
   $Header$
 */

#ifndef __DDD__
#define __DDD__

/* for size_t */
#include <cstddef>

#include "namespace.h"

START_UGDIM_NAMESPACE

#define DDD_VERSION    "1.9"


/****************************************************************************/
/*                                                                          */
/* settings for switching C_FRONTEND/CPP_FRONTEND                           */
/*                                                                          */
/****************************************************************************/

#ifdef DDD_FRONTEND_C
#define C_FRONTEND
#endif

#ifdef DDD_FRONTEND_CPP
#define CPP_FRONTEND
#endif


/* check FRONTEND-setting for plausibility */
#if defined(C_FRONTEND) && defined(CPP_FRONTEND)
#error DDD Configuration Error: C_FRONTEND and CPP_FRONTEND are set.
#endif




/* default frontend is C_FRONTEND */
#ifndef C_FRONTEND
 #ifndef CPP_FRONTEND
   #define C_FRONTEND
 #endif
#endif


/* helpful macros for FRONTEND switching, will be #undef'd when ddd.h ends */
#define _FPTR
#define _OBJREF   DDD_HDR


/* F77SYM(lsym,usym) macro is defined in compiler.h. 961127 KB */

/*
   #ifdef __cplusplus
   #ifndef CPP_FRONTEND
   extern "C" {
   #endif
   #endif
 */

/****************************************************************************/
/*                                                                          */
/* compile time constants defining static data size (i.e. arrays)           */
/* other constants                                                          */
/*                                                                          */
/****************************************************************************/



/* return types for DDD functions */
/* NOTE: changes must be also done in fddd.f */
typedef enum {
  DDD_RET_OK            = 0,            /* function was executed ok             */
  DDD_RET_ERROR_UNKNOWN = 1,            /* unknown error condition              */
  DDD_RET_ERROR_NOMEM   = 2             /* function aborted due to mem shortage */
} DDD_RET;


/* types of elements for StructRegister */
/* (use negative values for combination with positive DDD_TYPEs) */
/* NOTE: changes must be also done in fddd.f */
typedef enum {
  EL_DDDHDR   =  0,                     /* element type: DDD header             */
  EL_GDATA    = -1,                     /* element type: global data            */
  EL_LDATA    = -2,                     /* element type: local data             */
  EL_GBITS    = -3,                     /* element type: bitwise, 1=global      */
  EL_DATAPTR  = -4,                     /* element type: data pointer           */
  EL_OBJPTR   = -5,                     /* element type: object pointer         */
  EL_CONTINUE = -6,                     /* continued element definition list    */
  EL_END      = -7                      /* end of element definition list       */
} DDD_ELEM_TYPE;




/* options for DDD_SetOption */
/* NOTE: changes must be also done in fddd.f */
typedef enum {
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
  OPT_INFO_JOIN,                   /* display some statistical info during join */
  OPT_INFO_IF_WITH_ATTR,           /* display interfaces detailed (with attrs)  */

  OPT_XFER_PRUNE_DELETE,           /* prune del-cmd in del/xfercopy-combination */

  OPT_IF_REUSE_BUFFERS,            /* reuse interface buffs as long as possible */
  OPT_IF_CREATE_EXPLICIT,          /* dont (re-)create interfaces automatically */

  OPT_CPLMGR_USE_FREELIST,         /* use freelist for coupling-memory (default)*/

  OPT_END
} DDD_OPTION;


/* NOTE: changes must be also done in fddd.f */
enum OptConsts {
  OPT_OFF = 0,
  OPT_ON
};

enum OptConstIdent {
  IDMODE_LISTS = 1,         /* ordering of each identify-tupel is relevant      */
  IDMODE_SETS               /* ordering of each identify-tupel is not sensitive */
};

enum OptConstXfer {
  XFER_SHOW_NONE     = 0x0000,        /* show no statistical infos              */
  XFER_SHOW_OBSOLETE = 0x0001,        /* show #obsolete xfer-commands           */
  XFER_SHOW_MEMUSAGE = 0x0002,        /* show sizes of message buffers          */
  XFER_SHOW_MSGSALL  = 0x0004         /* show message contents by LowComm stats */
};

enum OptConstJoin {
  JOIN_SHOW_NONE     = 0x0000,        /* show no statistical infos              */
  JOIN_SHOW_OBSOLETE = 0x0001,        /* show #obsolete join-commands           */
  JOIN_SHOW_MEMUSAGE = 0x0002,        /* show sizes of message buffers          */
  JOIN_SHOW_MSGSALL  = 0x0004         /* show message contents by LowComm stats */
};



/* direction of interface communication (DDD_IFOneway) */
/* NOTE: changes must be also done in fddd.f */
typedef enum {
  IF_FORWARD  = 1,                     /* communicate from A to B               */
  IF_BACKWARD = 2                      /* communicate from B to A               */
} DDD_IF_DIR;


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
  XFER_NEW,


  /* return value for DDD_XferIsPrunedDelete */
  XFER_PRUNED_TRUE = 0x9100,
  XFER_PRUNED_FALSE,
  XFER_PRUNED_ERROR,


  /* return value for DDD_XferObjIsResent */
  XFER_RESENT_TRUE = 0x9200,
  XFER_RESENT_FALSE,
  XFER_RESENT_ERROR
};


/* several default modes for priority handling */
enum PrioMatrixDefaults {
  PRIOMERGE_MAXIMUM = 0,
  PRIOMERGE_MINIMUM
};



/* constants for management of temporary memory allocation/deletion */
enum TMemRequests {
  TMEM_ANY     = 0x0000,
  TMEM_MSG,
  TMEM_OBJLIST,
  TMEM_CPL,

  TMEM_XFER    = 0x1000,
  TMEM_LOWCOMM,

  TMEM_JOIN    = 0x2000,

  TMEM_CONS    = 0x3000,

  TMEM_IDENT   = 0x4000
};



/****************************************************************************/
/*                                                                          */
/* data structures and new types                                            */
/*                                                                          */
/****************************************************************************/


/*
        new DDD types, used during access of DDD functional interface
 */
#ifdef DDD_GID_DEBUG
struct ddd_gid_debug
{
  unsigned int val;
  /* ddd_gid_debug(unsigned int v) : val(v) {} */
  /* ddd_gid_debug() : val(0) {} */
  bool operator < (const ddd_gid_debug & other) { return val < other.val; }
  bool operator > (const ddd_gid_debug & other) { return val > other.val; }
  bool operator < (int i) { return val < i; }
  bool operator > (int i) { return val > i; }
  bool operator == (const ddd_gid_debug & other) { return val == other.val; }
  bool operator != (const ddd_gid_debug & other) { return val != other.val; }
  bool operator == (int i) { return val == i; }
  bool operator != (int i) { return val != i; }
  ddd_gid_debug operator ++ (int) { ddd_gid_debug x(*this); x.val++; return x; }
  ddd_gid_debug operator + (int i) { ddd_gid_debug x(*this); x.val = x.val + i; return x; }
  ddd_gid_debug operator - (int i) { ddd_gid_debug x(*this); x.val = x.val - i; return x; }
  ddd_gid_debug& operator= (unsigned int v) { val = v; return *this; }
  ddd_gid_debug& operator= (unsigned long int v) { val = v; return *this; }
  ddd_gid_debug& operator= (long int v) { val = v; return *this; }
  ddd_gid_debug& operator= (int v) { val = v; return *this; }
};
typedef ddd_gid_debug DDD_GID;
#define DDD_GID_TO_INT(A) A.val
#else
#ifdef DDD_GID_T
typedef DDD_GID_T DDD_GID;
#else
typedef unsigned long DDD_GID;
#endif
#define DDD_GID_TO_INT(A) (unsigned int) A
#endif
typedef unsigned int DDD_TYPE;
typedef unsigned int DDD_IF;
typedef unsigned int DDD_PROC;
typedef unsigned int DDD_PRIO;
typedef unsigned int DDD_ATTR;

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
  DDD_GID gid;            /* global id */

        #ifdef C_FRONTEND
  char empty[4];                 /* 4 unused bytes in current impl. */
        #endif
} DDD_HEADER;

#ifdef CPP_FRONTEND
typedef unsigned int DDD_INDEX;
typedef char           * DDD_OBJ;
class DDD_Object;  // forward declaration
typedef DDD_Object     * DDD_HDR;
#endif
#ifdef C_FRONTEND
typedef char           * DDD_OBJ;
typedef DDD_HEADER     * DDD_HDR;
#endif

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
  HANDLER_XFERCOPYMANIP,
  HANDLER_END=999
};


/* handler prototypes */

/* handlers related to certain DDD_TYPE (i.e., member functions) */
#if defined(C_FRONTEND) || \
  (defined(CPP_FRONTEND) && ! defined(WITH_VIRTUAL_HANDLERS))

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
#endif
#if defined(C_FRONTEND)
typedef void (*HandlerXFERCOPYMANIP)(DDD_OBJ _FPTR);
#endif



/* handlers not related to DDD_TYPE (i.e., global functions) */
typedef DDD_TYPE (*HandlerGetRefType)(DDD_OBJ _FPTR, DDD_OBJ _FPTR);



#if defined(C_FRONTEND)
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
/*#ifndef __cplusplus*/
#ifdef C_FRONTEND
#define DDD_InfoPriority(ddd_hdr)    ((ddd_hdr)->prio)
#define DDD_InfoGlobalId(ddd_hdr)    ((ddd_hdr)->gid)
#define DDD_InfoAttr(ddd_hdr)       ((ddd_hdr)->attr)
#define DDD_InfoType(ddd_hdr)       ((ddd_hdr)->typ)
#endif
/*#endif*/


/****************************************************************************/
/*                                                                          */
/* declaration of DDD functional interface                                  */
/*                                                                          */
/****************************************************************************/


#ifdef CPP_FRONTEND

/**
        DDD Type class.
        Each DDD object has a previously specified DDD\_Type.

        \todoTBC
 */

// currently not used!
class DDD_Type
{
public:
  DDD_Type (DDD_TYPE type)  {
    _dddtype = type;
  }

  void ChangeName (char*);
  void Display();

private:
  DDD_TYPE _dddtype;
};



/**
        DDD Library class.
        Construct one single instance of the class in order to use
        the functionality of the DDD library.

        \todoTBC
 */

class DDD_Library
{
public:
  DDD_Library (int *, char ***);
  ~DDD_Library ();

  /// DDD\_Library is a Singleton
  static DDD_Library* Instance();

  /// shows status of DDD library
  void Status (void);
  void SetOption (DDD_OPTION, int);

  DDD_PROC InfoMe (void);
  DDD_PROC InfoMaster (void);
  DDD_PROC InfoProcs (void);
  int IsMaster (void)  {
    return InfoMe()==InfoMaster();
  }

  /// redirect DDD stdout by registering callback function
  void LineOutRegister (void (*func)(const char *s));

  // from TypeManager
  DDD_TYPE TypeDeclareStruct (const char* n="");
  DDD_TYPE TypeDeclare (const char* n="")
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
  void TypeDisplay(DDD_TYPE);
  int InfoTypes (void);
  int InfoHdrOffset (DDD_TYPE);

                #ifndef WITH_VIRTUAL_HANDLERS
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
                #endif

  // from PrioManager
  void PrioMergeDefault (DDD_TYPE, int);
  void PrioMergeDefine (DDD_TYPE, DDD_PRIO, DDD_PRIO, DDD_PRIO);
  DDD_PRIO PrioMerge (DDD_TYPE, DDD_PRIO, DDD_PRIO);
  void PrioMergeDisplay (DDD_TYPE);

  // Identification
  void    IdentifyBegin (void);
  DDD_RET IdentifyEnd (void);

  // Transfer
  void XferBegin (void);
  DDD_RET XferEnd (void);

  // Prio
  void PrioBegin (void);
  DDD_RET PrioEnd (void);

  // Join
  void JoinBegin (void);
  DDD_RET JoinEnd (void);

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


/****************************************************************************/

/**
        DDD Object class.
        Each distributed object (\ie, DDD object) will have to inherit
        from #DDD_Object#.

        \todoTBC
 */

class DDD_Object
{
public:
  DDD_Object (DDD_TYPE, DDD_PRIO, DDD_ATTR a=0);
  void Init (DDD_TYPE, DDD_PRIO, DDD_ATTR a=0);
  DDD_Object ();
  ~DDD_Object ();


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


                #ifdef WITH_VIRTUAL_HANDLERS
  // DDD Handlers as virtual functions
  virtual void HandlerLDATACONSTRUCTOR (void) { }
  virtual void HandlerDESTRUCTOR       (void) { }
  virtual void HandlerDELETE           (void) { }
  virtual void HandlerUPDATE           (void) { }
  virtual void HandlerOBJMKCONS        (int) { }
  virtual void HandlerSETPRIORITY      (DDD_PRIO) { }
  virtual void HandlerXFERCOPY         (DDD_PROC, DDD_PRIO) { }
  virtual void HandlerXFERDELETE       (void) { }
  virtual void HandlerXFERGATHER       (int, DDD_TYPE, void *) { }
  virtual void HandlerXFERSCATTER      (int, DDD_TYPE, void *, int) { }
  virtual void HandlerXFERGATHERX      (int, DDD_TYPE, char **) { }
  virtual void HandlerXFERSCATTERX     (int, DDD_TYPE, char **, int) { }
                #endif

  //protected:
public:            // currently open to public, make protected later
  DDD_HEADER _hdr;
};


/****************************************************************************/

/**
        DDD ObjectOf template class.
        This template class simplifies the usage of the base class
   #DDD_Object#.

        \todoTBC
 */

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


/****************************************************************************/

/**
        DDD IndexObject class.
        For objects stored in arrays, it is more convenient to
        use the #DDD_IndexObject# class instead of the common
        base class #DDD_Object#.

        \todoTBC
 */

class DDD_IndexObject : public DDD_Object
{
public:
  DDD_IndexObject (DDD_TYPE, DDD_INDEX, DDD_PRIO, DDD_ATTR a=0);
  ~DDD_IndexObject() {}

  DDD_INDEX Index()    {
    return _index;
  }
  operator int()       {
    return _index;
  }

private:
  DDD_INDEX _index;
};



/**
        DDD IndexObjectOf template class.
        This template class simplifies the usage of its base class
   #DDD_IndexObject#.

        \todoTBC
 */

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



/****************************************************************************/

/**
        DDD Communicator class.
        Base class for all communicator classes.

        \todoTBC
 */

class DDD_Communicator
{
  // empty
};


/**
        DDD GatherScatter class.
        This communicator specifies two member functions for packing
        and unpacking data for one distributed DDD object inside a
        \ddd{interface}.

        \todoTBC
 */

class DDD_GatherScatter : public DDD_Communicator
{
public:
  virtual int Gather  (DDD_Object*, void*) = 0;
  virtual int Scatter (DDD_Object*, void*) = 0;
};


/**
        DDD GatherScatterX class.
        This communicator specifies two member functions for packing
        and unpacking data for one distributed DDD object inside a
        \ddd{interface}.
        The handler functions #Gather# and #Scatter# take eXtended arguments.

        \todoTBC
 */

class DDD_GatherScatterX : public DDD_Communicator
{
public:
  virtual int Gather  (DDD_Object*, void*, DDD_PROC, DDD_PRIO) = 0;
  virtual int Scatter (DDD_Object*, void*, DDD_PROC, DDD_PRIO) = 0;
};


/**
        DDD Exec class.
        This communicator is used for local execution of a single
        function for each object inside a \ddd{interface}.

        \todoTBC
 */

class DDD_Exec : public DDD_Communicator
{
public:
  virtual int Exec (DDD_Object*) = 0;
};


/**
        DDD ExecX class.
        This communicator is used for local execution of a single
        function for each object inside a \ddd{interface}.
        The member function #Exec# takes eXtended arguments.

        \todoTBC
 */

class DDD_ExecX : public DDD_Communicator
{
public:
  virtual int Exec (DDD_Object*, DDD_PROC, DDD_PRIO) = 0;
};


/****************************************************************************/

/**
        DDD Interface class.
        This is the abstraction of distributed-graph overlap at processor
        borders.

        \todoTBC
 */

class DDD_Interface
{
public:
  DDD_Interface (int, DDD_TYPE O[],
                 int, DDD_PRIO A[], int, DDD_PRIO B[], char* n="");
  DDD_Interface (DDD_TYPE, DDD_PRIO, DDD_PRIO, char* n="");
  void SetName (const char *);

  static void DisplayAll (void);
  void Display (void);

  static size_t InfoMemoryAll (void);
  size_t InfoMemory (void);


  void Exchange  (                    size_t, DDD_GatherScatter*);
  void Exchange  (                    size_t, DDD_GatherScatter&);
  void Oneway    (         DDD_IF_DIR,size_t, DDD_GatherScatter*);
  void Oneway    (         DDD_IF_DIR,size_t, DDD_GatherScatter&);
  /* TODO: NIY
                  void ExecLocal (                            DDD_Exec*);
                  void AExchange (DDD_ATTR,           size_t, DDD_GatherScatter*);
                  void AOneway   (DDD_ATTR,DDD_IF_DIR,size_t, DDD_GatherScatter*);
                  void AExecLocal(DDD_ATTR,                   DDD_Exec*);

                  void Exchange  (                    size_t, DDD_GatherScatterX*);
                  void Oneway    (         DDD_IF_DIR,size_t, DDD_GatherScatterX*);
                  void ExecLocal (                            DDD_ExecX*);
                  void AExchange (DDD_ATTR,           size_t, DDD_GatherScatterX*);
                  void AOneway   (DDD_ATTR,DDD_IF_DIR,size_t, DDD_GatherScatterX*);
                  void AExecLocal(DDD_ATTR,                   DDD_ExecX*);
   */

private:
  void Init (int,DDD_TYPE*, int,DDD_PRIO*, int,DDD_PRIO*, char*);

private:
  DDD_IF _id;
};

#endif


/****************************************************************************/


/*
        General DDD Module
 */
#if defined(C_FRONTEND)
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
void     DDD_LineOutRegister (void (*func)(const char *s));
#endif


/*
        Type Manager Module
 */

#ifdef C_FRONTEND
DDD_TYPE DDD_TypeDeclare (const char *name);
int      DDD_InfoHdrOffset (DDD_TYPE);
#endif
#if defined(C_FRONTEND)
void     DDD_TypeDefine (DDD_TYPE _FPTR, ...);
void     DDD_TypeDisplay (DDD_TYPE _FPTR);

/* oldstyle setting of DDD-handlers, will be removed in later versions */
void     DDD_HandlerRegister (DDD_TYPE _FPTR, ...);
#endif
int      DDD_InfoTypes (void);


#if defined(C_FRONTEND)
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
void     DDD_SetHandlerXFERCOPYMANIP   (DDD_TYPE _FPTR, HandlerXFERCOPYMANIP);
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



/*
        Identification Environment Module
 */

#if defined(C_FRONTEND)
void     DDD_IdentifyBegin (void);
DDD_RET  DDD_IdentifyEnd (void);
void     DDD_IdentifyNumber (_OBJREF, DDD_PROC _FPTR, int _FPTR);
void     DDD_IdentifyString (_OBJREF, DDD_PROC _FPTR, char *);
void     DDD_IdentifyObject (_OBJREF, DDD_PROC _FPTR, _OBJREF);
#endif


/*
        Interface Module
 */

#ifdef C_FRONTEND
DDD_IF   DDD_IFDefine (int, DDD_TYPE O[], int, DDD_PRIO A[], int, DDD_PRIO B[]);
void     DDD_IFSetName (DDD_IF, const char *);
#endif

#if defined(C_FRONTEND)
void     DDD_IFDisplayAll (void);
void     DDD_IFDisplay (DDD_IF _FPTR);
size_t   DDD_IFInfoMemoryAll (void);
size_t   DDD_IFInfoMemory (DDD_IF _FPTR);
void     DDD_IFRefreshAll (void);

void     DDD_IFExchange   (DDD_IF _FPTR,                                size_t _FPTR, ComProcPtr,ComProcPtr);
void     DDD_IFOneway     (DDD_IF _FPTR,               DDD_IF_DIR _FPTR,size_t _FPTR, ComProcPtr,ComProcPtr);
void     DDD_IFExecLocal  (DDD_IF _FPTR,                                              ExecProcPtr);
void     DDD_IFAExchange  (DDD_IF _FPTR,DDD_ATTR _FPTR,                 size_t _FPTR, ComProcPtr,ComProcPtr);
void     DDD_IFAOneway    (DDD_IF _FPTR,DDD_ATTR _FPTR,DDD_IF_DIR _FPTR,size_t _FPTR, ComProcPtr,ComProcPtr);
void     DDD_IFAExecLocal (DDD_IF _FPTR,DDD_ATTR _FPTR,                               ExecProcPtr);
void     DDD_IFExchangeX  (DDD_IF _FPTR,                                size_t _FPTR, ComProcXPtr,ComProcXPtr);
void     DDD_IFOnewayX    (DDD_IF _FPTR,               DDD_IF_DIR _FPTR,size_t _FPTR, ComProcXPtr,ComProcXPtr);
void     DDD_IFExecLocalX (DDD_IF _FPTR,                                              ExecProcXPtr);
void     DDD_IFAExchangeX (DDD_IF _FPTR,DDD_ATTR _FPTR,                 size_t _FPTR, ComProcXPtr,ComProcXPtr);
void     DDD_IFAOnewayX   (DDD_IF _FPTR,DDD_ATTR _FPTR,DDD_IF_DIR _FPTR,size_t _FPTR, ComProcXPtr,ComProcXPtr);
void     DDD_IFAExecLocalX(DDD_IF _FPTR,DDD_ATTR _FPTR,                               ExecProcXPtr);
#endif


/*
        Transfer Environment Module
 */
#ifdef C_FRONTEND
int      DDD_XferWithAddData (void);
void     DDD_XferAddData (int _FPTR, DDD_TYPE _FPTR);
void     DDD_XferAddDataX (int _FPTR, DDD_TYPE _FPTR, size_t sizes[]);
int      DDD_XferIsPrunedDelete (_OBJREF);
int      DDD_XferObjIsResent (_OBJREF);
#endif
#if defined(C_FRONTEND)
void     DDD_XferBegin (void);
DDD_RET  DDD_XferEnd (void);
void     DDD_XferCopyObj (_OBJREF, DDD_PROC _FPTR, DDD_PRIO _FPTR);
void     DDD_XferCopyObjX (_OBJREF, DDD_PROC _FPTR, DDD_PRIO _FPTR, size_t _FPTR);
void     DDD_XferDeleteObj (_OBJREF);
void     DDD_XferPrioChange (_OBJREF, DDD_PRIO _FPTR);
#endif


/*
        Prio Environment Module
 */
#ifdef C_FRONTEND
void     DDD_PrioBegin (void);
DDD_RET  DDD_PrioEnd (void);
void     DDD_PrioChange (_OBJREF, DDD_PRIO _FPTR);
#endif



/*
        Join Environment Module
 */
#ifdef C_FRONTEND
void     DDD_JoinBegin (void);
DDD_RET  DDD_JoinEnd (void);
void     DDD_JoinObj (_OBJREF, DDD_PROC _FPTR, DDD_GID _FPTR);
#endif


/*
        Object Manager
 */

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
        Maintainance & Debugging
 */

#if defined(C_FRONTEND)
int      DDD_ConsCheck (void);  /* returns total #errors since V1.6.6 */
void     DDD_ListLocalObjects (void);
DDD_HDR  DDD_SearchHdr (DDD_GID _FPTR);
#endif


/****************************************************************************/

#undef _FPTR
#undef _OBJREF

/*
   #ifdef __cplusplus
   #ifndef CPP_FRONTEND
   }
   #endif
   #endif
 */
/****************************************************************************/

END_UGDIM_NAMESPACE

#endif
