// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      dddi.h                                                        */
/*                                                                          */
/* Purpose:   internal header file for ddd module                           */
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
/*            95/10/05 kb  added casts to mem-management macros             */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/



/*
   #define CheckIFMEM
   #define CheckPMEM
   #define CheckCplMEM
   #define CheckTmpMEM
   #define CheckMsgMEM
 */


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __DDDI_H__
#define __DDDI_H__

#include <assert.h>

#ifndef __COMPILER__
#include "compiler.h"
#endif

#include "include/memmgr.h"
#include "ppif.h"
#include "ctrl/stat.h"

#include "include/ddd.h"
#include "include/dddio.h"


/****************************************************************************/

/*
        macro in order to turn usage of register variables on/off
 */
#define REGISTER register

/*
        macro for exiting program in case of a severe error condition
 */
#define HARD_EXIT  assert(0)
/* #define HARD_EXIT  exit(1) */


/*
        macro for output of downward compatibility warning messages
 */
#define OLDSTYLE(txt) \
  if (me==master && DDD_GetOption(OPT_WARNING_OLDSTYLE)==OPT_ON)  \
  { DDD_PrintError('W', 8888, txt); }



/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/*** DDD internal parameters ***/

#define MAX_ELEMDESC   64    /* max. number of elements per TYPE_DESC       */
#define MAX_TYPEDESC   32    /* max. number of TYPE_DESC                    */
#define MAX_PRIO       32    /* max. number of DDD_PRIO                     */

#define MAX_OBJ    200000    /* max. number of locally registered objects   */
#define MAX_CPL     77700    /* max. number of local objects with coupling  */

#define MAX_TRIES 5000000    /* max. number of tries til timeout in IF-comm */

#define MAX_PROCBITS_IN_GID 10  /* this allows 2^10 procs and 2^22 objects  */

/* use maximum as default, if no Priomerge-matrix is available */
#define PRIOMERGE_DEFAULT PRIOMERGE_MAXIMUM



/*** DDD internal constants ***/

/* maximum number of procs allowed (limited by GID construction) */
#define MAX_PROCS   (1<<MAX_PROCBITS_IN_GID)


#define GID_INVALID  -1            /* invalid global id                     */
#define PRIO_INVALID (MAX_PRIO+1)  /* invalid priority                      */
#define PROC_INVALID (MAX_PROCS+1) /* invalid processor number              */
#define ERROR        -1            /* standard error indicator              */


/* types of virtual channels (for ppif interface) */
enum VChanType {
  VC_IDENT   = 15,               /* channels used for identification module     */
  VC_IFCOMM  = 16,               /* channels used for interface module          */
  VC_TOPO    = 17                /* channels used for xfer module (topology)    */
};


/* results of an prio-merge operation. see mgr/prio.c for more details. */
enum PrioMergeVals {
  PRIO_ERROR = -1,
  PRIO_UNKNOWN,
  PRIO_FIRST,
  PRIO_SECOND
};




/****************************************************************************/


/*
        macros for providing information for executable information system
 */
#ifdef C_FRONTEND
        #define INFO_FRONTEND "C_FRONTEND"
#endif
#ifdef CPP_FRONTEND
        #define INFO_FRONTEND "CPP_FRONTEND"
#endif
#ifdef F_FRONTEND
        #define INFO_FRONTEND "F_FRONTEND"
#endif


/*
        macros for controlling the executable information system
        this is based on usage off the RCS 'ident' command.
 */
/*
   #define DDD_RCS_STRING DDDRCSSTRING(DDD_VERSION,INFO_FRONTEND,MAX_OBJ,MAX_CPL)
   #define DDDRCSSTRING(ddd_version,info_frontend,max_obj,max_cpl)\
                DDDRCSSTRINGAUX(ddd_version,info_frontend,max_obj,max_cpl)
   #define DDDRCSSTRINGAUX(ddd_version,info_frontend,max_obj,max_cpl)\
                "$""State: DDD_VERSION="ddd_version" DDD_FRONTEND="info_frontend" DDD_MAXOBJ="#max_obj" DDD_MAXCPL="#max_cpl" $"

   #define RCSID(header,module_rcs_string) RCSIDAUX(header,module_rcs_string)
   #define RCSIDAUX(header,module_rcs_string) static char rcsid[] = #header ## module_rcs_string;
 */
#define RCSID(a,b)


/* insert the following line literally into all source files */
/* remove the question marks!! */
/* RCSID????("$Header$",DDD_RCS_STRING) */




/****************************************************************************/
/*                                                                          */
/* data structures exported by the corresponding source file                */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/* COUPLING: record for coupling of local object with foreign objects       */
/****************************************************************************/

typedef struct obj_coupl
{
  struct obj_coupl *_next;
  unsigned short proc;
  unsigned char prio;
  unsigned char flags;
  DDD_HDR obj;
} COUPLING;


#define CPL_NEXT(cpl)   ((cpl)->_next)



/****************************************************************************/
/* ELEM_DESC: description of one element in DDD object structure description */
/****************************************************************************/

/* TODO, in CPP_FRONTEND only one of the versions for C_FRONTEND and
        F_FRONTEND is needed, depending on setting of desc->storage.
        this should be a union for memory efficiency reasons.
 */
typedef struct _ELEM_DESC
{
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  int offset;                         /* element offset from object address     */
  unsigned char *gbits;               /* ptr to gbits array, if type==EL_GBITS  */
#endif

#if defined(F_FRONTEND) || defined(CPP_FRONTEND)
  char     *array;                        /* pointer to the array of this element   */
  int msgoffset;                                          /* offset of this element in the message  */
#endif

  size_t size;                        /* size of this element                   */
  int type;                           /* type of element, one of EL_xxx         */
  DDD_TYPE reftype;                   /* if EL_OBJPTR, type of ref. destination */
} ELEM_DESC;



/****************************************************************************/
/* TYPE_DESC: single DDD object structure description                       */
/****************************************************************************/

typedef struct _TYPE_DESC
{
  int mode;                             /* current TypeMode (DECLARE/DEFINE)    */
  char     *name;                       /* textual object description           */
  int currTypeDefCall;                  /* number of current call to TypeDefine */

#ifdef CPP_FRONTEND
  int storage;                          /* STORAGE_ARRAY or STORAGE_STRUCT      */

  /* if storage==STORAGE_ARRAY */
  int arraySize;                                        /* number of elements in the arrays     */
  //int    nextFree;				/* next free object in arrays           */
  int elemHeader;                       /* which rec. type-elem contains hdr?
                                           (offsetHeader gives local offset)    */
                                        /* (hasHeader must be TRUE)             */
#endif

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  /* if C_FRONTEND or (CPP_FRONTEND and storage==STORAGE_STRUCT) */
  int hasHeader;                        /* flag: real ddd type (with header)?   */
  int offsetHeader;                     /* offset of header from begin of obj   */
#endif

#ifdef F_FRONTEND
  int arraySize;                                                /* number of elements in the arrays     */
  int nextFree;                                                 /* next free object in arrays           */
  DDD_HDR hdr;                                          /* headers for all elements             */
#endif

  ELEM_DESC element[MAX_ELEMDESC];       /* element description array           */
  int nElements;                        /* number of elements in object         */
  size_t size;                          /* size of object, correctly aligned    */


  /* pointer to handler functions */
  HandlerLDATACONSTRUCTOR handlerLDATACONSTRUCTOR;
  HandlerDESTRUCTOR handlerDESTRUCTOR;
  HandlerDELETE handlerDELETE;
  HandlerUPDATE handlerUPDATE;
  HandlerOBJMKCONS handlerOBJMKCONS;
  HandlerSETPRIORITY handlerSETPRIORITY;
  HandlerXFERCOPY handlerXFERCOPY;
  HandlerXFERDELETE handlerXFERDELETE;
  HandlerXFERGATHER handlerXFERGATHER;
  HandlerXFERSCATTER handlerXFERSCATTER;
  HandlerXFERGATHERX handlerXFERGATHERX;
  HandlerXFERSCATTERX handlerXFERSCATTERX;
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
  HandlerXFERCOPYMANIP handlerXFERCOPYMANIP;
#endif
#ifdef F_FRONTEND
  HandlerALLOCOBJ handlerALLOCOBJ;
  HandlerFREEOBJ handlerFREEOBJ;
#endif


  DDD_PRIO *prioMatrix;                 /* 2D matrix for comparing priorities   */
  int prioDefault;                      /* default mode for PrioMerge           */

  /* redundancy for efficiency */
  int nPointers;                        /* number of outside references         */
  unsigned char  *cmask;                /* mask for fast type-dependent copy    */
} TYPE_DESC;



/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/

extern TYPE_DESC theTypeDefs[MAX_TYPEDESC];

extern DDD_HDR theObj[MAX_OBJ];
extern int nObjs;

extern COUPLING   *theCpl[MAX_CPL];
extern int theCplN[MAX_CPL];
extern int nCpls;                        /* number of coupling lists */
extern int nCplItems;                    /* number of couplings      */

extern int        *iBuffer;
extern char       *cBuffer;

extern int theOptions[OPT_END];


/* from topo.c */
extern VChannelPtr *theTopology;



/****************************************************************************/
/*                                                                          */
/* definitions previously hold in misc.h                                    */
/*                                                                          */
/****************************************************************************/

#ifndef ABS
#define ABS(i) (((i)<0) ? (-(i)) : (i))
#endif

#ifndef MIN
#define MIN(x,y) (((x)<(y)) ? (x) : (y))
#endif

#ifndef MAX
#define MAX(x,y) (((x)>(y)) ? (x) : (y))
#endif


/* round up to next alignment border */
#ifndef CEIL
#define CEIL(n) ((n)+((ALIGNMENT-((n)&(ALIGNMENT-1)))&(ALIGNMENT-1)))
#endif

/* round down to next alignment border */
#ifndef FLOOR
#define FLOOR(n) ((n)&ALIGNMASK)
#endif

#define YES     1
#define ON      1
#ifndef TRUE
#define TRUE    1
#endif

#define NO      0
#define OFF     0
#ifndef FALSE
#define FALSE   0
#endif




/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/

/* usage of flags in DDD_HEADER */
/* -- none yet -- */

/* usage of flags in COUPLING */
/* usage of 0x03 while interface-building, temporarily */
#define MASKCPLDIR 0x00000003
#define CPLDIR(c) (((INT)((c)->flags))&MASKCPLDIR)
#define SETCPLDIR(c,d) ((c)->flags) = (((c)->flags)&(~MASKCPLDIR))|((d)&MASKCPLDIR)



/* convert DDD_OBJ to DDD_HDR and vice versa */

#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
#define OBJ2HDR(obj,desc)  ((DDD_HDR)(((char *)obj)+((desc)->offsetHeader)))
#define HDR2OBJ(hdr,desc)  ((DDD_OBJ)(((char *)hdr)-((desc)->offsetHeader)))
#define OBJ_OBJ(hdr)       ((DDD_OBJ)(((char *)hdr)- \
                                      (theTypeDefs[OBJ_TYPE(hdr)].offsetHeader)))
#endif
#ifdef F_FRONTEND
#define OBJ2HDR(obj,desc) ((DDD_HDR)((desc)->hdr+(obj)))
#define HDR2OBJ(hd,desc)  ((DDD_OBJ)(((char *)(hd)-(char*)((desc)->hdr))/ \
                                     sizeof(DDD_HEADER)))
#define OBJ_OBJ(hd)       HDR2OBJ(hd,&theTypeDefs[OBJ_TYPE(hd)])
#endif

/* internal access of DDD_HEADER members */

/* type of object */
#define OBJ_TYPE(o)     ((o)->typ)

/* priority of object */
#define OBJ_PRIO(o)     ((o)->prio)

/* attr of object */
#define OBJ_ATTR(o)    ((o)->attr)

/* global id of object */
#define OBJ_GID(o)      ((o)->gid)

/* get index into global object table */
#define OBJ_INDEX(o)    ((o)->myIndex)

/* internal flags of object */
#define OBJ_FLAGS(o)    ((o)->flags)




/* get boolean: does object have couplings? */
#define HAS_COUPLING(o) (OBJ_INDEX(o)<nCpls)

/* get #couplings per object */
#define NCOUPLINGS(o) (HAS_COUPLING(o) ? theCplN[OBJ_INDEX(o)] : 0)

/* get pointer to object's coupling list */
#define THECOUPLING(o) (HAS_COUPLING(o) ? theCpl[OBJ_INDEX(o)] : NULL)


/* get default VChan to processor p */
#define VCHAN_TO(p)   (theTopology[(p)])


/* DDD_HDR may be invalid */
#define MarkHdrInvalid(hdr)    OBJ_INDEX(hdr)=MAX_OBJ
#define IsHdrInvalid(hdr)      OBJ_INDEX(hdr)==MAX_OBJ




#if defined(CheckPMEM) || defined(CheckIFMEM) || defined(CheckCplMEM) || defined(CheckMsgMEM) || defined(CheckTmpMEM)

static void *dummy_ptr;

#ifndef SST
#define SST (CEIL(sizeof(size_t)))
#define GET_SSTVAL(adr)   *(size_t *)(((char *)adr)-SST)
#endif

#endif


/* memory management */

/*** mapping memory allocation calls to memmgr_ calls ***/

#define AllocObj(s,t,p,a) memmgr_AllocOMEM((size_t)s,(int)t,(int)p,(int)a)
#define AllocHeap(s)      memmgr_AllocHMEM((size_t)s)
#define AllocCom(s)       memmgr_AllocAMEM((size_t)s)

#ifdef CheckPMEM
#define AllocFix(s)       (dummy_ptr=SST+(char *)memmgr_AllocPMEM(SST+(size_t)s));\
  GET_SSTVAL(dummy_ptr) = s;                                    \
  printf("%4d: MALL PFix adr=%08x size=%ld file=%s line=%d\n",\
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocFix(s)       memmgr_AllocPMEM((size_t)s)
#endif


#ifdef CheckMsgMEM
#define AllocMsg(s)       (dummy_ptr=SST+(char *)memmgr_AllocTMEM(SST+(size_t)s));\
  GET_SSTVAL(dummy_ptr) = s;                                    \
  printf("%4d: MALL TMsg adr=%08x size=%ld file=%s line=%d\n",\
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocMsg(s)       memmgr_AllocTMEM((size_t)s)
#endif


#ifdef CheckTmpMEM
#define AllocTmp(s)       (dummy_ptr=SST+(char *)memmgr_AllocTMEM(SST+(size_t)s));\
  GET_SSTVAL(dummy_ptr) = s;                                    \
  printf("%4d: MALL TTmp adr=%08x size=%ld file=%s line=%d\n",\
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocTmp(s)       memmgr_AllocTMEM((size_t)s)
#endif


#ifdef CheckCplMEM
#define AllocCpl(s)       (dummy_ptr=SST+(char *)memmgr_AllocAMEM(SST+(size_t)s));\
  GET_SSTVAL(dummy_ptr) = s;                                    \
  printf("%4d: MALL ACpl adr=%08x size=%ld file=%s line=%d\n",  \
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocCpl(s)       memmgr_AllocAMEM((size_t)s)
#endif

#ifdef CheckIFMEM
#define AllocIF(s)        (dummy_ptr=SST+(char *)memmgr_AllocAMEM(SST+(size_t)s));\
  GET_SSTVAL(dummy_ptr) = s;                                    \
  printf("%4d: MALL AIF  adr=%08x size=%ld file=%s line=%d\n",  \
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocIF(s)        memmgr_AllocAMEM((size_t)s)
#endif

#ifdef F_FRONTEND
#define AllocHdr(s)       memmgr_AllocAMEM((size_t)s)
#endif


/*** mapping memory free calls to memmgr calls ***/

#define FreeObj(mem,s,t)  memmgr_FreeOMEM(mem,(size_t)s,(int)t)
#define FreeHeap(mem)     memmgr_FreeHMEM(mem)
#define FreeCom(mem)      memmgr_FreeAMEM(mem)

#ifdef CheckPMEM
#define FreeFix(mem)      {               \
    size_t s=GET_SSTVAL(mem); \
    memmgr_FreePMEM(((char *)mem)-SST);  \
    printf("%4d: FREE PFix adr=%08x size=%ld file=%s line=%d\n",\
           me,mem,s,__FILE__,__LINE__); }
#else
#define FreeFix(mem)      memmgr_FreePMEM(mem)
#endif

#ifdef CheckMsgMEM
#define FreeMsg(mem)      {   \
    size_t s=GET_SSTVAL(mem); \
    memmgr_FreeTMEM(((char *)mem)-SST);    \
    printf("%4d: FREE TMsg adr=%08x size=%ld file=%s line=%d\n",\
           me,mem,s,__FILE__,__LINE__); }
#else
#define FreeMsg(mem)      memmgr_FreeTMEM(mem)
#endif


#ifdef CheckTmpMEM
#define FreeTmp(mem)      {                  \
    size_t s=GET_SSTVAL(mem); \
    memmgr_FreeTMEM(((char *)mem)-SST);       \
    printf("%4d: FREE TTmp adr=%08x size=%ld file=%s line=%d\n",\
           me,mem,s,__FILE__,__LINE__); }
#else
#define FreeTmp(mem)      memmgr_FreeTMEM(mem)
#endif


#ifdef CheckCplMEM
#define FreeCpl(mem)      {                   \
    size_t s=GET_SSTVAL(mem); \
    memmgr_FreeAMEM(((char *)mem)-SST);     \
    printf("%4d: FREE ACpl adr=%08x size=%ld file=%s line=%d\n",\
           me,mem,s,__FILE__,__LINE__); }
#else
#define FreeCpl(mem)      memmgr_FreeAMEM(mem)
#endif

#ifdef CheckIFMEM
#define FreeIF(mem)       { \
    size_t s=GET_SSTVAL(mem); \
    memmgr_FreeAMEM(((char *)mem)-SST);    \
    printf("%4d: FREE AIF  adr=%08x size=%ld file=%s line=%d\n",\
           me,mem,s,__FILE__,__LINE__); }
#else
#define FreeIF(mem)       memmgr_FreeAMEM(mem)
#endif

#ifdef F_FRONTEND
#define FreeHdr(mem)      memmgr_FreeAMEM(mem)
#endif


/* mapping mark/release heap calls to memmgr calls */
#define MarkHeap          memmgr_MarkHMEM
#define ReleaseHeap       memmgr_ReleaseHMEM



/****************************************************************************/

/* macros for mapping internal usage of external functions */
#if defined(C_FRONTEND) || defined(F_FRONTEND)
/* not used yet */
#endif


/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/

/* ddd.c */
int DDD_GetOption (DDD_OPTION);


/* typemgr.c */
#if defined(C_FRONTEND) || defined(F_FRONTEND)
void      ddd_TypeMgrInit (void);
#endif
void      ddd_TypeMgrExit (void);
int       ddd_TypeDefined (TYPE_DESC *);


/* cplmgr.c */
void      ddd_CplMgrInit (void);
void      ddd_CplMgrExit (void);
COUPLING *AddCoupling (DDD_HDR, DDD_PROC, DDD_PRIO);
COUPLING *ModCoupling (DDD_HDR, DDD_PROC, DDD_PRIO);
void      DelCoupling (DDD_HDR, DDD_PROC);
void      DisposeCouplingList (COUPLING *);
void      DDD_InfoCoupling (DDD_HDR);


/* prio.c */
int PriorityMerge (TYPE_DESC *, DDD_PRIO, DDD_PRIO, DDD_PRIO *);


/* if.c */
void      ddd_IFInit (void);
void      ddd_IFExit (void);
void      IFAllFromScratch (void);
void      DDD_InfoIFImpl (DDD_IF);
void      IFInvalidateShortcuts (DDD_TYPE);
int       DDD_CheckInterfaces (void);


/* xfer.c */
void      ddd_XferInit (void);
void      ddd_XferExit (void);
int       XferActive (void);
void      XferRegisterDelete (DDD_HDR);


/* ident.c */
void      ddd_IdentInit (void);
void      ddd_IdentExit (void);


/* cons.c */
void      ddd_ConsInit (void);
void      ddd_ConsExit (void);


/* cmds.c */
void      DDD_XferPrioChange (DDD_HDR, DDD_PRIO);


/* topo.c */
void      ddd_TopoInit (void);
void      ddd_TopoExit (void);
DDD_PROC *DDD_ProcArray (void);
void      DDD_GetChannels (int);
void      DDD_DisplayTopo (void);



/* objmgr.c */
void      DDD_HdrConstructorCopy (DDD_HDR, DDD_PRIO);
void      ObjCopyGlobalData (TYPE_DESC *, DDD_OBJ, DDD_OBJ, size_t);
DDD_HDR  *LocalObjectsList (void);
DDD_HDR  *LocalCoupledObjectsList (void);
#ifdef CPP_FRONTEND
DDD_OBJ  DDD_ObjNew (size_t, DDD_TYPE, DDD_PRIO, DDD_ATTR);
#endif
#ifdef F_FRONTEND
void      DDD_HdrDestructor (DDD_HDR);
#endif

/* reduct.c */
int       ddd_GlobalSumInt  (int);
int       ddd_GlobalMaxInt  (int);
int       ddd_GlobalMinInt  (int);

#endif
