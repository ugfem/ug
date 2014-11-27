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

/* RCS_ID
   $Header$
 */


/*
   #define CheckIFMEM
   #define CheckPMEM
   #define CheckCplMEM
   #define CheckMsgMEM
   #define CheckTmpMEM
   #define CheckHeapMem
 */


/****************************************************************************/
/*                                                                          */
/* auto include mechanism and other include files                           */
/*                                                                          */
/****************************************************************************/

#ifndef __DDDI_H__
#define __DDDI_H__

#include <climits>


#include <cassert>

#include "ugtypes.h"
#include "architecture.h"

#include "include/memmgr.h"
#include "ppif.h"
#include "ctrl/stat.h"

#include "dddstr.h"
#include "include/ddd.h"
#include "include/dddio.h"

START_UGDIM_NAMESPACE

/*
        macro for exiting program in case of a severe error condition
 */
#define HARD_EXIT  assert(0)
/* #define HARD_EXIT  exit(1) */


/*
        macros for correct return or premature abort of a procedure
 */
#define RET_ON_OK      return (TRUE)
#define RET_ON_ERROR   return (FALSE)
#define IS_OK(p)       ((p)==TRUE)
typedef int RETCODE;


/*
        macro for output of downward compatibility warning messages
 */
#define OLDSTYLE(txt) \
  if (me==master && DDD_GetOption(OPT_WARNING_OLDSTYLE)==OPT_ON)  \
  { DDD_PrintError('W', 1080, txt); }



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

#define MAX_ELEMDESC   64     /* max. number of elements per TYPE_DESC      */
#define MAX_TYPEDESC   32     /* max. number of TYPE_DESC                   */
#define MAX_PRIO       32     /* max. number of DDD_PRIO                    */
#define MAX_CPL_START  65536  /* max. number of local objects with coupling */

#ifdef WithFullObjectTable
#define MAX_OBJ_START  262144 /* max. number of locally registered objects  */
#else
#define MAX_OBJ_START  MAX_CPL_START
#endif



#define MAX_TRIES  50000000  /* max. number of tries until timeout in IF-comm */

#ifdef DDD_MAX_PROCBITS_IN_GID
#define MAX_PROCBITS_IN_GID DDD_MAX_PROCBITS_IN_GID
#else
#define MAX_PROCBITS_IN_GID 24  /* this allows 2^24 procs and 2^40 objects  */
#endif

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

/* string constants */
#define STR_NOMEM  "out of memory"


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
  unsigned short _proc;
  unsigned char prio;
  unsigned char _flags;
  DDD_HDR obj;
} COUPLING;


#define CPL_NEXT(cpl)   ((cpl)->_next)
#define CPL_PROC(cpl)   ((cpl)->_proc)



/****************************************************************************/
/* ELEM_DESC: description of one element in DDD object structure description */
/****************************************************************************/

/* TODO, in CPP_FRONTEND only one of the versions for C_FRONTEND and
        F_FRONTEND is needed, depending on setting of desc->storage.
        this should be a union for memory efficiency reasons.
 */
typedef struct _ELEM_DESC
{
  int offset;                         /* element offset from object address     */
  unsigned char *gbits;               /* ptr to gbits array, if type==EL_GBITS  */

#if defined(CPP_FRONTEND)
  char     *array;                        /* pointer to the array of this element   */
  int msgoffset;                                          /* offset of this element in the message  */
#endif

  size_t size;                        /* size of this element                   */
  int type;                           /* type of element, one of EL_xxx         */


  /* if type==EL_OBJPTR, the following entries will be used               */
  DDD_TYPE _reftype;                        /* DDD_TYPE of ref. destination     */

  /* if reftype==DDD_TYPE_BY_HANDLER, we must use a handler for
     determining the reftype on-the-fly (additional parameter during
     TypeDefine with EL_OBJPTR)                                        */
  HandlerGetRefType reftypeHandler;
} ELEM_DESC;


/* macros for accessing ELEM_DESC */
#define EDESC_REFTYPE(ed)          ((ed)->_reftype)
#define EDESC_SET_REFTYPE(ed,rt)   (ed)->_reftype=(rt)


/****************************************************************************/
/* TYPE_DESC: single DDD object structure description                       */
/****************************************************************************/

typedef struct _TYPE_DESC
{
  int mode;                             /* current TypeMode (DECLARE/DEFINE)    */
  const char  *name;                       /* textual object description           */
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

  /* if C_FRONTEND or (CPP_FRONTEND and storage==STORAGE_STRUCT) */
  int hasHeader;                        /* flag: real ddd type (with header)?   */
  int offsetHeader;                     /* offset of header from begin of obj   */


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
#if defined(C_FRONTEND)
  HandlerXFERCOPYMANIP handlerXFERCOPYMANIP;
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

extern DDD_HDR   *ddd_ObjTable;
extern int ddd_ObjTabSize;
extern int ddd_nObjs;

extern COUPLING   **ddd_CplTable;
extern short      *ddd_NCplTable;
extern int ddd_CplTabSize;
extern int ddd_nCpls;                    /* number of coupling lists */
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


/* internal access of DDD_HEADER members */

#ifdef CPP_FRONTEND
#define ACCESS_HDR(o,c)   ((o)->_hdr.c)
#else
#define ACCESS_HDR(o,c)   ((o)->c)
#endif


/* type of object */
#define OBJ_TYPE(o)     ACCESS_HDR(o,typ)

/* priority of object */
#define OBJ_PRIO(o)     ACCESS_HDR(o,prio)

/* attr of object */
#define OBJ_ATTR(o)     ACCESS_HDR(o,attr)

/* global id of object */
#define OBJ_GID(o)      ACCESS_HDR(o,gid)

/* get index into global object table */
#define OBJ_INDEX(o)    ACCESS_HDR(o,myIndex)

/* internal flags of object */
#define OBJ_FLAGS(o)    ACCESS_HDR(o,flags)


/****************************************************************************/

/* usage of flags in DDD_HEADER */
#define MASK_OBJ_PRUNED 0x00000001
#define OBJ_PRUNED(c) (((int)(OBJ_FLAGS(c)))&MASK_OBJ_PRUNED)
#define SET_OBJ_PRUNED(c,d) (OBJ_FLAGS(c)) = ((OBJ_FLAGS(c))&(~MASK_OBJ_PRUNED))|((d)&MASK_OBJ_PRUNED)

#define MASK_OBJ_RESENT  0x00000002
#define SHIFT_OBJ_RESENT 1
#define OBJ_RESENT(c) ((((int)(OBJ_FLAGS(c)))&MASK_OBJ_RESENT)>>SHIFT_OBJ_RESENT)
#define SET_OBJ_RESENT(c,d) (OBJ_FLAGS(c)) = ((OBJ_FLAGS(c))&(~MASK_OBJ_RESENT))|(((d)<<SHIFT_OBJ_RESENT)&MASK_OBJ_RESENT)


/* usage of flags in COUPLING */
/* usage of 0x03 while interface-building, temporarily */
#define MASKCPLDIR 0x00000003
#define CPLDIR(c) (((int)((c)->_flags))&MASKCPLDIR)
#define SETCPLDIR(c,d) ((c)->_flags) = (((c)->_flags)&(~MASKCPLDIR))|((d)&MASKCPLDIR)

/* usage of 0x10 for remembering the memory origin for the COUPLING struct */
#define MASKCPLMEM 0x00000010
#define CPLMEM_EXTERNAL  0x00
#define CPLMEM_FREELIST  0x10
#define CPLMEM(c) (((int)((c)->_flags))&MASKCPLMEM)
#define SETCPLMEM_EXTERNAL(c) ((c)->_flags) = (((c)->_flags)&(~MASKCPLMEM))|(CPLMEM_EXTERNAL)
#define SETCPLMEM_FREELIST(c) ((c)->_flags) = (((c)->_flags)&(~MASKCPLMEM))|(CPLMEM_FREELIST)



/* convert DDD_OBJ to DDD_HDR and vice versa */

#define OBJ2HDR(obj,desc)  ((DDD_HDR)(((char *)obj)+((desc)->offsetHeader)))
#define HDR2OBJ(hdr,desc)  ((DDD_OBJ)(((char *)hdr)-((desc)->offsetHeader)))
#define OBJ_OBJ(hdr)       ((DDD_OBJ)(((char *)hdr)- \
                                      (theTypeDefs[OBJ_TYPE(hdr)].offsetHeader)))


/****************************************************************************/

/*
        macros for access of coupling tables
 */

/* get boolean: does object have couplings? */
#define ObjHasCpl(o)      (OBJ_INDEX(o)<ddd_nCpls)

/* get #couplings per object */
#define ObjNCpl(o)        (ObjHasCpl(o) ? ddd_NCplTable[OBJ_INDEX(o)] : 0)
#define IdxNCpl(i)        (ddd_NCplTable[i])

/* get pointer to object's coupling list */
#define ObjCplList(o)     (ObjHasCpl(o) ? ddd_CplTable[OBJ_INDEX(o)] : NULL)
#define IdxCplList(i)     (ddd_CplTable[i])


/* increment/decrement number of coupled objects */
/*
   #define NCpl_Increment    { ddd_nCpls++; printf("%4d: nCpls++ now %d, %s:%d\n",me,ddd_nCpls,__FILE__,__LINE__); }
   #define NCpl_Decrement    { ddd_nCpls--; printf("%4d: nCpls-- now %d, %s:%d\n",me,ddd_nCpls,__FILE__,__LINE__); }
 */
#define NCpl_Increment    ddd_nCpls++;
#define NCpl_Decrement    ddd_nCpls--;

#define NCpl_Get          ddd_nCpls
#define NCpl_Init         ddd_nCpls=0


/* DDD_HDR may be invalid */
#define MarkHdrInvalid(hdr)    OBJ_INDEX(hdr)=(INT_MAX-1)
#define IsHdrInvalid(hdr)      OBJ_INDEX(hdr)==(INT_MAX-1)

#ifndef WithFullObjectTable
#define MarkHdrLocal(hdr)      OBJ_INDEX(hdr)=(INT_MAX)
#define IsHdrLocal(hdr)        OBJ_INDEX(hdr)==(INT_MAX)
#endif

/****************************************************************************/

/* get default VChan to processor p */
#define VCHAN_TO(p)   (theTopology[(p)])


/****************************************************************************/

/* types for StdIf-communication functions (see if/ifstd.ct) */
typedef int (*ExecProcHdrPtr)(DDD_HDR);
typedef int (*ExecProcHdrXPtr)(DDD_HDR, DDD_PROC, DDD_PRIO);
typedef int (*ComProcHdrPtr)(DDD_HDR, void *);
typedef int (*ComProcHdrXPtr)(DDD_HDR, void *, DDD_PROC, DDD_PRIO);


/****************************************************************************/
/*                                                                          */
/* specialties for CPP_FRONTEND                                             */
/*                                                                          */
/****************************************************************************/

#ifdef CPP_FRONTEND

/* not used up to now
   class DDD_ObjPtr
   {
        public:

                // access to DDD_HEADER of DDD_Object
                DDD_HDR operator-> ()  { return &(_obj->_hdr); }

                // access to DDD_Object itself
                DDD_Object* operator* ()  { return _obj; }

        private:
                DDD_Object*  _obj;
   };

   #define HdrPtr   DDD_ObjPtr
 */

#define CallHandler(desc,hname)     (desc->handler ## hname)
#define HParam(obj)                 obj,
#define HParamOnly(obj)             obj

#endif


#if defined(C_FRONTEND)
/*
   #define HdrPtr   DDD_HDR
 */

#endif




/****************************************************************************/
/*                                                                          */
/* memory management                                                        */
/*                                                                          */
/****************************************************************************/


#if defined(CheckPMEM) || defined(CheckIFMEM) || defined(CheckCplMEM) || defined(CheckMsgMEM) || defined(CheckTmpMEM) || defined(CheckHeapMem)

static void *dummy_ptr;
static char *mem_ptr;

#ifndef SST
#define SST (CEIL(sizeof(size_t)))
#define GET_SSTVAL(adr)   *(size_t *)(((char *)adr)-SST)
#endif

#endif


/*** mapping memory allocation calls to memmgr_ calls ***/

#define AllocObj(s,t,p,a) memmgr_AllocOMEM((size_t)s,(int)t,(int)p,(int)a)
#define AllocCom(s)       memmgr_AllocAMEM((size_t)s)


#ifdef CheckHeapMem
#define AllocHeap(s,key)  \
  (dummy_ptr = (mem_ptr=(char *)memmgr_AllocHMEM(SST+(size_t)(s), key)) != NULL ? \
               mem_ptr+SST : NULL);                                     \
  if (mem_ptr!=NULL) GET_SSTVAL(dummy_ptr) = s;                            \
  printf("%4d: MALL Heap adr=%08x size=%ld file=%s line=%d\n",             \
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocHeap(s,key)    memmgr_AllocHMEM((size_t)(s), key)
#endif


#ifdef CheckPMEM
#define AllocFix(s)  \
  (dummy_ptr = (mem_ptr=(char *)memmgr_AllocPMEM(SST+(size_t)s)) != NULL ? \
               mem_ptr+SST : NULL);                                     \
  if (mem_ptr!=NULL) GET_SSTVAL(dummy_ptr) = s;                            \
  printf("%4d: MALL PFix adr=%08x size=%ld file=%s line=%d\n",             \
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocFix(s)       memmgr_AllocPMEM((size_t)s)
#endif


#ifdef CheckMsgMEM
#define AllocMsg(s)  \
  (dummy_ptr = (mem_ptr=(char *)memmgr_AllocTMEM(SST+(size_t)s, TMEM_MSG)) \
               != NULL ?                                                            \
               mem_ptr+SST : NULL);                                     \
  if (mem_ptr!=NULL) GET_SSTVAL(dummy_ptr) = s;                            \
  printf("%4d: MALL TMsg adr=%08x size=%ld file=%s line=%d\n",             \
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocMsg(s)       memmgr_AllocTMEM((size_t)s, TMEM_MSG)
#endif


#ifdef CheckTmpMEM
#define AllocTmp(s)  \
  (dummy_ptr = (mem_ptr=(char *)memmgr_AllocTMEM(SST+(size_t)s, TMEM_ANY)) \
               != NULL ?                                                            \
               mem_ptr+SST : NULL);                                     \
  if (mem_ptr!=NULL) GET_SSTVAL(dummy_ptr) = s;                            \
  printf("%4d: MALL TTmp adr=%08x size=%ld file=%s line=%d\n",             \
         me,dummy_ptr,s,__FILE__,__LINE__)

#define AllocTmpReq(s,r)  \
  (dummy_ptr = (mem_ptr=(char *)memmgr_AllocTMEM(SST+(size_t)s, r))        \
               != NULL ?                                                            \
               mem_ptr+SST : NULL);                                     \
  if (mem_ptr!=NULL) GET_SSTVAL(dummy_ptr) = s;                            \
  printf("%4d: MALL TTmp adr=%08x size=%ld kind=%d file=%s line=%d\n",     \
         me,dummy_ptr,s,r,__FILE__,__LINE__)
#else
#define AllocTmp(s)       memmgr_AllocTMEM((size_t)s, TMEM_ANY)
#define AllocTmpReq(s,r)  memmgr_AllocTMEM((size_t)s, r)
#endif


#ifdef CheckCplMEM
#define AllocCpl(s)  \
  (dummy_ptr = (mem_ptr=(char *)memmgr_AllocAMEM(SST+(size_t)s)) != NULL ? \
               mem_ptr+SST : NULL);                                     \
  if (mem_ptr!=NULL) GET_SSTVAL(dummy_ptr) = s;                            \
  printf("%4d: MALL ACpl adr=%08x size=%ld file=%s line=%d\n",             \
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocCpl(s)       memmgr_AllocAMEM((size_t)s)
#endif

#ifdef CheckIFMEM
#define AllocIF(s)  \
  (dummy_ptr = (mem_ptr=(char *)memmgr_AllocAMEM(SST+(size_t)s)) != NULL ? \
               mem_ptr+SST : NULL);                                     \
  if (mem_ptr!=NULL) GET_SSTVAL(dummy_ptr) = s;                            \
  printf("%4d: MALL AIF  adr=%08x size=%ld file=%s line=%d\n",             \
         me,dummy_ptr,s,__FILE__,__LINE__)
#else
#define AllocIF(s)        memmgr_AllocAMEM((size_t)s)
#endif



/*** mapping memory free calls to memmgr calls ***/

#define FreeObj(mem,s,t)  memmgr_FreeOMEM(mem,(size_t)s,(int)t)
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
#define FreeMsg(mem,size)      {   \
    size_t s=GET_SSTVAL(mem); \
    memmgr_FreeTMEM(((char *)mem)-SST, TMEM_MSG);    \
    printf("%4d: FREE TMsg adr=%08x size=%ld file=%s line=%d\n",\
           me,mem,s,__FILE__,__LINE__); }
#else
#define FreeMsg(mem,size)      memmgr_FreeTMEM(mem, TMEM_MSG)
#endif


#ifdef CheckTmpMEM
#define FreeTmp(mem,size)      {                  \
    size_t s=GET_SSTVAL(mem); \
    memmgr_FreeTMEM(((char *)mem)-SST,TMEM_ANY);       \
    printf("%4d: FREE TTmp adr=%08x size=%ld file=%s line=%d\n",\
           me,mem,s,__FILE__,__LINE__); }
#define FreeTmpReq(mem,size,r)      {                  \
    size_t s=GET_SSTVAL(mem); \
    memmgr_FreeTMEM(((char *)mem)-SST,r);       \
    printf("%4d: FREE TTmp adr=%08x size=%ld kind=%d file=%s line=%d\n",\
           me,mem,s,r,__FILE__,__LINE__); }
#else
#define FreeTmp(mem,size)      memmgr_FreeTMEM(mem,TMEM_ANY)
#define FreeTmpReq(mem,size,r) memmgr_FreeTMEM(mem,r)
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



/* mapping mark/release heap calls to memmgr calls */
#define MarkHeap          memmgr_MarkHMEM
#define ReleaseHeap       memmgr_ReleaseHMEM



/****************************************************************************/

/* macros for mapping internal usage of external functions */
#if defined(C_FRONTEND)
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
#if defined(C_FRONTEND)
void      ddd_TypeMgrInit (void);
#endif
void      ddd_TypeMgrExit (void);
int       ddd_TypeDefined (TYPE_DESC *);


/* objmgr.c */
void      ddd_ObjMgrInit (void);
void      ddd_ObjMgrExit (void);
void      ddd_EnsureObjTabSize(int);


/* cplmgr.c */
void      ddd_CplMgrInit (void);
void      ddd_CplMgrExit (void);
COUPLING *AddCoupling (DDD_HDR, DDD_PROC, DDD_PRIO);
COUPLING *ModCoupling (DDD_HDR, DDD_PROC, DDD_PRIO);
void      DelCoupling (DDD_HDR, DDD_PROC);
void      DisposeCouplingList (COUPLING *);
void      DDD_InfoCoupling (DDD_HDR);


/* mgr/prio.c */
enum PrioMergeVals PriorityMerge (TYPE_DESC *, DDD_PRIO, DDD_PRIO, DDD_PRIO *);


/* if/if.c */
void      ddd_IFInit (void);
void      ddd_IFExit (void);
void      IFAllFromScratch (void);
void      DDD_InfoIFImpl (DDD_IF);
void      IFInvalidateShortcuts (DDD_TYPE);
int       DDD_CheckInterfaces (void);

/* if/ifcmds.c */
void   ddd_StdIFExchange   (size_t, ComProcHdrPtr,ComProcHdrPtr);
void   ddd_StdIFExecLocal  (        ExecProcHdrPtr);
void   ddd_StdIFExchangeX  (size_t, ComProcHdrXPtr,ComProcHdrXPtr);
void   ddd_StdIFExecLocalX (        ExecProcHdrXPtr);



/* xfer/xfer.c */
void      ddd_XferInit (void);
void      ddd_XferExit (void);
int       ddd_XferActive (void);
void      ddd_XferRegisterDelete (DDD_HDR);


/* xfer/cmds.c */
void      DDD_XferPrioChange (DDD_HDR, DDD_PRIO);


/* prio/pcmds.c */
void      ddd_PrioInit (void);
void      ddd_PrioExit (void);
int       ddd_PrioActive (void);


/* join/join.c */
void      ddd_JoinInit (void);
void      ddd_JoinExit (void);
int       ddd_JoinActive (void);


/* ident/ident.c */
void      ddd_IdentInit (void);
void      ddd_IdentExit (void);


/* basic/cons.c */
void      ddd_ConsInit (void);
void      ddd_ConsExit (void);


/* basic/topo.c */
void      ddd_TopoInit (void);
void      ddd_TopoExit (void);
DDD_PROC *DDD_ProcArray (void);
RETCODE   DDD_GetChannels (int);
void      DDD_DisplayTopo (void);



/* mgr/objmgr.c */
void      DDD_HdrConstructorCopy (DDD_HDR, DDD_PRIO);
void      ObjCopyGlobalData (TYPE_DESC *, DDD_OBJ, DDD_OBJ, size_t);
DDD_HDR  *LocalObjectsList (void);
void      FreeLocalObjectsList (DDD_HDR *);
DDD_HDR  *LocalCoupledObjectsList (void);
void      FreeLocalCoupledObjectsList (DDD_HDR *);
#ifdef CPP_FRONTEND
DDD_OBJ  DDD_ObjNew (size_t, DDD_TYPE, DDD_PRIO, DDD_ATTR);
#endif

/* basic/reduct.c */
int       ddd_GlobalSumInt  (int);
int       ddd_GlobalMaxInt  (int);
int       ddd_GlobalMinInt  (int);


/* ctrl/stat.c */
void      ddd_StatInit (void);
void      ddd_StatExit (void);


END_UGDIM_NAMESPACE

#endif
