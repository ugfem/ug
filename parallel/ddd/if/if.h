// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      if.h                                                          */
/*                                                                          */
/* Purpose:   include file for DDD interface module                         */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/09/21 kb  extracted from if.c                              */
/*            96/01/16 kb  added DDD_OBJ shortcut to ifHead and ifAttr      */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/


/*
   #define CtrlTimeouts
   #define CtrlTimeoutsDetailed
 */


/* maximum number of interfaces */
#define MAX_IF                  20


/* maximum length of interface description string */
#define IF_NAMELEN      80



enum CplDir {
  DirAB  =  0x01,
  DirBA  =  0x02,
  DirABA =  DirAB|DirBA
};





/* macros for easier coding of replicated source code */

/* loop over one interface (all ifHeads) */
#define ForIF(id,iter)  for((iter)=theIF[(id)].ifHead;  \
                            (iter)!=NULL;  \
                            (iter)=(iter)->next)


#ifdef CPP_FRONTEND
typedef DDD_Object* IFObjPtr;
#else
typedef DDD_OBJ IFObjPtr;
#endif


/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/* IF_ATTR: single part of interface, all couplings have same attr          */
/****************************************************************************/

typedef struct if_attr
{
  struct if_attr *next;

  /* note: the cplXX resp. objXX arrays are NOT contiguous in memory */
  COUPLING   **cplAB, **cplBA, **cplABA;
  IFObjPtr   *objAB,  *objBA,  *objABA;       /* object shortcut */
  int nItems, nAB, nBA, nABA;
  DDD_ATTR attr;

} IF_ATTR;



/****************************************************************************/
/* IF_PROC: descriptor of message and its contents/buffers for IF-communic. */
/****************************************************************************/

typedef struct if_proc
{
  struct     if_proc *next;
  IF_ATTR    *ifAttr;
  int nAttrs;

  /* note: the cplXX resp. objXX arrays ARE contiguous in memory */
  COUPLING   **cpl, **cplAB, **cplBA, **cplABA;
  IFObjPtr   *obj,  *objAB,  *objBA,  *objABA;       /* object shortcut */
  int nItems, nAB, nBA, nABA;
  DDD_PROC proc;

  VChannelPtr vc;
  msgid msgIn;
  msgid msgOut;
  char       *msgBufIn;
  char       *msgBufOut;
  int lenBufIn;
  int lenBufOut;

} IF_PROC;



/****************************************************************************/
/* IF_DEF: descriptor for one single interface                              */
/****************************************************************************/

typedef struct if_def
{
  IF_PROC   *ifHead;
  COUPLING  **cpl;              /* list of couplings belonging to interface     */
  int nItems;                   /* overall number of items in this interface    */

  IFObjPtr  *obj;              /* shortcut: list of object addresses in interf */
  int objValid;                 /* flag: is obj-table valid?                    */

  int nIfHeads;

  int nObjStruct;
  int nPrioA;
  int nPrioB;
  DDD_TYPE O[16];
  DDD_PRIO A[16];
  DDD_PRIO B[16];

  /* data for efficiency tuning */
  int maskO;

  /* data for nice user interaction */
  char name[IF_NAMELEN+1];            /* string for interface identification */

} IF_DEF;



/****************************************************************************/
/*                                                                          */
/* definition of variables for corresponding module                         */
/*                                                                          */
/****************************************************************************/

extern IF_DEF theIF[MAX_IF];
extern int nIFs;



/****************************************************************************/
/*                                                                          */
/* function declarations                                                    */
/*                                                                          */
/****************************************************************************/


/* ifuse.c */
void    IFGetMem (IF_PROC *, size_t, int, int);
int     IFInitComm (DDD_IF);
void    IFExitComm (DDD_IF);
void    IFInitSend (IF_PROC *);
int     IFPollSend (DDD_IF);
#ifdef CPP_FRONTEND
char*  IFCommLoopObjGather (DDD_GatherScatter&, IFObjPtr*, char*, size_t, int);
char*  IFCommLoopObjScatter (DDD_GatherScatter&, IFObjPtr*, char*, size_t, int);
#else
char *  IFCommLoopObj (ComProcPtr, IFObjPtr *, char *, size_t, int);
#endif
#ifndef CPP_FRONTEND
char *  IFCommLoopCpl (ComProcPtr, COUPLING **, char *, size_t, int);
char *  IFCommLoopCplX (ComProcXPtr, COUPLING **, char *, size_t , int);
#endif
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
#ifndef CPP_FRONTEND
void    IFExecLoopObj (ExecProcPtr, IFObjPtr *, int);
void    IFExecLoopCplX (ExecProcXPtr, COUPLING **, int);
#endif
#endif


/* ifobjsc.c */
void IFCreateObjShortcut (DDD_IF);
void IFCheckShortcuts (DDD_IF);
