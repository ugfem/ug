// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      cplmgr.c                                                      */
/*                                                                          */
/* Purpose:   management of couplings                                       */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/30 kb  begin                                            */
/*            94/08/24 kb  added DDD_InfoProcList()                         */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* include files                                                            */
/*            system include files                                          */
/*            application include files                                     */
/*                                                                          */
/****************************************************************************/

/* standard C library */
#include <config.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>

#include "dddi.h"

USING_UG_NAMESPACES

/* PPIF namespace: */
  USING_PPIF_NAMESPACE

#define DebugCoupling 10  /* 10 is off */

START_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/

/* size of segment of couplings (for memory allocation) */
#define CPLSEGM_SIZE 512



/****************************************************************************/
/*                                                                          */
/* data types                                                               */
/*                                                                          */
/****************************************************************************/

/*
        the storage of COUPLING items is done with the following scheme:
        allocation in segments of couplings, freeing into freelist.

        ALLOC:  try first to get one item out of freelist (memlistCpl),
                if that's not possible, get one from current segment;
                alloc segments from MemMgr.

        FREE:   put coupling into freelist.
 */


/* segment of Cpls */
typedef struct _CplSegm
{
  struct _CplSegm *next;
  int nItems;

  COUPLING item[CPLSEGM_SIZE];
} CplSegm;



/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)



static CplSegm *segmCpl = NULL;
static COUPLING *memlistCpl = NULL;
static int *localIBuffer;
static int nCplSegms;


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


static CplSegm *NewCplSegm (void)
{
  CplSegm *segm;

  segm = (CplSegm *) AllocTmpReq(sizeof(CplSegm), TMEM_CPL);

  if (segm==NULL)
  {
    DDD_PrintError('F', 2550, STR_NOMEM " during NewCoupling()");
    HARD_EXIT;
  }

  segm->next   = segmCpl;
  segmCpl      = segm;
  segm->nItems = 0;
  nCplSegms++;

  return(segm);
}


static void FreeCplSegms (void)
{
  CplSegm *segm = segmCpl;
  CplSegm *next = NULL;

  while (segm!=NULL)
  {
    next = segm->next;
    FreeTmpReq(segm, sizeof(CplSegm), TMEM_CPL);

    segm = next;
  }

  segmCpl = NULL;
  nCplSegms = 0;
  memlistCpl = NULL;
}


/****************************************************************************/

/* auxiliary function to init Coupling memory and initial data */
static void InitNewCoupling (COUPLING* cpl)
{
  /* set coupling to initial value, in order to find any bugs lateron.
     (hint by C.Wieners) */
  memset(cpl, 0, sizeof(COUPLING));

  /* init private data */
  cpl->_flags = 0;
}


static COUPLING *NewCoupling (void)
{
  COUPLING *cpl;

  if (DDD_GetOption(OPT_CPLMGR_USE_FREELIST)==OPT_ON)
  {
    /* allocate coupling from segments (which are allocated
       from segment-freelists) */

    if (memlistCpl==NULL)
    {
      CplSegm *segm = segmCpl;

      if (segm==NULL || segm->nItems==CPLSEGM_SIZE)
      {
        segm = NewCplSegm();
      }

      cpl = &(segm->item[segm->nItems++]);
    }
    else
    {
      cpl = memlistCpl;
      memlistCpl = CPL_NEXT(cpl);
    }

    /* init coupling memory and its private data */
    InitNewCoupling(cpl);

    /* remember memory origin for later disposal */
    SETCPLMEM_FREELIST(cpl);
  }
  else
  {
    /* allocate coupling directly */
    cpl = (COUPLING *) AllocTmpReq(sizeof(COUPLING), TMEM_CPL);

    if (cpl==NULL)
    {
      DDD_PrintError('F', 2551, STR_NOMEM " during NewCoupling()");
      HARD_EXIT;
    }

    /* init coupling memory and its private data */
    InitNewCoupling(cpl);

    /* remember memory origin for later disposal */
    SETCPLMEM_EXTERNAL(cpl);
  }

  nCplItems++;

  return(cpl);
}


static void DisposeCoupling (COUPLING *cpl)
{
  if (CPLMEM(cpl)==CPLMEM_FREELIST)
  {
    CPL_NEXT(cpl) = memlistCpl;
    memlistCpl = cpl;
  }
  else
  {
    FreeTmpReq(cpl, sizeof(COUPLING), TMEM_CPL);
  }

  nCplItems--;
}


/****************************************************************************/


static void AllocCplTables (long n)
{
  /* allocate coupling table */
  ddd_CplTable = (COUPLING **) AllocTmp(sizeof(COUPLING *) * n);
  if (ddd_CplTable==NULL)
  {
    sprintf(cBuffer, STR_NOMEM " for coupling table of size %ld",
            n * sizeof(COUPLING *));
    DDD_PrintError('E', 2510, cBuffer);
    HARD_EXIT;
  }

  ddd_NCplTable = (short *) AllocTmp(sizeof(short) * n);
  if (ddd_NCplTable==NULL)
  {
    sprintf(cBuffer, STR_NOMEM " for cpl-sizes table of size %ld",
            n * (long)sizeof(short));
    DDD_PrintError('E', 2511, cBuffer);
    HARD_EXIT;
  }

  ddd_CplTabSize = n;
}


static void IncreaseCplTabSize (void)
{
  COUPLING **old_CplTable   = ddd_CplTable;
  short     *old_NCplTable  = ddd_NCplTable;
  int old_CplTabSize = ddd_CplTabSize;

  /* compute new size (currently: double size) */
  ddd_CplTabSize = old_CplTabSize * 2;

  /* allocate new coupling table */
  ddd_CplTable = (COUPLING **) AllocTmp(sizeof(COUPLING *) * ddd_CplTabSize);
  if (ddd_CplTable==NULL)
  {
    sprintf(cBuffer, STR_NOMEM " for coupling table of size %ld",
            ((long)ddd_CplTabSize) * sizeof(COUPLING *));
    DDD_PrintError('W', 2512, cBuffer);

    /* restore data and return without action */
    ddd_CplTabSize = old_CplTabSize;
    ddd_CplTable   = old_CplTable;
    return;
  }

  /* copy data from old cpl-table to new one, assuming the old one is full */
  memcpy(ddd_CplTable, old_CplTable, sizeof(COUPLING *) * old_CplTabSize);

  /* free old one */
  FreeTmp(old_CplTable,0);


  /* now we alloc new ncpl-table. the number of entries in ddd_CplTable and
     ddd_NCplTable are equal, but the size per entry (in byte) for ddd_NCplTable
     is half. therefore, the following allocation will get exactly the memory
     chunk the previous FreeTmp freed (if only the next layer can manage this
     adequately, i.e. the MemMgr). */
  ddd_NCplTable = (short *) AllocTmp(sizeof(short) * ddd_CplTabSize);
  if (ddd_NCplTable==NULL)
  {
    sprintf(cBuffer, STR_NOMEM " for cpl-sizes table of size %ld",
            ((long)ddd_CplTabSize) * sizeof(short));
    DDD_PrintError('E', 2513, cBuffer);
    HARD_EXIT;
  }

  /* again, copy data */
  memcpy(ddd_NCplTable, old_NCplTable, sizeof(short) * old_CplTabSize);

  /* again free old table */
  FreeTmp(old_NCplTable,0);


  /* issue a warning in order to inform user */
  sprintf(cBuffer, "increased coupling table, now %d entries", ddd_CplTabSize);
  DDD_PrintError('W', 2514, cBuffer);


  ddd_EnsureObjTabSize(ddd_CplTabSize);
}




/****************************************************************************/
/*                                                                          */
/* Function:  AddCoupling                                                   */
/*                                                                          */
/* Purpose:   get new coupling record and init contents                     */
/*            if coupling is already there, no additional coupling is       */
/*            created. priority is adapted in this case.                    */
/*                                                                          */
/* Input:     hdr: DDD-header of object with new coupling                   */
/*            proc: owner of copy to be registered                          */
/*            prio: priority of copy                                        */
/*                                                                          */
/* Output:    ptr to new cpl record (or old one, if existed before)         */
/*            NULL on error                                                 */
/*                                                                          */
/****************************************************************************/

COUPLING *AddCoupling (DDD_HDR hdr, DDD_PROC proc, DDD_PRIO prio)
{
  COUPLING        *cp, *cp2;
  DDD_HDR oldObj;
  int objIndex;
  int freeCplIdx = NCpl_Get;

  assert(proc!=me);

#       if DebugCoupling<=1
  sprintf(cBuffer, "%4d: AddCoupling %08x proc=%d prio=%d\n",
          me, OBJ_GID(hdr), proc, prio);
  DDD_PrintDebug(cBuffer);
#       endif

  /* find or free position in coupling array */
  objIndex = OBJ_INDEX(hdr);
  if (! ObjHasCpl(hdr))
  {
    if (freeCplIdx==ddd_CplTabSize)
    {
      /* try to make CplTables larger ... */
      IncreaseCplTabSize();

      if (freeCplIdx==ddd_CplTabSize)
      {
        /* didn't work, give up. */
        DDD_PrintError('E', 2520, "no more couplings in AddCoupling");
        HARD_EXIT;
        /*return(NULL);*/
      }
    }

                #ifdef WithFullObjectTable
    oldObj = ddd_ObjTable[freeCplIdx];

    /* exchange object without coupling and object with coupling */
    /* free position freeCplIdx, move corresponding hdr reference
       elsewhere. */
    ddd_ObjTable[objIndex] = oldObj;
    OBJ_INDEX(oldObj)      = objIndex;
                #else
    assert(IsHdrLocal(hdr));

    /* hdr has been local, therefore not known by DDD, we have
       to register it now. */
    ddd_nObjs++;
                #endif


    assert(freeCplIdx<ddd_ObjTabSize);
    ddd_ObjTable[freeCplIdx] = hdr;
    OBJ_INDEX(hdr)           = freeCplIdx;

    objIndex = freeCplIdx;
    IdxCplList(objIndex) = NULL;
    IdxNCpl(objIndex) = 0;

    NCpl_Increment;
  }
  else
  {
    for(cp2=IdxCplList(objIndex); cp2!=NULL; cp2=CPL_NEXT(cp2))
    {
      if (CPL_PROC(cp2)==proc)
      {
        if (cp2->prio!=prio)
        {
          /* coupling upgrades/downgrades, are they allowed?
                                                  printf("%4d: diff in cpl, %05x old %d-%d new %d-%d\n",
                                                          me,OBJ_GID(hdr),cp2->proc,cp2->prio, proc, prio);
           */
          cp2->prio = prio;
        }
        /*
                                        DDD_PrintError('W', 2600, "coupling already known in AddCoupling");
         */
        return(cp2);
      }
    }
  }

  /* create new coupling record */
  cp = NewCoupling();
  if (cp==NULL) {
    DDD_PrintError('E', 2500, STR_NOMEM " in AddCoupling");
    return(NULL);
  }

  /* init contents */
  cp->obj = hdr;
  CPL_PROC(cp) = proc;
  cp->prio = prio;

  /* insert into theCpl array */
  CPL_NEXT(cp) = IdxCplList(objIndex);
  IdxCplList(objIndex) = cp;
  IdxNCpl(objIndex)++;

  return(cp);
}





/****************************************************************************/
/*                                                                          */
/* Function:  ModCoupling                                                   */
/*                                                                          */
/* Purpose:   find existing coupling record and modify priority             */
/*            this function does coupling upgrade/downgrade without         */
/*            complaining.                                                  */
/*                                                                          */
/* Input:     hdr: DDD-header of object with new coupling                   */
/*            proc: owner of copy to be modified                            */
/*            prio: new priority of copy                                    */
/*                                                                          */
/* Output:    ptr to old cpl record                                         */
/*            NULL on error                                                 */
/*                                                                          */
/****************************************************************************/

COUPLING *ModCoupling (DDD_HDR hdr, DDD_PROC proc, DDD_PRIO prio)
{
  COUPLING        *cp2;
  int objIndex;

  assert(proc!=me);

#       if DebugCoupling<=1
  sprintf(cBuffer, "%4d: ModCoupling %08x proc=%d prio=%d\n",
          me, OBJ_GID(hdr), proc, prio);
  DDD_PrintDebug(cBuffer);
#       endif

  /* find or free position in coupling array */
  objIndex = OBJ_INDEX(hdr);
  if (! ObjHasCpl(hdr))
  {
    /* there are no couplings for this object! */
    sprintf(cBuffer, "no couplings for %08x in ModCoupling", OBJ_GID(hdr));
    DDD_PrintError('E', 2530, cBuffer);
    return(NULL);
  }
  else
  {
    /* look if coupling exists and change it */
    for(cp2=IdxCplList(objIndex); cp2!=NULL; cp2=CPL_NEXT(cp2))
    {
      if (CPL_PROC(cp2)==proc)
      {
        cp2->prio = prio;
        return(cp2);
      }
    }
  }

  /* coupling not found */
  sprintf(cBuffer, "no coupling from %d for %08x in ModCoupling",
          proc, OBJ_GID(hdr));
  DDD_PrintError('E', 2531, cBuffer);
  HARD_EXIT;

  return(NULL);         /* never reach this */
}




/****************************************************************************/
/*                                                                          */
/* Function:  DelCoupling                                                   */
/*                                                                          */
/* Purpose:   remove coupling record from object                            */
/*                                                                          */
/* Input:     hdr: DDD-header of object with old coupling                   */
/*            proc: owner of copy to be un-registered                       */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void DelCoupling (DDD_HDR hdr, DDD_PROC proc)
{
  COUPLING        *cpl, *cplLast;
  int objIndex;

  objIndex = OBJ_INDEX(hdr);

  if (objIndex<NCpl_Get)
  {
    for(cpl=IdxCplList(objIndex), cplLast=NULL; cpl!=NULL; cpl=CPL_NEXT(cpl))
    {
      if(CPL_PROC(cpl)==proc)
      {
        if (cplLast==NULL)
        {
          IdxCplList(objIndex) = CPL_NEXT(cpl);
        }
        else {
          CPL_NEXT(cplLast) = CPL_NEXT(cpl);
        }
#                               if DebugCoupling<=1
        sprintf(cBuffer,"%4d: DelCoupling %07x on proc=%d, now %d cpls\n",
                me, OBJ_GID(hdr), proc, IdxNCpl(objIndex)-1);
        DDD_PrintDebug(cBuffer);
#                               endif

        DisposeCoupling(cpl);

        IdxNCpl(objIndex)--;

        if (IdxNCpl(objIndex)==0)
        {
          NCpl_Decrement;

                                        #ifdef WithFullObjectTable
          OBJ_INDEX(hdr) = NCpl_Get;
          OBJ_INDEX(ddd_ObjTable[NCpl_Get]) = objIndex;
          ddd_ObjTable[objIndex] = ddd_ObjTable[NCpl_Get];
          ddd_ObjTable[NCpl_Get] = hdr;
                                        #else
          /* we will not register objects without coupling,
             so we have to forget about hdr and mark it as local. */
          ddd_nObjs--; assert(ddd_nObjs==NCpl_Get);

          ddd_ObjTable[objIndex] = ddd_ObjTable[NCpl_Get];
          OBJ_INDEX(ddd_ObjTable[NCpl_Get]) = objIndex;

          MarkHdrLocal(hdr);
                                        #endif

          IdxCplList(objIndex) = IdxCplList(NCpl_Get);
          IdxNCpl(objIndex) = IdxNCpl(NCpl_Get);
        }
        break;
      }
      cplLast = cpl;
    }
  }
}


/****************************************************************************/
/*                                                                          */
/* Function:  DisposeCouplingList                                           */
/*                                                                          */
/* Purpose:   dispose complete coupling list                                */
/*                                                                          */
/* Input:     cpl: first element of coupling list                           */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void DisposeCouplingList (COUPLING *cpl)
{
  COUPLING *c, *next;

  c = cpl;
  while (c!=NULL)
  {
    next = CPL_NEXT(c);
    DisposeCoupling(c);
    c = next;
  }
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_InfoProcList                                              */
/*                                                                          */
/* Purpose:   return list of couplings of certain object                    */
/*                                                                          */
/* Input:     hdr: DDD-header of object with coupling                       */
/*                                                                          */
/* Output:    pointer to localIBuffer, which has been filled with:          */
/*               1) id of calling processor                                 */
/*               2) priority of local object copy on calling processor      */
/*               3) id of processor which holds an object copy              */
/*               4) priority of copy on that processor                      */
/*               5) 3+4 repeated for each coupling                          */
/*               6) processor number = -1 as end mark                       */
/*                                                                          */
/****************************************************************************/

#ifdef C_FRONTEND
int *DDD_InfoProcList (DDD_HDR hdr)
{
#endif
#ifdef CPP_FRONTEND
int *DDD_Object::InfoProcList (void)
{
  DDD_HDR hdr = this;
#endif

COUPLING *cpl;
int i, objIndex = OBJ_INDEX(hdr);

/* insert description of own (i.e. local) copy */
localIBuffer[0] = me;
localIBuffer[1] = OBJ_PRIO(hdr);

i=2;

/* append descriptions of foreign copies */
if (objIndex<NCpl_Get)
{
  for(cpl=IdxCplList(objIndex); cpl!=NULL; cpl=CPL_NEXT(cpl), i+=2) {
    localIBuffer[i]   = CPL_PROC(cpl);
    localIBuffer[i+1] = cpl->prio;
  }
}

/* append end mark */
localIBuffer[i] = -1;

return(localIBuffer);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_InfoProcPrio                                              */
/*                                                                          */
/* Purpose:   return first processor number with a given priority           */
/*                                                                          */
/* Input:     hdr:  DDD-header of object with coupling                      */
/*            prio: priority to search for                                  */
/*                                                                          */
/* Output:    id of processor which holds the object copy with prio         */
/*            (or procs if no such copy exists)                             */
/*                                                                          */
/****************************************************************************/

DDD_PROC DDD_InfoProcPrio (DDD_HDR hdr, DDD_PRIO prio)
{
  COUPLING *cpl;
  int objIndex = OBJ_INDEX(hdr);

  /* append descriptions of foreign copies */
  if (objIndex<NCpl_Get)
  {
    for(cpl=IdxCplList(objIndex); cpl!=NULL; cpl=CPL_NEXT(cpl))
    {
      if (cpl->prio == prio)
        return(CPL_PROC(cpl));
    }
  }

  /* eventually local copy has priority we are looking for */
  if (OBJ_PRIO(hdr)==prio)
    return(me);

  return(procs);
}


int DDD_InfoIsLocal (DDD_HDR hdr)
{
  return(! ObjHasCpl(hdr));
}


int DDD_InfoNCopies (DDD_HDR hdr)
{
  /*
     COUPLING *cpl;
     int n = 0;

     if (ObjHasCpl(hdr))
     {
          for(cpl=IdxCplList(OBJ_INDEX(hdr)); cpl!=NULL; cpl=CPL_NEXT(cpl))
                  n++;
     }
   */

  return(ObjNCpl(hdr));
}


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_InfoCoupling                                              */
/*                                                                          */
/* Purpose:   displays list of coupling for certain object                  */
/*                                                                          */
/* Input:     hdr: DDD-header of object with coupling                       */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void DDD_InfoCoupling (DDD_HDR hdr)
{
  COUPLING *cpl;
  int objIndex = OBJ_INDEX(hdr);

  sprintf(cBuffer, "%4d: InfoCoupling for object %07x (%05d/%05d)\n",
          me, OBJ_GID(hdr), objIndex, NCpl_Get);
  DDD_PrintLine(cBuffer);

  if (objIndex<NCpl_Get)
  {
    for(cpl=IdxCplList(objIndex); cpl!=NULL; cpl=CPL_NEXT(cpl))
    {
      sprintf(cBuffer, "%4d:    cpl %08x proc=%4d prio=%4d\n",
              me, cpl, CPL_PROC(cpl), cpl->prio);
      DDD_PrintLine(cBuffer);
    }
  }
}




/****************************************************************************/
/*                                                                          */
/* Function:  DDD_InfoCplMemory                                             */
/*                                                                          */
/* Purpose:   returns number of bytes used for coupling data                */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    size of memory used for couplings                             */
/*                                                                          */
/****************************************************************************/

size_t DDD_InfoCplMemory (void)
{
  size_t sum = 0;

  sum += sizeof(CplSegm) * nCplSegms;

  return(sum);
}



/****************************************************************************/
/*                                                                          */
/* Function:  CplMgrInit and CplMgrExit                                     */
/*                                                                          */
/* Purpose:   init/exit coupling manager                                    */
/*                                                                          */
/* Input:     -                                                             */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void ddd_CplMgrInit (void)
{
  /* allocate first (smallest) coupling tables */
  AllocCplTables(MAX_CPL_START);


  localIBuffer = (int*)AllocFix((2*procs+1)*sizeof(int));
  if (localIBuffer==NULL)
  {
    DDD_PrintError('E', 2532, STR_NOMEM " for DDD_InfoProcList()");
    HARD_EXIT;
  }

  memlistCpl = NULL;
  segmCpl    = NULL;
  nCplSegms  = 0;
}


void ddd_CplMgrExit (void)
{
  FreeFix(localIBuffer);
  FreeCplSegms();

  FreeTmp(ddd_CplTable,0);
  FreeTmp(ddd_NCplTable,0);
}

/****************************************************************************/

END_UGDIM_NAMESPACE
