// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      unpack.c                                                      */
/*                                                                          */
/* Purpose:   receives and unpacks messages                                 */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   940201 kb  begin                                              */
/*            960508 kb  restructured completely. efficiency improvement.   */
/*            960718 kb  introduced lowcomm-layer (sets of messages)        */
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
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "dddi.h"
#include "xfer.h"



/* TODO kb 961210
   #define WANTED_NOCH_GENAUER_UNTERSUCHEN
 */

/* TODO kb 961210
   #define DEBUG_MERGE_MODE
 */
#define MERGE_MODE_IN_TESTZUSTAND


#define DebugUnpack  5  /* off: 5 */

/*#define DebugCouplingCons*/


/*
   #define AddCoupling(a,b,c)  printf("%4d: AC %d,%d/%d %d\n",me,__LINE__,b,c,(int) AddCoupling(a,b,c))
 */


/****************************************************************************/
/*                                                                          */
/* constant definitions                                                     */
/*                                                                          */
/****************************************************************************/


#ifdef C_FRONTEND
#       define SIZEOF_REF  sizeof(void *)
#       define NULL_REF    NULL
#else /* F_FRONTEND */
#       define SIZEOF_REF  sizeof(DDD_OBJ)
#       define NULL_REF    0
#endif



/****************************************************************************/
/*                                                                          */
/* definition of exported global variables                                  */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)




/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


static int sort_TENewCpl (const void *e1, const void *e2)
{
  TENewCpl   *ci1, *ci2;

  ci1 = (TENewCpl *)e1;
  ci2 = (TENewCpl *)e2;

  if (ci1->gid < ci2->gid) return(-1);
  if (ci1->gid > ci2->gid) return(1);

  if (ci1->dest < ci2->dest) return(-1);
  if (ci1->dest > ci2->dest) return(1);

  /* ascending priority */
  if (ci1->prio < ci2->prio) return(-1);
  if (ci1->prio > ci2->prio) return(1);

  return(0);
}


static int sort_ObjTabPtrs (const void *e1, const void *e2)
{
  OBJTAB_ENTRY   *ci1, *ci2;

  ci1 = *(OBJTAB_ENTRY **)e1;
  ci2 = *(OBJTAB_ENTRY **)e2;

  /* sort with ascending gid */
  if (ci1->gid < ci2->gid) return(-1);
  if (ci1->gid > ci2->gid) return(1);

  /* sort with decreasing priority */
  /* not necessary anymore. see first phase of
     AcceptReceivedObjects() for details. KB 970128
     if (ci1->prio < ci2->prio) return(1);
     if (ci1->prio > ci2->prio) return(-1);
   */

  return(0);
}



/****************************************************************************/


/*
        convert indices to symtab into references (pointers).
        the object objmem gets all references from template msgmem.
        msgmem and objmem may point to the same storage (this feature is
        used by PutDepData() ).
 */

static void LocalizeObject (int merge_mode, TYPE_DESC *desc,
                            char    *msgmem,
                            DDD_OBJ objmem,
                            SYMTAB_ENTRY *theSymTab)
{
  ELEM_DESC     *theElem;
  int e;
  DDD_OBJ obj = objmem;

  /*
     printf("%4d:    Localize {\n", me); fflush(stdout);
   */
  /* prepare map of structure elements */
  theElem = desc->element;

  /* loop over all pointers inside of object obj */
  for(e=0; e<desc->nElements; e++, theElem++)
  {
    if (theElem->type==EL_OBJPTR)
    {
      TYPE_DESC *refdesc = &theTypeDefs[theElem->reftype];
      int l;
#ifdef C_FRONTEND
      char      *msgrefarray = msgmem+theElem->offset;
      char      *objrefarray = objmem+theElem->offset;
      /*
         printf("%4d:    Localize e=%d typ=%s reftyp=%d size=%d\n",
              me,e,desc->name,theElem->reftype,theElem->size); fflush(stdout);
       */
#else
      char      *objrefarray = theElem->array + (theElem->size*obj);
#endif

      /* loop over single pointer array */
      for(l=0; l<theElem->size; l+=SIZEOF_REF)
      {
        INT stIdx;

        /* ref points to a reference inside objmem */
        DDD_OBJ *ref = (DDD_OBJ *) (objrefarray+l);


        /* reference had been replaced by SymTab-index */
#ifdef C_FRONTEND
        stIdx = (*(INT *)(msgrefarray+l)) - 1;
#else
        /* TODO: this is from V1_6_4_F77_3, not the actual version. */
        stIdx = ((INT)*ref) - 1;
#endif


        /* test for Localize execution in merge_mode */
        if (merge_mode && (*ref!=NULL_REF))
        {
          /* if we are in merge_mode, we do not update
             existing references. */
                                        #ifdef DEBUG_MERGE_MODE
          printf("%4d: loc-merge curr=%08x keep     e=%d l=%d\n",
                 me, OBJ_GID(OBJ2HDR(*ref,refdesc)), e,l);
                                        #endif

          /* it may happen here that different references
             are in incoming and existing object. this is implicitly
             resolved by using the existing reference and ignoring
             the incoming one. if the REF_COLLISION option is set,
                  we will issue a warning.
           */

          if (stIdx>=0 &&
              DDD_GetOption(OPT_WARNING_REF_COLLISION)==OPT_ON)
          {
            /* get corresponding symtab entry */
            if (theSymTab[stIdx].adr.hdr!=OBJ2HDR(*ref,refdesc))
            {
              sprintf(cBuffer,
                      "reference collision in %08x "
                      "(old=%08x, inc=%08x) in LocalizeObject\n",
                      OBJ_GID(OBJ2HDR(obj,desc)),
                      OBJ_GID(OBJ2HDR(*ref,refdesc)),
                      OBJ_GID(theSymTab[stIdx].adr.hdr));
              DDD_PrintError('W', 6540, cBuffer);
              /* assert(0);  ??? */
            }
          }

          continue;
        }


        /*
           printf("%4d:    Localize adr=%08x l=%d ref=%08x *ref=%08x stIdx=%d\n",
                me,objmem,l,ref,*ref,(int)stIdx); fflush(stdout);
         */
        if (stIdx>=0)
        {
          /* get corresponding symtab entry */
          SYMTAB_ENTRY *st = &(theSymTab[stIdx]);

          /*
                  convert reference from header to object itself
                  and replace index by pointer; if header==NULL,
                  referenced object is not known and *ref should
                  therefore be NULL, too!
           */

#ifdef MERGE_MODE_IN_TESTZUSTAND
          if (merge_mode)
          {
            if (st->adr.hdr!=NULL)
            {
                                                        #ifdef DEBUG_MERGE_MODE
              printf("%4d: loc-merge curr=%08x have_sym e=%d l=%d to %08x\n",
                     me, *ref, e,l,OBJ_GID(st->adr.hdr));
                                                        #endif

              *ref = HDR2OBJ(st->adr.hdr,refdesc);
            }
                                                #ifdef DEBUG_MERGE_MODE
            else
            {
              printf(
                "%4d: loc-merge curr=%08x have_sym e=%d l=%d to NULL_REF\n",
                me, *ref, e, l);
            }
                                                #endif
          }
          else
#endif
          {
            if (st->adr.hdr!=NULL)
              *ref = HDR2OBJ(st->adr.hdr,refdesc);
            else
              *ref = NULL_REF;
          }
        }
        else
        {
#ifdef MERGE_MODE_IN_TESTZUSTAND
          if (merge_mode)
          {
                                                #ifdef DEBUG_MERGE_MODE
            printf("%4d: loc-merge curr=%08x no_sym   e=%d l=%d\n",
                   me, *ref, e,l);
                                                #endif
          }
          else
#endif
          {
            *ref = NULL_REF;
          }
        }
      }
    }
  }
}



/****************************************************************************/



static void PutDepData (char *data,
                        TYPE_DESC *desc,
                        DDD_OBJ obj,
                        SYMTAB_ENTRY *theSymTab,
                        int newness)
{
  TYPE_DESC    *descDep;
  char         *chunk, *curr, *adr, **table;
  int i, j, chunks;
  int addCnt;
  DDD_TYPE addTyp;
  int          *depTable, depTableSize;


  /* get overall number of chunks */
  chunks = ((int *)data)[0];
  chunk  = data + CEIL(sizeof(int));


  /*
     printf("%4d: PutDepData (chunks=%d) {\n", me, chunks); fflush(stdout);
   */

  /* loop through all chunks */
  for(j=0; j<chunks; j++)
  {
    /* first entries of chunk are addCnt and addTyp */
    addCnt = ((int *)chunk)[0];
    addTyp = ((DDD_TYPE *)chunk)[1];
    chunk += CEIL(sizeof(int)+sizeof(DDD_TYPE));

    if (addCnt>=0)
    {
      if (addTyp!=DDD_USER_DATA)
      {
        /* convert pointers using SymTab */
        descDep = &theTypeDefs[addTyp];
        curr = chunk;
        for(i=0; i<addCnt; i++)
        {
          /* insert pointers into copy using SymTab */
          if (descDep->nPointers>0)
          {
            LocalizeObject(FALSE, descDep,
                           curr,
                           (DDD_OBJ)curr,
                           theSymTab);
          }
          curr += CEIL(descDep->size);

          /*
             printf("%4d: PutDepData   chunk %d, item %d/%d\n", me, j, i, addCnt); fflush(stdout);
           */
        }
      }
      else
      {
        /* addType==DDD_USER_DATA ->
              scatter stream of bytes with len addCnt */
        curr = chunk + CEIL(addCnt);
      }

      /* scatter data via handler */
      /*
         printf("%4d: PutDepData   XFERSCATTER ...\n", me); fflush(stdout);
       */
      if (desc->handler[HANDLER_XFERSCATTER]!=NULL)
        desc->handler[HANDLER_XFERSCATTER](obj,
                                           addCnt, addTyp, (void *)chunk, newness);
      /*
         printf("%4d: PutDepData   XFERSCATTER ok\n", me); fflush(stdout);
       */
    }
    else
    {
      /* variable sized chunks */
      addCnt *= -1;

      /* convert offset table into pointer table */
      descDep = &theTypeDefs[addTyp];
      table = (char **)chunk;
      chunk += CEIL(sizeof(int)*addCnt);
      for(i=0, adr=chunk; i<addCnt; i++)
      {
        table[i] = ((INT)table[i])+adr;

        /* insert pointers into copy using SymTab */
        if (addTyp!=DDD_USER_DATA)
        {
          curr = table[i];
          if (descDep->nPointers>0)
            LocalizeObject(FALSE, descDep,
                           curr,
                           (DDD_OBJ)curr,
                           theSymTab);
        }
      }

      /* scatter data via handler */
      if (desc->handler[HANDLER_XFERSCATTERX]!=NULL)
        desc->handler[HANDLER_XFERSCATTERX](obj,
                                            addCnt, addTyp, table, newness);
    }


    /*
       printf("%4d: PutDepData   finished chunk %d\n", me, j); fflush(stdout);
     */
    chunk = curr;
  }

  /*
     printf("%4d: PutDepData }\n", me, chunks); fflush(stdout);
   */
}




/****************************************************************************/


static void AcceptObjFromMsg (
  OBJTAB_ENTRY *theObjTab, int lenObjTab,
  char *theObjects,
  DDD_HDR *localCplObjs, int nLocalCplObjs)
{
  int i, j;

  for(i=0, j=0; i<lenObjTab; i++)
  {
    OBJTAB_ENTRY *ote = &theObjTab[i];
    TYPE_DESC    *desc = &theTypeDefs[ote->typ];

    if (ote->is_new == OTHERMSG)
    {
      /* object is in another message with higher priority */
      continue;
    }

    while ((j<nLocalCplObjs) && (OBJ_GID(localCplObjs[j]) < ote->gid))
      j++;

    if ((j<nLocalCplObjs) && (OBJ_GID(localCplObjs[j])==ote->gid))
    {
      /* object already here, compare priorities.
         this is the implementation of rule XFER-C3. */
      DDD_PRIO newprio;
      int ret = PriorityMerge(desc,
                              ote->prio, OBJ_PRIO(localCplObjs[j]), &newprio);

      if (ret==PRIO_FIRST || ret==PRIO_UNKNOWN)                    /* incoming is higher or equal */
      {
        DDD_OBJ copy;

#                       if DebugUnpack<=1
        sprintf(cBuffer, "%4d: NewPrio wins. %07x\n",me,
                ote->gid);
        DDD_PrintDebug(cBuffer);
#                       endif

        /* new priority wins -> recreate */
        /* all GDATA-parts are overwritten by contents of message */
        copy = (DDD_OBJ)(theObjects+ote->offset);
        ObjCopyGlobalData(desc,
                          HDR2OBJ(localCplObjs[j],desc), copy, ote->size);

        ote->is_new = PARTNEW;
      }
      else                    /* existing is higher than incoming */
      {
#                               if DebugUnpack<=1
        sprintf(cBuffer, "%4d: OldPrio wins. %07x\n",me,
                ote->gid);
        DDD_PrintDebug(cBuffer);
#                               endif

        /* new priority looses -> keep existing obj */
        ote->is_new = NOTNEW;
      }

      /* store pointer to local object */
      ote->hdr = localCplObjs[j];

      /* store old priority and set new one */
      ote->prio    = newprio;
      ote->oldprio = OBJ_PRIO(localCplObjs[j]);
      OBJ_PRIO(localCplObjs[j]) = newprio;
    }
    else
    {
      DDD_OBJ msgcopy, newcopy;

#                       if DebugUnpack<=1
      sprintf(cBuffer, "%4d: NewObject        %07x\n",me,
              ote->gid);
      DDD_PrintDebug(cBuffer);
#                       endif

      /* new object, create local copy */
      msgcopy = (DDD_OBJ)(theObjects+ote->offset);
      newcopy = DDD_ObjNew(ote->size, ote->typ, ote->prio, ote->attr);
      ote->hdr = OBJ2HDR(newcopy,desc);

      /* copy GDATA */
      ObjCopyGlobalData(desc, newcopy, msgcopy, ote->size);
      ote->is_new = TOTALNEW;

      /* construct HDR */
      DDD_HdrConstructorCopy(ote->hdr, ote->prio);

      /* construct LDATA */
      if (desc->handler[HANDLER_LDATACONSTRUCTOR]!=NULL)
#ifdef C_FRONTEND
        desc->handler[HANDLER_LDATACONSTRUCTOR](newcopy);
#else
        desc->handler[HANDLER_LDATACONSTRUCTOR](&newcopy);
#endif
    }
  }
}



static void AcceptReceivedObjects (
  LC_MSGHANDLE *theMsgs, int nRecvMsgs,
  OBJTAB_ENTRY **allRecObjs, int nRecObjs,
  DDD_HDR *localCplObjs, int nLocalCplObjs)
{
  /*
          allRecObjs is a pointer array to all OBJTAB_ENTRYs
          received in incoming messages. it is sorted according
          to (gid/ascending).

          1. collision detection for incoming objects with same
                  gid: accept object with merged priority. if several
                  such objects exist, an arbitrary one is chosen.
                  discard all other objects with same gid. (RULE XFER-C2)

          2. transfer objects from message into local memory.

          3. propagate hdr-pointer to all OBJTAB_ENTRYs with equal gid.
   */

  int i;

  if (nRecObjs==0)
    return;

  /* 1. collision detection */
  for(i=nRecObjs-1; i>0; i--)
  {
    if (allRecObjs[i]->gid != allRecObjs[i-1]->gid)
    {
      allRecObjs[i]->is_new = THISMSG;
    }
    else
    {
      DDD_PRIO newprio;
      int ret;

      ret = PriorityMerge(&theTypeDefs[allRecObjs[i]->typ],
                          allRecObjs[i]->prio, allRecObjs[i-1]->prio, &newprio);

      if (ret==PRIO_FIRST || ret==PRIO_UNKNOWN)
      {
        /* item i is winner */
        OBJTAB_ENTRY *tmp;

        allRecObjs[i]->prio = newprio;

        /* switch first item i in second position i-1 */
        /* item on first position will be discarded */
        tmp = allRecObjs[i];
        allRecObjs[i] = allRecObjs[i-1];
        allRecObjs[i-1] = tmp;
      }
      else
      {
        /* item i-1 is winner */
        allRecObjs[i-1]->prio = newprio;
      }

      /* mark item i invalid */
      allRecObjs[i]->is_new = OTHERMSG;
    }
  }
  allRecObjs[0]->is_new = THISMSG;

  /* now the first item in a series of items with equal gid is the
     THISMSG-item, the following are OTHERMSG-items */


  /* 2. transfer from message into local memory */
  for(i=0; i<nRecvMsgs; i++)
  {
    LC_MSGHANDLE xm = theMsgs[i];

    AcceptObjFromMsg(
      (OBJTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.objtab_id),
      LC_GetTableLen(xm, xferGlobals.objtab_id),
      (char *) LC_GetPtr(xm, xferGlobals.objmem_id),
      localCplObjs, nLocalCplObjs
      );
  }


  /* 3. propagate hdr-pointer */
  for(i=1; i<nRecObjs; i++)
  {
    if (allRecObjs[i]->is_new == OTHERMSG)
    {
      /* propagate hdr-pointer */
      allRecObjs[i]->hdr = allRecObjs[i-1]->hdr;
    }
  }
}



/****************************************************************************/

/*
        this function updates the couplings of local objects.
        the inputs for deciding which couplings have to be added are:

        for prev. existing objects:
                -  sending to new_owner-destinations
                -  incoming NewCpl-items for previously existing objects

        for new (incoming) objects:
                -  incoming NewCpl-items for new objects

        as a side effect this function sends XIAddCpl-items
        to new_owner-procs.
 */


static void UpdateCouplings (
  TENewCpl *itemsNC, int nNC,                     /* NewCpl  */
  OBJTAB_ENTRY **itemsO, int nO,                  /* Objects */
  DDD_HDR *itemsLCO, int nLCO,                   /* local objs with coupling */
  XIDelObj  **itemsDO, int nDO,                   /* XIDelObj */
  XICopyObj **itemsNO, int nNO)                   /* NewOwners */
{
  int iNC, iO, iDO, iNO, iLCO;

  /*
          each NewCpl either corresponds to an incoming object
          or to a local object, but not both.
   */

  /* loop for all incoming objects */
  for(iO=0, iNC=0, iDO=0; iO<nO; iO++)
  {
    REGISTER DDD_HDR hdr = itemsO[iO]->hdr;
    REGISTER DDD_GID gid = itemsO[iO]->gid;

    /* scan NewCpl-entries for given gid */
    while (iNC<nNC && itemsNC[iNC].gid < gid)
      iNC++;

    /* scan DelObj-entries for given gid */
    while (iDO<nDO && itemsDO[iDO]->gid < gid)
      iDO++;


    /* if there is a DelObj-item with same gid, then the object
       has been deleted and send by a remote proc afterwards.
       we must:
         - restore old couplings locally
             - invalidate XIDelCpl-items
     */
    if (iDO<nDO && itemsDO[iDO]->gid == gid)
    {
      XIDelCpl *dc = itemsDO[iDO]->delcpls;
      for( ; dc!=NULL; dc=dc->next)
      {
        /* restore previous coupling */
        if (dc->prio>=0)
          AddCoupling(hdr, dc->to, dc->prio);

        /* invalidate XIDelCpl-item */
        dc->to=procs;
      }

      /* restore only one time */
      itemsDO[iDO]->delcpls = NULL;
    }


    /* for all NewCpl-Items with same gid as incoming object */
    while (iNC<nNC && itemsNC[iNC].gid == gid)
    {
      /* there is a corresponding NewCpl-item */
      AddCoupling(hdr, itemsNC[iNC].dest, itemsNC[iNC].prio);

      {
        XIAddCpl *xc = NewXIAddCpl();
        xc->to      = itemsNC[iNC].dest;
        xc->te.gid  = gid;
        xc->te.proc = me;
        xc->te.prio = OBJ_PRIO(hdr);
      }

      iNC++;
    }
  }


  /* loop for previously existing objects */
  iNO=0; iNC=0; iLCO=0; iDO=0; iO=0;
  while (iNO<nNO || iNC<nNC)
  {
    DDD_GID gid;
    /*
       printf("%4d: NO-NC n %d-%d   i %d-%d\n", me, nNO, nNC, iNO, iNC);
     */

    /* scan all NewCpl-items with same (gid/dest), and take
       the last of each group (the one with highest prio). */
    while (iNC<nNC-1 && itemsNC[iNC+1].gid==itemsNC[iNC].gid &&
           itemsNC[iNC+1].dest==itemsNC[iNC].dest)
      iNC++;


    if (iNO>=nNO)
    {
      if (iNC<nNC)
      {
        TENewCpl *nc = &(itemsNC[iNC]);
        DDD_GID gidNC = nc->gid;

        /* scan local objects with couplings */
        while (iLCO<nLCO && OBJ_GID(itemsLCO[iLCO])<gidNC)
          iLCO++;

        if (iLCO<nLCO && OBJ_GID(itemsLCO[iLCO])==gidNC)
        {
          AddCoupling(itemsLCO[iLCO], nc->dest, nc->prio);
        }

        iNC++;
      }
      /* else: no moreNOs and no moreNCs, do nothing */
    }
    else
    {
      /* there is a new_owner-XICopyObj-item */
      XICopyObj *no = itemsNO[iNO];
      DDD_GID gidNO = no->gid;

      if (iNC>=nNC)
      {
        /* no more NewCpl-items */

        /* scan received objects */
        while (iO<nO && itemsO[iO]->gid < gidNO)
          iO++;

        /* check whether obj has been deleted during this xfer */
        while (iDO<nDO && itemsDO[iDO]->gid < gidNO)
          iDO++;

        if (! (iDO<nDO && itemsDO[iDO]->gid==gidNO))
        {
          /* there is no DelObj-item */
          AddCoupling(no->hdr, no->dest, no->prio);
        }
        else if (iO<nO && itemsO[iO]->gid==gidNO)
        {
          /* obj has been deleted and received again */
          AddCoupling(itemsO[iO]->hdr, no->dest, no->prio);
        }
        /*
                else: object has been deleted, and not been resent
                          ->hdr is invalid
                      there is no need for AddCoupling here.
         */

        iNO++;
      }
      else
      {
        /* moreNOs and moreNCs */
        TENewCpl *nc = &(itemsNC[iNC]);
        DDD_GID gidNC = nc->gid;
        DDD_HDR hdr   = NULL;

        /* scan local objects with couplings */
        while (iLCO<nLCO && OBJ_GID(itemsLCO[iLCO])<gidNC)
          iLCO++;

        /* check whether obj has been deleted during this xfer */
        while (iDO<nDO && itemsDO[iDO]->gid < gidNO)
          iDO++;

        /* scan received objects */
        while (iO<nO && itemsO[iO]->gid < gidNO)
          iO++;

        /* now the minimum of gids is relevant */
        if (gidNC==gidNO)
        {
          DDD_PROC destNC = nc->dest;
          DDD_PROC destNO = no->dest;


          if (! (iDO<nDO && itemsDO[iDO]->gid==gidNO))
          {
            hdr = no->hdr;
          }
          if (iLCO<nLCO && OBJ_GID(itemsLCO[iLCO])==gidNC)
          {
            hdr = itemsLCO[iLCO];
          }
          if (iO<nO && itemsO[iO]->gid==gidNO)
          {
            hdr = itemsO[iO]->hdr;
          }

          if (destNO < destNC)
          {
            if (hdr!=NULL)
              AddCoupling(hdr, destNO, itemsNO[iNO]->prio);
            iNO++;
          }
          else
          {
            if (destNO > destNC)
            {
              if (hdr!=NULL)
                AddCoupling(hdr, nc->dest, nc->prio);
#ifdef WANTED_NOCH_GENAUER_UNTERSUCHEN
              else { printf("%4d: WANTED 3  %d/%d\n",me,nc->dest,nc->prio); }
#endif
              iNC++;
            }
            else                                     /* destNO == destNC */
            {
              if (hdr!=NULL)
              {
                TYPE_DESC *desc = &theTypeDefs[OBJ_TYPE(hdr)];
                DDD_PRIO newprio;

                PriorityMerge(desc, no->prio, nc->prio, &newprio);
                AddCoupling(hdr, nc->dest, newprio);
              }
#ifdef WANTED_NOCH_GENAUER_UNTERSUCHEN
              else { printf("%4d: WANTED 4  %d/%d/%d\n",me,nc->dest,no->prio,nc->prio); }
#endif
              iNO++;
              iNC++;
            }
          }
        }
        else
        {
          if (gidNC<gidNO)
          {
            if (iLCO<nLCO && OBJ_GID(itemsLCO[iLCO])==gidNC)
            {
              AddCoupling(itemsLCO[iLCO], nc->dest, nc->prio);
            }
#ifdef WANTED_NOCH_GENAUER_UNTERSUCHEN
            else { printf("%4d: WANTED 5  %d/%d\n",me,nc->dest,nc->prio); }
#endif
            iNC++;
          }
          else
          {
            /* gidNC>gidNO */
            /* no more NewCpl-items */

            if (! (iDO<nDO && itemsDO[iDO]->gid==gidNO))
            {
              /* there is no DelObj-item */
              AddCoupling(no->hdr, no->dest, no->prio);
            }
#ifdef WANTED_NOCH_GENAUER_UNTERSUCHEN
            else { printf("%4d: WANTED 6  %d/%d\n",me,nc->dest,no->prio); }
#endif
            /*
                    else: object has been deleted, ->hdr is invalid
                    there is no need for AddCoupling here.
             */

            iNO++;
          }
        }
      }
    }
  }
}




/*
        this function handles local objects, which had been here before
        xfer. during xfer, another object with same gid was received
        and lead to a higher priority. this priority change must
        be communicated to all destination-procs, to which the local
        processor sent a copy during xfer.
 */
static void PropagateIncomings (
  XICopyObj **arrayNO, int nNO,
  OBJTAB_ENTRY **allRecObjs, int nRecObjs)
{
  int iRO, iNO;

  for(iRO=0, iNO=0; iRO<nRecObjs; iRO++)
  {
    int newness = allRecObjs[iRO]->is_new;

    if (newness==PARTNEW || newness==TOTALNEW)
    {
      COUPLING *cpl;

      /* object has been local before, but changed its prio */
      OBJTAB_ENTRY *ote = allRecObjs[iRO];

      /* scan received objects */
      while ((iNO<nNO) && (arrayNO[iNO]->gid < ote->gid))
        iNO++;

      /* communicate to all new_owner-destinations */
      while (iNO<nNO && arrayNO[iNO]->gid == ote->gid)
      {
        if (newness==PARTNEW)
        {
          XIModCpl *xc = NewXIModCpl();
          xc->to      = arrayNO[iNO]->dest;                               /* receiver of XIModCpl*/
          xc->te.gid  = ote->gid;                                         /* the object's gid    */
          xc->te.prio = OBJ_PRIO(ote->hdr);                               /* the obj's new prio  */
          xc->typ     = OBJ_TYPE(ote->hdr);                               /* the obj's ddd-type  */
        }

        iNO++;
      }

      /* communicate to all procs in coupling */
      for(cpl=THECOUPLING(ote->hdr); cpl!=NULL; cpl=cpl->next)
      {
        /*
                                        if (newness==PARTNEW)
                                        {
         */
        XIModCpl *xc = NewXIModCpl();
        xc->to      = cpl->proc;                                         /* receiver of XIModCpl*/
        xc->te.gid  = ote->gid;                                           /* the object's gid   */
        xc->te.prio = OBJ_PRIO(ote->hdr);                                 /* the obj's new prio */
        xc->typ     = OBJ_TYPE(ote->hdr);                                 /* the obj's ddd-type  */
        /*
                                        }
         */
      }
    }
  }
}


/****************************************************************************/


static void LocalizeObjects (LC_MSGHANDLE xm,
                             OBJTAB_ENTRY **allRecObjs, int nRecObjs,
                             DDD_HDR *localCplObjs, int nLocalCplObjs)
{
  SYMTAB_ENTRY *theSymTab;
  OBJTAB_ENTRY *theObjTab;
  char         *theObjects;
  int i, j;
  char         *data;
  int lenSymTab = (int) LC_GetTableLen(xm, xferGlobals.symtab_id);
  int lenObjTab = (int) LC_GetTableLen(xm, xferGlobals.objtab_id);

  STAT_RESET4;


  /* get table addresses inside message buffer */
  theSymTab = (SYMTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.symtab_id);
  theObjTab = (OBJTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.objtab_id);
  theObjects = (char *)        LC_GetPtr(xm, xferGlobals.objmem_id);


  /* insert pointers to known objects into SymTab */
  for(i=0, j=0; i<lenSymTab; i++)
  {
    while ((j<nLocalCplObjs) &&
           (OBJ_GID(localCplObjs[j]) < theSymTab[i].gid))
      j++;

    if (j==nLocalCplObjs)
    {
      /* no more valid local objects */
      theSymTab[i].adr.hdr = NULL;
    }
    else
    {
      if (OBJ_GID(localCplObjs[j]) == theSymTab[i].gid)
      {
        theSymTab[i].adr.hdr = localCplObjs[j];
      }
      else
      {
        theSymTab[i].adr.hdr = NULL;
      }
    }
  }


  /* insert new pointers in SymTab */
  for(i=0, j=0; i<lenSymTab; i++)
  {
    while ((j<nRecObjs) && (allRecObjs[j]->gid<theSymTab[i].gid))
      j++;

    if ((j<nRecObjs) && (allRecObjs[j]->gid==theSymTab[i].gid))
    {
      theSymTab[i].adr.hdr = allRecObjs[j]->hdr;
    }
  }


#       if DebugUnpack<=3
  sprintf(cBuffer, "%4d: converting pointers\n",me);
  DDD_PrintDebug(cBuffer);
#       endif

  /* convert pointers */
  for(i=0; i<lenObjTab; i++)               /* for all message items */
  {
    if (theObjTab[i].is_new==TOTALNEW)
    {
      TYPE_DESC *desc = &theTypeDefs[theObjTab[i].typ];
      DDD_OBJ obj   = HDR2OBJ(theObjTab[i].hdr, desc);

      if (desc->nPointers>0)
      {
        LocalizeObject(FALSE, desc,
                       (char *)(theObjects+theObjTab[i].offset),
                       obj,
                       theSymTab);
      }
    }


    /*
            TODO: hier geht das wissen aus etwaigen anderen kopien
            mit is_new==OTHERMSG verloren. referenzen aus diesen kopien,
            die in der entsprechenden kopie mit is_new==TOTAL_NEW
            nicht vorhanden waren, sind zwar in theSymTab der aktuellen
            msg bekannt, werden jedoch einfach ignoriert.
            hier muesste also strenggenommen eine art verschmelzung der
            information aus den verschiedenen symboltabellen hin.
            960509 KB
            implemented merge_mode for Localize. references from all copies
            will be merged into the local copy. 960813 KB
     */

                #ifdef MERGE_MODE_IN_TESTZUSTAND
    if (theObjTab[i].is_new==OTHERMSG || theObjTab[i].is_new==PARTNEW)
    {
      TYPE_DESC *desc = &theTypeDefs[theObjTab[i].typ];
      DDD_OBJ obj   = HDR2OBJ(theObjTab[i].hdr, desc);

      if (desc->nPointers>0)
      {
                                #ifdef DEBUG_MERGE_MODE
        printf("%4d: LocalizeObject in merge_mode, %08x prio %d\n",
               me, theObjTab[i].gid, theObjTab[i].prio);
                                #endif

        /* execute Localize in merge_mode */
        LocalizeObject(TRUE, desc,
                       (char *)(theObjects+theObjTab[i].offset),
                       obj,
                       theSymTab);
      }
    }
                #endif
  }

#       if DebugUnpack<=4
  sprintf(cBuffer, "%4d: UnpackSingleMessage, Phase 1 ready\n",me);
  DDD_PrintDebug(cBuffer);
  fflush(stdout);
#       endif

}




static void CallUpdateHandler (LC_MSGHANDLE xm)
{
  OBJTAB_ENTRY *theObjTab;
  int lenObjTab = (int) LC_GetTableLen(xm, xferGlobals.objtab_id);
  int i;

  /* get table addresses inside message buffer */
  theObjTab = (OBJTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.objtab_id);

  /* initialize new objects corresponding to application: update */
  for(i=0; i<lenObjTab; i++)               /* for all message items */
  {
    if (theObjTab[i].is_new == TOTALNEW)
    {
      TYPE_DESC *desc = &theTypeDefs[theObjTab[i].typ];
      DDD_OBJ obj   = HDR2OBJ(theObjTab[i].hdr, desc);

      /* call application handler for object updating */
      if (desc->handler[HANDLER_UPDATE]!=NULL)
        desc->handler[HANDLER_UPDATE](obj);
    }
  }
}



static void UnpackAddData (LC_MSGHANDLE xm)
{
  SYMTAB_ENTRY *theSymTab;
  OBJTAB_ENTRY *theObjTab;
  char         *theObjects;
  int i;
  char         *data;
  int lenSymTab = (int) LC_GetTableLen(xm, xferGlobals.symtab_id);
  int lenObjTab = (int) LC_GetTableLen(xm, xferGlobals.objtab_id);


  /* get table addresses inside message buffer */
  theSymTab = (SYMTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.symtab_id);
  theObjTab = (OBJTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.objtab_id);
  theObjects = (char *)        LC_GetPtr(xm, xferGlobals.objmem_id);


  /* scatter additional data via handler */
  for(i=0; i<lenObjTab; i++)               /* for all message items */
  {
    if (theObjTab[i].addLen>0)
    {
      TYPE_DESC *desc = &theTypeDefs[theObjTab[i].typ];
      DDD_OBJ obj   = HDR2OBJ(theObjTab[i].hdr, desc);
      int newness;

      /*
         printf("%4d: scatter %d/%d, addLen=%d, objadr=%08x gid=%08x\n",
              me,i,lenObjTab,theObjTab[i].addLen,obj,theObjTab[i].gid);
         fflush(stdout);
       */


      switch (theObjTab[i].is_new)
      {
      case OTHERMSG : newness=XFER_REJECT;   break;
      case NOTNEW :   newness=XFER_REJECT;   break;
      case PARTNEW :  newness=XFER_UPGRADE;  break;
      case TOTALNEW : newness=XFER_NEW;      break;
      }

      /*
              compute begin of data section. theObjTab[i].size is equal to
              desc->len for fixed sized objects and different for variable
              sized objects
       */
      data = (char *)(theObjects +
                      theObjTab[i].offset +
                      CEIL(theObjTab[i].size));

      PutDepData(data, desc, obj, theSymTab, newness);
    }
  }
}



/*
        in order to allow application reactions on a priority
        change, the SETPRIORITY-handler is called.

        TODO: is this really a reason for calling SETPRIORITY? or
        should there be a separate handler for this task?

        NOTE: due to the current implementation, the new priority
        has already been set in the local object's DDD_HEADER.
        but the SETPRIORITY-handler has to get the old priority inside
        the object and the new one as second argument. so we restore
        the old prio before calling the handler and set the newprio
        afterwards.
 */
static void CallSetPriorityHandler (LC_MSGHANDLE xm)
{
  OBJTAB_ENTRY *theObjTab;
  int lenObjTab = (int) LC_GetTableLen(xm, xferGlobals.objtab_id);
  int i;

  /* get table addresses inside message buffer */
  theObjTab = (OBJTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.objtab_id);

  for(i=0; i<lenObjTab; i++)               /* for all message items */
  {
    if ((theObjTab[i].is_new==NOTNEW || theObjTab[i].is_new==PARTNEW)
        && (theObjTab[i].oldprio != theObjTab[i].prio))
    {
      TYPE_DESC *desc = &theTypeDefs[theObjTab[i].typ];
      DDD_OBJ obj   = HDR2OBJ(theObjTab[i].hdr, desc);

      /* call application handler for object consistency */
      if (desc->handler[HANDLER_SETPRIORITY]!=NULL)
      {
        /* restore old priority in object */
        OBJ_PRIO(theObjTab[i].hdr) = theObjTab[i].oldprio;
        desc->handler[HANDLER_SETPRIORITY](obj, theObjTab[i].prio);

        /* restore new priority */
        OBJ_PRIO(theObjTab[i].hdr) = theObjTab[i].prio;
      }
    }
  }
}



static void CallObjMkConsHandler (LC_MSGHANDLE xm)
{
  OBJTAB_ENTRY *theObjTab;
  int lenObjTab = (int) LC_GetTableLen(xm, xferGlobals.objtab_id);
  int i;

  STAT_RESET4;

  /* get table addresses inside message buffer */
  theObjTab = (OBJTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.objtab_id);


  /* initialize new objects corresponding to application: consistency */
  for(i=0; i<lenObjTab; i++)               /* for all message items */
  {
    if (theObjTab[i].is_new==PARTNEW || theObjTab[i].is_new==TOTALNEW)
    {
      TYPE_DESC *desc = &theTypeDefs[theObjTab[i].typ];
      DDD_OBJ obj   = HDR2OBJ(theObjTab[i].hdr, desc);

      int newness = (theObjTab[i].is_new==PARTNEW) ?
                    XFER_UPGRADE : XFER_NEW;

      /* call application handler for object consistency */
      if (desc->handler[HANDLER_OBJMKCONS]!=NULL)
        desc->handler[HANDLER_OBJMKCONS](obj, newness);
    }
  }


  STAT_INCTIMER4(23);
}




/****************************************************************************/


/*
        unpack table of TEOldCpl-items.

        this function is called for each incoming message. for each
        incoming object which hasn't been here before, a set of old
        couplings is added as an estimate until the second xfer-
        communication gives more details.

        for OTHERMSG-objects there is always another object copy
        with TOTALNEW-flag. only the one with TOTALNEW-flag submits
        the set of TEOldCpl, the other sets will be equal (for consistent
        datasets before the Xfer) and are therefore redundant.
 */
static void UnpackOldCplTab (
  TEOldCpl *tabOC, int nOC,
  OBJTAB_ENTRY *tabO, int nO)
{
  int iO, iOC;

  iO = iOC = 0;
  while (iOC<nOC && iO<nO)
  {
    /* skip ObjTab-items until a TOTALNEW is found */
    while (iO<nO && tabO[iO].is_new!=TOTALNEW)
      iO++;

    if (iO<nO)
    {
      /* look for TEOldCpl-items with same gid.
         note: this relies on previous sorting via
         sort_XIOldCpl on sender side. */
      while (iOC<nOC && tabOC[iOC].gid<tabO[iO].gid)
        iOC++;

      /* found some TEOldCpl-items with same gid */
      /* add couplings now */
      while (iOC<nOC && tabOC[iOC].gid==tabO[iO].gid)
      {
        AddCoupling(tabO[iO].hdr,tabOC[iOC].proc,tabOC[iOC].prio);
        iOC++;
      }

      iO++;
    }
  }
}


/****************************************************************************/



/*
        main unpack procedure.
 */

void XferUnpack (LC_MSGHANDLE *theMsgs, int nRecvMsgs,
                 DDD_HDR *localCplObjs, int nLocalCplObjs,
                 XISetPrio **arraySP, int nSP,
                 XIDelObj **arrayDO, int nDO,
                 XICopyObj **arrayCO, int nCO,
                 XICopyObj **arrayNewOwners, int nNewOwners)
{
  TENewCpl     *allNewCpl;
  OBJTAB_ENTRY **unionObjTab;
  int lenObjTab, lenSymTab, nNewCpl;
  int i, pos1, pos2, len;


  lenObjTab=lenSymTab=nNewCpl=0;

  for(i=0; i<nRecvMsgs; i++)
  {
    LC_MSGHANDLE xm = theMsgs[i];
    lenObjTab += (int)LC_GetTableLen(xm, xferGlobals.objtab_id);
    lenSymTab += (int)LC_GetTableLen(xm, xferGlobals.symtab_id);
    nNewCpl += (int)LC_GetTableLen(xm, xferGlobals.newcpl_id);
  }

#       if DebugUnpack<=4
  sprintf(cBuffer, "%4d: SUM OF OBJ=%3d SYM=%3d NEW=%3d FROM %2d MSGS\n",
          me, lenObjTab, lenSymTab, nNewCpl, nRecvMsgs);
  DDD_PrintDebug(cBuffer);
#       endif

  STAT_RESET3;

  if (nNewCpl>0)
  {
    allNewCpl = (TENewCpl *) AllocTmp(sizeof(TENewCpl)*nNewCpl);

    if (allNewCpl==NULL) {
      DDD_PrintError('E', 6560, "not enough memory in XferUnpack");
      return;
    }
  } else {
    allNewCpl = NULL;
  }


  if (lenObjTab>0)
  {
    unionObjTab = (OBJTAB_ENTRY **)
                  AllocTmp(sizeof(OBJTAB_ENTRY *)*lenObjTab);

    if (unionObjTab==NULL) {
      DDD_PrintError('E', 6562, "not enough memory in XferUnpack");
      return;
    }
  } else {
    unionObjTab = NULL;
  }


  STAT_TIMER3(25); STAT_RESET3;


  /* create union tables: allNewCpl, unionObjTab */
  for(i=0, pos1=pos2=0; i<nRecvMsgs; i++)
  {
    LC_MSGHANDLE xm = theMsgs[i];

    len = LC_GetTableLen(xm, xferGlobals.newcpl_id);
    if (len>0)
    {
      memcpy(allNewCpl+pos1, LC_GetPtr(xm,xferGlobals.newcpl_id),
             sizeof(TENewCpl)*len);
      pos1 += len;
    }

    len = LC_GetTableLen(xm, xferGlobals.objtab_id);
    if (len>0)
    {
      OBJTAB_ENTRY *msg_ot = (OBJTAB_ENTRY *)
                             LC_GetPtr(xm,xferGlobals.objtab_id);
      OBJTAB_ENTRY **all_ot = unionObjTab+pos2;
      int oti;
      for(oti=0; oti<len; oti++, all_ot++, msg_ot++)
        *all_ot = msg_ot;

      pos2 += len;
    }
  }

  if (nNewCpl>0)
    qsort(allNewCpl, nNewCpl, sizeof(TENewCpl), sort_TENewCpl);

  if (lenObjTab>0)
    qsort(unionObjTab, lenObjTab,
          sizeof(OBJTAB_ENTRY *), sort_ObjTabPtrs);


#       if DebugUnpack<=2
  for(i=0; i<nNewCpl; i++)
  {
    sprintf(cBuffer, "%4d:   allNewCpl %08x on %4d/%d\n",me,
            allNewCpl[i].gid,allNewCpl[i].dest,allNewCpl[i].prio);
    DDD_PrintDebug(cBuffer);
  }
#       endif



  /* accept all received objects */
  if (nRecvMsgs>0)
  {
    AcceptReceivedObjects(
      theMsgs, nRecvMsgs,
      unionObjTab, lenObjTab,
      localCplObjs, nLocalCplObjs
      );
  }


  /*
          TODO: the following loops can be implemented more
          efficiently. in each loop, there is another loop
          across all objects inside the message. for each object,
          the TypeDesc is computed. the typedesc pointers should
          be computed once and stored somewhere. kb 970115
   */

  /* unpack all messages and update local topology */
  for(i=0; i<nRecvMsgs; i++)
    LocalizeObjects(theMsgs[i], unionObjTab, lenObjTab,
                    localCplObjs, nLocalCplObjs);

  /*
          at this point all new objects are established,
          their references point correctly to the neighbour objects.
          note: the references from neighbours to the new objects
          are not actualized yet! this has to be done via the
          application handler OBJMKCONS.
   */

  /* KB 941109
          the order of the next steps is crucial:
          1. update objects via handler
          2. add additional data items
          3. call set-prio handlers
          4. update consistency
   */

  /* for PARTNEW objects */
  for(i=0; i<nRecvMsgs; i++)
    CallSetPriorityHandler(theMsgs[i]);

  /* for TOTALNEW objects */
  for(i=0; i<nRecvMsgs; i++)
    CallUpdateHandler(theMsgs[i]);

  /* for all incoming objects */
  for(i=0; i<nRecvMsgs; i++)
    UnpackAddData(theMsgs[i]);

  /* for PARTNEW and TOTALNEW objects */
  for(i=0; i<nRecvMsgs; i++)
    CallObjMkConsHandler(theMsgs[i]);



#       if DebugXfer>1
  if (DDD_GetOption(OPT_DEBUG_XFERMESGS)==OPT_ON)
#       endif
  {
    for(i=0; i<nRecvMsgs; i++)
      XferDisplayMsg("OR", theMsgs[i]);
  }



  /* unpack all OldCpl-tabs */
  for(i=0; i<nRecvMsgs; i++)
  {
    LC_MSGHANDLE xm = theMsgs[i];
    UnpackOldCplTab(
      (TEOldCpl *)LC_GetPtr(xm,xferGlobals.oldcpl_id),
      LC_GetTableLen(xm, xferGlobals.oldcpl_id),
      (OBJTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.objtab_id),
      LC_GetTableLen(xm, xferGlobals.objtab_id) );
  }




  /* update couplings according to global cpl tab */
  UpdateCouplings(allNewCpl, nNewCpl,
                  unionObjTab, lenObjTab,
                  localCplObjs, nLocalCplObjs,
                  arrayDO, nDO,
                  arrayNewOwners, nNewOwners
                  );

  /* create new XI???Cpl-infos depending on allNewCpls for existing
     objects */
  PropagateCplInfos(arraySP, nSP, arrayDO, nDO, allNewCpl, nNewCpl);


  /* create some more XIModCpl-items due to incoming objects */
  PropagateIncomings(arrayNewOwners, nNewOwners, unionObjTab, lenObjTab);


  /* free temporary memory */
  if (allNewCpl!=NULL)
    FreeTmp(allNewCpl);
  if (unionObjTab!=NULL)
    FreeTmp(unionObjTab);
}



/****************************************************************************/
