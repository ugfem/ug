// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      pack.c                                                        */
/*                                                                          */
/* Purpose:   packs objects into messages                                   */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/01/31 kb  begin                                            */
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
#include <string.h>

#include "dddi.h"
#include "xfer.h"




/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/



/****************************************************************************/
/*                                                                          */
/* macros                                                                   */
/*                                                                          */
/****************************************************************************/

/* helpful macros for FRONTEND switching, will be #undef'd at EOF */
#ifdef F_FRONTEND
#define _FADR     &
#else
#define _FADR
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


#ifdef SORT_STATISTIK
static int n34=0, n35=0, n36=0;
#endif


static int sort_SymTabEntries (const void *e1, const void *e2)
{
  SYMTAB_ENTRY   *ci1, *ci2;

#ifdef SORT_STATISTIK
  n34++;
#endif

  ci1 = (SYMTAB_ENTRY *)e1;
  ci2 = (SYMTAB_ENTRY *)e2;

  if (ci1->gid < ci2->gid) return(-1);
  if (ci1->gid == ci2->gid) return(0);
  return(1);
}



static char *currentObjectMem;

static int sort_ObjTabEntries (const void *e1, const void *e2)
{
  DDD_GID g1, g2;

#ifdef SORT_STATISTIK
  n35++;
#endif

  g1 = OTE_GID(currentObjectMem, (OBJTAB_ENTRY *)e1);
  g2 = OTE_GID(currentObjectMem, (OBJTAB_ENTRY *)e2);

  /* sort with ascending gid */
  if (g1 < g2) return(-1);
  if (g1 > g2) return(1);

  return(0);
}


static int sort_MsgSize (const void *e1, const void *e2)
{
  XFERMSG   *xm1, *xm2;
  size_t s1, s2;

  xm1 = *((XFERMSG **)e1);
  xm2 = *((XFERMSG **)e2);

  s1 = LC_GetBufferSize(xm1->msg_h);
  s2 = LC_GetBufferSize(xm2->msg_h);

  /* sort with descending msg-size */
  if (s1 < s2) return(1);
  if (s1 > s2) return(-1);

  return(0);
}


/****************************************************************************/
/*                                                                          */
/* Function:  BuildSymTab                                                   */
/*                                                                          */
/* Purpose:   compute message SymTab entries for one single ddd-object.     */
/*                                                                          */
/* Input:     desc: descriptor of object                                    */
/*            obj:  DDD_OBJ ptr to ddd-object (in local mem) or NULL        */
/*            copy: copy of ddd-object (inside message buffer)              */
/*            theSymTab: actual portion of message SymTab                   */
/*                                                                          */
/* Output:    number of new entries into SymTab                             */
/*                                                                          */
/****************************************************************************/

static int BuildSymTab (TYPE_DESC *desc,
                        DDD_OBJ obj,
                        char *copy,
                        SYMTAB_ENTRY *theSymTab)
{
  ELEM_DESC   *theElem;
  int e, actSym;

  /*STAT_RESET4;*/

  /* reset local portion of SymTab */
  actSym = 0;

  /* prepare map of structure elements */
  theElem = desc->element;

  /* loop over all pointers inside of object obj */
  for(e=0; e<desc->nElements; e++, theElem++)
  {
    if (theElem->type==EL_OBJPTR)
    {
      TYPE_DESC *refdesc;
      int l;
      int rt_on_the_fly = (EDESC_REFTYPE(theElem)==DDD_TYPE_BY_HANDLER);

      /* determine reftype of this elem */
      if (! rt_on_the_fly)
      {
        /* we know the reftype of this element in advance */
        refdesc = &theTypeDefs[EDESC_REFTYPE(theElem)];
      }
      /* else: determine reftype on the fly by calling handler */


      /* loop over single pointer array */
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
      for(l=0; l<theElem->size; l+=sizeof(void *))
      {
        /* get address of outside reference */
        DDD_OBJ *ref = (DDD_OBJ *)(copy+theElem->offset+l);
#else
      for(l=0; l<theElem->size; l+=sizeof(DDD_OBJ))
      {
        /* F77TODO: DDD_OBJ* must be replaced by local objindex */
        /* get the index of the referenced object */
        DDD_OBJ *ref = (DDD_OBJ *)(copy+theElem->msgoffset);
#endif

        /* create symbol table entry */
        if (*ref!=NULL)
        {
          DDD_HDR refhdr;

          if (rt_on_the_fly)
          {
            DDD_TYPE rt;

            /* determine reftype on the fly by calling handler */
            assert(obj!=NULL);                                       /* we need a real object here */

            rt = theElem->reftypeHandler(obj, *ref);
            if (rt>=MAX_TYPEDESC)
            {
              DDD_PrintError('E', 6520,
                             "invalid referenced DDD_TYPE "
                             "returned by handler");
              HARD_EXIT;
            }
            refdesc = &theTypeDefs[rt];
          }

          /* get header of referenced object */
          refhdr = OBJ2HDR(*ref,refdesc);

          /* remember the GID of the referenced object */
          theSymTab[actSym].gid = OBJ_GID(refhdr);

          /* remember the address of the reference (in obj-copy) */
          theSymTab[actSym].adr.ref = ref;
          actSym++;
        }
      }
    }
  }

  /*STAT_INCTIMER4(33);*/

  /* return SymTab increment */
  return(actSym);
}



/****************************************************************************/
/*                                                                          */
/* Function:  GetDepData                                                    */
/*                                                                          */
/* Purpose:   fill object-dependent data into message. an appl. routine     */
/*            will be called to fill in the data actually. pointers are     */
/*            localized and the message SymTab is actualized.               */
/*                                                                          */
/* Input:     data: portion of message buffer reserved for dependent data   */
/*            desc: descriptor of object                                    */
/*            obj:  current ddd-object                                      */
/*            theSymTab: actual portion of message SymTab                   */
/*            xi:   single xferinfo for current ddd-object                  */
/*                                                                          */
/* Output:    number of new entries into SymTab                             */
/*                                                                          */
/****************************************************************************/

static int GetDepData (char *data,
                       TYPE_DESC *desc,
                       DDD_OBJ obj,
                       SYMTAB_ENTRY *theSymTab,
                       XICopyObj *xi)
{
  XFERADDDATA  *xa;
  TYPE_DESC    *descDep;
  char         *chunk, *adr, **table1, *next_chunk;
  int chunks, i, actSym, *table2;


  if (xi->addLen==0) return(0);

  chunks = 0;
  actSym = 0;


  /* first entry will be number of dependency chunks */
  chunk = data + CEIL(sizeof(int));


  /* loop through whole dependency data descriptor */
  for(xa=xi->add; xa!=NULL; xa=xa->next)
  {
    /* first entries of chunk are addCnt and addTyp */
    ((int *)chunk)[0]      = xa->addCnt;
    ((DDD_TYPE *)chunk)[1] = xa->addTyp;

    if (xa->sizes==NULL)
    {
      chunk += CEIL(sizeof(int)+sizeof(DDD_TYPE));

      /* then all records should be gathered via handler */
      if (desc->handlerXFERGATHER)
      {
                                #if defined(C_FRONTEND) || defined(F_FRONTEND)
        desc->handlerXFERGATHER(_FADR obj,
                                _FADR xa->addCnt, _FADR xa->addTyp, (void *)chunk);
                                #endif
                                #ifdef CPP_FRONTEND
        CallHandler(desc,XFERGATHER) (HParam(obj)
                                      xa->addCnt, xa->addTyp, (void *)chunk);
                                #endif
      }

      if (xa->addTyp<DDD_USER_DATA || xa->addTyp>DDD_USER_DATA_MAX)
      {
        /* insert pointers into symtab */
        descDep = &theTypeDefs[xa->addTyp];
        for(i=0; i<xa->addCnt; i++)
        {
          actSym += BuildSymTab(descDep, NULL,
                                chunk, &(theSymTab[actSym]));
          chunk += CEIL(descDep->size);
        }
      }
      else
      {
        /* no regular type -> send byte stream with length addCnt */
        chunk += CEIL(xa->addCnt);
      }
    }
    else
    {
      /* var-sized AddData items */
      ((int *)chunk)[0] *= -1;
      chunk += CEIL(sizeof(int)+sizeof(DDD_TYPE));

      /* create pointer array inside message */
      table1 = (char **)chunk;
      chunk += CEIL(sizeof(int)*xa->addCnt);
      for(i=0, adr=chunk; i<xa->addCnt; i++)
      {
        table1[i] = adr;
        adr += CEIL(xa->sizes[i]);
      }
      next_chunk = adr;

      /* then all records should be gathered via handler */
      if (desc->handlerXFERGATHERX)
      {
                                #if defined(C_FRONTEND) || defined(F_FRONTEND)
        desc->handlerXFERGATHERX(_FADR obj,
                                 _FADR xa->addCnt, _FADR xa->addTyp, table1);
                                #endif
                                #ifdef CPP_FRONTEND
        CallHandler(desc,XFERGATHERX) (HParam(obj)
                                       xa->addCnt, xa->addTyp, table1);
                                #endif
      }

      /* convert pointer table into offset table */
      table2 = (int *)table1;
      descDep = &theTypeDefs[xa->addTyp];
      adr = chunk;
      for(i=0; i<xa->addCnt; i++)
      {
        /* insert pointers into symtab */
        if (xa->addTyp<DDD_USER_DATA || xa->addTyp>DDD_USER_DATA_MAX)
        {
          actSym += BuildSymTab(descDep, NULL,
                                table1[i], &(theSymTab[actSym]));
        }

        table2[i] = (int)(table1[i]-adr);
      }

      chunk = next_chunk;
    }

    /* count chunks */
    chunks++;
  }


  /* remember number of chunks at the beginning of the deplist */
  ((int *)data)[0] = chunks;


  return(actSym);
}


#ifdef F_FRONTEND

/****************************************************************************/
/*                                                                          */
/* Function:  ObjToMsg                                                      */
/*                                                                          */
/* Purpose:   copy one fortran DDD_OBJ to the message buffer.				*/
/*                                                                          */
/* Input:     obj:  the DDD_OBJ												*/
/*            desc: the type description of the object						*/
/*			  msg:  the msg memory in the buffer							*/
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void ObjToMsg (DDD_OBJ obj, TYPE_DESC *desc, char *msg)

{
  int i, i2;

  /* Copy the header to the message */
  memcpy (msg, &(desc->hdr[obj]), sizeof(DDD_HEADER));
  msg += sizeof(DDD_HEADER);

  /* Now copy the non local object data to the message */
  for (i = 0; i < desc->nElements; i++)
  {
    ELEM_DESC *elem = &(desc->element[i]);
    char      *src  = elem->array + elem->size * obj;

    if (elem->type != EL_LDATA)
      for (i2 = elem->size; i2; i2--) *msg++ = *src++;
  }
}

#endif


/****************************************************************************/
/*                                                                          */
/* Function:  XferPackSingleMsgs                                            */
/*                                                                          */
/* Purpose:   build up one outgoing message completely, fill data into      */
/*            message buffer. objects and couplings will be packed, and     */
/*            several message tables will be constructed. pointers are      */
/*            localized inside the message. several applications handlers   */
/*            are called to fill in dependent data.                         */
/*                                                                          */
/* Input:     msg:  single message-send-info structure                      */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

static void XferPackSingleMsg (XFERMSG *msg)
{
  SYMTAB_ENTRY *theSymTab;
  OBJTAB_ENTRY *theObjTab;
  TENewCpl     *theNewCpl;
  TEOldCpl     *theOldCpl;
  char         *theObjects, *currObj;
  int i, actSym, actNewCpl, actOldCpl, actObj, recvProc;
  INT mi;


  /* recipient of this message */
  recvProc = msg->proc;

  /* get table addresses inside message */
  theSymTab = (SYMTAB_ENTRY *)LC_GetPtr(msg->msg_h, xferGlobals.symtab_id);
  theObjTab = (OBJTAB_ENTRY *)LC_GetPtr(msg->msg_h, xferGlobals.objtab_id);
  theNewCpl = (TENewCpl *)    LC_GetPtr(msg->msg_h, xferGlobals.newcpl_id);
  theOldCpl = (TEOldCpl *)    LC_GetPtr(msg->msg_h, xferGlobals.oldcpl_id);
  theObjects= (char *)LC_GetPtr(msg->msg_h, xferGlobals.objmem_id);


  /* build several tables inside message */
  actSym = actNewCpl = actOldCpl = actObj = 0;
  currObj = theObjects;
  for(i=0; i<msg->nObjItems; i++)          /* for all XICopyObj-items */
  {
    REGISTER XICopyObj *xi = msg->xferObjArray[i];
    REGISTER DDD_HDR hdr   = xi->hdr;
    TYPE_DESC *desc = &theTypeDefs[OBJ_TYPE(hdr)];
    DDD_OBJ obj   = HDR2OBJ(hdr,desc);
    /*COUPLING  *cpl;*/
    DDD_HDR copyhdr;

    /* build coupling table */
    /* skip cpl which describes object itself (receive proc) */
    /*
                    for(cpl=THECOUPLING(hdr); cpl!=NULL; cpl=cpl->next)
                    {
                            if (cpl->proc!=recvProc)
                            {
                                    theNewCpl[actNewCpl].gid  = OBJ_GID(hdr);
                                    theNewCpl[actNewCpl].proc = cpl->proc;
                                    theNewCpl[actNewCpl].prio = cpl->prio;
                                    actNewCpl++;
                            }
                    }
     */

    /* one coupling for object itself (send proc) */
    /*
                    theNewCpl[actNewCpl].gid  = OBJ_GID(hdr);
                    theNewCpl[actNewCpl].proc = me;
                    theNewCpl[actNewCpl].prio = OBJ_PRIO(hdr);
                    actNewCpl++;
     */


#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
    copyhdr = OBJ2HDR(currObj,desc);
#else
    copyhdr = (DDD_HDR) currObj;
#endif

    /* update object table */
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
    theObjTab[actObj].h_offset = (int)(((char *)copyhdr)-theObjects);
#else
    theObjTab[actObj].o_offset = (int)(currObj-theObjects);
    theObjTab[actObj].typ    = OBJ_TYPE(hdr);
    theObjTab[actObj].gid    = OBJ_GID(hdr);
    theObjTab[actObj].attr   = OBJ_ATTR(hdr);
    theObjTab[actObj].prio   = xi->prio;
#endif
    theObjTab[actObj].hdr      = NULL;
    theObjTab[actObj].addLen   = xi->addLen;
    theObjTab[actObj].size     = xi->size;              /* needed for variable-sized objects */
    actObj++;



    /*
            copy object into message. in the following xi->size
            equals desc->len for fixed-size objects.
     */
    /*STAT_RESET3;*/
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
    /* NOTE: object memory is copied _completely_, i.e., also LDATA-
       components are copied into message and sent to destination.
       then, on the receiving processor the data is sorted out... */
    memcpy(currObj, obj, xi->size);
#else
    ObjToMsg (obj, desc, currObj);
#endif
    /*STAT_INCTIMER3(32);*/

    /* insert priority into copy */
    OBJ_PRIO(copyhdr) = xi->prio;


    /* call application handler for direct manipulation */
    /* KB 941110:  moved from objmgr.c                  */
    /*
            Caution: this is a very, very dirty situation.
            HANDLER_XFERCOPYMANIP is able to manipulate the
            obj-copy inside the message. this handler should
            be removed in future DDD versions.
     */
                #if defined(C_FRONTEND)
    if (desc->handlerXFERCOPYMANIP)
    {
      /*
              NOTE: OBJ_TYPE could change during the
              execution of HANDLER_XFERCOPYMANIP. however,
              the position of DDD_HEADER inside the object
              should not change. therefore, we can remember
              the offsetHeader here and use it afterwards
              to adjust the desc.
       */
      int offset = desc->offsetHeader;

      /* now call handler */
      desc->handlerXFERCOPYMANIP(currObj);

      /* adjust new description according to new type */
      desc = &(theTypeDefs[OBJ_TYPE((DDD_HDR)(currObj+offset))]);
    }
                #endif

    /* build symbol table portion from object copy */
    actSym += BuildSymTab(desc, obj, (char *)currObj, &(theSymTab[actSym]));


    /* advance to next free object slot in message, c.f. alignment */
    currObj += CEIL(xi->size);


    /* gather additional data */
    if (xi->addLen>0)
    {
      actSym += GetDepData(currObj,
                           desc, obj, &(theSymTab[actSym]), xi);
      currObj += xi->addLen;
    }
  }


  /* for all XINewCpl items in this message */
  for(i=0; i<msg->nNewCpl; i++)
  {
    theNewCpl[actNewCpl]  = msg->xferNewCpl[i]->te;
    actNewCpl++;
  }

  /* for all XIOldCpl items in this message */
  for(i=0; i<msg->nOldCpl; i++)
  {
    theOldCpl[actOldCpl]  = msg->xferOldCpl[i]->te;
    actOldCpl++;
  }



  /* sort SymTab, ObjTab and CplTab */
  /*STAT_RESET3;*/
  qsort(theSymTab, actSym, sizeof(SYMTAB_ENTRY), sort_SymTabEntries);
  /*STAT_INCTIMER3(34); STAT_RESET3;*/

  /* sorting of objtab is necessary!! (see AcceptObjFromMsg) KB 960812 */
  currentObjectMem = theObjects;
  qsort(theObjTab, msg->nObjects, sizeof(OBJTAB_ENTRY), sort_ObjTabEntries);
  /*STAT_INCTIMER3(35); STAT_RESET3;*/


#ifdef SORT_STATISTIK
  sprintf(cBuffer, "ITEMS 34=%d, 35=%d, 36=%d\n", actSym,  msg->nObjects,  actNewCpl);
  DDD_PrintDebug(cBuffer);
  sprintf(cBuffer, "COMPS 34=%d, 35=%d, 36=%d\n", n34,  n35,  n36);
  DDD_PrintDebug(cBuffer);
#endif


  /* substitute all pointers by index into SymTab */
  /*TAT_RESET3;*/
  for(mi=0; mi<actSym; mi++)
  {
    /* patch SymTab index into reference location inside message */
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
    *(theSymTab[mi].adr.ref) = (DDD_OBJ)(mi+1);
#endif
#ifdef F_FRONTEND
    *(theSymTab[mi].adr.ref) = (mi+1);
#endif
  }
  /*STAT_INCTIMER3(37);*/


  /* TODO: theSymtab[].ref wird ab hier nicht mehr verwendet und muss nicht uebertragen werden! */


  /* set valid table entries */
  LC_SetTableLen(msg->msg_h, xferGlobals.symtab_id, actSym);
  LC_SetTableLen(msg->msg_h, xferGlobals.objtab_id, msg->nObjects);
  LC_SetTableLen(msg->msg_h, xferGlobals.newcpl_id, actNewCpl);
  LC_SetTableLen(msg->msg_h, xferGlobals.oldcpl_id, actOldCpl);


#if DebugXfer>1
  if (DDD_GetOption(OPT_DEBUG_XFERMESGS)==OPT_ON)
#endif
  XferDisplayMsg("OS", msg->msg_h);
}




/****************************************************************************/
/*                                                                          */
/* Function:  XferPackMsgs                                                  */
/*                                                                          */
/* Purpose:   allocate one message buffer for each outgoing message,        */
/*            fill buffer with message contents and initiate asynchronous   */
/*            send for each message.                                        */
/*                                                                          */
/* Input:     theMsgs: list of message-send-infos                           */
/*                                                                          */
/* Output:    -                                                             */
/*                                                                          */
/****************************************************************************/

void XferPackMsgs (XFERMSG *theMsgs)
{
  XFERMSG      *xm;

#if     DebugPack<=3
  sprintf(cBuffer, "%d: XferPackMsgs\n", me);
  DDD_PrintDebug(cBuffer);
  fflush(stdout);
#endif

  /* sort messages according to decreasing size. i.e., send
     biggest message first. LowComm will use this to handle
     situations with little memory ressources. */
  {
    int i, n;
    XFERMSG **xm_array;

    /* count number of messages */
    for(n=0, xm=theMsgs; xm!=NULL; xm=xm->next) n++;

    if (n>0)
    {
      /* alloc array of pointers to messages */
      xm_array = (XFERMSG **) OO_Allocate (sizeof(XFERMSG *) * n);
      if (xm_array!=NULL)
      {
        for(i=0, xm=theMsgs; i<n; xm=xm->next, i++) xm_array[i] = xm;

        /* sort array and relink list */
        qsort(xm_array, n, sizeof(XFERMSG *), sort_MsgSize);
        theMsgs = xm_array[0];
        for(i=0; i<n-1; i++) xm_array[i]->next = xm_array[i+1];
        if (n>1) xm_array[n-1]->next = NULL;

        /* free array */
        OO_Free (xm_array /*,0*/);
      }
      /* else
         {
              resorting msg-list is not possible due to memory shortage.
              simply don't do it.
         }
       */
    }
  }



  /* allocate buffer, pack messages and send away */
  for(xm=theMsgs; xm!=NULL; xm=xm->next)
  {
    if (! LC_MsgAlloc(xm->msg_h))
    {
      sprintf(cBuffer, STR_NOMEM " in XferPackMsgs (size=%ld)",
              (unsigned long) LC_GetBufferSize(xm->msg_h));
      DDD_PrintError('E', 6522, cBuffer);
      HARD_EXIT;
    }
    XferPackSingleMsg(xm);
    LC_MsgSend(xm->msg_h);
  }
}


/****************************************************************************/

#undef _FADR


/****************************************************************************/
