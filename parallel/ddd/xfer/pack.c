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

#include "dddi.h"
#include "xfer.h"


#define DebugPack   6  /* off: 6 */



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


static int n34=0, n35=0, n36=0;


static int sort_SymTabEntries (const void *e1, const void *e2)
{
  SYMTAB_ENTRY   *ci1, *ci2;
  /*n34++;*/

  ci1 = (SYMTAB_ENTRY *)e1;
  ci2 = (SYMTAB_ENTRY *)e2;

  if (ci1->gid < ci2->gid) return(-1);
  if (ci1->gid == ci2->gid) return(0);
  return(1);
}


static int sort_ObjTabEntries (const void *e1, const void *e2)
{
  OBJTAB_ENTRY   *ci1, *ci2;
  /*n35++;*/

  ci1 = (OBJTAB_ENTRY *)e1;
  ci2 = (OBJTAB_ENTRY *)e2;

  /* sort with ascending gid */
  if (ci1->gid < ci2->gid) return(-1);
  if (ci1->gid > ci2->gid) return(1);

  return(0);
}



/****************************************************************************/
/*                                                                          */
/* Function:  BuildSymTab                                                   */
/*                                                                          */
/* Purpose:   compute message SymTab entries for one single ddd-object.     */
/*                                                                          */
/* Input:     desc: descriptor of object                                    */
/*            copy: copy of ddd-object (inside message buffer)              */
/*            theSymTab: actual portion of message SymTab                   */
/*                                                                          */
/* Output:    number of new entries into SymTab                             */
/*                                                                          */
/****************************************************************************/

static int BuildSymTab (TYPE_DESC *desc, char *copy, SYMTAB_ENTRY *theSymTab)
{
  ELEM_DESC   *theElem;
  int e, actSym;

  STAT_RESET4;

  /* reset local portion of SymTab */
  actSym = 0;

  /* prepare map of structure elements */
  theElem = desc->element;

  /* loop over all pointers inside of object obj */
  for(e=0; e<desc->nElements; e++, theElem++)
  {
    if (theElem->type==EL_OBJPTR)
    {
      TYPE_DESC *refdesc = &theTypeDefs[theElem->reftype];
      int l;

      /* loop over single pointer array */
#ifdef C_FRONTEND
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
          /* get header of referenced object */
          DDD_HDR refhdr = OBJ2HDR(*ref,refdesc);

          /* remember the GID of the referenced object */
          theSymTab[actSym].gid = OBJ_GID(refhdr);

          /* remember the address of the reference (in obj-copy) */
          theSymTab[actSym].adr.ref = ref;
          actSym++;
        }
      }
    }
  }

  STAT_INCTIMER4(33);

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
      if (desc->handler[HANDLER_XFERGATHER]!=NULL)
        desc->handler[HANDLER_XFERGATHER](obj,
                                          xa->addCnt, xa->addTyp, (void *)chunk);

      if (xa->addTyp<DDD_USER_DATA || xa->addTyp>DDD_USER_DATA_MAX)
      {
        /* insert pointers into symtab */
        descDep = &theTypeDefs[xa->addTyp];
        for(i=0; i<xa->addCnt; i++)
        {
          actSym += BuildSymTab(descDep,
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
      if (desc->handler[HANDLER_XFERGATHERX]!=NULL)
        desc->handler[HANDLER_XFERGATHERX](obj,
                                           xa->addCnt, xa->addTyp, table1);

      /* convert pointer table into offset table */
      table2 = (int *)table1;
      descDep = &theTypeDefs[xa->addTyp];
      adr = chunk;
      for(i=0; i<xa->addCnt; i++)
      {
        /* insert pointers into symtab */
        if (xa->addTyp<DDD_USER_DATA || xa->addTyp>DDD_USER_DATA_MAX)
        {
          actSym += BuildSymTab(descDep,
                                table1[i], &(theSymTab[actSym]));
        }

        table2[i] = table1[i]-adr;
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
    COUPLING  *cpl;
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


    /* update object table */
    theObjTab[actObj].offset = currObj-theObjects;
    theObjTab[actObj].typ    = OBJ_TYPE(hdr);
    theObjTab[actObj].gid    = OBJ_GID(hdr);
    theObjTab[actObj].prio   = xi->prio;
    theObjTab[actObj].attr   = OBJ_ATTR(hdr);
    theObjTab[actObj].hdr    = NULL;
    theObjTab[actObj].addLen = xi->addLen;
    theObjTab[actObj].size   = xi->size;              /* needed for variable-sized objects */
    actObj++;



    /*
            copy object into message. in the following xi->size
            equals desc->len for fixed-size objects.
     */
    STAT_RESET3;
#ifdef C_FRONTEND
    memcpy(currObj, obj, xi->size);
    copyhdr = OBJ2HDR(currObj,desc);
#else
    ObjToMsg (obj, desc, currObj);
    copyhdr = (DDD_HDR) currObj;
#endif
    STAT_INCTIMER3(32);

    /* insert priority into copy */
    OBJ_PRIO(copyhdr) = xi->prio;


    /* call application handler for direct manipulation */
    /* KB 941110:  moved from objmgr.c                  */
    /*
            Caution: this is a very, very dirty situation.
            HANDLER_COPYMANIP is able to manipulate the
            obj-copy inside the message. this handler should
            be removed in future DDD versions.
     */
#ifdef C_FRONTEND
    if (desc->handler[HANDLER_COPYMANIP]!=NULL)
    {
      /*
              NOTE: OBJ_TYPE could change during the
              execution of HANDLER_COPYMANIP. however,
              the position of DDD_HEADER inside the object
              should not change. therefore, we can remember
              the offsetHeader here and use it afterwards
              to adjust the desc.
       */
      int offset = desc->offsetHeader;

      /* now call handler */
      desc->handler[HANDLER_COPYMANIP](currObj);

      /* adjust new description according to new type */
      desc = &(theTypeDefs[OBJ_TYPE((DDD_HDR)(currObj+offset))]);
    }
#endif

    /* build symbol table portion from object copy */
    actSym += BuildSymTab(desc, (char *)currObj, &(theSymTab[actSym]));


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
  STAT_RESET3;
  qsort(theSymTab, actSym, sizeof(SYMTAB_ENTRY), sort_SymTabEntries);
  STAT_INCTIMER3(34); STAT_RESET3;

  /* sorting of objtab is necessary!! (see AcceptObjFromMsg) KB 960812 */
  qsort(theObjTab, msg->nObjects, sizeof(OBJTAB_ENTRY), sort_ObjTabEntries);
  STAT_INCTIMER3(35); STAT_RESET3;


  /*
     sprintf(cBuffer, "ITEMS 34=%d, 35=%d, 36=%d\n", actSym,  msg->nObjects,  actNewCpl);
     DDD_PrintDebug(cBuffer);
     sprintf(cBuffer, "COMPS 34=%d, 35=%d, 36=%d\n", n34,  n35,  n36);
     DDD_PrintDebug(cBuffer);
   */


  /* substitute all pointers by index into SymTab */
  STAT_RESET3;
  for(mi=0; mi<actSym; mi++)
  {
    /* patch SymTab index into reference location inside message */
#ifdef C_FRONTEND
    *(theSymTab[mi].adr.ref) = (void *)(mi+1);
#else
    *(theSymTab[mi].adr.ref) = (mi+1);
#endif
  }
  STAT_INCTIMER3(37);


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

  /* pack messages and send away */
  for(xm=theMsgs; xm!=NULL; xm=xm->next)
  {
    XferPackSingleMsg(xm);
    LC_MsgSend(xm->msg_h);
  }
}
