// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ctrl.c                                                        */
/*                                                                          */
/* Purpose:   controls and displays messages, for debugging only            */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   94/02/16 kb  begin                                            */
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


/* #define DebugAllPointers */


/****************************************************************************/
/*                                                                          */
/* definition of static variables                                           */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


static void XferPtr (LC_MSGHANDLE xm, char *buf)
{
  SYMTAB_ENTRY *theSymTab;
  OBJTAB_ENTRY *theObjTab;
  char         *theObjects;
  int i;
  int lenSymTab = (int) LC_GetTableLen(xm, xferGlobals.symtab_id);
  int lenObjTab = (int) LC_GetTableLen(xm, xferGlobals.objtab_id);


  /* get table addresses inside message buffer */
  theSymTab = (SYMTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.symtab_id);
  theObjTab = (OBJTAB_ENTRY *) LC_GetPtr(xm, xferGlobals.objtab_id);
  theObjects = (char *)        LC_GetPtr(xm, xferGlobals.objmem_id);


  /* build symbol table */
  for(i=0; i<lenObjTab; i++)            /* for all objects in message */
  {
    DDD_HDR hdr   = (DDD_HDR)(theObjects+theObjTab[i].offset);
    TYPE_DESC *desc = &theTypeDefs[OBJ_TYPE(hdr)];
    DDD_OBJ obj   = HDR2OBJ(hdr,desc);
    int e;
    ELEM_DESC  *theElem = desc->element;

    /* loop over all pointers inside of object with DDD_HEADER hdr */
    for(e=0; e<desc->nElements; e++, theElem++)
    {
      if (theElem->type==EL_OBJPTR)
      {
        int l;

#ifdef C_FRONTEND
        for(l=0; l<theElem->size; l+=sizeof(void *))
        {
          /* ref points to a reference inside objmem */
          DDD_OBJ *ref = (DDD_OBJ *)(((char *)obj)+theElem->offset+l);
#else
        for(l=0; l<theElem->size; l+=sizeof(DDD_OBJ))
        {
          /* ref points to a reference inside objmem */
          DDD_OBJ *ref = (DDD_OBJ *)(((char *)obj)+theElem->msgoffset);
#endif
          /* reference had been replaced by SymTab-index */
          INT stIdx = ((int)*ref)-1;

          if (stIdx>=0)
          {
            /* get corresponding symtab entry */
            SYMTAB_ENTRY *st = &(theSymTab[stIdx]);

            sprintf(cBuffer, "%s 20        obj=%03d %03d st=%08x"
                    " gid=%08x (%08x==%08x)\n",
                    buf, theObjTab[i].offset, stIdx,
                    st, st->gid, st->adr.hdr, st->adr.ref);
            DDD_PrintDebug(cBuffer);
          }
        }
      }
    }
  }
}



void XferDisplayMsg (char *comment, LC_MSGHANDLE xm)
{
  SYMTAB_ENTRY *theSymTab;
  OBJTAB_ENTRY *theObjTab;
  TENewCpl     *theNewCpl;
  TEOldCpl     *theOldCpl;
  char         *theObjects, *currObj;
  char buf[30];
  int i, proc = LC_MsgGetProc(xm);
  int lenSymTab = (int) LC_GetTableLen(xm, xferGlobals.symtab_id);
  int lenObjTab = (int) LC_GetTableLen(xm, xferGlobals.objtab_id);
  int lenNewCpl = (int) LC_GetTableLen(xm, xferGlobals.newcpl_id);
  int lenOldCpl = (int) LC_GetTableLen(xm, xferGlobals.oldcpl_id);


  sprintf(buf, " %03d-%s-%03d ", me, comment, proc);

  /* get table addresses inside message */
  theSymTab = (SYMTAB_ENTRY *)LC_GetPtr(xm, xferGlobals.symtab_id);
  theObjTab = (OBJTAB_ENTRY *)LC_GetPtr(xm, xferGlobals.objtab_id);
  theNewCpl = (TENewCpl *)    LC_GetPtr(xm, xferGlobals.newcpl_id);
  theOldCpl = (TEOldCpl *)    LC_GetPtr(xm, xferGlobals.oldcpl_id);
  theObjects= (char *)LC_GetPtr(xm, xferGlobals.objmem_id);


  /* because of LC layer, this data can't be accessed anymore. KB 960718
          sprintf(cBuffer, "%s 00 MsgBuf=%08x\n", buf, xmdata);
          DDD_PrintDebug(cBuffer);
          sprintf(cBuffer, "%s 00 SymTab %04d\n", buf, theHeader->beginSymTab);
          DDD_PrintDebug(cBuffer);
          sprintf(cBuffer, "%s 01 ObjTab %04d\n", buf, theHeader->beginObjTab);
          DDD_PrintDebug(cBuffer);
          sprintf(cBuffer, "%s 02 CplTab %04d\n", buf, theHeader->beginCplTab);
          DDD_PrintDebug(cBuffer);
          sprintf(cBuffer, "%s 03 DelTab %04d\n", buf, theHeader->beginDelTab);
          DDD_PrintDebug(cBuffer);
          sprintf(cBuffer, "%s 04 ObjMem %04d\n", buf, theHeader->beginObjMem);
          DDD_PrintDebug(cBuffer);
   */


  sprintf(cBuffer, "%s 05 ObjTab.size=%05d\n", buf, lenObjTab);
  DDD_PrintDebug(cBuffer);
  sprintf(cBuffer, "%s 06 SymTab.size=%05d\n", buf, lenSymTab);
  DDD_PrintDebug(cBuffer);
  sprintf(cBuffer, "%s 07 NewCpl.size=%05d\n", buf, lenNewCpl);
  DDD_PrintDebug(cBuffer);
  sprintf(cBuffer, "%s 08 OldCpl.size=%05d\n", buf, lenOldCpl);
  DDD_PrintDebug(cBuffer);

  for(i=0; i<lenObjTab; i++)
  {
#ifdef C_FRONTEND
    DDD_OBJ obj = (DDD_OBJ)(theObjects + theObjTab[i].offset);

    sprintf(cBuffer, "%s 10 objtab    %06d typ=%1d gid=%08x "
            "hdr=%08x size=%05d add=%05d\n",
            buf, theObjTab[i].offset, theObjTab[i].typ,
            OBJ_GID(OBJ2HDR(obj,&theTypeDefs[theObjTab[i].typ])),
            theObjTab[i].hdr, theObjTab[i].size, theObjTab[i].addLen);
#else
    DDD_HDR hdr = (DDD_HDR)(theObjects + theObjTab[i].offset);

    sprintf(cBuffer, "%s 10 objtab    %06d typ=%1d gid=%08x "
            "hdr=%08x size=%05d add=%05d\n",
            buf, theObjTab[i].offset, theObjTab[i].typ,
            OBJ_GID(hdr),
            theObjTab[i].hdr, theObjTab[i].size, theObjTab[i].addLen);
#endif
    DDD_PrintDebug(cBuffer);
  }

  for(i=0; i<lenSymTab; i++)
  {
    sprintf(cBuffer, "%s 11 symtab %04d - %08x (%08x==%08x)\n",
            buf, i,
            theSymTab[i].gid, theSymTab[i].adr.hdr, theSymTab[i].adr.ref);
    DDD_PrintDebug(cBuffer);
  }

  for(i=0; i<lenNewCpl; i++)
  {
    sprintf(cBuffer, "%s 12 newcpl %04d - %08x %4d %4d\n",
            buf, i,
            theNewCpl[i].gid, theNewCpl[i].dest, theNewCpl[i].prio);
    DDD_PrintDebug(cBuffer);
  }

  for(i=0; i<lenOldCpl; i++)
  {
    sprintf(cBuffer, "%s 13 oldcpl %04d - %08x %4d %4d\n",
            buf, i,
            theOldCpl[i].gid, theOldCpl[i].proc, theOldCpl[i].prio);
    DDD_PrintDebug(cBuffer);
  }

#       ifdef DebugAllPointers
  XferPtr(xm, buf);
#       endif
}
