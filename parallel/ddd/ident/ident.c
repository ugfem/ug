// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ident.c                                                       */
/*                                                                          */
/* Purpose:   object identification for ddd module                          */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/30 kb  begin                                            */
/*            94/04/25 kb  major revision, all Identfy-functions impl.      */
/*            96/02/08 kb  fixed bug in vchannel handling (due to T3D)      */
/*            96/11/25 kb  added indirect identification                    */
/*                                                                          */
/* Remarks:                                                                 */
/*                                                                          */
/****************************************************************************/

#define DebugIdent 10  /* 10 is off */


/****************************************************************************/


/*
        in debuglevel DebugIdentCons, addititional data is sent with
        the identify-messages, in order to check the consistency of
        the identification tupels.

        this is for configuring the debug actions only.
 */
#define DebugIdentCons  8



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

#include "basic/notify.h"


#include "basic/oopp.h"    /* for object-orientated style via preprocessor */



/****************************************************************************/
/*                                                                          */
/* defines in the following order                                           */
/*                                                                          */
/*        compile time constants defining static data size (i.e. arrays)    */
/*        other constants                                                   */
/*        macros                                                            */
/*                                                                          */
/****************************************************************************/


/* types of IDENTINFO items, ID_OBJECT must be the smallest value!! */
#define ID_OBJECT   1
#define ID_NUMBER   2
#define ID_STRING   3


#define TUPEL_LEN(t)    ((int)((t)&0x3f))


/* overall mode of identification */
enum IdentMode {
  IMODE_IDLE = 0,          /* waiting for next DDD_IdentifyBegin() */
  IMODE_CMDS,              /* after DDD_IdentifyBegin(), before DDD_IdentifyEnd() */
  IMODE_BUSY               /* during DDD_IdentifyEnd() */
};



/***/

/* some macros for customizing oopp */
#define _NEWPARAMS
#define _NEWPARAMS_OR_VOID    void

#define __INDENT(n)   { int i; for(i=0; i<n; i++) fputs("   ",fp);}
#define _PRINTPARAMS  , int indent, FILE *fp
#define _PRINTPARAMS_DEFAULT  ,0,stdout
#define _INDENT       __INDENT(indent)
#define _INDENT1      __INDENT(indent+1)
#define _INDENT2      __INDENT(indent+2)
#define _PRINTNEXT    , indent+1, fp
#define _PRINTSAME    , indent, fp


/* map memory allocation calls */
#define OO_Allocate  ident_AllocTmp
#define OO_Free      ident_FreeTmp


/* extra prefix for all xfer-related data structures and/or typedefs */
/*#define ClassPrefix*/



/*
    NOTE: all container-classes from ooppcc.h are implemented in this
          source file by setting the following define.
 */
#define ContainerImplementation
#define _CHECKALLOC(ptr)   assert(ptr!=NULL)


/***/


static int TmpMem_kind = TMEM_ANY;

static void *ident_AllocTmp (size_t size)
{
  return AllocTmpReq(size, TmpMem_kind);
}

static void ident_FreeTmp (void *buffer)
{
  FreeTmpReq(buffer, 0, TmpMem_kind);
}




/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/


/****************************************************************************/
/* IDENTIFIER:                                                              */
/****************************************************************************/

typedef union {
  int number;
  char       *string;
  DDD_GID object;
} IDENTIFIER;



/****************************************************************************/
/* MSGITEM:                                                                 */
/****************************************************************************/

typedef struct {
  DDD_GID gid;
  DDD_PRIO prio;

#       if DebugIdent<=DebugIdentCons
  unsigned long tupel;         /* send tupel ID for checking consistency */
#       endif

} MSGITEM;



/****************************************************************************/
/* IDENTINFO:                                                               */
/****************************************************************************/


typedef struct _IDENTINFO {
  int typeId;
  int entry;
  IDENTIFIER id;


  MSGITEM msg;                     /* this item is sent to other procs */

  DDD_HDR hdr;

  struct _ID_TUPEL  *tupel;
} IDENTINFO;



/****************************************************************************/
/* ID_TUPEL:                                                                */
/****************************************************************************/

typedef struct _ID_REFDBY
{
  struct _IDENTINFO  *by;
  struct _ID_REFDBY  *next;
} ID_REFDBY;


typedef struct _ID_TUPEL {
  unsigned long tId;
  IDENTINFO     **infos;

  int nObjIds;                    /* number of entries with typeID==ID_OBJECT */

  int loi;                        /* level of indirection */
  ID_REFDBY     *refd;            /* list of referencing IdEntries */

} ID_TUPEL;


/****************************************************************************/
/* IdEntry:                                                                 */
/****************************************************************************/


typedef struct _IdEntry {
  IDENTINFO msg;
} IdEntry;


/* define container class */
#define SegmListOf   IdEntry
#define SegmSize     128
#include "basic/ooppcc.h"




/****************************************************************************/
/* ID_PLIST:                                                                */
/****************************************************************************/

typedef struct _ID_PLIST {
  DDD_PROC proc;
  int nEntries;
  int nIdentObjs;

  struct _ID_PLIST *next;
  IdEntrySegmList  *entries;

  IDENTINFO   **local_ids;
  ID_TUPEL    *indexmap;                  /* index-mapping of local_ids array */

  MSGITEM     *msgin, *msgout;
  msgid idin, idout;
} ID_PLIST;


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


static ID_PLIST   *thePLists;
static int cntIdents, nPLists;

static int identMode;



/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*
        management functions for IdentMode.

        these functions control the mode the ident-module is
        currently in. this is used for error detection, but
        also for correct detection of coupling inconsistencies
        and recovery.
 */

static char *IdentModeName (int mode)
{
  switch(mode)
  {
  case IMODE_IDLE : return "idle-mode";
  case IMODE_CMDS : return "commands-mode";
  case IMODE_BUSY : return "busy-mode";
  }
  return "unknown-mode";
}


static void IdentSetMode (int mode)
{
  identMode = mode;

#       if DebugIdent<=8
  sprintf(cBuffer, "%4d: IdentMode=%s.\n",
          me, IdentModeName(identMode));
  DDD_PrintDebug(cBuffer);
#       endif
}


static int IdentSuccMode (int mode)
{
  switch(mode)
  {
  case IMODE_IDLE : return IMODE_CMDS;
  case IMODE_CMDS : return IMODE_BUSY;
  case IMODE_BUSY : return IMODE_IDLE;
  }
  return IMODE_IDLE;
}



/*
   static int IdentMode (void)
   {
        return identMode;
   }
 */


static int IdentActive (void)
{
  return identMode!=IMODE_IDLE;
}


static int IdentStepMode (int old)
{
  if (identMode!=old)
  {
    sprintf(cBuffer, ERR_ID_WRONG_MODE,
            IdentModeName(identMode), IdentModeName(old));
    DDD_PrintError('E', 3070, cBuffer);
    return FALSE;
  }

  IdentSetMode(IdentSuccMode(identMode));
  return TRUE;
}


/****************************************************************************/


static void PrintPList (ID_PLIST *plist)
{
  sprintf(cBuffer, "%d: PList proc=%04d entries=%05d\n",
          me, plist->proc, plist->nEntries);
  DDD_PrintDebug(cBuffer);
}



/****************************************************************************/


static int compareId (IDENTINFO *el1, IDENTINFO *el2)
{
  int cmp;

  /* first compare id type (NUMBER, STRING or OBJECT) */
  if (el1->typeId < el2->typeId) return(-1);
  if (el1->typeId > el2->typeId) return(1);


  /* same typeIds, compare identification then */
  switch (el1->typeId) {
  case ID_NUMBER :
    if (el1->id.number < el2->id.number) return(-1);
    if (el1->id.number > el2->id.number) return(1);
    break;

  case ID_STRING :
    cmp = strcmp(el1->id.string, el2->id.string);
    if (cmp!=0) return(cmp);
    break;

  case ID_OBJECT :
    if (el1->id.object < el2->id.object) return(-1);
    if (el1->id.object > el2->id.object) return(1);
    break;
  }

  return(0);
}



/****************************************************************************/

/*
        two functions for sorting IdentifyXXX-requests into tupels:

        sort_intoTupelsLists keeps order of IdentifyXXX-issueing by
           application program, i.e., the ordering is relevant

        sort_intoTupelsSets reorders the IdentifyXXX-items inside each
           tupel itself;
           at this level the ordering is done only by typeId, where
           ID_OBJECT comes first. lateron the IdentifyObject-items will
           be sorted according to their gid (for objects with loi==0)
           or the index of the loi-1 object (for objects with loi>0)
 */


static int sort_intoTupelsLists (const void *e1, const void *e2)
{
  IDENTINFO       *el1, *el2;

  el1 = *((IDENTINFO **)e1);
  el2 = *((IDENTINFO **)e2);


  /* sort according to (old) global ids */
  if (el1->msg.gid < el2->msg.gid) return(-1);
  if (el1->msg.gid > el2->msg.gid) return(1);

  /* if equal, keep ordering of input  */
  if (el1->entry < el2->entry) return(-1);
  if (el1->entry > el2->entry) return(1);

  return(0);
}


static int sort_intoTupelsSets (const void *e1, const void *e2)
{
  IDENTINFO       *el1, *el2;

  el1 = *((IDENTINFO **)e1);
  el2 = *((IDENTINFO **)e2);


  /* sort according to (old) global ids */
  if (el1->msg.gid < el2->msg.gid) return(-1);
  if (el1->msg.gid > el2->msg.gid) return(1);

  /* if equal, sort according to identificator itself */
  return (compareId(el1,el2));
}




/****************************************************************************/


static int sort_loi (const void *e1, const void *e2)
{
  ID_TUPEL        *el1, *el2;

  el1 = (ID_TUPEL *) e1;
  el2 = (ID_TUPEL *) e2;

  /* sort according to level-of-indirection */
  if (el1->loi < el2->loi) return(-1);
  if (el1->loi > el2->loi) return(1);

  return(0);
}

static int sort_tupelOrder (const void *e1, const void *e2)
{
  ID_TUPEL *el1, *el2;
  int cmp, i, nIds;
  DDD_HDR el1hdr, el2hdr;

  el1 = (ID_TUPEL *) e1;
  el2 = (ID_TUPEL *) e2;


  /* sort according to tupel id */
  if (el1->tId < el2->tId) return(-1);
  if (el1->tId > el2->tId) return(1);

  /* ids are equal, sort according tupel value */

  /* recode tupel length from lowest 6 bits */
  nIds = TUPEL_LEN(el1->tId);


  /* compare until one tupel entry differs */
  for(i=0; i<nIds; i++) {
    if ((cmp=compareId(el1->infos[i], el2->infos[i])) != 0)
      return(cmp);
  }


  /* if tupels are equal by all means up to now, we
     sort according to DDD_TYPE of local object.
     hence, we can identify two pairs of local objects with
     the same tupel.
     this has to be ommitted if objects with different
     types should be identifiable. KB 960814
   */
  el1hdr = el1->infos[0]->hdr;
  el2hdr = el2->infos[0]->hdr;
  if (OBJ_TYPE(el1hdr) < OBJ_TYPE(el2hdr)) return(-1);
  if (OBJ_TYPE(el1hdr) > OBJ_TYPE(el2hdr)) return(1);


  if (el1hdr!=el2hdr)
  {
    sprintf(cBuffer, ERR_ID_SAME_TUPEL, OBJ_GID(el1hdr), OBJ_GID(el2hdr));
    DDD_PrintError('E', 3030, cBuffer);

    /*
       for(i=0; i<nIds; i++) {
            printf("%4d: tupel[%d]  %08x/%d  %08x/%d   (id/loi)\n",
                    me, i,
                    el1->infos[i]->id.object, el1->loi,
                    el2->infos[i]->id.object, el2->loi);
       }
     */

    HARD_EXIT;
  }

  return(0);
}


/****************************************************************************/

static void SetLOI (IDENTINFO *ii, int loi)
{
  ID_REFDBY *rby;
  ID_TUPEL  *tupel = ii->tupel;
  int i;

  /*
     printf("%4d: %08x SetLOI(%d, %d)\n", me, ii->msg.gid, loi, ii->loi);
   */

  /* set loi to maximum of current and new value */
  tupel->loi = MAX(loi, tupel->loi);

  /* primitive cycle detection */
  if (tupel->loi > 64)
  {
    sprintf(cBuffer, ERR_ID_OBJ_CYCLE, ii->msg.gid, ii->id.object);
    DDD_PrintError('E', 3310, cBuffer);
    HARD_EXIT;
  }


  for(rby=tupel->refd; rby!=NULL; rby=rby->next)
  {
    SetLOI(rby->by, loi+1);

    /* TODO detection of cycles */
  }
}



static int sort_refd_gid (const void *e1, const void *e2)
{
  DDD_GID g1, g2;

  g1 = (*(IDENTINFO **) e1)->id.object;
  g2 = (*(IDENTINFO **) e2)->id.object;

  /* sort according to global id */
  if (g1 < g2) return(-1);
  if (g1 > g2) return(1);

  return(0);
}



static void ResolveDependencies (
  ID_TUPEL  *tupels, int nTupels,
  IDENTINFO **id, int nIds, int nIdentObjs)
{
  IDENTINFO **refd;
  int i, j;

  if (nIdentObjs==0)
    return;

  refd = (IDENTINFO **) AllocTmp(sizeof(IDENTINFO *)*nIdentObjs);
  if (refd==NULL) {
    DDD_PrintError('E', 3300, ERR_ID_NOMEM_RESOLV);
    return;
  }

  /* build array of pointers to objects being used for identification */
  for(i=0, j=0; i<nIds; i++)
  {
    if (id[i]->typeId==ID_OBJECT)
    {
      refd[j] = id[i];
      j++;
    }
  }

  /* sort it according to GID of referenced objects */
  qsort(refd, nIdentObjs, sizeof(IDENTINFO *), sort_refd_gid);


  /*
     for(j=0; j<nIdentObjs; j++)
          printf("%4d: DepObj %08x  (from %08x)\n",
                  me, refd[j]->id.object, OBJ_GID(refd[j]->hdr));
   */


  for(i=0, j=0; i<nTupels; i++)
  {
    while (j<nIdentObjs &&
           refd[j]->id.object < tupels[i].infos[0]->msg.gid)
      j++;

    while (j<nIdentObjs &&
           refd[j]->id.object == tupels[i].infos[0]->msg.gid)
    {
      ID_REFDBY *rby = (ID_REFDBY *)AllocTmpReq(sizeof(ID_REFDBY),TMEM_IDENT);
      if (rby==NULL) {
        DDD_PrintError('E', 3301, ERR_ID_NOMEM_RESOLV);
        return;
      }

      /* remember that idp[i] is referenced by refd[j] */
      rby->by        = refd[j];
      rby->next      = tupels[i].refd;
      tupels[i].refd = rby;

      j++;
    }
  }

  FreeTmp(refd,0);


  for(i=0; i<nTupels; i++)
  {
    ID_REFDBY *rby;

    for(rby=tupels[i].refd; rby!=NULL; rby=rby->next)
    {
      /* if loi>0, this subtree has been loi-ed before */
      if (tupels[i].loi==0)
        SetLOI(rby->by, tupels[i].loi+1);
    }
  }


#       if DebugIdent<=2
  /* display */
  for(i=0; i<nTupels; i++)
  {
    ID_REFDBY *rby;

    printf("%4d: %08x has loi %d\n",
           me, tupels[i].infos[0]->msg.gid, tupels[i].loi);

    for(rby=tupels[i].refd; rby!=NULL; rby=rby->next)
    {
      printf("%4d: %08x referenced by %08x\n",
             me, tupels[i].infos[0]->msg.gid, rby->by->msg.gid);
    }
  }
#       endif
}



static void CleanupLOI (ID_TUPEL *tupels, int nTupels)
{
  int i;

  for(i=0; i<nTupels; i++)
  {
    ID_REFDBY *rby, *next=0;

    for(rby=tupels[i].refd; rby!=NULL; rby=next)
    {
      next = rby->next;

      /* TODO use freelists */
      FreeTmpReq(rby, sizeof(ID_REFDBY), TMEM_IDENT);
    }
  }
}


/****************************************************************************/


/*
        tupel-id doesn't contain information about the data in
        the tupel, it does only contain information about the
        structure of a tupel!
 */
static void TupelInit (ID_TUPEL *tupel, IDENTINFO **id, int nIds)
{
  int i, nObjIds;
  unsigned long tId;

  /* init tupel auxiliary data */
  tupel->loi  = 0;
  tupel->refd = NULL;


  /* compute tupel id */
  tId = 0;
  nObjIds = 0;
  for(i=0; i<nIds; i++)
  {
    tId = (tId<<2) | id[i]->typeId;

    /* count entries with ID_OBJECT */
    if (id[i]->typeId==ID_OBJECT)
      nObjIds++;
  }


  /* code length of tupel into lowest 6 bits */
  tId = (tId<<6) | nIds;

  /* printf("%4d: compute tupel id = %08x\n", me, tId); */


  /* insert tupel id, number of ID_OBJECT-entries, and link to
     array of pointers to the tupel's IDENTINFO structs */
  tupel->tId     = tId;
  tupel->nObjIds = nObjIds;
  tupel->infos   = id;


  /* set link from IDENTINFOs to tupel */
  for(i=0; i<nIds; i++)
  {
    id[i]->tupel = tupel;
  }
}



static int IdentifySort (IDENTINFO **id, int nIds,
                         int nIdentObjs, MSGITEM *items_out, ID_TUPEL **indexmap_out,
#ifdef CPP_FRONTEND
                         DDD_PROC /*dest*/
#else
                         DDD_PROC dest
#endif
                         )
{
  ID_TUPEL *tupels;
  int i, j, last, nTupels;
  int keep_order_inside_tupel;


  /* sort to recognize identification tupels */
  /* in case of IDMODE_LISTS, the original ordering
     inside each tupel is kept. for IDMODE_SETS, each tupel
     is sorted according to the identificators themselves. */
  STAT_RESET3;
  switch (DDD_GetOption(OPT_IDENTIFY_MODE))
  {
  case IDMODE_LISTS :
    qsort(id, nIds, sizeof(IDENTINFO *), sort_intoTupelsLists);
    keep_order_inside_tupel = TRUE;
    break;

  case IDMODE_SETS :
    qsort(id, nIds, sizeof(IDENTINFO *), sort_intoTupelsSets);
    keep_order_inside_tupel = FALSE;
    break;

  default :
    DDD_PrintError('E', 3330, ERR_ID_UNKNOWN_OPT);
    HARD_EXIT;
  }
  STAT_INCTIMER3(T_QSORT_TUPEL);


  /* compute number of tupels and allocate tupel array */
  for(i=0, last=0, nTupels=1; i<nIds; i++)
  {
    if (id[i]->msg.gid > id[last]->msg.gid)
    {
      nTupels++;
      last=i;
    }
  }
  tupels = (ID_TUPEL *) AllocTmp(sizeof(ID_TUPEL)*nTupels);
  if (tupels==NULL) {
    DDD_PrintError('E', 3000, ERR_ID_NOMEM_SORT);
    return(0);
  }

  /* init tupels (e.g., compute tupel ids) */
  for(i=0, last=0, j=0; i<nIds; i++)
  {
    if (id[i]->msg.gid > id[last]->msg.gid)
    {
      TupelInit(&(tupels[j]), &(id[last]), i-last);
      j++;
      last=i;
    }
  }
  TupelInit(&(tupels[j]), &(id[last]), nIds-last);


  /*
          now, 'tupels' is an array of identification tupels,
          sorted according the gid of the object the tupel has
          been specified for.

          i.e., in a more abstract way tupels is a list of object
          gids which will be identified.
   */

  /* resolve dependencies caused by IdentifyObject,
     and set level-of-indirection accordingly */
  STAT_RESET3;
  ResolveDependencies(tupels, nTupels, id, nIds, nIdentObjs);
  STAT_INCTIMER3(T_RESOLVE_DEP);


  /*
          the setting for loi is used for the next sorting procedure,
          first level of indirection comes first.
   */

  /* sort array for loi */
  STAT_RESET3;
  qsort(tupels, nTupels, sizeof(ID_TUPEL), sort_loi);
  STAT_INCTIMER3(T_QSORT_LOI);


  STAT_RESET3;
  i=0; j=0;
  do {
    while (j<nTupels && tupels[i].loi==tupels[j].loi)
    {
      /* reorder because of changes in id.object */
      if (! keep_order_inside_tupel)
      {
        qsort(tupels[j].infos,
              tupels[j].nObjIds,
              sizeof(IDENTINFO *),
              sort_intoTupelsSets);
      }
      j++;
    }

    /* sort sub-array for tupelId, tupelValue */
    if (j-i > 1)
      qsort(tupels+i, j-i, sizeof(ID_TUPEL), sort_tupelOrder);

    /* inherit index to tupels referencing this one */
    while (i<j)
    {
      ID_REFDBY *rby;

      for(rby=tupels[i].refd; rby!=NULL; rby=rby->next)
      {
        /* dont use gid of referenced object (because it
           will be known only after identification), but its
           position in the identification table instead! */

        /*
           printf("%4d: insertRef dest=%d loi=%d i=%d, %08x <- %08x\n",
           me, dest, tupels[i].loi, i, OBJ_GID(tupels[i].infos[0]->hdr),
           OBJ_GID(rby->by->hdr));
         */


        rby->by->id.object = i;

        /* if the ordering is not significant, we must reorder
           the tupel after this operation. (i.e., for IDMODE_SETS).
         */
      }

      i++;
    }
    /* now i==j */

  } while (i<nTupels);
  STAT_INCTIMER3(T_BUILD_GRAPH);


  /* construct array which will be sent actually */
  STAT_RESET3;
  for(j=0; j<nTupels; j++)
  {
    /*
       int k;
     */
#               if DebugIdent<=1
    printf("%4d: Ident dest=%d msg_idx[ %08x ] = %5d, loi=%d\n",
           me, dest, tupels[j].infos[0]->msg.gid, j, tupels[j].loi);
#               endif

    /*
       for(k=0;k<tupels[j].nObjIds;k++)
       {
            printf("%4d:               msg_idx %d %08x\n", me,
                    k, tupels[j].infos[k]->id.object);
       }
     */



    items_out[j] = tupels[j].infos[0]->msg;

#               if DebugIdent<=DebugIdentCons
    /* send additional data for cons-checking */
    items_out[j].tupel = tupels[j].tId;
#               endif
  }
  STAT_INCTIMER3(T_CONSTRUCT_ARRAY);


  CleanupLOI(tupels, nTupels);

  /* return indexmap table, in order to keep ordering of tupels */
  /* note: this array has to be freed in the calling function! */
  *indexmap_out = tupels;


  return(nTupels);
}



static int InitComm (int nPartners)
{
  ID_PLIST  *plist;
  int i, err;
  DDD_PROC  *partners = DDD_ProcArray();

  /* fill partner processor numbers into array */
  for(plist=thePLists, i=0; i<nPartners; i++, plist=plist->next)
    partners[i] = plist->proc;

  if (! IS_OK(DDD_GetChannels(nPartners)))
  {
    return(FALSE);
  }


  /* initiate asynchronous receives and sends */
  for(plist=thePLists; plist!=NULL; plist=plist->next)
  {
    long *len_adr;

    plist->idin = RecvASync(VCHAN_TO(plist->proc),
                            ((char *)plist->msgin) - sizeof(long),
                            sizeof(MSGITEM)*plist->nEntries + sizeof(long), &err);

    /* store number of entries at beginning of message */
    len_adr = (long *) (((char *)plist->msgout) - sizeof(long));
    *len_adr = plist->nEntries;
    plist->idout = SendASync(VCHAN_TO(plist->proc),
                             ((char *)plist->msgout) - sizeof(long),
                             sizeof(MSGITEM)*plist->nEntries + sizeof(long), &err);
  }

  return(TRUE);
}


/****************************************************************************/

/*
        this routine checks whether number of idents per proc is indeed
        pairwise consistent.
 */

static void idcons_CheckPairs (void)
{
  NOTIFY_DESC *msgs = DDD_NotifyBegin(nPLists);
  ID_PLIST        *plist;
  int i, j, nRecvs, err=FALSE;

  for(i=0, plist=thePLists; plist!=NULL; plist=plist->next, i++)
  {
    msgs[i].proc = plist->proc;
    msgs[i].size = plist->nEntries;
  }

  /* communicate */
  nRecvs = DDD_Notify();
  if (nRecvs==ERROR)
  {
    DDD_PrintError('E', 3907, ERR_ID_NOTIFY_FAILED);
    HARD_EXIT;
  }


  /* perform checking */
  for(plist=thePLists; plist!=NULL; plist=plist->next)
  {
    for(j=0; j<nRecvs; j++)
    {
      if (msgs[j].proc==plist->proc)
        break;
    }

    if (j==nRecvs)
    {
      sprintf(cBuffer, ERR_ID_DIFF_IDENT, plist->proc, plist->nEntries);
      DDD_PrintError('E', 3900, cBuffer);
      err=TRUE;
    }
    else
    {
      if (msgs[j].size!=plist->nEntries)
      {
        sprintf(cBuffer, ERR_ID_DIFF_N_IDENT,
                msgs[j].size, plist->proc, plist->nEntries);
        DDD_PrintError('E', 3901, cBuffer);
        err=TRUE;
      }
    }
  }

  DDD_NotifyEnd();

  if (err)
  {
    DDD_PrintError('E', 3908, ERR_ID_ERRORS);
    HARD_EXIT;
  }
  else
  {
    DDD_PrintError('W', 3909, ERR_ID_OK);
  }
}



/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* Function:  DDD_IdentifyEnd                                               */
/*                                                                          */
/****************************************************************************/

/**
        End of identication phase.
        This function starts the object identification process. After a call to
        this function (on all processors) all {\bf Identify}-commands since the
        last call to \funk{IdentifyBegin} are executed. This involves a set
        of local communications between the processors.
 */

#if defined(C_FRONTEND) || defined(F_FRONTEND)
void DDD_IdentifyEnd (void)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::IdentifyEnd (void)
#endif
{
  ID_PLIST        *plist, *pnext=NULL;
  IdEntry *id;
  int cnt, j;

  /* REMARK: dont use the id->msg.msg.prio fields until they
          are explicitly set at line L1!                         */

  STAT_SET_MODULE(DDD_MODULE_IDENT);
  STAT_ZEROALL;

#       if DebugIdent<=9
  printf("%4d: DDD_IdentifyEnd.\n", me);
  fflush(stdout);
#       endif

  /* step mode and check whether call to IdentifyEnd is valid */
  if (!IdentStepMode(IMODE_CMDS))
  {
    DDD_PrintError('E', 3071, ERR_ID_ABORT_END);
    HARD_EXIT;
  };


#   if DebugIdent<=9
  idcons_CheckPairs();
#       endif


  STAT_RESET1;

  /* for each id_plist entry */
  for(plist=thePLists, cnt=0; plist!=NULL; plist=plist->next, cnt++)
  {
    /* allocate message buffers */
    /* use one alloc for three buffers */
    plist->local_ids = (IDENTINFO **) AllocTmp(
      sizeof(IDENTINFO *)*plist->nEntries +                    /* for local id-infos  */
      sizeof(long) +                                           /* len of incoming msg */
      sizeof(MSGITEM)    *plist->nEntries +                    /* for incoming msg    */
      sizeof(long) +                                           /* len of outgoing msg */
      sizeof(MSGITEM)    *plist->nEntries                      /* for outgoing msg    */
      );

    if (plist->local_ids==NULL)
    {
      DDD_PrintError('F',3100, ERR_ID_NOMEM_IDENT_END);
      HARD_EXIT;
    }
    plist->msgin  = (MSGITEM *)
                    (((char *)&plist->local_ids[plist->nEntries]) + sizeof(long));
    plist->msgout = (MSGITEM *)
                    (((char *)&plist->msgin[plist->nEntries]) + sizeof(long));


    /* construct pointer array to IDENTINFO structs */
    /* AND: fill in current priority from object's header */
    {
      IdEntrySegm *li;
      int i = 0;

      for(li=plist->entries->first; li!=NULL; li=li->next)
      {
        int entry;
        for(entry=0; entry<li->nItems; entry++)
        {
          IdEntry *id = &(li->data[entry]);

          plist->local_ids[i] = &(id->msg);
          id->msg.msg.prio = OBJ_PRIO(id->msg.hdr);                               /* L1 */
          i++;
        }
      }
    }


    /* sort outgoing items */
    STAT_RESET2;
    plist->nEntries = IdentifySort(plist->local_ids, plist->nEntries,
                                   plist->nIdentObjs,
                                   plist->msgout,  /* output: msgbuffer outgoing */
                                   &plist->indexmap, /* output: mapping of indices to local_ids array */
                                   plist->proc);
    STAT_INCTIMER2(T_PREPARE_SORT);


#               if DebugIdent<=5
    PrintPList(plist);
#               endif
  }
  STAT_TIMER1(T_PREPARE);
  STAT_SETCOUNT(N_PARTNERS, cnt);

  /* initiate comm-channels and send/receive calls */
  STAT_RESET1;
  if (!InitComm(cnt))
  {
    DDD_PrintError('E', 3074, ERR_ID_ABORT_END);
    HARD_EXIT;
  }


  /*
          each pair of processors now has a plist with one copy on
          each side. the actual OBJ_GID is computed as the minimum
          of the two local object ids on each processor.
   */

#       if DebugIdent<=4
  printf("%4d: DDD_IdentifyEnd. PLists ready.\n", me); fflush(stdout);
#       endif


  /* poll receives */
  for(plist=thePLists, j=0; j<cnt; )
  {
    if (plist->msgin!=NULL)
    {
      int ret, i;

      if ((ret=InfoARecv(VCHAN_TO(plist->proc), plist->idin))==1)
      {
        /* process single plist */
        MSGITEM   *msgin  = plist->msgin;
        ID_TUPEL  *msgout = plist->indexmap;

        /* check control data */
        long *len_adr = (long *) (((char *)msgin) - sizeof(long));
        if (*len_adr != plist->nEntries)
        {
          sprintf(cBuffer, ERR_ID_DIFF_N_OBJECTS,
                  (int)*len_adr, plist->proc, plist->nEntries);
          DDD_PrintError('E', 3902, cBuffer);
          HARD_EXIT;
        }

        for(i=0; i<plist->nEntries; i++, msgin++, msgout++)
        {
#                                       if DebugIdent<=1
          printf("%4d: identifying %08x with %08x/%d to %08x\n", me,
                 OBJ_GID(msgout->infos[0]->hdr), msgin->gid,
                 plist->proc,
                 MIN(OBJ_GID(msgout->infos[0]->hdr), msgin->gid));
#                                       endif

#                                       if DebugIdent<=DebugIdentCons
          if (msgout->tId != msgin->tupel)
          {
            sprintf(cBuffer, ERR_ID_INCONS_TUPELS,
                    OBJ_GID(msgout->infos[0]->hdr), me,
                    msgin->gid, plist->proc);
            DDD_PrintError('E', 3920, cBuffer);
            HARD_EXIT;
          }
#                                       endif

          /* compute new GID from minimum of both current GIDs */
          OBJ_GID(msgout->infos[0]->hdr) =
            MIN(OBJ_GID(msgout->infos[0]->hdr), msgin->gid);

          /* add a coupling for new object copy */
          AddCoupling(msgout->infos[0]->hdr, plist->proc, msgin->prio);
        }

        /* free indexmap (=tupel) array */
        FreeTmp(plist->indexmap,0);

        /* mark plist as finished */
        plist->msgin=NULL;
        j++;
      }
      else
      {
        if (ret==-1)
        {
          sprintf(cBuffer, ERR_ID_CANT_RECV, plist->proc);
          DDD_PrintError('E', 3921, cBuffer);
          HARD_EXIT;
        }
      }
    }

    /* next plist, perhaps restart */
    plist=plist->next; if (plist==NULL) plist=thePLists;
  };
  STAT_TIMER1(T_COMM_AND_IDENT);

  /* poll sends */
  for(plist=thePLists; plist!=0; plist=pnext)
  {
    pnext = plist->next;

    /* wait for correct send and free buffer */
    while(InfoASend(VCHAN_TO(plist->proc), plist->idout)!=1)
      ;

    /* now, the plist->entries list isn't needed anymore, free */
    IdEntrySegmList_Free(plist->entries);

    FreeTmp(plist->local_ids,0);
    FreeTmpReq(plist, sizeof(ID_PLIST), TMEM_IDENT);
  };


#       if DebugIdent<=8
  printf("%4d: DDD_IdentifyEnd. Rebuilding interfaces.\n", me);
  fflush(stdout);
#       endif



  /* rebuild interfaces after topological change */
  STAT_RESET1;
  IFAllFromScratch();
  STAT_TIMER1(T_BUILD_IF);


#       if DebugIdent<=9
  printf("%4d: DDD_IdentifyEnd. Ready.\n", me); fflush(stdout);
#       endif

  IdentStepMode(IMODE_BUSY);
}



/****************************************************************************/


static IdEntry *IdentifyIdEntry (DDD_HDR hdr, DDD_PROC proc)
{
  IdEntry     *id;
  ID_PLIST        *plist;

  /* check whether Identify-call is valid */
  if (!IdentActive())
  {
    DDD_PrintError('E', 3072, ERR_ID_NO_BEGIN);
    HARD_EXIT;
  }

  if (proc==me)
  {
    sprintf(cBuffer, ERR_ID_NOT_WITH_ME, OBJ_GID(hdr));
    DDD_PrintError('E', 3060, cBuffer);
    HARD_EXIT;
  }

  if (proc>=procs)
  {
    sprintf(cBuffer, ERR_ID_NOT_WITH_PROC, OBJ_GID(hdr), proc);
    DDD_PrintError('E', 3061, cBuffer);
    HARD_EXIT;
  }



  /* search current plist entries */
  for(plist=thePLists; plist!=NULL; plist=plist->next) {
    if (plist->proc==proc)
      break;
  }

  if (plist==NULL)
  {
    /* get new id_plist record */
    plist = (ID_PLIST *) AllocTmpReq(sizeof(ID_PLIST),TMEM_IDENT);
    if (plist==NULL) {
      DDD_PrintError('F', 3210, ERR_ID_NOMEM_IDENTRY);
      return;
    }

    plist->proc = proc;
    plist->nEntries = 0;
    plist->entries = New_IdEntrySegmList();
    plist->nIdentObjs = 0;
    plist->next = thePLists;
    thePLists = plist;
    nPLists++;
  }


  /* insert into current plist */
  id =  IdEntrySegmList_NewItem(plist->entries);
  plist->nEntries++;
  if (id->msg.typeId==ID_OBJECT)
  {
    plist->nIdentObjs++;
  }

  id->msg.hdr      = hdr;
  id->msg.msg.gid  = OBJ_GID(hdr);

  /* NOTE: priority can change between Identify-command and IdentifyEnd!
     therefore, scan priorities at the beginning of IdentifyEnd, and not
     here. KB 970730.
          id->msg.msg.prio = OBJ_PRIO(hdr);
   */

  id->msg.entry = cntIdents++;

  return(id);
}


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_IdentifyNumber                                            */
/*                                                                          */
/****************************************************************************/

/**
        DDD Object identification via integer number.
        After an initial call to \funk{IdentifyBegin}, this function
        identifies two object copies on separate processors. It has to be
        called on both processors with the same identification value.
        The necessary actions (e.g. message transfer) are executed via the
        final call to \funk{IdentifyEnd}; therefore a whole set of
        \funk{Identify}-operations is accumulated.

        After the identification both objects have the same DDD global
        object ID, which is build using the minimum of both local object IDs.

        The identification specified here may be detailed even further by
        additional calls to {\bf Identify}-operations with the same
        local object. This will construct an identification tupel from
        all {\bf Identify}-commands for this local object.

   @param hdr   DDD local object which has to be identified with another object
   @param proc  Owner (processor number) of corresponding local object. This processor has to issue a corresponding call to this {\bf DDD\_Identify}-function.
   @param ident Identification value. This is an arbitrary number to identify two corresponding operations on different processors.
 */

#ifdef C_FRONTEND
void DDD_IdentifyNumber (DDD_HDR hdr, DDD_PROC proc, int ident)
{
#endif
#ifdef CPP_FRONTEND
void DDD_Object::IdentifyNumber (DDD_PROC proc, int ident)
{
  DDD_HDR hdr = this;
#endif
#ifdef F_FRONTEND
void orgDDD_IdentifyNumber (DDD_HDR hdr, DDD_PROC proc, int ident);

void DDD_IdentifyNumber (DDD_TYPE *type, DDD_OBJ *obj, DDD_PROC *proc, int *ident)
{
  DDD_HDR hdr = OBJ2HDR(*obj,&(theTypeDefs[*type]));

  orgDDD_IdentifyNumber (hdr, *proc, *ident);
}

void orgDDD_IdentifyNumber (DDD_HDR hdr, DDD_PROC proc, int ident)
{
#endif
IdEntry *id;

        #if DebugIdent<=2
printf("%4d: IdentifyIdEntry %08x %02d with %4d num %d\n", me,
       OBJ_GID(hdr), OBJ_TYPE(hdr), proc, id->msg.id.number);
        #endif

id = IdentifyIdEntry(hdr, proc);
if (id==NULL) {
  DDD_PrintError('F', 3200, ERR_ID_NOMEM_IDNUMBER);
  return;
}

id->msg.typeId = ID_NUMBER;
id->msg.id.number = ident;
}


/****************************************************************************/
/*                                                                          */
/* Function:  DDD_IdentifyString                                            */
/*                                                                          */
/****************************************************************************/

/**
        DDD Object identification via character string.
        After an initial call to \funk{IdentifyBegin}, this function
        identifies two object copies on separate processors. It has to be
        called on both processors with the same identification string.
        The necessary actions (e.g. message transfer) are executed via the
        final call to \funk{IdentifyEnd}; therefore a whole set of
        \funk{Identify}-operations is accumulated.

        After the identification both objects have the same DDD global
        object ID, which is build using the minimum of both local object IDs.

        The identification specified here may be detailed even further by
        additional calls to {\bf Identify}-operations with the same
        local object. This will construct an identification tupel from
        all {\bf Identify}-commands for this local object.

   @param hdr   DDD local object which has to be identified with another object
   @param proc  Owner (processor number) of corresponding local object. This processor has to issue a corresponding call to this {\bf DDD\_Identify}-function.
   @param ident Identification value. This is an arbitrary string to identify two corresponding operations on different processors.
 */

#ifdef C_FRONTEND
void DDD_IdentifyString (DDD_HDR hdr, DDD_PROC proc, char *ident)
{
#endif
#ifdef CPP_FRONTEND
void DDD_Object::IdentifyString (DDD_PROC proc, char *ident)
{
  DDD_HDR hdr = this;
#endif
#ifdef F_FRONTEND
void orgDDD_IdentifyString (DDD_HDR hdr, DDD_PROC proc, char *ident);

void DDD_IdentifyString (DDD_TYPE *type, DDD_OBJ *obj, DDD_PROC *proc, char *ident)
{
  DDD_HDR hdr = OBJ2HDR(*obj,&(theTypeDefs[*type]));

  orgDDD_IdentifyString (hdr, *proc, ident);
}

void orgDDD_IdentifyString (DDD_HDR hdr, DDD_PROC proc, char *ident)
{
#endif
IdEntry *id;

        #if DebugIdent<=2
printf("%4d: IdentifyIdEntry %08x %02d with %4d str %s\n", me,
       OBJ_GID(hdr), OBJ_TYPE(hdr), proc, id->msg.id.string);
        #endif

id = IdentifyIdEntry(hdr, proc);
if (id==NULL) {
  DDD_PrintError('F', 3201, ERR_ID_NOMEM_IDSTRING);
  return;
}

id->msg.typeId = ID_STRING;
id->msg.id.string = ident;
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_IdentifyObject                                            */
/*                                                                          */
/****************************************************************************/

/**
        DDD Object identification via another DDD Object.
        After an initial call to \funk{IdentifyBegin}, this function
        identifies two object copies on separate processors. It has to be
        called on both processors with the same identification object.
        The necessary actions (e.g. message transfer) are executed via the
        final call to \funk{IdentifyEnd}; therefore a whole set of
        \funk{Identify}-operations is accumulated.

        After the identification both objects have the same DDD global
        object ID, which is build using the minimum of both local object IDs.

        The identification object {\em ident} must be either a distributed
        object known to both processors issueing the \funk{IdentifyObject}-command
        or a local object which is not known to these two processors, but which
        will also be identified during the current {\bf Identify}-process.

        The identification specified here may be detailed even further by
        additional calls to {\bf Identify}-operations with the same
        local object. This will construct an identification tupel from
        all {\bf Identify}-commands for this local object.

   @param hdr   DDD local object which has to be identified with another object
   @param proc  Owner (processor number) of corresponding local object. This processor has to issue a corresponding call to this {\bf DDD\_Identify}-function.
   @param ident Identification object. This is an arbitrary global object which is known to both processors involved to identify the two corresponding operations on these processors.
 */

#ifdef C_FRONTEND
void DDD_IdentifyObject (DDD_HDR hdr, DDD_PROC proc, DDD_HDR ident)
{
#endif
#ifdef CPP_FRONTEND
void DDD_Object::IdentifyObject (DDD_PROC proc, DDD_Object* ident)
{
  DDD_HDR hdr   = this;
#endif
#ifdef F_FRONTEND
void orgDDD_IdentifyObject (DDD_HDR hdr, DDD_PROC proc, DDD_HDR ident);

void DDD_IdentifyObject (DDD_TYPE *type, DDD_OBJ *obj, DDD_PROC *proc, DDD_TYPE *type2, DDD_OBJ *obj2)
{
  DDD_HDR hdr = OBJ2HDR(*obj,&(theTypeDefs[*type]));
  DDD_HDR ident = OBJ2HDR(*obj2,&(theTypeDefs[*type2]));

  orgDDD_IdentifyObject (hdr, *proc, ident);
}

void orgDDD_IdentifyObject (DDD_HDR hdr, DDD_PROC proc, DDD_HDR ident)
{
#endif
IdEntry *id;

        #if DebugIdent<=2
printf("%4d: IdentifyIdEntry %08x %02d with %4d gid %08x\n", me,
       OBJ_GID(hdr), OBJ_TYPE(hdr), proc, id->msg.id.object);
        #endif

id = IdentifyIdEntry(hdr, proc);
if (id==NULL) {
  DDD_PrintError('F', 3202, ERR_ID_NOMEM_IDOBJ);
  return;
}

id->msg.typeId = ID_OBJECT;

/* use OBJ_GID as estimate for identification value, this estimate
   might be replaced when the corresponding object is identified
   itself. then its index in the identify-message will be used.
   remember identification value in order to replace above estimate,
   if necessary (i.e., remember ptr to ddd-hdr) */
id->msg.id.object = OBJ_GID(ident);
}



/****************************************************************************/
/*                                                                          */
/* Function:  DDD_IdentifyBegin                                             */
/*                                                                          */
/****************************************************************************/

/**
        Begin identification phase.
        A call to this function establishes a global identification operation.
        It should be issued on all processors. After this call an arbitrary
        series of {\bf Identify}-commands may be issued. The global
        identification operation is carried out via a \funk{IdentifyEnd}
        call on each processor.

        All identification commands given for one local object will be collected
        into an {\em identification tupel}. Thus, object identificators can be
        constructed from several simple identification calls. DDD option
   #IDENTIFY_MODE# may be set before the \funk{IdentifyEnd} call
        in order to specify how the order of simple identificators is
        handled for each complex identification tupel:

        \begin{description}
        \item[#IDMODE_LISTS#:]%
        The order of all identification commands for one local object is kept.
        Both processors with corresponding complex identificators must issue
        the identification commands in the same order.
        %
        \item[#IDMODE_SETS#:]%
        The order of all identification commands for one local object is
        not relevant. The DDD identification module sorts the commands
        inside each complex identificator. Both processors with corresponding
        identification tupels may issue the identification commands in any
        order.
        \end{description}
 */

#if defined(C_FRONTEND) || defined(F_FRONTEND)
void DDD_IdentifyBegin (void)
#endif
#ifdef CPP_FRONTEND
void DDD_Library::IdentifyBegin (void)
#endif
{
  /* step mode and check whether call to IdentifyBegin is valid */
  if (!IdentStepMode(IMODE_IDLE))
  {
    DDD_PrintError('E', 3073, ERR_ID_ABORT_BEGIN);
    HARD_EXIT;
  }

  thePLists    = NULL;
  nPLists      = 0;
  cntIdents    = 0;
}


/****************************************************************************/


void ddd_IdentInit (void)
{
  IdentSetMode(IMODE_IDLE);
}


void ddd_IdentExit (void)
{}



/****************************************************************************/
