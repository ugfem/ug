// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ifcreate.c                                                    */
/*                                                                          */
/* Purpose:   routines concerning interfaces between processors             */
/*            part 1: creating and maintaining interfaces                   */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   93/11/30 kb  begin                                            */
/*            94/03/03 kb  complete rewrite                                 */
/*            94/09/12 kb  IFExchange & IFOneway rewrite, two bugs fixed    */
/*            94/09/21 kb  created from if.c                                */
/*            95/01/13 kb  added range functionality                        */
/*            96/01/08 kb  renamed range to attr                            */
/*            96/01/16 kb  added DDD_OBJ shortcut to avoid indirect addr.   */
/*            96/05/07 kb  removed need for global const MAX_COUPLING       */
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
#include <assert.h>

#include "dddi.h"
#include "if.h"




/****************************************************************************/
/*                                                                          */
/* definition of exported variables                                         */
/*                                                                          */
/****************************************************************************/

IF_DEF theIF[MAX_IF];
int nIFs;


/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only                  */
/*                                                                          */
/****************************************************************************/


/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)



static IF_PROC *memlistIFHead;
static IF_ATTR *memlistIFAttr;


/****************************************************************************/
/*                                                                          */
/* routines                                                                 */
/*                                                                          */
/****************************************************************************/


static IF_PROC *NewIFHead (void)
{
  IF_PROC *ifh;

  if (memlistIFHead==NULL)
  {
    ifh = (IF_PROC *) AllocIF(sizeof(IF_PROC));
  }
  else
  {
    ifh = memlistIFHead;
    memlistIFHead = ifh->next;
  }

  return(ifh);
}


static void DisposeIFHead (IF_PROC *ifh)
{
  ifh->next = memlistIFHead;
  memlistIFHead = ifh;
}



static IF_ATTR *NewIFAttr (void)
{
  IF_ATTR *ifr;

  if (memlistIFAttr==NULL)
  {
    ifr = (IF_ATTR *) AllocIF(sizeof(IF_ATTR));
  }
  else
  {
    ifr = memlistIFAttr;
    memlistIFAttr = ifr->next;
  }

  return(ifr);
}


static void DisposeIFAttr (IF_ATTR *ifr)
{
  ifr->next = memlistIFAttr;
  memlistIFAttr = ifr;
}




static int sort_type (const void *e1, const void *e2)
{
  if (*(DDD_TYPE *)e1 < *(DDD_TYPE *)e2) return(-1);
  if (*(DDD_TYPE *)e1 == *(DDD_TYPE *)e2) return(0);
  return(1);
}


static int sort_prio (const void *e1, const void *e2)
{
  if (*(DDD_PRIO *)e1 < *(DDD_PRIO *)e2) return(-1);
  if (*(DDD_PRIO *)e1 == *(DDD_PRIO *)e2) return(0);
  return(1);
}


/****************************************************************************/
/*                                                                          */
/* Function:  sort_IFCouplings                                              */
/*                                                                          */
/* Purpose:   qsort procedure for sorting interface couplings.              */
/*            coupling list will be ordered according to:                   */
/*                1. processor number of represented object copy            */
/*                    (increasing order)                                    */
/*                2. direction of interface according to priorities         */
/*                    (increasing order)                                    */
/*                3. attr property of objects                               */
/*                    (decreasing order)                                    */
/*                4. global ids of objects                                  */
/*                    (increasing order)                                    */
/*                                                                          */
/* Input:     two couplings                                                 */
/*                                                                          */
/* Output:    1, 0, -1 depending on order of couplings                      */
/*                                                                          */
/****************************************************************************/

static int sort_IFCouplings (const void *e1, const void *e2)
{
  COUPLING  *cp1, *cp2;
  INT gid1, gid2;
  DDD_ATTR attr1, attr2;

  cp1 = *((COUPLING **)e1);
  cp2 = *((COUPLING **)e2);

  if (cp1->proc < cp2->proc) return(-1);
  if (cp1->proc > cp2->proc) return(1);

  if (CPLDIR(cp1) < CPLDIR(cp2)) return(-1);
  if (CPLDIR(cp1) > CPLDIR(cp2)) return(1);

  attr1 = OBJ_ATTR(cp1->obj);
  attr2 = OBJ_ATTR(cp2->obj);
  if (attr1 > attr2) return(-1);
  if (attr1 < attr2) return(1);

  gid1 = OBJ_GID(cp1->obj);
  gid2 = OBJ_GID(cp2->obj);
  if (gid1 < gid2) return(-1);
  if (gid1 == gid2) return(0);
  return(1);
}



void IFDeleteAll (DDD_IF ifId)
{
  IF_PROC  *ifh, *ifhNext;
  IF_ATTR *ifr, *ifrNext;

  /* free IF_PROC memory */
  ifh=theIF[ifId].ifHead;
  while (ifh!=NULL)
  {
    ifhNext = ifh->next;

    /* free IF_ATTR memory */
    ifr=ifh->ifAttr;
    while (ifr!=NULL)
    {
      ifrNext = ifr->next;

      DisposeIFAttr(ifr);
      ifr = ifrNext;
    }


    /* if there are msg-buffers, then we must free them here */
    BufferFree(ifh->bufIn);
    BufferFree(ifh->bufOut);

    DisposeIFHead(ifh);

    ifh = ifhNext;
  }

  /* free memory for coupling table */
  if (theIF[ifId].cpl!=NULL)
  {
    FreeIF(theIF[ifId].cpl);
    theIF[ifId].cpl=NULL;
  }

  /* free memory for shortcut object table */
  if (theIF[ifId].obj!=NULL)
  {
    FreeIF(theIF[ifId].obj);
    theIF[ifId].obj=NULL;
  }

  /* reset pointers */
  theIF[ifId].ifHead = NULL;
  theIF[ifId].nIfHeads = 0;
}




/* TODO  el-set relation, VERY inefficient! */
static int is_elem (DDD_PRIO el, int n, DDD_PRIO *set)
{
  REGISTER int i;

  for(i=0; i<n; i++)
    if (set[i]==el)
      return(TRUE);

  return(FALSE);
}



static void update_channels (DDD_IF ifId)
{
  IF_PROC *ifh;
  int i;
  DDD_PROC *partners = DDD_ProcArray();

  if (theIF[ifId].nIfHeads==0)
    return;

  /* MarkHeap(); */

  for(i=0, ifh=theIF[ifId].ifHead; ifh!=NULL; i++, ifh=ifh->next)
  {
    partners[i] = ifh->proc;
  }

  DDD_GetChannels(theIF[ifId].nIfHeads);

  for(ifh=theIF[ifId].ifHead; ifh!=NULL; ifh=ifh->next) {
    ifh->vc = VCHAN_TO(ifh->proc);
  }

  /* ReleaseHeap(); */
}


/****************************************************************************/

/* collect couplings into interface array, for standard interface */

static COUPLING ** IFCollectStdCouplings (void)
{
  COUPLING **cplarray;
  int index, n;

  if (nCplItems==0)
  {
    return(NULL);
  }

  /* get memory for couplings inside STD_IF */
  cplarray = (COUPLING **) AllocIF(sizeof(COUPLING *)*nCplItems);
  if (cplarray==NULL) {
    DDD_PrintError('E', 4000, STR_NOMEM " in IFCreateFromScratch");
    HARD_EXIT;
  }

  /* collect couplings */
  n=0;
  for(index=0; index<NCPL_GET; index++)
  {
    COUPLING  *cpl;

    for(cpl=IdxCplList(index); cpl!=NULL; cpl=CPL_NEXT(cpl))
    {
      cplarray[n] = cpl;
      n++;
    }
  }
  /*
     printf("%04d: n=%d, nCplItems=%d\n",me,n,nCplItems);
   */
  assert(n==nCplItems);
  return(cplarray);
}


/****************************************************************************/

static void IFCreateFromScratch (COUPLING **tmpcpl, DDD_IF ifId)
{
  IF_PROC     *ifHead, *lastIfHead;
  IF_ATTR    *ifAttr, *lastIfAttr;
  int n, i;
  DDD_PROC lastproc;
  int STAT_MOD;

  STAT_GET_MODULE(STAT_MOD);
  STAT_SET_MODULE(DDD_MODULE_IF);

  /* first delete possible old interface */
  IFDeleteAll(ifId);

  STAT_RESET1;
  if (ifId==STD_INTERFACE)
  {
    theIF[ifId].cpl = IFCollectStdCouplings();
    n = nCplItems;
  }
  else
  {
    int index;


    /* collect relevant couplings into tmpcpl array */
    n=0;
    for(index=0; index<NCPL_GET; index++)
    {
      /* determine whether object belongs to IF */
      if ((1<<OBJ_TYPE(ddd_ObjTable[index])) & theIF[ifId].maskO)
      {
        int objInA, objInB;

        objInA = is_elem(OBJ_PRIO(ddd_ObjTable[index]),
                         theIF[ifId].nPrioA, theIF[ifId].A);
        objInB = is_elem(OBJ_PRIO(ddd_ObjTable[index]),
                         theIF[ifId].nPrioB, theIF[ifId].B);

        if (objInA || objInB)
        {
          COUPLING  *cpl;

          /* test coupling list */
          for(cpl=IdxCplList(index); cpl!=NULL; cpl=CPL_NEXT(cpl))
          {
            int cplInA, cplInB, dir;

            cplInA = is_elem(cpl->prio,
                             theIF[ifId].nPrioA, theIF[ifId].A);
            cplInB = is_elem(cpl->prio,
                             theIF[ifId].nPrioB, theIF[ifId].B);

            /* compute possible IF directions */
            dir = ((objInA&&cplInB) ? DirAB : 0) |
                  ((objInB&&cplInA) ? DirBA : 0);

            if (dir > 0)
            {
              SETCPLDIR(cpl,dir);
              tmpcpl[n] = cpl;
              n++;
            }
          }
        }
      }
    }

    if (n>0)
    {
      /* re-alloc cpllist, now with correct size */
      theIF[ifId].cpl = (COUPLING **) AllocIF(sizeof(COUPLING *)*n);
      if (theIF[ifId].cpl==NULL)
      {
        sprintf(cBuffer, STR_NOMEM
                " for IF %02d in IFCreateFromScratch", ifId);
        DDD_PrintError('E', 4001, cBuffer);
        HARD_EXIT;
      }

      /* copy data from temporary array */
      memcpy((void *)theIF[ifId].cpl, (void *)tmpcpl, sizeof(COUPLING *)*n);
    }
    else
    {
      theIF[ifId].cpl = NULL;
    }
  }
  STAT_TIMER1(T_CREATE_COLLECT);



  /* sort IF couplings */
  STAT_RESET1;
  if (n>1)
    qsort(theIF[ifId].cpl, n, sizeof(COUPLING *), sort_IFCouplings);
  STAT_TIMER1(T_CREATE_SORT);


  /* create IF_PROCs */
  STAT_RESET1;
  lastproc = PROC_INVALID;
  lastIfHead  = NULL;
  theIF[ifId].nIfHeads = 0;
  for(i=0; i<n; i++)
  {
    COUPLING  **cplp = &(theIF[ifId].cpl[i]);
    COUPLING  *cpl = *cplp;
    DDD_ATTR attr = OBJ_ATTR(cpl->obj);

    if (cpl->proc != lastproc)
    {
      /* create new IfHead */
      theIF[ifId].nIfHeads++;
      ifHead = NewIFHead();
      ifHead->nItems   = 0;
      ifHead->cpl      = cplp;
      ifHead->obj      = NULL;
      ifHead->nAB      = ifHead->nBA   = ifHead->nABA   = 0;
      ifHead->cplAB    = ifHead->cplBA = ifHead->cplABA = NULL;
      ifHead->proc     = cpl->proc;
      ifHead->next     = lastIfHead;
      lastIfHead = ifHead;
      lastproc   = ifHead->proc;

      BufferInit(ifHead->bufIn);
      BufferInit(ifHead->bufOut);

      ifHead->nAttrs = 1;
      ifHead->ifAttr = ifAttr = NewIFAttr();
      ifAttr->attr   = attr;
      ifAttr->nItems = 0;
      ifAttr->nAB    = ifAttr->nBA   = ifAttr->nABA   = 0;
      ifAttr->cplAB  = ifAttr->cplBA = ifAttr->cplABA = NULL;
      ifAttr->next   = NULL;
      lastIfAttr = ifAttr;
    }

    /* count #items per processor */
    ifHead->nItems++;


    /* keep current ifAttr or find new one? */
    if (attr!=ifAttr->attr)
    {
      IF_ATTR *ifR;

      /* does ifAttr already exist? */
      ifR = ifHead->ifAttr;
      while ((ifR!=NULL) && (ifR->attr!=attr))
      {
        ifR=ifR->next;
      }

      if (ifR!=NULL)
      {
        /* reuse existing ifAttr */
        ifAttr = ifR;
      }
      else
      {
        /* create new ifAttr */
        ifHead->nAttrs++;
        ifAttr = NewIFAttr();
        ifAttr->attr   = attr;
        ifAttr->nItems = 0;
        ifAttr->nAB    = ifAttr->nBA   = ifAttr->nABA   = 0;
        ifAttr->cplAB  = ifAttr->cplBA = ifAttr->cplABA = NULL;
        ifAttr->next   = NULL;
        lastIfAttr->next = ifAttr;
        lastIfAttr = ifAttr;
      }
    }


    /* count #items per processor and attr */
    ifAttr->nItems++;


    /* count #items per directions AB, BA or ABA
            and set beginnings of AB/BA/ABA subarrays */
    if (ifId!=STD_INTERFACE)
    {
      switch (CPLDIR(cpl))
      {
      case DirAB :
        ifHead->nAB++;
        if (ifHead->cplAB==0) ifHead->cplAB = cplp;
        ifAttr->nAB++;
        if (ifAttr->cplAB==0) ifAttr->cplAB = cplp;
        break;

      case DirBA :
        ifHead->nBA++;
        if (ifHead->cplBA==0) ifHead->cplBA = cplp;
        ifAttr->nBA++;
        if (ifAttr->cplBA==0) ifAttr->cplBA = cplp;
        break;

      case DirABA :
        ifHead->nABA++;
        if (ifHead->cplABA==0) ifHead->cplABA = cplp;
        ifAttr->nABA++;
        if (ifAttr->cplABA==0) ifAttr->cplABA = cplp;
        break;
      }
    }
  }
  STAT_TIMER1(T_CREATE_BUILD);

  /* remember anchor of ifHead list */
  if (theIF[ifId].nIfHeads>0) {
    theIF[ifId].ifHead = ifHead;
  }

  /* store overall number of coupling items */
  theIF[ifId].nItems = n;


  /* establish obj-table as an addressing shortcut */
  STAT_RESET1;
  IFCreateObjShortcut(ifId);
  STAT_TIMER1(T_CREATE_SHORTCUT);


  STAT_RESET1;
  update_channels(ifId);
  STAT_TIMER1(T_CREATE_COMM);

  /* TODO das handling der VCs muss noch erheblich verbessert werden */
  /* TODO durch das is_elem suchen ist alles noch VERY inefficient */

  STAT_SET_MODULE(STAT_MOD);
}


/****************************************************************************/
/*                                                                          */
/*  DDD_IFDefine                                                            */
/*                                                                          */
/****************************************************************************/

/**
        Definition of a DDD Interface.

        This function defines a new \ddd{interface}. Its argument list contains
        three arrays: the first one specifies a subset of the global DDD object set,
        the second and third ones specify a subset of all DDD priorities.
        After the initial creation of the new interface its ID is returned.

        During all following DDD operations ({\bf Identify}- as well as
        {\bf Transfer}-operations) the interface will be kept consistent, and can
        be used for communication via \funk{IFExchange} and
        \funk{IFOneway} and analogous functions.

   @param nO  number of entries in array {\it O[\ ]} (as specified with second argument).
   @param O   array of DDD type IDs; the new interface will only contain objects with one of these types.
   @param nA  number of entries in array {\it A[\ ]} (as specified with forth argument).
   @param A   first array of DDD priorities; the object set A from the new interface will only contain objects with one of these priorities.
   @param nB  number of entries in array {\it B[\ ]} (as specified with sixth argument).
   @param B   second array of DDD priorities; the object set B from the new interface will only contain objects with one of these priorities.
 */

#ifdef C_FRONTEND
DDD_IF DDD_IFDefine (
  int nO, DDD_TYPE O[],
  int nA, DDD_PRIO A[],
  int nB, DDD_PRIO B[])
{
#endif
#ifdef CPP_FRONTEND
DDD_Interface::DDD_Interface (DDD_TYPE t, DDD_PRIO p1, DDD_PRIO p2, char* name)
{
  Init(1, &t, 1, &p1, 1, &p2, name);
}

DDD_Interface::DDD_Interface (
  int nO, DDD_TYPE O[],
  int nA, DDD_PRIO A[],
  int nB, DDD_PRIO B[],
  char *name)
{
  Init(nO, O, nA, A, nB, B, name);
}

void DDD_Interface::Init (
  int nO, DDD_TYPE O[],
  int nA, DDD_PRIO A[],
  int nB, DDD_PRIO B[],
  char *name)
{
  _id = nIFs;
#endif
#ifdef F_FRONTEND
DDD_IF orgDDD_IFDefine(int, DDD_TYPE O[],int, DDD_PRIO A[],int, DDD_PRIO B[]);

void DDD_IFDefine (
  int *nO, DDD_TYPE O[],
  int *nA, DDD_PRIO A[],
  int *nB, DDD_PRIO B[],
  DDD_IF *ddd_if)
{
  *ddd_if = orgDDD_IFDefine(*nO,O,*nA,A,*nB,B);
}

DDD_IF orgDDD_IFDefine (
  int nO, DDD_TYPE O[],
  int nA, DDD_PRIO A[],
  int nB, DDD_PRIO B[])
{
#endif
int i;
COUPLING **tmpcpl;

if (nIFs==MAX_IF) {
  DDD_PrintError('E', 4100, "no more interfaces in DDD_IFDefine");
                #ifdef CPP_FRONTEND
  HARD_EXIT;
                #else
  return(0);
                #endif
}

/* construct interface definition */
theIF[nIFs].nObjStruct = nO;
theIF[nIFs].nPrioA     = nA;
theIF[nIFs].nPrioB     = nB;
memcpy(theIF[nIFs].O, O, nO*sizeof(DDD_TYPE));
memcpy(theIF[nIFs].A, A, nA*sizeof(DDD_PRIO));
memcpy(theIF[nIFs].B, B, nB*sizeof(DDD_PRIO));
if (nO>1) qsort(theIF[nIFs].O, nO, sizeof(DDD_TYPE), sort_type);
if (nA>1) qsort(theIF[nIFs].A, nA, sizeof(DDD_PRIO), sort_prio);
if (nB>1) qsort(theIF[nIFs].B, nB, sizeof(DDD_PRIO), sort_prio);


#if defined(C_FRONTEND) || defined(F_FRONTEND)
/* reset name string */
theIF[nIFs].name[0] = 0;
#endif
#ifdef CPP_FRONTEND
SetName(name);
#endif


/* compute hash tables for fast access */
theIF[nIFs].maskO = 0;
for(i=0; i<nO; i++)
  theIF[nIFs].maskO |= (1<<(unsigned int)O[i]);


/* create initial interface state */
theIF[nIFs].ifHead = NULL;
if (nCplItems>0)
{
  /* allocate temporary cpl-list, this will be too large for
     average interfaces. */
  tmpcpl = (COUPLING **) AllocTmp(sizeof(COUPLING *)*nCplItems);
  if (tmpcpl==NULL) {
    DDD_PrintError('E', 4002, STR_NOMEM " in IFDefine");
    HARD_EXIT;
  }

  IFCreateFromScratch(tmpcpl, nIFs);

  /* free temporary array */
  FreeTmp(tmpcpl);
}
else
  IFCreateFromScratch(NULL, nIFs);

nIFs++;

        #ifndef CPP_FRONTEND
return(nIFs-1);
        #endif
}



void StdIFDefine (void)
{
  /* exception: no OBJSTRUCT or priority entries */
  theIF[STD_INTERFACE].nObjStruct = 0;
  theIF[STD_INTERFACE].nPrioA     = 0;
  theIF[STD_INTERFACE].nPrioB     = 0;

  theIF[STD_INTERFACE].maskO = 0xffff;


  /* reset name string */
  theIF[nIFs].name[0] = 0;


  /* create initial interface state */
  theIF[STD_INTERFACE].ifHead = NULL;
  IFCreateFromScratch(NULL, STD_INTERFACE);
}



#ifdef C_FRONTEND
void DDD_IFSetName (DDD_IF ifId, char *name)
{
#endif
#ifdef CPP_FRONTEND
void DDD_Interface::SetName (char *name)
{
  DDD_IF ifId = _id;
#endif
#if defined(C_FRONTEND) || defined(CPP_FRONTEND)
/* ensure maximum length */
if (strlen(name) > IF_NAMELEN-1)
{
  name[IF_NAMELEN-1] = 0;
}

/* copy name string */
strcpy(theIF[ifId].name, name);
}
#endif


/****************************************************************************/

void DDD_InfoIFImpl (DDD_IF ifId)
{
  IF_PROC    *ifh;

  sprintf(cBuffer, "|\n| DDD_IFInfoImpl for proc=%03d, IF %02d\n", me, ifId);
  DDD_PrintLine(cBuffer);

  sprintf(cBuffer, "|   cpl=%08x  nIfHeads=%03d first=%08x\n",
          theIF[ifId].cpl, theIF[ifId].nIfHeads, theIF[ifId].ifHead);
  DDD_PrintLine(cBuffer);

  for(ifh=theIF[ifId].ifHead; ifh!=NULL; ifh=ifh->next)
  {
    int i;

    sprintf(cBuffer, "|   head=%08x cpl=%08x p=%03d nItems=%05d nAttrs=%03d\n",
            ifh, ifh->cpl, ifh->proc, ifh->nItems, ifh->nAttrs);
    DDD_PrintLine(cBuffer);

    sprintf(cBuffer, "|      nAB= %05d\n", ifh->nAB);
    DDD_PrintLine(cBuffer);
    for(i=0; i<ifh->nAB; i++)
    {
      COUPLING *c = ifh->cplAB[i];
      sprintf(cBuffer, "|         gid=%08x proc=%04d prio=%02d "
              "osc=%08x/%08x\n",
              OBJ_GID(c->obj), c->proc, c->prio,
              ifh->objAB[i], OBJ_OBJ(c->obj)
              );
      DDD_PrintLine(cBuffer);
    }

    sprintf(cBuffer, "|      nBA= %05d\n", ifh->nBA);
    DDD_PrintLine(cBuffer);
    for(i=0; i<ifh->nBA; i++)
    {
      COUPLING *c = ifh->cplBA[i];
      sprintf(cBuffer, "|         gid=%08x proc=%04d prio=%02d "
              "osc=%08x/%08x\n",
              OBJ_GID(c->obj), c->proc, c->prio,
              ifh->objBA[i], OBJ_OBJ(c->obj)
              );
      DDD_PrintLine(cBuffer);
    }

    sprintf(cBuffer, "|      nABA=%05d\n", ifh->nABA);
    DDD_PrintLine(cBuffer);
    for(i=0; i<ifh->nABA; i++)
    {
      COUPLING *c = ifh->cplABA[i];
      sprintf(cBuffer, "|         gid=%08x proc=%04d prio=%02d "
              "osc=%08x/%08x\n",
              OBJ_GID(c->obj), c->proc, c->prio,
              ifh->objABA[i], OBJ_OBJ(c->obj)
              );
      DDD_PrintLine(cBuffer);
    }
  }
  DDD_PrintLine("|\n");
}



static void IFDisplay (DDD_IF i)
{
  IF_PROC    *ifh;
  IF_ATTR    *ifr;
  int j;
  char buf[50];

  sprintf(cBuffer, "| IF %02d ", i);
  if (i==STD_INTERFACE)
  {
    sprintf(buf, "including all (%08x)\n|       prio all to all\n",
            theIF[i].maskO);
    strcat(cBuffer, buf);
  }
  else
  {
    strcat(cBuffer, "including ");
    for(j=0; j<theIF[i].nObjStruct; j++)
    {
      sprintf(buf, "%s ", theTypeDefs[theIF[i].O[j]].name);
      strcat(cBuffer, buf);
    }
    sprintf(buf, "(%08x)\n|       prio ", theIF[i].maskO);
    strcat(cBuffer, buf);

    for(j=0; j<theIF[i].nPrioA; j++)
    {
      sprintf(buf, "%d ", theIF[i].A[j]);
      strcat(cBuffer, buf);
    }
    strcat(cBuffer, "to ");
    for(j=0; j<theIF[i].nPrioB; j++)
    {
      sprintf(buf, "%d ", theIF[i].B[j]);
      strcat(cBuffer, buf);
    }
    strcat(cBuffer, "\n");
  }
  DDD_PrintLine(cBuffer);

  if (theIF[i].name[0]!=0)
  {
    sprintf(cBuffer, "|       '%s'\n", theIF[i].name);
    DDD_PrintLine(cBuffer);
  }

  for(ifh=theIF[i].ifHead; ifh!=NULL; ifh=ifh->next)
  {
    if (DDD_GetOption(OPT_INFO_IF_WITH_ATTR)==OPT_OFF)
    {
      sprintf(cBuffer, "|        %3d=%3d,%3d,%3d - %02d\n",
              ifh->nItems, ifh->nAB, ifh->nBA, ifh->nABA, ifh->proc);
      DDD_PrintLine(cBuffer);
    }
    else
    {
      sprintf(cBuffer, "|        %3d=%3d,%3d,%3d - %02d - #a=%05d\n",
              ifh->nItems, ifh->nAB, ifh->nBA, ifh->nABA,
              ifh->proc, ifh->nAttrs);
      DDD_PrintLine(cBuffer);

      for (ifr=ifh->ifAttr; ifr!=NULL; ifr=ifr->next)
      {
        sprintf(cBuffer, "|      a %3d=%3d,%3d,%3d - %04d\n",
                ifr->nItems, ifr->nAB, ifr->nBA, ifr->nABA, ifr->attr);
        DDD_PrintLine(cBuffer);
      }
    }
  }
}


/****************************************************************************/
/*                                                                          */
/*  DDD_IFDisplay                                                           */
/*                                                                          */
/****************************************************************************/

/**
        Display overview of single DDD interface.
        This function displays an overview table for one DDD-interface,
        its definition parameters and the current number of constituing objects
        on the calling processor.

        For each neighbour processor corresponding
        to that interface a relation line is displayed, which contains
        the overall number of objects inside the interface, the number of
        oneway relations outwards, the number of oneway relations inwards,
        the number of exchange relations and the neighbour processor number.

   @param aIF  the \ddd{interface} ID.
 */

#ifdef C_FRONTEND
void DDD_IFDisplay (DDD_IF aIF)
{
#endif
#ifdef CPP_FRONTEND
void DDD_Interface::Display (void)
{
  DDD_IF aIF = _id;
#endif
#ifdef F_FRONTEND
void DDD_IFDisplay (DDD_IF *_aIF)
{
  DDD_IF aIF = *_aIF;
#endif

if (aIF>=nIFs)
{
  sprintf(cBuffer, "invalid IF %02d in DDD_IFDisplay", aIF);
  DDD_PrintError('W', 4050, cBuffer);
  return;
}


sprintf(cBuffer, "|\n| DDD_IF-Info for proc=%03d\n", me);
DDD_PrintLine(cBuffer);

IFDisplay(aIF);

DDD_PrintLine("|\n");
}



/****************************************************************************/
/*                                                                          */
/*  DDD_IFDisplayAll                                                        */
/*                                                                          */
/****************************************************************************/

/**
        Display overview of all DDD interfaces.
        This function displays an overview table containing all DDD-interfaces,
        their definition parameters and the current number of constituing objects
        on the calling processor.

        For each interface and each neighbour processor corresponding
        to that interface a relation line is displayed, which contains
        the overall number of objects inside the interface, the number of
        oneway relations outwards, the number of oneway relations inwards,
        the number of exchange relations and the neighbour processor number.
 */

#if defined(C_FRONTEND) || defined(F_FRONTEND)
void DDD_IFDisplayAll (void)
#endif
#ifdef CPP_FRONTEND
void DDD_Interface::DisplayAll (void)
#endif
{
  int i;

  sprintf(cBuffer, "|\n| DDD_IF-Info for proc=%03d (all)\n", me);
  DDD_PrintLine(cBuffer);

  for(i=0; i<nIFs; i++)
  {
    IFDisplay(i);
  }

  DDD_PrintLine("|\n");
}



static void IFRebuildAll (void)
{
  /*
     DDD_ConsCheck();
   */

  /* create standard interface */
  IFCreateFromScratch(NULL, STD_INTERFACE);


  if (nIFs>1)
  {
    int i;

    if (nCplItems>0)
    {
      COUPLING **tmpcpl;

      /* allocate temporary cpl-list, this will be too large for
              average interfaces. */
      tmpcpl = (COUPLING **) AllocTmp(sizeof(COUPLING *)*nCplItems);
      if (tmpcpl==NULL)
      {
        DDD_PrintError('E', 4000, STR_NOMEM " in IFAllFromScratch");
        HARD_EXIT;
      }

      /* TODO: ausnutzen, dass STD_IF obermenge von allen interfaces ist */
      for(i=1; i<nIFs; i++)
      {
        IFCreateFromScratch(tmpcpl, i);

        /*
           DDD_InfoIFImpl(i);
         */
      }

      /* free temporary array */
      FreeTmp(tmpcpl);
    }
    else
    {
      /* delete old interface structures */
      for(i=1; i<nIFs; i++)
      {
        /* delete possible old interface */
        IFDeleteAll(i);
      }
    }
  }
}


void IFAllFromScratch (void)
{
  if (DDD_GetOption(OPT_IF_CREATE_EXPLICIT)==OPT_ON)
  {
    /* interfaces must be created explicitly by calling
       DDD_IFRefreshAll(). This is for doing timings from
       application level. */
    return;
  }

  IFRebuildAll();
}



void DDD_IFRefreshAll (void)
{
  if (DDD_GetOption(OPT_IF_CREATE_EXPLICIT)==OPT_OFF)
  {
    /* if interfaces are not created explicit, then they
       are always kept consistent automatically. therefore,
       this function is senseless. */
    /* return;  nevertheless, dont return, create interfaces
                once more. just to be sure. */
  }

  IFRebuildAll();
}


/****************************************************************************/

void ddd_IFInit (void)
{
  /* init lists of unused items */
  memlistIFHead = NULL;
  memlistIFAttr = NULL;

  theIF[0].ifHead = NULL;
  theIF[0].cpl    = NULL;

  /* init standard interface */
  StdIFDefine();

  /* no other interfaces yet */
  nIFs = 1;
}


void ddd_IFExit (void)
{
  int i;

  for(i=0; i<nIFs; i++)
    IFDeleteAll(i);
}


/****************************************************************************/


static size_t IFInfoMemory (DDD_IF ifId)
{
  IF_PROC *ifp;
  size_t sum=0;

  sum += sizeof(IF_PROC)    * theIF[ifId].nIfHeads;         /* component ifHead */
  sum += sizeof(COUPLING *) * theIF[ifId].nItems;           /* component cpl    */
  sum += sizeof(IFObjPtr)   * theIF[ifId].nItems;           /* component obj    */

  for(ifp=theIF[ifId].ifHead; ifp!=NULL; ifp=ifp->next)
  {
    sum += sizeof(IF_ATTR) * ifp->nAttrs;              /* component ifAttr */
  }

  return(sum);
}



#ifdef C_FRONTEND
size_t DDD_IFInfoMemory (DDD_IF ifId)
{
#endif
#ifdef CPP_FRONTEND
size_t DDD_Interface::InfoMemory (void)
{
  DDD_IF ifId = _id;
#endif
#ifdef F_FRONTEND
size_t DDD_IFInfoMemory (DDD_IF *_ifId)
{
  DDD_IF ifId = *_ifId;
#endif


if (ifId>=nIFs)
{
  sprintf(cBuffer, "invalid IF %02d in DDD_IFInfoMemory", ifId);
  DDD_PrintError('W', 4051, cBuffer);
  HARD_EXIT;
}

return(IFInfoMemory(ifId));
}



#ifdef C_FRONTEND
size_t DDD_IFInfoMemoryAll (void)
#endif
#ifdef CPP_FRONTEND
size_t DDD_Interface::InfoMemoryAll (void)
#endif
#ifdef F_FRONTEND
size_t DDD_IFInfoMemoryAll (void)
#endif
{
  int i;
  size_t sum = 0;


  for(i=0; i<nIFs; i++)
  {
    sum += IFInfoMemory(i);
  }

  return(sum);
}


/****************************************************************************/
