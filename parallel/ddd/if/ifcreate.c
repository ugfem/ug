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

#include "dddi.h"
#include "if.h"


#define DebugShowAttr


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



static IF_ATTR *NewIFAttr ()
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




static int sort_int (const void *e1, const void *e2)
{
  if (*(int *)e1 < *(int *)e2) return(-1);
  if (*(int *)e1 == *(int *)e2) return(0);
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
/*                3. attr property of objects                              */
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
}



/* TODO  el-set relation, VERY inefficient! */
static int is_elem (unsigned int el, int n, unsigned int *set)
{
  int i;

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

  MarkHeap();

  for(i=0, ifh=theIF[ifId].ifHead; ifh!=NULL; i++, ifh=ifh->next)
  {
    partners[i] = ifh->proc;
  }

  DDD_GetChannels(theIF[ifId].nIfHeads);

  for(ifh=theIF[ifId].ifHead; ifh!=NULL; ifh=ifh->next) {
    ifh->vc = VCHAN_TO(ifh->proc);
  }

  ReleaseHeap();
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
    DDD_PrintError('E', 4000, "not enough memory in IFCreateFromScratch");
    exit(1);
  }

  /* collect couplings */
  n=0;
  for(index=0; index<nCpls; index++)
  {
    COUPLING  *cpl;

    for(cpl=theCpl[index]; cpl!=NULL; cpl=cpl->next)
    {
      cplarray[n] = cpl;
      n++;
    }
  }
  /*
     printf("%04d: n=%d, nCplItems=%d\n",me,n,nCplItems);
   */

  if (n!=nCplItems) {
    DDD_PrintError('F', 9999, "internal nCplItems-mismatch");
    exit(1);
  }

  return(cplarray);
}


/****************************************************************************/

void IFCreateFromScratch (DDD_IF ifId)
{
  IF_PROC     *ifHead, *lastIfHead;
  IF_ATTR    *ifAttr, *lastIfAttr;
  int n, i;
  DDD_PROC lastproc;


  /* first delete possible old interface */
  IFDeleteAll(ifId);

  if (ifId==STD_INTERFACE)
  {
    theIF[ifId].cpl = IFCollectStdCouplings();
    n = nCplItems;
  }
  else
  {
    int index;

    /* get memory for couplings inside IF */
    if (nCplItems>0)
    {
      theIF[ifId].cpl = (COUPLING **) AllocIF(sizeof(COUPLING *)*nCplItems);
      /* TODO: nCplItems will be too big for average interfaces! */
      if (theIF[ifId].cpl==NULL) {
        DDD_PrintError('E', 4000, "not enough memory in IFCreateFromScratch");
        exit(1);
      }
    }
    else
    {
      theIF[ifId].cpl = NULL;
    }


    /* collect relevant couplings */
    n=0;
    for(index=0; index<nCpls; index++)
    {
      /* determine whether object belongs to IF */
      if ((1<<OBJ_TYPE(theObj[index])) & theIF[ifId].maskO)
      {
        int objInA, objInB;

        objInA = is_elem(OBJ_PRIO(theObj[index]),
                         theIF[ifId].nPrioA, theIF[ifId].A);
        objInB = is_elem(OBJ_PRIO(theObj[index]),
                         theIF[ifId].nPrioB, theIF[ifId].B);

        if (objInA || objInB)
        {
          COUPLING  *cpl;

          /* test coupling list */
          for(cpl=theCpl[index]; cpl!=NULL; cpl=cpl->next)
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
              theIF[ifId].cpl[n] = cpl;
              n++;
            }
          }
        }
      }
    }
  }


  /* sort IF couplings */
  if (n>1)
    qsort(theIF[ifId].cpl, n, sizeof(COUPLING *), sort_IFCouplings);


  /* create IF_PROCs */
  lastproc = -1;
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
      ifHead->nItems = 0;
      ifHead->cpl    = cplp;
      ifHead->nAB    = ifHead->nBA   = ifHead->nABA   = 0;
      ifHead->cplAB  = ifHead->cplBA = ifHead->cplABA = NULL;
      ifHead->proc   = cpl->proc;
      ifHead->next   = lastIfHead;
      lastIfHead = ifHead;
      lastproc   = ifHead->proc;

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

  /* remember anchor of ifHead list */
  if (theIF[ifId].nIfHeads>0) {
    theIF[ifId].ifHead = ifHead;
  }

  /* store overall number of coupling items */
  theIF[ifId].nItems = n;


  /* TODO: an dieser stelle koennte das alte (zu grosse!)
          cpl-array gegen ein kleineres der groesse theIF[ifId].nItems
          ausgetauscht werden ... */


  /* establish obj-table as an addressing shortcut */
  IFCreateObjShortcut(ifId);


  update_channels(ifId);

  /* TODO das handling der VCs muss noch erheblich verbessert werden */
  /* TODO durch das is_elem suchen ist alles noch VERY inefficient */
}


#ifdef C_FRONTEND
DDD_IF DDD_IFDefine (
  int nO, DDD_TYPE O[],
  int nA, DDD_PRIO A[],
  int nB, DDD_PRIO B[])
#else
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
#endif
{
  int i;

  if (nIFs==MAX_IF) {
    DDD_PrintError('E', 4100, "no more interfaces in DDD_IFDefine");
    return(0);
  }

  /* construct interface definition */
  theIF[nIFs].nObjStruct = nO;
  theIF[nIFs].nPrioA     = nA;
  theIF[nIFs].nPrioB     = nB;
  memcpy(theIF[nIFs].O, O, nO*sizeof(DDD_TYPE));
  memcpy(theIF[nIFs].A, A, nA*sizeof(DDD_PRIO));
  memcpy(theIF[nIFs].B, B, nB*sizeof(DDD_PRIO));
  if (nO>1) qsort(theIF[nIFs].O, nO, sizeof(DDD_TYPE), sort_int);
  if (nA>1) qsort(theIF[nIFs].A, nA, sizeof(DDD_PRIO), sort_int);
  if (nB>1) qsort(theIF[nIFs].B, nB, sizeof(DDD_PRIO), sort_int);


  /* compute hash tables for fast access */
  theIF[nIFs].maskO = 0;
  for(i=0; i<nO; i++)
    theIF[nIFs].maskO |= (1<<(unsigned int)O[i]);


  /* create initial interface state */
  theIF[nIFs].ifHead = NULL;
  IFCreateFromScratch(nIFs);

  nIFs++;

  return(nIFs-1);
}


void StdIFDefine()
{
  /* exception: no OBJSTRUCT or priority entries */
  theIF[STD_INTERFACE].nObjStruct = 0;
  theIF[STD_INTERFACE].nPrioA     = 0;
  theIF[STD_INTERFACE].nPrioB     = 0;

  theIF[STD_INTERFACE].maskO = 0xffff;

  /* create initial interface state */
  theIF[STD_INTERFACE].ifHead = NULL;
  IFCreateFromScratch(STD_INTERFACE);
}


void DDD_InfoIFImpl (DDD_IF ifId)
{
  IF_PROC    *ifh;
  IF_ATTR   *ifr;

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



void DDD_IFDisplay (void)
{
  IF_PROC    *ifh;
  IF_ATTR    *ifr;
  int i, j;
  char buf[50];

  sprintf(cBuffer, "|\n| DDD_IFInfo for proc=%03d\n", me);
  DDD_PrintLine(cBuffer);

  for(i=0; i<nIFs; i++)
  {
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

    for(ifh=theIF[i].ifHead; ifh!=NULL; ifh=ifh->next)
    {
#               ifndef DebugShowAttr
      sprintf(cBuffer, "|        %3d=%3d,%3d,%3d - %02d\n",
              ifh->nItems, ifh->nAB, ifh->nBA, ifh->nABA, ifh->proc);
      DDD_PrintLine(cBuffer);

#               else
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
#               endif
    }

  }
  DDD_PrintLine("|\n");
}




void IFAllFromScratch (void)
{
  int i;
  /*
     DDD_ConsCheck();
   */

  /* TODO effizienter: ausnutzen, dass STD_IF obermenge von allen interfaces ist */
  for(i=0; i<nIFs; i++)
  {
    IFCreateFromScratch(i);
    /*
       DDD_InfoIFImpl(i);
     */
  }
}


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
