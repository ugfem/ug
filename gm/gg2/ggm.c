// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      ggm.c                                                         */
/*                                                                          */
/* Purpose:   manager for grid generator                                                        */
/*                                                                          */
/* Author:    Wolfgang Hoffmann, Henrik Renz-Reichert	                    */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart, Germany										*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   08.03.94 begin, ug version 2.2                                */
/*                15.10.95 implemented in ug31                                  */
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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "compiler.h"
#include "gm.h"
#include "heaps.h"
#include "misc.h"
#include "ugm.h"
#include "debug.h"
#include "general.h"

#include "ggm.h"

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
/* data structures used in this source file (exported data structures are   */
/*        in the corresponding include file!)                               */
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

static MULTIGRID *MG;
static INT ggMGUDid;
static MG_GGDATA *myMGdata;

static INT IflObj;
static INT FlObj;
static INT FcObj;

/* RCS string */
RCSID("$Header$",UG_RCS_STRING)

/****************************************************************************/
/*                                                                          */
/* forward declarations of functions used before they are defined           */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/* Function:  SetFlagsfortemporaryGGObjects                                     */
/*                                                                          */
/* Purpose :                                                                                                                        */
/*                                                                          */
/****************************************************************************/

INT SetFlagsfortemporaryGGObjects(INT IflObject,INT FlObject,INT FcObject)
{
  IflObj = IflObject;
  FlObj  = FlObject;
  FcObj  = FcObject;

  return(0);

}

/****************************************************************************/
/*                                                                          */
/* Function:  CreateIndepFrontList                                              */
/*                                                                          */
/* Purpose:   return pointer to a new independent front list structure      */
/*                                                                          */
/* Input:     GRID *theGrid: grid where vertex should be inserted               */
/*                                                                          */
/* Output:    INDEPFRONTLIST * : pointer to new independent front list      */
/*            NULL: if an error occured                                     */
/*                                                                          */
/****************************************************************************/

INDEPFRONTLIST *CreateIndepFrontList (GRID *theGrid)
{
  INDEPFRONTLIST *ipfl;

  ipfl = GetMemoryForObject(theGrid->mg,sizeof(INDEPFRONTLIST),IflObj);
  if (ipfl==NULL) return(NULL);

  /* initialize data */
  CTRL(ipfl) = 0;
  SETOBJT(ipfl,IflObj);
  STARTFL(ipfl) = NULL;
  LASTFL(ipfl)  = NULL;
  NFL(ipfl) = 0;
  MYGRID(ipfl) = theGrid;

  /* insert in independent front list list */

  SUCCIFL(ipfl) = STARTIFL(myMGdata);
  if (SUCCIFL(ipfl)!=NULL) PREDIFL(SUCCIFL(ipfl)) = ipfl;
  PREDIFL(ipfl) = NULL;
  STARTIFL(myMGdata) = ipfl;
  if (LASTIFL(myMGdata)==NULL)
    LASTIFL(myMGdata) = ipfl;
  /* counters */
  NIFL(myMGdata)++;

  return(ipfl);

}

/****************************************************************************/
/*                                                                          */
/* Function:  CreateFrontList                                                           */
/*                                                                          */
/* Purpose:   return pointer to a new front list structure                              */
/*                                                                          */
/* Input:     GRID *theGrid: grid where vertex should be inserted               */
/*                                                                          */
/* Output:    FRONTLIST *: pointer to new front list                        */
/*            NULL: if an error occured                                     */
/*                                                                          */
/****************************************************************************/

FRONTLIST *CreateFrontList (INDEPFRONTLIST *myIFL)
{
  FRONTLIST *pfl;

  pfl = GetMemoryForObject(MYGRID(myIFL)->mg,sizeof(FRONTLIST),FlObj);
  if (pfl==NULL) return(NULL);

  /* initialize data */
  CTRL(pfl) = 0;
  SETOBJT(pfl,FlObj);
  STARTFC(pfl) = NULL;
  LASTFC(pfl)  = NULL;
  NFC(pfl) = 0;
  MYGRID(pfl) = MYGRID(myIFL);
  MYIFL(pfl)   = myIFL;

  /* insert in myIFL at begin of myIFL */
  SUCCFL(pfl) = STARTFL(myIFL);
  if (SUCCFL(pfl)!=NULL) PREDFL(SUCCFL(pfl)) = pfl;
  PREDFL(pfl) = NULL;
  STARTFL(myIFL) = pfl;
  if (LASTFL(myIFL)==NULL)
    LASTFL(myIFL) = pfl;

  /* increment counter */
  NFL(myIFL)++;

  return(pfl);
}

/****************************************************************************/
/*                                                                          */
/* Function:  CreateFrontComp                                                           */
/*                                                                          */
/* Purpose:   return pointer to a new front component structure                         */
/*                                                                          */
/* Input:     GRID *theGrid: grid where vertex should be inserted               */
/*                                                                          */
/* Output:    FRONTCOMP *: pointer to new front component                   */
/*            NULL: if an error occured                                     */
/*                                                                          */
/****************************************************************************/

FRONTCOMP *CreateFrontComp (FRONTLIST *mylist, FRONTCOMP *after, int ncomp, NODE **NodeHandle)
{
  MULTIGRID *theMG;
  FRONTCOMP *pfc,*FChandle;
  int i;

  theMG = MYMG(MYGRID(mylist));

  if (mylist==NULL)
    return (NULL);

  if  (ncomp<1)
    return (NULL);

  if  ( ncomp == 1 )
  {
    pfc = GetMemoryForObject(theMG,sizeof(FRONTCOMP),FcObj);
    if (pfc==NULL) return(NULL);

    /* initialize data */
    SETOBJT(pfc,FcObj);
    FRONTN(pfc) = *NodeHandle;
    MYFL(pfc) = mylist;

    /* insert in (cyclic!) front comp list */
    if (after==NULL)
    {
      if (STARTFC(mylist)==NULL)
      {
        STARTFC(mylist) = pfc;
        LASTFC(mylist)  = pfc;
        SUCCFC(pfc) = pfc;
        PREDFC(pfc) = pfc;
      }
      else
      {
        SUCCFC(pfc) = STARTFC(mylist);
        PREDFC(pfc) = LASTFC(mylist);
        PREDFC(SUCCFC(pfc)) = pfc;
        SUCCFC(PREDFC(pfc)) = pfc;
        STARTFC(mylist) = pfc;
      }
    }
    else
    {
      SUCCFC(pfc) = SUCCFC(after);
      PREDFC(pfc) = after;
      PREDFC(SUCCFC(pfc)) = pfc;
      SUCCFC(after) = pfc;
      if (after==LASTFC(mylist))
        LASTFC(mylist) = pfc;
    }

    /* counters */
    NFC(mylist)++;

    return(pfc);
  }

  /* get storage for ncomp FCs */
  FChandle = (FRONTCOMP *) GetMem(theMG->theHeap,ncomp*sizeof(FRONTCOMP),FROM_BOTTOM);

  if (FChandle==NULL)
    return (NULL);

  /* init controlword */
  for (i=0; i<ncomp; i++)
  {
    FChandle[i].control = 0;
    SETOBJT(&(FChandle[i]),FcObj);
    FRONTN(&(FChandle[i])) = NodeHandle[i];
    MYFL(&(FChandle[i])) = mylist;
    FCNGB(&(FChandle[i])) = NULL;
    FCNGBS(&(FChandle[i])) = NULL;
  }

  /* create pointer connections */
  for (i=1; i<ncomp; i++)
  {
    SUCCFC(&(FChandle[i-1])) = &(FChandle[i]);
    PREDFC(&(FChandle[i]))   = &(FChandle[i-1]);
  }

  if (STARTFC(mylist)==NULL)
  {
    /* mylist is empty: make FC list cyclic */
    SUCCFC(&(FChandle[ncomp-1])) = &(FChandle[0]);
    PREDFC(&(FChandle[0]))        = &(FChandle[ncomp-1]);
    STARTFC(mylist) = &(FChandle[0]);
    LASTFC(mylist) = &(FChandle[ncomp-1]);
    NFC(mylist) = ncomp;
  }
  else
  {
    if (after==NULL)
    {
      SUCCFC(&(FChandle[ncomp-1])) = STARTFC(mylist);
      PREDFC(&(FChandle[0])) = LASTFC(mylist);
      PREDFC(SUCCFC(&(FChandle[ncomp-1]))) = &(FChandle[ncomp-1]);
      SUCCFC(PREDFC(&(FChandle[0]))) = &(FChandle[0]);
      PREDFC(STARTFC(mylist)) = &(FChandle[ncomp-1]);
      STARTFC(mylist) = &(FChandle[0]);
    }
    else
    {
      SUCCFC(&(FChandle[ncomp-1])) = SUCCFC(after);
      PREDFC(&(FChandle[0])) = after;
      PREDFC(SUCCFC(&(FChandle[ncomp-1]))) = &(FChandle[ncomp-1]);
      SUCCFC(after) = &(FChandle[0]);
      if (after==LASTFC(mylist))
        LASTFC(mylist) = &(FChandle[ncomp-1]);
    }
    NFC(mylist) += ncomp;
  }

  return (&(FChandle[ncomp-1]));
}

/****************************************************************************/
/*                                                                          */
/* Function:  DisposeADVfront		                                        */
/*                                                                          */
/* Purpose:   remove whole advancing front data structure from grid			*/
/*                                                                          */
/* Input:     GRID *theGrid: grid to remove from                                                */
/*            FRONTLIST *theFL: front list to remove                        */
/*                                                                          */
/* Output:    INT 0: ok                                                                                                         */
/*                                                                          */
/****************************************************************************/

INT DisposeADVfront (GRID *theGrid)
{
  INDEPFRONTLIST *theIFL,*nextIFL;

  /* remove all ADVfront lists */
  for (theIFL=LASTIFL(myMGdata); theIFL!=NULL; theIFL=nextIFL)
  {
    nextIFL = PREDIFL(theIFL);                          /* remember pred before disposing */
    DisposeIndepFrontList(theIFL);
  }

  STARTIFL(myMGdata) = NULL;
  LASTIFL(myMGdata)  = NULL;
  NIFL(myMGdata)     = 0;

  return (0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  DisposeIndepFrontList                                             */
/*                                                                          */
/* Purpose:   remove independent front list including its front components	*/
/*                                                                          */
/* Input:     GRID *theGrid: grid to remove from                                                */
/*            FRONTLIST *theFL: front list to remove                        */
/*                                                                          */
/* Output:    INT 0: ok                                                                                                         */
/*                                                                          */
/****************************************************************************/

INT DisposeIndepFrontList (INDEPFRONTLIST *theIFL)
{
  GRID *theGrid;
  FRONTLIST *theFL;

  HEAPFAULT(theIFL);

  theGrid = MYGRID(theIFL);

  /* remove front lists from theIFL */
  for (theFL=STARTFL(theIFL); theFL!=NULL; theFL=SUCCFL(theFL))
    if (DisposeFrontList(theFL)>0) return (1);


  /* remove independent front list from list */
  if (PREDIFL(theIFL)!=NULL)
    SUCCIFL(PREDIFL(theIFL)) = SUCCIFL(theIFL);
  else
  {
    STARTIFL(myMGdata) = SUCCIFL(theIFL);
  }
  if (SUCCIFL(theIFL)!=NULL)
    PREDIFL(SUCCIFL(theIFL)) = PREDIFL(theIFL);
  if (LASTIFL(myMGdata)==theIFL)
    LASTIFL(myMGdata) = PREDIFL(theIFL);
  NIFL(myMGdata)--;

  /* delete the independent front list itself */
  PutFreeObject(theGrid->mg,theIFL,sizeof(INDEPFRONTLIST),IflObj);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  DisposeFrontComp                                                          */
/*                                                                          */
/* Purpose:   remove front component                                                                        */
/*                                                                          */
/* Input:     GRID *theGrid: grid to remove from                                                */
/*            FRONTCOMP *theFC: front component to remove                   */
/*                                                                          */
/* Output:    INT 0: ok                                                                                                         */
/*                                                                          */
/****************************************************************************/

INT DisposeFrontComp (FRONTLIST *myList, FRONTCOMP *theFC)
{
  HEAPFAULT(theFC);

  if (STARTFC(myList) == LASTFC(myList))
  {
    DisposeFrontList (myList);
    return(0);
  }

  SUCCFC(PREDFC(theFC)) = SUCCFC(theFC);
  PREDFC(SUCCFC(theFC)) = PREDFC(theFC);
  if (theFC == STARTFC(myList))
    STARTFC(myList) = SUCCFC(theFC);
  else if (theFC == LASTFC(myList))
    LASTFC(myList) = PREDFC(theFC);

  PutFreeObject(MYGRID(myList)->mg,theFC,sizeof(FRONTCOMP),FcObj);

  NFC(myList)--;

  return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  DisposeFrontList                                                          */
/*                                                                          */
/* Purpose:   remove front list including its front components		        */
/*                                                                          */
/* Input:     GRID *theGrid: grid to remove from                                                */
/*            FRONTLIST *theFL: front list to remove                        */
/*                                                                          */
/* Output:    INT 0: ok                                                                                                         */
/*                                                                          */
/****************************************************************************/

INT DisposeFrontList (FRONTLIST *theFL)
{
  MULTIGRID *theMG;
  INDEPFRONTLIST *myIFL;

  HEAPFAULT(theFL);

  myIFL = MYIFL(theFL);
  theMG = MYMG(MYGRID(theFL));

  /* remove front components of theFL */
  while (STARTFC(theFL) != LASTFC(theFL))
    DisposeFrontComp(theFL,STARTFC(theFL));

  HEAPFAULT(STARTFC(theFL));
  PutFreeObject(theMG,STARTFC(theFL),sizeof(FRONTCOMP),FcObj);

  /* remove front list from list */
  if (PREDFL(theFL)!=NULL)
    SUCCFL(PREDFL(theFL)) = SUCCFL(theFL);
  else
    STARTFL(myIFL) = SUCCFL(theFL);
  if (SUCCFL(theFL)!=NULL)
    PREDFL(SUCCFL(theFL)) = PREDFL(theFL);
  if (LASTFL(myIFL)==theFL)
    LASTFL(myIFL) = PREDFL(theFL);

  NFL(myIFL)--;

  /* delete the front list itself */
  PutFreeObject(theMG,theFL,sizeof(FRONTLIST),FlObj);

  return(0);
}

/****************************************************************************/
/*                                                                          */
/* Function:  GetMGdataPointer                                                          */
/*                                                                          */
/* Purpose:   delivers a pointer to the Multigrid Data                                      */
/*                                                                          */
/****************************************************************************/

MG_GGDATA *GetMGdataPointer (MULTIGRID *theMG)
{
  myMGdata  = (MG_GGDATA*) GEN_MGUD_ADR(theMG,OFFSET_IN_MGUD(ggMGUDid));
  return(myMGdata);
}

/****************************************************************************/
/*                                                                          */
/* initialization for grid generator library			                                */
/*                                                                          */
/****************************************************************************/

INT InitGGManager ()
{
  /* allocate storage in general mg user data */
  ggMGUDid = GetNewBlockID();
  if (DefineMGUDBlock(ggMGUDid,sizeof(MG_GGDATA))!=GM_OK) return (__LINE__);

  return(0);
}
