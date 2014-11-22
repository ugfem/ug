// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  pamg.c                                                                                                        */
/*																			*/
/* Purpose:   aux. routines for parallel amg								*/
/*																			*/
/* Author:	  Christian Wrobel		                                                                */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   980812 cw  begin                                                                          */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


#ifdef ModelP

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

#include <config.h>
#include <assert.h>

#include "parallel.h"
#include "gm.h"
#include "evm.h"
#include "general.h"
#include "ugm.h"
#include "refine.h"
#include "algebra.h"
#include "debug.h"
#include "ugdevices.h"
#include "namespace.h"

/* UG namespaces: */
USING_UG_NAMESPACES

#ifdef USE_FAMG
/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* definition of exported global variables									*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);

START_UGDIM_NAMESPACE

#define NBARRAYSIZE     40      /* only temp. for pamg check */
static INT MaximumMatrices;     /* only temp. for pamg check */
static INT pamgerrors;  /* only temp. for pamg check */

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/

/* TODO:remove*/
void NS_DIM_PREFIX printall(GRID *grid, char *text)
{
  VECTOR *v;
  for( v=PFIRSTVECTOR(grid); v!=NULL; v = SUCCVC(v))
  {
    int *pl = DDD_InfoProcList(PARHDR(v));
    int n;

    if( pl==NULL)
    {
      printf("%d pamgDo %s: %x "VINDEX_FMTX " keine pl\n",me,text,v,VINDEX_PRTX(v));
    }
    else
    {
      for(n=0; pl[2*n]!=-1; n++)
        ;

      switch(n)
      {
      case 0 :
        printf("%d pamgDo %s: %x "VINDEX_FMTX " KEINE pl???\n",me,text,v,VINDEX_PRTX(v));
        break;
      case 1 :
        printf("%d pamgDo %s: %x "VINDEX_FMTX " %d %d\n",me,text,v,VINDEX_PRTX(v),pl[0],pl[1]);
        break;
      case 2 :
        printf("%d pamgDo %s: %x "VINDEX_FMTX " %d %d, %d %d\n",me,text,v,VINDEX_PRTX(v),pl[0],pl[1],pl[2],pl[3]);
        break;
      case 3 :
        printf("%d pamgDo %s: %x "VINDEX_FMTX " %d %d, %d %d, %d %d\n",me,text,v,VINDEX_PRTX(v),pl[0],pl[1],pl[2],pl[3],pl[4],pl[5]);
        break;
      default :
        printf("%d pamgDo %s: %x "VINDEX_FMTX " hat %d pl Eintraege????\n",me,text,v,VINDEX_PRTX(v),n);
      }
    }
  }
}

static int CountMatrices (DDD_OBJ obj)
{
  VECTOR *v = (VECTOR *)obj;
  MATRIX *m, *m2;
  int i=0;

  m=VSTART(v);
  if(m==NULL)
    return 0;                   /* here is no matrix */

  for( m=MNEXT(m); m!=NULL; m = MNEXT(m) )
    for( m2=VSTART(MDEST(m)); m2!=NULL; m2 = MNEXT(m2),i++ )
      ;

  MaximumMatrices = MAX(MaximumMatrices,i);

  return (0);
}

static int Gather_pamgCheck (DDD_OBJ obj, void *data)
{
  VECTOR  *v = (VECTOR *)obj;
  MATRIX *m, *m2;
  DDD_GID *buf = (DDD_GID *)data;
  INT i;

  i=1;          /* keep first free for length */
  buf[i++] = (DDD_GID)DDD_InfoGlobalId(PARHDR(v));
  buf[i++] = (DDD_GID)me;
  m=VSTART(v);
  if(m!=NULL)
    for( m=MNEXT(m); m!=NULL; m = MNEXT(m) )
      for( m2=VSTART(MDEST(m)); m2!=NULL; m2 = MNEXT(m2) )
        buf[i++] = DDD_InfoGlobalId(PARHDR(MDEST(m2)));
  buf[0] = (DDD_GID)i;
  ASSERT(i <= MaximumMatrices+3);

  return (0);
}

static INT isAlocalGID( DDD_GID search_gid, DDD_GID loc_gids[NBARRAYSIZE*NBARRAYSIZE], INT nr_loc_gids )
/* 1 if found, 0 else */
/* could be improved by a bisection search; then loc_gids must be sorted */
{
  for( ; --nr_loc_gids>=0; )
    if( search_gid == *loc_gids++ )
      return 1;
  return 0;
}


/* needed for qsort static int sort_Gids (const void *e1, const void *e2)
   {
    if ((*(DDD_GID *)e1) < (*(DDD_GID *)e2)) return(-1);
    if ((*(DDD_GID *)e1) > (*(DDD_GID *)e2)) return(1);

    return(0);
   }*/


/* TODO: remove unused variants */
static int Scatter_pamgCheckQQQQQQQQQQQQQQQQQQQQQQQQ (DDD_OBJ obj, void *data)
{
  VECTOR  *v = (VECTOR *)obj;
  MATRIX *m, *m2;
  DDD_GID *buf = (DDD_GID *)data;
  DDD_GID gid,loc_gids[NBARRAYSIZE*NBARRAYSIZE];
  INT i, nr_local_gids, sender_gid, sender_pe;

  /* the assertions must be valid only for masters */
  if( PRIO(v)!=PrioMaster )
    return 0;

  nr_local_gids=0;
  m=VSTART(v);
  if(m!=NULL)
    for( m=MNEXT(VSTART(v)); m!=NULL; m = MNEXT(m) )
      for( m2=VSTART(MDEST(m)); m2!=NULL; m2 = MNEXT(m2) )
        loc_gids[nr_local_gids++] = DDD_InfoGlobalId(PARHDR(MDEST(m2)));

  /*qsort( loc_gid, nr_local_gids, sizeof(DDD_GID), sort_Gids); makes only sense for bisection search */

  sender_gid = (INT)buf[1];
  sender_pe = (INT)buf[2];
  for( i=3; i<(INT)buf[0]; i++ )
    if( !isAlocalGID(buf[i],loc_gids,nr_local_gids) )
    {
      pamgerrors++;
      UserWriteF("\nERROR GID %08x on PE %d had NB with GID %08x", sender_gid, sender_pe, buf[i]);
    }

  return 0;
}

static int Scatter_pamgCheck (DDD_OBJ obj, void *data)
{
  VECTOR  *v = (VECTOR *)obj;
  MATRIX *m, *m2;
  DDD_GID *buf = (DDD_GID *)data;
  DDD_GID gid,loc_gids[NBARRAYSIZE*NBARRAYSIZE];
  INT i, nr_local_gids, sender_gid, sender_pe;
  /* INT in_overlap2;	not necessary any longer because scatter is called only for matsers if using BorderVectorIF and IF_FORWARD*/

  assert(PRIO(v)==PrioMaster);
  nr_local_gids=0;
  m=VSTART(v);
  /*in_overlap2 = (PRIO(v)==PrioMaster);*/
  if(m!=NULL)
    for( m=MNEXT(VSTART(v)); m!=NULL; m = MNEXT(m) )
      for( m2=VSTART(MDEST(m)); m2!=NULL; m2 = MNEXT(m2) )
      {
        loc_gids[nr_local_gids++] = DDD_InfoGlobalId(PARHDR(MDEST(m2)));
        /*in_overlap2 |= (PRIO(MDEST(m2))==PrioMaster);*/
      }

  /*if(!in_overlap2)*/
  /*	return 0;	 only vectors within overlap 2 must have the complete neighbourhood */

  /*qsort( loc_gid, nr_local_gids, sizeof(DDD_GID), sort_Gids); makes only sense for bisection search */

  sender_gid = (INT)buf[1];
  sender_pe = (INT)buf[2];
  for( i=3; i<(INT)buf[0]; i++ )
    if( !isAlocalGID(buf[i],loc_gids,nr_local_gids) )
    {
      pamgerrors++;
      UserWriteF("\nERROR GID %08x on PE %d had NB with GID %08x", sender_gid, sender_pe, buf[i]);
    }

  return 0;
}

INT NS_DIM_PREFIX pamgCheckDo( MULTIGRID *theMG, INT level )
{
  GRID *grid;
  size_t sizePerVector;

  for( ; level > BOTTOMLEVEL(theMG); level -- )
  {
    grid = GRID_ON_LEVEL(theMG,level);
    UserWriteF( "pamgCheckDo: checking level %d:", level );

    MaximumMatrices=0;
    DDD_IFAExecLocal(BorderVectorIF, GRID_ATTR(grid), CountMatrices);
    MaximumMatrices = UG_GlobalMaxINT(MaximumMatrices);
    sizePerVector = (MaximumMatrices+3) * sizeof(DDD_GID);              /* 3 for additional infos */

    pamgerrors = 0;
    /* Border -> Master */
    DDD_IFAOneway(BorderVectorIF, GRID_ATTR(grid),IF_FORWARD, sizePerVector,
                  Gather_pamgCheck, Scatter_pamgCheck);
    pamgerrors= UG_GlobalSumINT(pamgerrors);
    if( pamgerrors == 0 )
      UserWrite( " OK\n" );
    else
      UserWriteF( " %d ERRORS\n", pamgerrors );

    /*printmgrid( grid, 1 );*/
  }

  return 0;
}

END_UGDIM_NAMESPACE

#endif /* USE_FAMG */

#endif  /* ModelP */
