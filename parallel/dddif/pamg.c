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

#ifdef __OVERLAP2__
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

#define NBARRAYSIZE     40      /* only temp. for pamg check */
static INT MaximumMatrices;     /* only temp. for pamg check */
static INT pamgerrors;  /* only temp. for pamg check */

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/

/****************************************************************************/

INT pamgDo( MULTIGRID *theMG, INT level )
{
  GRID *grid = GRID_ON_LEVEL(theMG,level);
  VECTOR *v;

  DDD_XferBegin();
  for( v=PFIRSTVECTOR(grid); PRIO(v) != PrioBorder && PRIO(v) != PrioMaster; v = SUCCVC(v))
  {
    DDD_XferPrioChange( PARHDR((NODE*)VOBJECT(v)), PrioBorder );                     /* TODO: cancel this line; its only for beauty in checks */
    DDD_XferPrioChange( PARHDR(v), PrioBorder );
  }
  /* elements with old ghostprios cause errors in check; but that's ok */
  DDD_XferEnd();

  return 0;
}

static int CountMatrices (DDD_OBJ obj)
{
  VECTOR *v = (VECTOR *)obj;
  MATRIX *m, *m2;
  int i=0;

  for( m=MNEXT(VSTART(v)); m!=NULL; m = MNEXT(m) )
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
  for( m=MNEXT(VSTART(v)); m!=NULL; m = MNEXT(m) )
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

static int Scatter_pamgCheck (DDD_OBJ obj, void *data)
{
  VECTOR  *v = (VECTOR *)obj;
  MATRIX *m, *m2;
  DDD_GID *buf = (DDD_GID *)data;
  DDD_GID gid,loc_gids[NBARRAYSIZE*NBARRAYSIZE];
  INT i, nr_local_gids, sender_gid, sender_pe;

  nr_local_gids=0;
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
      UserWriteF("\nERROR GID %d on PE %d had NB with GID %d", sender_gid, sender_pe, buf[i]);
    }
}

INT pamgCheckDo( MULTIGRID *theMG, INT level )
{
  GRID *grid;
  size_t sizePerVector;

  for( ; level >= BOTTOMLEVEL(theMG); level -- )
  {
    grid = GRID_ON_LEVEL(theMG,level);
    UserWriteF( "pamgCheckDo: checking level %d:", level );

    MaximumMatrices=0;
    DDD_IFAExecLocal(BorderVectorIF, GRID_ATTR(grid), CountMatrices);
    MaximumMatrices = UG_GlobalMaxINT(MaximumMatrices);
    sizePerVector = (MaximumMatrices+3) * sizeof(DDD_GID);              /* 3 for additional infos */

    pamgerrors = 0;
    DDD_IFAOneway(BorderVectorIF, GRID_ATTR(grid),IF_FORWARD, sizePerVector,
                  Gather_pamgCheck, Scatter_pamgCheck);
    pamgerrors= UG_GlobalSumINT(pamgerrors);
    if( pamgerrors == 0 )
      UserWrite( " OK\n" );
    else
      UserWriteF( "%d ERRORS\n", pamgerrors );

    /*printmgrid( grid, 1 );*/
  }

  return 0;
}
#endif /* __OVERLAP2__ */

#endif  /* ModelP */
