// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  coviseif.c                                                                                                    */
/*																			*/
/* Purpose:   interface ug <-> covise                                                   */
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken										*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de								*/
/*																			*/
/* History:   10.12.97 begin                                                                            */
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/

/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

/* ug includes */
#include "gm.h"
#include "evm.h"
#include "general.h"

/* covise includes */
/*#include "covise.h"*/



/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/

#define MAX_STRLEN    64
#define MAX_SOLUTIONS 16


/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct
{
  INT n_components;
  char name[MAX_STR];
} SOLUTION_DESC;


typedef struct
{
  /* multigrid description */
  INT min_level;        /* coarsest grid level, may be <0 for algebraic methods */
  INT max_level;        /* ==TOPLEVEL in ug */

  /* surface grid description */
  INT n_nodes;       /* number of surfaces grid nodes (UG) or vertices (Covise) */
  INT n_elems;       /* number of surface grid elements */
  INT n_conns;       /* sum of corners over all surface elements */

  /* numerical data */
  INT num_solutions;       /* number of solutions (may be scalar or vector) */
  SOLUTION_DESC solutions[MAX_SOLUTIONS];       /* solution descriptors */

} COVISE_HEADER;


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

/* ug: */
/* not used currently
   static EVALUES *ElemEval;
   static EVECTOR *ElemVec;
   static PreprocessingProcPtr EvalPreProc;
   static ElementEvalProcPtr EvalPlotProc;
   static ElementVectorProcPtr EvecPlotProc;
   static int FirstEvec;
 */


/* covise: */


/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* functions
   /*																			*/
/****************************************************************************/


static INT ComputeSurfaceGridStats (COVISE_HEADER *covise)
{
  INT l, n_node, n_elem, n_conn;

  /* surface grid up to current level */
  n_node = n_elem = n_conn = 0;
  for (l=covise->min_level; l<=covise->max_level; l++)
  {
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);

    /* reset USED flags in all objects to be counted */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      SETUSED(theNode,0);
    }

    /* count geometric objects */
    for (theElement=FIRSTELEMENT(theGrid);
         theElement!=NULL; theElement=SUCCE(theElement))
    {
      if ((EstimateHere(theElement)) || (l==covise->max_level))
      {
        int i, coe = CORNERS_OF_ELEM(theElement);
        n_elem++;

        for (i=0; i<coe; i++)
        {
          NODE *theNode = CORNER(theElement,i);

          n_conn++;

          if (USED(theNode)) continue;
          SETUSED(theNode,1);

          if ((SONNODE(theNode)==NULL) || (l==covise->max_level))
          {
                                                #ifdef ModelP
            if (PRIO(theNode) == PrioMaster)
                                                #endif
            n_node++;
          }
        }
      }
    }
  }
        #ifdef ModelP
  n_node = UG_GlobalSumINT(n_node);
  n_elem = UG_GlobalSumINT(n_elem);
  n_conn = UG_GlobalSumINT(n_conn);
        #endif
  covise->n_nodes = n_node;
  covise->n_elems = n_elem;
  covise->n_conns = n_conn;

  return(0);
}



static INT FillCoviseHeader (MULTIGRID *theMG, COVISE_HEADER *covise)
{
  /* extract multigrid statistics */
  /* currently no negative levels, complete grid up to TOPLEVEL */
  covise->min_level = 0;
  covise->max_level = TOPLEVEL(theMG);

  /* extract surface grid statistics */
  ComputeSurfaceGridStats(covise);



  INT num_solutions;       /* number of solutions (may be scalar or vector) */
  SOLUTION_DESC solutions[MAX_SOLUTIONS];       /* solution descriptors */



  return(0);
}






/****************************************************************************/
