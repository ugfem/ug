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

#include <string.h>


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

#define MAX_STRLEN     64
#define MAX_SOLUTIONS  16
#define MAX_COMPONENTS 3


/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/

typedef struct
{
  INT n_components;

  /* external name, for user */
  char name[MAX_STRLEN];

  /* UG components */
  INT comps[MAX_COMPONENTS];

} SOLUTION_DESC;


typedef struct
{
  /* multigrid description */
  INT min_level;        /* coarsest grid level, may be <0 for algebraic methods */
  INT max_level;        /* ==TOPLEVEL in ug */

  /* surface grid description */
  INT n_vertices;       /* number of surface grid vertices */
  INT n_elems;          /* number of surface grid elements */
  INT n_conns;          /* sum of corners over all surface elements */

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
   static PreprocesingProcPtr EvalPreProc;
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


static INT ComputeSurfaceGridStats (MULTIGRID *theMG, COVISE_HEADER *covise)
{
  NODE *theNode;
  INT l, n_vertices, n_elem, n_conn;

  /* surface grid up to current level */
  n_vertices = n_elem = n_conn = 0;
  for (l=covise->min_level; l<=covise->max_level; l++)
  {
    ELEMENT *theElement;
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);

    /* reset USED flags in all vertices to be counted */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      SETUSED(MYVERTEX(theNode),0);
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
          theNode = CORNER(theElement,i);

          n_conn++;

          if (USED(MYVERTEX(theNode))) continue;
          SETUSED(MYVERTEX(theNode),1);

          if ((SONNODE(theNode)==NULL) || (l==covise->max_level))
          {
                                                #ifdef ModelP
            if (PRIO(theNode) == PrioMaster)
                                                #endif
            n_vertices++;
          }
        }
      }
    }
  }
        #ifdef ModelP
  n_vertices = UG_GlobalSumINT(n_vertices);
  n_elem     = UG_GlobalSumINT(n_elem);
  n_conn     = UG_GlobalSumINT(n_conn);
        #endif
  covise->n_vertices = n_vertices;
  covise->n_elems    = n_elem;
  covise->n_conns    = n_conn;

  return(0);
}


static INT GetSolutionDescs (MULTIGRID *theMG, COVISE_HEADER *covise)
{
  /* TODO replace by UI */
  covise->num_solutions = 2;

  strcpy(covise->solutions[0].name, "Concentration");
  covise->solutions[0].n_components = 1;
  covise->solutions[0].comps[0]      = 0;

  strcpy(covise->solutions[1].name, "Pressure");
  covise->solutions[1].n_components = 1;
  covise->solutions[1].comps[0]      = 1;
}


static INT FillCoviseHeader (MULTIGRID *theMG, COVISE_HEADER *covise)
{
  /* extract multigrid statistics */
  /* currently no negative levels, complete grid up to TOPLEVEL */
  covise->min_level = 0;
  covise->max_level = TOPLEVEL(theMG);

  /* extract surface grid statistics */
  ComputeSurfaceGridStats(theMG, covise);

  /* extract solution descriptions */
  GetSolutionDescs(theMG, covise);

  return(0);
}



/****************************************************************************/


static INT ResetVertexFlags (MULTIGRID *theMG, INT min_level, INT max_level)
{
  INT l;

  /* reset used flags in vertices */
  for (l=min_level; l<=max_level; l++)
  {
    NODE *theNode;
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);

    /* reset USED flags in all vertices */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      SETUSED(MYVERTEX(theNode),0);
    }
  }

}


static INT SendSurfaceGrid (MULTIGRID *theMG, INT min_level, INT max_level)
{
  INT l;

  /* send surface vertices, part1: set flags */
  ResetVertexFlags(theMG, min_level, max_level);

  /* send surface vertices, part2: send data */
  for (l=min_level; l<=max_level; l++)
  {
    NODE *theNode;
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);

    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      INT vid;
      DOUBLE *pos;

      if (USED(MYVERTEX(theNode))) continue;
      SETUSED(MYVERTEX(theNode),1);

      /* extract data from vertex */
      vid = ID(MYVERTEX(theNode));
      /* TODO write vid to buffer */
      /* TODO use VXGID in ModelP */

      pos = CVECT(MYVERTEX(theNode));
      /* TODO write pos[DIM] to buffer */
    }
  }



  /* send surface elems and connectivity */
  for (l=min_level; l<=max_level; l++)
  {
    ELEMENT *theElement;
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);

    for (theElement=FIRSTELEMENT(theGrid);
         theElement!=NULL; theElement=SUCCE(theElement))
    {
      if ((EstimateHere(theElement)) || (l==max_level))
      {
        int i, coe = CORNERS_OF_ELEM(theElement);

        /* TODO write coe to buffer (elem) */

        for (i=0; i<coe; i++)
        {
          NODE *theNode = CORNER(theElement,i);
          INT vid;

          vid = ID(MYVERTEX(theNode));
          /* TODO write vid to buffer (conn) */
          /* TODO use VXGID in ModelP */
        }
      }
    }
  }
}



static INT SendSolution (MULTIGRID *theMG, COVISE_HEADER *covise, INT idx_sol)
{
  INT l;

  /* reset vertex flags */
  ResetVertexFlags(theMG, covise->min_level, covise->max_level);


  /* extract data, loop from max_level to min_level! */
  /* TODO: special handling in ModelP */
  for (l=covise->max_level; l>=covise->min_level; l--)
  {
    NODE *theNode;
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);

    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      VECTOR *theVector;
      INT i;

      if (USED(MYVERTEX(theNode))) continue;
      SETUSED(MYVERTEX(theNode),1);

      theVector = NVECTOR(theNode);

      /* extract data from vector */
      for(i=0; i<covise->solutions[idx_sol].n_components; i++)
      {
        INT comp = covise->solutions[idx_sol].comps[i];
        DOUBLE val = VVALUE(theVector,comp);
        /* TODO write val to buffer */
      }
    }
  }
}



/****************************************************************************/


/****************************************************************************/
