// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  coviseif.cpp                                                                                                  */
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


extern "C" {

  /* ug includes */
#define ListElement UG_ListElement
#include "gm.h"
#include "evm.h"
#include "general.h"
#undef ListElement

} /* extern "C" */


/* covise includes */
#include "covise_connect.h"
#include "covise_msg.h"



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

#define MAX_ITEMS_SENT 4096


/* communication constants */
#define DFLT_PORT      31700


enum CommonConsts {
  MT_NONE,
  MT_UGHEADER
};


/****************************************************************************/
/*																			*/
/* data structures used in this source file (exported data structures are	*/
/*		  in the corresponding include file!)								*/
/*																			*/
/****************************************************************************/


/* common data types (for communication UG<->Covise) */
typedef int INT32;
typedef float FLOAT32;



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


/* covise: */
static ClientConnection* covise_connection=NULL;


/* RCS string */
static char RCS_ID("$Header$",UG_RCS_STRING);


/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* functions																*/
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

  return(0);
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


/* TODO renumber doesn't work in parallel */
static INT RenumberVertices (MULTIGRID *theMG, INT min_level, INT max_level)
{
  INT l, i;

  ResetVertexFlags(theMG, min_level, max_level);

  /* reset used flags in vertices */
  i = 0;
  for (l=min_level; l<=max_level; l++)
  {
    NODE *theNode;
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);

    /* reset USED flags in all vertices */
    for (theNode=FIRSTNODE(theGrid); theNode!=NULL; theNode=SUCCN(theNode))
    {
      if (USED(MYVERTEX(theNode))) continue;
      SETUSED(MYVERTEX(theNode),1);

      ID(MYVERTEX(theNode)) = i;
      i++;
    }
  }

  return(0);
}


static INT SendSurfaceGrid (MULTIGRID *theMG, COVISE_HEADER *covise)
{
  INT l, remaining, sent;
  TokenBuffer tb;
  Message* msg = new Message(tb);
  msg->type = (covise_msg_type)0;

  /* renumber vertex IDs */
  /* TODO doesn't work in ModelP */
  RenumberVertices(theMG, covise->min_level, covise->max_level);

  /* send surface vertices, part1: set flags */
  ResetVertexFlags(theMG, covise->min_level, covise->max_level);

  /* send surface vertices, part2: send data */
  sent = 0;

  /* start first buffer */
  tb.reset();
  remaining = MIN(covise->n_vertices,MAX_ITEMS_SENT);
  tb << remaining;

  for (l=covise->min_level; l<=covise->max_level; l++)
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
      /* TODO use VXGID in ModelP */
      vid = ID(MYVERTEX(theNode));
      pos = CVECT(MYVERTEX(theNode));

      tb << (INT32)vid;
      tb << (FLOAT32)pos[0];
      tb << (FLOAT32)pos[1];
      tb << (FLOAT32)pos[2];
      remaining--;
      sent++;

      if (remaining==0)
      {
        /* send this buffer */
        covise_connection->send_msg(msg);

        /* start next buffer */
        tb.reset();
        remaining = MIN(covise->n_vertices - sent, MAX_ITEMS_SENT);
        tb << remaining;
      }
    }
  }


  /* next buffer */
  tb.reset();
  remaining = MIN(covise->n_elems,MAX_ITEMS_SENT);
  tb << remaining;
  sent = 0;

  /* send surface elems and connectivity */
  for (l=covise->min_level; l<=covise->max_level; l++)
  {
    ELEMENT *theElement;
    GRID *theGrid = GRID_ON_LEVEL(theMG,l);

    for (theElement=FIRSTELEMENT(theGrid);
         theElement!=NULL; theElement=SUCCE(theElement))
    {
      if ((EstimateHere(theElement)) || (l==covise->max_level))
      {
        int i, coe = CORNERS_OF_ELEM(theElement);

        tb << (INT32)coe;

        for (i=0; i<coe; i++)
        {
          NODE *theNode = CORNER(theElement,i);
          INT vid;

          vid = ID(MYVERTEX(theNode));
          /* TODO use VXGID in ModelP */
          tb << (INT32)vid;
        }

        remaining--;
        sent++;

        if (remaining==0)
        {
          /* send this buffer */
          covise_connection->send_msg(msg);

          /* start next buffer */
          tb.reset();
          remaining = MIN(covise->n_elems - sent, MAX_ITEMS_SENT);
          tb << remaining;
        }
      }
    }
  }

  /* cleanup */
  delete msg;

  return(0);
}



static INT SendSolution (MULTIGRID *theMG, COVISE_HEADER *covise, INT idx_sol)
{
  INT l, remaining, sent, n_comp_sol;
  TokenBuffer tb;
  Message* msg = new Message(tb);
  msg->type = (covise_msg_type)0;

  n_comp_sol = covise->solutions[idx_sol].n_components;

  /* reset vertex flags */
  ResetVertexFlags(theMG, covise->min_level, covise->max_level);


  /* start first buffer */
  tb.reset();
  remaining = MIN(covise->n_vertices,MAX_ITEMS_SENT);
  tb << remaining;
  tb << idx_sol;
  sent = 0;

  for (l=covise->min_level; l<=covise->max_level; l++)

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
        for(i=0; i<n_comp_sol; i++)
        {
          INT comp = covise->solutions[idx_sol].comps[i];
          tb << (FLOAT32) VVALUE(theVector,comp);
        }

        remaining--;
        sent++;

        if (remaining==0)
        {
          /* send this buffer */
          covise_connection->send_msg(msg);

          /* start next buffer */
          tb.reset();
          remaining = MIN(covise->n_vertices - sent, MAX_ITEMS_SENT);
          tb << remaining;
          tb << idx_sol;
        }
      }
    }

  delete msg;
  return(0);
}



/****************************************************************************/


/****************************************************************************/



extern "C" INT ConnectCovise (MULTIGRID *theMG, char *hostname)
{
  COVISE_HEADER covise;
  Host *remotehost=new Host(hostname);

  printf("CoviseIF: trying to connect to %s\n", hostname);

  if (covise_connection==NULL)
  {
    covise_connection =
      new ClientConnection(remotehost, DFLT_PORT, 0, (sender_type)0, 0);
    printf("CoviseIF: new ClientConnection =%08x\n");

    if (! covise_connection->is_connected())
    {
      delete covise_connection;
      covise_connection=NULL;
    }
  }

  if (covise_connection!=NULL && covise_connection->is_connected())
  {
    printf("CoviseIF: ok, connected to %s\n", hostname);
    Message* msg = new Message();
    msg->type = (covise_msg_type)0;
    int i;

    if (covise_connection->check_for_input())
    {
      printf("CoviseIF: input from %s\n", hostname);
      covise_connection->recv_msg(msg);
      printf("CoviseIF: received msg from %s, type %d\n", hostname, msg->type);

      switch (msg->type)
      {
      case (covise_msg_type)0 :
      {
        TokenBuffer tb(msg);
        int headerRequest,gridRequest,numSolutions,solutions[MAX_SOLUTIONS];
        tb >> headerRequest;
        tb >> gridRequest;
        tb >> numSolutions;
        for(i=0; i<numSolutions; i++)
          tb >> solutions[i];

        if (headerRequest)
        {
          printf("CoviseIF: headerRequest\n");
          TokenBuffer rtb;
          FillCoviseHeader(theMG, &covise);
          rtb << MT_UGHEADER;
          rtb << covise.min_level;
          rtb << covise.max_level;
          rtb << covise.n_vertices;
          rtb << covise.n_elems;
          rtb << covise.n_conns;
          rtb << covise.num_solutions;
          for(i=0; i<numSolutions; i++)
          {
            rtb << covise.solutions[i].n_components;
            rtb << covise.solutions[i].name;
          }

          Message rmsg(rtb);
          rmsg.type = (covise_msg_type)0;
          covise_connection->send_msg(&rmsg);
        }

        if (gridRequest)
        {
          printf("CoviseIF: gridRequest\n");
          SendSurfaceGrid(theMG, &covise);
        }

        if (numSolutions>0)
        {
          printf("CoviseIF: solRequest\n");
          for(i=0; i<numSolutions; i++)
          {
            SendSolution(theMG,&covise,solutions[i]);
          }
        }
      }
      break;

      case SOCKET_CLOSED :
      case CLOSE_SOCKET :
      case EMPTY :
      {
        printf("CoviseIF: SOCKET_CLOSED\n");
        delete covise_connection;
        covise_connection = NULL;
      }
      break;

      case QUIT :
      {
        printf("CoviseIF: QUIT\n");
        delete covise_connection;
        covise_connection = NULL;
        /* TODO exit UG itself */
      }
      break;
      }
    }

    delete msg;
  }


  /* cleanup */
  delete remotehost;

  return(0);       /* no error */
}



extern "C" INT InitCoviseIF (void)
{
  return(0);       /* no error */
}


extern "C" INT ExitCoviseIF (void)
{
  if (covise_connection!=NULL)
  {
    delete covise_connection;
    covise_connection=NULL;
  }

  return(0);       /* no error */
}


/****************************************************************************/
