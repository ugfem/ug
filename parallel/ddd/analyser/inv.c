// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*                                                                          */
/* File:      inventor.c                                                    */
/*                                                                          */
/* Purpose:   visualization in IRIS Inventor format                         */
/*                                                                          */
/* Author:    Klaus Birken                                                  */
/*            Rechenzentrum Uni Stuttgart                                   */
/*            Universitaet Stuttgart                                        */
/*            Allmandring 30                                                */
/*            70550 Stuttgart                                               */
/*            internet: birken@rus.uni-stuttgart.de                         */
/*                                                                          */
/* History:   95/03/28 kb  begin                                            */
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
#include <stdio.h>
#include <math.h>


#include "dddi.h"




/****************************************************************************/
/*                                                                          */
/* data structures                                                          */
/*                                                                          */
/****************************************************************************/


/* graph of DDD_TYPEs, with directed edges according to references */

typedef struct _TYPE_NODE
{
  TYPE_DESC *def;                  /* pointer to DDD_TYPE */
  struct _TYPE_EDGE *refs;         /* linked list of referenced types */

} TYPE_NODE;


typedef struct _TYPE_EDGE
{
  DDD_TYPE reftype;        /* referenced type */
  int n;                   /* number of references */

  struct _TYPE_EDGE *next;        /* linked list */
} TYPE_EDGE;



/****************************************************************************/
/*                                                                          */
/* local variables                                                          */
/*                                                                          */
/****************************************************************************/

/* Revision Control System string */
RCSID("$Header$",DDD_RCS_STRING)


/* this array corresponds to theTypeDefs */
static TYPE_NODE theTypeNodes[MAX_TYPEDESC];

/* this array is sorted due to reference relations */
static TYPE_NODE *theTypeGraph[MAX_TYPEDESC];



/****************************************************************************/
/*                                                                          */
/* subroutines                                                              */
/*                                                                          */
/****************************************************************************/


static TYPE_EDGE *GetTypeEdge (TYPE_NODE *tn, DDD_TYPE reftype)
{
  TYPE_EDGE *te;

  for(te=tn->refs; te!=NULL && te->reftype!=reftype; te=te->next)
    ;

  if (te==NULL)
  {
    te = (TYPE_EDGE *)AllocTmp(sizeof(TYPE_EDGE));
    te->reftype = reftype;
    te->n = 0;

    te->next = tn->refs;
    tn->refs = te;
  }

  return(te);
}


static void AnalyseTypes (void)
{
  int i;

  /* create graph of DDD_TYPE ref-relations */
  for(i=0; i<DDD_InfoTypes(); i++)
  {
    TYPE_DESC *td = &(theTypeDefs[i]);
    TYPE_NODE *tn = &(theTypeNodes[i]);
    int e;

    tn->def = td;
    tn->refs = NULL;

    theTypeGraph[i] = tn;

    for(e=0; e<td->nElements; e++)
    {
      ELEM_DESC *el = &(td->element[e]);

      if (el->type==EL_OBJPTR)
      {
        TYPE_EDGE *te = GetTypeEdge(tn, EDESC_REFTYPE(el));
        te->n += (el->size / sizeof(void *));
      }
    }

    printf("%4d: type %s (%03d) refs:\n", me, theTypeDefs[i].name, i);
    {
      TYPE_EDGE *te;
      for(te=tn->refs; te!=NULL; te=te->next)
      {
        printf("         %s (%03d), n=%d\n",
               theTypeDefs[te->reftype].name, te->reftype, te->n);
      }
    }
  }

  /* sort type-graph nodes according to mutual references */
  /*
          qsort(theTypeGraph, DDD_InfoTypes
   */
}



/****************************************************************************/


void DDD_GraphicalAnalyser (char *filename)
{
  FILE *fp;

  fp = fopen(filename, "w");

  if (me==0)
  {
    AnalyseTypes();
  }

  fclose(fp);
}
