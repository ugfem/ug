// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initddd.c														*/
/*																			*/
/* Purpose:   register ug structs for handling them by ddd					*/
/*																			*/
/* Author:	  Stefan Lang													*/
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70550 Stuttgart												*/
/*			  email: stefan@ica3.uni-stuttgart.de							*/
/*			  phone: 0049-(0)711-685-7003									*/
/*			  fax  : 0049-(0)711-685-7000									*/
/*																			*/
/* History:   09.05.95 begin, ugp version 3.0								*/
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

#include <stdlib.h>

#include "debug.h"

#include "parallel.h"
#include "gm.h"

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

/* DDD objects */
DDD_TYPE TypeVector;
DDD_TYPE TypeIVertex, TypeBVertex;
DDD_TYPE TypeNode;

#ifdef __TWODIM__
DDD_TYPE TypeTrElem, TypeTrBElem,
         TypeQuElem, TypeQuBElem;
#endif

#ifdef __THREEDIM__
DDD_TYPE TypeTeElem, TypeTeBElem;
#endif

/* DDD data objects */
DDD_TYPE TypeConnection;
DDD_TYPE TypeVSegment;
DDD_TYPE TypeEdge;
DDD_TYPE TypeElementSide;

/* DDD interfaces needed for distributed computation */
INTERFACE ElementIF, NodeIF;

/****************************************************************************/
/*																			*/
/* definition of variables global to this source file only (static!)		*/
/*																			*/
/****************************************************************************/

/* data for CVS */
static char rcsid[] = "$Header$";

/****************************************************************************/
/*																			*/
/* forward declarations of functions used before they are defined			*/
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* Function:  ddd_structInit                                                                                            */
/*																			*/
/* Purpose:   register ug data structures, which should be handled by       */
/*            DDD, to DDD                                                   */
/*																			*/
/* Input:     void                                                              */
/*																			*/
/* Output:    void                                                                                                              */
/*																			*/
/*																			*/
/****************************************************************************/

void ddd_structInit (void)
{
  VECTOR v;
  NODE n;
  struct ivertex iv;
  struct bvertex bv;
  struct triangle tr;
  struct quadrilateral qu;
  struct tetrahedron te;
  MATRIX m;
  VSEGMENT vs;
  EDGE e;
  ELEMENTSIDE es;


  /****************************************************************************/
  /*																			*/
  /*  register DDD objects													*/
  /*																			*/
  /****************************************************************************/

  TypeVector = DDD_StructRegister("Vector", &v.ddd_hdr, &v,
                                  EL_GDATA,  &v.control,
                                  EL_OBJPTR, &v.object, 1, /* bisher: EL_LDATA! */
                                  EL_LDATA,  &v.pred,
                                  EL_GDATA,  &v.index,
                                  EL_LDATA,  &v.start,

                                  /* TODO: value wird noch ausgelassen. feld variabler laenge? */
                                  /* bei entscheidung 'value': kein weiteres feld
                                          bei ent. 'userdata *': EL_GDATA-feld        */

                                  EL_END,    &v+1
                                  );

  TypeIVertex = DDD_StructRegister("IVertex", &iv.ddd_hdr, &iv,
                                   EL_GDATA,  &iv.control,
                                   EL_GDATA,  &iv.id, /* TODO: muss umgerechnet werden! */
                                   EL_GDATA,  &iv.x[0],
                                   EL_LDATA,  &iv.pred,
                                   /*		EL_OBJPTR, &iv.father,	1,  bisher vergessen? oder doch EL_LDATA */
                                   EL_OBJPTR, &iv.topnode, 1, /* wirklich OBJPTR? */
                                   EL_END,    &iv+1
                                   );

  TypeBVertex = DDD_StructRegister("BVertex", &bv.ddd_hdr, &bv,
                                   EL_GDATA,  &bv.control,
                                   EL_GDATA,  &bv.id, /* TODO: muss umgerechnet werden! */
                                   EL_GDATA,  &bv.x[0],
                                   EL_LDATA,  &bv.pred,
                                   /*		EL_OBJPTR, &iv.father,	1,  bisher vergessen? oder doch EL_LDATA */
                                   EL_OBJPTR, &bv.topnode, 1, /* wirklich OBJPTR? */
                                   EL_LDATA,  &bv.vseg, /* different from IVertex */
                                   EL_END,    &bv+1
                                   );

  TypeNode = DDD_StructRegister("Node", &n.ddd_hdr, &n,
                                EL_GDATA,  &n.control,
                                EL_GDATA,  &n.id, /* TODO: muss umgerechnet werden! */
                                EL_LDATA,  &n.pred,
                                EL_OBJPTR, &n.father, 3,

                #ifdef __NODEDATA__
                                EL_OBJPTR, &n.vector, 1, /* associated vector                 */
                #endif

                                EL_END,    &n+1
                                );

  /* CAUTION: elements must be registered in this position     */
  /*          because of direct mapping form TAG to DDD_TYPE   */

        #ifdef __TWODIM__
  TypeTrElem = DDD_StructRegister("TrElem", &tr.ddd_hdr, &tr,
                                  EL_GDATA,  &tr.control,
                                  EL_GDATA,  &tr.id, /* TODO: muss umgerechnet werden! */
                                  EL_LDATA,  &tr.pred,
                                  EL_OBJPTR, &tr.n,               3,
                                  EL_OBJPTR, &tr.father,  1,
                                  EL_OBJPTR, &tr.sons,    4,
                                  EL_OBJPTR, &tr.nb,              3,

                #ifdef __ELEMDATA__
                                  EL_OBJPTR, &tr.vector,  1,
                #endif

                                  EL_END,    ((char *)&tr)+(sizeof(struct triangle)-3*sizeof(struct elementside *))
                                  );

  TypeTrBElem = DDD_StructRegister("TrBElem", &tr.ddd_hdr, &tr,
                                   EL_GDATA,  &tr.control,
                                   EL_GDATA,  &tr.id, /* vorsicht! "global" id?? */
                                   EL_LDATA,  &tr.pred,
                                   EL_OBJPTR, &tr.n,               3,
                                   EL_OBJPTR, &tr.father,  1,
                                   EL_OBJPTR, &tr.sons,    4,
                                   EL_OBJPTR, &tr.nb,              3,

                #ifdef __ELEMDATA__
                                   EL_OBJPTR, &tr.vector,  1,
                #endif

                                   EL_LDATA,  &tr.side, /* this is a boundary TrElem */
                                   EL_END,    &tr+1
                                   );

  TypeQuElem = DDD_StructRegister("QuElem", &qu.ddd_hdr, &qu,
                                  EL_GDATA,  &qu.control,
                                  EL_GDATA,  &qu.id, /* vorsicht! "global" id?? */
                                  EL_LDATA,  &qu.pred,
                                  EL_OBJPTR, &qu.n,               4,
                                  EL_OBJPTR, &qu.father,  1,
                                  EL_OBJPTR, &qu.sons,    4,
                                  EL_OBJPTR, &qu.nb,              4,

                #ifdef __ELEMDATA__
                                  EL_OBJPTR, &qu.vector,  1,
                #endif

                                  EL_END,    ((char *)&qu)+(sizeof(struct quadrilateral)-4*sizeof(struct elementside *))
                                  );

  TypeQuBElem = DDD_StructRegister("QuBElem", &qu.ddd_hdr, &qu,
                                   EL_GDATA,  &qu.control,
                                   EL_GDATA,  &qu.id, /* vorsicht! "global" id?? */
                                   EL_LDATA,  &qu.pred,
                                   EL_OBJPTR, &qu.n,               4,
                                   EL_OBJPTR, &qu.father,  1,
                                   EL_OBJPTR, &qu.sons,    4,
                                   EL_OBJPTR, &qu.nb,              4,

                #ifdef __ELEMDATA__
                                   EL_OBJPTR, &qu.vector,  1,
                #endif

                                   EL_LDATA,  &qu.side, /* this is a boundary QuElem */
                                   EL_END,    &qu+1
                                   );
        #endif /* TWODIM */

        #ifdef __THREEDIM__
  TypeTeElem = DDD_StructRegister("TeElem", &te.ddd_hdr, &te,
                                  EL_GDATA,  &te.control,
                                  EL_GDATA,  &te.id, /* vorsicht! "global" id?? */
                                  EL_LDATA,  &te.pred,
                                  EL_OBJPTR, &te.n,                       4,
                                  EL_OBJPTR, &te.father,          1,
                                  EL_OBJPTR, &te.sons,            1,
                                  EL_OBJPTR, &te.nb,                      4,

                #ifdef __ELEMDATA__
                                  EL_OBJPTR, &te.vector,          1,
                #endif

                #ifdef __SIDEDATA__
                                  EL_OBJPTR, &te.sidevector,      4,
                #endif

                                  EL_END,    ((char *)&te)+(sizeof(struct tetrahedron)-4*sizeof(struct elementside *))
                                  );

  TypeTeBElem = DDD_StructRegister("TeBElem", &te.ddd_hdr, &te,
                                   EL_GDATA,  &te.control,
                                   EL_GDATA,  &te.id, /* vorsicht! "global" id?? */
                                   EL_LDATA,  &te.pred,
                                   EL_OBJPTR, &te.n,                       4,
                                   EL_OBJPTR, &te.father,          1,
                                   EL_OBJPTR, &te.sons,            1,
                                   EL_OBJPTR, &te.nb,                      4,

                #ifdef __ELEMDATA__
                                   EL_OBJPTR, &te.vector,          1,
                #endif

                #ifdef __SIDEDATA__
                                   EL_OBJPTR,  &te.sidevector, 4,
                #endif

                                   EL_LDATA,  &te.side, /* this is a boundary TeElem */
                                   EL_END,    &te+1
                                   );
        #endif /* THREEDIM */


  /****************************************************************************/
  /*																			*/
  /*  register DDD data objects												*/
  /*																			*/
  /****************************************************************************/

  TypeConnection = DDD_StructRegister("Connection", NULL, &m,
                                      EL_GDATA,  &m.control,
                                      EL_LDATA,  &m.next,
                                      EL_OBJPTR, &m.vect,   1,

                                      /* value[1] wird definitiv weggelassen. */

                                      EL_END,    &m+1
                                      );


  TypeVSegment = DDD_StructRegister("VSegment", NULL, &vs,
                                    EL_GDATA,  &vs.control,
                                    EL_LDATA,  &vs.segdesc,
                                    EL_GDATA,  &vs.lambda,
                                    EL_LDATA,  &vs.next,
                                    EL_GDATA,  &vs+1, /* add a dummy containing segment id */
                                    EL_END,    ((char *)&vs)+(sizeof(VSEGMENT)+sizeof(INT))
                                    );


  TypeEdge = DDD_StructRegister("Edge", NULL, &e,
                                /* link 0 data */
                                EL_GDATA,  &e.links[0].control,
                                EL_LDATA,  &e.links[0].next,
                                EL_OBJPTR, &e.links[0].nbnode, 1,
                #ifdef __version23__
                                EL_LDATA,  &e.links[0].matelem,
                #endif

                                /* link 1 data */
                                EL_GDATA,  &e.links[1].control,
                                EL_LDATA,  &e.links[1].next,
                                EL_OBJPTR, &e.links[1].nbnode, 1,
                #ifdef __version23__
                                EL_LDATA,  &e.links[1].matelem,
                #endif

                #ifdef __MIDNODE__
                                EL_OBJPTR, &e.midnode,  1,
                #endif

                #ifdef __EDGEDATA__
                                EL_LDATA,  &e.vector,  1,
                #endif

                #ifdef __version23__
                                EL_LDATA,  &e.data,
                #endif

                                EL_END,    &e+1
                                );

  TypeElementSide = DDD_StructRegister("ElementSide", NULL, &es,
                                       EL_GDATA,  &es.control,
                                       EL_LDATA,  &es.segdesc,
                                       EL_GDATA,  &es.lambda[0][0],
                                       EL_LDATA,  &es.pred,
                                       EL_LDATA,  &es.succ,
                                       EL_GDATA,  &es+1, /* add a dummy containing segment id */
                                       EL_END,    ((char *)&es)+(sizeof(ELEMENTSIDE)+sizeof(INT))
                                       );
}


/****************************************************************************/
/*																			*/
/* Function:  ddd_ifInit                                                                                                        */
/*																			*/
/* Purpose:   defines the communication interfaces needed in ug for             */
/*		      management by DDD
   /*																			*/
/* Input:     void                                                              */
/*																			*/
/* Output:    void                                                                                                              */
/*																			*/
/*																			*/
/****************************************************************************/

#ifdef __TWODIM__
void ddd_ifInit (void)
{
  DDD_TYPE O[16];
  int A[16];
  int B[16];

  O[0] = TypeTrElem;
  O[1] = TypeTrElem;
  A[0] = 0; A[1] = 7;
  B[0] = 0; B[1] = 7;
  ElementIF = DDD_IFDefine(2,O,2,A,2,B);

  O[0] = TypeNode;
  A[0] = 0; A[1] = 7;
  B[0] = 0; B[1] = 7;
  NodeIF = DDD_IFDefine(1,O,2,A,2,B);
}
#endif


/****************************************************************************/
/*																			*/
/* Function:  InitParallel                                                                                              */
/*																			*/
/* Purpose:   initializes the ddd library by defining the ug internal       */
/*            issues:                                                       */
/*                      - format of handled structs                         */
/*                      - description of handlers                           */
/*						- definition of interfaces							*/
/*																			*/
/* Input:     int    argc: number of arguments                                  */
/*		      char **argv: list of argument pointers						*/
/*																			*/
/* Output:    int:   error value											*/
/*																			*/
/*																			*/
/****************************************************************************/


int InitParallel (int argc, char **argv)
{
  DDD_Init(argc, argv);

  /* register ddd data types and dependent types */
  ddd_structInit();

  /* display ddd types */
  IFDEBUG(dddif,1)
  DDD_StructDisplay(TypeVector);
  DDD_StructDisplay(TypeIVertex);
  DDD_StructDisplay(TypeBVertex);
  DDD_StructDisplay(TypeNode);

        #ifdef __TWODIM__
  DDD_StructDisplay(TypeTrElem);
  DDD_StructDisplay(TypeTrBElem);
  DDD_StructDisplay(TypeQuElem);
  DDD_StructDisplay(TypeQuBElem);
        #endif

        #ifdef __THREEDIM__
  DDD_StructDisplay(TypeTeElem);
  DDD_StructDisplay(TypeTeBElem);
        #endif

  /* display dependent types */
  DDD_StructDisplay(TypeConnection);
  DDD_StructDisplay(TypeVSegment);
  DDD_StructDisplay(TypeEdge);
  DDD_StructDisplay(TypeElementSide);
  ENDDEBUG

  ddd_HandlerInit();

        #ifdef __TWODIM__
  ddd_ifInit();
        #endif

  return 0;          /* no error */
}


/****************************************************************************/
/*																			*/
/* Function:  ExitParallel                                                                                              */
/*																			*/
/* Purpose:   exit the parallel application on ddd level                    */
/*																			*/
/* Input:     void                                                              */
/*																			*/
/* Output:    void                                                                                                              */
/*																			*/
/*																			*/
/****************************************************************************/

void ExitParallel (void)
{
  DDD_Exit();
}

#endif /* ModelP */
