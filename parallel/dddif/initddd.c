// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initddd.c														*/
/*																			*/
/* Purpose:   register ug structs for handling them by ddd					*/
/*																			*/
/* Author:	  Stefan Lang, Klaus Birken										*/
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


/****************************************************************************/
/*																			*/
/* defines in the following order											*/
/*																			*/
/*		  compile time constants defining static data size (i.e. arrays)	*/
/*		  other constants													*/
/*		  macros															*/
/*																			*/
/****************************************************************************/


/* macro for easier definition of DDD_TYPEs */
#define  ELDEF(comp)    &(comp),sizeof(comp)



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
DDD_TYPE TypeTrElem, TypeTrBElem;
DDD_TYPE TypeQuElem, TypeQuBElem;
#endif

#ifdef __THREEDIM__
DDD_TYPE TypeTeElem, TypeTeBElem;
DDD_TYPE TypePyElem, TypePyBElem;
DDD_TYPE TypeHeElem, TypeHeBElem;
#endif


/* DDD data objects */
DDD_TYPE TypeMatrix;
DDD_TYPE TypeVSegment;
DDD_TYPE TypeEdge;
DDD_TYPE TypeElementSide;

/* DDD interfaces needed for distributed computation */
DDD_IF ElementIF, NodeIF;


/* DDD global controls */
DDD_CTRL dddctrl;


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

enum ElemTypeFlag { Inside, Boundary };

static void ddd_InitGenericElement (INT tag, DDD_TYPE dddType, int etype)
{
  struct generic_element *ge=0;
  GENERAL_ELEMENT *desc = element_descriptors[tag];

  size_t ps  = sizeof(void *);
  void   **r = ge->refs;


  /* initialize base part (valid for all elements) */
  DDD_TypeDefine(dddType, ge,
                 EL_DDDHDR, &(ge->ddd),
                 EL_GDATA,  ELDEF(ge->control),

                 /* TODO: id muss umgerechnet werden! (?) */
                 EL_GDATA,  ELDEF(ge->id),
                 EL_GDATA,  ELDEF(ge->flag),
                 EL_GDATA,  ELDEF(ge->property),
                 EL_LDATA,  ELDEF(ge->pred),
                 EL_LDATA,  ELDEF(ge->succ),
                 EL_CONTINUE);


  /* initialize generic part */
  /* NOTE: references to 'union element' are denoted by the ref-type
     dddType (i.e., the DDD_TYPE of the currently defined element itself).
     TODO: this should be replaced by a more explicit TypeGenericElement in
     later versions (code would be more readable). */

  DDD_TypeDefine(dddType, ge,
                 EL_OBJPTR, r+n_offset[tag],       ps*desc->corners_of_elem, TypeNode,
                 EL_OBJPTR, r+father_offset[tag],  ps*1,                     dddType,
                 EL_OBJPTR, r+sons_offset[tag],    ps*desc->max_sons_of_elem,dddType,
                 EL_OBJPTR, r+nb_offset[tag],      ps*desc->sides_of_elem,   dddType,
                 EL_CONTINUE);


  /* optional components */

  if (dddctrl.elemData)
    DDD_TypeDefine(dddType, ge,
                   EL_OBJPTR, r+evector_offset[tag], ps*1,     TypeVector,
                   EL_CONTINUE);

  if (dddctrl.sideData)
    DDD_TypeDefine(dddType, ge,
                   EL_OBJPTR, r+svector_offset[tag], ps*desc->sides_of_elem, TypeVector,
                   EL_CONTINUE);

  if (etype==Inside)
  {
    DDD_TypeDefine(dddType, ge, EL_END, desc->inner_size);

    /* init type mapping arrays */
    dddctrl.ugtypes[dddType] = MAPPED_INNER_OBJT(tag);
    dddctrl.types[MAPPED_INNER_OBJT(tag)] = dddType;
    dddctrl.dddObj[MAPPED_INNER_OBJT(tag)] = TRUE;
  }
  else
  {
    DDD_TypeDefine(dddType, ge,
                   EL_LDATA, r+side_offset[tag],  ps*desc->sides_of_elem,
                   EL_END, desc->bnd_size);

    /* init type mapping arrays */
    dddctrl.ugtypes[dddType] = MAPPED_BND_OBJT(tag);
    dddctrl.types[MAPPED_BND_OBJT(tag)] = dddType;
    dddctrl.dddObj[MAPPED_BND_OBJT(tag)] = TRUE;
  }

}


/****************************************************************************/
/*																			*/
/* Function:  ddd_DeclareTypes												*/
/*																			*/
/* Purpose:   declare ug data structures as DDD_TYPES						*/
/*																			*/
/* Input:     void                                                              */
/*																			*/
/* Output:    void                                                                                                              */
/*																			*/
/*																			*/
/****************************************************************************/

static void ddd_DeclareTypes (void)
{
  /* NOTE: (960410 KB)

          - handling of Vector and Matrix types may be
            wrong, types array entries might be used differently
            by ug (e.g., during allocation/deallocation). TODO ueberpruefen!
            (beziehung algebra.c <-> ugm.c)

          - UGTypes for elements will be computed later. therefore the
            initialization of types[] entries for elements
            is done (after the UGTypes are known) during ddd_InitGenericElement.
            TODO remove this clumsy exception! how?

          - variables TypeXXXX containing the proper DDD_TYPEs may be
            superfluous. it would be possible to replace all occurences
            by macros like DDDType(VEOBJ), which would be implemented as
     #define DDDType(ugtype)  (dddctrl.types[ugtype])
            TODO check this!
            pros: no double information. currently, TypeXXX may differ
                  from corresponding dddctrl.types[] entry.
            cons: will this be compatible with alloc/dealloc and TypeDefine
                  of ug-general-elements?
   */


  /* 1. DDD objects (with DDD_HEADER) */

  TypeVector      = DDD_TypeDeclare("Vector");
  dddctrl.ugtypes[TypeVector] = VEOBJ;
  dddctrl.types[VEOBJ] = TypeVector;
  dddctrl.dddObj[VEOBJ] = TRUE;

  TypeIVertex     = DDD_TypeDeclare("IVertex");
  dddctrl.ugtypes[TypeIVertex] = IVOBJ;
  dddctrl.types[IVOBJ] = TypeIVertex;
  dddctrl.dddObj[IVOBJ] = TRUE;

  TypeBVertex     = DDD_TypeDeclare("BVertex");
  dddctrl.ugtypes[TypeBVertex] = BVOBJ;
  dddctrl.types[BVOBJ] = TypeBVertex;
  dddctrl.dddObj[BVOBJ] = TRUE;

  TypeNode        = DDD_TypeDeclare("Node");
  dddctrl.ugtypes[TypeNode] = NDOBJ;
  dddctrl.types[NDOBJ] = TypeNode;
  dddctrl.dddObj[NDOBJ] = TRUE;


        #ifdef __TWODIM__
  TypeTrElem      = DDD_TypeDeclare("TrElem");
  TypeTrBElem     = DDD_TypeDeclare("TrBElem");
  TypeQuElem      = DDD_TypeDeclare("QuElem");
  TypeQuBElem     = DDD_TypeDeclare("QuBElem");
        #endif /* TWODIM */

        #ifdef __THREEDIM__
  TypeTeElem      = DDD_TypeDeclare("TeElem");
  TypeTeBElem     = DDD_TypeDeclare("TeBElem");
  TypePyElem      = DDD_TypeDeclare("PyElem");
  TypePyBElem     = DDD_TypeDeclare("PyBElem");
  TypeHeElem      = DDD_TypeDeclare("HeElem");
  TypeHeBElem     = DDD_TypeDeclare("HeBElem");
        #endif /* THREEDIM */


  /* 2. DDD data objects (without DDD_HEADER) */

  TypeMatrix  = DDD_TypeDeclare("Matrix");
  dddctrl.ugtypes[TypeMatrix] = MAOBJ;
  dddctrl.types[MAOBJ] = TypeMatrix;

  TypeVSegment    = DDD_TypeDeclare("VSegment");
  dddctrl.ugtypes[TypeVSegment] = VSOBJ;
  dddctrl.types[VSOBJ] = TypeVSegment;

  TypeEdge        = DDD_TypeDeclare("Edge");
  dddctrl.ugtypes[TypeEdge] = EDOBJ;
  dddctrl.types[EDOBJ] = TypeEdge;

  TypeElementSide = DDD_TypeDeclare("ElementSide");
  dddctrl.ugtypes[TypeElementSide] = ESOBJ;
  dddctrl.types[ESOBJ] = TypeElementSide;
}


/****************************************************************************/
/*																			*/
/* Function:  ddd_DefineTypes												*/
/*																			*/
/* Purpose:   define previously declared DDD_TYPES							*/
/*            Note: this function depends on previous definition of all		*/
/*                  necessary ug-generic-elements.				                        */
/*																			*/
/* Input:     void                                                              */
/*																			*/
/* Output:    void                                                                                                              */
/*																			*/
/*																			*/
/****************************************************************************/

static void ddd_DefineTypes (void)
{
  VECTOR v;
  NODE n;
  struct ivertex iv;
  struct bvertex bv;

  MATRIX m;
  VSEGMENT vs;
  EDGE e;
  ELEMENTSIDE es;


  /* 1. DDD objects (with DDD_HEADER) */

  DDD_TypeDefine(TypeVector, &v,
                 EL_DDDHDR, &v.ddd,
                 EL_GDATA,  ELDEF(v.control),

                 /* object must be LDATA, because reftype may be a non-DDD-object */
                 /* (e.g., edge). therefore, 'object' must be updated by MKCONS-  */
                 /* handler of associated object. 960404 KB */
                 EL_LDATA, ELDEF(v.object),
                 EL_LDATA,  ELDEF(v.pred),
                 EL_LDATA,  ELDEF(v.succ),
                 EL_GDATA,  ELDEF(v.index),
                 EL_GDATA,  ELDEF(v.skip),
                 EL_LDATA,  ELDEF(v.start),

                 /* TODO is this LDATA or GDATA? */
                 EL_LDATA,  ELDEF(v.block_descr),

                #ifdef __INTERPOLATION_MATRIX__
                 EL_LDATA,  ELDEF(v.istart),
                #endif

                 /* TODO: value wird noch ausgelassen. feld variabler laenge? */
                 /* bei entscheidung 'value': kein weiteres feld
                         bei ent. 'userdata *': EL_GDATA-feld        */
                 EL_END,    &v+1
                 );


  DDD_TypeDefine(TypeIVertex, &iv,
                 EL_DDDHDR, &iv.ddd,
                 EL_GDATA,  ELDEF(iv.control),

                 /* TODO: muss umgerechnet werden! */
                 EL_GDATA,  ELDEF(iv.id),
                 EL_GDATA,  ELDEF(iv.x),
                 EL_GDATA,  ELDEF(iv.xi),
                 EL_LDATA,  ELDEF(iv.pred),
                 EL_LDATA,  ELDEF(iv.succ),
                 EL_LDATA,  ELDEF(iv.data),

                 /* TODO muss father LDATA oder OBJPTR sein? */
                 EL_LDATA,  ELDEF(iv.father),

                 /* TODO topnode wirklich OBJPTR? */
                 EL_OBJPTR, ELDEF(iv.topnode), TypeNode,
                 EL_END,    &iv+1
                 );


  DDD_TypeDefine(TypeBVertex, &bv,
                 EL_DDDHDR, &bv.ddd,
                 EL_GDATA,  ELDEF(bv.control),

                 /* TODO: muss umgerechnet werden! */
                 EL_GDATA,  ELDEF(bv.id),
                 EL_GDATA,  ELDEF(bv.x),
                 EL_GDATA,  ELDEF(bv.xi),
                 EL_LDATA,  ELDEF(bv.pred),
                 EL_LDATA,  ELDEF(bv.succ),
                 EL_LDATA,  ELDEF(bv.data),

                 /* TODO muss father LDATA oder OBJPTR sein? */
                 EL_LDATA,  ELDEF(bv.father),

                 /* TODO topnode wirklich OBJPTR? */
                 EL_OBJPTR, ELDEF(bv.topnode), TypeNode,
                 EL_LDATA,  ELDEF(bv.vseg),     /* different from IVertex */
                 EL_END,    &bv+1
                 );


  DDD_TypeDefine(TypeNode, &n,
                 EL_DDDHDR, &n.ddd,
                 EL_GDATA,  ELDEF(n.control),

                 /* TODO: muss umgerechnet werden! */
                 EL_GDATA,  ELDEF(n.id),
                 EL_GDATA,  ELDEF(n.index),
                 EL_LDATA,  ELDEF(n.pred),
                 EL_LDATA,  ELDEF(n.succ),

                 /* TODO was ist start? */
                 EL_LDATA,  ELDEF(n.start),
                 EL_OBJPTR, ELDEF(n.father),   TypeNode,
                 EL_OBJPTR, ELDEF(n.son),      TypeNode,

                 /* TODO: ref-typ muss eigentlich {TypeIVertex,TypeBVertex} sein! */
                 EL_OBJPTR, ELDEF(n.myvertex), TypeIVertex,
                 EL_CONTINUE);

  if (dddctrl.nodeData)
    DDD_TypeDefine(TypeNode, &n,
                   EL_OBJPTR, ELDEF(n.vector), TypeVector,
                   EL_CONTINUE);

  DDD_TypeDefine(TypeNode, &n,
                 EL_END,    &n+1);


        #ifdef __TWODIM__
  ddd_InitGenericElement(TRIANGLE,      TypeTrElem,  Inside);
  ddd_InitGenericElement(TRIANGLE,      TypeTrBElem, Boundary);
  ddd_InitGenericElement(QUADRILATERAL, TypeQuElem,  Inside);
  ddd_InitGenericElement(QUADRILATERAL, TypeQuBElem, Boundary);
        #endif /* TWODIM */


        #ifdef __THREEDIM__
  ddd_InitGenericElement(TETRAHEDRON, TypeTeElem,  Inside);
  ddd_InitGenericElement(TETRAHEDRON, TypeTeBElem, Boundary);
  ddd_InitGenericElement(PYRAMID,     TypePyElem,  Inside);
  ddd_InitGenericElement(PYRAMID,     TypePyBElem, Boundary);
  ddd_InitGenericElement(HEXAHEDRON,  TypeHeElem,  Inside);
  ddd_InitGenericElement(HEXAHEDRON,  TypeHeBElem, Boundary);
        #endif /* THREEDIM */



  /* 2. DDD data objects (without DDD_HEADER) */

  DDD_TypeDefine(TypeMatrix, &m,
                 EL_GDATA,  ELDEF(m.control),
                 EL_LDATA,  ELDEF(m.next),
                 EL_OBJPTR, ELDEF(m.vect),   TypeVector,
                 EL_LDATA,  ELDEF(m.value),
                 EL_END,    &m+1
                 );


  DDD_TypeDefine(TypeVSegment, &vs,
                 EL_GDATA,  ELDEF(vs.control),

                 /* TODO: thePatch LDATA oder OBJPTR? */
                 EL_LDATA,  ELDEF(vs.thePatch),
                 EL_GDATA,  ELDEF(vs.lambda),
                 EL_LDATA,  ELDEF(vs.next),
                 EL_END,    &vs+1

                 /* TODO: alte version. loeschen?
                                 EL_GDATA,  ELDEF(vs+1),	   add a dummy containing segment id
                                 EL_END,    ((char *)&vs)+(sizeof(VSEGMENT)+sizeof(INT))
                  */
                 );


  DDD_TypeDefine(TypeEdge, &e,
                 /* link 0 data */
                 EL_GDATA,  ELDEF(e.links[0].control),
                 EL_LDATA,  ELDEF(e.links[0].next),
                 EL_OBJPTR, ELDEF(e.links[0].nbnode), TypeNode,

                 /* link 1 data */
                 EL_GDATA,  ELDEF(e.links[1].control),
                 EL_LDATA,  ELDEF(e.links[1].next),
                 EL_OBJPTR, ELDEF(e.links[1].nbnode), TypeNode,

                 EL_OBJPTR, ELDEF(e.midnode),  TypeNode,
                 EL_CONTINUE);

  if (dddctrl.edgeData)
    DDD_TypeDefine(TypeEdge, &e,
                   EL_LDATA,  ELDEF(e.vector), TypeVector,
                   EL_CONTINUE);

  DDD_TypeDefine(TypeEdge, &e, EL_END, &e+1);


  DDD_TypeDefine(TypeElementSide, &es,
                 EL_GDATA,  ELDEF(es.control),
                 EL_LDATA,  ELDEF(es.thePatch),
                 EL_GDATA,  ELDEF(es.lambda),
                 EL_LDATA,  ELDEF(es.pred),
                 EL_LDATA,  ELDEF(es.succ),
                 EL_END,    &es+1

                 /* TODO: alte version. loeschen?
                                 EL_GDATA,  &es+1,	 add a dummy containing segment id
                                 EL_END,    ((char *)&es)+(sizeof(ELEMENTSIDE)+sizeof(INT))
                  */
                 );
}


/****************************************************************************/
/*																			*/
/* Function:  ddd_IfInit                                                                                                        */
/*																			*/
/* Purpose:   defines the communication interfaces needed in ug for             */
/*		      management by DDD                                             */
/*																			*/
/* Input:     void                                                              */
/*																			*/
/* Output:    void                                                                                                              */
/*																			*/
/*																			*/
/****************************************************************************/

#ifdef __TWODIM__
static void ddd_IfInit (void)
{
  DDD_TYPE O[8];
  int A[8];
  int B[8];

  O[0] = TypeTrElem;
  O[1] = TypeTrBElem;
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
/* Function:  InitDDDTypes                                                                                              */
/*																			*/
/* Purpose:   define DDD_TYPEs. this function must be called once before	*/
/*            creation of DDD-objects. It depends on correct and complete	*/
/*            initialization of all ug-generic-elements, therefore it must	*/
/*            must be called after completion of InitElementTypes().		*/
/*            As InitElementTypes() will be called whenever new Multigrids	*/
/*            are created/opened, an execution guard prevents this function */
/*            from multiple execution.										*/
/*																			*/
/* Input:     void															*/
/*																			*/
/* Output:    void															*/
/*																			*/
/****************************************************************************/


void InitDDDTypes (void)
{
  /* prevent from multiple execution */
  if (dddctrl.allTypesDefined)
    return;
  dddctrl.allTypesDefined = TRUE;

  ddd_DefineTypes();


  /* display ddd types */
  IFDEBUG(dddif,1)
  DDD_TypeDisplay(TypeVector);
  DDD_TypeDisplay(TypeIVertex);
  DDD_TypeDisplay(TypeBVertex);
  DDD_TypeDisplay(TypeNode);

        #ifdef __TWODIM__
  DDD_TypeDisplay(TypeTrElem);
  DDD_TypeDisplay(TypeTrBElem);
  DDD_TypeDisplay(TypeQuElem);
  DDD_TypeDisplay(TypeQuBElem);
        #endif

        #ifdef __THREEDIM__
  DDD_TypeDisplay(TypeTeElem);
  DDD_TypeDisplay(TypeTeBElem);
  DDD_TypeDisplay(TypePyElem);
  DDD_TypeDisplay(TypePyBElem);
  DDD_TypeDisplay(TypeHeElem);
  DDD_TypeDisplay(TypeHeBElem);
        #endif

  /* display dependent types */
  DDD_TypeDisplay(TypeMatrix);
  DDD_TypeDisplay(TypeVSegment);
  DDD_TypeDisplay(TypeEdge);
  DDD_TypeDisplay(TypeElementSide);
  ENDDEBUG

  ddd_HandlerInit();
}


/****************************************************************************/
/*																			*/
/* Function:  InitCurrMG													*/
/*																			*/
/* Purpose:   initialize the current multigrid which is handled by DDD		*/
/*																			*/
/* Input:	  MULTIGRID *MG:	the multigrid to handle						*/
/*																			*/
/* Output:	  void															*/
/*																			*/
/*																			*/
/****************************************************************************/

void InitCurrMG (MULTIGRID *MG)
{
  dddctrl.currMG = MG;

  dddctrl.nodeData = TYPE_DEF_IN_MG(dddctrl.currMG,NODEVECTOR);
  dddctrl.edgeData = TYPE_DEF_IN_MG(dddctrl.currMG,EDGEVECTOR);
  dddctrl.elemData = TYPE_DEF_IN_MG(dddctrl.currMG,ELEMVECTOR);
        #ifdef __THREEDIM__
  dddctrl.sideData = TYPE_DEF_IN_MG(dddctrl.currMG,SIDEVECTOR);
        #else
  dddctrl.sideData = FALSE;
        #endif
}



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
/* Input:     int    *argc: pointer to number of arguments                      */
/*		      char ***argv: pointer to list of argument pointers			*/
/*																			*/
/* Output:    int:   error value											*/
/*																			*/
/*																			*/
/****************************************************************************/


int InitParallel (int *argc, char ***argv)
{
  int i;

  DDD_Init(argc, argv);

  /* initialize context */
  /* TODO: malloc() should be replaced by HEAPs or ddd_memmgr */
  dddctrl.context = (int *)malloc(sizeof(int)*procs);

  /* initial context is master processor only */
  for(i=0; i<procs; i++)
    dddctrl.context[i] = 0;
  dddctrl.context[master] = 1;


  /* initialize type mapping arrays */
  for(i=0; i<MAXOBJECTS; i++)
  {
    dddctrl.types[i] = -1;
    dddctrl.dddObj[i] = FALSE;
  }
  for(i=0; i<MAXDDDTYPES; i++)
  {
    dddctrl.ugtypes[i] = -1;
  }


  /* declare DDD_TYPES, definition must be done later */
  ddd_DeclareTypes();
  dddctrl.allTypesDefined = FALSE;



        #ifdef __TWODIM__
  ddd_IfInit();
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
  /* free memory allocated by InitParallel */
  if (dddctrl.context!=0)
    free(dddctrl.context);

  DDD_Exit();
}

#endif /* ModelP */
